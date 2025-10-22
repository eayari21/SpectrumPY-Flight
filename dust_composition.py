"""Interactive dust composition analysis window.

This module provides a publication-quality Qt window for analysing
time-of-flight (TOF) dust spectra.  The tool can be invoked from the
main IDEX quicklook interface and offers the following capabilities:

* Inspect the three TOF gain stages individually or merge them into a
  single, unsaturated spectrum using the instrument's gain hierarchy.
* Interactively define a horizontal baseline and see the plot update in
  real time.
* Add, edit, and visualise exponentially modified Gaussian (EMG) fits
  for individual mass lines.  Fits may be adjusted numerically through
  a table interface and are persisted back to the HDF5 file.
* Display both the mass (bottom axis) and time (top axis) scales with
  editable stretch/shift parameters.
* Provide composition tables, relative abundance estimates, and a
  simple best-guess description of the analysed sample.

The window has been designed to be resilient in the face of partially
populated datasets—missing gain stages, absent analysis groups, or
read-only files are handled gracefully with clear status messages.
"""

from __future__ import annotations

import csv
import html
import math
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np

try:  # pragma: no cover - Qt import guard
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import (
        QAbstractItemView,
        QDialog,
        QDialogButtonBox,
        QDoubleSpinBox,
        QFormLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QSizePolicy,
        QSplitter,
        QStatusBar,
        QTableWidget,
        QTableWidgetItem,
        QVBoxLayout,
        QWidget,
    )
except Exception:  # pragma: no cover - fallback to PyQt6
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QAbstractItemView,
        QDialog,
        QDialogButtonBox,
        QDoubleSpinBox,
        QFormLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QSizePolicy,
        QSplitter,
        QStatusBar,
        QTableWidget,
        QTableWidgetItem,
        QVBoxLayout,
        QWidget,
    )

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

# --- ERFC shim compatible with NumPy 1.x/2.x (SciPy optional) ---
from typing import Union
import numpy as np

try:
    # Preferred on modern stacks
    from scipy.special import erfc as _erfc_impl  # type: ignore
except Exception:
    # Fallback: vectorize Python's math.erfc
    from math import erfc as _math_erfc

    def _erfc_impl(x: Union[np.ndarray, float, int]) -> np.ndarray:
        arr = np.asarray(x, dtype=float)
        if arr.ndim == 0:
            # Preserve scalar behaviour so callers see a 0-D array
            return np.array(_math_erfc(float(arr)))
        vec = np.vectorize(lambda v: _math_erfc(float(v)), otypes=[float])
        return vec(arr)
# --- end shim ---

def _erfc(values: np.ndarray) -> np.ndarray:
    """Return the complementary error function for *values*."""

    return np.asarray(_erfc_impl(values), dtype=float)

# ---------------------------------------------------------------------------
# Constants & helper utilities
# ---------------------------------------------------------------------------

GAIN_HIGH = 1600.0
GAIN_MEDIUM = 40.0
GAIN_LOW = 1.0

GAIN_MAP = {
    "TOF H": GAIN_HIGH,
    "TOF M": GAIN_MEDIUM,
    "TOF L": GAIN_LOW,
}

COMBINED_DATASET = "CombinedSignal"
COMBINED_TIME_DATASET = "CombinedTime"
ANALYSIS_GROUP = "Analysis"
DUST_GROUP = "DustComposition"
MASS_LINES_DATASET = "MassLines"


def _emg_model(time_values: np.ndarray, amplitude: float, mu: float, sigma: float, lam: float) -> np.ndarray:
    """Evaluate an exponentially modified Gaussian (EMG)."""

    arr = np.asarray(time_values, dtype=float)
    if arr.size == 0:
        return np.zeros(0, dtype=float)
    safe_sigma = sigma if abs(sigma) > 1.0e-15 else 1.0e-15
    safe_lambda = lam if abs(lam) > 1.0e-15 else 1.0e-15
    safe_amplitude = float(amplitude)
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        exponent = np.exp((safe_lambda / 2.0) * (2.0 * mu + safe_lambda * safe_sigma**2 - 2.0 * arr))
    argument = (mu + safe_lambda * safe_sigma**2 - arr) / (np.sqrt(2.0) * safe_sigma)
    return (safe_lambda * safe_amplitude / 2.0) * exponent * _erfc(argument)


def _emg_cdf(time_values: np.ndarray, mu: float, sigma: float, lam: float) -> np.ndarray:
    """Return the cumulative distribution of a unit-area EMG."""

    arr = np.asarray(time_values, dtype=float)
    if arr.size == 0:
        return np.zeros(0, dtype=float)
    safe_sigma = sigma if abs(sigma) > 1.0e-15 else 1.0e-15
    safe_lambda = lam if abs(lam) > 1.0e-15 else 1.0e-15
    z = (arr - mu) / safe_sigma
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        term1 = 0.5 * _erfc(-z / np.sqrt(2.0))
        exponent = np.exp(-safe_lambda * (arr - mu) + 0.5 * (safe_lambda * safe_sigma) ** 2)
        term2 = 0.5 * _erfc(-(z - safe_lambda * safe_sigma) / np.sqrt(2.0))
    cdf = term1 - exponent * term2
    return np.clip(cdf, 0.0, 1.0)


def _emg_area(
    amplitude: float,
    mu: float,
    sigma: float,
    lam: float,
    start: Optional[float] = None,
    end: Optional[float] = None,
) -> float:
    """Return the analytic EMG area between ``start`` and ``end``."""

    if not math.isfinite(amplitude):
        return 0.0
    safe_sigma = sigma if abs(sigma) > 1.0e-15 else 1.0e-15
    safe_lambda = lam if abs(lam) > 1.0e-15 else 1.0e-15
    start_cdf = 0.0
    end_cdf = 1.0
    if start is not None and math.isfinite(start):
        start_cdf = float(_emg_cdf(np.array([start], dtype=float), mu, safe_sigma, safe_lambda)[0])
    if end is not None and math.isfinite(end):
        end_cdf = float(_emg_cdf(np.array([end], dtype=float), mu, safe_sigma, safe_lambda)[0])
    area = float(amplitude) * (end_cdf - start_cdf)
    return area


def _estimate_amplitude_from_curve(
    time_axis: np.ndarray, fit_values: np.ndarray, mu: float, sigma: float, lam: float
) -> float:
    """Infer the EMG amplitude from sampled ``fit_values``."""

    if time_axis.size == 0 or fit_values.size == 0:
        return 0.0
    unit_model = _emg_model(time_axis, 1.0, mu, abs(sigma), abs(lam))
    with np.errstate(divide="ignore", invalid="ignore"):
        mask = np.isfinite(unit_model) & np.isfinite(fit_values) & (unit_model > 0)
        if not np.any(mask):
            return 0.0
        ratios = fit_values[mask] / unit_model[mask]
    ratios = ratios[np.isfinite(ratios)]
    if ratios.size == 0:
        return 0.0
    return float(np.median(ratios))


def _contiguous_mask(condition: np.ndarray, min_samples: int) -> np.ndarray:
    """Return a mask retaining only runs longer than ``min_samples``."""

    condition = np.asarray(condition, dtype=bool)
    if condition.size == 0:
        return np.zeros(0, dtype=bool)
    padded = np.concatenate(([False], condition, [False])).astype(int)
    diff = np.diff(padded)
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    mask = np.zeros_like(condition, dtype=bool)
    for start, end in zip(starts, ends):
        if end - start >= max(1, min_samples):
            mask[start:end] = True
    return mask


def detect_saturation(values: np.ndarray, times: np.ndarray) -> np.ndarray:
    """Heuristically detect saturated portions of a waveform."""

    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.zeros(0, dtype=bool)

    magnitude = np.nanmax(np.abs(arr))
    if not np.isfinite(magnitude) or magnitude == 0.0:
        return np.zeros_like(arr, dtype=bool)

    grad = np.abs(np.gradient(arr))
    derivative_threshold = 0.0025 * magnitude
    plateau = grad < derivative_threshold

    amplitude_threshold = np.nanpercentile(np.abs(arr), 99.7)
    high_amp = np.abs(arr) >= amplitude_threshold

    plateau_mask = plateau & high_amp

    if plateau_mask.size < 2:
        return plateau_mask

    times = np.asarray(times, dtype=float)
    if times.size >= 2:
        dt = float(np.nanmedian(np.diff(times)))
    else:
        dt = 0.0

    if not np.isfinite(dt) or dt <= 0.0:
        min_samples = 12
    else:
        min_samples = max(8, int(math.ceil(1.0 / max(dt, 1.0e-6))))

    return _contiguous_mask(plateau_mask, min_samples)


def _load_mass_reference() -> List[Tuple[float, str]]:
    """Return a list of reference masses and their labels."""

    reference: List[Tuple[float, str]] = []
    csv_path = Path(__file__).with_name("mass_comb.csv")
    if not csv_path.exists():
        return reference
    try:
        with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                try:
                    mass = float(row.get("Mass", "nan"))
                except Exception:
                    continue
                label = row.get("Name", "") or f"Mass {mass:.2f}"
                reference.append((mass, label))
    except Exception:
        return []
    return reference


MASS_REFERENCE = _load_mass_reference()


def nearest_mass_name(target: float) -> str:
    """Return the nearest reference species name for ``target`` mass."""

    if not MASS_REFERENCE or not np.isfinite(target):
        return f"m={target:.2f}"
    best_mass, best_name = MASS_REFERENCE[0]
    best_delta = abs(best_mass - target)
    for mass, name in MASS_REFERENCE[1:]:
        delta = abs(mass - target)
        if delta < best_delta:
            best_mass, best_name, best_delta = mass, name, delta
    return best_name


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class MassLineFit:
    """Container describing a single EMG mass-line fit."""

    line_id: int
    label: str
    mu: float
    sigma: float
    lam: float
    amplitude: float
    time_start: float
    time_end: float
    mass_guess: float
    abundance: float
    time_axis: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))
    fit_values: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))
    color: str = "#d62728"

    def parameters(self) -> Tuple[float, float, float, float]:
        return (self.amplitude, self.mu, self.sigma, self.lam)

    def as_row(self) -> Sequence[float | str]:
        return (
            self.label,
            f"{self.mass_guess:.3f}",
            f"{self.mu:.6f}",
            f"{self.sigma:.6f}",
            f"{self.lam:.6f}",
            f"{self.amplitude:.6f}",
            f"{self.abundance * 100.0:.2f}",
        )

    def window_area(self) -> float:
        area = _emg_area(self.amplitude, self.mu, abs(self.sigma), abs(self.lam), self.time_start, self.time_end)
        return max(area, 0.0)

    def total_area(self) -> float:
        area = _emg_area(self.amplitude, self.mu, abs(self.sigma), abs(self.lam))
        return max(area, 0.0)


# ---------------------------------------------------------------------------
# Dialogs
# ---------------------------------------------------------------------------


class ManualMassLineDialog(QDialog):
    """Dialog allowing users to manually create a mass-line fit."""

    def __init__(self, parent: QWidget, defaults: Dict[str, float | str]):
        super().__init__(parent)
        self.setWindowTitle("Add Manual Mass Line")
        self._result: Optional[Dict[str, float | str]] = None

        layout = QVBoxLayout(self)
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)

        self.label_edit = QLineEdit(self)
        self.label_edit.setText(str(defaults.get("label", "Manual line")))
        form.addRow("Label:", self.label_edit)

        self.mass_spin = QDoubleSpinBox(self)
        self.mass_spin.setDecimals(6)
        self.mass_spin.setRange(-1e6, 1e6)
        self.mass_spin.setValue(float(defaults.get("mass", 0.0)))
        form.addRow("Mass (amu):", self.mass_spin)

        self.mu_spin = QDoubleSpinBox(self)
        self.mu_spin.setDecimals(6)
        self.mu_spin.setRange(-1e6, 1e6)
        self.mu_spin.setValue(float(defaults.get("mu", 0.0)))
        form.addRow("μ (µs):", self.mu_spin)

        self.sigma_spin = QDoubleSpinBox(self)
        self.sigma_spin.setDecimals(6)
        self.sigma_spin.setRange(1e-9, 1e6)
        self.sigma_spin.setSingleStep(1e-3)
        self.sigma_spin.setValue(float(max(abs(float(defaults.get("sigma", 0.01))), 1e-9)))
        form.addRow("σ (µs):", self.sigma_spin)

        self.lambda_spin = QDoubleSpinBox(self)
        self.lambda_spin.setDecimals(6)
        self.lambda_spin.setRange(1e-9, 1e6)
        self.lambda_spin.setSingleStep(1e-3)
        self.lambda_spin.setValue(float(max(abs(float(defaults.get("lam", 1.0))), 1e-9)))
        form.addRow("λ (µs⁻¹):", self.lambda_spin)

        self.amplitude_spin = QDoubleSpinBox(self)
        self.amplitude_spin.setDecimals(6)
        self.amplitude_spin.setRange(0.0, 1e12)
        self.amplitude_spin.setSingleStep(1e-3)
        self.amplitude_spin.setValue(float(max(abs(float(defaults.get("amplitude", 1.0))), 1.0e-9)))
        form.addRow("A (DN·µs):", self.amplitude_spin)

        self.start_spin = QDoubleSpinBox(self)
        self.start_spin.setDecimals(6)
        self.start_spin.setRange(-1e6, 1e6)
        self.start_spin.setValue(float(defaults.get("time_start", 0.0)))
        form.addRow("Start time (µs):", self.start_spin)

        self.end_spin = QDoubleSpinBox(self)
        self.end_spin.setDecimals(6)
        self.end_spin.setRange(-1e6, 1e6)
        self.end_spin.setValue(float(defaults.get("time_end", 1.0)))
        form.addRow("End time (µs):", self.end_spin)

        layout.addLayout(form)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=self)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def accept(self) -> None:  # type: ignore[override]
        try:
            self._result = self._collect_values()
        except ValueError as exc:  # pragma: no cover - UI feedback
            QMessageBox.warning(self, "Invalid Input", str(exc))
            return
        super().accept()

    def collected_values(self) -> Optional[Dict[str, float | str]]:
        """Return the values entered by the user, if available."""

        return self._result

    def _collect_values(self) -> Dict[str, float | str]:
        label = self.label_edit.text().strip() or "Manual line"
        mass = float(self.mass_spin.value())
        mu = float(self.mu_spin.value())
        sigma = float(self.sigma_spin.value())
        lam = float(self.lambda_spin.value())
        amplitude = float(self.amplitude_spin.value())
        start = float(self.start_spin.value())
        end = float(self.end_spin.value())
        if sigma <= 0.0:
            raise ValueError("σ must be positive.")
        if lam <= 0.0:
            raise ValueError("λ must be positive.")
        if amplitude <= 0.0:
            raise ValueError("A must be positive.")
        if end <= start:
            raise ValueError("End time must be greater than start time.")
        return {
            "label": label,
            "mass": mass,
            "mu": mu,
            "sigma": sigma,
            "lam": lam,
            "amplitude": amplitude,
            "time_start": start,
            "time_end": end,
        }


class InspectMassLineDialog(QDialog):
    """Focused inspector for visualising and editing a single EMG mass line."""

    _EMG_HTML = (
        "<span style='font-size:16px;'>"
        "<b>f(t)</b> = <b>A</b> · <b>λ</b>/2 · exp\u207b((t − <b>μ</b>)·<b>λ</b> − (<b>λ</b>·<b>σ</b>)² / 2) · erfc\u208b((<b>μ</b> + (<b>λ</b>·<b>σ</b>)² − t) / (\u221a2 · <b>σ</b>))"
        "</span>"
    )

    def __init__(
        self,
        parent: QWidget,
        line: MassLineFit,
        *,
        time_axis: np.ndarray,
        signal: np.ndarray,
        baseline: float,
        source_name: str,
        mass_converter: Callable[[float], float],
    ):
        super().__init__(parent)

        self.setModal(True)
        self.setWindowTitle("Inspect Mass Line")
        self.setMinimumSize(720, 680)

        self._line = replace(line)
        self._time_axis = np.asarray(time_axis, dtype=float).ravel()
        self._signal = np.asarray(signal, dtype=float).ravel()
        self._baseline = float(baseline)
        self._source_name = source_name
        self._mass_converter = mass_converter
        self._result: Optional[Dict[str, float | str]] = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        self.header_label = QLabel("", self)
        self.header_label.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(self.header_label)

        self.formula_label = QLabel(self)
        self.formula_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.formula_label.setTextFormat(Qt.TextFormat.RichText)
        self.formula_label.setWordWrap(True)
        self.formula_label.setText(self._EMG_HTML)
        layout.addWidget(self.formula_label)

        source_label = QLabel(f"Signal source: {html.escape(self._source_name)}", self)
        source_label.setStyleSheet("color: #495057; font-size: 13px;")
        layout.addWidget(source_label)

        figure_container = QWidget(self)
        figure_layout = QVBoxLayout(figure_container)
        figure_layout.setContentsMargins(0, 0, 0, 0)
        figure_layout.setSpacing(6)

        self.figure = Figure(figsize=(6.5, 3.8), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        figure_layout.addWidget(self.canvas)

        layout.addWidget(figure_container, stretch=1)

        parameter_box = QGroupBox("Editable parameters", self)
        form = QFormLayout(parameter_box)
        form.setContentsMargins(12, 12, 12, 12)
        form.setSpacing(8)

        self.label_edit = QLineEdit(parameter_box)
        self.label_edit.setText(self._line.label)
        self.label_edit.textChanged.connect(self._on_label_changed)
        form.addRow("Label:", self.label_edit)

        time_min = float(np.nanmin(self._time_axis)) if self._time_axis.size else -1.0e6
        time_max = float(np.nanmax(self._time_axis)) if self._time_axis.size else 1.0e6

        self.mu_spin = QDoubleSpinBox(parameter_box)
        self.mu_spin.setDecimals(6)
        self.mu_spin.setRange(time_min, time_max)
        self.mu_spin.setValue(self._line.mu)
        self.mu_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow("μ (µs):", self.mu_spin)

        self.sigma_spin = QDoubleSpinBox(parameter_box)
        self.sigma_spin.setDecimals(6)
        self.sigma_spin.setRange(1.0e-9, 1.0e6)
        self.sigma_spin.setValue(max(abs(self._line.sigma), 1.0e-6))
        self.sigma_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow("σ (µs):", self.sigma_spin)

        self.lambda_spin = QDoubleSpinBox(parameter_box)
        self.lambda_spin.setDecimals(6)
        self.lambda_spin.setRange(1.0e-9, 1.0e6)
        self.lambda_spin.setValue(max(abs(self._line.lam), 1.0e-6))
        self.lambda_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow("λ (µs⁻¹):", self.lambda_spin)

        self.amplitude_spin = QDoubleSpinBox(parameter_box)
        self.amplitude_spin.setDecimals(6)
        self.amplitude_spin.setRange(1.0e-12, 1.0e12)
        self.amplitude_spin.setValue(max(abs(self._line.amplitude), 1.0e-9))
        self.amplitude_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow("A (DN·µs):", self.amplitude_spin)

        self.start_spin = QDoubleSpinBox(parameter_box)
        self.start_spin.setDecimals(6)
        self.start_spin.setRange(time_min, time_max)
        self.start_spin.setValue(self._line.time_start)
        self.start_spin.valueChanged.connect(self._on_window_changed)
        form.addRow("Start time (µs):", self.start_spin)

        self.end_spin = QDoubleSpinBox(parameter_box)
        self.end_spin.setDecimals(6)
        self.end_spin.setRange(time_min, time_max)
        self.end_spin.setValue(self._line.time_end)
        self.end_spin.valueChanged.connect(self._on_window_changed)
        form.addRow("End time (µs):", self.end_spin)

        layout.addWidget(parameter_box)

        self.mass_hint_label = QLabel("", self)
        self.mass_hint_label.setStyleSheet("font-size: 14px; font-style: italic; color: #495057;")
        layout.addWidget(self.mass_hint_label)

        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

        self._update_header()
        self._update_mass_hint()
        self._update_plot()

    def _update_header(self) -> None:
        label = self.label_edit.text().strip() or "Mass line"
        self.header_label.setText(f"Inspecting: {html.escape(label)}")

    def _update_mass_hint(self) -> None:
        mu = float(self.mu_spin.value())
        try:
            mass_value = float(self._mass_converter(mu))
        except Exception:
            mass_value = float("nan")
        if math.isfinite(mass_value):
            text = f"Estimated mass from μ: {mass_value:.3f} amu"
        else:
            text = "Estimated mass from μ: unavailable"
        self.mass_hint_label.setText(text)

    def _update_plot(self) -> None:
        amplitude = float(max(self.amplitude_spin.value(), 1.0e-12))
        mu = float(self.mu_spin.value())
        sigma = float(max(self.sigma_spin.value(), 1.0e-9))
        lam = float(max(self.lambda_spin.value(), 1.0e-9))
        start = float(self.start_spin.value())
        end = float(self.end_spin.value())

        if end <= start:
            end = start + 1.0e-6

        if self._time_axis.size == 0 or self._signal.size == 0:
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No waveform data available", ha="center", va="center")
            ax.set_axis_off()
            self.canvas.draw_idle()
            return

        padding = max((end - start) * 0.25, 1.0e-3)
        window_min = start - padding
        window_max = end + padding
        mask = (self._time_axis >= window_min) & (self._time_axis <= window_max)
        time_data = self._time_axis[mask]
        signal_data = self._signal[mask]
        if time_data.size == 0:
            time_data = self._time_axis
            signal_data = self._signal

        fit_time = np.linspace(start, end, 1200)
        fit_values = _emg_model(fit_time, amplitude, mu, sigma, lam)

        offset = 0.0
        finite_signal = signal_data[np.isfinite(signal_data)]
        if finite_signal.size:
            min_signal = float(np.nanmin(finite_signal))
            if min_signal <= 0.0:
                offset = abs(min_signal) + 1.0e-9
        signal_plot = signal_data + offset
        fit_plot = fit_values + offset

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.scatter(time_data, signal_plot, s=22, c="#1f77b4", alpha=0.75, label="Waveform")
        ax.plot(fit_time, fit_plot, color="#d62728", linewidth=2.2, label="EMG fit")
        ax.set_xlabel("Time (µs)")
        ax.set_ylabel("Signal (DN – baseline)")
        ax.set_yscale("symlog", linthresh=1.0e-3)
        ax.set_xlim(window_min, window_max)
        ax.set_title("Zoomed mass line view", fontsize=14, fontweight="bold")
        ax.legend(loc="best")
        ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.35)
        self.canvas.draw_idle()

    def _on_label_changed(self, _text: str) -> None:
        self._update_header()

    def _on_parameter_changed(self, _value: float) -> None:
        self._update_mass_hint()
        self._update_plot()

    def _on_window_changed(self, _value: float) -> None:
        if self.start_spin.value() >= self.end_spin.value():
            if _value == self.start_spin.value():
                self.end_spin.blockSignals(True)
                self.end_spin.setValue(self.start_spin.value() + 1.0e-6)
                self.end_spin.blockSignals(False)
            else:
                self.start_spin.blockSignals(True)
                self.start_spin.setValue(self.end_spin.value() - 1.0e-6)
                self.start_spin.blockSignals(False)
        self._update_plot()

    def accept(self) -> None:  # type: ignore[override]
        try:
            self._result = self._collect_values()
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid parameters", str(exc))
            return
        super().accept()

    def _collect_values(self) -> Dict[str, float | str]:
        label = self.label_edit.text().strip() or "Mass line"
        amplitude = float(max(self.amplitude_spin.value(), 1.0e-12))
        mu = float(self.mu_spin.value())
        sigma = float(max(self.sigma_spin.value(), 1.0e-9))
        lam = float(max(self.lambda_spin.value(), 1.0e-9))
        start = float(self.start_spin.value())
        end = float(self.end_spin.value())
        if end <= start:
            raise ValueError("The end time must be greater than the start time.")
        try:
            mass_guess = float(self._mass_converter(mu))
        except Exception:
            mass_guess = float("nan")
        return {
            "label": label,
            "amplitude": amplitude,
            "mu": mu,
            "sigma": sigma,
            "lam": lam,
            "time_start": start,
            "time_end": end,
            "mass_guess": mass_guess,
        }

    def collected_values(self) -> Optional[Dict[str, float | str]]:
        return self._result

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def _safe_delete(group: h5py.Group, name: str) -> None:
    try:
        if name in group:
            del group[name]
    except Exception:
        pass


def _write_dataset(group: h5py.Group, name: str, data: np.ndarray) -> None:
    if not isinstance(group, h5py.Group):
        return
    try:
        _safe_delete(group, name)
        group.create_dataset(name, data=data)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------


class DustCompositionWindow(QMainWindow):
    """Qt window that orchestrates the dust composition analysis workflow."""

    def __init__(self, h5: h5py.File, event_name: str, parent: Optional[QWidget] = None):
        super().__init__(parent)

        self._h5 = h5
        self._event = event_name
        self._group: Optional[h5py.Group] = None

        if self._h5 is not None:
            grp = self._h5.get(event_name)
            if isinstance(grp, h5py.Group):
                self._group = grp

        self.setWindowTitle(f"Dust Composition Analysis — Event {event_name}")
        self.resize(1320, 880)

        self._time_axis = np.zeros(0)
        self._waveforms: Dict[str, np.ndarray] = {}
        self._combined: Optional[np.ndarray] = None
        self._combined_cached_mass: Optional[np.ndarray] = None
        self._baseline = 0.0
        self._mass_params = {"stretch": 1.0, "shift": 0.0}
        self._mass_params_loaded = False
        self._mass_lines: List[MassLineFit] = []
        self._mass_line_counter = 0
        self._selected_line_id: Optional[int] = None

        self._combined_axis = None
        self._combined_time_axis = None
        self._baseline_artist = None
        self._span_selector: Optional[SpanSelector] = None
        self._in_baseline_mode = False
        self._block_table_signals = False
        self._block_selection_signals = False

        self._load_datasets()
        self._load_saved_state()

        central = QWidget(self)
        layout = QHBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(10)
        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        layout.addWidget(splitter)
        self.setCentralWidget(central)

        figure_container = QWidget(self)
        figure_layout = QVBoxLayout(figure_container)
        figure_layout.setContentsMargins(0, 0, 0, 0)
        figure_layout.setSpacing(0)

        self.figure = Figure(figsize=(8.2, 6.4), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.canvas.mpl_connect("button_press_event", self._on_canvas_click)

        self.toolbar = NavigationToolbar2QT(self.canvas, figure_container)

        figure_layout.addWidget(self.toolbar)
        figure_layout.addWidget(self.canvas)
        splitter.addWidget(figure_container)

        self.control_panel = QWidget(self)
        self.control_layout = QVBoxLayout(self.control_panel)
        self.control_layout.setContentsMargins(0, 0, 0, 0)
        self.control_layout.setSpacing(12)
        splitter.addWidget(self.control_panel)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        self._build_controls()

        bar = QStatusBar(self)
        self.setStatusBar(bar)

        self._refresh_plot(initial=True)
        self._update_tables()
        self._update_summary()

    # ---- Data loading ---------------------------------------------------
    def _load_datasets(self) -> None:
        if not self._group:
            return
        try:
            time_ds = self._group.get("Time (high sampling)")
            if time_ds is not None:
                self._time_axis = np.asarray(time_ds[()], dtype=float).ravel()
        except Exception:
            self._time_axis = np.zeros(0)
        for channel in ("TOF L", "TOF M", "TOF H"):
            try:
                dataset = self._group.get(channel)
                if dataset is not None:
                    self._waveforms[channel] = np.asarray(dataset[()], dtype=float).ravel()
            except Exception:
                self._waveforms[channel] = np.zeros(0)

    def _load_saved_state(self) -> None:
        if not self._group:
            return
        try:
            analysis = self._group.require_group(ANALYSIS_GROUP)
        except Exception:
            return
        dust_group = analysis.get(DUST_GROUP)
        if not isinstance(dust_group, h5py.Group):
            return
        try:
            if COMBINED_TIME_DATASET in dust_group:
                stored_time = np.asarray(dust_group[COMBINED_TIME_DATASET][()], dtype=float)
                if stored_time.size:
                    self._time_axis = stored_time.ravel()
            if COMBINED_DATASET in dust_group:
                combined = np.asarray(dust_group[COMBINED_DATASET][()], dtype=float)
                if combined.size:
                    self._combined = combined.ravel()
        except Exception:
            pass
        try:
            self._baseline = float(dust_group.attrs.get("Baseline", self._baseline))
        except Exception:
            self._baseline = 0.0
        try:
            self._mass_params["stretch"] = float(dust_group.attrs.get("MassStretch", 1.0))
            self._mass_params["shift"] = float(dust_group.attrs.get("MassShift", 0.0))
            self._mass_params_loaded = True
        except Exception:
            self._mass_params = {"stretch": 1.0, "shift": 0.0}
            self._mass_params_loaded = False
        table = None
        if MASS_LINES_DATASET in dust_group:
            try:
                table = dust_group[MASS_LINES_DATASET][()]
            except Exception:
                table = None
        if table is not None:
            self._load_mass_line_table(table)
        if not self._mass_lines:
            self._load_mass_lines_from_channels(analysis)
        try:
            fits_group = dust_group.get("Fits")
        except Exception:
            fits_group = None
        if isinstance(fits_group, h5py.Group):
            for line in self._mass_lines:
                name = f"line_{line.line_id}"
                if name in fits_group:
                    line_group = fits_group[name]
                    try:
                        line.time_axis = np.asarray(line_group["time"][()], dtype=float)
                        line.fit_values = np.asarray(line_group["values"][()], dtype=float)
                        if (not math.isfinite(line.amplitude) or line.amplitude <= 0.0) and line.time_axis.size:
                            inferred = _estimate_amplitude_from_curve(
                                line.time_axis,
                                line.fit_values,
                                line.mu,
                                line.sigma,
                                line.lam,
                            )
                            if math.isfinite(inferred) and inferred > 0.0:
                                line.amplitude = inferred
                            else:
                                area = float(
                                    np.trapz(np.clip(line.fit_values, 0.0, None), line.time_axis)
                                )
                                if math.isfinite(area) and area > 0.0:
                                    line.amplitude = area
                    except Exception:
                        line.time_axis = np.zeros(0)
                        line.fit_values = np.zeros(0)
        if self._mass_lines and self._selected_line_id is None:
            self._selected_line_id = self._mass_lines[0].line_id

    def _load_mass_line_table(self, table: np.ndarray) -> bool:
        dtype = getattr(table, "dtype", None)
        if dtype is None or not getattr(dtype, "names", None):
            return False

        def _get_field(entry: np.void, candidates: Sequence[str], default=None):
            for name in candidates:
                if name in entry.dtype.names:
                    value = entry[name]
                    if isinstance(value, bytes):
                        try:
                            return value.decode("utf-8")
                        except Exception:
                            return value.decode("latin-1", "ignore")
                    return value
            return default

        def _coerce_float(value, default: float = 0.0) -> float:
            try:
                result = float(value)
            except Exception:
                return default
            if math.isfinite(result):
                return result
            return default

        mu_mass_pairs: List[Tuple[float, float]] = []
        loaded_any = False
        entries = np.atleast_1d(table)
        for entry in entries:
            try:
                label_raw = _get_field(entry, ("label", "name", "species"), default="")
                label = str(label_raw) if label_raw not in (None, "", b"") else ""
                line_id = int(_coerce_float(_get_field(entry, ("id", "line_id", "peak_index"), default=self._mass_line_counter)))
                mu = _coerce_float(_get_field(entry, ("mu", "time_mu", "center"), default=0.0))
                sigma = _coerce_float(_get_field(entry, ("sigma", "width"), default=0.01))
                lam = _coerce_float(_get_field(entry, ("lam", "lambda", "tau"), default=1.0), default=1.0)
                time_start = _coerce_float(
                    _get_field(entry, ("time_start", "t_start"), default=mu - 3.0 * abs(sigma)),
                    default=mu - 3.0 * abs(sigma),
                )
                time_end = _coerce_float(
                    _get_field(entry, ("time_end", "t_end"), default=mu + 3.0 * abs(sigma)),
                    default=mu + 3.0 * abs(sigma),
                )
                if not math.isfinite(time_start):
                    time_start = mu - 3.0 * abs(sigma)
                if not math.isfinite(time_end) or time_end <= time_start:
                    time_end = time_start + max(abs(sigma) * 6.0, 1e-6)
                mass_guess = _coerce_float(_get_field(entry, ("mass", "mass_guess"), default=0.0))
                amplitude = _coerce_float(
                    _get_field(entry, ("amplitude", "A", "area", "signal_amplitude"), default=float("nan")),
                    default=float("nan"),
                )
                if not math.isfinite(amplitude):
                    amplitude = 0.0
                abundance = _coerce_float(
                    _get_field(entry, ("abundance", "relative_abundance"), default=0.0),
                    default=0.0,
                )
                if not label:
                    label = nearest_mass_name(mass_guess)
                line = MassLineFit(
                    line_id=line_id,
                    label=label,
                    mu=mu,
                    sigma=sigma,
                    lam=lam,
                    amplitude=amplitude,
                    time_start=time_start,
                    time_end=time_end,
                    mass_guess=mass_guess,
                    abundance=abundance,
                )
                self._mass_lines.append(line)
                self._mass_line_counter = max(self._mass_line_counter, line.line_id + 1)
                if math.isfinite(mu) and math.isfinite(mass_guess):
                    mu_mass_pairs.append((mu, mass_guess))
                loaded_any = True
            except Exception:
                continue

        if loaded_any and not self._mass_params_loaded:
            self._estimate_mass_axis(mu_mass_pairs)
        return loaded_any

    def _load_mass_lines_from_channels(self, analysis: h5py.Group) -> None:
        if not isinstance(analysis, h5py.Group):
            return
        for channel in ("TOF H", "TOF M", "TOF L"):
            try:
                channel_group = analysis.get(channel)
            except Exception:
                channel_group = None
            if not isinstance(channel_group, h5py.Group):
                continue
            dataset = channel_group.get(MASS_LINES_DATASET)
            if dataset is None:
                continue
            try:
                table = dataset[()]
            except Exception:
                continue
            if self._load_mass_line_table(table):
                break

    def _estimate_mass_axis(self, pairs: Iterable[Tuple[float, float]]) -> None:
        mu_values: List[float] = []
        mass_values: List[float] = []
        for mu, mass in pairs:
            if math.isfinite(mu) and math.isfinite(mass):
                mu_values.append(mu)
                mass_values.append(mass)
        if len(mu_values) < 2:
            return
        mu_arr = np.asarray(mu_values, dtype=float)
        mass_arr = np.asarray(mass_values, dtype=float)
        try:
            A = np.vstack([mu_arr, np.ones_like(mu_arr)]).T
            result, *_ = np.linalg.lstsq(A, mass_arr, rcond=None)
            stretch = float(result[0])
            intercept = float(result[1])
        except Exception:
            return
        if not math.isfinite(stretch) or abs(stretch) < 1.0e-12:
            return
        shift = -intercept / stretch
        self._mass_params["stretch"] = stretch
        self._mass_params["shift"] = shift
        self._mass_params_loaded = True
        self._combined_cached_mass = None

    # ---- UI construction ------------------------------------------------
    def _build_controls(self) -> None:
        self._build_action_buttons()
        self._build_baseline_controls()
        self._build_mass_axis_controls()
        self._build_mass_line_table()
        self._build_summary_section()
        self.control_layout.addStretch(1)

    def _build_action_buttons(self) -> None:
        box = QGroupBox("Waveform Modes", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(8)
        self.combine_button = QPushButton("Combine TOF", box)
        self.combine_button.setCheckable(True)
        self.combine_button.setToolTip(
            "Merge the three TOF gain stages into a single spectrum using the instrument gain ratios."
        )
        self.combine_button.clicked.connect(self._toggle_combine)
        layout.addWidget(self.combine_button)
        self.reset_view_button = QPushButton("Reset View", box)
        self.reset_view_button.setToolTip("Return to the individual gain-stage plots.")
        self.reset_view_button.clicked.connect(self._reset_view)
        layout.addWidget(self.reset_view_button)
        self.save_button = QPushButton("Save Analysis", box)
        self.save_button.setToolTip("Persist the current dust composition analysis back into the HDF5 file.")
        self.save_button.clicked.connect(self._save_to_file)
        layout.addWidget(self.save_button)
        self.control_layout.addWidget(box)

    def _build_baseline_controls(self) -> None:
        box = QGroupBox("Baseline", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(6)
        self.baseline_button = QPushButton("Select Baseline", box)
        self.baseline_button.setCheckable(True)
        self.baseline_button.setToolTip("Click to enter baseline selection mode, then click on the plot to place a horizontal baseline.")
        self.baseline_button.toggled.connect(self._toggle_baseline_mode)
        layout.addWidget(self.baseline_button)
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        self.baseline_spin = QDoubleSpinBox(box)
        self.baseline_spin.setDecimals(6)
        self.baseline_spin.setRange(-1e6, 1e6)
        self.baseline_spin.setValue(self._baseline)
        self.baseline_spin.valueChanged.connect(self._on_baseline_spin_changed)
        form.addRow("Baseline value:", self.baseline_spin)
        layout.addLayout(form)
        self.control_layout.addWidget(box)

    def _build_mass_axis_controls(self) -> None:
        box = QGroupBox("Mass Axis", self.control_panel)
        layout = QFormLayout(box)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)
        self.mass_stretch_spin = QDoubleSpinBox(box)
        self.mass_stretch_spin.setDecimals(6)
        self.mass_stretch_spin.setRange(1e-6, 1e6)
        self.mass_stretch_spin.setValue(self._mass_params["stretch"])
        self.mass_stretch_spin.valueChanged.connect(self._on_mass_params_changed)
        layout.addRow("Stretch:", self.mass_stretch_spin)
        self.mass_shift_spin = QDoubleSpinBox(box)
        self.mass_shift_spin.setDecimals(6)
        self.mass_shift_spin.setRange(-1e6, 1e6)
        self.mass_shift_spin.setValue(self._mass_params["shift"])
        self.mass_shift_spin.valueChanged.connect(self._on_mass_params_changed)
        layout.addRow("Shift:", self.mass_shift_spin)
        self.auto_mass_button = QPushButton("Auto-calc axis", box)
        self.auto_mass_button.setToolTip(
            "Estimate the mass-axis stretch and shift using existing mass lines."
        )
        self.auto_mass_button.clicked.connect(self._auto_calculate_mass_axis)
        layout.addRow("", self.auto_mass_button)
        self.add_mass_button = QPushButton("Add Mass Line", box)
        self.add_mass_button.setToolTip("Select a region on the combined plot to fit an EMG mass line.")
        self.add_mass_button.setCheckable(True)
        self.add_mass_button.toggled.connect(self._toggle_mass_line_mode)
        layout.addRow("", self.add_mass_button)
        self.control_layout.addWidget(box)

    def _build_mass_line_table(self) -> None:
        box = QGroupBox("Mass Line Fits", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(6)
        self.mass_table = QTableWidget(box)
        self.mass_table.setColumnCount(7)
        self.mass_table.setHorizontalHeaderLabels([
            "Label",
            "Mass (amu)",
            "μ (µs)",
            "σ (µs)",
            "λ (µs⁻¹)",
            "A (DN·µs)",
            "Abundance (%)",
        ])
        header = self.mass_table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setDefaultSectionSize(130)
        self.mass_table.verticalHeader().setVisible(False)
        self.mass_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.mass_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.mass_table.itemChanged.connect(self._on_mass_table_changed)
        self.mass_table.itemSelectionChanged.connect(self._on_mass_table_selection_changed)
        layout.addWidget(self.mass_table)
        self.inspect_mass_button = QPushButton("Inspect Selected", box)
        self.inspect_mass_button.setToolTip(
            "Open a zoomed-in view of the selected mass line with editable EMG parameters."
        )
        self.inspect_mass_button.setEnabled(False)
        self.inspect_mass_button.clicked.connect(self._inspect_selected_mass_line)
        layout.addWidget(self.inspect_mass_button)
        self.manual_mass_button = QPushButton("Add Manual Line", box)
        self.manual_mass_button.setToolTip("Enter EMG parameters directly without selecting a region on the plot.")
        self.manual_mass_button.clicked.connect(self._prompt_manual_mass_line)
        layout.addWidget(self.manual_mass_button)
        self.remove_mass_button = QPushButton("Remove Selected", box)
        self.remove_mass_button.clicked.connect(self._remove_selected_mass_line)
        self.remove_mass_button.setEnabled(False)
        layout.addWidget(self.remove_mass_button)
        self.control_layout.addWidget(box)

    def _build_summary_section(self) -> None:
        box = QGroupBox("Composition Summary", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(6)
        self.summary_table = QTableWidget(box)
        self.summary_table.setColumnCount(3)
        self.summary_table.setHorizontalHeaderLabels(["Label", "Mass (amu)", "Relative (%)"])
        self.summary_table.verticalHeader().setVisible(False)
        self.summary_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.summary_table.horizontalHeader()
        header.setStretchLastSection(True)
        layout.addWidget(self.summary_table)
        self.sample_guess_label = QLabel("Sample guess: –", box)
        self.sample_guess_label.setWordWrap(True)
        layout.addWidget(self.sample_guess_label)
        self.control_layout.addWidget(box)

    # ---- Event handlers -------------------------------------------------
    def _toggle_combine(self, checked: bool) -> None:
        if checked:
            if self._combined is None:
                self._combined = self._combine_waveforms()
            self._refresh_plot()
        else:
            self._refresh_plot(show_combined=False)

    def _reset_view(self) -> None:
        self.combine_button.setChecked(False)
        self._refresh_plot(show_combined=False)

    def _toggle_baseline_mode(self, checked: bool) -> None:
        self._in_baseline_mode = checked
        if checked:
            self.statusBar().showMessage("Baseline mode active — click on the plot to set the baseline.", 8000)
        else:
            self.statusBar().clearMessage()

    def _toggle_mass_line_mode(self, checked: bool) -> None:
        if checked and not self.combine_button.isChecked():
            self.combine_button.setChecked(True)
            if self._combined is None:
                self._combined = self._combine_waveforms()
            self._refresh_plot()
        self._set_span_selector_active(checked)
        if checked and self._span_selector is not None:
            self.statusBar().showMessage("Select a region to fit an EMG mass line.", 6000)
        else:
            self.statusBar().clearMessage()

    def _prompt_manual_mass_line(self) -> None:
        defaults = self._manual_mass_defaults()
        dialog = ManualMassLineDialog(self, defaults)
        result_code = dialog.exec()
        try:
            accepted = result_code == QDialog.DialogCode.Accepted  # type: ignore[attr-defined]
        except Exception:
            accepted = int(result_code) == int(QDialog.DialogCode.Accepted)
        if not accepted:
            return
        data = dialog.collected_values()
        if not data:
            return
        line = MassLineFit(
            line_id=self._mass_line_counter,
            label=str(data.get("label", "Manual line")),
            mu=float(data.get("mu", 0.0)),
            sigma=float(abs(float(data.get("sigma", 0.01)))),
            lam=float(abs(float(data.get("lam", 1.0)))),
            amplitude=float(max(abs(float(data.get("amplitude", 1.0))), 1.0e-9)),
            time_start=float(data.get("time_start", 0.0)),
            time_end=float(data.get("time_end", 1.0)),
            mass_guess=float(data.get("mass", 0.0)),
            abundance=0.0,
            color="#2ca02c",
        )
        dense_time = np.linspace(line.time_start, line.time_end, 800)
        line.time_axis = dense_time
        line.fit_values = _emg_model(dense_time, line.amplitude, line.mu, abs(line.sigma), abs(line.lam))
        self._mass_line_counter += 1
        self._mass_lines.append(line)
        self._selected_line_id = line.line_id
        self._update_mass_line_abundances()
        self._update_tables()
        self._update_summary()
        self._refresh_plot()

    def _inspect_selected_mass_line(self) -> None:
        line = self._current_mass_line()
        if line is None:
            QMessageBox.information(self, "No Selection", "Select a mass line to inspect.")
            return
        time_axis, signal, source = self._inspection_waveform()
        if time_axis.size == 0 or signal.size == 0:
            QMessageBox.warning(self, "No Data", "No waveform data are available for inspection.")
            return

        def _mass_from_time(value: float) -> float:
            arr = np.array([value], dtype=float)
            converted = self._time_to_mass(arr)
            if converted.size:
                return float(converted[0])
            return float("nan")

        dialog = InspectMassLineDialog(
            self,
            line,
            time_axis=time_axis,
            signal=signal,
            baseline=self._baseline,
            source_name=source,
            mass_converter=_mass_from_time,
        )
        result_code = dialog.exec()
        try:
            accepted = result_code == QDialog.DialogCode.Accepted  # type: ignore[attr-defined]
        except Exception:
            accepted = int(result_code) == int(QDialog.DialogCode.Accepted)
        if not accepted:
            return
        values = dialog.collected_values()
        if not values:
            return

        line.label = str(values.get("label", line.label))
        line.amplitude = float(values.get("amplitude", line.amplitude))
        line.mu = float(values.get("mu", line.mu))
        line.sigma = float(values.get("sigma", line.sigma))
        line.lam = float(values.get("lam", line.lam))
        line.time_start = float(values.get("time_start", line.time_start))
        line.time_end = float(values.get("time_end", line.time_end))
        mass_guess = float(values.get("mass_guess", line.mass_guess))
        if math.isfinite(mass_guess):
            line.mass_guess = mass_guess
        self._recompute_mass_line(line, preserve_label=True)
        self._selected_line_id = line.line_id
        self._update_tables()
        self._update_summary()
        self._refresh_plot()

    def _set_span_selector_active(self, active: bool) -> None:
        if self._span_selector is not None:
            self._span_selector.disconnect_events()
            self._span_selector = None
        if active and self._combined_axis is not None:
            self._span_selector = SpanSelector(
                self._combined_axis,
                self._on_mass_region_selected,
                direction="horizontal",
                useblit=True,
                props=dict(alpha=0.3, facecolor="#ffdd55"),
            )

    def _on_canvas_click(self, event) -> None:
        if not self._in_baseline_mode:
            return
        if event is None or event.ydata is None:
            return
        if self.combine_button.isChecked() and self._combined_axis:
            valid_axes = {self._combined_axis}
            if self._combined_time_axis is not None:
                valid_axes.add(self._combined_time_axis)
            if event.inaxes not in valid_axes:
                return
        self._set_baseline(float(event.ydata), from_user=True)

    def _on_mass_region_selected(self, xmin: float, xmax: float) -> None:
        if not self.combine_button.isChecked():
            return
        if self._combined is None:
            self._combined = self._combine_waveforms()
        if self._combined is None or self._combined.size == 0:
            QMessageBox.information(self, "No Data", "Combined TOF waveform is empty; cannot fit EMG.")
            return
        if abs(xmax - xmin) <= 0:
            return
        self.add_mass_line(xmin, xmax)

    def _on_baseline_spin_changed(self, value: float) -> None:
        self._set_baseline(float(value))

    def _auto_calculate_mass_axis(self) -> None:
        pairs: List[Tuple[float, float]] = []
        for line in self._mass_lines:
            mu = float(line.mu)
            mass = float(line.mass_guess)
            if math.isfinite(mu) and math.isfinite(mass):
                pairs.append((mu, mass))
        if len(pairs) < 2:
            QMessageBox.information(
                self,
                "Insufficient Mass Lines",
                "Add at least two mass lines with valid μ and mass values to estimate the axis.",
            )
            return
        previous = dict(self._mass_params)
        self._estimate_mass_axis(pairs)
        stretch = float(self._mass_params.get("stretch", 1.0))
        shift = float(self._mass_params.get("shift", 0.0))
        if hasattr(self, "mass_stretch_spin") and hasattr(self, "mass_shift_spin"):
            self.mass_stretch_spin.blockSignals(True)
            self.mass_shift_spin.blockSignals(True)
            self.mass_stretch_spin.setValue(stretch)
            self.mass_shift_spin.setValue(shift)
            self.mass_stretch_spin.blockSignals(False)
            self.mass_shift_spin.blockSignals(False)
        # If the estimate failed the parameters will remain unchanged; force a refresh anyway.
        if (
            not math.isclose(previous.get("stretch", 1.0), stretch, rel_tol=1e-9, abs_tol=1e-9)
            or not math.isclose(previous.get("shift", 0.0), shift, rel_tol=1e-9, abs_tol=1e-9)
        ):
            self._on_mass_params_changed()
        else:
            # No change in parameters, but refresh the plot to ensure the UI stays in sync.
            self._refresh_plot()
            self._update_tables()
            self._update_summary()

    def _on_mass_params_changed(self) -> None:
        self._mass_params["stretch"] = float(self.mass_stretch_spin.value())
        self._mass_params["shift"] = float(self.mass_shift_spin.value())
        self._combined_cached_mass = None
        self._refresh_plot()
        self._update_tables()
        self._update_summary()

    def _on_mass_table_changed(self, item: QTableWidgetItem) -> None:
        if self._block_table_signals or item is None:
            return
        row = item.row()
        column = item.column()
        if not (0 <= row < len(self._mass_lines)):
            return
        line = self._mass_lines[row]
        text = item.text().strip()
        try:
            if column == 0:
                line.label = text or line.label
            elif column == 1:
                line.mass_guess = float(text)
            elif column == 2:
                line.mu = float(text)
            elif column == 3:
                line.sigma = float(text)
            elif column == 4:
                line.lam = float(text)
            elif column == 5:
                line.amplitude = float(text)
        except ValueError:
            self._block_table_signals = True
            for col, value in enumerate(line.as_row()):
                self.mass_table.item(row, col).setText(str(value))
            self._block_table_signals = False
            return
        if column in (2, 3, 4, 5):
            self._recompute_mass_line(line)
        if column == 1:
            line.label = nearest_mass_name(line.mass_guess)
        self._update_tables()
        self._update_summary()
        self._refresh_plot()

    # ---- Plotting -------------------------------------------------------
    def _refresh_plot(self, show_combined: Optional[bool] = None, initial: bool = False) -> None:
        show_combined = self.combine_button.isChecked() if show_combined is None else show_combined
        self.figure.clear()
        self._combined_axis = None
        self._combined_time_axis = None
        self._baseline_artist = None
        if show_combined:
            if self._combined is None:
                self._combined = self._combine_waveforms()
            ax = self.figure.add_subplot(111)
            self._combined_axis = ax
            self._plot_combined(ax)
        else:
            self._plot_individual_axes()
        self.canvas.draw_idle()
        # Recreate the interactive span selector after the axes are rebuilt so
        # that click-and-drag mass selections continue to function.
        if hasattr(self, "add_mass_button"):
            self._set_span_selector_active(self.add_mass_button.isChecked())
        else:
            self._set_span_selector_active(False)
        if initial and self._combined is not None:
            self.combine_button.setChecked(True)
            self._refresh_plot(show_combined=True)

    def _plot_individual_axes(self) -> None:
        available = [k for k, v in self._waveforms.items() if v.size]
        if not available:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No TOF waveforms available.", ha="center", va="center", transform=ax.transAxes, fontsize=16)
            ax.axis("off")
            return
        grid = self.figure.add_gridspec(len(available), 1, hspace=0.35)
        for idx, channel in enumerate(available):
            ax = self.figure.add_subplot(grid[idx, 0])
            data = self._waveforms[channel]
            n = min(data.size, self._time_axis.size)
            x = self._time_axis[:n]
            y = data[:n]
            ax.plot(x, y, label=channel)
            ax.set_ylabel("Signal [pC/Δt]", fontsize=12)
            ax.set_title(channel, fontsize=14, fontweight="bold")
            ax.grid(True, alpha=0.3)
            if idx == len(available) - 1:
                ax.set_xlabel("Time [µs]", fontsize=12)
            else:
                ax.set_xticklabels([])

    def _plot_combined(self, ax) -> None:
        if self._combined is None or self._combined.size == 0:
            ax.text(0.5, 0.5, "Combined waveform unavailable.", ha="center", va="center", transform=ax.transAxes, fontsize=16)
            ax.axis("off")
            return
        mass_axis = self._combined_mass_axis()
        n = min(self._combined.size, mass_axis.size)
        ax.set_zorder(2)
        corrected = self._combined[:n] - self._baseline
        ax.plot(mass_axis[:n], corrected, color="#1f77b4", linewidth=1.6, label="Combined TOF")
        ax.set_facecolor("#f9fbff")
        ax.grid(True, alpha=0.35)
        ax.set_xlabel("Mass [amu]", fontsize=14)
        ax.set_ylabel("Signal [pC/Δt]", fontsize=14)
        ax.set_title("Combined TOF Spectrum", fontsize=16, fontweight="bold")
        if self._time_axis.size:
            ax_time = ax.twiny()
            ax_time.set_xlim(self._time_axis.min(), self._time_axis.max())
            ax_time.set_xlabel("Time [µs]", fontsize=13)
            ax_time.set_zorder(1)
            # Ensure the auxiliary time axis does not intercept mouse events that
            # are intended for the mass axis.  When the twinned axis remains
            # interactive Matplotlib reports the mouse events against it instead
            # of ``ax``, preventing the ``SpanSelector`` from ever firing.  By
            # disabling navigation and hiding the background patch we allow the
            # events to fall through to the combined axis where the selector is
            # attached.
            ax_time.set_navigate(False)
            ax_time.patch.set_visible(False)
            self._combined_time_axis = ax_time
        if (self._baseline or self.baseline_spin.value() != 0.0) and np.isfinite(self._baseline):
            self._baseline_artist = ax.axhline(0.0, color="#aa3377", linestyle="--", linewidth=1.2, label="Baseline")
        selected_id = self._selected_line_id
        plotted_any = False
        for line in self._mass_lines:
            if not (line.time_axis.size and line.fit_values.size):
                continue
            try:
                mass_axis = self._time_to_mass(line.time_axis)
            except Exception:
                mass_axis = line.time_axis
            y_values = line.fit_values
            try:
                rgba = to_rgba(line.color)
                base_color = line.color
            except Exception:
                rgba = to_rgba("#d62728")
                base_color = "#d62728"
            if line.line_id == selected_id:
                ax.plot(
                    mass_axis,
                    y_values,
                    linestyle="-",
                    linewidth=2.6,
                    color=base_color,
                    zorder=5,
                    label=line.label,
                )
            else:
                faded = (rgba[0], rgba[1], rgba[2], 0.45)
                ax.plot(
                    mass_axis,
                    y_values,
                    linestyle="-",
                    linewidth=1.2,
                    color=faded,
                    zorder=3,
                    label=line.label,
                )
            finite_values = np.asarray(y_values, dtype=float)
            finite_mask = np.isfinite(finite_values)
            if not np.any(finite_mask):
                plotted_any = True
                continue
            finite_mass = np.asarray(mass_axis, dtype=float)
            safe_values = np.where(finite_mask, finite_values, -np.inf)
            peak_idx = int(np.argmax(safe_values))
            if not (0 <= peak_idx < finite_mass.size) or not np.isfinite(safe_values[peak_idx]):
                plotted_any = True
                continue
            peak_x = float(finite_mass[peak_idx])
            peak_y = float(finite_values[peak_idx])
            label_mass = peak_x if np.isfinite(peak_x) else line.mass_guess
            if np.isfinite(label_mass):
                label_mass_text = f"{line.label}@{label_mass:.1f} amu"
            else:
                label_mass_text = line.label
            try:
                abundance = float(line.abundance)
            except Exception:
                abundance = 0.0
            abundance_text = ""
            if np.isfinite(abundance):
                abundance_text = f" ({max(abundance, 0.0) * 100.0:.1f}%)"
            annotation_text = f"{label_mass_text}{abundance_text}" if label_mass_text else line.label
            if not annotation_text:
                plotted_any = True
                continue
            annotation_kwargs = dict(
                xy=(peak_x, peak_y),
                xytext=(0, 14),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=12 if line.line_id == selected_id else 11,
                fontweight="bold" if line.line_id == selected_id else "medium",
                color="#111111",
                bbox=dict(
                    boxstyle="round,pad=0.35",
                    fc="white",
                    ec=base_color if line.line_id == selected_id else rgba,
                    lw=1.0,
                    alpha=0.9 if line.line_id == selected_id else 0.75,
                ),
                zorder=6 if line.line_id == selected_id else 4,
            )
            ax.annotate(annotation_text, **annotation_kwargs)
            plotted_any = True
        if plotted_any and len(ax.lines) > 1:
            legend = ax.legend(loc="best", fontsize=10)
            if legend is not None:
                legend.set_zorder(7)

    # ---- Combination logic ---------------------------------------------
    def _combine_waveforms(self) -> Optional[np.ndarray]:
        if self._time_axis.size == 0:
            return None
        high = self._waveforms.get("TOF H")
        medium = self._waveforms.get("TOF M")
        low = self._waveforms.get("TOF L")
        if high is None and medium is None and low is None:
            return None
        lengths = [self._time_axis.size]
        for arr in (high, medium, low):
            if arr is not None and arr.size:
                lengths.append(arr.size)
        length = min(lengths) if lengths else 0
        if length <= 0:
            return None
        combined = np.zeros(length, dtype=float)
        if high is not None and high.size:
            combined[:] = high[:length]
            saturated_high = detect_saturation(high[:length], self._time_axis[:length])
        else:
            saturated_high = np.ones(length, dtype=bool)
        high_gain = GAIN_MAP.get("TOF H", 1.0)
        medium_gain = GAIN_MAP.get("TOF M", 1.0)
        low_gain = GAIN_MAP.get("TOF L", 1.0)
        if medium is not None and medium.size:
            medium_scaled = medium[:length] * (high_gain / medium_gain)
            medium_mask = detect_saturation(medium[:length], self._time_axis[:length])
            replace_mask = saturated_high & np.isfinite(medium_scaled)
            combined[replace_mask] = medium_scaled[replace_mask]
        else:
            medium_mask = np.ones(length, dtype=bool)
        if low is not None and low.size:
            low_scaled = low[:length] * (high_gain / low_gain)
            low_mask = detect_saturation(low[:length], self._time_axis[:length])
            replace_mask = saturated_high & medium_mask & np.isfinite(low_scaled)
            combined[replace_mask] = low_scaled[replace_mask]
        return combined

    # ---- Baseline + mass conversions -----------------------------------
    def _set_baseline(self, value: float, from_user: bool = False) -> None:
        self._baseline = float(value)
        if hasattr(self, "baseline_spin"):
            self.baseline_spin.blockSignals(True)
            self.baseline_spin.setValue(self._baseline)
            self.baseline_spin.blockSignals(False)
        self._update_mass_line_abundances()
        self._update_tables()
        self._update_summary()
        if from_user and hasattr(self, "baseline_button"):
            self.baseline_button.setChecked(False)
        self._refresh_plot()

    def _combined_mass_axis(self) -> np.ndarray:
        if self._combined_cached_mass is not None:
            return self._combined_cached_mass
        mass = self._time_to_mass(self._time_axis)
        self._combined_cached_mass = mass
        return mass

    def _time_to_mass(self, time_values: np.ndarray) -> np.ndarray:
        stretch = self._mass_params.get("stretch", 1.0)
        shift = self._mass_params.get("shift", 0.0)
        return stretch * (np.asarray(time_values, dtype=float) - shift)

    def _mass_to_time(self, mass_values: np.ndarray) -> np.ndarray:
        stretch = self._mass_params.get("stretch", 1.0)
        shift = self._mass_params.get("shift", 0.0)
        if stretch == 0:
            stretch = 1.0
        return np.asarray(mass_values, dtype=float) / stretch + shift
    # ---- Mass line management ------------------------------------------
    def _manual_mass_defaults(self) -> Dict[str, float | str]:
        if self._combined is not None and self._combined.size and self._time_axis.size:
            try:
                corrected = self._combined - self._baseline
                idx = int(np.nanargmax(np.clip(corrected, 0.0, None)))
                idx = int(np.clip(idx, 0, self._time_axis.size - 1))
                mu = float(self._time_axis[idx])
            except Exception:
                mu = float(self._time_axis[self._time_axis.size // 2])
        elif self._time_axis.size:
            mu = float(self._time_axis[self._time_axis.size // 2])
        else:
            mu = 0.0
        if self._time_axis.size >= 2:
            try:
                dt = float(np.nanmedian(np.diff(self._time_axis)))
            except Exception:
                dt = 0.0
        else:
            dt = 0.0
        sigma = max(abs(float(dt * 6.0)), 1.0e-3) if dt > 0 else 0.01
        lam = 1.0 / max(sigma, 1.0e-6)
        span = max(sigma * 8.0, 1.0e-3)
        start = mu - span / 2.0
        end = mu + span / 2.0
        mass_guess = float(self._time_to_mass(mu))
        label = nearest_mass_name(mass_guess)
        amplitude = 1.0
        if self._combined is not None and self._combined.size and self._time_axis.size:
            try:
                mask = (self._time_axis >= start) & (self._time_axis <= end)
                if np.any(mask):
                    segment = np.clip(self._combined[mask] - self._baseline, 0.0, None)
                    amplitude = float(np.trapz(segment, self._time_axis[mask]))
            except Exception:
                amplitude = 1.0
        amplitude = max(amplitude, 1.0e-9)
        return {
            "label": label,
            "mass": mass_guess,
            "mu": mu,
            "sigma": sigma,
            "lam": lam,
            "amplitude": amplitude,
            "time_start": start,
            "time_end": end,
        }

    def _current_mass_line(self) -> Optional[MassLineFit]:
        if self._selected_line_id is None:
            return None
        for line in self._mass_lines:
            if line.line_id == self._selected_line_id:
                return line
        return None

    def _inspection_waveform(self) -> Tuple[np.ndarray, np.ndarray, str]:
        label = "Waveform"
        signal: Optional[np.ndarray] = None
        time_axis = np.asarray(self._time_axis, dtype=float)
        if self._combined is None or not getattr(self._combined, "size", 0):
            combined = self._combine_waveforms()
            if combined is not None and combined.size:
                self._combined = combined
        if self._combined is not None and getattr(self._combined, "size", 0):
            signal = np.asarray(self._combined, dtype=float)
            label = "Combined TOF"
        else:
            for key in ("TOF LG", "TOF L", "TOF Low", "TOF M", "TOF H"):
                candidate = self._waveforms.get(key)
                if candidate is not None and getattr(candidate, "size", 0):
                    signal = np.asarray(candidate, dtype=float)
                    label = key
                    break
        if signal is None or time_axis.size == 0:
            return np.zeros(0, dtype=float), np.zeros(0, dtype=float), label
        length = min(time_axis.size, signal.size)
        trimmed_time = time_axis[:length]
        trimmed_signal = signal[:length] - self._baseline
        return trimmed_time, trimmed_signal, label

    def _update_inspect_button_state(self) -> None:
        has_selection = self._current_mass_line() is not None
        if hasattr(self, "inspect_mass_button"):
            self.inspect_mass_button.setEnabled(has_selection)
        if hasattr(self, "remove_mass_button"):
            self.remove_mass_button.setEnabled(has_selection)

    def add_mass_line(self, x_min: float, x_max: float) -> None:
        if self._combined is None:
            return
        mass_min, mass_max = sorted((x_min, x_max))
        time_min, time_max = self._mass_to_time(np.array([mass_min, mass_max]))
        if time_max - time_min <= 0:
            QMessageBox.warning(self, "Invalid Selection", "Please select a wider region to fit.")
            return
        mask = (self._time_axis >= time_min) & (self._time_axis <= time_max)
        if np.count_nonzero(mask) < 6:
            QMessageBox.warning(self, "Selection Too Small", "Select a region containing more samples.")
            return
        x = self._time_axis[mask]
        y = self._combined[mask] - self._baseline
        if not np.any(np.isfinite(y)):
            QMessageBox.warning(self, "No Signal", "The selected region does not contain valid data.")
            return
        try:
            mu_guess = float(x[np.nanargmax(y)])
        except Exception:
            mu_guess = float(np.nanmean(x))
        sigma_guess = max(float(np.nanstd(x)), 1.0e-6)
        lam_guess = max(1.0 / max(sigma_guess, 1.0e-6), 1.0e-6)
        amplitude_guess = float(np.trapz(np.clip(y, 0.0, None), x))
        if not math.isfinite(amplitude_guess) or amplitude_guess <= 0.0:
            amplitude_guess = max(float(np.nanmax(np.clip(y, 0.0, None))) * max(time_max - time_min, 1.0e-6), 1.0e-6)
        try:
            from scipy.optimize import curve_fit  # type: ignore
        except Exception as exc:
            QMessageBox.warning(self, "Missing SciPy", f"SciPy is required for EMG fitting:\n{exc}")
            return

        try:
            params, _ = curve_fit(
                lambda t, amplitude, mu, sigma, lam: _emg_model(t, amplitude, mu, abs(sigma), abs(lam)),
                x,
                y,
                p0=(amplitude_guess, mu_guess, sigma_guess, lam_guess),
                maxfev=20000,
            )
            amplitude_fit, mu_fit, sigma_fit, lam_fit = params
            amplitude_fit = abs(float(amplitude_fit))
            sigma_fit = abs(float(sigma_fit))
            lam_fit = abs(float(lam_fit))
        except Exception:
            QMessageBox.warning(self, "Fit Failed", "Unable to fit an EMG curve to the selected region.")
            return
        dense_time = np.linspace(time_min, time_max, 800)
        fit_curve = _emg_model(dense_time, amplitude_fit, mu_fit, sigma_fit, lam_fit)
        mass_guess = float(self._time_to_mass(mu_fit))
        label = nearest_mass_name(mass_guess)
        line = MassLineFit(
            line_id=self._mass_line_counter,
            label=label,
            mu=float(mu_fit),
            sigma=float(sigma_fit),
            lam=float(lam_fit),
            amplitude=float(amplitude_fit),
            time_start=float(time_min),
            time_end=float(time_max),
            mass_guess=mass_guess,
            abundance=0.0,
            time_axis=dense_time,
            fit_values=fit_curve,
            color="#ff7f0e",
        )
        self._mass_line_counter += 1
        self._mass_lines.append(line)
        self._selected_line_id = line.line_id
        self._update_mass_line_abundances()
        self._update_tables()
        self._update_summary()
        self._refresh_plot()

    def _recompute_mass_line(self, line: MassLineFit, *, preserve_label: bool = False) -> None:
        previous_label = line.label
        previous_guess = line.mass_guess
        previous_auto = nearest_mass_name(previous_guess) if math.isfinite(previous_guess) else ""
        dense_time = np.linspace(line.time_start, line.time_end, 800)
        line.amplitude = max(abs(float(line.amplitude)), 0.0)
        fit_curve = _emg_model(dense_time, line.amplitude, line.mu, abs(line.sigma), abs(line.lam))
        line.time_axis = dense_time
        line.fit_values = fit_curve
        converted = self._time_to_mass(np.array([line.mu], dtype=float))
        if converted.size:
            line.mass_guess = float(converted[0])
        else:
            line.mass_guess = float("nan")
        new_auto = nearest_mass_name(line.mass_guess) if math.isfinite(line.mass_guess) else previous_auto
        if not preserve_label:
            label_clean = previous_label.strip() if isinstance(previous_label, str) else ""
            if not label_clean or (previous_auto and label_clean == previous_auto):
                line.label = new_auto
        self._update_mass_line_abundances()

    def _remove_selected_mass_line(self) -> None:
        selection = self.mass_table.selectionModel()
        if selection is None:
            return
        rows = selection.selectedRows()
        if not rows:
            return
        idx = rows[0].row()
        if 0 <= idx < len(self._mass_lines):
            removed = self._mass_lines.pop(idx)
            if self._selected_line_id == removed.line_id:
                if self._mass_lines:
                    next_idx = min(idx, len(self._mass_lines) - 1)
                    self._selected_line_id = self._mass_lines[next_idx].line_id
                else:
                    self._selected_line_id = None
            self._update_mass_line_abundances()
            self._update_tables()
            self._update_summary()
            self._refresh_plot()
            self._update_inspect_button_state()
    def _update_mass_line_abundances(self) -> None:
        if self._combined is None or self._combined.size == 0:
            for line in self._mass_lines:
                line.abundance = 0.0
            return
        baseline_corrected = self._combined - self._baseline
        positive = np.clip(baseline_corrected, 0.0, None)
        total_area = float(np.trapz(positive, self._time_axis)) if self._time_axis.size else 0.0
        for line in self._mass_lines:
            area = line.window_area() if line.time_end > line.time_start else 0.0
            if total_area > 0:
                line.abundance = float(max(area, 0.0) / total_area)
            else:
                line.abundance = 0.0

    def _update_tables(self) -> None:
        self._block_table_signals = True
        self.mass_table.setRowCount(len(self._mass_lines))
        for row, line in enumerate(self._mass_lines):
            for col, value in enumerate(line.as_row()):
                item = self.mass_table.item(row, col)
                if item is None:
                    item = QTableWidgetItem()
                    self.mass_table.setItem(row, col, item)
                item.setText(str(value))
                if col == 6:
                    flags = item.flags()
                    item.setFlags(flags & ~Qt.ItemFlag.ItemIsEditable)
        self.mass_table.resizeColumnsToContents()
        self._block_selection_signals = True
        if self._selected_line_id is not None:
            found = False
            for row, line in enumerate(self._mass_lines):
                if line.line_id == self._selected_line_id:
                    self.mass_table.selectRow(row)
                    found = True
                    break
            if not found:
                self.mass_table.clearSelection()
                self._selected_line_id = None
        else:
            self.mass_table.clearSelection()
        self._block_selection_signals = False
        self._block_table_signals = False
        self._update_inspect_button_state()

    def _on_mass_table_selection_changed(self) -> None:
        if self._block_selection_signals:
            return
        selection = self.mass_table.selectionModel()
        if selection is None:
            return
        rows = selection.selectedRows()
        if not rows:
            if self._selected_line_id is not None:
                self._selected_line_id = None
                self._refresh_plot()
            self._update_inspect_button_state()
            return
        row = rows[0].row()
        if not (0 <= row < len(self._mass_lines)):
            return
        line = self._mass_lines[row]
        if self._selected_line_id != line.line_id:
            self._selected_line_id = line.line_id
            self._refresh_plot()
        self._update_inspect_button_state()

    def _update_summary(self) -> None:
        entries = sorted(self._mass_lines, key=lambda ln: ln.abundance, reverse=True)
        self.summary_table.setRowCount(len(entries))
        for row, line in enumerate(entries):
            for col, value in enumerate((line.label, f"{line.mass_guess:.3f}", f"{line.abundance * 100.0:.2f}")):
                item = self.summary_table.item(row, col)
                if item is None:
                    item = QTableWidgetItem()
                    self.summary_table.setItem(row, col, item)
                item.setText(str(value))
        if entries:
            guesses = [f"{line.label} ({line.abundance * 100.0:.1f}%)" for line in entries[:3]]
            text = "Sample guess: " + ", ".join(guesses)
        else:
            text = "Sample guess: insufficient data"
        self.sample_guess_label.setText(text)

    # ---- Persistence ----------------------------------------------------
    def _save_to_file(self) -> None:
        if not self._group:
            QMessageBox.warning(self, "No File", "The current event is unavailable or the file is closed.")
            return
        try:
            analysis = self._group.require_group(ANALYSIS_GROUP)
            dust_group = analysis.require_group(DUST_GROUP)
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", f"Unable to create analysis group:\n{exc}")
            return
        try:
            _write_dataset(dust_group, COMBINED_TIME_DATASET, self._time_axis)
            if self._combined is None:
                combined = self._combine_waveforms()
            else:
                combined = self._combined
            if combined is not None:
                _write_dataset(dust_group, COMBINED_DATASET, combined)
            dust_group.attrs["Baseline"] = self._baseline
            dust_group.attrs["MassStretch"] = self._mass_params.get("stretch", 1.0)
            dust_group.attrs["MassShift"] = self._mass_params.get("shift", 0.0)
            if self._mass_lines:
                str_dtype = h5py.string_dtype(encoding="utf-8", length=120)
                table = np.zeros(len(self._mass_lines), dtype=[
                    ("id", "i4"),
                    ("label", str_dtype),
                    ("mu", "f8"),
                    ("sigma", "f8"),
                    ("lam", "f8"),
                    ("amplitude", "f8"),
                    ("time_start", "f8"),
                    ("time_end", "f8"),
                    ("mass", "f8"),
                    ("area", "f8"),
                    ("abundance", "f8"),
                ])
                for idx, line in enumerate(self._mass_lines):
                    table[idx] = (
                        line.line_id,
                        line.label,
                        line.mu,
                        line.sigma,
                        line.lam,
                        line.total_area(),
                        line.time_start,
                        line.time_end,
                        line.mass_guess,
                        line.window_area(),
                        line.abundance,
                    )
                _write_dataset(dust_group, MASS_LINES_DATASET, table)
                fits_group = dust_group.require_group("Fits")
                for key in list(fits_group.keys()):
                    del fits_group[key]
                for line in self._mass_lines:
                    line_group = fits_group.require_group(f"line_{line.line_id}")
                    _write_dataset(line_group, "time", line.time_axis)
                    _write_dataset(line_group, "values", line.fit_values)
            else:
                _safe_delete(dust_group, MASS_LINES_DATASET)
            QMessageBox.information(self, "Saved", "Dust composition analysis saved to file.")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", f"Failed to save analysis:\n{exc}")


def launch_dust_composition_window(h5: h5py.File, event_name: str, parent: Optional[QWidget] = None) -> DustCompositionWindow:
    """Convenience helper used by the main quicklook window."""

    return DustCompositionWindow(h5=h5, event_name=event_name, parent=parent)
