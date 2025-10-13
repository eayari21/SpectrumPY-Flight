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
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np

try:  # pragma: no cover - Qt import guard
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import (
        QAbstractItemView,
        QDoubleSpinBox,
        QFormLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
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
        QDoubleSpinBox,
        QFormLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
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

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

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


def _emg_model(time_values: np.ndarray, mu: float, sigma: float, lam: float) -> np.ndarray:
    """Evaluate an exponentially modified Gaussian (EMG)."""

    arr = np.asarray(time_values, dtype=float)
    if arr.size == 0:
        return np.zeros(0, dtype=float)
    safe_sigma = sigma if abs(sigma) > 1.0e-15 else 1.0e-15
    safe_lambda = lam if abs(lam) > 1.0e-15 else 1.0e-15
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        exponent = np.exp((safe_lambda / 2.0) * (2.0 * mu + safe_lambda * safe_sigma**2 - 2.0 * arr))
    argument = (mu + safe_lambda * safe_sigma**2 - arr) / (np.sqrt(2.0) * safe_sigma)
    return (safe_lambda / 2.0) * exponent * np.erfc(argument)


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
    time_start: float
    time_end: float
    mass_guess: float
    abundance: float
    time_axis: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))
    fit_values: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))
    color: str = "#d62728"

    def parameters(self) -> Tuple[float, float, float]:
        return (self.mu, self.sigma, self.lam)

    def as_row(self) -> Sequence[float | str]:
        return (
            self.label,
            f"{self.mass_guess:.3f}",
            f"{self.mu:.6f}",
            f"{self.sigma:.6f}",
            f"{self.lam:.6f}",
            f"{self.abundance * 100.0:.2f}",
        )


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
        self._mass_lines: List[MassLineFit] = []
        self._mass_line_counter = 0

        self._combined_axis = None
        self._baseline_artist = None
        self._span_selector: Optional[SpanSelector] = None
        self._in_baseline_mode = False
        self._block_table_signals = False

        self._load_datasets()
        self._load_saved_state()

        central = QWidget(self)
        layout = QHBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(10)
        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        layout.addWidget(splitter)
        self.setCentralWidget(central)

        self.figure = Figure(figsize=(8.2, 6.4), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.canvas.mpl_connect("button_press_event", self._on_canvas_click)
        splitter.addWidget(self.canvas)

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
        except Exception:
            self._mass_params = {"stretch": 1.0, "shift": 0.0}
        if MASS_LINES_DATASET in dust_group:
            try:
                table = dust_group[MASS_LINES_DATASET][()]
            except Exception:
                table = None
            if table is not None:
                for entry in table:
                    try:
                        label = entry["label"].decode("utf-8") if hasattr(entry["label"], "decode") else str(entry["label"])
                        mass_line = MassLineFit(
                            line_id=int(entry["id"]),
                            label=label,
                            mu=float(entry["mu"]),
                            sigma=float(entry["sigma"]),
                            lam=float(entry["lam"]),
                            time_start=float(entry["time_start"]),
                            time_end=float(entry["time_end"]),
                            mass_guess=float(entry["mass"]),
                            abundance=float(entry.get("abundance", 0.0)),
                        )
                        self._mass_lines.append(mass_line)
                        self._mass_line_counter = max(self._mass_line_counter, mass_line.line_id + 1)
                    except Exception:
                        continue
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
                    except Exception:
                        line.time_axis = np.zeros(0)
                        line.fit_values = np.zeros(0)

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
        self.mass_table.setColumnCount(6)
        self.mass_table.setHorizontalHeaderLabels([
            "Label",
            "Mass (amu)",
            "μ (µs)",
            "σ (µs)",
            "λ (µs⁻¹)",
            "Abundance (%)",
        ])
        header = self.mass_table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setDefaultSectionSize(130)
        self.mass_table.verticalHeader().setVisible(False)
        self.mass_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.mass_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.mass_table.itemChanged.connect(self._on_mass_table_changed)
        layout.addWidget(self.mass_table)
        self.remove_mass_button = QPushButton("Remove Selected", box)
        self.remove_mass_button.clicked.connect(self._remove_selected_mass_line)
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
        if self.combine_button.isChecked() and self._combined_axis and event.inaxes != self._combined_axis:
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
        except ValueError:
            self._block_table_signals = True
            for col, value in enumerate(line.as_row()):
                self.mass_table.item(row, col).setText(str(value))
            self._block_table_signals = False
            return
        if column in (2, 3, 4):
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
        ax.plot(mass_axis[:n], self._combined[:n], color="#1f77b4", linewidth=1.6, label="Combined TOF")
        ax.set_facecolor("#f9fbff")
        ax.grid(True, alpha=0.35)
        ax.set_xlabel("Mass [amu]", fontsize=14)
        ax.set_ylabel("Signal [pC/Δt]", fontsize=14)
        ax.set_title("Combined TOF Spectrum", fontsize=16, fontweight="bold")
        if self._time_axis.size:
            ax_time = ax.twiny()
            ax_time.set_xlim(self._time_axis.min(), self._time_axis.max())
            ax_time.set_xlabel("Time [µs]", fontsize=13)
        if self._baseline or self.baseline_spin.value() != 0.0:
            self._baseline_artist = ax.axhline(self._baseline, color="#aa3377", linestyle="--", linewidth=1.2, label="Baseline")
        for line in self._mass_lines:
            if line.time_axis.size and line.fit_values.size:
                try:
                    mass = self._time_to_mass(line.time_axis)
                except Exception:
                    mass = line.time_axis
                ax.plot(mass, line.fit_values + self._baseline, linestyle="-", linewidth=1.4, label=line.label)
        if len(ax.lines) > 1:
            ax.legend(loc="best", fontsize=10)

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
        try:
            from scipy.optimize import curve_fit  # type: ignore
        except Exception as exc:
            QMessageBox.warning(self, "Missing SciPy", f"SciPy is required for EMG fitting:\n{exc}")
            return

        try:
            params, _ = curve_fit(
                lambda t, mu, sigma, lam: _emg_model(t, mu, abs(sigma), abs(lam)),
                x,
                y,
                p0=(mu_guess, sigma_guess, lam_guess),
                maxfev=20000,
            )
            mu_fit, sigma_fit, lam_fit = params
            sigma_fit = abs(float(sigma_fit))
            lam_fit = abs(float(lam_fit))
        except Exception:
            QMessageBox.warning(self, "Fit Failed", "Unable to fit an EMG curve to the selected region.")
            return
        dense_time = np.linspace(time_min, time_max, 800)
        fit_curve = _emg_model(dense_time, mu_fit, sigma_fit, lam_fit)
        mass_guess = float(self._time_to_mass(mu_fit))
        label = nearest_mass_name(mass_guess)
        line = MassLineFit(
            line_id=self._mass_line_counter,
            label=label,
            mu=float(mu_fit),
            sigma=float(sigma_fit),
            lam=float(lam_fit),
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
        self._update_mass_line_abundances()
        self._update_tables()
        self._update_summary()
        self._refresh_plot()

    def _recompute_mass_line(self, line: MassLineFit) -> None:
        dense_time = np.linspace(line.time_start, line.time_end, 800)
        fit_curve = _emg_model(dense_time, line.mu, abs(line.sigma), abs(line.lam))
        line.time_axis = dense_time
        line.fit_values = fit_curve
        line.mass_guess = float(self._time_to_mass(line.mu))
        line.label = nearest_mass_name(line.mass_guess)
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
            del self._mass_lines[idx]
            self._update_mass_line_abundances()
            self._update_tables()
            self._update_summary()
            self._refresh_plot()
    def _update_mass_line_abundances(self) -> None:
        if self._combined is None or self._combined.size == 0:
            for line in self._mass_lines:
                line.abundance = 0.0
            return
        baseline_corrected = self._combined - self._baseline
        positive = np.clip(baseline_corrected, 0.0, None)
        total_area = float(np.trapz(positive, self._time_axis)) if self._time_axis.size else 0.0
        for line in self._mass_lines:
            if line.time_axis.size and line.fit_values.size:
                area = np.trapz(np.clip(line.fit_values, 0.0, None), line.time_axis)
            else:
                area = 0.0
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
                if col == 5:
                    flags = item.flags()
                    item.setFlags(flags & ~Qt.ItemFlag.ItemIsEditable)
        self.mass_table.resizeColumnsToContents()
        self._block_table_signals = False

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
                    ("time_start", "f8"),
                    ("time_end", "f8"),
                    ("mass", "f8"),
                    ("abundance", "f8"),
                ])
                for idx, line in enumerate(self._mass_lines):
                    table[idx] = (
                        line.line_id,
                        line.label,
                        line.mu,
                        line.sigma,
                        line.lam,
                        line.time_start,
                        line.time_end,
                        line.mass_guess,
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
