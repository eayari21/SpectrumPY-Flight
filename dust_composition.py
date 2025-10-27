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
import json
import math
import re
import unicodedata
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np

try:  # pragma: no cover - Qt import guard
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import (
        QAbstractItemView,
        QCheckBox,
        QComboBox,
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
        QScrollArea,
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
        QCheckBox,
        QComboBox,
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
        QScrollArea,
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

from line_shapes import (
    double_emg as _line_double_emg,
    emg as _line_emg,
    generalized_normal as _line_generalized_normal,
    gaussian as _line_gaussian,
    hyper_emg as _line_hyper_emg,
    lorentzian as _line_lorentzian,
    voigt as _line_voigt,
)

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


@dataclass(frozen=True)
class ExtraFieldSpec:
    key: str
    label: str
    widget: str  # "spin" or "text"
    decimals: int = 6
    minimum: float = -1.0e6
    maximum: float = 1.0e6
    step: float = 1.0e-3
    tooltip: str = ""


@dataclass(frozen=True)
class RSFPreset:
    key: str
    label: str
    rsf_values: Dict[str, float]


@dataclass(frozen=True)
class RelativeSensitivityResult:
    enabled: bool
    values: Dict[int, float]
    normalised: Dict[int, float]


RSF_PRESETS: Tuple[RSFPreset, ...] = (
    RSFPreset(
        key="puma_pia",
        label="PUMA/PIA (Krueger 1996)",
        rsf_values={"Mg": 3.1, "Si": 1.0, "Fe": 1.1},
    ),
    RSFPreset(
        key="lama_orthopyroxene",
        label="Orthopyroxene, LAMA (Sternglass 1971)",
        rsf_values={"Mg": 5.50, "Si": 1.0, "Fe": 1.12},
    ),
    RSFPreset(
        key="suda_olivine_fo87",
        label="Olivine Fo87, SUDA (Hillier et al. 2018)",
        rsf_values={"Mg": 4.93, "Si": 1.0, "Fe": 1.50},
    ),
    RSFPreset(
        key="hyperdust_olivine_fo91",
        label="Olivine Fo91, Hyperdust (this work)",
        rsf_values={"Mg": 4.97, "Si": 1.0, "Fe": 1.32},
    ),
    RSFPreset(
        key="tofsims",
        label="TOF SIMS (Stephan 2001)",
        rsf_values={"Mg": 5.10, "Si": 1.0, "Fe": 2.40},
    ),
)


@dataclass(frozen=True)
class LineShapeConfig:
    key: str
    display: str
    formula_html: str
    amplitude_label: str = "A (DN·µs)"
    amplitude_tooltip: str = "Integral of the line shape across all time."
    mu_label: str = "μ (µs)"
    mu_tooltip: str = "Location of the line centre in microseconds."
    sigma_label: str = "σ (µs)"
    sigma_tooltip: str = "Width parameter for the selected shape."
    sigma_minimum: float = 1.0e-9
    sigma_step: float = 1.0e-3
    lam_label: str = "λ (µs⁻¹)"
    lam_tooltip: str = "Rate parameter controlling the exponential tail."
    lam_minimum: float = 1.0e-9
    lam_step: float = 1.0e-3
    lam_visible: bool = True
    extras: Tuple[ExtraFieldSpec, ...] = ()
    legend_label: str = "Fit"


def _render_formula(text: str) -> str:
    return "<span style='font-size:16px;'>" + text + "</span>"


LINE_SHAPES: Dict[str, LineShapeConfig] = {
    "emg": LineShapeConfig(
        key="emg",
        display="Exponentially Modified Gaussian (EMG)",
        formula_html=_render_formula(
            "<b>f(t)</b> = <b>A</b>·<b>λ</b>/2 · exp((<b>μ</b> − t)·<b>λ</b> + (<b>λ</b>·<b>σ</b>)²/2) · erfc((<b>μ</b> + (<b>λ</b>·<b>σ</b>)² − t)/(√2·<b>σ</b>))"
        ),
        lam_tooltip="Inverse of the exponential tail time constant (λ = 1/τ).",
        legend_label="EMG fit",
    ),
    "gaussian": LineShapeConfig(
        key="gaussian",
        display="Gaussian",
        formula_html=_render_formula(
            "<b>f(t)</b> = <b>A</b> / (<b>σ</b>√{2π}) · exp(−(t − <b>μ</b>)² / (2<b>σ</b>²))"
        ),
        lam_visible=False,
        sigma_tooltip="Standard deviation of the Gaussian core.",
        legend_label="Gaussian fit",
    ),
    "lorentzian": LineShapeConfig(
        key="lorentzian",
        display="Lorentzian",
        formula_html=_render_formula(
            "<b>f(t)</b> = <b>A</b>/π · <b>γ</b>/((t − <b>μ</b>)² + <b>γ</b>²)"
        ),
        sigma_label="γ (µs)",
        sigma_tooltip="Half-width at half-maximum (HWHM).",
        lam_visible=False,
        legend_label="Lorentzian fit",
    ),
    "voigt": LineShapeConfig(
        key="voigt",
        display="Voigt",
        formula_html=_render_formula(
            "Voigt(<b>μ</b>, <b>σ</b>, <b>γ</b>) = Re[w((t − <b>μ</b> + i<b>γ</b>)/(√2<b>σ</b>))] · <b>A</b>/(<b>σ</b>√{2π})"
        ),
        sigma_tooltip="Gaussian component width (σ).",
        lam_label="γ (µs)",
        lam_tooltip="Lorentzian half-width parameter (γ).",
        lam_minimum=1.0e-9,
        legend_label="Voigt fit",
    ),
    "double_emg": LineShapeConfig(
        key="double_emg",
        display="Double EMG",
        formula_html=_render_formula(
            "<b>f(t)</b> = w₁·EMG(<b>τ₁</b>) + (1−w₁)·EMG(<b>τ₂</b>)"
        ),
        sigma_tooltip="Shared Gaussian width for both components.",
        lam_label="τ₁ (µs)",
        lam_tooltip="Time constant of the first exponential tail (signed).",
        lam_minimum=-1.0,
        lam_step=1.0e-3,
        extras=(
            ExtraFieldSpec(
                key="tau2",
                label="τ₂ (µs):",
                widget="spin",
                decimals=6,
                minimum=-1.0e6,
                maximum=1.0e6,
                step=1.0e-3,
                tooltip="Time constant of the second exponential tail (signed).",
            ),
            ExtraFieldSpec(
                key="weight",
                label="w₁ (0–1):",
                widget="spin",
                decimals=3,
                minimum=0.0,
                maximum=1.0,
                step=0.01,
                tooltip="Mixing weight applied to the first EMG component.",
            ),
        ),
        legend_label="Double EMG fit",
    ),
    "hyper_emg": LineShapeConfig(
        key="hyper_emg",
        display="HyperEMG",
        formula_html=_render_formula(
            "<b>f(t)</b> = (1−∑w)·G + ∑wₗ·EMG(τₗ) + ∑wᵣ·EMG(τᵣ)"
        ),
        sigma_tooltip="Core Gaussian width shared across all components.",
        lam_visible=False,
        extras=(
            ExtraFieldSpec(
                key="taus_left",
                label="Left τ (µs):",
                widget="text",
                tooltip="Comma separated list of negative tail constants (µs).",
            ),
            ExtraFieldSpec(
                key="weights_left",
                label="Left weights:",
                widget="text",
                tooltip="Comma separated weights for left tails (sum ≤ 1).",
            ),
            ExtraFieldSpec(
                key="taus_right",
                label="Right τ (µs):",
                widget="text",
                tooltip="Comma separated list of positive tail constants (µs).",
            ),
            ExtraFieldSpec(
                key="weights_right",
                label="Right weights:",
                widget="text",
                tooltip="Comma separated weights for right tails (sum ≤ 1).",
            ),
        ),
        legend_label="HyperEMG fit",
    ),
    "generalized_normal": LineShapeConfig(
        key="generalized_normal",
        display="Generalized Normal",
        formula_html=_render_formula(
            "<b>f(t)</b> = <b>A</b>·β/(2<b>α</b>Γ(1/β)) · exp(−| (t−<b>μ</b>) / <b>α</b> |^β)"
        ),
        sigma_label="α (µs)",
        sigma_tooltip="Scale parameter controlling the spread (α).",
        sigma_minimum=1.0e-6,
        lam_label="β",
        lam_tooltip="Shape exponent β controlling tail heaviness.",
        lam_minimum=0.1,
        lam_step=0.1,
        legend_label="Generalized normal fit",
    ),
}


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


def _coerce_float_list(value: float | List[float] | str | None) -> List[float]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, np.ndarray)):
        result: List[float] = []
        for item in value:
            try:
                result.append(float(item))
            except Exception:
                continue
        return result
    if isinstance(value, str):
        tokens = [tok.strip() for tok in value.replace(";", ",").split(",") if tok.strip()]
        result: List[float] = []
        for tok in tokens:
            try:
                result.append(float(tok))
            except Exception:
                continue
        return result
    try:
        return [float(value)]
    except Exception:
        return []


def _evaluate_line_shape(
    shape: str,
    time_values: np.ndarray,
    amplitude: float,
    mu: float,
    sigma: float,
    lam: float,
    extras: Optional[Dict[str, float | List[float]]] = None,
) -> np.ndarray:
    arr = np.asarray(time_values, dtype=float)
    if arr.size == 0:
        return np.zeros(0, dtype=float)
    if not math.isfinite(amplitude) or amplitude <= 0.0:
        return np.zeros_like(arr, dtype=float)

    extras = extras or {}
    key = (shape or "emg").lower()
    safe_sigma = max(abs(float(sigma)), 1.0e-12)

    if key == "gaussian":
        return _line_gaussian(arr, mu, safe_sigma, area=float(amplitude))
    if key == "lorentzian":
        gamma = max(abs(float(sigma)), 1.0e-12)
        return _line_lorentzian(arr, mu, gamma, area=float(amplitude))
    if key == "voigt":
        gamma = max(abs(float(lam)), 1.0e-12)
        return _line_voigt(arr, mu, safe_sigma, gamma, area=float(amplitude))
    if key == "double_emg":
        tau1 = float(lam)
        if not math.isfinite(tau1) or abs(tau1) < 1.0e-9:
            tau1 = 1.0e-9
        tau2_raw = extras.get("tau2")
        try:
            tau2 = float(tau2_raw) if tau2_raw is not None else -tau1
        except Exception:
            tau2 = -tau1
        if not math.isfinite(tau2) or abs(tau2) < 1.0e-9:
            tau2 = -math.copysign(1.0e-9, tau2 if tau2 != 0 else tau1)
        weight_raw = extras.get("weight")
        try:
            weight = float(weight_raw)
        except Exception:
            weight = 0.5
        weight = float(np.clip(weight, 0.0, 1.0))
        return _line_double_emg(arr, mu, safe_sigma, tau1, tau2, w1=weight, area=float(amplitude))
    if key == "hyper_emg":
        taus_left = [-abs(val) for val in _coerce_float_list(extras.get("taus_left"))]
        taus_right = [abs(val) for val in _coerce_float_list(extras.get("taus_right"))]
        weights_left = _coerce_float_list(extras.get("weights_left"))
        weights_right = _coerce_float_list(extras.get("weights_right"))
        return _line_hyper_emg(
            arr,
            mu,
            safe_sigma,
            taus_left=taus_left,
            taus_right=taus_right,
            weights_left=weights_left,
            weights_right=weights_right,
            area=float(amplitude),
        )
    if key == "generalized_normal":
        alpha = max(abs(float(sigma)), 1.0e-9)
        beta = abs(float(lam)) if math.isfinite(lam) else 1.0
        beta = max(beta, 0.1)
        return _line_generalized_normal(arr, mu, alpha, beta, area=float(amplitude))

    safe_lambda = abs(float(lam)) if math.isfinite(lam) else 1.0
    safe_lambda = max(safe_lambda, 1.0e-9)
    tau = 1.0 / safe_lambda
    return _line_emg(arr, mu, safe_sigma, tau, area=float(amplitude))


def _line_window_area(
    shape: str,
    amplitude: float,
    mu: float,
    sigma: float,
    lam: float,
    extras: Optional[Dict[str, float | List[float]]],
    start: Optional[float],
    end: Optional[float],
) -> float:
    if start is None or end is None or not (math.isfinite(start) and math.isfinite(end)) or end <= start:
        return 0.0
    sample = np.linspace(start, end, 800)
    values = _evaluate_line_shape(shape, sample, amplitude, mu, sigma, lam, extras)
    if values.size == 0:
        return 0.0
    with np.errstate(invalid="ignore"):
        return float(np.trapz(np.clip(values, 0.0, None), sample))


def _estimate_amplitude_from_curve(time_axis: np.ndarray, fit_values: np.ndarray, line: MassLineFit) -> float:
    """Infer the amplitude scaling for ``line`` from sampled ``fit_values``."""

    if time_axis.size == 0 or fit_values.size == 0:
        return 0.0
    unit_model = _evaluate_line_shape(
        line.shape,
        time_axis,
        1.0,
        line.mu,
        abs(line.sigma),
        abs(line.lam),
        line.extra_params,
    )
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


ELEMENT_PATTERN = re.compile(r"[A-Z][a-z]?")
KNOWN_ELEMENTS = {
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
}


def _normalise_formula(formula: str) -> str:
    simplified = unicodedata.normalize("NFKD", formula or "")
    ascii_text = simplified.encode("ascii", "ignore").decode("ascii")
    ascii_text = ascii_text.replace(" ", "")
    ascii_text = ascii_text.replace("-", "")
    ascii_text = ascii_text.replace("~", "")
    ascii_text = ascii_text.replace("_", "")
    ascii_text = ascii_text.replace(".", "")
    return ascii_text


def _coerce_element(symbol: str) -> Optional[str]:
    if not symbol:
        return None
    if symbol in KNOWN_ELEMENTS:
        return symbol
    if len(symbol) == 2 and symbol[1] in {"x", "y", "z"}:
        base = symbol[0]
        if base in KNOWN_ELEMENTS:
            return base
    if len(symbol) > 1:
        base = symbol[0]
        if base in KNOWN_ELEMENTS:
            return base
    return symbol if symbol in KNOWN_ELEMENTS else None


def _formula_to_elements(formula: str) -> Tuple[str, ...]:
    tokens = set()
    cleaned = _normalise_formula(formula)
    for match in ELEMENT_PATTERN.findall(cleaned):
        element = _coerce_element(match)
        if element:
            tokens.add(element)
    return tuple(sorted(tokens))


@dataclass(frozen=True)
class SampleDefinition:
    category: str
    name: str
    formula: str
    elements: Tuple[str, ...]


SAMPLE_DATA: Tuple[Tuple[str, str, str], ...] = (
    ("Silicate", "Forsterite (olivine endmember)", "Mg2SiO4"),
    ("Silicate", "Fayalite (olivine endmember)", "Fe2SiO4"),
    ("Silicate", "Olivine (solid solution)", "MgFeSiO"),
    ("Silicate", "Enstatite (OPX)", "MgSiO"),
    ("Silicate", "Ferrosilite (OPX)", "FeSiO"),
    ("Silicate", "Pigeonite (CPX)", "MgFeCaSiO"),
    ("Silicate", "Diopside (CPX)", "CaMgSiO"),
    ("Silicate", "Hedenbergite (CPX)", "CaFeSiO"),
    ("Silicate", "Augite (CPX)", "CaNaMgFeAlTiSiO"),
    ("Silicate", "Anorthite (feldspar)", "CaAlSiO"),
    ("Silicate", "Albite (feldspar)", "NaAlSiO"),
    ("Silicate", "Orthoclase (feldspar)", "KAlSiO"),
    ("Silicate", "Labradorite (feldspar SS)", "NaCaAlSiO"),
    ("Silicate", "Quartz / silica polymorphs", "SiO"),
    ("Silicate", "Amorphous silicate (GEMS-like)", "MgFeSiO"),
    ("Hydrated silicate", "Serpentine", "MgSiOH"),
    ("Hydrated silicate", "Saponite (smectite)", "CaMgFeSiAlOH"),
    ("Hydrated silicate", "Montmorillonite (smectite)", "NaCaAlMgSiOH"),
    ("Hydrated silicate", "Chlorite", "MgFeAlSiOH"),
    ("Oxide/Spinel", "Magnetite", "FeO"),
    ("Oxide/Spinel", "Hematite", "FeO"),
    ("Oxide/Spinel", "Wollastonite", "FeO"),
    ("Oxide/Spinel", "Corundum", "AlO"),
    ("Oxide/Spinel", "Spinel", "MgAlO"),
    ("Oxide/Spinel", "Chromite", "FeCrO"),
    ("Oxide/Spinel", "Ilmenite", "FeTiO"),
    ("Oxide/Spinel", "Rutile", "TiO"),
    ("Oxide/Spinel", "Periclase", "MgO"),
    ("Oxide/Spinel", "Hercynite", "FeAlO"),
    ("Sulfide", "Troilite", "FeS"),
    ("Sulfide", "Pyrrhotite", "FeS"),
    ("Sulfide", "Pentlandite", "FeNiS"),
    ("Sulfide", "Pyrite", "FeS"),
    ("Sulfide", "Oldhamite", "CaS"),
    ("Sulfide", "Niningerite", "MgFeS"),
    ("Sulfide", "Sphalerite", "ZnS"),
    ("Sulfide", "Galena", "PbS"),
    ("Sulfide", "Cubanite", "CuFeS"),
    ("Metal/Alloy", "Fe-Ni metal (kamacite/taenite)", "FeNi"),
    ("Metal/Alloy", "Kamacite", "FeNi"),
    ("Metal/Alloy", "Taenite", "FeNi"),
    ("Metal/Alloy", "Nickel", "Ni"),
    ("Metal/Alloy", "Cobalt", "Co"),
    ("Carbide/Nitride/Boride", "Silicon carbide (moissanite)", "SiC"),
    ("Carbide/Nitride/Boride", "Titanium carbide", "TiC"),
    ("Carbide/Nitride/Boride", "Graphite intergrowth with TiC (hibonite nodules context)", "CTi"),
    ("Carbide/Nitride/Boride", "Silicon nitride", "SiN"),
    ("Carbide/Nitride/Boride", "Titanium nitride", "TiN"),
    ("Carbide/Nitride/Boride", "Boron carbide (rare)", "BC"),
    ("Carbonaceous", "Graphite", "C"),
    ("Carbonaceous", "Nanodiamond / diamond", "C"),
    ("Carbonaceous", "Amorphous carbon", "C"),
    ("Chondrite", "CI carbonaceous chondrite matrix", "MgFeSiONC"),
    ("Chondrite", "CM carbonaceous chondrite matrix", "MgFeSiONCHS"),
    ("Chondrite", "CO/CV carbonaceous chondrite matrix", "MgFeSiAlCaONC"),
    ("Chondrite", "Ordinary chondrite (H/L)", "MgFeSiON"),
    ("Chondrule", "Type I (Mg-rich) chondrule melt", "MgSiON"),
    ("Chondrule", "Type II (Fe-rich) chondrule melt", "MgFeSiON"),
    ("Chondrule", "Al-rich CAI-like inclusion", "CaAlTiSiON"),
    ("Organic (simple)", "Polycyclic aromatic hydrocarbons (PAHs)", "CH"),
    ("Organic (simple)", "Formaldehyde", "CHO"),
    ("Organic (simple)", "Methanol", "CHO"),
    ("Organic (simple)", "Formamide", "CHNO"),
    ("Organic (simple)", "Hydrogen cyanide", "HCN"),
    ("Organic (simple)", "Acetaldehyde", "CHO"),
    ("Organic (simple)", "Acetonitrile", "CHN"),
    ("Ice/Volatile", "Water ice", "HO"),
    ("Ice/Volatile", "Carbon dioxide ice", "CO"),
    ("Ice/Volatile", "Carbon monoxide ice", "CO"),
    ("Ice/Volatile", "Methane ice", "CH"),
    ("Ice/Volatile", "Ammonia ice", "NH"),
    ("Ice/Volatile", "Molecular nitrogen ice", "N"),
    ("Ice/Volatile", "Hydrogen sulfide ice", "HS"),
    ("Ice/Volatile", "Sulfur dioxide ice", "SO"),
    ("Carbonate", "Calcite", "CaCO"),
    ("Carbonate", "Dolomite", "CaMgCO"),
    ("Carbonate", "Magnesite", "MgCO"),
    ("Carbonate", "Siderite", "FeCO"),
    ("Sulfate", "Gypsum", "CaSOH"),
    ("Sulfate", "Anhydrite", "CaSO"),
    ("Sulfate", "Jarosite", "KFeSOH"),
    ("Halide/Salt", "Halite", "NaCl"),
    ("Halide/Salt", "Sylvite", "KCl"),
    ("Halide/Salt", "Perchlorate (generic)", "ClO"),
    ("Phosphate", "Apatite (fluor/hydroxyl/chloroapatite)", "CaPOFClH"),
    ("Phosphate", "Whitlockite/Merrillite (meteoritic Ca-phosphates)", "CaMgFeNaPOH"),
    ("Phosphide", "Schreibersite", "FeNiP"),
    ("Presolar oxide", "Hibonite", "CaAlO"),
    ("Presolar oxide", "Perovskite", "CaTiO"),
    ("Presolar oxide", "Spinel (presolar)", "MgAlO"),
    ("Presolar SiC", "Beta-SiC (mainstream presolar SiC)", "SiC"),
    ("Presolar grain", "Graphite spherule", "C"),
    ("Presolar grain", "Amorphous silicate (GEMS-like presolar)", "MgFeSiO"),
    ("Presolar grain", "Nanodiamond aggregate", "C"),
)


SAMPLE_LIBRARY: Tuple[SampleDefinition, ...] = tuple(
    SampleDefinition(category=category, name=name, formula=formula, elements=_formula_to_elements(formula))
    for category, name, formula in SAMPLE_DATA
)

SAMPLE_BY_NAME: Dict[str, SampleDefinition] = {sample.name: sample for sample in SAMPLE_LIBRARY}


@dataclass(frozen=True)
class SampleMatch:
    sample: SampleDefinition
    coverage: float
    score: float


@dataclass(frozen=True)
class MixtureMatch:
    primary: SampleDefinition
    secondary: SampleDefinition
    fractions: Tuple[float, float]
    coverage: float
    score: float


def _element_weights_from_lines(mass_lines: Sequence["MassLineFit"]) -> Dict[str, float]:
    weights: Dict[str, float] = {}
    for line in mass_lines:
        if not math.isfinite(line.abundance) or line.abundance <= 0.0:
            continue
        label = line.label or ""
        elements = _formula_to_elements(label)
        if not elements and math.isfinite(line.mass_guess):
            elements = _formula_to_elements(nearest_mass_name(line.mass_guess))
        if not elements:
            continue
        share = line.abundance / max(len(elements), 1)
        for element in elements:
            weights[element] = weights.get(element, 0.0) + share
    if not weights:
        return {}
    normaliser = sum(weights.values())
    if normaliser <= 0.0:
        return weights
    for element in list(weights):
        weights[element] = max(weights[element] / normaliser, 0.0)
    return weights


def _rank_sample_matches(weights: Dict[str, float], limit: int = 8) -> List[SampleMatch]:
    if not weights:
        return []
    observed_elements = set(weights)
    matches: List[SampleMatch] = []
    for sample in SAMPLE_LIBRARY:
        sample_elements = set(sample.elements)
        if not sample_elements:
            continue
        overlap = sample_elements & observed_elements
        if not overlap:
            continue
        coverage = sum(weights.get(elem, 0.0) for elem in sample_elements)
        balance = len(overlap) / len(sample_elements)
        penalty = 0.02 * len(sample_elements - observed_elements)
        score = max(coverage * (0.5 + 0.5 * balance) - penalty, 0.0)
        if score <= 0.0:
            continue
        matches.append(SampleMatch(sample=sample, coverage=coverage, score=score))
    matches.sort(key=lambda match: (match.score, match.coverage), reverse=True)
    return matches[:limit]


def _best_mixture_match(weights: Dict[str, float], candidates: Sequence[SampleMatch]) -> Optional[MixtureMatch]:
    if not weights or len(candidates) < 2:
        return None
    observed_elements = set(weights)
    best: Optional[MixtureMatch] = None
    best_score = 0.0
    for idx, first in enumerate(candidates[:-1]):
        first_elements = set(first.sample.elements)
        for second in candidates[idx + 1 :]:
            second_elements = set(second.sample.elements)
            union = first_elements | second_elements
            overlap = union & observed_elements
            if not overlap:
                continue
            coverage = sum(weights.get(elem, 0.0) for elem in union)
            balance = len(overlap) / len(union)
            score = coverage * (0.5 + 0.5 * balance)
            if score <= best_score:
                continue
            first_contrib = sum(weights.get(elem, 0.0) for elem in first_elements)
            second_contrib = sum(weights.get(elem, 0.0) for elem in second_elements)
            total = first_contrib + second_contrib
            if total <= 0.0:
                fractions = (0.5, 0.5)
            else:
                fractions = (first_contrib / total, second_contrib / total)
            best = MixtureMatch(
                primary=first.sample,
                secondary=second.sample,
                fractions=fractions,
                coverage=coverage,
                score=score,
            )
            best_score = score
    return best


def analyse_sample_matches(mass_lines: Sequence["MassLineFit"]) -> Tuple[Optional[SampleMatch], Optional[MixtureMatch]]:
    weights = _element_weights_from_lines(mass_lines)
    matches = _rank_sample_matches(weights)
    best = matches[0] if matches else None
    mixture = _best_mixture_match(weights, matches[:6]) if matches else None
    return best, mixture
SPECIES_CHOICES: List[Tuple[str, float]] = [
    ("1H", 1.008),
    ("2H", 2.014),
    ("3He", 3.016),
    ("4He", 4.003),
    ("12C", 12.0),
    ("13C", 13.003),
    ("14N", 14.003),
    ("15N", 15.0),
    ("16O", 15.995),
    ("17O", 16.999),
    ("18O", 17.999),
    ("23Na", 22.99),
    ("24Mg", 23.985),
    ("25Mg", 24.986),
    ("26Mg", 25.983),
    ("27Al", 26.982),
    ("28Si", 27.977),
    ("29Si", 28.976),
    ("30Si", 29.974),
    ("31P", 30.974),
    ("32S", 31.972),
    ("33S", 32.971),
    ("34S", 33.968),
    ("36S", 35.967),
    ("35Cl", 34.969),
    ("37Cl", 36.966),
    ("39K", 38.964),
    ("41K", 40.962),
    ("40Ca", 39.963),
    ("44Ca", 43.955),
    ("45Sc", 44.956),
    ("47Ti", 46.952),
    ("48Ti", 47.948),
    ("49Ti", 48.948),
    ("50Ti", 49.944),
    ("50Cr", 49.946),
    ("52Cr", 51.941),
    ("53Cr", 52.941),
    ("54Cr", 53.939),
    ("55Mn", 54.938),
    ("54Fe", 53.94),
    ("56Fe", 55.935),
    ("57Fe", 56.935),
    ("58Fe", 57.933),
    ("58Ni", 57.935),
    ("60Ni", 59.931),
    ("61Ni", 60.931),
    ("62Ni", 61.928),
    ("64Ni", 63.928),
    ("63Cu", 62.93),
    ("65Cu", 64.928),
    ("64Zn", 63.929),
    ("66Zn", 65.926),
    ("67Zn", 66.927),
    ("68Zn", 67.924),
    ("70Zn", 69.925),
    ("H", 1.008),
    ("CH", 13.018),
    ("NH", 15.015),
    ("OH", 17.007),
    ("H2", 2.016),
    ("CH2", 14.026),
    ("NH2", 16.022),
    ("H2O", 18.015),
    ("CH3", 15.034),
    ("NH3", 17.03),
    ("CO", 28.011),
    ("N2", 28.014),
    ("NO", 30.006),
    ("O2", 31.998),
    ("CO2", 44.009),
    ("H2S", 34.076),
    ("SO", 48.06),
    ("SO2", 64.059),
    ("CS", 44.071),
    ("COS", 60.072),
    ("CS2", 76.131),
    ("C2", 24.021),
    ("C2H", 25.029),
    ("C2H2", 26.037),
    ("C2H3", 27.045),
    ("C2H4", 28.053),
    ("C2H5", 29.061),
    ("C2H6", 30.069),
    ("C3", 36.032),
    ("C3H", 37.04),
    ("C3H2", 38.047),
    ("C3H3", 39.055),
    ("C3H4", 40.063),
    ("C3H5", 41.071),
    ("C3H6", 42.079),
    ("C3H7", 43.086),
    ("C3H8", 44.094),
    ("C4", 48.042),
    ("C4H", 49.05),
    ("C4H2", 50.058),
    ("C4H3", 51.066),
    ("C4H4", 52.074),
    ("C4H5", 53.082),
    ("C4H6", 54.09),
    ("C4H7", 55.097),
    ("C4H8", 56.105),
    ("C4H9", 57.113),
    ("C5H5", 65.061),
    ("C5H7", 67.077),
    ("C5H9", 69.093),
    ("C6H5", 77.059),
    ("C6H6", 78.047),
    ("C6H7", 79.055),
    ("C7H7", 91.055),
    ("C7H8", 92.063),
    ("C8H7", 103.067),
    ("C9H7", 115.078),
    ("C10H8", 128.062),
    ("NaO", 38.989),
    ("NaOH", 39.997),
    ("MgO", 40.304),
    ("MgOH", 41.312),
    ("AlO", 42.98),
    ("SiH", 29.093),
    ("SiH2", 30.101),
    ("SiH3", 31.108),
    ("SiC", 40.096),
    ("SiN", 42.092),
    ("SiO", 44.084),
    ("SiO2", 60.083),
    ("Si2", 56.17),
    ("Si2O", 72.169),
    ("Si2O2", 88.168),
    ("Si2O3", 104.167),
    ("Si2O4", 120.166),
    ("Si3O4", 148.249),
    ("Si3O5", 164.25),
    ("Si3O6", 180.249),
    ("Si4O6", 208.334),
    ("Si4O7", 224.333),
    ("Si4O8", 240.332),
    ("FeO", 71.844),
    ("FeOH", 72.852),
    ("Fe2O", 127.535),
    ("Fe2O2", 143.534),
    ("Fe2O3", 159.533),
    ("MgSi", 52.39),
    ("MgSiO", 68.389),
    ("MgSiO2", 84.388),
    ("FeSi", 83.93),
    ("FeSiO", 99.929),
    ("FeSiO2", 115.928),
    ("CaO", 56.077),
    ("CaOH", 57.085),
    ("Ca2O", 120.156),
    ("CaSi", 68.163),
    ("CaSiO", 84.162),
    ("CaSiO2", 100.161),
    ("NaCl", 58.44),
    ("KCl", 74.548),
    ("MgCl", 59.365),
    ("CaCl", 75.528),
    ("FeCl", 91.905),
    ("S2", 64.12),
    ("FeS", 87.905),
    ("FeS2", 119.965),
    ("MgS", 56.365),
    ("CaS", 72.138),
    ("CuS", 95.606),
    ("ZnS", 97.44),
    ("Na2S", 78.04),
    ("MgSiO3", 100.387),
    ("FeSiO3", 131.927),
    ("CaSiO3", 116.16),
    ("Mg2SiO4", 140.691),
    ("Fe2SiO4", 203.771),
    ("MgFeSiO4", 172.231),
    ("NaAlSi3O8", 262.218),
    ("KAlSi3O8", 278.327),
    ("CaAl2Si2O8", 278.203),
]

SPECIES_BY_LABEL: Dict[str, float] = {name: mass for name, mass in SPECIES_CHOICES}


def _species_display(label: str, mass: float) -> str:
    return f"{label} ({mass:.3f} amu)"


def _species_for_label(label: str) -> Optional[Tuple[str, float]]:
    mass = SPECIES_BY_LABEL.get(label.strip())
    if mass is None:
        return None
    return label.strip(), mass


def _infer_element_symbol(label: str) -> Optional[str]:
    """Best-effort extraction of an elemental symbol from *label*."""

    if not label:
        return None
    cleaned = label.strip()
    # Prefer any explicitly assigned species name if available.
    token = cleaned.split()[0]
    token = re.sub(r"[^A-Za-z]", "", token)
    if not token:
        return None
    if len(token) == 1:
        return token.upper()
    return token[0].upper() + token[1:].lower()


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
    """Container describing a single analytical mass-line fit."""

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
    assigned_species: Optional[str] = None
    assigned_mass: Optional[float] = None
    shape: str = "emg"
    extra_params: Dict[str, float | List[float]] = field(default_factory=dict)

    def parameters(self) -> Tuple[float, float, float, float]:
        return (self.amplitude, self.mu, self.sigma, self.lam)

    def evaluate(self, time_values: np.ndarray) -> np.ndarray:
        return _evaluate_line_shape(
            self.shape,
            time_values,
            self.amplitude,
            self.mu,
            self.sigma,
            self.lam,
            self.extra_params,
        )

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
        return max(
            _line_window_area(
                self.shape,
                self.amplitude,
                self.mu,
                self.sigma,
                self.lam,
                self.extra_params,
                self.time_start,
                self.time_end,
            ),
            0.0,
        )

    def total_area(self) -> float:
        if not math.isfinite(self.amplitude):
            return 0.0
        return max(float(self.amplitude), 0.0)


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
    """Focused inspector for visualising and editing a single mass line."""

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
        if not isinstance(self._line.extra_params, dict):
            self._line.extra_params = {}
        else:
            self._line.extra_params = dict(self._line.extra_params)
        if self._line.shape not in LINE_SHAPES:
            self._line.shape = "emg"
        self._time_axis = np.asarray(time_axis, dtype=float).ravel()
        self._signal = np.asarray(signal, dtype=float).ravel()
        self._baseline = float(baseline)
        self._source_name = source_name
        self._mass_converter = mass_converter
        self._result: Optional[Dict[str, float | str]] = None
        self._assigned_species: Optional[str] = line.assigned_species
        self._assigned_mass: Optional[float] = line.assigned_mass
        self._block_species_signal = False

        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(18, 18, 18, 18)
        outer_layout.setSpacing(12)

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        content = QWidget(scroll_area)
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(12)

        scroll_area.setWidget(content)
        outer_layout.addWidget(scroll_area, stretch=1)

        self.header_label = QLabel("", content)
        self.header_label.setStyleSheet("font-size: 20px; font-weight: 600;")
        content_layout.addWidget(self.header_label)

        config = LINE_SHAPES.get(self._line.shape, LINE_SHAPES["emg"])

        self.formula_label = QLabel(content)
        self.formula_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.formula_label.setTextFormat(Qt.TextFormat.RichText)
        self.formula_label.setWordWrap(True)
        self.formula_label.setText(config.formula_html)
        content_layout.addWidget(self.formula_label)

        source_label = QLabel(f"Signal source: {html.escape(self._source_name)}", content)
        source_label.setStyleSheet("color: #495057; font-size: 13px;")
        content_layout.addWidget(source_label)

        figure_container = QWidget(content)
        figure_layout = QVBoxLayout(figure_container)
        figure_layout.setContentsMargins(0, 0, 0, 0)
        figure_layout.setSpacing(6)

        self.figure = Figure(figsize=(6.5, 3.8), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(900, 420)
        self.toolbar = NavigationToolbar2QT(self.canvas, figure_container)
        figure_layout.addWidget(self.toolbar)
        figure_layout.addWidget(self.canvas)

        content_layout.addWidget(figure_container, stretch=1)

        parameter_box = QGroupBox("Editable parameters", content)
        form = QFormLayout(parameter_box)
        form.setContentsMargins(12, 12, 12, 12)
        form.setSpacing(8)

        self.species_combo = QComboBox(parameter_box)
        self.species_combo.addItem("Custom / manual entry", userData=None)
        for name, mass in SPECIES_CHOICES:
            self.species_combo.addItem(_species_display(name, mass), (name, mass))
        self.species_combo.currentIndexChanged.connect(self._on_species_changed)
        form.addRow("Species:", self.species_combo)

        self._block_shape_signal = False
        self.shape_combo = QComboBox(parameter_box)
        for key, shape_cfg in LINE_SHAPES.items():
            self.shape_combo.addItem(shape_cfg.display, userData=key)
        self.shape_combo.currentIndexChanged.connect(self._on_shape_changed)
        form.addRow("Line shape:", self.shape_combo)

        self.label_edit = QLineEdit(parameter_box)
        self.label_edit.setText(self._line.label)
        self.label_edit.textChanged.connect(self._on_label_changed)
        form.addRow("Label:", self.label_edit)

        time_min = float(np.nanmin(self._time_axis)) if self._time_axis.size else -1.0e6
        time_max = float(np.nanmax(self._time_axis)) if self._time_axis.size else 1.0e6

        self.mu_label = QLabel("μ (µs):", parameter_box)
        self.mu_spin = QDoubleSpinBox(parameter_box)
        self.mu_spin.setDecimals(6)
        self.mu_spin.setRange(time_min, time_max)
        self.mu_spin.setValue(self._line.mu)
        self.mu_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow(self.mu_label, self.mu_spin)

        self.sigma_label = QLabel("σ (µs):", parameter_box)
        self.sigma_spin = QDoubleSpinBox(parameter_box)
        self.sigma_spin.setDecimals(6)
        self.sigma_spin.setRange(1.0e-9, 1.0e6)
        self.sigma_spin.setValue(max(abs(self._line.sigma), 1.0e-6))
        self.sigma_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow(self.sigma_label, self.sigma_spin)

        self.lambda_label = QLabel("λ (µs⁻¹):", parameter_box)
        self.lambda_spin = QDoubleSpinBox(parameter_box)
        self.lambda_spin.setDecimals(6)
        self.lambda_spin.setRange(-1.0e6, 1.0e6)
        lam_initial = float(self._line.lam) if math.isfinite(self._line.lam) else 1.0
        if abs(lam_initial) < 1.0e-9:
            lam_initial = 1.0e-6
        self.lambda_spin.setValue(lam_initial)
        self.lambda_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow(self.lambda_label, self.lambda_spin)

        self.amplitude_label = QLabel("A (DN·µs):", parameter_box)
        self.amplitude_spin = QDoubleSpinBox(parameter_box)
        self.amplitude_spin.setDecimals(6)
        self.amplitude_spin.setRange(1.0e-12, 1.0e12)
        self.amplitude_spin.setValue(max(abs(self._line.amplitude), 1.0e-9))
        self.amplitude_spin.valueChanged.connect(self._on_parameter_changed)
        form.addRow(self.amplitude_label, self.amplitude_spin)

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

        self.extra_box = QGroupBox("Additional parameters", parameter_box)
        self.extra_box.setVisible(False)
        self.extra_form = QFormLayout(self.extra_box)
        self.extra_form.setContentsMargins(8, 8, 8, 8)
        self.extra_form.setSpacing(6)
        self._extra_fields: Dict[str, QWidget] = {}
        form.addRow(self.extra_box)

        content_layout.addWidget(parameter_box)

        self._set_shape_index(self._line.shape)
        self._apply_shape_config(self._line.shape)

        self.mass_hint_label = QLabel("", content)
        self.mass_hint_label.setStyleSheet("font-size: 14px; font-style: italic; color: #495057;")
        content_layout.addWidget(self.mass_hint_label)

        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        outer_layout.addWidget(button_box)

        initial_index = 0
        if self._assigned_species:
            assigned_mass = self._assigned_mass
            tolerance = 5.0e-3
            for idx, (name, mass) in enumerate(SPECIES_CHOICES, start=1):
                if name == self._assigned_species and (
                    assigned_mass is None or not math.isfinite(assigned_mass) or abs(assigned_mass - mass) <= tolerance
                ):
                    initial_index = idx
                    if assigned_mass is None or not math.isfinite(assigned_mass):
                        self._assigned_mass = mass
                    break

        self._set_species_index(initial_index)
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
        parts: List[str] = []
        if self._assigned_mass is not None and math.isfinite(self._assigned_mass):
            species = self._assigned_species or (self.label_edit.text().strip() or "Species")
            parts.append(f"Assigned mass: {self._assigned_mass:.3f} amu ({species})")
        if math.isfinite(mass_value):
            parts.append(f"Estimated mass from μ: {mass_value:.3f} amu")
        else:
            parts.append("Estimated mass from μ: unavailable")
        self.mass_hint_label.setText(" | ".join(parts))

    def _update_plot(self) -> None:
        shape_key = self._current_shape_key()
        config = self._shape_config(shape_key)
        amplitude = float(max(self.amplitude_spin.value(), 1.0e-12))
        mu = float(self.mu_spin.value())
        sigma = float(max(abs(self.sigma_spin.value()), max(config.sigma_minimum, 1.0e-9)))
        lam_value = float(self.lambda_spin.value())
        if config.lam_visible:
            lam = lam_value
            if config.lam_minimum > 0.0:
                lam = max(lam, config.lam_minimum)
        else:
            lam = lam_value if math.isfinite(lam_value) else 1.0
        extras = self._collect_extra_parameters(validate=False)
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
        fit_values = _evaluate_line_shape(shape_key, fit_time, amplitude, mu, sigma, lam, extras)

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
        ax.plot(
            fit_time,
            fit_plot,
            color="#d62728",
            linewidth=2.2,
            label=config.legend_label,
        )
        ax.set_xlabel("Time (µs)")
        ax.set_ylabel("Signal (DN – baseline)")
        ax.set_yscale("symlog", linthresh=1.0e-3)
        ax.set_xlim(window_min, window_max)
        ax.set_title("Zoomed mass line view", fontsize=14, fontweight="bold")
        ax.legend(loc="best")
        ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.35)
        self.canvas.draw_idle()

    def _shape_config(self, key: str) -> LineShapeConfig:
        return LINE_SHAPES.get(key, LINE_SHAPES["emg"])

    def _current_shape_key(self) -> str:
        data = self.shape_combo.currentData()
        if isinstance(data, str) and data in LINE_SHAPES:
            return data
        return "emg"

    def _set_shape_index(self, key: str) -> None:
        target = key if key in LINE_SHAPES else "emg"
        for idx in range(self.shape_combo.count()):
            if self.shape_combo.itemData(idx) == target:
                self._block_shape_signal = True
                try:
                    self.shape_combo.setCurrentIndex(idx)
                finally:
                    self._block_shape_signal = False
                return
        if self.shape_combo.count():
            self._block_shape_signal = True
            try:
                self.shape_combo.setCurrentIndex(0)
            finally:
                self._block_shape_signal = False

    def _apply_shape_config(self, key: str) -> None:
        config = self._shape_config(key)
        self.formula_label.setText(config.formula_html)
        self.amplitude_label.setText(config.amplitude_label)
        self.amplitude_spin.setToolTip(config.amplitude_tooltip)
        self.mu_label.setText(config.mu_label)
        self.mu_spin.setToolTip(config.mu_tooltip)
        self.sigma_label.setText(config.sigma_label)
        self.sigma_spin.setToolTip(config.sigma_tooltip)
        sigma_min = max(config.sigma_minimum, 1.0e-9)
        self.sigma_spin.setRange(sigma_min, 1.0e6)
        self.sigma_spin.setSingleStep(config.sigma_step)

        self.lambda_label.setText(config.lam_label)
        self.lambda_spin.setToolTip(config.lam_tooltip)
        self.lambda_label.setVisible(config.lam_visible)
        self.lambda_spin.setVisible(config.lam_visible)
        if config.lam_visible:
            if config.lam_minimum < 0.0:
                self.lambda_spin.setRange(-1.0e6, 1.0e6)
            else:
                lam_min = max(config.lam_minimum, 1.0e-9)
                self.lambda_spin.setRange(lam_min, 1.0e6)
            self.lambda_spin.setSingleStep(config.lam_step)
            if config.lam_minimum >= 0.0 and self.lambda_spin.value() <= 0.0:
                self.lambda_spin.setValue(max(config.lam_minimum, 1.0e-9))

        # Rebuild extra parameter widgets
        for widget in self._extra_fields.values():
            if isinstance(widget, QWidget):
                self.extra_form.removeRow(widget)
                try:
                    widget.deleteLater()
                except RuntimeError:
                    # Qt may have already scheduled the widget for deletion.
                    # In that case we simply ignore the error so shape
                    # reconfiguration can proceed without aborting.
                    pass
        self._extra_fields.clear()

        if config.extras:
            self.extra_box.show()
            for spec in config.extras:
                if spec.widget == "spin":
                    field = QDoubleSpinBox(self.extra_box)
                    field.setDecimals(spec.decimals)
                    field.setRange(spec.minimum, spec.maximum)
                    field.setSingleStep(spec.step)
                    if spec.tooltip:
                        field.setToolTip(spec.tooltip)
                    base_value = self._line.extra_params.get(spec.key)
                    if isinstance(base_value, (int, float)):
                        field.setValue(float(base_value))
                    elif spec.key == "weight":
                        field.setValue(0.5)
                    elif spec.key.startswith("tau"):
                        tau_default = float(self.lambda_spin.value())
                        if spec.key == "tau2":
                            tau_default = -abs(tau_default)
                        field.setValue(tau_default)
                    self.extra_form.addRow(spec.label, field)
                    self._extra_fields[spec.key] = field
                else:
                    field = QLineEdit(self.extra_box)
                    if spec.tooltip:
                        field.setToolTip(spec.tooltip)
                    base_value = self._line.extra_params.get(spec.key)
                    if isinstance(base_value, str):
                        field.setText(base_value)
                    elif isinstance(base_value, (list, tuple, np.ndarray)):
                        field.setText(", ".join(f"{float(val):.6g}" for val in base_value))
                    self.extra_form.addRow(spec.label, field)
                    self._extra_fields[spec.key] = field
        else:
            self.extra_box.hide()

    def _collect_extra_parameters(self, *, validate: bool) -> Dict[str, float | List[float]]:
        extras: Dict[str, float | List[float]] = {}
        config = self._shape_config(self._current_shape_key())
        for spec in config.extras:
            widget = self._extra_fields.get(spec.key)
            if widget is None:
                continue
            if isinstance(widget, QDoubleSpinBox):
                value = float(widget.value())
                if validate and spec.key.startswith("tau") and abs(value) < 1.0e-9:
                    raise ValueError("Time constants must be non-zero.")
                if validate and spec.key == "weight" and not (0.0 <= value <= 1.0):
                    raise ValueError("w₁ must be between 0 and 1.")
                extras[spec.key] = value
            elif isinstance(widget, QLineEdit):
                text = widget.text().strip()
                if not text:
                    extras[spec.key] = []
                else:
                    values = _coerce_float_list(text)
                    if validate and not values:
                        raise ValueError(f"{spec.label} requires at least one value.")
                    extras[spec.key] = values
        return extras

    def _on_shape_changed(self, _index: int) -> None:
        if self._block_shape_signal:
            return
        key = self._current_shape_key()
        self._line.shape = key
        self._line.extra_params = {}
        self._apply_shape_config(key)
        self._update_plot()

    def _on_label_changed(self, _text: str) -> None:
        label = self.label_edit.text().strip()
        if self._assigned_species and label != self._assigned_species:
            self._assigned_species = None
            self._assigned_mass = None
            self._set_species_index(0)
        self._update_header()
        self._update_mass_hint()

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
        shape_key = self._current_shape_key()
        config = self._shape_config(shape_key)
        amplitude = float(max(self.amplitude_spin.value(), 1.0e-12))
        mu = float(self.mu_spin.value())
        sigma = float(max(abs(self.sigma_spin.value()), max(config.sigma_minimum, 1.0e-9)))
        lam_value = float(self.lambda_spin.value())
        if config.lam_visible:
            lam = lam_value
        else:
            lam = lam_value if math.isfinite(lam_value) else (self._line.lam or 1.0)
        if sigma <= 0.0:
            raise ValueError("The width parameter must be positive.")
        if config.key in {"emg", "voigt", "generalized_normal"} and lam <= 0.0:
            raise ValueError("The selected line shape requires a positive secondary parameter.")
        if config.key == "double_emg" and abs(lam) < 1.0e-9:
            raise ValueError("τ₁ must be non-zero.")
        start = float(self.start_spin.value())
        end = float(self.end_spin.value())
        if end <= start:
            raise ValueError("The end time must be greater than the start time.")
        if self._assigned_mass is not None and math.isfinite(self._assigned_mass):
            mass_guess = float(self._assigned_mass)
        else:
            try:
                mass_guess = float(self._mass_converter(mu))
            except Exception:
                mass_guess = float("nan")
        extras = self._collect_extra_parameters(validate=True)
        self._line.extra_params = extras
        return {
            "label": label,
            "amplitude": amplitude,
            "mu": mu,
            "sigma": sigma,
            "lam": lam,
            "time_start": start,
            "time_end": end,
            "mass_guess": mass_guess,
            "assigned_species": self._assigned_species or "",
            "assigned_mass": self._assigned_mass,
            "shape": shape_key,
            "extra_params": extras,
        }

    def _set_species_index(self, index: int) -> None:
        if 0 <= index < self.species_combo.count():
            self._block_species_signal = True
            try:
                self.species_combo.setCurrentIndex(index)
            finally:
                self._block_species_signal = False

    def _on_species_changed(self, index: int) -> None:
        if self._block_species_signal:
            return
        data = self.species_combo.itemData(index)
        if data is None:
            self._assigned_species = None
            self._assigned_mass = None
            self._update_mass_hint()
            return
        species_label, species_mass = data
        self._assigned_species = str(species_label)
        self._assigned_mass = float(species_mass)
        self.label_edit.blockSignals(True)
        try:
            self.label_edit.setText(self._assigned_species)
        finally:
            self.label_edit.blockSignals(False)
        self._update_header()
        self._update_mass_hint()

    def collected_values(self) -> Optional[Dict[str, float | str]]:
        return self._result

# ---------------------------------------------------------------------------
# Relative sensitivity factors dialog
# ---------------------------------------------------------------------------


class RelativeSensitivityDialog(QDialog):
    """Dialog allowing users to apply relative sensitivity factors (RSFs)."""

    def __init__(
        self,
        parent: QWidget,
        mass_lines: Sequence[MassLineFit],
        existing_values: Optional[Dict[int, float]] = None,
        enabled: bool = False,
    ):
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle("Relative Sensitivity Factors")
        self.setMinimumSize(720, 520)

        self._result: Optional[RelativeSensitivityResult] = None
        self._rows: List[Dict[str, object]] = []
        self._current_values: Dict[int, float] = {}
        self._current_normalised: Dict[int, float] = {}
        self._current_total: float = 0.0

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(12)

        title = QLabel("Balance relative abundances with RSFs", self)
        title.setStyleSheet("font-size: 18px; font-weight: 600;")
        layout.addWidget(title)

        description = QLabel(
            "Select the fitted mass lines to include, assign RSF weights, and preview the renormalised "
            "relative abundances. Presets from recent laboratory measurements are provided for "
            "convenience.",
            self,
        )
        description.setWordWrap(True)
        description.setStyleSheet("color: #495057;")
        layout.addWidget(description)

        self.enable_checkbox = QCheckBox("Enable RSF weighting", self)
        self.enable_checkbox.setChecked(bool(enabled))
        self.enable_checkbox.toggled.connect(self._on_enable_toggled)
        layout.addWidget(self.enable_checkbox)

        preset_layout = QHBoxLayout()
        preset_layout.setSpacing(8)
        preset_label = QLabel("Preset instrument:", self)
        preset_label.setStyleSheet("font-weight: 500;")
        preset_layout.addWidget(preset_label)
        self.preset_combo = QComboBox(self)
        self.preset_combo.addItem("Custom values", userData=None)
        for preset in RSF_PRESETS:
            self.preset_combo.addItem(preset.label, userData=preset)
        self.preset_combo.currentIndexChanged.connect(self._on_preset_changed)
        preset_layout.addWidget(self.preset_combo, stretch=1)
        layout.addLayout(preset_layout)

        table = QTableWidget(self)
        table.setColumnCount(7)
        table.setHorizontalHeaderLabels(
            [
                "Use",
                "Mass line",
                "Element",
                "Mass (amu)",
                "Measured (%)",
                "RSF",
                "RSF-normalised (%)",
            ]
        )
        header = table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setDefaultSectionSize(128)
        table.verticalHeader().setVisible(False)
        table.setAlternatingRowColors(True)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(table, stretch=1)
        self.table = table

        existing_values = existing_values or {}

        for row_index, line in enumerate(mass_lines):
            table.insertRow(row_index)
            checkbox = QCheckBox(table)
            default_checked = line.line_id in existing_values
            if not existing_values and max(float(line.abundance), 0.0) > 0.0:
                default_checked = True
            checkbox.setChecked(default_checked)
            checkbox.toggled.connect(self._update_adjusted_values)
            table.setCellWidget(row_index, 0, checkbox)

            label_item = QTableWidgetItem(line.label or f"Line {line.line_id}")
            label_item.setFlags(label_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row_index, 1, label_item)

            element_label = line.assigned_species or line.label
            element_symbol = _infer_element_symbol(element_label)
            element_item = QTableWidgetItem(element_symbol or "–")
            element_item.setFlags(element_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            element_item.setTextAlignment(int(Qt.AlignmentFlag.AlignCenter))
            table.setItem(row_index, 2, element_item)

            mass_item = QTableWidgetItem(f"{line.mass_guess:.3f}")
            mass_item.setFlags(mass_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            mass_item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
            table.setItem(row_index, 3, mass_item)

            measured_fraction = max(float(line.abundance), 0.0)
            measured_item = QTableWidgetItem(f"{measured_fraction * 100.0:.2f}")
            measured_item.setFlags(measured_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            measured_item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
            table.setItem(row_index, 4, measured_item)

            spin = QDoubleSpinBox(table)
            spin.setDecimals(3)
            spin.setRange(0.0, 20.0)
            spin.setSingleStep(0.05)
            spin.setValue(float(existing_values.get(line.line_id, 1.0)))
            spin.valueChanged.connect(self._update_adjusted_values)
            table.setCellWidget(row_index, 5, spin)

            adjusted_item = QTableWidgetItem("0.00")
            adjusted_item.setFlags(adjusted_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            adjusted_item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
            table.setItem(row_index, 6, adjusted_item)

            row_info: Dict[str, object] = {
                "line_id": int(line.line_id),
                "checkbox": checkbox,
                "spin": spin,
                "base_fraction": measured_fraction,
                "element": element_symbol,
                "adjust_item": adjusted_item,
            }
            self._rows.append(row_info)

        self.summary_label = QLabel("", self)
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-style: italic; color: #495057;")
        layout.addWidget(self.summary_label)

        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

        self._on_enable_toggled(bool(enabled))

    def _on_enable_toggled(self, checked: bool) -> None:
        for row in self._rows:
            checkbox = row["checkbox"]
            spin = row["spin"]
            if isinstance(checkbox, QCheckBox):
                checkbox.setEnabled(True)
            if isinstance(spin, QDoubleSpinBox):
                spin.setEnabled(checked and isinstance(checkbox, QCheckBox) and checkbox.isChecked())
        self.preset_combo.setEnabled(checked)
        self._update_adjusted_values()

    def _on_preset_changed(self, index: int) -> None:
        data = self.preset_combo.itemData(index)
        if not isinstance(data, RSFPreset):
            return
        for row in self._rows:
            element = row.get("element")
            spin = row.get("spin")
            if not isinstance(spin, QDoubleSpinBox):
                continue
            if not isinstance(element, str):
                continue
            value = data.rsf_values.get(element)
            if value is None:
                continue
            spin.blockSignals(True)
            spin.setValue(float(value))
            spin.blockSignals(False)
        self._update_adjusted_values()

    def _update_adjusted_values(self) -> None:
        enabled = self.enable_checkbox.isChecked()
        total = 0.0
        contributions: Dict[int, float] = {}
        selected_count = 0
        for row in self._rows:
            checkbox = row.get("checkbox")
            spin = row.get("spin")
            base_fraction = float(row.get("base_fraction", 0.0))
            line_id = int(row.get("line_id", -1))
            include = False
            if isinstance(checkbox, QCheckBox):
                include = checkbox.isChecked()
            if isinstance(spin, QDoubleSpinBox):
                spin.setEnabled(enabled and include)
            if not enabled or not include:
                contributions[line_id] = 0.0
                continue
            selected_count += 1
            rsf = float(spin.value()) if isinstance(spin, QDoubleSpinBox) else 1.0
            weight = max(base_fraction, 0.0) * max(rsf, 0.0)
            contributions[line_id] = weight
            total += weight

        for row in self._rows:
            line_id = int(row.get("line_id", -1))
            adjust_item = row.get("adjust_item")
            if not isinstance(adjust_item, QTableWidgetItem):
                continue
            if not enabled:
                percent = float(row.get("base_fraction", 0.0)) * 100.0
            else:
                weight = contributions.get(line_id, 0.0)
                percent = (weight / total * 100.0) if total > 0.0 else 0.0
            adjust_item.setText(f"{percent:.2f}")

        if not enabled:
            self.summary_label.setText("RSFs disabled. Displaying measured relative abundances.")
            self._current_values = {}
            self._current_normalised = {}
            self._current_total = 0.0
            return

        if selected_count == 0:
            self.summary_label.setText("Select at least one mass line to apply RSFs.")
            self._current_values = {}
            self._current_normalised = {}
            self._current_total = 0.0
            return

        if total <= 0.0:
            self.summary_label.setText(
                "Selected lines have zero total abundance; RSF weights cannot be normalised yet."
            )
        else:
            self.summary_label.setText(
                f"RSFs applied to {selected_count} line(s); values renormalised to 100%."
            )

        values: Dict[int, float] = {}
        normalised: Dict[int, float] = {}
        for row in self._rows:
            checkbox = row.get("checkbox")
            spin = row.get("spin")
            line_id = int(row.get("line_id", -1))
            if not isinstance(checkbox, QCheckBox) or not checkbox.isChecked():
                continue
            values[line_id] = float(spin.value()) if isinstance(spin, QDoubleSpinBox) else 1.0
            weight = contributions.get(line_id, 0.0)
            normalised[line_id] = (weight / total) if total > 0.0 else 0.0

        self._current_values = values
        self._current_normalised = normalised
        self._current_total = total

    def accept(self) -> None:  # type: ignore[override]
        if self.enable_checkbox.isChecked():
            if not self._current_values:
                QMessageBox.warning(
                    self,
                    "No Lines Selected",
                    "Select at least one mass line before enabling RSFs.",
                )
                return
            if self._current_total <= 0.0:
                QMessageBox.warning(
                    self,
                    "Zero Total Abundance",
                    "The selected lines have zero total abundance and cannot be normalised.",
                )
                return
        enabled = self.enable_checkbox.isChecked()
        self._result = RelativeSensitivityResult(
            enabled=enabled,
            values=dict(self._current_values) if enabled else {},
            normalised=dict(self._current_normalised) if enabled else {},
        )
        super().accept()

    def result_data(self) -> Optional[RelativeSensitivityResult]:
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

    def __init__(
        self,
        h5: h5py.File,
        event_name: str,
        *,
        event_names: Optional[Sequence[str]] = None,
        on_event_changed: Optional[Callable[[str], None]] = None,
        parent: Optional[QWidget] = None,
    ):
        super().__init__(parent)

        self._h5 = h5
        self._event = event_name
        self._group: Optional[h5py.Group] = None

        if self._h5 is not None:
            grp = self._h5.get(event_name)
            if isinstance(grp, h5py.Group):
                self._group = grp

        self._event_names: List[str] = list(dict.fromkeys(event_names or []))
        if event_name and event_name not in self._event_names:
            self._event_names.append(event_name)
        self._external_event_callback = on_event_changed
        self._event_selector: Optional[QComboBox] = None
        self._block_event_selector = False

        self.setWindowTitle(f"Dust Composition Analysis — Event {event_name}")
        self.resize(1320, 880)

        self._initialise_analysis_state()

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
        self.canvas.setMinimumSize(1100, 720)
        self.canvas.mpl_connect("button_press_event", self._on_canvas_click)

        self.toolbar = NavigationToolbar2QT(self.canvas, figure_container)

        figure_layout.addWidget(self.toolbar)
        figure_layout.addWidget(self.canvas)

        figure_scroll = QScrollArea(self)
        figure_scroll.setWidgetResizable(True)
        figure_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        figure_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        figure_scroll.setWidget(figure_container)
        splitter.addWidget(figure_scroll)

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
        self._update_mass_axis_lock()

    # ---- Data loading ---------------------------------------------------
    def _initialise_analysis_state(self) -> None:
        """Reset analysis state so a new event can be loaded cleanly."""

        self._time_axis = np.zeros(0)
        self._waveforms = {}
        self._combined = None
        self._combined_cached_mass = None
        self._baseline = 0.0
        self._mass_params = {"stretch": 1.0, "shift": 0.0}
        self._mass_params_loaded = False
        self._mass_axis_lock_level = 0
        self._mass_lines = []
        self._mass_line_counter = 0
        self._selected_line_id = None
        self._manual_sample_guess = None
        self._auto_sample_match = None
        self._auto_mixture_match = None
        self._block_sample_guess_signal = False
        self._anchor_displayed_line_id = None

        self._combined_axis = None
        self._combined_time_axis = None
        self._baseline_artist = None
        span_selector = getattr(self, "_span_selector", None)
        if span_selector is not None:
            try:
                span_selector.disconnect_events()
            except Exception:
                pass
        self._span_selector = None
        self._in_baseline_mode = False
        self._block_table_signals = False
        self._block_selection_signals = False
        self._block_species_assignment = False
        self._rsf_enabled = False
        self._rsf_values = {}
        self._rsf_normalised = {}

    def _ensure_event_listed(self, event_name: str) -> None:
        if not event_name:
            return
        if event_name not in self._event_names:
            self._event_names.append(event_name)
            if self._event_selector is not None:
                self._event_selector.addItem(event_name)
        if self._event_selector is not None:
            self._event_selector.setEnabled(len(self._event_names) > 1)

    def _sync_event_selector(self, event_name: Optional[str]) -> None:
        if self._event_selector is None or not event_name:
            return
        if event_name not in self._event_names:
            self._ensure_event_listed(event_name)
        try:
            index = self._event_names.index(event_name)
        except ValueError:
            return
        if self._event_selector.currentIndex() != index:
            self._block_event_selector = True
            try:
                self._event_selector.setCurrentIndex(index)
            finally:
                self._block_event_selector = False
        self._event_selector.setEnabled(len(self._event_names) > 1)

    def _apply_loaded_state_to_controls(self) -> None:
        if hasattr(self, "baseline_spin"):
            self.baseline_spin.blockSignals(True)
            self.baseline_spin.setValue(self._baseline)
            self.baseline_spin.blockSignals(False)
        if hasattr(self, "mass_stretch_spin"):
            self.mass_stretch_spin.blockSignals(True)
            self.mass_stretch_spin.setValue(self._mass_params.get("stretch", 1.0))
            self.mass_stretch_spin.blockSignals(False)
        if hasattr(self, "mass_shift_spin"):
            self.mass_shift_spin.blockSignals(True)
            self.mass_shift_spin.setValue(self._mass_params.get("shift", 0.0))
            self.mass_shift_spin.blockSignals(False)
        if hasattr(self, "combine_button"):
            self.combine_button.blockSignals(True)
            self.combine_button.setChecked(False)
            self.combine_button.blockSignals(False)
        if hasattr(self, "baseline_button"):
            self.baseline_button.blockSignals(True)
            self.baseline_button.setChecked(False)
            self.baseline_button.blockSignals(False)
        if hasattr(self, "add_mass_button"):
            self.add_mass_button.blockSignals(True)
            self.add_mass_button.setChecked(False)
            self.add_mass_button.blockSignals(False)
        if hasattr(self, "anchor_mass_spin"):
            self.anchor_mass_spin.setEnabled(False)
        if hasattr(self, "anchor_button"):
            self.anchor_button.setEnabled(False)
        self._set_span_selector_active(False)
        self._in_baseline_mode = False

    def _on_event_selector_changed(self, index: int) -> None:
        if self._block_event_selector:
            return
        if 0 <= index < len(self._event_names):
            self._switch_event(self._event_names[index], from_user=True)

    def _switch_event(self, event_name: str, *, from_user: bool) -> None:
        if not event_name:
            return
        self._ensure_event_listed(event_name)
        previous = self._event
        if event_name == previous:
            self._sync_event_selector(event_name)
            return

        self._event = event_name
        self._group = None
        if self._h5 is not None:
            grp = self._h5.get(event_name)
            if isinstance(grp, h5py.Group):
                self._group = grp

        self.setWindowTitle(f"Dust Composition Analysis — Event {event_name}")

        self._initialise_analysis_state()
        self._load_datasets()
        self._load_saved_state()
        self._apply_loaded_state_to_controls()
        self._refresh_plot(initial=True)
        self._update_tables()
        self._update_summary()
        self._update_mass_axis_lock()
        self._sync_event_selector(event_name)

        try:
            self.statusBar().showMessage(f"Loaded event {event_name}", 4000)
        except Exception:
            pass

        if from_user and self._external_event_callback is not None:
            try:
                self._external_event_callback(event_name)
            except Exception:
                pass

    def set_current_event(self, event_name: Optional[str]) -> None:
        """Synchronise the window with the parent quicklook view."""

        if not event_name:
            return
        self._switch_event(event_name, from_user=False)

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
                                line,
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
                assigned_label_raw = _get_field(entry, ("assigned_species",), default="")
                if assigned_label_raw in (None, "", b""):
                    assigned_label = ""
                else:
                    try:
                        assigned_label = str(assigned_label_raw)
                    except Exception:
                        assigned_label = ""
                assigned_mass = _coerce_float(
                    _get_field(entry, ("assigned_mass",), default=float("nan")),
                    default=float("nan"),
                )
                if not label:
                    label = nearest_mass_name(mass_guess)
                shape_raw = _get_field(entry, ("shape", "line_shape", "model"), default="emg")
                if isinstance(shape_raw, bytes):
                    try:
                        shape_raw = shape_raw.decode("utf-8")
                    except Exception:
                        shape_raw = shape_raw.decode("latin-1", "ignore")
                shape_key = str(shape_raw).strip().lower() or "emg"
                if shape_key not in LINE_SHAPES:
                    shape_key = "emg"
                extras_raw = _get_field(entry, ("extras", "extra", "extra_params"), default="")
                extra_params: Dict[str, float | List[float]] = {}
                if extras_raw not in (None, "", b""):
                    try:
                        if isinstance(extras_raw, bytes):
                            extras_raw = extras_raw.decode("utf-8")
                        parsed = json.loads(extras_raw)
                        if isinstance(parsed, dict):
                            extra_params = {
                                str(key): value for key, value in parsed.items() if isinstance(key, str)
                            }
                    except Exception:
                        extra_params = {}
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
                    shape=shape_key,
                    extra_params=extra_params,
                )
                if assigned_label:
                    line.label = assigned_label
                    line.assigned_species = assigned_label
                if math.isfinite(assigned_mass):
                    line.assigned_mass = assigned_mass
                    line.mass_guess = assigned_mass
                elif assigned_label and assigned_label in SPECIES_BY_LABEL:
                    species_mass = SPECIES_BY_LABEL[assigned_label]
                    line.assigned_mass = species_mass
                    line.mass_guess = species_mass
                elif assigned_label:
                    line.assigned_mass = None
                else:
                    species_match = _species_for_label(line.label)
                    if species_match is not None:
                        _, species_mass = species_match
                        if math.isfinite(line.mass_guess) and abs(line.mass_guess - species_mass) <= 5.0e-3:
                            line.assigned_species = species_match[0]
                            line.assigned_mass = species_mass
                            line.mass_guess = species_mass
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

    def _assigned_reference_lines(self) -> List[MassLineFit]:
        lines: List[MassLineFit] = []
        for line in self._mass_lines:
            mass = line.assigned_mass
            if mass is None or not math.isfinite(mass):
                continue
            try:
                mu = float(line.mu)
            except Exception:
                continue
            if not math.isfinite(mu):
                continue
            lines.append(line)
        return lines

    def _apply_mass_axis_lock(self) -> None:
        level = getattr(self, "_mass_axis_lock_level", 0)
        if hasattr(self, "mass_shift_spin"):
            self.mass_shift_spin.setEnabled(level == 0)
        if hasattr(self, "mass_stretch_spin"):
            self.mass_stretch_spin.setEnabled(level < 2)
        if hasattr(self, "auto_mass_button"):
            self.auto_mass_button.setEnabled(level == 0)
        # Anchor controls mirror the lock level: once reference lines are
        # assigned the shift should be governed solely by those references.
        self._refresh_anchor_controls()

    def _update_mass_axis_lock(self) -> None:
        assigned_lines = self._assigned_reference_lines()
        previous_level = getattr(self, "_mass_axis_lock_level", 0)
        new_level = 2 if len(assigned_lines) >= 2 else 1 if len(assigned_lines) == 1 else 0

        if new_level == 0:
            self._mass_axis_lock_level = 0
            self._apply_mass_axis_lock()
            return

        if new_level == 1:
            if previous_level < 1 and assigned_lines:
                line = assigned_lines[0]
                try:
                    stretch = float(self._mass_params.get("stretch", 1.0))
                    mu = float(line.mu)
                    mass = float(line.assigned_mass)
                except Exception:
                    stretch = float("nan")
                    mu = float("nan")
                    mass = float("nan")
                if (
                    math.isfinite(stretch)
                    and abs(stretch) >= 1.0e-12
                    and math.isfinite(mu)
                    and math.isfinite(mass)
                ):
                    shift = mu - mass / stretch
                    if math.isfinite(shift):
                        previous = dict(self._mass_params)
                        self._mass_params["shift"] = shift
                        self._apply_mass_params_update(previous)
            self._mass_axis_lock_level = 1
            self._apply_mass_axis_lock()
            return

        if new_level >= 2:
            if previous_level < 2:
                pairs: List[Tuple[float, float]] = []
                for line in assigned_lines:
                    try:
                        mu = float(line.mu)
                        mass = float(line.assigned_mass)
                    except Exception:
                        continue
                    if math.isfinite(mu) and math.isfinite(mass):
                        pairs.append((mu, mass))
                if len(pairs) >= 2:
                    previous = dict(self._mass_params)
                    self._estimate_mass_axis(pairs)
                    self._apply_mass_params_update(previous)
            self._mass_axis_lock_level = 2
            self._apply_mass_axis_lock()

    # ---- UI construction ------------------------------------------------
    def _build_controls(self) -> None:
        self._build_event_selector()
        self._build_action_buttons()
        self._build_baseline_controls()
        self._build_mass_axis_controls()
        self._build_mass_line_tools()
        self._build_mass_assignment_controls()
        self._build_mass_line_table()
        self._build_summary_section()
        self.control_layout.addStretch(1)

    def _build_event_selector(self) -> None:
        self._ensure_event_listed(self._event)
        if not self._event_names:
            return

        box = QGroupBox("Event", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(6)

        description = QLabel(
            "Choose another event to review without leaving the analysis window.",
            box,
        )
        description.setWordWrap(True)
        description.setStyleSheet("color: #495057; font-size: 13px;")
        layout.addWidget(description)

        combo = QComboBox(box)
        combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        combo.setMinimumHeight(40)
        combo.setStyleSheet("font-size: 15px;")
        for name in self._event_names:
            combo.addItem(name)
        index = 0
        if self._event and self._event in self._event_names:
            index = self._event_names.index(self._event)
        self._block_event_selector = True
        try:
            combo.setCurrentIndex(index)
        finally:
            self._block_event_selector = False
        combo.currentIndexChanged.connect(self._on_event_selector_changed)
        combo.setEnabled(len(self._event_names) > 1)
        layout.addWidget(combo)

        self._event_selector = combo
        self.control_layout.addWidget(box)

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

        anchor_row = QWidget(box)
        anchor_row_layout = QHBoxLayout(anchor_row)
        anchor_row_layout.setContentsMargins(0, 0, 0, 0)
        anchor_row_layout.setSpacing(6)

        self.anchor_mass_spin = QDoubleSpinBox(anchor_row)
        self.anchor_mass_spin.setDecimals(6)
        self.anchor_mass_spin.setRange(-1e6, 1e6)
        self.anchor_mass_spin.setValue(1.0)
        anchor_row_layout.addWidget(self.anchor_mass_spin)

        self.anchor_button = QPushButton("Anchor selected line", anchor_row)
        self.anchor_button.setToolTip(
            "Adjust the mass-axis shift so the selected mass line matches the anchor mass."
        )
        self.anchor_button.clicked.connect(self._anchor_selected_line)
        anchor_row_layout.addWidget(self.anchor_button)

        layout.addRow("Anchor mass:", anchor_row)

        self.anchor_mass_spin.setEnabled(False)
        self.anchor_button.setEnabled(False)
        self.auto_mass_button = QPushButton("Auto-calc axis", box)
        self.auto_mass_button.setToolTip(
            "Estimate the mass-axis stretch and shift using existing mass lines."
        )
        self.auto_mass_button.clicked.connect(self._auto_calculate_mass_axis)
        layout.addRow("", self.auto_mass_button)
        self.control_layout.addWidget(box)

    def _build_mass_line_tools(self) -> None:
        box = QGroupBox("Mass Line Tools", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(6)
        self.add_mass_button = QPushButton("Add Mass Line", box)
        self.add_mass_button.setToolTip("Select a region on the combined plot to fit an EMG mass line.")
        self.add_mass_button.setCheckable(True)
        self.add_mass_button.toggled.connect(self._toggle_mass_line_mode)
        layout.addWidget(self.add_mass_button)
        self.rsf_button = QPushButton("Relative Sensitivity Factors…", box)
        self.rsf_button.setToolTip(
            "Enable and adjust relative sensitivity factors for the currently fitted mass lines."
        )
        self.rsf_button.clicked.connect(self._open_relative_sensitivity_dialog)
        layout.addWidget(self.rsf_button)
        self.control_layout.addWidget(box)

    def _build_mass_assignment_controls(self) -> None:
        box = QGroupBox("Mass Assignment", self.control_panel)
        layout = QVBoxLayout(box)
        layout.setSpacing(6)
        description = QLabel(
            "Assign the selected mass line to a reference species to lock its mass and update the axis.",
            box,
        )
        description.setWordWrap(True)
        description.setStyleSheet("color: #495057; font-size: 12px;")
        layout.addWidget(description)
        self.mass_assignment_label = QLabel("No mass line selected.", box)
        self.mass_assignment_label.setWordWrap(True)
        self.mass_assignment_label.setStyleSheet("font-weight: 500;")
        layout.addWidget(self.mass_assignment_label)
        self.mass_species_combo = QComboBox(box)
        self.mass_species_combo.addItem("Assign selected line…", userData=None)
        for name, mass in SPECIES_CHOICES:
            self.mass_species_combo.addItem(_species_display(name, mass), (name, mass))
        self.mass_species_combo.currentIndexChanged.connect(self._on_mass_species_chosen)
        layout.addWidget(self.mass_species_combo)
        self.control_layout.addWidget(box)
        self._set_mass_assignment_index(0)

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
        self.rsf_status_label = QLabel("", box)
        self.rsf_status_label.setWordWrap(True)
        self.rsf_status_label.setStyleSheet("color: #0b7285; font-style: italic;")
        self.rsf_status_label.hide()
        layout.addWidget(self.rsf_status_label)
        self.sample_guess_label = QLabel("Sample guess: –", box)
        self.sample_guess_label.setWordWrap(True)
        layout.addWidget(self.sample_guess_label)
        self.sample_guess_combo = QComboBox(box)
        self.sample_guess_combo.setEditable(True)
        self.sample_guess_combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.sample_guess_combo.setToolTip("Select a material from the reference library or enter your own description.")
        self.sample_guess_combo.addItem("Auto (closest match)", userData=None)
        for sample in SAMPLE_LIBRARY:
            self.sample_guess_combo.addItem(f"{sample.category} — {sample.name}", userData=sample.name)
        self.sample_guess_combo.currentIndexChanged.connect(self._on_sample_guess_changed)
        self.sample_guess_combo.editTextChanged.connect(self._on_sample_guess_text_changed)
        layout.addWidget(self.sample_guess_combo)
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

    def _open_relative_sensitivity_dialog(self) -> None:
        if not self._mass_lines:
            QMessageBox.information(
                self,
                "No Mass Lines",
                "Fit one or more mass lines before applying relative sensitivity factors.",
            )
            return
        existing = dict(self._rsf_values) if self._rsf_enabled else {}
        dialog = RelativeSensitivityDialog(
            self,
            self._mass_lines,
            existing_values=existing,
            enabled=self._rsf_enabled,
        )
        result_code = dialog.exec()
        try:
            accepted = result_code == QDialog.DialogCode.Accepted  # type: ignore[attr-defined]
        except Exception:
            accepted = int(result_code) == int(QDialog.DialogCode.Accepted)
        if not accepted:
            return
        result = dialog.result_data()
        if result is None:
            return
        if result.enabled and result.values:
            self._rsf_enabled = True
            self._rsf_values = dict(result.values)
            self._rsf_normalised = dict(result.normalised)
            self._recalculate_rsf_normalised()
        else:
            self._rsf_enabled = False
            self._rsf_values.clear()
            self._rsf_normalised.clear()
        self._update_summary()

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
            shape="emg",
        )
        dense_time = np.linspace(line.time_start, line.time_end, 800)
        line.time_axis = dense_time
        line.fit_values = line.evaluate(dense_time)
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
        shape_value = values.get("shape", line.shape)
        if isinstance(shape_value, str) and shape_value in LINE_SHAPES:
            line.shape = shape_value
        else:
            line.shape = "emg"
        extras_value = values.get("extra_params")
        if isinstance(extras_value, dict):
            line.extra_params = {
                str(key): value for key, value in extras_value.items()
            }
        else:
            line.extra_params = {}
        line.time_start = float(values.get("time_start", line.time_start))
        line.time_end = float(values.get("time_end", line.time_end))
        assigned_species = str(values.get("assigned_species", "")).strip()
        assigned_mass_value = values.get("assigned_mass")
        try:
            assigned_mass = float(assigned_mass_value)
        except Exception:
            assigned_mass = None
        if assigned_mass is not None and math.isfinite(assigned_mass):
            line.assigned_species = assigned_species or line.label
            line.assigned_mass = assigned_mass
            line.mass_guess = assigned_mass
            if assigned_species:
                line.label = assigned_species
        else:
            line.assigned_species = None
            line.assigned_mass = None
        if line.assigned_mass is None:
            try:
                mass_guess = float(values.get("mass_guess", line.mass_guess))
            except Exception:
                mass_guess = float("nan")
            if math.isfinite(mass_guess):
                line.mass_guess = mass_guess
        self._recompute_mass_line(line, preserve_label=bool(line.assigned_species))
        self._selected_line_id = line.line_id
        self._update_tables()
        self._update_summary()
        self._refresh_plot()
        self._update_mass_axis_lock()

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
        self._update_mass_axis_from_lines(show_warning=True)

    def _on_mass_params_changed(self) -> None:
        self._mass_params["stretch"] = float(self.mass_stretch_spin.value())
        self._mass_params["shift"] = float(self.mass_shift_spin.value())
        self._combined_cached_mass = None
        self._refresh_plot()
        self._update_tables()
        self._update_summary()

    def _refresh_anchor_controls(self) -> None:
        if not hasattr(self, "anchor_mass_spin") or not hasattr(self, "anchor_button"):
            return
        if getattr(self, "_mass_axis_lock_level", 0) >= 1:
            self.anchor_mass_spin.setEnabled(False)
            self.anchor_button.setEnabled(False)
            self._anchor_displayed_line_id = None
            return
        line = self._current_mass_line()
        if line is None:
            self.anchor_mass_spin.setEnabled(False)
            self.anchor_button.setEnabled(False)
            self._anchor_displayed_line_id = None
            return
        self.anchor_mass_spin.setEnabled(True)
        self.anchor_button.setEnabled(True)
        if self._anchor_displayed_line_id != line.line_id:
            target_mass: Optional[float] = None
            if line.assigned_mass is not None and math.isfinite(line.assigned_mass):
                target_mass = float(line.assigned_mass)
            elif math.isfinite(line.mass_guess):
                target_mass = float(line.mass_guess)
            if target_mass is not None and math.isfinite(target_mass):
                self.anchor_mass_spin.blockSignals(True)
                self.anchor_mass_spin.setValue(target_mass)
                self.anchor_mass_spin.blockSignals(False)
            self._anchor_displayed_line_id = line.line_id

    def _update_mass_axis_from_lines(self, *, show_warning: bool) -> bool:
        if getattr(self, "_mass_axis_lock_level", 0) >= 1:
            return False
        pairs: List[Tuple[float, float]] = []
        for line in self._mass_lines:
            mu = float(line.mu)
            mass = float(line.mass_guess)
            if math.isfinite(mu) and math.isfinite(mass):
                pairs.append((mu, mass))
        if len(pairs) < 2:
            if show_warning:
                QMessageBox.information(
                    self,
                    "Insufficient Mass Lines",
                    "Add at least two mass lines with valid μ and mass values to estimate the axis.",
                )
            return False
        previous = dict(self._mass_params)
        self._estimate_mass_axis(pairs)
        self._apply_mass_params_update(previous)
        return True

    def _apply_mass_params_update(self, previous: Dict[str, float]) -> None:
        stretch = float(self._mass_params.get("stretch", 1.0))
        shift = float(self._mass_params.get("shift", 0.0))
        if hasattr(self, "mass_stretch_spin") and hasattr(self, "mass_shift_spin"):
            self.mass_stretch_spin.blockSignals(True)
            self.mass_shift_spin.blockSignals(True)
            self.mass_stretch_spin.setValue(stretch)
            self.mass_shift_spin.setValue(shift)
            self.mass_stretch_spin.blockSignals(False)
            self.mass_shift_spin.blockSignals(False)
        if (
            not math.isclose(previous.get("stretch", 1.0), stretch, rel_tol=1e-9, abs_tol=1e-9)
            or not math.isclose(previous.get("shift", 0.0), shift, rel_tol=1e-9, abs_tol=1e-9)
        ):
            self._on_mass_params_changed()
        else:
            self._refresh_plot()
            self._update_tables()
            self._update_summary()

    def _set_mass_assignment_index(self, index: int) -> None:
        if not hasattr(self, "mass_species_combo"):
            return
        if not (0 <= index < self.mass_species_combo.count()):
            index = 0
        self._block_species_assignment = True
        try:
            self.mass_species_combo.setCurrentIndex(index)
        finally:
            self._block_species_assignment = False

    def _on_mass_species_chosen(self, index: int) -> None:
        if self._block_species_assignment:
            return
        if not hasattr(self, "mass_species_combo"):
            return
        data = self.mass_species_combo.itemData(index)
        if data is None:
            return
        line = self._current_mass_line()
        if line is None:
            self._set_mass_assignment_index(0)
            return
        species_label, species_mass = data
        self._assign_species_to_line(line, str(species_label), float(species_mass))
        self._selected_line_id = line.line_id
        self._update_mass_line_abundances()
        self._update_tables()
        self._update_summary()
        self._refresh_plot()
        self._refresh_assignment_display()
        self._set_mass_assignment_index(0)

    def _assign_species_to_line(self, line: MassLineFit, species_label: str, species_mass: float) -> None:
        line.label = species_label
        line.mass_guess = species_mass
        line.assigned_species = species_label
        line.assigned_mass = species_mass
        self._update_mass_axis_lock()

    def _anchor_selected_line(self) -> None:
        line = self._current_mass_line()
        if line is None:
            QMessageBox.information(self, "No Selection", "Select a mass line to anchor.")
            return
        try:
            target_mass = float(self.anchor_mass_spin.value())
        except Exception:
            target_mass = float("nan")
        if not math.isfinite(target_mass):
            QMessageBox.warning(self, "Invalid Mass", "Enter a finite anchor mass value.")
            return
        mu = float(line.mu)
        if not math.isfinite(mu):
            QMessageBox.warning(
                self,
                "Invalid Mass Line",
                "The selected mass line does not have a valid centre time (μ).",
            )
            return
        stretch = float(self._mass_params.get("stretch", 1.0))
        if not math.isfinite(stretch) or abs(stretch) < 1.0e-12:
            QMessageBox.warning(
                self,
                "Invalid Stretch",
                "The mass-axis stretch must be finite and non-zero to anchor a line.",
            )
            return
        shift = mu - target_mass / stretch
        if not math.isfinite(shift):
            QMessageBox.warning(
                self,
                "Anchor Failed",
                "Unable to compute a valid shift from the selected line and anchor mass.",
            )
            return
        previous = dict(self._mass_params)
        self._mass_params["shift"] = shift
        self._anchor_displayed_line_id = line.line_id
        self._apply_mass_params_update(previous)

    def _refresh_assignment_display(self) -> None:
        line = self._current_mass_line()
        if hasattr(self, "mass_species_combo"):
            self.mass_species_combo.setEnabled(line is not None)
        if not hasattr(self, "mass_assignment_label"):
            return
        if line is None:
            self.mass_assignment_label.setText("No mass line selected.")
            return
        if line.assigned_species and math.isfinite(line.mass_guess):
            self.mass_assignment_label.setText(
                f"Assigned: {line.assigned_species} ({line.mass_guess:.3f} amu)"
            )
        elif math.isfinite(line.mass_guess):
            self.mass_assignment_label.setText(f"Estimated mass: {line.mass_guess:.3f} amu")
        else:
            self.mass_assignment_label.setText("Mass estimate unavailable.")

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
                new_label = text or line.label
                previous_assigned = line.assigned_species
                line.label = new_label
                if previous_assigned and new_label != previous_assigned:
                    line.assigned_species = None
                    line.assigned_mass = None
            elif column == 1:
                line.mass_guess = float(text)
                line.assigned_species = None
                line.assigned_mass = None
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
            self._recompute_mass_line(line, preserve_label=bool(line.assigned_species))
        if column == 1:
            line.label = nearest_mass_name(line.mass_guess)
        self._update_tables()
        self._update_summary()
        self._refresh_plot()
        self._update_mass_axis_lock()

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
        self._refresh_assignment_display()
        self._refresh_anchor_controls()

    def _on_sample_guess_changed(self, index: int) -> None:
        if self._block_sample_guess_signal:
            return
        data = self.sample_guess_combo.itemData(index)
        if isinstance(data, str) and data:
            self._manual_sample_guess = data
        elif index == 0:
            self._manual_sample_guess = None
        else:
            text = self.sample_guess_combo.currentText().strip()
            self._manual_sample_guess = text or None
        self._update_summary()

    def _on_sample_guess_text_changed(self, text: str) -> None:
        if self._block_sample_guess_signal:
            return
        cleaned = text.strip()
        if not cleaned or cleaned.lower() == "auto (closest match)":
            self._manual_sample_guess = None
        else:
            self._manual_sample_guess = cleaned
        self._update_summary()

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
            color="#ff7f0e",
            shape="emg",
        )
        line.fit_values = line.evaluate(dense_time)
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
        fit_curve = line.evaluate(dense_time)
        line.time_axis = dense_time
        line.fit_values = fit_curve
        preserve_label = preserve_label or bool(line.assigned_species)
        if line.assigned_mass is not None and math.isfinite(line.assigned_mass):
            line.mass_guess = float(line.assigned_mass)
        else:
            converted = self._time_to_mass(np.array([line.mu], dtype=float))
            if converted.size:
                line.mass_guess = float(converted[0])
            else:
                line.mass_guess = float("nan")
        new_auto = nearest_mass_name(line.mass_guess) if math.isfinite(line.mass_guess) else previous_auto
        if line.assigned_species:
            line.label = line.assigned_species
        elif not preserve_label:
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
            self._update_mass_axis_lock()
    def _update_mass_line_abundances(self) -> None:
        if self._combined is None or self._combined.size == 0:
            for line in self._mass_lines:
                line.abundance = 0.0
            self._recalculate_rsf_normalised()
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
        self._recalculate_rsf_normalised()

    def _recalculate_rsf_normalised(self) -> None:
        if not self._rsf_enabled:
            self._rsf_normalised = {}
            return
        valid_ids = {line.line_id for line in self._mass_lines}
        for line_id in list(self._rsf_values.keys()):
            if line_id not in valid_ids:
                del self._rsf_values[line_id]
        totals: Dict[int, float] = {}
        total_weight = 0.0
        for line in self._mass_lines:
            factor = self._rsf_values.get(line.line_id)
            if factor is None:
                continue
            base_fraction = max(line.abundance, 0.0)
            weight = base_fraction * max(float(factor), 0.0)
            totals[line.line_id] = weight
            total_weight += weight
        if not totals:
            self._rsf_enabled = False
            self._rsf_normalised = {}
            return
        if total_weight <= 0.0:
            self._rsf_normalised = {line_id: 0.0 for line_id in totals}
            return
        self._rsf_normalised = {line_id: weight / total_weight for line_id, weight in totals.items()}

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
        if self._rsf_enabled:
            self._recalculate_rsf_normalised()

        def _display_fraction(line: MassLineFit) -> float:
            if self._rsf_enabled and self._rsf_normalised:
                return self._rsf_normalised.get(line.line_id, line.abundance)
            return line.abundance

        entries = sorted(self._mass_lines, key=_display_fraction, reverse=True)
        self.summary_table.setRowCount(len(entries))
        for row, line in enumerate(entries):
            display_fraction = _display_fraction(line)
            for col, value in enumerate((line.label, f"{line.mass_guess:.3f}", f"{display_fraction * 100.0:.2f}")):
                item = self.summary_table.item(row, col)
                if item is None:
                    item = QTableWidgetItem()
                    self.summary_table.setItem(row, col, item)
                item.setText(str(value))
        if self._rsf_enabled and self._rsf_normalised:
            matching_lines = [replace(line, abundance=_display_fraction(line)) for line in self._mass_lines]
        else:
            matching_lines = list(self._mass_lines)
        best_match, mixture_match = analyse_sample_matches(matching_lines)
        self._auto_sample_match = best_match
        self._auto_mixture_match = mixture_match
        def _format_match(match: Optional[SampleMatch]) -> Optional[str]:
            if not match:
                return None
            return f"{match.sample.name} ({match.coverage * 100.0:.0f}% elemental match)"

        def _format_mixture(match: Optional[MixtureMatch]) -> Optional[str]:
            if not match:
                return None
            primary_frac = match.fractions[0] * 100.0
            secondary_frac = match.fractions[1] * 100.0
            parts = [
                f"{primary_frac:.0f}% {match.primary.name}",
                f"{secondary_frac:.0f}% {match.secondary.name}",
            ]
            description = " + ".join(parts)
            if match.coverage > 0.0:
                description += f" ({match.coverage * 100.0:.0f}% elemental match)"
            return description

        auto_text = _format_match(best_match)
        mixture_text = _format_mixture(mixture_match)
        if self._manual_sample_guess:
            lines = [f"Sample guess (manual): {self._manual_sample_guess}"]
            if auto_text:
                lines.append(f"Auto suggestion: {auto_text}")
            if mixture_text:
                lines.append(f"Mixture suggestion: {mixture_text}")
        else:
            if auto_text:
                lines = [f"Sample guess: {auto_text}"]
                if mixture_text:
                    lines.append(f"Mixture suggestion: {mixture_text}")
            else:
                lines = ["Sample guess: insufficient data"]
        self.sample_guess_label.setText("\n".join(lines))
        if hasattr(self, "rsf_status_label"):
            if self._rsf_enabled:
                if self._rsf_normalised:
                    active = [line_id for line_id, value in self._rsf_normalised.items() if value > 0.0]
                    selected = len(self._rsf_values)
                    if active:
                        message = (
                            f"Relative sensitivity factors applied to {selected} mass line(s); "
                            "selected lines are renormalised, others show measured abundances."
                        )
                    else:
                        message = (
                            "RSF weighting enabled — awaiting non-zero abundances to renormalise values."
                        )
                else:
                    message = "RSF weighting enabled — select lines and factors to apply."
                self.rsf_status_label.setText(message)
                self.rsf_status_label.show()
            else:
                self.rsf_status_label.hide()
        self._block_sample_guess_signal = True
        if self._manual_sample_guess:
            index = -1
            for row in range(1, self.sample_guess_combo.count()):
                data = self.sample_guess_combo.itemData(row)
                if isinstance(data, str) and data == self._manual_sample_guess:
                    index = row
                    break
            if index >= 0:
                self.sample_guess_combo.setCurrentIndex(index)
            else:
                self.sample_guess_combo.setEditText(self._manual_sample_guess)
        else:
            self.sample_guess_combo.setCurrentIndex(0)
        self._block_sample_guess_signal = False

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
                extras_dtype = h5py.string_dtype(encoding="utf-8", length=2048)
                table = np.zeros(len(self._mass_lines), dtype=[
                    ("id", "i4"),
                    ("label", str_dtype),
                    ("assigned_species", str_dtype),
                    ("mu", "f8"),
                    ("sigma", "f8"),
                    ("lam", "f8"),
                    ("amplitude", "f8"),
                    ("time_start", "f8"),
                    ("time_end", "f8"),
                    ("mass", "f8"),
                    ("assigned_mass", "f8"),
                    ("area", "f8"),
                    ("abundance", "f8"),
                    ("shape", str_dtype),
                    ("extras", extras_dtype),
                ])
                for idx, line in enumerate(self._mass_lines):
                    assigned_species = line.assigned_species or ""
                    if line.assigned_mass is not None and math.isfinite(line.assigned_mass):
                        assigned_mass = float(line.assigned_mass)
                    else:
                        assigned_mass = float("nan")
                    try:
                        extras_serialized = json.dumps(line.extra_params)
                    except Exception:
                        extras_serialized = "{}"
                    table[idx] = (
                        line.line_id,
                        line.label,
                        assigned_species,
                        line.mu,
                        line.sigma,
                        line.lam,
                        line.total_area(),
                        line.time_start,
                        line.time_end,
                        line.mass_guess,
                        assigned_mass,
                        line.window_area(),
                        line.abundance,
                        line.shape,
                        extras_serialized,
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


def launch_dust_composition_window(
    h5: h5py.File,
    event_name: str,
    *,
    event_names: Optional[Sequence[str]] = None,
    on_event_changed: Optional[Callable[[str], None]] = None,
    parent: Optional[QWidget] = None,
) -> DustCompositionWindow:
    """Convenience helper used by the main quicklook window."""

    return DustCompositionWindow(
        h5=h5,
        event_name=event_name,
        event_names=event_names,
        on_event_changed=on_event_changed,
        parent=parent,
    )
