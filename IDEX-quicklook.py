#!/opt/anaconda3/envs/idex/bin/python
# -*- coding: utf-8 -*-
from __future__ import annotations

# =====================================================================
# IDEX Quicklook — QtAgg canvas + Matplotlib Navigation Toolbar
# Interactive viewer with channel selection, overlays, and fit editing tools
# =====================================================================

import os
import sys
import argparse
import tempfile
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import h5py
import numpy as np

from HDF_View import launch_hdf_viewer
from dust_composition import launch_dust_composition_window

# Qt-backed Matplotlib so NavigationToolbar works
import matplotlib
from matplotlib import mathtext
matplotlib.use("QtAgg")

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

# --------- Qt binding-agnostic imports (prefer PySide6, fallback PyQt6) ---------
_QT = None
try:
    from PySide6.QtCore import Qt, QSize
    from PySide6.QtGui import QAction, QFont, QPixmap, QImage
    from PySide6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy, QDialog, QPushButton,
        QHBoxLayout, QGridLayout, QTableWidget, QTableWidgetItem, QHeaderView,
        QCheckBox, QDialogButtonBox
    )
    _QT = "PySide6"
except Exception:
    from PyQt6.QtCore import Qt, QSize
    from PyQt6.QtGui import QAction, QFont, QPixmap, QImage
    from PyQt6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy, QDialog, QPushButton,
        QHBoxLayout, QGridLayout, QTableWidget, QTableWidgetItem, QHeaderView,
        QCheckBox, QDialogButtonBox
    )
    _QT = "PyQt6"

print(f"[info] Qt binding: {_QT}, Matplotlib backend: {matplotlib.get_backend()}")

matplotlib.rcParams.update({
    "figure.autolayout": False,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "lines.linewidth": 1.8,
})

# --------- Small utils ---------
def non_native_open_dialog(parent: QWidget, start_dir: str | None = None) -> str | None:
    if start_dir is None:
        start_dir = os.path.join(os.getcwd(), "HDF5")
    options = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.ReadOnly
    filename, _ = QFileDialog.getOpenFileName(
        parent, "Open HDF5", start_dir,
        "HDF5 Files (*.h5 *.hdf5);;All files (*)", options=options
    )
    return filename or None

def list_event_groups(h5: h5py.File) -> list[str]:
    ev = [k for k, v in h5.items() if isinstance(v, h5py.Group)]
    try:
        return sorted(ev, key=lambda s: int(str(s)))
    except Exception:
        return sorted(ev)

def dset(h5: h5py.File, path: str):
    try:
        obj = h5[path]
        return obj[()] if isinstance(obj, h5py.Dataset) else None
    except Exception:
        return None

def _normalize_name(value: str) -> str:
    return "".join(ch.lower() for ch in value if ch.isalnum())

def _pair_key(path: str) -> str:
    base = _normalize_name(os.path.basename(path))
    for token in ("time", "result", "fit", "curve"):
        base = base.replace(token, "")
    return base

def _friendly_label(path: str) -> str:
    base = os.path.basename(path)
    replacements = {
        "FitResult": "Fit",
        "fitresult": "Fit",
        "Result": "Fit",
    }
    for needle, repl in replacements.items():
        base = base.replace(needle, repl)
    return base


def _label_from_param_path(path: str) -> str:
    base = os.path.basename(path)
    replacements = (
        ("FitParameters", "Fit"),
        ("FitParams", "Fit"),
        ("fitparameters", "Fit"),
        ("fitparams", "Fit"),
        ("Parameters", "Fit"),
        ("Params", "Fit"),
    )
    for needle, repl in replacements:
        if needle in base:
            base = base.replace(needle, repl)
    if "Fit" in base and " Fit" not in base:
        base = base.replace("Fit", " Fit")
    base = " ".join(base.split()).strip()
    if not base:
        base = "Fit"
    if "Fit" not in base:
        base = f"{base} Fit"
    return f"{base} (calc)"


_MATH_TEXT_PARSER = mathtext.MathTextParser("Macosx")
_LATEX_CACHE: Dict[str, Optional[QPixmap]] = {}


def _latex_to_pixmap(latex: str) -> Optional[QPixmap]:
    if not latex:
        return None
    cached = _LATEX_CACHE.get(latex)
    if cached is not None:
        return cached
    try:
        ftimage, _ = _MATH_TEXT_PARSER.to_rgba(f"${latex}$", dpi=180)
    except Exception:
        _LATEX_CACHE[latex] = None
        return None

    buffer = np.ascontiguousarray(np.clip(ftimage * 255.0, 0, 255).astype(np.uint8))
    height, width, _channels = buffer.shape
    image = QImage(buffer.data, width, height, 4 * width, QImage.Format_RGBA8888)
    pixmap = QPixmap.fromImage(image.copy())
    _LATEX_CACHE[latex] = pixmap
    return pixmap


def _mass_identifier_from_path(path: str) -> Optional[str]:
    if not path or "/Masses/" not in path:
        return None
    tail = path.split("/Masses/", 1)[-1]
    segment = tail.split("/", 1)[0]
    for suffix in ("FitParameters", "FitParams", "Parameters", "Params"):
        if segment.endswith(suffix):
            segment = segment[: -len(suffix)]
            break
    segment = segment.strip().strip("_-")
    return segment or None


def _dataset_sort_key(channel: str, dataset_path: str) -> Tuple[int, object]:
    mass_label = _mass_identifier_from_path(dataset_path)
    if mass_label is not None:
        try:
            mass_value: object = float(mass_label)
        except ValueError:
            mass_value = mass_label
        return (0, mass_value)
    return (1, os.path.basename(dataset_path))


def _coerce_parameter_values(values: np.ndarray) -> Optional[np.ndarray]:
    try:
        arr = np.asarray(values)
    except Exception:
        return None

    if arr.size == 0:
        return None

    if arr.dtype.names:  # structured array
        for field in ("value", "values", "val"):
            if field in arr.dtype.names:
                arr = np.asarray(arr[field])
                break
        else:
            pieces = []
            for name in arr.dtype.names:
                try:
                    pieces.append(np.asarray(arr[name]).ravel())
                except Exception:
                    return None
            if not pieces:
                return None
            arr = np.concatenate(pieces)
    elif arr.dtype == object:
        coerced: List[float] = []
        for item in arr.ravel():
            if item is None:
                coerced.append(np.nan)
                continue
            if isinstance(item, (float, int, np.floating, np.integer)):
                coerced.append(float(item))
                continue
            if hasattr(item, "value"):
                try:
                    coerced.append(float(item.value))
                    continue
                except Exception:
                    pass
            if isinstance(item, (list, tuple, np.ndarray)):
                sub = np.asarray(item).ravel()
                if sub.size:
                    try:
                        coerced.append(float(sub[0]))
                        continue
                    except Exception:
                        pass
            try:
                coerced.append(float(item))
            except Exception:
                return None
        arr = np.asarray(coerced, dtype=float)
    else:
        arr = np.asarray(arr, dtype=float)

    arr = np.asarray(arr, dtype=float).ravel()
    if arr.size == 0:
        return None
    return arr

def _to_1d(array: np.ndarray) -> np.ndarray:
    arr = np.asarray(array)
    if arr.ndim == 0:
        return arr.reshape(1)
    return arr.ravel()

# --------- Channel & fit metadata ---------
FAMILY_HIGH = "high"
FAMILY_LOW = "low"


@dataclass(frozen=True)
class ChannelDefinition:
    dataset: str
    time_dataset: str
    family: str


CHANNEL_DEFS: Dict[str, ChannelDefinition] = {
    "TOF L": ChannelDefinition(dataset="TOF L", time_dataset="Time (high sampling)", family=FAMILY_HIGH),
    "TOF M": ChannelDefinition(dataset="TOF M", time_dataset="Time (high sampling)", family=FAMILY_HIGH),
    "TOF H": ChannelDefinition(dataset="TOF H", time_dataset="Time (high sampling)", family=FAMILY_HIGH),
    "Ion Grid": ChannelDefinition(dataset="Ion Grid", time_dataset="Time (low sampling)", family=FAMILY_LOW),
    "Target L": ChannelDefinition(dataset="Target L", time_dataset="Time (low sampling)", family=FAMILY_LOW),
    "Target H": ChannelDefinition(dataset="Target H", time_dataset="Time (low sampling)", family=FAMILY_LOW),
}

CHANNEL_ORDER: List[str] = ["TOF L", "TOF M", "TOF H", "Ion Grid", "Target L", "Target H"]

FIT_ELIGIBLE_CHANNELS = {"Ion Grid", "Target L", "Target H", "TOF H"}

BASELINE_PRIMARY_WINDOW = (-7.0, -5.0)
BASELINE_FALLBACK_THRESHOLD = -2.0

FAMILY_TITLES = {
    FAMILY_HIGH: "TOF Channels (High Sampling)",
    FAMILY_LOW: "Ion Grid / Target Channels (Low Sampling)",
}

FAMILY_YLABELS = {
    FAMILY_HIGH: r"TOF Channels [pC/$\Delta t$]",
    FAMILY_LOW: r"Ion Grid / Target Channels [pC]",
}

TIME_AXIS_LABEL = r"Time [$\mu$s]"


FIT_SCALE_MULTIPLIERS: Dict[str, float] = {
    "Ion Grid": 7.46e-4,
    "Target H": 1.63e-1,
    "Target L": 1.58e1,
}


@dataclass(frozen=True)
class FitModelMeta:
    parameter_labels: Tuple[str, ...]
    latex: str
    evaluator: Callable[..., np.ndarray]


def _idex_ion_grid_model(x: np.ndarray, p0: float, p1: float, p4: float, p5: float, p6: float) -> np.ndarray:
    """Evaluate the standard IDEX ion grid / target response model."""

    arr = np.asarray(x, dtype=float)
    shift = arr - p0
    step = np.heaviside(shift, 0.0)
    safe_p5 = p5 if abs(p5) > 1.0e-12 else 1.0e-12
    safe_p6 = p6 if abs(p6) > 1.0e-12 else 1.0e-12
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        rise = 1.0 - np.exp(-shift / safe_p5)
        decay = np.exp(-shift / safe_p6)
    return p1 + step * (p4 * rise * decay)


def _emg_model(x: np.ndarray, mu: float, sigma: float, lam: float) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.size == 0:
        return arr
    if abs(lam) < 1.0e-15:
        return np.zeros_like(arr, dtype=float)
    safe_sigma = sigma if abs(sigma) > 1.0e-15 else 1.0e-15
    safe_lambda = lam
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        exponent = np.exp((safe_lambda / 2.0) * (2.0 * mu + safe_lambda * safe_sigma**2 - 2.0 * arr))
    argument = (mu + safe_lambda * safe_sigma**2 - arr) / (np.sqrt(2.0) * safe_sigma)
    return (safe_lambda / 2.0) * exponent * np.erfc(argument)


ION_GRID_FIT = FitModelMeta(
    parameter_labels=(
        "t0 (impact time)",
        "c (constant offset)",
        "A (amplitude)",
        "t1 (rise time)",
        "t2 (discharge time)",
    ),
    latex=r"f(t) = c + H(t - t_0)\,A\,\left(1 - e^{-(t - t_0)/t_1}\right)e^{-(t - t_0)/t_2}",
    evaluator=_idex_ion_grid_model,
)

EMG_FIT = FitModelMeta(
    parameter_labels=(
        "μ (location)",
        "σ (width)",
        "λ (decay)",
    ),
    latex=r"f(t) = \frac{\lambda}{2} \exp\left[\frac{\lambda}{2}(2\mu + \lambda\sigma^2 - 2t)\right] \operatorname{erfc}\left(\frac{\mu + \lambda\sigma^2 - t}{\sqrt{2}\,\sigma}\right)",
    evaluator=_emg_model,
)

FIT_MODEL_BY_CHANNEL: Dict[str, FitModelMeta] = {
    "Ion Grid": ION_GRID_FIT,
    "Target L": ION_GRID_FIT,
    "Target H": ION_GRID_FIT,
    "TOF H": EMG_FIT,
}


def _evaluate_fit_curve(channel: str, time_values: np.ndarray, params: np.ndarray) -> Optional[np.ndarray]:
    """Evaluate a fit curve for the supplied channel using the stored parameters."""

    model = FIT_MODEL_BY_CHANNEL.get(channel)
    if not model:
        return None

    param_array = _coerce_parameter_values(params)
    if param_array is None or param_array.size < len(model.parameter_labels):
        return None

    first = np.asarray(param_array[: len(model.parameter_labels)], dtype=float)
    if not np.all(np.isfinite(first)):
        return None

    try:
        return model.evaluator(time_values, *first)
    except Exception:
        return None


def _masked_mean(values: np.ndarray, mask: np.ndarray) -> Optional[float]:
    values = np.asarray(values)
    mask = np.asarray(mask, dtype=bool)
    if values.size == 0 or mask.size != values.size:
        return None
    if not np.any(mask):
        return None
    subset = values[mask]
    subset = subset[np.isfinite(subset)]
    if subset.size == 0:
        return None
    return float(subset.mean())


def _fit_paths_from_param(dataset_path: str) -> Optional[Tuple[str, str]]:
    """Return the fit result/time dataset paths corresponding to a parameter dataset."""

    if not dataset_path:
        return None

    for needle in ("FitParams", "FitParameters", "FitParam"):
        if needle in dataset_path:
            result_path = dataset_path.replace(needle, "FitResult")
            time_path = dataset_path.replace(needle, "FitTime")
            return result_path, time_path

    return None


@dataclass
class FitData:
    time_series: Dict[str, np.ndarray] = field(default_factory=dict)
    value_series: Dict[str, np.ndarray] = field(default_factory=dict)
    parameter_series: Dict[str, np.ndarray] = field(default_factory=dict)
    extras: Dict[str, np.ndarray] = field(default_factory=dict)

    def has_overlay(self) -> bool:
        if any(True for _ in self.iter_time_result_pairs()):
            return True
        return bool(self.parameter_series)

    def has_parameters(self) -> bool:
        return bool(self.parameter_series)

    def iter_time_result_pairs(self) -> Iterable[Tuple[str, str, str, np.ndarray, np.ndarray]]:
        seen_keys = set()
        for time_path, time_values in self.time_series.items():
            pair_key = _pair_key(time_path)
            if not pair_key:
                continue
            for value_path, value_values in self.value_series.items():
                if _pair_key(value_path) != pair_key:
                    continue
                dedupe_key = (time_path, value_path)
                if dedupe_key in seen_keys:
                    continue
                seen_keys.add(dedupe_key)
                yield (
                    _friendly_label(value_path),
                    time_path,
                    value_path,
                    _to_1d(time_values),
                    _to_1d(value_values),
                )

    def iter_parameter_items(self) -> Iterable[Tuple[str, np.ndarray]]:
        for path, values in self.parameter_series.items():
            yield path, np.asarray(values)


def y_label_with_units(channel_name: str) -> str:
    if channel_name.startswith("TOF"):
        return rf"{channel_name} [pC/$\Delta t$]"
    return rf"{channel_name} [pC]"

# --------- Main Window ---------
class MainWindow(QMainWindow):
    def __init__(self, filename: str | None = None, eventnumber: int | None = None):
        super().__init__()
        self.setWindowTitle("IDEX Quicklook — Interactive Viewer")
        self.setMinimumSize(1250, 820)
        self.setStatusBar(QStatusBar(self))

        self._h5: Optional[h5py.File] = None
        self._events: List[str] = []
        self._current_event: Optional[str] = None
        self._filename: Optional[str] = None
        self._tmpdir = tempfile.TemporaryDirectory(prefix="idex_quicklook_")
        self._fit_cache: Dict[Tuple[str, str], FitData] = {}
        self._fit_param_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._fit_result_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._fit_time_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._baseline_cache: Dict[Tuple[str, str], float] = {}
        self._show_fit: Dict[str, bool] = {name: False for name in FIT_ELIGIBLE_CHANNELS}
        self.selected_channels = set(CHANNEL_ORDER)
        self._child_windows: List[QWidget] = []

        central = QWidget(self)
        self.vbox = QVBoxLayout(central)
        self.vbox.setContentsMargins(10, 10, 10, 10)
        self.vbox.setSpacing(12)
        self.setCentralWidget(central)

        self._build_toolbar()
        self._build_controls()

        self.figure = Figure(figsize=(12, 8), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.nav_toolbar = NavigationToolbar(self.canvas, self)
        self.nav_toolbar.setStyleSheet("font-size: 14px; padding: 6px;")
        self.vbox.addWidget(self.nav_toolbar)
        self.vbox.addWidget(self.canvas)

        self.statusBar().showMessage("Select an HDF5 file to get started.")

        if filename:
            self.open_file(filename)
        else:
            chosen = non_native_open_dialog(self)
            if not chosen:
                self.close()
                return
            self.open_file(chosen)

        if eventnumber is not None and self._events and 1 <= eventnumber <= len(self._events):
            self.event_combo.setCurrentIndex(eventnumber - 1)

    # ---- UI construction -------------------------------------------------
    def _build_toolbar(self):
        tb = QToolBar("Main", self)
        tb.setIconSize(QSize(22, 22))
        tb.setMovable(False)
        self.addToolBar(tb)

        act_open = QAction("Open HDF5…", self)
        act_open.setShortcut("Ctrl+O")
        act_open.triggered.connect(self.action_open)
        tb.addAction(act_open)

        act_view = QAction("View HDF Contents", self)
        act_view.setShortcut("Ctrl+Shift+H")
        act_view.setToolTip("Inspect the currently loaded HDF5 file in a tree viewer.")
        act_view.triggered.connect(self.action_view_hdf_contents)
        tb.addAction(act_view)

        act_dust = QAction("Dust Composition…", self)
        act_dust.setShortcut("Ctrl+D")
        act_dust.setToolTip("Open the dust composition analysis window for the current event.")
        act_dust.triggered.connect(self.action_open_dust_composition)
        tb.addAction(act_dust)

        act_reload = QAction("Reload", self)
        act_reload.setShortcut("Ctrl+R")
        act_reload.triggered.connect(self.reload_current)
        tb.addAction(act_reload)

        tb.addSeparator()

        act_quit = QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        tb.addAction(act_quit)

        tb.addSeparator()
        label = QLabel("Event:", self)
        label.setStyleSheet("font-size: 15px; font-weight: bold; padding-right: 6px;")
        tb.addWidget(label)

        self.event_combo = QComboBox(self)
        self.event_combo.setMinimumWidth(220)
        self.event_combo.setStyleSheet("font-size: 15px; min-height: 36px;")
        self.event_combo.currentIndexChanged.connect(self.on_event_changed)
        tb.addWidget(self.event_combo)

    def _build_controls(self):
        panel = QWidget(self)
        panel_layout = QVBoxLayout(panel)
        panel_layout.setContentsMargins(0, 0, 0, 0)
        panel_layout.setSpacing(8)
        panel_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        heading = QLabel("Display Controls", self)
        heading.setStyleSheet("font-size: 18px; font-weight: bold;")
        panel_layout.addWidget(heading)

        sub_label = QLabel("Choose which channels and overlays are shown:", self)
        sub_label.setStyleSheet("font-size: 15px;")
        panel_layout.addWidget(sub_label)

        channel_widget = QWidget(self)
        grid = QGridLayout(channel_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(10)
        self.channel_buttons: Dict[str, QPushButton] = {}

        for idx, name in enumerate(CHANNEL_ORDER):
            btn = QPushButton(name, self)
            btn.setCheckable(True)
            btn.setChecked(True)
            btn.setMinimumHeight(50)
            btn.setStyleSheet("font-size: 16px; font-weight: bold; padding: 10px 18px;")
            btn.clicked.connect(lambda checked, channel=name: self.on_channel_toggled(channel, checked))
            self.channel_buttons[name] = btn
            grid.addWidget(btn, idx // 3, idx % 3)

        panel_layout.addWidget(channel_widget)

        toggle_row = QHBoxLayout()
        toggle_row.setSpacing(10)

        self.overlay_button = QPushButton("Overlay same time axis", self)
        self.overlay_button.setCheckable(True)
        self.overlay_button.setMinimumHeight(50)
        self.overlay_button.setStyleSheet("font-size: 16px; font-weight: bold; padding: 10px 18px;")
        self.overlay_button.setToolTip("When enabled, channels with the same time base are drawn together.")
        self.overlay_button.clicked.connect(self.on_overlay_toggled)
        toggle_row.addWidget(self.overlay_button)

        self.fit_buttons: Dict[str, QPushButton] = {}
        for channel in sorted(FIT_ELIGIBLE_CHANNELS):
            btn = QPushButton(f"Show {channel} Fit", self)
            btn.setCheckable(True)
            btn.setMinimumHeight(50)
            btn.setStyleSheet("font-size: 16px; padding: 10px 18px;")
            btn.setToolTip("Overlay fit curves when available.")
            btn.clicked.connect(lambda checked, chan=channel: self.on_fit_toggled(chan, checked))
            toggle_row.addWidget(btn)
            self.fit_buttons[channel] = btn

        self.edit_params_button = QPushButton("Edit Fit Parameters", self)
        self.edit_params_button.setMinimumHeight(50)
        self.edit_params_button.setStyleSheet("font-size: 16px; font-weight: bold; padding: 10px 18px;")
        self.edit_params_button.clicked.connect(self.open_fit_parameter_dialog)
        toggle_row.addWidget(self.edit_params_button)

        toggle_row.addStretch(1)
        panel_layout.addLayout(toggle_row)

        self.vbox.addWidget(panel)

    # ---- Actions ---------------------------------------------------------
    def action_open(self):
        chosen = non_native_open_dialog(self)
        if chosen:
            self.open_file(chosen)

    def action_view_hdf_contents(self):
        if not self._filename:
            QMessageBox.information(
                self,
                "No File Loaded",
                "Open an HDF5 file to browse its contents.",
            )
            return

        try:
            viewer = launch_hdf_viewer(self._filename, parent=self)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Viewer Error",
                f"Unable to launch the HDF viewer:\n{exc}",
            )
            return

        viewer.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        viewer.raise_()
        viewer.activateWindow()

        self._child_windows.append(viewer)

        def _cleanup(*_args):
            if viewer in self._child_windows:
                self._child_windows.remove(viewer)

        viewer.destroyed.connect(_cleanup)

    def action_open_dust_composition(self):
        if not self._h5 or not self._current_event:
            QMessageBox.information(
                self,
                "No Event",
                "Open an HDF5 file and select an event before analysing dust composition.",
            )
            return

        try:
            window = launch_dust_composition_window(self._h5, self._current_event, parent=self)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Dust Composition Error",
                f"Unable to launch the dust composition window:\n{exc}",
            )
            return

        window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        window.show()

        self._child_windows.append(window)

        def _cleanup(*_args):
            if window in self._child_windows:
                self._child_windows.remove(window)

        window.destroyed.connect(_cleanup)

    def open_file(self, path: str, preferred_event: Optional[str] = None):
        try:
            if self._h5 is not None:
                self._h5.close()
        except Exception:
            pass

        self._h5 = None
        self._events = []
        self._current_event = None
        self._filename = None
        self._fit_cache.clear()
        self._fit_param_overrides.clear()
        self._fit_result_overrides.clear()
        self._fit_time_overrides.clear()
        self._baseline_cache.clear()

        try:
            try:
                self._h5 = h5py.File(path, "r+")
            except OSError:
                self._h5 = h5py.File(path, "r")
            self._filename = path
        except Exception as exc:
            QMessageBox.critical(self, "Open Error", f"Failed to open file:\n{path}\n\n{exc}")
            self.plot_event(None)
            return

        events = list_event_groups(self._h5)
        if not events:
            QMessageBox.warning(self, "No Events", "No top-level event groups found in this file.")
        self._events = events

        self.event_combo.blockSignals(True)
        self.event_combo.clear()
        self.event_combo.addItems(self._events)
        self.event_combo.blockSignals(False)

        target_event = None
        if preferred_event and preferred_event in self._events:
            target_event = preferred_event
        elif self._events:
            target_event = self._events[0]

        if target_event:
            idx = self._events.index(target_event)
            self.event_combo.blockSignals(True)
            self.event_combo.setCurrentIndex(idx)
            self.event_combo.blockSignals(False)
            self.plot_event(target_event)
        else:
            self.plot_event(None)

    def reload_current(self):
        if self._filename:
            self.open_file(self._filename, preferred_event=self._current_event)

    def on_event_changed(self, idx: int):
        if 0 <= idx < len(self._events):
            self.plot_event(self._events[idx])

    def on_channel_toggled(self, channel: str, checked: bool):
        if checked:
            self.selected_channels.add(channel)
        else:
            self.selected_channels.discard(channel)
        self.plot_event(self._current_event)

    def on_overlay_toggled(self, checked: bool):
        self.plot_event(self._current_event)

    def on_fit_toggled(self, channel: str, checked: bool):
        if not self._current_event or not self._h5:
            self._show_fit[channel] = False
            if channel in self.fit_buttons:
                self.fit_buttons[channel].setChecked(False)
            return

        data = self.get_fit_data(self._current_event, channel)
        if checked and not data.has_overlay():
            QMessageBox.information(self, "No Fit", f"No fit curves were found for {channel} in this event.")
            self._show_fit[channel] = False
            if channel in self.fit_buttons:
                self.fit_buttons[channel].setChecked(False)
            return

        self._show_fit[channel] = checked
        self.plot_event(self._current_event)

    # ---- Plotting --------------------------------------------------------
    def plot_event(self, event_name: Optional[str]):
        self.figure.clear()
        self.figure.set_constrained_layout(True)

        if not self._h5 or not event_name:
            self._current_event = None
            ax = self.figure.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "Open an HDF5 file and choose an event to visualize.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=18,
            )
            ax.axis("off")
            self.canvas.draw_idle()
            self.refresh_fit_controls()
            self.update_status_text()
            return

        self._current_event = event_name
        self.refresh_fit_controls()

        selected = [name for name in CHANNEL_ORDER if self.channel_buttons.get(name) and self.channel_buttons[name].isChecked()]
        missing: List[str] = []

        self.figure.suptitle(
            f"{os.path.basename(self._filename or '')} — Event {event_name}",
            fontsize=20,
            fontweight="bold",
        )

        if not selected:
            ax = self.figure.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "No channels selected.\nUse the buttons above to enable channels.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=18,
            )
            ax.axis("off")
            self.canvas.draw_idle()
            self.update_status_text()
            return

        overlay_mode = self.overlay_button.isChecked()

        try:
            if overlay_mode:
                families: List[Tuple[str, List[str]]] = []
                high = [ch for ch in selected if CHANNEL_DEFS[ch].family == FAMILY_HIGH]
                low = [ch for ch in selected if CHANNEL_DEFS[ch].family == FAMILY_LOW]
                if high:
                    families.append((FAMILY_HIGH, high))
                if low:
                    families.append((FAMILY_LOW, low))

                if not families:
                    ax = self.figure.add_subplot(111)
                    ax.text(
                        0.5,
                        0.5,
                        "No compatible channels selected for overlay.",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=18,
                    )
                    ax.axis("off")
                else:
                    for idx, (family, channels) in enumerate(families, start=1):
                        ax = self.figure.add_subplot(len(families), 1, idx)
                        plotted_any = False
                        for channel in channels:
                            plotted_any |= self._plot_channel(ax, event_name, channel, overlay_mode=True, missing_channels=missing)
                        self._style_overlay_axis(ax, family, bottom=(idx == len(families)))
                        if plotted_any and (len(channels) > 1 or any(self._show_fit.get(ch, False) for ch in channels)):
                            ax.legend(loc="best")
            else:
                ordered = [ch for ch in CHANNEL_ORDER if ch in selected]
                if not ordered:
                    ax = self.figure.add_subplot(111)
                    ax.text(
                        0.5,
                        0.5,
                        "No channels available to plot.",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=18,
                    )
                    ax.axis("off")
                else:
                    grid = self.figure.add_gridspec(len(ordered), 1, hspace=0.42)
                    for idx, channel in enumerate(ordered):
                        ax = self.figure.add_subplot(grid[idx, 0])
                        plotted_any = self._plot_channel(ax, event_name, channel, overlay_mode=False, missing_channels=missing)
                        self._style_single_axis(ax, bottom=(idx == len(ordered) - 1))
                        if plotted_any and self._show_fit.get(channel) and len(ax.lines) > 1:
                            ax.legend(loc="best")
        except Exception as exc:
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Plot error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=16)
            ax.axis("off")

        self.canvas.draw_idle()
        self.update_status_text(missing)

    def _plot_channel(self, ax, event_name: str, channel: str, overlay_mode: bool, missing_channels: List[str]) -> bool:
        definition = CHANNEL_DEFS[channel]
        time_data = dset(self._h5, f"/{event_name}/{definition.time_dataset}")
        value_data = dset(self._h5, f"/{event_name}/{definition.dataset}")

        base_plotted = False
        reason: Optional[str] = None
        if time_data is None or value_data is None:
            reason = "Missing data"
        else:
            t = _to_1d(time_data)
            y = _to_1d(value_data)
            n = min(len(t), len(y))
            if n == 0:
                reason = "Empty dataset"
            else:
                label = channel if overlay_mode else None
                ax.plot(t[:n], y[:n], label=label)
                base_plotted = True

        fit_plotted = self._plot_fit(ax, event_name, channel, overlay_mode)

        if base_plotted or fit_plotted:
            if not overlay_mode:
                ax.set_title(channel, fontsize=17, fontweight="bold")
                ax.set_ylabel(y_label_with_units(channel), fontsize=16)
            return True

        if reason:
            if not overlay_mode:
                self._draw_missing_message(ax, channel, reason)
            missing_channels.append(channel)
        return False

    def _estimate_baseline(self, event_name: str, channel: str, reference_time: Optional[np.ndarray]) -> float:
        key = (event_name, channel)
        cached = self._baseline_cache.get(key)
        if cached is not None:
            return cached

        baseline_value = 0.0
        if not self._h5:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        definition = CHANNEL_DEFS.get(channel)
        if not definition:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        raw_path = f"/{event_name}/{definition.dataset}"
        raw_data = dset(self._h5, raw_path)
        if raw_data is None:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        raw_values = _to_1d(raw_data)
        if raw_values.size == 0:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        raw_values = np.asarray(raw_values, dtype=float)

        times: Optional[np.ndarray] = None
        if reference_time is not None:
            ref = _to_1d(reference_time)
            if ref.size == raw_values.size:
                times = np.asarray(ref, dtype=float)

        if times is None:
            time_path = f"/{event_name}/{definition.time_dataset}"
            stored_time = dset(self._h5, time_path)
            if stored_time is not None:
                stored = _to_1d(stored_time)
                if stored.size == raw_values.size:
                    times = np.asarray(stored, dtype=float)

        candidate: Optional[float] = None
        if times is not None and times.size == raw_values.size:
            primary_mask = (times >= BASELINE_PRIMARY_WINDOW[0]) & (times <= BASELINE_PRIMARY_WINDOW[1])
            candidate = _masked_mean(raw_values, primary_mask)
            if candidate is None:
                fallback_mask = times < BASELINE_FALLBACK_THRESHOLD
                candidate = _masked_mean(raw_values, fallback_mask)

        if candidate is None:
            count = min(max(1, raw_values.size // 10), raw_values.size)
            subset = raw_values[:count]
            finite = subset[np.isfinite(subset)]
            if finite.size:
                candidate = float(finite.mean())
            else:
                candidate = 0.0

        baseline_value = float(candidate if np.isfinite(candidate) else 0.0)
        self._baseline_cache[key] = baseline_value
        return baseline_value

    def _iter_fit_curves(self, event_name: str, channel: str, data: FitData):
        yielded = False
        for label, _time_path, _value_path, time_values, fit_values in data.iter_time_result_pairs():
            times = np.asarray(time_values, dtype=float)
            values = np.asarray(fit_values, dtype=float)
            baseline_offset = self._estimate_baseline(event_name, channel, times)
            if baseline_offset:
                values = values + baseline_offset
            n = min(len(times), len(values))
            if n == 0:
                continue
            yield label, times[:n], values[:n]
            yielded = True

        if yielded or not data.parameter_series or not self._h5:
            return

        definition = CHANNEL_DEFS.get(channel)
        if not definition:
            return
        base_time_data = dset(self._h5, f"/{event_name}/{definition.time_dataset}")
        if base_time_data is None:
            return
        base_time = np.asarray(_to_1d(base_time_data), dtype=float)
        if base_time.size == 0:
            return

        model = FIT_MODEL_BY_CHANNEL.get(channel)
        if not model:
            return

        scale = FIT_SCALE_MULTIPLIERS.get(channel, 1.0)
        for path, raw_params in data.iter_parameter_items():
            derived = _fit_paths_from_param(path)
            if derived:
                result_path, time_path = derived
                if result_path in data.value_series and time_path in data.time_series:
                    continue

            param_array = _coerce_parameter_values(raw_params)
            if param_array is None or param_array.size < len(model.parameter_labels):
                continue

            first = np.asarray(param_array[: len(model.parameter_labels)], dtype=float)
            if not np.all(np.isfinite(first)):
                continue

            try:
                curve = model.evaluator(base_time, *first)
            except Exception:
                continue

            fit_values = np.asarray(curve, dtype=float)
            if fit_values.size != base_time.size:
                continue
            if scale != 1.0:
                fit_values = fit_values * scale
            baseline_offset = self._estimate_baseline(event_name, channel, base_time)
            if baseline_offset:
                fit_values = fit_values + baseline_offset

            label = _label_from_param_path(path)
            yield label, base_time, fit_values

    def _plot_fit(self, ax, event_name: str, channel: str, overlay_mode: bool) -> bool:
        if not self._show_fit.get(channel):
            return False

        data = self.get_fit_data(event_name, channel)
        plotted = False
        for label, time_values, fit_values in self._iter_fit_curves(event_name, channel, data):
            n = min(len(time_values), len(fit_values))
            if n == 0:
                continue
            legend_label = label
            if overlay_mode and not label.lower().startswith(channel.lower()):
                legend_label = f"{channel} {label}"
            ax.plot(time_values[:n], fit_values[:n], linestyle="--", linewidth=2.2, label=legend_label)
            plotted = True
        return plotted

    def _style_overlay_axis(self, ax, family: str, bottom: bool):
        ax.set_facecolor("#f8f9fb")
        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=14, width=1.5, length=7)
        ax.set_title(FAMILY_TITLES.get(family, "Channels"), fontsize=18, fontweight="bold")
        ax.set_ylabel(FAMILY_YLABELS.get(family, "Channel"), fontsize=16)
        if bottom:
            ax.set_xlabel(TIME_AXIS_LABEL, fontsize=16)
        else:
            ax.set_xlabel("")

    def _style_single_axis(self, ax, bottom: bool):
        ax.set_facecolor("#f8f9fb")
        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=14, width=1.5, length=7)
        if bottom:
            ax.set_xlabel(TIME_AXIS_LABEL, fontsize=16)
        else:
            ax.set_xlabel("")

    def _draw_missing_message(self, ax, channel: str, reason: str):
        ax.set_facecolor("#f8f9fb")
        ax.text(0.5, 0.5, f"{channel}\n{reason}", ha="center", va="center", transform=ax.transAxes, fontsize=16)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    def refresh_fit_controls(self):
        if not self._h5 or not self._current_event:
            for btn in self.fit_buttons.values():
                btn.setEnabled(False)
            return

        for channel, btn in self.fit_buttons.items():
            data = self.get_fit_data(self._current_event, channel)
            available = any(True for _ in self._iter_fit_curves(self._current_event, channel, data))
            btn.setEnabled(available)
            if not available:
                btn.setChecked(False)
                self._show_fit[channel] = False
                btn.setToolTip("No fit curve found for this event.")
            else:
                btn.setToolTip("Overlay the available fit curve on the channel plot.")

    def update_status_text(self, missing: Optional[List[str]] = None):
        parts: List[str] = []
        event_label = self._current_event or "–"
        parts.append(f"Event {event_label}")

        channels = [ch for ch in CHANNEL_ORDER if ch in self.selected_channels and self.channel_buttons.get(ch) and self.channel_buttons[ch].isChecked()]
        parts.append("Channels: " + (", ".join(channels) if channels else "none"))
        parts.append(f"Overlay: {'ON' if self.overlay_button.isChecked() else 'OFF'}")

        fits_on = [ch for ch, state in self._show_fit.items() if state]
        if fits_on:
            parts.append("Fits: " + ", ".join(sorted(fits_on)))

        if missing:
            parts.append("Missing: " + ", ".join(sorted(set(missing))))

        edited = sorted({ch for (evt, ch, _path) in self._fit_param_overrides.keys() if evt == self._current_event})
        if edited:
            parts.append("Edited params: " + ", ".join(edited))

        self.statusBar().showMessage(" | ".join(parts))

    # ---- Fit data helpers -----------------------------------------------
    def get_fit_data(self, event_name: str, channel: str) -> FitData:
        key = (event_name, channel)
        if key in self._fit_cache:
            return self._fit_cache[key]

        data = FitData()
        if not self._h5:
            self._fit_cache[key] = data
            return data

        group = self._h5.get(event_name)
        if not isinstance(group, h5py.Group):
            self._fit_cache[key] = data
            return data

        base_norm = _normalize_name(channel)

        def visitor(name: str, obj):
            if not isinstance(obj, h5py.Dataset):
                return
            dataset_name = name.split("/")[-1]
            norm = _normalize_name(dataset_name)
            if "fit" not in norm or base_norm not in norm:
                return
            try:
                arr = np.array(obj[()], copy=True)
            except Exception:
                return
            full_path = f"{group.name}/{name}"
            if "time" in norm:
                data.time_series[full_path] = arr
            elif "result" in norm or "curve" in norm:
                data.value_series[full_path] = arr
            elif "param" in norm:
                data.parameter_series[full_path] = arr
            else:
                data.extras[full_path] = arr

        group.visititems(visitor)

        for path in list(data.time_series.keys()):
            override_time = self._fit_time_overrides.get((event_name, channel, path))
            if override_time is not None:
                data.time_series[path] = np.array(override_time, copy=True)

        for path in list(data.value_series.keys()):
            override_result = self._fit_result_overrides.get((event_name, channel, path))
            if override_result is not None:
                data.value_series[path] = np.array(override_result, copy=True)

        scale = FIT_SCALE_MULTIPLIERS.get(channel, 1.0)
        if scale != 1.0:
            for path, values in list(data.value_series.items()):
                data.value_series[path] = np.asarray(values, dtype=float) * scale

        self._fit_cache[key] = data
        return data

    def open_fit_parameter_dialog(self):
        if not self._current_event or not self._h5:
            QMessageBox.information(self, "No Event", "Open an event before editing fit parameters.")
            return

        available: Dict[str, FitData] = {}
        for channel in FIT_ELIGIBLE_CHANNELS:
            data = self.get_fit_data(self._current_event, channel)
            if data.has_parameters():
                available[channel] = data

        if not available:
            QMessageBox.information(self, "No Parameters", "No editable fit parameters were found for this event.")
            return

        dialog = FitParameterDialog(
            parent=self,
            event_name=self._current_event,
            channel_data=available,
            value_getter=self.get_parameter_values,
            save_callback=self.update_fit_override,
            reset_callback=self.clear_fit_override,
        )
        dialog.exec()

    def get_parameter_values(self, event: str, channel: str, dataset_path: str, base_array: np.ndarray) -> np.ndarray:
        key = (event, channel, dataset_path)
        override = self._fit_param_overrides.get(key)
        if override is not None:
            return np.array(override, copy=True)
        return np.array(base_array, copy=True)

    def update_fit_override(self, event: str, channel: str, dataset_path: str, values: np.ndarray):
        array = np.array(values, copy=True)
        self._fit_param_overrides[(event, channel, dataset_path)] = array
        self._remove_fit_override(event, channel, dataset_path)
        recomputed = self._recalculate_fit(event, channel, dataset_path, array)
        if recomputed:
            message = f"Updated fit parameters for {channel}; recomputed fit curve (in-memory only)."
        else:
            message = f"Updated fit parameters for {channel}; fit curve unchanged (missing data)."
        self.statusBar().showMessage(message, 6000)
        self._fit_cache.pop((event, channel), None)
        self.plot_event(self._current_event)

    def clear_fit_override(self, event: str, channel: str, dataset_path: str):
        key = (event, channel, dataset_path)
        removed_param = self._fit_param_overrides.pop(key, None) is not None
        removed_curve = self._remove_fit_override(event, channel, dataset_path)
        if removed_param or removed_curve:
            self.statusBar().showMessage(f"Reverted fit parameters for {channel}.", 6000)
            self._fit_cache.pop((event, channel), None)
            self.plot_event(self._current_event)

    def _recalculate_fit(self, event: str, channel: str, dataset_path: str, params: np.ndarray) -> bool:
        if not self._h5:
            return False
        derived = _fit_paths_from_param(dataset_path)
        if not derived:
            return False
        result_path, time_path = derived
        time_data = dset(self._h5, time_path)
        if time_data is None:
            definition = CHANNEL_DEFS.get(channel)
            if definition:
                time_data = dset(self._h5, f"/{event}/{definition.time_dataset}")
        if time_data is None:
            return False
        time_array = np.array(time_data, copy=True)
        flat_time = np.asarray(time_array, dtype=float).ravel()
        curve = _evaluate_fit_curve(channel, flat_time, params)
        if curve is None:
            return False
        curve_array = np.array(curve, copy=True)
        if curve_array.size != time_array.size:
            return False
        curve_array = curve_array.reshape(time_array.shape)
        self._fit_time_overrides[(event, channel, time_path)] = time_array
        self._fit_result_overrides[(event, channel, result_path)] = curve_array
        self._fit_cache.pop((event, channel), None)
        return True

    def _remove_fit_override(self, event: str, channel: str, dataset_path: str) -> bool:
        removed = False
        derived = _fit_paths_from_param(dataset_path)
        if derived:
            result_path, time_path = derived
            if (event, channel, result_path) in self._fit_result_overrides:
                del self._fit_result_overrides[(event, channel, result_path)]
                removed = True
            if (event, channel, time_path) in self._fit_time_overrides:
                del self._fit_time_overrides[(event, channel, time_path)]
                removed = True
        return removed

    # --- Cleanup ---------------------------------------------------------
    def closeEvent(self, event):
        try:
            if self._h5 is not None:
                self._h5.close()
        except Exception:
            pass
        try:
            self._tmpdir.cleanup()
        except Exception:
            pass
        event.accept()

# --------- Fit parameter editor dialog ---------
class FitParameterDialog(QDialog):
    def __init__(
        self,
        parent: QWidget,
        event_name: str,
        channel_data: Dict[str, FitData],
        value_getter: Callable[[str, str, str, np.ndarray], np.ndarray],
        save_callback: Callable[[str, str, str, np.ndarray], None],
        reset_callback: Callable[[str, str, str], None],
    ):
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle(f"Fit Parameters — Event {event_name}")
        self.setMinimumSize(580, 520)

        self._event_name = event_name
        self._channel_data = channel_data
        self._value_getter = value_getter
        self._save_callback = save_callback
        self._reset_callback = reset_callback

        self._param_arrays: Dict[Tuple[str, str], np.ndarray] = {}
        for channel, data in channel_data.items():
            for path, array in data.iter_parameter_items():
                self._param_arrays[(channel, path)] = np.asarray(array)

        self._current_channel: Optional[str] = None
        self._current_dataset: Optional[str] = None
        self._current_shape: Optional[Tuple[int, ...]] = None
        self._is_updating_table = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        header = QLabel(f"Editing parameters for event {event_name}", self)
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        self.formula_label = QLabel("", self)
        self.formula_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.formula_label.setStyleSheet("font-size: 16px;")
        self.formula_label.setWordWrap(True)
        self.formula_label.setMinimumHeight(70)
        self._default_formula_text = "Select a dataset to view its fit function."
        self.formula_label.setText(self._default_formula_text)
        layout.addWidget(self.formula_label)

        chooser_row = QHBoxLayout()
        chooser_row.setSpacing(10)

        self.channel_combo = QComboBox(self)
        self.channel_combo.setStyleSheet("font-size: 15px; min-height: 36px;")
        for channel in sorted(channel_data.keys()):
            self.channel_combo.addItem(channel)
        self.channel_combo.currentTextChanged.connect(self._on_channel_changed)
        chooser_row.addWidget(QLabel("Channel:", self))
        chooser_row.addWidget(self.channel_combo)

        self.dataset_combo = QComboBox(self)
        self.dataset_combo.setStyleSheet("font-size: 15px; min-height: 36px;")
        self.dataset_combo.currentIndexChanged.connect(self._on_dataset_changed)
        chooser_row.addWidget(QLabel("Dataset:", self))
        chooser_row.addWidget(self.dataset_combo)
        chooser_row.addStretch(1)
        layout.addLayout(chooser_row)

        self.table = QTableWidget(self)
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Fit parameter", "Value"])
        header_view = self.table.horizontalHeader()
        header_view.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header_view.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectItems)
        self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.table.setEditTriggers(
            QTableWidget.EditTrigger.DoubleClicked | QTableWidget.EditTrigger.SelectedClicked
        )
        self.table.setStyleSheet("font-size: 15px;")
        self.table.itemChanged.connect(self._on_table_item_changed)
        layout.addWidget(self.table)

        self.info_label = QLabel("", self)
        self.info_label.setStyleSheet("font-size: 14px; color: #555555;")
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)

        self.feedback_label = QLabel("", self)
        self.feedback_label.setStyleSheet("font-size: 14px; color: #146c43;")
        self.feedback_label.setWordWrap(True)
        layout.addWidget(self.feedback_label)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close, self)
        self.apply_button = QPushButton("Apply Changes", self)
        self.apply_button.setStyleSheet("font-size: 15px; font-weight: bold; padding: 8px 18px;")
        self.apply_button.clicked.connect(self.apply_changes)
        buttons.addButton(self.apply_button, QDialogButtonBox.ButtonRole.ActionRole)

        self.reset_button = QPushButton("Reset to File Values", self)
        self.reset_button.setStyleSheet("font-size: 15px; padding: 8px 18px;")
        self.reset_button.clicked.connect(self.reset_changes)
        buttons.addButton(self.reset_button, QDialogButtonBox.ButtonRole.ResetRole)

        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        if self.channel_combo.count() > 0:
            self._on_channel_changed(self.channel_combo.currentText())

    # ---- UI helpers -----------------------------------------------------
    def _on_channel_changed(self, channel: str):
        self.dataset_combo.blockSignals(True)
        self.dataset_combo.clear()
        entries = [path for (chan, path) in self._param_arrays.keys() if chan == channel]
        entries.sort(key=lambda value: _dataset_sort_key(channel, value))
        for path in entries:
            label = self._dataset_display_name(channel, path)
            self.dataset_combo.addItem(label, path)
        self.dataset_combo.blockSignals(False)

        if entries:
            self.dataset_combo.setCurrentIndex(0)
            self._display_dataset(channel, entries[0])
        else:
            self._display_dataset(channel, None)

    def _on_dataset_changed(self, index: int):
        channel = self.channel_combo.currentText()
        dataset_path = self.dataset_combo.itemData(index)
        self._display_dataset(channel, dataset_path)

    def _display_dataset(
        self,
        channel: Optional[str],
        dataset_path: Optional[str],
        *,
        preserve_feedback: bool = False,
    ):
        if not channel or dataset_path is None:
            self._update_formula_display(None)
            self._clear_table("No fit parameters are available for this channel.")
            return

        base_array = self._param_arrays.get((channel, dataset_path))
        if base_array is None:
            self._update_formula_display(None)
            self._clear_table("No fit parameters are available for this channel.")
            return

        self._update_formula_display(channel)

        values = self._value_getter(self._event_name, channel, dataset_path, base_array)
        array = np.asarray(values)
        self._current_channel = channel
        self._current_dataset = dataset_path
        self._current_shape = array.shape

        flat = array.ravel()
        labels = self._parameter_labels(channel, dataset_path, flat.size)

        self._is_updating_table = True
        try:
            self.table.setRowCount(flat.size)
            for idx, value in enumerate(flat):
                label = labels[idx] if idx < len(labels) else f"p{idx + 1}"
                idx_item = QTableWidgetItem(label)
                idx_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                value_item = QTableWidgetItem(self._format_value(value))
                value_item.setFlags(value_item.flags() | Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(idx, 0, idx_item)
                self.table.setItem(idx, 1, value_item)
        finally:
            self._is_updating_table = False

        self.table.resizeColumnsToContents()
        display_name = self._dataset_display_name(channel, dataset_path)
        self.info_label.setText(f"{channel} — {display_name}\n{dataset_path} — shape {array.shape}")
        if not preserve_feedback:
            self.feedback_label.clear()

    def _dataset_display_name(self, channel: str, dataset_path: str) -> str:
        mass_label = _mass_identifier_from_path(dataset_path)
        if mass_label:
            return f"Mass {mass_label}"
        base = os.path.basename(dataset_path)
        return base or dataset_path

    def _parameter_labels(self, channel: str, dataset_path: str, count: int) -> List[str]:
        model = FIT_MODEL_BY_CHANNEL.get(channel)
        if not model:
            return [f"p{i + 1}" for i in range(count)]
        labels = list(model.parameter_labels)
        if len(labels) >= count:
            return labels[:count]
        extras = [f"p{i + 1}" for i in range(len(labels), count)]
        return labels + extras

    def _update_formula_display(self, channel: Optional[str]):
        self.formula_label.setPixmap(QPixmap())
        self.formula_label.setToolTip("")
        if not channel:
            self.formula_label.setText(self._default_formula_text)
            return

        model = FIT_MODEL_BY_CHANNEL.get(channel)
        if not model or not model.latex:
            self.formula_label.setText("Fit function preview unavailable for this selection.")
            return

        pixmap = _latex_to_pixmap(model.latex)
        if pixmap:
            self.formula_label.setText("")
            self.formula_label.setPixmap(pixmap)
        else:
            self.formula_label.setText(model.latex)
        self.formula_label.setToolTip(model.latex)

    def _set_feedback(self, text: str, *, success: bool):
        if not text:
            self.feedback_label.clear()
            return
        color = "#146c43" if success else "#b02a37"
        self.feedback_label.setStyleSheet(f"font-size: 14px; color: {color};")
        self.feedback_label.setText(text)

    def _apply_current_values(self, *, auto: bool) -> bool:
        if not self._current_channel or not self._current_dataset:
            if not auto:
                QMessageBox.information(
                    self,
                    "No Selection",
                    "Select a parameter set before applying changes.",
                )
            return False

        values: List[float] = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 1)
            text = item.text().strip() if item else ""
            if not text:
                values.append(0.0)
                continue
            try:
                values.append(float(text))
            except ValueError:
                message = f"Row {row + 1} contains a non-numeric value: {text}"
                if auto:
                    self._set_feedback(message, success=False)
                else:
                    QMessageBox.warning(self, "Invalid value", message)
                return False

        if self._current_shape:
            expected = int(np.prod(self._current_shape))
        else:
            expected = len(values)

        if len(values) != expected:
            message = "The number of edited values does not match the original parameter shape."
            if auto:
                self._set_feedback(message, success=False)
            else:
                QMessageBox.warning(self, "Shape mismatch", message)
            return False

        array = np.array(values, dtype=float)
        if self._current_shape:
            array = array.reshape(self._current_shape)

        self._save_callback(self._event_name, self._current_channel, self._current_dataset, array)
        self._display_dataset(self._current_channel, self._current_dataset, preserve_feedback=True)

        if auto:
            self._set_feedback("Fit updated with latest parameters.", success=True)
        else:
            self._set_feedback("Changes stored in-memory and fit recomputed.", success=True)
        return True

    def _on_table_item_changed(self, item: QTableWidgetItem):
        if self._is_updating_table or item is None or item.column() != 1:
            return
        self._apply_current_values(auto=True)

    def _clear_table(self, message: str = ""):
        self.table.setRowCount(0)
        self.info_label.setText(message)
        self.feedback_label.clear()
        self._current_channel = None
        self._current_dataset = None
        self._current_shape = None
        self._is_updating_table = False

    def _format_value(self, value: object) -> str:
        try:
            if isinstance(value, (np.floating, float)):
                return f"{float(value):.6g}"
            if isinstance(value, (np.integer, int)):
                return str(int(value))
        except Exception:
            pass
        return str(value)

    # ---- Actions --------------------------------------------------------
    def apply_changes(self):
        self._apply_current_values(auto=False)

    def reset_changes(self):
        if not self._current_channel or not self._current_dataset:
            return
        self._reset_callback(self._event_name, self._current_channel, self._current_dataset)
        self._display_dataset(self._current_channel, self._current_dataset, preserve_feedback=True)
        self._set_feedback("Reverted to values from the file.", success=False)

# --------- CLI / main ---------
def main():
    parser = argparse.ArgumentParser(description="Run the IDEX Quicklook (QtAgg + Matplotlib toolbar).")
    parser.add_argument("--filename", nargs="?", default=None, help="Path to the HDF5 file.")
    parser.add_argument("--eventnumber", nargs="?", type=int, default=None, help="1-based event index.")
    args = parser.parse_args()
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    base_font = QFont(app.font())
    size = base_font.pointSize()
    if size <= 0:
        size = 12
    else:
        size = max(size, 12)
    base_font.setPointSize(size)
    app.setFont(base_font)

    w = MainWindow(filename=args.filename, eventnumber=args.eventnumber)
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
