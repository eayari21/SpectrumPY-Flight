#!/opt/anaconda3/envs/idex/bin/python
# -*- coding: utf-8 -*-

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

# Qt-backed Matplotlib so NavigationToolbar works
import matplotlib
matplotlib.use("QtAgg")

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

# --------- Qt binding-agnostic imports (prefer PySide6, fallback PyQt6) ---------
_QT = None
try:
    from PySide6.QtCore import Qt, QSize
    from PySide6.QtGui import QAction, QFont
    from PySide6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy, QDialog, QPushButton,
        QHBoxLayout, QGridLayout, QTableWidget, QTableWidgetItem, QHeaderView,
        QCheckBox, QDialogButtonBox
    )
    _QT = "PySide6"
except Exception:
    from PyQt6.QtCore import Qt, QSize
    from PyQt6.QtGui import QAction, QFont
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

FIT_ELIGIBLE_CHANNELS = {"Ion Grid", "Target L", "Target H"}

FAMILY_TITLES = {
    FAMILY_HIGH: "TOF Channels (High Sampling)",
    FAMILY_LOW: "Ion Grid / Target Channels (Low Sampling)",
}

FAMILY_YLABELS = {
    FAMILY_HIGH: r"TOF Channels [pC/$\Delta t$]",
    FAMILY_LOW: r"Ion Grid / Target Channels [pC]",
}

TIME_AXIS_LABEL = r"Time [$\mu$s]"


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


def _evaluate_fit_curve(channel: str, time_values: np.ndarray, params: np.ndarray) -> Optional[np.ndarray]:
    """Evaluate a fit curve for the supplied channel using the stored parameters."""

    if channel not in FIT_ELIGIBLE_CHANNELS:
        return None

    param_array = np.asarray(params, dtype=float).ravel()
    if param_array.size < 5 or not np.all(np.isfinite(param_array[:5])):
        return None

    return _idex_ion_grid_model(time_values, *param_array[:5])


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
        return any(True for _ in self.iter_time_result_pairs())

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
        self._fit_multiplier_cache: Dict[Tuple[str, str, str], float] = {}
        self._fit_param_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._fit_result_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._fit_time_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._show_fit: Dict[str, bool] = {name: False for name in FIT_ELIGIBLE_CHANNELS}
        self.selected_channels = set(CHANNEL_ORDER)

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
        self._fit_multiplier_cache.clear()
        self._fit_param_overrides.clear()
        self._fit_result_overrides.clear()
        self._fit_time_overrides.clear()

        try:
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
        base_time: Optional[np.ndarray] = None
        base_values: Optional[np.ndarray] = None
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
                base_time = t[:n]
                base_values = y[:n]
                label = channel if overlay_mode else None
                ax.plot(base_time, base_values, label=label)
                base_plotted = True

        fit_plotted = self._plot_fit(
            ax,
            event_name,
            channel,
            overlay_mode,
            base_time=base_time,
            base_values=base_values,
        )

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

    def _plot_fit(
        self,
        ax,
        event_name: str,
        channel: str,
        overlay_mode: bool,
        base_time: Optional[np.ndarray] = None,
        base_values: Optional[np.ndarray] = None,
    ) -> bool:
        if not self._show_fit.get(channel):
            return False

        data = self.get_fit_data(event_name, channel)
        plotted = False
        for label, time_path, value_path, time_values, fit_values in data.iter_time_result_pairs():
            n = min(len(time_values), len(fit_values))
            if n == 0:
                continue
            scaled_values = self._scale_fit_values(
                event_name,
                channel,
                time_values[:n],
                fit_values[:n],
                value_path,
                base_time=base_time,
                base_values=base_values,
            )
            legend_label = f"{channel} {label}" if overlay_mode else label
            ax.plot(time_values[:n], scaled_values, linestyle="--", linewidth=2.2, label=legend_label)
            plotted = True
        return plotted

    def _scale_fit_values(
        self,
        event: str,
        channel: str,
        time_values: np.ndarray,
        fit_values: np.ndarray,
        value_path: str,
        base_time: Optional[np.ndarray] = None,
        base_values: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        array = np.asarray(fit_values, dtype=float)
        if array.size == 0:
            return array

        key = (event, channel, value_path)
        multiplier = self._fit_multiplier_cache.get(key)
        if multiplier is None:
            multiplier = self._infer_fit_multiplier(
                event,
                channel,
                time_values,
                array,
                base_time=base_time,
                base_values=base_values,
            )
            if multiplier is not None:
                self._fit_multiplier_cache[key] = multiplier

        if multiplier is None or not np.isfinite(multiplier) or multiplier <= 0.0:
            return array
        return array * multiplier

    def _infer_fit_multiplier(
        self,
        event: str,
        channel: str,
        time_values: np.ndarray,
        fit_values: np.ndarray,
        base_time: Optional[np.ndarray] = None,
        base_values: Optional[np.ndarray] = None,
    ) -> Optional[float]:
        times = np.asarray(time_values, dtype=float).ravel()
        curve = np.asarray(fit_values, dtype=float).ravel()
        if times.size == 0 or curve.size == 0:
            return None

        if base_time is None or base_values is None:
            if not self._h5:
                return None
            definition = CHANNEL_DEFS.get(channel)
            if not definition:
                return None
            base_time_data = dset(self._h5, f"/{event}/{definition.time_dataset}")
            base_value_data = dset(self._h5, f"/{event}/{definition.dataset}")
            if base_time_data is None or base_value_data is None:
                return None
            base_time = _to_1d(base_time_data)
            base_values = _to_1d(base_value_data)
        else:
            base_time = _to_1d(base_time)
            base_values = _to_1d(base_values)

        n = min(base_time.size, base_values.size)
        if n == 0:
            return None
        base_time = base_time[:n]
        base_values = base_values[:n]
        if base_time.size < 5 or curve.size < 5:
            return None

        order = np.argsort(base_time)
        base_time = base_time[order]
        base_values = base_values[order]
        unique_time, unique_indices = np.unique(base_time, return_index=True)
        base_time = unique_time
        base_values = base_values[unique_indices]
        if base_time.size < 5:
            return None

        interpolated = np.interp(times, base_time, base_values, left=np.nan, right=np.nan)
        valid = np.isfinite(interpolated) & np.isfinite(curve)
        if valid.sum() < max(8, curve.size // 4):
            return None

        target = interpolated[valid]
        model = curve[valid]
        target_center = target - np.median(target)
        model_center = model - np.median(model)

        signal_level = np.percentile(np.abs(model_center), 95)
        if not np.isfinite(signal_level) or signal_level <= 0.0:
            return None

        signal_mask = np.abs(model_center) >= (0.2 * signal_level)
        if signal_mask.sum() < 5:
            signal_mask = np.abs(model_center) >= (0.1 * signal_level)
        if signal_mask.sum() < 5:
            signal_mask = np.ones_like(model_center, dtype=bool)

        target_signal = target_center[signal_mask]
        model_signal = model_center[signal_mask]
        if model_signal.size < 3:
            return None

        amp_model = np.percentile(np.abs(model_signal), 90)
        amp_target = np.percentile(np.abs(target_signal), 90)
        if amp_model <= 0.0 or amp_target <= 0.0:
            return None

        scale = amp_target / amp_model
        if not np.isfinite(scale) or scale <= 0.0:
            return None
        if scale < 1e-8 or scale > 1e8:
            return None
        return scale

    def _clear_fit_multiplier_cache(self, event: str, channel: str, value_path: Optional[str] = None):
        if value_path is not None:
            self._fit_multiplier_cache.pop((event, channel, value_path), None)
            return

        keys_to_remove = [key for key in self._fit_multiplier_cache if key[0] == event and key[1] == channel]
        for key in keys_to_remove:
            self._fit_multiplier_cache.pop(key, None)

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
            available = data.has_overlay()
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
        recomputed, result_path = self._recalculate_fit(event, channel, dataset_path, array)
        if recomputed:
            message = f"Updated fit parameters for {channel}; recomputed fit curve (in-memory only)."
        else:
            message = f"Updated fit parameters for {channel}; fit curve unchanged (missing data)."
        self.statusBar().showMessage(message, 6000)
        self._fit_cache.pop((event, channel), None)
        if recomputed and result_path:
            self._clear_fit_multiplier_cache(event, channel, result_path)
        else:
            self._clear_fit_multiplier_cache(event, channel)
        self.plot_event(self._current_event)

    def clear_fit_override(self, event: str, channel: str, dataset_path: str):
        key = (event, channel, dataset_path)
        removed_param = self._fit_param_overrides.pop(key, None) is not None
        removed_curve = self._remove_fit_override(event, channel, dataset_path)
        if removed_param or removed_curve:
            self.statusBar().showMessage(f"Reverted fit parameters for {channel}.", 6000)
            self._fit_cache.pop((event, channel), None)
            self._clear_fit_multiplier_cache(event, channel)
            self.plot_event(self._current_event)

    def _recalculate_fit(
        self,
        event: str,
        channel: str,
        dataset_path: str,
        params: np.ndarray,
    ) -> Tuple[bool, Optional[str]]:
        if not self._h5:
            return False, None
        derived = _fit_paths_from_param(dataset_path)
        if not derived:
            return False, None
        result_path, time_path = derived
        time_data = dset(self._h5, time_path)
        if time_data is None:
            definition = CHANNEL_DEFS.get(channel)
            if definition:
                time_data = dset(self._h5, f"/{event}/{definition.time_dataset}")
        if time_data is None:
            return False, None
        time_array = np.array(time_data, copy=True)
        flat_time = np.asarray(time_array, dtype=float).ravel()
        params_array = np.asarray(params, dtype=float)
        curve = _evaluate_fit_curve(channel, flat_time, params_array)
        if curve is None:
            return False, None
        curve_array = np.array(curve, copy=True)
        if curve_array.size != time_array.size:
            return False, None
        curve_array = curve_array.reshape(time_array.shape)
        self._fit_time_overrides[(event, channel, time_path)] = time_array
        self._fit_result_overrides[(event, channel, result_path)] = curve_array
        self._fit_cache.pop((event, channel), None)
        return True, result_path

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
            self._clear_fit_multiplier_cache(event, channel, result_path)
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

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        header = QLabel(f"Editing parameters for event {event_name}", self)
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

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
        self.table.setHorizontalHeaderLabels(["Index", "Value"])
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
        for path in sorted(entries):
            self.dataset_combo.addItem(os.path.basename(path), path)
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

    def _display_dataset(self, channel: Optional[str], dataset_path: Optional[str]):
        if not channel or dataset_path is None:
            self._clear_table("No fit parameters are available for this channel.")
            return

        base_array = self._param_arrays.get((channel, dataset_path))
        if base_array is None:
            self._clear_table("No fit parameters are available for this channel.")
            return

        values = self._value_getter(self._event_name, channel, dataset_path, base_array)
        array = np.asarray(values)
        self._current_channel = channel
        self._current_dataset = dataset_path
        self._current_shape = array.shape

        flat = array.ravel()
        self.table.setRowCount(flat.size)
        for idx, value in enumerate(flat):
            idx_item = QTableWidgetItem(str(idx))
            idx_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            value_item = QTableWidgetItem(self._format_value(value))
            value_item.setFlags(value_item.flags() | Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(idx, 0, idx_item)
            self.table.setItem(idx, 1, value_item)

        self.table.resizeColumnsToContents()
        self.info_label.setText(f"{dataset_path} — shape {array.shape}")
        self.feedback_label.clear()

    def _clear_table(self, message: str = ""):
        self.table.setRowCount(0)
        self.info_label.setText(message)
        self.feedback_label.clear()
        self._current_channel = None
        self._current_dataset = None
        self._current_shape = None

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
        if not self._current_channel or not self._current_dataset:
            QMessageBox.information(self, "No Selection", "Select a parameter set before applying changes.")
            return

        values: List[float] = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 1)
            if item is None:
                continue
            text = item.text().strip()
            if not text:
                values.append(0.0)
                continue
            try:
                values.append(float(text))
            except ValueError:
                QMessageBox.warning(self, "Invalid value", f"Row {row} contains a non-numeric value: {text}")
                return

        if self._current_shape:
            expected = int(np.prod(self._current_shape))
        else:
            expected = len(values)

        if len(values) != expected:
            QMessageBox.warning(
                self,
                "Shape mismatch",
                "The number of edited values does not match the original parameter shape.",
            )
            return

        array = np.array(values, dtype=float)
        if self._current_shape:
            array = array.reshape(self._current_shape)

        self._save_callback(self._event_name, self._current_channel, self._current_dataset, array)
        self.feedback_label.setStyleSheet("font-size: 14px; color: #146c43;")
        self.feedback_label.setText("Changes stored in-memory for this session.")

    def reset_changes(self):
        if not self._current_channel or not self._current_dataset:
            return
        self._reset_callback(self._event_name, self._current_channel, self._current_dataset)
        self.feedback_label.setStyleSheet("font-size: 14px; color: #b02a37;")
        self.feedback_label.setText("Reverted to values from the file.")
        self._display_dataset(self._current_channel, self._current_dataset)

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
