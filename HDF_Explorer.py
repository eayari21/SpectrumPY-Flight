"""Interactive dashboard for exploring IDEX HDF5 analysis products.

This module replaces the legacy :mod:`HDF_Plotter` tool with a richer
experience that focuses on per-event scalar summaries, ratio analysis and a
modernised time-series explorer.  The widgets are intentionally designed to
mirror the data model used by :mod:`IDEX-quicklook` so that selections made in
this dashboard seamlessly hand off to the Quicklook viewer when users click on
a point.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import importlib.util

import h5py  # type: ignore[import]
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.collections import PathCollection
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon, QPixmap
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

# --- Quicklook integration -------------------------------------------------

file_path = "./IDEX-quicklook.py"
module_name = "IDEX_quicklook"

spec = importlib.util.spec_from_file_location(module_name, file_path)
if spec is None or spec.loader is None:
    raise ImportError("Unable to locate IDEX-quicklook module")
module = importlib.util.module_from_spec(spec)
sys.modules[module_name] = module
spec.loader.exec_module(module)

from IDEX_quicklook import (  # type: ignore[attr-defined]
    FIT_MODEL_BY_CHANNEL,
    MainWindow as QuicklookWindow,
    _label_from_param_path,
    _mass_identifier_from_path,
    _normalize_name,
)

# --- Utility containers ----------------------------------------------------

IMAGES_DIR = Path(__file__).resolve().parent / "Images"
IMAP_LOGO_CANDIDATES = (
    "IMAP_logo.png",
    "IMAP_logo.jpg",
    "IMAP_logo.jpeg",
    "imap_logo.png",
    "IMAP.png",
    "IMAP.jpeg",
)

FIT_PARAM_SUFFIXES = ("fitparams", "fitparameters", "fitparam")
TIME_KEYWORDS = ("time", "epoch")

SCALAR_SENTINEL = math.nan

sns.set_theme(style="whitegrid")


@dataclass
class EventRecord:
    """Represents a single analysis event within an HDF5 file."""

    filename: str
    event_name: str
    epoch: Optional[float]


@dataclass
class DataStore:
    """In-memory representation of the extracted HDF5 analysis content."""

    events: List[EventRecord] = field(default_factory=list)
    scalars: Dict[str, List[float]] = field(default_factory=dict)
    arrays: Dict[str, List[Optional[np.ndarray]]] = field(default_factory=dict)
    array_time_keys: Dict[str, Optional[str]] = field(default_factory=dict)
    dataset_paths: Dict[str, str] = field(default_factory=dict)

    def append_event(self, record: EventRecord) -> int:
        index = len(self.events)
        self.events.append(record)
        for values in self.scalars.values():
            values.append(SCALAR_SENTINEL)
        for values in self.arrays.values():
            values.append(None)
        return index

    def ensure_scalar_entry(self, key: str) -> List[float]:
        if key not in self.scalars:
            length = len(self.events)
            self.scalars[key] = [SCALAR_SENTINEL] * length
        return self.scalars[key]

    def ensure_array_entry(self, key: str) -> List[Optional[np.ndarray]]:
        if key not in self.arrays:
            length = len(self.events)
            self.arrays[key] = [None] * length
        return self.arrays[key]


# --- Helper widgets --------------------------------------------------------


@dataclass
class AxisSelection:
    mode: str  # "value" or "ratio"
    primary: Optional[str]
    secondary: Optional[str]


class AxisSelector(QWidget):
    def __init__(self, label: str, parent: Optional[QWidget] = None):
        super().__init__(parent)

        self._label = QLabel(label)
        self._mode = QComboBox()
        self._mode.addItems(["Value", "Ratio"])
        self._primary = QComboBox()
        self._secondary = QComboBox()

        self._secondary.setVisible(False)

        self._mode.currentIndexChanged.connect(self._update_visibility)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        layout.addWidget(self._label)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(self._primary)

        self._secondary.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        row.addWidget(self._secondary)

        layout.addLayout(row)

    # ------------------------------------------------------------------
    def _update_visibility(self):
        is_ratio = self._mode.currentText().lower() == "ratio"
        self._secondary.setVisible(is_ratio)

    def set_items(self, labels: Sequence[str]):
        def _apply(combo: QComboBox):
            current = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(labels)
            if current and current in labels:
                combo.setCurrentText(current)
            combo.blockSignals(False)

        _apply(self._primary)
        _apply(self._secondary)
        self._update_visibility()

    def current(self) -> AxisSelection:
        mode = self._mode.currentText().lower()
        primary = self._primary.currentText() or None
        secondary = self._secondary.currentText() or None
        if mode != "ratio":
            secondary = None
        return AxisSelection(mode=mode, primary=primary, secondary=secondary)


# --- Core application ------------------------------------------------------


def _find_logo_path() -> Optional[Path]:
    for candidate in IMAP_LOGO_CANDIDATES:
        path = IMAGES_DIR / candidate
        if path.exists():
            return path
    return None


def _parameter_labels_for(channel: Optional[str], count: int) -> List[str]:
    if not channel:
        return [f"p{i + 1}" for i in range(count)]
    model = FIT_MODEL_BY_CHANNEL.get(channel)
    if not model:
        return [f"p{i + 1}" for i in range(count)]
    labels = list(model.parameter_labels)
    if len(labels) >= count:
        return labels[:count]
    extras = [f"p{i + 1}" for i in range(len(labels), count)]
    return labels + extras


def _friendly_scalar_label(base_label: str, channel: Optional[str], mass_label: Optional[str], param_label: Optional[str] = None) -> str:
    pieces: List[str] = []
    if channel:
        pieces.append(channel)
    if mass_label:
        pieces.append(f"Mass {mass_label}")
    if not pieces:
        pieces.append(base_label)
    if param_label:
        pieces.append(param_label)
    return " — ".join(pieces)


def _infer_channel_from_name(name: str) -> Optional[str]:
    target = _normalize_name(name)
    matches: List[Tuple[int, str]] = []
    for channel in FIT_MODEL_BY_CHANNEL:
        norm = _normalize_name(channel)
        if norm in target:
            matches.append((len(norm), channel))
    if not matches:
        return None
    matches.sort(reverse=True)
    return matches[0][1]


def _is_time_dataset(name: str) -> bool:
    lowered = name.lower()
    return any(token in lowered for token in TIME_KEYWORDS)


class HDFDataExplorer(QWidget):
    def __init__(self, source_path: str | os.PathLike[str] = "HDF5"):
        super().__init__()
        self.setWindowTitle("SpectrumPY — HDF Explorer")

        logo_path = _find_logo_path()
        if logo_path is not None:
            self.setWindowIcon(QIcon(str(logo_path)))

        resolved = Path(source_path)
        if not resolved.is_absolute():
            resolved = (Path(__file__).resolve().parent / resolved).resolve()

        if resolved.is_file():
            self.hdf5_folder = resolved.parent
        else:
            self.hdf5_folder = resolved

        if not self.hdf5_folder.exists():
            QMessageBox.critical(
                self,
                "Data location not found",
                f"The selected data location '{self.hdf5_folder}' does not exist.",
            )
            self.data_store = DataStore()
        else:
            self.data_store = self._load_datasets()

        self._quicklook_window: Optional[QuicklookWindow] = None
        self._scatter_artists: Dict[PathCollection, np.ndarray] = {}
        self._pick_cid: Optional[int] = None

        self._build_ui()

    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(12)

        logo_pixmap = self._load_logo_pixmap()
        if logo_pixmap is not None:
            logo_label = QLabel()
            logo_label.setPixmap(logo_pixmap)
            logo_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            header_layout.addWidget(logo_label)

        title_label = QLabel("SpectrumPY: Flight Addition — HDF Explorer")
        title_label.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)
        title_label.setStyleSheet("font-size: 18px; font-weight: 600;")
        header_layout.addWidget(title_label)
        header_layout.addStretch()

        layout.addWidget(header)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter, stretch=1)

        controls = QWidget()
        controls_layout = QVBoxLayout(controls)
        controls_layout.setContentsMargins(0, 0, 0, 0)
        controls_layout.setSpacing(12)

        self.axis_x1 = AxisSelector("Panel 1 — X axis:")
        self.axis_y1 = AxisSelector("Panel 1 — Y axis:")
        self.axis_c1 = AxisSelector("Panel 1 — Color:")

        self.axis_x2 = AxisSelector("Panel 2 — X axis:")
        self.axis_y2 = AxisSelector("Panel 2 — Y axis:")
        self.axis_c2 = AxisSelector("Panel 2 — Color:")

        for widget in (self.axis_x1, self.axis_y1, self.axis_c1, self.axis_x2, self.axis_y2, self.axis_c2):
            controls_layout.addWidget(widget)
            line = QFrame()
            line.setFrameShape(QFrame.Shape.HLine)
            line.setFrameShadow(QFrame.Shadow.Sunken)
            controls_layout.addWidget(line)

        self.plot_button = QPushButton("Render scatter panels")
        self.plot_button.setStyleSheet("font-size: 15px; font-weight: 600; padding: 8px 14px;")
        self.plot_button.clicked.connect(self.plot_data)
        controls_layout.addWidget(self.plot_button)

        controls_layout.addSpacing(12)

        timeseries_group = QGroupBox("Timeseries explorer")
        timeseries_layout = QGridLayout(timeseries_group)
        timeseries_layout.setContentsMargins(12, 12, 12, 12)
        timeseries_layout.setHorizontalSpacing(8)
        timeseries_layout.setVerticalSpacing(6)

        self.timeseries_combo = QComboBox()
        self.timeseries_combo.setEditable(False)
        self.timeseries_combo.currentIndexChanged.connect(self.update_timeseries_plot)

        self.timeseries_button = QPushButton("Update timeseries")
        self.timeseries_button.clicked.connect(self.update_timeseries_plot)

        timeseries_layout.addWidget(QLabel("Quantity"), 0, 0)
        timeseries_layout.addWidget(self.timeseries_combo, 0, 1)
        timeseries_layout.addWidget(self.timeseries_button, 1, 0, 1, 2)

        controls_layout.addWidget(timeseries_group)
        controls_layout.addStretch(1)

        splitter.addWidget(controls)

        plot_panel = QWidget()
        plot_layout = QVBoxLayout(plot_panel)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(6)

        self.figure = plt.figure(figsize=(10, 11), constrained_layout=False)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, stretch=1)

        splitter.addWidget(plot_panel)
        splitter.setSizes([320, 780])

        self._refresh_selectors()
        self.plot_data()

    # ------------------------------------------------------------------
    def _load_logo_pixmap(self, *, max_height: int = 72) -> Optional[QPixmap]:
        logo_path = _find_logo_path()
        if logo_path is None:
            return None
        pixmap = QPixmap(str(logo_path))
        if pixmap.isNull():
            return None
        return pixmap.scaledToHeight(
            max_height,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )

    # ------------------------------------------------------------------
    def _refresh_selectors(self) -> None:
        scalar_labels = sorted(self.data_store.scalars.keys())
        if not scalar_labels:
            placeholder = ["No scalar datasets detected"]
            for selector in (self.axis_x1, self.axis_y1, self.axis_c1, self.axis_x2, self.axis_y2, self.axis_c2):
                selector.set_items(placeholder)
            self.plot_button.setEnabled(False)
        else:
            for selector in (self.axis_x1, self.axis_y1, self.axis_c1, self.axis_x2, self.axis_y2, self.axis_c2):
                selector.set_items(scalar_labels)
            self.plot_button.setEnabled(True)

        timeseries_labels = sorted(self.data_store.scalars.keys()) + sorted(self.data_store.arrays.keys())
        if timeseries_labels:
            self.timeseries_combo.blockSignals(True)
            current = self.timeseries_combo.currentText()
            self.timeseries_combo.clear()
            self.timeseries_combo.addItems(timeseries_labels)
            if current in timeseries_labels:
                self.timeseries_combo.setCurrentText(current)
            self.timeseries_combo.blockSignals(False)
            self.timeseries_button.setEnabled(True)
        else:
            self.timeseries_combo.clear()
            self.timeseries_combo.addItem("No datasets available")
            self.timeseries_combo.setEnabled(False)
            self.timeseries_button.setEnabled(False)

    # ------------------------------------------------------------------
    def _load_datasets(self) -> DataStore:
        store = DataStore()

        for filename in sorted(os.listdir(self.hdf5_folder)):
            if not filename.endswith(".h5"):
                continue
            file_path = self.hdf5_folder / filename
            try:
                handle = h5py.File(str(file_path), "r")
            except OSError:
                continue

            with handle:
                for event_name, group in handle.items():
                    if not isinstance(group, h5py.Group):
                        continue

                    epoch = self._extract_epoch(group)
                    event_index = store.append_event(EventRecord(filename=filename, event_name=event_name, epoch=epoch))

                    analysis = group.get("Analysis")
                    if not isinstance(analysis, h5py.Group):
                        continue

                    datasets: Dict[str, np.ndarray] = {}
                    for name, dataset in analysis.items():
                        if not isinstance(dataset, h5py.Dataset):
                            continue
                        try:
                            data = np.array(dataset[()])
                        except Exception:
                            continue
                        datasets[name] = data

                    self._process_event_datasets(
                        store,
                        event_index,
                        group_name=group.name,
                        datasets=datasets,
                    )
        return store

    # ------------------------------------------------------------------
    def _extract_epoch(self, group: h5py.Group) -> Optional[float]:
        for name in ("Epoch", "epoch", "EventEpoch"):
            if name in group:
                obj = group[name]
                if isinstance(obj, h5py.Dataset):
                    try:
                        data = np.array(obj[()])
                    except Exception:
                        continue
                    if data.size == 0:
                        continue
                    value = float(np.ravel(data)[0])
                    if math.isfinite(value):
                        return value
        return None

    # ------------------------------------------------------------------
    def _process_event_datasets(
        self,
        store: DataStore,
        event_index: int,
        *,
        group_name: str,
        datasets: Dict[str, np.ndarray],
    ) -> None:
        time_series: Dict[str, np.ndarray] = {}
        value_series: Dict[str, np.ndarray] = {}

        for name, array in datasets.items():
            full_path = f"{group_name}/Analysis/{name}"
            store.dataset_paths.setdefault(name, full_path)

            if array.ndim == 0 or array.size == 1:
                series = store.ensure_scalar_entry(name)
                series[event_index] = float(np.ravel(array)[0])
                continue

            lowered = name.lower()
            if any(suffix in lowered for suffix in FIT_PARAM_SUFFIXES):
                self._ingest_fit_parameters(store, event_index, name, full_path, array)
                continue

            if _is_time_dataset(name):
                time_series[name] = np.ravel(array)
            else:
                value_series[name] = np.array(array)

        for name, arr in value_series.items():
            series = store.ensure_array_entry(name)
            series[event_index] = np.array(arr, copy=True)
            if name not in store.array_time_keys:
                store.array_time_keys[name] = self._guess_time_key(name, time_series)

        for name, arr in time_series.items():
            series = store.ensure_array_entry(name)
            series[event_index] = np.array(arr, copy=True)
            store.array_time_keys.setdefault(name, None)

    # ------------------------------------------------------------------
    def _guess_time_key(self, dataset_name: str, time_series: Dict[str, np.ndarray]) -> Optional[str]:
        if not time_series:
            return None
        norm_target = _normalize_name(dataset_name)
        for time_name in time_series:
            norm_time = _normalize_name(time_name)
            if norm_time == norm_target:
                return time_name
        for time_name in time_series:
            norm_time = _normalize_name(time_name)
            if any(token in norm_time for token in ("fit", "high", "low")) and any(token in norm_target for token in ("fit", "tof", "ion", "target")):
                return time_name
        return next(iter(time_series))

    # ------------------------------------------------------------------
    def _ingest_fit_parameters(
        self,
        store: DataStore,
        event_index: int,
        dataset_name: str,
        full_path: str,
        values: np.ndarray,
    ) -> None:
        flat = np.ravel(values)
        channel = _infer_channel_from_name(dataset_name)
        mass_label = _mass_identifier_from_path(full_path)
        base_label = _label_from_param_path(full_path)
        param_labels = _parameter_labels_for(channel, flat.size)

        for idx, value in enumerate(flat):
            if not np.isfinite(value):
                scalar_value = SCALAR_SENTINEL
            else:
                scalar_value = float(value)
            label = _friendly_scalar_label(base_label, channel, mass_label, param_labels[idx])
            series = store.ensure_scalar_entry(label)
            series[event_index] = scalar_value

    # ------------------------------------------------------------------
    def _resolve_scalar(self, key: Optional[str]) -> np.ndarray:
        if not key or key not in self.data_store.scalars:
            return np.array([])
        values = np.asarray(self.data_store.scalars[key], dtype=float)
        return values

    def _resolve_axis(self, selection: AxisSelection) -> np.ndarray:
        base = self._resolve_scalar(selection.primary)
        if selection.mode != "ratio" or not selection.secondary:
            return base
        denominator = self._resolve_scalar(selection.secondary)
        if base.size != denominator.size:
            return np.array([])
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.divide(base, denominator, out=np.full_like(base, np.nan), where=denominator != 0)
        return ratio

    # ------------------------------------------------------------------
    def plot_data(self) -> None:
        self.figure.clear()
        self._scatter_artists.clear()

        if not self.data_store.events:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No events found.", ha="center", va="center", fontsize=16, transform=ax.transAxes)
            ax.axis("off")
            self.canvas.draw_idle()
            return

        gs = self.figure.add_gridspec(3, 1, height_ratios=[3, 3, 2], hspace=0.6)

        ax_top = self.figure.add_subplot(gs[0, 0])
        ax_bottom = self.figure.add_subplot(gs[1, 0])
        self.ax_timeseries = self.figure.add_subplot(gs[2, 0])

        specs = [
            (self.axis_x1.current(), self.axis_y1.current(), self.axis_c1.current(), ax_top, "Panel 1"),
            (self.axis_x2.current(), self.axis_y2.current(), self.axis_c2.current(), ax_bottom, "Panel 2"),
        ]

        palette = "viridis"
        point_size = 32

        for x_sel, y_sel, c_sel, axis, title in specs:
            x_values = self._resolve_axis(x_sel)
            y_values = self._resolve_axis(y_sel)
            color_values = self._resolve_axis(c_sel)

            finite_mask = np.isfinite(x_values) & np.isfinite(y_values)
            if color_values.size == finite_mask.size:
                color_valid = np.isfinite(color_values)
                finite_mask &= color_valid
            indices = np.nonzero(finite_mask)[0]

            axis.clear()
            axis.set_title(f"{title}: {x_sel.primary or ''} vs {y_sel.primary or ''}", fontsize=14, fontweight="bold")
            axis.set_facecolor("#f8fafc")
            axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)

            if indices.size == 0:
                axis.text(0.5, 0.5, "No overlapping data", transform=axis.transAxes, ha="center", va="center", fontsize=13)
                continue

            scatter = axis.scatter(
                x_values[indices],
                y_values[indices],
                c=color_values[indices] if color_values.size else None,
                cmap=palette,
                s=point_size,
                alpha=0.85,
                edgecolor="white",
                linewidth=0.6,
                picker=True,
            )
            if color_values.size:
                cbar = self.figure.colorbar(scatter, ax=axis)
                cbar.set_label(c_sel.primary or "")

            axis.set_xlabel(self._axis_label_from_selection(x_sel))
            axis.set_ylabel(self._axis_label_from_selection(y_sel))

            self._scatter_artists[scatter] = indices

        if self._pick_cid is not None:
            self.canvas.mpl_disconnect(self._pick_cid)
        self._pick_cid = self.canvas.mpl_connect("pick_event", self._on_pick)
        self.canvas.draw_idle()
        self.update_timeseries_plot()

    # ------------------------------------------------------------------
    def _axis_label_from_selection(self, selection: AxisSelection) -> str:
        if selection.mode != "ratio" or not selection.secondary:
            return selection.primary or ""
        return f"{selection.primary} / {selection.secondary}"

    # ------------------------------------------------------------------
    def update_timeseries_plot(self):
        if not hasattr(self, "ax_timeseries"):
            return

        key = self.timeseries_combo.currentText()
        axis = self.ax_timeseries
        axis.clear()
        axis.set_facecolor("#fbfdff")
        axis.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)

        if not key or (key not in self.data_store.scalars and key not in self.data_store.arrays):
            axis.text(0.5, 0.5, "Select a dataset to plot", ha="center", va="center", fontsize=13, transform=axis.transAxes)
            self.canvas.draw_idle()
            return

        if key in self.data_store.scalars:
            epochs = np.array([event.epoch if event.epoch is not None else np.nan for event in self.data_store.events])
            values = np.asarray(self.data_store.scalars[key], dtype=float)
            mask = np.isfinite(epochs) & np.isfinite(values)
            if np.count_nonzero(mask) == 0:
                axis.text(0.5, 0.5, "No epoch/time metadata for selected scalars.", ha="center", va="center", transform=axis.transAxes)
            else:
                axis.plot(epochs[mask], values[mask], marker="o", markersize=5, linewidth=1.4, color="#2563eb")
                axis.set_xlabel("Epoch")
                axis.set_ylabel(key)
        else:
            values_series = self.data_store.arrays.get(key, [])
            time_key = self.data_store.array_time_keys.get(key)
            times_series = self.data_store.arrays.get(time_key, []) if time_key else None

            plotted = False
            for idx, value in enumerate(values_series):
                if value is None or value.size == 0:
                    continue
                times = None
                if times_series:
                    times = times_series[idx] if idx < len(times_series) else None
                    if times is not None and times.size != value.size:
                        times = None
                if times is not None:
                    axis.plot(times, value, linewidth=1.0, label=f"Event {self.data_store.events[idx].event_name}")
                else:
                    axis.plot(value, linewidth=1.0, label=f"Event {self.data_store.events[idx].event_name}")
                plotted = True

            if plotted:
                axis.legend(loc="upper right", fontsize=8)
                axis.set_xlabel("Time" if time_key else "Sample index")
                axis.set_ylabel(key)
            else:
                axis.text(0.5, 0.5, "No timeseries values found for selection.", ha="center", va="center", transform=axis.transAxes)

        axis.set_title("Timeseries view", fontsize=14, fontweight="bold")
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    def _on_pick(self, event):
        artist = event.artist
        if not isinstance(artist, PathCollection):
            return
        indices = self._scatter_artists.get(artist)
        if indices is None:
            return
        chosen = event.ind
        if not chosen:
            return
        idx = indices[chosen[0]]
        if idx >= len(self.data_store.events):
            return
        record = self.data_store.events[idx]
        file_path = str(self.hdf5_folder / record.filename)
        self._launch_quicklook(file_path, record.event_name)

    # ------------------------------------------------------------------
    def _launch_quicklook(self, file_path: str, event_name: str) -> None:
        window = self._quicklook_window
        if window is not None and window.isVisible():
            window.raise_()
            window.activateWindow()
            window.open_file(file_path, preferred_event=event_name)
            return

        window = QuicklookWindow(filename=file_path)
        window.show()
        window.open_file(file_path, preferred_event=event_name)
        self._quicklook_window = window

    # ------------------------------------------------------------------
    def closeEvent(self, event):
        if self._quicklook_window is not None:
            try:
                self._quicklook_window.close()
            except Exception:
                pass
        event.accept()


# --- Entrypoint ------------------------------------------------------------


def main() -> None:
    app = QApplication(sys.argv)
    logo_path = _find_logo_path()
    if logo_path is not None:
        app.setWindowIcon(QIcon(str(logo_path)))

    folder = Path(__file__).resolve().parent / "HDF5"
    window = HDFDataExplorer(folder)
    window.resize(1100, 900)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
