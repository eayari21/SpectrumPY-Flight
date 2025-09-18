#!/opt/anaconda3/envs/idex/bin/python
# -*- coding: utf-8 -*-

# =====================================================================
# IDEX Quicklook — QtAgg canvas + Matplotlib Navigation Toolbar
# Panel layout: TL/TR, ML/MR, LL/LR (per user mapping)
# =====================================================================

import os
import sys
import argparse
import tempfile
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
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy
    )
    _QT = "PySide6"
except Exception:
    from PyQt6.QtCore import Qt, QSize
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy
    )
    _QT = "PyQt6"

print(f"[info] Qt binding: {_QT}, Matplotlib backend: {matplotlib.get_backend()}")

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

# --------- Plotting ---------
def build_axes(fig: Figure):
    """
    Create a 3x2 grid and return a dict with slots:
    TL, TR, ML, MR, LL, LR
    """
    gs = fig.add_gridspec(3, 2)
    axes = {
        "TL": fig.add_subplot(gs[0, 0]),
        "TR": fig.add_subplot(gs[0, 1]),
        "ML": fig.add_subplot(gs[1, 0]),
        "MR": fig.add_subplot(gs[1, 1]),
        "LL": fig.add_subplot(gs[2, 0]),
        "LR": fig.add_subplot(gs[2, 1]),
    }
    return axes

# Channel → slot mapping per your request
CHANNEL_SLOT = {
    "TOF L":    "TL",
    "Ion Grid": "TR",
    "TOF M":    "ML",
    "Target L": "MR",
    "TOF H":    "LL",
    "Target H": "LR",
}

# Time dataset name by channel family
TIME_NAME = {
    "TOF L":    "Time (high sampling)",
    "TOF M":    "Time (high sampling)",
    "TOF H":    "Time (high sampling)",
    "Ion Grid": "Time (low sampling)",
    "Target L": "Time (low sampling)",
    "Target H": "Time (low sampling)",
}

BOTTOM_SLOTS = {"LL", "LR"}  # add x-label only on bottom row

def y_label_with_units(channel_name: str) -> str:
    """Append units in square brackets after the channel name."""
    if channel_name.startswith("TOF"):
        return rf"{channel_name} [pC/$\Delta t$]"
    else:
        return rf"{channel_name} [pC]"

def plot_event_on_axes(h5: h5py.File, event_name: str, fig: Figure, axmap: dict):
    """Render six channels directly to the provided axes map (TL/TR/ML/MR/LL/LR)."""
    fig.suptitle(
        f"{os.path.basename(getattr(h5, 'filename', '?'))} — Event {event_name}",
        fontsize=14, fontweight="bold"
    )

    # Clear all first (keeps grid consistent on errors)
    for a in axmap.values():
        a.clear()

    for chan, slot in CHANNEL_SLOT.items():
        a = axmap[slot]
        tname = TIME_NAME[chan]
        t = dset(h5, f"/{event_name}/{tname}")
        y = dset(h5, f"/{event_name}/{chan}")

        if t is None or y is None:
            a.text(0.5, 0.5, f"Missing: '{chan}' or '{tname}'",
                   ha="center", va="center", transform=a.transAxes)
            a.set_ylabel(y_label_with_units(chan))
            a.grid(True, alpha=0.3)
            continue

        t = np.asarray(t).ravel()
        y = np.asarray(y).ravel()
        n = min(len(t), len(y))

        if n == 0:
            a.text(0.5, 0.5, "Empty dataset", ha="center", va="center", transform=a.transAxes)
        else:
            a.plot(t[:n], y[:n], lw=1.0)

        a.set_ylabel(y_label_with_units(chan))
        a.grid(True, alpha=0.3)

    # X labels on bottom row only
    for s in BOTTOM_SLOTS:
        axmap[s].set_xlabel(r"Time [$\mu$s]")

# --------- Main Window ---------
class MainWindow(QMainWindow):
    def __init__(self, filename: str | None = None, eventnumber: int | None = None):
        super().__init__()
        self.setWindowTitle("IDEX Quicklook — QtAgg + Matplotlib Toolbar")
        self.setMinimumSize(1100, 750)
        self.setStatusBar(QStatusBar(self))

        self._h5: h5py.File | None = None
               # list of str
        self._events: list[str] = []
        self._current_event: str | None = None
        self._filename: str | None = None
        self._tmpdir = tempfile.TemporaryDirectory(prefix="idex_quicklook_")

        # Central widget/layout
        central = QWidget(self)
        self.vbox = QVBoxLayout(central)
        self.setCentralWidget(central)

        # ---- File/Actions toolbar ----
        tb = QToolBar("Main", self)
        tb.setIconSize(QSize(18, 18))
        self.addToolBar(tb)

        act_open = QAction("Open HDF5…", self)
        act_open.triggered.connect(self.action_open)
        tb.addAction(act_open)

        act_reload = QAction("Reload", self)
        act_reload.triggered.connect(self.reload_current)
        tb.addAction(act_reload)

        tb.addSeparator()

        act_quit = QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        tb.addAction(act_quit)

        # Event selector
        self.event_combo = QComboBox(self)
        self.event_combo.setMinimumWidth(220)
        self.event_combo.currentIndexChanged.connect(self.on_event_changed)
        tb.addSeparator()
        tb.addWidget(QLabel("Event:", self))
        tb.addWidget(self.event_combo)

        # ---- Matplotlib canvas + Navigation toolbar ----
        self.figure = Figure(figsize=(12, 8), layout="constrained")
        self.axes = build_axes(self.figure)          # dict of TL/TR/ML/MR/LL/LR
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self.nav_toolbar = NavigationToolbar(self.canvas, self)

        self.vbox.addWidget(self.nav_toolbar)
        self.vbox.addWidget(self.canvas)

        # Open file from CLI or prompt
        if filename:
            self.open_file(filename)
        else:
            chosen = non_native_open_dialog(self)
            if not chosen:
                self.close()
                return
            self.open_file(chosen)

        # Pick event
        if self._events:
            if eventnumber is not None and 1 <= eventnumber <= len(self._events):
                idx = eventnumber - 1
            else:
                idx = 0
            self.event_combo.setCurrentIndex(idx)
            self.plot_event(self._events[idx])

    # --- Actions ---
    def action_open(self):
        chosen = non_native_open_dialog(self)
        if chosen:
            self.open_file(chosen)

    def open_file(self, path: str):
        # Close previous file
        try:
            if self._h5 is not None:
                self._h5.close()
        except Exception:
            pass
        self._h5 = None
        self._events = []
               # reset filename
        self._filename = None

        try:
            self._h5 = h5py.File(path, "r")
            self._filename = path
        except Exception as e:
            QMessageBox.critical(self, "Open Error", f"Failed to open file:\n{path}\n\n{e}")
            return

        ev = list_event_groups(self._h5)
        if not ev:
            QMessageBox.warning(self, "No Events", "No top-level event groups found.")
        self._events = ev
        self.event_combo.blockSignals(True)
        self.event_combo.clear()
        self.event_combo.addItems(self._events)
        self.event_combo.blockSignals(False)
        self._current_event = self._events[0] if self._events else None
        self.plot_event(self._current_event)

    def reload_current(self):
        if self._filename:
            self.open_file(self._filename)

    def on_event_changed(self, idx: int):
        if 0 <= idx < len(self._events):
            self.plot_event(self._events[idx])

    # --- Rendering ---
    def plot_event(self, event_name: str | None):
        if not self._h5 or not event_name:
            return
        self._current_event = event_name
        try:
            # Rebuild axes if needed (e.g., after reload)
            if not isinstance(self.axes, dict) or set(self.axes.keys()) != {"TL","TR","ML","MR","LL","LR"}:
                self.figure.clear()
                self.axes = build_axes(self.figure)
            plot_event_on_axes(self._h5, event_name, self.figure, self.axes)
            self.canvas.draw_idle()
        except Exception as e:
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Plot error:\n{e}", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            self.canvas.draw_idle()

    # --- Cleanup ---
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

# --------- CLI / main ---------
def main():
    parser = argparse.ArgumentParser(description="Run the IDEX Quicklook (QtAgg + Matplotlib toolbar).")
    parser.add_argument("--filename", nargs="?", default=None, help="Path to the HDF5 file.")
    parser.add_argument("--eventnumber", nargs="?", type=int, default=None, help="1-based event index.")
    args = parser.parse_args()
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = MainWindow(filename=args.filename, eventnumber=args.eventnumber)
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
