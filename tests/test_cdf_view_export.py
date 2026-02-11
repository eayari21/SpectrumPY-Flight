from __future__ import annotations

import csv
import sys
import types
from pathlib import Path

import numpy as np


class _Dummy:
    def __init__(self, *args, **kwargs):
        pass


def _install_qt_stubs() -> None:
    if "PyQt6" in sys.modules:
        return

    qt_module = types.ModuleType("PyQt6")
    qtcore = types.ModuleType("PyQt6.QtCore")
    qtgui = types.ModuleType("PyQt6.QtGui")
    qtwidgets = types.ModuleType("PyQt6.QtWidgets")

    class _Qt:
        class ItemDataRole:
            UserRole = 0

        class ItemFlag:
            ItemIsEnabled = 1
            ItemIsSelectable = 2

        class Orientation:
            Horizontal = 1

    class _TreeWidget(_Dummy):
        class SelectionMode:
            ExtendedSelection = 1

    class _TableWidget(_Dummy):
        class EditTrigger:
            NoEditTriggers = 0

    class _HeaderView(_Dummy):
        class ResizeMode:
            Stretch = 0
            ResizeToContents = 1

    class _FileDialog(_Dummy):
        class Option:
            DontUseNativeDialog = 0

    qtcore.Qt = _Qt
    qtgui.QAction = _Dummy
    qtwidgets.QApplication = _Dummy
    qtwidgets.QFileDialog = _FileDialog
    qtwidgets.QHeaderView = _HeaderView
    qtwidgets.QLineEdit = _Dummy
    qtwidgets.QMainWindow = _Dummy
    qtwidgets.QMessageBox = _Dummy
    qtwidgets.QPlainTextEdit = _Dummy
    qtwidgets.QSplitter = _Dummy
    qtwidgets.QTableWidget = _TableWidget
    qtwidgets.QTableWidgetItem = _Dummy
    qtwidgets.QTabWidget = _Dummy
    qtwidgets.QToolBar = _Dummy
    qtwidgets.QTreeWidget = _TreeWidget
    qtwidgets.QTreeWidgetItem = _Dummy
    qtwidgets.QVBoxLayout = _Dummy
    qtwidgets.QWidget = _Dummy

    qt_module.QtCore = qtcore
    qt_module.QtGui = qtgui
    qt_module.QtWidgets = qtwidgets

    sys.modules["PyQt6"] = qt_module
    sys.modules["PyQt6.QtCore"] = qtcore
    sys.modules["PyQt6.QtGui"] = qtgui
    sys.modules["PyQt6.QtWidgets"] = qtwidgets


_install_qt_stubs()

from spectrumpy_flight.CDF_View import CDFViewWindow


class _FakeCDF:
    def __init__(self, variables: dict[str, object]):
        self._variables = variables

    def cdf_info(self):
        return {"zVariables": list(self._variables.keys()), "rVariables": []}

    def varget(self, name: str):
        return self._variables[name]


def _make_window(tmp_path: Path, variables: dict[str, object]) -> CDFViewWindow:
    window = CDFViewWindow.__new__(CDFViewWindow)
    window._cdf = _FakeCDF(variables)
    window._filename = str(tmp_path / "sample.cdf")
    return window


def test_find_epoch_variable_falls_back_to_case_insensitive_match(tmp_path: Path):
    window = _make_window(tmp_path, {"EpoCh": [1.0, 2.0, 3.0], "counts": [7, 8, 9]})

    epoch, name = window._find_epoch_variable()

    assert name == "EpoCh"
    assert np.array_equal(epoch, np.array([1.0, 2.0, 3.0]))


def test_export_dataset_to_csv_uses_epoch_column(tmp_path: Path):
    window = _make_window(tmp_path, {"epoch": [10, 11], "counts": [100, 200]})
    output = tmp_path / "counts.csv"

    window._export_dataset_to_csv("counts", np.array([10, 11]), str(output))

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows == [["epoch", "counts"], ["10", "100"], ["11", "200"]]


def test_export_multiple_datasets_to_csv_handles_duplicate_names(tmp_path: Path):
    window = _make_window(tmp_path, {"epoch": [1, 2], "x": [5, 6]})
    output = tmp_path / "multi.csv"

    window._export_multiple_datasets_to_csv(["x", "x"], np.array([1, 2]), str(output))

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows[0] == ["epoch", "x", "x (2)"]
    assert rows[1] == ["1", "5", "5"]
    assert rows[2] == ["2", "6", "6"]
