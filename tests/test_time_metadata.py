import importlib
import importlib.util
import sys
import types
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

SRC_DIR = Path(__file__).resolve().parents[1] / "src" / "spectrumpy_flight"


def _install_qt_stubs(monkeypatch: pytest.MonkeyPatch) -> None:
    if "PyQt6" in sys.modules:
        return

    qt_module = types.ModuleType("PyQt6")
    qtcore = types.ModuleType("PyQt6.QtCore")
    for name in ("Qt", "QSize", "QTimer", "QUrl", "QBuffer", "QIODevice", "QSignalBlocker"):
        setattr(qtcore, name, type(name, (), {}))
    qtgui = types.ModuleType("PyQt6.QtGui")
    for name in ("QAction", "QFont", "QIcon", "QPalette", "QPixmap", "QImage", "QTextCursor", "QTextDocument"):
        setattr(qtgui, name, type(name, (), {}))
    qtwidgets = types.ModuleType("PyQt6.QtWidgets")
    for name in (
        "QApplication",
        "QFileDialog",
        "QMainWindow",
        "QMessageBox",
        "QStatusBar",
        "QToolBar",
        "QVBoxLayout",
        "QWidget",
        "QComboBox",
        "QLabel",
        "QSizePolicy",
        "QDialog",
        "QPushButton",
        "QHBoxLayout",
        "QGridLayout",
        "QTableWidget",
        "QTableWidgetItem",
        "QHeaderView",
        "QCheckBox",
        "QDialogButtonBox",
        "QMenu",
        "QMenuBar",
        "QToolButton",
        "QTextBrowser",
        "QListWidget",
        "QListWidgetItem",
        "QLineEdit",
        "QWidgetAction",
        "QStyle",
        "QSplitter",
        "QSlider",
        "QScrollArea",
        "QFrame",
        "QGroupBox",
        "QDoubleSpinBox",
        "QSpinBox",
        "QPlainTextEdit",
        "QTabWidget",
        "QTreeWidget",
        "QTreeWidgetItem",
        "QAbstractItemView",
        "QFormLayout",
    ):
        setattr(qtwidgets, name, type(name, (), {}))

    monkeypatch.setitem(sys.modules, "PyQt6", qt_module)
    monkeypatch.setitem(sys.modules, "PyQt6.QtCore", qtcore)
    monkeypatch.setitem(sys.modules, "PyQt6.QtGui", qtgui)
    monkeypatch.setitem(sys.modules, "PyQt6.QtWidgets", qtwidgets)
    qt_module.QtCore = qtcore
    qt_module.QtGui = qtgui
    qt_module.QtWidgets = qtwidgets

    backend = types.ModuleType("matplotlib.backends.backend_qtagg")
    figure_canvas = type("FigureCanvasQTAgg", (), {"required_interactive_framework": None})
    backend.FigureCanvasQTAgg = figure_canvas
    backend.FigureCanvas = figure_canvas
    backend.NavigationToolbar2QT = type("NavigationToolbar2QT", (), {})
    monkeypatch.setitem(sys.modules, "matplotlib.backends.backend_qtagg", backend)


def _install_quicklook_dependency_stubs(monkeypatch: pytest.MonkeyPatch) -> None:
    modules: dict[str, dict[str, object]] = {
        "spectrumpy_flight.dust_composition": {
            "launch_dust_composition_window": lambda *args, **kwargs: None,
            "GAIN_HIGH": 1,
            "GAIN_LOW": 0,
        },
        "spectrumpy_flight.dust_estimator_gui": {
            "launch_dust_estimator_window": lambda *args, **kwargs: None,
        },
        "spectrumpy_flight.noise_analysis": {
            "ChannelMeta": type("ChannelMeta", (), {}),
            "launch_noise_analysis_window": lambda *args, **kwargs: None,
        },
        "spectrumpy_flight.plot_style": {
            "apply_plot_style": lambda *args, **kwargs: None,
        },
        "spectrumpy_flight.HDF_View": {
            "launch_hdf_viewer": lambda *args, **kwargs: None,
        },
        "spectrumpy_flight.CDF_View": {
            "launch_cdf_viewer": lambda *args, **kwargs: None,
        },
        "spectrumpy_flight.IDEX_Definitions_View": {
            "launch_variable_definitions_viewer": lambda *args, **kwargs: None,
        },
    }

    for module_name, attributes in modules.items():
        if module_name in sys.modules:
            continue
        module = types.ModuleType(module_name)
        for attr_name, value in attributes.items():
            setattr(module, attr_name, value)
        monkeypatch.setitem(sys.modules, module_name, module)


def _load_quicklook(monkeypatch: pytest.MonkeyPatch):
    _install_qt_stubs(monkeypatch)
    _install_quicklook_dependency_stubs(monkeypatch)
    import matplotlib
    monkeypatch.setattr(matplotlib, "use", lambda *args, **kwargs: None)
    spec = importlib.util.spec_from_file_location("IDEX_quicklook", SRC_DIR / "IDEX-quicklook.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load IDEX quicklook module for testing")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, "IDEX_quicklook", module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def quicklook_functions(monkeypatch: pytest.MonkeyPatch):
    module = _load_quicklook(monkeypatch)
    return module._first_finite_scalar, module._guess_event_timestamp_ms


@pytest.fixture
def hdf_explorer_class(monkeypatch: pytest.MonkeyPatch):
    _install_qt_stubs(monkeypatch)
    _install_quicklook_dependency_stubs(monkeypatch)
    import matplotlib

    monkeypatch.setattr(matplotlib, "use", lambda *args, **kwargs: None)
    module = importlib.import_module("spectrumpy_flight.HDF_Explorer")
    return module.HDFDataExplorer


class _DummySource:
    def __init__(self, datasets):
        self._datasets = datasets

    def get_dataset(self, event, dataset_name):
        return self._datasets.get(dataset_name)


def _utc_timestamp(year, month, day, hour=0, minute=0, second=0):
    return datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc).timestamp()


def test_first_finite_scalar_accepts_iso_timestamp(quicklook_functions):
    first_scalar, _ = quicklook_functions
    iso_value = "2023-07-25T12:34:56Z"
    expected = _utc_timestamp(2023, 7, 25, 12, 34, 56)
    assert first_scalar([iso_value]) == pytest.approx(expected)


def test_guess_event_timestamp_prefers_metadata_timestamp(quicklook_functions):
    _, guess_timestamp = quicklook_functions
    expected_seconds = _utc_timestamp(2023, 7, 25, 8, 30)
    dummy = _DummySource(
        {
            "Metadata/Timestamp": np.array([expected_seconds]),
            "Analysis/Trigger Time (ms)": np.array([2750.0]),
        }
    )
    result = guess_timestamp(dummy, "0002")
    assert result == pytest.approx(expected_seconds * 1000.0)


def test_metadata_epoch_from_timestamp_string(hdf_explorer_class):
    explorer = object.__new__(hdf_explorer_class)
    expected = _utc_timestamp(2023, 7, 25, 9, 15)
    datasets = {"Timestamp": np.array(["2023-07-25T09:15:00Z"], dtype=object)}
    assert explorer._metadata_epoch_from(datasets) == pytest.approx(expected)


def test_first_finite_value_accepts_numpy_datetime(hdf_explorer_class):
    explorer = object.__new__(hdf_explorer_class)
    expected = _utc_timestamp(2023, 7, 25, 10, 45)
    datasets = np.array([np.datetime64("2023-07-25T10:45:00")])
    assert explorer._first_finite_value(datasets) == pytest.approx(expected)


def test_metadata_epoch_from_numpy_datetime(hdf_explorer_class):
    explorer = object.__new__(hdf_explorer_class)
    expected = _utc_timestamp(2023, 7, 25, 11, 30)
    datasets = {"Timestamp": np.array([np.datetime64("2023-07-25T11:30:00")])}
    assert explorer._metadata_epoch_from(datasets) == pytest.approx(expected)
