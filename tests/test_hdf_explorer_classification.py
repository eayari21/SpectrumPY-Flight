"""Behavioural checks for :mod:`spectrumpy_flight.HDF_Explorer` classification helpers."""

from __future__ import annotations

import importlib
import sys
import types
from typing import Any

import pytest


def _install_qt_stubs(monkeypatch: pytest.MonkeyPatch) -> None:
    """Inject lightweight Qt stand-ins so :mod:`HDF_Explorer` can be imported.

    The actual widget toolkit is unavailable in the test environment.  We only
    need to satisfy attribute lookups that occur during module import; the
    stubs deliberately omit runtime behaviour because the tests exercise
    classification helpers in isolation.
    """

    qt_module = types.ModuleType("PyQt6")

    # ---- QtCore -----------------------------------------------------------------
    qtcore = types.ModuleType("PyQt6.QtCore")

    class _Qt:  # pragma: no cover - trivial data container
        class ItemDataRole:
            UserRole = 100
            FontRole = 101
            ForegroundRole = 102

        class AlignmentFlag:
            AlignLeft = 0x1
            AlignVCenter = 0x2

        class Orientation:
            Horizontal = 0

        class TransformationMode:
            SmoothTransformation = 0

    qtcore.Qt = _Qt
    monkeypatch.setitem(sys.modules, "PyQt6", qt_module)
    monkeypatch.setitem(sys.modules, "PyQt6.QtCore", qtcore)
    qt_module.QtCore = qtcore

    # ---- QtGui ------------------------------------------------------------------
    qtgui = types.ModuleType("PyQt6.QtGui")

    class QIcon:  # pragma: no cover - trivial placeholder
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

    class QPalette:  # pragma: no cover - trivial placeholder
        class ColorRole:
            Mid = 0

        def color(self, role: Any) -> Any:
            return None

    class QPixmap:  # pragma: no cover - trivial placeholder
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

        def isNull(self) -> bool:
            return True

        def scaledToHeight(self, *args: Any, **kwargs: Any) -> "QPixmap":
            return self

    qtgui.QIcon = QIcon
    qtgui.QPalette = QPalette
    qtgui.QPixmap = QPixmap
    monkeypatch.setitem(sys.modules, "PyQt6.QtGui", qtgui)
    qt_module.QtGui = qtgui

    # ---- QtWidgets --------------------------------------------------------------
    qtwidgets = types.ModuleType("PyQt6.QtWidgets")

    class _Widget:  # pragma: no cover - trivial placeholder
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

    class _Layout(_Widget):  # pragma: no cover - trivial placeholder
        def addWidget(self, *args: Any, **kwargs: Any) -> None:
            pass

        def addLayout(self, *args: Any, **kwargs: Any) -> None:
            pass

        def setContentsMargins(self, *args: Any, **kwargs: Any) -> None:
            pass

        def setSpacing(self, *args: Any, **kwargs: Any) -> None:
            pass

    class QSizePolicy:  # pragma: no cover - trivial placeholder
        class Policy:
            Expanding = object()
            Fixed = object()

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

    class QMessageBox:  # pragma: no cover - trivial placeholder
        @staticmethod
        def information(*args: Any, **kwargs: Any) -> None:
            pass

        @staticmethod
        def critical(*args: Any, **kwargs: Any) -> None:
            pass

    class QSplitter(_Widget):  # pragma: no cover - trivial placeholder
        def setChildrenCollapsible(self, *args: Any, **kwargs: Any) -> None:
            pass

        def setCollapsible(self, *args: Any, **kwargs: Any) -> None:
            pass

        def addWidget(self, *args: Any, **kwargs: Any) -> None:
            pass

        def setSizes(self, *args: Any, **kwargs: Any) -> None:
            pass

    qtwidgets.QApplication = _Widget
    qtwidgets.QComboBox = _Widget
    qtwidgets.QDialog = _Widget
    qtwidgets.QFrame = _Widget
    qtwidgets.QGridLayout = _Layout
    qtwidgets.QGroupBox = _Widget
    qtwidgets.QHBoxLayout = _Layout
    qtwidgets.QLabel = _Widget
    qtwidgets.QMessageBox = QMessageBox
    qtwidgets.QPushButton = _Widget
    qtwidgets.QSizePolicy = QSizePolicy
    qtwidgets.QSlider = _Widget
    qtwidgets.QSpinBox = _Widget
    qtwidgets.QSplitter = QSplitter
    qtwidgets.QVBoxLayout = _Layout
    qtwidgets.QWidget = _Widget
    monkeypatch.setitem(sys.modules, "PyQt6.QtWidgets", qtwidgets)
    qt_module.QtWidgets = qtwidgets

    # ---- Matplotlib Qt backend --------------------------------------------------
    backend = types.ModuleType("matplotlib.backends.backend_qtagg")

    class _Canvas:  # pragma: no cover - trivial placeholder
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

    class _Toolbar:  # pragma: no cover - trivial placeholder
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

    backend.FigureCanvasQTAgg = _Canvas
    backend.NavigationToolbar2QT = _Toolbar
    monkeypatch.setitem(sys.modules, "matplotlib.backends.backend_qtagg", backend)

    # ---- Quicklook shim --------------------------------------------------------
    quicklook = types.ModuleType("IDEX_quicklook")

    quicklook.FIT_MODEL_BY_CHANNEL = {}

    class _QuicklookWindow:  # pragma: no cover - trivial placeholder
        pass

    def _normalize_name(value: str) -> str:  # pragma: no cover - trivial placeholder
        return "".join(ch.lower() for ch in value if ch.isalnum())

    def _label_from_param_path(path: str) -> str:  # pragma: no cover - trivial placeholder
        return path

    def _mass_identifier_from_path(path: str) -> str:  # pragma: no cover - trivial placeholder
        return ""

    quicklook.MainWindow = _QuicklookWindow
    quicklook._normalize_name = _normalize_name
    quicklook._label_from_param_path = _label_from_param_path
    quicklook._mass_identifier_from_path = _mass_identifier_from_path
    monkeypatch.setitem(sys.modules, "IDEX_quicklook", quicklook)


@pytest.fixture(name="hdf_explorer")
def _hdf_explorer_fixture(monkeypatch: pytest.MonkeyPatch):
    """Return the imported :mod:`HDF_Explorer` module with Qt dependencies stubbed."""

    _install_qt_stubs(monkeypatch)
    return importlib.import_module("spectrumpy_flight.HDF_Explorer")


def _classify(module, label: str) -> str:
    instance = module.HDFDataExplorer.__new__(module.HDFDataExplorer)
    return module.HDFDataExplorer._classify_dataset_label(  # type: ignore[arg-type]
        instance,
        label=label,
        dataset_name=label,
        full_path=f"/Event/Analysis/{label}",
        category="analysis",
    )


def test_sql_match_series_grouped_with_accelerator_metadata(hdf_explorer):
    assert (
        _classify(hdf_explorer, "SQLMatch/VelocityKilometersPerSecond")
        == hdf_explorer.CATEGORY_ACCELERATOR
    )
    assert (
        _classify(hdf_explorer, "AcceleratorMetadata/VelocityKilometersPerSecond")
        == hdf_explorer.CATEGORY_ACCELERATOR
    )


def test_mass_line_products_stay_in_waveform_category(hdf_explorer):
    assert (
        _classify(hdf_explorer, "AcceleratorMetadata/TOF H Mass Lines")
        == hdf_explorer.CATEGORY_WAVEFORM
    )
    assert (
        _classify(hdf_explorer, "TOF H/MassLines")
        == hdf_explorer.CATEGORY_WAVEFORM
    )
