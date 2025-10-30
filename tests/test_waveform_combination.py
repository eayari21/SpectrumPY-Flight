import sys
import types

import pytest


def _install_qt_stubs():
    if "PySide6" in sys.modules or "PyQt6" in sys.modules:
        return

    qt_module = types.ModuleType("PySide6")
    qt_module.__version__ = "0.0"
    qt_core = types.ModuleType("PySide6.QtCore")
    qt_gui = types.ModuleType("PySide6.QtGui")
    qt_widgets = types.ModuleType("PySide6.QtWidgets")

    Qt = types.SimpleNamespace(
        KeyboardModifier=types.SimpleNamespace(ControlModifier=0x01000000),
        ItemDataRole=types.SimpleNamespace(UserRole=32, FontRole=6, ForegroundRole=9),
        AlignmentFlag=types.SimpleNamespace(AlignCenter=0x0004, AlignRight=0x0002, AlignVCenter=0x0020),
        TextFormat=types.SimpleNamespace(RichText=1),
        ScrollBarPolicy=types.SimpleNamespace(ScrollBarAsNeeded=0, ScrollBarAlwaysOff=1),
        Orientation=types.SimpleNamespace(Horizontal=0),
        ItemFlag=types.SimpleNamespace(ItemIsEditable=0x0002),
    )

    palette_cls = type("QPalette", (), {"ColorRole": types.SimpleNamespace(Mid=0)})

    qt_core.Qt = Qt
    qt_gui.QPalette = palette_cls

    widget_names = [
        "QAbstractItemView",
        "QCheckBox",
        "QComboBox",
        "QDialog",
        "QDialogButtonBox",
        "QDoubleSpinBox",
        "QFormLayout",
        "QGridLayout",
        "QGroupBox",
        "QHBoxLayout",
        "QLabel",
        "QLineEdit",
        "QMainWindow",
        "QMessageBox",
        "QPushButton",
        "QScrollArea",
        "QSizePolicy",
        "QSplitter",
        "QStatusBar",
        "QTableWidget",
        "QTableWidgetItem",
        "QVBoxLayout",
        "QWidget",
        "QToolButton",
    ]

    for name in widget_names:
        setattr(qt_widgets, name, type(name, (), {}))

    qt_widgets.QSizePolicy.Policy = types.SimpleNamespace(Expanding=0)

    qt_module.QtCore = qt_core
    qt_module.QtGui = qt_gui
    qt_module.QtWidgets = qt_widgets

    sys.modules["PySide6"] = qt_module
    sys.modules["PySide6.QtCore"] = qt_core
    sys.modules["PySide6.QtGui"] = qt_gui
    sys.modules["PySide6.QtWidgets"] = qt_widgets


_install_qt_stubs()

if "matplotlib.backends.backend_qtagg" not in sys.modules:
    backend_stub = types.ModuleType("matplotlib.backends.backend_qtagg")

    class _Canvas:
        def __init__(self, *args, **kwargs):
            pass

    class _Toolbar:
        def __init__(self, *args, **kwargs):
            pass

    backend_stub.FigureCanvasQTAgg = _Canvas
    backend_stub.NavigationToolbar2QT = _Toolbar

    import importlib

    backends = importlib.import_module("matplotlib.backends")
    setattr(backends, "backend_qtagg", backend_stub)
    sys.modules["matplotlib.backends.backend_qtagg"] = backend_stub

np = pytest.importorskip("numpy")

from spectrumpy_flight.dust_composition import (
    GAIN_HIGH,
    GAIN_LOW,
    GAIN_MEDIUM,
    combine_waveform_channels,
    detect_saturation,
)


def test_detect_saturation_flags_clipped_segments_with_jitter():
    times = np.linspace(0.0, 31.5, 4096)
    base = 0.24 * np.exp(-0.5 * ((times - 6.0) / 0.8) ** 2)
    clip_level = 0.19
    signal = base.copy()
    clipped = base >= clip_level

    signal[clipped] = clip_level

    mask = detect_saturation(signal, times)

    saturated_window = (times >= 5.2) & (times <= 6.5)
    assert saturated_window.any()
    core_window = (times >= 5.4) & (times <= 6.3)
    assert core_window.any()
    assert mask[core_window].mean() > 0.8

    unsaturated_window = (times <= 4.5) | (times >= 7.0)
    assert unsaturated_window.any()
    assert not np.any(mask[unsaturated_window])


def test_combine_waveform_channels_replaces_saturated_high_with_medium():
    times = np.linspace(0.0, 31.5, 4096)
    physical_signal = 0.5 * np.exp(-0.5 * ((times - 6.0) / 0.4) ** 2)

    high_baseline = 1200.0
    medium_baseline = 75.0

    high = high_baseline + physical_signal * GAIN_HIGH
    medium = medium_baseline + physical_signal * GAIN_MEDIUM

    saturation_mask = (times >= 5.4) & (times <= 6.6)
    high[saturation_mask] = high_baseline + 0.75 * GAIN_HIGH

    combined = combine_waveform_channels(times, high, medium, None)
    assert combined is not None

    baseline_mask = times <= (times[0] + 1.0)
    high_mean = np.mean(high[baseline_mask])
    medium_mean = np.mean(medium[baseline_mask])

    expected = high - high_mean
    expected[saturation_mask] = (medium[saturation_mask] - medium_mean) * (GAIN_HIGH / GAIN_MEDIUM)

    assert np.allclose(np.mean(combined[baseline_mask]), 0.0, atol=1e-6)

    assert np.allclose(combined, expected, rtol=1e-6, atol=1e-6)


def test_combine_waveform_channels_uses_low_when_high_and_medium_saturate():
    times = np.linspace(0.0, 31.5, 2048)
    physical_signal = 0.25 + 0.45 * np.exp(-((times - 7.5) / 0.35) ** 2)

    high_baseline = 980.0
    medium_baseline = 63.0
    low_baseline = -4.5

    high = high_baseline + physical_signal * GAIN_HIGH
    medium = medium_baseline + physical_signal * GAIN_MEDIUM
    low = low_baseline + physical_signal * GAIN_LOW

    clip_mask = (times >= 6.9) & (times <= 8.2)
    high[clip_mask] = high_baseline + 0.6 * GAIN_HIGH
    medium[clip_mask] = medium_baseline + 0.6 * GAIN_MEDIUM

    combined = combine_waveform_channels(times, high, medium, low)
    assert combined is not None

    baseline_mask = times <= (times[0] + 1.0)
    high_mean = np.mean(high[baseline_mask])
    medium_mean = np.mean(medium[baseline_mask])
    low_mean = np.mean(low[baseline_mask])

    expected = high - high_mean
    expected[clip_mask] = (medium[clip_mask] - medium_mean) * (GAIN_HIGH / GAIN_MEDIUM)

    medium_clip_mask = clip_mask
    expected[medium_clip_mask] = (low[medium_clip_mask] - low_mean) * (GAIN_HIGH / GAIN_LOW)

    assert np.allclose(np.mean(combined[baseline_mask]), 0.0, atol=1e-6)

    assert np.allclose(combined, expected, rtol=1e-6, atol=1e-6)


def test_combine_waveform_channels_prefers_high_when_unsaturated():
    times = np.linspace(0.0, 31.5, 4096)
    physical_signal = 0.42 * np.exp(-0.5 * ((times - 6.0) / 0.55) ** 2)

    high_baseline = 1120.0
    medium_baseline = 74.0

    high = high_baseline + physical_signal * GAIN_HIGH
    medium = medium_baseline + 1.12 * physical_signal * GAIN_MEDIUM

    combined = combine_waveform_channels(times, high, medium, None)
    assert combined is not None

    baseline_mask = times <= (times[0] + 1.0)
    expected = high - np.mean(high[baseline_mask])

    assert np.allclose(combined, expected, rtol=1e-6, atol=1e-6)


def test_combine_waveform_channels_respects_manual_selection():
    times = np.linspace(0.0, 31.5, 1024)
    physical_signal = 0.38 * np.exp(-0.5 * ((times - 5.2) / 0.6) ** 2)

    medium_baseline = 66.0
    medium = medium_baseline + physical_signal * GAIN_MEDIUM

    combined = combine_waveform_channels(
        times,
        None,
        medium,
        None,
        enabled_channels=("TOF M",),
    )

    assert combined is not None

    baseline_mask = times <= (times[0] + 1.0)
    expected = medium - np.mean(medium[baseline_mask])

    assert np.allclose(combined, expected, rtol=1e-6, atol=1e-6)
