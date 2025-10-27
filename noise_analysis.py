"""Noise analysis window for SpectrumPY."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

try:  # pragma: no cover - prefer PyQt6
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QFont
    from PyQt6.QtWidgets import (
        QComboBox,
        QFrame,
        QGridLayout,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QMessageBox,
        QScrollArea,
        QSizePolicy,
        QStatusBar,
        QVBoxLayout,
        QWidget,
    )
except Exception:  # pragma: no cover - fallback to PySide6
    from PySide6.QtCore import Qt  # type: ignore
    from PySide6.QtGui import QFont  # type: ignore
    from PySide6.QtWidgets import (  # type: ignore
        QComboBox,
        QFrame,
        QGridLayout,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QMessageBox,
        QScrollArea,
        QSizePolicy,
        QStatusBar,
        QVBoxLayout,
        QWidget,
    )


@dataclass(frozen=True)
class ChannelMeta:
    name: str
    dataset: str
    time_dataset: str
    unit: str = ""


ChannelLoader = Callable[[str, str], Tuple[Optional[np.ndarray], Optional[np.ndarray]]]


class NoiseAnalysisWindow(QMainWindow):
    """Interactive window that provides publication-ready noise diagnostics."""

    def __init__(
        self,
        *,
        event_name: str,
        channels: Iterable[ChannelMeta],
        loader: ChannelLoader,
        event_names: Optional[Sequence[str]] = None,
        on_event_changed: Optional[Callable[[str], None]] = None,
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"Noise Analysis — Event {event_name}")
        self.setMinimumSize(960, 720)
        self.setStatusBar(QStatusBar(self))

        self._event_name = event_name
        self._loader = loader
        self._channel_definitions = list(channels)
        self._event_names: List[str] = list(dict.fromkeys(event_names or []))
        if event_name and event_name not in self._event_names:
            self._event_names.append(event_name)
        self._external_event_callback = on_event_changed
        self._event_combo: Optional[QComboBox] = None
        self._block_event_combo = False
        self._data_cache: Dict[Tuple[str, str], Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = {}
        self._channel_units: Dict[str, str] = {
            meta.name: meta.unit.strip() for meta in self._channel_definitions
        }

        self._build_ui()
        self._populate_event_selector()
        self._populate_channels()

    # ------------------------------------------------------------------
    # UI assembly
    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        self.setStyleSheet(
            """
            QMainWindow {
                background-color: #edf1fb;
            }
            QFrame#SelectionCard, QFrame#SummaryCard {
                background-color: #ffffff;
                border-radius: 18px;
                border: 1px solid #d0d4e6;
            }
            QLabel#TitleLabel {
                font-size: 26px;
                font-weight: 700;
                color: #1f2240;
            }
            QLabel#SubtitleLabel {
                font-size: 15px;
                color: #495057;
            }
            QComboBox {
                font-size: 15px;
                padding: 10px 16px;
                border-radius: 12px;
                border: 1px solid #adb5d3;
                background-color: #f8f9ff;
                color: #1f2240;
            }
            QComboBox::drop-down {
                width: 26px;
                border-left: 1px solid #adb5d3;
            }
            QLabel#SummaryLabel {
                font-size: 15px;
                color: #495057;
            }
            QLabel#SummaryValue {
                font-size: 17px;
                font-weight: 700;
                color: #212529;
            }
            """
        )

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setCentralWidget(scroll_area)

        central = QWidget(scroll_area)
        scroll_area.setWidget(central)

        layout = QVBoxLayout(central)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(18)

        selection_card = QFrame()
        selection_card.setObjectName("SelectionCard")
        selection_card_layout = QVBoxLayout(selection_card)
        selection_card_layout.setContentsMargins(28, 28, 28, 28)
        selection_card_layout.setSpacing(14)

        header_layout = QVBoxLayout()
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(6)

        self._title_label = QLabel("Noise Analysis")
        self._title_label.setObjectName("TitleLabel")
        header_layout.addWidget(self._title_label)

        self._subtitle_label = QLabel(f"Event {self._event_name}")
        self._subtitle_label.setObjectName("SubtitleLabel")
        header_layout.addWidget(self._subtitle_label)

        selection_card_layout.addLayout(header_layout)

        event_layout = QHBoxLayout()
        event_layout.setContentsMargins(0, 0, 0, 0)
        event_layout.setSpacing(12)

        event_label = QLabel("Event")
        event_label.setObjectName("SummaryLabel")
        event_layout.addWidget(event_label)

        event_combo = QComboBox()
        event_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        event_combo.setMinimumHeight(36)
        event_combo.currentIndexChanged.connect(self._on_event_combo_changed)
        event_layout.addWidget(event_combo, 1)

        event_layout.addStretch()
        selection_card_layout.addLayout(event_layout)
        self._event_combo = event_combo

        combo_layout = QHBoxLayout()
        combo_layout.setContentsMargins(0, 0, 0, 0)
        combo_layout.setSpacing(12)

        combo_label = QLabel("Channel")
        combo_label.setObjectName("SummaryLabel")
        combo_layout.addWidget(combo_label)

        self.channel_combo = QComboBox()
        self.channel_combo.currentIndexChanged.connect(self._update_analysis)
        combo_layout.addWidget(self.channel_combo, 1)

        combo_layout.addStretch()
        selection_card_layout.addLayout(combo_layout)

        layout.addWidget(selection_card)

        self.figure = Figure(figsize=(11, 8), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        grid = self.figure.add_gridspec(3, 1, height_ratios=[1, 1, 1])
        self.axes_hist = self.figure.add_subplot(grid[0])
        self.axes_power = self.figure.add_subplot(grid[1])
        self.axes_autocorr = self.figure.add_subplot(grid[2])

        self.nav_toolbar = NavigationToolbar2QT(self.canvas, self)
        self.nav_toolbar.setStyleSheet("font-size: 13px; padding: 6px;")
        layout.addWidget(self.nav_toolbar)
        layout.addWidget(self.canvas, 1)

        summary_card = QFrame()
        summary_card.setObjectName("SummaryCard")
        summary_card_layout = QGridLayout(summary_card)
        summary_card_layout.setContentsMargins(24, 24, 24, 24)
        summary_card_layout.setHorizontalSpacing(30)
        summary_card_layout.setVerticalSpacing(12)

        summary_items = [
            ("Gaussian Mean", "mean"),
            ("Gaussian σ", "sigma"),
            ("Gaussian Amplitude", "amplitude"),
            ("Poisson Noise", "poisson"),
            ("RMS Noise", "rms"),
        ]

        self._summary_labels: Dict[str, QLabel] = {}
        self._summary_title_labels: Dict[str, QLabel] = {}
        self._summary_titles: Dict[str, str] = {}
        for idx, (title, key) in enumerate(summary_items):
            label = QLabel(title)
            label.setObjectName("SummaryLabel")
            value = QLabel("–")
            value.setObjectName("SummaryValue")
            summary_card_layout.addWidget(label, idx // 3, (idx % 3) * 2)
            summary_card_layout.addWidget(value, idx // 3, (idx % 3) * 2 + 1)
            self._summary_labels[key] = value
            self._summary_title_labels[key] = label
            self._summary_titles[key] = title

        layout.addWidget(summary_card)

    # ------------------------------------------------------------------
    # Channel management
    # ------------------------------------------------------------------
    def _populate_event_selector(self) -> None:
        combo = self._event_combo
        if combo is None:
            return
        combo.blockSignals(True)
        combo.clear()
        if self._event_name and self._event_name not in self._event_names:
            self._event_names.append(self._event_name)
        for name in self._event_names:
            combo.addItem(name)
        index = 0
        if self._event_name and self._event_name in self._event_names:
            index = self._event_names.index(self._event_name)
        combo.setCurrentIndex(index)
        combo.blockSignals(False)
        combo.setEnabled(combo.count() > 1)

    def _ensure_event_listed(self, event_name: str) -> None:
        if not event_name:
            return
        if event_name not in self._event_names:
            self._event_names.append(event_name)
            if self._event_combo is not None:
                self._event_combo.addItem(event_name)
        if self._event_combo is not None:
            self._event_combo.setEnabled(self._event_combo.count() > 1)

    def _sync_event_combo(self, event_name: Optional[str]) -> None:
        combo = self._event_combo
        if combo is None or not event_name:
            return
        if event_name not in self._event_names:
            self._ensure_event_listed(event_name)
        try:
            index = self._event_names.index(event_name)
        except ValueError:
            return
        if combo.currentIndex() != index:
            self._block_event_combo = True
            try:
                combo.setCurrentIndex(index)
            finally:
                self._block_event_combo = False
        combo.setEnabled(combo.count() > 1)

    def _on_event_combo_changed(self, index: int) -> None:
        if self._block_event_combo:
            return
        if 0 <= index < len(self._event_names):
            self._switch_event(self._event_names[index], from_user=True)

    def _switch_event(self, event_name: str, *, from_user: bool) -> None:
        if not event_name:
            return
        self._ensure_event_listed(event_name)
        previous = self._event_name
        if event_name == previous:
            self._sync_event_combo(event_name)
            return

        self._event_name = event_name
        self.setWindowTitle(f"Noise Analysis — Event {event_name}")
        self._subtitle_label.setText(f"Event {event_name}")
        self._sync_event_combo(event_name)
        previous_channel = self.channel_combo.currentText().strip()
        self._populate_channels(preferred=previous_channel, trigger_update=False)
        try:
            self.statusBar().showMessage(f"Viewing event {event_name}", 4000)
        except Exception:
            pass

        self._update_analysis()

        if from_user and self._external_event_callback is not None:
            try:
                self._external_event_callback(event_name)
            except Exception:
                pass

    def set_current_event(self, event_name: Optional[str]) -> None:
        if not event_name:
            return
        self._switch_event(event_name, from_user=False)

    def _populate_channels(
        self,
        *,
        preferred: Optional[str] = None,
        trigger_update: bool = True,
    ) -> None:
        self.channel_combo.blockSignals(True)
        self.channel_combo.clear()

        available = [meta.name for meta in self._channel_definitions if self._has_data(meta.name)]
        for name in available:
            self.channel_combo.addItem(name)
        self.channel_combo.blockSignals(False)

        self.channel_combo.setEnabled(bool(available))

        if not available:
            QMessageBox.information(
                self,
                "No Data",
                "None of the standard channels contain data for this event.",
            )
            self._clear_axes("No channel data available")
            return

        index = 0
        if preferred and preferred in available:
            index = available.index(preferred)
        self.channel_combo.setCurrentIndex(index)

        if trigger_update:
            self._update_analysis()

    def _has_data(self, channel: str) -> bool:
        values, _ = self._get_channel_arrays(channel)
        if values is None:
            return False
        flat = np.asarray(values, dtype=float).ravel()
        flat = flat[np.isfinite(flat)]
        return flat.size > 0

    def _get_channel_arrays(self, channel: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        key = (self._event_name, channel)
        if key not in self._data_cache:
            self._data_cache[key] = self._loader(self._event_name, channel)
        return self._data_cache[key]

    # ------------------------------------------------------------------
    # Analysis & plotting
    # ------------------------------------------------------------------
    def _update_analysis(self) -> None:
        channel = self.channel_combo.currentText().strip()
        if not channel:
            return

        values, times = self._get_channel_arrays(channel)
        if values is None:
            self._clear_axes(f"{channel}: data unavailable")
            return

        signal = np.asarray(values, dtype=float).ravel()
        signal = signal[np.isfinite(signal)]
        if signal.size == 0:
            self._clear_axes(f"{channel}: empty dataset")
            return

        mean = float(np.mean(signal))
        sigma = float(np.std(signal))
        total = signal.size

        # Histogram and Gaussian fit
        self.axes_hist.cla()
        bins = min(80, max(20, total // 20))
        counts, edges = np.histogram(signal, bins=bins)
        centers = (edges[:-1] + edges[1:]) / 2.0
        widths = np.diff(edges)
        if centers.size:
            self.axes_hist.bar(
                centers,
                counts,
                width=widths,
                align="center",
                color="#748ffc",
                edgecolor="#ffffff",
                alpha=0.85,
                label="Amplitude histogram",
            )

        amplitude = 0.0
        if sigma > 0 and centers.size:
            bin_width = float(np.mean(widths)) if widths.size else 1.0
            amplitude = float(total * bin_width / (sigma * np.sqrt(2.0 * np.pi)))
            fit_y = amplitude * np.exp(-0.5 * ((centers - mean) / sigma) ** 2)
            self.axes_hist.plot(
                centers,
                fit_y,
                color="#364fc7",
                linewidth=2.4,
                label="Gaussian fit",
            )
            self.axes_hist.legend(loc="upper right")

        unit = self._channel_units.get(channel, "")
        amplitude_label = "Amplitude"
        if unit:
            amplitude_label = f"{amplitude_label} [{unit}]"
        self.axes_hist.set_title(f"{channel} Amplitude Distribution", fontsize=16, fontweight="bold")
        self.axes_hist.set_xlabel(amplitude_label)
        self.axes_hist.set_ylabel("Counts")
        self.axes_hist.set_facecolor("#f8f9fb")
        self.axes_hist.grid(True, alpha=0.3)
        self.axes_hist.tick_params(axis="both", labelsize=12, width=1.2, length=6)

        # Power spectrum
        self.axes_power.cla()
        detrended = signal - mean
        fft = np.fft.rfft(detrended)
        magnitude = np.abs(fft)

        dt = self._infer_sample_spacing(times)
        freqs = np.fft.rfftfreq(signal.size, d=dt if dt and dt > 0 else 1.0)
        if freqs.size and magnitude.size:
            positive_mask = (freqs > 0) & (magnitude > 0)
            if np.count_nonzero(positive_mask) >= 1:
                freqs = freqs[positive_mask]
                magnitude = magnitude[positive_mask]
            else:
                freqs = np.array([])
                magnitude = np.array([])

        if freqs.size and magnitude.size:
            self.axes_power.plot(freqs, magnitude, color="#4263eb", linewidth=2.0)
            if magnitude.size > 1:
                self.axes_power.fill_between(
                    freqs,
                    magnitude,
                    np.full_like(magnitude, magnitude.min()),
                    color="#bac8ff",
                    alpha=0.45,
                )
            self.axes_power.set_xscale("log")
            self.axes_power.set_yscale("log")
        self.axes_power.set_title(f"{channel} Power Spectrum", fontsize=16, fontweight="bold")
        self.axes_power.set_xlabel("Frequency (1/Δt)")
        self.axes_power.set_ylabel("|FFT|")
        self.axes_power.set_facecolor("#f8f9fb")
        self.axes_power.grid(True, alpha=0.3)
        self.axes_power.tick_params(axis="both", labelsize=12, width=1.2, length=6)

        # Autocorrelation function
        self.axes_autocorr.cla()
        autocorr_full = np.correlate(detrended, detrended, mode="full")
        autocorr = autocorr_full[autocorr_full.size // 2 :]
        if autocorr.size:
            normalization = np.arange(signal.size, 0, -1, dtype=float)
            autocorr = autocorr[: normalization.size] / normalization
            variance = float(np.var(detrended))
            if variance > 0:
                autocorr = autocorr / variance

        if autocorr.size:
            if dt and dt > 0:
                lags = np.arange(autocorr.size) * dt
                xlabel = "Lag (time)"
            else:
                lags = np.arange(autocorr.size)
                xlabel = "Lag (samples)"
            self.axes_autocorr.plot(lags, autocorr, color="#2f9e44", linewidth=2.0)
            self.axes_autocorr.set_xlim(left=0)
        else:
            xlabel = "Lag"

        self.axes_autocorr.set_title(f"{channel} Autocorrelation", fontsize=16, fontweight="bold")
        self.axes_autocorr.set_xlabel(xlabel)
        self.axes_autocorr.set_ylabel("Autocorr")
        self.axes_autocorr.set_facecolor("#f8f9fb")
        self.axes_autocorr.grid(True, alpha=0.3)
        self.axes_autocorr.tick_params(axis="both", labelsize=12, width=1.2, length=6)

        poisson = float(np.sqrt(abs(np.mean(np.clip(signal, a_min=0.0, a_max=None)))))
        rms = float(np.sqrt(np.mean(detrended**2)))

        self._apply_summary_units(unit)
        self._set_summary_value("mean", mean)
        self._set_summary_value("sigma", sigma)
        self._set_summary_value("amplitude", amplitude)
        self._set_summary_value("poisson", poisson)
        self._set_summary_value("rms", rms)

        sample_text = f"{channel} • {total} samples"
        if dt and dt > 0:
            sample_text += f" • Δt≈{dt:.3g}"
        self.statusBar().showMessage(sample_text)

        self.canvas.draw_idle()

    def _infer_sample_spacing(self, times: Optional[np.ndarray]) -> Optional[float]:
        if times is None:
            return None
        flat = np.asarray(times, dtype=float).ravel()
        flat = flat[np.isfinite(flat)]
        if flat.size < 2:
            return None
        diffs = np.diff(np.sort(flat))
        diffs = diffs[diffs > 0]
        if diffs.size == 0:
            return None
        return float(np.median(diffs))

    def _set_summary_value(self, key: str, value: float) -> None:
        label = self._summary_labels.get(key)
        if not label:
            return
        if not np.isfinite(value):
            label.setText("–")
            return
        if abs(value) >= 1e4 or (abs(value) < 1e-2 and value != 0.0):
            label.setText(f"{value:.3e}")
        else:
            label.setText(f"{value:.4f}")

    def _clear_axes(self, message: str) -> None:
        self.axes_hist.cla()
        self.axes_power.cla()
        self.axes_autocorr.cla()

        for ax in (self.axes_hist, self.axes_power, self.axes_autocorr):
            ax.set_facecolor("#f8f9fb")
            ax.grid(False)
            ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, fontsize=15)
            ax.set_xticks([])
            ax.set_yticks([])

        for value in self._summary_labels.values():
            value.setText("–")

        self._apply_summary_units("")

        self.statusBar().showMessage(message)
        self.canvas.draw_idle()

    def _apply_summary_units(self, unit: str) -> None:
        for key, label in self._summary_title_labels.items():
            base = self._summary_titles.get(key, label.text())
            if unit:
                label.setText(f"{base} [{unit}]")
            else:
                label.setText(base)


def launch_noise_analysis_window(
    *,
    event_name: str,
    channels: Iterable[ChannelMeta],
    loader: ChannelLoader,
    event_names: Optional[Sequence[str]] = None,
    on_event_changed: Optional[Callable[[str], None]] = None,
    parent: Optional[QWidget] = None,
) -> NoiseAnalysisWindow:
    """Create a :class:`NoiseAnalysisWindow` ready for display."""

    return NoiseAnalysisWindow(
        event_name=event_name,
        channels=channels,
        loader=loader,
        event_names=event_names,
        on_event_changed=on_event_changed,
        parent=parent,
    )
