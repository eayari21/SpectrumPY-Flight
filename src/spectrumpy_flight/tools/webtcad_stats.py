#!/usr/bin/env python3
"""
IDEX WebTCAD downloader + viewer (PyQt6)

- Prompts for username / password
- Lets the user choose start and stop times (UTC)
- Downloads all requested time series from WebTCAD
- Saves a merged CSV (one time column + all series)
- Displays the time series in a multi-panel matplotlib view
  with associated MIN / MAX channels overlaid on the same axes.

Dependencies:
    pip install PyQt6 requests pandas matplotlib
"""

import io
import os
import sys
from pathlib import Path
from dataclasses import dataclass
from urllib.parse import quote
from typing import Dict, Tuple, List

import pandas as pd
import requests
import matplotlib

# Allow running directly from a source checkout without installing the package
REPO_SRC = Path(__file__).resolve().parents[3] / "src"
if REPO_SRC.is_dir() and str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

matplotlib.use("QtAgg")  # Qt backend for PyQt6

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import matplotlib.dates as mdates

from PyQt6.QtCore import Qt, QDateTime
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QLabel,
    QLineEdit,
    QDateTimeEdit,
    QPushButton,
    QMessageBox,
    QGroupBox,
    QCheckBox,
    QScrollArea,
)

from spectrumpy_flight.plot_style import apply_plot_style

# Apply publication-style defaults shared with IDEX Quicklook
apply_plot_style()


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

BASE_URL = (
    "https://lasp.colorado.edu/ops/imap/webtcad/latis/dap/"
    "AnalogTelemetryItem_SID1.csv"
)

@dataclass
class SeriesSpec:
    """Identifier for a single analog telemetry series."""
    query_param: str
    value: str


# Individual series: label -> SeriesSpec
def _load_series_specs() -> tuple[Dict[str, SeriesSpec], List[str]]:
    specs: Dict[str, SeriesSpec] = {
        "TOFHGMIN": SeriesSpec("KEY", "IDX_HW.SCITOFHGMIN"),
        "Det_Voltage": SeriesSpec("TMID", "38"),
        "Target_Low_Gain_Min": SeriesSpec("KEY", "IDX_HW.SCITLRMIN"),
        "Target_High_Gain_Min": SeriesSpec("KEY", "IDX_HW.SCITHRMIN"),
        "Target_Low_Gain_Max": SeriesSpec("KEY", "IDX_HW.SCITLRMX"),
        "Target_High_Gain_Max": SeriesSpec("KEY", "IDX_HW.SCITHRMX"),
        "TOFLGMIN": SeriesSpec("KEY", "IDX_HW.SCITOFLGMIN"),
        "TOFMGMIN": SeriesSpec("KEY", "IDX_HW.SCITOFMGMIN"),
        "TOFHGMAX": SeriesSpec("KEY", "IDX_HW.SCITOFHGMX"),
        "TOFLGMAX": SeriesSpec("KEY", "IDX_HW.SCITOFLGMX"),
        "TOFMGMAX": SeriesSpec("KEY", "IDX_HW.SCITOFMGMX"),
        "IONGRIDMIN": SeriesSpec("KEY", "IDX_HW.SCIIONMIN"),
        "IONGRIDMAX": SeriesSpec("KEY", "IDX_HW.SCIIONMX"),
        "TOFHG_THRESHOLD": SeriesSpec("KEY", "IDX_HW.TRGTOFHGTHOLD"),
        "TOFLG_THRESHOLD": SeriesSpec("KEY", "IDX_HW.TRGTOFLGTHOLD"),
        "TOFMG_THRESHOLD": SeriesSpec("KEY", "IDX_HW.TRGTOFMGTHOLD"),
        "ADC_MODE_TOFLG": SeriesSpec("KEY", "IDX_HW.SCIADCMODETOFLG"),
        "ADC_MODE_TOFHG": SeriesSpec("KEY", "IDX_HW.SCIADCMODETOFHG"),
        "SCI_EVENT_CNT": SeriesSpec("KEY", "IDX_SW.SCISEVT_CNT"),
        "SCI_AID": SeriesSpec("KEY", "IDX_SW.SCIAID"),
    }

    warnings: List[str] = []

    optional_env = {
        "IONGRIDMIN": ("WEBTCAD_IONGRID_MIN_TMID", "WEBTCAD_IONGRID_MIN_KEY"),
        "IONGRIDMAX": ("WEBTCAD_IONGRID_MAX_TMID", "WEBTCAD_IONGRID_MAX_KEY"),
    }

    for key, (tmid_var, key_var) in optional_env.items():
        raw_tmid = os.getenv(tmid_var)
        raw_key = os.getenv(key_var)

        if raw_tmid not in (None, ""):
            specs[key] = SeriesSpec("TMID", raw_tmid)
            continue

        if raw_key not in (None, ""):
            specs[key] = SeriesSpec("KEY", raw_key)
            continue

    return specs, warnings


SERIES_SPECS, SERIES_CONFIG_WARNINGS = _load_series_specs()

@dataclass(frozen=True)
class PlotSeries:
    key: str
    label: str
    linestyle: str = "-"
    linewidth: float = 0.8


@dataclass(frozen=True)
class PlotGroup:
    title: str
    series: List[PlotSeries]


PLOT_GROUPS: List[PlotGroup] = [
    PlotGroup(
        "TOF L",
        [
            PlotSeries("TOFLGMIN", "Min"),
            PlotSeries("TOFLGMAX", "Max", linestyle="--"),
        ],
    ),
    PlotGroup(
        "TOF M",
        [
            PlotSeries("TOFMGMIN", "Min"),
            PlotSeries("TOFMGMAX", "Max", linestyle="--"),
        ],
    ),
    PlotGroup(
        "TOF H",
        [
            PlotSeries("TOFHGMIN", "Min"),
            PlotSeries("TOFHGMAX", "Max", linestyle="--"),
        ],
    ),
    PlotGroup(
        "Target L",
        [
            PlotSeries("Target_Low_Gain_Min", "Min"),
            PlotSeries("Target_Low_Gain_Max", "Max", linestyle="--"),
        ],
    ),
    PlotGroup(
        "Target H",
        [
            PlotSeries("Target_High_Gain_Min", "Min"),
            PlotSeries("Target_High_Gain_Max", "Max", linestyle="--"),
        ],
    ),
    PlotGroup(
        "Ion Grid",
        [
            PlotSeries("IONGRIDMIN", "Min"),
            PlotSeries("IONGRIDMAX", "Max", linestyle="--"),
        ],
    ),
    PlotGroup("TOF H Threshold", [PlotSeries("TOFHG_THRESHOLD", "Threshold")]),
    PlotGroup("TOF L Threshold", [PlotSeries("TOFLG_THRESHOLD", "Threshold")]),
    PlotGroup("TOF M Threshold", [PlotSeries("TOFMG_THRESHOLD", "Threshold")]),
    PlotGroup(
        "ADC Mode",
        [
            PlotSeries("ADC_MODE_TOFLG", "TOF LG"),
            PlotSeries("ADC_MODE_TOFHG", "TOF HG", linestyle="--"),
        ],
    ),
    PlotGroup("SCI Event Count", [PlotSeries("SCI_EVENT_CNT", "Count")]),
    PlotGroup("SCI AID", [PlotSeries("SCI_AID", "AID")]),
]

DETECTOR_VOLT_KEY = "Det_Voltage"


@dataclass
class SeriesData:
    """Container for a single downloaded series."""
    name: str
    df: pd.DataFrame          # columns: ["time", value_name, "time_dt"]
    value_col: str


# ----------------------------------------------------------------------
# Networking + data helpers
# ----------------------------------------------------------------------

def _query_value(spec: SeriesSpec) -> str:
    """
    Prepare a query value for the series.

    String-based KEY filters need quoting to be parsed correctly by the
    WebTCAD service; numeric TMIDs can be passed as-is.
    """
    if spec.query_param.upper() != "KEY":
        return spec.value

    if spec.value.startswith('"') and spec.value.endswith('"'):
        return spec.value

    return f'"{spec.value}"'


def build_url(spec: SeriesSpec, start_iso: str, stop_iso: str) -> str:
    """
    Construct a URL matching the format you supplied, with
    time>=start_iso and time<=stop_iso, and the formatted time column.
    """
    params = {
        spec.query_param: _query_value(spec),
        "time>=": start_iso,
        "time<=": stop_iso,
    }

    req = requests.PreparedRequest()
    req.prepare_url(BASE_URL, params)

    format_param = quote(
        "format_time(time,yyyy-DDD'T'HH:mm:ss.SSS)",
        safe="(),",
    )
    return f"{req.url}&{format_param}"


def fetch_series(
    session: requests.Session,
    series_key: str,
    spec: SeriesSpec,
    start_iso: str,
    stop_iso: str,
) -> SeriesData:
    """Download one series from WebTCAD and return as SeriesData."""
    url = build_url(spec, start_iso, stop_iso)
    resp = session.get(url)
    resp.raise_for_status()

    buf = io.StringIO(resp.text)
    df_raw = pd.read_csv(buf)

    if df_raw.empty:
        raise ValueError(f"No data for {series_key} in given time range.")

    # Assume first column is time string, second column is the value
    time_col = df_raw.columns[0]
    value_col = df_raw.columns[1]

    df = df_raw[[time_col, value_col]].copy()
    df.rename(columns={time_col: "time"}, inplace=True)

    # Parse time for plotting (yyyy-DDD'T'HH:mm:ss.SSS)
    df["time_dt"] = pd.to_datetime(
        df["time"],
        format="%Y-%jT%H:%M:%S.%f",
        errors="coerce",
    )

    return SeriesData(name=series_key, df=df, value_col=value_col)


def merge_to_single_csv(series_map: Dict[str, SeriesData], out_path: str) -> pd.DataFrame:
    """
    Merge all series on string 'time' column (outer join) and save
    a CSV with one 'time' column and one column per series.
    """
    combined_df = None

    for key, s in series_map.items():
        df_merge = s.df[["time", s.value_col]].copy()
        df_merge.rename(columns={s.value_col: key}, inplace=True)

        if combined_df is None:
            combined_df = df_merge
        else:
            combined_df = pd.merge(
                combined_df,
                df_merge,
                on="time",
                how="outer",
            )

    if combined_df is None or combined_df.empty:
        raise ValueError("No data to merge.")

    # Sort by time using parsed datetime (fall back to lexicographic order)
    try:
        t = pd.to_datetime(
            combined_df["time"],
            format="%Y-%jT%H:%M:%S.%f",
            errors="coerce",
        )
        combined_df["_sort"] = t
        combined_df.sort_values("_sort", inplace=True)
        combined_df.drop(columns=["_sort"], inplace=True)
    except Exception:
        combined_df.sort_values("time", inplace=True)

    # Reorder columns: time first, then sorted series names for stability
    cols = ["time"] + sorted([c for c in combined_df.columns if c != "time"])
    combined_df = combined_df[cols]

    combined_df.to_csv(out_path, index=False)
    return combined_df


# ----------------------------------------------------------------------
# Matplotlib canvas
# ----------------------------------------------------------------------

class TimeSeriesCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(9, 6), tight_layout=False)
        super().__init__(fig)
        self.setParent(parent)
        self.figure = fig
        self.axes: list = []
        self.shared_xlim: Tuple[float, float] | None = None
        self._updating_xlim = False

        # Match the publication-style background used across Quicklook tools
        self.figure.set_facecolor(matplotlib.rcParams.get("figure.facecolor", "white"))

    def wheelEvent(self, event):
        """Route wheel events to the scroll area instead of zooming the plot."""

        # Let users hold Ctrl to retain matplotlib's zoom behavior, but otherwise
        # pass the wheel event up so the surrounding QScrollArea can scroll.
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            super().wheelEvent(event)
        else:
            event.ignore()

    def reset_shared_xlim(self) -> None:
        """Clear any stored zoom state so new plots start autoscaled."""

        self.shared_xlim = None

    def plot_series(
        self,
        series_map: Dict[str, SeriesData],
        visible_groups: Dict[str, bool],
    ) -> None:
        self.figure.clear()

        valid_groups = [
            grp
            for grp in PLOT_GROUPS
            if visible_groups.get(grp.title, True)
            and any(series.key in series_map for series in grp.series)
        ]

        detector_enabled = visible_groups.get("Detector Voltage", True)
        nrows = len(valid_groups)
        if detector_enabled and DETECTOR_VOLT_KEY in series_map:
            nrows += 1

        if nrows == 0:
            self.draw()
            return

        axes = self.figure.subplots(nrows, 1, sharex=True)
        if nrows == 1:
            axes = [axes]
        self.axes = axes

        ax_idx = 0

        # Plot all min/max pairs
        for group in valid_groups:
            ax = axes[ax_idx]
            ax_idx += 1

            plotted = 0
            for series in group.series:
                if series.key not in series_map:
                    continue
                s = series_map[series.key]
                df = s.df
                ax.plot(
                    df["time_dt"],
                    df[s.value_col],
                    linewidth=series.linewidth,
                    linestyle=series.linestyle,
                    label=series.label,
                )
                plotted += 1

            ax.set_ylabel(group.title)
            ax.grid(True, linestyle="-", linewidth=0.8, alpha=0.35)
            if plotted > 1:
                ax.legend(loc="upper right")
            ax.tick_params(axis="both", labelsize=11)

        # Detector voltage in the last panel
        if detector_enabled and DETECTOR_VOLT_KEY in series_map:
            ax = axes[ax_idx]
            s = series_map[DETECTOR_VOLT_KEY]
            df = s.df
            ax.plot(df["time_dt"], df[s.value_col], linewidth=1.2)
            ax.set_ylabel("Det V")
            ax.grid(True, linestyle="-", linewidth=0.8, alpha=0.35)
            ax.tick_params(axis="both", labelsize=11)

        # Shared x-label on the last axis
        axes[-1].set_xlabel("Time [UTC]", fontsize=13)
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n%H:%M:%S"))

        for ax in axes:
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            ax.margins(x=0)
            ax.tick_params(axis="x", rotation=0)

        # Apply any previous x-zoom (e.g., after toggling visibility) and keep
        # all panels in sync when the user zooms or pans.
        if self.shared_xlim is not None:
            for ax in axes:
                ax.set_xlim(self.shared_xlim)

        for ax in axes:
            ax.callbacks.connect("xlim_changed", self._on_xlim_changed)

        self.figure.tight_layout(pad=1.2)
        # Increase canvas height so scroll areas have room to render all panels
        self.setMinimumHeight(max(1, nrows) * 260)
        self.draw()

    def _on_xlim_changed(self, ax):
        if self._updating_xlim:
            return

        self.shared_xlim = ax.get_xlim()
        self._updating_xlim = True
        try:
            for other_ax in self.axes:
                if other_ax is not ax:
                    other_ax.set_xlim(self.shared_xlim)
        finally:
            self._updating_xlim = False


# ----------------------------------------------------------------------
# Main window
# ----------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("IDEX WebTCAD Waveform Viewer")
        self.resize(1100, 700)
        self.current_series_map: Dict[str, SeriesData] = {}
        self.current_warnings: List[str] = []

        central = QWidget()
        self.setCentralWidget(central)

        main_layout = QVBoxLayout(central)

        # --- Controls group ---------------------------------------------------
        controls_group = QGroupBox("Connection && Time Range")
        controls_layout = QVBoxLayout(controls_group)

        # Credentials + times in a form layout
        form = QFormLayout()

        self.user_edit = QLineEdit()
        self.user_edit.setPlaceholderText("username")

        self.pass_edit = QLineEdit()
        self.pass_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self.pass_edit.setPlaceholderText("password")

        self.start_edit = QDateTimeEdit()
        self.start_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss 'UTC'")
        self.start_edit.setCalendarPopup(True)
        self.start_edit.setTimeSpec(Qt.TimeSpec.UTC)

        self.stop_edit = QDateTimeEdit()
        self.stop_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss 'UTC'")
        self.stop_edit.setCalendarPopup(True)
        self.stop_edit.setTimeSpec(Qt.TimeSpec.UTC)

        # Defaults: last 24 hours in UTC
        now_utc = QDateTime.currentDateTimeUtc()
        self.stop_edit.setDateTime(now_utc)
        self.start_edit.setDateTime(now_utc.addDays(-1))

        form.addRow("Username:", self.user_edit)
        form.addRow("Password:", self.pass_edit)
        form.addRow("Start time (UTC):", self.start_edit)
        form.addRow("Stop time (UTC):", self.stop_edit)

        controls_layout.addLayout(form)

        # Run button
        button_row = QHBoxLayout()
        button_row.addStretch(1)

        self.run_btn = QPushButton("Download && Plot")
        self.run_btn.setFixedHeight(36)
        self.run_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #2563EB;
                color: white;
                border-radius: 6px;
                padding: 6px 18px;
                font-size: 14px;
                font-weight: 500;
            }
            QPushButton:hover {
                background-color: #1D4ED8;
            }
            QPushButton:pressed {
                background-color: #1E40AF;
            }
            """
        )
        self.run_btn.clicked.connect(self.on_run_clicked)

        button_row.addWidget(self.run_btn)
        button_row.addStretch(1)

        controls_layout.addLayout(button_row)
        main_layout.addWidget(controls_group)

        # --- Visibility toggles ----------------------------------------------
        visibility_group = QGroupBox("Plot visibility")
        visibility_layout = QHBoxLayout(visibility_group)

        self.visibility_checks: Dict[str, QCheckBox] = {}
        for group in PLOT_GROUPS:
            cb = QCheckBox(group.title)
            cb.setChecked(True)
            cb.stateChanged.connect(self.on_visibility_changed)
            self.visibility_checks[group.title] = cb
            visibility_layout.addWidget(cb)

        det_cb = QCheckBox("Detector Voltage")
        det_cb.setChecked(True)
        det_cb.stateChanged.connect(self.on_visibility_changed)
        self.visibility_checks["Detector Voltage"] = det_cb
        visibility_layout.addWidget(det_cb)

        visibility_layout.addStretch(1)
        main_layout.addWidget(visibility_group)

        # --- Plot canvas + toolbar -------------------------------------------
        self.canvas = TimeSeriesCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)

        canvas_container = QWidget()
        canvas_layout = QVBoxLayout(canvas_container)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        canvas_layout.addWidget(self.toolbar)
        canvas_layout.addWidget(self.canvas)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setWidget(canvas_container)

        plot_layout = QVBoxLayout()
        plot_layout.addWidget(scroll)

        main_layout.addLayout(plot_layout, stretch=1)

        # Status bar
        self.statusBar().showMessage("Ready")

    # ------------------------------------------------------------------ #
    # Event handlers
    # ------------------------------------------------------------------ #

    def on_run_clicked(self):
        username = self.user_edit.text().strip()
        password = self.pass_edit.text()

        if not username or not password:
            QMessageBox.warning(self, "Missing credentials", "Enter username and password.")
            return

        # Get UTC datetimes and format like 2025-11-04T18:55:16.000Z
        start_dt = self.start_edit.dateTime().toUTC()
        stop_dt = self.stop_edit.dateTime().toUTC()

        if start_dt >= stop_dt:
            QMessageBox.warning(self, "Invalid time range", "Start time must be before stop time.")
            return

        start_iso = start_dt.toString("yyyy-MM-dd'T'HH:mm:ss") + ".000Z"
        stop_iso = stop_dt.toString("yyyy-MM-dd'T'HH:mm:ss") + ".000Z"

        self.run_btn.setEnabled(False)
        self.statusBar().showMessage("Downloading data from WebTCAD…")

        try:
            series_map, warnings = self.download_all_series(username, password, start_iso, stop_iso)
        except Exception as e:
            self.run_btn.setEnabled(True)
            self.statusBar().showMessage("Error")
            QMessageBox.critical(self, "Download error", str(e))
            return

        # Save merged CSV
        try:
            csv_dir = Path(__file__).resolve().parents[3] / "CSV"
            csv_dir.mkdir(parents=True, exist_ok=True)
            out_path = (
                csv_dir
                / f"idex_waveforms_{start_dt.toString('yyyyMMdd_HHmmss')}"
                f"_to_{stop_dt.toString('yyyyMMdd_HHmmss')}.csv"
            )
            merge_to_single_csv(series_map, str(out_path))
        except Exception as e:
            # Non-fatal for plotting
            QMessageBox.warning(self, "CSV error", f"Could not write CSV:\n{e}")
            out_path = None

        # Plot
        self.current_series_map = series_map
        self.current_warnings = warnings
        self.canvas.reset_shared_xlim()
        self.refresh_plots()

        msg = "Finished"
        if out_path:
            msg += f" – saved merged CSV to {out_path}"
        self.statusBar().showMessage(msg)

        if self.current_warnings:
            warn_text = "\n".join(self.current_warnings)
            QMessageBox.warning(
                self,
                "Completed with warnings",
                warn_text,
            )

        self.run_btn.setEnabled(True)

    # ------------------------------------------------------------------ #

    def download_all_series(
        self,
        username: str,
        password: str,
        start_iso: str,
        stop_iso: str,
    ) -> Tuple[Dict[str, SeriesData], List[str]]:
        """Download configured series, skipping failures with warnings."""
        session = requests.Session()
        session.auth = (username, password)

        series_map: Dict[str, SeriesData] = {}
        warnings: List[str] = list(SERIES_CONFIG_WARNINGS)

        for key, spec in SERIES_SPECS.items():
            try:
                self.statusBar().showMessage(
                    f"Fetching {key} ({spec.query_param}={spec.value})…"
                )
                QApplication.processEvents()  # keep UI responsive
                s = fetch_series(session, key, spec, start_iso, stop_iso)
                series_map[key] = s
            except Exception as e:
                warnings.append(
                    f"Skipped {key} ({spec.query_param}={spec.value}): {e}"
                )

        if not series_map:
            raise RuntimeError("No series downloaded.")

        return series_map, warnings

    # ------------------------------------------------------------------ #

    def visible_flags(self) -> Dict[str, bool]:
        return {name: cb.isChecked() for name, cb in self.visibility_checks.items()}

    def refresh_plots(self) -> None:
        if not self.current_series_map:
            return
        self.canvas.plot_series(self.current_series_map, self.visible_flags())

    def on_visibility_changed(self) -> None:
        self.refresh_plots()


# ----------------------------------------------------------------------
# Entrypoint
# ----------------------------------------------------------------------

def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
