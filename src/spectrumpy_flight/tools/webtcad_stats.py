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
import sys
from dataclasses import dataclass
from typing import Dict, Tuple, List

import pandas as pd
import requests
import matplotlib

matplotlib.use("QtAgg")  # Qt backend for PyQt6

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

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
)


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

BASE_URL = (
    "https://lasp.colorado.edu/ops/imap/webtcad/latis/dap/"
    "AnalogTelemetryItem_SID1.csv"
)

# Individual series: label -> TMID
SERIES_TMIDS: Dict[str, int] = {
    "TOFHGMIN": 3323,
    "Det_Voltage": 38,
    "Target_Low_Gain_Min": 3356,
    "Target_High_Gain_Min": 3349,
    "Target_Low_Gain_Max": 3357,
    "Target_High_Gain_Max": 3350,
    "TOFLGMIN": 3331,
    "TOFMGMIN": 3341,
    "TOFHGMAX": 3325,
    "TOFLGMAX": 3333,
    "TOFMGMAX": 3343,
    # Add Ion Grid min/max TMIDs here if desired:
    # "IONGRIDMIN": <TMID>,
    # "IONGRIDMAX": <TMID>,
}

# How the series are grouped for plotting:
# (axis title, min_key, max_key)
PLOT_GROUPS: List[Tuple[str, str, str]] = [
    ("TOF L", "TOFLGMIN", "TOFLGMAX"),
    ("TOF M", "TOFMGMIN", "TOFMGMAX"),
    ("TOF H", "TOFHGMIN", "TOFHGMAX"),
    ("Target L", "Target_Low_Gain_Min", "Target_Low_Gain_Max"),
    ("Target H", "Target_High_Gain_Min", "Target_High_Gain_Max"),
    # ("Ion Grid", "IONGRIDMIN", "IONGRIDMAX"),
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

def build_url(tmid: int, start_iso: str, stop_iso: str) -> str:
    """
    Construct a URL matching the format you supplied, with
    time>=start_iso and time<=stop_iso, and the formatted time column.
    """
    return (
        f"{BASE_URL}?TMID={tmid}"
        f"&time%3E={start_iso}"
        f"&time%3C={stop_iso}"
        f"&format_time(yyyy-DDD%27T%27HH:mm:ss.SSS)"
    )


def fetch_series(
    session: requests.Session,
    series_key: str,
    tmid: int,
    start_iso: str,
    stop_iso: str,
) -> SeriesData:
    """Download one series from WebTCAD and return as SeriesData."""
    url = build_url(tmid, start_iso, stop_iso)
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
        fig = Figure(figsize=(8, 6), tight_layout=True)
        super().__init__(fig)
        self.setParent(parent)
        self.figure = fig
        self.axes = []

        # Light UI-style background
        self.figure.set_facecolor("#f7f7fb")

    def plot_series(self, series_map: Dict[str, SeriesData]):
        self.figure.clear()

        valid_groups = [
            grp for grp in PLOT_GROUPS
            if (grp[1] in series_map) or (grp[2] in series_map)
        ]

        nrows = len(valid_groups)
        if DETECTOR_VOLT_KEY in series_map:
            nrows += 1

        if nrows == 0:
            self.draw()
            return

        axes = self.figure.subplots(nrows, 1, sharex=False)
        if nrows == 1:
            axes = [axes]
        self.axes = axes

        ax_idx = 0

        # Plot all min/max pairs
        for title, min_key, max_key in valid_groups:
            ax = axes[ax_idx]
            ax_idx += 1

            if min_key in series_map:
                s_min = series_map[min_key]
                df = s_min.df
                ax.plot(
                    df["time_dt"],
                    df[s_min.value_col],
                    linewidth=0.8,
                    label="Min",
                )

            if max_key in series_map:
                s_max = series_map[max_key]
                df = s_max.df
                ax.plot(
                    df["time_dt"],
                    df[s_max.value_col],
                    linewidth=0.8,
                    linestyle="--",
                    label="Max",
                )

            ax.set_ylabel(title)
            ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
            ax.legend(loc="upper right", fontsize=8)
            ax.tick_params(axis="both", labelsize=8)

        # Detector voltage in the last panel
        if DETECTOR_VOLT_KEY in series_map:
            ax = axes[ax_idx]
            s = series_map[DETECTOR_VOLT_KEY]
            df = s.df
            ax.plot(df["time_dt"], df[s.value_col], linewidth=0.9)
            ax.set_ylabel("Det V")
            ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
            ax.tick_params(axis="both", labelsize=8)

        # Shared x-label on the last axis
        axes[-1].set_xlabel("Time [UTC]", fontsize=9)

        for ax in axes:
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

        self.figure.tight_layout()
        self.draw()


# ----------------------------------------------------------------------
# Main window
# ----------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("IDEX WebTCAD Waveform Viewer")
        self.resize(1100, 700)

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
        self.start_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.start_edit.setCalendarPopup(True)

        self.stop_edit = QDateTimeEdit()
        self.stop_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.stop_edit.setCalendarPopup(True)

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

        # --- Plot canvas + toolbar -------------------------------------------
        self.canvas = TimeSeriesCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)

        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, stretch=1)

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
            series_map = self.download_all_series(username, password, start_iso, stop_iso)
        except Exception as e:
            self.run_btn.setEnabled(True)
            self.statusBar().showMessage("Error")
            QMessageBox.critical(self, "Download error", str(e))
            return

        # Save merged CSV
        try:
            out_name = (
                f"idex_waveforms_{start_dt.toString('yyyyMMdd_HHmmss')}"
                f"_to_{stop_dt.toString('yyyyMMdd_HHmmss')}.csv"
            )
            merge_to_single_csv(series_map, out_name)
        except Exception as e:
            # Non-fatal for plotting
            QMessageBox.warning(self, "CSV error", f"Could not write CSV:\n{e}")
            out_name = None

        # Plot
        self.canvas.plot_series(series_map)

        msg = "Finished"
        if out_name:
            msg += f" – saved merged CSV to {out_name}"
        self.statusBar().showMessage(msg)
        self.run_btn.setEnabled(True)

    # ------------------------------------------------------------------ #

    def download_all_series(
        self,
        username: str,
        password: str,
        start_iso: str,
        stop_iso: str,
    ) -> Dict[str, SeriesData]:
        """Download all configured series and return a map: key -> SeriesData."""
        session = requests.Session()
        session.auth = (username, password)

        series_map: Dict[str, SeriesData] = {}

        for key, tmid in SERIES_TMIDS.items():
            try:
                self.statusBar().showMessage(f"Fetching {key} (TMID={tmid})…")
                QApplication.processEvents()  # keep UI responsive
                s = fetch_series(session, key, tmid, start_iso, stop_iso)
                series_map[key] = s
            except Exception as e:
                raise RuntimeError(f"Error fetching {key} (TMID={tmid}): {e}") from e

        if not series_map:
            raise RuntimeError("No series downloaded.")

        return series_map


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
