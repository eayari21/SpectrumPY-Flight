"""Launch panel for SpectrumPY: Flight Edition.

Provides a welcome screen that highlights the mission imagery and lets the
user choose which analysis environment to enter.
"""

from __future__ import annotations

import os
os.environ.pop("QT_DEBUG_PLUGINS", None)

import subprocess
import sys
from pathlib import Path
from typing import List, Optional

from PyQt6.QtCore import Qt, QDateTime
from PyQt6.QtGui import QFont, QIcon, QPixmap, QResizeEvent
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QDateTimeEdit,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget,
)

import requests
from datetime import datetime, timezone

from .HDF_Explorer import HDFDataExplorer
from IDEX_quicklook import MainWindow as QuicklookWindow

APP_TITLE = "SpectrumPY: Flight Edition"
APP_AUTHOR = "Ethan Ayari"
REPO_ROOT = Path(__file__).resolve().parent
IMAGES_DIR = REPO_ROOT / "Images"
DOWNLOADS_DIR = REPO_ROOT / "Data" / "api_downloads"

IMAP_LOGO_CANDIDATES = (
    "IMAP_logo.png",
    "IMAP_logo.jpg",
    "IMAP_logo.jpeg",
    "imap_logo.png",
    "IMAP.png",
    "IMAP.jpeg",
)

INSTRUMENT_IMAGE_CANDIDATES = (
    "IDEX.jpeg",
    "instrument.jpg",
    "instrument.png",
    "Instrument.jpeg",
    "IDEX.png",
)

SUPPORTED_DATA_EXTENSIONS = (".h5", ".hdf5", ".cdf", ".trc")


def _find_image(candidates: tuple[str, ...]) -> Optional[Path]:
    for candidate in candidates:
        candidate_path = IMAGES_DIR / candidate
        if candidate_path.exists():
            return candidate_path
    return None


def _load_scaled_pixmap(path: Path) -> Optional[QPixmap]:
    pixmap = QPixmap(str(path))
    if pixmap.isNull():
        return None
    return pixmap


class ResponsiveImageLabel(QLabel):
    """QLabel that scales its pixmap while maintaining aspect ratio."""

    def __init__(self, *args, minimum_height: int | None = None, maximum_height: int | None = None, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._original_pixmap: Optional[QPixmap] = None
        self.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        if minimum_height is not None:
            self.setMinimumHeight(minimum_height)
        if maximum_height is not None:
            self.setMaximumHeight(maximum_height)

    def setPixmap(self, pixmap: QPixmap | None) -> None:  # type: ignore[override]
        if pixmap is None or pixmap.isNull():
            self._original_pixmap = None
            super().setPixmap(QPixmap())
            return

        self._original_pixmap = pixmap
        self._update_scaled_pixmap()

    def resizeEvent(self, event: QResizeEvent) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        if self._original_pixmap is not None:
            self._update_scaled_pixmap()

    def _update_scaled_pixmap(self) -> None:
        if not self._original_pixmap:
            return

        available_width = max(1, self.width())
        available_height = max(1, self.height())
        scaled = self._original_pixmap.scaled(
            available_width,
            available_height,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        super().setPixmap(scaled)


def _qt_to_datetime_utc(value: QDateTime) -> datetime:
    return value.toUTC().toPyDateTime().replace(tzinfo=timezone.utc)


def _format_day_of_year(dt: datetime) -> str:
    return dt.strftime("%Y-%jT%H:%M:%S")


def _format_iso_seconds(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


class FetchDataDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Fetch data via APIs")
        self.setModal(True)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.source_combo = QComboBox()
        self.source_combo.addItem("WebPODA (LASP)", "webpoda")
        self.source_combo.addItem("SCD Download Portal", "scd")
        form.addRow("Data source", self.source_combo)

        self.username_edit = QLineEdit()
        self.username_edit.setPlaceholderText("Enter your username")
        form.addRow("Username", self.username_edit)

        self.password_edit = QLineEdit()
        self.password_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self.password_edit.setPlaceholderText("Enter your password")
        form.addRow("Password", self.password_edit)

        now_utc = QDateTime.currentDateTimeUtc()
        self.start_edit = QDateTimeEdit(now_utc)
        self.start_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss 'UTC'")
        self.start_edit.setCalendarPopup(True)
        self.start_edit.setTimeSpec(Qt.TimeSpec.UTC)
        form.addRow("Start time", self.start_edit)

        self.end_edit = QDateTimeEdit(now_utc)
        self.end_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss 'UTC'")
        self.end_edit.setCalendarPopup(True)
        self.end_edit.setTimeSpec(Qt.TimeSpec.UTC)
        form.addRow("Stop time", self.end_edit)

        layout.addLayout(form)

        notice = QLabel(
            "Credentials are used only to perform the download and are not stored."
        )
        notice.setWordWrap(True)
        notice.setStyleSheet("color: #475569; font-size: 12px;")
        layout.addWidget(notice)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def values(self) -> dict[str, object]:
        return {
            "source": self.source_combo.currentData(),
            "username": self.username_edit.text(),
            "password": self.password_edit.text(),
            "start": self.start_edit.dateTime(),
            "stop": self.end_edit.dateTime(),
        }


class LaunchWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(APP_TITLE)
        self.setMinimumSize(800, 680)

        self._child_windows: List[QWidget] = []
        self._selected_path: Optional[Path] = None

        logo_path = _find_image(IMAP_LOGO_CANDIDATES)
        if logo_path:
            app = QApplication.instance()
            if app:
                app.setWindowIcon(QIcon(str(logo_path)))
            self.setWindowIcon(QIcon(str(logo_path)))

        central = QWidget(self)
        central.setObjectName("centralWidget")
        self.setCentralWidget(central)

        outer_layout = QVBoxLayout(central)
        outer_layout.setSpacing(0)
        outer_layout.setContentsMargins(48, 48, 48, 48)

        card = QFrame()
        card.setObjectName("launchCard")
        card.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        card.setMaximumWidth(840)

        layout = QVBoxLayout(card)
        layout.setSpacing(24)
        layout.setContentsMargins(48, 48, 48, 48)

        if logo_path:
            logo_pixmap = _load_scaled_pixmap(logo_path)
            if logo_pixmap:
                logo_label = ResponsiveImageLabel(maximum_height=120)
                logo_label.setPixmap(logo_pixmap)
                layout.addWidget(logo_label)

        title_label = QLabel(APP_TITLE)
        title_font = QFont()
        title_font.setPointSize(26)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        title_label.setObjectName("launchTitle")
        layout.addWidget(title_label)

        author_label = QLabel(f"by {APP_AUTHOR}")
        author_font = QFont()
        author_font.setPointSize(14)
        author_label.setFont(author_font)
        author_label.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        author_label.setStyleSheet("color: #475569;")
        layout.addWidget(author_label)

        instrument_path = _find_image(INSTRUMENT_IMAGE_CANDIDATES)
        if instrument_path:
            instrument_pixmap = _load_scaled_pixmap(instrument_path)
            if instrument_pixmap:
                instrument_container = QFrame()
                instrument_container.setObjectName("instrumentFrame")
                instrument_container_layout = QVBoxLayout(instrument_container)
                instrument_container_layout.setContentsMargins(16, 16, 16, 16)
                instrument_container_layout.setSpacing(0)

                instrument_label = ResponsiveImageLabel(minimum_height=260)
                instrument_label.setPixmap(instrument_pixmap)
                instrument_container_layout.addWidget(instrument_label)

                layout.addWidget(instrument_container)

        description = QLabel(
            "Select a science data product to begin. Once a file is selected you can choose "
            "whether to explore it with the HDF Plotter or jump directly into the IDEX Quicklook interface."
        )
        description.setWordWrap(True)
        description.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        description.setStyleSheet(
            "font-size: 15px; line-height: 1.5; color: #334155; padding: 0 12px;"
        )
        layout.addWidget(description)

        file_row = QHBoxLayout()
        file_row.setSpacing(16)
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Choose a CDF/HDF5/trace or packet file to work with")
        self.file_edit.setReadOnly(True)
        self.file_edit.setMinimumHeight(42)
        self.file_edit.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        file_row.addWidget(self.file_edit)

        browse_button = QPushButton("Browse…")
        browse_button.setMinimumHeight(42)
        browse_button.clicked.connect(self.select_file)
        browse_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        file_row.addWidget(browse_button)

        layout.addLayout(file_row)

        fetch_button = QPushButton("Fetch data via APIs")
        fetch_button.setMinimumHeight(42)
        fetch_button.clicked.connect(self.fetch_data_via_api)
        layout.addWidget(fetch_button)

        button_row = QHBoxLayout()
        button_row.setSpacing(16)
        button_row.addStretch()

        self.hdf_button = QPushButton("Open in HDF Plotter")
        self.hdf_button.setEnabled(False)
        self.hdf_button.setMinimumSize(200, 48)
        self.hdf_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.hdf_button.clicked.connect(self.launch_hdf_explorer)
        button_row.addWidget(self.hdf_button)

        self.quicklook_button = QPushButton("Open in IDEX Quicklook")
        self.quicklook_button.setEnabled(False)
        self.quicklook_button.setMinimumSize(200, 48)
        self.quicklook_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.quicklook_button.clicked.connect(self.launch_quicklook)
        button_row.addWidget(self.quicklook_button)

        button_row.addStretch()
        layout.addLayout(button_row)

        layout.addItem(QSpacerItem(20, 20, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        footer = QLabel("IMAP • Integrated data exploration for SpectrumPY")
        footer.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        footer.setStyleSheet("color: #4a5568; font-size: 12px;")
        layout.addWidget(footer)

        outer_layout.addStretch(1)
        outer_layout.addWidget(card, alignment=Qt.AlignmentFlag.AlignHCenter)
        outer_layout.addStretch(1)

        self.setStyleSheet(
            """
            QWidget#centralWidget {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                            stop:0 #0f172a, stop:1 #1e293b);
            }
            QFrame#launchCard {
                background: rgba(255, 255, 255, 0.96);
                border-radius: 24px;
            }
            QFrame#instrumentFrame {
                background: #0f172a;
                border-radius: 18px;
            }
            QLabel {
                color: #1a202c;
            }
            QLabel#launchTitle {
                color: #0b1120;
            }
            QPushButton {
                background-color: #2563eb;
                color: white;
                border: none;
                border-radius: 12px;
                padding: 12px 24px;
                font-size: 15px;
                font-weight: 600;
            }
            QPushButton:disabled {
                background-color: #94a3b8;
                color: rgba(255, 255, 255, 0.85);
            }
            QPushButton:hover:!disabled {
                background-color: #1d4ed8;
            }
            QPushButton:pressed:!disabled {
                background-color: #1e40af;
            }
            QLineEdit {
                background: rgba(15, 23, 42, 0.04);
                border: 1px solid rgba(15, 23, 42, 0.12);
                border-radius: 12px;
                padding: 0 16px;
                font-size: 14px;
            }
            """
        )

    # ------------------------------------------------------------------
    def select_file(self) -> None:
        start_dir = self._selected_path.parent if self._selected_path else (REPO_ROOT / "HDF5")
        if not start_dir.exists():
            start_dir = REPO_ROOT

        file_dialog = QFileDialog(self, "Select data file", str(start_dir))
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilters([
            "Science Data (*.h5 *.hdf5 *.cdf *.trc *.bin)",
            "HDF5 (*.h5 *.hdf5)",
            "CDF (*.cdf)",
            "Trace (*.trc)",
            "Binary Packets (*.bin)",
            "All files (*)",
        ])

        if file_dialog.exec() == QFileDialog.DialogCode.Accepted:
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                self._update_selected_path(Path(selected_files[0]))

    def fetch_data_via_api(self) -> None:
        dialog = FetchDataDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        values = dialog.values()
        source = str(values.get("source") or "").strip()
        username = str(values.get("username") or "").strip()
        password = str(values.get("password") or "")
        start_qt = values.get("start")
        stop_qt = values.get("stop")

        if not username or not password:
            QMessageBox.warning(self, "Missing credentials", "Please enter both a username and password.")
            return

        if not isinstance(start_qt, QDateTime) or not isinstance(stop_qt, QDateTime):
            QMessageBox.critical(self, "Invalid input", "Unable to read the requested time range.")
            return

        start_dt = _qt_to_datetime_utc(start_qt)
        stop_dt = _qt_to_datetime_utc(stop_qt)
        if stop_dt <= start_dt:
            QMessageBox.warning(self, "Invalid time range", "The stop time must be after the start time.")
            return

        try:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            if source == "webpoda":
                downloaded_path = self._download_webpoda(username, password, start_dt, stop_dt)
            elif source == "scd":
                downloaded_path = self._download_scd(username, password, start_dt, stop_dt)
            else:
                QMessageBox.critical(self, "Unknown source", "Please choose a data source to continue.")
                return
        except requests.HTTPError as exc:
            QMessageBox.critical(
                self,
                "Download failed",
                f"The remote server returned an error: {exc}",
            )
            return
        except requests.RequestException as exc:
            QMessageBox.critical(
                self,
                "Download failed",
                f"Unable to fetch data from the selected API: {exc}",
            )
            return
        except OSError as exc:
            QMessageBox.critical(
                self,
                "Save failed",
                f"Unable to save the downloaded file: {exc}",
            )
            return
        finally:
            QApplication.restoreOverrideCursor()

        if not downloaded_path:
            QMessageBox.warning(
                self,
                "No data returned",
                "The selected API did not return any data for the requested interval.",
            )
            return

        self._update_selected_path(downloaded_path)
        QMessageBox.information(
            self,
            "Download complete",
            f"Data downloaded to:\n{downloaded_path}",
        )

    def _update_selected_path(self, path: Path) -> None:
        prepared_path = self._prepare_selected_path(path)
        if prepared_path is None:
            return

        self._selected_path = prepared_path
        self.file_edit.setText(str(prepared_path))
        enables_quicklook = prepared_path.suffix.lower() in SUPPORTED_DATA_EXTENSIONS
        enables_hdf = prepared_path.suffix.lower() in {".h5", ".hdf5"}
        self.quicklook_button.setEnabled(enables_quicklook)
        self.hdf_button.setEnabled(enables_hdf)

    def _prepare_selected_path(self, path: Path) -> Optional[Path]:
        suffix = path.suffix.lower()
        if not suffix or suffix == ".bin":
            converted = self._convert_packet_file(path)
            return converted
        return path

    def _download_webpoda(
        self,
        username: str,
        password: str,
        start: datetime,
        stop: datetime,
    ) -> Optional[Path]:
        DOWNLOADS_DIR.mkdir(parents=True, exist_ok=True)
        query_params = [
            ("time>=", _format_day_of_year(start)),
            ("time<", _format_day_of_year(stop)),
            ("project(time,packet)", ""),
            ('formatTime("yyyy-DDD\'T\'HH:mm:ss")', ""),
        ]
        url = "https://lasp.colorado.edu/ops/imap/poda/dap2/packets/SID1/IDX_SCI.asc"

        total_written = 0
        filename = f"IDX_SCI_{start.strftime('%Y%m%dT%H%M%S')}_{stop.strftime('%Y%m%dT%H%M%S')}.asc"
        target = DOWNLOADS_DIR / filename

        with requests.get(url, params=query_params, auth=(username, password), timeout=120, stream=True) as response:
            response.raise_for_status()
            with open(target, "wb") as handle:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        handle.write(chunk)
                        total_written += len(chunk)

        if total_written == 0:
            if target.exists():
                target.unlink()
            return None

        return target

    def _download_scd(
        self,
        username: str,
        password: str,
        start: datetime,
        stop: datetime,
    ) -> Optional[Path]:
        DOWNLOADS_DIR.mkdir(parents=True, exist_ok=True)
        params = {
            "start": _format_iso_seconds(start),
            "stop": _format_iso_seconds(stop),
        }
        url = "https://imap-mission.com/data/download"

        def _extract_filename(header_value: str | None) -> Optional[str]:
            if not header_value:
                return None
            parts = header_value.split("filename=")
            if len(parts) < 2:
                return None
            candidate = parts[1].strip()
            candidate = candidate.strip('"')
            candidate = candidate.strip("'")
            candidate = candidate.strip()
            return candidate or None

        fallback = f"imap_download_{start.strftime('%Y%m%dT%H%M%S')}_{stop.strftime('%Y%m%dT%H%M%S')}.bin"

        total_written = 0
        target: Optional[Path] = None

        with requests.get(url, params=params, auth=(username, password), timeout=120, stream=True) as response:
            response.raise_for_status()
            filename = _extract_filename(response.headers.get("Content-Disposition")) or fallback
            target = DOWNLOADS_DIR / filename
            with open(target, "wb") as handle:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        handle.write(chunk)
                        total_written += len(chunk)

        if not target:
            return None

        if total_written == 0 and target.exists():
            target.unlink()
            return None

        return target

    def _convert_packet_file(self, path: Path) -> Optional[Path]:
        output_dir = REPO_ROOT / "HDF5"
        output_dir.mkdir(parents=True, exist_ok=True)
        produced_path = output_dir / f"{path.name}.h5"

        command = [sys.executable, str(REPO_ROOT / "idex_packet.py"), "--file", str(path)]

        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            result = subprocess.run(
                command,
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            details = exc.stderr or exc.stdout or ""
            if details:
                details = details.strip()
                if len(details) > 1200:
                    details = details[:1200] + "…"
            if details:
                message = (
                    "idex_packet.py failed to convert the selected packet file.\n\n"
                    f"Details:\n{details}"
                )
            else:
                message = "idex_packet.py failed to convert the selected packet file."
            QMessageBox.critical(self, "Conversion failed", message)
            return None
        except OSError as exc:
            QMessageBox.critical(
                self,
                "Conversion failed",
                f"Unable to execute idex_packet.py: {exc}",
            )
            return None
        finally:
            QApplication.restoreOverrideCursor()

        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        if not produced_path.exists():
            QMessageBox.critical(
                self,
                "Conversion failed",
                "idex_packet.py finished without producing an HDF5 file.",
            )
            return None

        QMessageBox.information(
            self,
            "Packet converted",
            (
                "The selected packet was converted to HDF5.\n"
                f"Output file: {produced_path}"
            ),
        )
        return produced_path

    # ------------------------------------------------------------------
    def launch_hdf_explorer(self) -> None:
        if not self._selected_path:
            QMessageBox.information(self, "No data selected", "Please choose an HDF5 file first.")
            return

        if self._selected_path.suffix.lower() not in {".h5", ".hdf5"}:
            QMessageBox.warning(
                self,
                "Unsupported file",
                "The HDF Explorer works with HDF5 products. Select a *.h5 file to continue.",
            )
            return

        plotter = HDFDataExplorer(self._selected_path)
        plotter.show()
        self._register_child_window(plotter)

    def launch_quicklook(self) -> None:
        if not self._selected_path:
            QMessageBox.information(self, "No data selected", "Please choose a data file first.")
            return

        quicklook = QuicklookWindow(str(self._selected_path))
        quicklook.show()
        self._register_child_window(quicklook)

    def _register_child_window(self, window: QWidget) -> None:
        self._child_windows.append(window)
        window.destroyed.connect(lambda: self._child_windows.remove(window) if window in self._child_windows else None)


def main() -> int:
    app = QApplication(sys.argv)

    logo_path = _find_image(IMAP_LOGO_CANDIDATES)
    if logo_path:
        app.setWindowIcon(QIcon(str(logo_path)))

    window = LaunchWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
