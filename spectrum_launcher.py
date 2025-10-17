"""Launch panel for SpectrumPY: Flight Addition.

Provides a welcome screen that highlights the mission imagery and lets the
user choose which analysis environment to enter.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QIcon, QPixmap, QResizeEvent
from PyQt6.QtWidgets import (
    QApplication,
    QFileDialog,
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

from HDF_Explorer import HDFDataExplorer
from IDEX_quicklook import MainWindow as QuicklookWindow

APP_TITLE = "SpectrumPY: Flight Addition"
APP_AUTHOR = "Ethan Ayari"
REPO_ROOT = Path(__file__).resolve().parent
IMAGES_DIR = REPO_ROOT / "Images"

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
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setSpacing(18)
        layout.setContentsMargins(32, 32, 32, 32)

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
        layout.addWidget(title_label)

        author_label = QLabel(f"by {APP_AUTHOR}")
        author_font = QFont()
        author_font.setPointSize(14)
        author_label.setFont(author_font)
        author_label.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(author_label)

        instrument_path = _find_image(INSTRUMENT_IMAGE_CANDIDATES)
        if instrument_path:
            instrument_pixmap = _load_scaled_pixmap(instrument_path)
            if instrument_pixmap:
                instrument_label = ResponsiveImageLabel(minimum_height=260)
                instrument_label.setPixmap(instrument_pixmap)
                layout.addWidget(instrument_label)

        description = QLabel(
            "Select a science data product to begin. Once a file is selected you can choose "
            "whether to explore it with the HDF Plotter or jump directly into the IDEX Quicklook interface."
        )
        description.setWordWrap(True)
        description.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        description.setStyleSheet("font-size: 14px;")
        layout.addWidget(description)

        file_row = QHBoxLayout()
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Choose a CDF/HDF5/trace file to work with")
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

    # ------------------------------------------------------------------
    def select_file(self) -> None:
        start_dir = self._selected_path.parent if self._selected_path else (REPO_ROOT / "HDF5")
        if not start_dir.exists():
            start_dir = REPO_ROOT

        file_dialog = QFileDialog(self, "Select data file", str(start_dir))
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilters([
            "Science Data (*.h5 *.hdf5 *.cdf *.trc)",
            "HDF5 (*.h5 *.hdf5)",
            "CDF (*.cdf)",
            "Trace (*.trc)",
            "All files (*)",
        ])

        if file_dialog.exec() == QFileDialog.DialogCode.Accepted:
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                self._update_selected_path(Path(selected_files[0]))

    def _update_selected_path(self, path: Path) -> None:
        self._selected_path = path
        self.file_edit.setText(str(path))
        enables_quicklook = path.suffix.lower() in SUPPORTED_DATA_EXTENSIONS
        enables_hdf = path.suffix.lower() in {".h5", ".hdf5"}
        self.quicklook_button.setEnabled(enables_quicklook)
        self.hdf_button.setEnabled(enables_hdf)

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
