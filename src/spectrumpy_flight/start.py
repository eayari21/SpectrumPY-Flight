"""Convenience entry point for launching the SpectrumPY UI.

This entry point extends the default launcher with automatic processing of
directories that contain oscilloscope ``.trc`` files. When such a directory is
selected, the data is converted to an HDF5 product via :mod:`ImpactBook` and the
result is immediately opened in the quicklook interface.
"""

from __future__ import annotations

import os
import re
import sys
import traceback
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import Sequence

# ``python start.py`` executes the module as a script which clears ``__package__``.
# When that happens we manually add the parent directory of this file to
# ``sys.path`` so absolute imports resolve the installed package just like they do
# when launched via ``python -m`` or the console script entry point.
if __package__ in (None, ""):
    _this_file = Path(__file__).resolve()
    _package_root = _this_file.parent
    _src_root = _package_root.parent
    if str(_src_root) not in sys.path:
        sys.path.insert(0, str(_src_root))
    __package__ = "spectrumpy_flight"

os.environ.pop("QT_DEBUG_PLUGINS", None)

from spectrumpy_flight import package_path, spectrum_launcher
from spectrumpy_flight.readTrc import Trc
from spectrumpy_flight.spectrum_launcher import LaunchWindow

# --------- Qt binding-agnostic imports (prefer PySide6, fallback PyQt6) ---------
_QT = None
try:
    from PySide6.QtCore import QTimer, Qt
    from PySide6.QtGui import QCloseEvent, QIcon
    from PySide6.QtWidgets import (
        QApplication,
        QFileDialog,
        QDialog,
        QDialogButtonBox,
        QFormLayout,
        QInputDialog,
        QLineEdit,
        QMessageBox,
        QProgressDialog,
        QVBoxLayout,
    )
    _QT = "PySide6"
except Exception:
    from PyQt6.QtCore import QTimer, Qt
    from PyQt6.QtGui import QCloseEvent, QIcon
    from PyQt6.QtWidgets import (
        QApplication,
        QFileDialog,
        QDialog,
        QDialogButtonBox,
        QFormLayout,
        QInputDialog,
        QLineEdit,
        QMessageBox,
        QProgressDialog,
        QVBoxLayout,
    )
    _QT = "PyQt6"

REPO_ROOT = package_path()
TRC_SUFFIX = ".trc"


def _sanitize_experiment_name(name: str) -> str:
    sanitized = re.sub(r"[\\/:*?\"<>|]+", "_", name.strip())
    return sanitized or "impact_dataset"


class ChannelNamingDialog(QDialog):
    """Dialog that prompts the user to label each detected channel."""

    def __init__(self, suggestions: Sequence[str], parent: QDialog | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Channel names")
        self._names: list[str] | None = None

        layout = QVBoxLayout(self)
        form = QFormLayout()
        self._edits: list[QLineEdit] = []

        if not suggestions:
            suggestions = [f"Channel {index}" for index in range(1, 5)]

        for index, suggestion in enumerate(suggestions, start=1):
            edit = QLineEdit(self)
            edit.setText(suggestion)
            edit.setPlaceholderText(f"Channel {index}")
            form.addRow(f"Channel {index}", edit)
            self._edits.append(edit)

        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok,
            parent=self,
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _on_accept(self) -> None:
        names = [edit.text().strip() for edit in self._edits]
        if any(not name for name in names):
            QMessageBox.warning(self, "Missing channel name", "Please provide a name for each channel.")
            return
        self._names = names
        self.accept()

    def channel_names(self) -> list[str]:
        return list(self._names or [])


class AutoProcessingLaunchWindow(LaunchWindow):
    """Launch window that can automatically process directories of ``.trc`` files."""

    def __init__(self) -> None:
        super().__init__()
        self._executor = ThreadPoolExecutor(max_workers=1)
        self._active_future: Future[Path] | None = None
        self._progress_dialog: QProgressDialog | None = None

    # ------------------------------------------------------------------
    def select_file(self) -> None:  # type: ignore[override]
        start_dir = self._selected_path.parent if self._selected_path else (REPO_ROOT / "HDF5")
        if not start_dir.exists():
            start_dir = REPO_ROOT

        file_dialog = QFileDialog(self, "Select data file or directory", str(start_dir))
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setOption(QFileDialog.Option.DontUseNativeDialog, True)
        file_dialog.setOption(QFileDialog.Option.ShowDirsOnly, False)
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
            if not selected_files:
                return
            selected_path = Path(selected_files[0])
            if selected_path.is_dir():
                self._process_trc_directory(selected_path)
            else:
                self._update_selected_path(selected_path)

    # ------------------------------------------------------------------
    def _update_selected_path(self, path: Path) -> None:  # type: ignore[override]
        if path.is_dir():
            self._process_trc_directory(path)
            return

        if path.suffix.lower() == TRC_SUFFIX:
            self._process_trc_directory(path.parent)
            return

        super()._update_selected_path(path)

    # ------------------------------------------------------------------
    def _process_trc_directory(self, directory: Path) -> None:
        directory = directory.resolve()

        trc_files = sorted(directory.glob(f"*{TRC_SUFFIX}"))
        if not trc_files:
            QMessageBox.warning(
                self,
                "No trace files found",
                "The selected directory does not contain any .trc files.",
            )
            return

        if self._active_future is not None and not self._active_future.done():
            QMessageBox.information(
                self,
                "Processing in progress",
                "Please wait for the current ImpactBook run to finish before starting another.",
            )
            return

        suggestions = self._detect_channel_suggestions(trc_files)

        naming_dialog = ChannelNamingDialog(suggestions, parent=self)
        if naming_dialog.exec() != QDialog.Accepted:
            return

        channel_names = naming_dialog.channel_names()

        experiment_name, ok = QInputDialog.getText(
            self,
            "Experiment name",
            "Name for the generated HDF5 product:",
            text=_sanitize_experiment_name(directory.name),
        )
        if not ok:
            return

        experiment_name = _sanitize_experiment_name(experiment_name)
        if not experiment_name:
            QMessageBox.warning(self, "Invalid name", "Please provide a valid experiment name.")
            return

        progress = QProgressDialog("Processing trace directory…", None, 0, 0, self)
        progress.setWindowTitle("ImpactBook")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setCancelButton(None)
        progress.setRange(0, 0)
        progress.show()

        self._progress_dialog = progress

        def run_impactbook() -> Path:
            from spectrumpy_flight.ImpactBook import ImpactBook

            ImpactBook(channel_names, trcdir=str(directory), ExperimentName=experiment_name)
            return (REPO_ROOT / "HDF5" / f"{experiment_name}.h5").resolve()

        future = self._executor.submit(run_impactbook)
        self._active_future = future

        def _on_finished(done_future: Future[Path]) -> None:
            def finish() -> None:
                self._handle_processing_complete(done_future, directory, experiment_name)

            QTimer.singleShot(0, finish)

        future.add_done_callback(_on_finished)

    # ------------------------------------------------------------------
    def _handle_processing_complete(self, future: Future[Path], directory: Path, experiment_name: str) -> None:
        if self._progress_dialog is not None:
            self._progress_dialog.close()
            self._progress_dialog.deleteLater()
            self._progress_dialog = None

        self._active_future = None

        try:
            h5_path = future.result()
        except Exception as exc:  # pragma: no cover - GUI feedback path
            details = "".join(traceback.format_exception(exc))
            QMessageBox.critical(
                self,
                "ImpactBook failed",
                "ImpactBook was unable to process the selected directory.\n\n" + details,
            )
            return

        if not h5_path.exists():
            QMessageBox.warning(
                self,
                "Missing output",
                (
                    "ImpactBook completed without producing the expected HDF5 file.\n"
                    f"Expected location: {h5_path}"
                ),
            )
            return

        QMessageBox.information(
            self,
            "ImpactBook complete",
            (
                f"Directory '{directory}' has been processed into an HDF5 product.\n"
                f"Location: {h5_path}"
            ),
        )

        super()._update_selected_path(h5_path)
        self.launch_quicklook()

    # ------------------------------------------------------------------
    def _detect_channel_suggestions(self, trc_files: Sequence[Path]) -> list[str]:
        labels: list[str] = []
        seen: set[str] = set()
        trc_reader = Trc()

        for trc_path in trc_files:
            try:
                _, _, meta = trc_reader.open(str(trc_path))
            except OSError:
                continue

            label = (meta.get("TRACE_LABEL") or "").strip()
            if not label:
                label = trc_path.stem

            if label in seen:
                break

            seen.add(label)
            labels.append(label)

        if labels:
            return labels

        prefix_seen: set[str] = set()
        prefixes: list[str] = []
        for path in trc_files:
            stem = path.stem
            prefix = re.sub(r"\d+$", "", stem) or stem
            if prefix in prefix_seen:
                break
            prefix_seen.add(prefix)
            prefixes.append(prefix)

        return prefixes

    # ------------------------------------------------------------------
    def closeEvent(self, event: QCloseEvent) -> None:  # type: ignore[override]
        if self._executor is not None:
            self._executor.shutdown(wait=False, cancel_futures=True)
            self._executor = None  # type: ignore[assignment]
        super().closeEvent(event)


def main() -> int:
    app = QApplication(sys.argv)

    logo_path = spectrum_launcher._find_image(spectrum_launcher.IMAP_LOGO_CANDIDATES)
    if logo_path:
        app.setWindowIcon(QIcon(str(logo_path)))

    window = AutoProcessingLaunchWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
