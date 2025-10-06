"""Qt-based tree viewer for Common Data Format (CDF) files.

This module mirrors the ergonomics of the bundled :mod:`HDF_View` helper,
providing a lightweight browser for inspecting variable hierarchies, previewing
values, and reading attributes.  The implementation relies on :mod:`cdflib`
and therefore surfaces a clear error when the dependency is unavailable.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterable, Optional, Sequence

import numpy as np

try:  # pragma: no cover - optional dependency
    import cdflib  # type: ignore
except Exception as exc:  # pragma: no cover - import guard
    raise RuntimeError(
        "CDF viewer requires the optional 'cdflib' dependency."
    ) from exc

# ---------------------------------------------------------------------------
# Qt binding import (prefer PySide6, fall back to PyQt6)
# ---------------------------------------------------------------------------
_QT_API = None
try:  # pragma: no cover - import guard
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import (
        QApplication,
        QFileDialog,
        QHeaderView,
        QMainWindow,
        QMessageBox,
        QPlainTextEdit,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QTabWidget,
        QTreeWidget,
        QTreeWidgetItem,
        QVBoxLayout,
        QWidget,
    )
    _QT_API = "PySide6"
except Exception:  # pragma: no cover - import guard
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QApplication,
        QFileDialog,
        QHeaderView,
        QMainWindow,
        QMessageBox,
        QPlainTextEdit,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QTabWidget,
        QTreeWidget,
        QTreeWidgetItem,
        QVBoxLayout,
        QWidget,
    )
    _QT_API = "PyQt6"


def _format_scalar(value: object) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        try:
            return value.decode("utf-8")
        except Exception:
            return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return _format_scalar(value.item())
        return np.array2string(value, threshold=10)
    return str(value)


def _shape_to_text(shape: Sequence[int]) -> str:
    if not shape:
        return "Scalar"
    return "×".join(str(dim) for dim in shape)


class CDFViewWindow(QMainWindow):
    """Interactive viewer for browsing CDF variables and attributes."""

    MAX_PREVIEW_ROWS = 200
    MAX_PREVIEW_COLS = 60

    def __init__(self, filename: str, parent: Optional[QWidget] = None):
        super().__init__(parent)
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File not found: {filename}")

        try:
            self._cdf = cdflib.CDF(filename)
        except Exception as exc:
            raise RuntimeError(f"Failed to open CDF file: {filename}\n{exc}") from exc

        self._filename = filename
        self.setWindowTitle(f"CDF View — {os.path.basename(filename)}")
        self.resize(1100, 720)

        self._tree: QTreeWidget
        self._summary: QPlainTextEdit
        self._data_table: QTableWidget
        self._attr_table: QTableWidget
        self._build_ui()
        self._populate_tree()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        central = QWidget(self)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        layout.addWidget(splitter)

        self._tree = QTreeWidget(splitter)
        self._tree.setHeaderLabels(["Name", "Type", "Shape"])
        self._tree.setAlternatingRowColors(True)
        self._tree.header().setStretchLastSection(False)
        self._tree.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self._tree.header().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self._tree.header().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self._tree.currentItemChanged.connect(self._on_item_selected)

        tabs = QTabWidget(splitter)

        self._summary = QPlainTextEdit(tabs)
        self._summary.setReadOnly(True)
        tabs.addTab(self._summary, "Summary")

        self._data_table = QTableWidget(tabs)
        self._data_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._data_table.verticalHeader().setVisible(True)
        self._data_table.horizontalHeader().setVisible(True)
        self._data_table.horizontalHeader().setStretchLastSection(False)
        tabs.addTab(self._data_table, "Data Preview")

        self._attr_table = QTableWidget(tabs)
        self._attr_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._attr_table.setColumnCount(2)
        self._attr_table.setHorizontalHeaderLabels(["Attribute", "Value"])
        self._attr_table.horizontalHeader().setStretchLastSection(True)
        tabs.addTab(self._attr_table, "Attributes")

        self.setCentralWidget(central)

    # ------------------------------------------------------------------
    # Population helpers
    # ------------------------------------------------------------------
    def _populate_tree(self) -> None:
        self._tree.clear()
        root = QTreeWidgetItem([os.path.basename(self._filename) or self._filename, "File", ""])
        root.setData(0, Qt.ItemDataRole.UserRole, ("file", None))
        self._tree.addTopLevelItem(root)

        info = self._cdf.cdf_info()
        for category, label in (("zVariables", "Z Variables"), ("rVariables", "R Variables")):
            items: Iterable[str] = info.get(category, []) or []
            if not items:
                continue
            group_item = QTreeWidgetItem([label, "Group", ""])
            group_item.setData(0, Qt.ItemDataRole.UserRole, ("group", None))
            root.addChild(group_item)
            for varname in items:
                summary = self._variable_summary(varname)
                child = QTreeWidgetItem([varname, summary["type"], _shape_to_text(summary["shape"])])
                child.setData(0, Qt.ItemDataRole.UserRole, ("variable", varname))
                group_item.addChild(child)

        self._tree.expandItem(root)
        self._tree.setCurrentItem(root)

    def _variable_summary(self, name: str) -> Dict[str, object]:
        summary: Dict[str, object] = {"name": name, "shape": (), "type": "Variable"}
        try:
            info = self._cdf.varinq(name)
        except Exception as exc:
            summary["error"] = str(exc)
            return summary

        data_type = info.get("Data_Type_Description") or info.get("Data_Type")
        summary["type"] = str(data_type)
        summary["shape"] = tuple(info.get("Dim_Sizes", ()))
        summary["variance"] = info.get("Rec_Vary")
        summary["num_elem"] = info.get("Num_Elements")
        summary["last_rec"] = info.get("Last_Rec")
        return summary

    # ------------------------------------------------------------------
    # Selection handlers
    # ------------------------------------------------------------------
    def _on_item_selected(self, current: Optional[QTreeWidgetItem], previous: Optional[QTreeWidgetItem]) -> None:
        del previous  # unused
        if current is None:
            self._summary.clear()
            self._data_table.clear()
            self._attr_table.clear()
            return

        node_type, payload = current.data(0, Qt.ItemDataRole.UserRole) or (None, None)
        if node_type != "variable" or payload is None:
            self._summary.setPlainText("Select a variable to view its metadata.")
            self._data_table.clear()
            self._attr_table.clear()
            self._attr_table.setRowCount(0)
            self._data_table.setRowCount(0)
            return

        varname = str(payload)
        summary = self._variable_summary(varname)
        lines = [f"Variable: {varname}"]
        for key in ("type", "shape", "variance", "num_elem", "last_rec"):
            if key in summary and summary[key] is not None:
                value = summary[key]
                if key == "shape":
                    value = _shape_to_text(value if isinstance(value, Sequence) else ())
                lines.append(f"{key.replace('_', ' ').title()}: {value}")
        if "error" in summary:
            lines.append(f"Error: {summary['error']}")
        self._summary.setPlainText("\n".join(lines))

        self._populate_data_preview(varname)
        self._populate_attributes(varname)

    def _populate_data_preview(self, varname: str) -> None:
        self._data_table.clear()
        self._data_table.setRowCount(0)
        self._data_table.setColumnCount(0)

        try:
            data = self._cdf.varget(varname)
        except Exception as exc:
            self._data_table.setRowCount(1)
            self._data_table.setColumnCount(1)
            self._data_table.setItem(0, 0, QTableWidgetItem(f"Error reading data: {exc}"))
            return

        array = np.asarray(data)
        if array.ndim == 0:
            self._data_table.setRowCount(1)
            self._data_table.setColumnCount(1)
            self._data_table.setItem(0, 0, QTableWidgetItem(_format_scalar(array)))
            return

        if array.ndim == 1:
            array = array.reshape(array.shape[0], 1)

        rows = min(array.shape[0], self.MAX_PREVIEW_ROWS)
        cols = min(array.shape[1], self.MAX_PREVIEW_COLS)

        self._data_table.setRowCount(rows)
        self._data_table.setColumnCount(cols)
        for r in range(rows):
            for c in range(cols):
                try:
                    value = array[r, c]
                except Exception:
                    value = None
                self._data_table.setItem(r, c, QTableWidgetItem(_format_scalar(value)))

    def _populate_attributes(self, varname: str) -> None:
        self._attr_table.clearContents()
        self._attr_table.setRowCount(0)

        try:
            attrs = self._cdf.varattsget(varname, expand=True) or {}
        except Exception:
            attrs = {}

        rows = len(attrs)
        self._attr_table.setRowCount(rows)
        for idx, (key, value) in enumerate(sorted(attrs.items())):
            key_item = QTableWidgetItem(str(key))
            val_item = QTableWidgetItem(_format_scalar(value))
            key_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            val_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            self._attr_table.setItem(idx, 0, key_item)
            self._attr_table.setItem(idx, 1, val_item)

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------
    def closeEvent(self, event):  # pragma: no cover - GUI teardown
        try:
            self._cdf.close()
        except Exception:
            pass
        event.accept()


def launch_cdf_viewer(filename: str, parent: Optional[QWidget] = None) -> QWidget:
    """Launch the CDF viewer window and return it for bookkeeping."""

    window = CDFViewWindow(filename, parent=parent)
    window.show()
    return window


def _main() -> int:  # pragma: no cover - convenience CLI
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    filename, _ = QFileDialog.getOpenFileName(
        None,
        "Open CDF",
        os.getcwd(),
        "CDF files (*.cdf);;All files (*)",
    )
    if not filename:
        return 0

    try:
        window = launch_cdf_viewer(filename)
    except Exception as exc:
        QMessageBox.critical(None, "CDF Viewer", str(exc))
        return 1

    return app.exec()


if __name__ == "__main__":  # pragma: no cover - manual invocation
    sys.exit(_main())
