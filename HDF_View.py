#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lightweight HDF5 browser used by IDEX Quicklook.

This module provides a Qt-based tree viewer that emulates the core
capabilities of the HDFView application.  Users can explore the structure of
an HDF5 file, inspect attributes, and preview dataset values.  The viewer can
be launched standalone (``python HDF_View.py <file>``) or embedded within the
IDEX Quicklook GUI.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable, Optional, Sequence, Tuple

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# Qt binding import (prefer PySide6, fallback to PyQt6)
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
    """Return a human-readable representation of a scalar value."""

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


def _format_attribute(value: object) -> str:
    """Format attribute values for display in a table."""

    text = _format_scalar(value)
    if len(text) > 200:
        return text[:197] + "…"
    return text


def _shape_to_text(shape: Sequence[int]) -> str:
    return "×".join(str(dim) for dim in shape)


class HDFViewWindow(QMainWindow):
    """Simple HDF5 browser window."""

    MAX_PREVIEW_ROWS = 200
    MAX_PREVIEW_COLS = 50

    def __init__(self, filename: str, parent: Optional[QWidget] = None):
        super().__init__(parent)
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File not found: {filename}")

        try:
            self._h5 = h5py.File(filename, "r")
        except Exception as exc:
            raise RuntimeError(f"Failed to open HDF5 file: {filename}\n{exc}") from exc

        self._filename = filename

        self.setWindowTitle(f"HDF View — {os.path.basename(filename)}")
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

    def _populate_tree(self) -> None:
        self._tree.clear()
        root = QTreeWidgetItem([os.path.basename(self._filename) or self._filename, "File", ""])
        root.setData(0, Qt.ItemDataRole.UserRole, "/")
        root.setData(1, Qt.ItemDataRole.UserRole, "file")
        self._tree.addTopLevelItem(root)
        self._add_children(root, self._h5)
        self._tree.expandItem(root)
        self._tree.setCurrentItem(root)

    def _add_children(self, parent_item: QTreeWidgetItem, group: h5py.Group) -> None:
        keys = sorted(group.keys())
        for key in keys:
            obj = group[key]
            if isinstance(obj, h5py.Group):
                child = QTreeWidgetItem([key, "Group", ""])
                child.setData(0, Qt.ItemDataRole.UserRole, obj.name)
                child.setData(1, Qt.ItemDataRole.UserRole, "group")
                parent_item.addChild(child)
                self._add_children(child, obj)
            elif isinstance(obj, h5py.Dataset):
                shape_text = _shape_to_text(obj.shape) if obj.shape else "Scalar"
                child = QTreeWidgetItem([key, "Dataset", shape_text])
                child.setData(0, Qt.ItemDataRole.UserRole, obj.name)
                child.setData(1, Qt.ItemDataRole.UserRole, "dataset")
                parent_item.addChild(child)

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------
    def _on_item_selected(
        self,
        current: Optional[QTreeWidgetItem],
        previous: Optional[QTreeWidgetItem],
    ) -> None:
        del previous  # unused
        if current is None:
            return
        path = current.data(0, Qt.ItemDataRole.UserRole)
        node_type = current.data(1, Qt.ItemDataRole.UserRole)
        if node_type == "dataset":
            self._show_dataset(path)
        elif node_type in {"group", "file"}:
            self._show_group(path)
        else:
            self._summary.setPlainText("")
            self._clear_table(self._data_table)
            self._clear_table(self._attr_table)

    # ------------------------------------------------------------------
    # Display helpers
    # ------------------------------------------------------------------
    def _show_group(self, path: str) -> None:
        obj = self._h5[path] if path != "/" else self._h5
        lines = [f"Path: {path}", "Type: Group"]
        if isinstance(obj, h5py.File):
            lines[1] = "Type: File"
        lines.append(f"Children: {len(obj.keys())}")
        attr_items = sorted(obj.attrs.items(), key=lambda item: item[0])
        lines.append(f"Attributes: {len(attr_items)}")
        self._summary.setPlainText("\n".join(lines))
        self._clear_table(self._data_table)
        self._populate_attrs(attr_items)

    def _show_dataset(self, path: str) -> None:
        dset = self._h5[path]
        summary_lines = [
            f"Path: {path}",
            "Type: Dataset",
            f"Shape: {dset.shape if dset.shape else '()'}",
            f"Dtype: {dset.dtype}",
            f"Size: {dset.size}",
        ]
        if dset.chunks:
            summary_lines.append(f"Chunks: {dset.chunks}")
        if dset.compression:
            summary_lines.append(f"Compression: {dset.compression}")
        if dset.scaleoffset is not None:
            summary_lines.append(f"Scale/Offset: {dset.scaleoffset}")

        data = None
        preview_note = ""
        try:
            data = dset[()]
        except Exception as exc:
            preview_note = f"\n\nPreview unavailable: {exc}"
        self._summary.setPlainText("\n".join(summary_lines) + preview_note)

        if data is not None:
            self._populate_data_preview(np.asarray(data))
        else:
            self._clear_table(self._data_table)

        attr_items = sorted(dset.attrs.items(), key=lambda item: item[0])
        self._populate_attrs(attr_items)

    def _populate_data_preview(self, data: np.ndarray) -> None:
        if data.ndim == 0:
            self._data_table.setRowCount(1)
            self._data_table.setColumnCount(1)
            self._data_table.setHorizontalHeaderLabels(["Value"])
            self._data_table.setVerticalHeaderLabels([""])
            self._data_table.setItem(0, 0, QTableWidgetItem(_format_scalar(data.item())))
            return

        if data.ndim == 1:
            self._populate_vector_preview(data)
            return

        if data.ndim == 2:
            self._populate_matrix_preview(data)
            return

        # Fallback for higher dimensions – show flattened slice
        slice_data = data
        while slice_data.ndim > 2:
            slice_data = slice_data[0]
        message = (
            "Preview limited to the first slice along extra dimensions."
            if data.ndim > 2
            else ""
        )
        if message:
            current = self._summary.toPlainText()
            self._summary.setPlainText(current + f"\n\n{message}")
        if slice_data.ndim == 1:
            self._populate_vector_preview(slice_data)
        elif slice_data.ndim == 2:
            self._populate_matrix_preview(slice_data)
        else:
            self._clear_table(self._data_table)

    def _populate_vector_preview(self, data: np.ndarray) -> None:
        rows = min(data.shape[0], self.MAX_PREVIEW_ROWS)
        truncated = data.shape[0] > rows
        self._data_table.setColumnCount(2)
        self._data_table.setHorizontalHeaderLabels(["Index", "Value"])
        self._data_table.setRowCount(rows + (1 if truncated else 0))
        vertical_labels = []
        for row in range(rows):
            vertical_labels.append(str(row))
            self._data_table.setItem(row, 0, QTableWidgetItem(str(row)))
            self._data_table.setItem(row, 1, QTableWidgetItem(_format_scalar(data[row])))
        if truncated:
            ellipsis_row = rows
            vertical_labels.append("…")
            self._data_table.setItem(ellipsis_row, 0, QTableWidgetItem("…"))
            self._data_table.setItem(ellipsis_row, 1, QTableWidgetItem("…"))
        self._data_table.setVerticalHeaderLabels(vertical_labels)

    def _populate_matrix_preview(self, data: np.ndarray) -> None:
        rows = min(data.shape[0], self.MAX_PREVIEW_ROWS)
        cols = min(data.shape[1], self.MAX_PREVIEW_COLS)
        truncated_rows = data.shape[0] > rows
        truncated_cols = data.shape[1] > cols

        display_rows = rows + (1 if truncated_rows else 0)
        display_cols = cols + (1 if truncated_cols else 0)

        self._data_table.setRowCount(display_rows)
        self._data_table.setColumnCount(display_cols)

        horizontal_labels = [str(col) for col in range(cols)]
        if truncated_cols:
            horizontal_labels.append("…")
        self._data_table.setHorizontalHeaderLabels(horizontal_labels)

        vertical_labels = [str(row) for row in range(rows)]
        if truncated_rows:
            vertical_labels.append("…")
        self._data_table.setVerticalHeaderLabels(vertical_labels)

        for r in range(rows):
            for c in range(cols):
                self._data_table.setItem(r, c, QTableWidgetItem(_format_scalar(data[r, c])))

        if truncated_cols:
            for r in range(rows):
                self._data_table.setItem(r, cols, QTableWidgetItem("…"))
        if truncated_rows:
            for c in range(display_cols):
                self._data_table.setItem(rows, c, QTableWidgetItem("…"))

    def _populate_attrs(self, attrs: Iterable[Tuple[str, object]]) -> None:
        attrs = list(attrs)
        self._attr_table.setRowCount(len(attrs))
        if not attrs:
            self._attr_table.setVerticalHeaderLabels([])
            return
        for row, (name, value) in enumerate(attrs):
            self._attr_table.setItem(row, 0, QTableWidgetItem(str(name)))
            self._attr_table.setItem(row, 1, QTableWidgetItem(_format_attribute(value)))
        self._attr_table.resizeColumnToContents(0)

    def _clear_table(self, table: QTableWidget) -> None:
        table.clearContents()
        table.setRowCount(0)
        table.setColumnCount(0)

    # ------------------------------------------------------------------
    # Qt lifecycle
    # ------------------------------------------------------------------
    def closeEvent(self, event) -> None:  # noqa: D401 - Qt override
        """Ensure the backing HDF5 file handle is closed when the window exits."""

        try:
            if self._h5:
                self._h5.close()
        finally:
            super().closeEvent(event)


def launch_hdf_viewer(filename: str, parent: Optional[QWidget] = None) -> HDFViewWindow:
    """Create and show an :class:`HDFViewWindow` without managing the Qt event loop."""

    viewer = HDFViewWindow(filename, parent=parent)
    viewer.show()
    return viewer


def _choose_file_dialog(start_dir: Optional[str] = None) -> Optional[str]:
    options = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.ReadOnly
    directory = start_dir or os.getcwd()
    filename, _ = QFileDialog.getOpenFileName(
        None,
        "Open HDF5",
        directory,
        "HDF5 Files (*.h5 *.hdf5);;All files (*)",
        options=options,
    )
    return filename or None


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Browse the contents of an HDF5 file.")
    parser.add_argument("filename", nargs="?", help="Path to the HDF5 file to open.")
    args = parser.parse_args(argv)

    filename = args.filename
    if not filename:
        filename = _choose_file_dialog()
        if not filename:
            return 0

    app = QApplication.instance() or QApplication(sys.argv)
    try:
        viewer = HDFViewWindow(filename)
    except Exception as exc:  # pragma: no cover - UI error dialog
        QMessageBox.critical(None, "Open Error", str(exc))
        return 1
    viewer.show()
    return app.exec()


if __name__ == "__main__":  # pragma: no cover - manual execution entry
    sys.exit(main())
