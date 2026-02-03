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
import csv
import os
import sys
from typing import Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# Qt binding import (prefer PySide6, fallback to PyQt6)
# ---------------------------------------------------------------------------
_QT_API = None
try:  # pragma: no cover - import guard
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import (
        QApplication,
        QFileDialog,
        QHeaderView,
        QLineEdit,
        QMainWindow,
        QMessageBox,
        QPlainTextEdit,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QTabWidget,
        QToolBar,
        QTreeWidget,
        QTreeWidgetItem,
        QVBoxLayout,
        QWidget,
    )
    _QT_API = "PySide6"
except Exception:  # pragma: no cover - import guard
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import (
        QApplication,
        QFileDialog,
        QHeaderView,
        QLineEdit,
        QMainWindow,
        QMessageBox,
        QPlainTextEdit,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QTabWidget,
        QToolBar,
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
    TIME_DATASET_CANDIDATES = (
        "/Metadata/unpacked/utc_timestamp",
        "/Metadata/raw/utc_timestamp",
        "/Metadata/utc_timestamp",
    )

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

        toolbar = QToolBar("HDF Tools", self)
        toolbar.setMovable(False)
        self._search_edit = QLineEdit(self)
        self._search_edit.setPlaceholderText("Search datasets and groups…")
        self._search_edit.textChanged.connect(self._apply_filter)
        toolbar.addWidget(self._search_edit)
        export_action = QAction("Export Selected to CSV", self)
        export_action.setToolTip("Export selected datasets with utc_timestamp to CSV")
        export_action.triggered.connect(self._export_selected_to_csv)
        toolbar.addAction(export_action)
        self.addToolBar(toolbar)

        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        splitter.setChildrenCollapsible(False)
        if hasattr(splitter, "setCollapsible"):
            splitter.setCollapsible(0, False)
            splitter.setCollapsible(1, False)
        layout.addWidget(splitter)

        self._tree = QTreeWidget(splitter)
        self._tree.setHeaderLabels(["Name", "Type", "Shape"])
        self._tree.setAlternatingRowColors(True)
        self._tree.setSelectionMode(QTreeWidget.SelectionMode.ExtendedSelection)
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

    def _apply_filter(self, text: str) -> None:
        query = text.strip().lower()
        root = self._tree.topLevelItem(0)
        if root is None:
            return
        if not query:
            self._reset_filter(root)
            return
        self._filter_item(root, query)

    def _reset_filter(self, item: QTreeWidgetItem) -> None:
        item.setHidden(False)
        for index in range(item.childCount()):
            self._reset_filter(item.child(index))

    def _filter_item(self, item: QTreeWidgetItem, query: str) -> bool:
        name = item.text(0).lower()
        path = str(item.data(0, Qt.ItemDataRole.UserRole) or "").lower()
        match = query in name or (path and query in path)
        child_match = False
        for index in range(item.childCount()):
            child = item.child(index)
            child_match = self._filter_item(child, query) or child_match
        visible = match or child_match
        item.setHidden(not visible)
        if child_match:
            item.setExpanded(True)
        return visible

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

    def _find_time_dataset(self) -> Tuple[np.ndarray, str]:
        for path in self.TIME_DATASET_CANDIDATES:
            if path in self._h5:
                return np.asarray(self._h5[path][()]), path

        matches: List[str] = []

        def visitor(name: str, obj: h5py.Dataset) -> None:
            if isinstance(obj, h5py.Dataset):
                dataset_name = obj.name.split("/")[-1].lower()
                if dataset_name == "utc_timestamp":
                    matches.append(obj.name)

        self._h5.visititems(visitor)
        if not matches:
            raise ValueError("utc_timestamp dataset not found in the HDF5 file.")
        path = matches[0]
        return np.asarray(self._h5[path][()]), path

    def _coerce_vector(self, data: np.ndarray, label: str) -> np.ndarray:
        if data.ndim == 0:
            return np.array([data.item()])
        if data.ndim != 1:
            raise ValueError(f"{label} must be a 1D dataset to export.")
        return data

    def _selected_dataset_paths(self) -> List[str]:
        datasets: List[str] = []
        for item in self._tree.selectedItems():
            if item.data(1, Qt.ItemDataRole.UserRole) == "dataset":
                path = item.data(0, Qt.ItemDataRole.UserRole)
                if path:
                    datasets.append(path)
        return datasets

    def _export_selected_to_csv(self) -> None:
        datasets = self._selected_dataset_paths()
        if not datasets:
            QMessageBox.information(self, "Export to CSV", "Select one or more datasets to export.")
            return

        try:
            time_data, time_path = self._find_time_dataset()
            time_vector = self._coerce_vector(np.asarray(time_data), time_path)
        except Exception as exc:
            QMessageBox.critical(self, "Export Error", f"Unable to locate utc_timestamp.\n{exc}")
            return

        if len(datasets) == 1:
            dataset_path = datasets[0]
            default_name = os.path.basename(dataset_path.strip("/")) or "export"
            filename, _ = QFileDialog.getSaveFileName(
                self,
                "Export CSV",
                f"{default_name}.csv",
                "CSV Files (*.csv);;All files (*)",
                options=QFileDialog.Option.DontUseNativeDialog,
            )
            if not filename:
                return
            try:
                self._export_dataset_to_csv(dataset_path, time_vector, filename)
            except Exception as exc:
                QMessageBox.critical(self, "Export Error", str(exc))
                return
            QMessageBox.information(self, "Export Complete", f"Saved {filename}")
            return

        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Export Folder",
            os.getcwd(),
            options=QFileDialog.Option.DontUseNativeDialog,
        )
        if not directory:
            return
        exported: List[str] = []
        errors: List[str] = []
        for dataset_path in datasets:
            filename = os.path.join(directory, f"{os.path.basename(dataset_path.strip('/'))}.csv")
            try:
                self._export_dataset_to_csv(dataset_path, time_vector, filename)
                exported.append(filename)
            except Exception as exc:
                errors.append(f"{dataset_path}: {exc}")
        message = f"Exported {len(exported)} file(s)."
        if errors:
            message += "\n\nErrors:\n" + "\n".join(errors)
        QMessageBox.information(self, "Export Complete", message)

    def _export_dataset_to_csv(self, dataset_path: str, time_vector: np.ndarray, filename: str) -> None:
        if dataset_path not in self._h5:
            raise ValueError(f"Dataset not found: {dataset_path}")
        data = np.asarray(self._h5[dataset_path][()])
        values = self._coerce_vector(data, dataset_path)
        if values.shape[0] != time_vector.shape[0]:
            raise ValueError(
                "Dataset length does not match utc_timestamp length "
                f"({values.shape[0]} vs {time_vector.shape[0]})."
            )
        header_name = os.path.basename(dataset_path.strip("/")) or dataset_path
        with open(filename, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["utc_timestamp", header_name])
            for time_value, data_value in zip(time_vector, values):
                writer.writerow([_format_scalar(time_value), _format_scalar(data_value)])

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
        rows = data.shape[0]
        self._data_table.setColumnCount(2)
        self._data_table.setHorizontalHeaderLabels(["Index", "Value"])
        self._data_table.setRowCount(rows)
        vertical_labels = []
        for row in range(rows):
            vertical_labels.append(str(row))
            self._data_table.setItem(row, 0, QTableWidgetItem(str(row)))
            self._data_table.setItem(row, 1, QTableWidgetItem(_format_scalar(data[row])))
        self._data_table.setVerticalHeaderLabels(vertical_labels)

    def _populate_matrix_preview(self, data: np.ndarray) -> None:
        rows = data.shape[0]
        cols = data.shape[1]
        self._data_table.setRowCount(rows)
        self._data_table.setColumnCount(cols)

        horizontal_labels = [str(col) for col in range(cols)]
        self._data_table.setHorizontalHeaderLabels(horizontal_labels)

        vertical_labels = [str(row) for row in range(rows)]
        self._data_table.setVerticalHeaderLabels(vertical_labels)

        for r in range(rows):
            for c in range(cols):
                self._data_table.setItem(r, c, QTableWidgetItem(_format_scalar(data[r, c])))

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
    parser.add_argument(
        "path",
        nargs="?",
        help="Path to the HDF5 file to open.",
    )
    parser.add_argument(
        "-f",
        "--file",
        "--filename",
        dest="filename",
        help="Explicit path to the HDF5 file to open.",
    )
    args = parser.parse_args(argv)

    filename = args.filename or args.path
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
