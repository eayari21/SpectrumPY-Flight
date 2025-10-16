"""Qt viewer for the "IDEX CDF Variable Definitions" spreadsheet."""

from __future__ import annotations

from typing import Optional

try:  # pragma: no cover - Qt binding import guard
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import (
        QApplication,
        QLabel,
        QLineEdit,
        QListWidget,
        QListWidgetItem,
        QMainWindow,
        QPlainTextEdit,
        QSpinBox,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QVBoxLayout,
        QWidget,
    )
    QT_API = "PySide6"
except Exception:  # pragma: no cover - fallback binding
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QApplication,
        QLabel,
        QLineEdit,
        QListWidget,
        QListWidgetItem,
        QMainWindow,
        QPlainTextEdit,
        QSpinBox,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QVBoxLayout,
        QWidget,
    )
    QT_API = "PyQt6"

from idex_variable_definitions import (
    VariableDefinition,
    load_variable_definitions,
)

import numpy as np


_COEFFICIENT_HEADERS = ("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7")


class VariableDefinitionsWindow(QMainWindow):
    """Interactive browser for inspecting piecewise polynomial definitions."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        catalog = load_variable_definitions()
        self._catalog = catalog
        self._definitions = sorted(catalog.definitions, key=lambda d: d.cdf_varname.lower())
        self._current_definition: Optional[VariableDefinition] = None

        self.setWindowTitle("IDEX Variable Definitions")
        self.resize(1100, 720)

        self._build_ui()
        self._populate_list()

    def _build_ui(self) -> None:
        central = QWidget(self)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        self._search_field = QLineEdit(central)
        self._search_field.setPlaceholderText("Filter by variable name or packet field…")
        self._search_field.textChanged.connect(self._filter_definitions)
        layout.addWidget(self._search_field)

        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        layout.addWidget(splitter)

        self._list = QListWidget(splitter)
        self._list.currentItemChanged.connect(self._on_definition_selected)

        detail_widget = QWidget(splitter)
        detail_layout = QVBoxLayout(detail_widget)
        detail_layout.setContentsMargins(12, 12, 12, 12)
        detail_layout.setSpacing(10)

        self._title_label = QLabel("Select a variable to view its definition", detail_widget)
        self._title_label.setWordWrap(True)
        self._title_label.setStyleSheet("font-size: 20px; font-weight: 600;")
        detail_layout.addWidget(self._title_label)

        self._packet_label = QLabel(detail_widget)
        self._packet_label.setWordWrap(True)
        detail_layout.addWidget(self._packet_label)

        self._units_label = QLabel(detail_widget)
        detail_layout.addWidget(self._units_label)

        self._bits_label = QLabel(detail_widget)
        self._bits_label.setWordWrap(True)
        detail_layout.addWidget(self._bits_label)

        self._summary_edit = QPlainTextEdit(detail_widget)
        self._summary_edit.setReadOnly(True)
        self._summary_edit.setPlaceholderText("Plain language summary…")
        detail_layout.addWidget(self._summary_edit)

        self._notes_edit = QPlainTextEdit(detail_widget)
        self._notes_edit.setReadOnly(True)
        self._notes_edit.setPlaceholderText("Var_notes")
        detail_layout.addWidget(self._notes_edit)

        self._segments_table = QTableWidget(detail_widget)
        self._segments_table.setColumnCount(1 + len(_COEFFICIENT_HEADERS))
        headers = ["DN Range"] + list(_COEFFICIENT_HEADERS)
        self._segments_table.setHorizontalHeaderLabels(headers)
        self._segments_table.verticalHeader().setVisible(False)
        self._segments_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        detail_layout.addWidget(self._segments_table)

        self._dn_input = QSpinBox(detail_widget)
        self._dn_input.setRange(0, 4095)
        self._dn_input.setPrefix("DN: ")
        self._dn_input.valueChanged.connect(self._update_evaluated_value)
        detail_layout.addWidget(self._dn_input)

        self._evaluated_label = QLabel("Converted value: —", detail_widget)
        detail_layout.addWidget(self._evaluated_label)

        self.setCentralWidget(central)

    def _populate_list(self) -> None:
        self._list.clear()
        for definition in self._definitions:
            item = QListWidgetItem(definition.cdf_varname or definition.packet_variable)
            item.setData(Qt.ItemDataRole.UserRole, definition)
            item.setToolTip(f"Packet field: {definition.packet_variable}\nUnits: {definition.units or '—'}")
            self._list.addItem(item)
        if self._list.count() > 0:
            self._list.setCurrentRow(0)

    def _filter_definitions(self, text: str) -> None:
        text = text.strip().lower()
        for index in range(self._list.count()):
            item = self._list.item(index)
            definition = item.data(Qt.ItemDataRole.UserRole)
            haystack = " ".join(
                part for part in (
                    definition.cdf_varname,
                    definition.packet_variable,
                    definition.var_notes,
                )
                if part
            ).lower()
            item.setHidden(text not in haystack)
        self._select_first_visible()

    def _select_first_visible(self) -> None:
        for index in range(self._list.count()):
            if not self._list.item(index).isHidden():
                self._list.setCurrentRow(index)
                return
        self._current_definition = None
        self._clear_details()

    def _on_definition_selected(self, current: Optional[QListWidgetItem], _previous: Optional[QListWidgetItem]) -> None:
        if current is None or current.isHidden():
            self._current_definition = None
            self._clear_details()
            return
        definition = current.data(Qt.ItemDataRole.UserRole)
        if not isinstance(definition, VariableDefinition):
            self._current_definition = None
            self._clear_details()
            return

        self._current_definition = definition
        self._title_label.setText(f"{definition.cdf_varname} — {definition.units or 'dimensionless'}")
        packet = definition.packet_variable or "—"
        self._packet_label.setText(f"Packet field: <b>{packet}</b>")
        dn_span = definition.dn_range()
        if dn_span:
            span_text = f"DN range covered: {dn_span[0]} – {dn_span[1]}"
        else:
            span_text = "DN range covered: (not specified)"
        self._bits_label.setText(
            " | ".join(
                part
                for part in (
                    span_text,
                    f"Start bit: {definition.starting_bit}" if definition.starting_bit is not None else "",
                    f"Padding bits: {definition.padding_bits}" if definition.padding_bits is not None else "",
                    f"Unsigned bits: {definition.unsigned_bits}" if definition.unsigned_bits is not None else "",
                )
                if part
            )
        )
        self._units_label.setText(f"Units: {definition.units or '—'}")
        self._summary_edit.setPlainText(definition.summary or "")
        self._notes_edit.setPlainText(definition.var_notes or "")
        self._populate_segments(definition)
        self._configure_dn_input(definition)
        self._update_evaluated_value(self._dn_input.value())

    def _populate_segments(self, definition: VariableDefinition) -> None:
        self._segments_table.setRowCount(len(definition.segments))
        for row, segment in enumerate(definition.segments):
            dn_text = f"{segment.dn_range[0]} – {segment.dn_range[1]}"
            self._segments_table.setItem(row, 0, QTableWidgetItem(dn_text))
            for col, coeff in enumerate(segment.coefficients, start=1):
                display = f"{coeff:.6g}" if coeff else "0"
                self._segments_table.setItem(row, col, QTableWidgetItem(display))
        self._segments_table.resizeColumnsToContents()

    def _configure_dn_input(self, definition: VariableDefinition) -> None:
        dn_span = definition.dn_range()
        if dn_span:
            self._dn_input.blockSignals(True)
            self._dn_input.setRange(dn_span[0], dn_span[1])
            self._dn_input.setValue(dn_span[0])
            self._dn_input.blockSignals(False)
        else:
            self._dn_input.blockSignals(True)
            self._dn_input.setRange(0, 4095)
            self._dn_input.setValue(0)
            self._dn_input.blockSignals(False)

    def _update_evaluated_value(self, dn_value: int) -> None:
        if self._current_definition is None:
            self._evaluated_label.setText("Converted value: —")
            return
        converted = self._current_definition.evaluate(dn_value)
        if isinstance(converted, np.ndarray):
            converted = converted.item()
        units = self._current_definition.units or ""
        self._evaluated_label.setText(f"Converted value: {converted:.6g} {units}")

    def _clear_details(self) -> None:
        self._title_label.setText("Select a variable to view its definition")
        self._packet_label.clear()
        self._units_label.clear()
        self._bits_label.clear()
        self._summary_edit.clear()
        self._notes_edit.clear()
        self._segments_table.setRowCount(0)
        self._evaluated_label.setText("Converted value: —")


def launch_variable_definitions_viewer(parent: Optional[QWidget] = None) -> VariableDefinitionsWindow:
    """Instantiate and return the viewer window."""

    window = VariableDefinitionsWindow(parent=parent)
    window.show()
    return window


__all__ = ["VariableDefinitionsWindow", "launch_variable_definitions_viewer", "QT_API"]
