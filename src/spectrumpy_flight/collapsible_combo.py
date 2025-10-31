"""Reusable Qt combo box with collapsible category headings."""

from __future__ import annotations

from typing import Any, List, MutableMapping, Optional, Sequence, Tuple

_QT_AVAILABLE = True
try:  # pragma: no cover - prefer PySide6 when available
    from PySide6.QtCore import QEvent, QModelIndex, Qt
    from PySide6.QtGui import QFont, QPalette
    from PySide6.QtWidgets import QAbstractItemView, QComboBox, QListView
except Exception:  # pragma: no cover - fallback to PyQt6
    try:
        from PyQt6.QtCore import QEvent, QModelIndex, Qt
        from PyQt6.QtGui import QFont, QPalette
        from PyQt6.QtWidgets import QAbstractItemView, QComboBox, QListView
    except Exception:  # pragma: no cover - final fallback for test environments
        _QT_AVAILABLE = False


if _QT_AVAILABLE:

    class CollapsibleComboBox(QComboBox):
        """A ``QComboBox`` that supports collapsible category headings."""

        _HEADING_ROLE = int(Qt.ItemDataRole.UserRole) + 101
        _GROUP_ROLE = int(Qt.ItemDataRole.UserRole) + 102

        def __init__(self, parent: Optional[object] = None) -> None:
            super().__init__(parent)  # type: ignore[arg-type]
            self._reset_groups()

            view = QListView(self)
            view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            view.setUniformItemSizes(True)
            view.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
            view.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
            view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            view.viewport().installEventFilter(self)
            view.setStyleSheet(
                """
                QListView {
                    padding: 4px 0;
                }
                QListView::item {
                    padding: 6px 14px;
                }
                """
            )
            self.setView(view)
            self.setMaxVisibleItems(18)

        # ------------------------------------------------------------------
        def clear(self) -> None:  # type: ignore[override]
            super().clear()
            self._reset_groups()

        # ------------------------------------------------------------------
        def eventFilter(self, obj: object, event: QEvent) -> bool:
            if obj is self.view().viewport() and event.type() == QEvent.Type.MouseButtonPress:
                if hasattr(event, "position"):
                    pos = event.position().toPoint()  # type: ignore[attr-defined]
                else:  # pragma: no cover - compatibility guard
                    pos = getattr(event, "pos")()
                index = self.view().indexAt(pos)
                if index.isValid() and bool(index.data(self._HEADING_ROLE)):
                    group_id = index.data(self._GROUP_ROLE)
                    if isinstance(group_id, int):
                        self._toggle_group(group_id)
                        return True
            return super().eventFilter(obj, event)

        # ------------------------------------------------------------------
        def set_grouped_items(
            self,
            groups: Sequence[Tuple[str, Sequence[Any]]],
            *,
            placeholder: Optional[str] = None,
            placeholder_data: Any = None,
            preserve_selection: bool = True,
        ) -> None:
            current_text = self.currentText() if preserve_selection else ""
            self.blockSignals(True)
            super().clear()
            self._reset_groups()

            if placeholder is not None:
                super().addItem(placeholder, placeholder_data)
                self._placeholder_index = 0
            else:
                self._placeholder_index = None

            for group_index, (title, entries) in enumerate(groups):
                if group_index > 0 or self._placeholder_index is not None and self.count() > 0:
                    separator_row = self.count()
                    self.insertSeparator(separator_row)
                    self._separator_rows[group_index] = separator_row
                else:
                    self._separator_rows[group_index] = None

                self._add_heading(title, group_index)
                for entry in entries:
                    if isinstance(entry, (tuple, list)) and entry:
                        text = str(entry[0])
                        data = entry[1] if len(entry) > 1 else entry[0]
                    else:
                        text = str(entry)
                        data = entry
                    self._add_group_item(text, data, group_index)

            if preserve_selection and current_text:
                match_index = self.findText(current_text)
                if match_index != -1:
                    self.setCurrentIndex(match_index)
                    self._ensure_index_visible(match_index)
                elif self._placeholder_index is not None:
                    self.setCurrentIndex(self._placeholder_index)
            elif self._placeholder_index is not None and self.count() > 0:
                self.setCurrentIndex(self._placeholder_index)

            self.blockSignals(False)
            for group_id in list(self._group_children):
                self._refresh_group_visibility(group_id)
                self._update_heading_text(group_id)

        # ------------------------------------------------------------------
        def _reset_groups(self) -> None:
            self._group_children: MutableMapping[int, List[int]] = {}
            self._group_headings: MutableMapping[int, int] = {}
            self._group_titles: MutableMapping[int, str] = {}
            self._group_state: MutableMapping[int, bool] = {}
            self._separator_rows: MutableMapping[int, Optional[int]] = {}
            self._placeholder_index: Optional[int] = None

        # ------------------------------------------------------------------
        def _model_index(self, row: int) -> QModelIndex:
            return self.model().index(row, self.modelColumn(), self.rootModelIndex())

        # ------------------------------------------------------------------
        def _add_heading(self, title: str, group_id: int) -> None:
            row = self.count()
            super().addItem(title)
            index = self._model_index(row)
            font: QFont = self.font()
            font.setBold(True)
            self.model().setData(index, font, Qt.ItemDataRole.FontRole)
            muted = self.palette().color(QPalette.ColorRole.Mid)
            self.model().setData(index, muted, Qt.ItemDataRole.ForegroundRole)
            self.model().setData(index, True, self._HEADING_ROLE)
            self.model().setData(index, group_id, self._GROUP_ROLE)
            self.model().setData(index, 0, Qt.ItemDataRole.UserRole - 1)
            self._group_children[group_id] = []
            self._group_headings[group_id] = row
            self._group_titles[group_id] = title
            self._group_state[group_id] = True
            self._update_heading_text(group_id)

        # ------------------------------------------------------------------
        def _add_group_item(self, text: str, data: Any, group_id: int) -> None:
            super().addItem(text, data)
            row = self.count() - 1
            self._group_children.setdefault(group_id, []).append(row)
            index = self._model_index(row)
            self.model().setData(index, group_id, self._GROUP_ROLE)

        # ------------------------------------------------------------------
        def _toggle_group(self, group_id: int) -> None:
            expanded = self._group_state.get(group_id, True)
            self._group_state[group_id] = not expanded
            self._refresh_group_visibility(group_id)
            self._update_heading_text(group_id)

        # ------------------------------------------------------------------
        def _refresh_group_visibility(self, group_id: int) -> None:
            expanded = self._group_state.get(group_id, True)
            hidden = not expanded
            for row in self._group_children.get(group_id, []):
                self._set_row_hidden(row, hidden)
            separator_row = self._separator_rows.get(group_id)
            if separator_row is not None:
                self._set_row_hidden(separator_row, hidden)

        # ------------------------------------------------------------------
        def _update_heading_text(self, group_id: int) -> None:
            row = self._group_headings.get(group_id)
            if row is None:
                return
            index = self._model_index(row)
            title = self._group_titles.get(group_id, "")
            arrow = "▾" if self._group_state.get(group_id, True) else "▸"
            self.model().setData(index, f"{arrow} {title}", Qt.ItemDataRole.DisplayRole)

        # ------------------------------------------------------------------
        def _ensure_index_visible(self, index: int) -> None:
            model_index = self._model_index(index)
            group_id = model_index.data(self._GROUP_ROLE)
            if isinstance(group_id, int) and not self._group_state.get(group_id, True):
                self._group_state[group_id] = True
                self._refresh_group_visibility(group_id)
                self._update_heading_text(group_id)

        # ------------------------------------------------------------------
        def _set_row_hidden(self, row: int, hidden: bool) -> None:
            model = self.model()
            if model is None:
                return

            item = model.item(row)
            if item is not None and hasattr(item, "setHidden"):
                item.setHidden(hidden)
                return

            view = self.view()
            if hasattr(view, "setRowHidden"):
                view.setRowHidden(row, hidden)
                return

            index = self._model_index(row)
            role = Qt.ItemDataRole.UserRole - 1
            model.setData(index, 0 if hidden else 1, role)


else:  # pragma: no cover - fallback stub for environments without Qt bindings

    class CollapsibleComboBox:
        """Minimal stand-in used when Qt bindings are unavailable."""

        def __init__(self, *args: object, **kwargs: object) -> None:
            self._items: List[str] = []
            self._data: List[Any] = []
            self._current_index: int = -1
            self._editable = False

        def set_grouped_items(
            self,
            groups: Sequence[Tuple[str, Sequence[Any]]],
            *,
            placeholder: Optional[str] = None,
            placeholder_data: Any = None,
            preserve_selection: bool = True,
        ) -> None:
            previous = self.currentText() if preserve_selection else ""
            self.clear()
            if placeholder is not None:
                self.addItem(placeholder, placeholder_data)
                self._current_index = 0
            for _, entries in groups:
                for entry in entries:
                    if isinstance(entry, (tuple, list)) and entry:
                        text = str(entry[0])
                        data = entry[1] if len(entry) > 1 else entry[0]
                    else:
                        text = str(entry)
                        data = entry
                    self.addItem(text, data)
            if preserve_selection and previous:
                self.setCurrentText(previous)

        def addItem(self, text: str, userData: Any = None) -> None:
            self._items.append(str(text))
            self._data.append(userData)

        def clear(self) -> None:
            self._items.clear()
            self._data.clear()
            self._current_index = -1

        def count(self) -> int:
            return len(self._items)

        def currentText(self) -> str:
            if 0 <= self._current_index < len(self._items):
                return self._items[self._current_index]
            return ""

        def findText(self, text: str) -> int:
            try:
                return self._items.index(text)
            except ValueError:
                return -1

        def setCurrentIndex(self, index: int) -> None:
            if 0 <= index < len(self._items):
                self._current_index = index

        def setCurrentText(self, text: str) -> None:
            index = self.findText(text)
            if index != -1:
                self._current_index = index

        def itemData(self, index: int) -> Any:
            if 0 <= index < len(self._data):
                return self._data[index]
            return None

        # No-op compatibility methods -------------------------------------
        def setEditable(self, editable: bool) -> None:
            self._editable = editable

        def setInsertPolicy(self, *args: object, **kwargs: object) -> None:
            pass

        def setSizeAdjustPolicy(self, *args: object, **kwargs: object) -> None:
            pass

        def setEnabled(self, *args: object, **kwargs: object) -> None:
            pass

        def setStyleSheet(self, *args: object, **kwargs: object) -> None:
            pass

        def setToolTip(self, *args: object, **kwargs: object) -> None:
            pass

        def blockSignals(self, *args: object, **kwargs: object) -> None:
            pass

        def setMaxVisibleItems(self, *args: object, **kwargs: object) -> None:
            pass

        def setView(self, *args: object, **kwargs: object) -> None:
            pass
