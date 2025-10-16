#!/opt/anaconda3/envs/idex/bin/python
# -*- coding: utf-8 -*-
from __future__ import annotations

# =====================================================================
# IDEX Quicklook — QtAgg canvas + Matplotlib Navigation Toolbar
# Interactive viewer with channel selection, overlays, and fit editing tools
# =====================================================================

import os
import sys
import argparse
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

import textwrap

try:
    import h5py  # type: ignore
except Exception:  # pragma: no cover - optional dependency for environments without h5py
    h5py = None

import numpy as np

from HDF_View import launch_hdf_viewer
from dust_composition import launch_dust_composition_window
try:  # pragma: no cover - optional dependency, loaded lazily
    from HDF_View import launch_hdf_viewer
except Exception:  # pragma: no cover
    launch_hdf_viewer = None

try:  # pragma: no cover - optional dependency, loaded lazily
    from CDF_View import launch_cdf_viewer
except Exception:  # pragma: no cover
    launch_cdf_viewer = None

try:  # pragma: no cover - optional dependency
    from IDEX_Definitions_View import launch_variable_definitions_viewer
except Exception:  # pragma: no cover
    launch_variable_definitions_viewer = None

# Qt-backed Matplotlib so NavigationToolbar works
import matplotlib
from matplotlib import mathtext
matplotlib.use("QtAgg")

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

# --------- Qt binding-agnostic imports (prefer PySide6, fallback PyQt6) ---------
_QT = None
try:
    from PySide6.QtCore import Qt, QSize, QTimer, QUrl
    from PySide6.QtGui import QAction, QFont, QIcon, QPixmap, QImage, QTextCursor, QTextDocument
    from PySide6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy, QDialog, QPushButton,
        QHBoxLayout, QGridLayout, QTableWidget, QTableWidgetItem, QHeaderView,
        QCheckBox, QDialogButtonBox, QMenu, QMenuBar, QToolButton, QTextBrowser,
        QListWidget, QListWidgetItem, QLineEdit, QWidgetAction, QStyle, QSplitter
    )
    _QT = "PySide6"
except Exception:
    from PyQt6.QtCore import Qt, QSize, QTimer, QUrl
    from PyQt6.QtGui import QAction, QFont, QIcon, QPixmap, QImage, QTextCursor, QTextDocument
    from PyQt6.QtWidgets import (
        QApplication, QFileDialog, QMainWindow, QMessageBox, QStatusBar, QToolBar,
        QVBoxLayout, QWidget, QComboBox, QLabel, QSizePolicy, QDialog, QPushButton,
        QHBoxLayout, QGridLayout, QTableWidget, QTableWidgetItem, QHeaderView,
        QCheckBox, QDialogButtonBox, QMenu, QMenuBar, QToolButton, QTextBrowser,
        QListWidget, QListWidgetItem, QLineEdit, QWidgetAction, QStyle, QSplitter
    )
    _QT = "PyQt6"

print(f"[info] Qt binding: {_QT}, Matplotlib backend: {matplotlib.get_backend()}")

matplotlib.rcParams.update({
    "figure.autolayout": False,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "lines.linewidth": 1.8,
})

IMAGES_DIR = Path(__file__).resolve().parent / "Images"
IMAP_LOGO_CANDIDATES = (
    "IMAP_logo.png",
    "IMAP_logo.jpg",
    "IMAP_logo.jpeg",
    "imap_logo.png",
    "IMAP.png",
    "IMAP.jpeg",
)


def _find_brand_logo() -> Optional[Path]:
    for candidate in IMAP_LOGO_CANDIDATES:
        logo_path = IMAGES_DIR / candidate
        if logo_path.exists():
            return logo_path
    return None


def _load_brand_pixmap(*, max_height: int = 72) -> Optional[QPixmap]:
    logo_path = _find_brand_logo()
    if logo_path is None:
        return None

    pixmap = QPixmap(str(logo_path))
    if pixmap.isNull():
        return None

    # ``scaledToHeight`` already preserves the aspect ratio; in PyQt6/PySide6 it only
    # accepts the target height and an optional transformation mode.  Passing an
    # aspect-ratio flag raises a ``TypeError`` under PyQt6, so we only supply the
    # desired height and transformation mode here.
    return pixmap.scaledToHeight(
        max_height,
        Qt.TransformationMode.SmoothTransformation,
    )

# --------- Small utils ---------
def prompt_for_data_file(
    parent: QWidget,
    start_dir: str | None = None,
    *,
    preferred: Optional[str] = None,
) -> Optional[str]:
    """Show a non-native file dialog that accepts HDF5 and CDF payloads."""

    if start_dir is None:
        repo_root = Path(__file__).resolve().parent
        preferred_map = {
            "cdf": repo_root / "CDF",
            "hdf5": repo_root / "HDF5",
        }
        if preferred and (target_dir := preferred_map.get(preferred.lower())) and target_dir.exists():
            start_dir = str(target_dir)
        else:
            default_dir = repo_root / "HDF5"
            if not default_dir.exists():
                default_dir = repo_root
            start_dir = str(default_dir)

    options = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.ReadOnly

    filters = [
        "Data Files (*.h5 *.hdf5 *.cdf)",
        "HDF5 Files (*.h5 *.hdf5)",
        "CDF Files (*.cdf)",
        "All files (*)",
    ]

    if preferred == "hdf5":
        filters = [filters[1], filters[0], filters[2], filters[3]]
    elif preferred == "cdf":
        filters = [filters[2], filters[0], filters[1], filters[3]]

    filename, _ = QFileDialog.getOpenFileName(
        parent,
        "Open Data File",
        start_dir,
        ";;".join(filters),
        options=options,
    )
    return filename or None

def _normalize_name(value: str) -> str:
    return "".join(ch.lower() for ch in value if ch.isalnum())

def _pair_key(path: str) -> str:
    base = _normalize_name(os.path.basename(path))
    for token in ("time", "result", "fit", "curve"):
        base = base.replace(token, "")
    return base

def _friendly_label(path: str) -> str:
    base = os.path.basename(path)
    replacements = {
        "FitResult": "Fit",
        "fitresult": "Fit",
        "Result": "Fit",
    }
    for needle, repl in replacements.items():
        base = base.replace(needle, repl)
    return base


def _label_from_param_path(path: str) -> str:
    base = os.path.basename(path)
    replacements = (
        ("FitParameters", "Fit"),
        ("FitParams", "Fit"),
        ("fitparameters", "Fit"),
        ("fitparams", "Fit"),
        ("Parameters", "Fit"),
        ("Params", "Fit"),
    )
    for needle, repl in replacements:
        if needle in base:
            base = base.replace(needle, repl)
    if "Fit" in base and " Fit" not in base:
        base = base.replace("Fit", " Fit")
    base = " ".join(base.split()).strip()
    if not base:
        base = "Fit"
    if "Fit" not in base:
        base = f"{base} Fit"
    return f"{base} (calc)"


@dataclass(frozen=True)
class DocumentationEntry:
    path: Path
    title: str
    text: str
    lines: List[str]

    @property
    def display_name(self) -> str:
        return self.title or self.path.stem.replace("_", " ").title()

    @property
    def relative_path(self) -> str:
        try:
            return str(self.path.relative_to(Path(__file__).resolve().parent))
        except ValueError:
            return str(self.path)


def _iter_document_paths() -> List[Path]:
    root = Path(__file__).resolve().parent
    candidates: List[Path] = []
    readme = root / "README.md"
    if readme.exists():
        candidates.append(readme)
    docs_dir = root / "docs"
    if docs_dir.exists():
        candidates.extend(sorted(docs_dir.rglob("*.md")))
    return candidates


def _derive_title(path: Path, text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            return stripped.lstrip("# ").strip()
    return path.stem.replace("_", " ").title()


def _load_documentation_entries() -> List[DocumentationEntry]:
    entries: List[DocumentationEntry] = []
    for path in _iter_document_paths():
        try:
            text = path.read_text(encoding="utf-8")
        except Exception:
            continue
        lines = text.splitlines()
        title = _derive_title(path, text)
        entries.append(DocumentationEntry(path=path, title=title, text=text, lines=lines))
    return entries


_DOCUMENTATION_ENTRIES: List[DocumentationEntry] = _load_documentation_entries()


class DocumentationCenter(QDialog):
    """A lightweight reader with full-text search across project documentation."""

    def __init__(self, parent: Optional[QWidget] = None, *, initial_query: str = "") -> None:
        super().__init__(parent)
        self.setWindowTitle("SpectrumPY Flight Documentation Center")
        self.resize(960, 640)
        self.setModal(False)
        self._documents = _DOCUMENTATION_ENTRIES
        self._current_entry: Optional[DocumentationEntry] = None
        self._current_query: str = ""

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        heading = QLabel("SpectrumPY Documentation & Tutorials", self)
        heading.setObjectName("docHeading")
        heading.setStyleSheet("font-size: 22px; font-weight: 600;")
        layout.addWidget(heading)

        description = QLabel(
            "Search the bundled README and guides, or browse tutorials with a single click.",
            self,
        )
        description.setStyleSheet("font-size: 14px; color: #4a5568;")
        layout.addWidget(description)

        search_row = QHBoxLayout()
        search_row.setSpacing(10)

        self.search_field = QLineEdit(self)
        self.search_field.setPlaceholderText("Search documentation…")
        self.search_field.setClearButtonEnabled(True)
        self.search_field.setMinimumWidth(280)
        self.search_field.returnPressed.connect(self.perform_search)
        search_row.addWidget(self.search_field, stretch=1)

        search_button = QPushButton("Search", self)
        search_button.setMinimumHeight(32)
        search_button.clicked.connect(self.perform_search)
        search_row.addWidget(search_button)

        show_all_button = QPushButton("Show All", self)
        show_all_button.setMinimumHeight(32)
        show_all_button.clicked.connect(self.show_all_documents)
        search_row.addWidget(show_all_button)

        layout.addLayout(search_row)

        self.result_label = QLabel("Browse documentation", self)
        self.result_label.setStyleSheet("font-size: 13px; color: #4a5568;")
        layout.addWidget(self.result_label)

        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        splitter.setChildrenCollapsible(False)
        layout.addWidget(splitter, stretch=1)

        self.results = QListWidget(splitter)
        self.results.setAlternatingRowColors(True)
        self.results.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.results.setMinimumWidth(280)
        self.results.itemSelectionChanged.connect(self._on_result_selected)
        self.results.itemActivated.connect(self._on_result_activated)

        self.viewer = QTextBrowser(splitter)
        self.viewer.setOpenExternalLinks(True)
        self.viewer.setStyleSheet("font-size: 14px; line-height: 1.5em; background: #fefefe; padding: 12px; border-radius: 12px;")

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        self.show_all_documents()

        if initial_query:
            self.search_field.setText(initial_query)
            self.perform_search()

    # ------------------------------------------------------------------
    def focus_search(self) -> None:
        self.search_field.setFocus()
        self.search_field.selectAll()

    def set_query(self, value: str) -> None:
        self.search_field.setText(value)
        if value:
            self.perform_search()
        else:
            self.show_all_documents()

    def _entry_for_path(self, path_str: str) -> Optional[DocumentationEntry]:
        path = Path(path_str).resolve()
        for entry in self._documents:
            if entry.path.resolve() == path:
                return entry
        return None

    def show_all_documents(self) -> None:
        self.results.clear()
        self._current_query = ""
        for entry in self._documents:
            item = QListWidgetItem(f"{entry.display_name}")
            item.setData(Qt.ItemDataRole.UserRole, (str(entry.path), 0, ""))
            self.results.addItem(item)
        self.result_label.setText("Browse documentation")
        if self._documents:
            self.results.setCurrentRow(0)

    def perform_search(self) -> None:
        query = self.search_field.text().strip()
        self.results.clear()
        if not query:
            self.show_all_documents()
            return

        lowered = query.lower()
        matches: List[Tuple[DocumentationEntry, int, str]] = []
        for entry in self._documents:
            for idx, line in enumerate(entry.lines, start=1):
                if lowered in line.lower():
                    snippet = textwrap.shorten(line.strip(), width=100, placeholder="…")
                    matches.append((entry, idx, snippet))

        if not matches:
            item = QListWidgetItem("No matches found")
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            self.results.addItem(item)
            self.result_label.setText(f"No results for '{query}'.")
            self.viewer.setPlainText("")
            return

        self._current_query = query
        for entry, line_no, snippet in matches:
            display = f"{entry.display_name} — L{line_no}: {snippet}"
            item = QListWidgetItem(display)
            item.setData(Qt.ItemDataRole.UserRole, (str(entry.path), line_no, query))
            self.results.addItem(item)
        self.result_label.setText(f"{len(matches)} match(es) for '{query}'.")
        self.results.setCurrentRow(0)

    def _display_entry(self, entry: DocumentationEntry, *, line_no: int = 0, highlight: str = "") -> None:
        self._current_entry = entry
        document = self.viewer.document()
        if hasattr(document, "setMarkdownFeatures"):
            feature = getattr(QTextDocument.MarkdownFeature, "MarkdownDialectGitHub", None)
            if feature is None:
                feature = getattr(QTextDocument.MarkdownFeature, "MarkdownDialectCommonMark", None)
            if feature is not None:
                document.setMarkdownFeatures(feature)
        if hasattr(document, "setBaseUrl"):
            document.setBaseUrl(QUrl.fromLocalFile(str(entry.path.parent)))
        header = f"**{entry.display_name}**  \n`{entry.relative_path}`\n\n{entry.text}"
        if hasattr(self.viewer, "setMarkdown"):
            self.viewer.setMarkdown(header)
        else:
            plain = f"{entry.display_name}\n{entry.relative_path}\n\n{entry.text}"
            self.viewer.setPlainText(plain)
        if line_no > 0:
            cursor = self.viewer.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            for _ in range(line_no + 2):
                cursor.movePosition(QTextCursor.MoveOperation.Down)
            self.viewer.setTextCursor(cursor)
            self.viewer.ensureCursorVisible()
        if highlight:
            cursor = self.viewer.document().find(
                highlight,
                self.viewer.textCursor(),
                QTextDocument.FindFlag.FindCaseSensitively,
            )
            if cursor.isNull():
                cursor = self.viewer.document().find(
                    highlight,
                    0,
                    QTextDocument.FindFlag.FindCaseSensitively,
                )
            if cursor.isNull():
                cursor = self.viewer.document().find(highlight)
            if not cursor.isNull():
                self.viewer.setTextCursor(cursor)
                self.viewer.ensureCursorVisible()

    def _on_result_selected(self) -> None:
        item = self.results.currentItem()
        if item is None:
            return
        data = item.data(Qt.ItemDataRole.UserRole)
        if not data:
            return
        path_str, line_no, highlight = data
        entry = self._entry_for_path(path_str)
        if entry:
            self._display_entry(entry, line_no=line_no, highlight=highlight)

    def _on_result_activated(self, item: QListWidgetItem) -> None:
        data = item.data(Qt.ItemDataRole.UserRole)
        if not data:
            return
        path_str, line_no, highlight = data
        entry = self._entry_for_path(path_str)
        if entry:
            self._display_entry(entry, line_no=line_no, highlight=highlight)


_MATH_TEXT_PARSER = mathtext.MathTextParser("agg")
_LATEX_CACHE: Dict[str, Optional[QPixmap]] = {}


def _latex_to_pixmap(latex: str) -> Optional[QPixmap]:
    if not latex:
        return None
    cached = _LATEX_CACHE.get(latex)
    if cached is not None:
        return cached
    try:
        ftimage, _ = _MATH_TEXT_PARSER.to_rgba(f"${latex}$", dpi=180)
    except Exception as exc:
        print(f"[warn] Failed to render LaTeX '{latex}': {exc}")
        _LATEX_CACHE[latex] = None
        return None

    buffer = np.ascontiguousarray(np.clip(ftimage * 255.0, 0, 255).astype(np.uint8))
    height, width, _channels = buffer.shape
    image = QImage(buffer.data, width, height, 4 * width, QImage.Format_RGBA8888)
    pixmap = QPixmap.fromImage(image.copy())
    _LATEX_CACHE[latex] = pixmap
    return pixmap


def _mass_identifier_from_path(path: str) -> Optional[str]:
    if not path or "/Masses/" not in path:
        return None
    tail = path.split("/Masses/", 1)[-1]
    segment = tail.split("/", 1)[0]
    for suffix in ("FitParameters", "FitParams", "Parameters", "Params"):
        if segment.endswith(suffix):
            segment = segment[: -len(suffix)]
            break
    segment = segment.strip().strip("_-")
    return segment or None


def _dataset_sort_key(channel: str, dataset_path: str) -> Tuple[int, object]:
    mass_label = _mass_identifier_from_path(dataset_path)
    if mass_label is not None:
        try:
            mass_value: object = float(mass_label)
        except ValueError:
            mass_value = mass_label
        return (0, mass_value)
    return (1, os.path.basename(dataset_path))


def _coerce_parameter_values(values: np.ndarray) -> Optional[np.ndarray]:
    try:
        arr = np.asarray(values)
    except Exception:
        return None

    if arr.size == 0:
        return None

    if arr.dtype.names:  # structured array
        for field in ("value", "values", "val"):
            if field in arr.dtype.names:
                arr = np.asarray(arr[field])
                break
        else:
            pieces = []
            for name in arr.dtype.names:
                try:
                    pieces.append(np.asarray(arr[name]).ravel())
                except Exception:
                    return None
            if not pieces:
                return None
            arr = np.concatenate(pieces)
    elif arr.dtype == object:
        coerced: List[float] = []
        for item in arr.ravel():
            if item is None:
                coerced.append(np.nan)
                continue
            if isinstance(item, (float, int, np.floating, np.integer)):
                coerced.append(float(item))
                continue
            if hasattr(item, "value"):
                try:
                    coerced.append(float(item.value))
                    continue
                except Exception:
                    pass
            if isinstance(item, (list, tuple, np.ndarray)):
                sub = np.asarray(item).ravel()
                if sub.size:
                    try:
                        coerced.append(float(sub[0]))
                        continue
                    except Exception:
                        pass
            try:
                coerced.append(float(item))
            except Exception:
                return None
        arr = np.asarray(coerced, dtype=float)
    else:
        arr = np.asarray(arr, dtype=float)

    arr = np.asarray(arr, dtype=float).ravel()
    if arr.size == 0:
        return None
    return arr

def _to_1d(array: np.ndarray) -> np.ndarray:
    arr = np.asarray(array)
    if arr.ndim == 0:
        return arr.reshape(1)
    return arr.ravel()

# --------- Channel & fit metadata ---------
FAMILY_HIGH = "high"
FAMILY_LOW = "low"


@dataclass(frozen=True)
class ChannelDefinition:
    dataset: str
    time_dataset: str
    family: str


CHANNEL_DEFS: Dict[str, ChannelDefinition] = {
    "TOF L": ChannelDefinition(dataset="TOF L", time_dataset="Time (high sampling)", family=FAMILY_HIGH),
    "TOF M": ChannelDefinition(dataset="TOF M", time_dataset="Time (high sampling)", family=FAMILY_HIGH),
    "TOF H": ChannelDefinition(dataset="TOF H", time_dataset="Time (high sampling)", family=FAMILY_HIGH),
    "Ion Grid": ChannelDefinition(dataset="Ion Grid", time_dataset="Time (low sampling)", family=FAMILY_LOW),
    "Target L": ChannelDefinition(dataset="Target L", time_dataset="Time (low sampling)", family=FAMILY_LOW),
    "Target H": ChannelDefinition(dataset="Target H", time_dataset="Time (low sampling)", family=FAMILY_LOW),
}

CHANNEL_ORDER: List[str] = ["TOF L", "TOF M", "TOF H", "Ion Grid", "Target L", "Target H"]

FIT_ELIGIBLE_CHANNELS = {"Ion Grid", "Target L", "Target H", "TOF H"}

BASELINE_PRIMARY_WINDOW = (-7.0, -5.0)
BASELINE_FALLBACK_THRESHOLD = -2.0

FAMILY_YLABELS = {
    FAMILY_HIGH: r"$TOF$ [pC/\Delta t]",
    FAMILY_LOW: r"$Q$ [pC]",
}

TIME_AXIS_LABEL = r"Time [$\mu$s]"


FIT_SCALE_MULTIPLIERS: Dict[str, float] = {
    "Ion Grid": 7.46e-4,
    "Target H": 1.63e-1,
    "Target L": 1.58e1,
}


@dataclass(frozen=True)
class FitModelMeta:
    parameter_labels: Tuple[str, ...]
    latex: str
    evaluator: Callable[..., np.ndarray]


ION_GRID_PARAMETER_LABELS_FULL: Tuple[str, ...] = (
    "t0 (impact onset)",
    "C0 (baseline)",
    "C1 (image amplitude)",
    "C2 (impact amplitude)",
    "τ0 (image decay)",
    "τ1 (impact rise)",
    "τ2 (impact decay)",
)

ION_GRID_PARAMETER_LABELS_LEGACY: Tuple[str, ...] = (
    "t0 (impact onset)",
    "C0 (baseline)",
    "C2 (impact amplitude)",
    "τ1 (impact rise)",
    "τ2 (impact decay)",
)

ION_GRID_LABEL_OVERRIDES: Dict[int, Tuple[str, ...]] = {
    len(ION_GRID_PARAMETER_LABELS_LEGACY): ION_GRID_PARAMETER_LABELS_LEGACY,
}


def _idex_ion_grid_model(x: np.ndarray, *params: float) -> np.ndarray:
    """Evaluate the standard IDEX ion grid / target response model."""

    arr = np.asarray(x, dtype=float)
    if arr.size == 0:
        return arr

    if len(params) >= 7:
        t0, c0, c1, c2, tau0, tau1, tau2 = params[:7]
        shift = arr - t0
        step = np.heaviside(shift, 0.0)
        safe_tau0 = tau0 if abs(tau0) > 1.0e-12 else 1.0e-12
        safe_tau1 = tau1 if abs(tau1) > 1.0e-12 else 1.0e-12
        safe_tau2 = tau2 if abs(tau2) > 1.0e-12 else 1.0e-12
        with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
            image_decay = np.exp(-shift / safe_tau0)
            rise = 1.0 - np.exp(-shift / safe_tau1)
            decay = np.exp(-shift / safe_tau2)
        return c0 + step * (-c1 * image_decay + c2 * rise * decay - c1)

    if len(params) >= 5:
        p0, p1, p4, p5, p6 = params[:5]
        shift = arr - p0
        step = np.heaviside(shift, 0.0)
        safe_p5 = p5 if abs(p5) > 1.0e-12 else 1.0e-12
        safe_p6 = p6 if abs(p6) > 1.0e-12 else 1.0e-12
        with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
            rise = 1.0 - np.exp(-shift / safe_p5)
            decay = np.exp(-shift / safe_p6)
        return p1 + step * (p4 * rise * decay)

    if len(params) >= 1:
        baseline = params[1] if len(params) > 1 else params[0]
        return np.full_like(arr, float(baseline), dtype=float)

    return np.zeros_like(arr, dtype=float)


def _emg_model(x: np.ndarray, mu: float, sigma: float, lam: float) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.size == 0:
        return arr
    if abs(lam) < 1.0e-15:
        return np.zeros_like(arr, dtype=float)
    safe_sigma = sigma if abs(sigma) > 1.0e-15 else 1.0e-15
    safe_lambda = lam
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        exponent = np.exp((safe_lambda / 2.0) * (2.0 * mu + safe_lambda * safe_sigma**2 - 2.0 * arr))
    argument = (mu + safe_lambda * safe_sigma**2 - arr) / (np.sqrt(2.0) * safe_sigma)
    return (safe_lambda / 2.0) * exponent * np.erfc(argument)


ION_GRID_FIT = FitModelMeta(
    parameter_labels=ION_GRID_PARAMETER_LABELS_FULL,
    latex=(
        r"w(t; t_0, C_0, C_1, C_2, \tau_0, \tau_1, \tau_2) = C_0"
        r" - H(t - t_0)\,C_1\,e^{-(t - t_0)/\tau_0}"
        r" + H(t - t_0)\left[C_2\left(1 - e^{-(t - t_0)/\tau_1}\right)e^{-(t - t_0)/\tau_2} - C_1\right]"
    ),
    evaluator=_idex_ion_grid_model,
)

EMG_FIT = FitModelMeta(
    parameter_labels=(
        "μ (location)",
        "σ (width)",
        "λ (decay)",
    ),
    latex=(
        r"f(t) = \frac{\lambda}{2} \exp\left[\frac{\lambda}{2}(2\mu + \lambda\sigma^2 - 2t)\right]"
        r" \mathrm{erfc}\left(\frac{\mu + \lambda\sigma^2 - t}{\sqrt{2}\,\sigma}\right)"
    ),
    evaluator=_emg_model,
)

FIT_MODEL_BY_CHANNEL: Dict[str, FitModelMeta] = {
    "Ion Grid": ION_GRID_FIT,
    "Target L": ION_GRID_FIT,
    "Target H": ION_GRID_FIT,
    "TOF H": EMG_FIT,
}


def _evaluate_fit_curve(channel: str, time_values: np.ndarray, params: np.ndarray) -> Optional[np.ndarray]:
    """Evaluate a fit curve for the supplied channel using the stored parameters."""

    model = FIT_MODEL_BY_CHANNEL.get(channel)
    if not model:
        return None

    param_array = _coerce_parameter_values(params)
    if param_array is None or param_array.size < len(model.parameter_labels):
        return None

    first = np.asarray(param_array[: len(model.parameter_labels)], dtype=float)
    if not np.all(np.isfinite(first)):
        return None

    try:
        return model.evaluator(time_values, *first)
    except Exception:
        return None


def _masked_mean(values: np.ndarray, mask: np.ndarray) -> Optional[float]:
    values = np.asarray(values)
    mask = np.asarray(mask, dtype=bool)
    if values.size == 0 or mask.size != values.size:
        return None
    if not np.any(mask):
        return None
    subset = values[mask]
    subset = subset[np.isfinite(subset)]
    if subset.size == 0:
        return None
    return float(subset.mean())


def _fit_paths_from_param(dataset_path: str) -> Optional[Tuple[str, str]]:
    """Return the fit result/time dataset paths corresponding to a parameter dataset."""

    if not dataset_path:
        return None

    for needle in ("FitParams", "FitParameters", "FitParam"):
        if needle in dataset_path:
            result_path = dataset_path.replace(needle, "FitResult")
            time_path = dataset_path.replace(needle, "FitTime")
            return result_path, time_path

    return None


@dataclass
class FitData:
    time_series: Dict[str, np.ndarray] = field(default_factory=dict)
    value_series: Dict[str, np.ndarray] = field(default_factory=dict)
    parameter_series: Dict[str, np.ndarray] = field(default_factory=dict)
    extras: Dict[str, np.ndarray] = field(default_factory=dict)

    def has_overlay(self) -> bool:
        if any(True for _ in self.iter_time_result_pairs()):
            return True
        return bool(self.parameter_series)

    def has_parameters(self) -> bool:
        return bool(self.parameter_series)

    def iter_time_result_pairs(self) -> Iterable[Tuple[str, str, str, np.ndarray, np.ndarray]]:
        seen_keys = set()
        for time_path, time_values in self.time_series.items():
            pair_key = _pair_key(time_path)
            if not pair_key:
                continue
            for value_path, value_values in self.value_series.items():
                if _pair_key(value_path) != pair_key:
                    continue
                dedupe_key = (time_path, value_path)
                if dedupe_key in seen_keys:
                    continue
                seen_keys.add(dedupe_key)
                yield (
                    _friendly_label(value_path),
                    time_path,
                    value_path,
                    _to_1d(time_values),
                    _to_1d(value_values),
                )

    def iter_parameter_items(self) -> Iterable[Tuple[str, np.ndarray]]:
        for path, values in self.parameter_series.items():
            yield path, np.asarray(values)


class BaseDataSource(ABC):
    """Abstract data provider used by the quicklook controller."""

    def __init__(self, filename: str):
        self.filename = filename

    # ------------------------------------------------------------------
    # Lifecycle helpers
    # ------------------------------------------------------------------
    def close(self) -> None:
        """Close any open resources (default: no-op)."""

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------
    @abstractmethod
    def list_events(self) -> List[str]:
        """Return the ordered list of events contained in the file."""

    @abstractmethod
    def get_dataset(self, event: str, dataset_name: str) -> Optional[np.ndarray]:
        """Return the dataset for ``/{event}/{dataset_name}`` or ``None``."""

    @abstractmethod
    def gather_fit_data(self, event: str, channel: str) -> FitData:
        """Collect fit products for ``channel`` within ``event``."""

    def describe(self) -> str:
        return os.path.basename(self.filename)

    def supports_structure_viewer(self) -> bool:
        return False

    def launch_structure_viewer(self, parent: QWidget | None = None) -> QWidget:
        raise NotImplementedError("Structure viewer not implemented for this source")


class HDF5DataSource(BaseDataSource):
    """Adapter that exposes an HDF5 file via the quicklook datasource API."""

    def __init__(self, filename: str):
        if h5py is None:  # pragma: no cover - handled in runtime environment
            raise ImportError("h5py is required to open HDF5 files but is not installed.")
        super().__init__(filename)
        self._file = h5py.File(filename, "r")

    def close(self) -> None:
        try:
            self._file.close()
        except Exception:
            pass

    def list_events(self) -> List[str]:
        events: List[str] = []
        for key, value in self._file.items():
            if isinstance(value, h5py.Group):
                events.append(str(key))
        try:
            return sorted(events, key=lambda token: int(str(token)))
        except Exception:
            return sorted(events)

    def get_dataset(self, event: str, dataset_name: str) -> Optional[np.ndarray]:
        path = f"/{event}/{dataset_name}"
        try:
            obj = self._file[path]
        except Exception:
            return None
        if isinstance(obj, h5py.Dataset):
            try:
                return np.array(obj[()], copy=True)
            except Exception:
                return None
        return None

    def gather_fit_data(self, event: str, channel: str) -> FitData:
        data = FitData()
        group = self._file.get(event)
        if not isinstance(group, h5py.Group):
            return data

        base_norm = _normalize_name(channel)

        def visitor(name: str, obj):
            if not isinstance(obj, h5py.Dataset):
                return
            dataset_name = name.split("/")[-1]
            norm = _normalize_name(dataset_name)
            if "fit" not in norm or base_norm not in norm:
                return
            try:
                arr = np.array(obj[()], copy=True)
            except Exception:
                return
            full_path = f"{group.name}/{name}"
            if "time" in norm:
                data.time_series[full_path] = arr
            elif "result" in norm or "curve" in norm:
                data.value_series[full_path] = arr
            elif "param" in norm:
                data.parameter_series[full_path] = arr
            else:
                data.extras[full_path] = arr

        group.visititems(visitor)
        return data

    def supports_structure_viewer(self) -> bool:
        return True

    def launch_structure_viewer(self, parent: QWidget | None = None) -> QWidget:
        if launch_hdf_viewer is None:
            raise RuntimeError("HDF viewer is unavailable in this environment.")
        return launch_hdf_viewer(self.filename, parent=parent)


class CDFDataSource(BaseDataSource):
    """Adapter around ``cdflib`` to serve CDF content to the viewer."""

    #: Mapping from quicklook dataset names to CDF variable identifiers.
    DATASET_MAP: Dict[str, str] = {
        "Time (high sampling)": "time_high_sample_rate",
        "Time (low sampling)": "time_low_sample_rate",
        "TOF L": "tof_low",
        "TOF M": "tof_mid",
        "TOF H": "tof_high",
        "Ion Grid": "ion_grid",
        "Target L": "target_low",
        "Target H": "target_high",
    }

    #: Fit metadata required for quicklook overlays.
    FIT_VARIABLES: Dict[str, Dict[str, str]] = {
        "Ion Grid": {
            "result": "ion_grid_fit_results",
            "params": "ion_grid_fit_parameters",
        },
        "Target L": {
            "result": "target_low_fit_results",
            "params": "target_low_fit_parameters",
        },
        "Target H": {
            "result": "target_high_fit_results",
            "params": "target_high_fit_parameters",
        },
        "TOF H": {
            "params": "tof_peak_fit_parameters",
        },
    }

    def __init__(self, filename: str):
        try:
            import cdflib  # type: ignore
        except Exception as exc:  # pragma: no cover - dependency not installed
            raise ImportError(
                "cdflib is required to open CDF files but is not installed."
            ) from exc

        super().__init__(filename)
        self._cdflib = cdflib
        self._cdf = cdflib.CDF(filename)
        self._cache: Dict[str, np.ndarray] = {}
        self._event_count = self._resolve_event_count()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def close(self) -> None:
        try:
            self._cdf.close()
        except Exception:
            pass

    def _resolve_event_count(self) -> int:
        try:
            epoch = self._cdf.varget("epoch")
            return int(np.asarray(epoch).shape[0])
        except Exception:
            # Fallback: inspect a known waveform variable.
            for varname in ("tof_low", "tof_mid", "ion_grid"):
                try:
                    candidate = self._cdf.varget(varname)
                    return int(np.asarray(candidate).shape[0])
                except Exception:
                    continue
        return 0

    def _event_index(self, event: str) -> int:
        try:
            return max(0, int(str(event)) - 1)
        except Exception:
            return 0

    @lru_cache(maxsize=None)
    def _cached_variable(self, varname: str) -> np.ndarray:
        try:
            return np.asarray(self._cdf.varget(varname))
        except Exception:
            return np.asarray([])

    def list_events(self) -> List[str]:
        return [str(idx + 1) for idx in range(self._event_count)]

    def get_dataset(self, event: str, dataset_name: str) -> Optional[np.ndarray]:
        varname = self.DATASET_MAP.get(dataset_name)
        if not varname:
            return None
        data = self._cached_variable(varname)
        if data.size == 0:
            return None
        index = self._event_index(event)
        if index >= data.shape[0]:
            return None
        try:
            slice_data = data[index]
        except Exception:
            return None
        return np.array(slice_data, copy=True)

    def gather_fit_data(self, event: str, channel: str) -> FitData:
        data = FitData()
        mapping = self.FIT_VARIABLES.get(channel)
        if not mapping:
            return data

        index = self._event_index(event)
        event_prefix = f"/{event}/Analysis"

        # Time bases reuse the shared sampling axes.
        if channel in ("Ion Grid", "Target L", "Target H"):
            time_name = "Time (low sampling)"
        else:
            time_name = "Time (high sampling)"

        base_time = self.get_dataset(event, time_name)
        if base_time is not None:
            time_key = f"{event_prefix}/{channel} Fit Time"
            data.time_series[time_key] = np.array(base_time, copy=True)

        result_var = mapping.get("result")
        if result_var:
            result_data = self._cached_variable(result_var)
            if result_data.size and index < result_data.shape[0]:
                result_key = f"{event_prefix}/{channel} Fit Result"
                data.value_series[result_key] = np.array(result_data[index], copy=True)

        params_var = mapping.get("params")
        if params_var:
            params_data = self._cached_variable(params_var)
            if params_data.size and index < params_data.shape[0]:
                params_key = f"{event_prefix}/{channel} Fit Parameters"
                data.parameter_series[params_key] = np.array(params_data[index], copy=True)

        return data

    def supports_structure_viewer(self) -> bool:
        return launch_cdf_viewer is not None

    def launch_structure_viewer(self, parent: QWidget | None = None) -> QWidget:
        if launch_cdf_viewer is None:
            raise RuntimeError("CDF viewer is not available. Install optional GUI components.")
        return launch_cdf_viewer(self.filename, parent=parent)


def create_data_source(filename: str) -> BaseDataSource:
    """Instantiate the appropriate data source for ``filename``."""

    suffix = Path(filename).suffix.lower()
    if suffix in {".h5", ".hdf5"}:
        return HDF5DataSource(filename)
    if suffix == ".cdf":
        return CDFDataSource(filename)

    # Fallback: attempt HDF5 first, then CDF.
    try:
        return HDF5DataSource(filename)
    except Exception:
        return CDFDataSource(filename)


Y_AXIS_LABELS: Dict[str, str] = {
    "Target L": r"$Q_{TL}$ [pC]",
    "Target H": r"$Q_{TH}$ [pC]",
    "Ion Grid": r"$Q_{IG}$ [pC]",
    "TOF L": r"$TOF_{L}$ [pC/\Delta t]",
    "TOF M": r"$TOF_{M}$ [pC/\Delta t]",
    "TOF H": r"$TOF_{H}$ [pC/\Delta t]",
}


def y_label_with_units(channel_name: str) -> str:
    return Y_AXIS_LABELS.get(channel_name, channel_name)

# --------- Main Window ---------
class MainWindow(QMainWindow):
    def __init__(self, filename: str | None = None, eventnumber: int | None = None):
        super().__init__()
        self.setWindowTitle("IDEX Quicklook — Interactive Viewer")
        self.setMinimumSize(1250, 820)
        self.setStatusBar(QStatusBar(self))

        logo_path = _find_brand_logo()
        if logo_path is not None:
            self.setWindowIcon(QIcon(str(logo_path)))

        self._data_source: Optional[BaseDataSource] = None
        self._events: List[str] = []
        self._current_event: Optional[str] = None
        self._filename: Optional[str] = None
        self._h5: Optional[Any] = None
        self._tmpdir = tempfile.TemporaryDirectory(prefix="idex_quicklook_")
        self._fit_cache: Dict[Tuple[str, str], FitData] = {}
        self._fit_param_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._fit_result_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._fit_time_overrides: Dict[Tuple[str, str, str], np.ndarray] = {}
        self._baseline_cache: Dict[Tuple[str, str], float] = {}
        self._show_fit: Dict[str, bool] = {name: False for name in FIT_ELIGIBLE_CHANNELS}
        # ``refresh_fit_controls`` is invoked as soon as data are loaded, so make
        # sure ``fit_buttons`` always exists even before the control panel is
        # constructed.  Some startup paths (for example early failures while
        # building the control panel) were leaving the attribute undefined which
        # caused an ``AttributeError`` when plotting the first event.
        self.fit_buttons: Dict[str, QPushButton] = {}
        self.selected_channels = set(CHANNEL_ORDER)
        self._child_windows: List[QWidget] = []
        self._documentation_center: Optional[DocumentationCenter] = None

        self._create_actions()

        self._apply_modern_palette()

        central = QWidget(self)
        self.vbox = QVBoxLayout(central)
        self.vbox.setContentsMargins(10, 10, 10, 10)
        self.vbox.setSpacing(12)
        self.setCentralWidget(central)

        branding_widget = self._build_branding_banner()
        if branding_widget is not None:
            self.vbox.addWidget(branding_widget)

        self._build_menubar()
        self._build_toolbar()
        self._build_controls()

        self.figure = Figure(figsize=(12, 8), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.nav_toolbar = NavigationToolbar(self.canvas, self)
        self.nav_toolbar.setStyleSheet("font-size: 14px; padding: 6px;")
        self.vbox.addWidget(self.nav_toolbar)
        self.vbox.addWidget(self.canvas)

        self.statusBar().showMessage("Select a data file (HDF5 or CDF) to get started.")

        if filename:
            self.open_file(filename)
        else:
            chosen = prompt_for_data_file(self)
            if not chosen:
                self.close()
                return
            self.open_file(chosen)

        if eventnumber is not None and self._events and 1 <= eventnumber <= len(self._events):
            self.event_combo.setCurrentIndex(eventnumber - 1)

    # ---- UI construction -------------------------------------------------
    def _apply_modern_palette(self) -> None:
        self.setStyleSheet(
            """
            QMainWindow {
                background-color: #edf1fb;
            }
            QMenuBar {
                background-color: #ffffff;
                font-size: 14px;
            }
            QMenu {
                font-size: 14px;
                padding: 6px 10px;
            }
            QToolBar {
                background-color: #ffffff;
                padding: 6px;
                border: none;
            }
            QStatusBar {
                background-color: #f4f6fb;
                font-size: 13px;
            }
            QListWidget {
                font-size: 14px;
            }
            QTextBrowser {
                font-size: 14px;
            }
            """
        )
        status = self.statusBar()
        status.setStyleSheet("background-color: #f4f6fb; font-size: 13px; padding: 4px 8px;")
        status.setContentsMargins(6, 0, 6, 0)

    def _create_actions(self):
        self.open_any_action = QAction("Open…", self)
        self.open_any_action.setShortcut("Ctrl+O")
        self.open_any_action.triggered.connect(self.action_open)

        self.open_hdf5_action = QAction("Open HDF5…", self)
        self.open_hdf5_action.triggered.connect(lambda: self.action_open(preferred="hdf5"))

        self.open_cdf_action = QAction("Open CDF…", self)
        self.open_cdf_action.triggered.connect(lambda: self.action_open(preferred="cdf"))

        self.view_structure_action = QAction("Open Data Browser", self)
        self.view_structure_action.setShortcut("Ctrl+B")
        self.view_structure_action.triggered.connect(self.action_view_structure)

        self.open_variable_definitions_action = QAction("Variable Definitions…", self)
        self.open_variable_definitions_action.setEnabled(launch_variable_definitions_viewer is not None)
        self.open_variable_definitions_action.triggered.connect(self.action_open_variable_definitions)

        self.reload_action = QAction("Reload", self)
        self.reload_action.setShortcut("Ctrl+R")
        self.reload_action.triggered.connect(self.reload_current)

        self.close_file_action = QAction("Close File", self)
        self.close_file_action.triggered.connect(self.close_current_file)

        self.quit_action = QAction("Quit", self)
        self.quit_action.setShortcut("Ctrl+Q")
        self.quit_action.triggered.connect(self.close)

        self.save_png_action = QAction("Save Plot as PNG…", self)
        self.save_png_action.triggered.connect(lambda: self.save_plot("png"))

        self.save_pdf_action = QAction("Save Plot as PDF…", self)
        self.save_pdf_action.triggered.connect(lambda: self.save_plot("pdf"))

        self.save_svg_action = QAction("Save Plot as SVG…", self)
        self.save_svg_action.triggered.connect(lambda: self.save_plot("svg"))

        self.edit_fit_action = QAction("Edit Fit Parameters", self)
        self.edit_fit_action.setShortcut("Ctrl+E")
        self.edit_fit_action.triggered.connect(self.open_fit_parameter_dialog)

        self.reset_fit_action = QAction("Reset Fit Overrides", self)
        self.reset_fit_action.triggered.connect(self.reset_all_overrides)

        help_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_DialogHelpButton)
        self.help_action = QAction(help_icon, "Documentation Center", self)
        self.help_action.setShortcut("F1")
        self.help_action.setToolTip("Open the documentation center (F1)")
        self.help_action.setStatusTip("Open the searchable documentation center")
        self.help_action.setIconText("?")
        self.help_action.triggered.connect(self.open_documentation_center)

        self.search_docs_action = QAction("Search Documentation…", self)
        self.search_docs_action.setShortcut("Ctrl+F1")
        self.search_docs_action.setStatusTip("Jump straight to the documentation search panel")
        self.search_docs_action.triggered.connect(lambda: self.open_documentation_center(show_search=True))

    def _build_menubar(self):
        menubar = self.menuBar()
        if menubar is None:
            menubar = QMenuBar(self)
            self.setMenuBar(menubar)

        menubar.clear()
        menubar.setStyleSheet("font-size: 14px; background-color: #ffffff; padding: 4px;")

        file_menu = menubar.addMenu("&File")
        file_menu.addAction(self.open_any_action)
        file_menu.addAction(self.open_hdf5_action)
        file_menu.addAction(self.open_cdf_action)
        file_menu.addSeparator()

        save_menu = QMenu("Export Plot", self)
        save_menu.addAction(self.save_png_action)
        save_menu.addAction(self.save_pdf_action)
        save_menu.addAction(self.save_svg_action)
        file_menu.addMenu(save_menu)

        file_menu.addSeparator()
        file_menu.addAction(self.reload_action)
        file_menu.addAction(self.close_file_action)
        file_menu.addSeparator()
        file_menu.addAction(self.quit_action)

        edit_menu = menubar.addMenu("&Edit")
        edit_menu.addAction(self.edit_fit_action)
        edit_menu.addAction(self.reset_fit_action)

        view_menu = menubar.addMenu("&View")
        view_menu.addAction(self.view_structure_action)
        view_menu.addAction(self.open_variable_definitions_action)

        help_menu = menubar.addMenu("&Help")
        help_menu.addAction(self.help_action)
        help_menu.addAction(self.search_docs_action)
        help_menu.addSeparator()

        self.menu_search_field = QLineEdit(self)
        self.menu_search_field.setPlaceholderText("Search documentation…")
        self.menu_search_field.setClearButtonEnabled(True)
        self.menu_search_field.returnPressed.connect(self._trigger_menu_search)

        search_widget = QWidgetAction(self)
        search_widget.setDefaultWidget(self.menu_search_field)
        help_menu.addAction(search_widget)

    def _trigger_menu_search(self) -> None:
        query = self.menu_search_field.text().strip()
        self.open_documentation_center(initial_query=query, show_search=True)

    def _build_toolbar(self):
        tb = QToolBar("Main", self)
        tb.setIconSize(QSize(22, 22))
        tb.setMovable(False)
        tb.setStyleSheet("QToolBar { background-color: #ffffff; border: none; padding: 8px; spacing: 8px; }")
        self.addToolBar(tb)

        tb.addAction(self.open_any_action)
        tb.addAction(self.open_hdf5_action)
        tb.addAction(self.open_cdf_action)
        tb.addSeparator()
        tb.addAction(self.reload_action)

        export_button = QToolButton(self)
        export_button.setText("Export Plot")
        export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        export_menu = QMenu(export_button)
        export_menu.addAction(self.save_png_action)
        export_menu.addAction(self.save_pdf_action)
        export_menu.addAction(self.save_svg_action)
        export_button.setMenu(export_menu)
        export_button.setStyleSheet(
            """
            QToolButton {
                font-size: 15px;
                font-weight: 600;
                padding: 6px 14px;
                border-radius: 12px;
                background-color: #4263eb;
                color: #ffffff;
            }
            QToolButton:hover {
                background-color: #3b5bdb;
            }
            """
        )
        tb.addWidget(export_button)

        tb.addSeparator()

        act_dust = QAction("Dust Composition…", self)
        act_dust.setShortcut("Ctrl+D")
        act_dust.setToolTip("Open the dust composition analysis window for the current event.")
        act_dust.triggered.connect(self.action_open_dust_composition)
        self.addAction(act_dust)

        self.dust_button = QPushButton("Dust Composition", self)
        self.dust_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.dust_button.setMinimumHeight(46)
        self.dust_button.setStyleSheet(
            """
            QPushButton {
                font-size: 16px;
                font-weight: 700;
                padding: 10px 22px;
                border-radius: 14px;
                background-color: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                                 stop:0 #845ef7, stop:1 #5c7cfa);
                color: #ffffff;
                border: 1px solid #5f3dc4;
            }
            QPushButton:hover {
                background-color: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                                 stop:0 #7048e8, stop:1 #4c6ef5);
            }
            QPushButton:pressed {
                background-color: #364fc7;
            }
            """
        )
        self.dust_button.clicked.connect(act_dust.trigger)
        tb.addWidget(self.dust_button)

        act_reload = QAction("Reload", self)
        act_reload.setShortcut("Ctrl+R")
        act_reload.triggered.connect(self.reload_current)
        tb.addAction(act_reload)
        tb.addAction(self.view_structure_action)
        if launch_variable_definitions_viewer is not None:
            tb.addAction(self.open_variable_definitions_action)
        tb.addSeparator()

        tb.addAction(self.close_file_action)
        tb.addSeparator()

        tb.addAction(self.quit_action)
        tb.addSeparator()

        tb.addAction(self.help_action)
        tb.addSeparator()
        label = QLabel("Event:", self)
        label.setStyleSheet("font-size: 15px; font-weight: bold; padding-right: 6px;")
        tb.addWidget(label)

        self.event_combo = QComboBox(self)
        self.event_combo.setMinimumWidth(220)
        self.event_combo.setStyleSheet("font-size: 15px; min-height: 36px;")
        self.event_combo.currentIndexChanged.connect(self.on_event_changed)
        tb.addWidget(self.event_combo)

    def _build_branding_banner(self) -> Optional[QWidget]:
        pixmap = _load_brand_pixmap(max_height=64)
        if pixmap is None:
            return None

        container = QWidget(self)
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)

        logo_label = QLabel()
        logo_label.setPixmap(pixmap)
        logo_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        layout.addWidget(logo_label)

        title_label = QLabel("SpectrumPY: Flight Addition")
        title_label.setStyleSheet("font-size: 18px; font-weight: 600; color: #1f2937;")
        title_label.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(title_label)

        layout.addStretch()
        return container

    def _build_controls(self):
        panel = QWidget(self)
        panel.setObjectName("controlPanel")
        panel.setStyleSheet(
            """
            QWidget#controlPanel {
                background-color: #ffffff;
                border-radius: 16px;
                border: 1px solid #d6dfee;
                padding: 16px;
            }
            """
        )
        panel_layout = QVBoxLayout(panel)
        panel_layout.setContentsMargins(0, 0, 0, 0)
        panel_layout.setSpacing(8)
        panel_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        heading = QLabel("Display Controls", self)
        heading.setObjectName("controlHeading")
        heading.setStyleSheet("font-size: 19px; font-weight: 600;")
        panel_layout.addWidget(heading)

        sub_label = QLabel("Choose which channels and overlays are shown:", self)
        sub_label.setStyleSheet("font-size: 15px; color: #4a5568;")
        panel_layout.addWidget(sub_label)

        channel_widget = QWidget(self)
        grid = QGridLayout(channel_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(10)
        self.channel_buttons: Dict[str, QPushButton] = {}

        toggle_style = (
            """
            QPushButton {
                font-size: 16px;
                font-weight: 600;
                padding: 10px 18px;
                border-radius: 12px;
                background-color: #e8f0ff;
                border: 1px solid #c3d0ff;
            }
            QPushButton:hover {
                background-color: #dbe4ff;
            }
            QPushButton:checked {
                background-color: #4c6ef5;
                color: #ffffff;
                border: 1px solid #364fc7;
            }
            """
        )
        primary_style = (
            """
            QPushButton {
                font-size: 16px;
                font-weight: 600;
                padding: 10px 18px;
                border-radius: 12px;
                background-color: #4263eb;
                color: #ffffff;
            }
            QPushButton:hover {
                background-color: #3b5bdb;
            }
            """
        )

        self._primary_channel_buttons: List[QPushButton] = []
        for idx, name in enumerate(CHANNEL_ORDER):
            btn = QPushButton(name, self)
            btn.setCheckable(True)
            btn.setChecked(True)
            btn.setMinimumHeight(50)
            btn.setStyleSheet(toggle_style)
            btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            btn.clicked.connect(lambda checked, channel=name: self.on_channel_toggled(channel, checked))
            self.channel_buttons[name] = btn
            grid.addWidget(btn, idx // 3, idx % 3)
            self._primary_channel_buttons.append(btn)

        panel_layout.addWidget(channel_widget)

        toggle_row = QHBoxLayout()
        toggle_row.setSpacing(10)

        self.overlay_button = QPushButton("Overlay same time axis", self)
        self.overlay_button.setCheckable(True)
        self.overlay_button.setMinimumHeight(50)
        self.overlay_button.setStyleSheet(toggle_style)
        self.overlay_button.setToolTip("When enabled, channels with the same time base are drawn together.")
        self.overlay_button.clicked.connect(self.on_overlay_toggled)
        toggle_row.addWidget(self.overlay_button)

        self.fit_buttons: Dict[str, QPushButton] = {}
        for channel in sorted(FIT_ELIGIBLE_CHANNELS):
            btn = QPushButton(f"Show {channel} Fit", self)
            btn.setCheckable(True)
            btn.setMinimumHeight(50)
            btn.setStyleSheet(toggle_style)
            btn.setToolTip("Overlay fit curves when available.")
            btn.clicked.connect(lambda checked, chan=channel: self.on_fit_toggled(chan, checked))
            toggle_row.addWidget(btn)
            self.fit_buttons[channel] = btn

        self.edit_params_button = QPushButton("Edit Fit Parameters", self)
        self.edit_params_button.setMinimumHeight(50)
        self.edit_params_button.setStyleSheet(primary_style)
        self.edit_params_button.clicked.connect(self.open_fit_parameter_dialog)
        toggle_row.addWidget(self.edit_params_button)

        toggle_row.addStretch(1)
        panel_layout.addLayout(toggle_row)

        self._harmonize_primary_button_widths()
        QTimer.singleShot(0, self._harmonize_primary_button_widths)

        self.vbox.addWidget(panel)

    def _harmonize_primary_button_widths(self) -> None:
        """Ensure channel toggles share the width of the primary action button."""

        if not getattr(self, "_primary_channel_buttons", None):
            return

        reference_widgets = [self.edit_params_button, self.overlay_button, *self.fit_buttons.values()]
        reference_width = 0
        for widget in reference_widgets:
            if widget is None:
                continue
            reference_width = max(reference_width, widget.sizeHint().width())

        if reference_width <= 0:
            return

        reference_height = 0
        for widget in reference_widgets:
            if widget is None:
                continue
            reference_height = max(reference_height, widget.sizeHint().height())

        for btn in self._primary_channel_buttons:
            btn.setMinimumWidth(reference_width)
            if reference_height > 0:
                btn.setMinimumHeight(max(btn.minimumHeight(), reference_height))

    def _reset_layout_engine(self) -> None:
        """Re-enable Matplotlib's constrained layout after clearing the figure."""

        try:
            self.figure.set_constrained_layout(True)
        except Exception:
            try:
                self.figure.set_layout_engine("constrained")
            except Exception:
                pass

    def open_documentation_center(
        self,
        initial_query: Optional[str] = None,
        *,
        show_search: bool = False,
    ) -> None:
        query = initial_query or ""
        if self._documentation_center is None:
            self._documentation_center = DocumentationCenter(self, initial_query=query)
            self._documentation_center.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

            def _clear_reference(_object=None) -> None:
                self._documentation_center = None

            self._documentation_center.destroyed.connect(_clear_reference)
        else:
            if query:
                self._documentation_center.set_query(query)
            elif not self._documentation_center.isVisible():
                self._documentation_center.set_query("")

        if self._documentation_center is None:
            return

        self._documentation_center.show()
        self._documentation_center.raise_()
        self._documentation_center.activateWindow()

        if show_search:
            self._documentation_center.focus_search()

    # ---- Actions ---------------------------------------------------------
    def action_open(self, preferred: Optional[str] = None):
        chosen = prompt_for_data_file(self, preferred=preferred)
        if chosen:
            self.open_file(chosen)

    def action_view_structure(self):
        if not self._filename or not self._data_source:
            QMessageBox.information(
                self,
                "No File Loaded",
                "Open a data file to browse its structure.",
            )
            return

        if not self._data_source.supports_structure_viewer():
            QMessageBox.information(
                self,
                "Viewer Unavailable",
                "A structure viewer is not available for this data source.",
            )
            return

        try:
            viewer = self._data_source.launch_structure_viewer(parent=self)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Viewer Error",
                f"Unable to launch the data browser:\n{exc}",
            )
            return

        viewer.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        viewer.raise_()
        viewer.activateWindow()

        self._child_windows.append(viewer)

        def _cleanup(*_args):
            if viewer in self._child_windows:
                self._child_windows.remove(viewer)

        viewer.destroyed.connect(_cleanup)

    def action_open_variable_definitions(self):
        if launch_variable_definitions_viewer is None:
            QMessageBox.information(
                self,
                "Viewer Unavailable",
                "The variable definitions viewer is not available in this environment.",
            )
            return

        try:
            window = launch_variable_definitions_viewer(parent=self)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Viewer Error",
                f"Unable to launch the variable definitions viewer:\n{exc}",
            )
            return

        window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        window.raise_()
        window.activateWindow()

        self._child_windows.append(window)

        def _cleanup(*_args):
            if window in self._child_windows:
                self._child_windows.remove(window)

        window.destroyed.connect(_cleanup)

    def action_open_dust_composition(self):
        if not self._h5 or not self._current_event:
            QMessageBox.information(
                self,
                "No Event",
                "Open an HDF5 file and select an event before analysing dust composition.",
            )
            return

        try:
            window = launch_dust_composition_window(self._h5, self._current_event, parent=self)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Dust Composition Error",
                f"Unable to launch the dust composition window:\n{exc}",
            )
            return

        window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        window.show()

        self._child_windows.append(window)

        def _cleanup(*_args):
            if window in self._child_windows:
                self._child_windows.remove(window)

        window.destroyed.connect(_cleanup)

    def close_current_file(self):
        if not self._filename and not self._data_source:
            return
        self._reset_state()
        self.plot_event(None)
        self.statusBar().showMessage("No data file loaded.")

    def reset_all_overrides(self):
        if not (self._fit_param_overrides or self._fit_result_overrides or self._fit_time_overrides):
            QMessageBox.information(
                self,
                "No Overrides",
                "No in-memory fit overrides are currently applied.",
            )
            return

        self._fit_param_overrides.clear()
        self._fit_result_overrides.clear()
        self._fit_time_overrides.clear()
        self._fit_cache.clear()
        self._baseline_cache.clear()
        self.plot_event(self._current_event)
        self.statusBar().showMessage("Cleared all in-memory fit edits.", 6000)

    def save_plot(self, fmt: str):
        fmt = fmt.lower()
        if fmt not in {"png", "pdf", "svg"}:
            raise ValueError(f"Unsupported export format: {fmt}")

        if self.figure is None:
            return

        base_name = "idex_quicklook"
        if self._filename:
            base_name = Path(self._filename).stem
        if self._current_event:
            base_name += f"_event_{self._current_event}"

        suggested = Path(self._tmpdir.name) / f"{base_name}.{fmt}"

        path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export Plot ({fmt.upper()})",
            str(suggested),
            f"{fmt.upper()} files (*.{fmt});;All files (*)",
        )

        if not path:
            return

        export_path = Path(path)
        if export_path.suffix.lower() != f".{fmt}":
            export_path = export_path.with_suffix(f".{fmt}")

        try:
            self.figure.savefig(str(export_path), format=fmt, bbox_inches="tight")
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Export Failed",
                f"Could not export the plot:\n{exc}",
            )
            return

        self.statusBar().showMessage(f"Saved plot to {export_path}", 6000)

    # ---- Data helpers -----------------------------------------------------
    def _reset_state(self):
        if self._h5 is not None:
            try:
                self._h5.close()
            except Exception:
                pass
        self._h5 = None
        if self._data_source is not None:
            try:
                self._data_source.close()
            except Exception:
                pass
        self._data_source = None
        self._filename = None
        self._events = []
        self._current_event = None
        self._fit_cache.clear()
        self._fit_param_overrides.clear()
        self._fit_result_overrides.clear()
        self._fit_time_overrides.clear()
        self._baseline_cache.clear()
        self._show_fit = {name: False for name in FIT_ELIGIBLE_CHANNELS}

        self.event_combo.blockSignals(True)
        self.event_combo.clear()
        self.event_combo.blockSignals(False)

    def _get_dataset(self, event: str, dataset: str) -> Optional[np.ndarray]:
        if not self._data_source:
            return None
        return self._data_source.get_dataset(event, dataset)

    def _get_dataset_by_path(self, path: str) -> Optional[np.ndarray]:
        if not self._data_source:
            return None
        cleaned = path.strip("/")
        if not cleaned:
            return None
        parts = cleaned.split("/", 1)
        if len(parts) == 1:
            if not self._current_event:
                return None
            return self._data_source.get_dataset(self._current_event, parts[0])
        event, dataset = parts
        return self._data_source.get_dataset(event, dataset)

    def open_file(self, path: str, preferred_event: Optional[str] = None):
        self._reset_state()

        suffix = Path(path).suffix.lower()
        h5_handle = None

        try:
            if h5py is not None and suffix in {".h5", ".hdf5"}:
                try:
                    h5_handle = h5py.File(path, "r+")
                except OSError:
                    h5_handle = h5py.File(path, "r")
            source = create_data_source(path)
        except ImportError as exc:
            QMessageBox.critical(
                self,
                "Open Error",
                f"Failed to open file:\n{path}\n\n{exc}",
            )
            if h5_handle is not None:
                try:
                    h5_handle.close()
                except Exception:
                    pass
            self.plot_event(None)
            return
        except Exception as exc:
            if h5_handle is not None:
                try:
                    h5_handle.close()
                except Exception:
                    pass
            QMessageBox.critical(
                self,
                "Open Error",
                f"Failed to open file:\n{path}\n\n{exc}",
            )
            self.plot_event(None)
            return

        if h5_handle is None and h5py is not None and isinstance(source, HDF5DataSource):
            try:
                h5_handle = h5py.File(path, "r+")
            except OSError:
                try:
                    h5_handle = h5py.File(path, "r")
                except Exception:
                    h5_handle = None

        self._h5 = h5_handle
        self._data_source = source
        self._filename = path

        events = source.list_events()
        if not events:
            QMessageBox.warning(self, "No Events", "No top-level event groups found in this file.")
        self._events = events

        self.event_combo.blockSignals(True)
        self.event_combo.clear()
        self.event_combo.addItems(self._events)
        self.event_combo.blockSignals(False)

        target_event = None
        if preferred_event and preferred_event in self._events:
            target_event = preferred_event
        elif self._events:
            target_event = self._events[0]

        if target_event:
            idx = self._events.index(target_event)
            self.event_combo.blockSignals(True)
            self.event_combo.setCurrentIndex(idx)
            self.event_combo.blockSignals(False)
            self.plot_event(target_event)
        else:
            self.plot_event(None)

    def reload_current(self):
        if self._filename:
            self.open_file(self._filename, preferred_event=self._current_event)

    def on_event_changed(self, idx: int):
        if 0 <= idx < len(self._events):
            self.plot_event(self._events[idx])

    def on_channel_toggled(self, channel: str, checked: bool):
        if checked:
            self.selected_channels.add(channel)
        else:
            self.selected_channels.discard(channel)
        self.plot_event(self._current_event)

    def on_overlay_toggled(self, checked: bool):
        self.plot_event(self._current_event)

    def on_fit_toggled(self, channel: str, checked: bool):
        if not self._current_event or not self._data_source:
            self._show_fit[channel] = False
            if channel in self.fit_buttons:
                self.fit_buttons[channel].setChecked(False)
            return

        data = self.get_fit_data(self._current_event, channel)
        if checked and not data.has_overlay():
            QMessageBox.information(self, "No Fit", f"No fit curves were found for {channel} in this event.")
            self._show_fit[channel] = False
            if channel in self.fit_buttons:
                self.fit_buttons[channel].setChecked(False)
            return

        self._show_fit[channel] = checked
        self.plot_event(self._current_event)

    # ---- Plotting --------------------------------------------------------
    def plot_event(self, event_name: Optional[str]):
        self.figure.clear()
        self._reset_layout_engine()

        if not self._data_source or not event_name:
            self._current_event = None
            ax = self.figure.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "Open a data file and choose an event to visualize.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=18,
            )
            ax.axis("off")
            self.canvas.draw_idle()
            self.refresh_fit_controls()
            self.update_status_text()
            return

        self._current_event = event_name
        self.refresh_fit_controls()

        selected = [name for name in CHANNEL_ORDER if self.channel_buttons.get(name) and self.channel_buttons[name].isChecked()]
        missing: List[str] = []

        self.figure.suptitle(
            f"{os.path.basename(self._filename or '')} — Event {event_name}",
            fontsize=20,
            fontweight="bold",
        )

        if not selected:
            ax = self.figure.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "No channels selected.\nUse the buttons above to enable channels.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=18,
            )
            ax.axis("off")
            self.canvas.draw_idle()
            self.update_status_text()
            return

        overlay_mode = self.overlay_button.isChecked()

        try:
            if overlay_mode:
                families: List[Tuple[str, List[str]]] = []
                high = [ch for ch in selected if CHANNEL_DEFS[ch].family == FAMILY_HIGH]
                low = [ch for ch in selected if CHANNEL_DEFS[ch].family == FAMILY_LOW]
                if high:
                    families.append((FAMILY_HIGH, high))
                if low:
                    families.append((FAMILY_LOW, low))

                if not families:
                    ax = self.figure.add_subplot(111)
                    ax.text(
                        0.5,
                        0.5,
                        "No compatible channels selected for overlay.",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=18,
                    )
                    ax.axis("off")
                else:
                    for idx, (family, channels) in enumerate(families, start=1):
                        ax = self.figure.add_subplot(len(families), 1, idx)
                        plotted_any = False
                        for channel in channels:
                            plotted_any |= self._plot_channel(ax, event_name, channel, overlay_mode=True, missing_channels=missing)
                        self._style_overlay_axis(ax, family, bottom=(idx == len(families)))
                        if plotted_any and (len(channels) > 1 or any(self._show_fit.get(ch, False) for ch in channels)):
                            ax.legend(loc="best")
            else:
                ordered = [ch for ch in CHANNEL_ORDER if ch in selected]
                if not ordered:
                    ax = self.figure.add_subplot(111)
                    ax.text(
                        0.5,
                        0.5,
                        "No channels available to plot.",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=18,
                    )
                    ax.axis("off")
                else:
                    grid = self.figure.add_gridspec(len(ordered), 1)
                    for idx, channel in enumerate(ordered):
                        ax = self.figure.add_subplot(grid[idx, 0])
                        plotted_any = self._plot_channel(ax, event_name, channel, overlay_mode=False, missing_channels=missing)
                        next_family = None
                        if idx + 1 < len(ordered):
                            next_family = CHANNEL_DEFS[ordered[idx + 1]].family
                        self._style_single_axis(
                            ax,
                            channel=channel,
                            bottom=(idx == len(ordered) - 1),
                            next_family=next_family,
                        )
                        if plotted_any and self._show_fit.get(channel) and len(ax.lines) > 1:
                            ax.legend(loc="best")
        except Exception as exc:
            self.figure.clear()
            self._reset_layout_engine()
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Plot error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=16)
            ax.axis("off")

        self.canvas.draw_idle()
        self.update_status_text(missing)

    def _plot_channel(self, ax, event_name: str, channel: str, overlay_mode: bool, missing_channels: List[str]) -> bool:
        definition = CHANNEL_DEFS[channel]
        time_data = self._get_dataset(event_name, definition.time_dataset)
        value_data = self._get_dataset(event_name, definition.dataset)

        base_plotted = False
        reason: Optional[str] = None
        if time_data is None or value_data is None:
            reason = "Missing data"
        else:
            t = _to_1d(time_data)
            y = _to_1d(value_data)
            n = min(len(t), len(y))
            if n == 0:
                reason = "Empty dataset"
            else:
                label = channel if overlay_mode else None
                ax.plot(t[:n], y[:n], label=label)
                base_plotted = True

        fit_plotted = self._plot_fit(ax, event_name, channel, overlay_mode)

        if base_plotted or fit_plotted:
            if not overlay_mode:
                ax.set_ylabel(y_label_with_units(channel), fontsize=16)
            return True

        if reason:
            if not overlay_mode:
                self._draw_missing_message(ax, channel, reason)
            missing_channels.append(channel)
        return False

    def _estimate_baseline(self, event_name: str, channel: str, reference_time: Optional[np.ndarray]) -> float:
        key = (event_name, channel)
        cached = self._baseline_cache.get(key)
        if cached is not None:
            return cached

        baseline_value = 0.0
        if not self._data_source:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        definition = CHANNEL_DEFS.get(channel)
        if not definition:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        raw_data = self._get_dataset(event_name, definition.dataset)
        if raw_data is None:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        raw_values = _to_1d(raw_data)
        if raw_values.size == 0:
            self._baseline_cache[key] = baseline_value
            return baseline_value

        raw_values = np.asarray(raw_values, dtype=float)

        times: Optional[np.ndarray] = None
        if reference_time is not None:
            ref = _to_1d(reference_time)
            if ref.size == raw_values.size:
                times = np.asarray(ref, dtype=float)

        if times is None:
            stored_time = self._get_dataset(event_name, definition.time_dataset)
            if stored_time is not None:
                stored = _to_1d(stored_time)
                if stored.size == raw_values.size:
                    times = np.asarray(stored, dtype=float)

        candidate: Optional[float] = None
        if times is not None and times.size == raw_values.size:
            primary_mask = (times >= BASELINE_PRIMARY_WINDOW[0]) & (times <= BASELINE_PRIMARY_WINDOW[1])
            candidate = _masked_mean(raw_values, primary_mask)
            if candidate is None:
                fallback_mask = times < BASELINE_FALLBACK_THRESHOLD
                candidate = _masked_mean(raw_values, fallback_mask)

        if candidate is None:
            count = min(max(1, raw_values.size // 10), raw_values.size)
            subset = raw_values[:count]
            finite = subset[np.isfinite(subset)]
            if finite.size:
                candidate = float(finite.mean())
            else:
                candidate = 0.0

        baseline_value = float(candidate if np.isfinite(candidate) else 0.0)
        self._baseline_cache[key] = baseline_value
        return baseline_value

    def _iter_fit_curves(self, event_name: str, channel: str, data: FitData):
        yielded = False
        for label, _time_path, _value_path, time_values, fit_values in data.iter_time_result_pairs():
            times = np.asarray(time_values, dtype=float)
            values = np.asarray(fit_values, dtype=float)
            baseline_offset = self._estimate_baseline(event_name, channel, times)
            if baseline_offset:
                values = values + baseline_offset
            n = min(len(times), len(values))
            if n == 0:
                continue
            yield label, times[:n], values[:n]
            yielded = True

        if yielded or not data.parameter_series or not self._data_source:
            return

        definition = CHANNEL_DEFS.get(channel)
        if not definition:
            return
        base_time_data = self._get_dataset(event_name, definition.time_dataset)
        if base_time_data is None:
            return
        base_time = np.asarray(_to_1d(base_time_data), dtype=float)
        if base_time.size == 0:
            return

        model = FIT_MODEL_BY_CHANNEL.get(channel)
        if not model:
            return

        scale = FIT_SCALE_MULTIPLIERS.get(channel, 1.0)
        for path, raw_params in data.iter_parameter_items():
            derived = _fit_paths_from_param(path)
            if derived:
                result_path, time_path = derived
                if result_path in data.value_series and time_path in data.time_series:
                    continue

            param_array = _coerce_parameter_values(raw_params)
            if param_array is None or param_array.size < len(model.parameter_labels):
                continue

            first = np.asarray(param_array[: len(model.parameter_labels)], dtype=float)
            if not np.all(np.isfinite(first)):
                continue

            try:
                curve = model.evaluator(base_time, *first)
            except Exception:
                continue

            fit_values = np.asarray(curve, dtype=float)
            if fit_values.size != base_time.size:
                continue
            if scale != 1.0:
                fit_values = fit_values * scale
            baseline_offset = self._estimate_baseline(event_name, channel, base_time)
            if baseline_offset:
                fit_values = fit_values + baseline_offset

            label = _label_from_param_path(path)
            yield label, base_time, fit_values

    def _plot_fit(self, ax, event_name: str, channel: str, overlay_mode: bool) -> bool:
        if not self._show_fit.get(channel):
            return False

        data = self.get_fit_data(event_name, channel)
        plotted = False
        for label, time_values, fit_values in self._iter_fit_curves(event_name, channel, data):
            n = min(len(time_values), len(fit_values))
            if n == 0:
                continue
            legend_label = label
            if overlay_mode and not label.lower().startswith(channel.lower()):
                legend_label = f"{channel} {label}"
            ax.plot(time_values[:n], fit_values[:n], linestyle="--", linewidth=2.2, label=legend_label)
            plotted = True
        return plotted

    def _style_overlay_axis(self, ax, family: str, bottom: bool):
        ax.set_facecolor("#f8f9fb")
        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=14, width=1.5, length=7)
        ax.set_ylabel(FAMILY_YLABELS.get(family, "Channel"), fontsize=16)
        show_x = bottom or family == FAMILY_HIGH
        if show_x:
            ax.set_xlabel(TIME_AXIS_LABEL, fontsize=16)
        else:
            ax.set_xlabel("")

    def _style_single_axis(self, ax, channel: str, bottom: bool, next_family: Optional[str]):
        ax.set_facecolor("#f8f9fb")
        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=14, width=1.5, length=7)
        current_family = CHANNEL_DEFS.get(channel)
        family_name = current_family.family if current_family else None
        show_x = bottom or (next_family is not None and next_family != family_name)
        if show_x:
            ax.set_xlabel(TIME_AXIS_LABEL, fontsize=16)
        else:
            ax.set_xlabel("")

    def _draw_missing_message(self, ax, channel: str, reason: str):
        ax.set_facecolor("#f8f9fb")
        ax.text(0.5, 0.5, f"{channel}\n{reason}", ha="center", va="center", transform=ax.transAxes, fontsize=16)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    def refresh_fit_controls(self):
        if not self._data_source or not self._current_event:
            for btn in self.fit_buttons.values():
                btn.setEnabled(False)
            return

        for channel, btn in self.fit_buttons.items():
            data = self.get_fit_data(self._current_event, channel)
            available = any(True for _ in self._iter_fit_curves(self._current_event, channel, data))
            btn.setEnabled(available)
            if not available:
                btn.setChecked(False)
                self._show_fit[channel] = False
                btn.setToolTip("No fit curve found for this event.")
            else:
                btn.setToolTip("Overlay the available fit curve on the channel plot.")

    def update_status_text(self, missing: Optional[List[str]] = None):
        parts: List[str] = []
        event_label = self._current_event or "–"
        parts.append(f"Event {event_label}")

        channels = [ch for ch in CHANNEL_ORDER if ch in self.selected_channels and self.channel_buttons.get(ch) and self.channel_buttons[ch].isChecked()]
        parts.append("Channels: " + (", ".join(channels) if channels else "none"))
        parts.append(f"Overlay: {'ON' if self.overlay_button.isChecked() else 'OFF'}")

        fits_on = [ch for ch, state in self._show_fit.items() if state]
        if fits_on:
            parts.append("Fits: " + ", ".join(sorted(fits_on)))

        if missing:
            parts.append("Missing: " + ", ".join(sorted(set(missing))))

        edited = sorted({ch for (evt, ch, _path) in self._fit_param_overrides.keys() if evt == self._current_event})
        if edited:
            parts.append("Edited params: " + ", ".join(edited))

        self.statusBar().showMessage(" | ".join(parts))

    # ---- Fit data helpers -----------------------------------------------
    def get_fit_data(self, event_name: str, channel: str) -> FitData:
        key = (event_name, channel)
        if key in self._fit_cache:
            return self._fit_cache[key]

        data = FitData()
        if not self._data_source:
            self._fit_cache[key] = data
            return data

        base = self._data_source.gather_fit_data(event_name, channel)
        data.time_series = {path: np.array(values, copy=True) for path, values in base.time_series.items()}
        data.value_series = {path: np.array(values, copy=True) for path, values in base.value_series.items()}
        data.parameter_series = {path: np.array(values, copy=True) for path, values in base.parameter_series.items()}
        data.extras = {path: np.array(values, copy=True) for path, values in base.extras.items()}

        for path in list(data.time_series.keys()):
            override_time = self._fit_time_overrides.get((event_name, channel, path))
            if override_time is not None:
                data.time_series[path] = np.array(override_time, copy=True)

        for path in list(data.value_series.keys()):
            override_result = self._fit_result_overrides.get((event_name, channel, path))
            if override_result is not None:
                data.value_series[path] = np.array(override_result, copy=True)

        scale = FIT_SCALE_MULTIPLIERS.get(channel, 1.0)
        if scale != 1.0:
            for path, values in list(data.value_series.items()):
                data.value_series[path] = np.asarray(values, dtype=float) * scale

        self._fit_cache[key] = data
        return data

    def open_fit_parameter_dialog(self):
        if not self._current_event or not self._data_source:
            QMessageBox.information(self, "No Event", "Open an event before editing fit parameters.")
            return

        available: Dict[str, FitData] = {}
        for channel in FIT_ELIGIBLE_CHANNELS:
            data = self.get_fit_data(self._current_event, channel)
            if data.has_parameters():
                available[channel] = data

        if not available:
            QMessageBox.information(self, "No Parameters", "No editable fit parameters were found for this event.")
            return

        dialog = FitParameterDialog(
            parent=self,
            event_name=self._current_event,
            channel_data=available,
            value_getter=self.get_parameter_values,
            save_callback=self.update_fit_override,
            reset_callback=self.clear_fit_override,
        )
        dialog.exec()

    def get_parameter_values(self, event: str, channel: str, dataset_path: str, base_array: np.ndarray) -> np.ndarray:
        key = (event, channel, dataset_path)
        override = self._fit_param_overrides.get(key)
        if override is not None:
            return np.array(override, copy=True)
        return np.array(base_array, copy=True)

    def update_fit_override(self, event: str, channel: str, dataset_path: str, values: np.ndarray):
        array = np.array(values, copy=True)
        self._fit_param_overrides[(event, channel, dataset_path)] = array
        self._remove_fit_override(event, channel, dataset_path)
        recomputed = self._recalculate_fit(event, channel, dataset_path, array)
        if recomputed:
            message = f"Updated fit parameters for {channel}; recomputed fit curve (in-memory only)."
        else:
            message = f"Updated fit parameters for {channel}; fit curve unchanged (missing data)."
        self.statusBar().showMessage(message, 6000)
        self._fit_cache.pop((event, channel), None)
        self.plot_event(self._current_event)

    def clear_fit_override(self, event: str, channel: str, dataset_path: str):
        key = (event, channel, dataset_path)
        removed_param = self._fit_param_overrides.pop(key, None) is not None
        removed_curve = self._remove_fit_override(event, channel, dataset_path)
        if removed_param or removed_curve:
            self.statusBar().showMessage(f"Reverted fit parameters for {channel}.", 6000)
            self._fit_cache.pop((event, channel), None)
            self.plot_event(self._current_event)

    def _recalculate_fit(self, event: str, channel: str, dataset_path: str, params: np.ndarray) -> bool:
        if not self._data_source:
            return False
        derived = _fit_paths_from_param(dataset_path)
        if not derived:
            return False
        result_path, time_path = derived
        time_data = self._get_dataset_by_path(time_path)
        if time_data is None:
            definition = CHANNEL_DEFS.get(channel)
            if definition:
                time_data = self._get_dataset(event, definition.time_dataset)
        if time_data is None:
            return False
        time_array = np.array(time_data, copy=True)
        flat_time = np.asarray(time_array, dtype=float).ravel()
        curve = _evaluate_fit_curve(channel, flat_time, params)
        if curve is None:
            return False
        curve_array = np.array(curve, copy=True)
        if curve_array.size != time_array.size:
            return False
        curve_array = curve_array.reshape(time_array.shape)
        self._fit_time_overrides[(event, channel, time_path)] = time_array
        self._fit_result_overrides[(event, channel, result_path)] = curve_array
        self._fit_cache.pop((event, channel), None)
        return True

    def _remove_fit_override(self, event: str, channel: str, dataset_path: str) -> bool:
        removed = False
        derived = _fit_paths_from_param(dataset_path)
        if derived:
            result_path, time_path = derived
            if (event, channel, result_path) in self._fit_result_overrides:
                del self._fit_result_overrides[(event, channel, result_path)]
                removed = True
            if (event, channel, time_path) in self._fit_time_overrides:
                del self._fit_time_overrides[(event, channel, time_path)]
                removed = True
        return removed

    # --- Cleanup ---------------------------------------------------------
    def closeEvent(self, event):
        if self._data_source is not None:
            try:
                self._data_source.close()
            except Exception:
                pass
        try:
            self._tmpdir.cleanup()
        except Exception:
            pass
        event.accept()

# --------- Fit parameter editor dialog ---------
class FitParameterDialog(QDialog):
    def __init__(
        self,
        parent: QWidget,
        event_name: str,
        channel_data: Dict[str, FitData],
        value_getter: Callable[[str, str, str, np.ndarray], np.ndarray],
        save_callback: Callable[[str, str, str, np.ndarray], None],
        reset_callback: Callable[[str, str, str], None],
    ):
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle(f"Fit Parameters — Event {event_name}")
        self.setMinimumSize(580, 520)

        self._event_name = event_name
        self._channel_data = channel_data
        self._value_getter = value_getter
        self._save_callback = save_callback
        self._reset_callback = reset_callback

        self._param_arrays: Dict[Tuple[str, str], np.ndarray] = {}
        for channel, data in channel_data.items():
            for path, array in data.iter_parameter_items():
                self._param_arrays[(channel, path)] = np.asarray(array)

        self._current_channel: Optional[str] = None
        self._current_dataset: Optional[str] = None
        self._current_shape: Optional[Tuple[int, ...]] = None
        self._is_updating_table = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        header = QLabel(f"Editing parameters for event {event_name}", self)
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        self.formula_label = QLabel("", self)
        self.formula_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.formula_label.setStyleSheet("font-size: 16px;")
        self.formula_label.setWordWrap(True)
        self.formula_label.setMinimumHeight(70)
        self._default_formula_text = "Select a dataset to view its fit function."
        self.formula_label.setText(self._default_formula_text)
        layout.addWidget(self.formula_label)

        chooser_row = QHBoxLayout()
        chooser_row.setSpacing(10)

        self.channel_combo = QComboBox(self)
        self.channel_combo.setStyleSheet("font-size: 15px; min-height: 36px;")
        for channel in sorted(channel_data.keys()):
            self.channel_combo.addItem(channel)
        self.channel_combo.currentTextChanged.connect(self._on_channel_changed)
        chooser_row.addWidget(QLabel("Channel:", self))
        chooser_row.addWidget(self.channel_combo)

        self.dataset_combo = QComboBox(self)
        self.dataset_combo.setStyleSheet("font-size: 15px; min-height: 36px;")
        self.dataset_combo.currentIndexChanged.connect(self._on_dataset_changed)
        chooser_row.addWidget(QLabel("Dataset:", self))
        chooser_row.addWidget(self.dataset_combo)
        chooser_row.addStretch(1)
        layout.addLayout(chooser_row)

        self.table = QTableWidget(self)
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Fit parameter", "Value"])
        header_view = self.table.horizontalHeader()
        header_view.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header_view.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectItems)
        self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.table.setEditTriggers(
            QTableWidget.EditTrigger.DoubleClicked | QTableWidget.EditTrigger.SelectedClicked
        )
        self.table.setStyleSheet("font-size: 15px;")
        self.table.itemChanged.connect(self._on_table_item_changed)
        layout.addWidget(self.table)

        self.info_label = QLabel("", self)
        self.info_label.setStyleSheet("font-size: 14px; color: #555555;")
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)

        self.feedback_label = QLabel("", self)
        self.feedback_label.setStyleSheet("font-size: 14px; color: #146c43;")
        self.feedback_label.setWordWrap(True)
        layout.addWidget(self.feedback_label)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close, self)
        self.apply_button = QPushButton("Apply Changes", self)
        self.apply_button.setStyleSheet("font-size: 15px; font-weight: bold; padding: 8px 18px;")
        self.apply_button.clicked.connect(self.apply_changes)
        buttons.addButton(self.apply_button, QDialogButtonBox.ButtonRole.ActionRole)

        self.reset_button = QPushButton("Reset to File Values", self)
        self.reset_button.setStyleSheet("font-size: 15px; padding: 8px 18px;")
        self.reset_button.clicked.connect(self.reset_changes)
        buttons.addButton(self.reset_button, QDialogButtonBox.ButtonRole.ResetRole)

        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        if self.channel_combo.count() > 0:
            self._on_channel_changed(self.channel_combo.currentText())

    # ---- UI helpers -----------------------------------------------------
    def _on_channel_changed(self, channel: str):
        self.dataset_combo.blockSignals(True)
        self.dataset_combo.clear()
        entries = [path for (chan, path) in self._param_arrays.keys() if chan == channel]
        entries.sort(key=lambda value: _dataset_sort_key(channel, value))
        for path in entries:
            label = self._dataset_display_name(channel, path)
            self.dataset_combo.addItem(label, path)
        self.dataset_combo.blockSignals(False)

        if entries:
            self.dataset_combo.setCurrentIndex(0)
            self._display_dataset(channel, entries[0])
        else:
            self._display_dataset(channel, None)

    def _on_dataset_changed(self, index: int):
        channel = self.channel_combo.currentText()
        dataset_path = self.dataset_combo.itemData(index)
        self._display_dataset(channel, dataset_path)

    def _display_dataset(
        self,
        channel: Optional[str],
        dataset_path: Optional[str],
        *,
        preserve_feedback: bool = False,
    ):
        if not channel or dataset_path is None:
            self._update_formula_display(None)
            self._clear_table("No fit parameters are available for this channel.")
            return

        base_array = self._param_arrays.get((channel, dataset_path))
        if base_array is None:
            self._update_formula_display(None)
            self._clear_table("No fit parameters are available for this channel.")
            return

        self._update_formula_display(channel)

        values = self._value_getter(self._event_name, channel, dataset_path, base_array)
        array = np.asarray(values)
        self._current_channel = channel
        self._current_dataset = dataset_path
        self._current_shape = array.shape

        flat = array.ravel()
        labels = self._parameter_labels(channel, dataset_path, flat.size)

        self._is_updating_table = True
        try:
            self.table.setRowCount(flat.size)
            for idx, value in enumerate(flat):
                label = labels[idx] if idx < len(labels) else f"p{idx + 1}"
                idx_item = QTableWidgetItem(label)
                idx_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                value_item = QTableWidgetItem(self._format_value(value))
                value_item.setFlags(value_item.flags() | Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(idx, 0, idx_item)
                self.table.setItem(idx, 1, value_item)
        finally:
            self._is_updating_table = False

        self.table.resizeColumnsToContents()
        display_name = self._dataset_display_name(channel, dataset_path)
        self.info_label.setText(f"{channel} — {display_name}\n{dataset_path} — shape {array.shape}")
        if not preserve_feedback:
            self.feedback_label.clear()

    def _dataset_display_name(self, channel: str, dataset_path: str) -> str:
        mass_label = _mass_identifier_from_path(dataset_path)
        if mass_label:
            return f"Mass {mass_label}"
        base = os.path.basename(dataset_path)
        return base or dataset_path

    def _parameter_labels(self, channel: str, dataset_path: str, count: int) -> List[str]:
        model = FIT_MODEL_BY_CHANNEL.get(channel)
        if not model:
            return [f"p{i + 1}" for i in range(count)]
        if model is ION_GRID_FIT:
            override = ION_GRID_LABEL_OVERRIDES.get(count)
            if override:
                return list(override)
        labels = list(model.parameter_labels)
        if len(labels) >= count:
            return labels[:count]
        extras = [f"p{i + 1}" for i in range(len(labels), count)]
        return labels + extras

    def _update_formula_display(self, channel: Optional[str]):
        self.formula_label.setPixmap(QPixmap())
        self.formula_label.setToolTip("")
        if not channel:
            self.formula_label.setText(self._default_formula_text)
            return

        model = FIT_MODEL_BY_CHANNEL.get(channel)
        if not model or not model.latex:
            self.formula_label.setText("Fit function preview unavailable for this selection.")
            return

        pixmap = _latex_to_pixmap(model.latex)
        if pixmap:
            self.formula_label.setText("")
            self.formula_label.setPixmap(pixmap)
        else:
            self.formula_label.setText(model.latex)
        self.formula_label.setToolTip(model.latex)

    def _set_feedback(self, text: str, *, success: bool):
        if not text:
            self.feedback_label.clear()
            return
        color = "#146c43" if success else "#b02a37"
        self.feedback_label.setStyleSheet(f"font-size: 14px; color: {color};")
        self.feedback_label.setText(text)

    def _apply_current_values(self, *, auto: bool) -> bool:
        if not self._current_channel or not self._current_dataset:
            if not auto:
                QMessageBox.information(
                    self,
                    "No Selection",
                    "Select a parameter set before applying changes.",
                )
            return False

        values: List[float] = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 1)
            text = item.text().strip() if item else ""
            if not text:
                values.append(0.0)
                continue
            try:
                values.append(float(text))
            except ValueError:
                message = f"Row {row + 1} contains a non-numeric value: {text}"
                if auto:
                    self._set_feedback(message, success=False)
                else:
                    QMessageBox.warning(self, "Invalid value", message)
                return False

        if self._current_shape:
            expected = int(np.prod(self._current_shape))
        else:
            expected = len(values)

        if len(values) != expected:
            message = "The number of edited values does not match the original parameter shape."
            if auto:
                self._set_feedback(message, success=False)
            else:
                QMessageBox.warning(self, "Shape mismatch", message)
            return False

        array = np.array(values, dtype=float)
        if self._current_shape:
            array = array.reshape(self._current_shape)

        self._save_callback(self._event_name, self._current_channel, self._current_dataset, array)
        self._display_dataset(self._current_channel, self._current_dataset, preserve_feedback=True)

        if auto:
            self._set_feedback("Fit updated with latest parameters.", success=True)
        else:
            self._set_feedback("Changes stored in-memory and fit recomputed.", success=True)
        return True

    def _on_table_item_changed(self, item: QTableWidgetItem):
        if self._is_updating_table or item is None or item.column() != 1:
            return
        self._apply_current_values(auto=True)

    def _clear_table(self, message: str = ""):
        self.table.setRowCount(0)
        self.info_label.setText(message)
        self.feedback_label.clear()
        self._current_channel = None
        self._current_dataset = None
        self._current_shape = None
        self._is_updating_table = False

    def _format_value(self, value: object) -> str:
        try:
            if isinstance(value, (np.floating, float)):
                return f"{float(value):.6g}"
            if isinstance(value, (np.integer, int)):
                return str(int(value))
        except Exception:
            pass
        return str(value)

    # ---- Actions --------------------------------------------------------
    def apply_changes(self):
        self._apply_current_values(auto=False)

    def reset_changes(self):
        if not self._current_channel or not self._current_dataset:
            return
        self._reset_callback(self._event_name, self._current_channel, self._current_dataset)
        self._display_dataset(self._current_channel, self._current_dataset, preserve_feedback=True)
        self._set_feedback("Reverted to values from the file.", success=False)

# --------- CLI / main ---------
def main():
    parser = argparse.ArgumentParser(description="Run the IDEX Quicklook (QtAgg + Matplotlib toolbar).")
    parser.add_argument("--filename", nargs="?", default=None, help="Path to the data file (HDF5 or CDF).")
    parser.add_argument("--eventnumber", nargs="?", type=int, default=None, help="1-based event index.")
    args = parser.parse_args()
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    base_font = QFont(app.font())
    size = base_font.pointSize()
    if size <= 0:
        size = 12
    else:
        size = max(size, 12)
    base_font.setPointSize(size)
    app.setFont(base_font)

    w = MainWindow(filename=args.filename, eventnumber=args.eventnumber)
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
