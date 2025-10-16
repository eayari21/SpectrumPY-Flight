"""Helpers for parsing and inspecting the IDEX CDF variable definitions."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union
from zipfile import ZipFile

import numpy as np

_SPREADSHEET_NAME = "IDEX CDF Variable Definitions.xlsx"
_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_COEFFICIENT_COLUMNS = ("CO", "C1", "C2", "C3", "C4", "C5", "C6", "C7")

NumberLike = Union[int, float, np.ndarray]


def _normalize_text(value: str) -> str:
    """Return a case-insensitive, quote-normalised representation of ``value``."""

    return " ".join(value.replace("“", '"').replace("”", '"').strip().split()).lower()


def _parse_optional_int(text: str) -> Optional[int]:
    text = text.strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _parse_optional_float(text: str) -> Optional[float]:
    text = text.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _parse_dn_range(text: str) -> Optional[Tuple[int, int]]:
    text = text.strip()
    if not text:
        return None
    match = re.match(r"^(-?\d+)\s*-\s*(-?\d+)$", text)
    if not match:
        return None
    start, stop = match.groups()
    return int(start), int(stop)


def _read_spreadsheet_rows(path: Path) -> List[List[str]]:
    """Return the spreadsheet content as a list of row lists."""

    with ZipFile(path) as archive:
        shared_strings: List[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = archive.read("xl/sharedStrings.xml")
            shared_strings = _extract_shared_strings(shared_root)

        sheet_root = archive.read("xl/worksheets/sheet1.xml")
        return _extract_sheet_rows(sheet_root, shared_strings)


def _extract_shared_strings(xml_bytes: bytes) -> List[str]:
    import xml.etree.ElementTree as ET

    root = ET.fromstring(xml_bytes)
    values = []
    for node in root.findall(f"{_NS}si"):
        values.append("".join(node.itertext()))
    return values


def _extract_sheet_rows(xml_bytes: bytes, shared_strings: Sequence[str]) -> List[List[str]]:
    import xml.etree.ElementTree as ET

    root = ET.fromstring(xml_bytes)
    rows: List[List[str]] = []

    def col_to_index(column: str) -> int:
        index = 0
        for char in column:
            if char.isalpha():
                index = index * 26 + (ord(char.upper()) - ord("A") + 1)
        return index - 1

    for row in root.findall(f".//{_NS}row"):
        values: List[str] = []
        current_col = 0
        for cell in row.findall(f"{_NS}c"):
            ref = cell.get("r")
            if ref:
                col_letters = "".join(ch for ch in ref if ch.isalpha())
                target = col_to_index(col_letters)
                while current_col < target:
                    values.append("")
                    current_col += 1

            value_node = cell.find(f"{_NS}v")
            if value_node is None:
                value = ""
            elif cell.get("t") == "s":
                value = shared_strings[int(value_node.text)]
            else:
                value = value_node.text or ""

            values.append(value)
            current_col += 1

        rows.append(values)

    return rows


@dataclass
class PolynomialSegment:
    """Piecewise polynomial description for a DN range."""

    dn_range: Tuple[int, int]
    coefficients: Tuple[float, ...]

    def contains(self, dn: NumberLike) -> np.ndarray:
        values = np.asarray(dn)
        return (values >= self.dn_range[0]) & (values <= self.dn_range[1])

    def evaluate(self, dn: NumberLike) -> np.ndarray:
        values = np.asarray(dn, dtype=float)
        return np.polynomial.polynomial.polyval(values, self.coefficients)


@dataclass
class VariableDefinition:
    """Container describing an entry from the IDEX variable definitions sheet."""

    cdf_varname: str
    item_name: str
    packet_variable: str
    starting_bit: Optional[int]
    padding_bits: Optional[int]
    unsigned_bits: Optional[int]
    var_notes: str
    summary: str
    units: str
    segments: List[PolynomialSegment] = field(default_factory=list)

    def evaluate(self, dn: NumberLike) -> NumberLike:
        if not self.segments:
            return dn

        arr = np.asarray(dn, dtype=float)
        scalar = arr.ndim == 0
        flat = arr.reshape(-1)
        result = np.empty_like(flat, dtype=float)
        filled = np.zeros_like(flat, dtype=bool)

        for segment in self.segments:
            mask = (flat >= segment.dn_range[0]) & (flat <= segment.dn_range[1])
            if not np.any(mask):
                continue
            result[mask] = segment.evaluate(flat[mask])
            filled[mask] = True

        first, last = self.segments[0], self.segments[-1]
        below = flat < first.dn_range[0]
        if np.any(below):
            result[below] = first.evaluate(flat[below])
            filled[below] = True

        above = flat > last.dn_range[1]
        if np.any(above):
            result[above] = last.evaluate(flat[above])
            filled[above] = True

        if not filled.all():
            remaining = ~filled
            result[remaining] = flat[remaining]

        if scalar:
            return float(result[0])
        return result.reshape(arr.shape)

    def dn_range(self) -> Optional[Tuple[int, int]]:
        if not self.segments:
            return None
        start = min(seg.dn_range[0] for seg in self.segments)
        stop = max(seg.dn_range[1] for seg in self.segments)
        return start, stop


class VariableDefinitionsCatalog:
    """Catalog of :class:`VariableDefinition` entries keyed by helpful identifiers."""

    def __init__(self, definitions: Iterable[VariableDefinition]):
        self.definitions: List[VariableDefinition] = list(definitions)
        self._by_notes: Dict[str, VariableDefinition] = {}
        self._by_packet: Dict[str, VariableDefinition] = {}
        self._by_cdf_name: Dict[str, VariableDefinition] = {}

        for definition in self.definitions:
            if definition.var_notes:
                self._by_notes[_normalize_text(definition.var_notes)] = definition
            if definition.packet_variable:
                self._by_packet[_normalize_text(definition.packet_variable)] = definition
            if definition.cdf_varname:
                self._by_cdf_name[_normalize_text(definition.cdf_varname)] = definition

    def find_by_var_notes(self, var_notes: str) -> Optional[VariableDefinition]:
        return self._by_notes.get(_normalize_text(var_notes))

    def find_by_packet_variable(self, packet_name: str) -> Optional[VariableDefinition]:
        return self._by_packet.get(_normalize_text(packet_name))

    def find_by_cdf_varname(self, cdf_varname: str) -> Optional[VariableDefinition]:
        return self._by_cdf_name.get(_normalize_text(cdf_varname))


def _build_definitions_from_rows(rows: Sequence[Sequence[str]]) -> List[VariableDefinition]:
    if not rows:
        return []

    header = list(rows[0])
    column_index = {name: idx for idx, name in enumerate(header)}

    def get(row: Sequence[str], key: str) -> str:
        idx = column_index.get(key)
        if idx is None or idx >= len(row):
            return ""
        return row[idx]

    definitions: List[VariableDefinition] = []
    current: Optional[VariableDefinition] = None

    for raw_row in rows[1:]:
        if not any(value.strip() for value in raw_row):
            continue

        cdf_varname = get(raw_row, "CDF Varname").strip()
        if cdf_varname:
            current = VariableDefinition(
                cdf_varname=cdf_varname,
                item_name=get(raw_row, "Item Name (Algorithm Template)").strip(),
                packet_variable=get(raw_row, "Variable Name in Packet").strip(),
                starting_bit=_parse_optional_int(get(raw_row, "Starting Bit")),
                padding_bits=_parse_optional_int(get(raw_row, "Nbits padding before")),
                unsigned_bits=_parse_optional_int(get(raw_row, "Unsigned nbits")),
                var_notes=get(raw_row, "Var_notes").strip(),
                summary=get(raw_row, "Plain Language Summary").strip(),
                units=get(raw_row, "Units").strip(),
            )
            definitions.append(current)

        if current is None:
            continue

        dn_range = _parse_dn_range(get(raw_row, "DN Range"))
        if dn_range is None:
            continue

        coefficients: List[float] = []
        for key in _COEFFICIENT_COLUMNS:
            value = get(raw_row, key).strip()
            coeff = _parse_optional_float(value) or 0.0
            coefficients.append(coeff)

        current.segments.append(
            PolynomialSegment(dn_range=dn_range, coefficients=tuple(coefficients))
        )

    return definitions


@lru_cache(maxsize=1)
def load_variable_definitions(path: Optional[Path] = None) -> VariableDefinitionsCatalog:
    """Load and cache the spreadsheet contents as a catalog."""

    if path is None:
        path = Path(__file__).resolve().parent / _SPREADSHEET_NAME

    rows = _read_spreadsheet_rows(path)
    definitions = _build_definitions_from_rows(rows)
    return VariableDefinitionsCatalog(definitions)


def evaluate_by_var_notes(var_notes: str, dn: NumberLike, *, path: Optional[Path] = None) -> Optional[float]:
    catalog = load_variable_definitions(path)
    definition = catalog.find_by_var_notes(var_notes)
    if definition is None:
        return None
    return definition.evaluate(dn)


__all__ = [
    "PolynomialSegment",
    "VariableDefinition",
    "VariableDefinitionsCatalog",
    "evaluate_by_var_notes",
    "load_variable_definitions",
]
