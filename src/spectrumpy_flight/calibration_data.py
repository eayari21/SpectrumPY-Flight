"""Utilities for accelerator metadata and calibration matrix lookups."""

from __future__ import annotations

import csv
import math
import zipfile
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

try:  # Optional dependency used for live accelerator SQL queries
    from .idex_packet import SQLMatchCriteria, SQLMatchResult, query_dust_events
except Exception:  # pragma: no cover - optional at runtime
    SQLMatchCriteria = None  # type: ignore[assignment]
    SQLMatchResult = None  # type: ignore[assignment]
    query_dust_events = None  # type: ignore[assignment]

from . import package_path

__all__ = [
    "AcceleratorMatch",
    "AcceleratorMatchFinder",
    "CalibrationEntry",
    "CalibrationMatrix",
]


EXCEL_EPOCH = datetime(1899, 12, 30)


def _excel_serial_to_datetime(value: str | float | int | None) -> datetime | None:
    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    # Excel incorrectly treats 1900 as a leap year, so offset accordingly.
    return EXCEL_EPOCH + timedelta(days=numeric)


def _column_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha()).upper()
    value = 0
    for ch in letters:
        value = value * 26 + (ord(ch) - ord("A") + 1)
    return value


def _row_index(cell_ref: str) -> int:
    digits = "".join(ch for ch in cell_ref if ch.isdigit())
    return int(digits) if digits else 0


def _load_xlsx_sheet(path: Path, sheet_rel_path: str) -> list[tuple[int, dict[int, str]]]:
    """Return rows for the requested worksheet as column-indexed dictionaries."""

    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            from xml.etree import ElementTree as ET

            ns = {"t": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            xml = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for si in xml.findall(".//t:si", ns):
                text = "".join(node.text or "" for node in si.findall(".//t:t", ns))
                shared_strings.append(text)

        from xml.etree import ElementTree as ET

        ns = {"t": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        xml = ET.fromstring(archive.read(f"xl/{sheet_rel_path}"))
        rows: list[tuple[int, dict[int, str]]] = []
        for row in xml.findall(".//t:sheetData/t:row", ns):
            row_idx_text = row.attrib.get("r")
            if not row_idx_text:
                continue
            row_idx = int(row_idx_text)
            mapping: dict[int, str] = {}
            for cell in row.findall("t:c", ns):
                ref = cell.attrib.get("r")
                if not ref:
                    continue
                col_idx = _column_index(ref)
                cell_type = cell.attrib.get("t")
                value_elem = cell.find("t:v", ns)
                text_elem = cell.find("t:is/t:t", ns)
                if value_elem is None and text_elem is not None:
                    value = text_elem.text or ""
                elif value_elem is None:
                    continue
                else:
                    raw_value = value_elem.text or ""
                    if cell_type == "s":
                        try:
                            idx = int(raw_value)
                        except ValueError:
                            idx = -1
                        value = shared_strings[idx] if 0 <= idx < len(shared_strings) else ""
                    else:
                        value = raw_value
                mapping[col_idx] = value.strip()
            if mapping:
                rows.append((row_idx, mapping))
        return rows


def _sheet_reference(path: Path, sheet_name: str) -> str | None:
    with zipfile.ZipFile(path) as archive:
        from xml.etree import ElementTree as ET

        ns = {"t": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib.get("Id"): rel.attrib.get("Target")
            for rel in rels.findall(
                "{http://schemas.openxmlformats.org/package/2006/relationships}Relationship"
            )
        }
        for sheet in workbook.findall("t:sheets/t:sheet", ns):
            name = sheet.attrib.get("name")
            if name and name.strip() == sheet_name:
                rel_id = sheet.attrib.get(
                    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
                )
                if rel_id and rel_id in rel_map:
                    return rel_map[rel_id]
    return None


@dataclass(frozen=True)
class AcceleratorMatch:
    """Description of a matched accelerator record."""

    record_id: int | None
    timestamp_ms: float
    mass_kg: float
    velocity_mps: float
    charge_c: float
    radius_m: float
    estimate_quality: float
    source: str
    experiment_name: str | None = None
    experiment_description: str | None = None
    dust_type: str | None = None
    group_id: str | None = None


@dataclass(frozen=True)
class CalibrationEntry:
    """Metadata for a dust accelerator run extracted from the test matrix."""

    campaign: str
    date: datetime | None
    run_number: int | None
    material: str | None
    target_location: str | None
    azimuthal_location: str | None
    speed_range: str | None
    reference_voltage: float | None
    target_voltage: float | None
    detector_voltage: float | None
    rejection_voltage: float | None
    csas_used: str | None
    notes: str | None

    @property
    def experiment_key(self) -> str | None:
        if self.date is None or self.run_number is None:
            return None
        return f"{self.date:%Y_%m_%d}_run{self.run_number:d}"


class CalibrationMatrix:
    """Load and query the calibration test matrix workbook."""

    def __init__(self, path: Path | None = None) -> None:
        if path is None:
            path = package_path("lookup", "IDEX Calibration Test Matrix Post Env Dust.xlsx")
        self.path = path
        self.entries: dict[str, CalibrationEntry] = {}
        self.campaign_overview: list[str] = []
        self._load()

    def _load(self) -> None:
        if not self.path.exists():
            return
        fm_sheet = _sheet_reference(self.path, "FM Dust Testing")
        olivine_sheet = _sheet_reference(self.path, "Olivine 12-3 thru 12-15")
        aluminum_sheet = _sheet_reference(self.path, "Al 12-18 thru 12-20")
        if fm_sheet:
            self.campaign_overview = self._summarise_campaign(fm_sheet)
        for rel in (olivine_sheet, aluminum_sheet):
            if rel:
                self._ingest_runs(rel)

    def _summarise_campaign(self, rel_path: str) -> list[str]:
        rows = _load_xlsx_sheet(self.path, rel_path)
        if not rows:
            return []
        header_row = rows[0][1]
        column_lookup = {
            _column_index_ref: header_row[_column_index_ref]
            for _column_index_ref in sorted(header_row)
        }
        entries: list[str] = []
        for row_idx, mapping in rows[1:]:
            if row_idx <= 1:
                continue
            label = mapping.get(1)
            if not label:
                continue
            row_text = [label]
            for col_idx in sorted(mapping):
                if col_idx == 1:
                    continue
                header = column_lookup.get(col_idx)
                if not header:
                    continue
                value = mapping[col_idx]
                if value:
                    row_text.append(f"{header}: {value}")
            if len(row_text) > 1:
                entries.append("; ".join(row_text))
        return entries

    def _ingest_runs(self, rel_path: str) -> None:
        rows = _load_xlsx_sheet(self.path, rel_path)
        if not rows:
            return
        header_row = rows[0][1]
        column_names: dict[int, str] = {}
        for col_idx, name in header_row.items():
            column_names[col_idx] = name
        for _row_idx, mapping in rows[1:]:
            date_value = mapping.get(2)
            run_value = mapping.get(3)
            if run_value is None:
                continue
            try:
                run_number = int(float(run_value))
            except (TypeError, ValueError):
                run_number = None
            date = _excel_serial_to_datetime(date_value)
            material = mapping.get(7)
            target_location = mapping.get(8)
            azimuthal = mapping.get(9)
            speed_range = mapping.get(10)
            ref_voltage = _coerce_float(mapping.get(11))
            target_voltage = _coerce_float(mapping.get(12))
            detector_voltage = _coerce_float(mapping.get(13))
            rejection_voltage = _coerce_float(mapping.get(14))
            csas = mapping.get(15) or mapping.get(14)
            notes = mapping.get(16)
            campaign = mapping.get(1) or column_names.get(1, "Run")
            entry = CalibrationEntry(
                campaign=campaign.strip() if campaign else "Run",
                date=date,
                run_number=run_number,
                material=material or None,
                target_location=target_location or None,
                azimuthal_location=azimuthal or None,
                speed_range=speed_range or None,
                reference_voltage=ref_voltage,
                target_voltage=target_voltage,
                detector_voltage=detector_voltage,
                rejection_voltage=rejection_voltage,
                csas_used=(csas or None),
                notes=notes or None,
            )
            key = entry.experiment_key
            if key:
                self.entries[key] = entry

    def lookup(self, experiment_name: str | None) -> CalibrationEntry | None:
        if not experiment_name:
            return None
        key = experiment_name.strip()
        if not key:
            return None
        return self.entries.get(key)

    def overview_text(self) -> list[str]:
        return self.campaign_overview


def _coerce_float(value: str | float | int | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


class AcceleratorMatchFinder:
    """Locate accelerator metadata via SQL or local CSV fallbacks."""

    def __init__(
        self,
        csv_path: Path | None = None,
        *,
        time_tolerance_ms: float = 2_000.0,
        timezone_offsets_hours: Sequence[int] = (0, 6, -6, 7, -7),
    ) -> None:
        if csv_path is None:
            csv_path = package_path("lookup", "IDEX_FM_2023.csv")
        self.csv_path = csv_path
        self.time_tolerance_ms = float(time_tolerance_ms)
        self._timezone_offsets_ms = tuple(
            sorted({int(hours * 3_600_000) for hours in timezone_offsets_hours})
        )
        self._csv_records = self._load_csv_records()
        self._csv_timestamps = np.array(
            [record["timestamp_ms"] for record in self._csv_records], dtype=float
        )
        self._server_successes = 0

    def _load_csv_records(self) -> list[dict[str, float | str | None]]:
        records: list[dict[str, float | str | None]] = []
        if not self.csv_path.exists():
            return records
        with self.csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                try:
                    timestamp_ms = float(row.get("Integer Timestamp (ms)"))
                    velocity = float(row.get("Velocity (m/s)"))
                    mass = float(row.get("Mass (kg)"))
                    charge = float(row.get("Charge (C)"))
                    radius = float(row.get("Radius (m) "))
                except (TypeError, ValueError):
                    continue
                if not all(math.isfinite(v) for v in (timestamp_ms, velocity, mass, charge, radius)):
                    continue
                estimate_quality = _coerce_float(row.get("Estimate Quality")) or 0.0
                record_id = row.get("Dust Event ID")
                experiment_name = (row.get("Current Experiment Name") or "").strip()
                experiment_description = (row.get("Current Experiment Description") or "").strip()
                dust_type = (row.get("Dust Type") or "").strip()
                group_id = (row.get("Current Group ID") or "").strip()
                records.append(
                    {
                        "record_id": int(record_id) if record_id else None,
                        "timestamp_ms": float(timestamp_ms),
                        "velocity_mps": float(velocity),
                        "mass_kg": float(mass),
                        "charge_c": float(charge),
                        "radius_m": float(radius),
                        "estimate_quality": float(estimate_quality),
                        "experiment_name": experiment_name or None,
                        "experiment_description": experiment_description or None,
                        "dust_type": dust_type or None,
                        "group_id": group_id or None,
                    }
                )
        records.sort(key=lambda item: item["timestamp_ms"])
        return records

    @property
    def server_successes(self) -> int:
        return self._server_successes

    def find(
        self,
        instrument_timestamp_ms: float,
        *,
        timezone_offset_ms: float = 0.0,
        velocity_mps: float | None = None,
    ) -> AcceleratorMatch | None:
        match = self._query_server(instrument_timestamp_ms, timezone_offset_ms, velocity_mps)
        if match is not None:
            self._server_successes += 1
            return match
        return self._query_csv(instrument_timestamp_ms, timezone_offset_ms, velocity_mps)

    def _query_server(
        self,
        timestamp_ms: float,
        timezone_offset_ms: float,
        velocity_mps: float | None,
    ) -> AcceleratorMatch | None:
        if SQLMatchCriteria is None or query_dust_events is None:
            return None
        offsets_ms = self._candidate_offsets(timezone_offset_ms)
        for offset in offsets_ms:
            adjusted_time = timestamp_ms + offset
            criteria = SQLMatchCriteria(
                time_ms=adjusted_time,
                time_window_ms=self.time_tolerance_ms,
                velocity_kmps=(velocity_mps / 1000.0) if velocity_mps is not None else None,
                velocity_tolerance_kmps=0.3,
                min_quality=3,
                restrict_time=True,
                restrict_velocity=velocity_mps is not None,
                limit=10,
            )
            try:
                results, _sql, _params = query_dust_events(criteria)
            except Exception:
                continue
            match = self._choose_best_result(results, adjusted_time)
            if match is not None:
                return AcceleratorMatch(
                    record_id=match.record_id,
                    timestamp_ms=float(match.timestamp_ms or 0.0),
                    mass_kg=float(match.mass_kg or 0.0),
                    velocity_mps=float(match.velocity_mps or 0.0),
                    charge_c=float(match.charge_c or 0.0),
                    radius_m=float(match.radius_m or 0.0),
                    estimate_quality=float(match.estimate_quality or 0.0),
                    source="server",
                    experiment_name=(match.experiment_tag or match.metadata.get("ExperimentTag") if match.metadata else None),
                    experiment_description=(
                        match.experiment_description
                        or (match.metadata.get("ExperimentDescription") if match.metadata else None)
                    ),
                    dust_type=(match.metadata.get("DustSourceNotes") if match.metadata else None),
                    group_id=str(match.experiment_settings_id) if match.experiment_settings_id is not None else None,
                )
        return None

    def _candidate_offsets(self, initial_offset_ms: float) -> Iterable[float]:
        offsets = [initial_offset_ms]
        for candidate in self._timezone_offsets_ms:
            if candidate not in offsets:
                offsets.append(float(candidate))
            neg = -float(candidate)
            if neg not in offsets:
                offsets.append(neg)
        seen = set()
        for offset in offsets:
            if offset in seen:
                continue
            seen.add(offset)
            yield offset

    def _choose_best_result(
        self,
        results: Sequence[SQLMatchResult],
        reference_time_ms: float,
    ) -> SQLMatchResult | None:
        best: SQLMatchResult | None = None
        best_score = float("inf")
        for result in results:
            if result is None:
                continue
            timestamp = result.timestamp_ms
            if timestamp is None or not math.isfinite(timestamp):
                continue
            quality = result.estimate_quality if result.estimate_quality is not None else 0.0
            if quality < 3:
                continue
            delta = abs(float(timestamp) - reference_time_ms)
            if delta > self.time_tolerance_ms:
                continue
            score = delta - quality * 10.0
            if score < best_score:
                best_score = score
                best = result
        return best

    def _query_csv(
        self,
        timestamp_ms: float,
        timezone_offset_ms: float,
        velocity_mps: float | None,
    ) -> AcceleratorMatch | None:
        if not self._csv_records:
            return None
        best_record: dict[str, float | str | None] | None = None
        best_delta = float("inf")
        for offset in self._candidate_offsets(timezone_offset_ms):
            adjusted = float(timestamp_ms + offset)
            idx = int(np.searchsorted(self._csv_timestamps, adjusted))
            for candidate_idx in (idx - 1, idx, idx + 1):
                if candidate_idx < 0 or candidate_idx >= len(self._csv_records):
                    continue
                record = self._csv_records[candidate_idx]
                delta = abs(record["timestamp_ms"] - adjusted)
                if delta > self.time_tolerance_ms:
                    continue
                if velocity_mps is not None and math.isfinite(velocity_mps):
                    if not math.isfinite(record["velocity_mps"]):
                        continue
                    if abs(record["velocity_mps"] - velocity_mps) > 1_500.0:
                        continue
                if delta < best_delta:
                    best_delta = delta
                    best_record = record
        if best_record is None:
            return None
        return AcceleratorMatch(
            record_id=best_record["record_id"] if best_record["record_id"] is not None else None,
            timestamp_ms=float(best_record["timestamp_ms"]),
            mass_kg=float(best_record["mass_kg"]),
            velocity_mps=float(best_record["velocity_mps"]),
            charge_c=float(best_record["charge_c"]),
            radius_m=float(best_record["radius_m"]),
            estimate_quality=float(best_record["estimate_quality"]),
            source="csv",
            experiment_name=best_record.get("experiment_name") if isinstance(best_record.get("experiment_name"), str) else None,
            experiment_description=best_record.get("experiment_description") if isinstance(best_record.get("experiment_description"), str) else None,
            dust_type=best_record.get("dust_type") if isinstance(best_record.get("dust_type"), str) else None,
            group_id=best_record.get("group_id") if isinstance(best_record.get("group_id"), str) else None,
        )

