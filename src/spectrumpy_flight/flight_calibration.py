"""Flight calibration analysis and reporting utilities."""

from __future__ import annotations

import json
import math
import re
import sys
import textwrap
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Mapping, Sequence, cast

import numpy as np

try:  # pragma: no cover - optional dependency for runtime environments
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except Exception:  # pragma: no cover - matplotlib may be unavailable in some contexts
    plt = None  # type: ignore[assignment]
    PdfPages = None  # type: ignore[assignment]

from . import package_path
from .calibration_data import (
    AcceleratorMatch,
    AcceleratorMatchFinder,
    CalibrationEntry,
    CalibrationMatrix,
)
from .dust_schedule import DustScheduleEntry, load_dust_schedule
from .olivine_metrics import FE_ISOTOPES, MG_ISOTOPES, SI_ISOTOPES

__all__ = [
    "DustScheduleEntry",
    "load_dust_schedule",
    "FlightCalibrationAnalyzer",
    "generate_flight_calibration_report",
]




def _discover_data_files(root: Path) -> list[Path]:
    """Return data products sorted chronologically, preferring decoded HDF5."""

    preferred: dict[datetime, Path] = {}
    for path in sorted(root.rglob("ois_output_*")):
        if path.is_dir():
            continue
        timestamp = _parse_filename_timestamp(path)
        if timestamp is None:
            continue
        existing = preferred.get(timestamp)
        suffix = path.suffix.lower()
        if existing is None:
            preferred[timestamp] = path
            continue
        existing_suffix = existing.suffix.lower()
        if existing_suffix != ".h5" and suffix == ".h5":
            preferred[timestamp] = path
        elif existing_suffix == ".h5" and suffix != ".h5":
            continue
        elif suffix == existing_suffix == ".h5" and str(path) < str(existing):
            preferred[timestamp] = path
    return [preferred[key] for key in sorted(preferred)]


def _parse_filename_timestamp(path: Path) -> datetime | None:
    stem = path.stem if path.suffix else path.name
    parts = stem.split("ois_output_")
    if len(parts) != 2:
        return None
    timestamp = parts[1]
    try:
        dt = datetime.strptime(timestamp, "%Y%m%d_%H%M%S")
    except ValueError:
        return None
    return dt.replace(tzinfo=timezone.utc)


def _normalise_epoch_to_ms(value: float) -> float | None:
    if not np.isfinite(value):
        return None
    magnitude = abs(value)
    if magnitude >= 1e12:  # already in milliseconds
        return float(value)
    if magnitude >= 1e9:  # seconds since epoch
        return float(value) * 1_000.0
    if magnitude >= 1e6:  # microseconds
        return float(value) / 1_000.0
    if magnitude >= 1e3:  # milliseconds but small absolute
        return float(value)
    return float(value) * 1_000.0


def _read_scalar(dataset) -> float | None:
    if dataset is None:
        return None
    try:
        data = dataset[()]
    except Exception:
        return None
    arr = np.asarray(data).ravel()
    for value in arr:
        try:
            numeric = float(value)
        except Exception:
            continue
        if np.isfinite(numeric):
            return numeric
    return None


def _read_text(dataset) -> str | None:
    if dataset is None:
        return None
    try:
        value = dataset[()]
    except Exception:
        return None
    if isinstance(value, bytes):
        try:
            text = value.decode("utf-8")
        except Exception:
            text = value.decode("latin1", "ignore")
    else:
        text = str(value)
    text = text.strip()
    return text or None


def _ensure_matplotlib() -> None:
    if plt is None or PdfPages is None:  # pragma: no cover - runtime guard
        raise RuntimeError(
            "Matplotlib is required to generate the flight calibration report."
        )


_EXPERIMENT_PATTERN = re.compile(r"(20\d{2})[_-](\d{2})[_-](\d{2})[_\-]run(\d+)", re.IGNORECASE)
_RUN_PATTERN = re.compile(r"run\s*(\d+)", re.IGNORECASE)
_DATE_PATTERN = re.compile(r"(20\d{2})[-_/](\d{2})[-_/](\d{2})")


def _experiment_name_candidates(
    experiment_name: str | None, experiment_description: str | None
) -> list[str]:
    candidates: list[str] = []
    if experiment_name:
        text = experiment_name.strip()
        if text:
            candidates.append(text)
    if experiment_description:
        description = experiment_description.strip()
        if description:
            found = _EXPERIMENT_PATTERN.search(description)
            if found:
                year, month, day, run = found.groups()
                candidates.append(f"{year}_{month}_{day}_run{int(run):d}")
            else:
                date_match = _DATE_PATTERN.search(description)
                run_match = _RUN_PATTERN.search(description)
                if date_match and run_match:
                    year, month, day = date_match.groups()
                    run_value = int(run_match.group(1))
                    candidates.append(f"{year}_{month}_{day}_run{run_value:d}")
    return candidates


def _infer_experiment_name(match: AcceleratorMatch | None) -> str | None:
    if match is None:
        return None
    for candidate in _experiment_name_candidates(
        match.experiment_name, match.experiment_description
    ):
        if candidate:
            return candidate
    return None


_RSF_PRESETS: list[tuple[str, dict[str, float]]] = [
    ("Measured abundances", {"Mg": 1.0, "Si": 1.0, "Fe": 1.0}),
    ("PUMA/PIA (Krueger 1996)", {"Mg": 3.1, "Si": 1.0, "Fe": 1.1}),
    ("Orthopyroxene, LAMA (Sternglass 1971)", {"Mg": 5.50, "Si": 1.0, "Fe": 1.12}),
    ("Olivine Fo87, SUDA (Hillier et al. 2018)", {"Mg": 4.93, "Si": 1.0, "Fe": 1.50}),
    ("Olivine Fo91, Hyperdust (this work)", {"Mg": 4.97, "Si": 1.0, "Fe": 1.32}),
    ("TOF SIMS (Stephan 2001)", {"Mg": 5.10, "Si": 1.0, "Fe": 2.40}),
]


@dataclass
class EventRecord:
    file: Path
    event_id: str
    dust_type: str
    instrument_model: str
    timestamp: datetime
    accelerator_timestamp_ms: float | None
    accelerator_mass_kg: float | None
    accelerator_velocity_mps: float | None
    accelerator_charge_c: float | None
    accelerator_radius_m: float | None
    accelerator_source: str
    accelerator_estimate_quality: float | None
    has_accelerator_match: bool
    accelerator_experiment_name: str | None
    accelerator_experiment_description: str | None
    instrument_mass_estimate_kg: float | None
    instrument_velocity_estimate_mps: float | None
    target_rise_time_us: float | None
    target_velocity_fit_mps: float | None
    target_velocity_rise_mps: float | None
    target_velocity_ratio_mps: float | None
    collection_efficiency: float | None
    mass_resolutions: list[float]
    snr_by_channel: dict[str, float]
    chi_sq_by_channel: dict[str, float]
    reduced_chi_sq_by_channel: dict[str, float]
    ternary_point: tuple[float, float, float] | None
    mass_line_species: tuple[str, ...]
    mass_line_areas: dict[str, float]
    calibration_campaign: str | None
    calibration_material: str | None
    calibration_speed_range: str | None
    calibration_target_location: str | None
    calibration_azimuthal_location: str | None
    calibration_notes: str | None
    reference_voltage: float | None
    target_voltage: float | None
    detector_voltage: float | None
    rejection_voltage: float | None


def _mass_resolution(mass: float, sigma: float) -> float | None:
    if not np.isfinite(mass) or not np.isfinite(sigma):
        return None
    if sigma <= 0.0:
        return None
    return float(mass / (2.3548200450309493 * sigma))


def _decode_species(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "ignore").strip()
    return str(value).strip()


def _mass_lines_to_areas(table: np.ndarray) -> dict[str, float]:
    areas: dict[str, float] = {}
    for row in table:
        species = _decode_species(row["assigned_species"])
        if not species:
            continue
        area = float(row["area"])
        if not np.isfinite(area):
            continue
        areas[species] = areas.get(species, 0.0) + max(area, 0.0)
    return areas


def _ternary_from_areas(areas: Mapping[str, float]) -> tuple[float, float, float] | None:
    mg = sum(areas.get(species, 0.0) for species in MG_ISOTOPES)
    si = sum(areas.get(species, 0.0) for species in SI_ISOTOPES)
    fe = sum(areas.get(species, 0.0) for species in FE_ISOTOPES)
    if mg <= 0.0 or si <= 0.0 or fe <= 0.0:
        return None
    total = mg + si + fe
    if total <= 0.0:
        return None
    return float(mg / total), float(si / total), float(fe / total)


def _ternary_with_rsf(
    areas: Mapping[str, float],
    rsf: Mapping[str, float],
) -> tuple[float, float, float] | None:
    mg = sum(areas.get(species, 0.0) for species in MG_ISOTOPES)
    si = sum(areas.get(species, 0.0) for species in SI_ISOTOPES)
    fe = sum(areas.get(species, 0.0) for species in FE_ISOTOPES)
    mg_factor = float(rsf.get("Mg", 1.0) or 1.0)
    si_factor = float(rsf.get("Si", 1.0) or 1.0)
    fe_factor = float(rsf.get("Fe", 1.0) or 1.0)
    if mg_factor <= 0.0 or si_factor <= 0.0 or fe_factor <= 0.0:
        return None
    mg /= mg_factor
    si /= si_factor
    fe /= fe_factor
    if mg <= 0.0 or si <= 0.0 or fe <= 0.0:
        return None
    total = mg + si + fe
    if total <= 0.0:
        return None
    return float(mg / total), float(si / total), float(fe / total)


def _ternary_to_cartesian(point: tuple[float, float, float]) -> tuple[float, float]:
    mg, si, fe = point
    x = 0.5 * (2 * si + fe)
    y = (math.sqrt(3) / 2.0) * fe
    return float(x), float(y)


def _format_stats(values: Sequence[float]) -> dict[str, float | int]:
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return {"count": 0}
    return {
        "count": int(arr.size),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr, ddof=0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def _format_numeric(value: float | None, precision: int = 2, *, scale: float = 1.0) -> str:
    if value is None:
        return "—"
    try:
        numeric = float(value)
    except Exception:
        return "—"
    if not np.isfinite(numeric):
        return "—"
    numeric *= scale
    return f"{numeric:.{precision}f}"


_SPEED_RANGE_PATTERN = re.compile(r"\d+(?:\.\d+)?")


def _speed_label_to_value(label: str | None) -> float | None:
    if not label:
        return None
    matches = _SPEED_RANGE_PATTERN.findall(label)
    if not matches:
        return None
    try:
        numbers = [float(match) for match in matches]
    except ValueError:
        return None
    if not numbers:
        return None
    return float(sum(numbers) / len(numbers))


def _material_label(record: EventRecord) -> str:
    for candidate in (record.calibration_material, record.dust_type):
        if not candidate:
            continue
        label = str(candidate).strip()
        if label:
            return label
    return "Unknown"


@dataclass
class FlightCalibrationAnalyzer:
    schedule: Sequence[DustScheduleEntry]
    timestamp_tolerance_ms: float = 2_000.0
    material_filter: frozenset[str] | None = None
    match_finder: AcceleratorMatchFinder | None = None
    calibration_matrix: CalibrationMatrix | None = None
    events: list[EventRecord] = field(default_factory=list)
    skipped_files: list[Path] = field(default_factory=list)
    missing_hdf: list[Path] = field(default_factory=list)
    applied_timezone_offset_ms: float = 0.0
    _active_timezone_offset_ms: float = field(init=False, default=0.0, repr=False)

    def __post_init__(self) -> None:
        if self.material_filter:
            normalised = frozenset(
                material.strip().lower()
                for material in self.material_filter
                if material and material.strip()
            )
            self.material_filter = normalised or None

    def classify_timestamp(self, timestamp: datetime) -> DustScheduleEntry | None:
        for entry in self.schedule:
            if entry.contains(timestamp):
                return entry
        return None

    def collect(self, data_root: Path) -> None:
        files = _discover_data_files(data_root)

        base_events, base_skipped, base_missing = self._collect_with_offset(
            files, 0.0
        )

        best_events = base_events
        best_skipped = base_skipped
        best_missing = base_missing
        best_offset = 0.0
        best_count = len(base_events)

        if best_count < 100:
            timezone_offsets_ms = []
            for hours in (6, 7):
                delta = float(hours * 3_600_000)
                timezone_offsets_ms.extend((delta, -delta))

            for offset_ms in timezone_offsets_ms:
                events, skipped, missing = self._collect_with_offset(files, offset_ms)
                if len(events) > best_count:
                    best_events = events
                    best_skipped = skipped
                    best_missing = missing
                    best_count = len(events)
                    best_offset = offset_ms

        def _event_sort_key(record: EventRecord) -> tuple[str, str, float, str]:
            campaign = (record.calibration_campaign or "").strip().lower()
            material = _material_label(record).strip().lower()
            accel_time = record.accelerator_timestamp_ms
            if accel_time is None or not math.isfinite(accel_time):
                accel_time = record.timestamp.timestamp() * 1000.0
            return (campaign, material, float(accel_time), record.event_id)

        best_events.sort(key=_event_sort_key)

        self.events = best_events
        self.skipped_files = best_skipped
        self.missing_hdf = best_missing
        self.applied_timezone_offset_ms = best_offset

    def _collect_with_offset(
        self, files: Sequence[Path], offset_ms: float
    ) -> tuple[list[EventRecord], list[Path], list[Path]]:
        events: list[EventRecord] = []
        skipped_files: list[Path] = []
        missing_hdf: list[Path] = []
        previous_offset = self._active_timezone_offset_ms
        self._active_timezone_offset_ms = offset_ms
        try:
            for path in files:
                timestamp = _parse_filename_timestamp(path)
                if timestamp is None:
                    continue
                schedule_entry = self.classify_timestamp(timestamp)
                if self.material_filter is not None:
                    if schedule_entry is None:
                        continue
                    material_key = schedule_entry.material.strip().lower()
                    if material_key not in self.material_filter:
                        continue
                hdf_path: Path | None
                if path.suffix.lower() == ".h5":
                    hdf_path = path
                    if not hdf_path.exists():
                        missing_hdf.append(hdf_path)
                        continue
                else:
                    hdf_path = self._derive_hdf_path(path)
                    if not hdf_path.exists():
                        if not self._decode_raw_file(path, hdf_path):
                            missing_hdf.append(hdf_path)
                            continue
                try:
                    self._process_file(
                        hdf_path,
                        path,
                        timestamp,
                        schedule_entry.material if schedule_entry else "Unknown",
                        schedule_entry.instrument_model if schedule_entry else "Unknown",
                        events,
                    )
                except Exception:
                    skipped_files.append(hdf_path)
        finally:
            self._active_timezone_offset_ms = previous_offset

        return events, skipped_files, missing_hdf

    def _derive_hdf_path(self, raw_path: Path) -> Path:
        parts = list(raw_path.parts)
        for idx, part in enumerate(parts):
            if part.lower() == "data":
                parts[idx] = "HDF5"
        candidate = Path(parts[0])
        for segment in parts[1:]:
            candidate /= segment
        if candidate.suffix.lower() != ".h5":
            candidate = candidate.with_suffix(".h5")
        return candidate

    def _decode_raw_file(self, raw_path: Path, target_hdf: Path) -> bool:
        """Decode a raw ``ois_output`` capture into an HDF5 analysis file."""

        try:
            from .idex_packet import IDEXEvent  # type: ignore[import]
        except Exception:
            return False

        try:
            packets = IDEXEvent(str(raw_path))
        except Exception:
            return False

        try:
            packets.write_to_hdf5(packets.data, str(raw_path))
        except Exception:
            return False

        return target_hdf.exists()

    def _process_file(
        self,
        hdf_path: Path,
        raw_path: Path,
        timestamp: datetime,
        material: str,
        instrument_model: str,
        collector: list[EventRecord],
    ) -> None:
        import h5py

        with h5py.File(hdf_path, "r") as handle:
            for event_id in sorted(handle.keys()):
                event_group = handle.get(event_id)
                if event_group is None:
                    continue
                record = self._process_event(
                    event_group,
                    raw_path,
                    event_id,
                    timestamp,
                    material,
                    instrument_model,
                )
                if record is not None:
                    collector.append(record)

    def _match_from_group(self, accel_group) -> AcceleratorMatch | None:
        if accel_group is None:
            return None
        mass_kg = _read_scalar(accel_group.get("MassKilograms"))
        charge_c = _read_scalar(accel_group.get("ChargeCoulombs"))
        radius_m = _read_scalar(accel_group.get("RadiusMeters"))
        timestamp_ms = _read_scalar(accel_group.get("IntegerTimestamp"))
        estimate_quality = _read_scalar(accel_group.get("EstimateQuality"))
        if (
            estimate_quality is None
            or mass_kg is None
            or charge_c is None
            or radius_m is None
            or timestamp_ms is None
        ):
            return None
        if estimate_quality < 3 or mass_kg <= 0.0 or charge_c <= 0.0 or radius_m <= 0.0:
            return None
        velocity_mps = _read_scalar(accel_group.get("VelocityMetersPerSecond"))
        record_id_value = _read_scalar(accel_group.get("RecordID"))
        record_id = None
        if record_id_value is not None and np.isfinite(record_id_value):
            try:
                record_id = int(record_id_value)
            except Exception:
                record_id = None
        experiment_name = _read_text(accel_group.get("ExperimentTag"))
        experiment_description = _read_text(accel_group.get("ExperimentDescription"))
        if not experiment_name:
            experiment_name = _read_text(accel_group.get("ExperimentSettingsKey"))
        dust_type = _read_text(accel_group.get("DustSourceNotes"))
        if not dust_type:
            dust_type = _read_text(accel_group.get("DustTypeID"))
        group_id_text = _read_text(accel_group.get("ExperimentSettingsID"))
        if not group_id_text:
            group_id_value = _read_scalar(accel_group.get("ExperimentSettingsID"))
            if group_id_value is not None and np.isfinite(group_id_value):
                group_id_text = str(int(group_id_value))
        velocity_value = float(velocity_mps) if velocity_mps is not None and np.isfinite(velocity_mps) else float("nan")
        return AcceleratorMatch(
            record_id=record_id,
            timestamp_ms=float(timestamp_ms),
            mass_kg=float(mass_kg),
            velocity_mps=velocity_value,
            charge_c=float(charge_c),
            radius_m=float(radius_m),
            estimate_quality=float(estimate_quality),
            source="hdf",
            experiment_name=experiment_name,
            experiment_description=experiment_description,
            dust_type=dust_type,
            group_id=group_id_text,
        )

    def _resolve_calibration_entry(
        self, match: AcceleratorMatch, inferred_name: str | None
    ) -> CalibrationEntry | None:
        if self.calibration_matrix is None:
            return None
        candidates = []
        if inferred_name:
            candidates.append(inferred_name.strip())
        if match.experiment_name:
            candidates.append(match.experiment_name.strip())
        for candidate in candidates:
            if not candidate:
                continue
            entry = self.calibration_matrix.lookup(candidate)
            if entry is not None:
                return entry
        return None

    def _process_event(
        self,
        event_group,
        raw_file: Path,
        event_id: str,
        file_timestamp: datetime,
        material: str,
        instrument_model: str,
    ) -> EventRecord | None:
        analysis = event_group.get("Analysis")
        metadata = event_group.get("Metadata")
        if analysis is None:
            return None

        instrument_timestamp_ms: float | None = None
        if metadata is not None:
            epoch_dataset = metadata.get("Epoch")
            if epoch_dataset is not None:
                epoch_value = _read_scalar(epoch_dataset)
                if epoch_value is not None:
                    instrument_timestamp_ms = _normalise_epoch_to_ms(epoch_value)
        if instrument_timestamp_ms is None:
            return None

        channel_names = ("Ion Grid", "Target L", "Target H")
        snr_channels = channel_names + ("TOF H",)

        snr_by_channel: dict[str, float] = {}
        chi_sq_by_channel: dict[str, float] = {}
        reduced_chi_sq_by_channel: dict[str, float] = {}
        mass_resolutions: list[float] = []
        ternary_point: tuple[float, float, float] | None = None

        tof_group = analysis.get("TOF H")
        mass_line_species: set[str] = set()
        mass_line_areas: dict[str, float] = {}
        if tof_group is not None:
            mass_lines = tof_group.get("MassLines")
            if mass_lines is not None:
                try:
                    table = np.asarray(mass_lines[()])
                except Exception:
                    table = None
                if table is not None:
                    dtype_names = getattr(table, "dtype", None)
                    dtype_names = getattr(dtype_names, "names", None)
                    for row in table:
                        species = _decode_species(row["assigned_species"])
                        if not species and dtype_names and "label" in dtype_names:
                            species = _decode_species(row["label"])
                        if species:
                            mass_line_species.add(species)

                        assigned_mass = float(row["assigned_mass"])
                        if not np.isfinite(assigned_mass) or assigned_mass <= 0.0:
                            continue
                        sigma = float(row["sigma"])
                        resolution = _mass_resolution(assigned_mass, sigma)
                        if resolution is not None:
                            mass_resolutions.append(resolution)
                    areas = _mass_lines_to_areas(table)
                    if areas:
                        mass_line_areas = areas
                        ternary_point = _ternary_from_areas(areas)

        target_rise_time = _read_scalar(analysis.get("Target H 10-90 Risetime"))
        target_velocity_fit_kmps = _read_scalar(analysis.get("Target H Velocity Estimate"))
        target_velocity_rise_kmps = _read_scalar(analysis.get("Target H Velocity From Rise"))
        target_velocity_ratio_kmps = _read_scalar(analysis.get("Target H Velocity From Ratio"))

        def _kmps_to_mps(value: float | None) -> float | None:
            if value is None:
                return None
            try:
                numeric = float(value)
            except Exception:
                return None
            if not np.isfinite(numeric):
                return None
            return numeric * 1000.0

        instrument_mass: float | None = None
        instrument_velocity: float | None = None
        collection_efficiency: float | None = None

        for channel in snr_channels:
            snr_dataset = analysis.get(f"{channel} SNR")
            snr_value = _read_scalar(snr_dataset)
            if snr_value is not None:
                snr_by_channel[channel] = snr_value

        for channel in channel_names:
            chi_dataset = analysis.get(f"{channel} Chi Squared")
            chi_value = _read_scalar(chi_dataset)
            if chi_value is not None:
                chi_sq_by_channel[channel] = chi_value
            red_dataset = analysis.get(f"{channel} Reduced Chi Squared")
            red_value = _read_scalar(red_dataset)
            if red_value is not None:
                reduced_chi_sq_by_channel[channel] = red_value

            if instrument_mass is None:
                mass_dataset = analysis.get(f"{channel} Dust Mass Estimate")
                mass_value = _read_scalar(mass_dataset)
                if mass_value is not None and np.isfinite(mass_value) and mass_value > 0.0:
                    instrument_mass = float(mass_value)

            if instrument_velocity is None:
                velocity_dataset = analysis.get(f"{channel} Velocity Estimate")
                velocity_value = _read_scalar(velocity_dataset)
                converted = _kmps_to_mps(velocity_value)
                if converted is not None:
                    instrument_velocity = converted

            if collection_efficiency is None:
                ce_dataset = analysis.get(f"{channel} Collection Efficiency")
                ce_value = _read_scalar(ce_dataset)
                if ce_value is not None and np.isfinite(ce_value):
                    collection_efficiency = float(ce_value)

        target_velocity_fit_mps = _kmps_to_mps(target_velocity_fit_kmps)
        target_velocity_rise_mps = _kmps_to_mps(target_velocity_rise_kmps)
        target_velocity_ratio_mps = _kmps_to_mps(target_velocity_ratio_kmps)

        accel_group = analysis.get("AcceleratorMetadata")
        raw_timestamp_ms: float | None = None
        raw_velocity_mps: float | None = None
        raw_mass_kg: float | None = None
        raw_charge_c: float | None = None
        raw_radius_m: float | None = None
        raw_quality: float | None = None
        raw_experiment_name: str | None = None
        raw_experiment_description: str | None = None
        raw_dust_type: str | None = None

        if accel_group is not None:
            raw_timestamp_ms = _read_scalar(accel_group.get("IntegerTimestamp"))
            raw_velocity_mps = _read_scalar(accel_group.get("VelocityMetersPerSecond"))
            raw_mass_kg = _read_scalar(accel_group.get("MassKilograms"))
            raw_charge_c = _read_scalar(accel_group.get("ChargeCoulombs"))
            raw_radius_m = _read_scalar(accel_group.get("RadiusMeters"))
            raw_quality = _read_scalar(accel_group.get("EstimateQuality"))
            raw_experiment_name = _read_text(accel_group.get("ExperimentTag"))
            raw_experiment_description = _read_text(
                accel_group.get("ExperimentDescription")
            )
            if not raw_experiment_name:
                raw_experiment_name = _read_text(
                    accel_group.get("ExperimentSettingsKey")
                )
            raw_dust_type = _read_text(accel_group.get("DustSourceNotes"))
            if not raw_dust_type:
                raw_dust_type = _read_text(accel_group.get("DustTypeID"))

        match = self._match_from_group(accel_group)
        adjusted_timestamp_ms = instrument_timestamp_ms + self._active_timezone_offset_ms

        def _is_valid_match(candidate: AcceleratorMatch | None) -> bool:
            if candidate is None:
                return False
            if not math.isfinite(candidate.timestamp_ms):
                return False
            if abs(candidate.timestamp_ms - adjusted_timestamp_ms) > self.timestamp_tolerance_ms:
                return False
            if candidate.mass_kg <= 0.0 or candidate.charge_c <= 0.0 or candidate.radius_m <= 0.0:
                return False
            return True

        best_match = match if _is_valid_match(match) else None

        if best_match is None and self.match_finder is not None:
            velocity_hint = instrument_velocity
            if velocity_hint is None:
                velocity_hint = (
                    target_velocity_fit_mps
                    or target_velocity_rise_mps
                    or raw_velocity_mps
                )
            fallback = self.match_finder.find(
                instrument_timestamp_ms,
                timezone_offset_ms=self._active_timezone_offset_ms,
                velocity_mps=velocity_hint,
            )
            if _is_valid_match(fallback):
                best_match = fallback

        has_match = best_match is not None

        accelerator_timestamp_ms: float | None = None
        accelerator_mass_kg: float | None = None
        accelerator_velocity: float | None = None
        accelerator_charge_c: float | None = None
        accelerator_radius_m: float | None = None
        accelerator_source = "unmatched"
        accelerator_quality: float | None = None
        experiment_name: str | None = None
        experiment_description: str | None = None
        match_dust_type: str | None = None

        if best_match is not None:
            accelerator_timestamp_ms = float(best_match.timestamp_ms)
            accelerator_mass_kg = float(best_match.mass_kg)
            accelerator_velocity = (
                float(best_match.velocity_mps)
                if math.isfinite(best_match.velocity_mps)
                else None
            )
            accelerator_charge_c = float(best_match.charge_c)
            accelerator_radius_m = float(best_match.radius_m)
            accelerator_source = best_match.source
            accelerator_quality = float(best_match.estimate_quality)
            experiment_description = best_match.experiment_description
            match_dust_type = best_match.dust_type
            experiment_name = _infer_experiment_name(best_match)

        if experiment_name is None:
            for candidate in _experiment_name_candidates(
                raw_experiment_name, raw_experiment_description
            ):
                if candidate:
                    experiment_name = candidate
                    break

        if experiment_description is None:
            experiment_description = raw_experiment_description

        if accelerator_timestamp_ms is None and raw_timestamp_ms is not None:
            accelerator_timestamp_ms = float(raw_timestamp_ms)
        if accelerator_mass_kg is None and raw_mass_kg is not None and np.isfinite(raw_mass_kg):
            accelerator_mass_kg = float(raw_mass_kg)
        if accelerator_velocity is None and raw_velocity_mps is not None and np.isfinite(raw_velocity_mps):
            accelerator_velocity = float(raw_velocity_mps)
        if accelerator_charge_c is None and raw_charge_c is not None and np.isfinite(raw_charge_c):
            accelerator_charge_c = float(raw_charge_c)
        if accelerator_radius_m is None and raw_radius_m is not None and np.isfinite(raw_radius_m):
            accelerator_radius_m = float(raw_radius_m)
        if accelerator_quality is None and raw_quality is not None and np.isfinite(raw_quality):
            accelerator_quality = float(raw_quality)

        calibration_entry: CalibrationEntry | None = None
        if self.calibration_matrix is not None:
            candidate_names: list[str] = []
            if best_match is not None:
                candidate_names.extend(
                    _experiment_name_candidates(
                        best_match.experiment_name, best_match.experiment_description
                    )
                )
            candidate_names.extend(
                _experiment_name_candidates(
                    raw_experiment_name, raw_experiment_description
                )
            )
            seen: set[str] = set()
            for candidate in candidate_names:
                if not candidate:
                    continue
                if candidate in seen:
                    continue
                seen.add(candidate)
                entry = self.calibration_matrix.lookup(candidate)
                if entry is not None:
                    calibration_entry = entry
                    break

        calibration_material = (
            (calibration_entry.material if calibration_entry and calibration_entry.material else None)
            or match_dust_type
            or raw_dust_type
            or material
        )

        return EventRecord(
            file=raw_file,
            event_id=str(event_id),
            dust_type=material,
            instrument_model=instrument_model,
            timestamp=file_timestamp,
            accelerator_timestamp_ms=accelerator_timestamp_ms,
            accelerator_mass_kg=accelerator_mass_kg,
            accelerator_velocity_mps=accelerator_velocity,
            accelerator_charge_c=accelerator_charge_c,
            accelerator_radius_m=accelerator_radius_m,
            accelerator_source=accelerator_source,
            accelerator_estimate_quality=accelerator_quality,
            has_accelerator_match=has_match,
            accelerator_experiment_name=experiment_name,
            accelerator_experiment_description=experiment_description,
            instrument_mass_estimate_kg=instrument_mass,
            instrument_velocity_estimate_mps=instrument_velocity,
            target_rise_time_us=float(target_rise_time) if target_rise_time is not None and np.isfinite(target_rise_time) else None,
            target_velocity_fit_mps=target_velocity_fit_mps,
            target_velocity_rise_mps=target_velocity_rise_mps,
            target_velocity_ratio_mps=target_velocity_ratio_mps,
            collection_efficiency=collection_efficiency,
            mass_resolutions=mass_resolutions,
            snr_by_channel=snr_by_channel,
            chi_sq_by_channel=chi_sq_by_channel,
            reduced_chi_sq_by_channel=reduced_chi_sq_by_channel,
            ternary_point=ternary_point,
            mass_line_species=tuple(sorted(mass_line_species)),
            mass_line_areas=mass_line_areas,
            calibration_campaign=calibration_entry.campaign if calibration_entry else None,
            calibration_material=calibration_material,
            calibration_speed_range=calibration_entry.speed_range if calibration_entry else None,
            calibration_target_location=calibration_entry.target_location if calibration_entry else None,
            calibration_azimuthal_location=calibration_entry.azimuthal_location if calibration_entry else None,
            calibration_notes=calibration_entry.notes if calibration_entry else None,
            reference_voltage=calibration_entry.reference_voltage if calibration_entry else None,
            target_voltage=calibration_entry.target_voltage if calibration_entry else None,
            detector_voltage=calibration_entry.detector_voltage if calibration_entry else None,
            rejection_voltage=calibration_entry.rejection_voltage if calibration_entry else None,
        )

    def build_summary(self) -> dict[str, object]:
        by_material: dict[str, list[EventRecord]] = {}
        for record in self.events:
            label = _material_label(record)
            by_material.setdefault(label, []).append(record)

        material_stats: dict[str, dict[str, object]] = {}
        best_mass_resolution: tuple[str, float] | None = None
        best_collection_efficiency: tuple[str, float] | None = None

        for material, records in by_material.items():
            resolutions = [res for record in records for res in record.mass_resolutions]
            ce_values = [
                record.collection_efficiency
                for record in records
                if record.collection_efficiency is not None and np.isfinite(record.collection_efficiency)
            ]
            mass_ratio = [
                record.instrument_mass_estimate_kg / record.accelerator_mass_kg
                for record in records
                if record.instrument_mass_estimate_kg is not None
                and record.accelerator_mass_kg is not None
                and np.isfinite(record.instrument_mass_estimate_kg)
                and np.isfinite(record.accelerator_mass_kg)
                and record.accelerator_mass_kg > 0.0
            ]
            velocity_ratio = [
                record.instrument_velocity_estimate_mps / record.accelerator_velocity_mps
                for record in records
                if record.instrument_velocity_estimate_mps is not None
                and record.accelerator_velocity_mps is not None
                and np.isfinite(record.instrument_velocity_estimate_mps)
                and np.isfinite(record.accelerator_velocity_mps)
                and record.accelerator_velocity_mps > 0.0
            ]
            ternary_points = [
                record.ternary_point for record in records if record.ternary_point is not None
            ]

            resolution_stats = _format_stats(resolutions)
            ce_stats = _format_stats(ce_values)
            mass_ratio_stats = _format_stats(mass_ratio)
            velocity_ratio_stats = _format_stats(velocity_ratio)

            if resolution_stats.get("median") is not None:
                candidate = float(resolution_stats["median"])  # type: ignore[index]
                if best_mass_resolution is None or candidate > best_mass_resolution[1]:
                    best_mass_resolution = (material, candidate)

            if ce_stats.get("median") is not None:
                candidate = float(ce_stats["median"])  # type: ignore[index]
                if best_collection_efficiency is None or candidate > best_collection_efficiency[1]:
                    best_collection_efficiency = (material, candidate)

            material_stats[material] = {
                "event_count": len(records),
                "mass_resolution": resolution_stats,
                "collection_efficiency": ce_stats,
                "mass_ratio": mass_ratio_stats,
                "velocity_ratio": velocity_ratio_stats,
                "ternary_point_count": len(ternary_points),
            }

        matched_events = sum(1 for record in self.events if record.has_accelerator_match)

        summary = {
            "total_events": len(self.events),
            "matched_events": matched_events,
            "unmatched_events": len(self.events) - matched_events,
            "materials": material_stats,
            "best_mass_resolution": {
                "material": best_mass_resolution[0] if best_mass_resolution else None,
                "median": best_mass_resolution[1] if best_mass_resolution else None,
            },
            "best_collection_efficiency": {
                "material": best_collection_efficiency[0] if best_collection_efficiency else None,
                "median": best_collection_efficiency[1] if best_collection_efficiency else None,
            },
            "skipped_files": [str(path) for path in self.skipped_files],
            "missing_hdf": [str(path) for path in self.missing_hdf],
            "applied_timezone_offset_ms": self.applied_timezone_offset_ms,
        }
        return summary

    def write_summary(self, path: Path) -> None:
        path.write_text(json.dumps(self.build_summary(), indent=2), encoding="utf-8")

    def write_report(self, path: Path) -> None:
        _ensure_matplotlib()
        assert plt is not None and PdfPages is not None

        by_material: dict[str, list[EventRecord]] = {}
        for record in self.events:
            label = _material_label(record)
            by_material.setdefault(label, []).append(record)

        with PdfPages(path) as pdf:
            fig = plt.figure(figsize=(8.5, 11))
            fig.suptitle("Flight Calibration Summary", fontsize=18)
            ax = fig.add_subplot(111)
            ax.axis("off")

            summary = self.build_summary()
            best_res = summary["best_mass_resolution"]
            best_ce = summary["best_collection_efficiency"]

            source_counts = Counter(
                record.accelerator_source or "unknown" for record in self.events
            )
            source_text = ", ".join(
                f"{source}: {count}"
                for source, count in sorted(source_counts.items())
            )

            timezone_offset_ms = summary.get("applied_timezone_offset_ms", 0.0) or 0.0
            timezone_hours = timezone_offset_ms / 3_600_000.0

            matched_count = sum(1 for record in self.events if record.has_accelerator_match)
            unmatched_count = len(self.events) - matched_count

            lines = [
                f"Total events processed: {len(self.events)}",
                f"Events with accelerator matches: {matched_count}",
                f"Events without accelerator matches: {unmatched_count}",
                f"Files skipped (missing dependencies/errors): {len(self.skipped_files)}",
                f"Files missing decoded HDF5 products: {len(self.missing_hdf)}",
            ]
            if source_text:
                lines.append(f"Accelerator metadata sources → {source_text}")
            lines.append(
                "Timezone search offset applied to accelerator matching: "
                f"{timezone_hours:+.1f} h"
            )
            if best_res.get("material"):
                lines.append(
                    "Best median mass resolution: "
                    f"{best_res['material']} ({best_res['median']:.2f})"
                )
            if best_ce.get("material"):
                lines.append(
                    "Best median collection efficiency: "
                    f"{best_ce['material']} ({best_ce['median']:.3f})"
                )
            lines.extend(
                [
                    "",  # spacer line
                    "This summary consolidates accelerator SQL matches, CSV fallbacks,",
                    "and embedded HDF metadata. Events are sorted using the calibration",
                    "lookup tables so that campaigns share contiguous pages in the report.",
                    "The table below ranks each material by event count and lists the",
                    "median instrument/accelerator ratios for quick agreement checks.",
                ]
            )

            ax.text(0.01, 0.97, "\n".join(lines), va="top", ha="left", fontsize=11)

            def _format_median(stats: Mapping[str, object], precision: int) -> str:
                value = stats.get("median") if isinstance(stats, Mapping) else None
                if value is None:
                    return "—"
                try:
                    numeric = float(value)
                except Exception:
                    return "—"
                if not np.isfinite(numeric):
                    return "—"
                return f"{numeric:.{precision}f}"

            material_rows: list[list[str]] = []
            materials_summary = summary.get("materials", {})
            for material, stats in sorted(
                materials_summary.items(),
                key=lambda item: item[1].get("event_count", 0),
                reverse=True,
            ):
                event_count = stats.get("event_count", 0)
                mass_ratio = _format_median(stats.get("mass_ratio", {}), 2)
                velocity_ratio = _format_median(stats.get("velocity_ratio", {}), 2)
                collection_eff = _format_median(stats.get("collection_efficiency", {}), 3)
                material_rows.append(
                    [
                        material,
                        str(event_count),
                        mass_ratio,
                        velocity_ratio,
                        collection_eff,
                    ]
                )

            if material_rows:
                table = ax.table(
                    cellText=material_rows,
                    colLabels=[
                        "Material",
                        "Events",
                        "Median mass ratio",
                        "Median velocity ratio",
                        "Median collection eff.",
                    ],
                    loc="lower center",
                    cellLoc="center",
                    bbox=[0.0, 0.0, 1.0, 0.42],
                )
                table.auto_set_font_size(False)
                table.set_fontsize(9)
                table.scale(1.0, 1.2)

            pdf.savefig(fig)
            plt.close(fig)

            self._plot_data_inventory(pdf)
            self._plot_speed_efficiency_trends(pdf)
            if self.calibration_matrix is not None:
                self._plot_calibration_overview(pdf)
            self._plot_methodology_page(pdf)
            self._plot_combination_diagnostics_pages(pdf)
            self._plot_olivine_rsf_series(pdf)

            for material, records in by_material.items():
                self._plot_mass_resolution(pdf, material, records)
                self._plot_snr(pdf, material, records)
                self._plot_chi_squared(pdf, material, records)
                self._plot_collection_efficiency(pdf, material, records)
                self._plot_mass_velocity(pdf, material, records)
                self._plot_ternary(pdf, material, records)
                try:
                    self._plot_mass_line_probabilities(pdf, material, records)
                    self._plot_mass_line_yields(pdf, material, records)
                except Exception:
                    continue

    def _schedule_overview_rows(self) -> list[list[str]]:
        rows: list[list[str]] = []
        entries = sorted(self.schedule, key=lambda item: item.start)
        for entry in entries:
            count_text = f"{entry.count:,}" if entry.count else "—"
            instrument = entry.instrument_model or "Unknown"
            material = entry.material or "Unknown"
            rows.append([entry.label, instrument, material, count_text])
        return rows

    def _campaign_summary_rows(self) -> list[list[str]]:
        grouped: dict[str, dict[str, object]] = {}
        for record in self.events:
            campaign = record.calibration_campaign or "Uncatalogued campaign"
            bucket = grouped.setdefault(
                campaign,
                {
                    "materials": set(),
                    "speed_ranges": set(),
                    "targets": set(),
                    "charges": [],
                    "collection": [],
                    "count": 0,
                    "start": None,
                    "end": None,
                },
            )
            bucket["count"] = int(bucket.get("count", 0)) + 1
            materials = cast(set[str], bucket.setdefault("materials", set()))
            speed_ranges = cast(set[str], bucket.setdefault("speed_ranges", set()))
            targets = cast(set[str], bucket.setdefault("targets", set()))
            charges = cast(list[float], bucket.setdefault("charges", []))
            collection = cast(list[float], bucket.setdefault("collection", []))
            material_label = _material_label(record)
            if material_label:
                materials.add(material_label)
            if record.calibration_speed_range:
                speed_ranges.add(record.calibration_speed_range)
            if record.calibration_target_location:
                target_text = record.calibration_target_location
                if record.calibration_azimuthal_location:
                    target_text = (
                        f"{target_text} (az {record.calibration_azimuthal_location})"
                    )
                targets.add(target_text)
            charge = record.accelerator_charge_c
            if charge is not None and np.isfinite(charge) and charge > 0.0:
                charges.append(float(charge))
            ce = record.collection_efficiency
            if ce is not None and np.isfinite(ce):
                collection.append(float(ce))
            event_time = record.timestamp
            if event_time.tzinfo is None:
                event_time = event_time.replace(tzinfo=timezone.utc)
            else:
                event_time = event_time.astimezone(timezone.utc)
            start = bucket["start"]
            end = bucket["end"]
            if start is None or event_time < start:  # type: ignore[operator]
                bucket["start"] = event_time
            if end is None or event_time > end:  # type: ignore[operator]
                bucket["end"] = event_time

        rows: list[list[str]] = []
        for campaign, bucket in sorted(
            grouped.items(), key=lambda item: int(item[1]["count"]), reverse=True
        ):
            materials = sorted(bucket["materials"]) if bucket.get("materials") else []
            speed_ranges = sorted(bucket["speed_ranges"]) if bucket.get("speed_ranges") else []
            targets = sorted(bucket["targets"]) if bucket.get("targets") else []
            start = bucket.get("start")
            end = bucket.get("end")
            if isinstance(start, datetime) and isinstance(end, datetime):
                start_text = start.strftime("%Y-%m-%d %H:%M")
                end_text = end.strftime("%Y-%m-%d %H:%M")
                time_range = f"{start_text} – {end_text} UTC"
            else:
                time_range = "—"
            charge_stats = _format_stats(
                cast(Sequence[float], bucket.get("charges", []))
            )
            ce_stats = _format_stats(
                cast(Sequence[float], bucket.get("collection", []))
            )
            median_charge = _format_numeric(charge_stats.get("median"), 3, scale=1e12)
            median_ce = _format_numeric(ce_stats.get("median"), 3)
            rows.append(
                [
                    campaign,
                    ", ".join(materials) if materials else "—",
                    ", ".join(speed_ranges) if speed_ranges else "—",
                    ", ".join(targets) if targets else "—",
                    time_range,
                    str(bucket.get("count", 0)),
                    median_charge,
                    median_ce,
                ]
            )
        return rows

    def _target_speed_rows(self) -> list[list[str]]:
        grouped: dict[tuple[str, str], dict[str, object]] = {}
        for record in self.events:
            speed = record.calibration_speed_range or "Unspecified"
            target = record.calibration_target_location or "Unspecified target"
            key = (target, speed)
            bucket = grouped.setdefault(
                key,
                {"count": 0, "collection": [], "charges": []},
            )
            bucket["count"] = int(bucket.get("count", 0)) + 1
            charges = cast(list[float], bucket.setdefault("charges", []))
            collection = cast(list[float], bucket.setdefault("collection", []))
            ce = record.collection_efficiency
            if ce is not None and np.isfinite(ce):
                collection.append(float(ce))
            charge = record.accelerator_charge_c
            if charge is not None and np.isfinite(charge) and charge > 0.0:
                charges.append(float(charge))

        rows: list[list[str]] = []
        for (target, speed), bucket in sorted(
            grouped.items(), key=lambda item: (item[0][0], item[0][1])
        ):
            charge_stats = _format_stats(
                cast(Sequence[float], bucket.get("charges", []))
            )
            ce_stats = _format_stats(
                cast(Sequence[float], bucket.get("collection", []))
            )
            rows.append(
                [
                    target,
                    speed,
                    str(bucket.get("count", 0)),
                    _format_numeric(charge_stats.get("median"), 3, scale=1e12),
                    _format_numeric(ce_stats.get("median"), 3),
                ]
            )
        return rows

    def _speed_efficiency_points(self) -> list[tuple[float, float, str]]:
        points: list[tuple[float, float, str]] = []
        for record in self.events:
            ce = record.collection_efficiency
            if ce is None or not np.isfinite(ce):
                continue
            speed_value = _speed_label_to_value(record.calibration_speed_range)
            if speed_value is None or not np.isfinite(speed_value):
                continue
            location = record.calibration_target_location or "Unspecified target"
            points.append((float(speed_value), float(ce), location))
        return points

    def _plot_data_inventory(self, pdf: PdfPages) -> None:
        assert plt is not None
        schedule_rows = self._schedule_overview_rows()
        campaign_rows = self._campaign_summary_rows()
        if not schedule_rows and not campaign_rows:
            return
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle("Data inventory overview", fontsize=16)
        grid = fig.add_gridspec(2, 1, height_ratios=(0.4, 0.6), hspace=0.35)

        ax_schedule = fig.add_subplot(grid[0])
        ax_schedule.axis("off")
        ax_schedule.set_title(
            "Dust accelerator usage (lookup/IDEX_Dust_Testing.csv)",
            loc="left",
            fontsize=12,
            pad=10,
        )
        if schedule_rows:
            table = ax_schedule.table(
                cellText=schedule_rows,
                colLabels=["Window", "Instrument", "Material", "Planned shots"],
                loc="center",
                cellLoc="center",
                bbox=[0.0, 0.0, 1.0, 0.85],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1.0, 1.1)
        else:
            ax_schedule.text(
                0.02,
                0.5,
                "No dust schedule information was available.",
                va="center",
                ha="left",
                fontsize=11,
            )

        ax_campaign = fig.add_subplot(grid[1])
        ax_campaign.axis("off")
        ax_campaign.set_title(
            "Processed HDF5 coverage by calibration campaign",
            loc="left",
            fontsize=12,
            pad=10,
        )
        if campaign_rows:
            table = ax_campaign.table(
                cellText=campaign_rows,
                colLabels=[
                    "Campaign",
                    "Materials",
                    "Speed windows",
                    "Targets",
                    "UTC span",
                    "Events",
                    "Median charge (pC)",
                    "Median CE",
                ],
                loc="center",
                cellLoc="center",
                colLoc="center",
                bbox=[0.0, -0.05, 1.0, 1.05],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1.0, 1.15)
        else:
            ax_campaign.text(
                0.02,
                0.5,
                "No processed events were available to summarise.",
                va="center",
                ha="left",
                fontsize=11,
            )

        pdf.savefig(fig)
        plt.close(fig)

    def _plot_speed_efficiency_trends(self, pdf: PdfPages) -> None:
        assert plt is not None
        rows = self._target_speed_rows()
        points = self._speed_efficiency_points()
        if not rows and not points:
            return
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle("Collection efficiency trends", fontsize=16)
        grid = fig.add_gridspec(2, 1, height_ratios=(0.5, 0.5), hspace=0.35)

        ax_table = fig.add_subplot(grid[0])
        ax_table.axis("off")
        ax_table.set_title(
            "Target/speed combinations derived from the calibration matrix",
            loc="left",
            fontsize=12,
            pad=10,
        )
        if rows:
            table = ax_table.table(
                cellText=rows,
                colLabels=[
                    "Target location",
                    "Speed window",
                    "Events",
                    "Median charge (pC)",
                    "Median CE",
                ],
                loc="center",
                cellLoc="center",
                bbox=[0.0, 0.0, 1.0, 0.85],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1.0, 1.15)
        else:
            ax_table.text(
                0.02,
                0.5,
                "No calibration metadata was attached to the processed events.",
                va="center",
                ha="left",
                fontsize=11,
            )

        ax_plot = fig.add_subplot(grid[1])
        if points:
            locations = sorted({location for _, _, location in points})
            for location in locations:
                xs = [speed for speed, _, loc in points if loc == location]
                ys = [ce for _, ce, loc in points if loc == location]
                ax_plot.scatter(xs, ys, label=location, alpha=0.7, s=30)
            ax_plot.set_xlabel("Speed window centre (km/s)")
            ax_plot.set_ylabel("Collection efficiency")
            ax_plot.grid(True, alpha=0.3)
            if len(locations) <= 12:
                ax_plot.legend(loc="best", fontsize="small")
        else:
            ax_plot.axis("off")
            ax_plot.text(
                0.02,
                0.5,
                "Collection efficiency vs. speed trend could not be plotted",
                va="center",
                ha="left",
            )

        pdf.savefig(fig)
        plt.close(fig)

    def _plot_calibration_overview(self, pdf: PdfPages) -> None:
        assert plt is not None
        if self.calibration_matrix is None:
            return
        overview_lines = list(self.calibration_matrix.overview_text())
        grouped: dict[str, list[CalibrationEntry]] = {}
        for entry in self.calibration_matrix.entries.values():
            if entry.date is not None:
                key = entry.date.strftime("%Y-%m-%d")
            else:
                key = entry.campaign or "Undated campaign"
            grouped.setdefault(key, []).append(entry)

        blocks: list[str] = []
        if overview_lines:
            blocks.append(
                "Campaign objectives\n" + "\n".join(f"  • {line}" for line in overview_lines)
            )

        for date_key in sorted(grouped):
            entries = grouped[date_key]
            entries.sort(key=lambda item: (item.run_number if item.run_number is not None else 999))
            block_lines = [date_key]
            for entry in entries:
                parts: list[str] = []
                if entry.run_number is not None:
                    parts.append(f"Run {entry.run_number}")
                if entry.material:
                    parts.append(f"Material: {entry.material}")
                if entry.target_location:
                    if entry.azimuthal_location:
                        parts.append(
                            f"Target {entry.target_location} (azimuth {entry.azimuthal_location})"
                        )
                    else:
                        parts.append(f"Target {entry.target_location}")
                if entry.speed_range:
                    parts.append(f"Speed window: {entry.speed_range} km/s")
                voltage_terms: list[str] = []
                if entry.reference_voltage is not None:
                    voltage_terms.append(f"V_ref={entry.reference_voltage:.0f} V")
                if entry.target_voltage is not None:
                    voltage_terms.append(f"V_target={entry.target_voltage:.0f} V")
                if entry.detector_voltage is not None:
                    voltage_terms.append(f"V_det={entry.detector_voltage:.0f} V")
                if entry.rejection_voltage is not None:
                    voltage_terms.append(f"V_reject={entry.rejection_voltage:.0f} V")
                if voltage_terms:
                    parts.append(", ".join(voltage_terms))
                if entry.csas_used:
                    parts.append(f"CSAs: {entry.csas_used}")
                description = ", ".join(parts) if parts else "Configuration recorded"
                block_lines.append(f"  • {description}")
                if entry.notes:
                    wrapped = textwrap.fill(
                        f"Notes: {entry.notes}", width=80, subsequent_indent="      "
                    )
                    block_lines.append(f"    {wrapped}")
            blocks.append("\n".join(block_lines))

        if not blocks:
            return

        fig, ax = plt.subplots(figsize=(8.5, 11))
        fig.suptitle("Dust testing campaign summary", fontsize=16)
        ax.axis("off")
        ax.text(
            0.02,
            0.98,
            "\n\n".join(blocks),
            va="top",
            ha="left",
            fontsize=11,
            family="monospace",
        )
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_methodology_page(self, pdf: PdfPages) -> None:
        assert plt is not None
        lines = [
            "Data provenance:",
            "  • Events are matched against the accelerator SQL archive when available "
            "with quality ≥ 3.",
            "  • If fewer than ten SQL matches succeed, the offline lookup/IDEX_FM_2023.csv "
            "catalog is consulted.",
            "  • Instrument timestamps are tested against UTC as well as ±6 h and ±7 h offsets "
            "to accommodate MST/DST logging.",
            "  • All timestamps must agree within ±2 seconds (2,000 ms).",
            "",
            "Derived quantities:",
            r"  • Target rise time uses the 10/90 metric: $t_{rise} = t_{90\%} - t_{10\%}$.",
            r"  • Velocity estimates come from Target H waveform fits ($v_{fit}$) and the rise-time lookup "
            r"tables ($v_{rise}$) and are compared to accelerator speeds.",
            r"  • Target masses use the calibrated charge-to-mass curves stored in the analysis products and are "
            r"compared against accelerator shot masses.",
            r"  • Ternary abundances sum fitted mass-line areas for Mg, Si, and Fe. RSFs renormalise these areas via "
            r"$A_{RSF} = A_{meas} / \mathrm{RSF}$ prior to conversion to barycentric coordinates.",
            "",
            "Figure annotations:",
            "  • Diagnostic plots include 1:1 reference lines or linear fits to highlight agreement with accelerator data.",
            "  • Captions describe the processing steps and list the voltages recorded for each accelerator configuration.",
        ]
        fig, ax = plt.subplots(figsize=(8.5, 11))
        fig.suptitle("Methodology and calculations", fontsize=16)
        ax.axis("off")
        ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", fontsize=11)
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_combination_diagnostics_pages(self, pdf: PdfPages) -> None:
        assert plt is not None
        groups: dict[tuple[str, str, str], list[EventRecord]] = {}
        for record in self.events:
            composition = (record.calibration_material or record.dust_type or "Unknown").strip()
            speed_range = (record.calibration_speed_range or "Unknown").strip()
            target_location = (record.calibration_target_location or "Unknown").strip()
            key = (composition or "Unknown", speed_range or "Unknown", target_location or "Unknown")
            groups.setdefault(key, []).append(record)
        for key in sorted(groups):
            self._plot_combination_diagnostics(pdf, key, groups[key])

    def _plot_combination_diagnostics(
        self,
        pdf: PdfPages,
        combination: tuple[str, str, str],
        records: Sequence[EventRecord],
    ) -> None:
        assert plt is not None
        composition, speed_range, target_location = combination
        rise_points: list[tuple[float, float]] = []
        fit_points: list[tuple[float, float]] = []
        rise_velocity_points: list[tuple[float, float]] = []
        mass_points: list[tuple[float, float]] = []

        def _to_kmps(value: float | None) -> float | None:
            if value is None or not np.isfinite(value):
                return None
            return float(value) / 1000.0

        for record in records:
            accel_velocity = record.accelerator_velocity_mps
            if accel_velocity is None or not np.isfinite(accel_velocity):
                continue
            accel_kmps = _to_kmps(accel_velocity)
            if accel_kmps is None:
                continue
            if record.target_rise_time_us is not None and np.isfinite(record.target_rise_time_us):
                rise_points.append((accel_kmps, float(record.target_rise_time_us)))
            if record.target_velocity_fit_mps is not None and np.isfinite(record.target_velocity_fit_mps):
                fit_points.append((accel_kmps, float(record.target_velocity_fit_mps) / 1000.0))
            if record.target_velocity_rise_mps is not None and np.isfinite(record.target_velocity_rise_mps):
                rise_velocity_points.append(
                    (accel_kmps, float(record.target_velocity_rise_mps) / 1000.0)
                )
            if (
                record.instrument_mass_estimate_kg is not None
                and np.isfinite(record.instrument_mass_estimate_kg)
                and record.accelerator_mass_kg is not None
                and np.isfinite(record.accelerator_mass_kg)
                and record.accelerator_mass_kg > 0.0
                and record.instrument_mass_estimate_kg > 0.0
            ):
                mass_points.append(
                    (float(record.accelerator_mass_kg), float(record.instrument_mass_estimate_kg))
                )

        if not (rise_points or fit_points or rise_velocity_points or mass_points):
            return

        fig, axes = plt.subplots(3, 1, figsize=(8.5, 11), constrained_layout=True)
        fig.suptitle(
            f"{composition} — {speed_range} — Target {target_location} (n={len(records)})",
            fontsize=16,
        )

        ax_rise, ax_speed, ax_mass = axes

        if rise_points:
            rise_arr = np.asarray(rise_points, dtype=float)
            ax_rise.scatter(rise_arr[:, 0], rise_arr[:, 1], color="#0b5394", alpha=0.7)
            if rise_arr.shape[0] >= 2:
                x = rise_arr[:, 0]
                y = rise_arr[:, 1]
                coeffs = np.polyfit(x, y, deg=1)
                xs = np.linspace(float(np.min(x)), float(np.max(x)), 200)
                ax_rise.plot(xs, coeffs[0] * xs + coeffs[1], color="#990000", linestyle="--", label="Linear fit")
                ax_rise.legend(loc="best")
        else:
            ax_rise.text(0.5, 0.5, "No rise-time data", transform=ax_rise.transAxes, ha="center", va="center")
        ax_rise.set_xlabel("Accelerator speed (km/s)")
        ax_rise.set_ylabel("Target rise time (µs)")
        ax_rise.grid(True, alpha=0.3)
        ax_rise.text(0.02, 0.92, r"$t_{rise} = t_{90\%} - t_{10\%}$", transform=ax_rise.transAxes)

        plotted_any = False
        if fit_points:
            fit_arr = np.asarray(fit_points, dtype=float)
            ax_speed.scatter(fit_arr[:, 0], fit_arr[:, 1], label="$v_{fit}$", color="#1c4587", alpha=0.75)
            plotted_any = True
        if rise_velocity_points:
            rise_vel_arr = np.asarray(rise_velocity_points, dtype=float)
            ax_speed.scatter(rise_vel_arr[:, 0], rise_vel_arr[:, 1], label="$v_{rise}$", color="#6aa84f", alpha=0.75)
            plotted_any = True
        if plotted_any:
            combined = []
            if fit_points:
                combined.extend(fit_points)
            if rise_velocity_points:
                combined.extend(rise_velocity_points)
            combined_arr = np.asarray(combined, dtype=float)
            low = float(np.min(combined_arr[:, 0]))
            high = float(np.max(combined_arr[:, 0]))
            diag = np.linspace(low, high, 200)
            ax_speed.plot(diag, diag, color="#444444", linestyle=":", label="1:1 reference")
        else:
            ax_speed.text(0.5, 0.5, "No velocity estimates", transform=ax_speed.transAxes, ha="center", va="center")
        ax_speed.set_xlabel("Accelerator speed (km/s)")
        ax_speed.set_ylabel("Instrument speed (km/s)")
        ax_speed.grid(True, alpha=0.3)
        if plotted_any:
            ax_speed.legend(loc="best")

        if mass_points:
            mass_arr = np.asarray(mass_points, dtype=float)
            ax_mass.scatter(mass_arr[:, 0], mass_arr[:, 1], color="#38761d", alpha=0.75)
            ax_mass.set_xscale("log")
            ax_mass.set_yscale("log")
            diag_min = float(min(np.min(mass_arr[:, 0]), np.min(mass_arr[:, 1])))
            diag_max = float(max(np.max(mass_arr[:, 0]), np.max(mass_arr[:, 1])))
            diag = np.linspace(diag_min, diag_max, 200)
            ax_mass.plot(diag, diag, color="#444444", linestyle=":", label="1:1 reference")
            ax_mass.legend(loc="best")
        else:
            ax_mass.text(0.5, 0.5, "No mass estimates", transform=ax_mass.transAxes, ha="center", va="center")
        ax_mass.set_xlabel("Accelerator mass (kg)")
        ax_mass.set_ylabel("Target H mass estimate (kg)")
        ax_mass.grid(True, alpha=0.3, which="both")
        ax_mass.text(0.02, 0.92, "Masses derived from Target H charge calibration", transform=ax_mass.transAxes)

        fig.text(
            0.5,
            0.02,
            r"\textit{Figure: Rise time, velocity, and mass comparisons for the specified composition/speed/location combination.}",
            ha="center",
            fontsize=10,
        )
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_olivine_rsf_series(self, pdf: PdfPages) -> None:
        assert plt is not None
        olivine_records = [
            record
            for record in self.events
            if (record.calibration_material or record.dust_type or "").lower().find("olivine") >= 0
        ]
        if not olivine_records:
            return
        preset_points: list[tuple[str, list[tuple[float, float]]]] = []
        for label, rsf in _RSF_PRESETS:
            cartesian: list[tuple[float, float]] = []
            for record in olivine_records:
                if not record.mass_line_areas:
                    continue
                barycentric = _ternary_with_rsf(record.mass_line_areas, rsf)
                if barycentric is None:
                    continue
                cartesian.append(_ternary_to_cartesian(barycentric))
            if cartesian:
                preset_points.append((label, cartesian))
        if not preset_points:
            return

        cols = 3
        rows = math.ceil(len(preset_points) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(8.5, 11))
        axes = np.atleast_2d(axes)
        fig.suptitle("Olivine ternary scatter with RSF adjustments", fontsize=16)

        triangle = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3) / 2], [0, 0]], dtype=float)
        for idx, (label, points) in enumerate(preset_points):
            row = idx // cols
            col = idx % cols
            ax = axes[row, col]
            ax.plot(triangle[:, 0], triangle[:, 1], color="black")
            point_array = np.asarray(points, dtype=float)
            ax.scatter(point_array[:, 0], point_array[:, 1], alpha=0.7, color="#0b5394")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, math.sqrt(3) / 2 + 0.05)
            ax.set_title(label, fontsize=11)
            ax.text(0, -0.05, "Mg", ha="center")
            ax.text(1, -0.05, "Si", ha="center")
            ax.text(0.5, math.sqrt(3) / 2 + 0.03, "Fe", ha="center")
            factors = ", ".join(f"{element}: {value:.2f}" for element, value in rsf.items())
            ax.text(0.02, 0.92, f"RSF {factors}", transform=ax.transAxes, fontsize=9)

        total_axes = rows * cols
        if total_axes > len(preset_points):
            for idx in range(len(preset_points), total_axes):
                row = idx // cols
                col = idx % cols
                axes[row, col].axis("off")

        fig.text(
            0.5,
            0.02,
            r"\textit{Figure: Mg--Si--Fe distributions for olivine shots with different relative sensitivity factors applied.}",
            ha="center",
            fontsize=10,
        )
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_histogram(
        self,
        pdf: PdfPages,
        values: Sequence[float],
        title: str,
        xlabel: str,
        bins: int = 40,
        log: bool = False,
    ) -> None:
        assert plt is not None
        arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
        fig, ax = plt.subplots(figsize=(8, 6))
        if arr.size:
            ax.hist(arr, bins=bins, color="#0b5394", alpha=0.75)
            if log:
                ax.set_yscale("log")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Count")
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_mass_resolution(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        values = [value for record in records for value in record.mass_resolutions]
        if not values:
            return
        title = f"Mass resolution distribution — {material}"
        self._plot_histogram(pdf, values, title, "m/Δm")

    def _plot_snr(self, pdf: PdfPages, material: str, records: Sequence[EventRecord]) -> None:
        for channel in ("TOF H", "Ion Grid", "Target L", "Target H"):
            values = [
                record.snr_by_channel.get(channel)
                for record in records
                if channel in record.snr_by_channel
            ]
            values = [v for v in values if v is not None and np.isfinite(v)]
            if not values:
                continue
            title = f"{channel} SNR — {material}"
            self._plot_histogram(pdf, values, title, "SNR", bins=30, log=True)

    def _plot_chi_squared(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        for channel in ("Ion Grid", "Target L", "Target H"):
            chi_values = [
                record.chi_sq_by_channel.get(channel)
                for record in records
                if channel in record.chi_sq_by_channel
            ]
            chi_values = [v for v in chi_values if v is not None and np.isfinite(v)]
            if chi_values:
                title = f"{channel} χ² — {material}"
                self._plot_histogram(pdf, chi_values, title, "χ²", bins=30, log=True)

            red_values = [
                record.reduced_chi_sq_by_channel.get(channel)
                for record in records
                if channel in record.reduced_chi_sq_by_channel
            ]
            red_values = [v for v in red_values if v is not None and np.isfinite(v)]
            if red_values:
                title = f"{channel} reduced χ² — {material}"
                self._plot_histogram(pdf, red_values, title, "χ²_r", bins=30, log=True)

    def _plot_collection_efficiency(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        values = [
            record.collection_efficiency
            for record in records
            if record.collection_efficiency is not None and np.isfinite(record.collection_efficiency)
        ]
        if not values:
            return
        title = f"Collection efficiency — {material}"
        self._plot_histogram(pdf, values, title, "η", bins=30)

    def _plot_mass_velocity(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        accel_masses = [
            record.accelerator_mass_kg
            for record in records
            if record.accelerator_mass_kg is not None
            and np.isfinite(record.accelerator_mass_kg)
        ]
        instrument_masses = [
            record.instrument_mass_estimate_kg
            for record in records
            if record.instrument_mass_estimate_kg is not None
            and np.isfinite(record.instrument_mass_estimate_kg)
        ]
        accel_velocities = [
            record.accelerator_velocity_mps
            for record in records
            if record.accelerator_velocity_mps is not None
            and np.isfinite(record.accelerator_velocity_mps)
        ]
        instrument_velocities = [
            record.instrument_velocity_estimate_mps
            for record in records
            if record.instrument_velocity_estimate_mps is not None
            and np.isfinite(record.instrument_velocity_estimate_mps)
        ]

        if accel_masses:
            self._plot_histogram(
                pdf,
                accel_masses,
                f"Accelerator mass — {material}",
                "Mass (kg)",
                bins=30,
                log=True,
            )
        if instrument_masses:
            self._plot_histogram(
                pdf,
                instrument_masses,
                f"Instrument mass estimate — {material}",
                "Mass (kg)",
                bins=30,
                log=True,
            )
        if accel_velocities:
            self._plot_histogram(
                pdf,
                accel_velocities,
                f"Accelerator velocity — {material}",
                "Velocity (m/s)",
                bins=30,
            )
        if instrument_velocities:
            self._plot_histogram(
                pdf,
                instrument_velocities,
                f"Instrument velocity estimate — {material}",
                "Velocity (m/s)",
                bins=30,
            )

    def _plot_ternary(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        assert plt is not None
        points = [record.ternary_point for record in records if record.ternary_point is not None]
        if not points:
            return

        cartesian = np.array([_ternary_to_cartesian(point) for point in points])
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.set_title(f"Mg–Si–Fe ternary diagram — {material}")
        triangle = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3) / 2], [0, 0]])
        ax.plot(triangle[:, 0], triangle[:, 1], color="black")
        ax.scatter(cartesian[:, 0], cartesian[:, 1], alpha=0.6, color="#38761d")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0, -0.05, "Mg", ha="center")
        ax.text(1, -0.05, "Si", ha="center")
        ax.text(0.5, math.sqrt(3) / 2 + 0.03, "Fe", ha="center")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, math.sqrt(3) / 2 + 0.05)
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_mass_line_probabilities(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        assert plt is not None

        events: list[tuple[float, set[str]]] = []
        for record in records:
            velocity = record.accelerator_velocity_mps
            if velocity is None or not np.isfinite(velocity):
                velocity = record.instrument_velocity_estimate_mps
            if velocity is None or not np.isfinite(velocity):
                for candidate in (
                    record.target_velocity_fit_mps,
                    record.target_velocity_rise_mps,
                    record.target_velocity_ratio_mps,
                ):
                    if candidate is not None and np.isfinite(candidate):
                        velocity = float(candidate)
                        break
            if velocity is None or not np.isfinite(velocity):
                continue
            species = set(record.mass_line_species)
            events.append((float(velocity), species))

        if not events:
            return

        species_counts: Counter[str] = Counter()
        for _, species in events:
            species_counts.update(species)

        eligible_species = [name for name, count in species_counts.items() if count >= 10]
        if not eligible_species:
            return

        velocities = np.array([velocity for velocity, _ in events], dtype=float)
        if velocities.size == 0:
            return

        finite_mask = np.isfinite(velocities)
        if not np.any(finite_mask):
            return
        velocities = velocities[finite_mask]
        if velocities.size == 0:
            return

        min_velocity = float(np.min(velocities))
        max_velocity = float(np.max(velocities))

        if not np.isfinite(min_velocity) or not np.isfinite(max_velocity):
            return

        if math.isclose(min_velocity, max_velocity):
            bin_edges = np.array([min_velocity - 0.5, max_velocity + 0.5])
        else:
            base_bins = int(math.sqrt(len(events)))
            if base_bins < 1:
                base_bins = 1
            bin_count = min(20, max(1, base_bins))
            bin_edges = np.linspace(min_velocity, max_velocity, bin_count + 1)
            if np.unique(bin_edges).size <= 1:
                bin_edges = np.array([min_velocity - 0.5, max_velocity + 0.5])
        bin_count = max(1, len(bin_edges) - 1)

        counts_total = np.zeros(bin_count, dtype=float)
        counts_by_species: dict[str, np.ndarray] = {
            name: np.zeros(bin_count, dtype=float) for name in eligible_species
        }

        for velocity, species in events:
            idx = np.searchsorted(bin_edges, velocity, side="right") - 1
            if idx < 0:
                idx = 0
            elif idx >= bin_count:
                idx = bin_count - 1
            counts_total[idx] += 1.0
            for name in species:
                if name in counts_by_species:
                    counts_by_species[name][idx] += 1.0

        valid_bins = counts_total > 0
        if not np.any(valid_bins):
            return

        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        fig, ax = plt.subplots(figsize=(8, 6))
        plotted_any = False
        for name, counts in sorted(counts_by_species.items()):
            probabilities = np.divide(
                counts,
                counts_total,
                out=np.zeros_like(counts),
                where=counts_total > 0,
            )
            if not np.any(valid_bins & np.isfinite(probabilities)):
                continue
            ax.plot(
                bin_centers[valid_bins],
                probabilities[valid_bins],
                marker="o",
                label=name,
            )
            plotted_any = True

        if not plotted_any:
            plt.close(fig)
            return

        ax.set_title(
            f"Mass line appearance probability vs. impact velocity — {material}"
        )
        ax.set_xlabel("Impact velocity (m/s)")
        ax.set_ylabel("Appearance probability")
        ax.set_ylim(0.0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize="small")
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_mass_line_yields(
        self, pdf: PdfPages, material: str, records: Sequence[EventRecord]
    ) -> None:
        assert plt is not None

        entries_by_species: dict[str, list[tuple[float, float]]] = {}
        species_counts: Counter[str] = Counter()

        for record in records:
            velocity = record.accelerator_velocity_mps
            if velocity is None or not np.isfinite(velocity):
                velocity = record.instrument_velocity_estimate_mps
            if velocity is None or not np.isfinite(velocity):
                for candidate in (
                    record.target_velocity_fit_mps,
                    record.target_velocity_rise_mps,
                    record.target_velocity_ratio_mps,
                ):
                    if candidate is not None and np.isfinite(candidate):
                        velocity = float(candidate)
                        break
            if velocity is None or not np.isfinite(velocity):
                continue
            mass_kg = record.accelerator_mass_kg
            if mass_kg is None or not np.isfinite(mass_kg) or mass_kg <= 0.0:
                mass_kg = record.instrument_mass_estimate_kg
            if mass_kg is None or not np.isfinite(mass_kg) or mass_kg <= 0.0:
                continue
            for species, area in record.mass_line_areas.items():
                if not np.isfinite(area) or area <= 0.0:
                    continue
                yield_value = area / mass_kg
                if not np.isfinite(yield_value):
                    continue
                entries_by_species.setdefault(species, []).append(
                    (float(velocity), float(yield_value))
                )
                species_counts[species] += 1

        eligible_species = [
            name for name, count in species_counts.items() if count >= 10
        ]
        if not eligible_species:
            return

        fig, ax = plt.subplots(figsize=(8, 6))
        plotted_any = False
        for name in sorted(eligible_species):
            points = entries_by_species.get(name)
            if not points:
                continue
            velocities_arr = np.asarray([velocity for velocity, _ in points], dtype=float)
            yields_arr = np.asarray([yield_value for _, yield_value in points], dtype=float)
            valid = np.isfinite(velocities_arr) & np.isfinite(yields_arr)
            if not np.any(valid):
                continue
            ax.scatter(
                velocities_arr[valid],
                yields_arr[valid],
                label=name,
                alpha=0.7,
                s=30,
            )
            plotted_any = True

        if not plotted_any:
            plt.close(fig)
            return

        ax.set_title(f"Mass line ion yield vs. impact velocity — {material}")
        ax.set_xlabel("Impact velocity (m/s)")
        ax.set_ylabel("Ion yield (area / dust mass)")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize="small")
        pdf.savefig(fig)
        plt.close(fig)


def generate_flight_calibration_report(
    data_root: Path | str,
    output_dir: Path | str,
    report_name: str = "flight_calibration_report.pdf",
    schedule_path: Path | None = None,
    material_filter: Sequence[str] | None = None,
) -> Mapping[str, Path]:
    """Process IDEX flight data and produce a detailed calibration report."""

    schedule = load_dust_schedule(schedule_path)
    try:
        match_finder = AcceleratorMatchFinder()
    except Exception:
        match_finder = None
    try:
        calibration_matrix = CalibrationMatrix()
    except Exception:
        calibration_matrix = None
    analyzer = FlightCalibrationAnalyzer(
        schedule,
        material_filter=material_filter,
        match_finder=match_finder,
        calibration_matrix=calibration_matrix,
    )
    analyzer.collect(Path(data_root))

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    pdf_path = output / report_name
    summary_path = output / "flight_calibration_summary.json"

    analyzer.write_report(pdf_path)
    analyzer.write_summary(summary_path)

    return {
        "pdf": pdf_path,
        "summary": summary_path,
        "events": analyzer.events,
    }


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate the IDEX flight calibration report."
    )
    parser.add_argument(
        "data_root",
        help="Directory containing decoded HDF5 products or raw Data files.",
    )
    parser.add_argument(
        "--output-dir",
        default="flight_calibration",
        help="Directory where the report and summary JSON will be written.",
    )
    parser.add_argument(
        "--report-name",
        default="flight_calibration_report.pdf",
        help="Filename for the generated PDF report.",
    )
    parser.add_argument(
        "--schedule",
        default=None,
        help="Optional path to a custom dust testing schedule CSV file.",
    )
    parser.add_argument(
        "--material",
        dest="materials",
        action="append",
        help=(
            "Restrict analysis to schedule windows for the specified material. "
            "May be provided multiple times."
        ),
    )
    args = parser.parse_args(argv)

    schedule_path = Path(args.schedule) if args.schedule else None
    results = generate_flight_calibration_report(
        Path(args.data_root),
        Path(args.output_dir),
        report_name=args.report_name,
        schedule_path=schedule_path,
        material_filter=tuple(args.materials) if args.materials else None,
    )
    pdf_path: Path = results["pdf"]
    summary_path: Path = results["summary"]
    print(f"Report saved to {pdf_path}")
    print(f"Summary saved to {summary_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI hook
    raise SystemExit(main())
