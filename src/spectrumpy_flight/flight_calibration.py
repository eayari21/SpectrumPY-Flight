"""Flight calibration analysis and reporting utilities."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Mapping, Sequence

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
from .olivine_metrics import FE_ISOTOPES, MG_ISOTOPES, SI_ISOTOPES

__all__ = [
    "DustScheduleEntry",
    "load_dust_schedule",
    "FlightCalibrationAnalyzer",
    "generate_flight_calibration_report",
]


if sys.version_info >= (3, 10):

    @dataclass(slots=True)
    class DustScheduleEntry:
        """Describe a dust accelerator campaign window."""

        label: str
        start: datetime
        end: datetime
        instrument_model: str
        material: str
        count: int

        def contains(self, timestamp: datetime) -> bool:
            return self.start <= timestamp <= self.end

else:

    @dataclass
    class DustScheduleEntry:
        """Describe a dust accelerator campaign window."""

        __slots__ = (
            "label",
            "start",
            "end",
            "instrument_model",
            "material",
            "count",
        )

        label: str
        start: datetime
        end: datetime
        instrument_model: str
        material: str
        count: int

        def contains(self, timestamp: datetime) -> bool:
            return self.start <= timestamp <= self.end


def _parse_day(text: str) -> int:
    digits = "".join(ch for ch in text if ch.isdigit())
    if not digits:
        raise ValueError(f"Could not parse day from '{text}'.")
    return int(digits)


def _parse_date_range(text: str) -> tuple[datetime, datetime]:
    cleaned = text.strip().replace(",", "")
    parts = cleaned.split()
    if not parts:
        raise ValueError("Empty date string in dust schedule.")

    if len(parts) == 2:
        month_name, year_text = parts
        month = datetime.strptime(month_name, "%B").month
        year = int(year_text)
        start = datetime(year, month, 1, tzinfo=timezone.utc)
        if month == 12:
            end_month = 1
            end_year = year + 1
        else:
            end_month = month + 1
            end_year = year
        end = datetime(end_year, end_month, 1, tzinfo=timezone.utc) - timedelta(seconds=1)
        return start, end

    if len(parts) == 3 and "-" in parts[1]:
        month_name, day_range, year_text = parts
        start_day_text, end_day_text = day_range.split("-", 1)
        month = datetime.strptime(month_name, "%B").month
        year = int(year_text)
        start_day = _parse_day(start_day_text)
        end_day = _parse_day(end_day_text)
        start = datetime(year, month, start_day, tzinfo=timezone.utc)
        end = datetime(year, month, end_day, 23, 59, 59, tzinfo=timezone.utc)
        return start, end

    raise ValueError(f"Unsupported dust schedule date format: '{text}'.")


def load_dust_schedule(path: Path | None = None) -> list[DustScheduleEntry]:
    """Load the dust testing schedule from a CSV file."""

    if path is None:
        path = package_path("lookup", "IDEX_Dust_Testing.csv")

    entries: list[DustScheduleEntry] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            date_label = (row.get("Date") or "").strip()
            instrument_model = (row.get("Instrument Model") or "").strip()
            material = (row.get("Material") or "").strip()
            count_text = (row.get("Count") or "0").strip()
            if not date_label:
                continue
            start, end = _parse_date_range(date_label)
            try:
                count = int(count_text)
            except Exception:
                count = 0
            entries.append(
                DustScheduleEntry(
                    label=date_label,
                    start=start,
                    end=end,
                    instrument_model=instrument_model,
                    material=material,
                    count=count,
                )
            )
    entries.sort(key=lambda entry: entry.start)
    return entries


def _discover_data_files(root: Path) -> list[Path]:
    candidates: list[Path] = []
    for path in sorted(root.rglob("ois_output_*")):
        if path.is_dir():
            continue
        candidates.append(path)
    return candidates


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


def _ensure_matplotlib() -> None:
    if plt is None or PdfPages is None:  # pragma: no cover - runtime guard
        raise RuntimeError(
            "Matplotlib is required to generate the flight calibration report."
        )


@dataclass
class EventRecord:
    file: Path
    event_id: str
    dust_type: str
    instrument_model: str
    timestamp: datetime
    accelerator_mass_kg: float
    accelerator_velocity_mps: float | None
    accelerator_charge_c: float
    accelerator_radius_m: float
    instrument_mass_estimate_kg: float | None
    instrument_velocity_estimate_mps: float | None
    collection_efficiency: float | None
    mass_resolutions: list[float]
    snr_by_channel: dict[str, float]
    chi_sq_by_channel: dict[str, float]
    reduced_chi_sq_by_channel: dict[str, float]
    ternary_point: tuple[float, float, float] | None
    mass_line_species: tuple[str, ...]


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


@dataclass
class FlightCalibrationAnalyzer:
    schedule: Sequence[DustScheduleEntry]
    timestamp_tolerance_ms: float = 2_000.0
    events: list[EventRecord] = field(default_factory=list)
    skipped_files: list[Path] = field(default_factory=list)
    missing_hdf: list[Path] = field(default_factory=list)

    def classify_timestamp(self, timestamp: datetime) -> DustScheduleEntry | None:
        for entry in self.schedule:
            if entry.contains(timestamp):
                return entry
        return None

    def collect(self, data_root: Path) -> None:
        files = _discover_data_files(data_root)
        for path in files:
            timestamp = _parse_filename_timestamp(path)
            if timestamp is None:
                continue
            schedule_entry = self.classify_timestamp(timestamp)
            hdf_path = path if path.suffix.lower() == ".h5" else self._derive_hdf_path(path)
            if not hdf_path.exists():
                self.missing_hdf.append(hdf_path)
                continue
            try:
                self._process_file(
                    hdf_path,
                    path,
                    timestamp,
                    schedule_entry.material if schedule_entry else "Unknown",
                    schedule_entry.instrument_model if schedule_entry else "Unknown",
                )
            except Exception:
                self.skipped_files.append(hdf_path)

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

    def _process_file(
        self,
        hdf_path: Path,
        raw_path: Path,
        timestamp: datetime,
        material: str,
        instrument_model: str,
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
                    self.events.append(record)

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

        accel_group = analysis.get("AcceleratorMetadata") if analysis else None
        if accel_group is None:
            return None

        estimate_quality = _read_scalar(accel_group.get("EstimateQuality"))
        mass_kg = _read_scalar(accel_group.get("MassKilograms"))
        charge_c = _read_scalar(accel_group.get("ChargeCoulombs"))
        radius_m = _read_scalar(accel_group.get("RadiusMeters"))
        accel_timestamp_ms = _read_scalar(accel_group.get("IntegerTimestamp"))
        velocity_mps = _read_scalar(accel_group.get("VelocityMetersPerSecond"))

        if (
            estimate_quality is None
            or mass_kg is None
            or charge_c is None
            or radius_m is None
            or accel_timestamp_ms is None
        ):
            return None
        if estimate_quality < 3 or mass_kg <= 0.0 or charge_c <= 0.0 or radius_m <= 0.0:
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
        if abs(accel_timestamp_ms - instrument_timestamp_ms) > self.timestamp_tolerance_ms:
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
                    ternary_point = _ternary_from_areas(areas)

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
                    instrument_mass = mass_value

            if instrument_velocity is None:
                velocity_dataset = analysis.get(f"{channel} Velocity Estimate")
                velocity_value = _read_scalar(velocity_dataset)
                if velocity_value is not None and np.isfinite(velocity_value):
                    instrument_velocity = velocity_value

            if collection_efficiency is None:
                ce_dataset = analysis.get(f"{channel} Collection Efficiency")
                ce_value = _read_scalar(ce_dataset)
                if ce_value is not None and np.isfinite(ce_value):
                    collection_efficiency = ce_value

        return EventRecord(
            file=raw_file,
            event_id=str(event_id),
            dust_type=material,
            instrument_model=instrument_model,
            timestamp=file_timestamp,
            accelerator_mass_kg=float(mass_kg),
            accelerator_velocity_mps=None if velocity_mps is None else float(velocity_mps),
            accelerator_charge_c=float(charge_c),
            accelerator_radius_m=float(radius_m),
            instrument_mass_estimate_kg=instrument_mass,
            instrument_velocity_estimate_mps=instrument_velocity,
            collection_efficiency=collection_efficiency,
            mass_resolutions=mass_resolutions,
            snr_by_channel=snr_by_channel,
            chi_sq_by_channel=chi_sq_by_channel,
            reduced_chi_sq_by_channel=reduced_chi_sq_by_channel,
            ternary_point=ternary_point,
            mass_line_species=tuple(sorted(mass_line_species)),
        )

    def build_summary(self) -> dict[str, object]:
        by_material: dict[str, list[EventRecord]] = {}
        for record in self.events:
            by_material.setdefault(record.dust_type, []).append(record)

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
                if record.instrument_mass_estimate_kg
                and record.accelerator_mass_kg
                and np.isfinite(record.instrument_mass_estimate_kg)
                and np.isfinite(record.accelerator_mass_kg)
                and record.accelerator_mass_kg > 0.0
            ]
            velocity_ratio = [
                record.instrument_velocity_estimate_mps / record.accelerator_velocity_mps
                for record in records
                if record.instrument_velocity_estimate_mps
                and record.accelerator_velocity_mps
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

        summary = {
            "total_events": len(self.events),
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
        }
        return summary

    def write_summary(self, path: Path) -> None:
        path.write_text(json.dumps(self.build_summary(), indent=2), encoding="utf-8")

    def write_report(self, path: Path) -> None:
        _ensure_matplotlib()
        assert plt is not None and PdfPages is not None

        by_material: dict[str, list[EventRecord]] = {}
        for record in self.events:
            by_material.setdefault(record.dust_type, []).append(record)

        with PdfPages(path) as pdf:
            fig = plt.figure(figsize=(8.5, 11))
            fig.suptitle("Flight Calibration Summary", fontsize=18)
            ax = fig.add_subplot(111)
            ax.axis("off")

            lines = [
                f"Total events with accelerator matches: {len(self.events)}",
                f"Files skipped (missing dependencies/errors): {len(self.skipped_files)}",
                f"Files missing decoded HDF5 products: {len(self.missing_hdf)}",
            ]
            summary = self.build_summary()
            best_res = summary["best_mass_resolution"]
            best_ce = summary["best_collection_efficiency"]
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
            ax.text(0.01, 0.95, "\n".join(lines), va="top", ha="left", fontsize=12)
            pdf.savefig(fig)
            plt.close(fig)

            for material, records in by_material.items():
                self._plot_mass_resolution(pdf, material, records)
                self._plot_snr(pdf, material, records)
                self._plot_chi_squared(pdf, material, records)
                self._plot_collection_efficiency(pdf, material, records)
                self._plot_mass_velocity(pdf, material, records)
                self._plot_ternary(pdf, material, records)
                try:
                    self._plot_mass_line_probabilities(pdf, material, records)
                except Exception:
                    continue

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
        accel_masses = [record.accelerator_mass_kg for record in records]
        instrument_masses = [
            record.instrument_mass_estimate_kg
            for record in records
            if record.instrument_mass_estimate_kg is not None
        ]
        accel_velocities = [
            record.accelerator_velocity_mps
            for record in records
            if record.accelerator_velocity_mps is not None
        ]
        instrument_velocities = [
            record.instrument_velocity_estimate_mps
            for record in records
            if record.instrument_velocity_estimate_mps is not None
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


def generate_flight_calibration_report(
    data_root: Path | str,
    output_dir: Path | str,
    report_name: str = "flight_calibration_report.pdf",
    schedule_path: Path | None = None,
) -> Mapping[str, Path]:
    """Process IDEX flight data and produce a detailed calibration report."""

    schedule = load_dust_schedule(schedule_path)
    analyzer = FlightCalibrationAnalyzer(schedule)
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
    args = parser.parse_args(argv)

    schedule_path = Path(args.schedule) if args.schedule else None
    results = generate_flight_calibration_report(
        Path(args.data_root),
        Path(args.output_dir),
        report_name=args.report_name,
        schedule_path=schedule_path,
    )
    pdf_path: Path = results["pdf"]
    summary_path: Path = results["summary"]
    print(f"Report saved to {pdf_path}")
    print(f"Summary saved to {summary_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI hook
    raise SystemExit(main())
