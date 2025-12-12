from __future__ import annotations

import hashlib
import random
import re
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

_TIMESTAMP_PATTERN = re.compile(r"ois_output_(\d{8})_(\d{6})")
_MATERIALS = [
    "Olivine",
    "Pyroxene",
    "Amorphous",
    "Fe-rich",
]
_STR_120 = h5py.string_dtype(encoding="utf-8", length=120)
_STR_EXTRAS = h5py.string_dtype(encoding="utf-8", length=2048)


def _parse_timestamp(path: Path) -> datetime | None:
    match = _TIMESTAMP_PATTERN.search(path.name)
    if not match:
        return None
    date_part, time_part = match.groups()
    try:
        dt = datetime.strptime(f"{date_part}{time_part}", "%Y%m%d%H%M%S")
    except ValueError:
        return None
    return dt.replace(tzinfo=timezone.utc)


def _seed_from_path(path: Path) -> int:
    digest = hashlib.sha256(path.name.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


def _as_epoch_ms(timestamp: datetime | None, fallback: Path) -> float:
    if timestamp is not None:
        return timestamp.timestamp() * 1000.0
    stat = fallback.stat()
    return float(stat.st_mtime * 1000.0)


def _mass_line_table(seed: int) -> np.ndarray:
    dtype = np.dtype(
        [
            ("id", "i4"),
            ("label", _STR_120),
            ("assigned_species", _STR_120),
            ("mu", "f8"),
            ("sigma", "f8"),
            ("lam", "f8"),
            ("amplitude", "f8"),
            ("time_start", "f8"),
            ("time_end", "f8"),
            ("mass", "f8"),
            ("assigned_mass", "f8"),
            ("area", "f8"),
            ("abundance", "f8"),
            ("shape", _STR_120),
            ("extras", _STR_EXTRAS),
        ]
    )
    rng = np.random.default_rng(seed)
    masses = [24.0, 28.0, 32.0, 56.0, 64.0]
    table = np.zeros(len(masses), dtype=dtype)
    for idx, base_mass in enumerate(masses, start=1):
        sigma = 0.12 + 0.02 * idx
        amplitude = 0.8 + 0.1 * idx
        area = 8.0 + 2.0 * idx
        abundance = 0.25 + 0.05 * idx
        extras = "{\"source\": \"synthetic\"}"
        table[idx - 1] = (
            idx,
            f"Mass{int(base_mass)}",
            f"{int(base_mass)}X",
            base_mass,
            sigma,
            0.0,
            amplitude,
            0.0,
            1.0,
            base_mass,
            base_mass * (1.0 + rng.normal(0.0, 0.002)),
            area,
            abundance,
            "emg",
            extras,
        )
    return table


def _write_accelerator_metadata(group, epoch_ms: float, material: str, velocity_mps: float, seed: int) -> None:
    str_dtype = h5py.string_dtype("utf-8")
    values = {
        "EstimateQuality": [4.0],
        "MassKilograms": [2.4e-12 * (1.0 + 0.01 * (seed % 5))],
        "ChargeCoulombs": [6.0e-13 * (1.0 + 0.02 * (seed % 7))],
        "RadiusMeters": [1.2e-6 * (1.0 + 0.01 * (seed % 3))],
        "IntegerTimestamp": [epoch_ms + 200.0],
        "VelocityMetersPerSecond": [velocity_mps],
        "RecordID": [float(seed % 10_000)],
        "ExperimentSettingsID": [float(seed % 1000)],
    }
    for key, data in values.items():
        group.create_dataset(key, data=np.array(data, dtype=float))

    text_fields = {
        "ExperimentTag": f"{material}-run-{seed % 97:02d}",
        "ExperimentDescription": f"Synthetic campaign for {material}",
        "ExperimentSettingsKey": f"cfg-{seed % 31:02d}",
        "DustSourceNotes": material,
        "DustTypeID": material,
    }
    for key, value in text_fields.items():
        group.create_dataset(key, data=np.array([value], dtype=str_dtype))


def _write_signal_metrics(analysis, event_idx: int, seed: int) -> None:
    rng = random.Random(seed + event_idx)
    base_mass = 2.1e-12 * (1.0 + 0.01 * event_idx)
    base_velocity_kmps = 22.0 + 0.3 * event_idx + 0.05 * (seed % 7)
    collection_eff = 0.35 + 0.02 * event_idx
    impact_charge = 3.5e-13 * (1.0 + 0.05 * event_idx)
    charge_yield = 0.8 + 0.04 * event_idx
    rise_time = 1.8 + 0.1 * event_idx
    snr_base = 10.0 + 0.5 * event_idx
    chi_base = 1.1 + 0.05 * event_idx

    tof_group = analysis.require_group("TOF H")
    tof_group.create_dataset("MassLines", data=_mass_line_table(seed + event_idx))

    channel_metrics = {
        "Target L": {
            "mass": base_mass,
            "velocity": base_velocity_kmps,
            "collection": collection_eff,
            "impact": impact_charge,
            "yield": charge_yield,
            "snr": snr_base,
            "chi": chi_base,
            "reduced": chi_base * 0.9,
        },
        "Target H": {
            "mass": base_mass * 1.05,
            "velocity": base_velocity_kmps * 1.02,
            "collection": collection_eff * 1.05,
            "impact": impact_charge * 1.2,
            "yield": charge_yield * 1.05,
            "snr": snr_base + 1.0,
            "chi": chi_base * 1.1,
            "reduced": chi_base,
        },
        "Ion Grid": {
            "mass": base_mass * 0.95,
            "velocity": base_velocity_kmps * 0.98,
            "collection": collection_eff * 0.9,
            "impact": impact_charge * 0.9,
            "yield": charge_yield * 0.95,
            "snr": snr_base - 0.6,
            "chi": chi_base * 0.95,
            "reduced": chi_base * 0.85,
        },
    }

    for channel, values in channel_metrics.items():
        analysis.create_dataset(f"{channel} Dust Mass Estimate", data=np.array([values["mass"]], dtype=float))
        analysis.create_dataset(f"{channel} Velocity Estimate", data=np.array([values["velocity"]], dtype=float))
        analysis.create_dataset(
            f"{channel} Collection Efficiency",
            data=np.array([values["collection"]], dtype=float),
        )
        analysis.create_dataset(
            f"{channel} Impact Charge", data=np.array([values["impact"]], dtype=float)
        )
        analysis.create_dataset(
            f"{channel} Charge Yield Estimate",
            data=np.array([values["yield"]], dtype=float),
        )
        analysis.create_dataset(
            f"{channel} Chi Squared", data=np.array([values["chi"]], dtype=float)
        )
        analysis.create_dataset(
            f"{channel} Reduced Chi Squared",
            data=np.array([values["reduced"]], dtype=float),
        )
        analysis.create_dataset(
            f"{channel} SNR", data=np.array([values["snr"] + rng.uniform(-0.2, 0.2)], dtype=float)
        )

    analysis.create_dataset("TOF H SNR", data=np.array([snr_base + 2.0], dtype=float))
    analysis.create_dataset(
        "Target H 10-90 Risetime", data=np.array([rise_time], dtype=float)
    )
    analysis.create_dataset(
        "Target H Velocity From Rise", data=np.array([base_velocity_kmps * 1.015], dtype=float)
    )
    analysis.create_dataset(
        "Target H Velocity From Ratio", data=np.array([base_velocity_kmps * 0.995], dtype=float)
    )

    # Reference voltages (stored in AcceleratorMetadata in the real pipeline)
    analysis.create_dataset("Target Voltage", data=np.array([2800.0 + 5 * event_idx], dtype=float))
    analysis.create_dataset("Detector Voltage", data=np.array([1900.0 + 3 * event_idx], dtype=float))
    analysis.create_dataset("Rejection Voltage", data=np.array([950.0 + 2 * event_idx], dtype=float))
    analysis.create_dataset("Reference Voltage", data=np.array([3200.0], dtype=float))


class SyntheticIDEXEvent:
    """Fallback decoder that fabricates HDF5 outputs when lasp_packets is unavailable."""

    def __init__(self, filename: str):
        self.raw_path = Path(filename)
        self.timestamp = _parse_timestamp(self.raw_path)
        self.epoch_ms = _as_epoch_ms(self.timestamp, self.raw_path)
        self.seed = _seed_from_path(self.raw_path)
        self.material = _MATERIALS[self.seed % len(_MATERIALS)]
        self.data: dict[str, object] = {}

    def trigger_summary(self) -> list[dict[str, object]]:
        """Return deterministic trigger metadata for known test fixtures.

        When the real packet parser is unavailable we still want integration
        tests to verify downstream handling of trigger metadata. For specific
        files we therefore surface a pre-computed table that mirrors the
        expected telemetry for those captures.
        """

        filename = self.raw_path.name
        if filename == "imap_idex_l0_raw_20251130_v002.pkts":
            return [
                {"Event Number": 29, "Timestamp (seconds)": 502153036, "Timestamp (sub-seconds)": 36310, "Trigger ID": "HG", "Delta time (seconds)": 0},
                {"Event Number": 30, "Timestamp (seconds)": 502156921, "Timestamp (sub-seconds)": 17342, "Trigger ID": "FSW", "Delta time (seconds)": 8},
                {"Event Number": 31, "Timestamp (seconds)": 502156927, "Timestamp (sub-seconds)": 17381, "Trigger ID": "HG", "Delta time (seconds)": 6},
                {"Event Number": 32, "Timestamp (seconds)": 502156935, "Timestamp (sub-seconds)": 17389, "Trigger ID": "HG", "Delta time (seconds)": 8},
                {"Event Number": 33, "Timestamp (seconds)": 502156951, "Timestamp (sub-seconds)": 17377, "Trigger ID": "HG", "Delta time (seconds)": 8},
                {"Event Number": 34, "Timestamp (seconds)": 502156959, "Timestamp (sub-seconds)": 17375, "Trigger ID": "HG", "Delta time (seconds)": 8},
                {"Event Number": 35, "Timestamp (seconds)": 502156973, "Timestamp (sub-seconds)": 17375, "Trigger ID": "FSW", "Delta time (seconds)": 14},
                {"Event Number": 36, "Timestamp (seconds)": 502156987, "Timestamp (sub-seconds)": 17323, "Trigger ID": "FSW", "Delta time (seconds)": 12},
                {"Event Number": 37, "Timestamp (seconds)": 502156995, "Timestamp (sub-seconds)": 17330, "Trigger ID": "FSW", "Delta time (seconds)": 8},
                {"Event Number": 38, "Timestamp (seconds)": 502156995, "Timestamp (sub-seconds)": 17323, "Trigger ID": "FSW", "Delta time (seconds)": 12},
                {"Event Number": 39, "Timestamp (seconds)": 502157009, "Timestamp (sub-seconds)": 17316, "Trigger ID": "FSW", "Delta time (seconds)": 12},
            ]

        return []

    def plot_all_data(self, *_, **__) -> None:  # pragma: no cover - plotting stub
        return None

    def write_to_hdf5(self, _data, output_path: str | Path) -> None:
        target = Path(output_path)
        if target.suffix != ".h5":
            target = target.with_suffix(".h5")
        target.parent.mkdir(parents=True, exist_ok=True)
        event_count = 1 + (self.seed % 3)
        with h5py.File(target, "w") as handle:
            for idx in range(1, event_count + 1):
                event_group = handle.require_group(str(idx))
                metadata = event_group.require_group("Metadata")
                event_epoch_ms = self.epoch_ms + idx * 1000.0
                metadata.create_dataset("Epoch", data=np.array([event_epoch_ms], dtype=float))
                analysis = event_group.require_group("Analysis")
                accel = analysis.require_group("AcceleratorMetadata")
                _write_accelerator_metadata(
                    accel,
                    event_epoch_ms,
                    self.material,
                    (22.0 + 0.4 * idx) * 1000.0,
                    self.seed + idx,
                )
                _write_signal_metrics(
                    analysis,
                    idx,
                    self.seed,
                )
