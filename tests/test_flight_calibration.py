import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
h5py = pytest.importorskip("h5py")
pytest.importorskip("matplotlib")

from spectrumpy_flight.flight_calibration import (
    FlightCalibrationAnalyzer,
    generate_flight_calibration_report,
    load_dust_schedule,
)


def _create_mass_lines() -> np.ndarray:
    str_dtype = h5py.string_dtype(encoding="utf-8", length=120)
    extras_dtype = h5py.string_dtype(encoding="utf-8", length=2048)
    table = np.zeros(
        3,
        dtype=[
            ("id", "i4"),
            ("label", str_dtype),
            ("assigned_species", str_dtype),
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
            ("shape", str_dtype),
            ("extras", extras_dtype),
        ],
    )
    table[0] = (1, "Mg24", "24Mg", 24.0, 0.21, 0.0, 1.0, 0.0, 1.0, 24.0, 24.3, 12.0, 0.4, "emg", "{}")
    table[1] = (2, "Si28", "28Si", 28.0, 0.18, 0.0, 1.0, 0.0, 1.0, 28.0, 28.1, 8.0, 0.35, "emg", "{}")
    table[2] = (3, "Fe56", "56Fe", 56.0, 0.30, 0.0, 1.0, 0.0, 1.0, 56.0, 56.2, 5.0, 0.25, "emg", "{}")
    return table


def _write_test_hdf(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        event_group = handle.require_group("1")
        metadata = event_group.require_group("Metadata")
        epoch_ms = 1701777600123.0
        metadata.create_dataset("Epoch", data=np.array([epoch_ms]))

        analysis = event_group.require_group("Analysis")
        accel_group = analysis.require_group("AcceleratorMetadata")
        accel_group.create_dataset("EstimateQuality", data=np.array([4.0]))
        accel_group.create_dataset("MassKilograms", data=np.array([2.5e-12]))
        accel_group.create_dataset("ChargeCoulombs", data=np.array([6.5e-13]))
        accel_group.create_dataset("RadiusMeters", data=np.array([1.1e-6]))
        accel_group.create_dataset("IntegerTimestamp", data=np.array([epoch_ms + 10.0]))
        accel_group.create_dataset("VelocityMetersPerSecond", data=np.array([2250.0]))

        tof_group = analysis.require_group("TOF H")
        tof_group.create_dataset("MassLines", data=_create_mass_lines())

        for channel, snr_value in {
            "TOF H": 12.5,
            "Ion Grid": 8.4,
            "Target L": 15.0,
            "Target H": 9.6,
        }.items():
            analysis.create_dataset(f"{channel} SNR", data=np.array([snr_value]))

        analysis.create_dataset("Target L Chi Squared", data=np.array([1.2]))
        analysis.create_dataset("Target L Reduced Chi Squared", data=np.array([1.05]))
        analysis.create_dataset("Target H Chi Squared", data=np.array([1.4]))
        analysis.create_dataset("Target H Reduced Chi Squared", data=np.array([0.98]))
        analysis.create_dataset("Ion Grid Chi Squared", data=np.array([2.0]))
        analysis.create_dataset("Ion Grid Reduced Chi Squared", data=np.array([1.3]))

        analysis.create_dataset("Target L Dust Mass Estimate", data=np.array([2.7e-12]))
        analysis.create_dataset("Target L Velocity Estimate", data=np.array([2230.0]))
        analysis.create_dataset("Target L Collection Efficiency", data=np.array([0.42]))


def _resolve_reports_dir(tmp_path: Path) -> Path:
    """Return a persistent directory for calibration reports.

    Preference order:
    1. Directory provided via the SPECTRUMPY_FLIGHT_REPORT_DIR environment variable.
    2. A ``reports`` directory at the repository root.
    3. The temporary path provided by pytest (per-test sandbox).
    """

    env_path = os.getenv("SPECTRUMPY_FLIGHT_REPORT_DIR")
    candidates = []
    if env_path:
        candidates.append(Path(env_path).expanduser())

    # Repository root is one level above the tests/ directory.
    repo_root = Path(__file__).resolve().parents[1]
    candidates.append(repo_root / "reports")

    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        else:
            return candidate

    output_dir = tmp_path / "reports"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def test_generate_flight_calibration_report(tmp_path):
    data_root = tmp_path / "Data"
    data_root.mkdir()
    hdf_path = data_root / "ois_output_20231205_120000.h5"
    _write_test_hdf(hdf_path)

    output_dir = _resolve_reports_dir(tmp_path)
    results = generate_flight_calibration_report(data_root, output_dir)

    pdf_path = results["pdf"]
    summary_path = results["summary"]
    assert pdf_path.exists(), "Expected PDF report to be created."
    assert summary_path.exists(), "Expected JSON summary to be created."

    summary = json.loads(summary_path.read_text())
    assert summary["total_events"] == 1
    materials = summary["materials"]
    assert "Olivine" in materials
    olivine_stats = materials["Olivine"]
    assert olivine_stats["event_count"] == 1
    resolution_stats = olivine_stats["mass_resolution"]
    assert resolution_stats["count"] == 3
    assert resolution_stats["median"] > 0
    assert summary["best_mass_resolution"]["material"] == "Olivine"


def test_schedule_loader_defaults(tmp_path):
    custom_csv = tmp_path / "schedule.csv"
    custom_csv.write_text(
        "Date,Instrument Model,Material,Count\nDecember 2023,IDEX FM,Olivine,10\n",
        encoding="utf-8",
    )
    schedule = load_dust_schedule(custom_csv)
    assert len(schedule) == 1
    analyzer = FlightCalibrationAnalyzer(schedule)
    timestamp = datetime(2023, 12, 5, tzinfo=timezone.utc)
    entry = analyzer.classify_timestamp(timestamp)
    assert entry is not None
    assert entry.material == "Olivine"
