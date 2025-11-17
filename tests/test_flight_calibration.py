import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("matplotlib")

from spectrumpy_flight.flight_calibration import (
    FlightCalibrationAnalyzer,
    generate_flight_calibration_report,
    load_dust_schedule,
)
from spectrumpy_flight import drive_idex_packet


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def data_root(repo_root: Path) -> Path:
    return repo_root / "Data"


@pytest.fixture(scope="session")
def hdf5_root(repo_root: Path) -> Path:
    path = repo_root / "HDF5"
    path.mkdir(parents=True, exist_ok=True)
    return path


@pytest.fixture(scope="session")
def report_root(repo_root: Path) -> Path:
    path = repo_root / "reports"
    path.mkdir(parents=True, exist_ok=True)
    return path


@pytest.fixture(scope="session")
def populated_hdf5(data_root: Path, hdf5_root: Path, repo_root: Path):
    raw_files = sorted(drive_idex_packet.find_source_files([data_root]))
    missing = []
    for raw_file in raw_files:
        needs, _ = drive_idex_packet.needs_conversion(raw_file, hdf5_root)
        if needs:
            missing.append(raw_file)
    if missing:
        log_dir = repo_root / "drive_logs" / "tests"
        log_dir.mkdir(parents=True, exist_ok=True)
        argv = [
            "--inputs",
            str(data_root),
            "--out",
            str(hdf5_root),
            "--log-dir",
            str(log_dir),
            "--max-procs",
            "2",
            "--threads-per-proc",
            "1",
        ]
        drive_idex_packet.main(argv)
    hdf5_files = sorted(hdf5_root.glob("*.h5"))
    return {
        "raw_files": raw_files,
        "hdf5_dir": hdf5_root,
        "hdf5_files": hdf5_files,
    }


@pytest.fixture(scope="session")
def flight_report_summary(
    populated_hdf5,
    data_root: Path,
    report_root: Path,
):
    output_dir = report_root / "pytest_full_report"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = generate_flight_calibration_report(data_root, output_dir)
    summary = json.loads(results["summary"].read_text())
    return {"results": results, "summary": summary, "output_dir": output_dir}


def test_drive_idex_packet_generates_hdf_for_all_data(populated_hdf5):
    raw_files: list[Path] = populated_hdf5["raw_files"]
    hdf5_dir: Path = populated_hdf5["hdf5_dir"]
    assert raw_files, "Expected repository Data/ directory to provide raw captures."
    for raw_file in raw_files:
        hdf_path = hdf5_dir / f"{raw_file.name}.h5"
        assert hdf_path.exists(), f"Missing converted file for {raw_file.name}"
        with h5py.File(hdf_path, "r") as handle:
            assert handle.keys(), "Converted HDF5 file should expose at least one event group."
            first_event = next(iter(handle.keys()))
            event_group = handle[first_event]
            assert "Metadata" in event_group
            assert "Analysis" in event_group


def test_generate_flight_calibration_report(flight_report_summary, populated_hdf5):
    results = flight_report_summary["results"]
    summary = flight_report_summary["summary"]
    pdf_path: Path = results["pdf"]
    summary_path: Path = results["summary"]
    assert pdf_path.exists() and pdf_path.stat().st_size > 0
    assert summary_path.exists()
    assert summary["total_events"] == len(summary["events"])
    inventory = summary["hdf5_inventory"]
    expected = {
        f"{raw_file.name}.h5" for raw_file in populated_hdf5["raw_files"]
    }
    assert {
        Path(entry["hdf5_file"]).name for entry in inventory
    } == expected
    lookup_tables = summary["lookup_tables"]
    assert "IDEX_FM_2023.csv" in lookup_tables
    assert lookup_tables["IDEX_FM_2023.csv"].get("contents")


def test_material_filter_limits_data(
    flight_report_summary,
    data_root: Path,
    report_root: Path,
):
    summary = flight_report_summary["summary"]
    assert summary["materials"], "Expected at least one material classification."
    material_name = next(
        material
        for material, stats in summary["materials"].items()
        if stats.get("event_count", 0)
    )
    filtered_dir = report_root / "pytest_filtered_report"
    if filtered_dir.exists():
        shutil.rmtree(filtered_dir)
    filtered_dir.mkdir(parents=True, exist_ok=True)
    filtered = generate_flight_calibration_report(
        data_root,
        filtered_dir,
        material_filter=(material_name,),
    )
    filtered_summary = json.loads(filtered["summary"].read_text())
    assert set(filtered_summary["materials"].keys()) == {material_name}
    assert filtered_summary["total_events"] == summary["materials"][material_name]["event_count"]


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


def test_collect_decodes_missing_hdf(populated_hdf5, data_root: Path):
    analyzer = FlightCalibrationAnalyzer(load_dust_schedule())
    analyzer.collect(data_root)
    assert analyzer.events
    assert not analyzer.missing_hdf
