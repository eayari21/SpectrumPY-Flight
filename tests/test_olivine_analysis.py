import json
import os
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
h5py = pytest.importorskip("h5py")

from spectrumpy_flight import package_path
from spectrumpy_flight.olivine_metrics import EXPECTED_MASS_LINES, generate_olivine_metrics


_OLIVINE_ENV_VAR = "OLIVINE_TEST_DATA_ROOT"
_PACKAGE_ROOT = package_path()


def _h5_has_path(handle: h5py.File, path: str) -> bool:
    try:
        handle[path]
    except KeyError:
        return False
    return True


def _collect_olivine_inputs() -> list[Path]:
    """Locate decoded olivine HDF5 files for the regression test."""

    candidates: list[Path] = []
    env_value = os.environ.get(_OLIVINE_ENV_VAR)
    if env_value:
        candidates.append(Path(env_value))

    for default in (_PACKAGE_ROOT / "HDF5", _PACKAGE_ROOT / "Data"):
        candidates.append(default)

    discovered: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        if not candidate.exists():
            continue
        if candidate.is_file():
            if candidate.suffix.lower() == ".h5" and candidate.name.startswith("ois_output_"):
                resolved = candidate.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    discovered.append(candidate)
            continue

        for path in sorted(candidate.rglob("*.h5")):
            if not path.name.startswith("ois_output_"):
                continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            discovered.append(path)

    if not discovered:
        pytest.skip(
            "No decoded olivine HDF5 files were found. "
            "Populate the package HDF5 directory or set the OLIVINE_TEST_DATA_ROOT environment variable."
        )

    return discovered


def test_olivine_analysis_generates_mass_lines_and_flags(tmp_path):
    input_paths = _collect_olivine_inputs()

    events_seen = 0
    for input_path in input_paths:
        with h5py.File(input_path, "r") as handle:
            assert handle.keys(), f"Decoded file {input_path} contained no events."
            for event_id in handle.keys():
                events_seen += 1
                analysis_base = f"{event_id}/Analysis"
                assert _h5_has_path(handle, analysis_base), "Analysis group missing from event."

                flag_base = f"{analysis_base}/Flags"
                assert _h5_has_path(handle, f"{flag_base}/FailedFits"), "Failed fit flags missing."
            assert _h5_has_path(handle, f"{flag_base}/SaturatedChannels"), "Saturation flags missing."
            assert _h5_has_path(handle, f"{flag_base}/Notes"), "Notes flags missing."

            tof_mass_path = f"{analysis_base}/TOF H/MassLines"
            if _h5_has_path(handle, tof_mass_path):
                mass_lines = handle[tof_mass_path][()]
                assert mass_lines.dtype.names is not None, "Mass line dataset missing fields."
                assert {"mass", "abundance"}.issubset(mass_lines.dtype.names), (
                    "Mass line dataset missing mass or abundance information."
                )
                if mass_lines.size:
                    assert np.all(np.isfinite(mass_lines["mass"])), "Mass values should be finite when present."

            for channel in ("Ion Grid", "Target H", "Target L"):
                channel_dataset = f"{event_id}/{channel}"
                if _h5_has_path(handle, channel_dataset):
                    fit_path = f"{analysis_base}/{channel}FitParams"
                    assert _h5_has_path(handle, fit_path), f"Missing fit parameters for {channel}."
                    mass_estimate_path = f"{analysis_base}/{channel}MassEstimate"
                    charge_path = f"{analysis_base}/{channel}ImpactCharge"
                    assert _h5_has_path(handle, mass_estimate_path), f"Missing mass estimate for {channel}."
                    assert _h5_has_path(handle, charge_path), f"Missing impact charge for {channel}."

    assert events_seen > 0, "No events were processed from the decoded olivine inputs."

    report_dir = tmp_path / "metrics"
    report = generate_olivine_metrics(input_paths, report_dir)
    pdf_path = report["pdf"]
    summary_path = report["summary"]
    assert pdf_path.exists(), "Expected PDF metrics report was not created."
    assert summary_path.exists(), "Expected metrics JSON summary was not created."

    summary_data = json.loads(summary_path.read_text())
    assert summary_data["files_processed"] == len(input_paths), "Unexpected file count in summary JSON."
    assert summary_data["events_processed"] == report["events"], "Event count mismatch in summary JSON."
    assert "mass_analysis" in summary_data, "Mass analysis results missing from summary JSON."
    mass_analysis = summary_data["mass_analysis"]
    assert isinstance(mass_analysis.get("events"), list), "Mass analysis events should be a list."
    species_names = {name for name, _ in EXPECTED_MASS_LINES}
    if mass_analysis["events"]:
        first_event = mass_analysis["events"][0]
        assert "relative_abundances" in first_event, "Relative abundances missing from mass analysis event."
        assert species_names.issubset(set(first_event["relative_abundances"].keys()))
    stats = mass_analysis.get("relative_abundance_stats")
    assert isinstance(stats, dict), "Relative abundance statistics missing from mass analysis."
    assert species_names.issubset(set(stats.keys())), "Not all expected species present in abundance statistics."
    ternary_points = mass_analysis.get("ternary_points")
    assert isinstance(ternary_points, list), "Ternary composition points should be recorded as a list."
