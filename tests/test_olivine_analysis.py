import pytest
import json

np = pytest.importorskip("numpy")
h5py = pytest.importorskip("h5py")

from olivine_metrics import EXPECTED_MASS_LINES, generate_olivine_metrics


def _h5_has_path(handle: h5py.File, path: str) -> bool:
    try:
        handle[path]
    except KeyError:
        return False
    return True


def _build_synthetic_event(handle: h5py.File, name: str) -> None:
    event_group = handle.create_group(name)

    # Waveform datasets required by trigger interpolation helpers.
    time_low = np.linspace(0.0, 10.0, 50)
    time_high = np.linspace(0.0, 5.0, 100)
    pulse = np.exp(-((time_low - 2.0) ** 2) / 0.5)

    event_group.create_dataset("Time (low sampling)", data=time_low)
    event_group.create_dataset("Time (high sampling)", data=time_high)

    for channel, waveform in {
        "Ion Grid": pulse,
        "Target H": pulse * 0.8,
        "Target L": pulse * 1.2,
        "TOF H": np.exp(-((time_high - 1.5) ** 2) / 0.2),
    }.items():
        event_group.create_dataset(channel, data=waveform)

    analysis_group = event_group.create_group("Analysis")
    flags_group = analysis_group.create_group("Flags")
    flags_group.create_dataset("FailedFits", data=np.array([0], dtype=np.int8))
    flags_group.create_dataset("SaturatedChannels", data=np.array([0], dtype=np.int8))
    flags_group.create_dataset("Notes", data=np.array([0], dtype=np.int8))

    tof_group = analysis_group.create_group("TOF H")
    tof_group.create_dataset(
        "MassLines",
        data=np.array([(32.0, 0.75)], dtype=[("mass", "f8"), ("abundance", "f8")]),
    )

    metadata_group = event_group.create_group("Metadata")
    metadata_group.create_dataset("Ion Grid Saturated", data=np.array([0], dtype=np.int8))

    for channel in ("Ion Grid", "Target H", "Target L"):
        analysis_group.create_dataset(f"{channel}FitParams", data=np.array([0, 0, 0, 1, 5], dtype=float))
        analysis_group.create_dataset(f"{channel}MassEstimate", data=np.array([1.5], dtype=float))
        analysis_group.create_dataset(f"{channel}ImpactCharge", data=np.array([2.5], dtype=float))
        analysis_group.create_dataset(f"{channel}ChiSquared", data=np.array([1.2], dtype=float))
        analysis_group.create_dataset(f"{channel}ReducedChiSquared", data=np.array([0.8], dtype=float))
        analysis_group.create_dataset(f"{channel} SNR", data=np.array([10.0], dtype=float))

    analysis_group.create_dataset("TOF H SNR", data=np.array([12.0], dtype=float))


def test_olivine_analysis_generates_mass_lines_and_flags(tmp_path):
    output_path = tmp_path / "synthetic_event.h5"
    with h5py.File(output_path, "w") as handle:
        _build_synthetic_event(handle, "Event_0001")

    assert output_path.exists(), "Expected analysis HDF5 output was not created."

    with h5py.File(output_path, "r") as handle:
        for event_id in handle.keys():
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

    report_dir = tmp_path / "metrics"
    report = generate_olivine_metrics([output_path], report_dir)
    pdf_path = report["pdf"]
    summary_path = report["summary"]
    assert pdf_path.exists(), "Expected PDF metrics report was not created."
    assert summary_path.exists(), "Expected metrics JSON summary was not created."

    summary_data = json.loads(summary_path.read_text())
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
