import subprocess
from pathlib import Path

import h5py
import pandas as pd


def test_first_transmit_test_trigger_params(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    script_dir = repo_root / "src" / "spectrumpy_flight"

    cmd = [
        "python",
        "idex_packet.py",
        "-f",
        "Data/Flight/imap_idex_l0_raw_20251130_v002.pkts",
    ]
    result = subprocess.run(
        cmd,
        cwd=script_dir,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, f"Command failed: {result.stderr or result.stdout}"

    reports_dir = repo_root / "reports"
    actual_path = reports_dir / "first_transmit_trigger_params.csv"
    assert actual_path.exists(), "Trigger summary was not produced"

    expected_path = repo_root / "tests" / "data" / "first_transmit_expected_trigger_params.csv"
    expected_df = pd.read_csv(expected_path)
    actual_df = pd.read_csv(actual_path)

    pd.testing.assert_frame_equal(
        actual_df.reset_index(drop=True),
        expected_df.reset_index(drop=True),
        check_like=True,
    )

    comparison_report = actual_df.copy()
    for column in expected_df.columns:
        comparison_report[f"{column} matches"] = actual_df[column] == expected_df[column]

    report_path = reports_dir / "first_transmit_test_trigger_params_report.csv"
    comparison_report.to_csv(report_path, index=False)

    hdf5_path = script_dir / "Data" / "Flight" / "imap_idex_l0_raw_20251130_v002.h5"
    assert hdf5_path.exists(), "HDF5 output was not written"

    trigger_mode_map = {"HG": "HGThreshold"}

    with h5py.File(hdf5_path, "r") as handle:
        for _, row in expected_df.iterrows():
            event = str(int(row["Event Number"]))
            assert event in handle, f"Event {event} missing from HDF5"

            metadata = handle[event]["Metadata"]
            trigger_mode = metadata["TriggerMode"].asstr()[0]
            expected_mode = trigger_mode_map.get(row["Trigger ID"], row["Trigger ID"])
            assert trigger_mode == expected_mode

            assert metadata["TimestampSeconds"][0] == row["Timestamp (seconds)"]
            assert metadata["TimestampSubseconds"][0] == row["Timestamp (sub-seconds)"]
