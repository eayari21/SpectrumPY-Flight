import re
from datetime import datetime
from pathlib import Path
from typing import Optional

import pytest

np = pytest.importorskip("numpy")
h5py = pytest.importorskip("h5py")

from idex_packet import IDEXEvent
from olivine_metrics import generate_olivine_metrics


DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
START_DATE = datetime(2023, 12, 3)
END_DATE = datetime(2023, 12, 15)


def _extract_datetime_from_name(name: str) -> Optional[datetime]:
    match = re.search(r"(\d{8})", name)
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%m%d%Y")
    except ValueError:
        return None


def _olivine_files():
    for entry in sorted(DATA_DIR.glob("*")):
        if not entry.is_file():
            continue
        event_date = _extract_datetime_from_name(entry.name)
        if event_date is None:
            continue
        if START_DATE <= event_date <= END_DATE:
            yield entry


OLIVINE_FILES = list(_olivine_files())

if not OLIVINE_FILES:
    pytest.skip("No olivine data files found in Data/", allow_module_level=True)


def _h5_has_path(handle: h5py.File, path: str) -> bool:
    try:
        handle[path]
    except KeyError:
        return False
    return True


@pytest.mark.parametrize("data_path", OLIVINE_FILES)
def test_olivine_analysis_generates_mass_lines_and_flags(tmp_path, data_path):
    output_path = tmp_path / f"{data_path.name}.h5"
    packets = IDEXEvent(str(data_path))
    packets.write_to_hdf5(packets.data, str(output_path))

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
