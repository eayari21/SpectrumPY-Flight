"""Pytest configuration for the SpectrumPY-Flight test suite."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, Tuple

import pytest

SRC_DIR = Path(__file__).resolve().parent.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import spectrumpy_flight  # noqa: E402, F401


def _decode_strings(values: Iterable[object]) -> Iterable[str]:
    for value in values:
        if isinstance(value, bytes):
            text = value.decode("utf-8", errors="ignore")
        else:
            text = str(value)
        yield text.strip()


@pytest.fixture(scope="session")
def idex_mass_analysis() -> Dict[Tuple[str, int], dict]:
    pytest.importorskip("lasp_packets")
    np = pytest.importorskip("numpy")
    h5py = pytest.importorskip("h5py")

    from spectrumpy_flight import package_path
    from spectrumpy_flight.idex_packet import IDEXEvent, _resolve_output_path

    package_root = package_path()
    data_dir = package_root / "Data"
    raw_files = (
        "ois_output_12132023_223729",
        "ois_output_12182023_185430",
    )

    missing = [name for name in raw_files if not (data_dir / name).exists()]
    if missing:
        pytest.skip(
            "Required packet capture(s) missing: " + ", ".join(sorted(missing))
        )

    results: Dict[Tuple[str, int], dict] = {}
    for raw_name in raw_files:
        input_path = data_dir / raw_name
        event = IDEXEvent(str(input_path))
        event.write_to_hdf5(event.data, str(input_path))
        output_path = _resolve_output_path(str(input_path))
        try:
            with h5py.File(output_path, "r") as handle:
                for event_id in handle.keys():
                    tof_group = handle.get(f"{event_id}/Analysis/TOF H")
                    if tof_group is None or "MassLines" not in tof_group:
                        continue
                    table = tof_group["MassLines"][()]
                    species = list(_decode_strings(table["assigned_species"]))
                    abundances = np.asarray(table["abundance"], dtype=float)
                    if abundances.size == 0:
                        continue
                    mapping = {}
                    for name, value in zip(species, abundances):
                        if not name:
                            continue
                        mapping[name] = float(value * 100.0)
                    if not mapping:
                        continue
                    stretch = float(tof_group.attrs.get("MassStretch", np.nan))
                    results[(raw_name, int(event_id))] = {
                        "stretch_ns": stretch,
                        "abundances": mapping,
                    }
        finally:
            if output_path.exists():
                output_path.unlink()

    if not results:
        pytest.skip("No mass line results were generated from the sample packets.")
    return results
