from __future__ import annotations

import math

from spectrumpy_flight.calibration_data import AcceleratorMatchFinder, CalibrationMatrix


def test_accelerator_match_finder_csv_fallback() -> None:
    finder = AcceleratorMatchFinder()
    # Known timestamp from the lookup table (UTC epoch milliseconds).
    timestamp_ms = 1702495880712
    match = finder.find(timestamp_ms, timezone_offset_ms=0.0)
    assert match is not None
    assert match.source in {"csv", "server", "hdf"}
    # Velocity should agree with the CSV metadata (within tolerance).
    assert math.isclose(match.velocity_mps, 2237.1, rel_tol=1e-3)


def test_calibration_matrix_lookup() -> None:
    matrix = CalibrationMatrix()
    entry = matrix.lookup("2023_12_13_run0")
    assert entry is not None
    assert entry.material == "Olivine"
    assert entry.target_location == "5"
    assert entry.speed_range == "2-90"
