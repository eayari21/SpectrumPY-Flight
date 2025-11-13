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
    assert match.schedule_label is not None
    assert match.campaign is not None


def test_accelerator_match_finder_search_window() -> None:
    finder = AcceleratorMatchFinder()
    timestamp_ms = 1702495880712
    matches = finder.search(timestamp_ms, limit=5)
    assert matches
    assert 1 <= len(matches) <= 5
    assert all(match.schedule_label is not None for match in matches)
    assert all(match.campaign is not None for match in matches)


def test_calibration_matrix_lookup() -> None:
    matrix = CalibrationMatrix()
    entry = matrix.lookup("2023_12_13_run0")
    assert entry is not None
    assert entry.material == "Olivine"
    assert entry.target_location == "5"
    assert entry.speed_range == "2-90"
