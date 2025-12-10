from datetime import datetime, timezone

from spectrumpy_flight.idex_packet import _parse_filename_epoch


def test_parse_filename_epoch_year_first_with_time():
    filename = "idex_waveforms_20251128_074509_to_20251203_004509.txt"
    dt, has_time = _parse_filename_epoch(filename)

    assert has_time is True
    assert dt == datetime(2025, 11, 28, 7, 45, 9, tzinfo=timezone.utc)


def test_parse_filename_epoch_month_first_with_time():
    filename = "ois_output_12212023_181601"
    dt, has_time = _parse_filename_epoch(filename)

    assert has_time is True
    assert dt == datetime(2023, 12, 21, 18, 16, 1, tzinfo=timezone.utc)


def test_parse_filename_epoch_year_first_date_only():
    filename = "capture_20250214"
    dt, has_time = _parse_filename_epoch(filename)

    assert has_time is False
    assert dt == datetime(2025, 2, 14, tzinfo=timezone.utc)
