from spectrumpy_flight.idex_packet import parse_hs_waveform
from spectrumpy_flight.json_idex_header import parse_hs_waveform as parse_hs_waveform_json


def _build_hs_waveform(blocks: int) -> str:
    # Each block: 2 pad bits + 3x10-bit samples
    chunk = "00" + f"{1:010b}" + f"{2:010b}" + f"{3:010b}"
    return chunk * blocks


def test_parse_hs_waveform_returns_nominal_tof_channel_length():
    waveform_raw = _build_hs_waveform(2731)  # 2731 * 3 = 8193 decoded values

    parsed = parse_hs_waveform(waveform_raw)

    assert len(parsed) == 8192


def test_parse_hs_waveform_json_returns_nominal_tof_channel_length():
    waveform_raw = _build_hs_waveform(2731)  # 2731 * 3 = 8193 decoded values

    parsed = parse_hs_waveform_json(waveform_raw)

    assert len(parsed) == 8192
