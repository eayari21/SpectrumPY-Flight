import numpy as np
import pytest

from spectrumpy_flight.idex_packet import BaseIDEXEvent


def _build_event_with_header(header):
    event = object.__new__(BaseIDEXEvent)
    event.header = header
    event.hspretrigblocks = 0
    event.lspretrigblocks = 0
    event.hgdelay = 0
    event.mgdelay = 0
    event.lgdelay = 0
    return event


def test_threshold_trigger_aligns_to_zero_time():
    event_id = 1
    sample_delay = 12
    header = {
        (event_id, "HSPretriggerBlocks"): 0,
        (event_id, "TOFDelay_H"): sample_delay,
        (event_id, "TriggerOrigin"): "HG",
        (event_id, "TriggerMode"): "HGThreshold",
        (event_id, "TriggerLevel"): 0.5,
    }
    event = _build_event_with_header(header)
    trigger_origin = event._trigger_origin(event_id)
    trigger_mode = event.header[(event_id, "TriggerMode")]

    if trigger_origin == "HG":
        channel = "TOF H"
    elif trigger_origin == "MG":
        channel = "TOF M"
    else:
        channel = "TOF L"

    sample_count = 1024
    time_values = event._build_time_array(
        sample_count, high_rate=True, event_id=event_id, channel=channel
    )
    t0_index = 512 - 1
    trigger_index = t0_index + sample_delay
    waveform = np.zeros(sample_count, dtype=float)
    waveform[trigger_index:] = 0.6

    if trigger_mode.endswith("Threshold"):
        threshold = event.header[(event_id, "TriggerLevel")]
        crossing_index = int(np.argmax(waveform >= threshold))
        assert crossing_index == trigger_index
        assert time_values[crossing_index] == pytest.approx(0.0)
