"""Authority-backed descriptions for IDEX TXHDR packed/unpacked metadata fields."""

from __future__ import annotations

from typing import Optional

TXHDR_AUTHORITY_REFERENCE = (
    "LASP Document 168606 (SUDA FSW Science Packets White Paper, Rev E)"
)


PACKED_FIELD_DESCRIPTIONS = {
    "IDX__TXHDRSCIEVTLEN": "Science event packet length word.",
    "IDX__TXHDRTIMESEC1": "Upper spacecraft-clock seconds word (most-significant 16 bits).",
    "IDX__TXHDRTIMESEC2": "Lower spacecraft-clock seconds word (least-significant 16 bits).",
    "IDX__TXHDRTIMESUBS": "Sub-second count in 20 microsecond ticks.",
    "IDX__TXHDRTRIGOFFSET": "Trigger offset between trigger and stored waveform window.",
    "IDX__TXHDRTRIGID": "Trigger identity / trigger source register.",
    "IDX__TXHDREVTNUM": "Event number and trigger-origin bitfield.",
    "IDX__TXHDRBLOCKS": "Packed pre/post-trigger block counts for HS and LS waveform streams.",
    "IDX__TXHDRHGTRIGCTRL1": "High-gain trigger control register 1.",
    "IDX__TXHDRHGTRIGCTRL2": "High-gain trigger control register 2.",
    "IDX__TXHDRLGTRIGCTRL1": "Low-gain trigger control register 1.",
    "IDX__TXHDRLGTRIGCTRL2": "Low-gain trigger control register 2.",
    "IDX__TXHDRMGTRIGCTRL1": "Mid-gain trigger control register 1.",
    "IDX__TXHDRMGTRIGCTRL2": "Mid-gain trigger control register 2.",
    "IDX__TXHDRLSADC": "Low-speed ADC configuration/status register.",
    "IDX__TXHDRPOLSTAT": "Polarities/status register for triggering and processing.",
    "IDX__TXHDRPOLCTRL": "Polarity control register.",
    "IDX__TXHDRCOINENA": "Coincidence enable bitmask.",
    "IDX__TXHDRLSTRIGMODE": "Low-speed trigger mode selector.",
    "IDX__TXHDRMGTRIGMODE": "Mid-gain trigger mode selector.",
    "IDX__TXHDRLGTRIGMODE": "Low-gain trigger mode selector.",
    "IDX__TXHDRHGTRIGMODE": "High-gain trigger mode selector.",
    "IDX__TXHDRTOFMAX": "Packed TOF maximum values.",
    "IDX__TXHDRTOFMIN": "Packed TOF minimum values.",
    "IDX__TXHDRLS0MAXMIN": "Packed LS0 channel max/min values.",
    "IDX__TXHDRLS1MAXMIN": "Packed LS1 channel max/min values.",
    "IDX__TXHDRLS2MAXMIN": "Packed LS2 channel max/min values.",
    "IDX__TXHDRFIFODELAY": "FIFO delay register used in low-speed timing chain.",
    "IDX__TXHDRSAMPDELAY": "Packed high-speed per-channel sample-delay register.",
    "IDX__TXHDRTRANSCNT": "Transmit / transaction counter.",
    "IDX__TXHDRPROCHKCH01": "Processor housekeeping ADC channels 0 and 1 (packed).",
    "IDX__TXHDRPROCHKCH23": "Processor housekeeping ADC channels 2 and 3 (packed).",
    "IDX__TXHDRPROCHKCH45": "Processor housekeeping ADC channels 4 and 5 (packed).",
    "IDX__TXHDRPROCHKCH67": "Processor housekeeping ADC channels 6 and 7 (packed).",
    "IDX__TXHDRHVPSHKCH01": "HVPS housekeeping ADC channels 0 and 1 (packed).",
    "IDX__TXHDRHVPSHKCH23": "HVPS housekeeping ADC channels 2 and 3 (packed).",
    "IDX__TXHDRHVPSHKCH45": "HVPS housekeeping ADC channels 4 and 5 (packed).",
    "IDX__TXHDRHVPSHKCH67": "HVPS housekeeping ADC channels 6 and 7 (packed).",
    "IDX__TXHDRLVHK0CH01": "LVPS HK ADC0 channels 0 and 1 (packed).",
    "IDX__TXHDRLVHK0CH23": "LVPS HK ADC0 channels 2 and 3 (packed).",
    "IDX__TXHDRLVHK0CH45": "LVPS HK ADC0 channels 4 and 5 (packed).",
    "IDX__TXHDRLVHK0CH67": "LVPS HK ADC0 channels 6 and 7 (packed).",
    "IDX__TXHDRLVHK1CH01": "LVPS HK ADC1 channels 0 and 1 (packed).",
    "IDX__TXHDRLVHK1CH23": "LVPS HK ADC1 channels 2 and 3 (packed).",
    "IDX__TXHDRLVHK1CH45": "LVPS HK ADC1 channels 4 and 5 (packed).",
    "IDX__TXHDRLVHK1CH67": "LVPS HK ADC1 channels 6 and 7 (packed).",
    "IDX__TXHDRLVHK2CH01": "LVPS HK ADC2 channels 0 and 1 (packed).",
    "IDX__TXHDRLVHK2CH23": "LVPS HK ADC2 channels 2 and 3 (packed).",
    "IDX__TXHDRLVHK2CH45": "LVPS HK ADC2 channels 4 and 5 (packed).",
    "IDX__TXHDRLVHK2CH67": "LVPS HK ADC2 channels 6 and 7 (packed).",
    "IDX__TXHDRFSWAIDCOPY": "Copy of FSW application identifier.",
    "IDX__TXHDRFSWBINCOPY": "Copy of FSW build / binary identifier.",
    "IDX__TXHDRFSWMAJOR": "FSW major version number.",
    "IDX__TXHDRFSWMINOR": "FSW minor version number.",
    "IDX__TXHDRFSWPATCH": "FSW patch version number.",
    "IDX__TXHDRFSWHVSTAT": "FSW high-voltage status word.",
    "IDX__TXHDRFSWMEM0": "FSW memory status word 0.",
    "IDX__TXHDRFSWMEM1": "FSW memory status word 1.",
    "IDX__TXHDRFPGAVER": "FPGA version word.",
}

UNPACKED_FIELD_DESCRIPTIONS = {
    "HSPretriggerBlocks": "Decoded high-speed pre-trigger block count from IDX__TXHDRBLOCKS.",
    "HSPosttriggerBlocks": "Decoded high-speed post-trigger block count from IDX__TXHDRBLOCKS.",
    "LSPretriggerBlocks": "Decoded low-speed pre-trigger block count from IDX__TXHDRBLOCKS.",
    "LSPosttriggerBlocks": "Decoded low-speed post-trigger block count from IDX__TXHDRBLOCKS.",
    "TOFDelay_H": "Decoded TOF-H sample delay from IDX__TXHDRSAMPDELAY.",
    "TOFDelay_M": "Decoded TOF-M sample delay from IDX__TXHDRSAMPDELAY.",
    "TOFDelay_L": "Decoded TOF-L sample delay from IDX__TXHDRSAMPDELAY.",
    "TriggerOrigin": "Decoded trigger-origin labels from IDX__TXHDREVTNUM bitfield.",
    "TriggerOffsetMicroseconds": "Trigger offset converted to microseconds.",
    "HSPretriggerOffsetMicroseconds": "High-speed pre-trigger offset converted to microseconds.",
    "LSPretriggerOffsetMicroseconds": "Low-speed pre-trigger offset converted to microseconds.",
    "FIFODelayMicroseconds": "FIFO delay converted to microseconds.",
    "SampleDelayMicroseconds_TOF_H": "TOF-H sample delay converted to microseconds.",
    "SampleDelayMicroseconds_TOF_M": "TOF-M sample delay converted to microseconds.",
    "SampleDelayMicroseconds_TOF_L": "TOF-L sample delay converted to microseconds.",
    "epoch": "Combined spacecraft event epoch (seconds since spacecraft epoch).",
    "timestamp_utc": "ISO-8601 UTC timestamp derived from TXHDR time fields.",
}


def describe_field(name: str, *, packed: bool) -> Optional[str]:
    table = PACKED_FIELD_DESCRIPTIONS if packed else UNPACKED_FIELD_DESCRIPTIONS
    if name in table:
        return table[name]
    if name.startswith("IDX__TXHDRSP"):
        return "TXHDR spare/reserved word."
    if name.startswith("IDX__TXHDRFSWPT"):
        return "FSW pass-through telemetry word."
    if name.startswith("IDX__TXHDR"):
        return "TXHDR telemetry/header field from science packet."
    return None

