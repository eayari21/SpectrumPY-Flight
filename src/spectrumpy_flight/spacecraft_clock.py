"""Utilities for working with the spacecraft clock reference epoch."""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
from typing import Union

Number = Union[int, float]

SPACECRAFT_EPOCH = datetime(2010, 1, 1, tzinfo=timezone.utc)
"""Start of the spacecraft clock (Jan 1, 2010 UTC)."""


def combine_coarse_fine_seconds(
    coarse: Number,
    fine: Number,
    *,
    fine_resolution_seconds: float = 20e-6,
) -> float:
    """Combine coarse and fine counters into elapsed seconds.

    Parameters
    ----------
    coarse
        Integral seconds reported by the spacecraft clock.
    fine
        Fractional counter representing ``fine_resolution_seconds`` per count.
    fine_resolution_seconds
        Number of seconds represented by one fine counter increment.  The
        default mirrors the value used in the packet decoder.

    Returns
    -------
    float
        Total elapsed seconds relative to :data:`SPACECRAFT_EPOCH`.
    """

    return float(coarse) + float(fine) * float(fine_resolution_seconds)


def spacecraft_seconds_to_datetime(seconds: Number) -> datetime:
    """Convert elapsed spacecraft-clock seconds to a UTC ``datetime``."""

    return SPACECRAFT_EPOCH + timedelta(seconds=float(seconds))


def spacecraft_seconds_to_unix_seconds(seconds: Number) -> float:
    """Convert elapsed spacecraft-clock seconds to Unix epoch seconds."""

    return spacecraft_seconds_to_datetime(seconds).timestamp()


__all__ = [
    "SPACECRAFT_EPOCH",
    "combine_coarse_fine_seconds",
    "spacecraft_seconds_to_datetime",
    "spacecraft_seconds_to_unix_seconds",
]
