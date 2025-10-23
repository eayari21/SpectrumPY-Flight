"""Utility helpers shared between IDEX analysis scripts."""
from __future__ import annotations

from typing import Callable, Dict, Iterable

import numpy as np

__all__ = [
    "RISE_METRIC_SUFFIXES",
    "compute_rise_metrics",
]


RISE_METRIC_SUFFIXES: Dict[str, str] = {
    "t10": "10pct Time",
    "t90": "90pct Time",
    "v10": "10pct Value",
    "v90": "90pct Value",
    "rise": "10-90 Risetime",
}


def _interpolate_crossing(
    times: np.ndarray,
    values: np.ndarray,
    level: float,
    *,
    rising: bool,
) -> float:
    """Return the interpolated time where ``values`` first cross ``level``.

    The ``rising`` flag controls whether we look for the first sample on or above
    (``rising=True``) or on or below (``rising=False``) the requested ``level``.
    """
    if not np.isfinite(level):
        return float("nan")

    if rising:
        crossings = np.where(values >= level)[0]

        def compare(y: float) -> bool:
            return y < level

    else:
        crossings = np.where(values <= level)[0]

        def compare(y: float) -> bool:
            return y > level

    if crossings.size == 0:
        return float("nan")

    idx = int(crossings[0])
    if idx == 0:
        return float(times[0])

    prev_idx = idx - 1
    while prev_idx >= 0 and not compare(values[prev_idx]):
        prev_idx -= 1

    if prev_idx < 0:
        return float(times[idx])

    t0 = float(times[prev_idx])
    t1 = float(times[idx])
    y0 = float(values[prev_idx])
    y1 = float(values[idx])

    if not np.isfinite(y0) or not np.isfinite(y1):
        return float("nan")
    if y1 == y0:
        return float(t1)

    fraction = (level - y0) / (y1 - y0)
    return float(t0 + fraction * (t1 - t0))


def compute_rise_metrics(
    time: Iterable[float],
    fit_curve: Iterable[float],
    baseline: float | None = None,
) -> Dict[str, float]:
    """Compute 10/90 rise metrics for the provided fit curve.

    The threshold levels are defined as 10% and 90% of the span between the
    minimum and maximum values of the fit curve.  This matches the visual
    requirement in the UI of marking when the fit first reaches a given fraction
    of its full excursion.
    """
    times = np.asarray(list(time), dtype=float)
    values = np.asarray(list(fit_curve), dtype=float)

    metrics = {key: float("nan") for key in RISE_METRIC_SUFFIXES}

    if times.size == 0 or values.size == 0:
        return metrics

    mask = np.isfinite(times) & np.isfinite(values)
    if not np.any(mask):
        return metrics

    times = times[mask]
    values = values[mask]

    if times.size < 2:
        return metrics

    order = np.argsort(times)
    times = times[order]
    values = values[order]

    if baseline is None or not np.isfinite(baseline):
        baseline = 0.0

    min_idx = int(np.nanargmin(values))
    max_idx = int(np.nanargmax(values))

    min_val = float(values[min_idx])
    max_val = float(values[max_idx])

    if not (np.isfinite(min_val) and np.isfinite(max_val)):
        return metrics

    amplitude = max_val - min_val
    if amplitude <= 0 or not np.isfinite(amplitude):
        return metrics

    if min_idx <= max_idx:
        segment = slice(min_idx, max_idx + 1)
        low_val = min_val
        rising = True
    else:
        segment = slice(max_idx, min_idx + 1)
        low_val = min_val
        rising = False

    times_segment = times[segment]
    values_segment = values[segment]

    if times_segment.size < 2 or values_segment.size < 2:
        return metrics

    level_10 = low_val + 0.10 * amplitude
    level_90 = low_val + 0.90 * amplitude

    t10 = _interpolate_crossing(times_segment, values_segment, level_10, rising=rising)
    t90 = _interpolate_crossing(times_segment, values_segment, level_90, rising=rising)

    if np.isfinite(t10):
        metrics["t10"] = float(t10)
        metrics["v10"] = float(np.interp(t10, times, values))
    if np.isfinite(t90):
        metrics["t90"] = float(t90)
        metrics["v90"] = float(np.interp(t90, times, values))
    if np.isfinite(metrics["t10"]) and np.isfinite(metrics["t90"]):
        metrics["rise"] = float(metrics["t90"] - metrics["t10"])

    return metrics
