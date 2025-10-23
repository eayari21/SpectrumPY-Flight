"""Utility helpers shared between IDEX analysis scripts."""
from __future__ import annotations

from typing import Dict, Iterable

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


def _interpolate_crossing(times: np.ndarray, values: np.ndarray, level: float) -> float:
    """Return the interpolated time where ``values`` first cross ``level``."""
    if not np.isfinite(level):
        return float("nan")

    above = np.where(values >= level)[0]
    if above.size == 0:
        return float("nan")

    idx = int(above[0])
    if idx == 0:
        return float(times[0])

    t0 = float(times[idx - 1])
    t1 = float(times[idx])
    y0 = float(values[idx - 1])
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
    """Compute 10/90 rise metrics for the provided fit curve."""
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
        window = max(1, int(round(0.1 * values.size)))
        baseline = float(np.nanmedian(values[:window]))
    top = float(np.nanmax(values))
    amplitude = top - baseline

    if not np.isfinite(amplitude) or amplitude <= 0:
        return metrics

    level_10 = baseline + 0.10 * amplitude
    level_90 = baseline + 0.90 * amplitude

    t10 = _interpolate_crossing(times, values, level_10)
    t90 = _interpolate_crossing(times, values, level_90)

    if np.isfinite(t10):
        metrics["t10"] = float(t10)
        metrics["v10"] = float(np.interp(t10, times, values))
    if np.isfinite(t90):
        metrics["t90"] = float(t90)
        metrics["v90"] = float(np.interp(t90, times, values))
    if np.isfinite(metrics["t10"]) and np.isfinite(metrics["t90"]):
        metrics["rise"] = float(metrics["t90"] - metrics["t10"])

    return metrics
