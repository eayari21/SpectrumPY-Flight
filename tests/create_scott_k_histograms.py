#!/usr/bin/env python3
"""
Utility script to accumulate histogram inputs from decoded IDEX .h5 files.

For every HDF5 file under the given root (recursing by default) the script
collects:
    * Target charges (Target L/Target H impact charge datasets)
    * Ion grid charges
    * Integrated detector charge      = Sum((signal - baseline) * gain * dt)
    * Peak detector current           = (peak(signal) - baseline) * gain

The gain and sampling interval are taken from HDF5 attributes when present.
If either is missing, a conservative fallback of ``gain=1`` and ``dt=1`` is
used so that the waveforms still contribute to the histograms.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np

ANALYSIS_GROUP = "Analysis"
DEFAULT_ROOT = Path(__file__).resolve().parent.parent / "Data"
DEFAULT_GAINS_MA_PER_DN: Dict[str, float] = {
    "TOF H": 2.89e-4,
    "TOF M": 1.13e-2,
    "TOF L": 5.14e-4,
    "Ion Grid": 7.46e-4,
    "Target H": 1.63e-1,
    "Target L": 1.58e1,
}


def iter_h5_files(root: Path, recursive: bool = True) -> List[Path]:
    paths: List[Path] = []
    if recursive:
        for dirpath, _, filenames in os.walk(root):
            for fn in filenames:
                if fn.lower().endswith(".h5"):
                    paths.append(Path(dirpath, fn))
    else:
        for fn in root.iterdir():
            if fn.is_file() and fn.name.lower().endswith(".h5"):
                paths.append(fn)
    paths.sort()
    return paths


def _normalize_spaces(name: str) -> str:
    return " ".join(str(name).split())


def _get_child(group: h5py.Group, candidate: str):
    if candidate in group:
        return group[candidate]
    wanted = _normalize_spaces(candidate)
    for key in group.keys():
        if _normalize_spaces(key) == wanted:
            return group[key]
    return None


def _read_scalar(ds: h5py.Dataset) -> Optional[float]:
    try:
        value = ds[()]
    except Exception:
        return None
    arr = np.asarray(value)
    if arr.shape == ():
        return float(arr)
    if arr.size == 1:
        return float(arr.reshape(()))
    return None


def _charge_datasets(analysis_group: h5py.Group, channels: Sequence[str]) -> List[float]:
    charges: List[float] = []
    for channel in channels:
        for suffix in ("ImpactCharge", " Impact Charge"):
            ds = _get_child(analysis_group, f"{channel}{suffix}")
            if isinstance(ds, h5py.Dataset):
                value = _read_scalar(ds)
                if value is not None and math.isfinite(value):
                    charges.append(value)
    return charges


def _channel_baseline(analysis_group: h5py.Group, channel: str) -> float:
    channel_group = analysis_group.get(channel)
    if isinstance(channel_group, h5py.Group):
        baseline = channel_group.attrs.get("Baseline")
        try:
            baseline_value = float(baseline)
            if math.isfinite(baseline_value):
                return baseline_value
        except Exception:
            pass
    return 0.0


def _channel_gain(analysis_group: h5py.Group, channel: str) -> float:
    channel_group = analysis_group.get(channel)
    if isinstance(channel_group, h5py.Group):
        for key in ("Gain", "Gain[mA/DN]", "Gain_mA_per_DN"):
            gain = channel_group.attrs.get(key)
            try:
                gain_value = float(gain)
                if math.isfinite(gain_value) and gain_value > 0:
                    return gain_value
            except Exception:
                continue
    return DEFAULT_GAINS_MA_PER_DN.get(channel, 1.0)


def _delta_time(event_group: h5py.Group, channel: str) -> float:
    if channel.startswith("TOF"):
        time_ds = event_group.get("Time (high sampling)")
    else:
        time_ds = event_group.get("Time (low sampling)")
    if isinstance(time_ds, h5py.Dataset):
        try:
            arr = np.asarray(time_ds[()], dtype=float)
            if arr.size >= 2:
                delta = np.diff(arr)
                finite_delta = delta[np.isfinite(delta)]
                if finite_delta.size:
                    return float(np.mean(finite_delta))
        except Exception:
            pass
    return 1.0


def _waveform(event_group: h5py.Group, channel: str) -> Optional[np.ndarray]:
    ds = _get_child(event_group, channel)
    if isinstance(ds, h5py.Dataset):
        try:
            return np.asarray(ds[()], dtype=float)
        except Exception:
            return None
    return None


def _compute_waveform_metrics(event_group: h5py.Group, analysis_group: h5py.Group, channel: str) -> Tuple[Optional[float], Optional[float]]:
    data = _waveform(event_group, channel)
    if data is None or data.size == 0:
        return None, None
    baseline = _channel_baseline(analysis_group, channel)
    gain = _channel_gain(analysis_group, channel)
    dt = _delta_time(event_group, channel)

    corrected = data - baseline
    integrated = float(np.sum(corrected) * gain * dt)
    peak = float(np.max(corrected) * gain)
    return integrated, peak


def process_file(path: Path) -> Tuple[List[float], List[float], List[float], List[float]]:
    target_charges: List[float] = []
    ion_grid_charges: List[float] = []
    integrated_charges: List[float] = []
    peak_currents: List[float] = []

    with h5py.File(path, "r") as handle:
        for event_name, event_obj in handle.items():
            if not isinstance(event_obj, h5py.Group):
                continue
            analysis_group = _get_child(event_obj, ANALYSIS_GROUP)
            if not isinstance(analysis_group, h5py.Group):
                continue

            target_charges.extend(_charge_datasets(analysis_group, ("Target L", "Target H")))
            ion_grid_charges.extend(_charge_datasets(analysis_group, ("Ion Grid",)))

            for channel in ("Target L", "Target H", "TOF H", "TOF M", "TOF L"):
                integrated, peak = _compute_waveform_metrics(event_obj, analysis_group, channel)
                if integrated is not None and math.isfinite(integrated):
                    integrated_charges.append(integrated)
                if peak is not None and math.isfinite(peak):
                    peak_currents.append(peak)

    return target_charges, ion_grid_charges, integrated_charges, peak_currents


def summarise(values: Sequence[float], label: str, bins: int = 50) -> None:
    if not values:
        print(f"No values collected for {label}.")
        return
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    print(f"Collected {len(values)} {label} value(s); {finite.size} finite.")
    if finite.size:
        counts, edges = np.histogram(finite, bins=bins)
        print(f"{label} histogram (bin edges then counts):")
        print(edges)
        print(counts)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Create aggregated histograms from .h5 files.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="Root directory to search for .h5 files (default: Data)")
    parser.add_argument("--no-recursive", action="store_true", help="Only search the top-level directory for .h5 files")
    args = parser.parse_args(argv)

    files = iter_h5_files(args.root, recursive=not args.no_recursive)
    if not files:
        print(f"No .h5 files found under {args.root}")
        return 1

    target_charges: List[float] = []
    ion_grid_charges: List[float] = []
    integrated_charges: List[float] = []
    peak_currents: List[float] = []

    for h5_path in files:
        t, ig, integ, peak = process_file(h5_path)
        target_charges.extend(t)
        ion_grid_charges.extend(ig)
        integrated_charges.extend(integ)
        peak_currents.extend(peak)

    summarise(target_charges, "target charge")
    summarise(ion_grid_charges, "ion grid charge")
    summarise(integrated_charges, "integrated detector charge")
    summarise(peak_currents, "peak detector current")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
