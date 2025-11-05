#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
An algorithm to convert the time axis to mass for TOF spectra.
__author__      = Ethan Ayari, Mihály Horányi,
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""
import os
import sys
import h5py
import random
import csv
import math
import copy
from pathlib import Path
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns

if __package__ is None or __package__ == "":
    _MODULE_DIR = Path(__file__).resolve().parent
    _PACKAGE_ROOT = _MODULE_DIR.parent
    for _path in (_MODULE_DIR, _PACKAGE_ROOT):
        _path_str = str(_path)
        if _path_str not in sys.path:
            sys.path.append(_path_str)
    from plot_style import apply_plot_style
    from mass_calibration import TOFMassCal
else:
    from .plot_style import apply_plot_style
    from .mass_calibration import TOFMassCal

from matplotlib import colors
from scipy.signal import find_peaks

apply_plot_style()
plt.rcParams['agg.path.chunksize'] = 10_000

MASS_STRETCH_MIN_US = 1.3
MASS_STRETCH_MAX_US = 1.6
MAX_CALIBRATION_ORDER = 4

MIN_EXTRA_MASS_AMU = 0.75
MIN_MASS_LINE_SEPARATION_AMU = 0.2


@dataclass(frozen=True)
class MassReference:
    int_mass: int
    mass_value: float
    label: str
    species: str


_MASS_REFERENCE_SPEC: Sequence[Tuple[int, str, Sequence[str]]] = (
    (1, "H", ("H",)),
    (12, "C", ("C",)),
    (16, "O", ("O",)),
    (23, "Na", ("Na",)),
    (24, "Mg", ("Mg24",)),
    (25, "Mg", ("Mg25",)),
    (26, "Mg", ("Mg26",)),
    (28, "Si", ("Si28",)),
    (29, "Si", ("Si29",)),
    (30, "Si", ("Si30",)),
    (39, "K", ("K39",)),
    (40, "Ca", ("Ca40",)),
    (54, "Fe", ("Fe54",)),
    (56, "Fe", ("Fe56",)),
    (57, "Fe", ("Fe57",)),
    (58, "Fe", ("Fe58",)),
)


_LAST_ASSIGNMENTS: Optional[Dict[str, object]] = None


def _resolve_mass_row(table: pd.DataFrame, int_mass: int, candidates: Sequence[str]) -> Tuple[float, str]:
    mask = table["Name"].astype(str).isin(candidates)
    if mask.any():
        row = table.loc[mask].iloc[0]
    else:
        masses = table["Mass"].astype(float).to_numpy(copy=False)
        idx = int(np.argmin(np.abs(masses - float(int_mass))))
        row = table.iloc[idx]
    return float(row["Mass"]), str(row["Name"])


@lru_cache(maxsize=1)
def _load_reference_masses() -> Tuple[MassReference, ...]:
    path = Path(__file__).resolve().with_name("mass_comb.csv")
    try:
        table = pd.read_csv(path)
    except Exception as exc:
        raise RuntimeError(f"Unable to read mass combination table: {exc}") from exc

    references: List[MassReference] = []
    for int_mass, species, candidates in _MASS_REFERENCE_SPEC:
        try:
            mass_value, label = _resolve_mass_row(table, int_mass, candidates)
        except Exception:
            mass_value = float(int_mass)
            label = candidates[0] if candidates else str(int_mass)
        references.append(MassReference(
            int_mass=int_mass,
            mass_value=float(mass_value),
            label=str(label),
            species=species,
        ))
    references.sort(key=lambda ref: ref.mass_value)
    return tuple(references)


def _prepare_signal(tof: Sequence[float]) -> np.ndarray:
    data = np.asarray(tof, dtype=float)
    if data.size == 0:
        return data
    baseline_window = max(32, data.size // 10)
    baseline = float(np.nanmedian(data[:baseline_window]))
    adjusted = np.nan_to_num(data - baseline, nan=0.0, posinf=0.0, neginf=0.0)
    adjusted[np.isnan(adjusted)] = 0.0
    return np.clip(adjusted, 0.0, None)


def _smooth_signal(signal: np.ndarray, width: int = 5) -> np.ndarray:
    if signal.size < width or width <= 1:
        return signal
    kernel = np.ones(width, dtype=float) / float(width)
    return np.convolve(signal, kernel, mode="same")


def _robust_std(values: np.ndarray) -> float:
    if values.size == 0:
        return 0.0
    median = float(np.median(values))
    return 1.4826 * float(np.median(np.abs(values - median)))


def _compute_mass_axis(time_zero_us: np.ndarray, stretch_us: float, shift_us: float) -> np.ndarray:
    shifted = time_zero_us - shift_us
    masses = np.empty_like(time_zero_us, dtype=float)
    valid = shifted > 0
    masses[valid] = (shifted[valid] / stretch_us) ** 2
    masses[~valid] = 0.0
    return np.clip(masses, 0.0, 400.0)


def _select_evenly_spaced_indices(indices: Sequence[int], count: int) -> List[int]:
    """Return ``count`` indices that are as evenly spaced as possible."""

    unique_sorted = np.array(sorted({int(idx) for idx in indices}), dtype=int)
    if count <= 0 or unique_sorted.size == 0:
        return []
    if count >= unique_sorted.size:
        return unique_sorted.tolist()

    segments = np.array_split(unique_sorted, count)
    selected: List[int] = []
    for segment in segments:
        if segment.size == 0:
            continue
        selected.append(int(segment[segment.size // 2]))

    if len(selected) < count:
        remaining = [int(idx) for idx in unique_sorted if int(idx) not in selected]
        for idx in remaining:
            selected.append(idx)
            if len(selected) >= count:
                break

    return selected[:count]


def _assign_mass_lines(
    detection_signal: np.ndarray,
    fit_signal: np.ndarray,
    time_axis: np.ndarray,
    time_zero: np.ndarray,
    step_us: float,
    stretch_us: float,
    shift_us: float,
    mass_scale: np.ndarray,
    references: Sequence[MassReference],
    calibration: Optional[TOFMassCal] = None,
    *,
    max_remaining_peaks: Optional[int] = None,
) -> Dict[str, object]:
    if detection_signal.size == 0 or step_us <= 0.0:
        return {
            "peaks": np.array([], dtype=int),
            "mass_lines": [],
            "stretch": float(stretch_us),
            "shift": float(shift_us),
            "step": float(step_us),
        }

    if max_remaining_peaks is None:
        max_total_lines: Optional[int] = None
    else:
        try:
            max_total_lines = int(max_remaining_peaks)
        except (TypeError, ValueError):
            max_total_lines = None

    robust_noise = _robust_std(detection_signal)
    absolute_max = float(np.nanmax(detection_signal)) if detection_signal.size else 0.0
    prominence = max(absolute_max * 0.015, robust_noise * 4.5, 1e-6)
    distance = max(5, int(round(0.18 / max(step_us, 1e-6))))
    peaks, _ = find_peaks(detection_signal, prominence=prominence, distance=distance)
    peaks = peaks.astype(int, copy=False)

    if max_total_lines is not None and max_total_lines > 0:
        desired_peak_count = max(max_total_lines, len(references))
        adjusted_prominence = prominence
        for _ in range(5):
            if peaks.size <= desired_peak_count * 3:
                break
            adjusted_prominence *= 1.5
            new_peaks, _ = find_peaks(
                detection_signal,
                prominence=adjusted_prominence,
                distance=distance,
            )
            new_peaks = new_peaks.astype(int, copy=False)
            if new_peaks.size >= peaks.size:
                break
            peaks = new_peaks
        prominence = adjusted_prominence

    tolerance = max(6, int(round(0.22 / max(step_us, 1e-6))))
    half_window = max(8, int(round(0.35 / max(step_us, 1e-6))))

    available_peaks = list(peaks)
    used_indices = set()
    noise_floor = robust_noise
    min_peak_height = max(absolute_max * 0.0025, noise_floor * 3.5, 1e-6)
    tail_threshold = max(min_peak_height * 0.6, noise_floor * 2.5, 1e-9)
    mass_lines: List[Dict[str, object]] = []
    accepted_masses: List[float] = []
    line_id = 1
    for reference in references:
        if calibration is not None:
            expected_time_zero = float(calibration.mass_to_tof([reference.mass_value])[0])
        else:
            expected_time_zero = shift_us + stretch_us * math.sqrt(reference.mass_value)
        expected_index = int(round(expected_time_zero / step_us)) if step_us > 0 else 0
        if expected_index < 0 or expected_index >= detection_signal.size:
            continue
        peak_index: Optional[int] = None
        if available_peaks:
            distances = [abs(idx - expected_index) for idx in available_peaks]
            best_pos = int(np.argmin(distances))
            if distances[best_pos] <= tolerance:
                peak_index = int(available_peaks.pop(best_pos))
                used_indices.add(peak_index)

        if peak_index is None:
            start_idx = max(0, expected_index - tolerance)
            end_idx = min(detection_signal.size, expected_index + tolerance + 1)
            if end_idx - start_idx >= 3:
                local_window = detection_signal[start_idx:end_idx]
                local_idx = int(np.argmax(local_window))
                candidate_index = start_idx + local_idx
                if candidate_index not in used_indices:
                    candidate_height = float(detection_signal[candidate_index])
                    local_baseline = float(np.nanmedian(local_window))
                    if (
                        candidate_height >= min_peak_height
                        and candidate_height - local_baseline >= max(noise_floor, 1e-6)
                    ):
                        peak_index = candidate_index
                        used_indices.add(peak_index)
                        if peak_index in available_peaks:
                            available_peaks.remove(peak_index)

        if peak_index is None:
            continue
        if mass_scale.size <= peak_index:
            continue
        mass_at_peak = float(mass_scale[peak_index])
        if not math.isfinite(mass_at_peak) or abs(mass_at_peak - reference.int_mass) > 1.25:
            continue
        start_idx = peak_index
        while start_idx > 0 and detection_signal[start_idx] >= tail_threshold:
            start_idx -= 1
        end_idx = peak_index
        last_index = detection_signal.size - 1
        while end_idx < last_index and detection_signal[end_idx] >= tail_threshold:
            end_idx += 1
        start_idx = max(0, start_idx - 2)
        end_idx = min(detection_signal.size, end_idx + 3)
        if end_idx - start_idx < 6:
            start_idx = max(0, peak_index - half_window)
            end_idx = min(detection_signal.size, peak_index + half_window + 1)
        mass_lines.append({
            "line_id": line_id,
            "int_mass": reference.int_mass,
            "label": reference.label,
            "species": reference.species,
            "mass_reference": reference.mass_value,
            "peak_index": peak_index,
            "time": float(time_axis[peak_index]) if time_axis.size > peak_index else float("nan"),
            "time_offset": float(time_zero[peak_index]) if time_zero.size > peak_index else float("nan"),
            "expected_time": float(time_axis[0] + expected_time_zero) if time_axis.size else float(expected_time_zero),
            "window": (start_idx, end_idx),
            "peak_amplitude": float(fit_signal[peak_index]) if fit_signal.size > peak_index else 0.0,
            "mass_scale_value": mass_at_peak,
        })
        line_id += 1
        if math.isfinite(mass_at_peak):
            accepted_masses.append(float(mass_at_peak))

    if max_total_lines is None:
        remaining_slots: Optional[int] = None
    else:
        remaining_slots = max(0, max_total_lines - len(mass_lines))

    if available_peaks:
        remaining = sorted(set(available_peaks) - used_indices)
        if remaining_slots is not None:
            if remaining_slots == 0:
                remaining = []
            elif len(remaining) > remaining_slots:
                remaining = _select_evenly_spaced_indices(remaining, remaining_slots)
            else:
                remaining = _select_evenly_spaced_indices(remaining, len(remaining))
        else:
            remaining = _select_evenly_spaced_indices(remaining, len(remaining))
        for peak_index in remaining:
            if peak_index < 0 or peak_index >= detection_signal.size:
                continue
            height = float(detection_signal[peak_index])
            if not np.isfinite(height) or height < min_peak_height:
                continue
            if mass_scale.size > peak_index:
                mass_scale_value = float(mass_scale[peak_index])
            else:
                mass_scale_value = float("nan")
            if not math.isfinite(mass_scale_value):
                continue
            if mass_scale_value < MIN_EXTRA_MASS_AMU:
                continue
            if any(abs(mass_scale_value - existing) < MIN_MASS_LINE_SEPARATION_AMU for existing in accepted_masses):
                continue
            start_idx = peak_index
            while start_idx > 0 and detection_signal[start_idx] >= tail_threshold:
                start_idx -= 1
            end_idx = peak_index
            last_index = detection_signal.size - 1
            while end_idx < last_index and detection_signal[end_idx] >= tail_threshold:
                end_idx += 1
            start_idx = max(0, start_idx - 2)
            end_idx = min(detection_signal.size, end_idx + 3)
            if end_idx - start_idx < 6:
                start_idx = max(0, peak_index - half_window)
                end_idx = min(detection_signal.size, peak_index + half_window + 1)
            mass_lines.append({
                "line_id": line_id,
                "int_mass": int(round(mass_scale_value)) if math.isfinite(mass_scale_value) else 0,
                "label": f"Line{line_id}",
                "species": "",
                "mass_reference": float("nan"),
                "peak_index": peak_index,
                "time": float(time_axis[peak_index]) if time_axis.size > peak_index else float("nan"),
                "time_offset": float(time_zero[peak_index]) if time_zero.size > peak_index else float("nan"),
                "expected_time": float(time_axis[peak_index]) if time_axis.size > peak_index else float("nan"),
                "window": (start_idx, end_idx),
                "peak_amplitude": float(fit_signal[peak_index]) if fit_signal.size > peak_index else 0.0,
                "mass_scale_value": mass_scale_value,
            })
            line_id += 1
            accepted_masses.append(float(mass_scale_value))
            if remaining_slots is not None:
                remaining_slots -= 1
                if remaining_slots <= 0:
                    break

    return {
        "peaks": peaks,
        "mass_lines": mass_lines,
        "stretch": float(stretch_us),
        "shift": float(shift_us),
        "step": float(step_us),
        "calibration": calibration.to_dict() if calibration is not None else None,
    }


def get_last_mass_line_assignments() -> Dict[str, object]:
    if _LAST_ASSIGNMENTS is None:
        return {
            "peaks": np.array([], dtype=int),
            "mass_lines": [],
            "stretch": float("nan"),
            "shift": float("nan"),
            "step": float("nan"),
            "calibration": None,
            "origin": float("nan"),
        }
    # Return a deep copy to prevent callers from mutating cached data.
    result = copy.deepcopy(_LAST_ASSIGNMENTS)
    if "peaks" in result:
        result["peaks"] = np.asarray(result["peaks"], dtype=int)
    return result


def peak_time2mass(TOF, time):
    """Same as the original time2mass algorithm, but now we use peak locations 
    to create a simpler comb where every mass line has an amplitude of 1."""
    mass_comb = pd.read_csv("mass_comb.csv")

    masses = np.round(mass_comb["Mass"],1)

    time = time - time[0]
    # plt.plot(time, TOF, label="Time Scale")
    # plt.show()

    random_stretch = np.linspace(1400, 1500, 10)
    shift = []
    corr = []
    product = []

    # Find indices of all peaks (local maxima) and troughs (local minima)
    peaks, _ = find_peaks(TOF, prominence=20)

    # troughs, _ = find_peaks(-TOF)

    # Combine indices of peaks and troughs to get all extrema
    extrema_indices = np.sort(peaks)

    # || Plot the extrema
    plt.figure(figsize=(10, 6))
    plt.plot(time, TOF, label='TOF vs Time')
    plt.plot(time[extrema_indices], TOF[extrema_indices], 'rx', label='Extrema')

    plt.title('TOF vs Time with Extrema Marked')
    plt.xlabel('Time [s]')
    plt.ylabel('TOF [dN]')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Initialize TOF comb array
    TOF_comb = np.zeros(len(TOF))
    # Set extrema and neighboring indices to 1
    for idx in extrema_indices:
        TOF_comb[max(0, idx-20):min(len(TOF), idx+20)] = 10

    

    # || Plot the TOF comb
    plt.figure(figsize=(10, 6))
    plt.plot(time, TOF_comb)

    plt.title('TOF Correlation Comb')
    plt.xlabel('Time [s]')
    plt.ylabel('TOF Comb')
    plt.legend()
    plt.grid(True)
    plt.show()

    for stretch in random_stretch:

        # || Step 1): Make a vector with all zeros and a length of 8189, same as the TOF length of record: t_i
        t_i = np.zeros(len(TOF))

        # || 
        m_i = np.linspace(1, 200, 200)
        # Set elements in m_i to zero if they are not in masses
        # m_i = np.where(np.isin(m_i, masses), m_i, 0)
        # print(f"m_i = {m_i}, masses = {masses}")

        # || Step 3): Calculate the 200 times (T_i-s)  these masses should show up for a T_off = 0  and stretch A (best guess) T_i = T_off + A*sqrt(m_i) and set the closest t_i =1 to each of the calculated Ti-s, the rest stays zero.

        # Guess a shift of zero
        T_off = 0

        # Guess a stretch of 1462 ns (electronics)
        A = stretch

        T_i = T_off + A*np.sqrt(masses)

        # Loop through the elements of T_i and set corresponding elements in t_i to 1
        for idx in np.round(T_i).astype(int):
            # if(0 <= idx < len(t_i) and (t_i[idx] in masses)):
            if(0 <= idx < len(t_i)):
                # print(idx)
                t_i[idx] = 1

        # || Step 4): Calculate the cross-correlation with the original TOF, the max will give you the best lag (T_off) for a given A.

        # Cross-correlate T_i with TOF
        cross_correlation = np.correlate(t_i, TOF_comb, mode='full')

        lags = np.linspace(1, len(cross_correlation), len(cross_correlation))
        plt.scatter(lags, cross_correlation)
        plt.xlabel("Lags", fontsize=14, fontweight="bold")
        plt.ylabel("Correlation", fontsize=14, fontweight="bold")
        plt.title(f"Correlation for stretch = {stretch}")
        plt.show()

        # lags = np.linspace(1, len(TOF), len(TOF))
        # plt.scatter(lags, t_i*TOF)
        # plt.xlabel("Lags", fontsize=14, fontweight="bold")
        # plt.ylabel(r"Product $t_{i} \cdot TOF$", fontsize=14, fontweight="bold")
        # plt.title(f"Product for stretch = {stretch}")
        # plt.show()

        # Find the lag corresponding to the maximum correlation
        best_lag = np.argmax(cross_correlation)

        # Filter out lags outside the specified range
        # if 0 <= best_lag <= 8000:
        #     shift.append(best_lag)
        #     corr.append(cross_correlation[np.argmax(cross_correlation)])
    

        # t_i = t_i - best_lag

        # print(f"||===Best Lag = {best_lag}===||")
        shift.append(best_lag)
        # print(f"||===Correlation = {cross_correlation[best_lag]}===||")
        corr.append(cross_correlation[best_lag])

        product.append(np.sum(t_i*TOF))
        # print(f"||===STRETCH = {stretch}, LAG = {best_lag}, CORRELATION = {cross_correlation[best_lag]}")


    print(f"Best correlation = {max(corr)}, best stretch = {random_stretch[np.argmax(corr)]} best lag = {shift[np.argmax(corr)]}")

    # t = a*sqrt(m)

    # plt.scatter(random_stretch, corr, c="black", label=r"Stretch parameter $\alpha$")
    # plt.xlabel(r"Stretch parameter $\alpha \, [\frac{s^{2}}{amu}]$", fontsize=14, fontweight="bold")
    # plt.ylabel("Correlation", fontsize=14, fontweight="bold")
    # plt.show()


    # plt.scatter(shift, corr, c="red", label=r"Shift parameter $t_{0}$")
    # plt.xlabel(r"Shift parameter $t_{0}$ [s]", fontsize=14, fontweight="bold")
    # plt.ylabel("Correlation", fontsize=14, fontweight="bold")
    # plt.show()

    # plt.scatter(shift, product, c="red", label=r"Shift parameter $t_{0}$")
    # plt.xlabel(r"Shift parameter $t_{0}$ [s]", fontsize=14, fontweight="bold")
    # plt.ylabel("Product", fontsize=14, fontweight="bold")
    # plt.show()

    # samplesize = 1e-9*.75  # Oscilloscope sampling rate
    samplesize = 0.0038466235767167234e-6  # FM sampling rate (quartz oscillator)

    mass_scale = [((time-samplesize*shift[np.argmax(corr)])/(1e-9*random_stretch[np.argmax(corr)]))**2 for time in time]
    # mass_scale = [((time-samplesize*9416)/(1e-9*1444))**2 for time in time]
    # mass_scale = [((time*1e3)/(random_stretch[np.argmax(corr)]))**2 for time in time]

    # print("mass scale max/min = ",min(mass_scale), " , ", max(mass_scale))
    # print(f"&&&&&& mass scale = {mass_scale}")


    # best_lag is the lag that aligns T_i with TOF for the maximum correlation

    # || Step 5): Use an extrema finder to get to the best A, alternatively, one could try optimizing for the two parameters simultaneously.



    fig, ax1 = plt.subplots()

    # Plotting TOF against mass_scale on the main axis
    ax1.plot(mass_scale, TOF, 'b-')
    ax1.set_xlabel('Mass Scale [amu]', fontsize=14)
    ax1.set_ylabel('TOF [dN]', fontsize=14)

    # Create a secondary x-axis for time
    secax = ax1.secondary_xaxis('top')
    secax.set_xlabel('Time [s]', fontsize=14)

    # Define tick positions and compute corresponding time values
    best_shift = shift[np.argmax(corr)]*1e-9
    best_stretch = random_stretch[np.argmax(corr)]*1e-9

    # Define tick positions and compute corresponding time values
    tick_positions = np.linspace(min(mass_scale), max(mass_scale), num=4)  # Adjust num as needed
    tick_labels = [f'{np.sqrt(t * best_stretch**2 + best_shift**2):.2e}' for t in tick_positions]


    secax.set_xticks(tick_positions)
    secax.set_xticklabels(tick_labels)

    plt.title('TOF vs Mass Scale and Time')
    plt.grid(True)
    plt.show()

def time2mass(TOF, time, *, allow_out_of_range: bool = False, max_auto_lines: Optional[int] = None):
    """Return optimised mass-axis parameters for a TOF waveform."""

    global _LAST_ASSIGNMENTS

    references = _load_reference_masses()
    tof = np.asarray(TOF, dtype=float)
    time_axis = np.asarray(time, dtype=float)
    origin_us = float(time_axis[0]) if time_axis.size else 0.0

    if tof.size == 0 or time_axis.size == 0 or tof.size != time_axis.size:
        mass_scale = np.zeros_like(time_axis, dtype=float)
        _LAST_ASSIGNMENTS = {
            "peaks": np.array([], dtype=int),
            "mass_lines": [],
            "stretch": float("nan"),
            "shift": float("nan"),
            "step": float("nan"),
            "calibration": None,
            "origin": origin_us,
        }
        return float(MASS_STRETCH_MIN_US), 0.0, mass_scale

    time_zero = time_axis - time_axis[0]
    step_us = float(np.median(np.diff(time_zero))) if time_zero.size > 1 else 0.0
    if not np.isfinite(step_us) or step_us <= 0.0:
        mass_scale = np.zeros_like(time_axis, dtype=float)
        _LAST_ASSIGNMENTS = {
            "peaks": np.array([], dtype=int),
            "mass_lines": [],
            "stretch": float("nan"),
            "shift": float("nan"),
            "step": float(step_us),
            "calibration": None,
            "origin": origin_us,
        }
        return float(MASS_STRETCH_MIN_US), 0.0, mass_scale

    detection_signal = _prepare_signal(tof)
    smoothed = _smooth_signal(detection_signal)

    sqrt_masses = np.sqrt(np.array([ref.mass_value for ref in references], dtype=float))
    if sqrt_masses.size == 0:
        mass_scale = np.zeros_like(time_axis, dtype=float)
        _LAST_ASSIGNMENTS = {
            "peaks": np.array([], dtype=int),
            "mass_lines": [],
            "stretch": float("nan"),
            "shift": float("nan"),
            "step": float(step_us),
            "calibration": None,
            "origin": origin_us,
        }
        return float(MASS_STRETCH_MIN_US), 0.0, mass_scale

    stretch_min = MASS_STRETCH_MIN_US
    stretch_max = MASS_STRETCH_MAX_US
    if allow_out_of_range:
        stretch_min = min(stretch_min, 0.6)
        stretch_max = max(stretch_max, 2.5)

    weights = 1.0 / np.maximum(sqrt_masses, 1.0)
    best: Optional[Dict[str, float]] = None

    def evaluate(stretch_values: Iterable[float], shift_window: Optional[Tuple[int, int]] = None) -> None:
        nonlocal best
        for stretch in stretch_values:
            if stretch <= 0:
                continue
            expected_times = stretch * sqrt_masses
            base_indices = np.round(expected_times / step_us).astype(int)
            if base_indices.size == 0:
                continue
            shift_min_samples = int(-base_indices[0])
            shift_max_samples = int(smoothed.size - 1 - base_indices[-1])
            if shift_max_samples < shift_min_samples:
                continue
            if shift_window is not None:
                shift_min_clamped = max(shift_min_samples, shift_window[0])
                shift_max_clamped = min(shift_max_samples, shift_window[1])
            else:
                shift_min_clamped = shift_min_samples
                shift_max_clamped = shift_max_samples
            if shift_max_clamped < shift_min_clamped:
                continue
            shifts = np.arange(shift_min_clamped, shift_max_clamped + 1, dtype=int)
            if shifts.size == 0:
                continue
            index_matrix = base_indices[np.newaxis, :] + shifts[:, np.newaxis]
            index_matrix = np.clip(index_matrix, 0, smoothed.size - 1)
            values = smoothed[index_matrix]
            scores = (values * weights[np.newaxis, :]).sum(axis=1)
            best_idx = int(np.argmax(scores))
            score = float(scores[best_idx])
            shift_samples = int(shifts[best_idx])
            if best is None or score > best["score"]:
                best = {
                    "stretch": float(stretch),
                    "shift_samples": float(shift_samples),
                    "score": score,
                }

    coarse_stretches = np.linspace(stretch_min, stretch_max, 181)
    evaluate(coarse_stretches)

    if best is None:
        fallback = float(stretch_min)
        mass_scale = _compute_mass_axis(time_zero, fallback, 0.0)
        _LAST_ASSIGNMENTS = _assign_mass_lines(
            smoothed,
            tof,
            time_axis,
            time_zero,
            step_us,
            fallback,
            0.0,
            mass_scale,
            references,
            max_remaining_peaks=max_auto_lines,
        )
        _LAST_ASSIGNMENTS["origin"] = origin_us
        return fallback, 0.0, mass_scale

    refined_min = best["stretch"] - 0.01
    refined_max = best["stretch"] + 0.01
    refined_min = max(stretch_min, refined_min)
    refined_max = min(stretch_max, refined_max) if not allow_out_of_range else max(refined_min + 1e-6, refined_max)
    refined_stretches = np.linspace(refined_min, refined_max, 41)
    shift_center = int(round(best["shift_samples"]))
    shift_window = (shift_center - 60, shift_center + 60)
    evaluate(refined_stretches, shift_window=shift_window)

    shift_samples = int(round(best["shift_samples"]))
    shift_us = float(shift_samples * step_us)
    stretch_us = float(best["stretch"])
    mass_scale = _compute_mass_axis(time_zero, stretch_us, shift_us)
    assignments = _assign_mass_lines(
        smoothed,
        tof,
        time_axis,
        time_zero,
        step_us,
        stretch_us,
        shift_us,
        mass_scale,
        references,
        max_remaining_peaks=max_auto_lines,
    )
    assignments["origin"] = origin_us

    if assignments.get("mass_lines"):
        shift_candidates: List[float] = []
        for line in assignments["mass_lines"]:
            mass_reference = line.get("mass_reference")
            time_offset = line.get("time_offset")
            if mass_reference is None or time_offset is None:
                continue
            try:
                mass_reference_value = float(mass_reference)
                time_offset_value = float(time_offset)
            except Exception:
                continue
            if not (math.isfinite(mass_reference_value) and math.isfinite(time_offset_value)):
                continue
            if mass_reference_value <= 0.0:
                continue
            shift_candidates.append(time_offset_value - stretch_us * math.sqrt(mass_reference_value))
        if shift_candidates:
            refined_shift = float(np.median(shift_candidates))
            if math.isfinite(refined_shift):
                shift_us = refined_shift
                mass_scale = _compute_mass_axis(time_zero, stretch_us, shift_us)
                assignments = _assign_mass_lines(
                    smoothed,
                    tof,
                    time_axis,
                    time_zero,
                    step_us,
                    stretch_us,
                    shift_us,
                    mass_scale,
                    references,
                    max_remaining_peaks=max_auto_lines,
                )
                assignments["origin"] = origin_us

    calibration_model: Optional[TOFMassCal] = None
    reference_masses: List[float] = []
    reference_times: List[float] = []
    if assignments.get("mass_lines"):
        for line in assignments["mass_lines"]:
            try:
                mass_reference = float(line.get("mass_reference", float("nan")))
            except Exception:
                mass_reference = float("nan")
            try:
                time_offset = float(line.get("time_offset", float("nan")))
            except Exception:
                time_offset = float("nan")
            if not (np.isfinite(mass_reference) and np.isfinite(time_offset)):
                continue
            reference_masses.append(mass_reference)
            reference_times.append(time_offset)

    if len(reference_masses) >= 2:
        try:
            calibration_model = TOFMassCal.from_lines(
                reference_masses,
                reference_times,
                max_order=MAX_CALIBRATION_ORDER,
                mass_range=(0.0, 400.0),
                enforce_monotonic=True,
            )
        except Exception:
            calibration_model = None

    if calibration_model is not None:
        shift_us = float(calibration_model.coeffs[0])
        stretch_us = float(calibration_model.coeffs[1]) if calibration_model.coeffs.size > 1 else stretch_us
        mass_scale = calibration_model.tof_to_mass(time_zero)
        assignments = _assign_mass_lines(
            smoothed,
            tof,
            time_axis,
            time_zero,
            step_us,
            stretch_us,
            shift_us,
            mass_scale,
            references,
            calibration=calibration_model,
            max_remaining_peaks=max_auto_lines,
        )
        assignments["origin"] = origin_us
    else:
        assignments["calibration"] = None

    if not allow_out_of_range:
        stretch_us = float(np.clip(stretch_us, stretch_min, stretch_max))
    assignments["stretch"] = float(stretch_us)
    assignments["shift"] = float(shift_us)
    _LAST_ASSIGNMENTS = assignments

    return stretch_us, shift_us, mass_scale

# ||
# ||
# || Test code: Import file and run a single waveform
# %%
if __name__ == "__main__":
    # || Parameters of interest
    amplitude = []
    time = []
    mass = []

    # || Read in the Fe on EM scope .h5 file
    fname = "/Users/ethanayari/Dropbox/IDEX Pipeline/3_Codes/2_Level 0->1/Python_Packet_Crushing/IDEX_Decom/TL/1/TL1.h5"

    with h5py.File(fname, "r") as f:
        spectra_group = f["Spectra"]
        
        for key in spectra_group.keys():
            spectrum = spectra_group[key]
            
            if isinstance(spectrum, h5py.Group):
                if "Amplitude" in spectrum.keys():
                    amplitude.append(spectrum["Amplitude"][:])
                
                if "Time" in spectrum.keys():
                    time.append(spectrum["Time"][:])
                
                if "Mass" in spectrum.keys():
                    mass.append(spectrum["Mass"][:])

    # || Read in the aluminum FM .h5 file
    # fname = "/Users/impact/Dropbox/Mac (2)/Desktop/Post_Env/Aluminum/science_data/20231218/run0_aluminum_vel_gt_3_HGCSA/HDF5/ois_output_12182023_184030.h5"

    # || Read in the olivine FM .h5 file
    # fname = "/Users/impact/Dropbox/Mac (2)/Desktop/Post_Env/Olivine/20231214/run1_Det2500V_vel_gt8_pos4/ois_output_12152023_013217.h5"

    # with h5py.File(fname, "r") as f:

    #         for key in f.keys():
    #             # Access the datasets
    #             time_data = f[f"{key}/Time (high sampling)"][:]
    #             amplitude_data = f[f"{key}/TOF H"][:]
                
    #             # Convert the datasets to float
    #             time_data = np.array(time_data, dtype=float)
    #             amplitude_data = np.array(amplitude_data, dtype=float)
                
    #             # Append the converted data to the lists
    #             time.append(time_data*1e-6)
    #             amplitude.append(amplitude_data)

    # print(time[6])

    peak_time2mass(amplitude[0], time[0])

    # time2mass(amplitude[1], time[1])

    # time2mass(amplitude[2], time[2])
    # for i in range(len(TOF)):
        # time2mass(TOF[i], time[i])
