#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

"""
A Python object to store IDEX packets.
__author__      = Ethan Ayari & Gavin Medley, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""

# || Python libraries
import argparse
import json
import math
import os
import socket
import sys
import bitstring
import h5py
import shutil
import struct
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


def _float_or_nan(value: Optional[float]) -> float:
    """Return a floating point value or ``np.nan`` if conversion is not possible."""

    if value is None:
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")

if __package__ is None or __package__ == "":
    _MODULE_DIR = Path(__file__).resolve().parent
    _PACKAGE_ROOT = _MODULE_DIR.parent
    for _path in (_MODULE_DIR, _PACKAGE_ROOT):
        _path_str = str(_path)
        if _path_str not in sys.path:
            sys.path.append(_path_str)
    from plot_style import apply_plot_style
    from idex_analysis_utils import RISE_METRIC_SUFFIXES, compute_rise_metrics
    from rice_decode import idex_rice_Decode
    from time2mass import time2mass, get_last_mass_line_assignments
    from lookup.dust_estimator import estimate_particle, load_default_tables
else:
    from .plot_style import apply_plot_style
    from .idex_analysis_utils import RISE_METRIC_SUFFIXES, compute_rise_metrics
    from .rice_decode import idex_rice_Decode
    from .time2mass import time2mass, get_last_mass_line_assignments
    from .lookup.dust_estimator import estimate_particle, load_default_tables

apply_plot_style()
import numpy as np

MASS_STRETCH_MIN = 1.3
MASS_STRETCH_MAX = 1.6

COMBINED_SIGNAL_DATASET = "CombinedSignal"
COMBINED_TIME_DATASET = "CombinedTime"
DUST_ANALYSIS_GROUP = "Analysis/DustComposition"

try:
    import cupy as cp  # Optional GPU acceleration
    _HAS_CUPY = True
except Exception:  # pragma: no cover - cupy is optional
    cp = None
    _HAS_CUPY = False

from datetime import datetime, timedelta, timezone

from scipy.optimize import curve_fit
from scipy.signal import detrend, butter, filtfilt
from scipy.integrate import quad
from scipy.special import erfc


# || LASP software
from lasp_packets import xtcedef  # Gavin Medley's xtce UML implementation
from lasp_packets import parser  # Gavin Medley's constant bitstream implementation
import cdflib.cdfwrite as cdfwrite
import cdflib.cdfread as cdfread

# %%IDEX ION GRID FUNCTION DEFINITON
def IDEXIonGrid(x, P0, P1, P4, P5, P6):
    return P1 + np.heaviside(x-P0, 0) * ( P4 * (1.0 - np.exp(-(x-P0)/P5)) * np.exp( -(x-P0)/P6))

# Define the EMG function
def EMG(x, mu, sigma, lam):
    prefactor = lam / 2
    exponent = np.exp((lam / 2) * (2 * mu + lam * sigma**2 - 2 * x))
    erfc_part = erfc((mu + lam * sigma**2 - x) / (np.sqrt(2) * sigma))
    return prefactor * exponent * erfc_part

# Function to calculate the area under the EMG fit curve
def calculate_area_under_emg(x_slice, param):
    if(type(param) is not int):
        # Extract EMG fit parameters: mu, sigma, lam
        mu, sigma, lam = param

        # Perform numerical integration using scipy.integrate.quad
        area, error = quad(EMG, x_slice[0], x_slice[-1], args=(mu, sigma, lam))
        
        return area
    else:
        return 0.0

# Helper function to apply the polynomial transformation
def apply_polynomial(coeffs, X):
    # Compute the value using the polynomial formula
    return sum(coeffs[i] * (X ** i) for i in range(len(coeffs)))

# Helper function to create dataset if it doesn't exist
def create_dataset_if_not_exists(hdf5_file, dataset_path, data, *, dtype=None):
    if dataset_path in hdf5_file:
        print(f"Dataset '{dataset_path}' already exists. Skipping creation.")
        return hdf5_file[dataset_path]
    group_path = os.path.dirname(dataset_path)
    if group_path and group_path != '/':
        hdf5_file.require_group(group_path)
    if dtype is not None:
        return hdf5_file.create_dataset(dataset_path, data=data, dtype=dtype)
    return hdf5_file.create_dataset(dataset_path, data=data)


def _contiguous_mask(condition: np.ndarray, min_samples: int) -> np.ndarray:
    condition = np.asarray(condition, dtype=bool)
    if condition.size == 0:
        return np.zeros(0, dtype=bool)

    padded = np.concatenate(([False], condition, [False])).astype(int)
    diff = np.diff(padded)
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]

    mask = np.zeros_like(condition, dtype=bool)
    for start, end in zip(starts, ends):
        if end - start >= max(1, min_samples):
            mask[start:end] = True
    return mask


def _detect_saturation_segments(values: np.ndarray, times: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.zeros(0, dtype=bool)

    magnitude = np.nanmax(np.abs(arr))
    if not np.isfinite(magnitude) or magnitude == 0.0:
        return np.zeros_like(arr, dtype=bool)

    grad = np.abs(np.gradient(arr))
    derivative_threshold = 0.0025 * magnitude
    plateau = grad < derivative_threshold

    repeated = np.zeros_like(arr, dtype=bool)
    if arr.size >= 2:
        diffs = np.abs(np.diff(arr))
        repeat_tol = max(1.0e-9, 1.0e-4 * magnitude)
        repeats = diffs <= repeat_tol
        if repeats.any():
            repeated[1:] |= repeats
            repeated[:-1] |= repeats

    amplitude_threshold = np.nanpercentile(np.abs(arr), 99.7)
    high_amp = np.abs(arr) >= amplitude_threshold

    plateau_mask = (plateau | repeated) & high_amp

    extreme_mask = np.zeros_like(arr, dtype=bool)
    tolerance = 0.003 * magnitude + 1.0e-9
    max_val = float(np.nanmax(arr))
    min_val = float(np.nanmin(arr))
    if np.isfinite(max_val) and max_val > 0.0:
        extreme_mask |= (max_val - arr) <= tolerance
    if np.isfinite(min_val) and min_val < 0.0:
        extreme_mask |= (arr - min_val) <= tolerance
    plateau_mask |= extreme_mask & high_amp

    if plateau_mask.size < 2:
        return plateau_mask

    times = np.asarray(times, dtype=float)
    if times.size >= 2:
        dt = float(np.nanmedian(np.diff(times)))
    else:
        dt = 0.0

    if not np.isfinite(dt) or dt <= 0.0:
        min_samples = 12
    else:
        min_samples = max(8, int(math.ceil(1.0 / max(dt, 1.0e-6))))

    return _contiguous_mask(plateau_mask, min_samples)


def _first_microsecond_mean(values: Optional[np.ndarray], times: Optional[np.ndarray]) -> float:
    if values is None:
        return 0.0

    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0

    if times is None:
        sample_count = min(arr.size, 260)
        if sample_count == 0:
            return 0.0
        return float(np.nanmean(arr[:sample_count]))

    time_arr = np.asarray(times, dtype=float)
    if time_arr.size == 0:
        sample_count = min(arr.size, 260)
        if sample_count == 0:
            return 0.0
        return float(np.nanmean(arr[:sample_count]))

    length = min(arr.size, time_arr.size)
    if length == 0:
        return 0.0

    arr = arr[:length]
    time_arr = time_arr[:length]

    if length >= 2:
        dt = float(np.nanmedian(np.diff(time_arr)))
    else:
        dt = 0.0

    if np.isfinite(dt) and dt > 0.0 and dt < 1.0e-6:
        window = 1.0e-6
    else:
        window = 1.0

    start = float(time_arr[0])
    mask = (time_arr - start) <= window

    if not np.any(mask):
        if np.isfinite(dt) and dt > 0.0:
            samples = int(math.ceil(window / dt))
        else:
            samples = 260
        samples = min(max(samples, 1), length)
        mask = np.zeros(length, dtype=bool)
        mask[:samples] = True

    return float(np.nanmean(arr[mask]))


def _combine_waveform_channels(
    time_axis: np.ndarray,
    high: Optional[np.ndarray],
    medium: Optional[np.ndarray],
    low: Optional[np.ndarray],
    gain_map: Optional[Dict[str, float]] = None,
) -> Optional[np.ndarray]:
    del gain_map

    if time_axis is None:
        return None

    times = np.asarray(time_axis, dtype=float)
    if times.size == 0:
        return None

    candidates: List[Tuple[str, np.ndarray, int]] = []
    for label, channel in (("TOF H", high), ("TOF M", medium), ("TOF L", low)):
        if channel is None or not getattr(channel, "size", 0):
            continue
        arr = np.asarray(channel, dtype=float)
        length = min(arr.size, times.size)
        if length <= 0:
            continue
        candidates.append((label, arr[:length], length))

    if not candidates:
        return None

    for label, array, length in candidates:
        if label == "TOF H":
            sample_count = min(length, 50)
            if sample_count == 0:
                baseline = 0.0
            else:
                baseline = float(np.nanmean(array[:sample_count]))
            return array - baseline

    label, array, length = candidates[0]
    truncated_times = times[:length]
    baseline = _first_microsecond_mean(array, truncated_times)
    return array - baseline


def _estimate_baseline(time_array: np.ndarray, signal: np.ndarray) -> float:
    values = np.asarray(signal, dtype=float)
    if values.size == 0:
        return 0.0

    times = np.asarray(time_array, dtype=float)
    if times.size == values.size:
        mask = (times >= -7.0) & (times <= -5.0)
        if np.any(mask):
            return float(np.nanmedian(values[mask]))

    sample_count = min(values.size, 64)
    if sample_count == 0:
        return 0.0
    return float(np.nanmedian(values[:sample_count]))


def _serialise_mass_lines(group: h5py.Group, mass_lines: List[Dict[str, object]]) -> None:
    str_dtype = h5py.string_dtype(encoding='utf-8', length=120)
    extras_dtype = h5py.string_dtype(encoding='utf-8', length=2048)
    table = np.zeros(len(mass_lines), dtype=[
        ('id', 'i4'),
        ('label', str_dtype),
        ('assigned_species', str_dtype),
        ('mu', 'f8'),
        ('sigma', 'f8'),
        ('lam', 'f8'),
        ('amplitude', 'f8'),
        ('time_start', 'f8'),
        ('time_end', 'f8'),
        ('mass', 'f8'),
        ('assigned_mass', 'f8'),
        ('area', 'f8'),
        ('abundance', 'f8'),
        ('shape', str_dtype),
        ('extras', extras_dtype),
    ])
    for idx, record in enumerate(mass_lines):
        extras_serialized = "{}"
        try:
            extras_serialized = json.dumps(record.get('extras', {}))
        except Exception:
            extras_serialized = "{}"
        assigned_species = record.get('species', '') or ''
        assigned_mass = float(record.get('assigned_mass', np.nan))
        table[idx] = (
            int(record.get('line_id', idx + 1)),
            str(record.get('label', f"Line{idx + 1}")),
            str(assigned_species),
            float(record.get('mu', np.nan)),
            float(record.get('sigma', np.nan)),
            float(record.get('lam', np.nan)),
            float(record.get('amplitude', 0.0)),
            float(record.get('time_start', np.nan)),
            float(record.get('time_end', np.nan)),
            float(record.get('mass_guess', np.nan)),
            assigned_mass,
            float(record.get('area', 0.0)),
            float(record.get('abundance', 0.0)),
            str(record.get('shape', 'emg')),
            extras_serialized,
        )
    if 'MassLines' in group:
        del group['MassLines']
    group.create_dataset('MassLines', data=table)

    fits_group = group.require_group('Fits')
    for key in list(fits_group.keys()):
        del fits_group[key]
    for record in mass_lines:
        line_group = fits_group.require_group(f"line_{int(record.get('line_id', 0))}")
        for key in list(line_group.keys()):
            del line_group[key]
        line_group.create_dataset('time', data=np.asarray(record.get('fit_time', []), dtype=float))
        line_group.create_dataset('values', data=np.asarray(record.get('fit_curve', []), dtype=float))


def _analyse_mass_lines(signal: np.ndarray, time_axis: np.ndarray) -> Optional[Dict[str, object]]:
    signal = np.asarray(signal, dtype=float)
    time_axis = np.asarray(time_axis, dtype=float)
    if signal.size == 0 or signal.size != time_axis.size:
        return None

    stretch, shift, mass_scale = time2mass(signal, time_axis)
    stretch = float(np.clip(stretch, MASS_STRETCH_MIN, MASS_STRETCH_MAX))
    assignments = get_last_mass_line_assignments() or {}
    peaks = np.asarray(assignments.get('peaks', np.array([], dtype=int)), dtype=int)

    mass_line_records: List[Dict[str, object]] = []
    total_area = 0.0
    for line_info in assignments.get('mass_lines', []):
        peak_index = int(line_info.get('peak_index', 0))
        window = line_info.get('window', (peak_index - 10, peak_index + 10))
        start = max(0, int(window[0]))
        end = min(signal.size, int(window[1]))
        if end - start < 4:
            continue
        x_slice = np.asarray(time_axis[start:end], dtype=float)
        y_slice = np.asarray(signal[start:end], dtype=float)
        if x_slice.size == 0 or y_slice.size != x_slice.size:
            continue
        param, _param_cov, sig_amp, fitted_curve = FitEMG(x_slice, y_slice)
        if param is None:
            continue
        area = calculate_area_under_emg(x_slice, param)
        chi_sq, red_chi = calculate_chi_squared(y_slice, fitted_curve, len(param))
        mass_value = float(line_info.get('mass_reference', line_info.get('mass_scale_value', np.nan)))
        record = {
            'line_id': int(line_info.get('line_id', len(mass_line_records) + 1)),
            'label': str(line_info.get('label', f"Line{line_info.get('line_id', len(mass_line_records) + 1)}")),
            'species': str(line_info.get('species', '')),
            'mu': float(param[0]),
            'sigma': float(param[1]),
            'lam': float(param[2]),
            'amplitude': float(max(area, 0.0)),
            'time_start': float(x_slice[0]),
            'time_end': float(x_slice[-1]),
            'mass_guess': mass_value,
            'assigned_mass': mass_value,
            'area': float(max(area, 0.0)),
            'abundance': 0.0,
            'shape': 'emg',
            'extras': {},
            'fit_time': np.asarray(x_slice, dtype=float),
            'fit_curve': np.asarray(fitted_curve, dtype=float),
            'chi_sq': float(chi_sq),
            'red_chi': float(red_chi),
            'sig_amp': float(sig_amp),
            'peak_index': peak_index,
            'mass_scale_value': float(line_info.get('mass_scale_value', np.nan)),
        }
        total_area += record['area']
        mass_line_records.append(record)

    if total_area > 0.0:
        for record in mass_line_records:
            record['abundance'] = float(max(record['area'], 0.0) / total_area)
    else:
        for record in mass_line_records:
            record['abundance'] = 0.0

    valid_peaks = peaks[(peaks >= 0) & (peaks < len(mass_scale))]
    if valid_peaks.size:
        kappa = float(np.mean(mass_scale[valid_peaks] - np.round(mass_scale[valid_peaks])))
    else:
        kappa = np.nan

    return {
        'mass_scale': np.array(mass_scale, dtype=float),
        'peaks': valid_peaks,
        'kappa': float(kappa) if np.isfinite(kappa) else np.nan,
        'stretch': stretch,
        'shift': float(shift),
        'mass_lines': mass_line_records,
        'assignments': assignments,
        'total_area': float(total_area),
    }


def _compute_particle_estimate(
    charge_c: Optional[float],
    rise_time_us: Optional[float],
    ratio: Optional[float],
    *,
    rise_params,
    ratio_params,
    yield_params,
) -> Optional[object]:
    if charge_c is None or not np.isfinite(charge_c) or charge_c <= 0.0:
        return None
    if rise_params is None or yield_params is None:
        return None
    try:
        return estimate_particle(
            charge_c=charge_c,
            rise_time=rise_time_us,
            ion_to_target_ratio=ratio,
            rise_params=rise_params,
            ratio_params=ratio_params,
            yield_params=yield_params,
        )
    except Exception:
        return None


def _replace_data_dir(path: Path) -> Path:
    parts = list(path.parts)
    for idx, part in enumerate(parts):
        if part.lower() == 'data':
            parts[idx] = 'HDF5'
            base = Path(parts[0]) if parts else Path('HDF5')
            for segment in parts[1:]:
                base /= segment
            return base
    return path


def _resolve_output_path(filename: str) -> Path:
    input_path = Path(filename).expanduser()
    stem = input_path.stem if input_path.suffix else input_path.name
    parent = input_path.parent
    target_parent = _replace_data_dir(parent)
    if target_parent == parent and not target_parent.is_absolute():
        target_parent = Path(__file__).resolve().parent / "HDF5"
    target_parent.mkdir(parents=True, exist_ok=True)
    return target_parent / f"{stem}.h5"


def _bitstring_to_ints(waveform_raw: str, pad_bits: int, value_bits: int,
                       values_per_block: int, trim_tail: int = 0):
    if not waveform_raw:
        return []

    block_bits = pad_bits + value_bits * values_per_block
    ascii_bits = np.frombuffer(waveform_raw.encode('ascii'), dtype=np.uint8) - 48
    usable = (ascii_bits.size // block_bits) * block_bits
    if usable == 0:
        return []

    ascii_bits = ascii_bits[:usable].reshape(-1, block_bits)
    if pad_bits:
        ascii_bits = ascii_bits[:, pad_bits:]
    ascii_bits = ascii_bits.reshape(-1, value_bits).astype(np.int16, copy=False)

    powers = 1 << np.arange(value_bits - 1, -1, -1, dtype=np.int64)
    if _HAS_CUPY and ascii_bits.size > 32_000:
        try:
            ints_gpu = cp.asarray(ascii_bits)
            powers_gpu = cp.asarray(powers)
            ints = cp.asnumpy(ints_gpu.dot(powers_gpu))
        except Exception:
            ints = ascii_bits.dot(powers)
    else:
        ints = ascii_bits.dot(powers)

    if trim_tail and ints.size:
        ints = ints[:-trim_tail]
    return ints.tolist()


def calculate_chi_squared(observed, model, num_params):
    observed = np.asarray(observed, dtype=float)
    model = np.asarray(model, dtype=float)
    valid_mask = np.isfinite(observed) & np.isfinite(model)
    if not np.any(valid_mask):
        return np.nan, np.nan

    residuals = observed[valid_mask] - model[valid_mask]
    chi_sq = float(np.sum(residuals ** 2))
    dof = int(np.count_nonzero(valid_mask) - num_params)
    reduced_chi_sq = float(chi_sq / dof) if dof > 0 else np.nan
    return chi_sq, reduced_chi_sq


def calculate_snr(signal, time=None, baseline_range=(-7, -5)):
    signal = np.asarray(signal, dtype=float)
    if signal.size == 0:
        return np.nan

    if time is not None and len(time) == len(signal):
        time = np.asarray(time, dtype=float)
        baseline_mask = (time >= baseline_range[0]) & (time <= baseline_range[1])
    else:
        baseline_mask = np.zeros(signal.shape, dtype=bool)

    if not np.any(baseline_mask):
        baseline_mask = np.ones(signal.shape, dtype=bool)

    baseline_segment = signal[baseline_mask]
    if baseline_segment.size == 0:
        return np.nan

    baseline_mean = float(np.nanmean(baseline_segment))
    baseline_std = float(np.nanstd(baseline_segment))
    if baseline_std == 0 or np.isnan(baseline_std):
        return np.nan

    peak_amplitude = float(np.nanmax(signal) - baseline_mean)
    return peak_amplitude / baseline_std


def detect_saturation(signal, min_repeats=3):
    signal = np.asarray(signal)
    if signal.size < min_repeats:
        return False

    if _HAS_CUPY:
        try:
            signal_gpu = cp.asarray(signal)
            repeats = cp.diff(signal_gpu) == 0
            if not bool(cp.any(repeats)):
                return False
            boundaries = cp.nonzero(cp.diff(cp.concatenate((cp.array([False]), repeats, cp.array([False])))))[0]
            if boundaries.size == 0:
                return False
            run_lengths = cp.asnumpy(boundaries[1::2] - boundaries[::2])
            return run_lengths.size > 0 and (run_lengths.max() + 1) >= min_repeats
        except Exception:
            pass

    repeats = np.diff(signal) == 0
    if not np.any(repeats):
        return False
    changes = np.flatnonzero(np.diff(np.concatenate(([False], repeats, [False]))))
    if changes.size == 0:
        return False
    run_lengths = changes[1::2] - changes[::2]
    return run_lengths.size > 0 and (run_lengths.max() + 1) >= min_repeats

# Fit routine for EMG
def FitEMG(time, amplitude):
    x = np.asarray(time)
    y = np.asarray(amplitude)

    # || Initial Guess for the parameters of the EMG
    mu_guess = x[np.argmax(y)]  # Initial guess for the mean
    sigma_guess = np.std(x) / 10  # Initial guess for standard deviation
    lam_guess = 1 / (x[-1] - x[0])  # Initial guess for decay rate

    p0 = [mu_guess, sigma_guess, lam_guess]  # Initial parameter guesses

    # Fit the data using curve_fit
    try:
        param, param_cov = curve_fit(EMG, x, y, p0=p0, maxfev=100_000)

        # Generate the fitted curve
        result = EMG(x, *param)
        sig_amp = max(result) - np.mean(y)

        return param, param_cov, sig_amp, result
    except RuntimeError as e:
        print(f"Fit failed: {e}")
        return None, None, None, None

# %%Target Signal Fitting Routine %% #

# || Very noisy due to "microphonics", so we will:
# || 1) Remove a linear baseline (y = a*x + b), and 
# || 2) Remove a sinusoidal background (y = c*sin(d*x + e)

def _robust_sigma(x):
    """MAD-based robust std."""
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))

def _uniform_moving_avg(y, n):
    if n <= 1:
        return y.astype(float)
    c = np.cumsum(np.insert(y, 0, 0.0))
    out = (c[n:] - c[:-n]) / float(n)
    # pad to original length (reflect)
    pad_left = np.full(n // 2, out[0])
    pad_right = np.full(len(y) - len(out) - len(pad_left), out[-1])
    return np.concatenate([pad_left, out, pad_right]).astype(float)

def _find_onset_time(time, y, smooth_us=0.8, zthr=4.0):
    """
    Detect onset as the first index where smoothed dy/dt exceeds (median + zthr*robust_sigma).
    Returns np.nan if not found.
    """
    time = np.asarray(time, float)
    y = np.asarray(y, float)
    if len(time) < 5:
        return np.nan

    dt = np.median(np.diff(time))
    if not np.isfinite(dt) or dt <= 0:
        return np.nan

    # smooth over ~smooth_us
    n = max(3, int(round(abs(smooth_us / dt))))
    ys = _uniform_moving_avg(y, n)
    dy = np.gradient(ys, time)

    mu = np.median(dy)
    sig = _robust_sigma(dy)
    if sig <= 0 or not np.isfinite(sig):
        return np.nan

    z = (dy - mu) / sig
    idx = np.where(z > zthr)[0]
    return time[idx[0]] if idx.size else np.nan

def _quietest_window_mask(time, y, win_us=3.0, prefer_left_of=None):
    """
    Return a boolean mask selecting the quietest (lowest variance) sliding window of length win_us.
    If prefer_left_of is provided, prefer windows fully to the left of that time;
    fall back to global quietest if none exist.
    """
    time = np.asarray(time, float)
    y = np.asarray(y, float)
    n = len(time)
    if n == 0:
        return np.zeros(0, dtype=bool)

    dt = np.median(np.diff(time))
    k = max(3, int(round(abs(win_us / dt))))  # window length in samples
    if k >= n:
        m = np.zeros(n, dtype=bool)
        m[:] = True
        return m

    # compute rolling variance (simple, fast)
    # use cumulative sums for speed
    y2 = y * y
    c = np.cumsum(np.insert(y, 0, 0.0))
    c2 = np.cumsum(np.insert(y2, 0, 0.0))
    window_var = (c2[k:] - c2[:-k]) / k - ((c[k:] - c[:-k]) / k) ** 2

    # indices denote windows [i, i+k)
    candidates = np.arange(len(window_var))

    if prefer_left_of is not None:
        # only windows fully to the left of prefer_left_of
        left_mask = time[candidates + k - 1] < prefer_left_of
        left_candidates = candidates[left_mask]
        if left_candidates.size:
            i0 = left_candidates[np.argmin(window_var[left_mask])]
        else:
            i0 = candidates[np.argmin(window_var)]
    else:
        i0 = candidates[np.argmin(window_var)]

    m = np.zeros(n, dtype=bool)
    m[i0:i0 + k] = True
    return m

def FitTargetSignal(time, targetAmp,
                    pre_margin_us=2.0,   # left padding before onset for fit window
                    post_margin_us=60.0, # right padding after onset for fit window
                    baseline_win_us=3.0  # baseline window length
                    ):
    """
    Adaptive fit for the target signal. No hard-coded [-7, -5] µs gate.
    Returns (param, param_cov, sig_amp, time, filtered_full, fit_curve_full, chi_sq, red_chi)
    """

    # -- inputs as float arrays
    time = np.asarray(time, dtype=float)
    original_signal = np.asarray(targetAmp, dtype=float)
    signal = np.copy(original_signal)

    # Guard
    if time.size != signal.size or time.size == 0:
        # Return empty-like but consistent
        return {
            'params': np.array([]),
            'param_cov': np.empty((0, 0)),
            'signal_amplitude': np.nan,
            'time': time.astype(float),
            'filtered_signal': signal.astype(float),
            'fit_curve': np.full_like(signal, np.nan, dtype=float),
            'chi_sq': np.nan,
            'red_chi': np.nan,
            'rise_metrics': compute_rise_metrics([], []),
        }

    # Baseline guess (robust)
    baseline_guess = float(np.median(original_signal[:max(5, len(original_signal)//20)]))

    # -- Step 1: detect onset to guide masks
    t_onset = _find_onset_time(time, signal, smooth_us=0.8, zthr=4.0)

    # -- Step 2: choose a baseline/noise segment automatically
    # Prefer a quiet window *before* the onset; otherwise the global quietest window.
    if np.isfinite(t_onset):
        baseline_mask = _quietest_window_mask(time, signal, win_us=baseline_win_us, prefer_left_of=t_onset)
    else:
        baseline_mask = _quietest_window_mask(time, signal, win_us=baseline_win_us, prefer_left_of=None)

    baselineraw = signal[baseline_mask]
    baselinedomain = time[baseline_mask]

    # -- Step 3: remove linear background (fit only on baseline)
    try:
        slopeguess = 0.0
        # fit y = m*x + b on baseline, then detrend full signal using m,b
        (m_est, b_est), _ = curve_fit(LinearFit, baselinedomain, baselineraw,
                                      p0=[slopeguess, float(np.median(baselineraw))],
                                      maxfev=100_000)
        signal = signal - LinearFit(time, m_est, b_est)
    except Exception:
        # fallback: simple scipy detrend
        try:
            signal = detrend(signal)
        except Exception:
            # keep as-is
            pass

    # -- Step 4: optional sinusoidal background (fit on baseline region only)
    try:
        baseline_segment = signal[baseline_mask]
        if baseline_segment.size > 3:
            amp0 = float(np.ptp(baseline_segment)) if np.isfinite(np.ptp(baseline_segment)) else float(np.max(np.abs(baseline_segment)))
            amp0 = amp0 if np.isfinite(amp0) and amp0 > 0 else float(np.std(baseline_segment))
            p0 = [amp0, 1.0 / max(1e-6, np.median(np.diff(baselinedomain))), 0.0]  # [A, f, phi] crude init
            sineparam, _ = curve_fit(SineFit, baselinedomain, baseline_segment, p0=p0, maxfev=100_000)
            sinebase = SineFit(time, sineparam[0], sineparam[1], sineparam[2])
            signal = signal - sinebase
    except Exception:
        # ignore sinusoid removal if unstable
        pass

    # -- Step 5: low-pass filter (if available)
    try:
        signal = butter_lowpass_filter(signal, time)
    except Exception:
        pass

    filtered_full = np.asarray(signal, dtype=float)

    # -- Step 6: build an adaptive fit window around the onset
    if np.isfinite(t_onset):
        fit_left = t_onset - float(pre_margin_us)
        fit_right = t_onset + float(post_margin_us)
        fit_mask = (time >= fit_left) & (time <= fit_right)
        # safety: if the mask is tiny (e.g., onset near boundaries), expand
        if np.count_nonzero(fit_mask) < max(20, len(time)//50):
            fit_mask = np.ones_like(time, dtype=bool)
    else:
        # fallback to full domain
        fit_mask = np.ones_like(time, dtype=bool)

    fit_time = time[fit_mask]
    filtered_segment = filtered_full[fit_mask]

    # robust baseline within the fit window: use left-most 20% (or ≤ pre_margin_us) of the window
    if fit_time.size:
        left_span = min(pre_margin_us, 0.2 * (fit_time[-1] - fit_time[0]) if fit_time[-1] > fit_time[0] else pre_margin_us)
        base_mask_local = fit_time <= (fit_time[0] + left_span)
    else:
        base_mask_local = np.zeros(0, dtype=bool)

    yBaseline = np.where(base_mask_local, filtered_segment, np.nan)
    baseline_mean = float(np.nanmean(yBaseline)) if np.any(base_mask_local) else 0.0

    ionTime = fit_time.astype(float)
    ionAmp = filtered_segment.astype(float)

    # -- Step 7: parameter initial guesses for IDEXIonGrid
    # t0 near onset (if found) in the local coordinates
    t0 = float(ionTime[0]) if ionTime.size else 0.0
    if np.isfinite(t_onset):
        # set t0 close to onset but within window
        t0 = float(np.clip(t_onset, ionTime[0], ionTime[-1])) if ionTime.size else float(t_onset)

    # amplitude guess: difference between high percentile and baseline
    if ionAmp.size:
        hi = np.nanpercentile(ionAmp, 95)
        amplitude_guess = float(hi - baseline_mean)
        if not np.isfinite(amplitude_guess) or amplitude_guess <= 0:
            amplitude_guess = float(np.nanmax(ionAmp) - np.nanmin(ionAmp))
    else:
        amplitude_guess = 0.0

    # shape/time constants: keep your defaults but allow wider basin
    t1 = 3.71
    t2 = 37.1

    # baseline for model = baseline_mean (more stable than very first sample)
    baseline_for_model = baseline_mean if np.isfinite(baseline_mean) else baseline_guess

    # -- Step 8: fit
    try:
        param, param_cov = curve_fit(
            IDEXIonGrid,
            ionTime,
            ionAmp,
            p0=[t0, baseline_for_model, amplitude_guess, t1, t2],
            maxfev=100_000,
        )
    except Exception:
        # fall back: try without strict t0 (use window start) and smaller maxfev
        try:
            param, param_cov = curve_fit(
                IDEXIonGrid,
                ionTime,
                ionAmp,
                p0=[float(ionTime[0]) if ionTime.size else 0.0, baseline_for_model,
                    max(1e-6, amplitude_guess), t1, t2],
                maxfev=50_000,
            )
        except Exception:
            # give up gracefully with NaNs
            nanarr = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
            return {
                'params': nanarr,
                'param_cov': np.full((5, 5), np.nan),
                'signal_amplitude': np.nan,
                'time': time.astype(float),
                'filtered_signal': filtered_full,
                'fit_curve': np.full_like(filtered_full, np.nan, dtype=float),
                'chi_sq': np.nan,
                'red_chi': np.nan,
                'rise_metrics': compute_rise_metrics([], []),
            }

    # -- Step 9: compose full fit curve over the full time array (NaN outside fit window)
    fit_slice = IDEXIonGrid(ionTime, *param)
    fit_curve_full = np.full_like(filtered_full, np.nan, dtype=float)
    fit_curve_full[fit_mask] = fit_slice

    sig_amp = float(np.nanmax(fit_slice) - baseline_mean) if fit_slice.size else 0.0

    # goodness of fit on the fit window only (valid where model is defined)
    valid_mask = np.isfinite(fit_curve_full)
    residuals = filtered_full[valid_mask] - fit_curve_full[valid_mask]
    chi_sq = float(np.sum(residuals ** 2)) if residuals.size else np.nan
    dof = int(np.count_nonzero(valid_mask) - len(param))
    red_chi = float(chi_sq / dof) if dof > 0 and np.isfinite(chi_sq) else np.nan

    rise_metrics = compute_rise_metrics(time.astype(float), fit_curve_full, baseline_mean)
    return {
        'params': param,
        'param_cov': param_cov,
        'signal_amplitude': sig_amp,
        'time': time.astype(float),
        'filtered_signal': filtered_full,
        'fit_curve': fit_curve_full,
        'chi_sq': chi_sq,
        'red_chi': red_chi,
        'rise_metrics': rise_metrics,
    }

# ||
# ||
# || Generator object from LASP packets
# || to read in the data
class IDEXEvent:
    def __init__(self, filename: str):
        """Test parsing a real XTCE document"""
        # TODO: CHge location of xml definition
        module_root = Path(__file__).resolve().parent
        idex_xtce = module_root / "idex_combined_science_definition.xml"
        idex_definition = xtcedef.XtcePacketDefinition(xtce_document=str(idex_xtce))
        # assert isinstance(idex_definition, xtcedef.XtcePacketDefinition)


        idex_packet_file = filename
        print(f"Reading in data file {idex_packet_file}")
        idex_binary_data = bitstring.ConstBitStream(filename=idex_packet_file)
        print("Data import completed, writing packet structures.")

        idex_parser = parser.PacketParser(idex_definition)
        idex_packet_generator = idex_parser.generator(idex_binary_data,
                                                    # skip_header_bits=64,
                                                    skip_header_bits=32,  # For sciData
                                                    show_progress=True,
                                                    yield_unrecognized_packet_errors=True)
    

        print("Packet structures written.")
        idex_binary_data.pos = 0
        idex_packet_generator = idex_parser.generator(idex_binary_data)
        self.data = {}
        self.header={}
        self.lspretrigblocks = 0
        self.lsposttrigblocks = 0
        self.hspretrigblocks = 0
        self.hsposttrigblocks = 0
        self.hgdelay = 0
        self.hstime = np.array([], dtype=float)
        self.lstime = np.array([], dtype=float)
        evtnum = 0
        for pkt in idex_packet_generator:
            print(evtnum)
            if 'IDX__SCI0TYPE' in pkt.data:
                # print(evtnum)
                if pkt.data['IDX__SCI0TYPE'].raw_value == 1:
                    evtnum += 1
                    print(pkt.data)

                    # Iterate over all items in pkt.data and store them in the header
                    for key, item in pkt.data.items():
                        self.header[(evtnum, key)] = item.derived_value
                        print(f"{key} = {self.header[(evtnum, key)]}")
                    print(f"^*****Event header {evtnum}******^")

                    # sciEvtnum = bin(pkt.data['IDX__SCI0EVTNUM'].derived_value).replace('b', '')


                    # print(f"NBlocks = binary: {bin(pkt.data['IDX__TXHDRBLOCKS'].derived_value)} hex: {hex(pkt.data['IDX__TXHDRBLOCKS'].derived_value)}")
                    
                    # nBlocks = bin(pkt.data['IDX__TXHDRBLOCKS'].derived_value).replace('b', '')

                    # Extract the 17-22-bit integer (usually 8)
                    self.lspretrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value >> 16) &  0b1111

                    # Extract the next 4-bit integer (usually 8)
                    self.lsposttrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value >> 12) & 0b1111

                    # Extract the next 6 bits integer (usually 32)
                    self.hspretrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value >> 6) & 0b111111

                    # Extract the first 6 bits (usually 32)
                    self.hsposttrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value) & 0b111111


                    print("HS pre trig sampling blocks: ", self.hspretrigblocks)

                    print("LS pre trig sampling blocks: ", self.lspretrigblocks)

                    print("HS post trig sampling blocks: ", self.hsposttrigblocks)

                    print("LS post trig sampling blocks: ", self.lsposttrigblocks)

                    print(f"IDX__TXHDRHVPSHKCH01 = {pkt.data['IDX__TXHDRHVPSHKCH01'].derived_value}")

                    # Extract raw DN value for Voltage reading of Detector on HVPS Board (ADC CHnel 0)
                    self.header[(evtnum, 'detector_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH01'].derived_value) & 0b111111111111
                    print("Detector voltage = ", self.header[(evtnum, 'detector_voltage')])

                    # Extract raw DN value for Voltage reading of Sensor on HVPS Board (ADC CHnel 1)
                    self.header[(evtnum, 'sensor_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH01'].derived_value >> 16) & 0b111111111111
                    print("Sensor voltage = ", self.header[(evtnum, 'sensor_voltage')])

                    # HVPS Board signal "Target Voltage" (ADC CHnel 23)
                    self.header[(evtnum, 'target_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH23'].derived_value) & 0b111111111111
                    print("Target voltage = ", self.header[(evtnum, 'target_voltage')])

                    # HVPS Board signal "Reflectron Voltage" (ADC CHnel 23)
                    self.header[(evtnum, 'reflectron_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH23'].derived_value >> 16) & 0b111111111111
                    print("Reflectron voltage = ", self.header[(evtnum, 'reflectron_voltage')])

                    # HVPS Board signal "Rejection Voltage" (ADC CHnel 45)
                    self.header[(evtnum, 'rejection_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH45'].derived_value) & 0b111111111111
                    print("Rejection voltage = ", self.header[(evtnum, 'rejection_voltage')])

                    # HVPS Board signal "Current for the HVPS sensor" (ADC CHnel 45)
                    self.header[(evtnum, 'current_hvps_sensor')] = (pkt.data['IDX__TXHDRHVPSHKCH45'].derived_value >> 16) & 0b111111111111
                    print("Current for HVPS sensor = ", self.header[(evtnum, 'current_hvps_sensor')])

                    # HVPS Board signal "Positive current for the HVPS sensor" (ADC CHnel 67)
                    self.header[(evtnum, 'positive_current_hvps')] = (pkt.data['IDX__TXHDRHVPSHKCH67'].derived_value) & 0b111111111111
                    print("Positive current for HVPS sensor = ", self.header[(evtnum, 'positive_current_hvps')])

                    # HVPS Board signal "Negative current for the HVPS sensor" (ADC CHnel 67)
                    self.header[(evtnum, 'negative_current_hvps')] = (pkt.data['IDX__TXHDRHVPSHKCH67'].derived_value >> 16) & 0b111111111111
                    print("Negative current for HVPS sensor = ", self.header[(evtnum, 'negative_current_hvps')])

                    # LVPS Board signal "Voltage of +3.3V reference" (ADC CHnel 01)
                    self.header[(evtnum, 'voltage_3V3_ref')] = (pkt.data['IDX__TXHDRLVHK0CH01'].derived_value) & 0b111111111111
                    print("Voltage +3.3V reference = ", self.header[(evtnum, 'voltage_3V3_ref')])

                    # LVPS Board signal "Voltage of +3.3V operational reference" (ADC CHnel 01)
                    self.header[(evtnum, 'voltage_3V3_op_ref')] = (pkt.data['IDX__TXHDRLVHK0CH01'].derived_value >> 16) & 0b111111111111
                    print("Voltage +3.3V operational reference = ", self.header[(evtnum, 'voltage_3V3_op_ref')])

                    # LVPS Board signal "Voltage on -6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'voltage_neg6V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH23'].derived_value) & 0b111111111111
                    print("Voltage -6V bus = ", self.header[(evtnum, 'voltage_neg6V_bus')])

                    # LVPS Board signal "Voltage on +6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'voltage_pos6V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH23'].derived_value >> 16) & 0b111111111111
                    print("Voltage +6V bus = ", self.header[(evtnum, 'voltage_pos6V_bus')])

                    # LVPS Board signal "Voltage on +16V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'voltage_pos16V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH45'].derived_value) & 0b111111111111
                    print("Voltage +16V bus = ", self.header[(evtnum, 'voltage_pos16V_bus')])

                    # LVPS Board signal "Voltage on +3.3V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'voltage_pos3V3_bus')] = (pkt.data['IDX__TXHDRLVHK0CH45'].derived_value >> 16) & 0b111111111111
                    print("Voltage +3.3V bus = ", self.header[(evtnum, 'voltage_pos3V3_bus')])

                    # LVPS Board signal "Voltage on -5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_neg5V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH67'].derived_value) & 0b111111111111
                    print("Voltage -5V bus = ", self.header[(evtnum, 'voltage_neg5V_bus')])

                    # LVPS Board signal "Voltage on +5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_pos5V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH67'].derived_value >> 16) & 0b111111111111
                    print("Voltage +5V bus = ", self.header[(evtnum, 'voltage_pos5V_bus')])

                    # LVPS Board signal "Current on +3.3V bus" (ADC CHnel 01)
                    self.header[(evtnum, 'current_3V3_bus')] = (pkt.data['IDX__TXHDRLVHK1CH01'].derived_value) & 0b111111111111
                    print("Current +3.3V bus = ", self.header[(evtnum, 'current_3V3_bus')])

                    # LVPS Board signal "Current on +16V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'current_16V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH23'].derived_value >> 16) & 0b111111111111
                    print("Current +16V bus = ", self.header[(evtnum, 'current_16V_bus')])

                    # LVPS Board signal "Current on +6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'current_6V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH23'].derived_value) & 0b111111111111
                    print("Current +6V bus = ", self.header[(evtnum, 'current_6V_bus')])

                    # LVPS Board signal "Current on -6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'current_neg6V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH23'].derived_value >> 16) & 0b111111111111
                    print("Current -6V bus = ", self.header[(evtnum, 'current_neg6V_bus')])

                    # LVPS Board signal "Current on +5V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'current_5V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH45'].derived_value) & 0b111111111111
                    print("Current +5V bus = ", self.header[(evtnum, 'current_5V_bus')])

                    # LVPS Board signal "Current on -5V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'current_neg5V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH45'].derived_value >> 16) & 0b111111111111
                    print("Current -5V bus = ", self.header[(evtnum, 'current_neg5V_bus')])

                    # LVPS Board signal "Current on +2.5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'current_2V5_bus')] = (pkt.data['IDX__TXHDRLVHK1CH67'].derived_value) & 0b111111111111
                    print("Current +2.5V bus = ", self.header[(evtnum, 'current_2V5_bus')])

                    # LVPS Board signal "Current on -2.5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'current_neg2V5_bus')] = (pkt.data['IDX__TXHDRLVHK1CH67'].derived_value >> 16) & 0b111111111111
                    print("Current -2.5V bus = ", self.header[(evtnum, 'current_neg2V5_bus')])


                    # LVPS Board signal "Current on the 1V POL" (ADC CHnel 01)
                    self.header[(evtnum, 'current_1V_pol')] = (pkt.data['IDX__TXHDRPROCHKCH01'].derived_value) & 0b111111111111
                    print("Current on the 1V POL = ", self.header[(evtnum, 'current_1V_pol')])

                    # LVPS Board signal "Current on the 1.9V POL" (ADC CHnel 01)
                    self.header[(evtnum, 'current_1.9V_pol')] = (pkt.data['IDX__TXHDRPROCHKCH01'].derived_value >> 16) & 0b111111111111
                    print("Current on the 1.9V POL = ", self.header[(evtnum, 'current_1.9V_pol')])

                    # LVPS Board signal "ProcBd Temperature 1" (ADC CHnel 23)
                    self.header[(evtnum, 'temperature_1')] = (pkt.data['IDX__TXHDRPROCHKCH23'].derived_value) & 0b111111111111
                    print("ProcBd Temperature 1 = ", self.header[(evtnum, 'temperature_1')])

                    # LVPS Board signal "ProcBd Temperature 2" (ADC CHnel 23)
                    self.header[(evtnum, 'temperature_2')] = (pkt.data['IDX__TXHDRPROCHKCH23'].derived_value >> 16) & 0b111111111111
                    print("ProcBd Temperature 2 = ", self.header[(evtnum, 'temperature_2')])

                    # LVPS Board signal "Voltage on 1V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'voltage_1V_bus')] = (pkt.data['IDX__TXHDRPROCHKCH45'].derived_value) & 0b111111111111
                    print("Voltage on 1V bus = ", self.header[(evtnum, 'voltage_1V_bus')])

                    # LVPS Board signal "FPGA Temperature" (ADC CHnel 45)
                    self.header[(evtnum, 'fpga_temperature')] = (pkt.data['IDX__TXHDRPROCHKCH45'].derived_value >> 16) & 0b111111111111
                    print("FPGA Temperature = ", self.header[(evtnum, 'fpga_temperature')])

                    # LVPS Board signal "Voltage on 1.9V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_1.9V_bus')] = (pkt.data['IDX__TXHDRPROCHKCH67'].derived_value) & 0b111111111111
                    print("Voltage on 1.9V bus = ", self.header[(evtnum, 'voltage_1.9V_bus')])

                    # LVPS Board signal "Voltage on 3.3V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_3.3V_bus')] = (pkt.data['IDX__TXHDRPROCHKCH67'].derived_value >> 16) & 0b111111111111
                    print("Voltage on 3.3V bus = ", self.header[(evtnum, 'voltage_3.3V_bus')])






                    # Define the coefficients that are the same for every variable
                    coefficients = ['CO', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

                    mapping_dict = {
                        'detector_voltage': 'Last measurement in raw DN for HVPS Board signal “Detector Voltage”',
                        'sensor_voltage': 'Last measurement in raw DN for HVPS Board signal “Sensor Voltage"',
                        'target_voltage': 'Last measurement in raw DN for HVPS Board signal “Target Voltage”',
                        'reflectron_voltage': 'Last measurement in raw DN for HVPS Board signal “Reflectron Voltage”',
                        'rejection_voltage': 'Last measurement in raw DN for HVPS Board signal “Rejection Voltage”',
                        'current_hvps_sensor': 'Last measurement in raw DN for HVPS Board signal “Detector Current”',
                        'positive_current_hvps': 'Last measurement in raw DN for HVPS Board signal “Sensor IP”',
                        'negative_current_hvps': 'Last measurement in raw DN for HVPS Board signal “Sensor IN”',
                        'voltage_3V3_ref': 'Last measurement in raw DN for LVPS Board signal “P3.3VREF_HK”',
                        'voltage_3V3_op_ref': 'Last measurement in raw DN for LVPS Board signal “P3.3VREF_OP”',
                        'voltage_neg6V_bus': 'Last measurement in raw DN for LVPS Board signal “N6V”',
                        'voltage_pos6V_bus': 'Last measurement in raw DN for LVPS Board signal “P6V”',
                        'voltage_pos16V_bus': 'Last measurement in raw DN for LVPS Board signal “P16V”',
                        'voltage_pos3V3_bus': 'Last measurement in raw DN for LVPS Board signal “P3.3V”',
                        'voltage_neg5V_bus': 'Last measurement in raw DN for LVPS Board signal “N5V”',
                        'voltage_pos5V_bus': 'Last measurement in raw DN for LVPS Board signal “P5V”',
                        'current_3V3_bus': 'Last measurement in raw DN for LVPS Board signal “P3.3_IMON”',
                        'current_16V_bus': 'Last measurement in raw DN for LVPS Board signal “P16V_IMON”',
                        'current_6V_bus': 'Last measurement in raw DN for LVPS Board signal “P6V_IMON”',
                        'current_neg6V_bus': 'Last measurement in raw DN for LVPS Board signal “N6V_IMON”',
                        'current_5V_bus': 'Last measurement in raw DN for LVPS Board signal “P5V_IMON”',
                        'current_neg5V_bus': 'Last measurement in raw DN for LVPS Board signal “N5V_IMON”',
                        'current_2V5_bus': 'Last measurement in raw DN for LVPS Board signal “P2.5V_IMON”',
                        'current_neg2V5_bus': 'Last measurement in raw DN for LVPS Board signal “N2.5V_IMON”',
                        'spare_signal': 'Last measurement in raw DN for LVPS Board signal “Spare”',
                        'current_1V_pol':'Last measurement in raw DN for Processor Board signal “1V POL Current”',
                        'current_1.9V_pol':'Last measurement in raw DN for Processor Board signal “1.9V POL Current”',
                        'temperature_1': 'Last measurement in raw DN for Processor Board signal “ProcBd Temp1”',
                        'temperature_2': 'Last measurement in raw DN for Processor Board signal “ProcBd Temp2”',
                        'voltage_1V_bus': 'Last measurement in raw DN for Processor Board signal “1V Voltage”',
                        'fpga_temperature': 'Last measurement in raw DN for Processor Board signal “FPGA Temp”',
                        'voltage_1.9V_bus': 'Last measurement in raw DN for Processor Board signal “1.9V Voltage”',
                        'voltage_3.3V_bus': 'Last measurement in raw DN for Processor Board signal “3.3V Voltage”'
                    }

                    # Create a reverse dictionary
                    reverse_mapping_dict = {value: key for key, value in mapping_dict.items()}

                    # Read in Scott K's instrument settings conversions
                    settings_df = pd.read_excel("IDEX CDF Variable Definitions.xlsx")
                    # Normalize Var_notes by converting curly quotes to straight quotes and removing NaN values
                    # settings_df['Var_notes'] = settings_df['Var_notes'].replace(np.nan, '', regex=True)  # Replace NaN with empty strings
                    # settings_df['Var_notes'] = settings_df['Var_notes'].str.replace('“', '"').str.replace('”', '"')  # Replace curly quotes with straight quotes
                    # settings_df['Var_notes'] = settings_df['Var_notes'].str.strip()  # Remove any leading/trailing spaces

                    # print("Var notes: ", [note for note in settings_df["Var_notes"]])

                    # print("Coefficients: ", settings_df.columns)

                    # Step 1: Create the mapping from Var_notes to row indices
                    var_to_row = {var_note: idx for idx, var_note in enumerate(settings_df['Var_notes'])}

                    print(f"var_to_row = {var_to_row}")

                    print(f"\n \n ***** var_name being converted for each instrument setting ***** \n \n")

                    for var_name, row_idx in var_to_row.items():
                        try:
                            print("Matching row ", len(settings_df.iloc[[row_idx]].columns))  # Print the entire row for the problematic variable

                            # Extract the polynomial coefficients from the spreadsheet
                            coeffs = settings_df.iloc[row_idx][coefficients].values  # Access the row by index and then columns by labels
                            print(f"coeffs for {var_name}= {coeffs}")
                            target_value = var_name
                            print(f"var_name = {var_name}")


                            # Get the corresponding key
                            var_name = reverse_mapping_dict.get(target_value)
                            print(f"var_name = {var_name}")
                            
                            # Get the raw DN value for this variable (from your script)
                            X = self.header[(evtnum, var_name)]

                            print(f"Header info = {X}")
                            
                            # Apply the polynomial transformation
                            transformed_value = apply_polynomial(coeffs, X)
                            self.header[(evtnum, var_name)] = transformed_value
                            
                            print(f"Transformed value for {var_name} = {transformed_value}")
                        
                        except KeyError as e:
                            print(f"KeyError: {e} - Could not find entry for {var_name}, {row_idx}. Please check the Var_notes or variable mapping.")
                        except Exception as e:
                            print(f"An error occurred: {e}")


                     # Account for HS trigger delay
                    self.TOFdelay = pkt.data['IDX__TXHDRSAMPDELAY'].derived_value  # Last two bits are padding

                    # Mask to extract 10-bit values
                    mask = 0b1111111111

                    self.lgdelay = (self.TOFdelay) & mask # First 10 bits (0-9)
                    self.mgdelay = (self.TOFdelay >> 10) & mask # Next 10 bits (10-19)
                    self.hgdelay = (self.TOFdelay >> 20) & mask # Next 10 bits (20-29)
                    print(f"High gain delay = {self.hgdelay} samples.")
                    print(f"Mid gain delay = {self.mgdelay} samples.")
                    print(f"Low gain delay = {self.lgdelay} samples.")

                    if(pkt.data['IDX__TXHDRLSTRIGMODE'].derived_value!='DIS'):  # If this was a LS (Target Low Gain) trigger (DIS=disabled)
                        print(f"Low sampling trigger mode = {pkt.data['IDX__TXHDRLSTRIGMODE'].derived_value}")
                        self.Triggerorigin = 'LS' 
                        print("Low sampling trigger mode enabled.")
                        # Check the first 25th-bit integer for a coincidence trigger
                        # coincidence = (trigmode >> 24) &  0b1
                        # if(coincidence==1):
                            # self.Triggermode = 'Coincidence'
                        # else:
                            # self.Triggermode = 'Threshold'
                    print("High sampling trigger mode = ", pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value, pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value, pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value)
                    # Mask for extracting 11-bit and 10-bit values
                    mask_11_bit = 0b11111111111  # 11-bit mask
                    mask_10_bit = 0b1111111111   # 10-bit mask
                    if(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value!=0):
                        print("Low gain TOF trigger mode enabled.")
                        self.Triggerorigin = 'LG'
                        # Extract the first 11 bits (bits 21-31)
                        minsamples = pkt.data['IDX__TXHDRLGTRIGCTRL1'].derived_value & mask_11_bit
                        # Extract the second 11 bits (bits 10-20)
                        maxsamples = (pkt.data['IDX__TXHDRLGTRIGCTRL1'].derived_value >> 11) & mask_11_bit
                        # Extract the last 10 bits (bits 0-9)
                        self.header[(evtnum, 'TriggerLevel')] = 2.89e-4*((pkt.data['IDX__TXHDRLGTRIGCTRL1'].derived_value >> 22) & mask_10_bit)
                        print(f"Trigger level = {self.header[(evtnum, 'TriggerLevel')]}")

                        if(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value==1):
                            print("Threshold trigger mode enabled for low gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "LGThreshold"
                        elif(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value==2):
                            print("Single pulse mode enabled for low gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "LGSinglePulse"
                        else:
                            print("Double pulse mode enabled for low gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "LGDoublePulse"

                    if(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value!=0):
                        print("Mid gain TOF trigger mode enabled.")
                        self.Triggerorigin = 'MG'
                        # Extract the first 11 bits (bits 21-31)
                        minsamples = pkt.data['IDX__TXHDRMGTRIGCTRL1'].derived_value & mask_11_bit
                        # Extract the second 11 bits (bits 10-20)
                        maxsamples = (pkt.data['IDX__TXHDRMGTRIGCTRL1'].derived_value  >> 11) & mask_11_bit
                        # Extract the last 10 bits (bits 0-9)
                        self.header[(evtnum, 'TriggerLevel')] = 2.89e-4*((pkt.data['IDX__TXHDRMGTRIGCTRL1'].derived_value >> 22) & mask_10_bit)
                        print(f"Trigger level = {self.header[(evtnum, 'TriggerLevel')]}")
                        if(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value==1):
                            print("Threshold trigger mode enabled for mid gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "MGThreshold"
                        elif(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value==2):
                            print("Single pulse mode enabled for mid gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "MGSinglePulse"
                        else:
                            print("Double pulse mode enabled for mid gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "MGDoublePulse"

                    if(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value!=0):
                        print("High gain trigger mode enabled.")
                        self.Triggerorigin = 'HG'
                        # Extract the first 11 bits (bits 21-31)
                        minsamples = pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value & mask_11_bit
                        # Extract the second 11 bits (bits 10-20)
                        maxsamples = (pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value  >> 11) & mask_11_bit
                        # Extract the last 10 bits (bits 0-9)
                        self.header[(evtnum, 'TriggerLevel')] = 2.89e-4*((pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value >> 22) & mask_10_bit)
                        print(f"For {pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value}, HG Trigger level = {self.header[(evtnum, 'TriggerLevel')]}, sample settings = {minsamples}, {maxsamples}")

                        if(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value==1):
                            print("Threshold trigger mode enabled for high gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "HGThreshold"
                        elif(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value==2):
                            print("Single pulse mode enabled for high gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "HGSinglePulse"
                        else:
                            print("Double pulse mode enabled for high gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "HGDoublePulse"

                    print(f"AID = {pkt.data['IDX__SCI0AID'].derived_value}")  # Instrument event number
                    print(f"Event number = {pkt.data['IDX__SCI0EVTNUM'].raw_value}")  # Event number out of how many events constitute the file
                    # print(f"Time = {pkt.data['IDX__SCI0TIME32'].derived_value}")  # Time in 20 ns intervals


                    print(f"Rice compression enabled = {bool(pkt.data['IDX__SCI0COMP'].raw_value)}")
                    compressed = bool(pkt.data['IDX__SCI0COMP'].raw_value)  # If we need to decompress the data


                    # self.header[evtnum][f"TimeIntervals"] = pkt.data['IDX__SCI0TIME32'].derived_value  # Store the number of 20 us intervals in the respective CDF "Time" variables
                    self.header[(evtnum, 'Timestamp')] = pkt.data['SHCOARSE'].derived_value + 20*(10**(-6))*pkt.data['SHFINE'].derived_value # Use this as the CDF epoch
                    print(f"Timestamp = {self.header[(evtnum, 'Timestamp')]} seconds since epoch (Midnight January 1st, 2012)")

                    # Convert to MST (UTC-7)
                    utc_time = datetime(2010, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=self.header[(evtnum, 'Timestamp')])
                    # mst_offset = timedelta(hours=-7)
                    # mst_time = utc_time + mst_offset
                    print(f"Trigger time = {utc_time}")
                    self.header[(evtnum, 'Timestamp')] = utc_time.timestamp()


                if pkt.data['IDX__SCI0TYPE'].raw_value in [2, 4, 8, 16, 32, 64]:
                    # print(self.data.keys())

                    if (evtnum, pkt.data['IDX__SCI0TYPE'].raw_value) not in self.data.keys():  # If this is a new entry,
                        self.data.update({(evtnum, pkt.data['IDX__SCI0TYPE'].raw_value): pkt.data['IDX__SCI0RAW'].raw_value})
                    else:
                        self.data[(evtnum, pkt.data['IDX__SCI0TYPE'].raw_value)] += pkt.data['IDX__SCI0RAW'].raw_value


        # Parse the waveforms according to the scitype present (high gain and low gain CHnels encode waveform data differently).
        i = 1
        for scitype, waveform in self.data.items():
            if(compressed):  # If we need to decompress the data
                        print(waveform)
                        compressedFile = "test_compressed.txt"
                        dataFile = open(compressedFile, "wb")
                        index = 0
                        # print(f"||===waveform = {waveform}===||")
                        while index < len(waveform):
                            # Get 4 bytes (32 bits) from the 'waveform' binary string
                            data = waveform[index: index + 32]

                            # Convert the binary string to bytes using 'int' and 'to_bytes'
                            uint32 = int(data, 2).to_bytes(4, byteorder='big')

                            # Write the bytes to the file
                            dataFile.write(uint32)

                            index = index + 32

                        dataFile.close()
                        # print(waveform)
                        decompressor = RiceGolombDecompressor(waveform)
                        
                        if(scitype[1] < 12):  # LS
                            nsamples = 8*(self.lspretrigblocks + self.lsposttrigblocks)
                            # copy = gpt_rice_Decode(waveform, True, nsamples)
                            # copy = rice_Decode(compressedFile, f"test.txt", False, nsamples)
                            copy = decompressor.decompress(10)
                            waveform = copy
                        else:  # HS
                            nsamples = 512*(self.hspretrigblocks + self.hsposttrigblocks) # - pkt.data['IDX__TXHDRSAMPDELAY']
                            # copy = gpt_rice_Decode(waveform, True, nsamples)
                            # copy = rice_Decode(compressedFile, f"test.txt", True, nsamples)
                            copy = decompressor.decompress(12)
                            waveform = copy

            self.data[scitype] = parse_waveform_data(waveform, scitype[1])
        
        names = {2: "TOF H", 4: "TOF L", 8: "TOF M", 16: "Target L", 32: "Target H", 64: "Ion Grid"}
        datastore = {}
        for scitype, waveform in self.data.items():
            datastore[(scitype[0], names[scitype[1]])] = self.data[(scitype[0], scitype[1])]
        self.data = datastore
        self.numevents = evtnum
        # print(self.data.keys())

    # ||
    # ||
    # || Gather all of the events
    # || and plot them
    def _high_trigger_offset(self) -> float:
        pre_blocks = getattr(self, "lspretrigblocks", 0)
        delay = getattr(self, "hgdelay", 0)
        return 8 * (1.0 / 4.0625) * (pre_blocks + 1) - (1.0 / 260.0) * delay

    def _low_trigger_offset(self) -> float:
        pre_blocks = getattr(self, "hspretrigblocks", 0)
        return 512 * (1.0 / 260.0) * (pre_blocks + 1)

    def _build_time_array(self, sample_count: int, *, high_rate: bool) -> np.ndarray:
        if sample_count <= 0:
            return np.array([], dtype=float)
        spacing = 1.0 / 260.0 if high_rate else 1.0 / 4.0625
        offset = self._high_trigger_offset() if high_rate else self._low_trigger_offset()
        time_values = np.arange(sample_count, dtype=float) * spacing
        return time_values - offset

    def plot_all_data(self, packets, fname: str):
        fname = os.path.split(fname)[-1]
        # Create a folder to store the plots
        PlotFolder = os.path.join(os.getcwd(), f"Plots/{fname}")
        if os.path.exists(PlotFolder):  # If it exists, remove it
            shutil.rmtree(PlotFolder)
        os.makedirs(PlotFolder)

        # print("Number of packet items = ", len(packets.items()))
        fig, ax = plt.subplots(nrows=6)  # Make this general
        fig.set_size_inches(18.5, 10.5)
        for i, (k, v) in enumerate(packets.items()):  # k[0] = Event number, k[1] = CHnel name, v=waveform data
            # fig = plt.figure(figsize=(17,12)) 
            # print(i%6)
            i=i%6  # We take modulo 6 so it is the same for each event
            x = np.linspace(0, len(v), len(v))  # Number of samples
            # Scale number of samples by ~4 ns (high rate) or ~250 ns (low rate) to get to time.
            if(i<=2):
                x *= 1/260  # HS
                self.hstime = x
            else:
                x *= 1/4.0625  # LS
                self.lstime = x

            # print("array length = ", len(v))

            
            print(f"Length of Channel {k[1]} = {len(v)}")
            # ax[i].fill_between(x, v, color='r')
            ax[i].set_xlabel(r"Time [$\mu$ s]", font="Times New Roman", fontsize=30, fontweight='bold')
            if(i<=2):
                self.lstriggertime = 8*(1/4.0625)*(self.lspretrigblocks+1) - (1/260)*self.hgdelay
                print(f"High sampling trigger time = {self.lstriggertime} microseconds.")
                self.hstime = self.hstime - self.lstriggertime
                ax[i].axvline(min(self.hstime)+self.lstriggertime, c="red", lw=2)

            else:
                self.hstriggertime = 512*(1/260)*(self.hspretrigblocks+1)  #  - (1/260)*self.hgdelay
                print(f"Low sampling trigger time = {self.hstriggertime} microseconds.")
                self.lstime = self.lstime-self.hstriggertime
                ax[i].axvline(min(self.lstime)+self.hstriggertime, c="red", lw=2)
                
            plt.subplots_adjust(bottom=0.2)


            plt.suptitle(f"{fname} Event {k[0]}", font="Times New Roman", fontsize=30, fontweight='bold')
            # plt.tight_layout()

            if i==5:  #  End of the event, lets free up some memory
                ax[0].plot(self.hstime, packets[(k[0], "TOF L")])
                ax[0].set_ylabel("TOF L", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "TOF L")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "TOF L")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "TOF L")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "TOF L")])} [dN]'
                ax[0].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[0].transAxes)
                # ax[0].set_xlim([0, 31.5])
                
                ax[1].plot(self.hstime, packets[(k[0], "TOF M")])
                ax[1].set_ylabel("TOF M", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "TOF M")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "TOF M")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "TOF M")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "TOF M")])} [dN]'
                ax[1].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[1].transAxes)
                # ax[1].set_xlim([0, 31.5])
                
                ax[2].plot(self.hstime, packets[(k[0], "TOF H")])
                ax[2].set_ylabel("TOF H", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "TOF H")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "TOF H")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "TOF H")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "TOF H")])} [dN]'
                ax[2].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[2].transAxes)
                # ax[2].set_xlim([0, 31.5])

                ax[3].plot(self.lstime, packets[(k[0], "Ion Grid")])
                ax[3].set_ylabel("Ion Grid", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "Ion Grid")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Ion Grid")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Ion Grid")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Ion Grid")])} [dN]'
                ax[3].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[3].transAxes)
                # ax[3].set_xlim([0, 126.5])
                
                if(self.header[(k[0], 'Timestamp')] < 494_733_600):  # If we are before September 27th, 2023 then we use the old definitions
                
                    ax[4].plot(self.lstime, packets[(k[0], "Target L")])
                    ax[4].set_ylabel("Target LG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target L")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target L")])} [dN]'
                    ax[4].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[4].transAxes)
                    # ax[4].set_xlim([0, 126.5])
                    
                    ax[5].plot(self.lstime, packets[(k[0], "Target H")])
                    ax[5].set_ylabel("Target HG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target H")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target H")])} [dN]'
                    ax[5].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[5].transAxes)
                    # ax[5].set_xlim([0, 126.5])

                else:
                    ax[4].plot(self.lstime, packets[(k[0], "Target H")])
                    ax[4].set_ylabel("Target HG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target H")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target H")])} [dN]'
                    ax[4].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[4].transAxes)
                    # ax[4].set_xlim([0, 126.5])
                    
                    ax[5].plot(self.lstime, packets[(k[0], "Target L")])
                    ax[5].set_ylabel("Target LG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target L")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target L")])} [dN]'
                    ax[5].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[5].transAxes)
                    # ax[5].set_xlim([0, 126.5])


                plt.savefig(os.path.join(PlotFolder, f"{fname}_Event_{k[0]}.png"), dpi=100)
                # plt.show()
                plt.close()
                del fig, ax
                fig, ax = plt.subplots(nrows=6)  # Make this general
                fig.set_size_inches(17.5, 10.5)

                
    
    # ||
    # ||
    # || Write the waveform data 
    # || to an HDF5 file
    def write_to_hdf5(self, waveforms: dict, filename: str):
        output_path = _resolve_output_path(filename)
        if output_path.exists():
            output_path.unlink()

        try:
            rise_table, ratio_table, yield_table = load_default_tables()
            rise_params = rise_table.get('combined')
            ratio_params = ratio_table.get('combined') if ratio_table else None
            yield_params = yield_table.get('combined')
        except Exception as exc:
            print(f"Warning: unable to load dust estimator tables: {exc}")
            rise_params = None
            ratio_params = None
            yield_params = None

        conversion_factors = {
            'TOF H': 2.89e-4,
            'TOF M': 1.13e-2,
            'TOF L': 5.14e-4,
            'Ion Grid': 7.46e-4,
            'Target H': 1.63e-1,
            'Target L': 1.58e1,
        }
        target_channels = {'Target L', 'Target H', 'Ion Grid'}

        waveform_items = list(waveforms.items())

        high_sample_count = next(
            (len(data) for (_, channel), data in waveform_items if channel in {'TOF H', 'TOF M', 'TOF L'}),
            0,
        )
        low_sample_count = next(
            (len(data) for (_, channel), data in waveform_items if channel in {'Ion Grid', 'Target H', 'Target L'}),
            0,
        )

        self.hstime = self._build_time_array(high_sample_count, high_rate=True)
        self.lstime = self._build_time_array(low_sample_count, high_rate=False)

        event_saturation_flags: Dict[str, bool] = {}
        flags_by_event: Dict[str, Dict[str, List[str]]] = {}
        event_results: Dict[str, Dict[str, object]] = {}

        def _prepare_waveform(item):
            (event_id, channel), data = item
            event_key = str(event_id)
            analysis = {
                'key': (event_key, channel),
                'data': np.asarray(data),
                'transformed': None,
                'channel_saturated': False,
                'snr': None,
                'mass': None,
                'target_fit': None,
                'logs': [],
                'time_array': np.array([], dtype=float),
                'fit_failures': [],
                'notes': [],
                'baseline': 0.0,
                'charge_c': None,
                'rise_time': None,
            }

            if channel in conversion_factors:
                transformed_data = analysis['data'].astype(float) * conversion_factors[channel]
                analysis['transformed'] = transformed_data
                analysis['channel_saturated'] = detect_saturation(analysis['data'])

                if channel in {'TOF H', 'TOF M', 'TOF L'}:
                    time_array = self._build_time_array(len(transformed_data), high_rate=True)
                elif channel in target_channels:
                    time_array = self._build_time_array(len(transformed_data), high_rate=False)
                else:
                    time_array = np.arange(len(transformed_data), dtype=float)
                analysis['time_array'] = np.asarray(time_array, dtype=float)

                if channel == 'TOF H':
                    sample_count = min(transformed_data.size, 50)
                    if sample_count == 0:
                        baseline_value = 0.0
                    else:
                        baseline_value = float(np.nanmean(transformed_data[:sample_count]))
                else:
                    baseline_value = _estimate_baseline(analysis['time_array'], transformed_data)
                analysis['baseline'] = baseline_value

                if analysis['time_array'].size:
                    min_len = min(len(analysis['time_array']), len(transformed_data))
                    if min_len > 0:
                        snr_value = calculate_snr(
                            transformed_data[:min_len],
                            analysis['time_array'][:min_len],
                        )
                    else:
                        snr_value = calculate_snr(transformed_data)
                else:
                    snr_value = calculate_snr(transformed_data)
                analysis['snr'] = snr_value

                if channel == 'TOF H':
                    baseline_corrected = transformed_data - baseline_value
                    mass_data = _analyse_mass_lines(baseline_corrected, analysis['time_array'])
                    if mass_data is not None:
                        assignments = mass_data.get('assignments', {})
                        peaks = np.asarray(assignments.get('peaks', np.array([], dtype=int)), dtype=int)
                        analysis['logs'].append(
                            f"Auto-assigned {len(assignments.get('mass_lines', []))} mass line(s) with {peaks.size} detected peaks"
                        )
                    analysis['mass'] = mass_data
                else:
                    analysis['mass'] = None
            else:
                analysis['baseline'] = 0.0
                analysis['snr'] = calculate_snr(analysis['data'])

            if channel in target_channels:
                signal_for_fit = analysis['transformed'] if analysis['transformed'] is not None else analysis['data']
                if analysis['time_array'].size == signal_for_fit.size:
                    target_time = analysis['time_array']
                else:
                    target_time = self._build_time_array(signal_for_fit.size, high_rate=False)
                target_time = np.asarray(target_time, dtype=float)
                target_fit = FitTargetSignal(target_time, signal_for_fit)
                analysis['target_fit'] = target_fit
                params = np.asarray(target_fit.get('params', np.array([])), dtype=float)
                if params.size == 0 or not np.all(np.isfinite(params)):
                    analysis['fit_failures'].append(f"{channel}Fit")
                    analysis['notes'].append(f"Target fit for {channel} returned invalid parameters")
                analysis['time_array'] = target_time
                sig_amp = target_fit.get('signal_amplitude', np.nan)
                charge_c = float(sig_amp) * 1e-12 if np.isfinite(sig_amp) else None
                analysis['charge_c'] = charge_c
                rise_metrics = target_fit.get('rise_metrics', {})
                analysis['rise_time'] = rise_metrics.get('rise')
            else:
                analysis['target_fit'] = None

            return analysis

        max_workers = max(1, min(32, (os.cpu_count() or 1)))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            analyses = list(executor.map(_prepare_waveform, waveform_items))

        with h5py.File(output_path, 'w') as h:
            for (evtnum, key), value in self.header.items():
                if "IDX__" in key:
                    print(f"Skipping header field {key}")
                    continue
                dataset_path = f"/{evtnum}/Metadata/{key}"
                if isinstance(value, str):
                    dtype = h5py.string_dtype(encoding='utf-8')
                    create_dataset_if_not_exists(h, dataset_path, data=value, dtype=dtype)
                else:
                    create_dataset_if_not_exists(h, dataset_path, data=np.array(value))
                print(f"Created dataset: {dataset_path} with value: {value}")

            for analysis in analyses:
                event_key, channel = analysis['key']
                data = analysis['data']

                for log in analysis['logs']:
                    print(log)

                event_saturation_flags.setdefault(event_key, False)
                event_flags = flags_by_event.setdefault(
                    event_key,
                    {
                        'failed_fits': [],
                        'saturated_channels': [],
                        'notes': [],
                    },
                )

                event_info = event_results.setdefault(event_key, {'channels': {}, 'spice_written': False})
                event_info['channels'][channel] = analysis

                ra_values = np.deg2rad(np.random.uniform(0, 15, size=1))
                dec_values = np.deg2rad(np.random.uniform(-5, 5, size=1))
                roll_values = np.deg2rad(np.random.uniform(-0.5, 0.5, size=1))
                pitch_values = np.deg2rad(np.random.uniform(-0.5, 0.5, size=1))
                yaw_values = np.deg2rad(np.random.uniform(-0.5, 0.5, size=1))
                position_x = np.random.uniform(1.5e6, 1.52e6, size=1)
                position_y = np.random.uniform(-2000, 2000, size=1)
                position_z = np.random.uniform(-2000, 2000, size=1)
                velocity_x = np.random.uniform(-1.0, 1.0, size=1)
                velocity_y = np.random.uniform(-1.0, 1.0, size=1)
                velocity_z = np.random.uniform(-0.5, 0.5, size=1)

                if not event_info['spice_written']:
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/RightAscension", ra_values)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Declination", dec_values)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Attitude/Roll", roll_values)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Attitude/Pitch", pitch_values)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Attitude/Yaw", yaw_values)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Ephemeris/PositionX", position_x)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Ephemeris/PositionY", position_y)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Ephemeris/PositionZ", position_z)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Ephemeris/VelocityX", velocity_x)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Ephemeris/VelocityY", velocity_y)
                    create_dataset_if_not_exists(h, f"/{event_key}/SpiceData/Ephemeris/VelocityZ", velocity_z)
                    event_info['spice_written'] = True

                metadata_path = f"/{event_key}/Metadata/Epoch"
                if metadata_path not in h:
                    timestamp = self.header.get((int(event_key), 'Timestamp'), 0)
                    create_dataset_if_not_exists(h, metadata_path, data=np.array(timestamp))

                transformed = analysis['transformed']
                if transformed is not None:
                    dataset_path = f"/{event_key}/{channel}"
                    if dataset_path in h:
                        del h[dataset_path]
                    h.create_dataset(dataset_path, data=transformed)
                    channel_saturated = analysis['channel_saturated']
                    create_dataset_if_not_exists(
                        h,
                        f"/{event_key}/Metadata/{channel} Saturated",
                        data=np.array([int(channel_saturated)], dtype=np.int8),
                    )
                    if channel_saturated:
                        event_saturation_flags[event_key] = True
                        event_flags['saturated_channels'].append(channel)
                        event_flags['notes'].append(f"{channel} channel saturation detected")

                    if analysis['fit_failures']:
                        event_flags['failed_fits'].extend(analysis['fit_failures'])

                    if analysis['notes']:
                        event_flags['notes'].extend(analysis['notes'])

                    if analysis['snr'] is not None:
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel} SNR",
                            data=np.array([analysis['snr']], dtype=float),
                        )

                    if channel == 'TOF L':
                        time_data = analysis.get('time_array', self._build_time_array(len(data), high_rate=True))
                        h.create_dataset(f"/{event_key}/Time (high sampling)", data=time_data)
                    if channel == 'Ion Grid':
                        time_data = analysis.get('time_array', self._build_time_array(len(data), high_rate=False))
                        h.create_dataset(f"/{event_key}/Time (low sampling)", data=time_data)

                target_fit = analysis['target_fit']
                if target_fit is not None:
                    param = target_fit.get('params', np.array([]))
                    fit_time = target_fit.get('time', np.array([]))
                    fit_curve = target_fit.get('fit_curve', np.array([]))
                    filtered_signal = target_fit.get('filtered_signal', np.array([]))
                    sig_amp = target_fit.get('signal_amplitude', np.nan)
                    chi_sq = target_fit.get('chi_sq', np.nan)
                    red_chi = target_fit.get('red_chi', np.nan)

                    create_dataset_if_not_exists(h, f"/{event_key}/Analysis/{channel} Fit Parameters", data=np.asarray(param, dtype=float))
                    create_dataset_if_not_exists(h, f"/{event_key}/Analysis/{channel} Fit Time", data=np.asarray(fit_time, dtype=float))
                    create_dataset_if_not_exists(h, f"/{event_key}/Analysis/{channel} Fit Result", data=np.asarray(fit_curve, dtype=float))
                    create_dataset_if_not_exists(h, f"/{event_key}/Analysis/{channel} Filtered Signal", data=np.asarray(filtered_signal, dtype=float))
                    create_dataset_if_not_exists(h, f"/{event_key}/Analysis/{channel} Chi Squared", data=np.array([chi_sq], dtype=float))
                    create_dataset_if_not_exists(h, f"/{event_key}/Analysis/{channel} Reduced Chi Squared", data=np.array([red_chi], dtype=float))

                    rise_metrics = target_fit.get('rise_metrics', {})
                    for key, suffix in RISE_METRIC_SUFFIXES.items():
                        value = rise_metrics.get(key)
                        if value is None or not np.isfinite(value):
                            continue
                        dataset_path = f"/{event_key}/Analysis/{channel} {suffix}"
                        create_dataset_if_not_exists(h, dataset_path, data=np.array([value], dtype=float))

            for event_key, saturated in event_saturation_flags.items():
                create_dataset_if_not_exists(
                    h,
                    f"/{event_key}/Metadata/AnyChannelSaturated",
                    data=np.array([int(saturated)], dtype=np.int8),
                )

            string_dtype = h5py.string_dtype(encoding='utf-8')
            for event_key, flag_values in flags_by_event.items():
                flag_base = f"/{event_key}/Analysis/Flags"
                flag_group = h.require_group(flag_base)

                failed = sorted(set(flag_values['failed_fits']))
                saturated = sorted(set(flag_values['saturated_channels']))
                notes = sorted(set(flag_values['notes']))

                datasets = {
                    'FailedFits': failed,
                    'SaturatedChannels': saturated,
                    'Notes': notes,
                }

                for name, entries in datasets.items():
                    dataset_path = f"{flag_base}/{name}"
                    if dataset_path in h:
                        del h[dataset_path]
                    if entries:
                        data = np.array(entries, dtype=object)
                        flag_group.create_dataset(name, data=data, dtype=string_dtype)
                    else:
                        flag_group.create_dataset(name, shape=(0,), dtype=string_dtype)

            for event_key, info in event_results.items():
                channels = info.get('channels', {})
                high_analysis = channels.get('TOF H')
                medium_analysis = channels.get('TOF M')
                low_analysis = channels.get('TOF L')

                time_axis = None
                if high_analysis and high_analysis['time_array'].size:
                    time_axis = np.asarray(high_analysis['time_array'], dtype=float)
                elif medium_analysis and medium_analysis['time_array'].size:
                    time_axis = np.asarray(medium_analysis['time_array'], dtype=float)
                elif low_analysis and low_analysis['time_array'].size:
                    time_axis = np.asarray(low_analysis['time_array'], dtype=float)

                combined = _combine_waveform_channels(
                    time_axis,
                    high_analysis['transformed'] if high_analysis else None,
                    medium_analysis['transformed'] if medium_analysis else None,
                    low_analysis['transformed'] if low_analysis else None,
                )

                if combined is not None and time_axis is not None and time_axis.size:
                    combined = np.asarray(combined, dtype=float)
                    length = min(combined.size, time_axis.size)
                    if length == 0:
                        continue
                    combined = combined[:length]
                    combined_time = np.asarray(time_axis[:length], dtype=float)

                    dust_group = h.require_group(f"/{event_key}/{DUST_ANALYSIS_GROUP}")
                    if COMBINED_SIGNAL_DATASET in dust_group:
                        del dust_group[COMBINED_SIGNAL_DATASET]
                    dust_group.create_dataset(COMBINED_SIGNAL_DATASET, data=combined)
                    if COMBINED_TIME_DATASET in dust_group:
                        del dust_group[COMBINED_TIME_DATASET]
                    dust_group.create_dataset(COMBINED_TIME_DATASET, data=combined_time)

                    baseline_candidate = None
                    for candidate_channel in ('TOF H', 'TOF M', 'TOF L'):
                        candidate = channels.get(candidate_channel)
                        if candidate is None:
                            continue
                        baseline_value = candidate.get('baseline')
                        if baseline_value is not None and np.isfinite(baseline_value):
                            baseline_candidate = float(baseline_value)
                            break
                    if baseline_candidate is None or not np.isfinite(baseline_candidate):
                        baseline_candidate = 0.0
                    dust_group.attrs['Baseline'] = baseline_candidate

                    mass_data = _analyse_mass_lines(combined, combined_time)
                    if mass_data is not None:
                        dust_group.attrs['MassStretch'] = float(mass_data.get('stretch', np.nan))
                        dust_group.attrs['MassShift'] = float(mass_data.get('shift', np.nan))
                        mass_lines = mass_data.get('mass_lines', [])
                        if mass_lines:
                            _serialise_mass_lines(dust_group, mass_lines)
                        else:
                            if 'MassLines' in dust_group:
                                del dust_group['MassLines']
                            if 'Fits' in dust_group:
                                del dust_group['Fits']
                        event_group = h.require_group(event_key)
                        if 'Mass' in event_group:
                            del event_group['Mass']
                        event_group.create_dataset('Mass', data=mass_data['mass_scale'])
                    else:
                        dust_group.attrs['MassStretch'] = float('nan')
                        dust_group.attrs['MassShift'] = float('nan')

                ion_charge = None
                ion_analysis = channels.get('Ion Grid')
                if ion_analysis is not None:
                    ion_charge = ion_analysis.get('charge_c')

                for channel_name in ('Ion Grid', 'Target L', 'Target H'):
                    channel_analysis = channels.get(channel_name)
                    if channel_analysis is None:
                        continue
                    charge_c = channel_analysis.get('charge_c')
                    impact_value = float(charge_c) if charge_c is not None and np.isfinite(charge_c) else np.nan
                    create_dataset_if_not_exists(
                        h,
                        f"/{event_key}/Analysis/{channel_name} Impact Charge",
                        data=np.array([impact_value], dtype=float),
                    )

                    if channel_name in {'Target L', 'Target H'}:
                        rise_time = channel_analysis.get('rise_time')
                        ratio = None
                        if ion_charge is not None and charge_c is not None and charge_c > 0.0:
                            ratio = ion_charge / charge_c if charge_c else None
                        estimate = _compute_particle_estimate(
                            impact_value if np.isfinite(impact_value) and impact_value > 0.0 else None,
                            rise_time,
                            ratio,
                            rise_params=rise_params,
                            ratio_params=ratio_params,
                            yield_params=yield_params,
                        )
                        if estimate is not None:
                            mass_value = _float_or_nan(getattr(estimate, 'mass_kg', np.nan))
                            velocity_value = _float_or_nan(getattr(estimate, 'velocity_kms', np.nan))
                            yield_value = _float_or_nan(getattr(estimate, 'yield_c_per_kg', np.nan))
                            velocity_details = getattr(estimate, 'velocity_details', None)
                            if velocity_details is not None:
                                rise_velocity = _float_or_nan(getattr(velocity_details, 'rise_time', np.nan))
                                ratio_velocity = _float_or_nan(
                                    getattr(velocity_details, 'collection_efficiency', np.nan)
                                )
                                velocity_source = getattr(velocity_details, 'source', '') or ''
                            else:
                                rise_velocity = np.nan
                                ratio_velocity = np.nan
                                velocity_source = ''
                        else:
                            mass_value = np.nan
                            velocity_value = np.nan
                            yield_value = np.nan
                            rise_velocity = np.nan
                            ratio_velocity = np.nan
                            velocity_source = ''
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Dust Mass Estimate",
                            data=np.array([mass_value], dtype=float),
                        )
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Velocity Estimate",
                            data=np.array([velocity_value], dtype=float),
                        )
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Charge Yield Estimate",
                            data=np.array([yield_value * 1.0e12], dtype=float),
                        )
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Velocity From Rise",
                            data=np.array([rise_velocity], dtype=float),
                        )
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Velocity From Ratio",
                            data=np.array([ratio_velocity], dtype=float),
                        )
                        if velocity_source:
                            source_dataset = [velocity_source]
                            create_dataset_if_not_exists(
                                h,
                                f"/{event_key}/Analysis/{channel_name} Velocity Source",
                                data=source_dataset,
                                dtype=h5py.string_dtype(encoding='utf-8'),
                            )
                    else:
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Dust Mass Estimate",
                            data=np.array([np.nan], dtype=float),
                        )
                        create_dataset_if_not_exists(
                            h,
                            f"/{event_key}/Analysis/{channel_name} Velocity Estimate",
                            data=np.array([np.nan], dtype=float),
                        )

# ||
# ||
# || Parse the high sampling rate data, this
# || should be 10-bit blocks
def parse_hs_waveform(waveform_raw: str):
    """Parse a binary string representing a high gain waveform"""
    ints = _bitstring_to_ints(waveform_raw, pad_bits=2, value_bits=10, values_per_block=3, trim_tail=4)
    print(len(ints))
    return ints

# ||
# ||
# || Parse the low sampling rate data, this
# || should be 12-bit blocks
def parse_ls_waveform(waveform_raw: str):
    """Parse a binary string representing a low gain waveform"""
    ints = _bitstring_to_ints(waveform_raw, pad_bits=8, value_bits=12, values_per_block=2)
    print(len(ints))
    return ints

# ||
# ||
# || Use the SciType flag to determine the sampling rate of
# || the data we are trying to parse
def parse_waveform_data(waveform: str, scitype: int):
    """Parse the binary string that represents a waveform"""
    print(f'Parsing waveform for scitype={scitype}')
    if scitype in (2, 4, 8):
        return parse_hs_waveform(waveform)
    else:
        return parse_ls_waveform(waveform)

# ||
# ||
# || Write the waveform data 
# || to CDF files
def write_to_cdf(packets):
    
    cdf_master = cdfread.CDF('imap_idex_l0-raw_0000000_v01.cdf')
    if (cdf_master.file != None):
    # Get the cdf's specification
        info=cdf_master.cdf_info()
        cdf_file=cdfwrite.CDF('./IDEX_SSIM.cdf',cdf_spec=info,delete=True)
    # if (cdf_file.file == None):
    #     cdf_master.close()
    #     raise OSError('Problem writing file.... Stop')

    # Get the global attributes
    globalaAttrs=cdf_master.globalattsget(expand=True)
    # Write the global attributes
    cdf_file.write_globalattrs(globalaAttrs)
    zvars=info['zVariables']
    print('no of zvars=',len(zvars))
    # Loop thru all the zVariables --> What are zvars vs rvars?
    for x in range (0, len(zvars)):
        # Get the variable's specification
        varinfo=cdf_master.varinq(zvars[x])
        print('Z =============>',x,': ', varinfo['Variable'])


# Z =============> 0 :  Epoch
# Z =============> 1 :  IDEX_Trigger
# Z =============> 2 :  TOF_Low
# Z =============> 3 :  TOF_Mid
# Z =============> 4 :  TOF_High
# Z =============> 5 :  Time_Low_SR
# Z =============> 6 :  Time_High_SR
# Z =============> 7 :  Target_Low
# Z =============> 8 :  Target_High
# Z =============> 9 :  Ion_Grid

        if(varinfo['Variable']=="Epoch"):
            vardata = None
        if(varinfo['Variable']=="IDEX_Trigger"):
            vardata = packets.header[(1,"Timestamp")]
        if(varinfo['Variable']=="TOF_Low"):
            print(len(np.array(packets.data[(1,"TOF L")])))
            vardata = np.array(packets.data[(1,"TOF L")], float)
        if(varinfo['Variable']=="TOF_Mid"):
            vardata = np.array(packets.data[(1,"TOF M")])
        if(varinfo['Variable']=="TOF_High"):
            vardata = np.array(packets.data[(1,"TOF H")])
        if(varinfo['Variable']=="Time_Low_SR"):
            vardata = np.linspace(0, len(packets.data[(1,"Target L")]), len(len(packets.data[(1,"Target L")])))
        if(varinfo['Variable']=="Time_High_SR"):
            vardata = np.linspace(0, len(packets.data[(1,"TOF L")]), len(len(packets.data[(1,"Target L")])))
        if(varinfo['Variable']=="Target_Low"):
            vardata = np.array(packets.data[(1,"Target L")])
        if(varinfo['Variable']=="Target_High"):
            vardata = np.array(packets.data[(1,"Target H")])
        if(varinfo['Variable']=="Ion_Grid"):
            vardata = np.array(packets.data[(1,"Ion Grid")])
        # Get the variable's attributes
        varattrs=cdf_master.varattsget(zvars[x], expand=True)
        if (varinfo['Sparse'].lower() == 'no_sparse'):
            # A variable with no sparse records... get the variable data
            # vardata= None
            # Create the zVariable, write out the attributes and data
            cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
        else:
            # A variable with sparse records...
            # data is in this form [physical_record_numbers, data_values]
            # physical_record_numbers (0-based) contains the real record
            # numbers. For example, a variable has only 3 physical records
            # at [0, 5, 10]:
            varrecs=[0,5,10]

            # vardata=None  # np.asarray([.,.,.,..])
            # Create the zVariable, and optionally write out the attributes
            # and data
            cdf_file.write_var(varinfo, var_attrs=varattrs,
                        var_data=[varrecs,vardata])
    rvars=info['rVariables']
    print('no of rvars=',len(rvars))
    # Loop thru all the rVariables
    for x in range (0, len(rvars)):
        varinfo=cdf_master.varinq(rvars[x])
        print('R =============>',x,': ', varinfo['Variable'])
        varattrs=cdf_master.varattsget(rvars[x], expand=True)
        if (varinfo['Sparse'].lower() == 'no_sparse'):
            vardata=None
            # Create the rVariable, write out the attributes and data
            cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
        else:
            varrecs= None  # [.,.,.,..]
            vardata= None  # np.asarray([.,.,.,..])
            cdf_file.write_var(varinfo, var_attrs=varattrs,
                        var_data=[varrecs,vardata])
    cdf_master.close()
    cdf_file.close()

# || Test code: Import file and write the relevant data to an hdf5 file
if __name__ == "__main__":
    # Initalize parsing object to pass filename
    aparser = argparse.ArgumentParser()
    aparser.add_argument("--file", "-f", type=str, required=True)
    args = aparser.parse_args()

    packets = IDEXEvent(args.file)
    # print(packets.data.keys())
    try:
        packets.plot_all_data(packets.data, args.file)
    except Exception as e:
        print(e)
    packets.write_to_hdf5(packets.data, args.file)
    # write_to_cdf(packets)
