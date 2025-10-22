#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
An algorithm to convert the time axis to mass for TOF spectra.
__author__      = Ethan Ayari, Mihály Horányi,
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""
import os
import h5py
import random
import csv

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns

from pathlib import Path
from matplotlib import colors
from scipy.optimize import curve_fit, minimize
from scipy.signal import find_peaks
# plt.style.use('seaborn-v0_8-pastel')
plt.style.use('seaborn-pastel')
# ||
# ||
# ||

# plt.rcParams['font.family'] = 'DejaVuSerif-BoldItalic'
# plt.rcParams['font.serif'] = 'DejaVuSerif-BoldItalic'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['agg.path.chunksize'] = 10_000

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

def time2mass(TOF, time):
    r"""Return an optimized mass axis for a TOF waveform.

    The routine estimates the stretch (``A``) and temporal offset (``t0``)
    parameters in the classic

    .. math:: t = t_0 + A \sqrt{m}

    relationship by aligning prominent peaks with integer mass numbers from
    :datafile:`mass_comb.csv`. The optimization keeps the resulting mass scale
    inside a physically realistic 1–300 amu window (allowing a small buffer up
    to 400 amu for the heaviest species).
    """

    mass_table = pd.read_csv(Path(__file__).resolve().with_name("mass_comb.csv"))
    available_masses = np.array(sorted({int(round(m)) for m in mass_table["Mass"] if 0 < m < 400}))

    if len(available_masses) == 0:
        raise ValueError("No valid masses found in mass_comb.csv")

    tof = np.asarray(TOF, dtype=float)
    time = np.asarray(time, dtype=float)

    time_offset = time - time[0]
    step = np.median(np.diff(time_offset)) if len(time_offset) > 1 else 0.0

    # Determine the native time units (seconds or microseconds).
    if step > 1e-6:
        scale_to_seconds = 1e-6
    else:
        scale_to_seconds = 1.0

    time_seconds = time_offset * scale_to_seconds
    time_ns = time_seconds * 1e9

    # Identify the most prominent peaks in the spectrum.
    prominence = max(np.ptp(tof) * 0.05, 1.0)
    peaks, _ = find_peaks(tof, prominence=prominence)

    if len(peaks) < 3:
        # Attempt a secondary search with a lower threshold.
        peaks, _ = find_peaks(tof)

    if len(peaks) == 0:
        # Fallback to a default mapping if no peaks are detected.
        default_stretch_ns = 1450.0
        mass_scale = _compute_mass_axis(time_ns, default_stretch_ns, 0.0)
        return default_stretch_ns, 0.0, mass_scale

    # Focus on the most intense peaks to drive the calibration.
    peak_order = np.argsort(tof[peaks])[-min(len(peaks), 20):]
    peak_times_ns = np.sort(time_ns[peaks][peak_order])

    def objective(params):
        stretch_ns, shift_ns = params
        if stretch_ns <= 0:
            return 1e9

        shifted = peak_times_ns - shift_ns
        if np.any(shifted <= 0):
            return 1e9

        masses = (shifted / stretch_ns) ** 2

        # Penalize masses outside the desired operating range.
        penalty = 0.0
        penalty += np.sum((1.0 - masses[masses < 1.0]) ** 2)
        penalty += np.sum((masses[masses > 350.0] - 350.0) ** 2)

        # Match peaks to the nearest known (integer) mass numbers.
        diffs = masses[:, None] - available_masses[None, :]
        nearest = available_masses[np.argmin(np.abs(diffs), axis=1)]
        residual = masses - nearest

        # Encourage a realistic mass span across the entire waveform.
        full_shifted = time_ns - shift_ns
        if np.any(full_shifted > 0):
            mass_full = (full_shifted[full_shifted > 0] / stretch_ns) ** 2
            penalty += max(0.0, mass_full.max() - 400.0) ** 2
            penalty += max(0.0, 1.0 - mass_full.min()) ** 2
        else:
            penalty += 1e6

        return float(np.mean(residual ** 2) + 0.01 * penalty)

    stretch_guess = 1450.0
    shift_guess = max(0.0, peak_times_ns[0] - stretch_guess * np.sqrt(max(available_masses[0], 1)))

    bounds = [(800.0, 4000.0), (-500.0, max(time_ns[-1], 1000.0))]
    result = minimize(objective, x0=[stretch_guess, shift_guess], bounds=bounds, method="L-BFGS-B")

    if not result.success:
        best_stretch_ns, best_shift_ns = stretch_guess, shift_guess
    else:
        best_stretch_ns, best_shift_ns = result.x

    mass_scale = _compute_mass_axis(time_ns, best_stretch_ns, best_shift_ns)

    # Convert the temporal shift to seconds for downstream consumers.
    shift_seconds = best_shift_ns * 1e-9

    return float(best_stretch_ns), float(shift_seconds), mass_scale


def _compute_mass_axis(time_ns: np.ndarray, stretch_ns: float, shift_ns: float) -> np.ndarray:
    """Calculate the mass scale ensuring it remains in the physical range."""

    shifted = time_ns - shift_ns
    masses = np.empty_like(time_ns, dtype=float)
    valid = shifted > 0
    masses[valid] = (shifted[valid] / stretch_ns) ** 2
    masses[~valid] = 0.0

    masses = np.clip(masses, 0.0, 400.0)
    return masses

    # ||
    # ||
    # ||
    # %%

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
