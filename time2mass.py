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

from matplotlib import colors
from scipy.optimize import curve_fit
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
    """After the baseline correction:

    1) Make a vector with all zeros and a length of 8189, same as the TOF length of record: t_i
    2) Make a vector for masses m_i: 1, 2, 3, 4, …… 200 (200 length, and the elements are the integers 1,2, 3,  ... 200)
    3) Calculate the 200 times (T_i-s)  these masses should show up for a T_off = 0  and stretch A (best guess)
              T_i = T_off + A*sqrt(m_i) and set the closest t_i =1 to each of the calculated Ti-s, the rest stays zero.
    4) Calculate the cross-correlation with the original TOF, the max will give you the best lag (T_off) for a given A.
    5) Use an extrema finder to get to the best A, alternatively, one could try optimizing for the two parameters simultaneously. """

    print(os.listdir())  # Try to find the mass comb file
    mass_comb = pd.read_csv("../mass_comb.csv")

    masses = np.round(mass_comb["Mass"],1)

    time = time - time[0]
    # plt.plot(time, TOF, label="Time Scale")
    # plt.show()

    random_stretch = np.linspace(1400, 1500, 10)
    shift = []
    corr = []
    product = []

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
        cross_correlation = np.correlate(t_i, TOF, mode='full')

        # lags = np.linspace(1, len(cross_correlation), len(cross_correlation))
        # plt.scatter(lags, cross_correlation)
        # plt.xlabel("Lags", fontsize=14, fontweight="bold")
        # plt.ylabel("Correlation", fontsize=14, fontweight="bold")
        # plt.title(f"Correlation for stretch = {stretch}")
        # plt.show()

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

    mass_scale = [((time*1e-6-samplesize*shift[np.argmax(corr)])/(1e-9*random_stretch[np.argmax(corr)]))**2 for time in time]
    # mass_scale = [((time-samplesize*9416)/(1e-9*1444))**2 for time in time]
    # mass_scale = [((time*1e3)/(random_stretch[np.argmax(corr)]))**2 for time in time]

    # print("mass scale max/min = ",min(mass_scale), " , ", max(mass_scale))
    # print(f"&&&&&& mass scale = {mass_scale}")


    # best_lag is the lag that aligns T_i with TOF for the maximum correlation

    # || Step 5): Use an extrema finder to get to the best A, alternatively, one could try optimizing for the two parameters simultaneously.



    # fig, ax1 = plt.subplots()

    # # Plotting TOF against mass_scale on the main axis
    # ax1.plot(mass_scale, TOF, 'b-')
    # ax1.set_xlabel('Mass Scale [amu]', fontsize=14)
    # ax1.set_ylabel('TOF [dN]', fontsize=14)

    # # Create a secondary x-axis for time
    # secax = ax1.secondary_xaxis('top')
    # secax.set_xlabel('Time [s]', fontsize=14)

    # # Define tick positions and compute corresponding time values
    # best_shift = shift[np.argmax(corr)]*1e-9
    # best_stretch = random_stretch[np.argmax(corr)]*1e-9

    # # Define tick positions and compute corresponding time values
    # tick_positions = np.linspace(min(mass_scale), max(mass_scale), num=4)  # Adjust num as needed
    # tick_labels = [f'{np.sqrt(t * best_stretch**2 + best_shift**2):.2e}' for t in tick_positions]


    # secax.set_xticks(tick_positions)
    # secax.set_xticklabels(tick_labels)

    # plt.title('TOF vs Mass Scale and Time')
    # plt.grid(True)
    # plt.show()

    # || Heatap for stretch and shift

    # print('\n'.join(['*'] * 20))
    # print(f"stretch = {random_stretch}, \n shift = {shift}, \n corr = {corr}")
    random_stretch = np.array(random_stretch)
    shift = np.array(shift)
    # corr_matrix = np.array(corr)
    
    # Initialize corr_values array
    corr_values = np.zeros((len(random_stretch), len(shift)))

    # Compute correlation values
    for i in range(len(random_stretch)):
        for j in range(len(shift)):

            # Convert lists to arrays for calculations
            temp_mass_scale = np.array([(time - samplesize * shift[j]) / (1e-9 * random_stretch[i])**2 for time in time])
            
            # Ensure mass_scale is also a NumPy array
            mass_scale_array = np.array(mass_scale)
            
            # Check that temp_mass_scale and mass_scale have the same length
            if len(temp_mass_scale) != len(mass_scale_array):
                raise ValueError(f"Length mismatch: temp_mass_scale has length {len(temp_mass_scale)}, but mass_scale has length {len(mass_scale_array)}")

            # Compute the correlation value
            corr_values[i, j] = np.sum(temp_mass_scale - mass_scale_array)


    # Plotting
    # plt.figure(figsize=(10, 8))
    # plt.imshow(corr_values, extent=[shift.min(), shift.max(), random_stretch.min(), random_stretch.max()],
    #         origin='lower', cmap='magma', aspect='auto')
    # plt.colorbar(label='Correlation Value')

    # plt.ylabel(r'Stretch $\alpha$ [ns]', fontsize=20, fontweight="bold")
    # plt.xlabel('Shift $t_{0}$ [ns]', fontsize=20, fontweight="bold")
    # plt.title('Correlation Map', fontsize=20, fontweight="bold")
    # plt.show()
    return random_stretch[np.argmax(corr)], samplesize*shift[np.argmax(corr)], mass_scale_array

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
