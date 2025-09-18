import h5py
import numpy as np
from time2mass import time2mass
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-pastel")

# Function to extract TOF H values where Time (High Sampling) is between -7 and -5
def extract_snr(hdf_file):
    TOFmaxes = []
    sigmas = []
    kappas = []
    stretches = []
    shifts = []

    # Open the HDF5 file
    with h5py.File(hdf_file, 'r') as f:
        # Loop through all the numbered groups
        for group_key in f.keys():
            # print(group_key)
            # Access each group (event number)
            group = f[group_key]
            # print(group.keys())
            # Check if "TOF H" and "Time (High Sampling)" datasets are present
            if "TOF H" in group.keys() and "Time (high sampling)" in group.keys():
                print(f"Extracting baselines for event {group_key}")
                # Read the "TOF H" and "Time (High Sampling)" datasets
                tof_h = group["TOF H"][:]
                time_high_sampling = group["Time (high sampling)"][:]

                # Calculate the mass scale for the TOF spectrum
                stretch, shift, mass_scale = time2mass(tof_h, time_high_sampling)
                print(f"Time scale: max={max(time_high_sampling)}, min={min(time_high_sampling)}")
                print(f"Mass scale: max={max(mass_scale)}, min={min(mass_scale)}")
                print(f"Sample size = {time_high_sampling[1] - time_high_sampling[0]} microseconds")

                # Find indices of all peaks (local maxima) and troughs (local minima)
                peaks, _ = find_peaks(tof_h, prominence=20)
                print(f"peaks = {peaks}")

                kappa = np.mean([mass_scale[peak]-np.round(mass_scale[peak], 1) for peak in peaks])

                # kappa = np.mean(np.abs([mass_scale[peak]-np.round(mass_scale[peak], 1) for peak in peaks]))
                print(f"Kappa = {kappa}")


                # Find indices where Time (High Sampling) is between -7 and -5
                mask = np.logical_and(time_high_sampling >= -7, time_high_sampling <= -5)

                # Append the corresponding TOF H values to the baseline array
                TOFmaxes.append(max(tof_h) - np.mean(tof_h[mask]))
                sigmas.append(np.std(tof_h[mask]))
                kappas.append(kappa)
                stretches.append(stretch)
                shifts.append(shift)

                # # Create subplots for the histograms
                # fig, ax = plt.subplots(1, 1, figsize=(8, 12))

                # # Plot histogram for SNRs
                # ax.hist(tof_h[mask], bins=20, color='b', alpha=0.7)
                # ax.set_xlabel("Noise amplitude [dN]", fontsize=20, font="Times New Roman", fontweight="bold")
                # ax.set_ylabel("Count", fontsize=20, font="Times New Roman", fontweight="bold")
                # ax.set_title("Noise distribution", fontsize=22, font="Times New Roman", fontweight="bold")

    return np.array(TOFmaxes), np.array(sigmas), np.array(kappas), np.array(shifts), np.array(stretches)

# hdf_file = 'HDF5/ois_output_12182023_205615.h5'
hdf_file = 'HDF5/ois_output_12182023_184030.h5'

maxes, sigmas, kappas, shifts, stretches = extract_snr(hdf_file)
print(f"sigmas {len(sigmas)}, maxes = {len(maxes)}")

SNRs = [int(m/s) for m, s in zip(maxes, sigmas)]
print(len(SNRs))

plt.scatter(SNRs, kappas)
plt.xlabel("SNR", fontsize = 20, font="Times New roman", fontweight="bold")
plt.ylabel(r"Kappa $\kappa$", fontsize = 20, font="Times New roman", fontweight="bold")
plt.show()

# Define bin size and calculate number of bins
bin_size = 40
num_bins = len(SNRs) // bin_size

# Sort SNRs and kappas together for binning
sorted_data = sorted(zip(SNRs, kappas), key=lambda x: x[0])
sorted_SNRs, sorted_kappas = zip(*sorted_data)

# Create lists to hold bin means and standard deviations
bin_centers = []
kappa_means = []
kappa_stds = []

# Bin the data
for i in range(num_bins):
    bin_SNRs = sorted_SNRs[i * bin_size : (i + 1) * bin_size]
    bin_kappas = sorted_kappas[i * bin_size : (i + 1) * bin_size]
    
    # Calculate mean SNR for bin center, and mean & std of kappa values in bin
    bin_centers.append(np.mean(bin_SNRs))
    kappa_means.append(np.mean(bin_kappas))
    kappa_stds.append(np.std(bin_kappas))

# Plot the binned data with error bars
plt.errorbar(bin_centers, kappa_means, yerr=kappa_stds, fmt='o', color="skyblue", 
             ecolor="gray", elinewidth=2, capsize=.1)
plt.xlabel("SNR", fontsize=20, font="Times New Roman", fontweight="bold")
plt.ylabel(r"Kappa $\kappa$", fontsize=20, font="Times New Roman", fontweight="bold")
plt.title("Binned SNR vs. Kappa with Error Bars", fontsize=20, font="Times New Roman", fontweight="bold")
plt.show()

# Create subplots for the histograms
fig, axs = plt.subplots(3, 1, figsize=(8, 12))

# Plot histogram for SNRs
axs[0].hist(SNRs, bins=20, color='b', alpha=0.7)
axs[0].set_xlabel("SNR", fontsize=20, font="Times New Roman", fontweight="bold")
axs[0].set_ylabel("Count", fontsize=20, font="Times New Roman", fontweight="bold")
axs[0].set_title("SNR Distribution", fontsize=22, font="Times New Roman", fontweight="bold")

# Plot histogram for Stretches
axs[1].hist(stretches, bins=20, color='g', alpha=0.7)
axs[1].set_xlabel("Stretch", fontsize=20, font="Times New Roman", fontweight="bold")
axs[1].set_ylabel("Count", fontsize=20, font="Times New Roman", fontweight="bold")
axs[1].set_title("Stretch Distribution", fontsize=22, font="Times New Roman", fontweight="bold")

# Plot histogram for Shifts
axs[2].hist(shifts, bins=20, color='r', alpha=0.7)
axs[2].set_xlabel("Shift", fontsize=20, font="Times New Roman", fontweight="bold")
axs[2].set_ylabel("Count", fontsize=20, font="Times New Roman", fontweight="bold")
axs[2].set_title("Shift Distribution", fontsize=22, font="Times New Roman", fontweight="bold")

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()



