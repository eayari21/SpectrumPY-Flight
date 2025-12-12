#!/usr/bin/env python3
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# ------------------------------------------------------------
# User Settings
# ------------------------------------------------------------
INPUT_FILE = "HDF5/Flight/imap_idex_l0_raw_20251130_v001.h5"
EVENT = "29"

OUT_PLOT = "Plots/evt29_mass_spectrum_linear.png"

# Known IDEX mass lines (µs) — EDIT THESE FOR YOUR INSTRUMENT
# You can refine these from calibration tables.
mass_lines = {
    "H⁺":  1.90,
    "C⁺":  2.80,
    "O⁺":  3.20,
    "Mg⁺": 4.10,
    "Si⁺": 4.50,
    "Fe⁺": 5.10,
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def read_ds(grp, name):
    return np.array(grp[name]) if name in grp else None

def subtract_baseline(x, npre=200):
    baseline = np.median(x[:npre])
    return x - baseline, baseline

# ------------------------------------------------------------
# Load event data
# ------------------------------------------------------------
with h5py.File(INPUT_FILE, "r") as f:
    grp = f[EVENT]

    raw_h = read_ds(grp, "TOF H")
    raw_m = read_ds(grp, "TOF M")
    raw_l = read_ds(grp, "TOF L")

    time_hi = read_ds(grp, "Time (high sampling)")

# Convert channels to charge units
tof_h = (1/4.67) * raw_h
tof_m = 2.93e-3 * raw_m
tof_l = raw_l

# ------------------------------------------------------------
# Baseline subtract MID + LOW
# ------------------------------------------------------------
tof_m_bs, base_m = subtract_baseline(tof_m)
tof_l_bs, base_l = subtract_baseline(tof_l)

print("Baselines:")
print("  MID =", base_m)
print("  LOW =", base_l)

# ------------------------------------------------------------
# Gain-match LOW to MID
# ------------------------------------------------------------
# Derived gain ratio (event-by-event)
mask_pos = time_hi > 0
peak_m = np.max(tof_m_bs[mask_pos])
peak_l = np.max(tof_l_bs[mask_pos])
gain_l_over_m = peak_m / peak_l

# OR use your ~40 hardware gain
# gain_l_over_m = 1/40

print("\nGain correction factor (LOW → MID scale):", gain_l_over_m)

tof_l_corr = tof_l_bs * gain_l_over_m
tof_m_corr = tof_m_bs                         # MID is reference scale

# ------------------------------------------------------------
# Build the combined spectrum (overlay, no summing)
# ------------------------------------------------------------
plt.figure(figsize=(14, 6))
plt.plot(time_hi, tof_m_corr, label="TOF Mid (gain-corrected)", linewidth=1.2)
plt.plot(time_hi, tof_l_corr, label="TOF Low (gain-corrected)", linewidth=1.0)

plt.xlabel("Time (µs)")
plt.ylabel("Signal (baseline-subtracted, gain-corrected)")
plt.title("Combined TOF Spectrum (Mid + Low)\nFlat Baseline with Mass Lines")
plt.grid(True)
plt.ylim(bottom=-0.002)   # small negative to show baseline
plt.legend()

# ------------------------------------------------------------
# Peak detection for annotation
# ------------------------------------------------------------
# Combine MID and LOW for peak detection
combined = tof_m_corr + tof_l_corr

# Detect peaks over a threshold
peaks, _ = find_peaks(combined, height=np.max(combined)*0.1, distance=20)

# Annotate detected peaks
for idx in peaks:
    t = time_hi[idx]
    amp = combined[idx]
    plt.plot(t, amp, "ro", markersize=4)
    plt.text(t, amp + 0.0005, f"{t:.2f} µs", fontsize=8, rotation=45)

# ------------------------------------------------------------
# Label expected mass lines (from lookup table)
# ------------------------------------------------------------
for label, t_line in mass_lines.items():
    # find nearest index
    idx = np.argmin(np.abs(time_hi - t_line))
    amp = combined[idx]

    plt.axvline(t_line, color="k", linestyle="--", alpha=0.5)
    plt.text(
        t_line + 0.05,
        amp + 0.001,
        label,
        fontsize=10,
        color="k"
    )

# ------------------------------------------------------------
# Save + show
# ------------------------------------------------------------
plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=300)
plt.show()

print("\nSaved annotated mid+low mass spectrum:", OUT_PLOT)
