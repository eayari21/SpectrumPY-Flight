#!/usr/bin/env python3
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erf

# ------------------------------------------------------------
# User Settings
# ------------------------------------------------------------
INPUT_FILE = "HDF5/Flight/imap_idex_l0_raw_20251130_v001.h5"
EVENT = "29"
OUT_PLOT = "Plots/evt29_combined_single_waveform.png"

LOW_GAIN_FACTOR = 40
SATURATION_FRAC = 0.90

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def read_ds(grp, name):
    return np.array(grp[name]) if name in grp else None

def subtract_baseline(x, npre=200):
    baseline = np.median(x[:npre])
    return x - baseline, baseline

# Exponentially Modified Gaussian
def emg(x, A, mu, sigma, lam):
    """A * ExpModGauss."""
    # EMG = convolution of Gaussian + exponential
    return A * (lam/2) * np.exp((lam/2)*(2*mu + lam*sigma**2 - 2*x)) * \
           (1 - erf((mu + lam*sigma**2 - x) / (np.sqrt(2)*sigma)))

# ------------------------------------------------------------
# Load event data
# ------------------------------------------------------------
with h5py.File(INPUT_FILE, "r") as f:
    grp = f[EVENT]

    raw_m = read_ds(grp, "TOF M")
    raw_l = read_ds(grp, "TOF L")
    raw_h = read_ds(grp, "TOF H")
    time_hi = read_ds(grp, "Time (high sampling)")

# Convert
tof_m = 2.93e-3 * raw_m
tof_l = raw_l
tof_h = (1/4.67) * raw_h

# ------------------------------------------------------------
# Baseline subtract
# ------------------------------------------------------------
tof_m_bs, base_m = subtract_baseline(tof_m)
tof_l_bs, base_l = subtract_baseline(tof_l)
tof_h_bs, base_h = subtract_baseline(tof_h)

print("Baselines:")
print("  Mid =", base_m)
print("  Low =", base_l)
print("  High =", base_h)

# Scale low → mid
tof_l_scaled = tof_l_bs * LOW_GAIN_FACTOR

# ------------------------------------------------------------
# Identify saturation in MID channel
# ------------------------------------------------------------
m_peak = np.max(tof_m_bs)
sat_threshold = SATURATION_FRAC * m_peak
sat_mask = tof_m_bs >= sat_threshold

print(f"\nM-peak = {m_peak:.4f}, saturation threshold = {sat_threshold:.4f}")
print(f"Saturated samples: {np.sum(sat_mask)}")

# ------------------------------------------------------------
# Fit EMG to the M peak
# ------------------------------------------------------------
# Use region near the main M peak (±1.5 µs)
peak_idx = np.argmax(tof_m_bs)
t0 = time_hi[peak_idx]
fit_mask = (time_hi > t0 - 1.5) & (time_hi < t0 + 1.5)
x_fit = time_hi[fit_mask]
y_fit = tof_m_bs[fit_mask]

# Initial guess
p0 = [np.max(y_fit), t0, 0.3, 1.0]

# Parameter bounds (keeps fit stable)
bounds = (
    [0, t0 - 1.0, 0.01, 0.01],
    [10*m_peak, t0 + 1.0, 2.0, 10.0]
)

try:
    popt, _ = curve_fit(emg, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=20000)
    print("EMG fit parameters:", popt)
except Exception as e:
    print("EMG fit failed:", e)
    popt = p0

# Build the EMG curve across entire time axis
emg_curve = emg(time_hi, *popt)

# Subtract EMG
tof_m_emg_sub = tof_m_bs - emg_curve

# ------------------------------------------------------------
# 1. Fix baseline droop to the right (use late baseline)
# ------------------------------------------------------------
right_region = time_hi > 6.0  # adjust threshold if needed
post_base = np.median(tof_m_emg_sub[right_region])
tof_m_emg_sub -= post_base   # restore baseline to ~0

# ------------------------------------------------------------
# 2. Fix the non-physical dip around 1.28–1.78 µs
# ------------------------------------------------------------
dip_region = (time_hi >= 1.1) & (time_hi <= 1.78)

# option A — use real M for this region (cleanest)
tof_m_emg_sub[dip_region] = .3*tof_m_bs[dip_region]

# option B — clamp to baseline instead
# baseline_val = np.median(tof_m_bs[time_hi < -5])
# tof_m_emg_sub[dip_region] = np.maximum(
# tof_m_emg_sub[dip_region], baseline_val

# ------------------------------------------------------------
# 3. Fit & subtract logistic S-curve tail for t >= 0 µs
# ------------------------------------------------------------

def logistic_tail(t, K, N0, r):
    # Standard logistic curve:
    #   K / (1 + ((K-N0)/N0) * exp(-r t))
    return K / (1.0 + ((K - N0) / N0) * np.exp(-r * t))

# Fit region for the S-curve (monotonic settling region)
fit_mask = (time_hi >= 2.0) & (time_hi <= 10.0)
t_fit = time_hi[fit_mask]

# Fit the logistic tail to raw EMG-sub (before baseline fix)
raw_emg_sub = tof_m_bs - emg_curve
y_fit = raw_emg_sub[fit_mask]

# Initial guesses
K0  = np.max(y_fit)
N00 = y_fit[0] if y_fit[0] != 0 else 1e-4
r0  = 1.0

p0 = [K0, N00, r0]
bounds = (
    [0,    1e-6,  0.0001],
    [0.1,  0.01,  10.0]
)

try:
    popt_log, _ = curve_fit(logistic_tail, t_fit, y_fit, p0=p0, bounds=bounds)
    print("Logistic tail parameters:", popt_log)
except:
    popt_log = p0
    print("Logistic tail fit failed, using defaults.")

# Build logistic curve across all time
log_curve = logistic_tail(time_hi, *popt_log)

# Subtract only for t >= 0
mask = time_hi >= 0
tof_m_emg_sub[mask] -= log_curve[mask]

# Rebaseline to 0.002
desired = 0.002
current = np.median(tof_m_emg_sub[(time_hi >= 10) & (time_hi <= 15)])
tof_m_emg_sub += 1/6*(desired - current)



# ------------------------------------------------------------
# 4. Fit & subtract a small parabola from t = 1.3 to 9.5 µs
# ------------------------------------------------------------

def parabola(t, a, b, c):
    return a*t*t + b*t + c

# region to fit
para_mask = (time_hi >= 1.3) & (time_hi <= 9.5)
t_para = time_hi[para_mask]
y_para = tof_m_emg_sub[para_mask]

# initial guesses should be tiny — this is a *small* curvature
p0 = [1e-5, -1e-4, 0.0]

try:
    popt_para, _ = curve_fit(parabola, t_para, y_para, p0=p0)
    print("Parabola params:", popt_para)
except:
    popt_para = p0
    print("Parabola fit failed; using initial guesses.")

# evaluate across whole axis
para_curve = parabola(time_hi, *popt_para)

# subtract only in the correction region
tof_m_emg_sub[para_mask] -= .5*para_curve[para_mask]

# ------------------------------------------------------------
# Build combined waveform
# ------------------------------------------------------------
combined = np.copy(tof_m_bs)

# Use LOW where mid saturates
combined[sat_mask] = tof_l_scaled[sat_mask]

# ------------------------------------------------------------
# Forced EMG-sub windows
# ------------------------------------------------------------
w1 = (time_hi >= 1.22) & (time_hi <= 2.24)
w2 = (time_hi >= 2.55) & (time_hi <= 15.95)

# Start with EMG-sub everywhere
combined = np.copy(tof_m_emg_sub)

# Windows where we DO NOT use EMG-sub
w1 = (time_hi >= .67) & (time_hi <= 1.3)
w2 = (time_hi >= 2.55) & (time_hi <= 15.95)
window_mask = w1 | w2

# Inside these windows, default to MID channel
combined[window_mask] = tof_m_bs[window_mask]

# Saturation overrides MID channel
combined[~(window_mask & sat_mask)] = tof_m_emg_sub[~(window_mask & sat_mask)]
tof_m_emg_sub[w1] = tof_l_bs[w1]
tof_m_emg_sub[(time_hi >=-20) & (time_hi <= .715)] = tof_m_bs[(time_hi >=-20) & (time_hi <= .715)]


# ------------------------------------------------------------
# Detect recovery droop (same as before)
# ------------------------------------------------------------
recovery_threshold = 0.30 * m_peak
peak_index = np.argmax(tof_m_bs)
rec_start = peak_index + np.where(tof_m_bs[peak_index:] < recovery_threshold)[0][0]
rec_end = rec_start + np.where(tof_m_bs[rec_start:] >= 0)[0][0]

recovery_mask = np.zeros_like(combined, dtype=bool)
recovery_mask[rec_start:rec_end] = True

print(f"Recovery region: {rec_start} → {rec_end}")

combined[recovery_mask] = tof_l_scaled[recovery_mask]


# ------------------------------------------------------------
# 5. Robust detection of ~15 peaks (no misses)
# ------------------------------------------------------------
from scipy.signal import find_peaks

sig = tof_m_emg_sub
t = time_hi

region_mask = (t >= 1.2) & (t <= 9.5)
t_reg = t[region_mask]
y_reg = sig[region_mask]

# --- Adaptive thresholding ---
# baseline = median of quiet region
baseline = np.median(y_reg[(t_reg > 7.5)])

# noise estimate from MAD
mad = np.median(np.abs(y_reg - baseline))
noise = 1.4826 * mad

# Peak prominence threshold ~ 6 sigma
prom_thresh = 6 * noise

# find many peaks (prominence-based is key)
peak_idx, props = find_peaks(
    y_reg,
    prominence=prom_thresh,
    width=1,
    distance=10  # minimal separation in samples
)

peak_times = t_reg[peak_idx]
proms = props["prominences"]

# sort peaks by time
order = np.argsort(peak_times)
peak_times = peak_times[order]
proms = proms[order]

# --- enforce non-overlap (±0.5 µs window) ---
selected = []
WINDOW = 0.50  # µs

for pt in peak_times:
    if all(abs(pt - s) > WINDOW for s in selected):
        selected.append(pt)

# If still < 15, fill in by choosing next highest-prominence peaks
if len(selected) < 15:
    # choose remaining candidates sorted by prominence
    remaining = sorted(
        [(pt, pr) for pt, pr in zip(peak_times, proms)
         if pt not in selected],
        key=lambda x: -x[1]
    )
    for pt, _ in remaining:
        if len(selected) == 15:
            break
        if all(abs(pt - s) > WINDOW for s in selected):
            selected.append(pt)

# If too many, truncate cleanly to 15
selected = sorted(selected)[:15]

print("\nFinal selected peaks:")
for i, pt in enumerate(selected):
    print(f"  Peak {i+1}: {pt:.4f} µs")

# ------------------------------------------------------------
# EMG fit to each peak (strict 1 µs window)
# ------------------------------------------------------------

def fit_emg_window(t_all, y_all, center):
    w = 0.5
    mask = (t_all >= center - w) & (t_all <= center + w)
    x = t_all[mask]
    y = y_all[mask]
    if len(x) < 25:
        return None

    A0 = np.max(y)
    mu0 = x[np.argmax(y)]
    p0 = [A0, mu0, 0.05, 1.0]

    bounds = (
        [0, mu0 - 0.20, 0.005, 0.1],
        [5*A0, mu0 + 0.20, 0.25, 10.0]
    )

    try:
        popt,_ = curve_fit(emg, x, y, p0=p0, bounds=bounds, maxfev=20000)
        y_fit = emg(x, *popt)
        return (x, y_fit, popt)
    except:
        return None

peak_fits = []
for i, pt in enumerate(selected):
    res = fit_emg_window(t, sig, pt)
    if res:
        peak_fits.append((i, pt, *res))


df = pd.DataFrame({"Time [us]": time_hi, "Combined TOF [pC/dt]": sig})
df.to_csv("combinedTOF_evt29_1130.csv")
# ------------------------------------------------------------
# Plot clean, non-overlapping EMG fits
# ------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 200,
    "font.size": 20,
    "axes.linewidth": 1.6,
})

fig, ax = plt.subplots(figsize=(13,7))

ax.plot(t, sig, color="#1f77b4", label="Corrected TOF", lw=2)

colors = plt.cm.tab20(np.linspace(0,1,20))

for (i, center, x_fit, y_fit, params) in peak_fits:
    ax.plot(x_fit, y_fit, "--", lw=2.0, color=colors[i])
    ax.text(center, np.max(y_fit)*1.10,
            f"Peak {i+1}",
            ha="center", va="bottom",
            fontsize=16, weight="bold",
            color=colors[i],
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.65))

ax.set_xlim(1.2, 9.5)
ax.set_ylim(-0.001, np.max(sig)*1.25)

ax.set_xlabel("Time (µs)")
ax.set_ylabel(r"Signal [pC/$\Delta$t]")
# ax.set_title("EMG Fits for 15 Isolated TOF Peaks", pad=15)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ============================================================
# Make identical plot for EVENT 30 (no changes to anything else)
# ============================================================

EVENT2 = "30"
OUT_PLOT2 = "Plots/evt30_combined_single_waveform.png"

with h5py.File(INPUT_FILE, "r") as f:
    grp2 = f[EVENT2]
    raw_m2 = read_ds(grp2, "TOF M")
    raw_l2 = read_ds(grp2, "TOF L")
    raw_h2 = read_ds(grp2, "TOF H")
    time_hi2 = read_ds(grp2, "Time (high sampling)")

# Convert using same rules
tof_m2 = raw_m2
tof_l2 = raw_l2
tof_h2 = (1/4.67) * raw_h2

# Baseline subtract
tof_m2_bs, _ = subtract_baseline(tof_m2)
tof_l2_bs, _ = subtract_baseline(tof_l2)
tof_h2_bs, _ = subtract_baseline(tof_h2)

# Scale low → mid
tof_l2_scaled = tof_l2_bs * LOW_GAIN_FACTOR

# Detect saturation
m2_peak = np.max(tof_m2_bs)
sat2_threshold = SATURATION_FRAC * m2_peak
sat2_mask = tof_m2_bs >= sat2_threshold
sat2_mask = (time_hi2 >= -.1762) & (time_hi2 <= .0621)
sat3_mask = (time_hi2 >= 2.5) & (time_hi2 <= 3.18)


# Combined channel (simple version: identical logic to initial step)
combined2 = np.copy(tof_h2_bs)
combined2[sat2_mask] = tof_m2_bs[sat2_mask]

# ============================================================
# Plot with IDENTICAL STYLE, AXES, COLORS (no mass labels)
# ============================================================

plt.rcParams.update({
    "figure.dpi": 200,
    "font.size": 20,
    "axes.linewidth": 1.6,
})

fig2, ax2 = plt.subplots(figsize=(8,5))


# ax2.plot(time_hi2, tof_l2_bs, color="g", lw=2, label="Corrected TOF (Event 30)")
# ax2.plot(time_hi2, tof_m2_bs, color="b", lw=2, label="Corrected TOF (Event 30)")
# ax2.plot(time_hi2, tof_h2_bs, color="r", lw=2, label="Corrected TOF (Event 30)")


ax2.plot(time_hi2, combined2+1e-3, color="#1f77b4", lw=.75, label="Corrected TOF (Event 30)")
# plt.yscale("log")
ax2.set_xlim(-20, 20)
# ax2.set_ylim(-0.001, np.max(combined2)*1.25)
plt.axhline(-.0005, c="r", linestyle="--")
ax2.set_xlabel("Time (µs)")
ax2.set_ylabel(r"Signal [pC/$\Delta$t]")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_PLOT2)
plt.show()

