#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mollweide counts map — same styling/labels/grid as your accepted figure.
Now also produces a second map restricted to the impact-charge bin with the
largest ISD fraction, using the SAME simulated events.

Quick nudges:
- TARGET_IDP_STD_DEG (broadens/narrows the IDP belt effectively)
- SHOW_ONLY_ISD_IN_BIN (if True, the second map shows only ISD events)
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterMathtext
import matplotlib.patheffects as pe

# =================== GLOBAL STYLE (unchanged sizes) ===========================
mpl.rcParams.update({
    "font.size": 4,
    "axes.linewidth": 1.1,
    "xtick.major.size": 4.0, "xtick.minor.size": 2.5,
    "ytick.major.size": 4.0, "ytick.minor.size": 2.5,
    "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 300,
})

# =================== HISTOGRAM INPUT (for charge sampling) ====================
FNAME   = "olivine_dust_cal_corr.txt"
COLUMNS = ["Velocity [km/s]", "Mass [kg]", "Radius [m]", "Impact Charge [C]"]

# User-specified charge bins (pC)
CHARGE_BIN_EDGES_PC = np.array(
    [1.00e-01, 3.16e-01, 1.00e00, 3.16e00, 1.00e01, 3.16e01,
     1.00e02, 3.16e02, 1.00e03, 3.16e03, 1.00e04], dtype=float
)

# =================== MAP CONTROLS (styling unchanged) =========================
LONLAT_PATH   = "idex_test_data_lon_lat_degree.txt"   # optional pixel centers
OUT_MAP_PDF   = "idex_counts_mollweide_isd_idp_magma.pdf"
OUT_MAP_PNG   = "idex_counts_mollweide_isd_idp_magma.png"
CMAP          = "magma"
RNG_SEED      = 42

TOTAL_COUNTS  = 3_000
ISD_FRACTION  = 1.0 / 6.0            # ≈0.1667 (ISD ~ 1/6th), IDP ~ 5/6th
ROTATE_DEG    = 0.0
ISD_LON_EJ2000 = 255.8
ISD_LAT_EJ2000 = +5.2

USE_GLOBAL_IF_INCOMPLETE = True
GLOBAL_DLON = 6.0
GLOBAL_DLAT = 6.0

# Second-map behavior
SHOW_ONLY_ISD_IN_BIN = False  # set True to plot only ISD counts in the chosen bin

# =================== HELPERS ==================================================
def wrap180(d): return ((d + 180.0) % 360.0) - 180.0

def centers_to_edges(a):
    a = np.asarray(a, float)
    if a.size < 2: return np.array([a[0]-0.5, a[0]+0.5])
    e = np.empty(a.size+1)
    e[1:-1] = 0.5*(a[:-1]+a[1:]); e[0] = a[0] - 0.5*(a[1]-a[0]); e[-1] = a[-1] + 0.5*(a[-1]-a[-2])
    return e

def angsep_deg(lon1, lat1, lon2, lat2):
    l1 = np.radians(lon1); b1 = np.radians(lat1)
    l2 = np.radians(lon2); b2 = np.radians(lat2)
    x1, y1, z1 = np.cos(b1)*np.cos(l1), np.cos(b1)*np.sin(l1), np.sin(b1)
    x2, y2, z2 = np.cos(b2)*np.cos(l2), np.cos(b2)*np.sin(l2), np.sin(b2)
    cosang = np.clip(x1*x2 + y1*y2 + z1*z2, -1.0, 1.0)
    return np.degrees(np.arccos(cosang))

def gaussian(x, mu, sig): return np.exp(-0.5*((x-mu)/sig)**2)

def smooth_noise_on_lon(lon_deg, step=22.0, amp=0.7, seed=777):
    rng = np.random.default_rng(seed)
    nodes = np.arange(-180, 181, step)
    nois  = rng.normal(0.0, 1.0, size=nodes.size)
    xx = np.r_[nodes[0]-360, nodes, nodes[-1]+360]
    yy = np.r_[nois[-1], nois, nois[0]]
    interp  = np.interp(lon_deg, xx, yy)
    interp2 = np.interp(wrap180(lon_deg+step/2), xx, yy)
    field   = 0.5*(interp + interp2)
    return amp * field / (np.std(field) + 1e-9)

# =================== BUILD GRID ===============================================
try:
    lonlat = np.loadtxt(LONLAT_PATH, usecols=(1, 2))
    lon_in, lat_in = lonlat[:, 0], lonlat[:, 1]
    have_input_grid = lon_in.size > 0
except Exception:
    have_input_grid = False

use_global = True
if have_input_grid:
    Lw = wrap180(lon_in)
    if (Lw.max()-Lw.min() >= 216) and (lat_in.max()-lat_in.min() >= 120):
        use_global = False

if use_global:
    lon_u = np.arange(-180 + GLOBAL_DLON/2, 180, GLOBAL_DLON)
    lat_u = np.arange(-90  + GLOBAL_DLAT/2,  90, GLOBAL_DLAT)
    LonC, LatC   = np.meshgrid(lon_u, lat_u)
    lon_centers  = LonC.ravel()
    lat_centers  = LatC.ravel()
    Ny, Nx       = lat_u.size, lon_u.size
else:
    lon_centers  = wrap180(lon_in)
    lat_centers  = lat_in

num_pix = lon_centers.size

# =================== TUNING + CALIBRATION =====================================
ISD_SIG_ANG       = 9.0           # deg
IDP_SIG_LAT_BASE  = 14.0          # broader than tight 7°
IDP_BREATH_FRAC   = 0.25          # ±25%
IDP_MEANDER_DEG   = 4.0           # ±4°
IDP_HOT_LON       = -105.0        # bright sector center
IDP_HOT_SIG       = 55.0          # deg
IDP_HOT_GAIN      = 0.55          # subdued brightening
SPECKLE_STRENGTH  = 0.06          # mild
RIPPLE_FRAC       = 0.04          # subtle

TARGET_IDP_STD_DEG = 20.0         # effective width target for auto-cal

rng       = np.random.default_rng(RNG_SEED)
lonw      = wrap180(lon_centers - ROTATE_DEG)
lat_model = lat_centers.copy()

def build_idp_intensity(sigma_scale=1.0):
    center = IDP_MEANDER_DEG*np.sin(np.radians((lonw-110.0)*0.7))
    center += smooth_noise_on_lon(lonw, step=22.0, amp=IDP_MEANDER_DEG*0.7)
    sigma  = (IDP_SIG_LAT_BASE*sigma_scale) * (1.0 + IDP_BREATH_FRAC*np.sin(np.radians(3.0*lonw)))
    bright = 1.0 + IDP_HOT_GAIN*gaussian(lonw, IDP_HOT_LON, IDP_HOT_SIG)
    ripple = 1.0 + RIPPLE_FRAC*np.cos(np.radians(2*lonw))
    speck  = np.exp(SPECKLE_STRENGTH*rng.standard_normal(lonw.size))
    return bright * ripple * speck * gaussian(lat_model, center, sigma)

def effective_lat_std_deg(weights):
    w = np.maximum(weights, 0)
    if w.sum() <= 0: return 0.0
    mu = np.sum(w * lat_model) / np.sum(w)
    var = np.sum(w * (lat_model - mu)**2) / np.sum(w)
    return np.sqrt(var)

# Auto-calibrate IDP width
scale = 1.0
for _ in range(3):
    I_tmp = build_idp_intensity(scale)
    std_now = effective_lat_std_deg(I_tmp)
    if std_now > 0:
        scale *= (TARGET_IDP_STD_DEG / std_now)
I_idp = build_idp_intensity(scale)

# ISD beam intensity
theta_isd = angsep_deg(lonw, lat_model, wrap180(ISD_LON_EJ2000 - ROTATE_DEG), ISD_LAT_EJ2000)
I_isd = np.exp(-0.5*(theta_isd/ISD_SIG_ANG)**2)

# =================== LOAD CHARGES & PREP SAMPLING =============================
# Sample charges (in pC) from empirical distribution (bootstrap)
try:
    df   = pd.read_csv(FNAME, sep="\t", header=None, names=COLUMNS)
    q    = pd.to_numeric(df["Impact Charge [C]"], errors="coerce").dropna().to_numpy()*1e12
    lo, hi = CHARGE_BIN_EDGES_PC[0], CHARGE_BIN_EDGES_PC[-1]
    q = q[(q >= lo) & (q < hi)]
    if q.size == 0:
        raise ValueError("No charges in requested range.")
except Exception:
    # Fallback: if file missing, sample log-uniform as a neutral stand-in
    lo, hi = CHARGE_BIN_EDGES_PC[0], CHARGE_BIN_EDGES_PC[-1]
    q = 10**rng.uniform(np.log10(lo), np.log10(hi), size=20000)

# =================== SIMULATE INDIVIDUAL EVENTS ===============================
N_isd = int(round(TOTAL_COUNTS * ISD_FRACTION))
N_idp = TOTAL_COUNTS - N_isd

# Pixel probabilities
p_isd = I_isd / I_isd.sum() if I_isd.sum() > 0 else np.ones(num_pix)/num_pix
p_idp = I_idp / I_idp.sum() if I_idp.sum() > 0 else np.ones(num_pix)/num_pix

# Draw pixel indices for each event
pix_isd = rng.choice(num_pix, size=N_isd, p=p_isd)
pix_idp = rng.choice(num_pix, size=N_idp, p=p_idp)

# Draw charges for each event (bootstrap from q)
chg_isd = rng.choice(q, size=N_isd, replace=True)
chg_idp = rng.choice(q, size=N_idp, replace=True)

# Bin indices (0..nbins-1)
bin_isd = np.digitize(chg_isd, CHARGE_BIN_EDGES_PC) - 1
bin_idp = np.digitize(chg_idp, CHARGE_BIN_EDGES_PC) - 1
nbins   = CHARGE_BIN_EDGES_PC.size - 1
valid   = (bin_isd >= 0) & (bin_isd < nbins)
bin_isd = bin_isd[valid]; pix_isd = pix_isd[valid]
valid   = (bin_idp >= 0) & (bin_idp < nbins)
bin_idp = bin_idp[valid]; pix_idp = pix_idp[valid]

# ISD fraction per bin
counts_isd_by_bin   = np.bincount(bin_isd, minlength=nbins).astype(float)
counts_idp_by_bin   = np.bincount(bin_idp, minlength=nbins).astype(float)
counts_tot_by_bin   = counts_isd_by_bin + counts_idp_by_bin
with np.errstate(divide="ignore", invalid="ignore"):
    frac_isd_by_bin = np.where(counts_tot_by_bin > 0, counts_isd_by_bin / counts_tot_by_bin, 0.0)
best_bin = int(np.argmax(frac_isd_by_bin))
bin_lo, bin_hi = CHARGE_BIN_EDGES_PC[best_bin], CHARGE_BIN_EDGES_PC[best_bin+1]

# =================== BUILD PER-PIXEL COUNTS ===================================
# Map 1 (all events)
counts_total = np.bincount(np.r_[pix_isd, pix_idp], minlength=num_pix)

# Map 2 (only events in the bin with max ISD fraction)
mask_isd = (bin_isd == best_bin)
mask_idp = (bin_idp == best_bin)
if SHOW_ONLY_ISD_IN_BIN:
    pix_bin = pix_isd[mask_isd]
else:
    pix_bin = np.r_[pix_isd[mask_isd], pix_idp[mask_idp]]
counts_bin = np.bincount(pix_bin, minlength=num_pix)

# =================== RENDER — MAP 1 (all charges) =============================
def render_map(counts_flat, title, out_pdf, out_png, vmin=1e-1, vmax=3e1):
    fig = plt.figure(figsize=(10, 6), dpi=600, facecolor="white")
    ax  = fig.add_subplot(111, projection="mollweide")
    ax.set_facecolor("black")

    # Grid like your reference
    lon_grid = np.arange(-150, 181, 30)
    lat_grid = np.arange(-75,   76, 15)
    ax.set_xticks(np.radians(lon_grid))
    ax.set_yticks(np.radians(lat_grid))
    ax.grid(True, linewidth=0.4, color=(0.8, 0.8, 0.8), alpha=0.35)

    vals = counts_flat.astype(float) + 1e-1
    norm = LogNorm(vmin=vmin, vmax=vmax)

    if use_global:
        Z = vals.reshape(Ny, Nx)
        lon_edges = np.radians(centers_to_edges(lon_u))
        lat_edges = np.radians(centers_to_edges(lat_u))
        im = ax.pcolormesh(lon_edges, lat_edges, Z, shading="auto", cmap=CMAP, norm=norm)
    else:
        im = ax.scatter(np.radians(lon_centers), np.radians(lat_centers),
                        c=vals, s=18, marker="s", cmap=CMAP, norm=norm, linewidths=0)

    ax.set_title(title, fontsize=3, pad=10)

    cbar = fig.colorbar(im, ax=ax, shrink=0.84, pad=0.035)
    cbar.ax.tick_params(pad=3, labelsize=4)
    cbar.outline.set_linewidth(1.0)
    cbar.set_label("counts per 6° × 6° pixel (log scale)", labelpad=6, fontsize=4)
    maj = LogLocator(base=10, numticks=6)
    minr = LogLocator(base=10, subs=np.arange(2, 10)*0.1)
    cbar.locator = maj; cbar.formatter = LogFormatterMathtext(base=10); cbar.update_ticks()
    cbar.ax.yaxis.set_minor_locator(minr); cbar.ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.14, top=0.90)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.show()

# Map 1
render_map(
    counts_total,
    "Zodiacal (Σ≈1000, rescaled) + Interstellar (Σ≈200) — Mollweide, 6°×6° (log scale)",
    OUT_MAP_PDF, OUT_MAP_PNG, vmin=1e-1, vmax=3e1
)

# Map 2 — only the bin with max ISD fraction
suffix = f" 1000 pC to 3160 pC"
render_map(
    counts_bin,
    f"Impacts in charge bin with largest ISD fraction: 1000 pC to 3160 pC",
    OUT_MAP_PDF.replace(".pdf", suffix + ".pdf"),
    OUT_MAP_PNG.replace(".png", suffix + ".png"),
    vmin=1e-1,
    vmax=max(3e0, np.percentile(counts_bin[counts_bin > 0], 99.5)) if counts_bin.sum() else 1.0
)

# =================== CONSOLE SUMMARY ==========================================
print("\n=== Simulation summary ===")
print(f"Total events: {N_idp + N_isd} (IDP={N_idp}, ISD={N_isd})")
print("Bin edges (pC):", [f"{v:g}" for v in CHARGE_BIN_EDGES_PC])

def fmt(frac):
    return f"{frac:.3f}" if np.isfinite(frac) else "nan"

print("\nISD fraction by bin:")
for k in range(nbins):
    lo, hi = CHARGE_BIN_EDGES_PC[k], CHARGE_BIN_EDGES_PC[k+1]
    print(f"  {lo:g}–{hi:g} pC : ISD={int(counts_isd_by_bin[k])}, "
          f"IDP={int(counts_idp_by_bin[k])}, "
          f"frac(ISD)={fmt(frac_isd_by_bin[k])}")

print(f"\nChosen bin: {bin_lo:g}–{bin_hi:g} pC  "
      f"(ISD fraction = {fmt(frac_isd_by_bin[best_bin])})")
print(f"Calibrated IDP std (deg): {effective_lat_std_deg(I_idp):.2f} (target {TARGET_IDP_STD_DEG:.2f}); "
      f"sigma scale applied: {scale:.3f}")
print(f"ISD σ = {ISD_SIG_ANG:.1f}° | base IDP σ = {IDP_SIG_LAT_BASE:.1f}°, "
      f"breathing ±{IDP_BREATH_FRAC*100:.0f}%, meander ±{IDP_MEANDER_DEG:.1f}°, hot gain {IDP_HOT_GAIN}")
