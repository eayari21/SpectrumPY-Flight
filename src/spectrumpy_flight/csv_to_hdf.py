import h5py
import numpy as np
import pandas as pd
import os

# ------------------------------------------------------------
# User paths
# ------------------------------------------------------------
CSV_FILE = "combinedTOF_evt29_1130.csv"
H5_FILE  = "HDF5/Flight/imap_idex_l0_raw_20251130_v001.h5"   # <-- update if needed
EVENT    = "29"   # event number

# ------------------------------------------------------------
# Load CSV (expects 2 columns: time, signal)
# ------------------------------------------------------------
df = pd.read_csv(CSV_FILE)

# support headers or no headers
if df.shape[1] != 2:
    raise ValueError("CSV must have exactly 2 columns: time, signal")

time_arr   = df.iloc[:,0].to_numpy()
signal_arr = df.iloc[:,1].to_numpy()

print(f"Loaded CSV: {len(time_arr)} samples")

# ------------------------------------------------------------
# HDF5 write
# ------------------------------------------------------------
with h5py.File(H5_FILE, "r+") as f:
    base = f[f"{EVENT}/Analysis/DustComposition"]

    # --- CombinedSignal dataset ---
    if "CombinedSignal" in base:
        del base["CombinedSignal"]
    base.create_dataset("CombinedSignal", data=signal_arr)

    # --- CombinedTime dataset ---
    if "CombinedTime" in base:
        del base["CombinedTime"]
    base.create_dataset("CombinedTime", data=time_arr)

print("✓ Finished writing CombinedSignal and CombinedTime to HDF5")
