#!/usr/bin/env python3
"""
Collect scalar values from /<event>/Analysis/Target LImpactCharge across .h5 files
using STRICT h5py extraction: ds[()].

- Handles common scalar encodings:
    * true scalar (0-d) numeric
    * 1-element numeric arrays
    * scalar bytes/str holding a number (e.g., b"1.23")
- Prints a brief debug line for the first few datasets it reads.

Usage:
  python collect_target_LImpactCharge.py
  python collect_target_LImpactCharge.py --root "/path/to/Flight_Packets_Analysis_3" --csv limap_charges.csv
"""

import os
import sys
import argparse
from typing import List, Tuple
import h5py
import csv
import numpy as np

DEFAULT_ROOT  = "/Users/etay8828/Desktop/IDEX_Cal/idex-decom/HDF5/Flight_Packets_Analysis_3"
ANALYSIS_GRP  = "Analysis"
DATASET_NAME  = "Target LImpactCharge"   # exact visual name
DEBUG_SHOW    = 8                        # show first N debug lines

def iter_h5_files(root: str, recursive: bool = True) -> List[str]:
    root = os.path.expanduser(root)
    files = []
    if recursive:
        for dirpath, _, filenames in os.walk(root):
            for fn in filenames:
                if fn.lower().endswith(".h5"):
                    files.append(os.path.join(dirpath, fn))
    else:
        for fn in os.listdir(root):
            if fn.lower().endswith(".h5"):
                files.append(os.path.join(root, fn))
    files.sort()
    return files

def _normalize_spaces(s: str) -> str:
    # collapse all unicode spaces to a single ASCII space
    return " ".join(str(s).split())

def _get_dataset_by_name(grp: h5py.Group, wanted: str):
    """
    Return grp[wanted] if present; otherwise try a forgiving space-normalized match
    (to survive weird non-breaking spaces). Still returns the object as-is; reading is via [()].
    """
    if wanted in grp:
        return grp[wanted]
    wanted_norm = _normalize_spaces(wanted)
    for k in grp.keys():
        if _normalize_spaces(k) == wanted_norm:
            return grp[k]
    return None

def _strict_scalar_from_ds(ds) -> float:
    """
    STRICT read via ds[()], with minimal conversions for common scalar encodings.
    No attribute lookups, no heuristics beyond interpreting a scalar number.
    """
    val = ds[()]  # <-- the only data access
    # bytes/str that actually hold a number
    if isinstance(val, (bytes, bytearray)):
        s = val.decode("utf-8", "ignore").strip()
        return float(s)
    if isinstance(val, str):
        return float(val.strip())
    # numpy scalar or array
    arr = np.asarray(val)
    if arr.shape == ():  # true 0-d scalar
        return float(arr)
    if arr.size == 1:    # 1-element array
        return float(arr.reshape(()))
    # More than one element is unexpected for a scalar dataset
    raise ValueError(f"Dataset is not scalar (shape={arr.shape}, dtype={arr.dtype})")

def extract_LImpactCharge_from_file(path: str, dbg_left: list) -> List[Tuple[str, float]]:
    """
    Return list of (event_name, value) for all /<event>/Analysis/Target LImpactCharge in this file.
    """
    out: List[Tuple[str, float]] = []
    with h5py.File(path, "r") as f:
        # Fast path: event groups at root
        for key, obj in f.items():
            if not isinstance(obj, h5py.Group):
                continue
            if ANALYSIS_GRP in obj:
                ag = obj[ANALYSIS_GRP]
                ds = _get_dataset_by_name(ag, DATASET_NAME)
                if isinstance(ag, h5py.Group) and isinstance(ds, h5py.Dataset):
                    try:
                        val = _strict_scalar_from_ds(ds)
                    except Exception:
                        continue
                    out.append((str(key), val))
                    # Debug print for the first few only
                    if dbg_left[0] > 0:
                        _dbg_print(path, key, ds, val)
                        dbg_left[0] -= 1

        # Fallback traversal if nothing found above (structure variance)
        if not out:
            def visitor(name, obj):
                if isinstance(obj, h5py.Dataset) and name.endswith(f"{ANALYSIS_GRP}/{DATASET_NAME}"):
                    parts = name.split("/")
                    event = parts[-3] if len(parts) >= 3 else "(unknown)"
                    try:
                        val = _strict_scalar_from_ds(obj)
                    except Exception:
                        return
                    out.append((event, val))
                    if dbg_left[0] > 0:
                        _dbg_print(path, event, obj, val, full_name=name)
                        dbg_left[0] -= 1
            f.visititems(lambda n, o: visitor(n, o))
    return out

def _dbg_print(file_path, event, ds, val, full_name=None):
    try:
        dtype = ds.dtype
        shape = ds.shape
        raw = ds[()]
        # best-effort attribute peek without assumptions
        attrs = {k: ds.attrs.get(k) for k in ("_FillValue", "scale_factor", "add_offset")}
    except Exception:
        dtype = shape = raw = attrs = "<?>"
    name = full_name if full_name else f"/{event}/{ANALYSIS_GRP}/{ds.name.split('/')[-1]}"
    print(f"[DEBUG] {os.path.basename(file_path)} :: {name} :: dtype={dtype}, shape={shape}, "
          f"raw={repr(raw)} -> parsed={val}, isfinite={np.isfinite(val)}")

def write_csv(rows: List[Tuple[str, str, float]], csv_path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "event", "Target LImpactCharge"])
        for r in rows:
            w.writerow(r)

def main():
    p = argparse.ArgumentParser(description="Collect Analysis/Target LImpactCharge from .h5 files (strict ds[()] read).")
    p.add_argument("--root", default=DEFAULT_ROOT, help="Root folder containing .h5 files")
    p.add_argument("--no-recursive", action="store_true", help="Do not descend into subdirectories")
    p.add_argument("--csv", default="", help="Optional path to write a CSV of [file,event,value]")
    args = p.parse_args()

    files = iter_h5_files(args.root, recursive=not args.no_recursive)
    if not files:
        print(f"No .h5 files found under: {args.root}", file=sys.stderr)
        sys.exit(2)

    charges: List[float] = []
    rows: List[Tuple[str, str, float]] = []
    total_events = 0
    missing_in_files = 0
    dbg_left = [DEBUG_SHOW]

    for fp in files:
        found = extract_LImpactCharge_from_file(fp, dbg_left)
        if not found:
            missing_in_files += 1
        for event, val in found:
            charges.append(val)
            rows.append((os.path.basename(fp), event, val))
        total_events += len(found)

    # Summary
    print(f"Scanned {len(files)} file(s) under: {args.root}")
    print(f"Found {total_events} event(s) with '{ANALYSIS_GRP}/{DATASET_NAME}'.")
    if missing_in_files:
        print(f"{missing_in_files} file(s) contained no matching dataset.")

    # Preview
    if charges:
        preview_n = min(10, len(charges))
        print(f"First {preview_n} values:", charges[:preview_n])
        finite = np.isfinite(charges)
        if finite.any():
            finite_vals = np.array(charges, float)[finite]
            print(f"Finite stats → min={finite_vals.min()}, max={finite_vals.max()}, count={finite_vals.size}/{len(charges)} finite")
        else:
            print("All parsed values are NaN/inf.")

    # Optional CSV
    if args.csv:
        write_csv(rows, args.csv)
        print(f"Wrote CSV: {args.csv}")

    return charges, rows

if __name__ == "__main__":
    main()
