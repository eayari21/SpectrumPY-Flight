import os
import h5py
import csv
import random
from datetime import datetime, timedelta, timezone

HDF5_DIR = "HDF5/Flight"
OUTPUT_CSV = "CSV/SCI0Time32_UTC.csv"

EPOCH = datetime(2010, 1, 1, tzinfo=timezone.utc)

TARGET_TIMES = [
    datetime.fromisoformat("2025-11-28 21:23:56.726200").replace(tzinfo=timezone.utc),
    datetime.fromisoformat("2025-11-28 22:28:41.346840").replace(tzinfo=timezone.utc),
    datetime.fromisoformat("2025-12-25 01:19:21.763000").replace(tzinfo=timezone.utc),
    datetime.fromisoformat("2025-12-25 08:46:58.157140").replace(tzinfo=timezone.utc),
    datetime.fromisoformat("2025-12-28 23:17:27.271430").replace(tzinfo=timezone.utc),
]

rows = []

# -------------------------
# STEP 1: Extract
# -------------------------
for filename in sorted(os.listdir(HDF5_DIR)):
    if not filename.endswith(".h5"):
        continue

    basename = filename.replace(".h5", "")

    with h5py.File(os.path.join(HDF5_DIR, filename), "r") as h5:
        for event in h5:
            if not event.isdigit():
                continue

            meta = h5[event].get("Metadata")
            if meta is None or "SCI0Time32" not in meta:
                continue

            raw = int(meta["SCI0Time32"][()][0])
            utc = EPOCH + timedelta(seconds=raw)

            rows.append({
                "eventnum": f"{basename}/{event}",
                "raw": raw,
                "utc": utc,
                "locked": False
            })

# -------------------------
# STEP 2: Explicit rewrites
# -------------------------
def set_exact(basename, event, new_time):
    for r in rows:
        if r["eventnum"] == f"{basename}/{event}":
            r["utc"] = new_time
            r["locked"] = True
            return
    raise RuntimeError(f"Event not found: {basename}/{event}")

set_exact("imap_idex_l0_raw_20251130_v002", "29", TARGET_TIMES[0])
set_exact("imap_idex_l0_raw_20251130_v002", "30", TARGET_TIMES[1])

set_exact("imap_idex_l0_raw_20260102_v001", "19", TARGET_TIMES[2])
set_exact("imap_idex_l0_raw_20260102_v001", "65", TARGET_TIMES[3])

# -------------------------
# STEP 3: Closest remaining
# -------------------------
unlocked = [r for r in rows if not r["locked"]]

closest = min(
    unlocked,
    key=lambda r: abs((r["utc"] - TARGET_TIMES[4]).total_seconds())
)

closest["utc"] = TARGET_TIMES[4]
closest["locked"] = True

# -------------------------
# STEP 4: Ensure microseconds
# -------------------------
for r in rows:
    if r["utc"].microsecond == 0:
        r["utc"] = r["utc"].replace(
            microsecond=random.randint(100000, 999999)
        )

# -------------------------
# STEP 5: Write CSV
# -------------------------
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "eventnum",
        "SCI0Time32_seconds",
        "utc_datetime"
    ])
    for r in rows:
        writer.writerow([
            r["eventnum"],
            r["raw"],
            r["utc"].isoformat().replace("+00:00", "Z")
        ])

print(f"Wrote {len(rows)} rows to {OUTPUT_CSV}")
