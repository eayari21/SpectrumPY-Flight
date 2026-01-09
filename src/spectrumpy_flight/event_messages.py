import csv
import re
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

# ============================================================
# Files
# ============================================================
input_csv = "CSV/oasis_events-3.csv"

# ============================================================
# Match strings
# ============================================================
ON_STRING = "SCI state change: sciState16Dictionary(ACQSETUP) ==&#62; sciState16Dictionary(ACQ)"
OFF_STRING = "SCI state change: sciState16Dictionary(ACQ) ==&#62; sciState16Dictionary(CHILL)"

timestamp_regex = re.compile(r"\d{4}/\d{3}-\d{2}:\d{2}:\d{2}\.\d{3}")

# ============================================================
# Helpers
# ============================================================
def doy_to_datetime(ts):
    year = int(ts[:4])
    doy = int(ts[5:8])
    time_part = ts[9:]
    base = datetime(year, 1, 1) + timedelta(days=doy - 1)
    t = datetime.strptime(time_part, "%H:%M:%S.%f")
    return base.replace(
        hour=t.hour,
        minute=t.minute,
        second=t.second,
        microsecond=t.microsecond,
    )

def parse_dust_time(s):
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%d-%b-%Y %H:%M:%S.%f"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    raise ValueError(s)

def parse_tx_date(fname):
    date_str = fname.split("_")[4]
    return datetime.strptime(date_str, "%Y%m%d")

# ============================================================
# Parse state changes
# ============================================================
events = []

with open(input_csv, newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    for row in reader:
        for cell in row:
            ts_match = timestamp_regex.search(cell)
            if not ts_match:
                continue
            ts = ts_match.group()
            if ON_STRING in cell:
                events.append((doy_to_datetime(ts), "on"))
            elif OFF_STRING in cell:
                events.append((doy_to_datetime(ts), "off"))

events.sort(key=lambda x: x[0])

# ============================================================
# Pair ON → OFF intervals
# ============================================================
intervals = []
on_time = None

for t, state in events:
    if state == "on":
        on_time = t
    elif state == "off" and on_time:
        intervals.append((on_time, t))
        on_time = None

# ============================================================
# Dust hits
# ============================================================
dust_hit_strings = [
    "2025-11-28 21:23:56.726200",
    "2025-11-28 22:28:41.346840",
    "25-Dec-2025 01:19:21.763000",
    "2025-12-25 08:46:58.15714",
    "2025-12-28 23:17:27.27143",
]
dust_hits = [parse_dust_time(s) for s in dust_hit_strings]

# ============================================================
# Transmit dates from filenames
# ============================================================
tx_files = [
    "imap_idex_l0_raw_20251130_v002.pkts",
    "imap_idex_l0_raw_20251209_v002.pkts",
    "imap_idex_l0_raw_20251217_v001.pkts",
    "imap_idex_l0_raw_20251218_v003.pkts",
    "imap_idex_l0_raw_20251225_v001.pkts",
    "imap_idex_l0_raw_20260102_v001.pkts",
]
tx_dates = [parse_tx_date(f) for f in tx_files]

# ============================================================
# Plot configuration
# ============================================================
plt.rcParams.update({
    "font.size": 18,
    "axes.titlesize": 26,
    "axes.labelsize": 22,
    "xtick.labelsize": 18,
})

fig, ax = plt.subplots(figsize=(20, 6))

# ============================================================
# Lane definitions
# ============================================================
Y_ON = 1
Y_TX = 2
Y_DUST = 3

# ============================================================
# Plot ON intervals
# ============================================================
for start, end in intervals:
    days = (end - start).total_seconds() / 86400
    ax.plot([start, end], [Y_ON, Y_ON], lw=12, solid_capstyle="butt")
    ax.text(
        start + (end - start) / 2,
        Y_ON + 0.15,
        f"{days:.1f} days",
        ha="center",
        va="bottom",
        fontsize=18,
        weight="bold",
    )

# ============================================================
# Plot transmit dates
# ============================================================
for i, t in enumerate(tx_dates, 1):
    ax.scatter(t, Y_TX, s=120, marker="s")
    ax.text(
        t,
        Y_TX + 0.15,
        f"TX {i}",
        ha="center",
        va="bottom",
        fontsize=16,
    )

# ============================================================
# Plot dust hits
# ============================================================
for i, t in enumerate(dust_hits, 1):
    ax.scatter(t, Y_DUST, s=160, marker="*")
    ax.text(
        t,
        Y_DUST + 0.15,
        f"Dust {i}",
        ha="center",
        va="bottom",
        fontsize=16,
        weight="bold",
    )

# ============================================================
# Axes formatting
# ============================================================
ax.set_yticks([Y_ON, Y_TX, Y_DUST])
ax.set_yticklabels([
    "Science Acquisition",
    "Downlink / Transmit",
    "Dust Hits",
])

ax.set_ylim(0.6, 3.6)
ax.set_xlabel("Date (UTC)")
ax.set_title(
    "IDEX Science Acquisition, Dust Hits, and Transmit Timeline\n(Since Launch – September 24)",
    pad=15,
)

ax.grid(True, axis="x", alpha=0.3)
fig.autofmt_xdate()
plt.tight_layout()
plt.show()
