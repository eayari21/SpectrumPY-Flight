import csv
import re
from datetime import datetime, timedelta
from collections import defaultdict

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
# Transmit dates
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
# Determine full time span
# ============================================================
t_min = min(t for t, _ in events)
t_max = max(t for t, _ in events)

# ============================================================
# Plot OFF background (red)
# ============================================================
ax.plot([t_min, t_max], [Y_ON, Y_ON],
        lw=12, color="tab:red", solid_capstyle="butt")

# ============================================================
# Plot ON intervals (blue)
# ============================================================
for i, (start, end) in enumerate(intervals):
    if i == len(intervals) - 1:
        start = start - timedelta(days=2)

    days = (end - start).total_seconds() / 86400

    ax.plot([start, end], [Y_ON, Y_ON],
            lw=12, color="tab:blue", solid_capstyle="butt")

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
# Plot transmit dates (closer labels, stronger staggering)
# ============================================================
tx_colors = plt.cm.tab10.colors
TX_LABEL_STEP = 0.22
TX_LABEL_OFFSET = 0.10

for i, t in enumerate(tx_dates, 1):
    c = tx_colors[(i - 1) % len(tx_colors)]
    ax.scatter(t, Y_TX, s=120, marker="s", color=c)
    ax.text(
        t,
        Y_TX + TX_LABEL_OFFSET + TX_LABEL_STEP * ((i - 1) % 3),
        f"TX {i}",
        ha="center",
        va="bottom",
        fontsize=16,
        color=c,
    )

# ============================================================
# Plot dust hits (closer labels, stronger vertical staggering)
# ============================================================
dust_colors = plt.cm.tab10.colors
dust_offsets = defaultdict(int)
time_bin = timedelta(hours=12)

DUST_STEP = 0.22
DUST_LABEL_OFFSET = 0.08

for i, t in enumerate(sorted(dust_hits), 1):
    key = None
    for prev_t in dust_offsets:
        if abs(t - prev_t) < time_bin:
            key = prev_t
            break

    if key is None:
        dust_offsets[t] = 0
        offset = 0
    else:
        dust_offsets[key] += 1
        offset = dust_offsets[key]

    c = dust_colors[(i - 1) % len(dust_colors)]
    y = Y_DUST + DUST_STEP * offset

    ax.scatter(t, y, s=180, marker="*", color=c)
    ax.text(
        t,
        y + DUST_LABEL_OFFSET,
        f"Dust {i}",
        ha="center",
        va="bottom",
        fontsize=16,
        weight="bold",
        color=c,
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

ax.set_ylim(0.6, 4.1)
ax.set_xlabel("Date (UTC)")
ax.set_title(
    "IDEX Science Acquisition, Dust Hits, and Transmit Timeline\n(Since Launch – September 24)",
    pad=15,
)

ax.grid(True, axis="x", alpha=0.3)
fig.autofmt_xdate()
plt.tight_layout()
plt.show()

# ============================================================
# Daily duty cycle stats
# ============================================================
on_seconds_by_day = defaultdict(float)

for start, end in intervals:
    current = start
    while current.date() <= end.date():
        day_start = datetime.combine(current.date(), datetime.min.time())
        day_end = day_start + timedelta(days=1)

        seg_start = max(start, day_start)
        seg_end = min(end, day_end)

        if seg_start < seg_end:
            on_seconds_by_day[current.date()] += (
                seg_end - seg_start
            ).total_seconds()

        current = day_start + timedelta(days=1)

# ============================================================
# Print results
# ============================================================
for day in sorted(on_seconds_by_day):
    pct = 100 * on_seconds_by_day[day] / 86400
    doy = day.timetuple().tm_yday
    print(f"{day.month:02d}/{day.day:02d} (DOY {doy}), {pct:.1f}%")
