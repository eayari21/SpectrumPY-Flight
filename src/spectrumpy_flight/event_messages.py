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
    for fmt in (
        "%Y-%m-%dT%H:%M:%S.%fZ",
        "%Y-%m-%d %H:%M:%S.%f",
        "%d-%b-%Y %H:%M:%S.%f",
    ):
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
# Dust hits (FULL, UNIQUE)
# ============================================================
dust_hit_strings = [
    "2025-11-26T12:13:06.552647Z",
    "2025-11-28T21:23:56.726200Z",
    "2025-11-28T22:28:41.346840Z",
    "2025-11-29T12:08:49.962823Z",
    "2025-12-08T22:56:21.604477Z",
    "2025-12-15T08:08:12.610856Z",
    "2025-12-24T17:33:17.770637Z",
    "2025-12-25T01:19:21.763000Z",
    "2025-12-25T08:46:58.157140Z",
]

dust_hits = sorted({parse_dust_time(s) for s in dust_hit_strings})

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
# Time span
# ============================================================
t_min = min(t for t, _ in events)
t_max = max(t for t, _ in events)

# ============================================================
# OFF background
# ============================================================
ax.plot([t_min, t_max], [Y_ON, Y_ON],
        lw=12, color="tab:red", solid_capstyle="butt")

# ============================================================
# ON intervals
# ============================================================
for i, (start, end) in enumerate(intervals):
    if i == len(intervals) - 1:
        start -= timedelta(days=2)

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
# Transmits
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
# Dust hits — ABSOLUTE NON-OVERLAP (FORCED STACK)
# ============================================================
dust_colors = plt.cm.tab10.colors

DUST_STEP = 0.35
DUST_LABEL_OFFSET = 0.10

for i, t in enumerate(dust_hits, 1):
    c = dust_colors[(i - 1) % len(dust_colors)]

    # force every dust hit onto its own row
    y = Y_DUST + (i - 1) * DUST_STEP

    ax.scatter(
        t, y,
        s=180,
        marker="*",
        color=c,
        zorder=3,
    )

    ax.text(
        t,
        y + DUST_LABEL_OFFSET,
        f"Dust {i}",
        ha="center",
        va="bottom",
        fontsize=16,
        weight="bold",
        color=c,
        zorder=4,
    )

# expand y-limits to fit all dust hits
ax.set_ylim(0.6, Y_DUST + DUST_STEP * (len(dust_hits) + 1))



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
