from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
import textwrap
from typing import Dict, Iterable, List, Tuple

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pytest

matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from spectrumpy_flight import idex_packet


CSV_TIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"
TIME_TOLERANCE_SECONDS = 6.0
REPORTS_DIR = Path("reports")


def _read_csv_counts(path: Path) -> Dict[int, float]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        next(reader)  # header
        mapping: Dict[int, float] = {}
        for row in reader:
            if len(row) < 2:
                continue
            time_str, value_str = row[0].strip(), row[1].strip()
            try:
                timestamp = datetime.strptime(time_str, CSV_TIME_FORMAT).replace(
                    tzinfo=timezone.utc
                )
                value = int(float(value_str))
            except ValueError:
                continue
            mapping[value] = timestamp.timestamp()
    return mapping


def _extract_event_records(pkts_path: Path) -> List[Tuple[float, int, int]]:
    packets = idex_packet.IDEXEvent(str(pkts_path))
    records: List[Tuple[float, int, int]] = []
    event_ids: Iterable[int] = sorted({key[0] for key in packets.header})
    for event_id in event_ids:
        timestamp = packets.header.get((event_id, "Timestamp"))
        scievt_count = packets.header.get((event_id, "IDX__TXHDREVTNUM"))
        sent_unproc = packets.header.get((event_id, "IDX__TXHDRTRANSCNT"))
        if timestamp is None or scievt_count is None or sent_unproc is None:
            continue
        records.append((float(timestamp), int(scievt_count), int(sent_unproc)))
    return records


def _render_pdf_report(
    pdf_path: Path,
    *,
    matched_events: int,
    max_delta: float,
    tolerance: float,
    plot_fig: matplotlib.figure.Figure,
) -> None:
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    background = (
        "This report documents alignment between timestamps provided by IDEX LASP "
        "packet data and spacecraft counters exported to CSV. Maintaining "
        "alignment ensures science events are synchronized with spacecraft "
        "telemetry for downstream analysis."
    )
    methods = (
        "Event records were read from raw .pkts files using idex_packet. Each event "
        "timestamp was matched against reference science and transmission counters "
        "from IDX_SW.SCIEVT_CNT.csv and IDX_SW.SCISENTUNPROCCNT.csv. Differences "
        "between packet timestamps and reference counters were compared against the "
        f"{tolerance:.1f}s tolerance threshold."
    )
    results = (
        f"A total of {matched_events} events were compared. The maximum absolute "
        f"timestamp difference observed was {max_delta:.2f} seconds, which is "
        f"within the {tolerance:.1f}-second tolerance threshold."
    )

    with PdfPages(pdf_path) as pdf:
        text_fig = plt.figure(figsize=(8.5, 11))
        text_fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        sections = [
            ("Background", background),
            ("Methods", methods),
            ("Results", results),
        ]

        y = 0.9
        for heading, body in sections:
            text_fig.text(0.1, y, heading, fontsize=16, fontweight="bold")
            y -= 0.05
            for line in textwrap.wrap(body, 85):
                text_fig.text(0.1, y, line, fontsize=12, wrap=True)
                y -= 0.03
            y -= 0.04

        pdf.savefig(text_fig, bbox_inches="tight")
        plt.close(text_fig)

        pdf.savefig(plot_fig, bbox_inches="tight")


@pytest.mark.skipif(
    not idex_packet._HAS_LASP_PACKETS,  # type: ignore[attr-defined]
    reason="idex_packet timestamp alignment requires lasp_packets dependency",
)
def test_idex_packet_timestamps_match_science_counts(tmp_path: Path) -> None:
    pkts_files = [
        Path("src/spectrumpy_flight/Data/Flight/imap_idex_l0_raw_20251130_v002.pkts"),
        Path("src/spectrumpy_flight/Data/Flight/imap_idex_l0_raw_20251209_v002.pkts"),
    ]
    csv_counts = {
        "scievt": _read_csv_counts(
            Path("tests/timestamp_test/IDX_SW.SCIEVT_CNT.csv")
        ),
        "sent_unproc": _read_csv_counts(
            Path("tests/timestamp_test/IDX_SW.SCISENTUNPROCCNT.csv")
        ),
    }

    all_records: List[Tuple[float, int, int]] = []
    for pkts_path in pkts_files:
        all_records.extend(_extract_event_records(pkts_path))

    assert all_records, "Expected idex_packet to provide timestamped science events."

    deltas: List[float] = []
    scievt_times: List[Tuple[float, int]] = []
    sent_unproc_times: List[Tuple[float, int]] = []

    for timestamp, scievt_count, sent_unproc in all_records:
        if scievt_count in csv_counts["scievt"]:
            delta = abs(timestamp - csv_counts["scievt"][scievt_count])
            deltas.append(delta)
            scievt_times.append((timestamp, scievt_count))
        if sent_unproc in csv_counts["sent_unproc"]:
            delta = abs(timestamp - csv_counts["sent_unproc"][sent_unproc])
            deltas.append(delta)
            sent_unproc_times.append((timestamp, sent_unproc))

    assert deltas, "No matching counts between idex_packet output and CSV references."
    assert np.all(np.array(deltas) <= TIME_TOLERANCE_SECONDS)

    matched_events = len(deltas)
    max_delta = max(deltas)

    # Report
    report_path = tmp_path / "idex_timestamp_alignment.txt"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("IDX timestamp alignment summary\n")
        handle.write(f"Total events compared: {matched_events}\n")
        handle.write(
            f"Max |timestamp - reference|: {max_delta:.2f} seconds (tolerance {TIME_TOLERANCE_SECONDS:.1f}s)\n"
        )

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    scievt_csv_times = sorted(csv_counts["scievt"].items())
    sent_csv_times = sorted(csv_counts["sent_unproc"].items())
    if scievt_csv_times:
        ax.plot(
            [datetime.fromtimestamp(t, tz=timezone.utc) for t in dict(scievt_csv_times).values()],
            [count for count, _ in scievt_csv_times],
            label="SCIEVT_CNT (CSV)",
        )
    if sent_csv_times:
        ax.plot(
            [datetime.fromtimestamp(t, tz=timezone.utc) for t in dict(sent_csv_times).values()],
            [count for count, _ in sent_csv_times],
            label="SCISENTUNPROCCNT (CSV)",
        )
    if scievt_times:
        scievt_times.sort()
        ax.scatter(
            [datetime.fromtimestamp(t, tz=timezone.utc) for t, _ in scievt_times],
            [count for _, count in scievt_times],
            color="black",
            marker="x",
            label="idex_packet (SCIEVT)",
        )
    if sent_unproc_times:
        sent_unproc_times.sort()
        ax.scatter(
            [datetime.fromtimestamp(t, tz=timezone.utc) for t, _ in sent_unproc_times],
            [count for _, count in sent_unproc_times],
            color="orange",
            marker="+",
            label="idex_packet (SCISENTUNPROCCNT)",
        )

    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("Count")
    ax.set_title("IDEX timestamp alignment with spacecraft counters")
    ax.legend()
    fig.autofmt_xdate()
    plot_path = tmp_path / "idex_timestamp_alignment.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")

    pdf_path = REPORTS_DIR / "idex_timestamp_alignment_report.pdf"
    _render_pdf_report(
        pdf_path,
        matched_events=matched_events,
        max_delta=max_delta,
        tolerance=TIME_TOLERANCE_SECONDS,
        plot_fig=fig,
    )
    plt.close(fig)

    assert report_path.exists()
    assert plot_path.exists()
    assert pdf_path.exists()
