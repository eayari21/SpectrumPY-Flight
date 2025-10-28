from __future__ import annotations

import argparse
import json
import math
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Mapping, MutableMapping, Sequence

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


SNR_CHANNELS: Sequence[str] = (
    "TOF L",
    "TOF M",
    "TOF H",
    "Ion Grid",
    "Target H",
    "Target L",
)

TARGET_CHANNELS: Sequence[str] = ("Ion Grid", "Target H", "Target L")

TRIGGER_PAIRS: Mapping[str, tuple[str, str]] = {
    "Ion Grid vs TOF H": ("Ion Grid", "TOF H"),
    "TOF H vs Target": ("TOF H", "Target"),
    "Target vs Ion Grid": ("Target", "Ion Grid"),
}


def _finite_values(values: Iterable[float] | np.ndarray | float) -> np.ndarray:
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        return arr
    mask = np.isfinite(arr)
    return arr[mask]


def _hist_figure(
    data: Iterable[float],
    *,
    title: str,
    xlabel: str,
    caption: str,
    requirement_lines: Sequence[tuple[float | None, str]] | None = None,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 6))
    finite = _finite_values(data)

    if finite.size:
        bins = min(40, max(10, int(round(math.sqrt(finite.size)))))
        ax.hist(
            finite,
            bins=bins,
            color="#2a6ea6",
            edgecolor="white",
        )
        ax.set_ylabel("Count")
        if requirement_lines:
            for value, label in requirement_lines:
                if value is None:
                    continue
                ax.axvline(value, color="#c94b32", linestyle="--", linewidth=1.5)
                ax.text(
                    value,
                    ax.get_ylim()[1] * 0.9,
                    label,
                    rotation=90,
                    color="#c94b32",
                    va="top",
                    ha="left",
                    fontsize=9,
                )
    else:
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=14,
        )
        ax.set_axis_off()

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    fig.text(0.5, 0.02, caption, ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


def _summary_figure(summary_text: Sequence[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")

    y = 0.94
    ax.text(0.05, y, "Olivine Test Metrics Summary", fontsize=18, fontweight="bold")
    y -= 0.06

    for paragraph in summary_text:
        wrapped = textwrap.fill(paragraph, 90)
        ax.text(0.05, y, wrapped, fontsize=11, va="top")
        y -= 0.05 * (wrapped.count("\n") + 1)

    return fig


def _dataset_scalar(dataset: h5py.Dataset | None) -> float | None:
    if dataset is None:
        return None
    data = np.asarray(dataset[()]).astype(float).ravel()
    if data.size == 0:
        return None
    value = float(data.flat[0])
    if not np.isfinite(value):
        return None
    return value


def _interpolate_first_crossing(time: np.ndarray, signal: np.ndarray, level: float) -> float | None:
    indices = np.where(signal >= level)[0]
    if indices.size == 0:
        return None
    idx = int(indices[0])
    if idx == 0:
        return float(time[0])
    prev_idx = idx - 1
    y0 = float(signal[prev_idx])
    y1 = float(signal[idx])
    t0 = float(time[prev_idx])
    t1 = float(time[idx])
    if not (np.isfinite(y0) and np.isfinite(y1) and np.isfinite(t0) and np.isfinite(t1)):
        return None
    if y1 == y0:
        return float(t1)
    frac = (level - y0) / (y1 - y0)
    return float(t0 + frac * (t1 - t0))


def _infer_trigger_time(event_group: h5py.Group, analysis_group: h5py.Group | None, channel: str) -> float | None:
    if analysis_group is not None:
        metric_dataset = analysis_group.get(f"{channel} 10pct Time")
        value = _dataset_scalar(metric_dataset)
        if value is not None:
            return value

        params_dataset = analysis_group.get(f"{channel}FitParams")
        if params_dataset is not None:
            params = np.asarray(params_dataset[()]).astype(float).ravel()
            if params.size and np.isfinite(params[0]):
                return float(params[0])

    waveform_dataset = event_group.get(channel)
    if waveform_dataset is None:
        return None

    data = np.asarray(waveform_dataset[()], dtype=float).ravel()
    if data.size == 0:
        return None

    if channel.startswith("TOF"):
        time_dataset_name = "Time (high sampling)"
    else:
        time_dataset_name = "Time (low sampling)"

    time_dataset = event_group.get(time_dataset_name)
    if time_dataset is None:
        return None

    time = np.asarray(time_dataset[()], dtype=float).ravel()
    n = min(time.size, data.size)
    if n == 0:
        return None
    time = time[:n]
    data = data[:n]

    baseline_region = max(5, data.size // 20)
    baseline = float(np.nanmedian(data[:baseline_region]))
    amplitude = float(np.nanmax(data) - baseline)
    if not np.isfinite(amplitude) or amplitude <= 0:
        return None

    level = baseline + 0.10 * amplitude
    return _interpolate_first_crossing(time, data, level)


@dataclass
class MetricCollector:
    snr: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    rise: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    decay: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    chi_sq: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    reduced_chi_sq: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    trigger_deltas: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    saturation_counts: Counter[str] = field(default_factory=Counter)
    target_usage: Counter[str] = field(default_factory=Counter)
    files_processed: int = 0
    events_processed: int = 0

    def consume_file(self, path: Path) -> None:
        with h5py.File(path, "r") as handle:
            self.files_processed += 1
            for event_id in handle.keys():
                event_group = handle[event_id]
                self.events_processed += 1
                analysis_group = event_group.get("Analysis")
                metadata_group = event_group.get("Metadata")

                if analysis_group is not None:
                    for dataset_name, dataset in analysis_group.items():
                        if dataset_name.endswith("SNR"):
                            channel = dataset_name[:-3]
                            values = _finite_values(dataset[()])
                            if values.size:
                                self.snr[channel].extend(values.tolist())

                    for channel in TARGET_CHANNELS:
                        params_dataset = analysis_group.get(f"{channel}FitParams")
                        if params_dataset is not None:
                            params = np.asarray(params_dataset[()], dtype=float).ravel()
                            if params.size >= 5:
                                rise_value = float(params[3])
                                decay_value = float(params[4])
                                if np.isfinite(rise_value):
                                    self.rise[channel].append(rise_value)
                                if np.isfinite(decay_value):
                                    self.decay[channel].append(decay_value)

                        chi_dataset = analysis_group.get(f"{channel}ChiSquared")
                        if chi_dataset is not None:
                            chi_values = _finite_values(chi_dataset[()])
                            if chi_values.size:
                                self.chi_sq[channel].extend(chi_values.tolist())

                        red_dataset = analysis_group.get(
                            f"{channel}ReducedChiSquared"
                        )
                        if red_dataset is not None:
                            red_values = _finite_values(red_dataset[()])
                            if red_values.size:
                                self.reduced_chi_sq[channel].extend(
                                    red_values.tolist()
                                )

                if metadata_group is not None:
                    for dataset_name, dataset in metadata_group.items():
                        if dataset_name.endswith(" Saturated"):
                            channel = dataset_name.replace(" Saturated", "")
                            values = np.asarray(dataset[()]).astype(int).ravel()
                            count = int(np.count_nonzero(values))
                            if count:
                                self.saturation_counts[channel] += count

                trigger_times: dict[str, float] = {}
                if analysis_group is not None:
                    for channel in TARGET_CHANNELS + ("Target",):
                        trigger_times[channel] = float("nan")

                for channel in {"TOF H", "Ion Grid", "Target H", "Target L"}:
                    value = _infer_trigger_time(event_group, analysis_group, channel)
                    if value is not None and np.isfinite(value):
                        trigger_times[channel] = value

                target_channel = None
                if "Target H" in trigger_times and np.isfinite(trigger_times["Target H"]):
                    target_channel = "Target H"
                elif "Target L" in trigger_times and np.isfinite(trigger_times["Target L"]):
                    target_channel = "Target L"

                if target_channel is not None:
                    self.target_usage[target_channel] += 1
                    trigger_times["Target"] = trigger_times[target_channel]

                for pair_name, (channel_a, channel_b) in TRIGGER_PAIRS.items():
                    value_a = trigger_times.get(channel_a)
                    value_b = trigger_times.get(channel_b)
                    if value_a is None or value_b is None:
                        continue
                    if not (np.isfinite(value_a) and np.isfinite(value_b)):
                        continue
                    delta = abs(float(value_a) - float(value_b))
                    self.trigger_deltas[pair_name].append(delta)

    def _summary_lines(self) -> list[str]:
        lines = [
            f"Processed {self.events_processed} events across {self.files_processed} file(s).",
        ]

        if self.saturation_counts:
            saturation_bits = [
                f"{channel}: {count}" for channel, count in sorted(self.saturation_counts.items())
            ]
            lines.append(
                "Saturation detections — " + ", ".join(saturation_bits)
            )
        else:
            lines.append("No saturation flags were detected in the analysed files.")

        if self.target_usage:
            usage_bits = [
                f"{channel} ({count} event{'s' if count != 1 else ''})"
                for channel, count in sorted(self.target_usage.items())
            ]
            lines.append(
                "Trigger deltas paired with: " + ", ".join(usage_bits)
            )
        else:
            lines.append("No target trigger times were available for delta calculations.")

        lines.append(
            "All histograms include the instrument requirements from the signal handbook when thresholds are defined."
        )
        return lines

    def write_report(self, pdf_path: Path) -> None:
        summary_fig = _summary_figure(self._summary_lines())

        with PdfPages(pdf_path) as pdf:
            pdf.savefig(summary_fig)
            plt.close(summary_fig)

            for channel in SNR_CHANNELS:
                requirement = []
                if channel.startswith("TOF"):
                    requirement.append((3.0, "SNR ≥ 3"))
                else:
                    requirement.append((6.0, "SNR ≥ 6"))
                fig = _hist_figure(
                    self.snr.get(channel, []),
                    title=f"{channel} Signal-to-Noise Ratio",
                    xlabel="SNR",
                    caption=f"Histogram of {channel} SNR across all analysed olivine events.",
                    requirement_lines=requirement,
                )
                pdf.savefig(fig)
                plt.close(fig)

            for channel in TARGET_CHANNELS:
                fig = _hist_figure(
                    self.rise.get(channel, []),
                    title=f"{channel} Rise Time (τ_rise)",
                    xlabel="Rise time (µs)",
                    caption=f"Distribution of fitted rise constants for the {channel} waveform fits.",
                    requirement_lines=((0.1, "Spec lower"), (100.0, "Spec upper")),
                )
                pdf.savefig(fig)
                plt.close(fig)

            for channel in TARGET_CHANNELS:
                fig = _hist_figure(
                    self.decay.get(channel, []),
                    title=f"{channel} Decay Time (τ_decay)",
                    xlabel="Decay time (µs)",
                    caption=f"Distribution of fitted decay constants for the {channel} waveform fits.",
                    requirement_lines=((1e-5, "Spec lower"), (10.0, "Spec upper")),
                )
                pdf.savefig(fig)
                plt.close(fig)

            for pair_name in ("Ion Grid vs TOF H", "TOF H vs Target", "Target vs Ion Grid"):
                fig = _hist_figure(
                    self.trigger_deltas.get(pair_name, []),
                    title=f"Trigger Delta: {pair_name}",
                    xlabel="|Δt| (µs)",
                    caption="Absolute timing offset between trigger references for the channel pair.",
                    requirement_lines=((20.0, "Spec upper"),),
                )
                pdf.savefig(fig)
                plt.close(fig)

            for channel in TARGET_CHANNELS:
                fig = _hist_figure(
                    self.chi_sq.get(channel, []),
                    title=f"{channel} χ²",
                    xlabel="χ²",
                    caption=f"Goodness-of-fit χ² statistic for the {channel} model fits.",
                    requirement_lines=(),
                )
                pdf.savefig(fig)
                plt.close(fig)

            for channel in TARGET_CHANNELS:
                fig = _hist_figure(
                    self.reduced_chi_sq.get(channel, []),
                    title=f"{channel} Reduced χ²",
                    xlabel="Reduced χ²",
                    caption=f"Reduced χ² values for the {channel} waveform fits.",
                    requirement_lines=((1.0, "Ideal"),),
                )
                pdf.savefig(fig)
                plt.close(fig)

    def write_summary_json(self, path: Path) -> None:
        def stats(values: Iterable[float]) -> Mapping[str, float | int]:
            arr = _finite_values(values)
            if arr.size == 0:
                return {"count": 0}
            return {
                "count": int(arr.size),
                "mean": float(np.mean(arr)),
                "median": float(np.median(arr)),
                "std": float(np.std(arr, ddof=0)),
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
            }

        summary = {
            "files_processed": self.files_processed,
            "events_processed": self.events_processed,
            "snr": {channel: stats(values) for channel, values in self.snr.items()},
            "rise": {channel: stats(values) for channel, values in self.rise.items()},
            "decay": {channel: stats(values) for channel, values in self.decay.items()},
            "chi_squared": {channel: stats(values) for channel, values in self.chi_sq.items()},
            "reduced_chi_squared": {
                channel: stats(values) for channel, values in self.reduced_chi_sq.items()
            },
            "trigger_deltas": {
                pair: stats(values) for pair, values in self.trigger_deltas.items()
            },
            "saturation_counts": dict(self.saturation_counts),
            "target_usage": dict(self.target_usage),
        }

        path.write_text(json.dumps(summary, indent=2))


def _iter_input_paths(inputs: Iterable[Path | str]) -> Iterator[Path]:
    for raw_path in inputs:
        path = Path(raw_path)
        if not path.exists():
            continue
        if path.is_dir():
            yield from sorted(path.rglob("*.h5"))
        else:
            yield path


def generate_olivine_metrics(
    inputs: Iterable[Path | str],
    output_dir: Path | str,
    report_name: str = "olivine_metrics_report.pdf",
) -> Mapping[str, Path]:
    collector = MetricCollector()
    consumed_paths = []
    for path in _iter_input_paths(inputs):
        collector.consume_file(path)
        consumed_paths.append(path)

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    pdf_path = output / report_name
    summary_path = output / "olivine_metrics_summary.json"

    collector.write_report(pdf_path)
    collector.write_summary_json(summary_path)

    return {
        "pdf": pdf_path,
        "summary": summary_path,
        "inputs": consumed_paths,
        "events": collector.events_processed,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate olivine calibration histograms and summary report."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more HDF5 files or directories containing decoded olivine events.",
    )
    parser.add_argument(
        "--output-dir",
        default="olivine_metrics",
        help="Destination directory for the generated report and summary JSON.",
    )
    parser.add_argument(
        "--report-name",
        default="olivine_metrics_report.pdf",
        help="Filename for the generated PDF report.",
    )
    args = parser.parse_args(argv)

    results = generate_olivine_metrics(args.inputs, args.output_dir, args.report_name)
    pdf_path: Path = results["pdf"]
    summary_path: Path = results["summary"]
    events: int = results["events"]
    print(f"Processed {events} event(s).")
    print(f"Report saved to {pdf_path}")
    print(f"Summary saved to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
