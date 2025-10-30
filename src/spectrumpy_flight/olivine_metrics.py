from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, MutableMapping, Sequence

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages

from .line_shapes import emg as _emg_profile
from .time2mass import time2mass as optimise_time2mass
from .lookup import dust_estimator as _dust_estimator

try:  # pragma: no cover - SciPy is optional in some environments
    from scipy.optimize import curve_fit as _curve_fit  # type: ignore

    _HAVE_CURVE_FIT = True
except Exception:  # pragma: no cover - gracefully degrade without SciPy
    _HAVE_CURVE_FIT = False
    _curve_fit = None


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

EXPECTED_MASS_LINES: Sequence[tuple[str, float]] = (
    ("1H", 1.0),
    ("12C", 12.0),
    ("16O", 16.0),
    ("23Na", 23.0),
    ("24Mg", 24.0),
    ("25Mg", 25.0),
    ("26Mg", 26.0),
    ("28Si", 28.0),
    ("29Si", 29.0),
    ("30Si", 30.0),
    ("38K", 38.0),
    ("54Fe", 54.0),
    ("56Fe", 56.0),
    ("57Fe", 57.0),
    ("58Fe", 58.0),
)

MG_ISOTOPES = {"24Mg", "25Mg", "26Mg"}
SI_ISOTOPES = {"28Si", "29Si", "30Si"}
FE_ISOTOPES = {"54Fe", "56Fe", "57Fe", "58Fe"}

MASS_FIT_WINDOW = 0.6


@dataclass
class ParticleEstimateRecord:
    event_id: str
    channel: str
    charge_c: float
    velocity_kms: float
    mass_kg: float
    yield_c_per_kg: float
    rise_time_us: float | None
    ion_to_target_ratio: float | None
    velocity_source: str


@dataclass
class FigureSpec:
    stem: str
    caption: str
    label: str
    builder: Callable[[], plt.Figure]


@dataclass
class FigureAsset:
    path: Path
    caption: str
    label: str


def _resolve_default_coefficients():
    try:
        rise_table, ratio_table, yield_table = _dust_estimator.load_default_tables()
    except Exception:
        return None, None, None

    def _choose(table, preferred="combined"):
        if table is None:
            return None
        try:
            return table.get(preferred)
        except KeyError:
            labels = list(table.labels())
            if not labels:
                raise
            return table.get(labels[0])

    try:
        rise_params = _choose(rise_table)
        ratio_params = _choose(ratio_table)
        yield_params = _choose(yield_table)
    except Exception:
        return None, None, None
    return rise_params, ratio_params, yield_params


_RISE_PARAMS, _RATIO_PARAMS, _YIELD_PARAMS = _resolve_default_coefficients()


def _finite_values(values: Iterable[float] | np.ndarray | float) -> np.ndarray:
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        return arr
    mask = np.isfinite(arr)
    return arr[mask]


def _approximate_mode(values: np.ndarray) -> float | None:
    if values.size == 0:
        return None
    unique = np.unique(values)
    if unique.size == 1:
        return float(unique[0])
    bins = int(max(10, min(50, math.sqrt(values.size))))
    counts, edges = np.histogram(values, bins=bins)
    if not np.any(counts):
        return float(values[0])
    idx = int(np.argmax(counts))
    upper_index = min(idx + 1, edges.size - 1)
    return float((edges[idx] + edges[upper_index]) / 2.0)


def _descriptive_stats(values: Iterable[float]) -> dict[str, float | int | None]:
    arr = _finite_values(values)
    if arr.size == 0:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "mode": None,
            "std": None,
        }

    stats = {
        "count": int(arr.size),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr, ddof=0)),
    }
    mode = _approximate_mode(arr)
    stats["mode"] = None if mode is None else float(mode)
    return stats


def _format_value(value: float | int | None) -> str:
    if value is None:
        return "--"
    try:
        numeric = float(value)
    except Exception:
        return "--"
    if not math.isfinite(numeric):
        return "--"
    abs_value = abs(numeric)
    if abs_value and (abs_value >= 1_000 or abs_value < 1e-2):
        return f"{numeric:.3e}"
    return f"{numeric:.3f}"


def _maybe_log_axis(ax: plt.Axes, axis: str, values: Iterable[float]) -> None:
    arr = _finite_values(values)
    if arr.size == 0:
        return
    positive = arr[arr > 0]
    if positive.size == 0:
        return
    span = float(np.max(positive) / np.min(positive))
    if span >= 100.0:
        if axis == "x":
            ax.set_xscale("log")
        else:
            ax.set_yscale("log")


def _latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    escaped = str(text)
    for original, replacement in replacements.items():
        escaped = escaped.replace(original, replacement)
    return escaped


def _render_stats_table(
    title: str,
    label: str,
    stats_map: Mapping[str, Mapping[str, float | int | None]],
) -> str:
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{_latex_escape(title)}}}",
        rf"\label{{{label}}}",
        r"\begin{tabular}{lrrrrrrr}",
        r"\toprule",
        "Metric & Count & Min & Max & Mean & Median & Mode & Std. Dev. \\\\",
        r"\midrule",
    ]
    if not stats_map:
        lines.append(r"\multicolumn{8}{c}{No data available.} \\")
    else:
        for name in sorted(stats_map):
            stats = stats_map[name]
            row = " & ".join(
                [
                    _latex_escape(name),
                    str(stats.get("count", 0)),
                    _format_value(stats.get("min")),
                    _format_value(stats.get("max")),
                    _format_value(stats.get("mean")),
                    _format_value(stats.get("median")),
                    _format_value(stats.get("mode")),
                    _format_value(stats.get("std")),
                ]
            )
            lines.append(f"{row} \\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    return "\n".join(lines)


def _render_count_table(title: str, label: str, counts: Mapping[str, int]) -> str:
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{_latex_escape(title)}}}",
        rf"\label{{{label}}}",
        r"\begin{tabular}{lr}",
        r"\toprule",
        "Channel & Occurrences \\\\",
        r"\midrule",
    ]
    if counts:
        for name, count in sorted(counts.items()):
            lines.append(f"{_latex_escape(name)} & {count} \\")
    else:
        lines.append(r"\multicolumn{2}{c}{No data available.} \\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    return "\n".join(lines)


def _render_longtable(
    caption: str,
    label: str,
    column_spec: str,
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
) -> str:
    header_line = " & ".join(headers) + " \\\\"
    lines = [
        rf"\begin{{longtable}}{{{column_spec}}}",
        rf"\caption{{{_latex_escape(caption)}}}\label{{{label}}}\\",
        r"\toprule",
        header_line.rstrip(),
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        header_line.rstrip(),
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
    ]
    if rows:
        for row in rows:
            lines.append(" & ".join(row) + r" \\")
    else:
        lines.append(
            rf"\multicolumn{{{len(headers)}}}{{c}}{{No data available.}} \\")
    lines.append(r"\end{longtable}")
    return "\n".join(lines)


def _compile_latex(tex_path: Path) -> Path | None:
    engines = [engine for engine in ("tectonic", "pdflatex", "xelatex", "lualatex") if shutil.which(engine)]
    for engine in engines:
        if engine == "tectonic":
            command = [engine, "--outdir", str(tex_path.parent), tex_path.name]
        else:
            command = [engine, "-interaction=nonstopmode", "-halt-on-error", tex_path.name]
        result = subprocess.run(
            command,
            cwd=tex_path.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode == 0:
            candidate = tex_path.with_suffix(".pdf")
            if candidate.exists():
                return candidate
    return None


def _mass_vs_velocity_figure(records: Sequence[ParticleEstimateRecord]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("white")
    if not records:
        ax.axis("off")
        ax.text(0.5, 0.5, "No velocity estimates available.", ha="center", va="center")
        return fig

    speeds = np.array([rec.velocity_kms for rec in records], dtype=float)
    masses = np.array([rec.mass_kg for rec in records], dtype=float)
    channels = [rec.channel for rec in records]
    color_map = {
        "Target H": "#1f77b4",
        "Target L": "#d62728",
        "Ion Grid": "#2ca02c",
    }
    colors = [color_map.get(channel, "#555555") for channel in channels]

    ax.scatter(
        speeds,
        masses,
        c=colors,
        edgecolor="white",
        linewidth=0.7,
        alpha=0.9,
    )
    handles = []
    for channel, color in color_map.items():
        if channel in channels:
            handles.append(
                ax.scatter([], [], c=color, label=channel, edgecolor="white", linewidth=0.7)
            )
    if handles:
        ax.legend(handles=handles, frameon=False, loc="lower right")
    ax.set_xlabel("Velocity (km/s)", fontsize=12)
    ax.set_ylabel("Mass (kg)", fontsize=12)
    ax.set_title("Mass versus velocity estimates", fontsize=14)
    _maybe_log_axis(ax, "x", speeds)
    _maybe_log_axis(ax, "y", masses)
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)
    fig.tight_layout()
    return fig


def _yield_vs_velocity_figure(records: Sequence[ParticleEstimateRecord]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("white")
    if not records:
        ax.axis("off")
        ax.text(0.5, 0.5, "No velocity estimates available.", ha="center", va="center")
        return fig

    speeds = np.array([rec.velocity_kms for rec in records], dtype=float)
    yields = np.array([rec.yield_c_per_kg for rec in records], dtype=float)
    channels = [rec.channel for rec in records]
    color_map = {
        "Target H": "#1f77b4",
        "Target L": "#d62728",
        "Ion Grid": "#2ca02c",
    }

    for channel, color in color_map.items():
        mask = [idx for idx, rec in enumerate(records) if rec.channel == channel]
        if not mask:
            continue
        ax.scatter(
            speeds[mask],
            yields[mask],
            c=color,
            edgecolor="white",
            linewidth=0.7,
            alpha=0.9,
            label=channel,
        )

    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, loc="lower right")
    ax.set_xlabel("Velocity (km/s)", fontsize=12)
    ax.set_ylabel("Charge yield (C/kg)", fontsize=12)
    ax.set_title("Charge yield as a function of velocity", fontsize=14)
    _maybe_log_axis(ax, "x", speeds)
    _maybe_log_axis(ax, "y", yields)
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)
    fig.tight_layout()
    return fig


def _abundance_stack_figure(
    mass_results: Sequence[MassAnalysisResult],
    estimates: Mapping[str, ParticleEstimateRecord],
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.patch.set_facecolor("white")
    results_by_event = {result.event_id: result for result in mass_results}
    paired: list[tuple[float, MassAnalysisResult]] = []
    for event_id, estimate in estimates.items():
        result = results_by_event.get(event_id)
        if result is None or result.total_area <= 0.0:
            continue
        paired.append((estimate.velocity_kms, result))

    if not paired:
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No paired mass spectra and velocity estimates.",
            ha="center",
            va="center",
        )
        return fig

    species_order = [name for name, _ in EXPECTED_MASS_LINES]
    speeds = np.array([item[0] for item in paired], dtype=float)
    bins = max(4, min(12, int(math.sqrt(len(paired)))))
    if float(np.max(speeds)) == float(np.min(speeds)):
        bin_edges = np.linspace(float(np.min(speeds)) * 0.9, float(np.max(speeds)) * 1.1 + 1e-6, bins + 1)
    else:
        bin_edges = np.linspace(float(np.min(speeds)), float(np.max(speeds)), bins + 1)
    bin_indices = np.digitize(speeds, bin_edges, right=False) - 1
    bin_indices = np.clip(bin_indices, 0, bins - 1)

    stack_values = {species: [0.0] * bins for species in species_order}
    counts = [0] * bins
    for idx, (_, result) in zip(bin_indices, paired):
        counts[idx] += 1
        for species in species_order:
            stack_values[species][idx] += float(result.relative_abundances.get(species, 0.0))

    for idx, count in enumerate(counts):
        if count > 0:
            for species in species_order:
                stack_values[species][idx] /= count

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    y_values = [stack_values[species] for species in species_order]
    ax.stackplot(bin_centers, y_values, labels=species_order, alpha=0.9)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Velocity (km/s)", fontsize=12)
    ax.set_ylabel("Relative abundance", fontsize=12)
    ax.set_title("Stacked olivine line abundances vs. velocity", fontsize=14)
    _maybe_log_axis(ax, "x", bin_centers)
    ax.legend(loc="upper right", ncol=2, fontsize=8, frameon=False)
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)
    fig.tight_layout()
    return fig


def _mass_line_probability_figure(
    mass_results: Sequence[MassAnalysisResult],
    estimates: Mapping[str, ParticleEstimateRecord],
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.patch.set_facecolor("white")
    results_by_event = {result.event_id: result for result in mass_results}
    paired: list[tuple[float, MassAnalysisResult]] = []
    for event_id, estimate in estimates.items():
        result = results_by_event.get(event_id)
        if result is None or result.total_area <= 0.0:
            continue
        paired.append((estimate.velocity_kms, result))

    if not paired:
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No paired mass spectra and velocity estimates.",
            ha="center",
            va="center",
        )
        return fig

    species_order = [name for name, _ in EXPECTED_MASS_LINES]
    speeds = np.array([item[0] for item in paired], dtype=float)
    bins = max(4, min(12, int(math.sqrt(len(paired)))))
    if float(np.max(speeds)) == float(np.min(speeds)):
        bin_edges = np.linspace(float(np.min(speeds)) * 0.9, float(np.max(speeds)) * 1.1 + 1e-6, bins + 1)
    else:
        bin_edges = np.linspace(float(np.min(speeds)), float(np.max(speeds)), bins + 1)
    bin_indices = np.digitize(speeds, bin_edges, right=False) - 1
    bin_indices = np.clip(bin_indices, 0, bins - 1)

    probabilities = {species: [0.0] * bins for species in species_order}
    counts = [0] * bins
    for idx, (_, result) in zip(bin_indices, paired):
        counts[idx] += 1
        for species in species_order:
            abundance = float(result.relative_abundances.get(species, 0.0))
            if abundance > 0.0:
                probabilities[species][idx] += 1.0

    for idx, count in enumerate(counts):
        if count > 0:
            for species in species_order:
                probabilities[species][idx] /= count

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    for species in species_order:
        ax.plot(
            bin_centers,
            probabilities[species],
            label=species,
            linewidth=1.6,
        )

    ax.set_ylim(0, 1)
    ax.set_xlabel("Velocity (km/s)", fontsize=12)
    ax.set_ylabel("Appearance probability", fontsize=12)
    ax.set_title("Mass-line detection probability vs. velocity", fontsize=14)
    _maybe_log_axis(ax, "x", bin_centers)
    ax.legend(loc="upper right", ncol=2, fontsize=8, frameon=False)
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)
    fig.tight_layout()
    return fig


def _saturation_bar_figure(saturation_counts: Mapping[str, int]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("white")
    if not saturation_counts:
        ax.axis("off")
        ax.text(0.5, 0.5, "No saturation flags were recorded.", ha="center", va="center")
        return fig

    items = sorted(saturation_counts.items())
    channels = [item[0] for item in items]
    counts = np.array([item[1] for item in items], dtype=float)
    positions = np.arange(len(channels))
    bars = ax.bar(positions, counts, color="#2a6ea6", edgecolor="white")
    ax.set_xticks(positions, [channel.replace(" ", "\n") for channel in channels])
    ax.set_ylabel("Occurrences", fontsize=12)
    ax.set_title("Channel saturation events", fontsize=14)
    for rect, count in zip(bars, counts):
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            rect.get_height() + max(counts) * 0.02,
            f"{int(count)}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    fig.tight_layout()
    return fig


def _hist_figure(
    data: Iterable[float],
    *,
    title: str,
    xlabel: str,
    caption: str,
    requirement_lines: Sequence[tuple[float | None, str]] | None = None,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("white")
    finite = _finite_values(data)

    if finite.size:
        positive = finite[finite > 0]
        nonpositive_count = int(finite.size - positive.size)

        if positive.size:
            min_positive = float(np.min(positive))
            max_positive = float(np.max(positive))
            if np.isclose(min_positive, max_positive):
                min_positive *= 0.9
                max_positive *= 1.1

            log_min = float(np.log10(min_positive))
            log_max = float(np.log10(max_positive))
            if np.isclose(log_min, log_max):
                log_min -= 0.05
                log_max += 0.05

            bin_edges = np.logspace(log_min, log_max, num=41)
            ax.hist(
                positive,
                bins=bin_edges,
                color="#2a6ea6",
                edgecolor="white",
                linewidth=0.7,
            )
            ax.set_xscale("log")
            ax.set_ylabel("Count", fontsize=12)
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.set_major_locator(LogLocator(base=10))
            ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
            ax.xaxis.set_major_formatter(LogFormatterSciNotation())
            ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.5)
            ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.3)
            ax.set_axisbelow(True)
            ax.tick_params(axis="both", labelsize=10)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

            if requirement_lines:
                y_max = ax.get_ylim()[1]
                for value, label in requirement_lines:
                    if value is None or value <= 0:
                        continue
                    ax.axvline(value, color="#c94b32", linestyle="--", linewidth=1.5)
                    ax.annotate(
                        label,
                        xy=(value, y_max),
                        xytext=(5, -10),
                        textcoords="offset points",
                        rotation=90,
                        color="#c94b32",
                        va="top",
                        ha="left",
                        fontsize=9,
                    )

            if nonpositive_count:
                ax.text(
                    0.98,
                    0.95,
                    f"{nonpositive_count} values ≤ 0 omitted",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=9,
                    color="#6c757d",
                )
        else:
            ax.text(
                0.5,
                0.5,
                "Data must be positive for logarithmic histogramming.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
            )
            ax.set_axis_off()
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

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=12)
    fig.text(0.5, 0.02, caption, ha="center", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


def _mass_abundance_figure(results: Sequence["MassAnalysisResult"]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(10, 6))
    species_order = [name for name, _ in EXPECTED_MASS_LINES]

    if results and any(r.total_area > 0 for r in results):
        medians = []
        for species in species_order:
            values = [r.relative_abundances.get(species, 0.0) for r in results]
            arr = _finite_values(values)
            medians.append(float(np.median(arr)) if arr.size else 0.0)

        positions = np.arange(len(species_order))
        ax.bar(positions, medians, color="#2a6ea6", edgecolor="white")
        ax.set_xticks(positions, species_order, rotation=45, ha="right")
        ax.set_ylabel("Median relative abundance")
        ax.set_ylim(0.0, 1.0)
        caption = "Median per-event relative abundance for the expected olivine mass lines."
    else:
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No mass spectral fits were available to compute relative abundances.",
            ha="center",
            va="center",
            fontsize=12,
        )
        caption = (
            "Mass spectrum analysis requires valid TOF H waveforms and SciPy's curve_fit."
        )

    fig.suptitle("Olivine Mass Line Relative Abundances")
    fig.text(0.5, 0.02, caption, ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


def _ternary_to_cartesian(mg: float, si: float, fe: float) -> tuple[float, float]:
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0) / 2.0]])
    bary = np.array([mg, si, fe], dtype=float)
    total = float(bary.sum())
    if total <= 0.0:
        bary = np.array([1.0, 0.0, 0.0])
        total = 1.0
    bary = bary / total
    coord = bary @ vertices
    return float(coord[0]), float(coord[1])


def _ternary_figure(points: Sequence[tuple[float, float, float]]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_axis_off()

    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0) / 2.0]])
    loop = np.vstack([vertices, vertices[0]])
    ax.plot(loop[:, 0], loop[:, 1], color="#444444", linewidth=1.5)
    ax.text(-0.05, -0.05, "Mg", fontsize=12, fontweight="bold")
    ax.text(1.05, -0.05, "Si", fontsize=12, fontweight="bold", ha="right")
    ax.text(0.5, np.sqrt(3.0) / 2.0 + 0.05, "Fe", fontsize=12, fontweight="bold", ha="center")

    if points:
        coords = np.array([_ternary_to_cartesian(*p) for p in points])
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c="#2a6ea6",
            edgecolor="white",
            linewidth=0.75,
            s=60,
        )
    else:
        ax.text(
            0.5,
            0.4,
            "No events provided Mg–Si–Fe abundances for the ternary diagram.",
            ha="center",
            va="center",
            fontsize=12,
        )

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, np.sqrt(3.0) / 2.0 + 0.1)
    fig.suptitle("Mg–Si–Fe Ternary Composition")
    fig.text(
        0.5,
        0.02,
        "Points are normalised per event using the fitted Mg, Si, and Fe line areas.",
        ha="center",
        fontsize=10,
    )
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
class MassLineFit:
    species: str
    target_mass: float
    fit_mass: float | None
    sigma: float | None
    tau: float | None
    area: float
    success: bool

    def to_dict(self) -> dict[str, float | bool | None]:
        return {
            "species": self.species,
            "target_mass": float(self.target_mass),
            "fit_mass": None if self.fit_mass is None else float(self.fit_mass),
            "sigma": None if self.sigma is None else float(self.sigma),
            "tau": None if self.tau is None else float(self.tau),
            "area": float(self.area),
            "success": bool(self.success),
        }


@dataclass
class MassAnalysisResult:
    event_id: str
    baseline: float
    stretch_us: float
    shift_us: float
    fits: list[MassLineFit]
    relative_abundances: dict[str, float]
    raw_areas: dict[str, float]

    def __post_init__(self) -> None:
        self.baseline = float(self.baseline)
        self.stretch_us = float(self.stretch_us)
        self.shift_us = float(self.shift_us)

    @property
    def total_area(self) -> float:
        return float(sum(self.raw_areas.values()))

    def ternary_point(self) -> tuple[float, float, float] | None:
        mg = sum(self.raw_areas.get(species, 0.0) for species in MG_ISOTOPES)
        si = sum(self.raw_areas.get(species, 0.0) for species in SI_ISOTOPES)
        fe = sum(self.raw_areas.get(species, 0.0) for species in FE_ISOTOPES)
        if mg <= 0.0 or si <= 0.0 or fe <= 0.0:
            return None
        total = mg + si + fe
        if total <= 0.0:
            return None
        return float(mg / total), float(si / total), float(fe / total)

    def to_dict(self) -> dict[str, object]:
        return {
            "event_id": self.event_id,
            "baseline": float(self.baseline),
            "stretch_us": float(self.stretch_us),
            "shift_us": float(self.shift_us),
            "relative_abundances": {k: float(v) for k, v in self.relative_abundances.items()},
            "raw_areas": {k: float(v) for k, v in self.raw_areas.items()},
            "fits": [fit.to_dict() for fit in self.fits],
        }


def _fit_mass_line(
    mass_axis: np.ndarray,
    signal: np.ndarray,
    species: str,
    target_mass: float,
) -> MassLineFit:
    mask = (mass_axis >= target_mass - MASS_FIT_WINDOW) & (
        mass_axis <= target_mass + MASS_FIT_WINDOW
    )
    if not np.any(mask):
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    mass_window = np.asarray(mass_axis[mask], dtype=float)
    intensity_window = np.asarray(signal[mask], dtype=float)
    if mass_window.size < 5 or not np.any(np.isfinite(intensity_window)):
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    finite_mask = np.isfinite(intensity_window)
    if not np.any(finite_mask):
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    intensity_window = intensity_window[finite_mask]
    mass_window = mass_window[finite_mask]

    if mass_window.size < 5:
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    positive_mask = intensity_window > 0.0
    area_guess = float(
        np.trapz(np.maximum(intensity_window, 0.0), mass_window)
    )
    if area_guess <= 0.0 and not np.any(positive_mask):
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    peak_idx = int(np.argmax(intensity_window))
    mu_guess = float(mass_window[peak_idx])
    sigma_guess = 0.25
    tau_guess = 0.2

    area_upper = max(
        area_guess * 10.0,
        float(np.abs(intensity_window).max() * (mass_window[-1] - mass_window[0]) * 5.0),
        1.0,
    )

    if not _HAVE_CURVE_FIT:
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    try:
        popt, _ = _curve_fit(
            lambda m, mu, sigma, tau, area: _emg_profile(
                m, mu, sigma, tau, area=area
            ),
            mass_window,
            intensity_window,
            p0=[mu_guess, sigma_guess, tau_guess, max(area_guess, 1.0e-6)],
            bounds=(
                [target_mass - MASS_FIT_WINDOW, 0.05, 0.01, 0.0],
                [target_mass + MASS_FIT_WINDOW, 1.5, 10.0, area_upper],
            ),
            maxfev=10000,
        )
    except Exception:
        return MassLineFit(species, target_mass, None, None, None, 0.0, False)

    fit_mass, sigma, tau, area = map(float, popt)
    area = max(area, 0.0)
    return MassLineFit(species, target_mass, fit_mass, sigma, tau, area, True)


def _analyse_mass_spectrum(event_id: str, event_group: h5py.Group) -> MassAnalysisResult | None:
    waveform = event_group.get("TOF H")
    time_dataset = event_group.get("Time (high sampling)")
    if waveform is None or time_dataset is None:
        return None

    tof = np.asarray(waveform[()], dtype=float).ravel()
    time = np.asarray(time_dataset[()], dtype=float).ravel()
    n = min(tof.size, time.size)
    if n < 10:
        return None
    tof = tof[:n]
    time = time[:n]

    try:
        stretch_us, shift_us, mass_axis = optimise_time2mass(tof, time)
    except Exception:
        return None

    baseline_region = max(10, n // 20)
    baseline = float(np.nanmedian(tof[:baseline_region])) if baseline_region else float(np.nanmedian(tof))
    signal = np.asarray(tof - baseline, dtype=float)

    fits: list[MassLineFit] = []
    raw_areas: dict[str, float] = {}
    for species, target_mass in EXPECTED_MASS_LINES:
        fit = _fit_mass_line(mass_axis, signal, species, target_mass)
        fits.append(fit)
        raw_areas[species] = float(fit.area if fit.success else 0.0)

    total_area = float(sum(raw_areas.values()))
    if total_area > 0.0:
        relative = {species: area / total_area for species, area in raw_areas.items()}
    else:
        relative = {species: 0.0 for species in raw_areas.keys()}

    return MassAnalysisResult(
        event_id=event_id,
        baseline=baseline,
        stretch_us=stretch_us,
        shift_us=shift_us,
        fits=fits,
        relative_abundances=relative,
        raw_areas=raw_areas,
    )


def _relative_stats(
    results: Sequence[MassAnalysisResult],
) -> dict[str, dict[str, float | int | None]]:
    stats: dict[str, dict[str, float | int | None]] = {}
    for species, _ in EXPECTED_MASS_LINES:
        values = [r.relative_abundances.get(species, 0.0) for r in results if r.total_area > 0.0]
        arr = _finite_values(values)
        descriptor = _descriptive_stats(arr)
        stats[species] = {
            "count": descriptor["count"],
            "min": None if descriptor["min"] is None else float(descriptor["min"]),
            "max": None if descriptor["max"] is None else float(descriptor["max"]),
            "mean": None if descriptor["mean"] is None else float(descriptor["mean"]),
            "median": None if descriptor["median"] is None else float(descriptor["median"]),
            "mode": None if descriptor["mode"] is None else float(descriptor["mode"]),
            "std": None if descriptor["std"] is None else float(descriptor["std"]),
        }
    return stats


def _aggregate_elemental_fractions(
    results: Sequence[MassAnalysisResult],
) -> tuple[float, float, float] | None:
    mg_values = [
        sum(r.relative_abundances.get(species, 0.0) for species in MG_ISOTOPES)
        for r in results
        if r.total_area > 0.0
    ]
    si_values = [
        sum(r.relative_abundances.get(species, 0.0) for species in SI_ISOTOPES)
        for r in results
        if r.total_area > 0.0
    ]
    fe_values = [
        sum(r.relative_abundances.get(species, 0.0) for species in FE_ISOTOPES)
        for r in results
        if r.total_area > 0.0
    ]

    mg_arr = _finite_values(mg_values)
    si_arr = _finite_values(si_values)
    fe_arr = _finite_values(fe_values)
    if not (mg_arr.size and si_arr.size and fe_arr.size):
        return None

    return (
        float(np.median(mg_arr) * 100.0),
        float(np.median(si_arr) * 100.0),
        float(np.median(fe_arr) * 100.0),
    )


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
    mass_results: list[MassAnalysisResult] = field(default_factory=list)
    _ternary_points: list[tuple[float, float, float]] = field(default_factory=list)
    impact_charge: MutableMapping[str, list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    particle_estimates: list[ParticleEstimateRecord] = field(default_factory=list)
    event_estimates: dict[str, ParticleEstimateRecord] = field(default_factory=dict)

    def consume_file(self, path: Path) -> None:
        with h5py.File(path, "r") as handle:
            self.files_processed += 1
            for event_id in handle.keys():
                event_group = handle[event_id]
                self.events_processed += 1
                analysis_group = event_group.get("Analysis")
                metadata_group = event_group.get("Metadata")

                event_rise_times: dict[str, float | None] = {
                    channel: None for channel in TARGET_CHANNELS
                }
                event_impact_charges: dict[str, float] = {}

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
                                    event_rise_times[channel] = rise_value
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

                    for channel in TARGET_CHANNELS + ("Ion Grid",):
                        charge_dataset = analysis_group.get(f"{channel}ImpactCharge")
                        if charge_dataset is None:
                            continue
                        charges = _finite_values(charge_dataset[()])
                        if charges.size:
                            self.impact_charge[channel].extend(charges.tolist())
                            event_impact_charges[channel] = float(np.mean(charges))

                if metadata_group is not None:
                    for dataset_name, dataset in metadata_group.items():
                        if dataset_name.endswith(" Saturated"):
                            channel = dataset_name.replace(" Saturated", "")
                            values = np.asarray(dataset[()]).astype(int).ravel()
                            count = int(np.count_nonzero(values))
                            if count:
                                self.saturation_counts[channel] += count

                mass_result = _analyse_mass_spectrum(event_id, event_group)
                if mass_result is not None:
                    self.mass_results.append(mass_result)
                    ternary = mass_result.ternary_point()
                    if ternary is not None:
                        self._ternary_points.append(ternary)

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

                if (
                    target_channel is not None
                    and _RISE_PARAMS is not None
                    and _YIELD_PARAMS is not None
                ):
                    target_charge = event_impact_charges.get(target_channel)
                    if target_charge is not None and target_charge > 0.0:
                        rise_time = event_rise_times.get(target_channel)
                        ion_charge = event_impact_charges.get("Ion Grid")
                        ratio = None
                        if (
                            ion_charge is not None
                            and ion_charge > 0.0
                            and target_charge > 0.0
                        ):
                            try:
                                ratio = ion_charge / target_charge
                            except ZeroDivisionError:
                                ratio = None
                        try:
                            estimate = _dust_estimator.estimate_particle(
                                charge_c=target_charge,
                                rise_time=rise_time,
                                ion_to_target_ratio=ratio,
                                rise_params=_RISE_PARAMS,
                                ratio_params=_RATIO_PARAMS,
                                yield_params=_YIELD_PARAMS,
                                velocity_range=(1.0, 100.0),
                            )
                        except Exception:
                            estimate = None
                        if estimate is not None:
                            record = ParticleEstimateRecord(
                                event_id=event_id,
                                channel=target_channel,
                                charge_c=float(estimate.charge_c),
                                velocity_kms=float(estimate.velocity_kms),
                                mass_kg=float(estimate.mass_kg),
                                yield_c_per_kg=float(estimate.yield_c_per_kg),
                                rise_time_us=None if rise_time is None else float(rise_time),
                                ion_to_target_ratio=None if ratio is None else float(ratio),
                                velocity_source=estimate.velocity_details.source,
                            )
                            self.particle_estimates.append(record)
                            self.event_estimates[event_id] = record

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

        if self.mass_results:
            elemental = _aggregate_elemental_fractions(self.mass_results)
            if elemental is not None:
                mg_percent, si_percent, fe_percent = elemental
                lines.append(
                    "Median elemental fractions — "
                    + ", ".join(
                        [
                            f"Mg: {mg_percent:.1f}%",
                            f"Si: {si_percent:.1f}%",
                            f"Fe: {fe_percent:.1f}%",
                        ]
                    )
                )

            species_stats = _relative_stats(self.mass_results)
            if species_stats:
                formatted = ", ".join(
                    f"{species}: {(values.get('median') or 0.0) * 100.0:.1f}%"
                    for species, values in species_stats.items()
                )
                lines.append("Median per-line abundances — " + formatted)
        else:
            lines.append(
                "No mass spectral waveforms were available for the olivine abundance analysis."
            )
        return lines

    def write_report(self, pdf_path: Path) -> None:
        pdf_path = Path(pdf_path)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        summary_lines = self._summary_lines()

        figure_specs: list[FigureSpec] = []

        figure_specs.append(
            FigureSpec(
                stem="summary",
                caption="Executive summary of the analysed olivine events.",
                label="fig:summary",
                builder=lambda summary_lines=summary_lines: _summary_figure(summary_lines),
            )
        )

        for channel in SNR_CHANNELS:
            requirement: list[tuple[float | None, str]] = []
            if channel.startswith("TOF"):
                requirement.append((3.0, "SNR ≥ 3"))
            else:
                requirement.append((6.0, "SNR ≥ 6"))
            data = list(self.snr.get(channel, []))
            stem = f"snr_{channel.replace(' ', '_').replace('/', '_').lower()}"
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"Distribution of signal-to-noise ratios for {channel}.",
                    label=f"fig:{stem}",
                    builder=lambda data=data, channel=channel, requirement=requirement: _hist_figure(
                        data,
                        title=f"{channel} Signal-to-Noise Ratio",
                        xlabel="SNR",
                        caption=f"Histogram of {channel} SNR across all analysed olivine events.",
                        requirement_lines=requirement,
                    ),
                )
            )

        for channel in TARGET_CHANNELS:
            stem = f"rise_{channel.replace(' ', '_').lower()}"
            data = list(self.rise.get(channel, []))
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"Fitted rise-time constants for the {channel} channel.",
                    label=f"fig:{stem}",
                    builder=lambda data=data, channel=channel: _hist_figure(
                        data,
                        title=f"{channel} Rise Time (τ_rise)",
                        xlabel="Rise time (µs)",
                        caption=f"Distribution of fitted rise constants for the {channel} waveform fits.",
                        requirement_lines=((0.1, "Spec lower"), (100.0, "Spec upper")),
                    ),
                )
            )

        for channel in TARGET_CHANNELS:
            stem = f"decay_{channel.replace(' ', '_').lower()}"
            data = list(self.decay.get(channel, []))
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"Fitted decay-time constants for the {channel} channel.",
                    label=f"fig:{stem}",
                    builder=lambda data=data, channel=channel: _hist_figure(
                        data,
                        title=f"{channel} Decay Time (τ_decay)",
                        xlabel="Decay time (µs)",
                        caption=f"Distribution of fitted decay constants for the {channel} waveform fits.",
                        requirement_lines=((1e-5, "Spec lower"), (10.0, "Spec upper")),
                    ),
                )
            )

        for pair_name in ("Ion Grid vs TOF H", "TOF H vs Target", "Target vs Ion Grid"):
            stem = f"trigger_{pair_name.replace(' ', '_').replace('/', '_').lower()}"
            deltas = list(self.trigger_deltas.get(pair_name, []))
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"Absolute trigger offset distribution for {pair_name} events.",
                    label=f"fig:{stem}",
                    builder=lambda deltas=deltas, pair_name=pair_name: _hist_figure(
                        deltas,
                        title=f"Trigger Delta: {pair_name}",
                        xlabel="|Δt| (µs)",
                        caption="Absolute timing offset between trigger references for the channel pair.",
                        requirement_lines=((20.0, "Spec upper"),),
                    ),
                )
            )

        for channel in TARGET_CHANNELS:
            stem = f"chi_{channel.replace(' ', '_').lower()}"
            data = list(self.chi_sq.get(channel, []))
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"χ² goodness-of-fit statistics for the {channel} waveform model fits.",
                    label=f"fig:{stem}",
                    builder=lambda data=data, channel=channel: _hist_figure(
                        data,
                        title=f"{channel} χ²",
                        xlabel="χ²",
                        caption=f"Goodness-of-fit χ² statistic for the {channel} model fits.",
                        requirement_lines=(),
                    ),
                )
            )

        for channel in TARGET_CHANNELS:
            stem = f"reduced_chi_{channel.replace(' ', '_').lower()}"
            data = list(self.reduced_chi_sq.get(channel, []))
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"Reduced χ² values for the {channel} waveform fits.",
                    label=f"fig:{stem}",
                    builder=lambda data=data, channel=channel: _hist_figure(
                        data,
                        title=f"{channel} Reduced χ²",
                        xlabel="Reduced χ²",
                        caption=f"Reduced χ² values for the {channel} waveform fits.",
                        requirement_lines=((1.0, "Ideal"),),
                    ),
                )
            )

        for channel in TARGET_CHANNELS + ("Ion Grid",):
            stem = f"impact_charge_{channel.replace(' ', '_').lower()}"
            data = list(self.impact_charge.get(channel, []))
            figure_specs.append(
                FigureSpec(
                    stem=stem,
                    caption=f"Impact charge distribution derived for the {channel} channel.",
                    label=f"fig:{stem}",
                    builder=lambda data=data, channel=channel: _hist_figure(
                        data,
                        title=f"{channel} Impact Charge",
                        xlabel="Charge (C)",
                        caption=f"Impact charge estimates for {channel} obtained from fitted waveforms.",
                        requirement_lines=None,
                    ),
                )
            )

        figure_specs.append(
            FigureSpec(
                stem="mass_abundance",
                caption="Median relative abundances for the expected olivine mass lines.",
                label="fig:mass-abundance",
                builder=lambda results=list(self.mass_results): _mass_abundance_figure(results),
            )
        )

        figure_specs.append(
            FigureSpec(
                stem="ternary",
                caption="Mg–Si–Fe ternary composition derived from mass-line fits.",
                label="fig:ternary",
                builder=lambda points=list(self._ternary_points): _ternary_figure(points),
            )
        )

        figure_specs.append(
            FigureSpec(
                stem="mass_vs_velocity",
                caption="Particle mass as a function of inferred impact velocity.",
                label="fig:mass-vs-velocity",
                builder=lambda records=list(self.particle_estimates): _mass_vs_velocity_figure(records),
            )
        )

        figure_specs.append(
            FigureSpec(
                stem="yield_vs_velocity",
                caption="Charge yield per unit mass as a function of velocity.",
                label="fig:yield-vs-velocity",
                builder=lambda records=list(self.particle_estimates): _yield_vs_velocity_figure(records),
            )
        )

        figure_specs.append(
            FigureSpec(
                stem="abundance_stack",
                caption="Stacked relative abundances of olivine mass lines vs. impact velocity.",
                label="fig:abundance-stack",
                builder=lambda results=list(self.mass_results), estimates=dict(self.event_estimates): _abundance_stack_figure(
                    results,
                    estimates,
                ),
            )
        )

        figure_specs.append(
            FigureSpec(
                stem="mass_line_probability",
                caption="Detection probability for each mass line as a function of velocity.",
                label="fig:mass-line-probability",
                builder=lambda results=list(self.mass_results), estimates=dict(self.event_estimates): _mass_line_probability_figure(
                    results,
                    estimates,
                ),
            )
        )

        figure_specs.append(
            FigureSpec(
                stem="saturation_counts",
                caption="Frequency of saturation flags by channel.",
                label="fig:saturation",
                builder=lambda counts=dict(self.saturation_counts): _saturation_bar_figure(counts),
            )
        )

        figures_dir = pdf_path.parent / f"{pdf_path.stem}_figures"
        figures_dir.mkdir(parents=True, exist_ok=True)

        figure_assets: list[FigureAsset] = []
        for spec in figure_specs:
            fig = spec.builder()
            asset_path = figures_dir / f"{spec.stem}.pdf"
            fig.savefig(asset_path, format="pdf", bbox_inches="tight")
            plt.close(fig)
            figure_assets.append(
                FigureAsset(path=asset_path, caption=spec.caption, label=spec.label)
            )

        def _ordered_stats(
            source: Mapping[str, Sequence[float]],
            preferred: Sequence[str] | None = None,
        ) -> dict[str, Mapping[str, float | int | None]]:
            stats_map: dict[str, Mapping[str, float | int | None]] = {}
            seen: set[str] = set()
            if preferred is not None:
                for key in preferred:
                    values = source.get(key, [])
                    stats_map[key] = _descriptive_stats(values)
                    seen.add(key)
            for key in sorted(source):
                if key in seen:
                    continue
                stats_map[key] = _descriptive_stats(source[key])
            return stats_map

        snr_stats = _ordered_stats(self.snr, SNR_CHANNELS)
        rise_stats = _ordered_stats(self.rise, TARGET_CHANNELS)
        decay_stats = _ordered_stats(self.decay, TARGET_CHANNELS)
        chi_stats = _ordered_stats(self.chi_sq, TARGET_CHANNELS)
        reduced_stats = _ordered_stats(self.reduced_chi_sq, TARGET_CHANNELS)
        impact_stats = _ordered_stats(self.impact_charge, TARGET_CHANNELS + ("Ion Grid",))
        trigger_stats = {
            pair: _descriptive_stats(self.trigger_deltas.get(pair, []))
            for pair in TRIGGER_PAIRS
        }

        particle_metrics: dict[str, Mapping[str, float | int | None]] = {
            "Velocity (km/s)": _descriptive_stats(
                [record.velocity_kms for record in self.particle_estimates]
            ),
            "Mass (kg)": _descriptive_stats(
                [record.mass_kg for record in self.particle_estimates]
            ),
            "Charge Yield (C/kg)": _descriptive_stats(
                [record.yield_c_per_kg for record in self.particle_estimates]
            ),
        }
        rise_samples = [record.rise_time_us for record in self.particle_estimates if record.rise_time_us is not None]
        if rise_samples:
            particle_metrics["Rise Time (µs)"] = _descriptive_stats(rise_samples)
        ratio_samples = [
            record.ion_to_target_ratio
            for record in self.particle_estimates
            if record.ion_to_target_ratio is not None
        ]
        if ratio_samples:
            particle_metrics["Ion/Target Charge Ratio"] = _descriptive_stats(ratio_samples)

        mass_line_stats = _relative_stats(self.mass_results)

        mass_fit_rows: list[list[str]] = []
        for result in self.mass_results:
            for fit in result.fits:
                if not fit.success:
                    continue
                rel = result.relative_abundances.get(fit.species, 0.0)
                rel_display = f"{rel * 100.0:.2f}\\%"
                mass_fit_rows.append(
                    [
                        _latex_escape(result.event_id),
                        _latex_escape(fit.species),
                        _format_value(fit.target_mass),
                        _format_value(fit.fit_mass),
                        _format_value(fit.sigma),
                        _format_value(fit.tau),
                        _format_value(fit.area),
                        rel_display,
                    ]
                )

        particle_rows: list[list[str]] = []
        for record in self.particle_estimates:
            particle_rows.append(
                [
                    _latex_escape(record.event_id),
                    _latex_escape(record.channel),
                    _format_value(record.charge_c),
                    _format_value(record.velocity_kms),
                    _format_value(record.mass_kg),
                    _format_value(record.yield_c_per_kg),
                    _format_value(record.rise_time_us),
                    _format_value(record.ion_to_target_ratio),
                    _latex_escape(record.velocity_source),
                ]
            )

        tex_path = pdf_path.with_suffix(".tex")
        latex_lines: list[str] = [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{graphicx}",
            r"\usepackage{booktabs}",
            r"\usepackage{longtable}",
            r"\usepackage{caption}",
            r"\usepackage{hyperref}",
            r"\hypersetup{colorlinks=true, linkcolor=blue, urlcolor=blue}",
            r"\begin{document}",
            r"\title{Olivine Metrics Scientific Report}",
            r"\author{SpectrumPY Flight Toolkit}",
            r"\date{\today}",
            r"\maketitle",
            r"\tableofcontents",
            r"\clearpage",
        ]

        latex_lines.append(r"\section{Overview}")
        latex_lines.append(r"\begin{itemize}")
        for line in summary_lines:
            latex_lines.append(rf"\item {_latex_escape(line)}")
        latex_lines.append(r"\end{itemize}")

        latex_lines.append(r"\section{Operational Metrics}")
        latex_lines.append(
            _render_count_table(
                "Saturation flag counts",
                "tab:saturation-counts",
                dict(self.saturation_counts),
            )
        )
        latex_lines.append(
            _render_count_table(
                "Preferred target channel usage",
                "tab:target-usage",
                dict(self.target_usage),
            )
        )

        latex_lines.append(r"\section{Channel Statistics}")
        latex_lines.append(
            _render_stats_table(
                "Signal-to-noise ratio statistics",
                "tab:snr",
                snr_stats,
            )
        )
        latex_lines.append(
            _render_stats_table(
                "Rise-time fit parameters",
                "tab:rise",
                rise_stats,
            )
        )
        latex_lines.append(
            _render_stats_table(
                "Decay-time fit parameters",
                "tab:decay",
                decay_stats,
            )
        )
        latex_lines.append(
            _render_stats_table(
                "χ² statistics",
                "tab:chi",
                chi_stats,
            )
        )
        latex_lines.append(
            _render_stats_table(
                "Reduced χ² statistics",
                "tab:reduced-chi",
                reduced_stats,
            )
        )
        latex_lines.append(
            _render_stats_table(
                "Trigger delta statistics",
                "tab:trigger-delta",
                trigger_stats,
            )
        )

        latex_lines.append(r"\section{Impact Charge and Particle Estimates}")
        latex_lines.append(
            _render_stats_table(
                "Impact charge statistics by channel",
                "tab:impact-charge",
                impact_stats,
            )
        )
        latex_lines.append(
            _render_stats_table(
                "Particle estimate summary",
                "tab:particle-summary",
                particle_metrics,
            )
        )
        latex_lines.append(
            _render_longtable(
                "Per-event particle estimates",
                "tab:particle-events",
                "llrrrrrrl",
                [
                    "Event",
                    "Channel",
                    "Charge (C)",
                    "Velocity (km/s)",
                    "Mass (kg)",
                    "Yield (C/kg)",
                    "Rise (µs)",
                    "Ion/Target",
                    "Velocity Source",
                ],
                particle_rows,
            )
        )

        latex_lines.append(r"\section{Mass Line Analysis}")
        latex_lines.append(
            _render_stats_table(
                "Relative abundance statistics for mass lines (fractional units)",
                "tab:mass-stats",
                mass_line_stats,
            )
        )
        latex_lines.append(
            _render_longtable(
                "Successful mass-line fits by event",
                "tab:mass-fits",
                "llrrrrrr",
                [
                    "Event",
                    "Species",
                    "Target Mass",
                    "Fit Mass",
                    "σ",
                    "τ",
                    "Area",
                    "Rel. Abundance",
                ],
                mass_fit_rows,
            )
        )

        latex_lines.append(r"\section{Figures}")
        for asset in figure_assets:
            rel_path = Path(os.path.relpath(asset.path, tex_path.parent)).as_posix()
            latex_lines.extend(
                [
                    r"\begin{figure}[htbp]",
                    r"\centering",
                    rf"\includegraphics[width=0.85\linewidth]{{{rel_path}}}",
                    rf"\caption{{{_latex_escape(asset.caption)}}}",
                    rf"\label{{{asset.label}}}",
                    r"\end{figure}",
                ]
            )

        latex_lines.append(r"\end{document}")

        tex_path.write_text("\n".join(latex_lines), encoding="utf-8")

        compiled_pdf = _compile_latex(tex_path)
        if compiled_pdf is not None:
            if compiled_pdf != pdf_path:
                shutil.move(compiled_pdf, pdf_path)
            for extension in (".aux", ".log", ".out"):
                aux = tex_path.with_suffix(extension)
                if aux.exists():
                    try:
                        aux.unlink()
                    except Exception:
                        pass
        else:
            print(
                "Warning: LaTeX engine not found; falling back to Matplotlib-generated PDF report.",
            )
            with PdfPages(pdf_path) as pdf:
                summary_fig = _summary_figure(summary_lines)
                pdf.savefig(summary_fig)
                plt.close(summary_fig)
                for spec in figure_specs:
                    fig = spec.builder()
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
            "mass_analysis": {
                "events": [result.to_dict() for result in self.mass_results],
                "relative_abundance_stats": {
                    species: {
                        "count": stat_values["count"],
                        "median": stat_values["median"],
                        "mean": stat_values["mean"],
                        "std": stat_values["std"],
                    }
                    for species, stat_values in _relative_stats(self.mass_results).items()
                },
                "ternary_points": [
                    {"Mg": float(mg), "Si": float(si), "Fe": float(fe)}
                    for mg, si, fe in self._ternary_points
                ],
                "median_elemental_percent": (
                    {
                        "Mg": elemental[0],
                        "Si": elemental[1],
                        "Fe": elemental[2],
                    }
                    if (elemental := _aggregate_elemental_fractions(self.mass_results)) is not None
                    else None
                ),
            },
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
