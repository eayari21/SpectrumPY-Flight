#!/usr/bin/env python3
"""Generate a human-readable summary of selected variables inside an IDEX L2A CDF.

The script relies on the external ``strings`` utility (available in the base
container) to avoid a hard dependency on ``cdflib``.  It is intentionally
light-weight so it can run in restricted environments where installing the
NASA CDF runtime is not possible.

Example usage::

    python tools/cdf_strings_report.py CDF/IDEX_Pre_Launch_CDFs/imap_idex_l2a_sci-1week_20231219_v999.cdf \
        --variables tof_peak_fit_parameters target_low_fit_results

The output is Markdown that can be redirected to a file."""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Sequence

ATTRIBUTE_KEYS: tuple[str, ...] = (
    "CATDESC",
    "DEPEND_0",
    "DEPEND_1",
    "DEPEND_2",
    "DEPEND_3",
    "DEPEND_4",
    "DICT_KEY",
    "DISPLAY_TYPE",
    "FIELDNAM",
    "FORMAT",
    "FILLVAL",
    "LABLAXIS",
    "SCALETYP",
    "UNITS",
    "VAR_NOTES",
)

DEFAULT_VARIABLES: tuple[str, ...] = (
    "epoch",
    "time_high_sample_rate",
    "time_low_sample_rate",
    "time_high_sample_rate_index",
    "time_low_sample_rate_index",
    "tof_high",
    "tof_mid",
    "tof_low",
    "tof_peak_fit_parameters",
    "tof_peak_area_under_fit",
    "tof_peak_chi_square",
    "tof_peak_reduced_chi_square",
    "tof_peak_kappa",
    "tof_snr",
    "target_fit_parameter_index",
    "target_fit_parameter_labels",
    "target_low",
    "target_low_fit_parameters",
    "target_low_fit_results",
    "target_low_impact_charge",
    "target_low_velocity_estimate",
    "target_low_dust_mass_estimate",
    "target_low_chi_squared",
    "target_low_reduced_chi_squared",
    "target_high",
    "target_high_fit_parameters",
    "target_high_fit_results",
    "target_high_impact_charge",
    "target_high_velocity_estimate",
    "target_high_dust_mass_estimate",
    "target_high_chi_squared",
    "target_high_reduced_chi_squared",
    "ion_grid",
    "ion_grid_fit_parameters",
    "ion_grid_fit_results",
    "ion_grid_impact_charge",
    "ion_grid_velocity_estimate",
    "ion_grid_dust_mass_estimate",
    "ion_grid_chi_squared",
    "ion_grid_reduced_chi_squared",
)


def _iter_strings(cdf_path: Path) -> Iterable[str]:
    if shutil.which("strings") is None:
        raise SystemExit("The 'strings' utility is required but not available.")
    proc = subprocess.Popen(
        ["strings", "-n", "1", str(cdf_path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    try:
        for raw_line in proc.stdout:
            yield raw_line.rstrip("\n")
    finally:
        proc.stdout.close()
        proc.wait()


def _collect_strings(cdf_path: Path) -> list[str]:
    return [line.strip() for line in _iter_strings(cdf_path)]


def _first_index(haystack: Sequence[str], needle: str) -> int | None:
    needle_lower = needle.lower()
    for idx, value in enumerate(haystack):
        if needle_lower in value.lower():
            return idx
    return None


def _is_metadata_value(token: str) -> bool:
    token = token.strip()
    if not token or token in ATTRIBUTE_KEYS:
        return False
    printable = ''.join(ch for ch in token if ch.isprintable())
    if printable != token:
        return False
    if token.count(' ') > 16:
        return False
    if any(ch.isalpha() for ch in token):
        return True
    numeric = token.replace('.', '').replace('-', '').replace('E', '').replace('+', '')
    if numeric.isdigit() and len(numeric) > 1:
        return True
    return False


def _clean_context_token(token: str) -> str | None:
    token = token.strip()
    if not token:
        return None
    if len(token) == 1 and not token.isalnum():
        return None
    printable = ''.join(ch for ch in token if ch.isprintable())
    if printable != token:
        return None
    if token in ATTRIBUTE_KEYS:
        return None
    return token


def _extract(strings_dump: Sequence[str], variable: str, window: int = 200, context_size: int = 12) -> tuple[dict[str, str], list[str]]:
    start = _first_index(strings_dump, variable)
    if start is None:
        return {}, []
    meta: dict[str, str] = {}
    context: list[str] = []
    slice_end = min(len(strings_dump), start + window)
    window_slice = strings_dump[start:slice_end]
    for idx, token in enumerate(window_slice):
        if token in ATTRIBUTE_KEYS:
            for candidate in window_slice[idx + 1:]:
                if _is_metadata_value(candidate):
                    meta.setdefault(token, candidate.strip())
                    break
        elif len(context) < context_size:
            cleaned = _clean_context_token(token)
            if cleaned:
                context.append(cleaned)
    return meta, context


def _format_section(name: str, metadata: dict[str, str], context: list[str]) -> str:
    lines = [f"### {name}"]
    if metadata:
        for key in ATTRIBUTE_KEYS:
            if key in metadata:
                lines.append(f"* **{key}**: {metadata[key]}")
    if context:
        lines.append("*Context sample:* " + ", ".join(context[:12]))
    if len(lines) == 1:
        lines.append("*No metadata located in the ASCII scan.*")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cdf", type=Path, help="Path to the CDF file to inspect")
    parser.add_argument(
        "--variables",
        nargs="*",
        default=DEFAULT_VARIABLES,
        help="Variable names to summarize",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=220,
        help="Number of strings to examine after the variable name",
    )
    args = parser.parse_args()

    if not args.cdf.exists():
        raise SystemExit(f"CDF file not found: {args.cdf}")

    strings_dump = _collect_strings(args.cdf)
    seen: set[str] = set()
    print(f"## Metadata summary for {args.cdf.name}\n")
    for variable in args.variables:
        if variable in seen:
            continue
        seen.add(variable)
        metadata, context = _extract(strings_dump, variable, window=args.window)
        print(_format_section(variable, metadata, context))
        print()


if __name__ == "__main__":
    main()
