"""Utilities to estimate dust particle velocity and mass.

This module implements the piecewise power-law fits published for IDEX
calibrations (E. Grun et al., 1992) to recover dust impact velocities and
masses from the target rise time, the collection efficiency (ion grid/target
ratio), and the measured target charge amplitude.

The calibration coefficients are distributed alongside the code under the
``lookup`` directory.  They can be extended at runtime by providing custom CSV
files that follow the layout of ``yield.csv`` and ``t_rise.csv``.  The helper
functions in this file remain agnostic to the coefficient source, accepting any
mapping that resolves the seven parameters of the fit in
Eq. 9 of the attached documentation.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Union

import csv
import math

Number = Union[int, float]
ScalarOrSequence = Union[Number, Sequence[Number]]


@dataclass(frozen=True)
class CurveParameters:
    """Parameterisation of the piecewise power-law curve used in the fits."""

    A: float
    a1: float
    a2: float
    a3: float
    vb: float
    vc: float
    k: float

    def evaluate(self, x: ScalarOrSequence) -> Union[float, List[float]]:
        """Evaluate the calibration curve at ``x``.

        Parameters
        ----------
        x:
            Scalar or sequence of positive values.  Scalars return a float,
            sequences return a list of floats.
        """

        def _eval_scalar(value: float) -> float:
            if value <= 0:
                raise ValueError("Curve evaluation is only defined for positive inputs.")
            base = self.A * value ** self.a1
            term_b = 1.0 + (value / self.vb) ** self.k if self.vb else 1.0
            term_c = 1.0 + (value / self.vc) ** self.k if self.vc else 1.0
            exp_b = (self.a2 - self.a1) / self.k
            exp_c = (self.a3 - self.a2) / self.k
            return base * term_b ** exp_b * term_c ** exp_c

        if isinstance(x, (int, float)):
            return _eval_scalar(float(x))
        return [_eval_scalar(float(value)) for value in x]


class CoefficientTable:
    """Load and serve coefficient sets from CSV files."""

    REQUIRED_FIELDS: Tuple[str, ...] = ("A", "a1", "a2", "a3", "vb", "vc", "k")

    def __init__(self, rows: Mapping[str, CurveParameters]):
        self._rows = dict(rows)

    @classmethod
    def from_csv(cls, path: Union[str, Path], default_label: str = "combined") -> "CoefficientTable":
        with Path(path).open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise ValueError(f"CSV file {path} has no header.")
            label_column = cls._identify_label_column(reader.fieldnames)
            rows: MutableMapping[str, CurveParameters] = {}
            auto_index = 0
            for raw_row in reader:
                if all((value or "").strip() == "" for value in raw_row.values()):
                    continue
                label = (raw_row.get(label_column, "") or "").strip() if label_column else ""
                if not label:
                    label = default_label if not rows else f"{default_label}_{auto_index}"
                    auto_index += 1
                params = cls._row_to_parameters(raw_row)
                rows[label.lower()] = params
        return cls(rows)

    @staticmethod
    def _identify_label_column(columns: Sequence[str]) -> Optional[str]:
        normalised = [CoefficientTable._normalise_header(name) for name in columns]
        for original, norm in zip(columns, normalised):
            if norm in {"material", "label", "name"}:
                return original
        return None

    @staticmethod
    def _normalise_header(name: str) -> str:
        normalised = (
            name.lower()
            .replace(" ", "")
            .replace("μ", "u")
            .replace("[", "")
            .replace("]", "")
            .replace("/", "")
            .replace("Δ", "delta")
            .replace("σ", "sigma")
            .replace("₁", "1")
            .replace("₂", "2")
            .replace("₃", "3")
            .replace("_", "")
        )
        if normalised.startswith("vb"):
            return "vb"
        if normalised.startswith("vc"):
            return "vc"
        return normalised

    @classmethod
    def _row_to_parameters(cls, row: Mapping[str, str]) -> CurveParameters:
        extracted: Dict[str, float] = {}
        for key, value in row.items():
            if value is None or value.strip() == "":
                continue
            normalised = cls._normalise_header(key)
            if normalised in {"a", "a0"}:
                extracted["A"] = float(value)
                continue
            if normalised in cls.REQUIRED_FIELDS:
                extracted[normalised] = float(value)
        missing = [field for field in cls.REQUIRED_FIELDS if field not in extracted]
        if missing:
            raise ValueError(f"Coefficient row missing required fields: {missing}")
        return CurveParameters(
            A=extracted["A"],
            a1=extracted["a1"],
            a2=extracted["a2"],
            a3=extracted["a3"],
            vb=extracted["vb"],
            vc=extracted["vc"],
            k=extracted["k"],
        )

    def get(self, label: str) -> CurveParameters:
        try:
            return self._rows[label.lower()]
        except KeyError as exc:
            raise KeyError(f"Unknown coefficient set '{label}'. Available: {sorted(self._rows)}") from exc

    def labels(self) -> Sequence[str]:
        return sorted(self._rows)


@dataclass
class VelocityEstimates:
    rise_time: Optional[float]
    collection_efficiency: Optional[float]
    selected: float
    source: str


@dataclass
class ParticleEstimate:
    charge_c: float
    velocity_kms: float
    mass_kg: float
    yield_c_per_kg: float
    velocity_details: VelocityEstimates


def compute_velocity_from_rise_time(
    rise_time: Optional[Number],
    params: CurveParameters,
) -> Optional[float]:
    if rise_time is None:
        return None
    return float(params.evaluate(float(rise_time)))


def compute_velocity_from_collection_efficiency(
    ratio: Optional[Number],
    params: CurveParameters,
) -> Optional[float]:
    if ratio is None:
        return None
    return float(params.evaluate(float(ratio)))


def select_velocity(
    rise_velocity: Optional[float],
    ratio_velocity: Optional[float],
    valid_range: Tuple[float, float] = (1.0, 100.0),
) -> Tuple[float, str]:
    candidates: List[Tuple[str, float]] = []
    for source, value in ("rise_time", rise_velocity), ("collection_efficiency", ratio_velocity):
        if value is None or not math.isfinite(value):
            continue
        if valid_range[0] <= value <= valid_range[1]:
            candidates.append((source, value))
    if len(candidates) == 2:
        avg = sum(value for _, value in candidates) / 2.0
        return avg, "average"
    if candidates:
        source, value = candidates[0]
        return value, source
    for source, value in ("rise_time", rise_velocity), ("collection_efficiency", ratio_velocity):
        if value is None or not math.isfinite(value):
            continue
        return max(min(value, valid_range[1]), valid_range[0]), f"clamped_{source}"
    raise ValueError("No valid velocity estimates were provided.")


def compute_mass_from_charge(
    charge_c: Number,
    velocity_kms: float,
    yield_params: CurveParameters,
) -> Tuple[float, float]:
    if charge_c <= 0:
        raise ValueError("Charge must be positive to compute mass.")
    yield_value = float(yield_params.evaluate(float(velocity_kms)))
    if yield_value <= 0:
        raise ValueError("Yield evaluation returned a non-positive value.")
    mass = float(charge_c) / yield_value
    return mass, yield_value


def estimate_particle(
    *,
    charge_c: Number,
    rise_time: Optional[Number],
    ion_to_target_ratio: Optional[Number],
    rise_params: CurveParameters,
    ratio_params: Optional[CurveParameters],
    yield_params: CurveParameters,
    velocity_range: Tuple[float, float] = (1.0, 100.0),
) -> ParticleEstimate:
    rise_velocity = compute_velocity_from_rise_time(rise_time, rise_params)
    ratio_velocity = (
        compute_velocity_from_collection_efficiency(ion_to_target_ratio, ratio_params)
        if ratio_params is not None
        else None
    )
    selected_velocity, source = select_velocity(rise_velocity, ratio_velocity, velocity_range)
    mass, charge_yield = compute_mass_from_charge(charge_c, selected_velocity, yield_params)
    return ParticleEstimate(
        charge_c=float(charge_c),
        velocity_kms=selected_velocity,
        mass_kg=mass,
        yield_c_per_kg=charge_yield,
        velocity_details=VelocityEstimates(
            rise_time=rise_velocity,
            collection_efficiency=ratio_velocity,
            selected=selected_velocity,
            source=source,
        ),
    )


def load_default_tables() -> Tuple[CoefficientTable, Optional[CoefficientTable], CoefficientTable]:
    base = Path(__file__).parent
    rise_table = CoefficientTable.from_csv(base / "t_rise.csv")
    ratio_path = base / "collection_efficiency.csv"
    ratio_table = CoefficientTable.from_csv(ratio_path) if ratio_path.exists() else None
    yield_table = CoefficientTable.from_csv(base / "yield.csv")
    return rise_table, ratio_table, yield_table


def _parse_optional_float(value: str) -> Optional[float]:
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    return float(value)


def estimate_from_csv(
    input_rows: Iterable[Mapping[str, str]],
    *,
    rise_params: CurveParameters,
    ratio_params: Optional[CurveParameters],
    yield_params: CurveParameters,
    velocity_range: Tuple[float, float] = (1.0, 100.0),
) -> Iterator[Tuple[Mapping[str, str], ParticleEstimate]]:
    for row in input_rows:
        charge = float(row["target_charge_c"])
        rise_time = _parse_optional_float(row.get("target_rise_time_us"))
        ratio = _parse_optional_float(row.get("ion_to_target_ratio"))
        estimate = estimate_particle(
            charge_c=charge,
            rise_time=rise_time,
            ion_to_target_ratio=ratio,
            rise_params=rise_params,
            ratio_params=ratio_params,
            yield_params=yield_params,
            velocity_range=velocity_range,
        )
        yield row, estimate


def write_estimates_csv(
    rows_with_estimates: Iterable[Tuple[Mapping[str, str], ParticleEstimate]],
    *,
    output_path: Union[str, Path],
) -> None:
    rows_with_estimates = list(rows_with_estimates)
    if not rows_with_estimates:
        raise ValueError("No rows to write.")
    original_headers = list(rows_with_estimates[0][0].keys())
    extra_headers = [
        "velocity_km_s",
        "mass_kg",
        "charge_yield_c_per_kg",
        "velocity_rise_time_km_s",
        "velocity_collection_efficiency_km_s",
        "velocity_selection_source",
    ]
    with Path(output_path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=original_headers + extra_headers)
        writer.writeheader()
        for original_row, estimate in rows_with_estimates:
            output = dict(original_row)
            output.update(
                velocity_km_s=f"{estimate.velocity_kms:.6g}",
                mass_kg=f"{estimate.mass_kg:.6g}",
                charge_yield_c_per_kg=f"{estimate.yield_c_per_kg:.6g}",
                velocity_rise_time_km_s=(
                    f"{estimate.velocity_details.rise_time:.6g}"
                    if estimate.velocity_details.rise_time is not None
                    else ""
                ),
                velocity_collection_efficiency_km_s=(
                    f"{estimate.velocity_details.collection_efficiency:.6g}"
                    if estimate.velocity_details.collection_efficiency is not None
                    else ""
                ),
                velocity_selection_source=estimate.velocity_details.source,
            )
            writer.writerow(output)


def _load_input_csv(path: Union[str, Path]) -> List[Mapping[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if not rows:
            raise ValueError("Input CSV is empty.")
        required = {"target_charge_c", "target_rise_time_us", "ion_to_target_ratio"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input CSV missing required columns: {sorted(missing)}")
        return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Estimate dust particle velocity and mass.")
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        help="CSV with target charge, rise time, and ion/target ratio",
    )
    parser.add_argument(
        "-f",
        "--file",
        "--filename",
        dest="input_file",
        type=Path,
        help="Explicit CSV containing target charge, rise time, and ion/target ratio",
    )
    parser.add_argument("output", type=Path, help="Destination CSV with estimates appended")
    parser.add_argument("--rise-coeff", default="combined", help="Coefficient label for the rise-time curve")
    parser.add_argument(
        "--ratio-coeff",
        default="combined",
        help="Coefficient label for the collection-efficiency curve (ignored if unavailable)",
    )
    parser.add_argument("--yield-coeff", default="combined", help="Coefficient label for the charge-yield curve")
    parser.add_argument("--rise-table", type=Path, default=Path(__file__).with_name("t_rise.csv"))
    parser.add_argument(
        "--ratio-table",
        type=Path,
        default=Path(__file__).with_name("collection_efficiency.csv"),
        help="CSV file with ion/target ratio coefficients (optional)",
    )
    parser.add_argument("--yield-table", type=Path, default=Path(__file__).with_name("yield.csv"))
    parser.add_argument(
        "--min-velocity",
        type=float,
        default=1.0,
        help="Minimum acceptable velocity before clamping [km/s]",
    )
    parser.add_argument(
        "--max-velocity",
        type=float,
        default=100.0,
        help="Maximum acceptable velocity before clamping [km/s]",
    )
    args = parser.parse_args(argv)

    rise_table = CoefficientTable.from_csv(args.rise_table, default_label="combined")
    rise_params = rise_table.get(args.rise_coeff)

    ratio_params = None
    if args.ratio_table.exists():
        ratio_table = CoefficientTable.from_csv(args.ratio_table, default_label="combined")
        if args.ratio_coeff.lower() in ratio_table.labels():
            ratio_params = ratio_table.get(args.ratio_coeff)
        else:
            raise KeyError(
                f"Coefficient '{args.ratio_coeff}' not found in {args.ratio_table}. "
                f"Available: {', '.join(ratio_table.labels())}"
            )

    yield_table = CoefficientTable.from_csv(args.yield_table, default_label="combined")
    yield_params = yield_table.get(args.yield_coeff)

    input_path = args.input_file or args.input
    if input_path is None:
        raise SystemExit("An input CSV must be supplied using -f/--file/--filename or as a positional argument.")

    rows = _load_input_csv(input_path)
    estimates = list(
        estimate_from_csv(
            rows,
            rise_params=rise_params,
            ratio_params=ratio_params,
            yield_params=yield_params,
            velocity_range=(args.min_velocity, args.max_velocity),
        )
    )
    write_estimates_csv(estimates, output_path=args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
