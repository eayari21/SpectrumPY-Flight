"""Helper utilities for annotating Mg–Si–Fe ternary diagrams with reference lines."""

from __future__ import annotations

from typing import Iterable, Mapping, Sequence, Tuple


# Each entry contains (label, start_composition, end_composition) where the
# compositions are provided as relative atom counts for Mg, Si, and Fe before
# normalisation to barycentric space.
_MG_SI_FE_MINERAL_SERIES: Tuple[
    Tuple[str, Mapping[str, float], Mapping[str, float]],
    ...,
] = (
    (
        "Olivine (Fo–Fa)",
        {"Mg": 2.0, "Si": 1.0, "Fe": 0.0},
        {"Mg": 0.0, "Si": 1.0, "Fe": 2.0},
    ),
    (
        "Orthopyroxene (En–Fs)",
        {"Mg": 1.0, "Si": 1.0, "Fe": 0.0},
        {"Mg": 0.0, "Si": 1.0, "Fe": 1.0},
    ),
    (
        "Fo–En tie line",
        {"Mg": 2.0, "Si": 1.0, "Fe": 0.0},
        {"Mg": 1.0, "Si": 1.0, "Fe": 0.0},
    ),
)


def _normalise_endpoint(
    endpoint: Mapping[str, float],
    element_order: Sequence[str],
) -> Tuple[float, float, float]:
    values = [max(float(endpoint.get(symbol, 0.0)), 0.0) for symbol in element_order]
    total = sum(values)
    if total <= 0.0:
        return (0.0, 0.0, 0.0)
    return (values[0] / total, values[1] / total, values[2] / total)


def iter_mg_si_fe_series(
    element_order: Sequence[str],
) -> Iterable[Tuple[str, Tuple[float, float, float], Tuple[float, float, float]]]:
    """Yield Mg–Si–Fe mineral reference lines in the requested element order.

    Parameters
    ----------
    element_order:
        A three-element sequence giving the symbol mapped to each barycentric
        axis (e.g. ("Mg", "Si", "Fe")).  The returned tuples are normalised to
        this order so they can be passed directly to barycentric conversion
        helpers.
    """

    if set(element_order) != {"Mg", "Si", "Fe"}:
        return

    for label, start, end in _MG_SI_FE_MINERAL_SERIES:
        start_triplet = _normalise_endpoint(start, element_order)
        end_triplet = _normalise_endpoint(end, element_order)
        if not any(start_triplet) or not any(end_triplet):
            continue
        yield label, start_triplet, end_triplet


__all__ = ["iter_mg_si_fe_series"]

