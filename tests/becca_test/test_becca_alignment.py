from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import pytest

# The tabulated reference data compiled by Becca for the same events.
REFERENCE_TABLE_PATH = Path(__file__).with_name(
    "Spectrum_vs_SpectrumPY_SPECTRUM_VALUES.txt"
)
assert REFERENCE_TABLE_PATH.exists(), "Becca's reference table is missing."

EXPECTED_SPECIES = (
    "H",
    "H2",
    "C",
    "CH",
    "O",
    "Na",
    "24Mg",
    "25Mg",
    "26Mg",
    "C2H2",
    "Al",
    "Si",
    "C2H5",
    "39K",
    "Ca",
    "41K",
    "Mass 48",
    "Al2",
    "56Fe",
    "Mass 64",
    "Al2O",
)

@pytest.fixture(scope="module")
def reference_results(idex_mass_analysis) -> Dict[Tuple[str, int], dict]:
    return idex_mass_analysis


@pytest.fixture(scope="module")
def becca_reference() -> Dict[Tuple[str, int], dict]:
    lines = [line for line in REFERENCE_TABLE_PATH.read_text().splitlines() if line.strip()]
    assert lines, "Becca's reference table is empty."

    data: Dict[Tuple[str, int], dict] = {}
    for line in lines[1:]:
        parts = line.split()
        expected_length = 4 + len(EXPECTED_SPECIES)
        assert len(parts) == expected_length, (
            "Reference row has an unexpected column count; "
            "did the table structure change?"
        )
        data_file, dust, event_id, stretch = parts[:4]
        key = (data_file, int(event_id))
        species_values = [float(value) for value in parts[4:]]
        data[key] = {
            "dust": dust,
            "stretch_ns": float(stretch),
            "abundances": dict(zip(EXPECTED_SPECIES, species_values)),
        }
    return data


@pytest.mark.parametrize(
    "event_key",
    [
        ("ois_output_12132023_223729", 6),
        ("ois_output_12132023_223729", 14),
        ("ois_output_12132023_223729", 22),
        ("ois_output_12182023_185430", 4),
        ("ois_output_12182023_185430", 7),
        ("ois_output_12182023_185430", 18),
    ],
)
def test_regression_matches_becca_reference(event_key, reference_results, becca_reference):
    assert event_key in reference_results, f"Regression results missing {event_key}."
    assert event_key in becca_reference, f"Becca reference missing {event_key}."

    regression = reference_results[event_key]
    becca = becca_reference[event_key]

    assert regression["stretch_ns"] == pytest.approx(
        becca["stretch_ns"], rel=1e-6, abs=1e-3
    )

    regression_abundances = {
        species: float(regression["abundances"].get(species, 0.0))
        for species in EXPECTED_SPECIES
    }
    becca_abundances = {
        species: float(becca["abundances"].get(species, 0.0))
        for species in EXPECTED_SPECIES
    }

    max_difference = 0.0
    for species in EXPECTED_SPECIES:
        difference = abs(regression_abundances[species] - becca_abundances[species])
        max_difference = max(max_difference, difference)
        assert regression_abundances[species] == pytest.approx(
            becca_abundances[species], rel=1e-3, abs=5e-3
        )

    # The automated pipeline reproduces Becca's reference table to within numerical
    # precision, so highlight the achieved maximum deviation in case the regression
    # drifts in the future.
    assert max_difference == pytest.approx(0.0, abs=5e-4)
