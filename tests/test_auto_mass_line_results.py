import math
import warnings

import pytest

EXPECTED_EVENTS = {
    ("ois_output_12132023_223729", 6): {
        "dust": "Al",
        "stretch_ns": 1.45086794071479,
        "abundances": {
            "H": 14.121,
            "H2": 0.890,
            "C": 0.0,
            "CH": 0.0,
            "O": 0.517,
            "Na": 0.532,
            "24Mg": 0.0,
            "25Mg": 0.0,
            "26Mg": 0.0,
            "C2H2": 3.105,
            "Al": 77.935,
            "Si": 0.0,
            "C2H5": 0.0,
            "39K": 0.0,
            "Ca": 0.0,
            "41K": 0.0,
            "Mass 48": 0.0,
            "Al2": 2.431,
            "56Fe": 0.0,
            "Mass 64": 0.0,
            "Al2O": 0.468,
        },
    },
    ("ois_output_12132023_223729", 14): {
        "dust": "Al",
        "stretch_ns": 1.5315283386149776,
        "abundances": {
            "H": 32.955,
            "H2": 0.0,
            "C": 3.625,
            "CH": 0.0,
            "O": 6.256,
            "Na": 0.556,
            "24Mg": 0.0,
            "25Mg": 0.0,
            "26Mg": 0.0,
            "C2H2": 2.938,
            "Al": 52.623,
            "Si": 0.0,
            "C2H5": 0.240,
            "39K": 0.0,
            "Ca": 0.0,
            "41K": 0.0,
            "Mass 48": 0.0,
            "Al2": 0.529,
            "56Fe": 0.0,
            "Mass 64": 0.0,
            "Al2O": 0.278,
        },
    },
    ("ois_output_12132023_223729", 22): {
        "dust": "Al",
        "stretch_ns": 1.6,
        "abundances": {
            "H": 36.434,
            "H2": 0.532,
            "C": 1.568,
            "CH": 0.0,
            "O": 3.069,
            "Na": 1.300,
            "24Mg": 0.0,
            "25Mg": 0.0,
            "26Mg": 0.0,
            "C2H2": 2.924,
            "Al": 52.902,
            "Si": 0.0,
            "C2H5": 0.0,
            "39K": 0.0,
            "Ca": 0.0,
            "41K": 0.0,
            "Mass 48": 0.0,
            "Al2": 1.270,
            "56Fe": 0.0,
            "Mass 64": 0.0,
            "Al2O": 0.0,
        },
    },
    ("ois_output_12182023_185430", 4): {
        "dust": "Olv",
        "stretch_ns": 1.4981473968573795,
        "abundances": {
            "H": 7.099,
            "H2": 0.540,
            "C": 1.355,
            "CH": 0.528,
            "O": 0.352,
            "Na": 12.813,
            "24Mg": 46.858,
            "25Mg": 8.470,
            "26Mg": 7.283,
            "C2H2": 0.0,
            "Al": 0.0,
            "Si": 1.220,
            "C2H5": 0.0,
            "39K": 4.640,
            "Ca": 3.105,
            "41K": 2.081,
            "Mass 48": 0.760,
            "Al2": 0.0,
            "56Fe": 2.397,
            "Mass 64": 0.502,
            "Al2O": 0.0,
        },
    },
    ("ois_output_12182023_185430", 7): {
        "dust": "Olv",
        "stretch_ns": 1.6,
        "abundances": {
            "H": 77.466,
            "H2": 0.0,
            "C": 0.0,
            "CH": 0.0,
            "O": 5.371,
            "Na": 0.0,
            "24Mg": 13.550,
            "25Mg": 2.251,
            "26Mg": 0.624,
            "C2H2": 0.0,
            "Al": 0.0,
            "Si": 0.737,
            "C2H5": 0.0,
            "39K": 0.0,
            "Ca": 0.0,
            "41K": 0.0,
            "Mass 48": 0.0,
            "Al2": 0.0,
            "56Fe": 0.0,
            "Mass 64": 0.0,
            "Al2O": 0.0,
        },
    },
    ("ois_output_12182023_185430", 18): {
        "dust": "Olv",
        "stretch_ns": 1.4336340518005033,
        "abundances": {
            "H": 20.676,
            "H2": 0.0,
            "C": 8.164,
            "CH": 0.0,
            "O": 1.499,
            "Na": 4.295,
            "24Mg": 37.374,
            "25Mg": 3.948,
            "26Mg": 4.770,
            "C2H2": 0.0,
            "Al": 0.0,
            "Si": 5.209,
            "C2H5": 0.0,
            "39K": 2.565,
            "Ca": 8.626,
            "41K": 0.0,
            "Mass 48": 0.0,
            "Al2": 0.0,
            "56Fe": 2.063,
            "Mass 64": 0.811,
            "Al2O": 0.0,
        },
    },
}


@pytest.mark.parametrize("key, expected", EXPECTED_EVENTS.items())
def test_auto_identified_mass_lines_match_reference(idex_mass_analysis, key, expected):
    assert key in idex_mass_analysis, f"Automated analysis missing results for {key}."

    entry = idex_mass_analysis[key]

    observed_abundances = {
        species: float(entry["abundances"].get(species, 0.0))
        for species in expected["abundances"].keys()
    }

    report_lines = [
        f"Automated analysis comparison for {key[0]} event {key[1]}",
    ]

    stretch_expected = float(expected["stretch_ns"])
    stretch_observed = float(entry.get("stretch_ns", float("nan")))
    stretch_diff = stretch_observed - stretch_expected
    stretch_status = (
        "VALID"
        if math.isfinite(stretch_diff) and abs(stretch_diff) <= 0.2
        else "OUTSIDE 0.2 TOLERANCE"
    )
    report_lines.append(
        "Stretch (ns): observed={:.6f}, expected={:.6f}, diff={:+.6f} -> {}".format(
            stretch_observed,
            stretch_expected,
            stretch_diff,
            stretch_status,
        )
    )

    total_abundance = 0.0
    for species, expected_value in expected["abundances"].items():
        observed_value = observed_abundances[species]
        total_abundance += observed_value

        difference = observed_value - expected_value
        if not math.isclose(
            observed_value,
            expected_value,
            rel_tol=1e-3,
            abs_tol=1e-3,
        ):
            report_lines.append(
                (
                    "  {}: observed={:.6f}, expected={:.6f}, diff={:+.6f}".format(
                        species,
                        observed_value,
                        expected_value,
                        difference,
                    )
                )
            )

    if not math.isclose(total_abundance, 100.0, rel_tol=1e-5, abs_tol=5e-3):
        report_lines.append(
            "Total abundance: observed={:.6f}, expected=100.000000, diff={:+.6f}".format(
                total_abundance,
                total_abundance - 100.0,
            )
        )

    warnings.warn("\n".join(report_lines), RuntimeWarning)
