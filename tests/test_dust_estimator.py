import math
import pytest

from spectrumpy_flight.lookup.dust_estimator import (
    CurveParameters,
    ParticleEstimate,
    compute_mass_from_charge,
    estimate_particle,
    select_velocity,
)


def make_flat_curve(value: float) -> CurveParameters:
    return CurveParameters(A=value, a1=0.0, a2=0.0, a3=0.0, vb=1e9, vc=1e9, k=1.0)


def test_curve_parameters_evaluate_scalar_and_sequence():
    params = CurveParameters(A=2.0, a1=1.0, a2=1.0, a3=1.0, vb=1e9, vc=1e9, k=1.0)
    assert math.isclose(params.evaluate(3.0), 6.0)
    result = params.evaluate([1.0, 2.0, 3.0])
    assert result == [2.0, 4.0, 6.0]


def test_select_velocity_prefers_average_when_both_valid():
    value, source = select_velocity(12.0, 8.0)
    assert source == "average"
    assert math.isclose(value, 10.0)


def test_select_velocity_clamps_when_outside_range():
    value, source = select_velocity(0.2, None)
    assert source == "clamped_rise_time"
    assert value == 1.0


def test_estimate_particle_averages_and_computes_mass():
    rise_params = make_flat_curve(10.0)
    ratio_params = make_flat_curve(12.0)
    yield_params = make_flat_curve(1.0)

    estimate = estimate_particle(
        charge_c=2.0,
        rise_time=1.0,
        ion_to_target_ratio=1.0,
        rise_params=rise_params,
        ratio_params=ratio_params,
        yield_params=yield_params,
    )

    assert isinstance(estimate, ParticleEstimate)
    assert math.isclose(estimate.velocity_kms, 11.0)
    assert math.isclose(estimate.mass_kg, 2.0)
    assert math.isclose(estimate.yield_c_per_kg, 1.0)
    assert estimate.velocity_details.source == "average"
    assert math.isclose(estimate.velocity_details.rise_time, 10.0)
    assert math.isclose(estimate.velocity_details.collection_efficiency, 12.0)


def test_compute_mass_from_charge_rejects_non_positive_charge():
    with pytest.raises(ValueError):
        compute_mass_from_charge(0.0, 1.0, make_flat_curve(1.0))
