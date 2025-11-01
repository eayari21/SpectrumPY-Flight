import numpy as np

from spectrumpy_flight.mass_calibration import TOFMassCal, fit_tof_series, invert_tof_to_mass


def test_fit_tof_series_linear_model():
    masses = np.array([16.0, 64.0])
    s = np.sqrt(masses)
    t0 = 0.75
    a1 = 1.5
    times = t0 + a1 * s

    coeffs = fit_tof_series(masses, times)

    assert coeffs.shape == (2,)
    assert np.isclose(coeffs[0], t0, rtol=1e-9, atol=1e-9)
    assert np.isclose(coeffs[1], a1, rtol=1e-9, atol=1e-9)


def test_tof_mass_calibration_quadratic():
    masses = np.array([16.0, 40.0, 100.0])
    s = np.sqrt(masses)
    coeffs_true = np.array([0.25, 1.45, 0.015])
    times = coeffs_true[0] + coeffs_true[1] * s + coeffs_true[2] * s**2

    calibration = TOFMassCal.from_lines(masses, times, max_order=3)
    assert calibration is not None
    assert calibration.order == 2

    recovered_masses = calibration.tof_to_mass(times)
    assert np.allclose(recovered_masses, masses, rtol=1e-9, atol=1e-9)


def test_invert_tof_to_mass_vectorised():
    coeffs = np.array([0.1, 1.4])
    times = np.linspace(0.1, 5.0, 50)
    masses = invert_tof_to_mass(times, coeffs)
    expected = np.maximum((times - coeffs[0]) / coeffs[1], 0.0) ** 2
    assert np.allclose(masses, expected)
