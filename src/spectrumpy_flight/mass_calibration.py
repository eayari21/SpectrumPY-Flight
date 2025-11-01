"""Utilities for fitting and inverting time-of-flight mass calibration series."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np


def fit_tof_series(
    known_masses_u: Sequence[float],
    measured_times_us: Sequence[float],
    max_order: Optional[int] = None,
    weights: Optional[Sequence[float]] = None,
) -> np.ndarray:
    """Fit a polynomial series in ``s = sqrt(m)`` to TOF calibration data.

    Parameters
    ----------
    known_masses_u:
        Reference masses in atomic mass units.
    measured_times_us:
        Measured time-of-flight values (typically microseconds).
    max_order:
        Highest power of ``s`` to include.  If ``None`` the order is ``N-1``
        where ``N`` is the number of calibration lines.
    weights:
        Optional weights for each calibration line.  The weights are applied
        as ``sqrt(weights)`` inside the least-squares solve so callers can pass
        ``1/sigma**2`` style inputs.

    Returns
    -------
    numpy.ndarray
        Coefficient vector ``[t0, a1, ..., aK]`` of length ``K+1``.
    """

    masses = np.asarray(known_masses_u, dtype=float)
    times = np.asarray(measured_times_us, dtype=float)
    if masses.size == 0 or masses.size != times.size:
        raise ValueError("Masses and times must have the same non-zero length")

    order = masses.size - 1 if max_order is None else min(max_order, masses.size - 1)
    if order < 1:
        raise ValueError("At least two calibration lines are required to fit the series")

    s = np.sqrt(np.clip(masses, 0.0, None))
    design = np.vstack([np.power(s, k) for k in range(order + 1)]).T

    if weights is not None:
        w = np.sqrt(np.asarray(weights, dtype=float))
        if w.shape != masses.shape:
            raise ValueError("Weights must match the shape of the calibration masses")
        design = design * w[:, None]
        times = times * w

    coeffs, *_ = np.linalg.lstsq(design, times, rcond=None)
    return coeffs


def _evaluate_series(coeffs: np.ndarray, s: np.ndarray) -> np.ndarray:
    result = np.zeros_like(s, dtype=float) + coeffs[0]
    for k in range(1, coeffs.size):
        result += coeffs[k] * np.power(s, k)
    return result


def _evaluate_series_derivative(coeffs: np.ndarray, s: np.ndarray) -> np.ndarray:
    derivative = np.zeros_like(s, dtype=float)
    for k in range(1, coeffs.size):
        derivative += k * coeffs[k] * np.power(s, k - 1)
    return derivative


def invert_tof_to_mass(
    tof_us: Sequence[float],
    coeffs: Sequence[float],
    *,
    m_min_u: float = 0.0,
    m_max_u: float = 400.0,
    newton_iters: int = 6,
    tol: float = 1e-9,
) -> np.ndarray:
    """Invert ``t = f(s)`` to obtain masses for arbitrary TOF samples."""

    t = np.asarray(tof_us, dtype=float)
    c = np.asarray(coeffs, dtype=float)
    if c.size < 2:
        raise ValueError("At least two coefficients are required for inversion")

    t0 = float(c[0])
    a1 = float(c[1]) if c.size > 1 else 0.0
    if not np.isfinite(a1) or abs(a1) < 1e-12:
        raise ValueError("First-order coefficient must be finite and non-zero")

    s = np.maximum((t - t0) / a1, 0.0)
    s_max = float(np.sqrt(max(m_max_u, 0.0)))
    if not np.isfinite(s_max) or s_max <= 0.0:
        s_max = 1.0

    for _ in range(max(1, newton_iters)):
        f = _evaluate_series(c, s) - t
        df = _evaluate_series_derivative(c, s)
        with np.errstate(divide="ignore", invalid="ignore"):
            step = np.where(np.abs(df) > 1e-18, f / df, 0.0)
        s_new = s - step
        # Ensure non-negative, finite, and within bounds
        s_new = np.where(np.isfinite(s_new) & (s_new >= 0.0), s_new, 0.5 * s)
        s_new = np.clip(s_new, 0.0, s_max)
        if np.all(np.abs(s_new - s) < tol):
            s = s_new
            break
        s = s_new

    masses = np.square(s)
    masses = np.clip(masses, max(0.0, m_min_u), max(m_min_u, m_max_u))
    return masses


@dataclass
class TOFMassCal:
    """Encapsulate a TOF mass calibration series."""

    coeffs: np.ndarray
    mass_range_u: Tuple[float, float] = (0.0, 400.0)

    def __post_init__(self) -> None:
        coeffs = np.asarray(self.coeffs, dtype=float)
        if coeffs.size < 2:
            raise ValueError("Calibration requires at least t0 and a1 coefficients")
        object.__setattr__(self, "coeffs", coeffs)
        m_min, m_max = self.mass_range_u
        m_min = float(np.clip(m_min, 0.0, None))
        m_max = float(max(m_min, m_max))
        object.__setattr__(self, "mass_range_u", (m_min, m_max))

    @property
    def order(self) -> int:
        return self.coeffs.size - 1

    def mass_to_tof(self, masses_u: Sequence[float]) -> np.ndarray:
        masses = np.asarray(masses_u, dtype=float)
        s = np.sqrt(np.clip(masses, 0.0, None))
        return _evaluate_series(self.coeffs, s)

    def tof_to_mass(self, tof_us: Sequence[float], **kwargs) -> np.ndarray:
        return invert_tof_to_mass(
            tof_us,
            self.coeffs,
            m_min_u=self.mass_range_u[0],
            m_max_u=self.mass_range_u[1],
            **kwargs,
        )

    def derivative(self, masses_u: Sequence[float]) -> np.ndarray:
        masses = np.asarray(masses_u, dtype=float)
        s = np.sqrt(np.clip(masses, 0.0, None))
        return _evaluate_series_derivative(self.coeffs, s)

    def is_monotonic(self, samples: int = 256) -> bool:
        m_min, m_max = self.mass_range_u
        s_min = float(np.sqrt(max(m_min, 0.0)))
        s_max = float(np.sqrt(max(m_max, 0.0)))
        if not np.isfinite(s_max) or s_max <= s_min:
            s_max = s_min + 1.0
        s_vals = np.linspace(s_min, s_max, max(samples, 16))
        derivative = _evaluate_series_derivative(self.coeffs, s_vals)
        return np.all(derivative > 0.0)

    def to_dict(self) -> dict:
        return {
            "series_type": "sqrtm_poly",
            "coeffs": self.coeffs.tolist(),
            "mass_range_u": [float(self.mass_range_u[0]), float(self.mass_range_u[1])],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "TOFMassCal":
        coeffs = data.get("coeffs", [])
        mass_range = tuple(data.get("mass_range_u", (0.0, 400.0)))
        return cls(np.asarray(coeffs, dtype=float), mass_range_u=mass_range)

    @classmethod
    def from_lines(
        cls,
        masses_u: Sequence[float],
        times_us: Sequence[float],
        *,
        max_order: Optional[int] = None,
        weights: Optional[Sequence[float]] = None,
        mass_range: Tuple[float, float] = (0.0, 400.0),
        enforce_monotonic: bool = True,
    ) -> Optional["TOFMassCal"]:
        masses = np.asarray(masses_u, dtype=float)
        times = np.asarray(times_us, dtype=float)
        if masses.size < 2 or masses.size != times.size:
            return None

        limit_order = masses.size - 1 if max_order is None else min(max_order, masses.size - 1)
        for order in range(limit_order, 0, -1):
            coeffs = fit_tof_series(masses, times, max_order=order, weights=weights)
            cal = cls(coeffs, mass_range_u=mass_range)
            if not enforce_monotonic or cal.is_monotonic():
                return cal
        return None

    def evaluate(self, masses_u: Sequence[float]) -> np.ndarray:
        return self.mass_to_tof(masses_u)
