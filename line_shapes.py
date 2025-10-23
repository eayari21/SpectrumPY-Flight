"""Common analytical line-shape definitions used by the dust composition UI."""

from __future__ import annotations

from math import pi, sqrt
from typing import Iterable, List
from scipy.special import erfc as comperror

import numpy as np

try:  # pragma: no cover - NumPy 2 compatibility
    erfc = getattr(np, "erfc")
except AttributeError:  # pragma: no cover - fallback for very old NumPy
    erfc = comperror

try:  # pragma: no cover - SciPy is optional at runtime
    from scipy.special import wofz as _wofz  # type: ignore

    _HAVE_SCIPY = True
except Exception:  # pragma: no cover - gracefully degrade without SciPy
    _HAVE_SCIPY = False


def gaussian(x: np.ndarray | Iterable[float], mu: float, sigma: float, *, area: float = 1.0) -> np.ndarray:
    """Return a Gaussian centred at ``mu`` with standard deviation ``sigma``."""

    arr = np.asarray(x, dtype=float)
    safe_sigma = sigma if abs(sigma) > 1.0e-18 else 1.0e-18
    prefactor = area / (safe_sigma * sqrt(2.0 * pi))
    exponent = -0.5 * ((arr - mu) / safe_sigma) ** 2
    return prefactor * np.exp(exponent)


def lorentzian(x: np.ndarray | Iterable[float], mu: float, gamma: float, *, area: float = 1.0) -> np.ndarray:
    """Return a Lorentzian line shape with half-width ``gamma``."""

    arr = np.asarray(x, dtype=float)
    safe_gamma = abs(gamma) if abs(gamma) > 1.0e-18 else 1.0e-18
    return area * (safe_gamma / pi) / ((arr - mu) ** 2 + safe_gamma**2)


def voigt(
    x: np.ndarray | Iterable[float],
    mu: float,
    sigma: float,
    gamma: float,
    *,
    area: float = 1.0,
) -> np.ndarray:
    """Return a Voigt profile using SciPy when available, or a pseudo-Voigt approximation."""

    arr = np.asarray(x, dtype=float)
    safe_sigma = abs(sigma) if abs(sigma) > 1.0e-18 else 1.0e-18
    safe_gamma = abs(gamma) if abs(gamma) > 1.0e-18 else 1.0e-18
    if _HAVE_SCIPY:
        z = ((arr - mu) + 1j * safe_gamma) / (safe_sigma * sqrt(2.0))
        return (area / (safe_sigma * sqrt(2.0 * pi))) * np.real(_wofz(z))

    # Pseudo-Voigt approximation (Kielkopf-like)
    f_g = 2.0 * sqrt(2.0 * np.log(2.0)) * safe_sigma
    f_l = 2.0 * safe_gamma
    f = (
        f_l**5
        + 2.69269 * f_l**4 * f_g
        + 2.42843 * f_l**3 * f_g**2
        + 4.47163 * f_l**2 * f_g**3
        + 0.07842 * f_l * f_g**4
        + f_g**5
    ) ** (1.0 / 5.0)
    eta = 1.36603 * (f_l / f) - 0.47719 * (f_l / f) ** 2 + 0.11116 * (f_l / f) ** 3
    eta = np.clip(eta, 0.0, 1.0)
    gamma_eff = f / 2.0
    sigma_eff = f / (2.0 * sqrt(2.0 * np.log(2.0)))
    return eta * lorentzian(arr, mu, gamma_eff, area=area) + (1.0 - eta) * gaussian(arr, mu, sigma_eff, area=area)


def emg(
    x: np.ndarray | Iterable[float],
    mu: float,
    sigma: float,
    tau: float,
    *,
    area: float = 1.0,
) -> np.ndarray:
    """Return an exponentially modified Gaussian (EMG). Positive ``tau`` gives a right tail."""

    arr = np.asarray(x, dtype=float)
    safe_sigma = abs(sigma) if abs(sigma) > 1.0e-18 else 1.0e-18
    if tau == 0.0:
        return gaussian(arr, mu, safe_sigma, area=area)
    lam = 1.0 / tau
    xeff = arr if tau > 0.0 else (2.0 * mu - arr)
    arg = (mu + lam * safe_sigma**2 - xeff) / (sqrt(2.0) * safe_sigma)
    pref = (lam * area) / 2.0
    expo = np.exp(0.5 * lam * (2.0 * mu + lam * safe_sigma**2 - 2.0 * xeff))
    return pref * expo * erfc(arg)


def double_emg(
    x: np.ndarray | Iterable[float],
    mu: float,
    sigma: float,
    tau1: float,
    tau2: float,
    *,
    w1: float = 0.5,
    area: float = 1.0,
) -> np.ndarray:
    """Return a convex combination of two EMGs sharing ``mu`` and ``sigma``."""

    arr = np.asarray(x, dtype=float)
    weight = float(np.clip(w1, 0.0, 1.0))
    return emg(arr, mu, sigma, tau1, area=area * weight) + emg(arr, mu, sigma, tau2, area=area * (1.0 - weight))


def hyper_emg(
    x: np.ndarray | Iterable[float],
    mu: float,
    sigma: float,
    *,
    taus_left: List[float] | None = None,
    taus_right: List[float] | None = None,
    weights_left: List[float] | None = None,
    weights_right: List[float] | None = None,
    area: float = 1.0,
) -> np.ndarray:
    """Return a HyperEMG (Gaussian core plus optional EMG tails)."""

    arr = np.asarray(x, dtype=float)
    taus_left = [] if taus_left is None else list(taus_left)
    taus_right = [] if taus_right is None else list(taus_right)

    wl = np.array([] if weights_left is None else weights_left, dtype=float)
    wr = np.array([] if weights_right is None else weights_right, dtype=float)

    wl = np.clip(wl, 0.0, None)
    wr = np.clip(wr, 0.0, None)

    total_weights = wl.sum() + wr.sum()
    if total_weights > 1.0:
        wl = wl / (total_weights + 1.0e-18)
        wr = wr / (total_weights + 1.0e-18)
        total_weights = wl.sum() + wr.sum()

    core_weight = max(0.0, 1.0 - total_weights)
    result = np.zeros_like(arr, dtype=float)
    if core_weight > 0.0:
        result += gaussian(arr, mu, sigma, area=area * core_weight)

    for tau, weight in zip(taus_left, wl):
        if weight <= 0.0:
            continue
        result += emg(arr, mu, sigma, -abs(tau), area=area * weight)

    for tau, weight in zip(taus_right, wr):
        if weight <= 0.0:
            continue
        result += emg(arr, mu, sigma, abs(tau), area=area * weight)

    return result


def generalized_normal(
    x: np.ndarray | Iterable[float],
    mu: float,
    alpha: float,
    beta: float,
    *,
    area: float = 1.0,
) -> np.ndarray:
    """Return the Subbotin (generalised normal) distribution."""

    from math import gamma as _gamma

    arr = np.asarray(x, dtype=float)
    safe_alpha = abs(alpha) if abs(alpha) > 1.0e-18 else 1.0e-18
    safe_beta = abs(beta) if abs(beta) > 1.0e-18 else 1.0e-18
    coef = safe_beta / (2.0 * safe_alpha * _gamma(1.0 / safe_beta))
    return area * coef * np.exp(-np.abs((arr - mu) / safe_alpha) ** safe_beta)


__all__ = [
    "gaussian",
    "lorentzian",
    "voigt",
    "emg",
    "double_emg",
    "hyper_emg",
    "generalized_normal",
]
