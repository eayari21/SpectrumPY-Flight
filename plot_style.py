"""Project-wide Matplotlib styling helpers.

The goal of this module is to keep figures visually consistent across the
SpectrumPY tools regardless of which entry point is used.  The helper below
configures fonts, tick visibility, and background colours so that exported
figures match the high-contrast examples provided by science operations.
"""

from __future__ import annotations

import os
from functools import lru_cache
from typing import Mapping

from matplotlib import cycler, pyplot as plt

__all__ = ["apply_plot_style", "available_color_schemes"]

_BASE_RC: Mapping[str, object] = {
    "font.family": "DejaVu Sans",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "Nimbus Sans", "Liberation Sans"],
    "font.size": 13,
    "axes.labelsize": 15,
    "axes.labelweight": "semibold",
    "axes.titlesize": 17,
    "axes.titleweight": "semibold",
    "axes.linewidth": 1.4,
    "axes.labelpad": 6.0,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 6.5,
    "xtick.major.width": 1.25,
    "ytick.major.size": 6.5,
    "ytick.major.width": 1.25,
    "xtick.minor.size": 3.5,
    "xtick.minor.width": 1.0,
    "ytick.minor.size": 3.5,
    "ytick.minor.width": 1.0,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "legend.frameon": False,
    "legend.fontsize": 12,
    "figure.dpi": 110,
    "figure.autolayout": False,
    "savefig.bbox": "tight",
}

_LIGHT_RC: Mapping[str, object] = {
    "figure.facecolor": "#ffffff",
    "axes.facecolor": "#ffffff",
    "axes.edgecolor": "#121212",
    "axes.labelcolor": "#121212",
    "axes.prop_cycle": cycler(color=[
        "#0c5da5",  # blue
        "#ff8100",  # orange
        "#7a68a6",  # purple
        "#8fb339",  # green
        "#c62828",  # red
        "#0096c7",  # cyan
    ]),
    "xtick.color": "#121212",
    "ytick.color": "#121212",
    "grid.color": "#d4d4d4",
    "grid.linestyle": "-",
    "grid.linewidth": 0.7,
    "text.color": "#121212",
    "savefig.facecolor": "#ffffff",
    "savefig.edgecolor": "#ffffff",
}

_DARK_RC: Mapping[str, object] = {
    "figure.facecolor": "#0b0d16",
    "axes.facecolor": "#0b0d16",
    "axes.edgecolor": "#f5f5f5",
    "axes.labelcolor": "#f5f5f5",
    "axes.prop_cycle": cycler(color=[
        "#ffbf3d",  # golden yellow
        "#ff6f91",  # pink
        "#7bf1a8",  # mint
        "#70a0ff",  # light blue
        "#f2a541",  # orange
        "#d0d0ff",  # lavender
    ]),
    "xtick.color": "#f5f5f5",
    "ytick.color": "#f5f5f5",
    "grid.color": "#2a2d3a",
    "grid.linestyle": "-",
    "grid.linewidth": 0.7,
    "text.color": "#f5f5f5",
    "savefig.facecolor": "#0b0d16",
    "savefig.edgecolor": "#0b0d16",
}

_COLOR_SCHEMES: Mapping[str, Mapping[str, object]] = {
    "light": _LIGHT_RC,
    "dark": _DARK_RC,
}


def available_color_schemes() -> tuple[str, ...]:
    """Return the known colour schemes in alphabetical order."""

    return tuple(sorted(_COLOR_SCHEMES))


@lru_cache(maxsize=None)
def _resolve_color_scheme(name: str | None) -> Mapping[str, object]:
    if not name:
        return _COLOR_SCHEMES["light"]

    normalized = name.strip().lower()
    if normalized not in _COLOR_SCHEMES:
        raise ValueError(
            f"Unknown plot color scheme '{name}'. Available options: {', '.join(available_color_schemes())}."
        )
    return _COLOR_SCHEMES[normalized]


def apply_plot_style(color_scheme: str | None = None) -> None:
    """Apply the project Matplotlib style.

    Parameters
    ----------
    color_scheme:
        Optional colour scheme to apply (``"light"`` or ``"dark"``).  When not
        provided, the function checks the ``IDEX_PLOT_STYLE`` environment
        variable; if it is unset the light theme is used.
    """

    requested_scheme = color_scheme or os.environ.get("IDEX_PLOT_STYLE", "light")

    plt.style.use("default")
    plt.rcParams.update(_BASE_RC)
    plt.rcParams.update(_resolve_color_scheme(requested_scheme))

    # Tighten legend spacing and tick label padding to match the provided
    # reference figures.  These values are set after the general rcParams update
    # so they are applied regardless of downstream customisation.
    plt.rcParams["legend.borderaxespad"] = 0.8
    plt.rcParams["xtick.major.pad"] = 6
    plt.rcParams["ytick.major.pad"] = 6

