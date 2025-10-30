"""SpectrumPY-Flight science and visualization toolkit."""

from __future__ import annotations

from importlib import metadata


try:
    __version__ = metadata.version("spectrumpy-flight")
except metadata.PackageNotFoundError:  # pragma: no cover - fallback during local dev
    __version__ = "0.0.0"

__all__ = ["__version__"]