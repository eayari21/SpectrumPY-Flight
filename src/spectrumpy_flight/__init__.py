"""SpectrumPY-Flight package metadata and helpers."""

from __future__ import annotations

from pathlib import Path

__all__ = ["__version__", "package_path", "default_hdf5_dir", "hdf5_search_paths"]

__version__ = "1.1.0"

_PACKAGE_ROOT = Path(__file__).resolve().parent


def package_path(*parts: str) -> Path:
    """Return the absolute path to a resource within the package."""

    return _PACKAGE_ROOT.joinpath(*parts)


from .paths import default_hdf5_dir, hdf5_search_paths  # noqa: E402
