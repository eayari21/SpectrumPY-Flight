"""Shared filesystem helpers for SpectrumPY-Flight."""

from __future__ import annotations

import os
from pathlib import Path

from . import package_path

__all__ = ["DEFAULT_HDF5_ENV", "default_hdf5_dir", "hdf5_search_paths"]

# Environment variable that overrides the default HDF5 directory.
DEFAULT_HDF5_ENV = "SPECTRUMPY_HDF5_DIR"


def _cwd_hdf5_dir() -> Path:
    """Return the HDF5 directory relative to the current working tree."""

    return Path.cwd() / "HDF5"


def default_hdf5_dir(create: bool = False) -> Path:
    """Return the preferred HDF5 directory for generated products.

    The lookup order is:

    1. :data:`DEFAULT_HDF5_ENV` when set.  The directory is expanded with
       ``Path.expanduser`` so ``"~/data"`` style values work.
    2. ``./HDF5`` relative to the current working directory.  This mirrors the
       historical layout used by the shell helpers and existing pipelines.
    3. The package's bundled ``HDF5`` directory.  This is primarily useful when
       reading shipped resources (for example, within the test suite) rather
       than writing new analysis products.

    Parameters
    ----------
    create:
        When ``True``, ensure the resolved directory exists.
    """

    env_value = os.environ.get(DEFAULT_HDF5_ENV)
    if env_value:
        target = Path(env_value).expanduser()
        if create:
            target.mkdir(parents=True, exist_ok=True)
        return target

    cwd_target = _cwd_hdf5_dir()
    if cwd_target.exists() or create:
        if create:
            cwd_target.mkdir(parents=True, exist_ok=True)
        return cwd_target

    package_target = package_path("HDF5")
    if package_target.exists():
        if create:
            package_target.mkdir(parents=True, exist_ok=True)
        return package_target

    # Fall back to the working tree location even if it does not currently
    # exist so callers have a predictable path to work with.
    if create:
        cwd_target.mkdir(parents=True, exist_ok=True)
    return cwd_target


def hdf5_search_paths() -> list[Path]:
    """Return candidate directories that may contain analysis HDF5 files."""

    paths: list[Path] = []

    env_value = os.environ.get(DEFAULT_HDF5_ENV)
    if env_value:
        paths.append(Path(env_value).expanduser())

    paths.append(_cwd_hdf5_dir())
    paths.append(package_path("HDF5"))

    # Remove duplicates while preserving order.
    unique_paths: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        unique_paths.append(path)

    return unique_paths
