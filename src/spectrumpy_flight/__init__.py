"""SpectrumPY-Flight package metadata and helpers."""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["__version__", "package_path", "default_hdf5_dir", "hdf5_search_paths"]

__version__ = "1.1.0"

_PACKAGE_ROOT = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# HDF5 file locking
# ---------------------------------------------------------------------------
#
# Many of our regression tests open decoded olivine HDF5 files that live on
# shared or network filesystems.  The default HDF5 behaviour is to lock files
# when they are opened, which can raise ``BlockingIOError`` on filesystems that
# do not fully support POSIX advisory locks.  ``HDF5_USE_FILE_LOCKING`` can be
# set to ``FALSE`` to disable the locking logic, which is safe for our use case
# because the datasets are read-only during tests.  Set the default here so all
# consumers of the package automatically benefit without having to remember to
# tweak their environment.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def package_path(*parts: str) -> Path:
    """Return the absolute path to a resource within the package."""

    return _PACKAGE_ROOT.joinpath(*parts)


from .paths import default_hdf5_dir, hdf5_search_paths  # noqa: E402
