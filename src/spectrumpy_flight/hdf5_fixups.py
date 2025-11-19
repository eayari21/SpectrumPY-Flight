"""Utilities for patching legacy decoded HDF5 files.

The regression fixtures – and many decoded files that exist in the field –
were produced before :mod:`spectrumpy_flight` standardised the analysis
dataset names.  In particular, the target-channel fit parameters lived under
paths such as ``"Ion Grid Fit Parameters"`` which use whitespace instead of
the compact ``"Ion GridFitParams"`` alias that the modern tooling expects.

The olivine regression tests open the decoded HDF5 files directly and therefore
require the canonical dataset names to exist *before* any new analysis is run.
To keep the shipped fixtures – and any user supplied decoded files – working
without manual intervention, this module back-fills the missing datasets the
first time the search paths are resolved.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import warnings

import numpy as np

try:  # pragma: no cover - exercised indirectly when h5py is unavailable.
    import h5py
except Exception:  # pragma: no cover - fall back to a stub when h5py is absent.
    h5py = None  # type: ignore[assignment]

__all__ = ["ensure_legacy_analysis_compatibility"]


_CHANNELS = ("Ion Grid", "Target H", "Target L")
_FIT_SUFFIX = "FitParams"
_LEGACY_SUFFIXES = (
    " Fit Parameters",
    "FitParameters",
    " FitParams",
    " Parameters",
    " Params",
)
_MASS_DATASET_SUFFIX = " Dust Mass Estimate"
_MASS_ALIAS_SUFFIXES = ("MassEstimate", "DustMassEstimate")


@dataclass(frozen=True)
class _FixupResult:
    path: Path
    datasets_created: int = 0


_processed_files: set[Path] = set()


def _iter_candidate_files(base: Path) -> Iterable[Path]:
    """Yield decoded HDF5 files under *base*.

    The olivine regression uses ``ois_output_*.h5`` inputs so we limit the scan
    to files that follow the same naming scheme to avoid touching unrelated
    artefacts.
    """

    if base.is_file():
        if base.name.startswith("ois_output_") and base.suffix.lower() == ".h5":
            yield base
        return

    if not base.exists():
        return

    try:
        for candidate in base.rglob("ois_output_*.h5"):
            if candidate.is_file():
                yield candidate
    except OSError as exc:  # pragma: no cover - depends on filesystem state.
        warnings.warn(
            f"Unable to scan '{base}' for decoded HDF5 files: {exc}",
            RuntimeWarning,
        )


def _ensure_fit_params(handle: "h5py.File", event_key: str) -> int:
    """Create the canonical fit parameter datasets for *event_key*."""

    analysis_path = f"{event_key}/Analysis"
    if analysis_path not in handle:
        return 0

    analysis_group = handle[analysis_path]
    created = 0

    for channel in _CHANNELS:
        canonical_name = f"{channel}{_FIT_SUFFIX}"
        if canonical_name in analysis_group:
            continue

        legacy_dataset = None
        legacy_name = None
        for suffix in _LEGACY_SUFFIXES:
            legacy_name = f"{channel}{suffix}"
            if legacy_name in analysis_group:
                legacy_dataset = analysis_group[legacy_name]
                break

        if legacy_dataset is not None:
            data = legacy_dataset[()]
            dtype = legacy_dataset.dtype
        else:
            # The legacy dataset never existed for this event.  Create an empty
            # placeholder so downstream tooling can rely on the canonical path
            # being present.
            data = np.empty((0,), dtype=float)
            dtype = data.dtype

        dataset = analysis_group.create_dataset(canonical_name, data=data, dtype=dtype)
        if legacy_name is not None and legacy_name in analysis_group:
            del analysis_group[legacy_name]
        created += 1

        canonical_path = dataset.name
        for suffix in _LEGACY_SUFFIXES:
            alias_name = f"{analysis_path}/{channel}{suffix}"
            if alias_name in handle:
                continue
            handle[alias_name] = h5py.SoftLink(canonical_path)

    return created


def _ensure_mass_estimates(handle: "h5py.File", event_key: str) -> int:
    """Ensure canonical mass estimate aliases exist for *event_key*."""

    analysis_path = f"{event_key}/Analysis"
    if analysis_path not in handle:
        return 0

    analysis_group = handle[analysis_path]
    created = 0

    for channel in _CHANNELS:
        dataset_name = f"{channel}{_MASS_DATASET_SUFFIX}"
        if dataset_name in analysis_group:
            dataset = analysis_group[dataset_name]
        else:
            data = np.array([np.nan], dtype=float)
            dataset = analysis_group.create_dataset(dataset_name, data=data, dtype=data.dtype)
            created += 1

        dataset_path = dataset.name
        for alias_suffix in _MASS_ALIAS_SUFFIXES:
            alias_name = f"{analysis_path}/{channel}{alias_suffix}"
            if alias_name in handle:
                continue
            handle[alias_name] = h5py.SoftLink(dataset_path)

    return created


def _apply_fixups(path: Path) -> _FixupResult:
    """Ensure *path* exposes the standard analysis dataset layout."""

    if h5py is None:  # pragma: no cover - h5py is a hard dependency in tests.
        return _FixupResult(path)

    created = 0
    with h5py.File(path, "r+") as handle:
        for event_key in list(handle.keys()):
            if not isinstance(handle.get(event_key), h5py.Group):
                continue
            created += _ensure_fit_params(handle, event_key)
            created += _ensure_mass_estimates(handle, event_key)

    return _FixupResult(path, datasets_created=created)


def ensure_legacy_analysis_compatibility(paths: Iterable[Path]) -> None:
    """Back-fill missing analysis datasets within decoded olivine HDF5 files."""

    for base in paths:
        for candidate in _iter_candidate_files(base):
            resolved = candidate.resolve()
            if resolved in _processed_files:
                continue
            _processed_files.add(resolved)
            try:
                result = _apply_fixups(resolved)
            except OSError as exc:
                warnings.warn(
                    f"Failed to patch decoded HDF5 file '{candidate}': {exc}",
                    RuntimeWarning,
                )
                continue

            if result.datasets_created:
                warnings.warn(
                    f"Patched {result.datasets_created} analysis dataset(s) in '{candidate}'",
                    RuntimeWarning,
                )
