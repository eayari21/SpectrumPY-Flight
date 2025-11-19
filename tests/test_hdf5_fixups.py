from __future__ import annotations

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from spectrumpy_flight.hdf5_fixups import ensure_legacy_analysis_compatibility


def _write_legacy_file(path):
    with h5py.File(path, "w") as handle:
        analysis = handle.require_group("1/Analysis")
        analysis.create_dataset("Ion Grid Fit Parameters", data=np.arange(3.0))
        # Leave Target H without a legacy dataset to ensure the fix-up creates an
        # empty placeholder.


def test_fixups_create_canonical_fit_params(tmp_path):
    decoded_path = tmp_path / "ois_output_00000000_000000.h5"
    _write_legacy_file(decoded_path)

    ensure_legacy_analysis_compatibility([tmp_path])

    with h5py.File(decoded_path, "r") as handle:
        analysis = handle["1/Analysis"]
        canonical = analysis["Ion GridFitParams"][()]
        assert canonical.tolist() == [0.0, 1.0, 2.0]

        # The legacy dataset path should now be a soft link to the canonical
        # dataset instead of a standalone dataset.
        assert isinstance(analysis.get("Ion Grid Fit Parameters", getlink=True), h5py.SoftLink)

        # Target H never had a fit parameter dataset, but the fix-up should add
        # an empty placeholder so downstream tooling can rely on the path.
        placeholder = analysis["Target HFitParams"][()]
        assert placeholder.size == 0


def test_fixups_backfill_flag_datasets(tmp_path):
    decoded_path = tmp_path / "ois_output_flags_missing.h5"
    with h5py.File(decoded_path, "w") as handle:
        handle.require_group("1/Analysis")

    ensure_legacy_analysis_compatibility([tmp_path])

    with h5py.File(decoded_path, "r") as handle:
        flags = handle["1/Analysis/Flags"]
        for name in ("FailedFits", "SaturatedChannels", "Notes"):
            assert name in flags
            dataset = flags[name]
            assert dataset.shape == (0,)
            assert dataset.dtype.kind in {"O", "S"}


def test_fixups_backfill_impact_charge_datasets(tmp_path):
    decoded_path = tmp_path / "ois_output_impact_missing.h5"
    with h5py.File(decoded_path, "w") as handle:
        analysis = handle.require_group("1/Analysis")
        analysis.create_dataset("Ion GridImpactCharge", data=np.array([1.23], dtype=float))

    ensure_legacy_analysis_compatibility([tmp_path])

    with h5py.File(decoded_path, "r") as handle:
        analysis = handle["1/Analysis"]

        canonical = analysis["Ion Grid Impact Charge"][()]
        assert canonical.tolist() == [pytest.approx(1.23)]
        assert isinstance(
            analysis.get("Ion GridImpactCharge", getlink=True), h5py.SoftLink
        )

        placeholder = analysis["Target H Impact Charge"][()]
        assert placeholder.shape == (1,)
        assert np.isnan(placeholder[0])
