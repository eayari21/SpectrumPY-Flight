import json

import numpy as np

from spectrumpy_flight.idex_packet import _normalise_hdf5_data_for_cdf


def test_normalise_hdf5_data_for_cdf_serializes_structured_array() -> None:
    dtype = np.dtype([
        ("id", "<i4"),
        ("label", "S16"),
        ("assigned_mass", "<f8"),
    ])
    values = np.array([(1, b"Fe", 55.845), (2, b"Mg", 24.305)], dtype=dtype)

    normalized, data_type, elements = _normalise_hdf5_data_for_cdf(values)

    assert data_type == 51
    assert elements >= 1
    assert normalized.shape == (1,)

    payload = json.loads(str(normalized[0]))
    assert payload[0]["id"] == 1
    assert payload[0]["label"] == "Fe"
    assert payload[1]["assigned_mass"] == 24.305
