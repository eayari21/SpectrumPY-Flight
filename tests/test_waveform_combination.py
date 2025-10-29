import pytest

np = pytest.importorskip("numpy")

from dust_composition import detect_saturation


def test_detect_saturation_flags_clipped_segments_with_jitter():
    times = np.linspace(0.0, 1.0, 2048)
    base = 0.24 * np.exp(-0.5 * ((times - 0.5) / 0.02) ** 2)
    clip_level = 0.19
    signal = base.copy()
    clipped = base >= clip_level

    # Simulate ADC clipping that introduces jitter while hugging the clip level
    # rather than producing a perfectly flat plateau.
    jitter = 2.5e-4 * np.sin(np.linspace(0.0, 18.0 * np.pi, clipped.sum()))
    signal[clipped] = clip_level - 5.0e-4 + jitter

    mask = detect_saturation(signal, times)

    saturated_window = (times >= 0.48) & (times <= 0.52)
    assert saturated_window.any()
    assert np.all(mask[saturated_window])

    unsaturated_window = (times <= 0.45) | (times >= 0.55)
    assert unsaturated_window.any()
    assert not np.any(mask[unsaturated_window])
