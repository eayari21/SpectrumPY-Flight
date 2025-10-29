import pytest

np = pytest.importorskip("numpy")

from dust_composition import combine_waveform_channels, detect_saturation, GAIN_HIGH, GAIN_LOW, GAIN_MEDIUM


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


def test_combine_waveform_channels_replaces_saturated_high_with_medium():
    times = np.linspace(0.0, 1.0, 1024)
    physical_signal = 0.5 * np.exp(-0.5 * ((times - 0.5) / 0.05) ** 2)

    high_baseline = 1200.0
    medium_baseline = 75.0

    high = high_baseline + physical_signal * GAIN_HIGH
    medium = medium_baseline + physical_signal * GAIN_MEDIUM

    saturation_mask = (times >= 0.45) & (times <= 0.55)
    high[saturation_mask] = high_baseline + 0.75 * GAIN_HIGH

    combined = combine_waveform_channels(times, high, medium, None)
    assert combined is not None

    expected = high.copy()
    expected[saturation_mask] = high_baseline + (medium[saturation_mask] - medium_baseline) * (
        GAIN_HIGH / GAIN_MEDIUM
    )

    assert np.allclose(combined, expected, rtol=1e-6, atol=1e-6)


def test_combine_waveform_channels_uses_low_when_high_and_medium_saturate():
    times = np.linspace(0.0, 1.0, 512)
    physical_signal = 0.25 + 0.45 * np.exp(-((times - 0.52) / 0.04) ** 2)

    high_baseline = 980.0
    medium_baseline = 63.0
    low_baseline = -4.5

    high = high_baseline + physical_signal * GAIN_HIGH
    medium = medium_baseline + physical_signal * GAIN_MEDIUM
    low = low_baseline + physical_signal * GAIN_LOW

    clip_mask = (times >= 0.48) & (times <= 0.54)
    high[clip_mask] = high_baseline + 0.6 * GAIN_HIGH
    medium[clip_mask] = medium_baseline + 0.6 * GAIN_MEDIUM

    combined = combine_waveform_channels(times, high, medium, low)
    assert combined is not None

    expected = high.copy()
    expected[clip_mask] = high_baseline + (low[clip_mask] - low_baseline) * (
        GAIN_HIGH / GAIN_LOW
    )

    assert np.allclose(combined, expected, rtol=1e-6, atol=1e-6)
