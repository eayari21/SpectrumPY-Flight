# Time-to-Mass Calibration Theory and Practice

Time-of-flight (TOF) mass spectrometers record ion arrival times that must be mapped to mass-to-charge ratios before the data can be interpreted.  The `spectrumpy_flight` package provides utilities in [`mass_calibration.py`](../mass_calibration.py) and the [`time2mass.py`](../time2mass.py) application for estimating calibration polynomials, performing robust inversions, and applying the resulting time-to-mass transforms to spectra and quicklook products.  This document summarizes the theory, algorithms, and recommended workflows for mass calibration within the project.

## Physical model

Under the standard constant-field TOF approximation the travel time for an ion with mass-to-charge ratio $m / q$ is well approximated by

$$ t(m) = t_0 + \sum_{k=1}^{K} a_k \left(\sqrt{m}\right)^k, $$

where $t_0$ is the electronic delay, $a_1$ scales with the inverse extraction field strength, and the higher-order coefficients capture detector timing offsets and energy spreads.  Introducing the reduced coordinate $s = \sqrt{m}$ yields a polynomial series in $s$

$$ t(s) = t_0 + a_1 s + a_2 s^2 + \cdots + a_K s^K. $$

The coefficients $a_k$ are determined from calibration data by solving a weighted least-squares system.  Given $N$ calibration lines ($N \ge 2$), the design matrix $\mathbf{A} \in \mathbb{R}^{N \times (K+1)}$ has entries

$$ A_{i,k} = s_i^k, \qquad s_i = \sqrt{m_i}, $$

and the optimal coefficients solve the normal equations

$$ (\mathbf{A}^T \mathbf{W} \mathbf{A})\,\mathbf{c} = \mathbf{A}^T \mathbf{W} \mathbf{t}. $$

Here $\mathbf{W}$ is a diagonal weight matrix with elements $w_i = 1/\sigma_i^2$ when measurement uncertainties $\sigma_i$ are known.  The `fit_tof_series` routine directly computes the least-squares solution via `numpy.linalg.lstsq`, which internally performs an SVD for numerical stability.

## Inversion by Newton iteration

To convert a measured TOF sample back to a mass coordinate we must invert $t(s)$.  The `invert_tof_to_mass` helper and the [`TOFMassCal`](../mass_calibration.py) data class implement a safeguarded Newton iteration.  Starting from the affine approximation

$$ s^{(0)} = \max\left(\frac{t - t_0}{a_1}, 0\right), $$

the iteration refines the solution using

$$ s^{(n+1)} = \min\!\left(\max\!\left(s^{(n)} - \frac{f\big(s^{(n)}\big)}{f'\big(s^{(n)}\big)},\; 0\right),\; s_{\max}\right) $$

with

$$ f(s) = t(s) - t, \qquad f'(s) = \frac{\mathrm{d} t(s)}{\mathrm{d}s} = \sum_{k=1}^{K} k\,a_k\,s^{k-1}. $$

The iteration terminates when $|s^{(n+1)} - s^{(n)}| < \varepsilon$ (default $10^{-9}$) or when a fixed number of iterations is reached.  After convergence the squared reduced coordinate gives the calibrated mass,

$$ m = \left(s^{(n)}\right)^2, $$

which is finally clipped to the instrument-supported range $[m_{\min}, m_{\max}]$.

## Stretch-and-shift initialization

The interactive `time2mass.py` workflow provides an intuitive stretch-and-shift step to initialize the polynomial fit before loading higher-order corrections.  The algorithm computes a "time-zero" axis and performs the following operations:

1. **Stretch search.**  Estimate the global scale factor $\alpha$ by scanning $\alpha \in [1.3, 1.6]\,\mu\text{s}$ (configurable via `MASS_STRETCH_MIN_US`/`MAX_US`).  For each candidate the software rescales the time axis as $s = (t - \Delta)/\alpha$ and evaluates the agreement with reference lines.
2. **Shift detection.**  Apply a median-based baseline removal to the detection signal, smooth with a moving average kernel, and locate peaks using SciPy's `find_peaks`.  The shift $\Delta$ maximizing the reference alignment is selected.
3. **Polynomial bootstrapping.**  Convert the stretch-and-shift estimate into synthetic calibration masses, invoke `fit_tof_series`, and optionally merge with previously saved coefficients.

This initialization yields a near-physical prior for the Newton inversion, improving robustness on faint or noisy spectra.

## Algorithmic recipes

### Fitting a calibration series

1. Collect at least two reference lines with known masses $m_i$ and measured TOF values $t_i$.
2. Choose the maximum polynomial order $K$.  Empirically, $K = 3$ or $4$ suffices for `spectrumpy_flight` datasets.
3. Invoke
   ```python
   from spectrumpy_flight.mass_calibration import fit_tof_series
   coeffs = fit_tof_series(known_masses_u=masses, measured_times_us=times, max_order=K)
   ```
4. Create a calibration object to apply the solution:
   ```python
   from spectrumpy_flight.mass_calibration import TOFMassCal
   cal = TOFMassCal(coeffs, mass_range_u=(0.0, 400.0))
   ```
5. Use `cal.mass_to_tof` or `cal.tof_to_mass` to transform arrays.

### Applying the Newton inversion to a spectrum

```python
from spectrumpy_flight.mass_calibration import TOFMassCal
from spectrumpy_flight.time2mass import _prepare_signal

cal = TOFMassCal(coeffs, mass_range_u=(1.0, 350.0))
corrected_signal = _prepare_signal(raw_counts)
mass_axis = cal.tof_to_mass(time_axis_us)
```

### Saving and reusing coefficients

The `time2mass.py` interface stores calibration metadata alongside spectra.  To reuse coefficients:

```python
from spectrumpy_flight.mass_calibration import TOFMassCal
import numpy as np

coeffs = np.loadtxt("calibration_coeffs.txt")
cal = TOFMassCal(coeffs, mass_range_u=(1.0, 350.0))
# Later, evaluate arrival times for candidate ions
candidate_masses = np.array([23.0, 24.0, 39.0])
predicted_tof = cal.mass_to_tof(candidate_masses)
```

## Quality control and diagnostics

* **Residual analysis.**  Plot $t_i - t(m_i)$ versus $\sqrt{m_i}$ to detect systematic deviations.  Random residuals with RMS $< 10$ ns indicate a reliable calibration.
* **Derivative monitoring.**  The derivative $f'(s)$ must remain strictly positive across the mass range.  Negative derivatives cause the inversion to diverge; if encountered, reduce $K$ or constrain the fit with additional references.
* **Clipping behavior.**  The implementation clips Newton steps to the instrument bounds to prevent runaway solutions.  If many samples saturate at $m_{\max}$ consider expanding the range or revisiting the initial stretch.

## Integration with quicklook pipelines

* The quicklook tools (`IDEX-quicklook.py`, `qd_quicklook.py`, and related scripts) call into `time2mass.py` for automatic calibration when mass-aligned overlays are requested.
* The [`mass_calibration.py`](../mass_calibration.py) utilities are designed to be dependency-light so that headless workflows—such as the `run_all.py` batch processor—can precompute mass axes.
* Output products annotate headers with the coefficients `t0`, `a1`, ..., `aK` to support reproducibility and offline verification.

## References

* Dawson, P. H. *Quadrupole Mass Spectrometry and Its Applications*. Elsevier, 1976.
* Wiley, W. C., and H. H. McLaren. "Time-of-flight mass spectrometer with improved resolution." *Review of Scientific Instruments* 26, no. 12 (1955): 1150-1157.
