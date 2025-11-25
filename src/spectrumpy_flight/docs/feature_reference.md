# SpectrumPY-Flight Feature and Algorithm Reference

This reference surveys every bundled application, command-line script, and helper shipped with the SpectrumPY-Flight distribution.  It is designed to complement the tutorials by cataloguing inputs, outputs, mathematical models, and CLI options in one place.  When a tool exposes numerical processing steps, the relevant equations are written explicitly so analysts can reproduce the computations outside the GUI.

---

## Desktop entry points

### `spectrumpy_flight.start`

*Launch command*: `python -m spectrumpy_flight.start`

The launcher enumerates the packaged applications, documentation, and example data.  From the **Quicklook** tile it spawns `spectrum_launcher.py`, which in turn spins up a Qt event loop, injects the shared plot palette via `plot_style.apply_plot_style`, and forwards the selected target (`IDEX-quicklook.py`, `Scope-IDEX-quicklook.py`, or `qd_quicklook.py`).【F:src/spectrumpy_flight/start.py†L1-L116】【F:src/spectrumpy_flight/spectrum_launcher.py†L167-L289】  The launcher maintains a single-process lifecycle so that secondary windows inherit QApplication state (fonts, palettes, logging destinations).

### `IDEX-quicklook.py`

*Launch command*: `python -m spectrumpy_flight.IDEX-quicklook [-f DATAFILE] [--eventnumber N]`

The flagship Qt GUI for browsing decoded telemetry.  Major algorithmic components include:

* **Event loader** – Delegates to `paths.py` to resolve the active data directory, then uses `h5py`/`cdflib` readers.  When multiple events are present a `QListView` enumerates them; selecting an entry triggers asynchronous dataset hydration and view updates.【F:src/spectrumpy_flight/IDEX-quicklook.py†L328-L653】【F:src/spectrumpy_flight/paths.py†L1-L136】
* **Mass alignment** – Calls into `time2mass.py` to stretch-and-shift the TOF axis, then overlays polynomial mass calibrations (`TOFMassCal`) so the quicklook plots can display physical mass units.  The workflow follows the Newton inversion discussed in [Time-to-Mass Calibration Theory and Practice](time_to_mass_calibration.md).【F:src/spectrumpy_flight/IDEX-quicklook.py†L1770-L2035】【F:src/spectrumpy_flight/time2mass.py†L1-L240】
* **Noise workspace** – Computes Gaussian fits, FFT magnitudes, and autocorrelation coefficients.  Noise power spectra use the Welch estimator `P_xx(f) = (1/K) * sum_{k=1..K} |F{w_k[n] x[n]}|^2`, where Hann-windowed segments `w_k[n]` provide leakage control before averaging.【F:src/spectrumpy_flight/IDEX-quicklook.py†L2027-L2188】【F:src/spectrumpy_flight/noise_analysis.py†L260-L429】
* **Dust composition fitting** – Interfaces with `dust_composition.py` to model peaks using exponentially modified Gaussian (EMG) profiles of the form `EMG(t; mu, sigma, lambda) = (lambda/2) * exp((lambda/2) * (2*mu + lambda*sigma^2 - 2*t)) * erfc((mu + lambda*sigma^2 - t) / (sqrt(2) * sigma))`.  Optimisation relies on lmfit’s Levenberg–Marquardt backend with bounds derived from the currently highlighted window.【F:src/spectrumpy_flight/dust_composition.py†L1188-L1478】

### `Scope-IDEX-quicklook.py`

*Launch command*: `python -m spectrumpy_flight.Scope-IDEX-quicklook [--filename FILE] [--eventnumber N]`

A variant of the main quicklook tuned for oscilloscope captures (`.h5` converted from `.trc`).  It shares the same calibration and plotting stack as the flagship viewer but defaults to raw scope channel naming and emphasises trigger block navigation.【F:src/spectrumpy_flight/Scope-IDEX-quicklook.py†L1-L828】【F:src/spectrumpy_flight/Scope-IDEX-quicklook.py†L829-L858】

### `qd_quicklook.py`

*Launch command*: `python -m spectrumpy_flight.qd_quicklook --sourcefolder DIR --targetfolder LABEL`

Generates quick diagnostic PNGs and summary CSVs for charge detector (`QD`) campaigns.  The script extracts `.trc` waveforms, decimates them according to sampling metadata, and applies a cubic-spline interpolation before integrating charge packets.  Charge is computed as `Q = integral_{t0}^{t1} I(t) dt ≈ sum_{i=0}^{N-1} I_i * Δt`, with integration limits driven by automatic threshold crossing plus user overrides stored in the resulting metadata JSON.【F:src/spectrumpy_flight/qd_quicklook.py†L1-L903】

### `science_tool.py`

*Launch command*: `python -m spectrumpy_flight.science_tool`

An operator console for GS-OASIS packet pulls.  It constructs a generator over XTCE-defined packets, decodes Rice-compressed telemetry via `rice_decode.idex_rice_Decode`, and writes high-, mid-, and low-gain TOF arrays into a consolidated HDF5 structure.  Trigger delays are extracted bitwise from telemetry headers; for example the high-gain delay uses the mask `hgdelay = (TOFdelay >> 20) & 0b1111111111`, mirroring the firmware packing scheme.【F:src/spectrumpy_flight/science_tool.py†L1-L233】  The tool then launches interactive Matplotlib dashboards for each channel.

### `run_all.py`

*Launch command*: `python -m spectrumpy_flight.run_all`

Automates a full ingest by driving the packet decoder, quicklook exports, and documentation bundling.  The orchestrator inspects the `Data/` directory, calls `drive_idex_packet.py` to parallelise conversions, runs `IDEX-quicklook` in headless export mode, and archives the outputs.  Environment variables can override the default data paths; refer to `paths.py` for precedence rules.【F:src/spectrumpy_flight/run_all.py†L1-L213】【F:src/spectrumpy_flight/paths.py†L1-L136】

---

## Command-line utilities

### `idex_packet.py`

*Launch command*: `python -m spectrumpy_flight.idex_packet -f RAWFILE`

Required argument:

| Option | Description |
| --- | --- |
| `-f`, `--file`, `--filename` | Input oscilloscope or telemetry capture (extensionless or `.bin`). |

Although the CLI is intentionally minimal, the module exposes numerous toggles through top-level constants (GPU enablement, worker pool sizing, CDF export) so packaged wrappers can customise behaviour without destabilising the operator workflow.【F:src/spectrumpy_flight/idex_packet.py†L1-L1295】  Core algorithms include:

* **Rice/Golomb decompression** using a hybrid unary–binary decoder that enforces the instrument’s `k` parameter per channel.【F:src/spectrumpy_flight/rice_decode.py†L1-L274】
* **Time axis reconstruction** from trigger metadata, producing arrays `t_i = (i - N_pre) * Δt + t_offset`, where `N_pre` and `Δt` are derived from the XTCE header words.【F:src/spectrumpy_flight/idex_packet.py†L1002-L1290】
* **SNR estimation** via sliding baselines: `SNR = (max(x) - mu_noise) / sigma_noise`, with `mu_noise` and `sigma_noise` computed across the first 2 µs of pre-trigger samples.【F:src/spectrumpy_flight/SNR_Calculator.py†L12-L63】

### `drive_idex_packet.py`

*Launch command*: `python -m spectrumpy_flight.drive_idex_packet --inputs DIR [DIR ...] --out OUTPUTDIR [options]`

Key options:

| Option | Description |
| --- | --- |
| `--inputs DIR [DIR ...]` | One or more directories to scan for raw captures. |
| `--out OUTPUTDIR` | Destination folder for generated `.h5` files. |
| `--idex-script PATH` | Optional explicit path to `idex_packet.py`. |
| `--python EXE` | Interpreter used to spawn child processes (defaults to `sys.executable`). |
| `--max-procs N` | Maximum concurrent conversions. |
| `--threads-per-proc N` | Hard cap on BLAS/OpenMP threads inside each worker. |
| `--log-dir DIR` | Directory for per-job stdout/stderr captures. |
| `--nice N` | Apply POSIX niceness to reduce scheduling priority. |
| `--dry-run` | Print the planned workload without executing conversions. |

The driver recursively scans candidate trees, filters out already-converted captures, and launches multiple `idex_packet.py` subprocesses.  Jobs inherit a sanitized environment that caps BLAS threads:
`\text{threads}_\mathrm{BLAS} = \min\bigl(\text{CPU count}, \text{user limit}\bigr).`
The script writes stdout/stderr per job (`job_XXXX.out`, `job_XXXX.err`) to simplify post-run auditing.【F:src/spectrumpy_flight/drive_idex_packet.py†L24-L203】

### `combine_target_signals.py`

*Launch command*: `python -m spectrumpy_flight.combine_target_signals [--root DIR] [--csv OUTFILE] [--no-recursive]`

Scans decoded HDF5 directories and extracts the scalar dataset `/EVENT/Analysis/Target LImpactCharge` using a strict `ds[()]` read path.  The script normalises whitespace quirks in dataset names, supports recursive traversal, and optionally exports a CSV summary.

Outputs include per-file debug previews of the first few scalar reads plus statistics of finite values across all events.【F:src/spectrumpy_flight/combine_target_signals.py†L1-L220】  The parsing logic enforces scalar shapes, raising an exception if a dataset unexpectedly contains arrays larger than one element.

### `dust_composition.py`

*Launch command*: `python -m spectrumpy_flight.dust_composition [HDF5FILE]`

Provides a standalone window for fitting EMG and HyperEMG peaks, computing relative sensitivity corrections, and producing ternary abundance plots.  HyperEMG tails follow
`\mathrm{HyperEMG}(t) = \sum_{j} w_j \, \mathrm{EMG}(t; \mu, \sigma, \lambda_j), \qquad \sum_j w_j = 1.`
Outputs include CSV tables of fitted parameters, RSF-adjusted counts, and PDF figures.【F:src/spectrumpy_flight/dust_composition.py†L1-L4023】

### `dust_estimator_gui.py`

*Launch command*: `python -m spectrumpy_flight.dust_estimator_gui`

Estimates impact charge by applying calibration polynomials to TOF peaks.  The GUI accepts manual peak picks, performs linearised error propagation
`\sigma_m \approx \left|\frac{\partial m}{\partial t}\right| \sigma_t = \left|\frac{2 (t - t_0)}{a_1^2}\right| \sigma_t,`
assuming the first-order stretch dominates.  Results export as CSV snapshots for mission logs.【F:src/spectrumpy_flight/dust_estimator_gui.py†L1-L611】

### `ImpactBook.py`

*Launch command*: `python -m spectrumpy_flight.ImpactBook [-i INPUTDIR] [-e EXPERIMENT]`

Creates PDF “impact books” summarising events with stacked plots, channel notes, and metadata pages.  It leverages Matplotlib’s `PdfPages` to collate figures and writes interactive annotations describing trigger timing, saturation metrics, and SNR thresholds.【F:src/spectrumpy_flight/ImpactBook.py†L1-L820】

### `olivine_metrics.py`

*Launch command*: `python -m spectrumpy_flight.olivine_metrics [-f HDF5FILE] [--roi MASS1 MASS2]`

Runs a regression suite tailored to olivine calibration shots.  Gaussian line fits extract centroid drift, peak width, and asymmetry.  The algorithm minimises
`\chi^2 = \sum_{i} \frac{(y_i - y(t_i; A, \mu, \sigma))^2}{\sigma_i^2},`
with `lmfit.Model.fit`.  Output directories contain CSV summary tables, Matplotlib figures, and JSON parameter dumps.【F:src/spectrumpy_flight/olivine_metrics.py†L1-L2663】

### `noise_analysis.py`

*Launch command*: `python -m spectrumpy_flight.noise_analysis --file HDF5FILE --channel TOF_H`

Computes time-domain noise histograms, FFTs, and autocorrelation functions.  The autocorrelation coefficient at lag $\tau$ is
`\rho(\tau) = \frac{1}{(N-\tau) \sigma^2} \sum_{i=0}^{N-\tau-1} (x_i - \mu)(x_{i+\tau} - \mu),`
and confidence intervals follow Bartlett’s approximation to flag statistically significant oscillations.【F:src/spectrumpy_flight/noise_analysis.py†L1-L429】

### `time2mass.py`

*Launch command*: `python -m spectrumpy_flight.time2mass --input HDF5FILE`

Offers an interactive calibration assistant.  Users select peaks, adjust stretch/shift sliders, and export coefficient files.  Internally the script reuses the Newton inversion from [`mass_calibration.py`](../mass_calibration.py); full derivations appear in [Time-to-Mass Calibration Theory and Practice](time_to_mass_calibration.md).【F:src/spectrumpy_flight/time2mass.py†L1-L399】

### `read_from_ois.py`

*Launch command*: `python -m spectrumpy_flight.read_from_ois --input OISFILE --output DIR`

Reads GS-OASIS packet dumps, applies Rice decompression, and writes HDF5 snapshots.  The script shares most helper functions with `science_tool.py` but can run unattended for overnight ingest.  Output names follow the timestamp embedded in the filename (`ois_output_%Y%m%d_%H%M%S.h5`).【F:src/spectrumpy_flight/read_from_ois.py†L1-L203】

### `rice_decode.py`

*Launch command*: `python -m spectrumpy_flight.rice_decode -f RAWFILE [-o OUTFILE]`

Decodes standalone Rice-compressed products.  The decoder implements the canonical quotient–remainder split:
`\text{value} = q \cdot 2^k + r,`
where $q$ is extracted from unary-coded bits and $r$ from the following $k$ binary bits.  Optional `--interleave` toggles between telemetry and oscilloscope packing orders.【F:src/spectrumpy_flight/rice_decode.py†L1-L274】

---

## Data browsers and utilities

### `HDF_View.py`

*Launch command*: `python -m spectrumpy_flight.HDF_View FILE`

Provides a hierarchical tree explorer for HDF5 products with dataset previews, Matplotlib quick plots, and attribute inspection.  Large datasets stream lazily; slicing leverages NumPy-style indexing entered in the GUI.  Value histograms compute
`p(v) = \frac{1}{N}\sum_{i=1}^{N} \mathbf{1}[v_i \in \text{bin}(v)],`
enabling quick sanity checks on trigger distributions.【F:src/spectrumpy_flight/HDF_View.py†L1-L513】

### `HDF_Explorer.py`

Desktop browser for walking HDF5 and CDF trees simultaneously.  The viewer issues batched reads using `h5py.Dataset.read_direct` to avoid repeated locking overhead and displays attribute metadata alongside dataset previews.【F:src/spectrumpy_flight/HDF_Explorer.py†L1-L380】

### `CDF_View.py`

Specialised NASA CDF inspector with time conversion helpers.  The tool converts TT2000 timestamps to UTC strings via
`\text{UTC}(i) = \text{datetime64}(\text{nanoseconds since J2000} + \Delta_{\text{leap}}).`
Interactive widgets allow filtering by variable name, dimensionality, and attribute presence.【F:src/spectrumpy_flight/CDF_View.py†L1-L481】

### `SNR_Calculator.py`

Experimental batch SNR extraction for calibration campaigns.  For each event the script calls `time2mass.time2mass` to obtain the stretch/shift parameters, detects peaks with `scipy.signal.find_peaks`, and stores
`\kappa = \frac{1}{M} \sum_{m \in \text{peaks}} \left( m_\text{scale}(m) - \mathrm{round}\bigl(m_\text{scale}(m), 1\bigr) \right)`
as a fractional alignment metric.  Outputs include scatter plots of SNR versus $\kappa$ plus histograms of stretch and shift distributions.【F:src/spectrumpy_flight/SNR_Calculator.py†L1-L120】

### `line_shapes.py`

Defines Gaussian, Lorentzian, EMG, and HyperEMG profiles used across the GUI.  Functions expose gradients so lmfit can compute analytic Jacobians during optimisation.  For example, the Gaussian profile is
`G(t; A, \mu, \sigma) = A \exp\left(-\frac{(t-\mu)^2}{2\sigma^2}\right),`
with derivative $\partial G/\partial \mu = G(t)\,(t-\mu)/\sigma^2$.【F:src/spectrumpy_flight/line_shapes.py†L1-L181】

### `bitstream.py`

Handles low-level bit slicing for Rice/Golomb decoding.  The module implements constant-time masks using NumPy vector operations and optional CuPy acceleration.  Given a packed bit array, `_chunk_to_uint32` reconstructs integers as
`u_j = \sum_{b=0}^{31} s_{32j+b} 2^b,`
where $s_k$ are individual bits.【F:src/spectrumpy_flight/bitstream.py†L1-L213】

### `mass_calibration.py`

Library module defining `fit_tof_series`, `invert_tof_to_mass`, and `TOFMassCal`.  The polynomial fit solves a weighted least-squares problem in the reduced coordinate $s = \sqrt{m}$; see [Time-to-Mass Calibration Theory and Practice](time_to_mass_calibration.md) for the full derivation.  The class exposes `mass_to_tof`, `tof_to_mass`, and `derivative` helpers used throughout the quicklook stack.【F:src/spectrumpy_flight/mass_calibration.py†L1-L144】

---

## Shell helpers

| Script | Purpose |
| --- | --- |
| `process_packets.sh` | Converts raw captures into HDF5/CDF files, skipping already processed outputs.  Accepts `-f/--file/--filename` or directory arguments and defaults to scanning `Data/`.【F:process_packets.sh†L1-L56】 |
| `process_EM_Data.sh` | Mirror of `process_packets.sh` tuned for electromagnetic calibration directories; includes retry loops and logging to `logs/`.【F:process_EM_Data.sh†L1-L37】 |
| `run_science.sh` | Wraps `python -m spectrumpy_flight.run_all` after exporting `PYTHONPATH` so the repository tree is importable without installation.【F:run_science.sh†L1-L9】 |
| `analyze_velocities.sh` | Batch-processes `.trc` oscilloscope folders by invoking `qd_quicklook.py -s DIR -t LABEL` for every subdirectory under the configured `BASE_DIR`.【F:analyze_velocities.sh†L1-L22】 |

Each helper preserves the original exit codes of the underlying Python modules, allowing shell pipelines to detect failures reliably.

---

## Data resources

The `src/spectrumpy_flight/Data/` and `src/spectrumpy_flight/HDF5/` directories hold sample packets and decoded outputs.  Analysts can replicate tutorials by opening these references inside the quicklook GUIs.  These assets live only in the source repository to keep the PyPI wheels lightweight—fetch them from GitHub when you need practice data.  The `CDF/` tree contains mission-specific variable definitions used by `CDF_View.py` and the NASA CDF export path in `idex_packet.py`.

---

## Further reading

* [Time-to-Mass Calibration Theory and Practice](time_to_mass_calibration.md) – Detailed derivations and recipes for TOF mass calibration.
* [`quicklook_tutorial.md`](quicklook_tutorial.md) – Hands-on walkthrough of the GUI.
* [`fitting_reference.md`](fitting_reference.md) – Data Browser workflow, variable definitions, and target overlays.
* [`olivine_metrics.md`](olivine_metrics.md) – Specialist calibration campaign guidance.
