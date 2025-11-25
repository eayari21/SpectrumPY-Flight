# SpectrumPY Flight Signal Fitting & Operations Guide

## 1. Fitting models and evaluation flow

The automated fitting stack is concentrated in the packet decoder and quicklook toolchain. Raw waveform slices are preconditioned, model parameters are optimised with SciPy's `curve_fit`, and the resulting curves are persisted alongside derived scalars. Unless stated otherwise, time-like parameters are expressed in microseconds (µs), exponential rate coefficients are in inverse microseconds (µs⁻¹), and amplitudes are in detector charge units (pC/Δt for TOF channels or µA-equivalent for grid currents). The core models are:

### 1.1 Exponentially modified Gaussian (EMG)

Targeted at charge pick-up channels, the EMG definition inside `FitEMG` is

`E(x; \mu, \sigma, \lambda) = \frac{\lambda}{2} \exp\Bigl(\tfrac{\lambda}{2} (2\mu + \lambda \sigma^2 - 2x)\Bigr) \mathrm{erfc}\left(\frac{\mu + \lambda \sigma^2 - x}{\sqrt{2}\,\sigma}\right),`

with $(\mu, \sigma, \lambda)$ initialised from the waveform maximum, its spread, and the acquisition span respectively. The helper integrates the fitted pulse to recover charge by calling `quad(E, x_{\min}, x_{\max})`. Failed fits return sentinel zeros so that downstream aggregation can skip the record.【F:idex_packet.py†L58-L119】【F:idex_packet.py†L120-L186】

### 1.2 Baseline removal and noise suppression

Before target fits run, `FitTargetSignal` extracts a noise-only window, detrends it with a linear model,

`L(t; a, b) = a t + b,`

where $a$ is the slope (pC/Δt per µs) and $b$ is the intercept (pC/Δt). It optionally subtracts a coherent microphonic tone

`S(t; c, d, e) = c \sin(d t + e),`

with amplitude $c$ (pC/Δt), angular frequency $d$ (rad·µs⁻¹), and phase $e$ (rad). Residuals are low-pass filtered via `butter_lowpass_filter` before restricting the fit domain to the pre-trigger region that contains the physical response.【F:idex_packet.py†L188-L282】 Linear and sinusoidal primitives are shared with the oscilloscope-based tools so manual analyses can reproduce the automated cleaning.【F:qd_quicklook.py†L31-L70】【F:ImpactBook.py†L25-L64】 Saturation streaks are masked via `detect_saturation`, which flags plateaus longer than the configurable repeat threshold so gain-substituted traces fall back to medium/low channels when necessary.【F:idex_packet.py†L165-L199】

### 1.3 Ion grid response

The instrument response model couples a Heaviside rise with exponential growth and decay to capture collector charging:

`G(t; P_0, P_1, P_4, P_5, P_6) = P_1 + H(t-P_0)\, P_4 \bigl(1 - e^{-(t-P_0)/P_5}\bigr) e^{-(t-P_0)/P_6}.`

`FitTargetSignal` seeds the optimiser with nominal rise/decay constants and solves for $(P_0, P_1, P_4, P_5, P_6)$ on the filtered residual. The resulting waveform and peak-to-baseline amplitude are archived for the quicklook overlays.【F:idex_packet.py†L60-L116】【F:idex_packet.py†L232-L276】 Identical physics is available in the oscilloscope processing path, keeping lab-calibrated fits consistent with flight decoding.【F:qd_quicklook.py†L71-L120】

### 1.4 Pickup tube trajectory model

Velocity estimates in the oscilloscope (`ImpactBook` / `qd_quicklook`) path use a piecewise pickup tube model that integrates constant acceleration through the tube apertures and an RC decay when the particle leaves the field of view. The piecewise construction is defined as

`Q(t; t_0, q, v) = { q/Δt1 * (t - t0) for t0 < t < t0 + Δt1; q * exp(-(t - (t0 + Δt1))/τ) for t0 + Δt1 < t < t0 + Δt1 + Δt2; linear exit and exponential tail elsewhere }`

where $\Delta t_1 = d / v$ is the transit across an aperture, $\Delta t_2 = (\ell - 2d)/v$ is the interior flight time, and $\tau$ encodes the PUT RC constant. The full implementation enumerates each interval explicitly and feeds the composite signal into `curve_fit`. This mirrors the flight pipeline so velocity histograms agree once traces are exported to CSV.【F:qd_quicklook.py†L51-L123】

### 1.5 Dataset harvesting in Quicklook

Within the Qt viewer, each `HDF5DataSource`/`CDFDataSource` instance scans `/event/Analysis` groups for arrays whose normalised names contain `fit`, splitting the products into time bases, best-fit values, and saved parameter vectors. Any manual override applied in the UI is written back via `_recalculate_fit` so follow-up loads replay the revised solution.【F:IDEX-quicklook.py†L1009-L1220】【F:IDEX-quicklook.py†L2303-L2831】 Rise metrics are kept in sync with the waveform cache, ensuring the `10%`/`90%` crossing times and amplitudes stay coupled to the active parameter set.【F:idex_analysis_utils.py†L5-L101】【F:IDEX-quicklook.py†L2705-L2831】

### 1.6 Mass line shape library

Mass composition workflows rely on a shared catalogue of analytic line shapes implemented in `line_shapes.py`. The module provides Gaussian, Lorentzian, Voigt, EMG, double-EMG, HyperEMG, and generalised-normal profiles with consistent area-normalisation so integrated abundances remain comparable across species.【F:line_shapes.py†L17-L126】 The dust composition UI layers additional metadata on top of those primitives—each line shape advertises the parameter labels, tooltips, and any auxiliary controls (for example secondary time constants or weight lists) needed to collect user input.【F:dust_composition.py†L1387-L1478】

When you inspect a mass line, the dialog evaluates the selected shape against a dense time grid, annotates the waveform overlay, and updates abundance estimates by integrating the analytic curve within the chosen time window. Abundances are normalised against the baseline-corrected TOF energy so changing line shapes or tail parameters directly feeds into the percentage summaries and exported tables.【F:dust_composition.py†L3208-L3295】【F:dust_composition.py†L3336-L3484】

### 1.7 Signal-to-noise and rise metric evaluation

The packet decoder computes a per-channel signal-to-noise ratio by comparing the peak excursion of the filtered trace to the standard deviation of a quiet-time segment:

`\mathrm{SNR} = \frac{\max_t y(t) - \bar{y}_\mathrm{baseline}}{\sigma_\mathrm{baseline}},`

where the baseline statistics are drawn from the automatically selected window described in §1.2.【F:idex_packet.py†L205-L246】 Quicklook supplements this with 10–90% rise metrics that interpolate when the fit curve crosses each threshold, storing $t_{10}$, $t_{90}$, and their corresponding amplitudes in `/Analysis/<channel> 10pct Time`, etc.【F:idex_analysis_utils.py†L5-L101】【F:IDEX-quicklook.py†L2705-L2831】 These rise times (µs) and amplitudes (pC/Δt) are displayed alongside the fitted traces to confirm slew-rate expectations.

### 1.8 Fit editing pipeline

Launching **Edit Fit Parameters** (toolbar button or `Ctrl+E`) opens the `FitParameterDialog`, which enumerates every dataset carrying editable coefficients, renders the analytic form, and exposes per-parameter editors. The redesigned dialog groups datasets by channel, highlights the active waveform with a matplotlib preview, and lets you place baselines, amplitudes, and timing anchors directly on the plot via selector buttons before recomputing the analytic curve.【F:IDEX-quicklook.py†L3520-L3927】 When you apply changes, `_recalculate_fit` re-evaluates the curve against the stored time axis, caches the override in memory, and optionally writes the updated parameters plus regenerated fit arrays back into the HDF5 analysis group. Status messages confirm whether the save succeeded, while resets purge both parameter and curve overrides so subsequent loads return to the file-backed solution.【F:IDEX-quicklook.py†L3927-L4172】

## 2. Command-line entry points

Every script that exposes an `argparse` interface is listed below. Commands assume you are in the repository root and your environment satisfies GUI or packet dependencies as needed.

| Script | Purpose | Invocation |
| --- | --- | --- |
| `IDEX-quicklook.py` | Launch the Qt Quicklook with optional file/event preselection. | `python IDEX-quicklook.py --filename <HDF/CDF> --eventnumber <index>` (flags optional).【F:IDEX-quicklook.py†L3180-L3201】 |
| `Scope-IDEX-quicklook.py` | Alternate Quicklook build for scope datasets. | `python Scope-IDEX-quicklook.py --filename <HDF5> --eventnumber <index>`; both flags optional.【F:Scope-IDEX-quicklook.py†L827-L844】 |
| `imfpy/gui/main.py` | Unified launcher with `launcher` and `quicklook` subcommands. | `python -m imfpy.gui.main launcher` for the welcome screen or `python -m imfpy.gui.main quicklook --filename <path> --eventnumber <index>`.【F:imfpy/gui/main.py†L1-L129】 |
| `ImpactBook.py` | Batch-fit oscilloscope TRC directories into HDF5/CSV outputs. | `python ImpactBook.py --trcdir <directory> --experimentname <label>` (both required).【F:ImpactBook.py†L794-L819】 |
| `qd_quicklook.py` | Legacy oscilloscope pipeline (calls `ImpactBook`). | `python qd_quicklook.py --sourcefolder <trc_dir> --targetfolder <experiment>`; required flags pass through to `ImpactBook`.【F:qd_quicklook.py†L900-L912】 |
| `idex_packet.py` | Decode raw L0 telemetry, compute fits, emit HDF5/CDF. | `python idex_packet.py --file <Data/packet.bin>` (required).【F:idex_packet.py†L1397-L1424】 |
| `rice_decode.py` | Standalone Rice decompressor for science frames. | `python rice_decode.py --sourcefile <compressed.bin> --targetfile <decoded.bin>`; optional `sampleCount` constant is set in-code.【F:rice_decode.py†L275-L298】 |
| `combine_target_signals.py` | Scan HDF5 trees for `Analysis/Target LImpactCharge`. | `python combine_target_signals.py --root <path> [--no-recursive] [--csv output.csv]`.【F:combine_target_signals.py†L146-L200】 |
| `tools/cdf_strings_report.py` | Summarise metadata strings extracted from a CDF ASCII dump. | `python tools/cdf_strings_report.py <file.cdf> [--variables V1 V2 …] [--window N]`. Use `-h/--help` to review every optional flag and usage example before running it against mission data.【F:tools/cdf_strings_report.py†L160-L214】 |
| `HDF_View.py` | Lightweight HDF5 browser used by Quicklook. | `python HDF_View.py <file.h5>`; omitting the path opens a file dialog. Returns Qt exit code.【F:HDF_View.py†L1-L197】【F:HDF_View.py†L400-L419】 |
| `ImpactBook.py` / `qd_quicklook.py` CSV outputs | Use GUI to visualise histograms after CLI execution. | Run the commands above; outputs are stored next to the working directory. |
| `imfpy/gui/main.py` packaged usage | Package entry point respects the same flags via `python -m`. | As above; import guard injects the repo root when run from source.【F:imfpy/gui/main.py†L1-L75】 |

If a script imports `argparse` but lacks a parser (for example, `science_tool.py` or `read_from_ois.py`), it is intended for interactive execution or environment-specific workflows and should be customised before use.

### 2.1 Detailed CLI behaviour

* **`idex_packet.py`** reads binary frames, decompresses optional Rice blocks, fans out waveform preparation across up to `os.cpu_count()` worker threads, and emits HDF5/CDF products. Use `--file` for a single packet or integrate with the driver in §7. Parallel waveform preparation reuses the same fitting routines described in §1.【F:idex_packet.py†L1189-L1245】【F:idex_packet.py†L1397-L1424】
* **`drive_idex_packet.py`** orchestrates batch conversions. `--inputs` accepts one or more roots with extensionless packet files, `--out` selects the HDF5 destination, `--max-procs` caps concurrent decoder launches, and `--threads-per-proc` injects BLAS/OpenMP limits so nested parallelism stays bounded.【F:drive_idex_packet.py†L1-L112】
* **`tools/process_packets_parallel.py`** powers the shell wrappers. It discovers raw packet files, skips events that already have HDF5 output, prints a resource summary, and launches subprocess decoders with a worker count derived from CPU and memory availability (override with `--max-workers` or `--memory-per-job`).【F:tools/process_packets_parallel.py†L1-L118】【F:tools/resource_utils.py†L1-L75】
* **`ImpactBook.py`** processes oscilloscope TRC trees. Pair `--trcdir` and `--experimentname` to generate HDF5 and CSV outputs that mirror the packet decoder’s structure.【F:ImpactBook.py†L794-L819】
* **`combine_target_signals.py`** traverses analysis trees, optionally writing combined CSV summaries. `--no-recursive` restricts the scan to a single directory, while `--csv` selects an export path.【F:combine_target_signals.py†L146-L200】

### 2.2 Quicklook GUI shortcuts

The Qt viewer exposes its tooling both via the toolbar and keyboard accelerators:

* **Data Browser** — `View → Open Data Browser` or `Ctrl+B` opens the HDF tree explorer described in §5.【F:IDEX-quicklook.py†L1992-L2055】
* **Variable Definitions** — `View → Variable Definitions…` loads the calibration polynomial viewer when available.【F:IDEX-quicklook.py†L1996-L2055】
* **Dust Composition** — Toolbar button or `Ctrl+D` launches the mass inference window for the active event.【F:IDEX-quicklook.py†L2139-L2196】【F:IDEX-quicklook.py†L2508-L2579】
* **Noise Analysis** — Toolbar button or `Ctrl+N` spawns the diagnostics window detailed in §4.【F:IDEX-quicklook.py†L2027-L2188】【F:IDEX-quicklook.py†L2579-L2656】
* **Edit Fit Parameters** — Toolbar button or menu action opens the parameter dialog covered in §1.8.【F:IDEX-quicklook.py†L2358-L2377】【F:IDEX-quicklook.py†L3450-L3641】

Additional channel toggles, overlay controls, and export buttons are grouped on the left-hand panel so the same interactions work whether the viewer is launched from the CLI scripts above or from the `imfpy` launcher module.【F:IDEX-quicklook.py†L2267-L2385】

## 3. Dust composition and mass inference

Opening **Dust Composition** builds a combined TOF waveform by stitching saturated high-gain samples with rescaled medium- and low-gain data using the saturation masks computed in §1.2. The combination honours gain factors stored in `GAIN_MAP` so all segments share the high-gain scale.【F:dust_composition.py†L3208-L3295】 Users can:

1. Toggle individual channels or the combined spectrum. The tool overlays both mass (amu) and time (µs) axes; a twinned Matplotlib axis keeps the time grid active for baseline edits.【F:dust_composition.py†L3208-L3295】【F:dust_composition.py†L3336-L3421】
2. Adjust baselines, stretch, and shift parameters. Baseline edits are applied immediately to abundance calculations, while stretch/shift reparameterise the time-to-mass conversion via $m = \text{stretch}\,(t-\text{shift})$.【F:dust_composition.py†L3388-L3484】
3. Add or edit mass lines. Default seeds pull the strongest peak from the combined waveform, estimate $\mu$ (µs), $\sigma$ (µs), $\lambda$ (µs⁻¹), and the integration window, then integrate the selected line shape to compute abundance fractions.【F:dust_composition.py†L3392-L3484】 Manual edits trigger `_recompute_mass_line`, update the abundance tables, and refresh the preview plot.【F:dust_composition.py†L3208-L3335】
4. Apply relative sensitivity factors (RSFs) and compare against the built-in material library. The RSF dialog lets you enable/disable lines, assign per-line factors, and renormalise contributions before exporting, while the sample guess picker ranks catalogued materials and ternary mixtures that best explain the observed abundances.【F:dust_composition.py†L3008-L3525】【F:dust_composition.py†L3526-L4023】【F:dust_composition.py†L4843-L5160】 Ternary plots, RSF-adjusted tables, and mixture diagnostics update instantly so you can iterate on calibrations and reference selections in one session.【F:dust_composition.py†L2139-L2764】【F:dust_composition.py†L3828-L3953】
5. Persist results back to HDF5. Saving writes the combined waveform, baseline, mass calibration parameters, RSF selections, and every line’s fitted curve under `/Analysis/DustComposition`, allowing subsequent Quicklook sessions to reload manual adjustments.【F:dust_composition.py†L3919-L4023】

## 4. Noise diagnostics and baseline studies

The **Noise Analysis** window renders three complementary views for any channel with available data:

1. **Amplitude histogram** — Builds a histogram with adaptive bins and overlays a Gaussian with mean $\mu$ and standard deviation $\sigma$ to quantify the baseline distribution in detector units.【F:noise_analysis.py†L260-L328】
2. **Power spectrum** — Computes the FFT of the detrended signal, plots $|\mathcal{F}(f)|$ on log–log axes, and shades the spectral magnitude to highlight coherent tones.【F:noise_analysis.py†L328-L346】
3. **Autocorrelation** — Normalises the autocorrelation by the variance and presents it against sample lag or time lag depending on the inferred sampling cadence.【F:noise_analysis.py†L346-L386】

Summary cards report the Gaussian amplitude (counts), Poisson-equivalent shot noise $\sqrt{\bar{y}}$, and RMS noise, all labelled with the appropriate channel units. Status messages annotate the number of samples analysed and the median Δt recovered from the time array.【F:noise_analysis.py†L386-L429】 Channels, units, and event lists are pre-populated from the Quicklook definitions, so switching channels or hopping to the next event updates axes, legend text, and summary annotations without reloading the window.【F:noise_analysis.py†L55-L221】【F:IDEX-quicklook.py†L2579-L2656】

## 5. Correlating the Data Browser with variable definitions

1. Open a science product in Quicklook and press **Ctrl +B** or choose `View → Open Data Browser`. The viewer delegates to the `HDF_View` tree, exposing `/event/Analysis` datasets, attributes, and previews in a split-pane layout.【F:IDEX-quicklook.py†L1559-L1627】【F:HDF_View.py†L1-L197】
2. Select any dataset to populate the summary and preview tables. The browser shows the dataset path, shape, dtype, and quick statistics so you can pinpoint the exact DN source that feeds a control or housekeeping register.【F:HDF_View.py†L200-L299】
3. Launch the variable definitions window via `View → Variable Definitions…`. The companion Qt tool loads the Excel translation sheet, filters by CDF or packet field, and renders piecewise polynomial coefficients for every DN span.【F:IDEX-quicklook.py†L1559-L1627】【F:IDEX_Definitions_View.py†L1-L199】
4. Enter the DN reported in the data browser into the `DN:` spin box. The viewer evaluates the appropriate polynomial and classifies the reading as *low*, *high*, *below*, or *exceeding* the calibrated range, letting you triage instrument states at a glance.【F:IDEX_Definitions_View.py†L200-L260】
5. Cross-reference the packet field name shown in the definitions window with the dataset path from the browser to confirm whether you are looking at commanded set-points, measured telemetry, or derived metrics. Because both tools share the same Excel-backed catalog, the conversions match the pipeline used during file generation.

## 6. Incorporating images and rich media

The documentation center preloads every Markdown file under `docs/` and renders them with Qt’s GitHub/CommonMark dialect. Relative links resolve against the source directory thanks to `QTextDocument.setBaseUrl`, so static assets placed under `docs/media/` or adjacent folders are available at runtime.【F:IDEX-quicklook.py†L240-L336】【F:IDEX-quicklook.py†L400-L470】 To include an illustration:

1. Save the image under `docs/media/` (for example, `docs/media/emg_fit.png`).
2. Reference it with standard Markdown syntax, `![EMG fit overlay](media/emg_fit.png)`. The viewer loads the resource when the documentation is opened.
3. For animated GIFs or video walkthroughs, export the media to a shareable format and embed it using Markdown image syntax (GIF) or a hyperlink to an external hosting service for large MP4 files. Qt’s viewer follows the link in the default browser.

Because Markdown is rendered directly in the GUI and GitHub, the same instructions apply to repository README updates and in-application help.

## 7. Analytical line shapes

The repository standardises line shapes in `line_shapes.py`; each profile assumes an area-normalised amplitude so integrated abundance (pC) corresponds to the `area` argument.【F:line_shapes.py†L1-L130】 Parameter units follow the same conventions as §1.

* **Gaussian** —
`G(x; \mu, \sigma, A) = \frac{A}{\sigma \sqrt{2\pi}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right),`
  where $\mu$ is the peak time (µs) and $\sigma$ is the standard deviation (µs).【F:line_shapes.py†L13-L24】

* **Lorentzian** —
`L(x; \mu, \gamma, A) = \frac{A}{\pi} \frac{\gamma}{(x-\mu)^2 + \gamma^2},`
  with half-width at half-maximum $\gamma$ (µs).【F:line_shapes.py†L27-L37】

* **Voigt** — A convolution of Gaussian and Lorentzian components using SciPy’s Faddeeva implementation when available; otherwise a pseudo-Voigt approximation blends both shapes by weight $\eta$ determined from the effective full-width at half-maximum.【F:line_shapes.py†L39-L73】

* **Exponentially modified Gaussian (EMG)** —
`\mathrm{EMG}(x; \mu, \sigma, \tau, A) = \frac{A}{2\tau} \exp\left(\frac{\tau^{-1}}{2} \left(2\mu + \tau^{-1}\sigma^2 - 2x\right)\right) \mathrm{erfc}\left(\frac{\mu + \tau^{-1}\sigma^2 - x}{\sqrt{2}\sigma}\right),`
  where $\tau$ is the exponential tail time constant (µs). A negative $\tau$ models a left tail.【F:line_shapes.py†L75-L101】

* **Double EMG** — Convex combination of two EMGs with shared $\mu$ and $\sigma$, separate $\tau_1$/$\tau_2$, and mixing weight $w_1$.【F:line_shapes.py†L103-L115】

* **HyperEMG** — Gaussian core plus optional left/right EMG tails with user-supplied time constants and weights that renormalise automatically.【F:line_shapes.py†L117-L162】

* **Generalised normal (Subbotin)** —
`N(x; \mu, \alpha, \beta, A) = \frac{A\beta}{2\alpha\,\Gamma(1/\beta)} \exp\left(-\left|\frac{x-\mu}{\alpha}\right|^{\beta}\right),`
  where $\alpha$ (µs) controls scale and $\beta$ sets the shape exponent.【F:line_shapes.py†L164-L187】

These profiles appear throughout the dust composition UI and the quicklook fit dialog, ensuring on-screen documentation and exported CSVs reference identical definitions.【F:dust_composition.py†L1387-L1478】【F:IDEX-quicklook.py†L3450-L3568】

## 8. Parallel decoding and optimisation strategy

The decoder combines intra-event threading with optional GPU acceleration:

* **Thread-level parallelism** — `_prepare_waveform` runs inside a `ThreadPoolExecutor`, allowing independent channel fits (TOF high/mid/low, target, grid) to execute concurrently per event.【F:idex_packet.py†L1189-L1245】 The driver script and the new `tools/process_packets_parallel.py` helper wrap this with process-level concurrency so multiple packet files are decoded simultaneously while respecting CPU and memory limits.【F:drive_idex_packet.py†L58-L112】【F:tools/process_packets_parallel.py†L1-L118】
* **Thread affinity for libraries** — `drive_idex_packet.py` injects `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and related environment variables so BLAS-backed FFT and fitting calls stay within the requested per-process budgets, avoiding oversubscription when dozens of events are decoded in parallel.【F:drive_idex_packet.py†L63-L87】 The parallel helper uses the same environment and exposes `--max-workers` plus `--memory-per-job` switches for tighter control.【F:tools/process_packets_parallel.py†L58-L115】
* **GPU acceleration** — The ASCII waveform unpacker attempts to offload bit-matrix dot-products to CuPy when arrays exceed 32 k samples, converting them back to NumPy once decoded.【F:idex_packet.py†L38-L84】 If CuPy is unavailable or the transfer fails, the code falls back to NumPy seamlessly.
* **Streaming safeguards and profiling** — `detect_saturation` opportunistically executes on CuPy arrays while falling back to NumPy, and `calculate_snr` reuses the same baseline windows selected during fitting so gain substitution logic and quality metrics remain identical regardless of hardware.【F:idex_packet.py†L165-L246】 Before launching workers the parallel helper prints CPU counts, available memory, and GPU detections, making it easier to validate resource assumptions for large backfills.【F:tools/process_packets_parallel.py†L70-L118】【F:tools/resource_utils.py†L1-L75】

These layers keep the science pipeline responsive on laptops (single-process decoding with NumPy) and on workstations (multi-process driver with GPU offload).

## 9. End-to-end workflow: decoding, fitting, and interpreting an event

1. **Decode Level-0 telemetry.** Run `python idex_packet.py --file Data/<packet.bin>` to parse XTCE definitions, decompress Rice frames, and emit `/event/*` HDF5 groups populated with raw waveforms and analysis scaffolding.【F:idex_packet.py†L1397-L1424】
2. **Inspect raw fits.** Open the generated `.h5` file in Quicklook. Each channel toggle can display both the original trace and the stored fit curves harvested via `gather_fit_data`, ensuring you see the exact model saved on disk.【F:IDEX-quicklook.py†L1009-L1220】【F:IDEX-quicklook.py†L2554-L2768】
3. **Adjust parameters when necessary.** Use `Edit Fit Parameters` to tweak EMG or ion-grid coefficients. Overrides call `_recalculate_fit`, persist the new arrays under `/Analysis/Fits`, and mark the event so future sessions pick up the change.【F:IDEX-quicklook.py†L2705-L2831】
4. **Assign mass axes and baselines.** Baseline-subtracted fits combine the linear and sinusoidal estimates from `FitTargetSignal` with the physical grid model, so toggling baseline overlays clarifies how the detector responds relative to a zero-charge reference.【F:idex_packet.py†L188-L276】【F:IDEX-quicklook.py†L2323-L2670】 The shared `time2mass` conversion utilities can be invoked from the analysis notebook layer to relabel axes once the trigger time is confirmed.【F:time2mass.py†L1-L120】
5. **Guess compositions.** With cleaned fits visible, open the dust composition tools or export CSVs for offline analysis. Composition helpers (see `dust_composition.py`) read the same fitted amplitudes and time-of-flight metrics, so the workflow from Quicklook to composition inference remains consistent.【F:dust_composition.py†L1260-L1300】
6. **Document findings.** Capture screenshots or attach annotated plots via the export menu. When writing mission logs, embed the figures in Markdown using the approach from Section 4 so the documentation center and GitHub present identical views.【F:IDEX-quicklook.py†L1559-L1659】

Following this pipeline keeps Level-0 decoding, fit refinement, and interpretation aligned across command-line utilities, GUI tooling, and downstream documentation.
