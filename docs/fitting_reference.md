# SpectrumPY Flight Signal Fitting & Operations Guide

## 1. Fitting models and evaluation flow

The automated fitting stack is concentrated in the packet decoder and quicklook toolchain. Raw waveform slices are preconditioned, model parameters are optimised with SciPy's `curve_fit`, and the resulting curves are persisted alongside derived scalars. The core models are:

### 1.1 Exponentially modified Gaussian (EMG)

Targeted at charge pick-up channels, the EMG definition inside `FitEMG` is

$$
E(x; \mu, \sigma, \lambda) = \frac{\lambda}{2} \exp\Bigl(\tfrac{\lambda}{2} (2\mu + \lambda \sigma^2 - 2x)\Bigr) \operatorname{erfc}\left(\frac{\mu + \lambda \sigma^2 - x}{\sqrt{2}\,\sigma}\right),
$$

with $(\mu, \sigma, \lambda)$ initialised from the waveform maximum, its spread, and the acquisition span respectively. The helper integrates the fitted pulse to recover charge by calling `quad(E, x_{\min}, x_{\max})`. Failed fits return sentinel zeros so that downstream aggregation can skip the record.【F:idex_packet.py†L58-L119】【F:idex_packet.py†L120-L186】

### 1.2 Baseline removal and noise suppression

Before target fits run, `FitTargetSignal` extracts a noise-only window, detrends it with a linear model,

$$
L(t; a, b) = a t + b,
$$

and optionally subtracts a coherent microphonic tone

$$
S(t; c, d, e) = c \sin(d t + e).
$$

Residuals are low-pass filtered via `butter_lowpass_filter` before restricting the fit domain to the pre-trigger region that contains the physical response.【F:idex_packet.py†L188-L282】 Linear and sinusoidal primitives are shared with the oscilloscope-based tools so manual analyses can reproduce the automated cleaning.【F:qd_quicklook.py†L31-L70】【F:ImpactBook.py†L25-L64】

### 1.3 Ion grid response

The instrument response model couples a Heaviside rise with exponential growth and decay to capture collector charging:

$$
G(t; P_0, P_1, P_4, P_5, P_6) = P_1 + H(t-P_0)\, P_4 \bigl(1 - e^{-(t-P_0)/P_5}\bigr) e^{-(t-P_0)/P_6}.
$$

`FitTargetSignal` seeds the optimiser with nominal rise/decay constants and solves for $(P_0, P_1, P_4, P_5, P_6)$ on the filtered residual. The resulting waveform and peak-to-baseline amplitude are archived for the quicklook overlays.【F:idex_packet.py†L60-L116】【F:idex_packet.py†L232-L276】 Identical physics is available in the oscilloscope processing path, keeping lab-calibrated fits consistent with flight decoding.【F:qd_quicklook.py†L71-L120】

### 1.4 Pickup tube trajectory model

Velocity estimates in the oscilloscope (`ImpactBook` / `qd_quicklook`) path use a piecewise pickup tube model that integrates constant acceleration through the tube apertures and an RC decay when the particle leaves the field of view. The piecewise construction is defined as

$$
Q(t; t_0, q, v) =
\begin{cases}
\frac{q}{\Delta t_1} (t - t_0) & t_0 < t < t_0 + \Delta t_1, \\
q\,e^{-\frac{t-(t_0+\Delta t_1)}{\tau}} & t_0 + \Delta t_1 < t < t_0 + \Delta t_1 + \Delta t_2, \\
\text{linear exit and exponential tail} & \text{elsewhere},
\end{cases}
$$

where $\Delta t_1 = d / v$ is the transit across an aperture, $\Delta t_2 = (\ell - 2d)/v$ is the interior flight time, and $\tau$ encodes the PUT RC constant. The full implementation enumerates each interval explicitly and feeds the composite signal into `curve_fit`. This mirrors the flight pipeline so velocity histograms agree once traces are exported to CSV.【F:qd_quicklook.py†L51-L123】

### 1.5 Dataset harvesting in Quicklook

Within the Qt viewer, each `HDF5DataSource`/`CDFDataSource` instance scans `/event/Analysis` groups for arrays whose normalised names contain `fit`, splitting the products into time bases, best-fit values, and saved parameter vectors. Any manual override applied in the UI is written back via `_recalculate_fit` so follow-up loads replay the revised solution.【F:IDEX-quicklook.py†L1009-L1220】【F:IDEX-quicklook.py†L2303-L2831】

### 1.6 Mass line shape library

Mass composition workflows rely on a shared catalogue of analytic line shapes implemented in `line_shapes.py`. The module provides Gaussian, Lorentzian, Voigt, EMG, double-EMG, HyperEMG, and generalised-normal profiles with consistent area-normalisation so integrated abundances remain comparable across species.【F:line_shapes.py†L17-L126】 The dust composition UI layers additional metadata on top of those primitives—each line shape advertises the parameter labels, tooltips, and any auxiliary controls (for example secondary time constants or weight lists) needed to collect user input.【F:dust_composition.py†L1387-L1478】

When you inspect a mass line, the dialog evaluates the selected shape against a dense time grid, annotates the waveform overlay, and updates abundance estimates by integrating the analytic curve within the chosen time window. Abundances are normalised against the baseline-corrected TOF energy so changing line shapes or tail parameters directly feeds into the percentage summaries and exported tables.【F:dust_composition.py†L3208-L3295】【F:dust_composition.py†L3336-L3484】

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

## 3. Correlating the Data Browser with variable definitions

1. Open a science product in Quicklook and press **Ctrl +B** or choose `View → Open Data Browser`. The viewer delegates to the `HDF_View` tree, exposing `/event/Analysis` datasets, attributes, and previews in a split-pane layout.【F:IDEX-quicklook.py†L1559-L1627】【F:HDF_View.py†L1-L197】
2. Select any dataset to populate the summary and preview tables. The browser shows the dataset path, shape, dtype, and quick statistics so you can pinpoint the exact DN source that feeds a control or housekeeping register.【F:HDF_View.py†L200-L299】
3. Launch the variable definitions window via `View → Variable Definitions…`. The companion Qt tool loads the Excel translation sheet, filters by CDF or packet field, and renders piecewise polynomial coefficients for every DN span.【F:IDEX-quicklook.py†L1559-L1627】【F:IDEX_Definitions_View.py†L1-L199】
4. Enter the DN reported in the data browser into the `DN:` spin box. The viewer evaluates the appropriate polynomial and classifies the reading as *low*, *high*, *below*, or *exceeding* the calibrated range, letting you triage instrument states at a glance.【F:IDEX_Definitions_View.py†L200-L260】
5. Cross-reference the packet field name shown in the definitions window with the dataset path from the browser to confirm whether you are looking at commanded set-points, measured telemetry, or derived metrics. Because both tools share the same Excel-backed catalog, the conversions match the pipeline used during file generation.

## 4. Incorporating images and rich media

The documentation center preloads every Markdown file under `docs/` and renders them with Qt’s GitHub/CommonMark dialect. Relative links resolve against the source directory thanks to `QTextDocument.setBaseUrl`, so static assets placed under `docs/media/` or adjacent folders are available at runtime.【F:IDEX-quicklook.py†L240-L336】【F:IDEX-quicklook.py†L400-L470】 To include an illustration:

1. Save the image under `docs/media/` (for example, `docs/media/emg_fit.png`).
2. Reference it with standard Markdown syntax, `![EMG fit overlay](media/emg_fit.png)`. The viewer loads the resource when the documentation is opened.
3. For animated GIFs or video walkthroughs, export the media to a shareable format and embed it using Markdown image syntax (GIF) or a hyperlink to an external hosting service for large MP4 files. Qt’s viewer follows the link in the default browser.

Because Markdown is rendered directly in the GUI and GitHub, the same instructions apply to repository README updates and in-application help.

## 5. End-to-end workflow: decoding, fitting, and interpreting an event

1. **Decode Level-0 telemetry.** Run `python idex_packet.py --file Data/<packet.bin>` to parse XTCE definitions, decompress Rice frames, and emit `/event/*` HDF5 groups populated with raw waveforms and analysis scaffolding.【F:idex_packet.py†L1397-L1424】
2. **Inspect raw fits.** Open the generated `.h5` file in Quicklook. Each channel toggle can display both the original trace and the stored fit curves harvested via `gather_fit_data`, ensuring you see the exact model saved on disk.【F:IDEX-quicklook.py†L1009-L1220】【F:IDEX-quicklook.py†L2554-L2768】
3. **Adjust parameters when necessary.** Use `Edit Fit Parameters` to tweak EMG or ion-grid coefficients. Overrides call `_recalculate_fit`, persist the new arrays under `/Analysis/Fits`, and mark the event so future sessions pick up the change.【F:IDEX-quicklook.py†L2705-L2831】
4. **Assign mass axes and baselines.** Baseline-subtracted fits combine the linear and sinusoidal estimates from `FitTargetSignal` with the physical grid model, so toggling baseline overlays clarifies how the detector responds relative to a zero-charge reference.【F:idex_packet.py†L188-L276】【F:IDEX-quicklook.py†L2323-L2670】 The shared `time2mass` conversion utilities can be invoked from the analysis notebook layer to relabel axes once the trigger time is confirmed.【F:time2mass.py†L1-L120】
5. **Guess compositions.** With cleaned fits visible, open the dust composition tools or export CSVs for offline analysis. Composition helpers (see `dust_composition.py`) read the same fitted amplitudes and time-of-flight metrics, so the workflow from Quicklook to composition inference remains consistent.【F:dust_composition.py†L1260-L1300】
6. **Document findings.** Capture screenshots or attach annotated plots via the export menu. When writing mission logs, embed the figures in Markdown using the approach from Section 4 so the documentation center and GitHub present identical views.【F:IDEX-quicklook.py†L1559-L1659】

Following this pipeline keeps Level-0 decoding, fit refinement, and interpretation aligned across command-line utilities, GUI tooling, and downstream documentation.
