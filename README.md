# SpectrumPY-Flight

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14948734.svg)](https://doi.org/10.5281/zenodo.14948734)

> Ground-system decoding, fitting, and visualization tools for the LASP/IMPACT IDEX campaign.

![SpectrumPY overview](src/spectrumpy_flight/docs/media/quicklook_overview.svg)

![Quicklook windows](src/spectrumpy_flight/docs/media/quicklook_windows.svg)

### Quicklook highlights

* **Auto ingest** – Drag in `.trc` scope captures or browse to existing telemetry; the launcher hands the file off to the right viewer without extra CLI flags.【F:src/spectrumpy_flight/spectrum_launcher.py†L196-L279】
* **Batch processing** – Shell helpers wrap `idex_packet.py` and driver scripts so you can queue directories of packets for decoding while the GUI stays available for inspection.【F:process_packets.sh†L1-L27】【F:drive_idex_packet.py†L24-L134】
* **Fit adjustment** – Reload events and tweak EMG parameters directly inside Quicklook; overrides persist to disk so later sessions start from your edits.【F:IDEX-quicklook.py†L1056-L1349】【F:dust_composition.py†L1129-L1478】
* **TOF mass spectrometry** – Combine high-/mid-/low-gain traces, fit analytic line shapes, and export abundance tables from the Dust Composition workspace.【F:dust_composition.py†L3008-L3525】

---

## Overview

SpectrumPY-Flight wraps packet decoding, fitting, and visualization into a
single workflow that covers the full time-of-flight (TOF) dust analysis loop:

* **Automatic TOF mass spectrometry** – Load raw telemetry or oscilloscope
  captures and generate calibrated mass axes with the Newton-stabilized series
  expansion implemented in `mass_calibration.py` and `time2mass.py`.
  Interactive stretch-and-shift bootstrapping keeps solutions stable even on
  faint spectra.【F:src/spectrumpy_flight/docs/time_to_mass_calibration.md†L1-L111】
* **Batch processing** – Launch the provided shell helpers to ingest new
  downlinks or laboratory campaigns without touching the GUI, then hand the
  derived products back to the quicklook tools for inspection.【F:process_packets.sh†L1-L27】【F:process_EM_Data.sh†L1-L37】
* **Fit adjustments on demand** – Recompute exponential-modified Gaussian
  (EMG) overlays, tweak lmfit parameters, and persist manual overrides directly
  inside the Dust Composition workspace.【F:IDEX-quicklook.py†L1056-L1349】【F:dust_composition.py†L1129-L1478】
* **Composition identification** – Compare detected species to reference
  libraries, apply relative sensitivity factors, and summarise abundances with
  ternary diagrams that highlight cometary dust classes and olivine content in
  real time.【F:dust_composition.py†L3008-L3484】【F:dust_composition.py†L3526-L4023】

---

## Launch in 10 seconds

```
python -m spectrumpy_flight.start
```

The Spectrum Launcher opens immediately and hands off to the Quicklook viewer so you can load telemetry, browse documentation, and export science plots without touching any other scripts.【F:src/spectrumpy_flight/start.py†L1-L7】【F:src/spectrumpy_flight/spectrum_launcher.py†L167-L289】

## Installation

Install the published package from PyPI when you want to run SpectrumPY outside
the repository checkout:

```
python -m pip install spectrumpy
```

The default install remains lightweight so packet analysis scripts such as
`idex_packet.py` and `drive_idex_packet.py` work on headless hosts without
Qt. Add the GUI and SQL helpers in a single step when you want to run the
Quicklook tooling:

```
python -m pip install "spectrumpy[quicklook]"
```

macOS-specific preparation steps are covered in
[`src/spectrumpy_flight/docs/macos_pip_install.md`](src/spectrumpy_flight/docs/macos_pip_install.md). The document explains
how to prepare Homebrew Python environments, optional extras, and Qt security
prompts for first-time launches on Ventura and newer systems.【F:src/spectrumpy_flight/docs/macos_pip_install.md†L1-L83】

Optional extras match the bundled feature sets declared in `pyproject.toml`.
Combine extras when needed (`spectrumpy[quicklook,pyqt,gpu,database]`) to pull in
all secondary dependencies listed in the metadata.【F:pyproject.toml†L32-L52】

| Extra | Included dependencies | When to use |
| --- | --- | --- |
| `quicklook` | PySide6 + qtawesome for the Qt GUI, SQLAlchemy + PyMySQL for accelerator matching | Enables the Spectrum Launcher, Quicklook viewer, and other GUI entry points on a clean environment. |
| `pyqt` | PyQt6 | Swap to the PyQt bindings instead of PySide6 (can be combined with `quicklook`). |
| `gpu` | CuPy | Accelerate packet decoding on CUDA-capable hosts. |
| `database` | SQLAlchemy, PyMySQL, MySQL connector | Install database adapters without the GUI stack (used for legacy telemetry archives). |

SpectrumPY ships lean wheels so `pip install spectrumpy` stays well under
PyPI's size limits. Sample telemetry captures and pre-generated CDF products
remain available in the source repository; download them separately when you
need tutorial data.

### Command-line entry points

The PyPI package exposes convenience launchers so the most common workflows
require no manual path management after installation:

* `idex-packet` – decode raw packets and emit HDF5 analysis products.
* `drive-idex-packet` – orchestrate parallel conversions across directories of
  captures.
* `idex-quicklook` – open the Qt quicklook viewer directly.
* `spectrumpy-start` – launch the SpectrumPY welcome screen with automatic
  oscilloscope ingestion.

---

## Core capabilities

* **Packet decommutation** – Decode binary telemetry with the XTCE mission definition and publish structured HDF5/CDF products for analysis.【F:idex_packet.py†L203-L440】【F:science_tool.py†L60-L377】
* **Batch processing** – Automate pipeline runs for new downlinks or oscilloscope campaigns with the provided shell helpers.【F:process_packets.sh†L1-L27】【F:process_EM_Data.sh†L1-L37】
* **Fit adjustment** – Refine lmfit overlays inside the GUI, reset overrides, persist overrides to disk, and compare solutions in real time.【F:IDEX-quicklook.py†L1056-L1349】【F:IDEX-quicklook.py†L3520-L3722】
* **Noise diagnostics** – Launch the dedicated analysis window for Gaussian histograms, FFT spectra, and autocorrelation studies per channel.【F:IDEX-quicklook.py†L2027-L2188】【F:noise_analysis.py†L260-L429】
* **Parallel decoding** – Scale packet processing across CPU cores with automatic worker discovery, per-job memory limits, and live progress reporting.【F:tools/process_packets_parallel.py†L1-L120】【F:tools/resource_utils.py†L1-L75】
* **TOF mass spectrometry** – Quantify dust populations with the composition workspace, apply relative sensitivity factors, match library samples, and inspect EMG mass lines in detail.【F:dust_composition.py†L3008-L3484】【F:dust_composition.py†L3526-L4023】

### Time-to-mass conversion details

The calibration workflow models the arrival-time series as a polynomial in the
reduced coordinate `s = sqrt(m)`:

```
t(s) = t_0 + a_1 s + a_2 s^2 + \cdots + a_K s^K.
```

Least-squares fitting recovers the coefficients `a_k` from calibration lines
and the safeguarded Newton inversion maps measured times back to masses while
clipping to the instrument range.【F:src/spectrumpy_flight/docs/time_to_mass_calibration.md†L13-L91】

After convergence the squared coordinate yields the final mass axis:

```
m = (s^(n))^2.
```

These steps are exposed through `TOFMassCal` and the `time2mass.py` CLI so that
both automated batch runs and interactive GUI sessions share the same
well-conditioned calibration pipeline.【F:src/spectrumpy_flight/docs/time_to_mass_calibration.md†L1-L146】

---

## Performance and parallelization details

### idex_packet.py

* **GPU-assisted bit unpacking** – `_bitstring_to_ints` converts ASCII bit streams into integers with a fully vectorised NumPy path and an automatic CuPy offload whenever more than 32 000 values are present, falling back seamlessly if GPU execution is unavailable or fails mid-run.【F:idex_packet.py†L96-L125】 The helper keeps allocations contiguous (`astype(..., copy=False)`) and trims padding without extra passes so the downstream waveform parsing never re-materialises the raw bitstring.【F:idex_packet.py†L96-L125】
* **Shared kernels for saturation and SNR detection** – Signal-quality metrics (SNR, saturation runs) operate on NumPy arrays with optional GPU execution that mirrors the CPU algorithm: `detect_saturation` lifts the differencing and run-length detection to CuPy when possible, then reuses the same boolean boundaries to avoid Python loops, while `calculate_snr` slices baseline windows in-place to minimise temporary buffers.【F:idex_packet.py†L128-L195】
* **Rice/Golomb decompression in-stream** – RiceGolombDecompressor instances emit fully decompressed traces directly into `self.data`, eliminating intermediate files for both low- and high-speed telemetry and reusing instrument trigger block counts to size the buffers exactly.【F:idex_packet.py†L941-L979】
* **Time-base reuse and deterministic scaling** – `_build_time_array` precomputes high-rate and low-rate time axes once per batch using cached trigger offsets, so waveform analyses (mass fitting, target modelling) only rescale existing arrays instead of rebuilding per channel.【F:idex_packet.py†L1002-L1290】
* **Threaded waveform analysis with capped fan-out** – `write_to_hdf5` packages every `(event, channel)` pair and submits them to a `ThreadPoolExecutor` whose worker count is bounded by both the host CPU count and a hard cap of 32, preventing oversubscription when many events are present.【F:idex_packet.py†L1182-L1295】 Each worker returns a structured dictionary that the coordinator consumes sequentially to materialise datasets, guaranteeing that HDF5 writes remain serial while expensive per-channel transforms (mass scaling, EMG fitting, target curve estimation) execute in parallel.【F:idex_packet.py†L1208-L1336】

### drive_idex_packet.py

* **Parallel process orchestration** – The driver discovers extensionless oscilloscope captures, filters out already-converted targets, and farms the remaining workload to a thread-backed execution pool that launches separate `idex_packet.py` processes so every core handles an independent capture.【F:drive_idex_packet.py†L24-L134】 Logging is performed per job (`job_XXXX.(out|err)`), enabling high-velocity runs without interleaved stdout noise.【F:drive_idex_packet.py†L51-L71】
* **Oversubscription controls** – Operator-provided `--max-procs` bounds the executor fan-out, while `build_env` injects BLAS/OpenMP thread caps (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, and Accelerate’s `VECLIB_MAXIMUM_THREADS`) so the heavy SciPy/Numpy kernels inside each packet conversion cannot oversubscribe cores when dozens of jobs launch at once.【F:drive_idex_packet.py†L39-L134】 Optional niceness further de-prioritises the child processes on POSIX hosts for cooperative scheduling in shared environments.【F:drive_idex_packet.py†L60-L71】
* **Resumable workflows** – By evaluating destination `.h5` existence up front (`needs_conversion`), the orchestrator skips already processed captures, enabling incremental reruns that immediately converge on the remaining work without repeated scans or redundant CPU time.【F:drive_idex_packet.py†L35-L110】

---

## Quicklook essentials

### Primary windows

* **Data Browser** – Press `Ctrl+B` or choose **View → Open Data Browser** to launch the structure tree for the active event. The viewer exposes every dataset, attribute, and preview panel inside the HDF5/CDF product.【F:IDEX-quicklook.py†L1888-L1964】【F:HDF_View.py†L60-L398】
* **Variable Definitions** – Use **View → Variable Definitions…** to open the spreadsheet companion that maps engineering units, calibration spans, and coefficients straight from the official Excel reference.【F:IDEX-quicklook.py†L1888-L1964】【F:IDEX_Definitions_View.py†L1-L199】
* **Dust Composition** – Click the **Dust Composition** toolbar button or press `Ctrl+D` to launch the EMG fitting studio for the current event.【F:IDEX-quicklook.py†L2036-L2083】【F:dust_composition.py†L1-L215】
* **Noise Analysis** – Use the toolbar button or `Ctrl+N` to open the diagnostics workspace with histogram, FFT, and autocorrelation summaries for the active channel.【F:IDEX-quicklook.py†L2027-L2188】【F:noise_analysis.py†L260-L429】
* **Documentation Center** – Hit `F1` (or the `?` button) to browse the bundled Markdown library without leaving the viewer; `Ctrl+F1` jumps directly to the global search panel.【F:IDEX-quicklook.py†L1888-L1964】

### Inspect individual mass lines

1. Open the Dust Composition window (`Ctrl+D`) from the main toolbar.【F:IDEX-quicklook.py†L2036-L2083】
2. Combine the TOF traces if needed, select a region on the stacked plot, and click **Add Mass Line** to fit a new EMG peak; manual entry is also available through **Add Manual Line**.【F:dust_composition.py†L1129-L1244】
3. Highlight a row in the **Mass Line Fits** table and press **Inspect Selected** to open the zoomed EMG editor with amplitude, μ, σ, λ, and window controls, plus optional HyperEMG/Voigt tails.【F:dust_composition.py†L1188-L1478】
4. Apply relative sensitivity factors (RSFs) or sample library templates to renormalise abundances, adjust parameters, and confirm—tables and ternary summaries refresh immediately.【F:dust_composition.py†L3008-L3484】【F:dust_composition.py†L3526-L4023】

---

## Packaging quick reference

### Shared setup

All platforms use the same PyInstaller spec. Start from a clean Python ≥3.10 environment and install the packaging requirements:

```
python -m pip install --upgrade pip
pip install -r packaging/packaging-requirements.txt
```

Then build the application bundle:

```
pyinstaller --noconfirm packaging/idex_quicklook.spec
```

【F:packaging/README.md†L26-L55】【F:packaging/idex_quicklook.spec†L1-L89】

### macOS

```
./packaging/macos/create_dmg.sh dist/IDEX-Quicklook.app "IDEX-Quicklook-$(git describe --tags --always)-$(uname -m).dmg"
```

The script wraps the generated `.app` into a signed-ready disk image; finish by running your usual `codesign` and notarization steps.【F:packaging/README.md†L57-L76】

### Windows

PyInstaller emits `dist/IDEX-Quicklook/IDEX-Quicklook.exe`. Compress that folder into a ZIP for delivery, or feed it into your MSI tool of choice, for example:

```
powershell -Command "Compress-Archive -Path dist/IDEX-Quicklook -DestinationPath IDEX-Quicklook-windows.zip"
```

Users unzip the archive and launch `IDEX-Quicklook.exe` directly.【F:packaging/README.md†L78-L88】

### Linux

Package the frozen binary as a tarball ready for distribution:

```
tar -C dist -czf IDEX-Quicklook-linux.tar.gz IDEX-Quicklook
```

The executable runs without root privileges and extracts its shared libraries on first launch.【F:packaging/README.md†L90-L108】

---

## Documentation & tutorials

* **Quicklook tutorial** – Guided walkthrough of the GUI workflow, shortcuts, and recommended analysis order.【F:src/spectrumpy_flight/docs/quicklook_tutorial.md†L1-L123】
* **Fitting reference** – Deep dive into baseline handling, Data Browser correlations, and variable-definition lookups.【F:src/spectrumpy_flight/docs/fitting_reference.md†L1-L115】
* **Packaging playbook** – Extended release checklists, validation steps, and troubleshooting strategies for desktop builds.【F:src/spectrumpy_flight/docs/packaging_tutorial.md†L1-L198】
* **Olivine metrics** – Instructions for running the olivine regression and locating the generated SNR, timing, and saturation reports.【F:src/spectrumpy_flight/docs/olivine_metrics.md†L1-L64】
* **Time-to-mass calibration** – Polynomial fitting, Newton inversion safeguards, and quick recipes for turning TOF traces into calibrated mass axes.【F:src/spectrumpy_flight/docs/time_to_mass_calibration.md†L1-L146】
* **Feature reference** – Exhaustive catalogue of launchers, CLIs, shell helpers, and their underlying algorithms with LaTeX formulas and usage examples.【F:src/spectrumpy_flight/docs/feature_reference.md†L1-L195】

Launch the in-app documentation center anytime (`F1`) to search and read these guides inside the Quicklook viewer.【F:IDEX-quicklook.py†L1888-L1964】

<!-- policy:begin -->
## License

This repository is released under the **Ayari Public No-Derivatives License (APND) v1.0**.
You may **download and use the software as-is**. You may **not** modify it or distribute
modified versions without written permission. See [`LICENSE`](./LICENSE).

## Citation & DOI

If you use this software in research or a product, please **cite it**. We archive
releases on **Zenodo** to mint a DOI. You can also download machine-readable metadata from [`CITATION.cff`](./CITATION.cff).

**How to cite (example):**
> Ayari, E. (2025). SpectrumPY-Flight (v1.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.14948734

For convenience, you can embed the DOI badge anywhere documentation is published:

```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14948734.svg)](https://doi.org/10.5281/zenodo.14948734)
```

## Contributing

This is a **read-only** public release. We do **not** accept external pull requests or patches.
Please open an **issue** for bugs or feature requests. For derivative-use exceptions, contact the author.

## Repro/Use

- Clone or download a release tarball.
- Use **unmodified** sources/binaries as described in the docs.
- Do **not** redistribute modified versions.
<!-- policy:end -->
