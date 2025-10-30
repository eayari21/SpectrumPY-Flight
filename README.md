# SpectrumPY-Flight

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14948734.svg)](https://doi.org/10.5281/zenodo.14948734)

> Ground-system decoding, fitting, and visualization tools for the LASP/IMPACT IDEX campaign.

![SpectrumPY overview](src/spectrumpy_flight/docs/media/quicklook_overview.svg)

---

## Launch in 10 seconds

```
spectrumpy-launcher
```

The Spectrum Launcher opens immediately and hands off to the Quicklook viewer so you can load telemetry, browse documentation, and export science plots without touching any other scripts.【F:src/spectrumpy_flight/start.py†L1-L7】【F:src/spectrumpy_flight/spectrum_launcher.py†L167-L289】

---

## Core capabilities

* **Packet decommutation** – Decode binary telemetry with the XTCE mission definition and publish structured HDF5/CDF products for analysis.【F:src/spectrumpy_flight/idex_packet.py†L203-L440】【F:src/spectrumpy_flight/science_tool.py†L60-L377】
* **Batch processing** – Automate pipeline runs for new downlinks or oscilloscope campaigns with the provided shell helpers.【F:process_packets.sh†L1-L27】【F:process_EM_Data.sh†L1-L37】
* **Fit adjustment** – Refine lmfit overlays inside the GUI, reset overrides, persist overrides to disk, and compare solutions in real time.【F:src/spectrumpy_flight/idex_quicklook.py†L1056-L1349】【F:src/spectrumpy_flight/idex_quicklook.py†L3520-L3722】
* **Noise diagnostics** – Launch the dedicated analysis window for Gaussian histograms, FFT spectra, and autocorrelation studies per channel.【F:src/spectrumpy_flight/idex_quicklook.py†L2027-L2188】【F:src/spectrumpy_flight/noise_analysis.py†L260-L429】
* **Parallel decoding** – Scale packet processing across CPU cores with automatic worker discovery, per-job memory limits, and live progress reporting.【F:src/spectrumpy_flight/tools/process_packets_parallel.py†L1-L120】【F:src/spectrumpy_flight/tools/resource_utils.py†L1-L75】
* **TOF mass spectrometry** – Quantify dust populations with the composition workspace, apply relative sensitivity factors, match library samples, and inspect EMG mass lines in detail.【F:src/spectrumpy_flight/dust_composition.py†L3008-L3484】【F:src/spectrumpy_flight/dust_composition.py†L3526-L4023】

---

## Quicklook essentials

### Primary windows

* **Data Browser** – Press `Ctrl+B` or choose **View → Open Data Browser** to launch the structure tree for the active event. The viewer exposes every dataset, attribute, and preview panel inside the HDF5/CDF product.【F:src/spectrumpy_flight/idex_quicklook.py†L1888-L1964】【F:src/spectrumpy_flight/HDF_View.py†L60-L398】
* **Variable Definitions** – Use **View → Variable Definitions…** to open the spreadsheet companion that maps engineering units, calibration spans, and coefficients straight from the official Excel reference.【F:src/spectrumpy_flight/idex_quicklook.py†L1888-L1964】【F:src/spectrumpy_flight/IDEX_Definitions_View.py†L1-L199】
* **Dust Composition** – Click the **Dust Composition** toolbar button or press `Ctrl+D` to launch the EMG fitting studio for the current event.【F:src/spectrumpy_flight/idex_quicklook.py†L2036-L2083】【F:src/spectrumpy_flight/dust_composition.py†L1-L215】
* **Noise Analysis** – Use the toolbar button or `Ctrl+N` to open the diagnostics workspace with histogram, FFT, and autocorrelation summaries for the active channel.【F:src/spectrumpy_flight/idex_quicklook.py†L2027-L2188】【F:src/spectrumpy_flight/noise_analysis.py†L260-L429】
* **Documentation Center** – Hit `F1` (or the `?` button) to browse the bundled Markdown library without leaving the viewer; `Ctrl+F1` jumps directly to the global search panel.【F:src/spectrumpy_flight/idex_quicklook.py†L1888-L1964】

### Inspect individual mass lines

1. Open the Dust Composition window (`Ctrl+D`) from the main toolbar.【F:src/spectrumpy_flight/idex_quicklook.py†L2036-L2083】
2. Combine the TOF traces if needed, select a region on the stacked plot, and click **Add Mass Line** to fit a new EMG peak; manual entry is also available through **Add Manual Line**.【F:src/spectrumpy_flight/dust_composition.py†L1129-L1244】
3. Highlight a row in the **Mass Line Fits** table and press **Inspect Selected** to open the zoomed EMG editor with amplitude, μ, σ, λ, and window controls, plus optional HyperEMG/Voigt tails.【F:src/spectrumpy_flight/dust_composition.py†L1188-L1478】
4. Apply relative sensitivity factors (RSFs) or sample library templates to renormalise abundances, adjust parameters, and confirm—tables and ternary summaries refresh immediately.【F:src/spectrumpy_flight/dust_composition.py†L3008-L3484】【F:src/spectrumpy_flight/dust_composition.py†L3526-L4023】

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

## Publishing a pip release

Follow this checklist to cut a new wheel and publish it to PyPI from macOS. All commands run in a terminal from the repository root after activating your virtual environment.

1. **Update metadata.** Set the new semantic version in `pyproject.toml` and commit documentation updates that describe the release scope.【F:pyproject.toml†L6-L35】
2. **Install build tooling.**
   ```bash
   python -m pip install --upgrade pip build twine
   ```
3. **Build source and wheel artifacts.**
   ```bash
   python -m build
   ```
   The distributions land in `dist/` as `spectrumpy_flight-<version>.tar.gz` and `.whl` archives.
4. **Smoke-test the wheel.** Create a throwaway virtual environment, install the freshly built wheel, and run the CLI entry points (`spectrumpy-launcher`, `spectrumpy-quicklook`, etc.) to confirm Qt launches and bundled data loads correctly.【F:pyproject.toml†L38-L44】
5. **Upload to PyPI.** If the tests pass, publish the artifacts with Twine:
   ```bash
   python -m twine upload dist/*
   ```
   For TestPyPI dry runs, append `--repository testpypi` to the command.
6. **Tag the release.** Record the version in git so the DOI archive and downstream users can track the exact build:
   ```bash
   git tag -a v<version> -m "Release v<version>"
   git push --tags
   ```

---

## Documentation & tutorials

* **Quicklook tutorial** – Guided walkthrough of the GUI workflow, shortcuts, and recommended analysis order.【F:src/spectrumpy_flight/docs/quicklook_tutorial.md†L1-L123】
* **Fitting reference** – Deep dive into baseline handling, Data Browser correlations, and variable-definition lookups.【F:src/spectrumpy_flight/docs/fitting_reference.md†L1-L115】
* **Packaging playbook** – Extended release checklists, validation steps, and troubleshooting strategies for desktop builds.【F:src/spectrumpy_flight/docs/packaging_tutorial.md†L1-L198】
* **Olivine metrics** – Instructions for running the olivine regression and locating the generated SNR, timing, and saturation reports.【F:src/spectrumpy_flight/docs/olivine_metrics.md†L1-L64】

Launch the in-app documentation center anytime (`F1`) to search and read these guides inside the Quicklook viewer.【F:src/spectrumpy_flight/idex_quicklook.py†L1888-L1964】

<!-- policy:begin -->
## License

This repository is released under the **Ayari Public No-Derivatives License (APND) v1.0**.
You may **download and use the software as-is**. You may **not** modify it or distribute
modified versions without written permission. See [`LICENSE`](./LICENSE).

## Citation & DOI

If you use this software in research or a product, please **cite it**. We archive
releases on **Zenodo** to mint a DOI. You can also download machine-readable metadata from [`CITATION.cff`](./CITATION.cff).

**How to cite (example):**
> Ayari, E. (2025). SpectrumPY-Flight (v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.14948734

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
