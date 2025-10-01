# SpectrumPY-Flight

SpectrumPY-Flight is a LASP/IMPACT toolkit used during the IDEX flight-model integration and test campaign to ingest ground-system telemetry, decode packets, fit science waveforms, and generate quicklook data products. The repository bundles Python and shell utilities for live acquisition, Rice decompression, curve fitting, HDF5/CDF export, oscilloscope calibration processing, map generation, and GUI inspection workflows.

## Getting started

### Prerequisites

* Python 3.10 (the historic flight environment used Anaconda-based 3.10.4 interpreters).【F:science_tool.py†L1-L58】【F:lmfit_idex_packet.py†L1-L44】
* System packages for Qt-based GUIs (PyQt6) and scientific Python stack (NumPy, SciPy, matplotlib, pandas, seaborn, h5py, lmfit, bitstring, cdflib).【F:science_tool.py†L24-L58】【F:lmfit_idex_packet.py†L12-L44】【F:HDF_Plotter.py†L1-L158】【F:HDF_View.py†L60-L163】
* Mission-specific packet definition file `idex_combined_science_definition.xml` colocated in the repository root.【F:science_tool.py†L63-L79】【F:idex_packet.py†L205-L222】

### Environment setup

1. Create and activate a Python environment with the dependencies above (`conda create -n spectrum_py python=3.10`, then `pip install bitstring numpy scipy matplotlib pandas seaborn h5py lmfit cdflib pyqt6`).
2. Ensure the repository root is on your `PYTHONPATH` so helper modules (e.g., `time2mass`, `rice_decode`) resolve.
3. Create `Data/`, `HDF5/`, `Plots/`, `IonGridFits/`, `TargetFits/`, and `QDFits/` folders if you plan to run the automation scripts outside the lab environment (the code generates experiment-specific subfolders inside these directories).【F:ImpactBook.py†L124-L136】【F:combine_target_signals.py†L30-L144】

### Quickstart workflow

1. **Capture telemetry** from the GS-OASIS socket using `read_from_ois.py`, which writes raw binary dumps into `output_from_ois/` when interrupted (`Ctrl+C`).【F:read_from_ois.py†L1-L43】
2. **Decode packets** by feeding the raw file into `lmfit_idex_packet.py` or `science_tool.py` to build HDF5 products with waveform and analysis groups.【F:lmfit_idex_packet.py†L46-L200】【F:idex_packet.py†L203-L376】
3. **Inspect results** through the PyQt GUIs (`IDEX-quicklook.py`, `HDF_View.py`, `HDF_Plotter.py`, `Scope-IDEX-quicklook.py`) or produce tabular CSV summaries with `combine_target_signals.py` and analytics scripts such as `SNR_Calculator.py` and `histograms_maps.py`.【F:IDEX-quicklook.py†L66-L1588】【F:HDF_View.py†L60-L398】【F:HDF_Plotter.py†L24-L406】【F:combine_target_signals.py†L30-L149】【F:SNR_Calculator.py†L1-L81】【F:histograms_maps.py†L63-L303】
4. **Automate batch processing** with the provided shell helpers (`process_packets.sh`, `run_science.sh`, `process_EM_Data.sh`, `analyze_velocities.sh`).【F:process_packets.sh†L1-L27】【F:run_science.sh†L1-L8】【F:process_EM_Data.sh†L1-L37】【F:analyze_velocities.sh†L1-L18】

## Repository layout

| Path | Purpose |
| ---- | ------- |
| `Data/` | Raw telemetry staging area scanned by `process_packets.sh` before decoding. |
| `HDF5/` | Output folder for decoded events, quicklook exports, and GUI inputs. |
| `Plots/` | Destination for generated figures (e.g., quicklook PNGs, histogram maps). |
| `CDF/`, `HDF/`, `HDF5`, `CDF Variable Definitions.xlsx` | Reference data products and metadata. |
| `mass_comb.csv` | Mass comb definition consumed by `time2mass.py`. |
| `idex_combined_science_definition.xml` | XTCE definition for LASP packet parser. |
| Shell scripts (`*.sh`) | Batch automation for science processing and calibration analysis. |
| Python scripts (`*.py`) | Telemetry ingest, decoding, fitting, visualization, and analytics modules (detailed below). |

## Python script reference

### Telemetry ingest and packet decoding

#### `read_from_ois.py`
* Connects to the GS-OASIS SSIM socket (port 7514), continuously receives 1 KiB chunks, and appends raw bytes until interrupted.【F:read_from_ois.py†L16-L36】
* On `KeyboardInterrupt`, timestamps the capture, stores it in `output_from_ois/`, and can forward the file to downstream decoders (commented hooks for baseline runs).【F:read_from_ois.py†L30-L43】

#### `science_tool.py`
* Provides the historical “Science Tool” end-to-end decoder. The `IDEXEvent` class loads the XTCE definition, iterates packets, extracts event headers, handles Rice-compressed waveform segments, and populates `self.data`/`self.header` dictionaries for each event.【F:science_tool.py†L60-L203】
* Header extraction includes trigger delays, LVPS/HVPS telemetry words, and timestamp conversion to UTC epochs for each event.【F:science_tool.py†L104-L203】
* `plot_all_data` renders six-panel quicklook plots of TOF, target, and ion-grid channels with metadata overlays and writes PNGs per event.【F:science_tool.py†L241-L346】
* `write_to_hdf5` archives waveforms and derived time arrays into per-event HDF5 groups (including metadata epochs).【F:science_tool.py†L350-L377】
* `parse_hs_waveform`, `parse_ls_waveform`, and `parse_waveform_data` decode 10-bit high-sampling and 12-bit low-sampling bitstreams using `bitstring`, stripping padding and returning integer arrays keyed by `SCITYPE`.【F:science_tool.py†L383-L416】

#### `idex_packet.py`
* Houses a richer packet object with SciPy fitting primitives. Global utilities implement the IDEX ion-grid model (`IDEXIonGrid`), exponentially modified Gaussian (`EMG`), EMG area integrator, generic polynomial evaluator, dataset creation helper, EMG fit routine, and target-signal fit pipeline with baseline detrending and sine suppression.【F:idex_packet.py†L34-L201】
* The `IDEXEvent` class mirrors the science tool but stores every telemetry word, translates HVPS/LVPS house-keeping via lookup tables, and writes waveform/metadata datasets into HDF5 with paired time axes.【F:idex_packet.py†L203-L376】
* Waveform parsing and CDF export helpers at the bottom mirror the science tool implementations.【F:idex_packet.py†L968-L1008】

#### `lmfit_idex_packet.py`
* Extends `idex_packet.py` by swapping SciPy curve fits for `lmfit` models, enabling parameter bounds and residual minimization. The shared helpers (`IDEXIonGrid`, `EMG`, `calculate_area_under_emg`, `apply_polynomial`, `create_dataset_if_not_exists`, `FitEMG`, `FitTargetSignal`) are rewritten with lmfit parameter objects and low-pass filtered residual handling.【F:lmfit_idex_packet.py†L46-L200】
* The event parser builds HDF5 science and analysis groups, writes fitted parameters, target/ion-grid fit diagnostics, TOF-derived mass scales, and optional CDF outputs. The class exposes methods to decompress Rice packets, fit signals per event, compute statistics, and persist results similar to the SciPy version (see file body for the complete workflow).【F:lmfit_idex_packet.py†L201-L1074】
* Duplicated waveform parsing helpers ensure compatibility with `science_tool.py`.【F:lmfit_idex_packet.py†L1034-L1074】

#### `rice_decode.py`
* Implements Corrine’s Rice decompression algorithm for both CLI use (`rice_Decode`) and programmatic ingestion (`idex_rice_Decode`). The functions read ASCII hex streams, rebuild byte arrays, step through predictor-select states, Golomb-decode samples, and reconstruct 10-bit/12-bit frames with overflow checks and compression statistics.【F:rice_decode.py†L16-L157】

#### `bitstream.py` and `david_rice_decode.py`
* Provide a minimal byte-oriented bitstream reader with `read`, `read_byte`, `read_signed`, end-of-stream detection, pointer accessors, and flush/close semantics reused by legacy Rice decoders.【F:bitstream.py†L1-L88】【F:david_rice_decode.py†L1-L117】
* `david_rice_decode.py` also exposes a command-line `rice_Decode` entry point that wraps the bitstream utility for offline decoding. 【F:david_rice_decode.py†L123-L133】

#### `time2mass.py`
* Supplies lookup-based conversions from TOF sample index to mass, using `peak_time2mass` to align measured peaks with reference combs and `time2mass` to compute stretch/shift parameters alongside the calibrated mass scale. 【F:time2mass.py†L42-L257】

### GUI quicklook and plotting tools

#### `IDEX-quicklook.py`
* Qt application for browsing HDF5 science products. Utility functions handle file dialogs, dataset access, name normalization, LaTeX math rendering, dataset sorting, and numeric coercion (`non_native_open_dialog`, `list_event_groups`, `dset`, `_normalize_name`, `_pair_key`, `_friendly_label`, `_label_from_param_path`, `_latex_to_pixmap`, `_mass_identifier_from_path`, `_dataset_sort_key`, `_coerce_parameter_values`, `_to_1d`).【F:IDEX-quicklook.py†L66-L243】
* Analytical helpers produce fit curves for ion-grid (`_idex_ion_grid_model`) and EMG waveforms (`_emg_model`), combine time arrays with parameter vectors (`_evaluate_fit_curve`, `_masked_mean`, `_fit_paths_from_param`), and label axes with units (`y_label_with_units`).【F:IDEX-quicklook.py†L304-L462】
* `FitData` dataclass-like container (via methods on `FitData` inner class) supplies overlay/parameter iteration (`has_overlay`, `has_parameters`, `iter_time_result_pairs`, `iter_parameter_items`).【F:IDEX-quicklook.py†L419-L448】
* `MainWindow` orchestrates the UI: constructor loads files/events, `_build_toolbar` and `_build_controls` create Qt widgets, `action_*` slots open files or auxiliary viewers, `open_file`/`reload_current` manage the active dataset, and event/channel toggles update the multi-axis matplotlib canvas (`on_event_changed`, `on_channel_toggled`, `on_overlay_toggled`, `on_fit_toggled`, `plot_event`). Plotting helpers (`_plot_channel`, `_estimate_baseline`, `_iter_fit_curves`, `_plot_fit`, `_style_overlay_axis`, `_style_single_axis`, `_draw_missing_message`) render raw, overlay, and fit data, while `refresh_fit_controls` and `update_status_text` synchronize controls and status bars.【F:IDEX-quicklook.py†L460-L1079】
* Fit management routines (`get_fit_data`, `open_fit_parameter_dialog`, `get_parameter_values`, `update_fit_override`, `clear_fit_override`, `_recalculate_fit`, `_remove_fit_override`) let analysts tweak parameters in a dedicated dialog class (`FitParameterDialog`), whose methods manage channel selection, dataset lists, value formatting, feedback, table editing, and apply/reset operations.【F:IDEX-quicklook.py†L1102-L1577】
* `main()` wires up `QApplication`, parses optional CLI arguments (`--filename`, `--eventnumber`), and launches the window.【F:IDEX-quicklook.py†L1585-L1595】

#### `HDF_View.py`
* Lightweight HDF5 browser: `_format_scalar`, `_format_attribute`, and `_shape_to_text` convert NumPy/scalar values for display. `HDFViewWindow` loads files, builds a tree/table UI (`_build_ui`, `_populate_tree`, `_add_children`), reacts to selections (`_on_item_selected`), and populates summary/data/attribute tabs (`_show_summary`, `_show_dataset`, `_show_group`, `_populate_attr_table`). File dialog helpers (`_choose_file_dialog`) and `launch_hdf_viewer` integrate with other tools, while `main()` provides a CLI entry point.【F:HDF_View.py†L60-L398】

#### `HDF_Plotter.py`
* Dashboard for multi-file scatter matrix exploration. `HDF5PlotterApp` constructor loads all `Analysis` datasets from an HDF5 folder, storing values with source-file provenance (`load_datasets`). `init_ui` builds dual-plot selectors and matplotlib canvases, and `plot_data` renders stacked scatter plots with interactive pick/hover callbacks that open the detailed quicklook window for selected points.【F:HDF_Plotter.py†L24-L406】
* Main guard bootstraps a `QApplication` against the `HDF5/` directory.【F:HDF_Plotter.py†L398-L410】

#### `Scope-IDEX-quicklook.py`
* PyQt GUI tailored for oscilloscope `.trc` campaigns. Custom formatters (`MyFormatter`) and plot canvases (`MplCanvas`) support navigation. `MainWindow` wires menu actions, TRC import/export, SQL logging, channel selection, waveform fitting (`fitTarget`, `fitIon`, `fitQD`), and quick navigation across traces (`upTrace`, `downTrace`, `chooseTrace`, `changeChannel`). Auxiliary widgets (`SQLWindow`, `ChooseTrace`) provide dialogs for SQL metadata and trace selection. Helper functions `read_hdf` and `displayTRC` load HDF5 quicklook exports and display TRC arrays directly.【F:Scope-IDEX-quicklook.py†L114-L773】

#### `ImpactBook.py`
* End-to-end oscilloscope processing pipeline. Global fits (`LinearFit`, `SineFit`, `QDFit`, `IDEXIonGrid`) capture background trends, pickup-tube physics, and ion-grid dynamics.【F:ImpactBook.py†L37-L118】
* `ImpactBook` orchestrates experiment ingestion: runs `generalReadTRC`, initializes output folders, and instantiates `ImpactEvent` for each trace.【F:ImpactBook.py†L120-L143】
* `ImpactEvent` builds pandas data frames, writes per-channel HDF5 datasets, performs target, ion-grid, and QD fits, derives amplitudes/timings, fetches accelerator metadata via SQL, and saves plots/CSV/HDF5 products (`FitTargetSignal`, `FitIonGridSignal`, `FitQDSignal`, `RetrieveSQL`).【F:ImpactBook.py†L145-L609】
* Module-level utilities `generalReadTRC`, `butter_lowpass_filter`, and `find_nearest` parse LeCroy TRC directories, denoise signals, and provide nearest-neighbor helpers.【F:ImpactBook.py†L636-L771】

#### `qd_quicklook.py`
* Streamlined version of `ImpactBook` for quick target/QD assessments. Shares the same fitting helpers and ImpactEvent workflow but includes optimized TRC reading (`optimizedReadTRC`), dual SQL retrieval paths, and additional smoothing utilities (`butter_lowpass_filter`, `find_nearest`).【F:qd_quicklook.py†L43-L878】

#### `readTrc.py`
* Implements a `Trc` class capable of opening LeCroy `.trc` files, parsing the binary header (`_readX`, `_readS`), sampling arrays (`_readSamples`), and constructing timestamps (`_getTimeStamp`).【F:readTrc.py†L43-L223】

### Analysis utilities

#### `combine_target_signals.py`
* Iterates over HDF5 products (`iter_h5_files`), retrieves scalar analysis values by dataset name (`_get_dataset_by_name`, `_strict_scalar_from_ds`), logs debug mismatches (`_dbg_print`), and writes consolidated CSV tables (`write_csv`) enumerating `LImpactCharge` metrics per event. CLI `main()` wires argument parsing and orchestrates extraction.【F:combine_target_signals.py†L30-L149】

#### `SNR_Calculator.py`
* `extract_snr` scans each event’s high-sampling TOF waveform, derives baseline statistics over a pre-trigger window, converts time to mass using `time2mass`, extracts peak prominence with `find_peaks`, computes stretch/shift/kappa diagnostics, and generates scatter plus histogram visualizations for SNR studies.【F:SNR_Calculator.py†L1-L103】

#### `histograms_maps.py`
* Generates synthetic longitude/latitude counts maps. Helpers wrap angles, build histogram edges, compute angular separations, synthesize Gaussian noise (`wrap180`, `centers_to_edges`, `angsep_deg`, `gaussian`, `smooth_noise_on_lon`). `build_idp_intensity` and `effective_lat_std_deg` craft intensity distributions. `render_map` renders Mollweide projections with colorbars, and `fmt` formats fractional labels for tick annotations.【F:histograms_maps.py†L63-L303】

#### `run_all.py`
* Convenience launcher that sequentially runs `science_tool.py` and a backup routine via `subprocess` (`run_science_tool`, `run_backup_script`).【F:run_all.py†L1-L53】

### Shell automation

* `process_packets.sh` walks the `Data/` tree for extensionless telemetry captures, skips already-decoded files, and invokes `lmfit_idex_packet.py` for each new dataset, writing HDF5 into `HDF5/`.【F:process_packets.sh†L1-L27】
* `run_science.sh` mirrors `run_all.py` by executing `science_tool.py` followed by `backup_science.py`.【F:run_science.sh†L1-L8】
* `process_EM_Data.sh` recursively processes calibration directories with `.trc` files, launching `ImpactBook.py` per experiment and logging failures.【F:process_EM_Data.sh†L1-L37】
* `analyze_velocities.sh` loops over accelerator campaign subdirectories and runs `qd_quicklook.py` with directory-specific titles to batch-generate analyses.【F:analyze_velocities.sh†L1-L18】

## Data products and downstream outputs

* HDF5 archives contain per-event waveform datasets (`TOF`, `Target`, `Ion Grid`, `Time (high sampling)`, `Time (low sampling)`), metadata groups (epoch timestamps), and analysis groups with fit parameters/curves when using the lmfit pipeline.【F:science_tool.py†L350-L377】【F:lmfit_idex_packet.py†L201-L1074】
* Quicklook GUIs support overlaying fit results, exporting plots, and opening HDF5 summaries for collaboration via Qt dialogs.【F:IDEX-quicklook.py†L460-L1588】【F:HDF_Plotter.py†L24-L406】【F:HDF_View.py†L60-L398】
* CSV consolidations (e.g., `combine_target_signals.py`) and SNR studies produce mission-readiness metrics for mass/charge calibration and detection performance.【F:combine_target_signals.py†L30-L149】【F:SNR_Calculator.py†L1-L103】

## Operational notes

* The scripts assume the working directory is the repository root so relative paths to `HDF5/`, `Data/`, and XML definitions resolve correctly.【F:science_tool.py†L63-L377】【F:process_packets.sh†L1-L27】
* Rice decompression requires ASCII-hex input; ensure raw telemetry is converted before running `idex_rice_Decode` directly. For automated workflows, rely on `lmfit_idex_packet.py` or `science_tool.py`, which invoke the decoder internally.【F:lmfit_idex_packet.py†L41-L200】【F:rice_decode.py†L16-L157】
* GUI applications require a display server; on headless systems use `xvfb-run` or similar wrappers when launching PyQt-based tools.【F:IDEX-quicklook.py†L460-L1595】【F:HDF_View.py†L60-L398】【F:HDF_Plotter.py†L24-L410】【F:Scope-IDEX-quicklook.py†L114-L773】

