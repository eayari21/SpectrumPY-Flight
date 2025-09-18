# SpectrumPY-Flight

SpectrumPY-Flight is a LASP/IMPACT toolkit built for the IDEX flight-model integration and test campaign. It bundles Python utilities that ingest, decode, analyze, and visualize telemetry from the GS-OASIS ground system.

## Flight telemetry ingest & decoding
- The primary `science_tool.py` script parses binary packets defined in `idex_combined_science_definition.xml`, using LASP’s XTCE parser to walk telemetered events, apply trigger logic, and harvest waveform payloads for subsequent analysis.
- Supporting modules such as `rice_decode.py` implement the mission’s Rice decompression algorithm, while `idex_packet.py` wraps parsed packets with signal-fitting routines (EMG, ion-grid, target fits) and exports results to CDF/HDF structures for downstream science products.
- Batch scripts automate telemetry handling: `process_packets.sh` scans acquisition directories, converts raw files with `lmfit_idex_packet.py`, and builds HDF5 outputs, and `bitstream.py` supplies low-level bitstream reading primitives used by the decoders.

## Data capture and orchestration
- `read_from_ois.py` listens to the OASIS socket server, captures live SSIM telemetry, and writes timestamped raw files that can be fed back into the decoding pipeline.
- Convenience launchers such as `run_all.py` and `run_science.sh` chain the main science tool with a backup routine, streamlining routine data runs during test operations.

## Visualization and quicklook tools
- PyQt-based applications (`IDEX-quicklook.py`, `HDF_Plotter.py`) provide interactive waveform and parameter inspection across HDF5 products, including multi-panel plotting, event selection, and hooks back into the quicklook window for detailed review.
- Additional GUIs such as `Scope-IDEX-quicklook.py` mirror the quicklook workflow for oscilloscope-derived data, emphasizing consistent interfaces across lab and flight contexts.

## Accelerator and scope data processing
- Laboratory waveforms from LeCroy oscilloscopes are ingested via `readTrc.py` and processed with `ImpactBook.py`/`qd_quicklook.py`, which fit pickup-tube, target, and ion-grid signals, archive results to HDF5/CSV, and visualize velocity distributions; automation scripts (`process_EM_Data.sh`, `analyze_velocities.sh`) batch these conversions over entire calibration campaigns.

## Analysis utilities
- Data-mining scripts include `combine_target_signals.py` (strict scalar extraction across HDF5 archives), `SNR_Calculator.py` (noise/statistics over TOF spectra), `time2mass.py` (correlating time-of-flight with mass combs), and `histograms_maps.py` (synthetic Mollweide counts maps), extending the toolkit into science-ready exploratory analysis.

