# IDEX L2A CDF structure vs. existing HDF5 quicklook format

This note captures the findings from scanning the new Level-2A CDF delivery and
compares the variable layout with the event-centric HDF5 files used by the
current `IDEX-quicklook` tool.  The summary is based on the ASCII metadata scan
produced by `tools/cdf_strings_report.py` (`docs/l2a_cdf_ascii_scan.md`) and the
existing HDF5 writer logic in `ImpactBook.py`.

## L2A CDF highlights

* **Primary dimensions** – The CDF advertises an `epoch` record and two time
  indices (`time_high_sample_rate_index`, `time_low_sample_rate_index`) that are
  reused across waveform variables.  The companion time axes document the
  sampling steps: the high-rate axis is tied to the TOF channels (~1/260 µs
  steps) while the low-rate axis covers Ion Grid and Target sensors (~1/4.025 µs
  steps).【F:docs/l2a_cdf_ascii_scan.md†L15-L38】
* **TOF family** – The quicklook TOF panels correspond to `tof_high`,
  `tof_mid`, and `tof_low`, with associated fit products (EMG fit parameters,
  area, χ², reduced χ², converted mass axis, and SNR).  The metadata shows those
  arrays depend on `epoch`, `mass_index`, and (for parameters) a
  `peak_fit_parameter_index`, mirroring the peak-fitting outputs already plotted
  in the HDF workflow.【F:docs/l2a_cdf_ascii_scan.md†L40-L68】
* **Target and Ion Grid fits** – For each low-rate channel the CDF includes the
  fitted waveform (`*_fit_results` spectrograms on
  `time_low_sample_rate_index`), impact charge, velocity, mass, and χ² quality
  metrics.  The context strings list the same physical descriptions used in the
  quicklook analysis (e.g., “Estimated Velocity (Target Low)” and “Charge from
  the dust impact on Ion Grid channel…”).【F:docs/l2a_cdf_ascii_scan.md†L75-L148】
* **Fit parameter dictionary** – `target_fit_parameter_index` and
  `target_fit_parameter_labels` supply the label sequence for the per-channel
  fit parameters (impact time, offset, amplitude, rise/decay constants).【F:docs/l2a_cdf_ascii_scan.md†L69-L80】
* **Derived axis labelling** – Although the ASCII scan yields raw tokens, the
  context snippets preserve unit annotations such as `km/s`, `Kg`, and `pC`
  alongside the derived quantities.  The script is therefore sufficient to
  recover the human-readable axis labels even before a full CDF reader is in
  place.【F:docs/l2a_cdf_ascii_scan.md†L85-L148】

## HDF5 quicklook expectations

* Each impact is stored in its own group (`/{Number}/…`) containing the raw
  waveforms and the aligned time vector.  The column names written by
  `ImpactBook` match the quicklook channel labels (`TOF L/M/H`, `Ion Grid`,
  `Target L/H`, etc.) and reuse the same time column for all tracks.【F:ImpactBook.py†L163-L172】
* Derived products (fit arrays, scalar diagnostics, accelerator solution, etc.)
  live in an `Analysis/` sub-group under the event.  Every derived value is
  persisted as an individual dataset whose name matches the quicklook sidebar
  labels (e.g., “Target Fit Result”, “Ion Grid Fit Time”, “Particle Velocity
  (kmps, QD Fit)”).【F:ImpactBook.py†L228-L270】

## Mapping considerations

1. **Event granularity** – The CDF is organized as per-impact records indexed
   by `epoch`.  Quicklook expects to populate an in-memory event dictionary per
   shot.  Using `cdflib` (or the ASCII scan as a temporary aid) we can iterate
   over the CDF record dimension, assembling a structure that mirrors the HDF
   event group (`waveforms`, `time axes`, `analysis metrics`).
2. **Waveform channels** – The TOF low/mid/high, Target low/high, and Ion Grid
   spectra should be read into arrays shaped `(epoch, sample_index)` with the
   corresponding `time_*_sample_rate` axis.  This maps directly onto the HDF5
   datasets created in quicklook (`/{Number}/TOF H`, `/Target H`, etc.).【F:docs/l2a_cdf_ascii_scan.md†L31-L133】【F:ImpactBook.py†L163-L172】
3. **Fit metadata** – The CDF fit products already provide both parameter
   indices and evaluation grids; they can populate the quicklook “Analysis”
   datasets without recomputation.  For example, `target_low_fit_results`
   corresponds to `Analysis/Target Fit Result`, while the scalar charge and χ²
   entries map to their similarly named datasets in the HDF writer.【F:docs/l2a_cdf_ascii_scan.md†L82-L148】【F:ImpactBook.py†L228-L270】
4. **Axis and unit labels** – The context strings expose unit hints (`pC`,
   `km/s`, `Kg`).  When wiring the quicklook GUI to CDF, these labels can seed the
   plot annotations exactly as today; once a full CDF reader is available we can
  capture the precise metadata (LABLAXIS/UNITS) rather than relying on the ASCII
   scan.

## Recommended next steps for CDF ingestion

1. **Adopt `cdflib` in the environment** – The ASCII approach is helpful for
   reconnaissance, but reading the actual arrays will require the `cdflib` (or
   NASA CDF) runtime so that the quicklook loader can request variables by name
   and respect record variance.
2. **Build a CDF adapter module** – Implement a loader that returns a
   quicklook-compatible dictionary:
   * Retrieve the shared axes (`time_high_sample_rate`,
     `time_low_sample_rate`) and index variables once per file.
   * For each `epoch` record, slice out the waveform variables and derived
     metrics, normalizing names to the existing quicklook keys.
   * Preserve the fit arrays and scalar quantities exactly as stored in the
     CDF; the HDF quicklook pipeline already expects these values and units.
3. **Extend the GUI detection logic** – Update the file-open dialog and
   quicklook controller to recognize `.cdf` inputs, calling the new adapter
   instead of the HDF reader while reusing the downstream plotting code.
4. **Validate with side-by-side plots** – Once the adapter is in place, compare
   the CDF-driven quicklook output with an HDF event to ensure waveform shapes,
   axes, and derived metrics agree within numerical tolerance.  The ASCII scan
   (`docs/l2a_cdf_ascii_scan.md`) can be used as a checklist to confirm that all
   expected variables were mapped.
