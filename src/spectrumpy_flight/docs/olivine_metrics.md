# Olivine Test Metrics

The olivine regression test now produces a complete calibration package each time
it is executed. The workflow was designed to satisfy the SNR, timing, and fit
quality requirements highlighted in the instrument signal handbook (see the
specification snapshot provided with the task request).

## Running the pytest harness

1. Place the calibrated olivine binary frames under `Data/` using the same
   naming convention as the existing fixtures (`MMDDYYYY` date tokens are used
   to detect the observation window).
2. Launch the regression from the repository root:

   ```bash
   pytest tests/test_olivine_analysis.py
   ```

   Pytest reports the temporary directory that contains the generated HDF5
   products and the metrics bundle. To persist the artefacts, supply a base
   directory, for example `pytest --basetemp=/tmp/olivine-run`.

3. Inspect the outputs located under `<tmp>/metrics/`:

   - `olivine_metrics_report.pdf` combines 21 histograms (SNR, rise/decay time,
     trigger deltas, χ², and reduced χ²) with a summary page that also includes
     saturation counts.
   - `olivine_metrics_summary.json` captures the numerical aggregates used to
     build the plots so the trends can be post-processed without re-running the
     test.

The regression continues to assert that the analysis HDF5 tree contains the
expected fit products while simultaneously generating the documentation bundle.

## Standalone CLI usage

The reporting machinery is accessible outside of pytest for ad-hoc analyses or
backfilling historical datasets.

```bash
python -m olivine_metrics Data/decoded_runs --output-dir reports/olivine
```

Any mix of HDF5 files or directories can be supplied; the script crawls each
directory for `*.h5` files, aggregates the metrics, and writes the PDF plus JSON
summary to the requested output directory.

## Contents of the report

Each histogram page includes axis labels, a descriptive caption, and the
relevant requirement threshold when available:

- **SNR** — One histogram per channel (TOF L/M/H, Ion Grid, Target H/L) with the
  minimum acceptable SNR annotated.
- **Rise/decay constants** — Rise (`τ_rise`) and decay (`τ_decay`) parameters for
  Ion Grid, Target H, and Target L fits with the permitted ranges shaded as
  reference lines.
- **Trigger offsets** — Absolute timing differences for the Ion Grid↔TOF H,
  TOF H↔Target, and Target↔Ion Grid pairs (20 µs upper bound).
- **Fit quality** — χ² and reduced χ² distributions for each low-rate channel.
- **Saturation statistics** — Reported on the summary page together with the
  number of processed events.

The JSON file stores counts, extrema, and central-tendency statistics for every
metric so mission analysts can rapidly query compliance across campaigns.
