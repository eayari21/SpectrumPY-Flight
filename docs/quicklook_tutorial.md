# SpectrumPY Quicklook Tutorial

The refreshed IDEX Quicklook application is designed to take you from raw HDF5/CDF science files to annotated plots in a few clicks. This guide walks through every major surface of the UI, from opening data to exploring the in-app documentation center.

## 1. Launching the viewer

1. Activate your Python environment with the GUI dependencies (`PySide6`/`PyQt6`, `matplotlib`, `numpy`, `h5py`, `cdflib`).
2. From the repository root, run `python IDEX-quicklook.py`. The welcome dialog lets you choose an HDF5 or CDF product the moment the window opens.
3. If you call the script with `--filename` and optionally `--eventnumber`, the viewer will skip the chooser and load that event immediately.

## 2. Opening HDF5 and CDF products

* **Open…** (toolbar or `File → Open…`) provides a unified dialog for HDF5/CDF inputs. The helper automatically defaults to the repository’s `HDF5/` directory the first time you launch it.
* **Open CDF…** starts the dialog in the `CDF/` directory so you can jump straight to exported Common Data Format products without re-navigating.
* Use **Reload** after running a new decode or fit cycle; the viewer will reread the active file and refresh plots without rebuilding your channel selections.

## 3. Navigating events

* The event selector lives on the right side of the toolbar. Scroll or use the arrow keys to move through events; the status bar will update with the event number and timestamp metadata when available.
* Use the **Open Data Browser** action to launch the contextual CDF/HDF structure viewers in parallel windows when you need to inspect attributes or raw tables.

## 4. Channel toggles and overlays

* The toggle panel groups every waveform by detector. Buttons glow indigo when active, giving you immediate feedback on what is plotted.
* Enable **Overlay same time axis** to draw the selected channels against the same x-axis when their sampling rate aligns.
* Toggle **Show Ion Grid Fit**, **Show Target L Fit**, and related buttons to display calculated EMG/ion-grid fits alongside raw traces.
* `Edit Fit Parameters` opens an interactive dialog where you can tweak fit coefficients, apply overrides, and revert to the on-disk solution.

## 5. Plotting canvas

* Each channel is rendered on its own axis within a synchronized matplotlib canvas. Interact with the embedded navigation toolbar to zoom, pan, or export views.
* The **Export Plot** button exposes PNG/PDF/SVG exports. File names incorporate the event number, making it easy to archive comparisons.

## 6. Built-in help & documentation search

* Click the question-mark **Help** button or press `F1` to launch the new Documentation Center. It opens a split-pane window that lists every bundled guide (README, tutorials, technical memos) and renders their content in-line.
* Use the search bar at the top of the help window to search across all Markdown files. Results show the file name, line number, and a short snippet; selecting a match jumps the reader to that location.
* The menu bar also includes a quick search widget under `Help` so you can start typing without leaving the main window.

## 7. Keyboard shortcuts

| Shortcut | Action |
| -------- | ------ |
| `Ctrl+O` | Open any supported file |
| `Ctrl+R` | Reload the current file |
| `Ctrl+B` | Open the contextual data browser |
| `Ctrl+E` | Launch the fit parameter editor |
| `Ctrl+F1` | Focus the documentation search bar |
| `F1` | Open the documentation center |
| `Ctrl+Q` | Quit the application |

## 8. Troubleshooting tips

* If the viewer reports that `cdflib` or `h5py` is missing, install the optional dependency in your environment (`pip install cdflib h5py`).
* The UI requires a display server. On headless Linux nodes use `xvfb-run python IDEX-quicklook.py` to virtualize an X server.
* When the fit overlay diverges from expectations, use `Reset Fit Overrides` to discard cached tweaks. The status bar confirms when overrides have been cleared.

---

For a deeper dive into the science pipelines, packet decoders, and automation scripts, see the project [README](../README.md) and the documentation index surfaced by the in-app help center.
