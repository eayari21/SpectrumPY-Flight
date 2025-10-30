# SpectrumPY Packaging & Distribution Playbook

This playbook translates the bare-bones build checklist into a comprehensive,
step-by-step packaging tutorial for every supported platform. Follow it when you
need to create a new PyInstaller bundle, ship an update, or validate that the
frozen desktop experience still mirrors the developer environment.

## Table of contents

1. [Overview](#overview)
2. [Pre-flight checks](#pre-flight-checks)
3. [Shared build workflow](#shared-build-workflow)
4. [Platform-specific packaging](#platform-specific-packaging)
5. [Running and validating the bundles](#running-and-validating-the-bundles)
6. [Shipping updates](#shipping-updates)
7. [Integrating with the Spectrum Launcher](#integrating-with-the-spectrum-launcher)
8. [When PyInstaller refuses to cooperate](#when-pyinstaller-refuses-to-cooperate)
9. [Troubleshooting reference](#troubleshooting-reference)

## Overview

SpectrumPY uses a single PyInstaller specification (`packaging/idex_quicklook.spec`)
to freeze the Quicklook GUI and ship a self-contained desktop application on
macOS, Windows, and Linux. The spec explicitly pulls the in-repo documentation,
Qt assets, and Matplotlib data files so the bundled application works offline
and renders correctly on first launch.【F:packaging/idex_quicklook.spec†L1-L89】

The packaging directory supplements the spec with helper scripts (such as the
macOS DMG creator) and a focused requirements file that installs PyInstaller and
the build-time hooks.【F:packaging/README.md†L1-L48】

## Pre-flight checks

Before invoking PyInstaller, complete these sanity checks to eliminate most
runtime surprises:

1. **Clean environment.** Create a fresh Python ≥3.10 virtual environment on the
   target platform. Mixing dependencies from the development environment can
   cause frozen executables to miss compiled extensions at runtime.
2. **Install build requirements.** Upgrade `pip` and install the packaging
   requirements file. This pulls the exact PyInstaller release and helper
   libraries used by the GitHub Actions workflow.【F:packaging/README.md†L10-L32】
3. **Verify documentation assets.** Confirm that `README.md`, the `docs/`
   directory, and key imagery in `Images/` exist. The spec copies these files so
   the in-app help center has full context when offline.【F:packaging/idex_quicklook.spec†L21-L47】
4. **Check Qt bindings.** Ensure `PySide6` works by launching `python -c "from
   PySide6.QtWidgets import QApplication"`. PyInstaller embeds Qt, but build-time
   imports must succeed.
5. **Run the Quicklook GUI once.** Execute `python IDEX-quicklook.py` and open a
   sample dataset to verify local fonts, Matplotlib backends, and mission assets.
   Troubleshoot here rather than after the freeze step.

## Shared build workflow

Every platform follows the same high-level build recipe:

1. **Activate the clean environment.**
2. **Install the mission runtime** (the same dependencies developers use) and
   the packaging extras:
   ```bash
   python -m pip install --upgrade pip
   pip install -r packaging/packaging-requirements.txt
   pip install -r requirements.txt  # or the project-specific dependency list
   ```
3. **Invoke PyInstaller with the shared spec:**
   ```bash
   pyinstaller --noconfirm packaging/idex_quicklook.spec
   ```
4. **Collect the bundle** from `dist/IDEX-Quicklook` (Windows/Linux) or
   `dist/IDEX-Quicklook.app` (macOS). The folder already contains the Spectrum
   Launcher, Quicklook viewer, documentation, and runtime dependencies.
5. **Archive the output** in a platform-friendly container (ZIP, DMG, or tar.gz)
   if you plan to distribute it to teammates.

The spec dynamically includes Matplotlib, seaborn, NumPy, SciPy, and Qt hidden
imports so you rarely need to hand-edit module lists.【F:packaging/idex_quicklook.spec†L49-L86】

## Platform-specific packaging

### macOS (Intel & Apple silicon)

1. **Match architectures.** Build on the same architecture you plan to ship to
   (Intel builds on macOS 13, Apple silicon on macOS 14 when using CI runners).
2. **Run the shared build workflow.** You should end up with
   `dist/IDEX-Quicklook.app`.
3. **Codesign (optional but recommended).**
   ```bash
   codesign --deep --force --options runtime \
     --sign "Developer ID Application: Your Name" \
     dist/IDEX-Quicklook.app
   ```
4. **Notarize (optional).** Submit the signed app to Apple with
   `xcrun notarytool submit` and staple the ticket once approved.
5. **Create a DMG.** Use the helper script to generate a distributable disk
   image:
   ```bash
   ./packaging/macos/create_dmg.sh dist/IDEX-Quicklook.app \
     "IDEX-Quicklook-$(git describe --tags --always)-$(uname -m).dmg"
   ```
   The script arranges the app bundle on a canvas and preserves the signed
   binaries.【F:packaging/README.md†L20-L40】
6. **Smoke-test Gatekeeper.** Mount the DMG on a clean macOS machine, drag the
   app to `/Applications`, and launch it to confirm codesigning/notarization
   succeeded.

### Windows

1. **Run the shared build workflow** in an elevated `cmd` or PowerShell session.
   PyInstaller emits `dist/IDEX-Quicklook/IDEX-Quicklook.exe`.
2. **Bundle for distribution.** Zip the `IDEX-Quicklook` folder (which includes
   all DLLs, Qt plugins, and assets) before sharing:
   ```powershell
   Compress-Archive -Path dist/IDEX-Quicklook -DestinationPath IDEX-Quicklook.zip
   ```
3. **Optional installer.** Feed the folder into WiX, Inno Setup, or another MSI
   creator if your team prefers a guided installer. Ensure the installer places
   the files in a writable directory (e.g., `%ProgramFiles%\SpectrumPY`).
4. **First-launch verification.** From a non-development account, double-click
   `IDEX-Quicklook.exe` and ensure the Spectrum Launcher opens with mission
   imagery, file selection, and Quicklook/HDF buttons enabled after choosing an
   HDF5 or CDF file.【F:spectrum_launcher.py†L167-L247】

### Linux

1. **Run the shared build workflow.** The output lives at
   `dist/IDEX-Quicklook/IDEX-Quicklook`.
2. **Archive the bundle** for distribution:
   ```bash
   tar -C dist -czf IDEX-Quicklook-linux.tar.gz IDEX-Quicklook
   ```
3. **Library compatibility.** Verify the target systems provide a modern glibc
   (Ubuntu 20.04+ works without extra steps). If you hit loader errors, rebuild
   on the distribution you plan to support.
4. **Desktop integration.** Optionally craft a `.desktop` file that calls the
   launcher binary so users can pin it to their environment.

## Running and validating the bundles

After packaging, validate that the frozen environment mirrors the source tree:

1. **Launch through the Spectrum Launcher.** Double-click the packaged binary or
   run `./IDEX-Quicklook` (Linux), `IDEX-Quicklook.exe` (Windows), or open the
   `.app` bundle (macOS). The launcher should:
   * Display the IMAP/IDEX imagery and metadata.【F:spectrum_launcher.py†L124-L169】
   * Allow browsing for `.h5`, `.hdf5`, `.cdf`, or `.trc` files and enable the
     appropriate buttons once a file is chosen.【F:spectrum_launcher.py†L176-L247】
   * Open the HDF Plotter or IDEX Quicklook in separate windows with the selected
     dataset.【F:spectrum_launcher.py†L249-L289】
2. **Use the welcome screen shortcuts.** Test both buttons to ensure the child
   windows spawn and remain registered so closing them does not terminate the
   launcher unexpectedly.【F:spectrum_launcher.py†L249-L279】
3. **Verify documentation availability.** Open the Quicklook viewer, press `F1`,
   and confirm the documentation center lists every Markdown guide (including
   this playbook) thanks to the spec bundling the `docs/` directory.【F:packaging/idex_quicklook.spec†L21-L47】
4. **Exercise plotting workflows.** Load an event, toggle channels, and export a
   plot to make sure Matplotlib/seaborn assets and fonts load correctly.

## Shipping updates

SpectrumPY relies on Git tags and the `desktop-builds.yml` workflow to automate
bundle creation. When you publish a new release, the workflow regenerates the
platform artifacts and attaches them to the release page for download.【F:packaging/README.md†L49-L72】

To deliver an update manually or via automation:

1. **Merge the new code** into the main branch and update documentation as
   needed.
2. **Tag the release:** `git tag -a vX.Y.Z -m "SpectrumPY vX.Y.Z" && git push --tags`.
3. **Monitor the workflow run** and download the generated `.dmg`, `.zip`, and
   `.tar.gz` artifacts once they finish.
4. **Notify users** to replace their previous installation with the new bundle.
   Because configuration lives outside the distribution folder, updates are
   non-destructive.【F:packaging/README.md†L64-L79】
5. **Document release notes.** Highlight user-visible changes, migration steps,
   and any new dependencies so analysts know when to upgrade.

## Integrating with the Spectrum Launcher

Every bundle includes `spectrum_launcher.py` and the minimalist `start.py`
wrapper. Encourage users to launch via the Spectrum Launcher because it provides
contextual imagery, guided file selection, and one-click access to the HDF
Plotter or full Quicklook interface.【F:start.py†L1-L7】【F:spectrum_launcher.py†L167-L289】

* **Default entry point.** The macOS `.app` and Windows executable call the
  launcher's `main()` function on startup, mirroring the behavior of running
  `python start.py` in a source checkout.
* **File chooser defaults.** The launcher opens in the `HDF5/` directory if it
  exists, falling back to the repository root. Users can switch to CDFs or trace
  files without restarting.【F:spectrum_launcher.py†L220-L247】
* **Multiple sessions.** Each viewer window registers itself with the launcher
  so closing Quicklook or the HDF Plotter leaves the welcome screen running for
  extra files.【F:spectrum_launcher.py†L249-L279】

## When PyInstaller refuses to cooperate

Local packaging environments can fail for a host of reasons—from mismatched
Qt binaries to stale build artifacts. Before throwing out the workflow, walk
through these guardrails:

1. **Reset the build directory.** Delete the `build/` and `dist/` folders and
   rerun PyInstaller with the `--clean` flag to discard cached metadata:
   ```bash
   rm -rf build dist
   pyinstaller --noconfirm --clean packaging/idex_quicklook.spec
   ```
2. **Recreate the environment.** Use Python 3.10 and install the pinned
   packaging requirements to guarantee a known-good PyInstaller/Qt pairing:
   ```bash
   python -m pip install --upgrade pip
   pip install -r packaging/packaging-requirements.txt
   ```
   This pulls the same versions that power the automated builds.【F:packaging/packaging-requirements.txt†L1-L13】
3. **Compare against CI.** If the build still fails, trigger the
   [desktop-builds workflow](../.github/workflows/desktop-builds.yml) or grab
   the latest release artifacts. Each run emits notarization-ready `.dmg`
   images plus Windows and Linux bundles, so you can validate fixes without a
   local freeze.【F:.github/workflows/desktop-builds.yml†L1-L82】

Once you have a working artifact from CI, you can diff the PyInstaller logs or
the embedded Python runtime to isolate machine-specific quirks before retrying
locally.

## Troubleshooting reference

* **Missing Qt platform plugins.** Ensure `PySide6` was installed in the build
  environment before running PyInstaller. The spec already collects the plugins,
  but build-time discovery depends on a healthy installation.
* **Documentation center shows empty list.** Confirm the `docs/` directory was
  present when you froze the app. If you renamed or moved documentation, update
  the spec's `datas` section accordingly.【F:packaging/idex_quicklook.spec†L21-L47】
* **Launch buttons stay disabled.** The launcher only enables the HDF button for
  `.h5/.hdf5` files and enables Quicklook for `.h5`, `.hdf5`, `.cdf`, or `.trc`.
  Pick a file with the appropriate extension or adjust
  `SUPPORTED_DATA_EXTENSIONS` in `spectrum_launcher.py` if you add new formats.
* **Linux glibc errors.** Rebuild the bundle on a distribution that matches your
  target environment or use a containerized build image pinned to the correct
  glibc release.

Keeping this playbook close by should make preparing, validating, and shipping
SpectrumPY desktop bundles a predictable process.
