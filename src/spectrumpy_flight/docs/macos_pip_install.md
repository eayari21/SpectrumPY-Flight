# macOS installation guide for the pip distribution

This guide walks through preparing a macOS host to install and run the
`Spectrumpy-Flight` pip package. It assumes macOS 13 (Ventura) or later with an
administrator account.

## 1. Install command-line prerequisites

1. Install the Xcode command-line tools (provides compilers and SDK headers):
   ```bash
   xcode-select --install
   ```
   Accept the prompt and wait for the installation to finish.
2. Install [Homebrew](https://brew.sh/) if it is not already on the system:
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
   Follow the instructions that Homebrew prints at the end to add it to your
   `PATH` (usually by updating `~/.zprofile`).
3. Use Homebrew to install an up-to-date Python interpreter and helpful build
   dependencies:
   ```bash
   brew install python@3.11 cmake pkg-config libffi
   ```
   The formula places `python3` in `/opt/homebrew/bin` on Apple Silicon or
   `/usr/local/bin` on Intel machines.

## 2. Create an isolated Python environment

Working inside a virtual environment keeps the IDEX tooling and its
dependencies isolated from the system interpreter.

1. Create a project folder for the virtual environment:
   ```bash
   mkdir -p ~/Envs
   cd ~/Envs
   python3 -m venv spectrumpy
   ```
2. Activate the environment whenever you plan to install or run the package:
   ```bash
   source ~/Envs/spectrumpy/bin/activate
   ```
3. Upgrade `pip`, `setuptools`, and `wheel` before installing any packages:
   ```bash
   python -m pip install --upgrade pip setuptools wheel
   ```

## 3. Install Spectrumpy-Flight from PyPI

With the environment activated, install the published package. The default
installation pulls in PySide6 for the Qt GUI bindings. If you prefer the PyQt6
bindings, use the optional extra shown below.

```bash
python -m pip install spectrumpy
# or, to install with the PyQt6 bindings instead of PySide6
python -m pip install "spectrumpy[pyqt]"
```

### Optional dependencies

* GPU acceleration for packet processing:
  ```bash
  python -m pip install "spectrumpy[gpu]"
  ```
  Choose the CuPy build that matches your CUDA toolkit if you install it
  separately.
* MySQL-compatible database connectors for legacy telemetry archives:
  ```bash
  python -m pip install "spectrumpy[database]"
  ```

## 4. Verify Qt integration

PySide6 bundles its own Qt frameworks. macOS occasionally blocks unsigned
frameworks the first time they are loaded. If the GUI refuses to start, open the
System Settings → *Privacy & Security* panel and allow the blocked Qt
components. Relaunch the app after granting access.

To confirm that the Qt platform plugins are discoverable, run the launcher from
Terminal while the virtual environment is active:

```bash
spectrumpy-gui launcher
```

You should see the SpectrumPY welcome window. To skip directly to the Quicklook
viewer, supply the `quicklook` subcommand:

```bash
spectrumpy-gui quicklook --filename /path/to/product.h5
```

## 5. Keep the installation up to date

Periodically update the package and its dependencies to pick up bug fixes and
feature improvements:

```bash
python -m pip install --upgrade spectrumpy
```

When you are done working, deactivate the environment:

```bash
deactivate
```

Reactivate it later with the `source ~/Envs/spectrumpy/bin/activate` command
before launching the tools again.
