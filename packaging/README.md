# Desktop packaging guide

This folder bundles the scripts and documentation required to publish the
`IDEX-quicklook.py` GUI as a standalone desktop application for macOS (Intel and
Apple silicon), Windows, and Linux. The build process relies on
[PyInstaller](https://pyinstaller.org/) to freeze the Python environment and
captures the Markdown documentation so the in-app help browser works offline.

If you prefer a narrative, end-to-end walkthrough, read the companion
[**SpectrumPY Packaging & Distribution Playbook**](../docs/packaging_tutorial.md)
for annotated checklists, validation steps, and troubleshooting recipes.

## Build overview

1. Create a clean Python 3.10+ environment on the target platform.
2. Install the mission dependencies alongside the packaging helpers:
   ```bash
   python -m pip install --upgrade pip
   pip install -r packaging/packaging-requirements.txt
   ```
3. Invoke PyInstaller with the shared spec:
   ```bash
   pyinstaller --noconfirm packaging/idex_quicklook.spec
   ```
4. Retrieve the bundle from the `dist/` directory and wrap it in a platform
   specific container (DMG, installer, or tarball).

The `packaging/idex_quicklook.spec` definition collects the documentation,
Qt/Matplotlib assets, and scientific dependencies that are dynamically imported
at runtime. The same spec works across macOS, Windows, and Linux, so you only
need to maintain one configuration.【F:packaging/idex_quicklook.spec†L1-L89】

## Platform notes

### macOS (Intel & Apple silicon)

* Run the build on the matching architecture (`macos-13` for Intel, `macos-14`
  for Apple silicon when using GitHub Actions). PyInstaller must run on the same
  architecture that will execute the app.
* After the build completes you will have `dist/IDEX-Quicklook.app`. Create a
  signed disk image with:
  ```bash
  ./packaging/macos/create_dmg.sh dist/IDEX-Quicklook.app "IDEX-Quicklook-$(git describe --tags --always)-$(uname -m).dmg"
  ```
* Codesign and notarize the `.app`/`.dmg` before distributing it outside your
  organization. The script leaves hooks for the usual `codesign` + `xcrun
  notarytool` commands.

### Windows

* PyInstaller produces `dist/IDEX-Quicklook/IDEX-Quicklook.exe`. Compress the
  folder into a ZIP file for side-loading or feed it into an MSI creator such as
  WiX if you need an installer.
* Executables inherit the embedded Python runtime, so users only need to unzip
  the folder and run `IDEX-Quicklook.exe`.

### Linux

* PyInstaller emits `dist/IDEX-Quicklook/IDEX-Quicklook`. Package it as a tar.gz
  for distribution:
  ```bash
  tar -C dist -czf IDEX-Quicklook-linux.tar.gz IDEX-Quicklook
  ```
* On first launch the executable will create a temporary folder for extracted
  shared libraries; no root privileges are required. Ensure the target systems
  provide a compatible glibc (Ubuntu 20.04+ works out of the box).

## Automated builds with GitHub Actions

The repository includes `.github/workflows/desktop-builds.yml`, which builds the
bundles for macOS (Intel/Apple silicon), Windows, and Linux whenever a release
is created or the workflow is run manually. Each job publishes zipped artifacts
(`.dmg`, `.zip`, `.tar.gz`) that you can attach to a GitHub release. Triggering a
new build is as easy as pushing a tag or using the **Run workflow** button.

## Shipping updates to users

When you push changes to the repository and tag a release, the GitHub Actions
workflow rebuilds the executables automatically. Users can update by downloading
and installing the latest artifact from the release page:

1. Visit the project’s **Releases** page and grab the newest platform package.
2. Replace the previous installation (overwrite the `.app`, unzip the fresh
   Windows folder, or untar the Linux bundle).
3. Launch the application – configuration and generated data live outside the
   bundle, so updates are non-destructive.

This approach keeps distribution simple: you publish updated archives, and users
opt-in to the latest tools whenever new science capabilities land in the repo.
