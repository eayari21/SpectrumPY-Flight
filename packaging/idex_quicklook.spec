# -*- mode: python ; coding: utf-8 -*-

"""PyInstaller specification for the SpectrumPY IDEX Quicklook desktop bundles.

This spec is intentionally shared across macOS, Windows, and Linux builds. It
collects project documentation, Qt assets, and Matplotlib data so the GUI works
out of the box in frozen bundles.
"""

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules


REPO_ROOT = Path(__file__).resolve().parent.parent
ENTRYPOINT = REPO_ROOT / "IDEX-quicklook.py"

if not ENTRYPOINT.exists():
    raise SystemExit(f"Unable to locate GUI entry point: {ENTRYPOINT}")

# ---------------------------------------------------------------------------
# Data files bundled with the executable.
# ---------------------------------------------------------------------------
datas = []

# Ship the top-level README and Markdown docs so the in-app documentation
# browser can render content without requiring network access.
readme_path = REPO_ROOT / "README.md"
if readme_path.exists():
    datas.append((str(readme_path), "."))

docs_dir = REPO_ROOT / "docs"
if docs_dir.exists():
    datas.append((str(docs_dir), "docs"))

# Matplotlib, seaborn, and Qt all load additional data files at runtime.
# `collect_data_files` makes sure the frozen bundle contains fonts, stylesheets,
# and plugin libraries that would otherwise be missing.
datas.extend(collect_data_files("matplotlib", include_py_files=False))
datas.extend(collect_data_files("seaborn", include_py_files=False))
datas.extend(collect_data_files("PySide6"))

# ---------------------------------------------------------------------------
# Hidden imports required by scientific packages that rely on lazy loading.
# ---------------------------------------------------------------------------
hiddenimports = []
hiddenimports.extend(collect_submodules("matplotlib"))
hiddenimports.extend(collect_submodules("mpl_toolkits"))
hiddenimports.extend(collect_submodules("numpy"))
hiddenimports.extend(collect_submodules("scipy"))
hiddenimports.extend(collect_submodules("PySide6"))

# ---------------------------------------------------------------------------
# PyInstaller build graph
# ---------------------------------------------------------------------------
block_cipher = None


a = Analysis(
    [str(ENTRYPOINT)],
    pathex=[str(REPO_ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="IDEX-Quicklook",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

if sys.platform == "darwin":
    app = BUNDLE(
        exe,
        name="IDEX-Quicklook.app",
        icon=None,
        bundle_identifier="edu.colorado.lasp.idex.quicklook",
        info_plist={"NSHighResolutionCapable": True},
    )
elif sys.platform == "win32":
    coll = COLLECT(
        exe,
        a.binaries,
        a.zipfiles,
        a.datas,
        strip=False,
        upx=True,
        upx_exclude=[],
        name="IDEX-Quicklook",
    )
else:
    coll = COLLECT(
        exe,
        a.binaries,
        a.zipfiles,
        a.datas,
        strip=False,
        upx=True,
        upx_exclude=[],
        name="IDEX-Quicklook",
    )
