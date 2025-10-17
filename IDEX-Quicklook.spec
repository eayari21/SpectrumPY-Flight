# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.utils.hooks import collect_dynamic_libs
from PyInstaller.utils.hooks import collect_submodules
from PyInstaller.utils.hooks import collect_all

datas = []
binaries = []
hiddenimports = ['importlib_resources.trees', 'matplotlib.backends.backend_qtagg', 'PyQt6.sip']
datas += collect_data_files('PyQt6')
binaries += collect_dynamic_libs('PyQt6')
hiddenimports += collect_submodules('matplotlib.backends')
hiddenimports += collect_submodules('matplotlib.backends.backend_qtagg')
tmp_ret = collect_all('PyQt6')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]
tmp_ret = collect_all('matplotlib')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]


a = Analysis(
    ['IDEX-quicklook.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['numpy.tests', 'numpy.f2py.tests', 'scipy.tests', 'matplotlib.tests', 'mpl_toolkits.tests'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='IDEX-Quicklook',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch='arm64',
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='IDEX-Quicklook',
)
app = BUNDLE(
    coll,
    name='IDEX-Quicklook.app',
    icon=None,
    bundle_identifier='com.lasp.idex.quicklook',
)
