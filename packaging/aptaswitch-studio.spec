# -*- mode: python ; coding: utf-8 -*-

import os
import sys

from PyInstaller.utils.hooks import collect_all, collect_submodules


# SPECPATH is injected by PyInstaller and points at this file's directory.
ROOT = os.path.abspath(os.path.join(SPECPATH, os.pardir))
STREAMLIT_DATAS, STREAMLIT_BINARIES, STREAMLIT_HIDDEN = collect_all("streamlit")
RESVG_DATAS, RESVG_BINARIES, RESVG_HIDDEN = collect_all("resvg_py")

if sys.platform == "darwin":
    ICON = os.path.join(ROOT, "packaging", "icons", "AppIcon.icns")
elif sys.platform == "win32":
    ICON = os.path.join(ROOT, "packaging", "icons", "AppIcon.ico")
else:
    ICON = None

a = Analysis(
    [os.path.join(ROOT, "scripts", "run_app.py")],
    pathex=[ROOT, os.path.join(ROOT, "src")],
    binaries=STREAMLIT_BINARIES + RESVG_BINARIES,
    datas=STREAMLIT_DATAS
    + RESVG_DATAS
    + [
        (
            os.path.join(ROOT, "src", "aptaswitch_studio", "web_app.py"),
            "aptaswitch_studio",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_studio", "assets"),
            "aptaswitch_studio/assets",
        ),
    ],
    hiddenimports=STREAMLIT_HIDDEN
    + RESVG_HIDDEN
    + collect_submodules("aptaswitch_core")
    + collect_submodules("aptaswitch_studio")
    + [
        "openpyxl",
        "numpy",
        "pandas",
        "scipy.sparse",
        "scipy.interpolate",
        "scipy.optimize",
        "yaml",
        "jinja2",
        "packaging.tags",
        "packaging.utils",
        "streamlit.testing.v1",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["nupack", "multistrand", "PySide6"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Apta2Switch-Studio",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    icon=ICON,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="Apta2Switch-Studio",
)

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="Apta2Switch-Studio.app",
        icon=ICON,
        bundle_identifier="com.apta2switch-studio.app",
        version="0.2.1",
        info_plist={
            "CFBundleName": "Δpta2Switch-Studio",
            "CFBundleDisplayName": "Δpta2Switch-Studio",
            "CFBundleShortVersionString": "0.2.1",
            "CFBundleVersion": "0.2.1b1",
            "NSHighResolutionCapable": True,
            "LSMultipleInstancesProhibited": True,
            "NSHumanReadableCopyright": "Apta2Switch-Studio contributors",
        },
    )
