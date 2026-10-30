# -*- mode: python ; coding: utf-8 -*-

import os
import sys
import runpy
import shutil
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_submodules, copy_metadata
from PyInstaller.building.datastruct import Tree


# SPECPATH is injected by PyInstaller and points at this file's directory.
ROOT = os.path.abspath(os.path.join(SPECPATH, os.pardir))
THERMODYNAMICS_RUNTIME = os.path.join(ROOT, ".local_deps", "thermodynamics")
THERMODYNAMICS_LIBRARY = "libaptaswitch_hybrid.dylib" if sys.platform == "darwin" else "libaptaswitch_hybrid.so"
for filename in (THERMODYNAMICS_LIBRARY, THERMODYNAMICS_LIBRARY + ".json",
                 "ViennaRNA-COPYING.txt", "THIRD-PARTY-NOTICES.txt"):
    if not os.path.isfile(os.path.join(THERMODYNAMICS_RUNTIME, filename)):
        raise SystemExit("Prepare the complete hybrid runtime with scripts/prepare_thermodynamics_runtime.py first.")
OXDNA_RUNTIME = os.path.abspath(os.environ.get(
    "APTASWITCH_OXDNA_RUNTIME", os.path.join(ROOT, ".local_deps", "oxdna-runtime")))
validate_science = runpy.run_path(os.path.join(ROOT, "scripts", "prepare_bundled_science.py"))["validate_payload"]
for directory, kind in ((OXDNA_RUNTIME, "oxdna"),):
    try:
        validate_science(Path(directory), kind)
    except (ValueError, OSError) as exc:
        raise SystemExit(f"Prepare the bundled scientific runtimes with scripts/prepare_bundled_science.py: {exc}")
STREAMLIT_DATAS, STREAMLIT_BINARIES, STREAMLIT_HIDDEN = collect_all("streamlit")
RESVG_DATAS, RESVG_BINARIES, RESVG_HIDDEN = collect_all("resvg_py")
OSTIR_DATAS, OSTIR_BINARIES, OSTIR_HIDDEN = collect_all("ostir")
RNA_DATAS, RNA_BINARIES, RNA_HIDDEN = collect_all("RNA")
# OSTIR reads .par files next to its Python package; RNA loads its native module.
# Explicit distribution metadata also preserves upstream copyright and licences.
OSTIR_DATAS += copy_metadata("OSTIR") + copy_metadata("ViennaRNA")
OSTIR_SOURCES = runpy.run_path(os.path.join(ROOT, "scripts", "prepare_ostir_sources.py"))["prepare_sources"](Path(ROOT))

if sys.platform == "darwin":
    ICON = os.path.join(ROOT, "packaging", "icons", "AppIcon.icns")
else:
    ICON = None

a = Analysis(
    [os.path.join(ROOT, "scripts", "run_app.py")],
    pathex=[ROOT, os.path.join(ROOT, "src")],
    binaries=STREAMLIT_BINARIES + RESVG_BINARIES + OSTIR_BINARIES + RNA_BINARIES,
    datas=STREAMLIT_DATAS
    + RESVG_DATAS
    + OSTIR_DATAS
    + RNA_DATAS
    + [
        (
            os.path.join(ROOT, "src", "aptaswitch_studio", "web_app.py"),
            "aptaswitch_studio",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_ffs_worker.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_kinetics_worker.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_native_ffs.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_kinetics.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_initialization.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_guided_worker.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_geometry.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_loop_geometry.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_box.py"),
            "aptaswitch_core",
        ),
        (
            os.path.join(ROOT, "src", "aptaswitch_studio", "assets"),
            "aptaswitch_studio/assets",
        ),
        (os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_free_md_worker.py"), "aptaswitch_core"),
        (os.path.join(ROOT, "src", "aptaswitch_core", "oxdna_free_md_metrics.py"), "aptaswitch_core"),
    ],
    hiddenimports=STREAMLIT_HIDDEN
    + RESVG_HIDDEN
    + OSTIR_HIDDEN
    + RNA_HIDDEN
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
    # oxpy runs exclusively in its authenticated standalone worker closure.
    excludes=["multistrand", "PySide6", "oxpy"],
    noarchive=False,
)

# Keep the isolated oxpy worker's relative library layout intact. macOS leaf
# signing refreshes its complete manifest before the enclosing app is sealed.
if sys.platform != "darwin":
    a.datas += Tree(OXDNA_RUNTIME, prefix="oxdna_runtime", typecode="DATA")
# The native loader verifies exact bytes against the manifest. BINARY entries
# can be rewritten by PyInstaller's dependency scanner; keep this self-contained
# library intact, then reseal its hash after macOS leaf signing, before app seal.
a.datas += Tree(THERMODYNAMICS_RUNTIME, prefix="thermodynamics", typecode="DATA")
a.datas += Tree(str(OSTIR_SOURCES), prefix="ostir_sources", typecode="DATA")

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
            "CFBundleVersion": "0.2.1b2",
            "NSHighResolutionCapable": True,
            "LSMultipleInstancesProhibited": True,
            "NSHumanReadableCopyright": "Apta2Switch-Studio contributors",
        },
    )
    # BUNDLE's DATA rewriting dereferences the worker's Python.framework
    # symlinks. Copy the already frozen closure intact after bundle assembly;
    # scripts/sign_macos_app.py signs it and reseals its manifest next.
    worker_destination = Path(app.name) / "Contents" / "Resources" / "oxdna_runtime"
    if worker_destination.exists():
        shutil.rmtree(worker_destination)
    shutil.copytree(OXDNA_RUNTIME, worker_destination, symlinks=True)
