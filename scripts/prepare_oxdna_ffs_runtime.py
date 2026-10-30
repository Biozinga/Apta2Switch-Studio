#!/usr/bin/env python3
"""Repair and rebuild the supplied oxDNA CPU FFS backend when needed.

The upstream MD force-reset API changed without the same initialization being
carried into some FFS_MD builds. Apply the *same* reset as MD_CPUBackend before
the existing native FFS force calculation; no force field, order parameter,
or stopping rule is modified. A source backup and a patch manifest are kept.
The bundled runtime builder performs the behavioral force/stop check afterwards.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def repair_source(source: Path) -> bool:
    backend = source / "src/Backends/FFS_MD_CPUBackend.cpp"
    baseline = source / "src/Backends/MD_CPUBackend.cpp"
    old = "void FFS_MD_CPUBackend::_ffs_compute_forces(void) {\n\t_interaction->begin_energy_computation();"
    new = ("void FFS_MD_CPUBackend::_ffs_compute_forces(void) {\n"
           "\t_interaction->begin_energy_and_force_computation();\n"
           "\t// Keep force initialization identical to MD_CPUBackend.\n"
           "\tfor(auto p : _particles) {\n"
           "\t\tp->set_initial_forces(current_step(), _box.get());\n"
           "\t}")
    text = backend.read_text()
    if new in text:
        return False
    if old not in text:
        if ("begin_energy_and_force_computation()" in text
                and "set_initial_forces(current_step(), _box.get())" in text):
            return False  # An upstream version already contains the fix.
        raise ValueError("Unrecognized FFS_MD backend; refusing an ambiguous source patch.")
    md = baseline.read_text()
    if not ("begin_energy_and_force_computation()" in md
            and "set_initial_forces(current_step(), _box.get())" in md):
        raise ValueError("The MD backend does not use the expected force-reset API.")
    backup = backend.with_suffix(".cpp.before_force_reset_fix")
    if backup.exists() and backup.read_text() != text:
        raise ValueError("An unrelated previous source backup exists; refusing to overwrite it.")
    backup.write_text(text)
    patched = text.replace(old, new, 1)
    backend.write_text(patched)
    (source / "aptaswitch_ffs_force_reset.json").write_text(json.dumps({
        "patch": "cpu_ffs_force_reset_v1",
        "file": "src/Backends/FFS_MD_CPUBackend.cpp",
        "original_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
        "reason": "Initialize per-step forces and torques exactly as MD_CPUBackend does.",
    }, indent=2) + "\n")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="The oxDNA source directory containing src/ and build/.")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--python", type=Path, help="External oxDNA virtual environment interpreter (required for build).")
    parser.add_argument("--patch-only", action="store_true")
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if not args.patch_only and (args.python is None or not args.python.expanduser().is_file()):
        parser.error("--python must name the external oxDNA Python interpreter")
    changed = repair_source(source)
    print("FFS CPU force initialization repaired." if changed else "FFS CPU force initialization already matches MD.")
    if not args.patch_only:
        python = str(args.python.expanduser().absolute())
        # The user's installation may have been moved; regenerate rather than
        # trusting absolute paths retained in an old CMake cache.
        subprocess.run(["cmake", "--fresh", "-S", str(source), "-B", str(source / "build"),
                        "-DPython=ON", "-DOxpySystemInstall=ON", "-DCMAKE_BUILD_TYPE=Release",
                        f"-DDEFAULT_PYTHON={python}"], check=True)
        subprocess.run(["cmake", "--build", str(source / "build"), "--target", "oxDNA", "core",
                        "--parallel", str(args.jobs)], check=True)
        subprocess.run([python, "-m", "pip", "install", "--no-deps", "--no-build-isolation",
                        "--force-reinstall", str(source / "build/python")], check=True)
        print("Build complete. Run the bundled runtime force/stop check before packaging.")


if __name__ == "__main__":
    main()
