#!/usr/bin/env python3
"""Rebuild oxDNA/oxpy solely from the flat sources directory shipped with the app.

Requires Python 3.11, CMake and a C++ compiler plus the packages recorded in
BUILDING.txt. Does not import Apta2Switch or download scientific source files.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import runpy
import shutil
import subprocess
import sys
import tarfile
import tempfile

REVISION = "75ae96f063954b2f29524bb897c274aee7ba4466"
ARCHIVE_SHA256 = "29e7780829ee967046168e55851a1e51017c6f37daeed3b8550b8d855704c878"
WORKERS = ("oxdna_ffs_worker", "oxdna_kinetics_worker", "oxdna_native_ffs", "oxdna_kinetics",
           "oxdna_initialization", "oxdna_guided_worker", "oxdna_geometry", "oxdna_loop_geometry", "oxdna_box", "oxdna_free_md_worker", "oxdna_free_md_metrics")
PARAMETERS = dict(zip(("DNA", "RNA", "NA"), ("oxDNA2_sequence_dependent_parameters.txt",
                 "rna_sequence_dependent_parameters.txt", "NA_sequence_dependent_parameters.txt")))


def verify_sources(directory: Path) -> None:
    required = ("oxDNA-source.tar.gz", "prepare_oxdna_ffs_runtime.py", "prepare_bundled_science.py",
                "run_oxdna_worker.py", *(name + ".py" for name in WORKERS))
    for name in required:
        if not (directory / name).is_file():
            raise RuntimeError(f"Missing distributed source: {name}")
    with (directory / "oxDNA-source.tar.gz").open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != ARCHIVE_SHA256:
            raise RuntimeError("The original oxDNA source archive checksum differs")


def rebuild(directory: Path, output: Path, jobs: int) -> None:
    verify_sources(directory)
    if output.exists():
        raise ValueError("Choose a new output directory; an existing runtime is never overwritten")
    # Only standard-library imports execute when loading this shipped helper.
    # These two helpers do not depend on an application checkout.
    helpers = runpy.run_path(str(directory / "prepare_bundled_science.py"))
    with tempfile.TemporaryDirectory(prefix="oxdna-source-rebuild-") as temporary:
        stage = Path(temporary)
        with tarfile.open(directory / "oxDNA-source.tar.gz") as archive:
            archive.extractall(stage, filter="data")
        source = stage / ("oxDNA-" + REVISION)
        runpy.run_path(str(directory / "prepare_oxdna_ffs_runtime.py"))["repair_source"](source)
        subprocess.run(["cmake", "-S", str(source), "-B", str(source / "build"), "-DPython=ON",
                        "-DOxpySystemInstall=ON", "-DCUDA=OFF", "-DNATIVE_COMPILATION=OFF",
                        "-DCMAKE_BUILD_TYPE=Release", f"-DDEFAULT_PYTHON={sys.executable}"], check=True)
        subprocess.run(["cmake", "--build", str(source / "build"), "--target", "oxDNA", "core",
                        "--parallel", str(jobs)], check=True)
        subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-build-isolation",
                        "--force-reinstall", str(source / "build/python")], check=True)
        command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir",
                   "--name", "aptaswitch-oxdna-worker", "--distpath", str(stage / "dist"),
                   "--workpath", str(stage / "work"), "--specpath", str(stage),
                   "--paths", str(directory), "--hidden-import", "oxpy",
                   "--collect-all", "oxpy", "--collect-all", "numpy",
                   "--add-binary", str(source / "build/bin/oxDNA") + ":bin"]
        for name in WORKERS:
            command += ["--hidden-import", name]
        subprocess.run([*command, str(directory / "run_oxdna_worker.py")], check=True,
                       env=dict(os.environ, PYINSTALLER_CONFIG_DIR=str(stage / "cache")))
        output.mkdir(parents=True)
        shutil.copytree(stage / "dist/aptaswitch-oxdna-worker", output / "worker", symlinks=True)
        (output / "parameters").mkdir()
        for name in PARAMETERS.values():
            shutil.copyfile(source / name, output / "parameters" / name)
        shutil.copyfile(source / "LICENSE", output / "LICENSE.txt")
        (output / "sources").mkdir()
        for path in directory.iterdir():
            if path.is_file():
                shutil.copyfile(path, output / "sources" / path.name)
        helpers["dependency_notices"](output)
    metadata = helpers["manifest"](output, "oxdna", worker="worker/aptaswitch-oxdna-worker",
        binary="worker/_internal/bin/oxDNA", parameters={key: "parameters/" + name for key, name in PARAMETERS.items()},
        source_revision=REVISION, source_archive_sha256=ARCHIVE_SHA256, license="GPL-3.0",
        system=platform.system(), machine=platform.machine())
    runtime = {"python": str(output / metadata["worker"]), "binary": str(output / metadata["binary"]),
               **{"seq_dep_file_" + key: str(output / value) for key, value in metadata["parameters"].items()}}
    subprocess.run([runtime["python"], str(output / "sources/oxdna_native_ffs.py"), "--probe", json.dumps(runtime)],
                   check=True, cwd=output, timeout=60)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--verify-sources", action="store_true")
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    if args.verify_sources:
        verify_sources(directory)
        print("Complete local source set and original archive verified")
        return
    if sys.version_info[:2] != (3, 11) or sys.platform not in ("darwin", "linux"):
        parser.error("Use Python 3.11 on macOS or Linux")
    if args.output is None or args.jobs < 1:
        parser.error("--output and a positive --jobs are required")
    rebuild(directory, args.output.expanduser().resolve(), args.jobs)


if __name__ == "__main__":
    main()
