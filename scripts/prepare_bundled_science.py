#!/usr/bin/env python3
"""Maintainer build of offline oxDNA/oxpy payloads.

Downloads happen only here, never in the application. The GPL oxDNA source,
our exact repair/build/worker sources, and upstream notices accompany its binary.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OXDNA_REVISION = "75ae96f063954b2f29524bb897c274aee7ba4466"
OXDNA_ARCHIVE_SHA256 = "29e7780829ee967046168e55851a1e51017c6f37daeed3b8550b8d855704c878"
PARAMETERS = ("oxDNA2_sequence_dependent_parameters.txt", "rna_sequence_dependent_parameters.txt", "NA_sequence_dependent_parameters.txt")
WORKERS = ("oxdna_ffs_worker", "oxdna_kinetics_worker", "oxdna_native_ffs", "oxdna_kinetics",
           "oxdna_initialization", "oxdna_guided_worker", "oxdna_geometry", "oxdna_loop_geometry", "oxdna_box", "oxdna_free_md_worker", "oxdna_free_md_metrics")


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def download(url, path, expected):
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        pending = path.with_suffix(path.suffix + ".download")
        with urllib.request.urlopen(url, timeout=60) as src, pending.open("wb") as out:
            shutil.copyfileobj(src, out)
        if digest(pending) != expected:
            pending.unlink()
            raise RuntimeError(f"Upstream checksum mismatch: {path.name}")
        pending.replace(path)
    if digest(path) != expected:
        raise RuntimeError(f"Cached checksum mismatch: {path.name}")
    return path


def manifest(root, kind, **metadata):
    files = {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob("*"))
             if p.is_file() and p.name != "runtime.json"}
    value = dict(schema_version=1, kind=kind, files=files, **metadata)
    (root / "runtime.json").write_text(json.dumps(value, indent=2) + "\n")
    return value


def validate_payload(root, kind=None):
    sys.path.insert(0, str(ROOT / "src"))
    from aptaswitch_core.bundled_science import payload_manifest
    kind = kind or "oxdna"
    if kind != "oxdna":
        raise ValueError("Only the oxDNA payload is supported")
    value = payload_manifest(Path(root), kind)
    required = {"worker", "binary", "parameters", "source_revision"}
    if not required.issubset(value):
        raise ValueError("Bundled scientific runtime metadata is incomplete")
    return value


def python_license_file():
    """Locate Python's full notice in framework/source and Debian installs."""
    version = f"python{sys.version_info.major}.{sys.version_info.minor}"
    candidates = (
        Path(sys.base_prefix) / "lib" / version / "LICENSE.txt",
        Path(sysconfig.get_path("stdlib")) / "LICENSE.txt",
        Path("/usr/share/doc") / version / "copyright",
        Path("/usr/share/doc") / ("lib" + version + "-stdlib") / "copyright",
    )
    for path in candidates:
        if path.is_file():
            return path
    raise RuntimeError("The build Python must provide its complete license; checked: "
                       + ", ".join(str(path) for path in candidates))


def dependency_notices(output):
    """Keep the actual interpreter/packager and native dependency notices."""
    notices = output / "licenses"
    notices.mkdir(exist_ok=True)
    shutil.copyfile(python_license_file(), notices / "Python-LICENSE.txt")
    for name in ("pyinstaller", "numpy", "setuptools"):
        distribution = importlib.metadata.distribution(name)
        found = False
        for index, entry in enumerate(sorted(distribution.files or ())):
            if ".dist-info/" in entry.as_posix() and any(word in entry.name.lower() for word in ("license", "copying", "notice")):
                source = Path(distribution.locate_file(entry))
                if source.is_file():
                    shutil.copyfile(source, notices / f"{name}-{index}-{entry.name}")
                    found = True
        if not found:
            raise RuntimeError(f"Missing installed license for bundled {name}")
    libraries = {path.name for path in (output / "worker").rglob("*") if path.is_file()}
    dependencies = (
        ("libssl", "openssl@3", "openssl"), ("libcrypto", "openssl@3", "openssl"),
        ("liblzma", "xz", "liblzma5"), ("libmpdec", "mpdecimal", "libmpdec3"),
        ("libffi", "libffi", "libffi8"), ("libexpat", "expat", "libexpat1"),
        ("libbz2", "bzip2", "libbz2-1.0"), ("libsqlite3", "sqlite", "libsqlite3-0"),
        ("libz.", "zlib", "zlib1g"),
    )
    for prefix, formula, debian in dependencies:
        if not any(name.startswith(prefix) for name in libraries):
            continue
        if sys.platform == "darwin":
            base = Path(subprocess.check_output(["brew", "--prefix", formula], text=True).strip())
            choices = sorted(path for path in base.rglob("*") if path.is_file()
                       and any(word in path.name.upper() for word in ("LICENSE", "COPYING", "COPYRIGHT", "NOTICE"))
                       and "man" not in path.parts)
        else:
            choices = [Path("/usr/share/doc") / debian / "copyright"]
        choices = [path for path in choices if path.is_file()]
        if not choices:
            raise RuntimeError(f"Missing license notice for bundled {prefix}")
        for index, source in enumerate(choices):
            shutil.copyfile(source, notices / f"{formula}-{index}-{source.name}")




def oxdna_source_archive(cache):
    """Use the authenticated source shipped in the wiki archive when present."""
    vendored = ROOT / "vendor/sources" / f"oxDNA-{OXDNA_REVISION}.tar.gz"
    if vendored.is_file():
        if digest(vendored) != OXDNA_ARCHIVE_SHA256:
            raise RuntimeError("Vendored oxDNA source archive checksum mismatch")
        return vendored
    return download(f"https://codeload.github.com/lorenzo-rovigatti/oxDNA/tar.gz/{OXDNA_REVISION}",
                    cache / f"oxdna-{OXDNA_REVISION}.tar.gz", OXDNA_ARCHIVE_SHA256)


def prepare_oxdna(cache, output, jobs):
    archive = oxdna_source_archive(cache)
    build = ROOT / ".local_deps/oxdna-build"
    build.mkdir(parents=True, exist_ok=True)
    source = build / f"oxDNA-{OXDNA_REVISION}"
    if not source.is_dir():
        with tarfile.open(archive) as package:
            package.extractall(build, filter="data")
    import runpy
    runpy.run_path(str(ROOT / "scripts/prepare_oxdna_ffs_runtime.py"))["repair_source"](source)
    subprocess.run(["cmake", "-S", str(source), "-B", str(source / "build"), "-DPython=ON",
                    "-DOxpySystemInstall=ON", "-DCUDA=OFF", "-DNATIVE_COMPILATION=OFF", "-DCMAKE_BUILD_TYPE=Release",
                    f"-DDEFAULT_PYTHON={sys.executable}"], check=True)
    subprocess.run(["cmake", "--build", str(source / "build"), "--target", "oxDNA", "core", "--parallel", str(jobs)], check=True)
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-build-isolation", "--force-reinstall", str(source / "build/python")], check=True)
    with tempfile.TemporaryDirectory(prefix="aptaswitch-oxdna-freeze-") as directory:
        stage = Path(directory)
        command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir",
                   "--name", "aptaswitch-oxdna-worker", "--distpath", str(stage / "dist"),
                   "--workpath", str(stage / "work"), "--specpath", str(stage),
                   "--paths", str(ROOT / "src/aptaswitch_core"), "--hidden-import", "oxpy",
                   "--collect-all", "oxpy", "--collect-all", "numpy",
                   "--add-binary", str(source / "build/bin/oxDNA") + ":bin"]
        for worker in WORKERS:
            command += ["--hidden-import", worker]
        freeze_environment = dict(os.environ, PYINSTALLER_CONFIG_DIR=str(stage / "cache"))
        subprocess.run([*command, str(ROOT / "scripts/run_oxdna_worker.py")], check=True,
                       env=freeze_environment)
        output.mkdir(parents=True, exist_ok=True)
        if (output / "worker").exists():
            shutil.rmtree(output / "worker")
        shutil.copytree(stage / "dist/aptaswitch-oxdna-worker", output / "worker", symlinks=True)
    (output / "parameters").mkdir(exist_ok=True)
    for name in PARAMETERS:
        shutil.copyfile(source / name, output / "parameters" / name)
    (output / "sources").mkdir(exist_ok=True)
    shutil.copyfile(archive, output / "sources/oxDNA-source.tar.gz")
    shutil.copyfile(source / "LICENSE", output / "LICENSE.txt")
    for name in ("prepare_bundled_science.py", "prepare_oxdna_ffs_runtime.py", "run_oxdna_worker.py", "rebuild_oxdna_worker.py"):
        shutil.copyfile(ROOT / "scripts" / name, output / "sources" / name)
    for name in WORKERS:
        shutil.copyfile(ROOT / "src/aptaswitch_core" / (name + ".py"), output / "sources" / (name + ".py"))
    (output / "sources/BUILDING.txt").write_text(
        "oxDNA is GPLv3. This payload includes its exact upstream source archive,\n"
        "Apta2Switch's force-reset repair, build script and Python worker sources.\n"
        "Copy this entire sources folder to a writable location. No application checkout is needed.\n"
        "Install Python 3.11, CMake, and a C++ compiler (Xcode command-line tools on macOS, GCC on Linux).\n"
        "In the copied sources folder, create a disposable build environment:\n"
        "python3.11 -m venv rebuild-env\n"
        "rebuild-env/bin/python -m pip install numpy==1.26.4 pyinstaller==6.21.0 setuptools==80.10.2\n"
        "rebuild-env/bin/python rebuild_oxdna_worker.py --verify-sources\n"
        "rebuild-env/bin/python rebuild_oxdna_worker.py --output ../rebuilt-oxdna --jobs 2\n"
        "The script uses only the shipped original source archive, repair and worker/helper files.\n"
        "It enables Python/oxpy, disables CUDA and NATIVE_COMPILATION, and uses Release CPU builds.\n"
        "It builds a standalone worker, copies parameters/notices/sources, writes hashes, then runs the NA FFS probe.\n"
        "Package installation may need Internet; no scientific source/model download is performed.\n"
        "Code-signing is platform-specific and is not reproduced by this source rebuild.\n"
        f"Upstream revision: {OXDNA_REVISION}\n")
    patch = source / "aptaswitch_ffs_force_reset.json"
    if patch.exists():
        shutil.copyfile(patch, output / "sources/aptaswitch_ffs_force_reset.json")
    dependency_notices(output)
    manifest(output, "oxdna", worker="worker/aptaswitch-oxdna-worker", binary="worker/_internal/bin/oxDNA",
        parameters={kind: "parameters/" + name for kind, name in zip(("DNA", "RNA", "NA"), PARAMETERS)},
        source_revision=OXDNA_REVISION, source_archive_sha256=OXDNA_ARCHIVE_SHA256,
        license="GPL-3.0", system=platform.system(), machine=platform.machine())
    validate_payload(output, "oxdna")
    from aptaswitch_core.oxdna_setup import validate_oxdna_runtime
    validate_oxdna_runtime(str(output))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--component", choices=("all", "oxdna"), default="all")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if sys.platform not in ("darwin", "linux"):
        raise SystemExit("Build native payloads on macOS or Linux; Windows uses the complete Linux payload.")
    cache = ROOT / ".local_deps/bundled-science-downloads"
    for component in (("oxdna",) if args.component == "all" else (args.component,)):
        output = ROOT / ".local_deps" / (component + "-runtime")
        if args.check:
            validate_payload(output, component)
        else:
            prepare_oxdna(cache, output, args.jobs)


if __name__ == "__main__":
    main()
