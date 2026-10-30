#!/usr/bin/env python3
"""Maintainer build of the offline hybrid runtime, with pinned Vienna sources.

End users receive the resulting library; no compiler or download is needed at
application launch. Windows uses this same Linux payload through its WSL runtime.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
VERSION = "2.7.2"
URL = f"https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-{VERSION}.tar.gz"
ARCHIVE_SHA256 = "1ab5f4a4f76fc85a2243546088e45f5d85f2d7a56cc656e969b005cce9bfab5f"
ORIGINAL_PARAMS = "a6da5b5038e07c6b4b0b843429115db84b6bcb8b48e51074f3180d2842358e7a"
PATCHED_PARAMS = "0163c59586c5389acc0b1284a97fda7fd0652716b1d26548466e5df735956cf9"
NOTICES = """Apta2Switch experimental hybrid thermodynamic runtime

ViennaRNA Package 2.7.2 — Institute for Theoretical Chemistry,
University of Vienna. Authors and contributors: https://www.tbi.univie.ac.at/RNA/
Its license is reproduced in ViennaRNA-COPYING.txt. The library is a derived work;
ViennaRNA's redistribution conditions apply independently of the application's
custom LicenseRef-Apta2Switch-Source license; ViennaRNA itself is not relicensed
by that grant. That custom license does not override any third-party GPL terms.
In particular, distribution for a fee requires consulting the
ViennaRNA authors. Source versions and parameter hashes are recorded below.

The isolated motif parameter provider uses ViennaRNA's Turner 2004 RNA and
Mathews 2004 DNA tables, and the public Sugimoto et al. (1995) RNA/DNA hybrid
stacking enthalpies and entropies (Biochemistry 34:11211–11216,
doi:10.1021/bi00035a029). Hybrid loop energies are an experimental composition
of homogeneous loop models; they are not an official ViennaRNA hybrid model.
No dangles or coaxial stacking are modeled. See the source model documentation.

Single RNA strands use RNAlib 2.7.2 for partition functions, pair probabilities,
joint samples and minimum free-energy structures, with the same RNA parameters
and declared sodium policy. RNAlib minima use centikcal energies and a rounded
per-pair sodium correction; hybrid minima retain the native PF-weight convention.
These solver and minimum-energy conventions are recorded in the build manifest.

Rebuild: scripts/prepare_thermodynamics_runtime.py. Source archive and hashes:
https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.2.tar.gz
The only Vienna source adjustment initializes SaltLoop index MAXLOOP+1 in
params.c (four loops). Exact source, parameter, static-library and runtime hashes
are recorded in the accompanying build manifest.
"""


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def third_party_notices() -> str:
    """Keep the incorporated Biopython table's readable binary notice intact."""
    source = ROOT / "src/aptaswitch_core/thermodynamics/hybrid_motifs.py"
    documentation = ast.get_docstring(ast.parse(source.read_text()))
    marker = "The small H/S table below is from Biopython's MeltingTemp.py:"
    if documentation is None or marker not in documentation:
        raise RuntimeError("Missing Biopython attribution for incorporated hybrid parameters")
    return NOTICES + "\n" + documentation[documentation.index(marker):] + "\n"


def validate_payload(root: Path) -> dict:
    root = Path(root)
    candidates = [root / name for name in ("libaptaswitch_hybrid.dylib", "libaptaswitch_hybrid.so", "aptaswitch_hybrid.dll")
                  if (root / name).is_file()]
    if len(candidates) != 1:
        raise RuntimeError("Exactly one platform-matching thermodynamic library is required")
    library = candidates[0]
    manifest = json.loads(library.with_suffix(library.suffix + ".json").read_text())
    if manifest.get("native_library_sha256") != sha256(library) or manifest.get("status") != "built":
        raise RuntimeError("Thermodynamic library does not match its build manifest")
    for name in ("ViennaRNA-COPYING.txt", "THIRD-PARTY-NOTICES.txt"):
        if not (root / name).is_file() or not (root / name).read_text().strip():
            raise RuntimeError(f"Missing thermodynamics license notice: {name}")
    return manifest


def prepare_source(cache: Path, archive: Path | None = None) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    # The wiki source aggregate carries the exact unmodified upstream archive.
    # Retain a verified cache copy for the OSTIR corresponding-source bundle.
    vendored = ROOT / "vendor/sources" / f"ViennaRNA-{VERSION}.tar.gz"
    if archive is None and vendored.is_file():
        archive = vendored
    if archive is not None:
        if sha256(archive) != ARCHIVE_SHA256:
            raise RuntimeError("ViennaRNA source archive checksum mismatch")
        cached = cache / f"ViennaRNA-{VERSION}.tar.gz"
        if Path(archive).resolve() != cached.resolve():
            shutil.copyfile(archive, cached)
        archive = cached
    if archive is None:
        archive = cache / f"ViennaRNA-{VERSION}.tar.gz"
        if not archive.is_file():
            pending = archive.with_suffix(".download")
            with urllib.request.urlopen(URL, timeout=60) as response, pending.open("wb") as output:
                shutil.copyfileobj(response, output)
            if sha256(pending) != ARCHIVE_SHA256:
                pending.unlink()
                raise RuntimeError("ViennaRNA source archive checksum mismatch")
            pending.replace(archive)
    if sha256(archive) != ARCHIVE_SHA256:
        raise RuntimeError("ViennaRNA source archive checksum mismatch")
    source = cache / f"ViennaRNA-{VERSION}"
    if not source.exists():
        with tarfile.open(archive, "r:gz") as package:
            if any(member.name.split('/')[0] != source.name for member in package.getmembers()):
                raise RuntimeError("Unexpected source archive layout")
            package.extractall(cache, filter="data")
    params = source / "src/ViennaRNA/params/params.c"
    if sha256(params) == ORIGINAL_PARAMS:
        value = params.read_text()
        needle = "for (i = 31; i <= MAXLOOP; i++) {"
        if value.count(needle) != 4:
            raise RuntimeError("Unexpected Vienna salt-loop source")
        params.write_text(value.replace(needle, "for (i = 31; i <= MAXLOOP + 1; i++) {"))
    if sha256(params) != PATCHED_PARAMS:
        raise RuntimeError("Vienna parameter source differs from the audited correction")
    configured = cache / "native-configured.json"
    build_flags = dict(CC=os.environ.get("CC", "cc"), CXX=os.environ.get("CXX", "c++"),
                       CFLAGS="-O2 -fPIC", CXXFLAGS="-O2 -fPIC")
    # Only the C library is used here; the Python API is supplied by its pinned
    # wheel. Disabling individual bindings alone still requires SWIG >= 4.3.
    configure_options = ["--without-swig", "--without-python", "--without-perl", "--without-doc",
                         "--without-kinfold", "--without-rnalocmin", "--without-rnaxplorer",
                         "--without-svm", "--without-gsl", "--disable-openmp", "--disable-lto",
                         "--disable-simd", "--disable-dependency-tracking"]
    identity = {**build_flags, "configure_options": configure_options,
                "params": PATCHED_PARAMS, "machine": platform.machine(), "platform": sys.platform}
    if not configured.exists() or json.loads(configured.read_text()) != identity:
        # Make does not know that compiler flags changed. Do not reuse objects
        # produced with a different toolchain/configuration under a new recipe.
        if (source / "src/ViennaRNA/Makefile").exists():
            subprocess.run(["make", "-C", "src/ViennaRNA", "clean"], cwd=source, check=True)
        subprocess.run(["./configure", *configure_options], cwd=source,
                       env={**os.environ, **build_flags}, check=True)
        configured.write_text(json.dumps(identity))
    subprocess.run(["make", "-C", "src/ViennaRNA", "-j2"], cwd=source, check=True)
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vienna-source", type=Path, help="Already compiled, audited Vienna 2.7.2 source tree")
    parser.add_argument("--archive", type=Path, help="Pinned source archive for an offline maintainer build")
    parser.add_argument("--output", type=Path, default=ROOT / ".local_deps/thermodynamics")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        validate_payload(args.output)
        print("Thermodynamic payload verified")
        return
    if sys.platform not in ("darwin", "linux"):
        raise SystemExit("Build on macOS or Linux; Windows ships the Linux runtime through WSL.")
    source = args.vienna_source or prepare_source(ROOT / ".local_deps/thermodynamics-build", args.archive)
    library = args.output / ("libaptaswitch_hybrid.dylib" if sys.platform == "darwin" else "libaptaswitch_hybrid.so")
    subprocess.run([sys.executable, str(ROOT / "scripts/build_hybrid_native.py"),
                    "--vienna-source", str(source), "--output", str(library),
                    "--compiler", os.environ.get("CXX", "c++")], check=True)
    shutil.copyfile(source / "COPYING", args.output / "ViennaRNA-COPYING.txt")
    (args.output / "THIRD-PARTY-NOTICES.txt").write_text(third_party_notices())
    validate_payload(args.output)
    sys.path.insert(0, str(ROOT / "src"))
    from aptaswitch_core.thermodynamics.setup import validate_runtime
    checked = validate_runtime(library)
    if not checked.ok:
        raise RuntimeError(checked.message)
    print(checked.message)


if __name__ == "__main__":
    main()
