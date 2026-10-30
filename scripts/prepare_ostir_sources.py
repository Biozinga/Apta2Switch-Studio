#!/usr/bin/env python3
"""Attach exact OSTIR and application sources to every frozen local build.

Downloads occur only at build time. This records corresponding source and
notices; it is not a determination of GPL/ViennaRNA licence compatibility.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import tarfile
import urllib.request

VERSION = "1.1.3"
SHA256 = "6fa38a4e867c8d38aab16b59addc07ce2de03ee766ef6018c24c7a78ff163e14"
URL = f"https://files.pythonhosted.org/packages/source/o/ostir/ostir-{VERSION}.tar.gz"
SOURCE_DIRS = ("src", "scripts", "packaging", "tests", ".github", "LICENSES", "vendor/licenses")
SOURCE_FILES = ("LICENSE", "THIRD_PARTY_NOTICES.md", "README.md", "pyproject.toml", "MANIFEST.in",
                ".gitignore", ".gitattributes", "requirements-constraints.txt", "UPSTREAM-SOURCES.json",
                ".streamlit/config.toml")
EXCLUDED = {"__pycache__", ".DS_Store", "build", "dist", "results", ".env", ".aws", ".ssh",
            ".venv", ".local_deps", "credentials", "secrets.toml", "validation-results", ".pytest_cache"}
NOTICE = """OSTIR 1.1.3 — Cameron Roots, Alexandra Lukasiewicz, Jeffrey Barrick and contributors
https://github.com/barricklab/ostir
GNU GPL version 3 or later; full text in OSTIR-LICENSE.txt.

OSTIR itself is unmodified. Apta2Switch supplies an exploratory constrained-RNA
ON/OFF adaptation; its source is included in Apta2Switch-source.tar.gz.
Apta2Switch-Studio source uses its custom Apta2Switch-Studio Source License
(LicenseRef-Apta2Switch-Source): use and modification are permitted, including
commercial use; redistribution must be free apart from physical-media costs,
with retained attribution and notices. Its full LICENSE accompanies the source,
with the MIT notice from earlier releases preserved in LICENSES/.
Third-party works keep their own licenses; this change does not relicense them.

ViennaRNA 2.7.2 has separate redistribution terms, including a fee restriction.
Including sources and these notices does NOT resolve the compatibility of those
terms with the GPL. Our own custom license does not override OSTIR GPL terms
or resolve that compatibility. Local use is the scope of this integration.
Distribution of this combined build requires resolving that compatibility
with the licensors;
neither a separate file nor a separate process is assumed to grant an exemption.

Source rebuild: unpack Apta2Switch-source.tar.gz, create a Python 3.11 virtual
environment, install -e '.[build]' and use scripts/build_macos.sh or build_linux.sh.
The exact upstream OSTIR and ViennaRNA sources are included alongside this file;
other scientific runtime sources accompany their own payloads. Build scripts
record additional dependencies and their source locations. Application source
contains no designs, user runs, local settings or credentials.
"""


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_application_source(root: Path, target: Path) -> None:
    """Snapshot actual build inputs, including uncommitted source changes."""
    files = [root / name for name in SOURCE_FILES
             if (root / name).is_file() and not (root / name).is_symlink()]
    for name in SOURCE_DIRS:
        for path in (root / name).rglob("*"):
            relative = path.relative_to(root)
            if (not path.is_file() or path.is_symlink()
                    or any(part in EXCLUDED or part.endswith(".egg-info")
                           or part.startswith(".env.") for part in relative.parts)
                    or path.suffix in {".pyc", ".pyo"}):
                continue
            files.append(path)
    with tarfile.open(target, "w:gz") as archive:
        for path in sorted(files):
            info = archive.gettarinfo(str(path), "Apta2Switch-Studio/" + path.relative_to(root).as_posix())
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            with path.open("rb") as stream:
                archive.addfile(info, stream)


def prepare_sources(root: Path, archive: Path | None = None) -> Path:
    root = Path(root)
    for name, version in (("OSTIR", VERSION), ("ViennaRNA", "2.7.2")):
        if importlib.metadata.version(name) != version:
            raise RuntimeError(f"Install the pinned {name}=={version} before freezing")
    output = root / ".local_deps/ostir-sources"
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"ostir-{VERSION}.tar.gz"
    if archive is None:
        vendored = root / "vendor/sources" / target.name
        if vendored.is_file():
            archive = vendored
    if archive is not None:
        if digest(Path(archive)) != SHA256:
            raise RuntimeError("OSTIR source checksum mismatch")
        if Path(archive).resolve() != target.resolve():
            shutil.copyfile(archive, target)
    elif not target.is_file():
        pending = target.with_suffix(".download")
        with urllib.request.urlopen(URL, timeout=60) as response, pending.open("wb") as stream:
            shutil.copyfileobj(response, stream)
        if digest(pending) != SHA256:
            pending.unlink()
            raise RuntimeError("OSTIR download checksum mismatch")
        pending.replace(target)
    if digest(target) != SHA256:
        raise RuntimeError("OSTIR source checksum mismatch")
    distribution = importlib.metadata.distribution("OSTIR")
    with tarfile.open(target) as package:
        license_bytes = package.extractfile(f"ostir-{VERSION}/LICENSE").read()
        # A pinned version alone would not detect local modifications to OSTIR.
        for entry in distribution.files or ():
            if entry.parts[0] == "ostir" and entry.suffix in {".py", ".par"}:
                upstream = package.extractfile(f"ostir-{VERSION}/{entry.as_posix()}").read()
                if upstream != Path(distribution.locate_file(entry)).read_bytes():
                    raise RuntimeError(f"Installed OSTIR differs from its source archive: {entry}")
    (output / "OSTIR-LICENSE.txt").write_bytes(license_bytes)
    # Reuse the authenticated source already required by the native engine build.
    import runpy
    vienna = runpy.run_path(str(root / "scripts/prepare_thermodynamics_runtime.py"))
    vienna_archive = root / ".local_deps/thermodynamics-build/ViennaRNA-2.7.2.tar.gz"
    if not vienna_archive.is_file() or digest(vienna_archive) != vienna["ARCHIVE_SHA256"]:
        raise RuntimeError("Prepare the authenticated ViennaRNA source before freezing OSTIR")
    shutil.copyfile(vienna_archive, output / vienna_archive.name)
    (output / "NOTICE.txt").write_text(NOTICE, encoding="utf-8")
    write_application_source(root, output / "Apta2Switch-source.tar.gz")
    manifest = {
        "ostir_version": VERSION, "viennarna_version": "2.7.2", "license": "GPL-3.0-or-later",
        "source_url": URL, "source_sha256": SHA256,
        "files": {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != "manifest.json"},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    print(prepare_sources(Path(__file__).resolve().parents[1], args.archive))
