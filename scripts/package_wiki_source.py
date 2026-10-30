#!/usr/bin/env python3
"""Build and verify the wiki source download; never include local user data.

The upstream archives are separate source packages with their original terms.
This is a source aggregate, not a redistributable combined executable or an
offline Python installer. Run with Python 3.11 or newer; no installed app needed.
"""
from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
APPLICATION_LICENSE = "LicenseRef-Apta2Switch-Source"
DIRECTORIES = ("src", "scripts", "packaging", "tests", "LICENSES", ".github", ".streamlit")
FILES = ("LICENSE", "THIRD_PARTY_NOTICES.md", "README.md", "pyproject.toml", "MANIFEST.in",
         ".gitignore", ".gitattributes", "requirements-constraints.txt", "UPSTREAM-SOURCES.json")
EXCLUDED = {"__pycache__", "build", "dist", "runs", "designs", ".git", ".local_deps",
            ".venv", ".build-venv", ".aws", ".ssh", ".DS_Store", "credentials", "secrets.toml",
            "validation-results", ".pytest_cache"}
SUFFIXES = {".py", ".sh", ".ps1", ".spec", ".md", ".txt", ".toml", ".json", ".jsonl",
            ".yml", ".yaml", ".patch", ".S", ".c", ".cpp", ".h", ".hpp", ".in", ".par",
            ".js", ".svg", ".png", ".jpg", ".jpeg", ".icns", ".ico", ".plist", ".html", ".css", ".csv"}
REQUIRED_FILES = {
    "LICENSE", "LICENSES/MIT-legacy.txt", "THIRD_PARTY_NOTICES.md", "README.md", "pyproject.toml",
    "scripts/build_macos.sh", "scripts/build_linux.sh",
    "scripts/build_environment.sh", "scripts/build_hybrid_native.py", "scripts/package_wiki_source.py",
    "src/aptaswitch_core/thermodynamics/hybrid_motifs.py",
    "tests/native_reference/reference_checks.py", "tests/fixtures/extension_native_fixture.json",
    *(f"src/aptaswitch_core/thermodynamics/native/{name}" for name in (
        "connected_dp.cpp", "connected_dp.hpp", "hybrid_api.cpp", "hybrid_api.hpp",
        "hybrid_provider.cpp", "hybrid_provider.hpp")),
}
UPSTREAM = (
    dict(name="ViennaRNA", version="2.7.2", filename="ViennaRNA-2.7.2.tar.gz",
         sha256="1ab5f4a4f76fc85a2243546088e45f5d85f2d7a56cc656e969b005cce9bfab5f",
         url="https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.2.tar.gz",
         license="ViennaRNA custom license", notice="ViennaRNA-2.7.2/COPYING",
         cache=(".local_deps/thermodynamics-build/ViennaRNA-2.7.2.tar.gz", ".local_deps/ostir-sources/ViennaRNA-2.7.2.tar.gz")),
    dict(name="oxDNA", version="75ae96f063954b2f29524bb897c274aee7ba4466",
         filename="oxDNA-75ae96f063954b2f29524bb897c274aee7ba4466.tar.gz",
         sha256="29e7780829ee967046168e55851a1e51017c6f37daeed3b8550b8d855704c878",
         url="https://codeload.github.com/lorenzo-rovigatti/oxDNA/tar.gz/75ae96f063954b2f29524bb897c274aee7ba4466",
         license="GPL-3.0", notice="oxDNA-75ae96f063954b2f29524bb897c274aee7ba4466/LICENSE",
         cache=(".local_deps/bundled-science-downloads/oxdna-75ae96f063954b2f29524bb897c274aee7ba4466.tar.gz",
                ".local_deps/oxdna-runtime-beta2-md-20261014/sources/oxDNA-source.tar.gz",
                ".local_deps/oxdna-runtime/sources/oxDNA-source.tar.gz")),
    dict(name="OSTIR", version="1.1.3", filename="ostir-1.1.3.tar.gz",
         sha256="6fa38a4e867c8d38aab16b59addc07ce2de03ee766ef6018c24c7a78ff163e14",
         url="https://files.pythonhosted.org/packages/source/o/ostir/ostir-1.1.3.tar.gz",
         license="GPL-3.0-or-later", notice="ostir-1.1.3/LICENSE",
         cache=(".local_deps/ostir-sources/ostir-1.1.3.tar.gz",)),
)


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def version(root: Path) -> str:
    import tomllib
    return tomllib.loads((root / "pyproject.toml").read_text())["project"]["version"]


def check_build_pins(root: Path) -> None:
    checks = (("scripts/prepare_thermodynamics_runtime.py", {"VERSION": UPSTREAM[0]["version"], "ARCHIVE_SHA256": UPSTREAM[0]["sha256"]}),
              ("scripts/prepare_bundled_science.py", {"OXDNA_REVISION": UPSTREAM[1]["version"], "OXDNA_ARCHIVE_SHA256": UPSTREAM[1]["sha256"]}),
              ("scripts/prepare_ostir_sources.py", {"VERSION": UPSTREAM[2]["version"], "SHA256": UPSTREAM[2]["sha256"]}))
    for filename, expected in checks:
        values = {}
        for node in ast.parse((root / filename).read_text()).body:
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        values[target.id] = node.value.value
        if any(values.get(key) != value for key, value in expected.items()):
            raise RuntimeError(f"Source-release pins differ from the actual builder: {filename}")


def selected_files(root: Path) -> list[Path]:
    paths = [root / name for name in FILES if (root / name).is_file() and not (root / name).is_symlink()]
    for directory in DIRECTORIES:
        if (root / directory).is_symlink():
            continue
        for base, directories, names in os.walk(root / directory, followlinks=False):
            directories[:] = sorted(name for name in directories if name not in EXCLUDED
                                    and not name.endswith(".egg-info") and not (Path(base) / name).is_symlink()
                                    and name != "results")
            for name in sorted(names):
                path = Path(base) / name
                if (name in EXCLUDED or name.startswith(".env") or path.is_symlink()
                        or path.suffix not in SUFFIXES):
                    continue
                paths.append(path)
    # macOS iCloud placeholders can block reads indefinitely. Do not silently
    # omit them or replace potentially newer source with an old release copy.
    placeholders = [str(p.relative_to(root)) for p in paths if getattr(p.stat(), "st_flags", 0) & 0x40000000]
    if placeholders:
        raise RuntimeError("Download these iCloud source files before packaging:\n" + "\n".join(placeholders))
    return sorted(paths)


def check_required_files(files) -> None:
    missing = REQUIRED_FILES - set(files)
    if missing:
        raise RuntimeError("Source release is missing build inputs or license notices: "
                           + ", ".join(sorted(missing)))


def upstream_archive(root: Path, entry: dict, offline: bool) -> Path:
    candidates = [root / "vendor/sources" / entry["filename"], *(root / name for name in entry["cache"])]
    target = root / ".local_deps/wiki-source-cache" / entry["filename"]
    candidates.append(target)
    for path in candidates:
        if path.is_file():
            if digest(path) != entry["sha256"]:
                raise RuntimeError(f"Source checksum mismatch: {path}")
            return path
    if offline:
        raise RuntimeError(f"Missing upstream source in offline mode: {entry['filename']}")
    target.parent.mkdir(parents=True, exist_ok=True)
    pending = target.with_suffix(".download")
    try:
        with urllib.request.urlopen(entry["url"], timeout=60) as response, pending.open("wb") as output:
            shutil.copyfileobj(response, output)
        if digest(pending) != entry["sha256"]:
            raise RuntimeError(f"Downloaded source checksum mismatch: {entry['filename']}")
        pending.replace(target)
    finally:
        pending.unlink(missing_ok=True)
    return target


def add_bytes(archive: tarfile.TarFile, name: str, value: bytes, mode: int = 0o644) -> str:
    info = tarfile.TarInfo(name)
    info.size, info.mode = len(value), mode
    archive.addfile(info, io.BytesIO(value))
    return hashlib.sha256(value).hexdigest()


def build(root: Path, output: Path, *, offline: bool = False) -> dict:
    root = root.resolve()
    check_build_pins(root)
    paths = selected_files(root)
    check_required_files(path.relative_to(root).as_posix() for path in paths)
    app_version = version(root)
    prefix = f"Apta2Switch-Studio-{app_version}"
    upstream = [(entry, upstream_archive(root, entry, offline)) for entry in UPSTREAM]
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema_version=1, project="Apta2Switch-Studio", version=app_version,
                    application_license=APPLICATION_LICENSE, delivery="source-aggregate",
                    build_targets=["macOS arm64", "macOS x86_64", "Linux x86_64"],
                    python="3.11",
                    feature_status={
                        "oxdna_ffs": {"complete_run_completed": False,
                                      "end_to_end_tested": False,
                                      "reason": "No complete run was performed due to time constraints."},
                        "full_reporter_context": {"functional": False, "end_to_end_tested": False},
                    }, upstream=[], files={})
    with tempfile.NamedTemporaryFile(dir=output.parent, suffix=".pending", delete=False) as temporary:
        pending = Path(temporary.name)
    try:
        with pending.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed, tarfile.open(fileobj=compressed, mode="w|") as archive:
            for path in paths:
                relative = path.relative_to(root).as_posix()
                before = path.stat()
                data = path.read_bytes()
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise RuntimeError(f"Source changed while packaging: {relative}")
                mode = 0o755 if before.st_mode & 0o111 else 0o644
                manifest["files"][relative] = add_bytes(archive, prefix + "/" + relative, data, mode)
            for entry, path in upstream:
                relative = "vendor/sources/" + entry["filename"]
                with path.open("rb") as source:
                    info = tarfile.TarInfo(prefix + "/" + relative)
                    info.size, info.mode = path.stat().st_size, 0o644
                    archive.addfile(info, source)
                manifest["files"][relative] = entry["sha256"]
                with tarfile.open(path) as package:
                    stream = package.extractfile(entry["notice"])
                    if stream is None:
                        raise RuntimeError(f"Missing upstream license: {entry['name']}")
                    notice = stream.read()
                license_path = "vendor/licenses/" + entry["name"] + "-LICENSE.txt"
                manifest["files"][license_path] = add_bytes(archive, prefix + "/" + license_path, notice)
                manifest["upstream"].append({key: entry[key] for key in ("name", "version", "url", "sha256", "license")} | {"path": relative, "license_path": license_path})
            add_bytes(archive, prefix + "/RELEASE-MANIFEST.json", (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
        check(pending)
        pending.chmod(0o644)
        pending.replace(output)
    finally:
        pending.unlink(missing_ok=True)
    output.with_name(output.name + ".sha256").write_text(digest(output) + "  " + output.name + "\n")
    return manifest


def check(path: Path) -> dict:
    with tarfile.open(path) as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        if len(names) != len(set(names)):
            raise RuntimeError("Duplicate source archive path")
        for member in members:
            parts = PurePosixPath(member.name).parts
            if (not member.isfile() or member.name.startswith("/") or "\\" in member.name
                    or PurePosixPath(member.name).as_posix() != member.name
                    or ".." in parts or len(parts) < 2 or not parts[0].startswith("Apta2Switch-Studio-")):
                raise RuntimeError(f"Unexpected source archive entry: {member.name}")
            if any(part in EXCLUDED or part.endswith(".egg-info") or part.startswith(".env") for part in parts):
                raise RuntimeError(f"Private/generated content in source archive: {member.name}")
            if "results" in parts[1:]:
                raise RuntimeError(f"Unexpected results directory: {member.name}")
        roots = {PurePosixPath(name).parts[0] for name in names}
        if len(roots) != 1:
            raise RuntimeError("Expected exactly one source directory")
        prefix = roots.pop() + "/"
        manifest = json.load(archive.extractfile(prefix + "RELEASE-MANIFEST.json"))
        if manifest.get("application_license") != APPLICATION_LICENSE or manifest.get("delivery") != "source-aggregate":
            raise RuntimeError("Invalid source release metadata")
        files = manifest["files"]
        if set(names) != {prefix + name for name in files} | {prefix + "RELEASE-MANIFEST.json"}:
            raise RuntimeError("Archive and manifest paths differ")
        for name, expected in files.items():
            stream = archive.extractfile(prefix + name)
            if hashlib.file_digest(stream, "sha256").hexdigest() != expected:
                raise RuntimeError(f"Source file differs from manifest: {name}")
        check_required_files(files)
        expected_upstream = {(entry["name"], entry["sha256"]) for entry in UPSTREAM}
        if {(entry["name"], entry["sha256"]) for entry in manifest["upstream"]} != expected_upstream:
            raise RuntimeError("Upstream source pins differ from the audited release")
        for entry in manifest["upstream"]:
            if files.get(entry["path"]) != entry["sha256"] or entry["license_path"] not in files:
                raise RuntimeError("Missing upstream source or license notice")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--offline", action="store_true", help="Require existing, hash-verified upstream archives")
    parser.add_argument("--check", type=Path, metavar="ARCHIVE", help="Verify an existing source download")
    args = parser.parse_args()
    if args.check:
        manifest = check(args.check)
        print(f"Verified {len(manifest['files'])} source files: {args.check}")
        return
    output = args.output or ROOT / "dist" / f"Apta2Switch-Studio-{version(ROOT)}-wiki-source.tar.gz"
    manifest = build(ROOT, output, offline=args.offline)
    print(f"Source archive: {output}\nSHA-256: {digest(output)}\nFiles: {len(manifest['files'])}")


if __name__ == "__main__":
    main()
