"""Read-only discovery and integrity checks for shipped scientific payloads."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys


def runtime_directory(name: str) -> Path:
    if getattr(sys, "frozen", False):
        root = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
        candidates = (root / name, root.parent / "Resources" / name)
    else:
        candidates = (Path(__file__).resolve().parents[2] / ".local_deps" / name.replace("_", "-"),)
    return next((path.resolve() for path in candidates if path.is_dir()), candidates[0].resolve())


def payload_manifest(root: Path, expected_kind: str) -> dict:
    root = Path(root).resolve()
    try:
        manifest = json.loads((root / "runtime.json").read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f"Bundled {expected_kind} runtime is missing or damaged; reinstall the complete application") from error
    if (not isinstance(manifest, dict) or manifest.get("schema_version") != 1
            or manifest.get("kind") != expected_kind or not isinstance(manifest.get("files"), dict)
            or not manifest["files"]):
        raise ValueError(f"Invalid bundled {expected_kind} manifest")
    for name, expected in manifest["files"].items():
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or "\\" in name:
            raise ValueError("Unsafe bundled runtime path")
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError(f"Missing bundled runtime file: {name}")
        with path.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if actual != expected:
            raise ValueError(f"Bundled runtime checksum mismatch: {name}")
    return manifest
