"""Automatic discovery and internal validation of bundled oxDNA/oxpy."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from functools import lru_cache

from .localization import tr
from .bundled_science import payload_manifest, runtime_directory
from .scientific_process import external_process_environment


PARAMETERS = {"DNA": "oxDNA2_sequence_dependent_parameters.txt",
              "RNA": "rna_sequence_dependent_parameters.txt",
              "NA": "NA_sequence_dependent_parameters.txt"}


def default_oxdna_directory() -> Path:
    return runtime_directory("oxdna_runtime")


def configured_oxdna_home() -> str:
    return str(default_oxdna_directory())


def _probe_oxpy(python: str) -> str:
    code = ("import oxpy; "
            "assert all(hasattr(oxpy, name) for name in ('Context', 'InputFile', 'OxpyManager')); "
            "print(getattr(oxpy, '__version__', 'available'))")
    try:
        result = subprocess.run([python, "-c", code], capture_output=True, text=True,
                                timeout=15, env=external_process_environment())
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValueError(tr("Le Python local oxDNA ne peut pas charger oxpy.")) from exc
    if result.returncode:
        raise ValueError(tr("Le Python local oxDNA ne peut pas charger oxpy."))
    return result.stdout.strip().splitlines()[-1] if result.stdout.strip() else "available"


def _probe_native_ffs(runtime: dict) -> dict:
    """Validate forces and native stop conditions, not just the Python import."""
    script = Path(__file__).with_name("oxdna_native_ffs.py")
    try:
        with tempfile.TemporaryDirectory(prefix="aptaswitch_oxdna_check_") as directory:
            result = subprocess.run(
                [runtime["python"], str(script), "--probe", json.dumps(runtime)],
                capture_output=True, text=True, timeout=30, cwd=directory,
                env=external_process_environment())
        if result.returncode:
            raise ValueError(result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "FFS_MD probe failed")
        checked = json.loads(result.stdout.strip().splitlines()[-1])
        if not checked.get("force_equivalence") or not checked.get("native_first_crossing"):
            raise ValueError("FFS_MD probe did not validate physical forces and native stops")
        return checked
    except (OSError, subprocess.TimeoutExpired, ValueError, IndexError) as exc:
        raise ValueError("The bundled oxDNA FFS engine failed its force/stop integrity check; reinstall the complete application.") from exc


def validate_oxdna_runtime(folder: str) -> dict:
    """Internal build/smoke check: executables, oxpy and all three NA models."""
    if not isinstance(folder, str) or not folder.strip():
        raise ValueError(tr("Indiquez le dossier oxDNA à vérifier."))
    root = Path(folder).expanduser().resolve()
    manifest = payload_manifest(root, "oxdna")
    for key in ("worker", "binary"):
        if manifest.get(key) not in manifest["files"]:
            raise ValueError("Bundled oxDNA executable provenance is incomplete")
    python, binary = root / manifest["worker"], root / manifest["binary"]
    if not python.is_file() or not binary.is_file():
        raise ValueError(tr("Installation oxDNA incomplète : oxpy et oxDNA sont nécessaires."))
    if os.name != "nt" and any(not os.access(path, os.X_OK) for path in (python, binary)):
        raise ValueError(tr("Le Python local et le moteur oxDNA doivent avoir la permission d’exécution."))
    runtime = {"home": str(root), "python": str(python), "binary": str(binary),
               "source_revision": manifest.get("source_revision"), "license": "GPL-3.0", "bundled": True}
    hashes = {}
    for kind, filename in PARAMETERS.items():
        relative = manifest.get("parameters", {}).get(kind)
        if relative not in manifest["files"]:
            raise ValueError("Bundled oxDNA parameter provenance is incomplete")
        path = root / relative
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(tr("Paramètres oxNA dépendants de la séquence absents : {file}.", file=filename))
        runtime["seq_dep_file_" + kind] = str(path)
        hashes[kind] = hashlib.sha256(path.read_bytes()).hexdigest()
    runtime["parameter_sha256"] = hashes
    runtime["oxpy_version"] = _probe_oxpy(str(python))
    runtime["native_ffs"] = _probe_native_ffs(runtime)
    return runtime


@lru_cache(maxsize=4)
def _checked_runtime(folder, manifest_stamp):
    return validate_oxdna_runtime(folder)


def detect_oxdna_runtime() -> dict | None:
    folder = configured_oxdna_home()
    try:
        info = (Path(folder) / "runtime.json").stat()
    except FileNotFoundError:
        return None
    return dict(_checked_runtime(folder, (info.st_mtime_ns, info.st_size, info.st_ino)))
