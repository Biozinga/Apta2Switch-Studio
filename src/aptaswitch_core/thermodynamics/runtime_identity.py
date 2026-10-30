"""Record the actual ViennaRNA binary, including locally audited corrections.

The native __version__ stays 2.7.2 after a source patch. A version string alone
therefore cannot distinguish results from the released and corrected builds.
Binary identity is always recorded; a supplied build manifest additionally
binds a known source/patch recipe to that exact loaded binary. A manifest is
provenance evidence, not a substitute for numerical model validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path


_SALT_LOOP_FIX = {
    "upstream_version": "2.7.2",
    "patch_id": "aptaswitch-salt-loop-31-v1",
    "upstream_archive_sha256": "1ab5f4a4f76fc85a2243546088e45f5d85f2d7a56cc656e969b005cce9bfab5f",
    "original_params_c_sha256": "a6da5b5038e07c6b4b0b843429115db84b6bcb8b48e51074f3180d2842358e7a",
    "patch_sha256": "700a7e70eb0aa2bc1b1e7778f827bd46ecaabba9ce4ab299b1919485ea424a7e",
    "expected_patched_params_c_sha256": "0163c59586c5389acc0b1284a97fda7fd0652716b1d26548466e5df735956cf9",
    "status": "built",
}


def _stamp(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino


@lru_cache(maxsize=16)
def _digest(path: Path, stamp: tuple[int, int, int, int]) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    if _stamp(path) != stamp:
        raise RuntimeError("ViennaRNA's native binary changed while its identity was being read")
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class NativeRuntimeIdentity:
    upstream_version: str
    module_path: str
    module_sha256: str
    patch_id: str | None = None
    manifest_sha256: str | None = None

    @property
    def version_label(self) -> str:
        return (f"{self.upstream_version};native-sha256={self.module_sha256};"
                f"patch={self.patch_id or 'unrecorded'}")


def identify_vienna_runtime(RNA, manifest_path: str | Path | None = None) -> NativeRuntimeIdentity:
    """Identify the already loaded module; never load/reload a different build.

    An absent manifest is explicitly recorded as unrecorded, never guessed to
    be the unmodified release. A supplied manifest must describe the audited
    salt-loop correction and match the binary actually imported by the caller.
    This also catches accidentally cached stock bindings in a patched worker.
    """
    if getattr(RNA, "__version__", None) != "2.7.2":
        raise RuntimeError("Only the audited ViennaRNA 2.7.2 runtime is supported")
    module_file = getattr(getattr(RNA, "_RNA", None), "__file__", None)
    if not isinstance(module_file, (str, Path)):
        raise RuntimeError("ViennaRNA's loaded native binary cannot be identified")
    try:
        path = Path(module_file).resolve(strict=True)
        if not path.is_file():
            raise ValueError("Native module path is not a file")
        module_digest = _digest(path, _stamp(path))
    except (OSError, ValueError) as error:
        raise RuntimeError("ViennaRNA's loaded native binary cannot be read") from error
    patch_id = manifest_digest = None
    if manifest_path is not None:
        try:
            raw = Path(manifest_path).read_bytes()
            manifest = json.loads(raw)
        except (OSError, TypeError, ValueError) as error:
            raise ValueError("Cannot read the ViennaRNA build manifest") from error
        if not isinstance(manifest, dict) or any(
                manifest.get(key) != value for key, value in _SALT_LOOP_FIX.items()):
            raise ValueError("Build manifest does not describe the audited ViennaRNA salt-loop fix")
        if manifest.get("native_module_sha256") != module_digest:
            raise ValueError("Build manifest does not match the loaded ViennaRNA binary")
        patch_id = manifest["patch_id"]
        manifest_digest = hashlib.sha256(raw).hexdigest()
    return NativeRuntimeIdentity(RNA.__version__, str(path), module_digest,
                                 patch_id, manifest_digest)
