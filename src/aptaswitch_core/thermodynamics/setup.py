"""Locate and validate the application's explicit, identified hybrid runtime.

The library and its build manifest travel together. No fallback imports an
unrelated folding engine; missing or invalid native runtimes fail explicitly.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from functools import lru_cache
import os
from pathlib import Path
import sys

from .native_backend import NativeHybridBackend
from .types import ComplexSpec, Conditions, Strand

ENGINE_ID = "native-hybrid"
MODEL_ID = "experimental-Sugimoto95-Vienna04-loopmean"
LIBRARY_BASENAME = ("aptaswitch_hybrid.dll" if sys.platform == "win32" else
                    "libaptaswitch_hybrid.dylib" if sys.platform == "darwin" else "libaptaswitch_hybrid.so")


@dataclass(frozen=True)
class RuntimeValidation:
    ok: bool
    path: str
    message: str
    provenance: dict = field(default_factory=dict)


def _library_path(path) -> Path:
    target = Path(path).expanduser()
    if target.is_dir():
        target = target / LIBRARY_BASENAME
    return target.resolve()


def configured_library_path() -> str | None:
    """Discover the bundled runtime; an environment override is for maintainers.

    An explicitly configured but missing path is returned so validation reports
    the error; it must never silently select another library instead.
    """
    explicit = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")
    if explicit:
        return str(_library_path(explicit))
    from ..run_storage import application_directory
    roots = [Path(__file__).resolve().parents[3], application_directory()]
    if getattr(sys, "frozen", False):
        bundle = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
        roots[:0] = [bundle, bundle.parent / "Resources"]
    for root in dict.fromkeys(roots):
        for suffix in ("thermodynamics", ".local_deps/thermodynamics", "aptaswitch_core/thermodynamics/lib"):
            candidate = root / suffix / LIBRARY_BASENAME
            if candidate.is_file():
                return str(candidate.resolve())
    return None


def _stamp(path: Path):
    info = path.stat()
    return info.st_mtime_ns, info.st_size, info.st_ino


@lru_cache(maxsize=8)
def _backend(path: str, binary_stamp, manifest_stamp):
    # NativeHybridBackend validates the manifest and binary when loading it.
    return NativeHybridBackend(path)


def get_backend(library_path=None) -> NativeHybridBackend:
    target = library_path or configured_library_path()
    if not target:
        raise RuntimeError("Le moteur thermodynamique intégré est absent. Réinstallez l’application complète.")
    path = _library_path(target)
    manifest = path.with_suffix(path.suffix + ".json")
    try:
        return _backend(str(path), _stamp(path), _stamp(manifest))
    except OSError as exc:
        raise RuntimeError(f"Le moteur hybride et son manifeste doivent être présents : {path}") from exc


def backend_for_request(request, *, library_path=None) -> NativeHybridBackend:
    """A newly requested calculation uses the configured native model.

    Historical request provenance remains attached to its stored results; new
    analyses record their own native identity rather than rewriting that data.
    """
    if getattr(request, "engine", None) != ENGINE_ID:
        raise ValueError("This calculation requires the native-hybrid engine.")
    backend = get_backend(library_path if library_path is not None else
                          getattr(request, "thermodynamics_library_path", None))
    if request.thermo.material != backend.provenance.model:
        raise ValueError("Le modèle demandé diffère du modèle du moteur hybride sélectionné.")
    return backend


def validate_runtime(path) -> RuntimeValidation:
    resolved = str(_library_path(path))
    try:
        backend = get_backend(resolved)
        spec = ComplexSpec((Strand("GCAAA", "rna", "switch"), Strand("TTTGC", "dna", "trigger")))
        conditions = Conditions(37, .12, 0)
        ensemble = backend.ensemble(spec, conditions, 2, seed=17)
        ensemble.analysis.partition.require_nonempty()
        minimum = backend.minimum(spec, conditions)
        probability = backend.structure_probability(spec, conditions, minimum.structure,
                                                   partition=ensemble.analysis.partition)
        if not 0 < probability <= 1 or len(ensemble.samples.structures) != 2:
            raise ValueError("Le test de probabilités du moteur hybride a échoué.")
        off = ComplexSpec((Strand("GCGAAAGC", "rna", "Switch"),))
        off_ensemble = backend.ensemble(off, conditions, 4, seed=17)
        off_minimum = backend.minimum(off, conditions)
        if (off_minimum.context.solver != "ViennaRNA-RNAlib-2.7.2"
                or off_ensemble.samples.context.solver != off_minimum.context.solver
                or len(off_ensemble.samples.structures) != 4
                or not 0 < backend.structure_probability(off, conditions, off_minimum.structure,
                    partition=off_ensemble.analysis.partition) <= 1):
            raise ValueError("Le test OFF ViennaRNA a échoué.")
        return RuntimeValidation(True, resolved, "Moteur hybride validé : ARN/ADN, probabilités et échantillonnage.",
                                 asdict(backend.provenance))
    except Exception as exc:
        return RuntimeValidation(False, resolved, str(exc))
