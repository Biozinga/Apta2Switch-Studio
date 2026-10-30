"""Recover a retained candidate's original design and validated ensemble caches.

This module never computes structures. The immutable snapshot supplies the MFE
and scientific settings; existing family artifacts are detached from the source
run so they remain usable when that run is moved or deleted.
"""

from __future__ import annotations

import json
import os
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

from aptaswitch_core.models import DesignCandidate, RunResult
from aptaswitch_core.run_storage import default_export_dir, legacy_export_dir, run_output_directory
from aptaswitch_core.structure_ensemble import ensemble_inputs, ensemble_key, ensemble_path, ensemble_sample_variant_keys, load_ensemble, recorded_native_ensemble_inputs, save_ensemble

from .run_library import _load_design
from .saved_switch_library import LIBRARY_FOLDER, SavedSwitch, _checked_path, _library_lock, assert_switch_current


def _saved_result(snapshot: dict, output_dir: Path) -> RunResult:
    return _load_design(deepcopy({
        **snapshot["source"], "request": snapshot["request"],
        "results": [snapshot["candidate"]],
    }), output_dir / "snapshot.json", retained_candidate=True)


def _comparison_source(analysis: dict) -> str | None:
    """Only an explicit comparison without redesign links to an earlier run."""
    if not isinstance(analysis, dict) or analysis.get("switch_design_repeated") is not False:
        return None
    run_id = analysis.get("source_run_id")
    if (not isinstance(run_id, str) or not run_id or run_id in {".", ".."}
            or any(character in run_id for character in "/\\:")
            or Path(run_id).name != run_id):
        return None
    return run_id


@lru_cache(maxsize=32)
def _recorded_comparison_source(path: str, modified_ns: int, size: int) -> str | None:
    # Retain only lineage metadata; probability matrices and complete results
    # must not accumulate in the UI cache. Stat arguments invalidate on change.
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("run_id") != Path(path).parent.name:
        return None
    return _comparison_source(payload.get("inverse_extension_analysis", {}))


def design_ensemble_sources(result, *, source_output_dir: Path | None = None) -> tuple[Path, ...]:
    """Find explicit comparison ancestors beside a trusted run directory.

    Recorded absolute paths are never followed. Each discovered ensemble still
    requires exact candidate/conditions and complete scientific validation.
    """
    # The displayed run directory is already selected by the application.
    # Canonicalize that boundary (including macOS /var -> /private/var), then
    # reject links in every discovered descendant. Explicit copy paths retain
    # the saved-library's stricter checks before any files are published.
    source = (_checked_path(source_output_dir) if source_output_dir is not None
              else Path(result.output_dir).expanduser().resolve())
    parent_id = _comparison_source(getattr(result, "inverse_extension_analysis", {}))
    sources, visited = [], set()
    for _ in range(64):
        if source in visited:
            break
        visited.add(source)
        try:
            _checked_path(source / "structure_ensembles")
            sources.append(source)
        except ValueError:
            pass
        if parent_id is None:
            break
        try:
            source = _checked_path(source.parent / parent_id)
            # The explicitly named parent can contain valid caches even when
            # its result JSON is absent. Only a valid record extends the chain.
            record = _checked_path(source / f"{parent_id}_results.json")
            details = record.stat()
            parent_id = (_recorded_comparison_source(str(record), details.st_mtime_ns, details.st_size)
                         if details.st_size <= 128 * 1024 * 1024 else None)
        except (OSError, ValueError, TypeError, RecursionError):
            parent_id = None
    return tuple(sources)


def copy_saved_design_ensembles(snapshot: dict, source_output_dir: Path, output_dir: Path) -> int:
    """Copy only complete caches matching this exact candidate and conditions.

    The source path is provided by a trusted caller, never taken from saved
    request paths. Filenames discover possible sample variants; the shared
    load_ensemble validator verifies their key, inputs, matrices and structures.
    Each artifact is published atomically, including when recovering old saves.
    """
    source = _checked_path(source_output_dir)
    destination = _checked_path(output_dir)
    if source == destination:
        return 0
    _checked_path(source / "structure_ensembles")
    result = _saved_result(snapshot, destination)
    return sum(_copy_design_ensemble_cache(result, directory, destination)
               for directory in design_ensemble_sources(result, source_output_dir=source))


def _copy_design_ensemble_cache(result, source: Path, destination: Path) -> int:
    cache = _checked_path(source / "structure_ensembles")
    if not cache.is_dir():
        return 0
    candidate = result.candidates[0]
    inputs = ensemble_inputs(result.request, candidate)
    variants = ensemble_sample_variant_keys(json.dumps(inputs, sort_keys=True, allow_nan=False))
    copied = 0
    for path in sorted(cache.glob("*.json")):
        if (path.is_symlink() or not path.is_file() or len(path.stem) != 64
                or set(path.stem) - set("0123456789abcdef")):
            continue
        target = _checked_path(ensemble_path(destination, path.stem))
        if target.is_file():
            # A detached artifact is immutable here. The shared explorer still
            # validates it on reading and reports corrupt local files directly.
            continue
        try:
            details = path.stat()
            selected_inputs = ({**inputs, "num_sample": variants[path.stem]} if path.stem in variants
                               else recorded_native_ensemble_inputs(path, result.request, candidate))
            if selected_inputs is None:
                continue
            artifact = load_ensemble(source, selected_inputs)
            if artifact is None:
                continue
        except (OSError, ValueError, TypeError, KeyError, RecursionError):
            # A missing, unrelated or incomplete cache cannot change the saved
            # candidate or prevent recovery of other valid sample variants.
            continue
        save_ensemble(destination, artifact)
        os.utime(target, ns=(details.st_atime_ns, details.st_mtime_ns))
        copied += 1
    return copied


def _source_directories(saved: SavedSwitch, run_id: str):
    """Resolve known run layouts without trusting absolute JSON paths."""
    roots = []
    if saved.path.parent.name == LIBRARY_FOLDER:
        base = saved.path.parent.parent
        roots.extend((base / "calculations", base))
    roots.extend((default_export_dir(), legacy_export_dir()))
    seen = set()
    for root in roots:
        try:
            source = _checked_path(run_output_directory(root, "design") / run_id)
        except ValueError:
            continue
        if source not in seen:
            seen.add(source)
            yield source


def load_saved_design(saved: SavedSwitch) -> tuple[RunResult, DesignCandidate]:
    """Restore the original MFE and any existing, independently stored families."""
    with _library_lock(saved.path.parent):
        assert_switch_current(saved)
        output_dir = _checked_path(saved.path / "calculations" / "design")
        cache = _checked_path(output_dir / "structure_ensembles")
        for path in cache.glob("*.json"):
            _checked_path(path)
        result = _saved_result(saved.snapshot, output_dir)
        for source in _source_directories(saved, result.run_id):
            try:
                _checked_path(source / "structure_ensembles")
            except ValueError:
                continue
            copy_saved_design_ensembles(saved.snapshot, source, output_dir)
        return result, result.candidates[0]
