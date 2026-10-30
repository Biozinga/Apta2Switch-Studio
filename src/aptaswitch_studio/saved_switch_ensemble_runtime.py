"""Candidate-owned structure families for a frozen full-reporter analysis.

Reuse the design ensemble worker, but publish under the original context run
and guard candidate identity exactly as for other saved-candidate artifacts.
"""
from __future__ import annotations

from aptaswitch_core.localization import tr

import json
from pathlib import Path
import re
import threading

from aptaswitch_core.saved_switch_ensemble import context_ensemble_source
from aptaswitch_core.structure_ensemble import ensemble_inputs, ensemble_key, ensemble_path
from aptaswitch_core.thermodynamics.setup import ENGINE_ID, get_backend
from .ensemble_runtime import EnsembleJob
from .saved_switch_library import _library_lock, assert_switch_current
from .saved_switch_runtime import _atomic_json, _inside


# Like the calculation thread, this registry survives Streamlit module reloads.
_JOBS, _LOCK = vars(threading).setdefault("_aptaswitch_context_ensemble_registry", ({}, threading.RLock()))


def _source(saved, analysis, *, backend=None):
    assert_switch_current(saved)
    run_id = analysis.get("run_id") if isinstance(analysis, dict) else None
    if (not isinstance(run_id, str) or re.fullmatch(r"[A-Za-z0-9_-]{1,100}", run_id) is None
            or analysis.get("snapshot_id") != saved.snapshot["snapshot_id"]):
        raise ValueError("Invalid context analysis for this saved candidate")
    folder = _inside(saved, saved.analyses_dir / run_id)
    path = _inside(saved, folder / "result.json")
    if path.stat().st_size > 32 * 1024 * 1024:
        raise ValueError("Saved context analysis is too large")
    data = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(data, dict) or data.get("snapshot_id") != saved.snapshot["snapshot_id"]
            or data.get("run_id") != run_id or data.get("status") not in {"completed", "cancelled", "failed"}):
        raise ValueError("Structure families require a saved context reference for this candidate")
    cache = _inside(saved, folder / "structure_ensembles")
    if cache.is_dir():
        for entry in cache.glob("*.json"):
            _inside(saved, entry)
    return context_ensemble_source(saved.snapshot, data, folder, backend=backend)


def load_context_ensemble_source(saved, analysis, *, backend=None):
    """Use the current candidate's local artifact, never a stored absolute path."""
    with _library_lock(saved.path.parent):
        return _source(saved, analysis, backend=backend)


def active_context_ensemble(saved_or_path):
    path = Path(getattr(saved_or_path, "path", saved_or_path)).resolve()
    with _LOCK:
        job = _JOBS.get(str(path))
        return job if job and job.snapshot().state == "running" else None


def active_context_ensemble_for_design(design_path):
    root = Path(design_path).expanduser().resolve()
    with _LOCK:
        jobs = tuple(_JOBS.values())
    return next((job for job in jobs if job.saved.path.resolve().is_relative_to(root)
                 and job.snapshot().state == "running"), None)


def _native_backend(result, import_path, backend):
    if getattr(result.request, "engine", "") != ENGINE_ID:
        raise ValueError("This calculation requires the native-hybrid engine.")
    result.request.thermo.validate_model()
    return backend if backend is not None else get_backend(import_path)


def create_context_ensemble_job(saved, analysis, result, candidate, import_path=None, *, num_sample=1000, backend=None):
    """Prepare a detached worker; the common renderer starts it explicitly."""
    backend = _native_backend(result, import_path, backend)
    canonical_result, canonical_candidate = load_context_ensemble_source(saved, analysis, backend=backend)
    inputs = ensemble_inputs(canonical_result.request, canonical_candidate, num_sample, backend=backend)
    if (ensemble_inputs(result.request, candidate, num_sample, backend=backend) != inputs
            or Path(result.output_dir).absolute() != canonical_result.output_dir.absolute()):
        raise ValueError("Structure-family inputs do not match the saved context analysis")
    return SavedSwitchEnsembleJob(saved, analysis, canonical_result, canonical_candidate,
                                  import_path, num_sample=num_sample, backend=backend)


class SavedSwitchEnsembleJob(EnsembleJob):
    def __init__(self, saved, analysis, result, candidate, import_path=None, *, num_sample=1000, backend=None):
        backend = _native_backend(result, import_path, backend)
        self.saved = saved
        self.analysis_identity = {"snapshot_id": analysis["snapshot_id"], "run_id": analysis["run_id"]}
        self.inputs = ensemble_inputs(result.request, candidate, num_sample, backend=backend)
        self.key = f"context_ensemble_{saved.snapshot['snapshot_id']}_{analysis['run_id']}_" + ensemble_key(self.inputs)
        super().__init__(result, candidate, import_path, num_sample=num_sample, save_result=self._publish, backend=backend)

    def start(self):
        # Lock order agrees with candidate replacement: library, then registry.
        with _library_lock(self.saved.path.parent), _LOCK:
            assert_switch_current(self.saved)
            from .saved_switch_free_md_runtime import active_free_md
            if active_free_md(self.saved):
                raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
            current_result, current_candidate = _source(self.saved, self.analysis_identity, backend=self._backend)
            if ensemble_inputs(current_result.request, current_candidate, self._num_sample, backend=self._backend) != self.inputs:
                raise ValueError("The saved context analysis changed before sampling started")
            active = active_context_ensemble(self.saved)
            if active is self:
                return
            if active is not None:
                raise ValueError("A structure-family analysis is already running for this candidate")
            key = str(self.saved.path.resolve())
            _JOBS[key] = self
        # Registration already blocks replacement; do not hold a filesystem
        # lock while calling worker code (also permits synchronous test runners).
        try:
            super().start()
        except Exception:
            with _LOCK:
                if _JOBS.get(key) is self:
                    _JOBS.pop(key, None)
            raise

    def _progress(self, current, total, message):
        assert_switch_current(self.saved)
        super()._progress(current, total, message)

    def _publish(self, result):
        if result.get("inputs") != self.inputs or result.get("key") != ensemble_key(self.inputs):
            raise ValueError("Structure-family output does not match the requested context")
        path = ensemble_path(self._result_source.output_dir, result["key"])
        _atomic_json(self.saved, path, {**result, "context_source": self.analysis_identity})

    def _run(self):
        try:
            super()._run()
        finally:
            with _LOCK:
                key = str(self.saved.path.resolve())
                if _JOBS.get(key) is self:
                    _JOBS.pop(key, None)
