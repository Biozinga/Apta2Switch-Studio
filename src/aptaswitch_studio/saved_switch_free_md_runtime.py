"""Candidate-owned, resumable DNA-DNA free-MD comparisons.

The immutable saved switch and frozen request own each experiment. A paused run
keeps its directory while another simulation can run; resuming never reads the
current form or draws a new starting configuration.
"""
from __future__ import annotations

from contextvars import copy_context
from datetime import datetime, timezone
import json
from pathlib import Path
import threading
import time
import uuid

from aptaswitch_core.localization import tr
from .saved_switch_library import (
    LIBRARY_FOLDER, _library_lock, load_switch,
    assert_switch_current as _current, saved_switch_identity as _identity,
)
from .saved_switch_runtime import _atomic_json, _inside
from .web_runtime import JobSnapshot

_JOBS, _LOCK = vars(threading).setdefault("_aptaswitch_free_md_registry", ({}, threading.RLock()))
_TERMINAL = {"completed", "cancelled", "failed"}
_METHOD = "FREE_MD_DNA_DNA"


def prepare_free_md_request(*args, **kwargs):
    from aptaswitch_core.oxdna_free_md import prepare_free_md_request as prepare
    return prepare(*args, **kwargs)


def run_free_md(*args, **kwargs):
    from aptaswitch_core.oxdna_free_md import run_free_md as run
    return run(*args, **kwargs)


def _copy(value):
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def _json(path):
    try:
        if path.is_symlink() or path.stat().st_size > 32 * 1024 * 1024:
            return {}
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _native_running(folder):
    from aptaswitch_core.oxdna_free_md import free_md_worker_running
    return free_md_worker_running(folder)


def _run_records(saved):
    root = _inside(saved, saved.path / "calculations" / "oxdna_free_md")
    for path in sorted(root.glob("*/request.json"), reverse=True):
        try:
            _inside(saved, path)
            metadata = _json(path)
            if (metadata.get("method") == _METHOD
                    and metadata.get("snapshot_id") == saved.snapshot["snapshot_id"]
                    and metadata.get("candidate_identity") == _identity(saved)):
                yield path.parent, metadata
        except (OSError, ValueError, KeyError):
            continue


def active_free_md(saved_or_path):
    path = Path(getattr(saved_or_path, "path", saved_or_path)).resolve()
    with _LOCK:
        job = _JOBS.get(str(path))
        if job and job.snapshot().state == "running":
            return job
        try:
            saved = saved_or_path if hasattr(saved_or_path, "snapshot") else load_switch(path)
            for folder, metadata in _run_records(saved):
                if _native_running(folder):
                    job = DetachedFreeMDJob(saved, folder, metadata)
                    _JOBS[str(path)] = job
                    return job
        except (OSError, ValueError):
            pass
        return None


def active_free_md_for_design(design_path):
    root = Path(design_path).expanduser().resolve()
    with _LOCK:
        jobs = tuple(_JOBS.values())
    running = next((job for job in jobs if job.saved.path.resolve().is_relative_to(root)
                    and job.snapshot().state == "running"), None)
    if running:
        return running
    for manifest in (root / LIBRARY_FOLDER).glob("*/snapshot.json"):
        if job := active_free_md(manifest.parent):
            return job
    return None


def _other_candidate_job(saved):
    from .saved_switch_ffs_runtime import active_ffs
    from .saved_switch_runtime import active_analysis
    from .saved_switch_ensemble_runtime import active_context_ensemble
    return active_analysis(saved) or active_ffs(saved) or active_context_ensemble(saved)


def _locked_candidate_job(saved):
    from .saved_switch_runtime import active_analysis
    # Registration is serialized by the candidate-library lock.
    from .saved_switch_ffs_runtime import active_ffs
    from .saved_switch_ensemble_runtime import active_context_ensemble
    return active_analysis(saved) or active_ffs(saved) or active_context_ensemble(saved)


def start_free_md(saved, runtime: dict, options: dict):
    if _other_candidate_job(saved):
        raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
    # Validation precedes publication. It does not compute native structures.
    request = prepare_free_md_request(saved.snapshot, runtime, options)
    with _library_lock(saved.path.parent), _LOCK:
        _current(saved)
        if existing := active_free_md(saved):
            return existing, False
        if _locked_candidate_job(saved):
            raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
        job = SavedSwitchFreeMDJob(saved, request)
        _JOBS[str(saved.path.resolve())] = job
    try:
        job.start()
    except Exception:
        with _LOCK:
            if _JOBS.get(str(saved.path.resolve())) is job:
                _JOBS.pop(str(saved.path.resolve()), None)
        raise
    return job, True


def resume_free_md(saved, run_dir):
    from aptaswitch_core.oxdna_free_md import load_free_md_resume_request
    if _other_candidate_job(saved):
        raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
    folder = _inside(saved, Path(run_dir))
    root = _inside(saved, saved.path / "calculations" / "oxdna_free_md")
    if folder.parent.resolve() != root.resolve():
        raise ValueError("Invalid free MD calculation folder")
    with _library_lock(saved.path.parent), _LOCK:
        _current(saved)
        if existing := active_free_md(saved):
            return existing, False
        if _locked_candidate_job(saved):
            raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
        result = next((row for row in saved_free_md_results(saved) if row["run_dir"] == str(folder)), None)
        if not result or not result.get("resumable"):
            raise ValueError(tr("Cette simulation ne dispose pas d’une sauvegarde permettant la reprise."))
        request = load_free_md_resume_request(folder)
        job = SavedSwitchFreeMDJob(saved, request, resume_dir=folder)
        _JOBS[str(saved.path.resolve())] = job
    try:
        job.start()
    except Exception:
        with _LOCK:
            if _JOBS.get(str(saved.path.resolve())) is job:
                _JOBS.pop(str(saved.path.resolve()), None)
        raise
    return job, True


def _saved_result(saved, folder, metadata):
    data = _json(folder / "result.json")
    if data.get("format") and (data.get("format") != "aptaswitch-oxdna-free-md-result"
            or data.get("version") != 1 or data.get("snapshot_id") != saved.snapshot["snapshot_id"]
            or data.get("candidate_identity") != _identity(saved)):
        return None
    if not data or data.get("status") not in _TERMINAL:
        data = {**_json(folder / "progress.json"), "status": "interrupted", "sampling_complete": False}
    request = _json(folder / "worker_request.json") or metadata.get("request") or {}
    data.update(requested_parameters=request, run_id=folder.name, run_dir=str(folder), output_dir=str(folder))
    data.setdefault("started_at", metadata.get("started_at", ""))
    data.setdefault("completed_at", metadata.get("started_at", ""))
    data["resumable"] = data.get("status") != "completed" and (folder / "resume_state.json").is_file()
    return data


def saved_free_md_results(saved):
    with _LOCK:
        job = _JOBS.get(str(saved.path.resolve()))
        active_folder = job.folder if job and job.snapshot().state == "running" else None
    results = []
    for folder, metadata in _run_records(saved):
        if folder == active_folder or _native_running(folder):
            continue
        data = _saved_result(saved, folder, metadata)
        if data is not None:
            results.append(data)
    return sorted(results, key=lambda row: str(row.get("completed_at", "")), reverse=True)


class DetachedFreeMDJob:
    def __init__(self, saved, folder, metadata):
        self.saved, self.folder, self.metadata = saved, folder, metadata
        self.run_id = folder.name
        self.request = _json(folder / "worker_request.json") or metadata.get("request") or {}
        self._finished, self._consumed = None, False

    def cancel(self):
        _current(self.saved)
        _inside(self.saved, self.folder / "CANCEL").touch()

    def snapshot(self):
        progress = _json(self.folder / "progress.json")
        running = _native_running(self.folder)
        result = None if running else _saved_result(self.saved, self.folder, self.metadata)
        if not running:
            self._finished = result
        message = progress.get("message", tr("MD libre ADN–ADN en cours…"))
        return JobSnapshot(state="running" if running else "done", current=progress.get("current", 0),
                           total=progress.get("total", 1), message=message, result=result,
                           error=(result or {}).get("error"), elapsed_seconds=float(progress.get("elapsed_seconds", 0)),
                           log=(message,))

    def consume_once(self):
        if self._consumed or self._finished is None:
            return False
        self._consumed = True
        return True


class SavedSwitchFreeMDJob:
    def __init__(self, saved, request, *, resume_dir=None):
        self.saved, self.request = saved, _copy(request)
        self.resuming = resume_dir is not None
        self.run_id = (Path(resume_dir).name if self.resuming else
                       datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_") + uuid.uuid4().hex[:10])
        self.folder = _inside(saved, saved.path / "calculations" / "oxdna_free_md" / self.run_id)
        self.folder.mkdir(parents=True, exist_ok=self.resuming)
        self._lock, self._cancel = threading.Lock(), threading.Event()
        self._state, self._result, self._error = "running", None, None
        self._current, self._total = 0, 1
        self._message = tr("Préparation de la MD libre ADN–ADN…")
        self._log = []
        self._started, self._finished, self._consumed = time.monotonic(), None, False
        old = _json(self.folder / "request.json") if self.resuming else {}
        self.started_at = old.get("started_at") or datetime.now(timezone.utc).isoformat()
        self.resumed_at = datetime.now(timezone.utc).isoformat() if self.resuming else None
        self._thread = threading.Thread(target=copy_context().run, args=(self._run,), daemon=True)

    def start(self):
        self._thread.start()

    def cancel(self):
        self._cancel.set()
        with self._lock:
            self._message = tr("Mise en pause demandée ; sauvegarde de la trajectoire…")

    def _progress(self, current, total, message):
        _current(self.saved)
        with self._lock:
            self._current, self._total, self._message = current, total, str(message)
            self._log = [*self._log[-199:], str(message)]

    def _run(self):
        try:
            _atomic_json(self.saved, self.folder / "request.json", {
                "format": "aptaswitch-oxdna-free-md-request", "version": 1, "method": _METHOD,
                "snapshot_id": self.saved.snapshot["snapshot_id"], "candidate_identity": _identity(self.saved),
                "run_id": self.run_id, "started_at": self.started_at, "request": self.request,
            })
            _current(self.saved)
            result = run_free_md(self.request, output_dir=self.folder,
                                 progress_callback=self._progress, cancel_check=self._cancel.is_set)
            if not isinstance(result, dict) or result.get("status") not in _TERMINAL:
                raise ValueError(tr("Le moteur oxDNA n’a pas renvoyé de résultat valide."))
        except Exception as exc:
            result = {"status": "cancelled" if self._cancel.is_set() else "failed", "error": str(exc)}
        result.update(format="aptaswitch-oxdna-free-md-result", version=1, method=_METHOD,
                      snapshot_id=self.saved.snapshot["snapshot_id"], candidate_identity=_identity(self.saved),
                      run_id=self.run_id, output_dir=str(self.folder), requested_parameters=self.request,
                      started_at=self.started_at, resumed_at=self.resumed_at,
                      completed_at=datetime.now(timezone.utc).isoformat())
        try:
            _atomic_json(self.saved, self.folder / "result.json", result)
        except Exception as exc:
            result.update(status="failed", error=str(exc))
        with self._lock:
            self._result, self._error = result, result.get("error")
            self._state = "failed" if result["status"] == "failed" else "done"
            self._finished = time.monotonic()
        with _LOCK:
            if _JOBS.get(str(self.saved.path.resolve())) is self:
                _JOBS.pop(str(self.saved.path.resolve()), None)

    def snapshot(self):
        with self._lock:
            return JobSnapshot(state=self._state, current=self._current, total=self._total,
                               message=self._message, result=self._result, error=self._error,
                               elapsed_seconds=(self._finished or time.monotonic()) - self._started,
                               log=tuple(self._log))

    def consume_once(self):
        with self._lock:
            if self._state == "running" or self._consumed:
                return False
            self._consumed = True
            return True
