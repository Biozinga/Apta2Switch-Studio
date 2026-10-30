"""Candidate-owned, cancellable oxDNA FFS jobs and portable result loading."""
from __future__ import annotations

from contextvars import copy_context
from datetime import datetime, timezone
import json
from pathlib import Path
import threading
import time
import uuid

from aptaswitch_core.localization import tr
from aptaswitch_core.oxdna_kinetics import (
    OPENING_PROTOCOL_VERSION, PROXIMITY_PROTOCOL_VERSION, validate_kinetic_protocol,
)
from .saved_switch_library import (
    _library_lock, load_switch, assert_switch_current as _current,
    saved_switch_identity as _identity,
)
from .saved_switch_runtime import _atomic_json, _inside
from .web_runtime import JobSnapshot


# Streamlit clears local modules after a source edit. The standard-library
# anchor survives that reload, just as the calculation threads do.
_JOBS, _LOCK = vars(threading).setdefault("_aptaswitch_ffs_registry", ({}, threading.RLock()))
_TERMINAL = {"completed", "cancelled", "failed"}
_ARTIFACTS = ("movie_file", "trajectory_file", "topology_file", "prepared_movie_file", "prepared_topology_file")


def prepare_ffs_request(*args, **kwargs):
    from aptaswitch_core.oxdna_ffs import prepare_ffs_request as prepare

    return prepare(*args, **kwargs)


def run_oxdna_ffs(*args, **kwargs):
    from aptaswitch_core.oxdna_ffs import run_oxdna_ffs as run

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
    from aptaswitch_core.oxdna_ffs import ffs_worker_running
    return ffs_worker_running(folder)


def _run_records(saved):
    root = _inside(saved, saved.path / "calculations" / "oxdna_ffs")
    for path in sorted(root.glob("*/request.json"), reverse=True):
        try:
            _inside(saved, path)
            metadata = _json(path)
            request = metadata.get("request")
            is_kinetics = (metadata.get("method") == "NATIVE_FFS"
                           or (isinstance(request, dict) and (request.get("kinetics") or {}).get("version")
                               in (OPENING_PROTOCOL_VERSION, PROXIMITY_PROTOCOL_VERSION)))
            if (is_kinetics and metadata.get("snapshot_id") == saved.snapshot["snapshot_id"]
                    and metadata.get("candidate_identity") == _identity(saved)):
                yield path.parent, metadata
        except (OSError, ValueError, KeyError):
            continue


def active_ffs(saved_or_path):
    path = Path(getattr(saved_or_path, "path", saved_or_path)).resolve()
    with _LOCK:
        job = _JOBS.get(str(path))
        if job and job.snapshot().state == "running":
            return job
        try:
            saved = saved_or_path if hasattr(saved_or_path, "snapshot") else load_switch(path)
            for folder, metadata in _run_records(saved):
                if _native_running(folder):
                    job = DetachedFFSJob(saved, folder, metadata)
                    _JOBS[str(path)] = job
                    return job
        except (OSError, ValueError):
            pass
        return None


def active_ffs_for_design(design_path):
    root = Path(design_path).expanduser().resolve()
    with _LOCK:
        jobs = tuple(_JOBS.values())
    running = next((job for job in jobs if job.saved.path.resolve().is_relative_to(root)
                    and job.snapshot().state == "running"), None)
    if running:
        return running
    from .saved_switch_library import LIBRARY_FOLDER
    for candidate in (root / LIBRARY_FOLDER).glob("*/snapshot.json"):
        if job := active_ffs(candidate.parent):
            return job
    return None


def start_ffs(saved, runtime: dict, options: dict):
    """Register under the same lock that protects candidate replacement."""
    from .saved_switch_runtime import active_analysis
    from .saved_switch_free_md_runtime import active_free_md

    # Context analysis startup holds its registry lock while publishing its
    # context sidecar. Query it before taking the library lock to avoid the
    # reverse lock order during simultaneous requests from two browser tabs.
    if active_analysis(saved) or active_free_md(saved):
        raise ValueError(tr("Attendez la fin du calcul de ce candidat avant de lancer la simulation."))
    with _library_lock(saved.path.parent), _LOCK:
        _current(saved)
        if active_free_md(saved):
            raise ValueError(tr("Attendez la fin du calcul de ce candidat avant de lancer la simulation."))
        existing = active_ffs(saved)
        if existing:
            return existing, False
        # Thermodynamic guide calculation belongs to the worker, outside the UI and
        # candidate-library lock. Keep an immutable copy of all chosen inputs.
        request = _copy({"runtime": runtime, "options": options})
        job = SavedSwitchFFSJob(saved, request)
        _JOBS[str(saved.path.resolve())] = job
        try:
            job.start()
        except Exception:
            _JOBS.pop(str(saved.path.resolve()), None)
            raise
        return job, True


def resume_ffs(saved, run_dir):
    """Resume this run in place, with the saved experiment and native banks."""
    from aptaswitch_core.oxdna_ffs import load_ffs_resume_request
    from .saved_switch_runtime import active_analysis
    from .saved_switch_free_md_runtime import active_free_md
    if active_analysis(saved) or active_free_md(saved):
        raise ValueError(tr("Attendez la fin du calcul de ce candidat avant de lancer la simulation."))
    folder = _inside(saved, Path(run_dir))
    root = _inside(saved, saved.path / "calculations" / "oxdna_ffs")
    if folder.parent.resolve() != root.resolve():
        raise ValueError("Invalid FFS calculation folder")
    with _library_lock(saved.path.parent), _LOCK:
        _current(saved)
        if active_free_md(saved):
            raise ValueError(tr("Attendez la fin du calcul de ce candidat avant de lancer la simulation."))
        if existing := active_ffs(saved):
            return existing, False
        result = next((item for item in saved_ffs_results(saved) if item["run_dir"] == str(folder)), None)
        if not result or not result.get("resumable"):
            raise ValueError(tr("Cette simulation ne dispose pas d’une sauvegarde permettant la reprise."))
        request = load_ffs_resume_request(folder)
        job = SavedSwitchFFSJob(saved, request, resume_dir=folder)
        _JOBS[str(saved.path.resolve())] = job
        try:
            job.start()
        except Exception:
            _JOBS.pop(str(saved.path.resolve()), None)
            raise
        return job, True


def _relative_artifact(value, folder: Path, old_folder) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        path = Path(value)
        if path.is_absolute():
            path = path.relative_to(Path(old_folder))
        if ".." in path.parts:
            return None
        target = folder / path
        if any(item.is_symlink() for item in (target, *target.parents) if item != folder.parent):
            return None
        if target.is_file() and target.resolve().is_relative_to(folder.resolve()):
            return path.as_posix()
    except (OSError, TypeError, ValueError):
        pass
    return None


def _saved_result(saved, folder, metadata):
    path = folder / "result.json"
    data = _json(path)
    if data.get("format") and (data.get("format") != "aptaswitch-oxdna-ffs-result"
            or data.get("version") != 2 or data.get("snapshot_id") != saved.snapshot["snapshot_id"]
            or data.get("candidate_identity") != _identity(saved)):
        return None
    frozen = _json(folder / "worker_request.json")
    requested = (frozen or metadata.get("request") or data.get("requested_parameters")
                 or {"kinetics": data.get("kinetics") or {"mode": data.get("simulation_mode", "opening")}})
    protocol = requested.get("kinetics") or {}
    mode = protocol.get("mode") or (requested.get("options") or {}).get("kinetics_mode", "opening")
    if mode not in {"opening", "proximity"}:
        # A first-contact-only result is not a rate for the complete pathway.
        return None
    if not data or data.get("status") not in _TERMINAL:
        data = {**_json(folder / "progress.json"), "status": "interrupted",
                "sampling_complete": False}
    data.setdefault("started_at", metadata.get("started_at", ""))
    data.setdefault("completed_at", metadata.get("started_at", ""))
    data.setdefault("run_id", folder.name)
    data.setdefault("method", "NATIVE_FFS")
    data["simulation_mode"] = mode
    data["requested_parameters"] = requested
    if mode == "proximity":
        data["rate_estimate"] = None
        data["kinetics"] = {**(data.get("kinetics") or {}), "mode": mode}
        data["kinetics"].pop("rate", None)
    data.setdefault("preparation", frozen.get("prepared"))
    for field in ("plan", "sampling", "box_size", "box_policy", "concentration_each_m", "visualization"):
        data.setdefault(field, frozen.get(field))
    data.setdefault("topology_file", "preparation/system.top")
    if not data.get("movie_file"):
        for name in ("trajectory.dat", "partial_trajectory.dat", "live.dat"):
            if (folder / name).is_file():
                data["movie_file"] = name
                break
    for field in _ARTIFACTS:
        if field in data:
            data[field] = _relative_artifact(data[field], folder, data.get("output_dir"))
    try:
        validate_kinetic_protocol(frozen.get("kinetics") or {})
        valid_protocol = True
    except (ValueError, TypeError):
        valid_protocol = False
    data["resumable"] = bool(isinstance(frozen.get("kinetics"), dict) and valid_protocol
                             and data.get("status") != "completed"
                             and (folder / "resume_state.json").is_file())
    data.update(output_dir=str(folder), run_dir=str(folder),
                result_path=str(path) if path.is_file() else None)
    return data


def saved_ffs_results(saved) -> list[dict]:
    """Include interrupted runs even if a crash prevented their final report."""
    with _LOCK:
        job = _JOBS.get(str(saved.path.resolve()))
        preparing_or_running = job.folder if job and job.snapshot().state == "running" else None
    results = []
    for folder, metadata in _run_records(saved):
        if folder == preparing_or_running or _native_running(folder):
            continue
        data = _saved_result(saved, folder, metadata)
        if data is not None:
            results.append(data)
    return sorted(results, key=lambda item: str(item.get("completed_at", "")), reverse=True)


class DetachedFFSJob:
    """Monitor a native worker recovered after a browser/server/module reload."""
    def __init__(self, saved, folder, metadata):
        self.saved, self.folder, self.run_id = saved, folder, folder.name
        self.metadata = metadata
        self._consumed = False
        self._finished = None

    def cancel(self):
        _current(self.saved)
        _inside(self.saved, self.folder / "CANCEL").touch()

    def snapshot(self):
        progress = _json(self.folder / "progress.json")
        running = _native_running(self.folder)
        result = None if running else _saved_result(self.saved, self.folder, self.metadata)
        if not running:
            self._finished = result
        message = progress.get("message", tr("Simulation oxDNA en cours…"))
        kinetics = progress.get("kinetics") or {}
        total = kinetics.get("replicas", progress.get("total", 1))
        current = max(0, kinetics["replica"] - 1) if running and "replica" in kinetics else progress.get("current", 0)
        return JobSnapshot(state="running" if running else "done",
                           current=current, total=total,
                           message=message, result=result, error=(result or {}).get("error"),
                           elapsed_seconds=float(progress.get("elapsed_seconds", 0)), log=(message,))

    def consume_once(self):
        if self._consumed or self._finished is None:
            return False
        self._consumed = True
        return True


class SavedSwitchFFSJob:
    def __init__(self, saved, request: dict, *, resume_dir=None):
        self.resuming = resume_dir is not None
        self.saved, self.request = saved, _copy(request)
        self.run_id = (Path(resume_dir).name if self.resuming else
                       datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_") + uuid.uuid4().hex[:10])
        self.folder = _inside(saved, saved.path / "calculations" / "oxdna_ffs" / self.run_id)
        self.folder.mkdir(parents=True, exist_ok=self.resuming)
        self._lock, self._cancel = threading.Lock(), threading.Event()
        self._state, self._result, self._error = "running", None, None
        self._current, self._total = 0, 1
        self._message = tr("Préparation de la simulation oxDNA…")
        self._log = []
        self._started, self._finished, self._consumed = time.monotonic(), None, False
        original = _json(self.folder / "request.json") if self.resuming else {}
        self.started_at = original.get("started_at") or datetime.now(timezone.utc).isoformat()
        self.resumed_at = datetime.now(timezone.utc).isoformat() if self.resuming else None
        self._prior_elapsed = float(_json(self.folder / "progress.json").get("elapsed_seconds", 0)) if self.resuming else 0.
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)

    def start(self):
        self._thread.start()

    def cancel(self):
        self._cancel.set()
        with self._lock:
            self._message = tr("Arrêt de la simulation demandé…")

    def _progress(self, current, total, message):
        _current(self.saved)
        with self._lock:
            self._current, self._total, self._message = current, total, str(message)
            self._log = [*self._log[-199:], str(message)]
        _atomic_json(self.saved, self.folder / "job_progress.json", {
            "snapshot_id": self.saved.snapshot["snapshot_id"], "status": "running",
            "current": current, "total": total, "message": str(message),
        })

    def _run(self):
        try:
            _current(self.saved)
            _atomic_json(self.saved, self.folder / "request.json", {
                "format": "aptaswitch-oxdna-ffs-request", "version": 2, "method": "NATIVE_FFS",
                "snapshot_id": self.saved.snapshot["snapshot_id"], "candidate_identity": _identity(self.saved),
                "run_id": self.run_id, "started_at": self.started_at, "request": self.request,
            })
            if self._cancel.is_set():
                result = {"status": "cancelled"}
            else:
                if self.resuming:
                    self._progress(0, 1, tr("Reprise de la simulation depuis les configurations sauvegardées…"))
                else:
                    self._progress(0, 1, tr("Calcul des guides 2D du switch OFF et du trigger isolés…"))
                    self.request = prepare_ffs_request(self.saved.snapshot, self.request["runtime"], self.request["options"])
                _current(self.saved)
                _atomic_json(self.saved, self.folder / "request.json", {
                    "format": "aptaswitch-oxdna-ffs-request", "version": 2, "method": "NATIVE_FFS",
                    "snapshot_id": self.saved.snapshot["snapshot_id"], "candidate_identity": _identity(self.saved),
                    "run_id": self.run_id, "started_at": self.started_at, "request": self.request,
                })
                if self._cancel.is_set():
                    result = {"status": "cancelled"}
                else:
                    result = run_oxdna_ffs(self.request, output_dir=self.folder,
                                          progress_callback=self._progress, cancel_check=self._cancel.is_set)
            if not isinstance(result, dict) or result.get("status") not in _TERMINAL:
                raise ValueError(tr("Le moteur oxDNA n’a pas renvoyé de résultat valide."))
        except Exception as exc:
            result = {"status": "cancelled" if self._cancel.is_set() else "failed", "error": str(exc)}
        mode = ((self.request.get("kinetics") or {}).get("mode")
                or (self.request.get("options") or {}).get("kinetics_mode", "opening"))
        result.update(format="aptaswitch-oxdna-ffs-result", version=2, method="NATIVE_FFS",
                      simulation_mode=mode,
                      snapshot_id=self.saved.snapshot["snapshot_id"], candidate_identity=_identity(self.saved),
                      run_id=self.run_id, output_dir=str(self.folder), requested_parameters=self.request,
                      started_at=self.started_at, resumed_at=self.resumed_at,
                      completed_at=datetime.now(timezone.utc).isoformat())
        if mode == "proximity":
            result["rate_estimate"] = None
            result["kinetics"] = {**(result.get("kinetics") or {}), "mode": mode}
            result["kinetics"].pop("rate", None)
        try:
            _atomic_json(self.saved, self.folder / "result.json", result)
        except Exception as exc:
            result.update(status="failed", error=str(exc))
        else:
            result["result_path"] = str(self.folder / "result.json")
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
                               elapsed_seconds=self._prior_elapsed + (self._finished or time.monotonic()) - self._started,
                               log=tuple(self._log))

    def consume_once(self):
        with self._lock:
            if self._state == "running" or self._consumed:
                return False
            self._consumed = True
            return True
