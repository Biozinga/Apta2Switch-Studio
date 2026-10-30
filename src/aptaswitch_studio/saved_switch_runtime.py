"""Persistent context analyses for named, retained switch candidates."""
from __future__ import annotations

from aptaswitch_core.localization import tr

from contextvars import copy_context
from datetime import datetime, timezone
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import tempfile
import threading
import time
import uuid

from aptaswitch_core.saved_switch_analysis import analyze_saved_switch, build_saved_switch_context, validate_analysis_config
from .saved_switch_library import assert_switch_current as _current, _library_lock
from aptaswitch_core.thermodynamics.setup import get_backend
from .web_runtime import JobSnapshot

_JOBS = {}
_LOCK = threading.RLock()


def _inside(saved, path: Path) -> Path:
    path = Path(path)
    if not path.resolve().is_relative_to(saved.path.resolve()):
        raise ValueError("Analysis path is outside the saved candidate folder")
    for parent in (path, *path.parents):
        if parent == saved.path:
            break
        if parent.is_symlink():
            raise ValueError("Saved candidate artifacts cannot be symbolic links")
    return path


def _atomic_json(saved, path: Path, data: dict) -> None:
    # Serialize identity verification and publication with candidate replacement.
    with _library_lock(saved.path.parent):
        _write_json(saved, path, data)


def _write_json(saved, path: Path, data: dict) -> None:
    _current(saved)
    _inside(saved, path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".analysis-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as target:
            json.dump(data, target, ensure_ascii=False, indent=2, allow_nan=False)
            target.flush()
            os.fsync(target.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def load_context(saved) -> dict:
    path = _inside(saved, saved.context_path)
    if not path.is_file():
        return {}
    if path.stat().st_size > 4 * 1024 * 1024:
        raise ValueError("Saved context is too large")
    data = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(data, dict) or data.get("snapshot_id") != saved.snapshot["snapshot_id"]
            or data.get("format") != "aptaswitch-saved-switch-context" or data.get("version") != 1
            or not isinstance(data.get("config"), dict)):
        raise ValueError("Invalid context for this saved candidate")
    return data


def _context_record(saved, config: dict, backend) -> dict:
    config = validate_analysis_config(saved.snapshot, config, backend=backend)
    context = build_saved_switch_context(saved.snapshot, config, backend=backend)
    data = {"format": "aptaswitch-saved-switch-context", "version": 1,
            "snapshot_id": saved.snapshot["snapshot_id"], "config": config, "context": context,
            "saved_at": datetime.now(timezone.utc).isoformat(),
            "thermodynamics": asdict(backend.provenance)}
    return data


def save_context(saved, config: dict, *, backend=None) -> dict:
    backend = backend if backend is not None else get_backend()
    data = _context_record(saved, config, backend)
    _atomic_json(saved, saved.context_path, data)
    return data


def configuration_key(saved, config: dict, *, backend=None) -> str:
    backend = backend if backend is not None else get_backend()
    return hashlib.sha256(json.dumps([saved.snapshot["snapshot_id"], config, asdict(backend.provenance)],
                                    sort_keys=True, allow_nan=False).encode()).hexdigest()


def saved_analyses(saved) -> list[dict]:
    root = _inside(saved, saved.analyses_dir)
    if not root.is_dir():
        return []
    records = []
    for path in sorted(root.glob("*/result.json"), reverse=True):
        try:
            _inside(saved, path)
            if path.stat().st_size > 32 * 1024 * 1024:
                continue
            data = json.loads(path.read_text(encoding="utf-8"))
            if (not isinstance(data, dict) or data.get("snapshot_id") != saved.snapshot["snapshot_id"]
                    or data.get("status") not in {"completed", "cancelled", "failed"}):
                continue
            data["result_path"] = str(path)
            records.append(data)
        except (OSError, ValueError, TypeError):
            continue
    # Run folders contain a second-resolution time followed by a random UUID;
    # their lexical order is not chronological when two runs finish quickly.
    return sorted(records, key=lambda row: (str(row.get("completed_at") or ""),
                                            str(row.get("run_id") or "")), reverse=True)


def active_analysis(saved):
    with _LOCK:
        job = _JOBS.get(str(saved.path.resolve()))
        return job if job and job.snapshot().state == "running" else None


def active_saved_analysis_for_design(design_path: Path | str):
    """Find a current design's worker after browser reload without disk scans."""
    root = Path(design_path).expanduser().resolve()
    with _LOCK:
        jobs = tuple(_JOBS.values())
    for job in jobs:
        if (job.saved.path.resolve().is_relative_to(root)
                and job.snapshot().state == "running"):
            return job
    return None


def start_analysis(saved, config: dict, import_path: str | None = None, *, backend=None):
    from .saved_switch_free_md_runtime import active_free_md
    _current(saved)
    if active_free_md(saved):
        raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
    backend = backend if backend is not None else get_backend(import_path)
    existing = active_analysis(saved)
    if existing:
        return existing, False
    context = _context_record(saved, config, backend)
    config = context["config"]
    with _library_lock(saved.path.parent), _LOCK:
        _current(saved)
        if active_free_md(saved):
            raise ValueError(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
        existing = active_analysis(saved)
        if existing:
            return existing, False
        # Publish under the library lock already held; do not nest flock.
        _write_json(saved, saved.context_path, context)
        job = SavedSwitchAnalysisJob(saved, config, import_path, backend=backend)
        _JOBS[str(saved.path.resolve())] = job
    try:
        job.start()
    except Exception:
        with _LOCK:
            if _JOBS.get(str(saved.path.resolve())) is job:
                _JOBS.pop(str(saved.path.resolve()), None)
        raise
    return job, True


class SavedSwitchAnalysisJob:
    def __init__(self, saved, config: dict, import_path: str | None = None, *, backend=None):
        _current(saved)
        self._backend = backend if backend is not None else get_backend(import_path)
        self.saved = saved
        self.config = json.loads(json.dumps(config, allow_nan=False))
        self.import_path = import_path
        self.key = configuration_key(saved, config, backend=self._backend)
        self.run_id = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_") + uuid.uuid4().hex[:10]
        self.folder = _inside(saved, saved.analyses_dir / self.run_id)
        self.folder.mkdir(parents=True, exist_ok=False)
        self._lock, self._cancel = threading.Lock(), threading.Event()
        self._state, self._result, self._error = "running", None, None
        self._current, self._total, self._message = 0, 1, ""
        self._log = []
        self._start, self._end, self._consumed = time.monotonic(), None, False
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)

    def start(self):
        self._thread.start()

    def cancel(self):
        self._cancel.set()

    def _progress(self, current, total, message):
        _current(self.saved)
        with self._lock:
            self._current, self._total, self._message = current, total, message
            self._log = [*self._log[-99:], str(message)]
        _atomic_json(self.saved, self.folder / "progress.json", {
            "snapshot_id": self.saved.snapshot["snapshot_id"], "status": "running",
            "current": current, "total": total, "message": message,
        })

    def _run(self):
        try:
            result = analyze_saved_switch(self.saved.snapshot, self.import_path, self.config,
                                         progress_callback=self._progress, cancel_check=self._cancel.is_set,
                                         backend=self._backend)
        except Exception as exc:
            result = {"status": "failed", "error": str(exc)}
        result.setdefault("config", self.config)
        result.setdefault("thermodynamics", asdict(self._backend.provenance))
        result.update(snapshot_id=self.saved.snapshot["snapshot_id"], configuration_key=self.key,
                      run_id=self.run_id, requested_config=self.config, completed_at=datetime.now(timezone.utc).isoformat())
        try:
            _atomic_json(self.saved, self.folder / "result.json", result)
        except Exception as exc:
            result.update(status="failed", error=str(exc))
        else:
            result["result_path"] = str(self.folder / "result.json")
        with self._lock:
            self._result, self._error = result, result.get("error")
            self._state = "done" if result.get("status") == "completed" else result.get("status", "failed")
            self._end = time.monotonic()
        with _LOCK:
            if _JOBS.get(str(self.saved.path.resolve())) is self:
                _JOBS.pop(str(self.saved.path.resolve()), None)

    def snapshot(self):
        with self._lock:
            return JobSnapshot(state=self._state, current=self._current, total=self._total,
                               message=self._message, result=self._result, error=self._error,
                               elapsed_seconds=(self._end or time.monotonic()) - self._start,
                               log=tuple(self._log))

    def consume_once(self):
        with self._lock:
            if self._state == "running" or self._consumed:
                return False
            self._consumed = True
            return True
