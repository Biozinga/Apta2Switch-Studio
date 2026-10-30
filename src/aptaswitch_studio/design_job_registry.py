"""Process-surviving design workers and strictly loaded recovery checkpoints."""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import os
import threading

from aptaswitch_core.models import RunResult
from aptaswitch_core.run_storage import run_output_directory
from aptaswitch_core.localization import tr


_JOBS, _RUN_JOBS, _LOCK = vars(threading).setdefault(
    '_aptaswitch_design_registry', ({}, {}, threading.RLock()))


def output_root(path: Path | str) -> Path:
    return run_output_directory(Path(path).expanduser().resolve(), 'design').parent


class DesignFileLock:
    """Nonblocking advisory lock, released by the OS if the app exits."""
    def __init__(self, path: Path):
        self.path = path
        self.handle = None

    def acquire(self, *, create=True) -> bool:
        if create:
            self.path.parent.mkdir(parents=True, exist_ok=True)
        elif not self.path.is_file():
            return True
        handle = self.path.open('a+b' if create else 'r+b')
        try:
            if os.name == 'nt':
                import msvcrt
                if self.path.stat().st_size == 0:
                    handle.write(b'\0')
                    handle.flush()
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except (BlockingIOError, OSError):
            handle.close()
            return False
        self.handle = handle
        return True

    def close(self):
        handle, self.handle = self.handle, None
        if handle is not None:
            try:
                if os.name == 'nt':
                    import msvcrt
                    handle.seek(0)
                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(handle, fcntl.LOCK_UN)
            finally:
                handle.close()


def find_design_job(path: Path | str):
    """Find the process-owned job by output root, category, run, or checkpoint."""
    canonical = Path(path).expanduser().resolve()
    if canonical.name == 'design_checkpoint.json':
        canonical = canonical.parent
    with _LOCK:
        return _RUN_JOBS.get(str(canonical)) or _JOBS.get(str(output_root(canonical)))


def register_design_job(job):
    """Atomically reserve one worker per output root; caller starts its thread."""
    root = output_root(job.request.output_dir)
    with _LOCK:
        previous = _JOBS.get(str(root))
        if previous is job:
            return False
        if previous is not None and previous.snapshot().state == 'running':
            raise ValueError(tr('Un calcul de conception est déjà en cours pour ce design.'))
        file_lock = DesignFileLock(root/'.design.lock')
        if not file_lock.acquire():
            raise ValueError(tr('Ce calcul de conception est déjà ouvert dans une autre instance de l’application.'))
        job._design_file_locks.append(file_lock)
        if previous is not None:
            for key, value in tuple(_RUN_JOBS.items()):
                if value is previous:
                    _RUN_JOBS.pop(key)
        _JOBS[str(root)] = job
        return True


def register_design_run(job, run_path: Path):
    canonical = Path(run_path).expanduser().resolve()
    expected_parent = run_output_directory(job.request.output_dir, 'design').resolve()
    if canonical.parent != expected_parent:
        raise ValueError(tr('Le checkpoint ne se trouve pas dans le dossier de résultats de ce design.'))
    with _LOCK:
        existing = _RUN_JOBS.get(str(canonical))
        if existing is job:
            return
        if existing is not None and existing.snapshot().state == 'running':
            raise ValueError(tr('Un calcul de conception utilise déjà ce checkpoint.'))
        file_lock = DesignFileLock(canonical/'design.worker.lock')
        if not file_lock.acquire():
            raise ValueError(tr('Ce checkpoint est déjà utilisé par une autre instance de l’application.'))
        job._design_file_locks.append(file_lock)
        _RUN_JOBS[str(canonical)] = job


def release_design_job(job, *, remove=False):
    with _LOCK:
        if remove:
            for registry in (_JOBS, _RUN_JOBS):
                for key, value in tuple(registry.items()):
                    if value is job:
                        registry.pop(key)
        for file_lock in reversed(job._design_file_locks):
            file_lock.close()
        job._design_file_locks.clear()


def completed_trials(result: RunResult) -> int:
    if getattr(result.request, "fixed_switch_sequence", None) is not None:
        count = result.optimization.get('evaluations', 0)
        return count if type(count) is int and count in (0, 1) else 0
    checkpoint = result.optimization.get('checkpoint', {})
    if isinstance(checkpoint, dict) and type(checkpoint.get('completed_trials')) is int:
        return checkpoint['completed_trials']
    return len(result.optimization.get('trials', []))


@dataclass(frozen=True)
class DesignCheckpoint:
    path: Path
    result: RunResult
    state: str

    @property
    def completed_trials(self) -> int:
        return completed_trials(self.result)

    @property
    def total_trials(self) -> int:
        return self.result.request.trials


def load_design_checkpoint(path: Path | str) -> DesignCheckpoint:
    from .run_library import load_run
    canonical = Path(path).expanduser().resolve()
    if canonical.is_dir():
        canonical = canonical/'design_checkpoint.json'
    if canonical.name != 'design_checkpoint.json':
        raise ValueError(tr('Ouvrez le checkpoint de conception à reprendre.'))
    loaded = load_run(canonical)
    if loaded.kind != 'design' or not isinstance(loaded.result, RunResult):
        raise ValueError(tr('Ce fichier ne contient pas un checkpoint de conception.'))
    result = loaded.result
    metadata = result.optimization.get('checkpoint')
    if (not isinstance(metadata, dict) or metadata.get('phase') not in
            {'design', 'inverse_extension', 'exports'} or result.status not in
            {'paused', 'running', 'completed', 'stopped', 'failed'}):
        raise ValueError(tr('Seul un calcul en pause ou interrompu peut être repris.'))
    # The reader validates scientific content; recovery also fixes ownership to
    # the actual file rather than trusting a serialized absolute output path.
    if Path(result.output_dir).resolve() != canonical.parent:
        raise ValueError(tr('Le checkpoint ne correspond pas à son dossier de calcul.'))
    if not isinstance(result.optimization.get('checkpoint'), dict):
        raise ValueError(tr('Ce résultat ne contient pas de checkpoint de conception reprenable.'))
    return DesignCheckpoint(canonical, result, 'paused' if result.status == 'paused' else 'interrupted')


def _locked_elsewhere(path: Path) -> bool:
    for lock_path in (path.parent/'design.worker.lock', path.parent.parent.parent/'.design.lock'):
        lock = DesignFileLock(lock_path)
        try:
            if not lock.acquire(create=False):
                return True
        except OSError:
            return True
        finally:
            lock.close()
    return False


def list_design_checkpoints(root: Path | str) -> list[DesignCheckpoint]:
    """Read valid paused/interrupted runs, excluding all actively owned files."""
    canonical = Path(root).expanduser().resolve()
    if canonical.name == 'design_checkpoint.json':
        sources = [canonical]
    elif (canonical/'design_checkpoint.json').is_file():
        sources = [canonical/'design_checkpoint.json']
    else:
        sources = run_output_directory(canonical, 'design').glob('*/design_checkpoint.json')
    found = []
    for source in sources:
        if source.is_symlink() or source.parent.is_symlink():
            continue
        current = find_design_job(source)
        if current is not None and current.snapshot().state == 'running':
            continue
        try:
            if _locked_elsewhere(source):
                continue
            found.append(load_design_checkpoint(source))
        except (OSError, ValueError, RecursionError):
            continue
    return sorted(found, key=lambda item: (item.result.started_at, item.result.run_id), reverse=True)


def resume_design_job(path: Path | str, *, finish_early=False, library_path=None):
    """Start from the persisted request; UI edits cannot change a continuation."""
    from .web_runtime import DesignRunJob
    checkpoint = load_design_checkpoint(path)
    result = checkpoint.result
    if library_path is not None:
        result.request = replace(result.request, thermodynamics_library_path=str(library_path))
    job = DesignRunJob(resume_result=result)
    job._resume_library_path = library_path
    if finish_early:
        job.finish_early()
    job.start()
    return job
