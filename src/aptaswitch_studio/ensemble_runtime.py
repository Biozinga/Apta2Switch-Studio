"""On-demand ensemble analysis, saved even if the user leaves Results."""
from __future__ import annotations

from contextvars import copy_context
import threading
import time

from aptaswitch_core.structure_ensemble import EnsembleCancelled, analyze_structure_ensemble, save_ensemble
from aptaswitch_core.thermodynamics.setup import ENGINE_ID, backend_for_request
from .web_runtime import JobSnapshot


class EnsembleJob:
    def __init__(self, result, candidate, import_path=None, *, num_sample=1000, save_result=None, backend=None):
        if getattr(result.request, "engine", "") != ENGINE_ID:
            raise ValueError("This calculation requires the native-hybrid engine.")
        result.request.thermo.validate_model()
        self._result_source, self._candidate = result, candidate
        self._backend = backend if backend is not None else backend_for_request(
            result.request, library_path=import_path or getattr(result.request, "thermodynamics_library_path", None))
        self._import_path, self._num_sample = import_path, num_sample
        self._save_result = save_result
        self._lock, self._cancel = threading.Lock(), threading.Event()
        self._state, self._result, self._error = "running", None, None
        self._current, self._total, self._message = 0, 4, ""
        self._start, self._end = time.monotonic(), None
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)

    def start(self):
        self._thread.start()

    def cancel(self):
        self._cancel.set()

    def _progress(self, current, total, message):
        with self._lock:
            self._current, self._total, self._message = current, total, message

    def _run(self):
        try:
            result = analyze_structure_ensemble(
                self._result_source.request, self._candidate, self._import_path, num_sample=self._num_sample,
                cancel_check=self._cancel.is_set, progress_callback=self._progress, backend=self._backend,
            )
            if self._cancel.is_set():
                raise EnsembleCancelled()
            if self._save_result is None:
                save_ensemble(self._result_source.output_dir, result)
            else:
                self._save_result(result)
        except EnsembleCancelled:
            with self._lock:
                self._state, self._end = "cancelled", time.monotonic()
        except Exception as exc:
            with self._lock:
                self._state, self._error, self._end = "failed", str(exc), time.monotonic()
        else:
            with self._lock:
                self._state, self._result, self._end = "done", result, time.monotonic()

    def snapshot(self):
        with self._lock:
            return JobSnapshot(state=self._state, current=self._current, total=self._total,
                               message=self._message, result=self._result, error=self._error,
                               elapsed_seconds=(self._end or time.monotonic()) - self._start)
