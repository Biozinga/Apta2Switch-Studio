"""Cancellable initial-system equilibrium analysis, independent of UI reruns."""

from __future__ import annotations

from contextvars import copy_context
from datetime import datetime
import threading
import time

from aptaswitch_core.localization import get_language, language_context, tr
from aptaswitch_core.tube_analysis import analyze_parasite_tubes
from aptaswitch_core.thermodynamics.setup import get_backend

from .web_runtime import JobSnapshot


class SystemTubeJob:
    def __init__(self, signature: tuple, *, backend=None):
        self.signature = signature
        self._backend = backend if backend is not None else get_backend(signature[2])
        self._lock = threading.Lock()
        self._cancel = threading.Event()
        self._language = get_language()
        self._state = "running"
        self._result = None
        self._error = None
        self._current, self._total = 0, 2
        self._message = tr("Préparation des deux tubes thermodynamiques…")
        self._log = []
        self._started = time.monotonic()
        self._finished = None
        self._consumed = False
        self._progress(0, 2, self._message)
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)

    def start(self) -> None:
        self._thread.start()

    def cancel(self) -> None:
        if not self._cancel.is_set():
            self._cancel.set()
            with self._lock, language_context(self._language):
                if self._state == "running":
                    self._message = tr("Annulation demandée — fin du calcul natif en cours…")
                    self._log.append(self._message)

    def _progress(self, current, total, message) -> None:
        with self._lock:
            self._current, self._total, self._message = current, total, message
            self._log.append(f"[{datetime.now().astimezone():%H:%M:%S}] {message}")

    def _run(self) -> None:
        aptamer, trigger, import_path, thermo, conditions = self.signature[:5]
        try:
            result = analyze_parasite_tubes(
                import_path, aptamer, trigger, thermo, conditions,
                cancel_check=self._cancel.is_set, progress_callback=self._progress, backend=self._backend,
            )
        except Exception as exc:
            with self._lock:
                self._state, self._error = "failed", str(exc)
                self._finished = time.monotonic()
            return
        with self._lock:
            self._result, self._state = result, "done"
            self._finished = time.monotonic()

    def snapshot(self) -> JobSnapshot:
        with self._lock:
            return JobSnapshot(
                state=self._state, current=self._current, total=self._total,
                message=self._message, result=self._result, error=self._error,
                elapsed_seconds=(self._finished or time.monotonic()) - self._started,
                log=tuple(self._log),
            )

    def consume_once(self) -> bool:
        with self._lock:
            if self._consumed or self._state == "running":
                return False
            self._consumed = True
            return True
