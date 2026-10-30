"""Background jobs and small helpers used by the local Streamlit interface."""

from __future__ import annotations

from aptaswitch_core.localization import get_language, language_context, tr

import io
import tempfile
import threading
import time
import zipfile
from contextvars import copy_context
from copy import deepcopy
from dataclasses import dataclass, replace
from datetime import datetime
from decimal import Decimal
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

from aptaswitch_core.run_storage import default_export_dir, legacy_export_dir
from aptaswitch_core.workflow import run_and_export
from .design_job_registry import (
    completed_trials, find_design_job, list_design_checkpoints, load_design_checkpoint,
    register_design_job, register_design_run, release_design_job, resume_design_job,
)


CONCENTRATION_UNITS = ("nM", "µM", "mM", "M")
_CONCENTRATION_FACTORS = {
    "nM": 1e-9,
    "µM": 1e-6,
    "uM": 1e-6,
    "μM": 1e-6,
    "mM": 1e-3,
    "M": 1.0,
}


def concentration_to_m(value: float, unit: str) -> float:
    """Convert a displayed concentration to mol/L, consistently across units."""

    try:
        factor = _CONCENTRATION_FACTORS[unit]
    except KeyError as exc:
        raise ValueError(tr('Unité de concentration inconnue : {v0}', v0=unit)) from exc
    # Decimal input keeps equivalent unit choices on the same float value
    # (e.g. 7 mM == 7000 µM), including keys for cached scientific analyses.
    return float(Decimal(str(float(value))) * Decimal(str(factor)))


def concentration_from_m(value_m: float, unit: str) -> float:
    """Convert mol/L to a selected display unit."""

    try:
        factor = _CONCENTRATION_FACTORS[unit]
    except KeyError as exc:
        raise ValueError(tr('Unité de concentration inconnue : {v0}', v0=unit)) from exc
    # Limit harmless binary floating-point tails before placing the value back
    # into a visible number input (e.g. show 5000, not 4999.999999999999).
    return float(f"{float(value_m) / factor:.12g}")


def preferred_concentration_unit(value_m: float) -> str:
    """Choose a readable unit for a molar concentration."""

    magnitude = abs(float(value_m))
    if magnitude >= 1.0:
        return "M"
    if magnitude >= 1e-3:
        return "mM"
    if magnitude >= 1e-6:
        return "µM"
    return "nM"


@lru_cache(maxsize=32)
def svg_to_png_bytes(svg: str, zoom: float = 2.0) -> bytes:
    """Rasterize one generated SVG into a high-resolution PNG."""

    import resvg_py

    png = resvg_py.svg_to_bytes(
        svg_string=svg,
        background="#ffffff",
        zoom=max(0.1, float(zoom)),
    )
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(tr('La conversion de la structure en PNG a échoué.'))
    return png


def normalized_export_dir(path_value: str | Path, app_dir: Path | None = None) -> Path:
    """Resolve blank/relative paths within the application's shared run library."""

    base = default_export_dir(app_dir)
    text = str(path_value).strip()
    if not text:
        return base
    path = Path(text).expanduser()
    return path if path.is_absolute() else base / path


def ensure_writable_export_dir(
    requested: Path,
    *,
    fallback: Path | None = None,
) -> tuple[Path, bool]:
    """Probe the requested root, then the app library and a writable user fallback."""

    candidates = [Path(requested)]
    fallbacks = [Path(fallback)] if fallback is not None else [default_export_dir(), legacy_export_dir()]
    for candidate in fallbacks:
        if candidate not in candidates:
            candidates.append(candidate)
    last_error: OSError | None = None
    for index, candidate in enumerate(candidates):
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(prefix=".aptaswitch-write-test-", dir=candidate):
                pass
        except OSError as exc:
            last_error = exc
            continue
        return candidate, index > 0
    if last_error is not None:
        raise last_error
    raise OSError(tr("Aucun dossier d'export disponible"))


@dataclass(frozen=True)
class JobSnapshot:
    """Immutable view of a background job, safe to read during a rerun."""

    state: str
    current: int
    total: int
    message: str
    result: Any = None
    exports: dict[str, str] | None = None
    error: str | None = None
    elapsed_seconds: float = 0.0
    eta_seconds: float | None = None
    log: tuple[str, ...] = ()
    control_requested: str | None = None
    phase: str = ''
    can_request_control: bool = False


class DesignRunJob:
    """Run a design without blocking Streamlit's single local UI session."""

    def __init__(self, request=None, *, analysis_source=None, resume_result=None):
        if resume_result is not None:
            if analysis_source is not None:
                raise ValueError('Inverse reanalysis cannot resume switch-design trials.')
            if request is not None and request.to_dict() != resume_result.request.to_dict():
                raise ValueError('A continuation must use its saved design request.')
            request = resume_result.request
        if request is None:
            raise ValueError('A design request is required.')
        self.request = deepcopy(request)
        self._analysis_source = analysis_source
        self._resume_result = deepcopy(resume_result)
        self._resume_library_path = None
        self._design_file_locks = []
        self._thread_started = False
        self._control_requested = None
        self._phase = (resume_result.optimization.get('checkpoint', {}).get('phase', 'design')
                       if resume_result is not None else
                       'design' if analysis_source is None else 'inverse_extension')
        self._checkpoint_path = (Path(resume_result.output_dir).resolve()/'design_checkpoint.json'
                                 if resume_result is not None else None)
        self._language = get_language()
        self._lock = threading.Lock()
        self._cancel = threading.Event()
        self._state = "running"
        self._current = completed_trials(resume_result) if resume_result is not None else 0
        self._total = request.trials
        self._message = tr('Préparation du calcul…')
        self._started_monotonic = time.monotonic()
        self._phase_started_monotonic = self._started_monotonic
        self._phase_start_current = self._current
        self._finished_monotonic: float | None = None
        self._log: list[str] = []
        self._result = None
        self._exports: dict[str, str] | None = None
        self._error: str | None = None
        self._consumed = False
        # Preserve the language selected when this job was created.
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)
        self._append_log_locked(self._message)

    def _append_log_locked(self, message: str) -> None:
        timestamp = datetime.now().astimezone().strftime("%H:%M:%S")
        self._log.append(f"[{timestamp}] {message}")
        # Bound the console history for long design runs.
        if len(self._log) > 500:
            self._log = self._log[-500:]

    def start(self) -> None:
        if self._thread_started:
            return
        if not register_design_job(self):
            return
        try:
            if self._resume_result is not None:
                register_design_run(self, self._checkpoint_path.parent)
                # Read again after both interprocess locks are acquired. Another
                # app may have advanced or completed the file since the click.
                persisted = load_design_checkpoint(self._checkpoint_path).result
                if self._resume_library_path is not None:
                    persisted.request = replace(persisted.request,
                        thermodynamics_library_path=str(self._resume_library_path))
                self._resume_result = deepcopy(persisted)
                self.request = deepcopy(persisted.request)
                self._current = completed_trials(persisted)
                self._total = self.request.trials
                self._phase = persisted.optimization.get('checkpoint', {}).get('phase', 'design')
                self._phase_start_current = self._current
            self._thread_started = True
            self._thread.start()
        except Exception:
            self._thread_started = False
            release_design_job(self, remove=True)
            raise

    @property
    def supports_continuation(self) -> bool:
        with self._lock:
            return self._analysis_source is None

    @property
    def control_requested(self) -> str | None:
        with self._lock:
            return self._control_requested

    @property
    def checkpoint_path(self) -> Path | None:
        with self._lock:
            return self._checkpoint_path

    def _request_control(self, action: str) -> None:
        with self._lock, language_context(self._language):
            if (self._state != 'running' or self._analysis_source is not None
                    or (self._thread_started and self._phase != 'design') or self._cancel.is_set()):
                return
            # Finishing is an explicit terminal request and wins over pause.
            if self._control_requested == 'finish':
                return
            self._control_requested = action
            self._message = (tr('Pause demandée — l’essai en cours sera terminé et sauvegardé.')
                             if action == 'pause' else
                             tr('Fin demandée — l’essai en cours sera terminé et sauvegardé.'))
            self._append_log_locked(self._message)

    def pause(self) -> None:
        self._request_control('pause')

    def finish_early(self) -> None:
        self._request_control('finish')

    def resume(self, *, finish_early=False, library_path=None):
        path = self.checkpoint_path
        if path is None:
            raise ValueError(tr('Aucun checkpoint de conception n’est disponible.'))
        return resume_design_job(path, finish_early=finish_early, library_path=library_path)

    def _control_check(self):
        with self._lock:
            return self._control_requested

    def _checkpoint(self, result) -> None:
        register_design_run(self, result.output_dir)
        with self._lock:
            self._checkpoint_path = Path(result.output_dir).resolve()/'design_checkpoint.json'
            self._result = deepcopy(result)
            self._current = completed_trials(result)
            self._total = result.request.trials
            self._phase = result.optimization.get('checkpoint', {}).get('phase', 'design')

    def cancel(self) -> None:
        self._cancel.set()
        with self._lock, language_context(self._language):
            if self._state == "running":
                self._message = tr('Annulation demandée…')
                self._append_log_locked(self._message)

    def _progress(self, current: int, total: int, message: str) -> None:
        with self._lock:
            if current < self._current or total != self._total:
                self._phase_started_monotonic = time.monotonic()
                self._phase_start_current = current
            self._current = current
            self._total = total
            self._message = message
            self._append_log_locked(message)

    def _run(self) -> None:
        try:
            self._run_calculation()
        finally:
            release_design_job(self)

    def _run_calculation(self) -> None:
        try:
            if self._analysis_source is None:
                result, exports = run_and_export(
                    self.request, progress_callback=self._progress, cancel_check=self._cancel.is_set,
                    control_check=self._control_check, resume_result=self._resume_result,
                    checkpoint_callback=self._checkpoint,
                )
            else:
                from aptaswitch_core.workflow import reanalyze_inverse_and_export
                result, exports = reanalyze_inverse_and_export(
                    self._analysis_source, self.request,
                    progress_callback=self._progress, cancel_check=self._cancel.is_set,
                )
        except Exception as exc:  # pragma: no cover - defensive thread boundary
            with self._lock:
                self._state = "failed"
                self._error = str(exc)
                self._message = tr('Le calcul a échoué.')
                self._finished_monotonic = time.monotonic()
                self._append_log_locked(tr('ERREUR · {v0}', v0=exc))
            return
        with self._lock:
            self._state = "paused" if result.status == "paused" else "done"
            self._result = result
            self._exports = exports
            self._current = (completed_trials(result) if self._analysis_source is None else len(result.candidates))
            self._total = self.request.trials
            self._message = (
                tr('Calcul annulé — résultats partiels exportés.')
                if result.status == "cancelled"
                else tr('Calcul en pause — les essais terminés sont sauvegardés.')
                if result.status == 'paused'
                else tr('Calcul terminé avec les essais déjà réalisés.')
                if result.status == 'stopped'
                else tr('Calcul terminé.')
            )
            self._finished_monotonic = time.monotonic()
            self._append_log_locked(self._message)

    def snapshot(self) -> JobSnapshot:
        with self._lock:
            stopped_at = self._finished_monotonic or time.monotonic()
            elapsed = max(0.0, stopped_at - self._started_monotonic)
            eta = None
            phase_progress = self._current - self._phase_start_current
            if self._state == "running" and phase_progress > 0 and self._current < self._total:
                phase_elapsed = max(0.0, stopped_at - self._phase_started_monotonic)
                eta = max(0.0, phase_elapsed / phase_progress * (self._total - self._current))
            elif self._state != "running" or self._current >= self._total:
                eta = 0.0
            return JobSnapshot(
                state=self._state,
                current=self._current,
                total=self._total,
                message=self._message,
                result=self._result,
                exports=dict(self._exports) if self._exports else None,
                error=self._error,
                elapsed_seconds=elapsed,
                eta_seconds=eta,
                log=tuple(self._log),
                control_requested=self._control_requested,
                phase=self._phase,
                can_request_control=(self._state == 'running' and self._analysis_source is None
                                     and self._phase == 'design' and not self._cancel.is_set()),
            )

    def consume_once(self) -> bool:
        """Return true once when the completed result should enter session state."""
        with self._lock:
            if self._consumed or self._state not in {"done", "paused", "failed"}:
                return False
            self._consumed = True
            return True


def exports_zip_bytes(paths: Mapping[str, str]) -> bytes:
    """Bundle files produced by a run for one browser download."""
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path_text in paths.values():
            path = Path(path_text)
            if path.is_file():
                archive.write(path, arcname=path.name)
    return payload.getvalue()
