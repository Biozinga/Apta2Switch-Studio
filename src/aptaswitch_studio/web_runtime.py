"""Background jobs and small helpers used by the local Streamlit interface."""

from __future__ import annotations

from aptaswitch_core.localization import get_language, language_context, tr

import io
import tempfile
import threading
import time
import zipfile
from contextvars import copy_context
from dataclasses import dataclass
from datetime import datetime
from decimal import Decimal
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

from aptaswitch_core.run_storage import default_export_dir, legacy_export_dir
from aptaswitch_core.workflow import run_and_export


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


class DesignRunJob:
    """Run a design without blocking Streamlit's single local UI session."""

    def __init__(self, request):
        self.request = request
        self._language = get_language()
        self._lock = threading.Lock()
        self._cancel = threading.Event()
        self._state = "running"
        self._current = 0
        self._total = request.trials
        self._message = tr('Préparation du calcul…')
        self._started_monotonic = time.monotonic()
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
        self._thread.start()

    def cancel(self) -> None:
        self._cancel.set()
        with self._lock, language_context(self._language):
            if self._state == "running":
                self._message = tr('Annulation demandée…')
                self._append_log_locked(self._message)

    def _progress(self, current: int, total: int, message: str) -> None:
        with self._lock:
            self._current = current
            self._total = total
            self._message = message
            self._append_log_locked(message)

    def _run(self) -> None:
        try:
            result, exports = run_and_export(
                self.request,
                progress_callback=self._progress,
                cancel_check=self._cancel.is_set,
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
            self._state = "done"
            self._result = result
            self._exports = exports
            self._current = len(result.candidates)
            self._message = (
                tr('Calcul annulé — résultats partiels exportés.')
                if result.status == "cancelled"
                else tr('Calcul terminé.')
            )
            self._finished_monotonic = time.monotonic()
            self._append_log_locked(self._message)

    def snapshot(self) -> JobSnapshot:
        with self._lock:
            stopped_at = self._finished_monotonic or time.monotonic()
            elapsed = max(0.0, stopped_at - self._started_monotonic)
            eta = None
            if self._state == "running" and 0 < self._current < self._total:
                eta = max(0.0, elapsed / self._current * (self._total - self._current))
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
            )

    def consume_once(self) -> bool:
        """Return true once when the completed result should enter session state."""
        with self._lock:
            if self._consumed or self._state not in {"done", "failed"}:
                return False
            self._consumed = True
            return True


def contiguous_ranges(indices) -> list[tuple[int, int]]:
    """Collapse zero-based indices into inclusive contiguous ranges."""
    ranges: list[list[int]] = []
    for value in sorted(set(indices)):
        if ranges and value == ranges[-1][1] + 1:
            ranges[-1][1] = value
        else:
            ranges.append([value, value])
    return [(start, end) for start, end in ranges]


def largest_unbound_window(length: int, bound_indices) -> tuple[int, int] | None:
    """Return the longest aptamer-free trigger window as (start, length)."""
    bound = set(bound_indices)
    best_start, best_length = 0, 0
    cursor = 0
    while cursor < length:
        if cursor in bound:
            cursor += 1
            continue
        end = cursor
        while end < length and end not in bound:
            end += 1
        if end - cursor > best_length:
            best_start, best_length = cursor, end - cursor
        cursor = end
    if best_length < 2:
        return None
    return best_start, best_length


def exports_zip_bytes(paths: Mapping[str, str]) -> bytes:
    """Bundle files produced by a run for one browser download."""
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path_text in paths.values():
            path = Path(path_text)
            if path.is_file():
                archive.write(path, arcname=path.name)
    return payload.getvalue()
