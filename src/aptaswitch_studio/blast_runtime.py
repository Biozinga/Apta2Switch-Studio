"""Explicit, cancellable BLAST jobs and portable local search records.

Constructing a job, loading a report and reopening a design never contacts NCBI.
Only ``start`` (or the explicit run helper) executes the client's request.
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
import threading
import time
import uuid
import zipfile
from contextvars import copy_context
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from aptaswitch_core.blast_specificity import BlastClient, BlastScreenRequest, Taxon
from aptaswitch_core.localization import get_language, language_context, tr

from .web_runtime import JobSnapshot


FORMAT = "aptaswitch-blast-screen"
FORMAT_VERSION = 1
_FINAL_STATUSES = frozenset({"completed", "partial", "failed", "cancelled"})
_MAX_REPORT_BYTES = 64 * 1024 * 1024
_HIT_COLUMNS = (
    "accession", "accessions", "subject_id", "title", "taxid", "scientific_name", "subject_length",
    "hit_number", "hsp_number", "query_from", "query_to", "hit_from", "hit_to",
    "query_strand", "hit_strand", "alignment_length", "identities", "identity_pct",
    "query_coverage_pct", "mismatches", "gaps", "evalue", "bit_score", "score",
    "longest_contiguous_match", "qseq", "hseq", "midline", "descriptions",
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value: Any) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=".blast-", suffix=".json", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, ensure_ascii=False, allow_nan=False, indent=2)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _empty_result(request: BlastScreenRequest, status: str, error: str | None = None) -> dict:
    result = {
        "status": status,
        "runs": [{"taxon": {"taxid": taxon.taxid, "name": taxon.name}, "scope": scope,
                  "status": status, "hits": [], **({"error": error} if error else {})}
                 for taxon in request.taxa for scope in request.scopes],
    }
    if error:
        result["error"] = error
    return result


def _write_table(result: dict, path: Path) -> None:
    rows = []
    for index, run in enumerate(result.get("runs", []), 1):
        taxon = run.get("taxon", {})
        hits = run.get("hits", [])
        base = {
            "search_index": index, "requested_taxid": taxon.get("taxid"),
            "requested_organism": taxon.get("name"), "scope": run.get("scope"),
            "search_status": run.get("status"), "rid": run.get("rid"), "database": run.get("database"),
            "query_length": run.get("query_length", len(result["request"]["sequence"])),
            "subject_hit_count": run.get("subject_hit_count"), "reported_hsp_count": len(hits),
            "search_error": run.get("error"),
        }
        # Empty rows explicitly carry completed/failed/cancelled status. A
        # failed request must never look like a successful zero-hit search.
        rows.extend({**base, **hit} for hit in hits)
        if not hits:
            rows.append(base)
    metadata = ("search_index", "requested_taxid", "requested_organism", "scope", "search_status", "rid",
                "database", "query_length", "subject_hit_count", "reported_hsp_count", "search_error")
    columns = list(dict.fromkeys([*metadata, *_HIT_COLUMNS, *(key for row in rows for key in row)]))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value
                             for key, value in row.items()})


def run_blast_and_export(
    request: BlastScreenRequest, output_dir: Path, *, client: Any = None,
    progress_callback: Callable[[int, int, str], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    """Execute an explicitly authorized screen and retain its complete evidence."""
    request = request.normalized()
    base = Path(output_dir).expanduser()
    if base.is_symlink() or (base / "blast").is_symlink():
        raise ValueError(tr("Le dossier des recherches BLAST ne peut pas être un lien symbolique."))
    root = base.resolve() / "blast"
    started = _now()
    run_id = f"{datetime.now(timezone.utc):%Y%m%d_%H%M%S}_blast_{uuid.uuid4().hex[:8]}"
    directory = root / run_id
    directory.mkdir(parents=True, exist_ok=False)
    _write_json(directory / "request.json", {
        "format": FORMAT, "format_version": FORMAT_VERSION, "run_id": run_id,
        "started_at": started, "request": request.to_dict(), "fingerprint": request.fingerprint(),
    })
    logs = []
    log_path = directory / "blast.log"

    def progress(current: int, total: int, message: str) -> None:
        line = f"[{_now()}] {message}"
        logs.append(line)
        with log_path.open("a", encoding="utf-8") as stream:
            stream.write(line + "\n")
        if progress_callback:
            progress_callback(current, total, message)

    total = len(request.taxa) * len(request.scopes)
    progress(0, total, tr("Préparation du contrôle BLAST…"))
    if cancel_check and cancel_check():
        result = _empty_result(request, "cancelled")
    else:
        try:
            result = (client if client is not None else BlastClient()).run(
                request, progress_callback=progress, cancel_check=cancel_check,
            )
            if not isinstance(result, dict) or result.get("status") not in _FINAL_STATUSES:
                raise ValueError(tr("Le service BLAST a renvoyé un état de recherche invalide."))
            runs = result.get("runs")
            if (not isinstance(runs, list)
                    or any(not isinstance(run, dict) or run.get("status") not in _FINAL_STATUSES
                           or not isinstance(run.get("hits"), list)
                           or any(not isinstance(hit, dict) for hit in run["hits"]) for run in runs)
                    or result["status"] == "completed" and (
                        len(runs) != total or any(run["status"] != "completed" for run in runs))):
                raise ValueError(tr("Le service BLAST a renvoyé un état de recherche invalide."))
            # Detach result data from a client or test fixture before adding
            # local metadata and replacing raw reports with relative filenames.
            result = json.loads(json.dumps(result, ensure_ascii=False, allow_nan=False))
        except Exception as exc:
            result = _empty_result(request, "failed", str(exc))
    result.update({
        "format": FORMAT, "format_version": FORMAT_VERSION,
        "run_id": run_id, "started_at": started, "completed_at": _now(),
        "request": request.to_dict(), "fingerprint": request.fingerprint(),
        "output_dir": str(directory),
    })
    exports = {"request": str(directory / "request.json"), "log": str(log_path)}
    for index, run in enumerate(result.get("runs", []), 1):
        raw = run.pop("raw_report", None)
        if raw is None:
            continue
        extension = {"json": "json", "xml": "xml", "text": "txt", "txt": "txt"}.get(run.get("raw_format"), "txt")
        filename = f"search_{index:03d}.{extension}"
        if not isinstance(raw, str):
            raw = json.dumps(raw, ensure_ascii=False, indent=2)
        (directory / filename).write_text(raw, encoding="utf-8")
        run["raw_report_file"] = filename
        exports[f"raw_{index}"] = str(directory / filename)
    final_message = {
        "completed": tr("Contrôle BLAST terminé."),
        "partial": tr("Contrôle BLAST incomplet : certaines recherches ont échoué."),
        "failed": tr("Le contrôle BLAST a échoué."),
        "cancelled": tr("Contrôle BLAST annulé."),
    }[result["status"]]
    completed = sum(run.get("status") == "completed" for run in result.get("runs", []))
    progress(completed, total, final_message)
    if result.get("error"):
        progress(completed, total, str(result["error"]))
    result["log"] = logs
    _write_json(directory / "results.json", result)
    exports["json"] = str(directory / "results.json")
    _write_table(result, directory / "hits.csv")
    exports["csv"] = str(directory / "hits.csv")
    archive_path = directory / f"{run_id}.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in exports.values():
            archive.write(path, arcname=Path(path).name)
    exports["zip"] = str(archive_path)
    return result, exports


def load_blast_result(path: Path) -> tuple[dict, dict[str, str]]:
    """Read a local terminal report and relocate its exports without any network."""
    source = Path(path).expanduser()
    if source.is_dir():
        source = source / "results.json"
    if source.is_symlink() or not source.is_file() or source.stat().st_size > _MAX_REPORT_BYTES:
        raise ValueError(tr("Le rapport BLAST enregistré est absent ou invalide."))
    source = source.resolve()
    try:
        result = json.loads(source.read_text(encoding="utf-8"))
        if (not isinstance(result, dict) or result.get("format") != FORMAT
                or result.get("format_version") != FORMAT_VERSION or result.get("status") not in _FINAL_STATUSES
                or not isinstance(result.get("runs"), list)):
            raise ValueError
        recorded = result["request"]
        request = BlastScreenRequest(
            sequence=recorded["sequence"], taxa=tuple(Taxon(**row) for row in recorded["taxa"]),
            scopes=tuple(recorded["scopes"]), max_hits=recorded["max_hits"],
        ).normalized()
        if result.get("fingerprint") != request.fingerprint():
            raise ValueError
        for key in ("started_at", "completed_at"):
            if not isinstance(result.get(key), str) or datetime.fromisoformat(result[key]).tzinfo is None:
                raise ValueError
        for run in result["runs"]:
            if (not isinstance(run, dict) or run.get("status") not in _FINAL_STATUSES
                    or not isinstance(run.get("hits"), list)
                    or any(not isinstance(hit, dict) for hit in run["hits"])):
                raise ValueError
        if result["status"] == "completed" and (
            len(result["runs"]) != len(request.taxa) * len(request.scopes)
            or any(run["status"] != "completed" for run in result["runs"])
        ):
            raise ValueError
    except (KeyError, TypeError, ValueError, RecursionError) as exc:
        raise ValueError(tr("Le rapport BLAST enregistré est absent ou invalide.")) from exc
    directory = source.parent
    result["output_dir"] = str(directory)
    exports = {"json": str(source)}
    filenames = {"request": "request.json", "log": "blast.log", "csv": "hits.csv"}
    run_id = result.get("run_id")
    if isinstance(run_id, str) and Path(run_id).name == run_id and run_id not in {".", ".."}:
        filenames["zip"] = f"{run_id}.zip"
    for index, run in enumerate(result["runs"], 1):
        filename = run.get("raw_report_file")
        if isinstance(filename, str) and Path(filename).name == filename and filename not in {".", ".."}:
            filenames[f"raw_{index}"] = filename
    for key, filename in filenames.items():
        artifact = directory / filename
        if artifact.is_file() and not artifact.is_symlink() and artifact.resolve().parent == directory:
            exports[key] = str(artifact)
    return result, exports


class BlastScreenJob:
    """Run a submitted screen in the background; cancellation stops future calls."""

    def __init__(self, request: BlastScreenRequest, output_dir: Path, *,
                 client_factory: Callable[[], Any] | None = None):
        self.request, self.output_dir = request.normalized(), Path(output_dir)
        self._client_factory = client_factory or BlastClient
        self._language = get_language()
        self._lock, self._cancel = threading.Lock(), threading.Event()
        self._state, self._current = "running", 0
        self._total = len(self.request.taxa) * len(self.request.scopes)
        self._message = tr("Préparation du contrôle BLAST…")
        self._started, self._finished = time.monotonic(), None
        self._result, self._exports, self._error = None, None, None
        self._log: list[str] = []
        self._consumed = False
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)

    def start(self) -> None:
        self._thread.start()

    def cancel(self) -> None:
        self._cancel.set()
        with self._lock, language_context(self._language):
            if self._state == "running":
                self._message = tr("Annulation BLAST demandée ; aucune nouvelle requête ne sera envoyée.")
                self._log.append(self._message)

    def _progress(self, current: int, total: int, message: str) -> None:
        with self._lock:
            self._current, self._total, self._message = current, total, message
            self._log.append(f"[{datetime.now().astimezone():%H:%M:%S}] {message}")
            self._log = self._log[-500:]

    def _run(self) -> None:
        try:
            # Even the client's construction is postponed until explicit start.
            client = None if self._cancel.is_set() else self._client_factory()
            result, exports = run_blast_and_export(
                self.request, self.output_dir, client=client,
                progress_callback=self._progress, cancel_check=self._cancel.is_set,
            )
        except Exception as exc:
            result, exports = _empty_result(self.request, "failed", str(exc)), None
            self._progress(self._current, self._total, str(exc))
        with self._lock:
            self._result, self._exports = result, exports
            self._error = result.get("error")
            self._state = "failed" if result["status"] == "failed" else "done"
            self._finished = time.monotonic()

    def snapshot(self) -> JobSnapshot:
        with self._lock:
            return JobSnapshot(
                state=self._state, current=self._current, total=self._total, message=self._message,
                result=self._result, exports=dict(self._exports) if self._exports else None,
                error=self._error, elapsed_seconds=max(0., (self._finished or time.monotonic()) - self._started),
                eta_seconds=0. if self._finished is not None else None, log=tuple(self._log),
            )

    def consume_once(self) -> bool:
        with self._lock:
            if self._consumed or self._state == "running":
                return False
            self._consumed = True
            return True
