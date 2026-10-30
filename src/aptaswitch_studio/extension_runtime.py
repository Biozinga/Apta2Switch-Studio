"""Background execution and reproducible exports for trigger extension searches."""

from __future__ import annotations

from aptaswitch_core.localization import get_language, language_context, tr

import csv
import json
import threading
import time
import zipfile
from contextvars import copy_context
from dataclasses import asdict, dataclass, field
from datetime import datetime
from numbers import Integral
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from aptaswitch_core.models import ProgressCallback, ThermoConditions
from aptaswitch_core.extension_thermodynamics import extension_model_settings
from aptaswitch_core.thermodynamics.setup import ENGINE_ID, get_backend
from aptaswitch_core.run_storage import new_run_directory
from aptaswitch_core.sequences import normalize_dna
from aptaswitch_core.extension_preview import ligand_positions
from aptaswitch_core.extension_visuals import export_extension_visuals
from aptaswitch_core.tube_analysis import TubeAnalysisConditions, analyze_parasite_tubes
from aptaswitch_core.trigger_extension import (
    design_trigger_extensions,
    estimate_extension_search,
    extension_split,
    EXTENSION_SIDES,
)

from .web_runtime import JobSnapshot
from .scientific_engines_ui import result_engine_label


@dataclass(frozen=True)
class ExtensionRunRequest:
    """Inputs to the thermodynamic extension search at one or both ends of the core."""

    aptamer_dna: str
    trigger_dna: str
    target_lengths: tuple[int, ...]
    thermo: ThermoConditions
    import_path: str
    output_dir: Path
    molecule_name: str = "Système personnalisé"
    top: int = 100
    max_ensemble_candidates: int = 50000
    extension_side: str = "5prime"
    ligand_aptamer_positions: tuple[int, ...] = ()
    parasite_top: int = 5
    tube_conditions: TubeAnalysisConditions = field(default_factory=TubeAnalysisConditions)
    thermodynamics_library_path: str | None = None
    engine: str = ENGINE_ID

    def __post_init__(self) -> None:
        if self.engine != ENGINE_ID:
            raise ValueError("This calculation requires the native-hybrid engine.")
        self.thermo.validate_model()
        self.thermo.validate_system_chemistry()
        if isinstance(self.parasite_top, bool) or not isinstance(self.parasite_top, Integral) or self.parasite_top < 0:
            raise ValueError(tr('Le nombre de candidats à analyser par longueur doit être un entier positif ou nul.'))
        object.__setattr__(self, "parasite_top", int(self.parasite_top))

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["target_lengths"] = list(self.target_lengths)
        data["ligand_aptamer_positions"] = list(self.ligand_aptamer_positions)
        data["output_dir"] = str(self.output_dir)
        return data


def request_from_extension_result(result: dict, *, import_path: str = "", output_dir: Path,
                                  thermodynamics_library_path: str | None = None) -> ExtensionRunRequest:
    """Resume saved sequences/physical inputs with the current catalog ceiling."""
    data = dict(result["request"])
    if data.get("engine") != ENGINE_ID:
        raise ValueError("This calculation requires the native-hybrid engine.")
    thermo = dict(data["thermo"])
    for field in ("screening_max_complex_size", "validation_max_complex_size"):
        if type(thermo.get(field)) is int and thermo[field] in (4, 5):
            thermo[field] = 3
    data["thermo"] = ThermoConditions(**thermo)
    tube_conditions = dict(data.get("tube_conditions", {}))
    if type(tube_conditions.get("max_complex_size")) is int and tube_conditions["max_complex_size"] in (4, 5):
        tube_conditions["max_complex_size"] = 3
    data["tube_conditions"] = TubeAnalysisConditions(**tube_conditions)
    data["target_lengths"] = tuple(data["target_lengths"])
    data["ligand_aptamer_positions"] = tuple(data.get("ligand_aptamer_positions", ()))
    data.update(import_path="", output_dir=output_dir,
                thermodynamics_library_path=thermodynamics_library_path)
    return ExtensionRunRequest(**data)


_CANDIDATE_COLUMNS = (
    "rank", "rank_within_length", "final_score", "trigger_length", "extension_length",
    "extension_side", "extension_5prime", "extension_3prime", "extension_5prime_length", "extension_3prime_length",
    "extended_trigger", "mfe_structure", "mfe_energy_kcal_mol",
    "core_pair_probability_rmse", "core_pair_probability_mae",
    "reference_pair_probability_loss", "extension_pair_probability_sum",
    "extension_to_aptamer_probability_sum", "extension_to_trigger_core_probability_sum",
    "trigger_only_parasite_pct", "mixed_trigger_parasite_pct", "mixed_trigger_target_pct",
)


def _timestamped(message: str) -> str:
    return f"[{datetime.now().astimezone().strftime('%H:%M:%S')}] {message}"


def _candidate_table(candidates: Sequence[Mapping[str, Any]]) -> tuple[list[str], list[list[Any]]]:
    """Keep all computed metrics, including nested probability data, in exports."""
    keys = set().union(*(candidate.keys() for candidate in candidates)) if candidates else set()
    columns = list(_CANDIDATE_COLUMNS)
    columns.extend(sorted(keys.difference(columns)))
    rows = [[_cell_value(candidate.get(column)) for column in columns] for candidate in candidates]
    return columns, rows


def _cell_value(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return value


def _tube_composition_table(result: Mapping[str, Any]) -> tuple[list[str], list[list[Any]]]:
    """Flatten every saved population with its sample and experimental conditions."""
    records: list[dict[str, Any]] = []
    reference = result.get("reference", {})
    samples = [("original", reference, result.get("input", {}).get("trigger_dna", ""))]
    samples.extend(("candidate", candidate, candidate.get("extended_trigger", "")) for candidate in result.get("candidates", []))
    for kind, sample, sequence in samples:
        analysis = sample.get("tube_analysis")
        if not isinstance(analysis, Mapping):
            continue
        common = {
            "sample": kind, "rank": sample.get("rank"),
            "rank_within_length": sample.get("rank_within_length"),
            "trigger_length": sample.get("trigger_length", len(sequence)),
            "trigger_sequence": analysis.get("trigger_sequence", sequence),
            "aptamer_sequence": analysis.get("aptamer_sequence", result.get("input", {}).get("aptamer_dna", "")),
            "analysis_status": analysis.get("status"), "error": analysis.get("error"),
            **analysis.get("conditions", {}), **analysis.get("model", {}),
        }
        tubes = analysis.get("tubes", {})
        if not tubes:
            records.append(common)
        for name, tube in tubes.items():
            for population in tube.get("complexes", []):
                records.append({
                    **common, "tube": name, **population,
                    **{f"summary_{key}": value for key, value in tube.get("summary", {}).items()},
                })
    columns = list(dict.fromkeys(key for record in records for key in record))
    return columns, [[_cell_value(record.get(column)) for column in columns] for record in records]


def _write_xlsx(result: Mapping[str, Any], path: Path, headers: list[str], rows: list[list[Any]]) -> None:
    import openpyxl
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    workbook = openpyxl.Workbook()
    long_values: list[list[Any]] = []

    def append(sheet, values: Sequence[Any]) -> None:
        prepared = []
        for column, value in enumerate(values, start=1):
            value = _cell_value(value)
            if isinstance(value, str) and len(value) > 32767:
                # Excel silently truncates longer cells. Preserve every character
                # in numbered chunks instead, while the JSON also remains complete.
                address = f"{get_column_letter(column)}{sheet.max_row + 1}"
                parts = [value[start:start + 30000] for start in range(0, len(value), 30000)]
                long_values.extend([sheet.title, address, number, part] for number, part in enumerate(parts, start=1))
                value = tr('Voir Données longues · {v0}!{v1} · {v2} parties', v0=sheet.title, v1=address, v2=len(parts))
            prepared.append(value)
        sheet.append(prepared)

    candidates = workbook.active
    candidates.title = tr('Candidats étendus')
    append(candidates, headers)
    for row in rows:
        append(candidates, row)
    candidates.auto_filter.ref = candidates.dimensions
    tube_headers, tube_rows = _tube_composition_table(result)
    if tube_headers:
        sheet = workbook.create_sheet(tr('Composition des tubes'))
        append(sheet, tube_headers)
        for row in tube_rows:
            append(sheet, row)
        sheet.auto_filter.ref = sheet.dimensions
    for name, values in (
        (tr('Paramètres'), result.get("request", {})),
        (tr('Référence'), result.get("reference", {})),
        (tr('Recherche'), result.get("mfe_screen", {})),
    ):
        sheet = workbook.create_sheet(name)
        append(sheet, [tr('Champ'), tr('Valeur')])
        if isinstance(values, Mapping):
            for key, value in values.items():
                append(sheet, [key, value])
        else:
            append(sheet, [tr('Données'), values])
    curve_names = {
        "mfe_progress": tr('Progression MFE'),
        "ensemble_progress": tr('Progression ensemble'),
        "ranking": tr('Classement (courbe)'),
    }
    for key, records in result.get("curves", {}).items():
        sheet = workbook.create_sheet(curve_names.get(key, tr('Courbe')))
        if isinstance(records, list) and all(isinstance(record, Mapping) for record in records):
            columns = list(dict.fromkeys(field for record in records for field in record))
            if columns:
                append(sheet, columns)
                for record in records:
                    append(sheet, [record.get(column) for column in columns])
                sheet.auto_filter.ref = sheet.dimensions
            else:
                append(sheet, [tr('Aucune donnée')])
        else:
            append(sheet, [tr('Champ'), tr('Valeur')])
            append(sheet, [key, records])
    metadata = workbook.create_sheet(tr('Exécution'))
    append(metadata, [tr('Champ'), tr('Valeur')])
    for key in ("run_id", "status", "started_at", "completed_at", "warnings", "parasite_analysis", "thermodynamics"):
        append(metadata, [key, result.get(key)])
    if long_values:
        sheet = workbook.create_sheet(tr('Données longues'))
        sheet.append([tr('Feuille'), tr('Cellule'), tr('Partie'), tr('Valeur')])
        for row in long_values:
            sheet.append(row)

    for sheet in workbook:
        sheet.freeze_panes = "A2"
        for cell in sheet[1]:
            cell.fill = PatternFill("solid", fgColor="193B35")
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        sheet.row_dimensions[1].height = 32
        for index, column in enumerate(sheet.columns, start=1):
            length = max((len(str(cell.value or "")) for cell in column), default=12)
            sheet.column_dimensions[get_column_letter(index)].width = min(72, max(14, length + 2))
        for row in sheet.iter_rows(min_row=2):
            for cell in row:
                if isinstance(cell.value, str):
                    # Metadata is data, including arbitrary user-provided target names.
                    cell.data_type = "s"
                cell.alignment = Alignment(vertical="top", wrap_text=True)
    workbook.save(path)


def _write_data_exports(result: Mapping[str, Any]) -> dict[str, str]:
    output_dir = Path(result["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    run_id = result["run_id"]
    json_path = output_dir / f"{run_id}_results.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    csv_path = output_dir / f"{run_id}_ranked.csv"
    headers, rows = _candidate_table(result.get("candidates", []))
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)
    xlsx_path = output_dir / f"{run_id}_ranked.xlsx"
    _write_xlsx(result, xlsx_path, headers, rows)
    paths = {"json": str(json_path), "csv": str(csv_path), "xlsx": str(xlsx_path)}
    tube_headers, tube_rows = _tube_composition_table(result)
    if tube_headers:
        tube_path = output_dir / f"{run_id}_tube_complexes.csv"
        with tube_path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(tube_headers)
            writer.writerows(tube_rows)
        paths["tube_csv"] = str(tube_path)
    paths.update(export_extension_visuals(result, output_dir))
    return paths


def _finish_exports(result: Mapping[str, Any], paths: dict[str, str], log: Sequence[str]) -> dict[str, str]:
    output_dir = Path(result["output_dir"])
    log_path = output_dir / f"{result['run_id']}.log"
    engine = result_engine_label(dict(result))
    model = result.get("thermodynamics", {}).get("model")
    if model:
        engine += f" · {model}"
    lines = [
        tr('Identifiant du run : {v0}', v0=result['run_id']),
        tr('Statut : {v0}', v0=tr(result['status'])),
        tr('Début : {v0}', v0=result['started_at']),
        tr('Fin : {v0}', v0=result['completed_at']),
        tr('Moteur : {engine} · recherche exhaustive d’extension ({side})',
           engine=engine,
           side=tr(EXTENSION_SIDES[result.get('input', {}).get('extension_side', '5prime')])),
        tr('Candidats exportés : {v0}', v0=len(result.get('candidates', []))),
        "",
        *log,
    ]
    warnings = result.get("warnings", [])
    if warnings:
        lines.extend(["", tr('Avertissements :'), *(f"- {warning}" for warning in warnings)])
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    paths["log"] = str(log_path)
    archive_path = output_dir / f"{result['run_id']}_complete.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path_value in paths.values():
            path = Path(path_value)
            archive.write(path, arcname=path.relative_to(output_dir))
    paths["zip"] = str(archive_path)
    return paths


def _analyze_extension_tubes(
    result: dict[str, Any], request: ExtensionRunRequest, *,
    progress: ProgressCallback, cancel_check: Callable[[], bool] | None, backend=None,
) -> None:
    """Analyze the original and the best candidates after the unchanged ranking."""
    selected: list[dict[str, Any]] = []
    candidates = result.get("candidates", [])
    for length in sorted(set(request.target_lengths)):
        ranked = sorted(
            (candidate for candidate in candidates if candidate.get("trigger_length") == length),
            key=lambda candidate: (candidate.get("rank_within_length", candidate.get("rank", float("inf"))), candidate.get("rank", float("inf"))),
        )
        selected.extend(ranked[:request.parasite_top])
    metadata = {
        "status": "pending" if request.parasite_top else "disabled",
        "top_per_length": request.parasite_top,
        "conditions": asdict(request.tube_conditions),
        "selected_count": len(selected), "completed_count": 0,
        "failed_count": 0, "original_status": "not_analyzed",
        "selection": [{
            "rank": candidate.get("rank"), "rank_within_length": candidate.get("rank_within_length"),
            "trigger_length": candidate.get("trigger_length"),
        } for candidate in selected],
    }
    result["parasite_analysis"] = metadata
    if not request.parasite_top:
        return
    if result.get("status") == "cancelled" or (cancel_check and cancel_check()):
        metadata["status"] = "cancelled"
        result["status"] = "cancelled"
        return
    samples = [(result.setdefault("reference", {}), request.trigger_dna, tr('Trigger original'))]
    samples.extend((candidate, candidate["extended_trigger"], tr(
        'Candidat {rank} · {length} nt', rank=candidate.get("rank_within_length", candidate.get("rank")),
        length=candidate.get("trigger_length"),
    )) for candidate in selected)
    total = 2 * len(samples)
    progress(0, total, tr(
        'Complexes parasites · original et {count} candidats ; jusqu’à {size} brins par complexe.',
        count=len(selected), size=request.tube_conditions.max_complex_size,
    ))
    for index, (sample, sequence, label) in enumerate(samples):
        if cancel_check and cancel_check():
            metadata["status"] = "cancelled"
            result["status"] = "cancelled"
            break

        def tube_progress(done: int, expected: int, message: str) -> None:
            progress(index * 2 + (2 * done // max(1, expected)), total, f"{label} · {message}")

        progress(index * 2, total, tr('Analyse des tubes · {sample}', sample=label))
        try:
            from aptaswitch_core.thermodynamics.application_analysis import analyze_parasite_tubes as analyze_native_tubes
            native = backend if backend is not None else get_backend(request.thermodynamics_library_path)
            analysis = analyze_native_tubes(
                native, request.aptamer_dna, sequence, request.thermo, request.tube_conditions,
                standard_state=native.standard_state, cancel_check=cancel_check,
                progress_callback=tube_progress,
            )
        except Exception as exc:
            analysis = {
                "status": "failed", "error": str(exc), "conditions": asdict(request.tube_conditions),
                "aptamer_sequence": request.aptamer_dna, "trigger_sequence": sequence, "tubes": {},
            }
            warning = tr('Analyse des tubes impossible pour {sample} : {error}', sample=label, error=exc)
            result.setdefault("warnings", []).append(warning)
            progress((index + 1) * 2, total, warning)
        sample["tube_analysis"] = analysis
        if index == 0:
            metadata["original_status"] = analysis.get("status")
        elif analysis.get("status") == "completed":
            metadata["completed_count"] += 1
        if analysis.get("status") == "failed":
            metadata["failed_count"] += 1
        for tube_name, summary_key, column in (
            ("trigger_only", "trigger_parasite_pct", "trigger_only_parasite_pct"),
            ("aptamer_trigger", "trigger_parasite_pct", "mixed_trigger_parasite_pct"),
            ("aptamer_trigger", "trigger_target_pct", "mixed_trigger_target_pct"),
        ):
            value = analysis.get("tubes", {}).get(tube_name, {}).get("summary", {}).get(summary_key)
            if value is not None:
                sample[column] = value
        if analysis.get("status") == "cancelled" or (cancel_check and cancel_check()):
            metadata["status"] = "cancelled"
            result["status"] = "cancelled"
            break
    else:
        metadata["status"] = "completed_with_errors" if metadata["failed_count"] else "completed"
        progress(total, total, tr(
            'Complexes parasites · {count}/{selected} candidats analysés ; {failed} analyse(s) en échec.',
            count=metadata["completed_count"], selected=metadata["selected_count"], failed=metadata["failed_count"],
        ))


def run_extension_and_export(
    request: ExtensionRunRequest,
    *,
    progress_callback: ProgressCallback | None = None,
    cancel_check: Callable[[], bool] | None = None,
    advance_check: Callable[[], bool] | None = None,
    phase_callback: Callable[[str], None] | None = None,
    resume_result: dict | None = None,
    backend=None,
) -> tuple[dict[str, Any], dict[str, str]]:
    """Execute the scientific search and preserve its exact result and provenance."""
    backend = backend if backend is not None else get_backend(request.thermodynamics_library_path)
    settings = extension_model_settings(request.thermo, backend)
    if resume_result is not None and (
            resume_result.get("thermodynamics") != asdict(backend.provenance)
            or resume_result.get("model") != settings):
        raise ValueError("The extension checkpoint uses a different thermodynamic model or conditions")
    saved_request = request.to_dict()
    if backend is not None:
        # Persist the effective model so a later request can validate its identity.
        saved_request["thermo"]["material"] = settings["material"]
        saved_request["import_path"] = ""
        saved_request["engine"] = ENGINE_ID
    started = datetime.now().astimezone()
    run_id, output_dir = new_run_directory(request.output_dir, "extension", started)
    log: list[str] = []
    current = 0
    trigger_length = len(normalize_dna(request.trigger_dna, field_name="trigger"))
    selected = ligand_positions(request.ligand_aptamer_positions, len(normalize_dna(request.aptamer_dna, field_name="aptamère")))
    estimate = estimate_extension_search(trigger_length, request.target_lengths)
    total = estimate["total_candidates_expected"]

    def checkpoint(partial: dict) -> None:
        # Atomic replacement preserves the previous checkpoint on a crash.
        output_dir.mkdir(parents=True, exist_ok=True)
        data = {**partial, "run_id": run_id, "request": saved_request,
                "started_at": started.isoformat(timespec="seconds"),
                "output_dir": str(output_dir), "checkpoint_version": 1}
        path = output_dir / "mfe_checkpoint.json"
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
        temporary.replace(path)

    def progress(done: int, expected: int, message: str) -> None:
        nonlocal current, total
        current, total = done, expected
        log.append(_timestamped(message))
        if progress_callback:
            progress_callback(done, expected, message)

    if resume_result is None:
        progress(0, total, tr('Recherche d’extensions · {v0:,} extensions à examiner.', v0=total).replace(",", " "))
    else:
        progress(0, resume_result["mfe_screen"].get("mfe_survivors_for_ensemble", 0),
                 tr("Reprise du groupe MFE retenu · vérification de la référence thermodynamique…"))
    result = dict(design_trigger_extensions(
        aptamer_dna=request.aptamer_dna,
        trigger_dna=request.trigger_dna,
        target_lengths=request.target_lengths,
        thermo=request.thermo,
        import_path=request.import_path,
        top=request.top,
        max_ensemble_candidates=request.max_ensemble_candidates,
        extension_side=request.extension_side,
        progress_callback=progress,
        cancel_check=cancel_check,
        advance_check=advance_check,
        phase_callback=phase_callback,
        checkpoint_callback=checkpoint,
        resume_result=resume_result,
        backend=backend,
    ))
    result["input"]["ligand_aptamer_positions"] = list(selected)
    if phase_callback:
        phase_callback("tubes")
    _analyze_extension_tubes(result, request, progress=progress, cancel_check=cancel_check, backend=backend)
    result.update({
        "run_id": run_id,
        "started_at": started.isoformat(timespec="seconds"),
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "request": saved_request,
        "output_dir": str(output_dir),
    })
    progress(current, total, tr('Écriture des séquences, métriques, structures 2D et courbes…'))
    if phase_callback:
        phase_callback("exports")
    exports = _write_data_exports(result)
    message = (
        tr('Extension annulée — résultats partiels exportés.')
        if result["status"] == "cancelled"
        else tr('Extension terminée — résultats exportés.')
    )
    progress(current, total, message)
    exports = _finish_exports(result, exports, log)
    return result, exports


class ExtensionRunJob:
    """A cancellable extension run with the same UI lifecycle as DesignRunJob."""

    def __init__(self, request: ExtensionRunRequest, *, resume_result: dict | None = None, backend=None):
        extension_split(0, request.extension_side)
        self.request = request
        self._backend = backend if backend is not None else get_backend(request.thermodynamics_library_path)
        self._language = get_language()
        self._lock = threading.Lock()
        self._cancel = threading.Event()
        self._advance = threading.Event()
        self._phase = "reference"
        self._resume_result = resume_result
        self._state = "running"
        self._current = 0
        trigger_length = len(normalize_dna(request.trigger_dna, field_name="trigger"))
        self._total = estimate_extension_search(trigger_length, request.target_lengths)["total_candidates_expected"]
        self._message = tr('Préparation de l’extension…')
        self._started_monotonic = time.monotonic()
        self._phase_started_monotonic = self._started_monotonic
        self._phase_start_current = 0
        self._finished_monotonic: float | None = None
        self._log = [_timestamped(self._message)]
        self._result: dict[str, Any] | None = None
        self._exports: dict[str, str] | None = None
        self._error: str | None = None
        self._consumed = False
        # Preserve the language selected when this job was created.
        context = copy_context()
        self._thread = threading.Thread(target=context.run, args=(self._run,), daemon=True)

    def start(self) -> None:
        self._thread.start()

    def cancel(self) -> None:
        self._cancel.set()
        with self._lock, language_context(self._language):
            if self._state == "running":
                self._message = tr('Annulation demandée — fin du calcul en cours…')
                self._append_log_locked(self._message)

    def advance_to_ensemble(self) -> None:
        with self._lock, language_context(self._language):
            if self._state == "running" and self._phase == "mfe" and not self._cancel.is_set():
                self._advance.set()
                self._message = tr("Passage à l’analyse des candidats retenus après le calcul en cours…")
                self._append_log_locked(self._message)

    @property
    def can_advance(self) -> bool:
        with self._lock:
            return (self._state == "running" and self._phase == "mfe"
                    and self._current > 0 and not self._advance.is_set() and not self._cancel.is_set())

    def _set_phase(self, phase: str) -> None:
        with self._lock:
            self._phase = phase
            self._phase_started_monotonic = time.monotonic()
            self._phase_start_current = 0

    def _append_log_locked(self, message: str) -> None:
        self._log.append(_timestamped(message))
        # The complete log is retained in the exported journal.
        if len(self._log) > 500:
            self._log = self._log[-500:]

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
            result, exports = run_extension_and_export(
                self.request,
                progress_callback=self._progress,
                cancel_check=self._cancel.is_set,
                advance_check=self._advance.is_set,
                phase_callback=self._set_phase,
                resume_result=self._resume_result,
                backend=self._backend,
            )
        except Exception as exc:
            with self._lock:
                self._state = "failed"
                self._error = str(exc)
                self._message = tr('L’extension a échoué.')
                self._finished_monotonic = time.monotonic()
                self._append_log_locked(tr('ERREUR · {v0}', v0=exc))
            return
        with self._lock:
            self._state = "done"
            self._result = result
            self._exports = exports
            if result["status"] == "completed":
                self._current = self._total
            self._message = (
                tr('Extension annulée — résultats partiels exportés.')
                if result["status"] == "cancelled"
                else tr('Extension terminée.')
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
            )

    def consume_once(self) -> bool:
        with self._lock:
            if self._consumed or self._state not in {"done", "failed"}:
                return False
            self._consumed = True
            return True
