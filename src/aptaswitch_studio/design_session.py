"""Portable scientific state for a named design, independent of Streamlit.

Only explicit form fields are stored. Worker objects, user settings and local
dependency locations never become project settings. Calculations remain the
source of truth for exported results and are reloaded relative to the project.
"""

from __future__ import annotations

import copy
import json
import math
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, MutableMapping

from aptaswitch_core.localization import tr
from aptaswitch_core.models import ReporterContext, ThermoConditions
from aptaswitch_core.thermodynamics.setup import configured_library_path
from aptaswitch_core.reporter_library import all_reporters
from aptaswitch_core.sequences import normalize_dna
from aptaswitch_core.tube_analysis import TubeAnalysisConditions

from .run_library import discover_runs, load_run
from .web_runtime import JobSnapshot, concentration_from_m


FORM_FIELDS = frozenset({
    "design_mode", "reading_switch_input", "reading_use_full_reporter",
    "molecule_name", "aptamer_input", "trigger_input", "t7_leader", "rbs_loop", "linker_auto",
    "frame_linker", "trials", "optimization_max_evaluations",
    "optimization_population_size", "optimization_lineages", "optimization_local_tolerance",
    "exclude_stop_candidates", "exclude_non_linear_candidates", "temperature_c",
    "trigger_concentration_value", "trigger_concentration_unit",
    "switch_concentration_value", "switch_concentration_unit",
    "sodium_concentration_value", "sodium_concentration_unit",
    "magnesium_concentration_value", "magnesium_concentration_unit",
    "reporter_key", "last_reporter_key", "weight_on", "weight_leak", "weight_defect",
    "weight_functional_off_toehold", "weight_functional_off_stem",
    "weight_functional_on_rbs_aug", "weight_functional_on_rbs_linker",
    "weight_ddg", "weight_rbs", "preview_state", "results_kind", "selected_trial",
    "weight_trigger_linearity", "weight_rbs_calculator",
    "candidate_state", "extension_enabled", "extension_side", "extension_mode", "extension_target_length",
    "extension_allow_long", "extension_top", "extension_ensemble_cap", "extension_parasite_top",
    "extension_temperature", "extension_sodium", "extension_sodium_unit",
    "extension_magnesium", "extension_magnesium_unit",
    "extension_ligand_positions", "extension_ligand_snapshot", "aptamer_indices",
    "system_aptamer_concentration", "system_aptamer_concentration_unit",
    "system_trigger_concentration", "system_trigger_concentration_unit", "system_tube_volume_ul",
    "inverse_extension_enabled", "inverse_extension_target_length", "inverse_extension_side",
    "inverse_extension_top", "inverse_extension_performance_weight",
})
# Old project forms remain importable, but their geometry controls are no
# longer editable inputs. Clear any such keys left in an existing session;
# recorded result/checkpoint architectures are independent and stay intact.
_RETIRED_FORM_FIELDS = frozenset({
    "customize_tswitch", "toehold_length", "rbs_aug_distance", "binding_start", "binding_end",
})
JOB_KEYS = ("design_job", "extension_job", "system_tube_job", "blast_job", "ensemble_job", "saved_switch_analysis_job", "saved_switch_ffs_job", "saved_switch_free_md_job")
_RESET_FIELDS = FORM_FIELDS | _RETIRED_FORM_FIELDS | frozenset({
    "active_design_path", "active_design_name", "output_dir", "trigger_snapshot",
    "last_result", "last_exports", "run_error", "run_notice", "sequence_action_error",
    "long_run_confirm", "saved_run_path", "run_open_path", "run_open_error", "run_library_extra_roots", "switch_results_section",
    "opened_design_source", "opened_extension_source", "design_original_aptamer",
    "design_original_trigger", "design_using_extension", "design_reporter", "design_save_error",
    "design_restore_warnings", "_design_saved_state", "_design_captured_tubes",
    "_design_captured_analysis", "reporter_return_to_design", "design_section",
    *JOB_KEYS,
})
_RESET_PREFIXES = ("extension_", "system_", "initial_", "switch_preview",
                   "switch_candidate", "system_tubes", "candidate_export", "blast_", "inverse_extension_", "ensemble_", "context_ensemble_", "saved_switch_", "results_table_", "results_candidate_")


def _json_copy(value: Any) -> Any:
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def validate_design_state(state: Any) -> None:
    """Reject malformed editable state before an import replaces another design."""
    def invalid(field: str) -> None:
        raise ValueError(tr("Les paramètres enregistrés du design sont invalides : {field}.", field=field))

    if not isinstance(state, dict):
        invalid("state")
    for key in ("form", "original_sequences", "selected_runs", "initial_analysis", "reporter", "blast"):
        if key in state and not isinstance(state[key], dict):
            invalid(key)
    if "using_extension" in state and not isinstance(state["using_extension"], bool):
        invalid("using_extension")
    blast = state.get("blast", {})
    if "taxa" in blast:
        taxa = blast["taxa"]
        if (not isinstance(taxa, list) or len(taxa) > 20
                or any(not isinstance(taxon, dict) or type(taxon.get("taxid")) is not int
                       or not 1 <= taxon["taxid"] <= 2 ** 31 - 1
                       or not isinstance(taxon.get("name"), str)
                       or not 1 <= len(taxon["name"].strip()) <= 300 for taxon in taxa)):
            invalid("blast.taxa")
        if len({taxon["taxid"] for taxon in taxa}) != len(taxa):
            invalid("blast.taxa")
    if "scopes" in blast and (not isinstance(blast["scopes"], list)
                              or any(not isinstance(scope, str) or scope not in {"genomic", "transcript"} for scope in blast["scopes"])
                              or len(set(blast["scopes"])) != len(blast["scopes"])):
        invalid("blast.scopes")
    if "max_hits" in blast and (type(blast["max_hits"]) is not int or not 1 <= blast["max_hits"] <= 500):
        invalid("blast.max_hits")
    bool_fields = {"linker_auto", "exclude_stop_candidates",
                   "exclude_non_linear_candidates", "extension_allow_long", "inverse_extension_enabled", "extension_enabled",
                   "reading_use_full_reporter"}
    integer_fields = {"trials", "optimization_max_evaluations",
                      "optimization_population_size", "optimization_lineages",
                      "selected_trial", "extension_target_length", "extension_top", "extension_ensemble_cap",
                      "extension_parasite_top",
                      "inverse_extension_target_length", "inverse_extension_top"}
    list_fields = {"extension_ligand_positions", "aptamer_indices"}
    numeric_fields = {key for key in FORM_FIELDS if key.startswith("weight_") or key.endswith("_concentration_value")}
    numeric_fields.update({"temperature_c", "extension_temperature", "extension_sodium", "extension_magnesium",
                           "system_aptamer_concentration", "system_trigger_concentration", "system_tube_volume_ul",
                           "inverse_extension_performance_weight", "optimization_local_tolerance"})
    for key, value in state.get("form", {}).items():
        if key in bool_fields:
            valid = isinstance(value, bool)
        elif key in integer_fields or key.startswith("extension_selected_"):
            minimum = 0 if key == "extension_parasite_top" or key.startswith("extension_selected_") else 1
            valid = type(value) is int and minimum <= value <= 2 ** 53 - 1
        elif key in numeric_fields:
            valid = type(value) in (int, float) and math.isfinite(value)
        elif key in list_fields:
            valid = isinstance(value, list) and all(type(item) is int and item >= 0 for item in value)
        elif key in FORM_FIELDS:
            valid = isinstance(value, str)
        else:
            continue
        if not valid:
            invalid(key)
    choices = {
        "design_mode": {"design", "reading"},
        "extension_side": {"5prime", "3prime", "both"},
        "inverse_extension_side": {"5prime", "3prime", "both"},
        "extension_mode": {"Deux versions : 24 et 30 nt", "Longueur personnalisée"},
        "results_kind": {"Design de switches", "Extension de triggers"},
        "preview_state": {"OFF", "ON + trigger"}, "candidate_state": {"OFF", "ON + trigger"},
    }
    form = state.get("form", {})
    for key, allowed in choices.items():
        if key in form and form[key] not in allowed:
            invalid(key)
    for key, value in form.items():
        if key.endswith("_unit") and key in FORM_FIELDS and value not in {"nM", "µM", "mM", "M"}:
            invalid(key)
    # Mirror bounded form controls before a staged import can replace a project.
    # Drafts may contain empty or partial sequence inputs. Retired geometry
    # controls in historical manifests are ignored rather than restored.
    bounds = {
        "trials": (1, 100000), "optimization_max_evaluations": (1, 100000), "temperature_c": (0, 100),
        "optimization_population_size": (1, 512), "optimization_lineages": (1, 64), "optimization_local_tolerance": (0, 1),
        "extension_temperature": (0, 100), "extension_top": (1, 1000),
        "extension_parasite_top": (0, 1000),
        "inverse_extension_target_length": (2, 1000), "inverse_extension_top": (1, 1000),
        "inverse_extension_performance_weight": (0, 100),
        "system_tube_volume_ul": (0.001, None),
        "extension_target_length": (2, None if form.get("extension_allow_long", False) else 30),
        **{key: (0, 20) for key in numeric_fields if key.startswith("weight_")},
    }
    for key, (minimum, maximum) in bounds.items():
        if key in form and (form[key] < minimum or (maximum is not None and form[key] > maximum)):
            invalid(key)
    if form.get("optimization_lineages", 1) > form.get("optimization_population_size", 512):
        invalid("optimization_lineages")
    # Concentration limits belong to a physical quantity, regardless of the
    # chosen display unit. Defaults match initialization when a field is absent.
    concentrations = (
        ("trigger_concentration_value", "trigger_concentration_unit", 5., "µM", 1e-9, 1e-2),
        ("switch_concentration_value", "switch_concentration_unit", 5., "µM", 1e-9, 1e-2),
        ("sodium_concentration_value", "sodium_concentration_unit", 120., "mM", .12, 1.),
        ("magnesium_concentration_value", "magnesium_concentration_unit", 0., "mM", 0., 0.),
        ("system_aptamer_concentration", "system_aptamer_concentration_unit", 5., "µM", 1e-12, None),
        ("system_trigger_concentration", "system_trigger_concentration_unit", 5., "µM", 1e-12, None),
        ("extension_sodium", "extension_sodium_unit", 120., "mM", .12, 1.),
        ("extension_magnesium", "extension_magnesium_unit", 0., "mM", 0., 0.),
    )
    for value_key, unit_key, default_value, default_unit, minimum, maximum in concentrations:
        if value_key not in form and unit_key not in form:
            continue
        value, unit = form.get(value_key, default_value), form.get(unit_key, default_unit)
        if (value < concentration_from_m(minimum, unit)
                or (maximum is not None and value > concentration_from_m(maximum, unit))):
            invalid(value_key)
    for key, value in state.get("original_sequences", {}).items():
        if key in {"aptamer", "trigger"} and not isinstance(value, str):
            invalid(key)
    for key, value in state.get("selected_runs", {}).items():
        if not isinstance(value, str) or Path(value).is_absolute() or ".." in Path(value).parts:
            invalid(key)
    analysis = state.get("initial_analysis", {})
    for key in ("structure", "tubes"):
        if key in analysis and not isinstance(analysis[key], dict):
            invalid(f"initial_analysis.{key}")
    if "structure" in analysis:
        signature = analysis["structure"].get("signature")
        if (not isinstance(signature, list) or len(signature) not in {7, 8}
                or not all(isinstance(item, str) for item in signature[:4])
                or not all(type(item) in (int, float) and math.isfinite(item) for item in signature[4:7])
                or (len(signature) == 8 and signature[7] is not None and not isinstance(signature[7], str))
                or not isinstance(analysis["structure"].get("result"), dict)):
            invalid("initial_analysis.structure")
    if "tubes" in analysis:
        tube = analysis["tubes"]
        if (not isinstance(tube.get("aptamer"), str) or not isinstance(tube.get("trigger"), str)
                or not all(isinstance(tube.get(key), dict) for key in ("thermo", "conditions", "result"))
                or not isinstance(tube.get("log", []), list)):
            invalid("initial_analysis.tubes")
    reporter = state.get("reporter", form.get("design_reporter"))
    if reporter is not None:
        required = {"key", "label", "sequence_after_start_rna"}
        optional = {"full_sequence_after_start_rna", "source_accession", "source_url"}
        if (not isinstance(reporter, dict) or not required.issubset(reporter)
                or set(reporter) - required - optional
                or any(not isinstance(reporter[name], str) for name in required)
                or any(value is not None and not isinstance(value, str) for name, value in reporter.items() if name in optional)
                or not reporter["key"]):
            invalid("reporter")
        try:
            ReporterContext(**reporter).validate()
        except (TypeError, ValueError):
            invalid("reporter")


def attach_saved_switch_jobs(session: Mapping[str, Any], *, ignore_ffs: bool = False) -> bool:
    """Recover current-design workers lost by a browser refresh.

    Registries outlive a browser session. Read-only mappings can still query
    their busy state; mutable sessions also regain navigation/progress handles.
    """
    active = session.get("active_design_path")
    if not active:
        return False
    from .saved_switch_runtime import active_saved_analysis_for_design
    from .saved_switch_ffs_runtime import active_ffs_for_design
    from .saved_switch_free_md_runtime import active_free_md_for_design
    from .saved_switch_ensemble_runtime import active_context_ensemble_for_design

    jobs = (
        ("saved_switch_analysis_job", active_saved_analysis_for_design(active)),
        ("saved_switch_ffs_job", active_ffs_for_design(active)),
        ("saved_switch_free_md_job", active_free_md_for_design(active)),
        ("ensemble_job", active_context_ensemble_for_design(active)),
    )
    if isinstance(session, MutableMapping):
        for key, job in jobs:
            if job is not None:
                session[key] = job
                if key == "ensemble_job":
                    session["ensemble_job_key"] = job.key
    return any(job is not None for key, job in jobs
               if not (ignore_ffs and key == "saved_switch_ffs_job"))


def active_job_running(session: Mapping[str, Any], *, ignore_ffs: bool = False) -> bool:
    """Keep exclusive actions guarded; context thermodynamics may share with FFS.

    FFS owns a detached worker and immutable inputs. Callers opting in still
    recover its monitor, while every other running job retains its guard.
    """
    recovered = attach_saved_switch_jobs(session, ignore_ffs=ignore_ffs)
    from .web_runtime import find_design_job
    root = (Path(session["active_design_path"]) / "calculations"
            if session.get("active_design_path") else session.get("output_dir"))
    job = find_design_job(root) if root else None
    design_running = job is not None and job.snapshot().state == "running"
    if design_running and isinstance(session, MutableMapping):
        session["design_job"] = job
    recovered = recovered or design_running
    return recovered or any(session.get(key) is not None and session[key].snapshot().state == "running"
                            for key in JOB_KEYS if not (ignore_ffs and key == "saved_switch_ffs_job"))


def reset_design_session(session: MutableMapping[str, Any]) -> None:
    """Discard project state only, after making sure no worker can write into it."""
    if active_job_running(session):
        raise ValueError(tr("Attendez la fin du calcul en cours avant de changer de design."))
    for key in list(session):
        if key in _RESET_FIELDS or key.startswith(_RESET_PREFIXES):
            del session[key]


def remember_original_sequences(session: MutableMapping[str, Any]) -> None:
    """Follow edits to original inputs until an extension is deliberately applied."""
    if session.get("design_using_extension"):
        return
    for target, field in (("aptamer", "aptamer_input"), ("trigger", "trigger_input")):
        if field in session:
            session[f"design_original_{target}"] = str(session[field])


def project_reporters(session: Mapping[str, Any]) -> dict[str, ReporterContext]:
    """Carry the chosen reporter with a project without modifying user settings."""
    reporters = dict(all_reporters())
    saved = session.get("design_reporter")
    if isinstance(saved, dict):
        try:
            reporter = ReporterContext(**saved)
            if (isinstance(reporter.key, str) and reporter.key
                    and isinstance(reporter.label, str)
                    and isinstance(reporter.sequence_after_start_rna, str)
                    and not set(reporter.sequence_after_start_rna.upper()) - set("ACGTU")):
                reporters[reporter.key] = reporter
        except TypeError:
            pass
    return reporters


def _relative_file(path: Any, root: Path) -> str | None:
    if not path:
        return None
    try:
        candidate = Path(path).resolve()
        relative = candidate.relative_to(root.resolve())
        return relative.as_posix() if candidate.is_file() else None
    except (OSError, TypeError, ValueError):
        return None


def _finished_results(session: MutableMapping[str, Any]) -> None:
    """A navigation callback can run before a finished worker's fragment reruns."""
    for key, result_key, exports_key in (
        ("design_job", "last_result", "last_exports"),
        ("extension_job", "extension_result", "extension_exports"),
        ("blast_job", "blast_result", "blast_exports"),
    ):
        job = session.get(key)
        if job is not None:
            snapshot = job.snapshot()
            if (snapshot.state in {"done", "paused"} or key == "blast_job" and snapshot.state == "failed") and snapshot.result is not None:
                session[result_key], session[exports_key] = snapshot.result, snapshot.exports or {}


def _capture_initial_analysis(session: Mapping[str, Any]) -> dict[str, Any]:
    captured = {}
    signature = session.get("system_analysis_signature")
    if signature and session.get("system_analysis_result") is not None:
        captured["structure"] = {
            "signature": [signature[0], signature[1], *signature[3:]],
            "result": session["system_analysis_result"],
        }
    job = session.get("system_tube_job")
    if job is not None:
        snapshot = job.snapshot()
        if snapshot.state != "running" and snapshot.result is not None:
            aptamer, trigger, _import_path, thermo, conditions = job.signature[:5]
            captured["tubes"] = {
                "aptamer": aptamer, "trigger": trigger, "thermo": asdict(thermo),
                "conditions": asdict(conditions), "result": snapshot.result,
                "backend_identity": job.signature[5] if len(job.signature) > 5 else None,
                "error": snapshot.error, "log": list(snapshot.log),
            }
    return _json_copy(captured)


def save_active_design(session: MutableMapping[str, Any]):
    """Save forms and portable references; return None when no design is open."""
    from .design_library import load_design, save_design

    active = session.get("active_design_path")
    if not active:
        return None
    root = Path(active).expanduser().resolve()
    for key in _RETIRED_FORM_FIELDS:
        session.pop(key, None)
    _finished_results(session)
    remember_original_sequences(session)
    form = {key: _json_copy(session[key]) for key in FORM_FIELDS if key in session}
    for key in session:
        if key.startswith("extension_selected_"):
            form[key] = _json_copy(session[key])
    chosen = project_reporters(session).get(session.get("reporter_key"))
    if chosen is not None:
        session["design_reporter"] = asdict(chosen)
    selected = {}
    for kind, exports_key in (("design", "last_exports"), ("extension", "extension_exports"), ("blast", "blast_exports")):
        exports = session.get(exports_key) or {}
        source = exports.get("checkpoint") if kind == "design" else None
        path = _relative_file(source or exports.get("json"), root)
        if path:
            selected[kind] = path
    state = {
        "form": form,
        "original_sequences": {target: session.get(f"design_original_{target}", "")
                               for target in ("aptamer", "trigger")},
        "using_extension": bool(session.get("design_using_extension")),
        "selected_runs": selected,
        "initial_analysis": _capture_initial_analysis(session),
    }
    if session.get("design_reporter"):
        state["reporter"] = _json_copy(session["design_reporter"])
    blast = {key: _json_copy(session[f"blast_{key}"]) for key in ("taxa", "scopes", "max_hits")
             if f"blast_{key}" in session}
    if blast:
        state["blast"] = blast
    # Avoid rewriting the manifest on every timed redraw or unchanged rerun.
    if session.get("_design_saved_state") == state:
        return load_design(root)
    project = save_design(root, state)
    _write_sequence_file(root, "original.json", state["original_sequences"])
    _write_sequence_file(root, "current.json", {
        "aptamer": form.get("aptamer_input", ""), "trigger": form.get("trigger_input", ""),
        "using_extension": state["using_extension"],
    })
    session["_design_saved_state"] = copy.deepcopy(state)
    session.pop("design_save_error", None)
    return project


def _write_sequence_file(root: Path, name: str, value: dict) -> None:
    folder = root / "sequences"
    if folder.is_symlink():
        raise ValueError(tr("Le dossier de séquences ne peut pas être un lien symbolique."))
    folder.mkdir(exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".sequences-", suffix=".json", dir=folder)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        os.replace(temporary, folder / name)
    finally:
        Path(temporary).unlink(missing_ok=True)


def save_active_design_safely(session: MutableMapping[str, Any]) -> None:
    """Fragment completion should still display results if saving is unavailable."""
    try:
        save_active_design(session)
    except (OSError, ValueError) as exc:
        session["design_save_error"] = str(exc)


def _safe_project_file(root: Path, relative: Any) -> Path | None:
    if not isinstance(relative, str) or Path(relative).is_absolute():
        return None
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError:
        return None
    return path if path.is_file() else None


def _restore_runs(root: Path, saved: Mapping[str, Any], session: MutableMapping[str, Any]) -> None:
    runs = discover_runs([root / "calculations"])
    for kind, result_key, exports_key in (
        ("design", "last_result", "last_exports"),
        ("extension", "extension_result", "extension_exports"),
    ):
        selected = _safe_project_file(root, saved.get(kind))
        candidates = ([selected] if selected else []) + [run.path for run in runs if run.kind == kind]
        for path in candidates:
            try:
                # A resumed design can have a newer checkpoint alongside the
                # exports saved at its previous pause.
                latest = path.parent if kind == "design" and (path.parent / "design_checkpoint.json").is_file() else path
                loaded = load_run(latest)
            except (OSError, ValueError):
                continue
            if loaded.kind != kind:
                continue
            session[result_key], session[exports_key] = loaded.result, loaded.exports
            session[f"opened_{kind}_source"] = str(loaded.source)
            break




def _restore_blast(root: Path, saved: Mapping[str, Any], session: MutableMapping[str, Any]) -> None:
    from datetime import datetime

    selected = _safe_project_file(root, saved.get("blast"))
    sources = sorted((root / "calculations" / "blast").glob("*/results.json"), reverse=True)
    if not selected and not sources:
        return
    from .blast_runtime import load_blast_result

    newest = None
    for source in ([selected] if selected else []) + sources:
        if not _relative_file(source, root):
            continue
        try:
            result, exports = load_blast_result(source)
        except (OSError, ValueError):
            continue
        timestamp = datetime.fromisoformat(result["completed_at"]).timestamp()
        if selected == source or newest is None or timestamp > newest[0]:
            newest = timestamp, result, exports
        if selected == source:
            break
    if newest:
        _, session["blast_result"], session["blast_exports"] = newest
        # Restored RIDs are provenance only; no job or consent is re-created.




class _StoredTubeJob:
    """Read-only completed snapshot accepted by the existing tube display."""
    restored = True

    def __init__(self, signature: tuple, snapshot: JobSnapshot):
        self.signature, self._snapshot = signature, snapshot

    def snapshot(self) -> JobSnapshot:
        return self._snapshot

    def consume_once(self) -> bool:
        return False


@dataclass(frozen=True)
class _StoredTubeConditions:
    """Historical signature only; new calculations use TubeAnalysisConditions."""
    aptamer_concentration_m: float
    trigger_concentration_m: float
    volume_l: float
    max_complex_size: int


def _restore_initial_analysis(saved: Mapping[str, Any], session: MutableMapping[str, Any]) -> None:
    import_path = configured_library_path()
    structure = saved.get("structure")
    if isinstance(structure, dict) and len(structure.get("signature", ())) in {7, 8}:
        signature = structure["signature"]
        session["system_analysis_signature"] = (*signature[:2], import_path, *signature[2:7], signature[7] if len(signature) == 8 else None)
        session["system_analysis_result"] = structure.get("result")
        session["system_analysis_error"] = ""
        # A portable result belongs to the installation that produced it. Its
        # original local path is deliberately not saved; replacing that path
        # with the current one must not certify a cache hit on this engine.
        # Keep the snapshot available offline, and refresh before reuse when
        # the native runtime is available (even if all physical inputs still match).
        session["system_analysis_requires_refresh"] = True
    tubes = saved.get("tubes")
    if isinstance(tubes, dict):
        try:
            values = dict(tubes["conditions"])
            maximum = values.get("max_complex_size", 3)
            if type(maximum) is int and 3 < maximum <= 5:
                # Validate physical values with the current input contract,
                # then retain the original scope only in a read-only signature.
                checked = TubeAnalysisConditions(**{**values, "max_complex_size": 3})
                conditions = _StoredTubeConditions(**{**asdict(checked), "max_complex_size": maximum})
            else:
                conditions = TubeAnalysisConditions(**values)
            signature = (tubes["aptamer"], tubes["trigger"], import_path,
                         ThermoConditions(**tubes["thermo"]), conditions, tubes.get("backend_identity"))
            snapshot = JobSnapshot(state="done", current=2, total=2, message="",
                                   result=tubes["result"], error=tubes.get("error"), log=tuple(tubes.get("log", [])))
            session["system_tube_signature"] = signature
            session["system_tube_job"] = _StoredTubeJob(signature, snapshot)
        except (TypeError, KeyError, ValueError):
            pass


def open_design_into_session(path: str | Path, session: MutableMapping[str, Any]):
    """Validate first, then switch the entire scientific state as one project."""
    from .design_library import load_design

    if active_job_running(session):
        raise ValueError(tr("Attendez la fin du calcul en cours avant de changer de design."))
    project = load_design(Path(path))
    validate_design_state(project.state)
    if session.get("active_design_path") and Path(session["active_design_path"]).resolve() == project.path.resolve():
        return save_active_design(session)
    save_active_design(session)
    reset_design_session(session)
    saved = project.state
    form = saved.get("form", {})
    for key, value in form.items():
        if key in FORM_FIELDS or key.startswith("extension_selected_"):
            session[key] = copy.deepcopy(value)
    try:
        session["trigger_snapshot"] = normalize_dna(session.get("trigger_input", ""), field_name="trigger")
    except ValueError:
        session["trigger_snapshot"] = None
    originals = saved.get("original_sequences", {})
    for target in ("aptamer", "trigger"):
        session[f"design_original_{target}"] = originals.get(target, session.get(f"{target}_input", ""))
    session["design_using_extension"] = bool(saved.get("using_extension"))
    reporter = saved.get("reporter", form.get("design_reporter"))
    if reporter:
        session["design_reporter"] = copy.deepcopy(reporter)
    session["active_design_path"], session["active_design_name"] = str(project.path), project.name
    session["output_dir"] = str(project.path / "calculations")
    _restore_runs(project.path, saved.get("selected_runs", {}), session)
    _restore_initial_analysis(saved.get("initial_analysis", {}), session)
    for key in ("taxa", "scopes", "max_hits"):
        if key in saved.get("blast", {}):
            session[f"blast_{key}"] = copy.deepcopy(saved["blast"][key])
    _restore_blast(project.path, saved.get("selected_runs", {}), session)
    session["nav"] = "Conception"
    return project
