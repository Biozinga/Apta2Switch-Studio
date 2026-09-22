"""Open saved calculations in the same result views as a newly finished run."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import streamlit as st

from aptaswitch_core.localization import set_language, tr

from aptaswitch_core.trigger_extension import EXTENSION_SIDES

from .run_library import discover_runs, load_run
from .desktop import external_process_environment
from .web_runtime import default_export_dir, legacy_export_dir, normalized_export_dir


RESULT_KINDS = {"design": "Design de switches", "extension": "Extension de triggers"}


def _running_job() -> bool:
    return any(
        job is not None and job.snapshot().state == "running"
        for job in (st.session_state.get("design_job"), st.session_state.get("extension_job"))
    )


def _open_run(path: str) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.run_open_error = ""
    if _running_job():
        st.session_state.run_open_error = tr("Attendez la fin du calcul en cours avant d’ouvrir un autre run.")
        return
    if not path.strip():
        st.session_state.run_open_error = tr("Choisissez un run ou indiquez son chemin.")
        return
    try:
        loaded = load_run(Path(path.strip()).expanduser())
    except (OSError, ValueError) as exc:
        st.session_state.run_open_error = str(exc)
        return

    # Only replace a finished job after the complete file has been validated.
    # Otherwise its monitor could put the previous result back on screen.
    if loaded.kind == "design":
        st.session_state.last_result = loaded.result
        st.session_state.last_exports = loaded.exports
        st.session_state.design_job = None
        st.session_state.run_error = ""
        st.session_state.run_notice = ""
        for key in ("selected_trial", "candidate_state"):
            st.session_state.pop(key, None)
    else:
        st.session_state.extension_result = loaded.result
        st.session_state.extension_exports = loaded.exports
        st.session_state.extension_job = None
        st.session_state.extension_error = ""
        st.session_state.extension_notice = ""
        for key in list(st.session_state):
            if key.startswith("extension_selected_"):
                del st.session_state[key]
    st.session_state.results_kind = RESULT_KINDS[loaded.kind]
    st.session_state[f"opened_{loaded.kind}_source"] = str(loaded.source)
    # Remember external locations during this session without moving any files.
    roots = st.session_state.get("run_library_extra_roots", [])
    st.session_state.run_library_extra_roots = list(dict.fromkeys([*roots, str(loaded.source.parent)]))
    st.session_state.saved_run_path = str(loaded.source)


def _open_selected() -> None:
    _open_run(st.session_state.get("saved_run_path") or "")


def _open_path() -> None:
    _open_run(st.session_state.get("run_open_path", ""))


def _show_folder(path: Path) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    try:
        path.mkdir(parents=True, exist_ok=True)
        if sys.platform == "darwin":
            subprocess.run(["open", str(path)], check=True, capture_output=True,
                           env=external_process_environment())
        elif sys.platform == "win32":
            os.startfile(str(path))
        else:
            subprocess.run(["xdg-open", str(path)], check=True, capture_output=True,
                           env=external_process_environment())
    except (OSError, subprocess.SubprocessError) as exc:
        st.session_state.run_open_error = tr("Impossible d’ouvrir le dossier {path} : {error}", path=path, error=exc)


def render_run_library() -> None:
    root = normalized_export_dir(st.session_state.get("output_dir", ""))
    roots = [root, default_export_dir(), legacy_export_dir()]
    roots.extend(Path(path) for path in st.session_state.get("run_library_extra_roots", []))
    runs = discover_runs(roots)
    by_path = {str(run.path): run for run in runs}
    if st.session_state.get("saved_run_path") not in by_path:
        st.session_state.saved_run_path = next(iter(by_path), None)

    def label(path: str) -> str:
        run = by_path[path]
        date = run.started_at.replace("T", " ")[:19] if run.started_at else run.run_id
        return tr("{date} · {kind} · {name} · {count} candidat(s) · {run_id}", date=date, kind=tr(RESULT_KINDS[run.kind]), name=run.molecule_name, count=run.candidate_count, run_id=run.run_id)

    run_labels = {path: label(path) for path in by_path}
    empty = not st.session_state.get("last_result") and not st.session_state.get("extension_result")
    with st.expander(tr("Ouvrir un run enregistré"), expanded=empty):
        st.caption(tr("Dossier des runs : {path}", path=root))
        st.caption(tr("Les nouveaux runs sont rangés dans designs/ ou extensions/. Les anciens dossiers sont aussi recherchés."))
        st.selectbox(
            tr("Run à ouvrir"), options=list(by_path), format_func=run_labels.__getitem__,
            key="saved_run_path", placeholder=tr("Aucun run enregistré trouvé"), disabled=not runs,
        )
        controls = st.columns(3)
        controls[0].button(
            tr("Ouvrir ce run"), key="open_saved_run", on_click=_open_selected,
            disabled=not runs or _running_job(), width="stretch", type="primary",
        )
        controls[1].button(tr("Actualiser la liste"), key="refresh_run_library", width="stretch")
        controls[2].button(tr("Afficher le dossier"), on_click=_show_folder, args=(root,), width="stretch")
        if st.session_state.get("saved_run_path"):
            st.caption(tr("Fichier sélectionné : {path}", path=st.session_state.saved_run_path))
        with st.popover(tr("Ouvrir depuis un chemin")):
            st.text_input(
                tr("Dossier du run ou fichier de résultats"), key="run_open_path",
                placeholder=tr("/chemin/vers/le/run"),
                help=tr("Accepte le dossier d’un run, son fichier *_results.json ou son projet .aptaswitch.json."),
            )
            st.button(tr("Ouvrir"), key="open_run_path", on_click=_open_path, disabled=_running_job())
        if _running_job():
            st.caption(tr("L’ouverture d’un autre run sera disponible à la fin du calcul en cours."))
    if st.session_state.get("run_open_error"):
        st.error(tr(st.session_state.run_open_error))


def render_saved_run_details(kind: str) -> None:
    """Show the preserved parameters and engine log, never re-run the engine."""
    source = st.session_state.get(f"opened_{kind}_source")
    job_key, result_key, export_key = (
        ("design_job", "last_result", "last_exports") if kind == "design"
        else ("extension_job", "extension_result", "extension_exports")
    )
    result = st.session_state.get(result_key)
    job = st.session_state.get(job_key)
    snapshot = job.snapshot() if kind == "extension" and job is not None else None
    if result is None:
        return
    if kind == "design" and (not source or job is not None):
        return
    if snapshot is not None and snapshot.state == "running":
        return
    if source and job is None:
        st.caption(tr("Run enregistré ouvert : {source} · résultats conservés, sans nouveau calcul.", source=source))
    request = result.request.to_dict() if kind == "design" else result.get("request", {})
    with st.expander(tr("Paramètres enregistrés et console du run")):
        if kind == "extension":
            inputs = result.get("input", {})
            model = result.get("model", {})
            side = inputs.get("extension_side", request.get("extension_side", "5prime"))
            st.caption(tr("Extrémité choisie pour ce run : {side}", side=tr(EXTENSION_SIDES.get(side, side))))
            if model:
                st.caption(
                    f"{model.get('material', '')} · {model.get('temperature_c', '')} °C · "
                    f"Na⁺ {model.get('sodium_m', '')} M · Mg²⁺ {model.get('magnesium_m', '')} M"
                )
            for warning in result.get("warnings", []):
                st.caption(tr(warning))
            st.json({
                tr("paramètres"): request,
                tr("entrée"): inputs,
                tr("modèle"): model,
                tr("criblage"): result.get("mfe_screen", {}),
            }, expanded=False)
        elif request:
            st.json(request, expanded=False)
        log = st.session_state.get(export_key, {}).get("log")
        log_text = None
        if log and Path(log).is_file():
            try:
                log_text = Path(log).read_text(encoding="utf-8", errors="replace")
            except OSError as exc:
                st.info(tr("Le journal n’est plus accessible : {error}", error=exc))
        if log_text is None and snapshot is not None and snapshot.log:
            log_text = "\n".join(snapshot.log)
        if log_text is not None:
            st.code("\n".join(tr(line) for line in log_text.splitlines()), language=None, wrap_lines=True, height=240)
        elif not log or not Path(log).is_file():
            st.caption(tr("Aucun journal enregistré à côté de ce fichier de résultats."))
