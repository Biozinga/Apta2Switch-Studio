"""Controls for durable design checkpoints, separate from inverse reanalysis."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from aptaswitch_core.design_checkpoint import CHECKPOINT_NAME
from aptaswitch_core.localization import set_language, tr
from aptaswitch_core.thermodynamics.setup import configured_library_path

from .design_session import active_job_running
from . import web_runtime


def _control_design(action: str) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    job = st.session_state.get("design_job")
    if job is None or job.snapshot().state != "running":
        return
    try:
        getattr(job, {"pause": "pause", "finish": "finish_early", "cancel": "cancel"}[action])()
    except (OSError, ValueError, RuntimeError) as exc:
        st.session_state.run_error = str(exc)


def render_design_run_controls(job, snapshot) -> None:
    """Keep the requested action visible while the current work finishes."""
    if getattr(job.request, "fixed_switch_sequence", None) is not None:
        st.button(tr("Mettre en pause"), key="pause_design_run", width="stretch",
                  disabled=bool(snapshot.control_requested) or not snapshot.can_request_control,
                  on_click=_control_design, args=("pause",))
        st.caption(tr("La pause conserve l’analyse terminée et ses paramètres. Un calcul natif en cours se termine avant l’enregistrement."))
        return
    if not job.supports_continuation:
        st.button(tr("Annuler le calcul"), key="cancel_design_reanalysis",
                  on_click=_control_design, args=("cancel",))
        return
    pending = snapshot.control_requested
    finalizing = not snapshot.can_request_control
    controls = st.columns(2)
    controls[0].button(tr("Terminer maintenant"), key="finish_design_run", width="stretch",
                       disabled=bool(pending) or finalizing, on_click=_control_design, args=("finish",))
    controls[1].button(tr("Mettre en pause"), key="pause_design_run", width="stretch",
                       disabled=bool(pending) or finalizing, on_click=_control_design, args=("pause",))
    if finalizing:
        st.info(tr("Finalisation en cours : comparaison des extensions si prévue, classement et exports. Aucun nouvel essai ne sera lancé."))
    elif pending == "pause":
        st.info(tr("Pause demandée : l’essai en cours se termine, puis les résultats et les paramètres sont enregistrés sur disque. Vous pourrez reprendre ce calcul après avoir fermé l’application."))
    elif pending == "finish":
        st.info(tr("Fin demandée : l’essai en cours se termine, puis les candidats obtenus sont classés, la comparaison des extensions est finalisée si elle est prévue, et les exports sont créés."))
    else:
        st.caption(tr("Ces actions attendent la fin de l’essai en cours. Terminer prépare le classement et les exports ; mettre en pause enregistre le calcul pour le reprendre plus tard."))


def _resume_design(path: str, finish: bool) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.run_error = ""
    try:
        if active_job_running(st.session_state):
            raise ValueError(tr("Un calcul est déjà en cours. Suivez-le dans Résultats."))
        job = web_runtime.resume_design_job(
            Path(path), finish_early=finish,
            library_path=configured_library_path(),
        )
        st.session_state.design_job = job
        st.session_state.last_exports = {}
        st.session_state.run_notice = ""
    except (OSError, ValueError, RuntimeError) as exc:
        st.session_state.run_error = str(exc)


def render_design_continuation(result) -> bool:
    """Show only resumable checkpoints; never rebuild a request from form fields."""
    metadata = result.optimization.get("checkpoint", {})
    unfinished = metadata.get("phase") in {"design", "inverse_extension", "exports"}
    if not unfinished:
        if result.status == "stopped":
            st.info(tr("Calcul terminé à votre demande · {completed}/{total} essais terminés.",
                       completed=len(result.optimization.get("trials", [])), total=result.request.trials))
        return False
    if getattr(result.request, "fixed_switch_sequence", None) is not None:
        completed = result.optimization.get("evaluations", 0)
        job = st.session_state.get("design_job")
        if job is not None and job.snapshot().state == "running":
            st.caption(tr("Lecture de séquences · {completed}/1 analyse terminée. Le calcul se poursuit.",
                          completed=completed))
            return True
        with st.container(border=True):
            st.markdown(tr("**Lecture de séquences · {completed}/1 analyse terminée**", completed=completed))
            st.caption(tr("La reprise conserve les séquences, les conditions et les paramètres de ce calcul. Les changements du formulaire ne s’appliquent pas à cette reprise."))
            path = Path(result.output_dir) / CHECKPOINT_NAME
            st.button(tr("Reprendre le calcul"), key="resume_design_run", type="primary",
                      width="stretch", disabled=not path.is_file(),
                      on_click=_resume_design, args=(str(path), False))
            if not path.is_file():
                st.error(tr("Le point de reprise est introuvable. Rouvrez le calcul depuis son dossier d’origine."))
        return True
    completed = metadata.get("completed_trials", len(result.optimization.get("trials", [])))
    job = st.session_state.get("design_job")
    if job is not None and job.snapshot().state == "running":
        st.caption(tr("Résultats partiels enregistrés · {completed}/{total} essais terminés. Le calcul se poursuit.",
                      completed=completed, total=result.request.trials))
        return True
    with st.container(border=True):
        label = "Calcul en pause" if result.status == "paused" else "Calcul interrompu"
        st.markdown(tr("**{state} · {completed}/{total} essais terminés**", state=tr(label),
                       completed=completed, total=result.request.trials))
        st.caption(tr("La reprise conserve les séquences, les conditions et les paramètres de ce calcul. Les changements du formulaire ne s’appliquent pas à cette reprise."))
        if result.status != "paused":
            st.caption(tr("Après une fermeture inattendue, seul un essai non encore enregistré peut être recommencé."))
        path = Path(result.output_dir) / CHECKPOINT_NAME
        missing = not path.is_file()
        controls = st.columns(2)
        controls[0].button(tr("Reprendre le calcul"), key="resume_design_run", type="primary",
                           width="stretch", disabled=missing, on_click=_resume_design, args=(str(path), False))
        controls[1].button(tr("Terminer avec les résultats actuels"), key="finish_saved_design_run",
                           width="stretch", disabled=missing, on_click=_resume_design, args=(str(path), True))
        st.caption(tr("Terminer avec les résultats actuels finalise la comparaison des extensions prévue et les exports, sans créer de nouvel essai."))
        if missing:
            st.error(tr("Le point de reprise est introuvable. Rouvrez le calcul depuis son dossier d’origine."))
    return True
