"""Named designs: the local home screen and project navigation."""

from __future__ import annotations

from datetime import datetime
from html import escape
from pathlib import Path

import streamlit as st

from aptaswitch_core.localization import tr
from .design_library import (
    DesignExistsError, create_design, default_design_root, export_design,
    import_design, list_designs, load_design,
)
from .design_session import (
    active_job_running, open_design_into_session,
    save_active_design,
)
from .run_library_ui import _show_folder


def persist_active_design() -> bool:
    """Keep recoverable disk failures visible without discarding the form."""
    try:
        save_active_design(st.session_state)
        st.session_state.pop("design_save_error", None)
        return True
    except (OSError, ValueError) as exc:
        st.session_state.design_save_error = str(exc)
        return False


def _home() -> None:
    persist_active_design()
    st.session_state.nav = "Mes designs"


def _activate(path: Path) -> None:
    open_design_into_session(path, st.session_state)
    st.session_state.nav = "Conception"
    st.session_state.pop("design_library_error", None)
    st.session_state.pop("design_pending_overwrite", None)
    st.session_state.pop("design_export_bytes", None)


def _open_selected() -> None:
    try:
        _activate(Path(st.session_state.design_to_open))
    except (OSError, ValueError) as exc:
        st.session_state.design_library_error = str(exc)


def _perform(action: dict, *, overwrite: bool = False) -> None:
    if active_job_running(st.session_state):
        st.session_state.design_library_error = tr("Attendez la fin du calcul avant de changer de design.")
        return
    if not persist_active_design():
        return
    root = default_design_root()
    try:
        if action["kind"] == "create":
            project = create_design(root, action["name"], overwrite=overwrite,
                                    state={"form": {"molecule_name": action["name"],
                                                    "design_mode": action.get("mode", "design")}})
        else:
            project = import_design(action["source"], root, name=action.get("name") or None,
                                    overwrite=overwrite)
        # The previous project was saved above. Do not re-save its state over
        # a newly confirmed replacement of that same folder when opening it.
        st.session_state.pop("active_design_path", None)
        _activate(project.path)
    except DesignExistsError as exc:
        st.session_state.design_pending_overwrite = action
        st.session_state.design_collision_message = str(exc)
        st.session_state.pop("design_library_error", None)
    except (OSError, ValueError) as exc:
        st.session_state.design_library_error = str(exc)


def _create() -> None:
    _perform({"kind": "create", "name": st.session_state.get("design_new_name", "").strip(),
              "mode": st.session_state.get("design_new_mode", "design")})


def _import_path() -> None:
    source = st.session_state.get("design_import_path", "").strip()
    if not source:
        st.session_state.design_library_error = tr("Indiquez le dossier ou l’archive du design à importer.")
        return
    _perform({"kind": "import", "source": Path(source).expanduser(),
              "name": st.session_state.get("design_import_name", "").strip()})


def _import_upload() -> None:
    uploaded = st.session_state.get("design_import_file")
    if uploaded is not None:
        _perform({"kind": "import", "source": uploaded.getvalue(),
                  "name": st.session_state.get("design_import_name", "").strip()})


def _confirm_overwrite() -> None:
    action = st.session_state.get("design_pending_overwrite")
    if action:
        _perform(action, overwrite=True)


def _cancel_overwrite() -> None:
    st.session_state.pop("design_pending_overwrite", None)
    st.session_state.pop("design_collision_message", None)


def _prepare_export() -> None:
    if not persist_active_design():
        return
    try:
        project = load_design(Path(st.session_state.active_design_path))
        st.session_state.design_export_bytes = export_design(project.path)
        st.session_state.design_export_revision = project.updated_at
        st.session_state.design_export_name = project.name + ".zip"
    except (OSError, ValueError) as exc:
        st.session_state.design_library_error = str(exc)


def render_design_toolbar() -> None:
    if not st.session_state.get("active_design_path"):
        return
    name = st.session_state.get("active_design_name", "")
    with st.container(key="active_design_banner", border=False):
        title, home, actions = st.columns([5, 2, 2], vertical_alignment="center")
        title.html(
            '<div class="active-design-identity" title="'
            + escape(tr("Sauvegarde automatique dans le dossier de ce design."), quote=True)
            + '"><div class="active-design-label">' + escape(tr("Design actif"))
            + '</div><div class="active-design-name">' + escape(str(name)) + '</div></div>'
        )
        home.button(tr("Mes designs"), key="design_home", on_click=_home, width="stretch")
        with actions.popover(tr("Dossier et export"), width="stretch"):
            st.caption(st.session_state.active_design_path)
            st.button(tr("Afficher le dossier"), key="show_design_folder", on_click=_show_folder,
                      args=(Path(st.session_state.active_design_path),), width="stretch")
            running = active_job_running(st.session_state)
            st.button(tr("Préparer l’export du design"), key="prepare_design_export",
                      on_click=_prepare_export, disabled=running, width="stretch")
            data = st.session_state.get("design_export_bytes")
            if data:
                revision = load_design(Path(st.session_state.active_design_path)).updated_at
                if revision == st.session_state.get("design_export_revision"):
                    st.download_button(tr("Télécharger le design (.zip)"), data,
                                       file_name=st.session_state.design_export_name,
                                       mime="application/zip", key="export_design", width="stretch")
                else:
                    st.caption(tr("Le design a changé. Préparez un nouvel export."))
    if st.session_state.get("design_restore_warnings"):
        for warning in st.session_state.design_restore_warnings:
            st.warning(str(warning))


def render_design_library() -> None:
    st.title(tr("Mes designs"), anchor=False,
             help=tr("Retrouvez un design, créez-en un nouveau ou importez un dossier complet. Séquences, extensions, calculs et structures 3D restent réunis."))
    root = default_design_root()
    running = active_job_running(st.session_state)
    if running:
        st.info(tr("Un calcul est en cours. Vous pourrez changer de design lorsqu’il sera terminé."))
        st.button(tr("Revenir au design en cours"), key="resume_active_design",
                  on_click=lambda: st.session_state.update(nav="Conception"), width="stretch")
    try:
        projects = list_designs(root)
    except (OSError, ValueError) as exc:
        st.error(str(exc))
        projects = []
    with st.container(border=True):
        st.subheader(tr("Designs enregistrés"))
        if projects:
            by_path = {str(project.path): project for project in projects}
            if st.session_state.get("design_to_open") not in by_path:
                st.session_state.design_to_open = next(iter(by_path))
            st.selectbox(tr("Choisir un design"), list(by_path), key="design_to_open",
                         format_func=lambda path: by_path[path].name)
            selected = by_path[st.session_state.design_to_open]
            form = selected.state.get("form", {})
            if not isinstance(form, dict):
                form = {}
            st.caption(tr("Modifié le {date} · Aptamère : {aptamer} nt · Trigger : {trigger} nt",
                          date=datetime.fromisoformat(selected.updated_at).astimezone().strftime("%Y-%m-%d %H:%M"),
                          aptamer=len("".join(str(form.get("aptamer_input", "")).split())),
                          trigger=len("".join(str(form.get("trigger_input", "")).split()))))
            if form.get("design_mode") == "reading":
                st.caption(tr("Lecture de séquences · Switch : {switch} nt",
                              switch=len("".join(str(form.get("reading_switch_input", "")).split()))))
            st.button(tr("Ouvrir ce design"), key="open_design", on_click=_open_selected,
                      disabled=running, type="primary", width="stretch")
        else:
            st.info(tr("Aucun design enregistré pour le moment. Créez votre premier design ci-dessous."))
    create_tab, import_tab = st.tabs([tr("Nouveau design"), tr("Importer un design")])
    with create_tab:
        st.text_input(tr("Nom du design"), key="design_new_name", placeholder=tr("Ex. : Adénosine — essai 1"))
        modes = {"design": tr("Conception"), "reading": tr("Lecture de séquences")}
        st.radio(tr("Mode du design"), tuple(modes), key="design_new_mode",
                 format_func=modes.__getitem__, horizontal=True, disabled=running)
        st.button(tr("Créer le design"), key="create_design", on_click=_create,
                  disabled=running, type="primary", width="stretch")
    with import_tab:
        st.text_input(tr("Nom après import (facultatif)"), key="design_import_name")
        st.file_uploader(tr("Archive du design (.zip)"), type=["zip"], key="design_import_file")
        st.button(tr("Importer l’archive"), key="import_design_upload", on_click=_import_upload,
                  disabled=running or st.session_state.get("design_import_file") is None,
                  type="primary", width="stretch")
        with st.expander(tr("Importer depuis un dossier local")):
            st.text_input(tr("Dossier du design, design.json ou archive ZIP"), key="design_import_path")
            st.button(tr("Importer ce design"), key="import_design_path", on_click=_import_path,
                      disabled=running, type="primary", width="stretch")
    if st.session_state.get("design_pending_overwrite"):
        st.warning(tr("Attention, ce design existe déjà. Voulez-vous l’écraser ? Une copie de l’ancien dossier sera conservée."))
        st.caption(st.session_state.get("design_collision_message", ""))
        yes, no = st.columns(2)
        yes.button(tr("Oui, écraser ce design"), key="overwrite_design_confirm",
                   on_click=_confirm_overwrite, disabled=running, type="primary", width="stretch")
        no.button(tr("Annuler"), key="overwrite_design_cancel", on_click=_cancel_overwrite, width="stretch")
    if st.session_state.get("design_library_error"):
        st.error(st.session_state.design_library_error)
    st.caption(tr("Dossier des designs : {path}", path=root))
