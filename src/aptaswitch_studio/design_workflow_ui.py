"""Lazy navigation through the initial system, optional extension and switch design."""

import streamlit as st

from aptaswitch_core.localization import tr
from .extension_ui import render_extension_controls
from .inverse_extension_ui import render_inverse_extension_controls


def _toggle_extension() -> None:
    enabled = st.session_state.extension_enabled
    if enabled:
        st.session_state.design_section = "extension"
    else:
        st.session_state.inverse_extension_enabled = False
        if st.session_state.design_section == "extension":
            st.session_state.design_section = "switch"


def select_design_section(section: str) -> None:
    st.session_state.design_section = section


def render_design_navigation() -> str:
    """Select exactly one design page without running hidden analyses."""
    st.checkbox(tr("Activer l’extension du trigger"), key="extension_enabled",
                on_change=_toggle_extension)
    if not st.session_state.extension_enabled:
        st.session_state.inverse_extension_enabled = False
    labels = {"system": tr("1 · Aptamère et trigger")}
    if st.session_state.extension_enabled:
        labels["extension"] = tr("2 · Extension du trigger")
    labels["switch"] = tr("{step} · Conception du toehold switch", step=len(labels) + 1)
    if st.session_state.design_section not in labels:
        st.session_state.design_section = "system"
    st.radio(tr("Étape de conception"), tuple(labels), key="design_section",
             format_func=labels.__getitem__, horizontal=True, label_visibility="collapsed")
    if st.session_state.design_section == "switch" and st.session_state.inverse_extension_enabled:
        st.caption(tr("Extension non exhaustive activée · trigger cible : {length} nt.",
                      length=st.session_state.inverse_extension_target_length))
    return st.session_state.design_section


def render_extension_page(aptamer: str, trigger: str) -> None:
    st.subheader(tr("Extension du trigger"),
                 help=tr("Les séquences et les conditions sont celles du système initial."))
    with st.expander(tr("Séquences utilisées")):
        for label, sequence in (("Aptamère", aptamer), ("Trigger", trigger)):
            st.caption(f"{tr(label)} · {len(sequence)} nt")
            if sequence:
                st.code(sequence, language=None, wrap_lines=True)
        st.button(tr("Modifier les séquences et les conditions"),
                  on_click=select_design_section, args=("system",), width="stretch")
    if not aptamer or not trigger:
        st.info(tr("Renseignez des séquences valides dans la conception avant de lancer une extension."))
    methods = {False: tr("Exhaustive · avant la conception"),
               True: tr("Non exhaustive · avec la conception")}
    st.radio(tr("Méthode d’extension"), tuple(methods), key="inverse_extension_enabled",
             format_func=methods.__getitem__, horizontal=True)
    if st.session_state.inverse_extension_enabled:
        render_inverse_extension_controls(trigger)
        st.button(tr("Continuer vers la conception du switch"),
                  on_click=select_design_section, args=("switch",), width="stretch")
    else:
        render_extension_controls(aptamer, trigger)
