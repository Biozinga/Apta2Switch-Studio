"""Shared concentration inputs for the initial system and switch design."""

from __future__ import annotations

import streamlit as st

from aptaswitch_core.localization import tr

from .web_runtime import CONCENTRATION_UNITS, concentration_from_m, concentration_to_m


def convert_concentration_display_unit(value_key: str, unit_key: str, previous_unit_key: str) -> None:
    """Change the displayed unit while preserving the physical concentration."""
    new_unit = st.session_state[unit_key]
    previous_unit = st.session_state.get(previous_unit_key, new_unit)
    if previous_unit != new_unit:
        value_m = concentration_to_m(st.session_state[value_key], previous_unit)
        st.session_state[value_key] = concentration_from_m(value_m, new_unit)
    st.session_state[previous_unit_key] = new_unit


def render_concentration_control(
    label: str, *, value_key: str, unit_key: str,
    min_m: float, max_m: float | None,
) -> float:
    """Align a value and unit selector, returning the concentration in mol/L."""
    if st.session_state.get(unit_key) not in CONCENTRATION_UNITS:
        st.session_state[unit_key] = "µM"
    unit = st.session_state[unit_key]
    previous_unit_key = f"{unit_key}_previous"
    if previous_unit_key not in st.session_state:
        st.session_state[previous_unit_key] = unit
    value_col, unit_col = st.columns([0.68, 0.32], vertical_alignment="bottom")
    value = value_col.number_input(
        tr(label),
        min_value=float(concentration_from_m(min_m, unit)),
        max_value=float(concentration_from_m(max_m, unit)) if max_m is not None else None,
        step={"nM": 1.0, "µM": 0.1, "mM": 0.1, "M": 0.001}[unit],
        format="%.6g", key=value_key,
    )
    unit_col.selectbox(
        tr("Unité"), options=CONCENTRATION_UNITS, key=unit_key,
        on_change=convert_concentration_display_unit,
        args=(value_key, unit_key, previous_unit_key),
        help=tr("Unité utilisée pour {label}.", label=tr(label).lower()),
    )
    value_m = concentration_to_m(value, st.session_state[unit_key])
    st.caption(tr("Équivalent : {value:.6g} M", value=value_m))
    return value_m
