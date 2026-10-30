"""Shared concentration inputs for the initial system and switch design."""

from __future__ import annotations

import streamlit as st

from aptaswitch_core.localization import tr

from .web_runtime import CONCENTRATION_UNITS, concentration_from_m, concentration_to_m


def convert_concentration_display_unit(value_key: str | tuple[str, ...], unit_key: str, previous_unit_key: str) -> None:
    """Change the displayed unit while preserving the physical concentration."""
    new_unit = st.session_state[unit_key]
    previous_unit = st.session_state.get(previous_unit_key, new_unit)
    if previous_unit != new_unit:
        for key in (value_key,) if isinstance(value_key, str) else value_key:
            value_m = concentration_to_m(st.session_state[key], previous_unit)
            st.session_state[key] = concentration_from_m(value_m, new_unit)
    st.session_state[previous_unit_key] = new_unit


def _display_unit(unit_key: str, default_unit: str) -> tuple[str, str]:
    previous_unit_key = f"{unit_key}_previous"
    if st.session_state.get(unit_key) not in CONCENTRATION_UNITS:
        previous = st.session_state.get(previous_unit_key)
        st.session_state[unit_key] = previous if previous in CONCENTRATION_UNITS else default_unit
    unit = st.session_state[unit_key]
    st.session_state[previous_unit_key] = unit
    return unit, previous_unit_key


def _number(label, value_key, unit, min_m, max_m, *, disabled=False, help=None):
    return st.number_input(
        tr(label), min_value=float(concentration_from_m(min_m, unit)),
        max_value=float(concentration_from_m(max_m, unit)) if max_m is not None else None,
        step={"nM": 1.0, "µM": 0.1, "mM": 0.1, "M": 0.001}[unit],
        format="%.6g", key=value_key, disabled=disabled, help=tr(help) if help else None,
    )


def render_concentration_control(
    label: str, *, value_key: str, unit_key: str,
    min_m: float, max_m: float | None,
    disabled: bool = False,
    value_m: float | None = None, default_unit: str = "µM", help: str | None = None,
) -> float:
    """Align a value and unit selector, returning the concentration in mol/L."""
    unit, previous_unit_key = _display_unit(unit_key, default_unit)
    if value_key not in st.session_state and value_m is not None:
        st.session_state[value_key] = concentration_from_m(value_m, unit)
    value_col, unit_col = st.columns([0.68, 0.32], vertical_alignment="bottom")
    with value_col:
        value = _number(label, value_key, unit, min_m, max_m, disabled=disabled, help=help)
    unit_col.selectbox(
        tr("Unité"), options=CONCENTRATION_UNITS, key=unit_key,
        on_change=convert_concentration_display_unit,
        args=(value_key, unit_key, previous_unit_key),
        help=tr("Unité utilisée pour {label}.", label=tr(label).lower()),
        disabled=disabled,
    )
    return concentration_to_m(value, st.session_state[unit_key])


def render_concentration_range(
    label: str, *, key: str, start_m: float, stop_m: float,
    min_m: float, max_m: float, default_unit: str = "mM", help: str | None = None,
) -> tuple[float, float]:
    """Two bounds with one shared unit; changing units preserves both bounds."""
    unit_key = key + "_unit"
    unit, previous_unit_key = _display_unit(unit_key, default_unit)
    keys = (key + "_min", key + "_max")
    for value_key, value_m in zip(keys, (start_m, stop_m)):
        if value_key not in st.session_state:
            st.session_state[value_key] = concentration_from_m(value_m, unit)
    start_col, stop_col, unit_col = st.columns([.38, .38, .24], vertical_alignment="bottom")
    with start_col:
        start = _number("Début", keys[0], unit, min_m, max_m, help=help)
    with stop_col:
        stop = _number("Fin", keys[1], unit, min_m, max_m, help=help)
    unit_col.selectbox(
        tr("Unité"), CONCENTRATION_UNITS, key=unit_key,
        on_change=convert_concentration_display_unit, args=(keys, unit_key, previous_unit_key),
        help=tr("Unité utilisée pour {label}.", label=tr(label).lower()),
    )
    return concentration_to_m(start, unit), concentration_to_m(stop, unit)
