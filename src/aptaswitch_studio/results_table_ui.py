"""Compact ranking, expandable to the complete saved result without rescoring."""

from __future__ import annotations

import hashlib
import json
import math
from numbers import Real

import pandas as pd
import streamlit as st

from aptaswitch_core.localization import set_language, tr
from aptaswitch_core.scoring import candidate_score_100
from .table_ui import render_readonly_table
from .scoring_help import ASSOCIATION_WITH_APTAMER_LABEL


_PENDING_SELECTION = "results_table_pending_trial"
_COMPACT_COLUMNS = ("Rang", "Essai", "Score", "ON %", ASSOCIATION_WITH_APTAMER_LABEL, "Défaut OFF/ON", "Defect global")
_PAGE_SIZE = 200


def candidate_status(candidate) -> str:
    if not candidate.excluded:
        return tr("Évalué — STOP prématuré" if candidate.has_stop_codon else "Évalué")
    reason = (candidate.exclusion_reason or "critère de criblage").partition(" — ")[0]
    if reason.startswith("Codon STOP prématuré"):
        reason = "STOP prématuré"
    return tr("Exclu — {reason}", reason=tr(reason))


def best_candidates(result) -> list:
    """Keep the recorded order and round 5% of eligible candidates upward."""
    eligible = [candidate for candidate in result.candidates if not candidate.excluded]
    return eligible[: (len(eligible) + 19) // 20]


def _number(value):
    return float(value) if isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value) else None


def result_rows(result, *, compact: bool = False, start: int = 0, stop: int | None = None) -> list[dict]:
    """Preserve recorded numeric values; missing measurements stay blank."""
    ranks = {
        candidate.trial: rank
        for rank, candidate in enumerate((c for c in result.candidates if not c.excluded), start=1)
    }
    rows = []
    candidates = best_candidates(result) if compact else result.candidates
    for candidate in candidates[start:stop]:
        regional = candidate.regional_defects or {}
        accessibility = getattr(candidate, "rbs_accessibility", {})
        rbs = getattr(candidate, "rbs_calculator", {})
        values = {
            "Rang": ranks.get(candidate.trial),
            "Essai": candidate.trial,
            "Statut": candidate_status(candidate),
            "Score": candidate_score_100(candidate, result.request.score_scale),
            "Defect global": _number(candidate.defect),
            "Défaut structure (global)": _number(candidate.defect_structural),
            "Défaut concentration (global)": _number(candidate.defect_concentration),
            "Défaut OFF": _number(regional.get("off")),
            "Défaut ON": _number(regional.get("on")),
            "Défaut OFF/ON": _number(regional.get("balanced")),
            "ON %": _number(candidate.on_yield_pct),
            ASSOCIATION_WITH_APTAMER_LABEL: _number(candidate.leak_pct),
            "ΔG OFF": _number(candidate.delta_g_off),
            "ΔG ON": _number(candidate.delta_g_on),
            "ΔΔG": _number(candidate.ddg_activation),
            "ΔG RBS–linker": _number(candidate.delta_g_rbs_linker),
            "Indice RBS (/100)": 100 * rbs["merit"] if rbs else None,
            "OFF · référence ARN": _number(rbs.get("off_rate")),
            "ON · ARN contraint": _number(rbs.get("on_rate")),
            "log10 ON/OFF · OSTIR": _number(rbs.get("log10_on_off")),
            "RBS → AUG entièrement libre ON (%)": (
                100 * accessibility["probability"] if accessibility else None),
            "ΔG ouverture RBS → AUG ON": _number(accessibility.get("opening_free_energy_kcal_mol")),
            "STOP": tr("oui" if candidate.has_stop_codon else "non"),
            "Premier STOP": candidate.first_stop if candidate.has_stop_codon else None,
        }
        columns = _COMPACT_COLUMNS if compact else values
        rows.append({tr(column): values[column] for column in columns})
    return rows


def _table_identity(result) -> str:
    payload = [result.run_id, str(result.output_dir), [(c.trial, c.excluded) for c in result.candidates]]
    return hashlib.sha256(json.dumps(payload).encode()).hexdigest()[:16]


def _apply_pending_selection(result, identity: str) -> None:
    pending = st.session_state.pop(_PENDING_SELECTION, None)
    if (isinstance(pending, dict) and pending.get("identity") == identity
            and pending.get("trial") in {candidate.trial for candidate in result.candidates}):
        # Apply the explicit dialog choice before the main candidate selector.
        st.session_state.selected_trial = pending["trial"]


def _table_display(result, *, compact: bool, start: int, stop: int):
    frame = pd.DataFrame(result_rows(result, compact=compact, start=start, stop=stop))
    formats = {tr(name): "{:.2f}" for name in ("Score", "ON %", ASSOCIATION_WITH_APTAMER_LABEL)}
    formats.update({tr(name): "{:.0f}" for name in ("Rang", "Essai", "Premier STOP")})
    return frame.style.format({name: value for name, value in formats.items() if name in frame},
                              precision=3, na_rep="—")


def _render_paginated_table(result, *, compact: bool, height: int, key: str) -> None:
    """Bound the static table size without changing the ranking or candidate choice."""
    total = len(best_candidates(result)) if compact else len(result.candidates)
    page = 1
    if total > _PAGE_SIZE:
        page = st.number_input(tr("Page du tableau"), min_value=1,
                               max_value=math.ceil(total / _PAGE_SIZE), value=1, step=1,
                               key=f"{key}_page")
    start = (page - 1) * _PAGE_SIZE
    stop = min(start + _PAGE_SIZE, total)
    if total > _PAGE_SIZE:
        st.caption(tr("Lignes {first} à {last} sur {total}.", first=start + 1, last=stop, total=total))
    render_readonly_table(_table_display(result, compact=compact, start=start, stop=stop),
                          height=height, key=key)


def _open_full_table(result, identity: str) -> None:
    @st.dialog(tr("Tableau complet des candidats"), width="large")
    def show() -> None:
        set_language(st.session_state.get("ui_language", "en"))
        st.caption(tr("Tous les candidats du classement enregistré, y compris les exclus."))
        _render_paginated_table(result, compact=False, height=550,
                                key=f"results_table_full_{identity}")
        trials = [candidate.trial for candidate in result.candidates]
        current = st.session_state.get("selected_trial")
        trial = st.selectbox(tr("Candidat à explorer"), trials,
                             index=trials.index(current) if current in trials else 0,
                             format_func=lambda value: tr("Essai {trial}", trial=value),
                             key=f"results_table_choose_candidate_{identity}")
        if st.button(tr("Explorer ce candidat"), width="stretch",
                     key=f"results_table_open_candidate_{identity}"):
            st.session_state[_PENDING_SELECTION] = {"identity": identity, "trial": trial}
            st.rerun(scope="app")
    show()


def render_results_table(result) -> None:
    """Show a static preview and apply an explicit full-table dialog choice."""
    identity = _table_identity(result)
    _apply_pending_selection(result, identity)
    eligible_count = sum(not candidate.excluded for candidate in result.candidates)
    best = best_candidates(result)
    if best:
        st.caption(tr("Les 5 % meilleurs : {shown} sur {total} candidats classés.",
                      shown=len(best), total=eligible_count))
        key = f"results_table_best_{identity}"
        _render_paginated_table(result, compact=True,
                                height=min(350, 38 + 35 * len(best)), key=key)
    if result.candidates and st.button(tr("Agrandir le tableau complet"),
                                       key=f"results_table_expand_{identity}", width="stretch"):
        _open_full_table(result, identity)
