"""Display bounded scores with their recorded comparison sum."""
from __future__ import annotations

import math
from numbers import Real

import streamlit as st

from aptaswitch_core.localization import tr
from aptaswitch_core.scoring import candidate_score_100


def _display(value, *, percent=False):
    if not isinstance(value, Real) or isinstance(value, bool) or not math.isfinite(value):
        return "—"
    return f"{100*value:.2f} %" if percent else f"{value:.4f}"


def render_candidate_score(result, candidate):
    if candidate.excluded:
        return
    scale = result.request.score_scale
    raw = getattr(candidate, "score_raw", None)
    with st.expander(tr("Détail du score /100")):
        columns = st.columns(2)
        columns[0].metric(tr("Score /100"), _display(candidate_score_100(candidate, scale)),
                          help=tr("La note /100 convertit la somme pondérée avec une échelle fixe, indépendante des autres candidats. Elle ne mesure pas une probabilité de fonctionnement."))
        columns[1].metric(tr("Score brut"), _display(raw))
        if result.request.inverse_extension is None:
            return
        weight = getattr(result.request.scoring, "trigger_linearity", 0.)
        st.caption(tr("Poids enregistré du malus de linéarité du trigger seul : {weight:g}.", weight=weight))
        diagnostics = getattr(candidate, "trigger_linearity", {}) or {}
        if not diagnostics:
            st.caption(tr("La linéarité du trigger seul n’a pas été enregistrée pour ce candidat."))
            return
        columns = st.columns(2)
        columns[0].metric(tr("Trigger seul · bases libres"), _display(diagnostics.get("unpaired_fraction"), percent=True),
                          help=tr("Fraction moyenne de bases libres dans l’ensemble thermodynamique du trigger isolé, aux conditions du design du switch. La structure ci-dessous est sa MFE enregistrée."))
        columns[1].metric(tr("Trigger seul · défaut de linéarité"), _display(diagnostics.get("normalized_defect")))
        if diagnostics.get("structure_mfe"):
            st.code(diagnostics["structure_mfe"], language=None, wrap_lines=True)
