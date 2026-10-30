"""Read-only joint opening diagnostic shared by results and saved switches."""
from __future__ import annotations

import streamlit as st

from aptaswitch_core.localization import tr


def render_rbs_accessibility(candidate):
    # A live session can retain candidates created before this optional diagnostic.
    measurement = getattr(candidate, "rbs_accessibility", {})
    if not measurement:
        return
    with st.expander(tr("Ouverture simultanée du RBS → AUG en ON")):
        probability, energy = st.columns(2)
        probability.metric(tr("Région entièrement libre (%)"), f"{100 * measurement['probability']:.3g}")
        opening = measurement['opening_free_energy_kcal_mol']
        energy.metric(tr("Coût d’ouverture (kcal/mol)"), "∞" if opening is None else f"{opening:.3f}")
        trigger_length = len(candidate.trigger_seq)
        st.caption(tr("Bases {start} à {end} du switch, du préfixe RBS à la fin d’AUG.",
                      start=measurement['start'] - trigger_length + 1,
                      end=measurement['end'] - trigger_length))
        st.caption(tr("Probabilité que toutes ces bases soient libres simultanément, dans le complexe ADN–ARN déjà formé. Ce n’est ni le produit des probabilités par base, ni le rendement d’association, ni une prédiction d’expression. Mesure informative, sans malus supplémentaire dans le score."))
        if measurement['constrained_ensemble_empty']:
            st.caption(tr("Aucune structure du modèle ne satisfait cette ouverture dans le complexe formé ; son coût est infini."))
        elif measurement['probability'] == 0:
            st.caption(tr("Probabilité trop petite pour être représentée numériquement ; le coût d’ouverture reste calculé à partir de son logarithme."))
