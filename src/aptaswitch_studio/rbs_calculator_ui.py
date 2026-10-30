"""Recorded OSTIR comparison; viewing never starts or rescales a calculation."""
from __future__ import annotations

import streamlit as st

from aptaswitch_core.localization import tr


def render_rbs_calculator(candidate):
    data = getattr(candidate, "rbs_calculator", {})
    if not data:
        return
    with st.container(border=True):
        st.subheader(tr("Initiation de traduction · OSTIR"), anchor=False)
        st.caption(tr("Comparaison exploratoire : ARN seul OFF / ARN contraint ON. Les valeurs ne sont pas des mesures de production de protéine."))
        off, on, ratio, merit = st.columns(4)
        off.metric(tr("OFF · référence ARN"), f"{data['off_rate']:.3g}")
        on.metric(tr("ON · ARN contraint"), f"{data['on_rate']:.3g}")
        logarithm = data['log10_on_off']
        ratio.metric(tr("Rapport ON/OFF exploratoire"),
                     f"{10 ** logarithm:.3g}×" if abs(logarithm) < 6 else f"10^{logarithm:.2f}×")
        merit.metric(tr("Indice RBS (/100)"), f"{100 * data['merit']:.2f}", help=tr(
            "Plus élevé est meilleur. L’indice combine l’activité ON et le contraste ON/OFF ; un ratio élevé entre deux états presque inactifs reste peu favorisé. Ce n’est pas un pourcentage de fonctionnement."))
        with st.expander(tr("Hypothèses et paramètres OSTIR")):
            st.write(tr("L’ARN du candidat est conservé. En ON, seules les bases censées être occupées par le trigger sont empêchées de former d’autres paires. Le RBS reste disponible pour le ribosome. Le calcul suppose le trigger fixé : il ne calcule ni l’énergie de liaison ADN–ARN, ni la concentration du complexe."))
            st.caption(tr("Référence distincte des conditions du run : E. coli · 37 °C · Na⁺ 1,021 M · OSTIR 1.1.3 / ViennaRNA 2.7.2. Taux en unités arbitraires OSTIR."))
            st.caption(tr("Fenêtre du modèle : bases {start} à {end} du switch, autour de l’AUG en position {aug}. Les interactions avec les bases hors de cette fenêtre ne sont pas évaluées.",
                start=data['window_start'] + 1, end=data['window_end'], aug=data['aug_start'] + 1))
            st.latex(r"M_{RBS}=\frac{ON}{ON+K}\max\!\left(0,\frac{ON-OFF}{ON+OFF}\right),\qquad K=e^{7.279194329}")
            st.caption(tr("K est un repère de normalisation du modèle, pas un objectif expérimental. Le score brut reçoit +100 × poids RBS ON/OFF × cet indice entre 0 et 1."))
