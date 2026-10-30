"""Read-only presentation of global and functional ensemble defects."""

from __future__ import annotations

import math
from numbers import Real

import streamlit as st

from aptaswitch_core.localization import tr
from .table_ui import render_readonly_table


FUNCTIONAL_ZONES = (
    ("off_toehold", "OFF", "Toehold libre"),
    ("off_stem", "OFF", "Tiges avant/après AUG"),
    ("on_rbs_aug", "ON", "Préfixe RBS → AUG libre"),
    ("on_rbs_linker", "ON", "Préfixe RBS → fin linker libre"),
)
FUNCTIONAL_DEFECT_DEFAULTS = {f"weight_functional_{key}": .75 if key == "on_rbs_linker" else 1.0
                              for key, _, _ in FUNCTIONAL_ZONES}


def _weight_fractions(values: dict, zones=FUNCTIONAL_ZONES) -> dict[str, float | None]:
    numbers = {key: _number(values.get(key)) for key, _, _ in zones}
    if any(value is None or value < 0 for value in numbers.values()):
        return {key: None for key in numbers}
    maximum = max(numbers.values())
    if maximum == 0:
        return {key: None for key in numbers}
    scaled = {key: value / maximum for key, value in numbers.items()}
    total = math.fsum(scaled.values())
    return {key: value / total for key, value in scaled.items()}


def _saved_fractions(regional: dict) -> dict[str, float | None]:
    """Use the weights recorded with these functional measurements."""
    return _weight_fractions(regional.get("normalized_weights", {}))


def _render_weight_summary(fractions: dict, *, saved: bool) -> None:
    if any(value is None for value in fractions.values()):
        return
    off = math.fsum(value for key, value in fractions.items() if key.startswith("off_"))
    on = math.fsum(value for key, value in fractions.items() if key.startswith("on_"))
    if saved:
        st.caption(tr("Poids effectifs enregistrés : OFF {off:g} % · ON {on:g} %.", off=100 * off, on=100 * on))
    else:
        st.caption(tr("Poids effectifs : OFF {off:g} % · ON {on:g} %.", off=100 * off, on=100 * on))
    for state, weight in (("OFF", off), ("ON", on)):
        if weight == 0:
            st.caption(tr("{state} est désactivé dans le compromis (poids nul) : son défaut agrégé vaut « — ». Les défauts locaux restent calculés.", state=state))


def render_regional_defect_controls() -> None:
    """Relative region priorities are separate from the overall score coefficient."""
    st.markdown(tr("#### Répartition du défaut OFF/ON"), help="\n\n".join((
        tr("Par défaut : 26,7 % pour chacun des trois premiers critères et 20 % pour RBS → fin du linker (poids 1 / 1 / 1 / 0,75). Le coefficient « Défaut OFF/ON · malus » règle séparément l’amplitude totale du malus."),
        tr("Les deux critères ON se recouvrent : la région préfixe RBS → AUG contribue aux deux. Cela renforce son importance ; leurs parts du défaut ne doivent pas être additionnées."),
    )))
    for key, default in FUNCTIONAL_DEFECT_DEFAULTS.items():
        st.session_state.setdefault(key, default)
    columns = st.columns(4)
    for column, (key, state, label) in zip(columns, FUNCTIONAL_ZONES):
        column.number_input(
            f"{state} · {tr(label)}", min_value=0.0, max_value=20.0, step=0.1,
            key=f"weight_functional_{key}",
            help=tr("Poids relatif de cette zone. 0 retire sa contribution au compromis ; son défaut local reste calculé."),
        )
    values = {key: st.session_state[f"weight_functional_{key}"] for key, _, _ in FUNCTIONAL_ZONES}
    fractions = _weight_fractions(values, FUNCTIONAL_ZONES)
    for column, (key, _, _) in zip(columns, FUNCTIONAL_ZONES):
        column.caption(tr("Part du compromis : {part}", part=_display(fractions[key], percent=True)))
    if any(value is None for value in fractions.values()):
        st.error(tr("Au moins un poids de zone doit être supérieur à 0. Pour désactiver tout le malus, réglez « Défaut OFF/ON · malus » à 0 en conservant des poids de zones valides."))
    _render_weight_summary(fractions, saved=False)


def _number(value) -> float | None:
    if isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value):
        return float(value)
    return None


def _display(value, *, percent: bool = False) -> str:
    number = _number(value)
    if number is None:
        return "—"
    return f"{100 * number:.4g} %" if percent else f"{number:.6g}"


def render_defect_method(result) -> None:
    """Explain the metric actually used by this saved run, without rescoring it."""
    st.markdown("**" + tr("Défaut OFF/ON") + "**", help="\n\n".join((
        tr("Le malus utilise quatre défauts locaux : toehold libre et tiges encadrant AUG en OFF ; préfixe RBS → fin AUG et préfixe RBS → fin linker libres en ON. Les poids sont ceux enregistrés avec ce calcul."),
        tr(
            "Les défauts sont des fractions de 0 à 1 ; plus faible est meilleur. Le défaut de structure "
            "global et le défaut de concentration global s’additionnent pour donner le Defect global. "
            "Le compromis régional est une mesure indépendante. « — » signifie non calculé, ou état "
            "désactivé lorsque ses deux poids sont nuls."
        ),
    )))
    _render_weight_summary(result.request.functional_defect_weights.fractions(), saved=True)


def _positions(indices) -> str:
    if not isinstance(indices, (list, tuple)) or not indices:
        return "—"
    if any(type(index) is not int or index < 0 for index in indices):
        return "—"
    values = sorted(set(index + 1 for index in indices))
    spans = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        spans.append(str(first) if first == previous else f"{first}–{previous}")
        first = previous = value
    spans.append(str(first) if first == previous else f"{first}–{previous}")
    return ", ".join(spans)


def render_candidate_defects(candidate) -> None:
    """Show saved region diagnostics; never derive ensemble values from the MFE."""
    regional = getattr(candidate, "regional_defects", None) or {}
    zones = FUNCTIONAL_ZONES
    regions = regional.get("regions") or {}
    with st.expander(tr("Détail des défauts OFF/ON")):
        if not regional:
            st.caption(tr("Ce résultat ne contient pas de défauts régionaux enregistrés ; leurs valeurs restent « — » dans le tableau."))
            return
        st.markdown("**" + tr("Défaut d’ensemble (0–1)") + "**", help="\n\n".join((
            tr("Ces défauts structuraux sont calculés sur l’ensemble thermodynamique des conformations, "
               "pas uniquement sur la structure MFE affichée. Les parts du compromis ci-dessous "
               "reprennent les poids enregistrés ; elles ne dépendent pas de la longueur des zones."),
            tr("La part du défaut rapporte les erreurs de la zone aux erreurs de tout le switch "
               "dans le même état, rapporteur inclus. Elle ne représente pas son poids dans le score ; "
               "elle vaut « — » lorsque le défaut total du switch est nul."),
            tr("Les deux critères ON se recouvrent : la région préfixe RBS → AUG contribue aux deux. Cela renforce son importance ; leurs parts du défaut ne doivent pas être additionnées."),
        )))
        fractions = _saved_fractions(regional)
        _render_weight_summary(fractions, saved=True)
        rows = []
        for key, state, label in zones:
            detail = regions.get(key) or {}
            rows.append({
                tr("État"): state,
                tr("Zone"): tr(label),
                tr("Défaut d’ensemble (0–1)"): _display(regional.get(key, detail.get("normalized_defect"))),
                tr("Part du compromis (%)"): _display(fractions[key], percent=True),
                tr("Nucléotides"): _display(detail.get("nucleotide_count")),
                tr("Nucléotides incorrects attendus"): _display(detail.get("expected_incorrect_nt")),
                tr("Part du défaut du switch (%)"): _display(detail.get("share_of_switch_defect"), percent=True),
                tr("Positions du switch (à partir de 1)"): _positions(detail.get("indices_switch")),
            })
        render_readonly_table(rows)
        parts = regional.get("off_stem_parts") or {}
        st.markdown("**" + tr("Tiges OFF") + "**", help="\n\n".join((
            tr("Détail des tiges OFF : 6 paires avant AUG et 9 après AUG dans le classique 30 nt. Les longueurs suivent l’architecture personnalisée."),
            tr("Le défaut des tiges est la moyenne sur leurs paires cibles : contributions 6/15 et 9/15 dans le classique. AUG et le bulge restent non appariés dans la cible OFF ; le maintien des tiges mesure leur séquestration indirecte."),
        )))
        render_readonly_table([{
            tr("Tige OFF"): tr(label),
            tr("Paires cibles"): _display((parts.get(key) or {}).get("pair_count")),
            tr("Défaut d’ensemble (0–1)"): _display((parts.get(key) or {}).get("normalized_defect")),
            tr("Positions du switch (à partir de 1)"): _positions((parts.get(key) or {}).get("indices_switch")),
        } for key, label in (("before_aug", "Avant AUG"), ("after_aug", "Après AUG"))])
