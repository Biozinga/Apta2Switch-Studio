"""Shared display of saved thermodynamic tube populations, without new calculations."""

from __future__ import annotations

import math

import altair as alt
import pandas as pd
import streamlit as st

from aptaswitch_core.localization import tr
from .table_ui import render_readonly_table


_CATEGORY_LABELS = {
    "free": "Libre",
    "target": "Complexe A + T",
    "parasite": "Complexes parasites",
}
_CATEGORY_COLORS = {"free": "#94A3B8", "target": "#14B8A6", "parasite": "#F59E0B"}
_SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _format_number(value: float) -> str:
    """Keep very small populations visible, without implying integer precision."""
    value = float(value)
    if not math.isfinite(value):
        return "—"
    if value == 0:
        return "0"
    if abs(value) >= 1e6 or abs(value) < 0.001:
        mantissa, exponent = f"{value:.2e}".split("e")
        return f"{mantissa} × 10{str(int(exponent)).translate(_SUPERSCRIPTS)}"
    if abs(value) >= 1000:
        rounded = float(f"{value:.3g}")
        return f"{rounded:,.0f}".replace(",", " ")
    return f"{value:.3g}"


def _format_pct(value: float) -> str:
    value = float(value)
    if not math.isfinite(value):
        return "—"
    if 0 < value < 0.001:
        return "< 0.001 %"
    return f"{value:.3f}".rstrip("0").rstrip(".") + " %"


def _composition_rows(tube: dict, *, mixed: bool) -> list[dict]:
    """Plot fractions of the input strands, never fractions of particles."""
    summary = tube["summary"]
    rows = []
    for species, name in (("trigger", "Trigger"), ("aptamer", tr("Aptamère"))):
        if species == "aptamer" and not mixed:
            continue
        for order, category in enumerate(_CATEGORY_LABELS):
            if category == "target" and not mixed:
                continue
            rows.append({
                "species": name,
                "category": tr(_CATEGORY_LABELS[category]),
                "fraction_pct": float(summary[f"{species}_{category}_pct"]),
                "order": order,
            })
    return rows


def _composition_chart(tube: dict, *, mixed: bool):
    categories = [category for category in _CATEGORY_LABELS if mixed or category != "target"]
    return (
        alt.Chart(pd.DataFrame(_composition_rows(tube, mixed=mixed)))
        .mark_bar(size=26, cornerRadius=3)
        .encode(
            x=alt.X("fraction_pct:Q", stack="zero", scale=alt.Scale(domain=[0, 100]),
                    title=tr("Part des brins introduits (%)"), axis=alt.Axis(values=[0, 25, 50, 75, 100])),
            y=alt.Y("species:N", sort=["Trigger", tr("Aptamère")], title=None,
                    axis=alt.Axis(labelFontSize=13)),
            color=alt.Color("category:N", title=None,
                            scale=alt.Scale(domain=[tr(_CATEGORY_LABELS[item]) for item in categories],
                                            range=[_CATEGORY_COLORS[item] for item in categories]),
                            legend=alt.Legend(orient="bottom")),
            order=alt.Order("order:Q"),
            tooltip=[alt.Tooltip("species:N", title=tr("Brin")),
                     alt.Tooltip("category:N", title=tr("Forme")),
                     alt.Tooltip("fraction_pct:Q", title="%", format=".4g")],
        )
        .properties(height=92 if mixed else 52)
        .configure_view(stroke=None)
    )


def _reference_delta(reference: dict | None, tube_key: str, current: float) -> str | None:
    if not reference or reference.get("status") != "completed":
        return None
    previous = reference.get("tubes", {}).get(tube_key, {}).get("summary", {}).get("trigger_parasite_pct")
    if previous is None:
        return None
    difference = current - float(previous)
    return tr("{value} points vs original", value=f"{difference:+.3f}")


def _reference_is_compatible(analysis: dict, reference: dict) -> bool:
    """Compare equilibrium fractions only under identical simulation settings.

    Tube volume only converts concentrations into expected particle counts; it
    does not change the percentages used by this comparison.
    """
    if analysis.get("model", {}) != reference.get("model", {}):
        return False
    conditions, previous = analysis.get("conditions", {}), reference.get("conditions", {})
    return all(
        conditions.get(field) is not None and conditions[field] == previous.get(field)
        for field in ("aptamer_concentration_m", "trigger_concentration_m", "max_complex_size")
    )


def _render_tube(tube: dict, *, mixed: bool, key: str, delta: str | None) -> None:
    summary = tube["summary"]
    metrics = st.columns(3)
    metrics[0].metric(
        tr("Trigger dans des complexes parasites"), _format_pct(summary["trigger_parasite_pct"]),
        delta=delta, delta_color="inverse",
        help="\n\n".join((
            tr("Pourcentage des brins de trigger introduits qui appartiennent à des complexes parasites. Un dimère T + T compte deux brins de trigger."),
            tr("Prédiction à l’équilibre, sans ligand. Les pourcentages portent sur les brins introduits ; T + T contient deux triggers. Dans le tube mixte, A + T est le complexe attendu ; les autres assemblages de plusieurs brins sont classés comme parasites."),
        )),
    )
    if mixed:
        metrics[1].metric(tr("Trigger dans A + T"), _format_pct(summary["trigger_target_pct"]),
                          help=tr("Complexe contenant exactement un aptamère et un trigger. Cette composition ne garantit ni la liaison au ligand ni un appariement correct."))
    else:
        metrics[1].metric(tr("Complexes parasites estimés"), _format_number(summary["parasite_complex_molecule_count"]),
                          help=tr("Nombre attendu de complexes parasites dans le volume indiqué. Un dimère compte ici comme un seul complexe."))
    metrics[2].metric(tr("Trigger libre"), _format_pct(summary["trigger_free_pct"]))

    st.altair_chart(_composition_chart(tube, mixed=mixed), width="stretch", key=f"{key}_composition")
    st.caption(tr(
        "Dans ce volume : ≈ {total} entités moléculaires, dont ≈ {parasites} complexes parasites.",
        total=_format_number(summary["total_complex_molecule_count"]),
        parasites=_format_number(summary["parasite_complex_molecule_count"]),
    ))
    with st.expander(tr("Composition détaillée du tube")):
        st.caption(tr("A = aptamère · T = trigger. Les formes ayant la même composition sont regroupées. Les nombres attendus sont calculés à partir des concentrations et du volume du tube."))
        records = []
        for row in tube["complexes"]:
            record = {
                tr("Forme"): row["label"],
                tr("Classe"): tr(_CATEGORY_LABELS[row["category"]]),
                tr("Concentration (µM)"): _format_number(row["concentration_m"] * 1e6),
                tr("Nombre attendu"): _format_number(row["complex_molecule_count"]),
                tr("% des triggers"): _format_pct(row["trigger_fraction_pct"]),
            }
            if mixed:
                record[tr("% des aptamères")] = _format_pct(row["aptamer_fraction_pct"])
            record[tr("% des entités moléculaires")] = _format_pct(row["particle_fraction_pct"])
            records.append(record)
        render_readonly_table(pd.DataFrame(records))


def render_tube_analysis(analysis: dict, *, key: str, reference: dict | None = None) -> None:
    """Render original or extended trigger populations using identical semantics."""
    if not analysis:
        st.caption(tr("Analyse des tubes non effectuée pour ce candidat."))
        return
    status = analysis.get("status")
    if status == "failed":
        st.warning(tr("L’analyse des tubes a échoué : {error}", error=tr(str(analysis.get("error", "—")))))
        return
    tubes = analysis.get("tubes", {})
    if status == "cancelled":
        st.caption(tr("Analyse des tubes interrompue. Seuls les calculs terminés sont affichés."))
    if not tubes:
        if status != "cancelled":
            st.caption(tr("Analyse des tubes non effectuée pour ce candidat."))
        return
    conditions = analysis["conditions"]
    st.caption(tr(
        "Aptamère (tube mixte) : {aptamer} µM · Trigger (chaque tube) : {trigger} µM · Volume : {volume} µL · Jusqu’à {size} brins par complexe",
        aptamer=_format_number(conditions["aptamer_concentration_m"] * 1e6),
        trigger=_format_number(conditions["trigger_concentration_m"] * 1e6),
        volume=_format_number(conditions["volume_l"] * 1e6), size=conditions["max_complex_size"],
    ))
    if reference and reference.get("status") == "completed" and not _reference_is_compatible(analysis, reference):
        st.caption(tr("Comparaison avec l’originale indisponible : les conditions de calcul diffèrent."))
        reference = None
    tabs = st.tabs([tr("Trigger seul"), tr("Aptamère + trigger")])
    for tab, tube_key in zip(tabs, ("trigger_only", "aptamer_trigger")):
        with tab:
            tube = tubes.get(tube_key)
            if not tube:
                st.caption(tr("Ce tube n’a pas été calculé."))
                continue
            _render_tube(tube, mixed=tube_key == "aptamer_trigger", key=f"{key}_{tube_key}",
                         delta=_reference_delta(reference, tube_key, tube["summary"]["trigger_parasite_pct"]))
