"""Controls and saved-result comparison for extension during switch design."""

from __future__ import annotations

import csv
from functools import partial
import io
import math

import streamlit as st
from aptaswitch_core.thermodynamics.setup import MODEL_ID

from aptaswitch_core.extension_visuals import render_extension_structure_svg
from aptaswitch_core.inverse_extension import ranked_inverse_candidates
from aptaswitch_core.localization import tr
from aptaswitch_core.sequences import normalize_dna
from aptaswitch_core.trigger_extension import EXTENSION_SIDES, extension_split
from aptaswitch_studio.design_session import active_job_running
from aptaswitch_studio.figure_ui import ExportOption, render_diagram, render_export_menu
from .table_ui import render_readonly_table


INVERSE_EXTENSION_DEFAULTS = {
    "inverse_extension_enabled": False,
    "inverse_extension_target_length": 30,
    "inverse_extension_side": "5prime",
    "inverse_extension_top": 20,
    "inverse_extension_performance_weight": 50,
}


def _initialize_controls() -> None:
    for key, value in INVERSE_EXTENSION_DEFAULTS.items():
        st.session_state.setdefault(key, value)


def render_inverse_extension_controls(trigger: str) -> None:
    """Describe the alternative workflow without modifying the original input."""
    _initialize_controls()
    st.subheader(tr("Extension non exhaustive"), help="\n\n".join((
        tr("Le trigger original reste entièrement fixe. Les bases X ajoutées aux extrémités "
           "sont choisies pendant le design du switch. Les meilleurs candidats sont ensuite "
           "comparés au complexe aptamère–trigger original."),
        tr("Comparaison aux conditions du système initial : même modèle ADN/ARN, température, Na⁺ et Mg²⁺."),
        tr("Ce mode n’explore pas toutes les extensions possibles. Une bonne performance "
           "calculée du switch ne garantit pas la conservation du complexe ni son efficacité expérimentale."),
    )))
    columns = st.columns(2)
    with columns[0]:
        st.number_input(
            tr("Longueur finale du trigger (nt)"), min_value=2, max_value=1000,
            step=1, key="inverse_extension_target_length",
        )
    with columns[1]:
        st.number_input(
            tr("Meilleurs candidats à comparer (top X)"), min_value=1, max_value=1000,
            step=1, key="inverse_extension_top",
            help=tr("Les X premiers candidats admissibles sont comparés automatiquement après le design. Les résultats sont enregistrés avec le calcul."),
        )
    st.radio(
        tr("Extrémité à étendre"), options=list(EXTENSION_SIDES),
        format_func=lambda side: tr(EXTENSION_SIDES[side]),
        key="inverse_extension_side", horizontal=True,
    )
    try:
        original = normalize_dna(trigger, field_name="trigger")
    except ValueError:
        st.info(tr("Saisissez un trigger original composé de bases A, C, G, T ou U pour voir les régions X."))
        return
    added = int(st.session_state.inverse_extension_target_length) - len(original)
    if added <= 0:
        st.error(tr("La longueur finale doit dépasser les {length} nt du trigger original.", length=len(original)))
        return
    prefix, suffix = extension_split(added, st.session_state.inverse_extension_side)
    st.code("X" * prefix + original + "X" * suffix, language=None, wrap_lines=True)
    st.caption(tr("Motif 5′ → 3′ : {prefix} X en 5′, {length} bases originales fixes, {suffix} X en 3′.",
                  prefix=prefix, length=len(original), suffix=suffix))


def _number(value: object) -> float | None:
    if type(value) in (int, float) and math.isfinite(value):
        return float(value)
    return None


def _render_model(model: dict) -> None:
    temperature, sodium, magnesium = (_number(model.get(key)) for key in (
        "temperature_c", "sodium_m", "magnesium_m",
    ))
    if all(value is not None for value in (temperature, sodium, magnesium)):
        st.caption(tr(
            "Conditions enregistrées du système initial : {material} · {temperature:g} °C · Na⁺ {sodium:g} M · Mg²⁺ {magnesium:g} M.",
            material=model.get("material", ""), temperature=temperature,
            sodium=sodium, magnesium=magnesium,
        ))


def _render_panel(aptamer: str, trigger: str, data: dict, *, title: str,
                  key: str, added: int = 0, suffix: int = 0) -> None:
    st.markdown(f"**{title}**")
    energy = _number(data.get("mfe_energy_kcal_mol"))
    energy_text = f"ΔG MFE = {energy:.2f} kcal/mol" if energy is not None else ""
    structure = data.get("mfe_structure", "")
    diagram = partial(
        render_extension_structure_svg, aptamer, trigger, structure,
        extension_length=added, extension_3prime_length=suffix,
        title=title, subtitle=energy_text, width=640, height=620,
        nucleotide_probabilities=data.get("nucleotide_probabilities"),
    )
    try:
        render_diagram(diagram, height=430, filename=key, key=key, structure=structure)
    except (ValueError, TypeError, KeyError) as exc:
        st.warning(tr("La structure 2D enregistrée ne peut pas être affichée : {error}", error=str(exc)))
    if trigger:
        st.code(trigger, language=None, wrap_lines=True)
    if energy_text:
        st.caption(energy_text)


def _render_reanalysis_control(result, on_reanalyze) -> None:
    settings = getattr(getattr(result, "request", None), "inverse_extension", None)
    if settings is None or on_reanalyze is None:
        return
    source = getattr(result, "run_id", "")
    if st.session_state.get("inverse_extension_reanalyze_source") != source:
        st.session_state.inverse_extension_reanalyze_source = source
        st.session_state.inverse_extension_reanalyze_top = min(1000, max(1, int(settings.top_candidates)))
    requested = st.number_input(
        tr("Nombre de meilleurs switches à comparer"), min_value=1, max_value=1000,
        step=1, key="inverse_extension_reanalyze_top",
    )
    running = active_job_running(st.session_state)
    recorded_model = getattr(getattr(settings, "thermo", None), "material", MODEL_ID)
    if recorded_model != MODEL_ID:
        st.caption(tr(
            "Nouvelle comparaison : moteur hybride expérimental Vienna + Sugimoto, avec les chimies et les conditions "
            "numériques enregistrées. Le résultat d’origine conserve son modèle {material}.",
            material=recorded_model,
        ))
    st.button(
        tr("Lancer une nouvelle comparaison"), on_click=on_reanalyze,
        args=(result, int(requested)), key="inverse_extension_reanalyze", disabled=running, type="primary", width="stretch",
        help=tr("Compare les switches déjà calculés, avec les séquences et les conditions enregistrées dans ce résultat."),
    )


def _render_trigger_linearity(candidate) -> None:
    """Display saved ensemble measurements; opening results never starts a calculation."""
    measurements = getattr(candidate, "trigger_linearity", {}) or {}
    with st.expander(tr("Trigger isolé · linéarité"), expanded=False):
        unpaired = _number(measurements.get("unpaired_fraction"))
        defect = _number(measurements.get("normalized_defect"))
        if unpaired is None or defect is None:
            st.caption(tr("La linéarité du trigger isolé n’a pas été calculée dans ce résultat."))
            return
        left, right = st.columns(2)
        left.metric(tr("Bases libres du trigger (%)"), f"{100 * unpaired:.2f}", help=tr(
            "La fraction de bases libres est une moyenne sur l’ensemble des conformations du trigger isolé, "
            "au modèle et aux conditions du design du switch. Elle n’est pas la probabilité d’un trigger entièrement linéaire. "
            "Le défaut vaut 1 moins cette fraction et contribue déjà au score du switch ; il n’est pas ajouté au compromis une seconde fois."
        ))
        right.metric(tr("Défaut du trigger isolé (0–1)"), f"{defect:.4f}")
        structure = measurements.get("structure_mfe")
        if structure:
            st.caption(tr("Structure MFE enregistrée du trigger isolé"))
            sequence = getattr(candidate, "trigger_seq", "")
            st.code("\n".join(value for value in (sequence, structure) if value), language=None, wrap_lines=True)


def render_inverse_extension_results(result, *, on_reanalyze=None) -> None:
    """Re-rank saved metrics only; do not calculate folds or change switch scores."""
    analysis = getattr(result, "inverse_extension_analysis", None)
    has_controls = on_reanalyze is not None and getattr(getattr(result, "request", None), "inverse_extension", None) is not None
    if not analysis and not has_controls:
        return
    analysis = analysis if isinstance(analysis, dict) else {}
    _initialize_controls()
    st.subheader(tr("Extension conçue avec le switch · comparaison"), help="\n\n".join((
        tr("Un impact structurel plus bas indique une meilleure conservation du complexe aptamère–trigger original."),
        tr("Le score du switch est affiché sur 100 ; il reste un indice comparatif, pas une probabilité de fonctionnement. "
           "Le classement conserve les différences du score brut, même lorsque les valeurs affichées sont arrondies."),
    )))
    if has_controls:
        if any((getattr(candidate, "inverse_extension_metrics", {}) or {}).get("status") == "completed"
               for candidate in result.candidates):
            with st.expander(tr("Nouvelle comparaison (facultatif)")):
                _render_reanalysis_control(result, on_reanalyze)
        else:
            _render_reanalysis_control(result, on_reanalyze)
    if not analysis:
        st.info(tr("Aucun candidat ne dispose encore d’une comparaison structurale complète."))
        return
    _render_model(analysis.get("model") or {})
    selected = analysis.get("selected_count", 0)
    analyzed = analysis.get("analyzed_count", 0)
    st.caption(tr("Comparaison enregistrée : {analyzed}/{selected} candidats comparés au système original. Le classement ci-dessous porte uniquement sur les candidats analysés.",
                  analyzed=analyzed, selected=selected))
    status = analysis.get("status")
    if status == "cancelled":
        st.warning(tr("Comparaison interrompue : seuls les résultats déjà calculés sont affichés."))
    elif status in {"failed", "partial"}:
        st.warning(tr("La comparaison est incomplète. Une analyse manquante ne signifie pas que la structure est conservée."))
    if analysis.get("error"):
        st.warning(tr("Détail de la comparaison : {error}", error=tr(str(analysis["error"]))))
    for warning in dict.fromkeys(analysis.get("warnings", [])):
        st.warning(tr(warning))
    weight = st.slider(
        tr("Priorité à la performance calculée du switch (%)"),
        min_value=0, max_value=100, step=5, key="inverse_extension_performance_weight",
        help="\n\n".join((
            tr("0 % privilégie la conservation structurale ; 100 % privilégie le score du switch."),
            tr("Le compromis est la moyenne pondérée des rangs normalisés parmi les candidats "
               "comparés, et non des scores bruts. Plus il est bas, meilleur est le compromis. "
               "Les scores du switch restent inchangés ; ce classement ne mesure pas un rendement expérimental."),
        )),
    ) / 100
    ranked = ranked_inverse_candidates(result, performance_weight=weight)
    if not ranked:
        st.info(tr("Aucun candidat ne dispose encore d’une comparaison structurale complète."))
        return
    rows = [{
        tr("Essai"): row["trial"],
        tr("Score du switch"): row["switch_score"],
        tr("Bases libres du trigger (%)"): (
            100 * value if (value := _number((getattr(row["candidate"], "trigger_linearity", {}) or {}).get("unpaired_fraction"))) is not None else None
        ),
        tr("Défaut du trigger isolé (0–1)"): _number(
            (getattr(row["candidate"], "trigger_linearity", {}) or {}).get("normalized_defect"),
        ),
        tr("Impact structurel"): row["conservation_score"],
        tr("Rang switch initial"): row["original_switch_rank"],
        tr("Rang performance (comparés)"): row["performance_rank"],
        tr("Rang structure (comparés)"): row["structure_rank"],
        tr("Compromis"): row["compromise_score"],
    } for row in ranked]
    render_readonly_table(rows)
    buffer = io.StringIO()
    exported_rows = [{**row, tr("Priorité performance (%)"): round(weight * 100)} for row in rows]
    writer = csv.DictWriter(buffer, fieldnames=list(exported_rows[0]))
    writer.writeheader()
    writer.writerows(exported_rows)
    render_export_menu([ExportOption(
        "CSV", buffer.getvalue().encode("utf-8-sig"),
        f"{getattr(result, 'run_id', 'design')}_inverse_extension_comparison.csv",
        "text/csv", key="inverse_extension_comparison_csv",
    )], label=tr("Exporter le classement comparatif"))
    by_trial = {row["trial"]: row for row in ranked}
    key = "inverse_extension_selected_trial"
    if st.session_state.get(key) not in by_trial:
        st.session_state[key] = ranked[0]["trial"]
    labels = {value: tr("Essai {trial} · compromis {score:.3f}",
                       trial=value, score=row["compromise_score"]) for value, row in by_trial.items()}
    trial = st.selectbox(
        tr("Candidat à comparer en 2D"), options=list(by_trial), key=key,
        format_func=labels.__getitem__,
        help=tr("Les deux vues utilisent les séquences et les structures MFE enregistrées lors du calcul, aux mêmes conditions du système initial."),
    )
    _render_trigger_linearity(by_trial[trial]["candidate"])
    _render_comparison_detail(analysis, by_trial[trial]["metrics"], trial, key_prefix="inverse_extension")


def render_inverse_extension_candidate(result, candidate, *, key_prefix: str) -> None:
    """Read the retained comparison without reranking a single candidate as first."""
    analysis = result.inverse_extension_analysis or {}
    metrics = candidate.inverse_extension_metrics or {}
    if metrics.get("status") != "completed":
        st.info(tr("Aucune comparaison d’extension n’a été enregistrée pour ce candidat."))
        if candidate.excluded:
            st.caption(tr("Les candidats exclus du design ne sont pas comparés."))
        elif analysis.get("status") == "completed":
            st.caption(tr("La comparaison du run portait sur ses {count} meilleurs candidats admissibles.",
                          count=analysis.get("selected_count", 0)))
        elif analysis.get("status") in {"failed", "cancelled", "partial"}:
            st.caption(tr("La comparaison du run source n’a pas été terminée."))
        return
    st.caption(tr("Comparaison enregistrée avec ce candidat, aux conditions du système aptamère–trigger original."))
    impact, rank = st.columns(2)
    value = _number(metrics.get("final_score"))
    impact.metric(tr("Impact structurel"), "—" if value is None else f"{value:.4g}",
                  help=tr("Un impact structurel plus bas indique une meilleure conservation du complexe aptamère–trigger original."))
    rank.metric(tr("Rang switch initial"), metrics.get("original_switch_rank", "—"))
    _render_comparison_detail(analysis, metrics, candidate.trial, key_prefix=key_prefix)
    with st.expander(tr("Conditions de la comparaison")):
        _render_model(analysis.get("model") or {})


def _render_comparison_detail(analysis, metrics, trial, *, key_prefix: str) -> None:
    """Share the exact saved structures between results and retained candidates."""
    preservation = st.columns(3)
    for column, label, field in zip(preservation,
                                   ("Paires MFE originales perdues", "Nouvelles paires du noyau", "Appariements impliquant l’extension"),
                                   ("missing_reference_mfe_pairs", "new_core_mfe_pairs", "extension_mfe_pairs")):
        column.metric(tr(label), metrics.get(field, "—"))
    reference = analysis.get("reference") or {}
    aptamer, original = reference.get("aptamer_sequence", ""), reference.get("trigger_sequence", "")
    extended = metrics.get("extended_trigger", "")
    if not aptamer or not original or not extended:
        st.warning(tr("Les séquences de référence ou de l’extension manquent dans ce résultat enregistré."))
        return
    columns = st.columns(2)
    with columns[0]:
        _render_panel(aptamer, original, reference, title=tr("Complexe original"), key=key_prefix + "_original")
    with columns[1]:
        _render_panel(
            aptamer, extended, metrics, title=tr("Complexe avec l’extension · essai {trial}", trial=trial),
            key=f"{key_prefix}_candidate_{trial}",
            added=int(metrics.get("extension_length", len(extended) - len(original))),
            suffix=len(metrics.get("extension_3prime", "")),
        )
