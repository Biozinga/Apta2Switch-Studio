"""Read-only summary of a saved candidate's original design measurements."""

from __future__ import annotations

import math
from numbers import Real

import streamlit as st

from aptaswitch_core.localization import tr
from aptaswitch_core.scoring import candidate_score_100
from .defect_ui import FUNCTIONAL_ZONES
from .table_ui import render_readonly_table
from .scientific_engines_ui import result_engine_label
from .accessibility_ui import render_rbs_accessibility
from .scoring_help import (
    ASSOCIATION_WITH_APTAMER_LABEL, ASSOCIATION_WITH_APTAMER_WEIGHT_LABEL,
    association_with_aptamer_help,
)


def _value(value, *, decimals=3):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        return "—"
    return f"{value:.{decimals}f}"


def _table(rows, *, label="Critère", value="Valeur"):
    render_readonly_table([{tr(label): tr(name), tr(value): content} for name, content in rows])


def render_saved_candidate_summary(saved, result, candidate):
    """Display recorded values; never rescore or use the editable full context."""
    request = result.request
    regional = candidate.regional_defects or {}
    with st.container(border=True):
        st.subheader(tr("Métriques du candidat · essai {trial}", trial=candidate.trial), anchor=False)
        score, association, leak, defect = st.columns(4)
        score.metric(tr("Score net (/100)"), _value(candidate_score_100(candidate, request.score_scale), decimals=2),
                     help=tr("La note /100 convertit la somme pondérée avec une échelle fixe, indépendante des autres candidats. Elle ne mesure pas une probabilité de fonctionnement."))
        association.metric(tr("Association ON (%)"), _value(candidate.on_yield_pct, decimals=2),
                           help=tr("Fraction des switches liés au trigger à l’équilibre, dans le tube switch + trigger."))
        leak.metric(tr(ASSOCIATION_WITH_APTAMER_LABEL), _value(candidate.leak_pct, decimals=2),
                    help=association_with_aptamer_help(request.thermo,
                        weight=saved.snapshot["request"].get("scoring", {}).get("leak")))
        defect.metric(tr("Défaut OFF/ON (0–1)"), _value(regional.get("balanced")),
                      help=tr("Compromis des quatre défauts locaux, avec les poids du run. Plus faible est meilleur."))
        energies, defects = st.columns(2)
        with energies:
            st.markdown("**" + tr("Énergies · kcal/mol") + "**", help=tr(
                "Énergies MFE enregistrées. ΔΔG = ΔG ON − ΔG OFF − ΔG du trigger seul ; ce n’est pas une barrière cinétique. ΔG RBS–linker décrit le repliement de cette région isolée."))
            _table([(label, _value(getattr(candidate, field))) for label, field in (
                ("ΔG OFF", "delta_g_off"), ("ΔG ON", "delta_g_on"),
                ("ΔΔG", "ddg_activation"), ("ΔG RBS–linker", "delta_g_rbs_linker"))])
        with defects:
            st.markdown("**" + tr("Défauts locaux · 0–1") + "**", help=tr(
                "Plus faible est meilleur. Les deux régions ON se recouvrent ; leurs défauts ne s’additionnent pas."))
            rows = [(f"{state} · {tr(label)}", _value(regional.get(key))) for key, state, label in FUNCTIONAL_ZONES]
            if request.inverse_extension is not None:
                rows.append(("Trigger seul · défaut de linéarité", _value((candidate.trigger_linearity or {}).get("normalized_defect"))))
            _table(rows)
        st.caption(tr("Défaut global : {total} · écart structural : {structural} · manque de complexes : {concentration} (0–1, informatif).",
                      total=_value(candidate.defect), structural=_value(candidate.defect_structural),
                      concentration=_value(candidate.defect_concentration)))
    from .rbs_calculator_ui import render_rbs_calculator
    render_rbs_calculator(candidate)
    render_rbs_accessibility(candidate)
    _render_run_settings(saved, result, candidate)


def _render_run_settings(saved, result, candidate):
    request = result.request
    thermo = request.thermo
    with st.expander(tr("Conditions du run et réglages du score")):
        st.caption(tr("Run d’origine : {run} · {count} essais demandés.", run=result.run_id, count=request.trials))
        source = saved.snapshot["source"]
        if source.get("started_at"):
            st.caption(tr("Démarré le {date}", date=source["started_at"]))
        st.markdown("**" + tr("Conditions du design") + "**")
        chemistry = {"rna": tr("ARN"), "dna": tr("ADN")}
        _table([
            ("Moteur du résultat", result_engine_label({"thermodynamics": result.thermodynamics,
                                                       "request": {"engine": request.engine}})),
            ("Modèle thermodynamique", thermo.material),
            ("Chimie des brins", tr("Switch : {switch} · trigger : {trigger}", switch=chemistry[thermo.switch_material], trigger=chemistry[thermo.trigger_material])),
            ("Température (°C)", _value(thermo.temperature_c, decimals=2)),
            ("Na⁺ (mM)", _value(thermo.sodium_m * 1000)), ("Mg²⁺ (mM)", _value(thermo.magnesium_m * 1000)),
            ("Switch (µM)", _value(thermo.resolved_switch_concentration_m * 1e6)),
            ("Trigger (µM)", _value(thermo.resolved_trigger_concentration_m * 1e6)),
            ("Rapporteur utilisé pour le design", tr("{name} · {length} nt", name=request.reporter.label, length=len(request.reporter.sequence_after_start_rna))),
            ("Complexes · criblage / validation", tr("{screening} / {validation} brins maximum", screening=thermo.screening_max_complex_size, validation=thermo.validation_max_complex_size)),
            ("Exclusion des STOP prématurés", tr("oui" if request.exclude_stop_candidates else "non")),
            ("Exclusion du RBS–linker non linéaire en ON", tr("oui" if request.exclude_non_linear_candidates else "non")),
        ], label="Paramètre")
        if result.thermodynamics:
            st.json(result.thermodynamics, expanded=False)
        if result.optimization:
            st.markdown("**" + tr("Optimisation enregistrée") + "**")
            trial = next((item for item in result.optimization.get("trials", [])
                          if item.get("trial") == candidate.trial), None)
            if trial and trial.get("selected_objective") and candidate.defect is not None:
                objective = trial["selected_objective"]
                weights = request.functional_defect_weights.fractions().values()
                local = sum(value * weight for value, weight in zip(objective["local_defects"], weights))
                _table([("Défaut local du candidat (%)", _value(100 * local)),
                        ("Défaut global du candidat (%)", _value(100 * objective["global_defect"])),
                        ("Évaluations", _value(trial.get("evaluations"), decimals=0)),
                        ("Candidats contrôlés", _value(trial.get("candidate_screened"), decimals=0))])
            st.json({"algorithm": result.optimization.get("algorithm"),
                     "objective": result.optimization.get("objective"),
                     "functional_weights": result.optimization.get("functional_weights"),
                     **({"rbs_stability": result.optimization["rbs_stability"]}
                        if result.optimization.get("rbs_stability") else {}),
                     "options": result.optimization.get("options", {}),
                     "trial": trial}, expanded=False)
        st.markdown("**" + tr("Poids du classement") + "**")
        weights = saved.snapshot["request"].get("scoring", {})
        terms = [("on_yield", "ON % · bonus"), ("leak", ASSOCIATION_WITH_APTAMER_WEIGHT_LABEL), ("defect", "Défaut OFF/ON · malus"),
                 ("ddg", "ΔΔG · liaison"), ("rbs", "ΔG RBS–linker")]
        if getattr(candidate, "rbs_calculator", {}):
            terms.append(("rbs_calculator", "RBS ON/OFF · bonus"))
        if request.inverse_extension is not None:
            terms.append(("trigger_linearity", "Trigger seul · défaut de linéarité"))
        _table([(label, _value(weights.get(key))) for key, label in terms], value="Poids enregistré")
        st.markdown("**" + tr("Répartition du défaut OFF/ON") + "**")
        fractions = request.functional_defect_weights.fractions()
        render_readonly_table([{tr("Zone"): f"{state} · {tr(label)}",
                               tr("Poids enregistré"): _value(getattr(request.functional_defect_weights, key)),
                               tr("Part du compromis (%)"): _value(100 * fractions[key], decimals=1)}
                              for key, state, label in FUNCTIONAL_ZONES])
        st.caption(tr("Score brut enregistré : {score}. Le score net utilise une conversion fixe sur 100, sans recalcul des mesures arrondies.",
                      score=_value(candidate.score_raw, decimals=4)))
        st.caption(tr("Somme pondérée : + wON × ON − wA × association en présence d’aptamère − 100 × wdéfaut × défaut OFF/ON − wΔΔG × ΔΔG + wRBS × ΔG RBS–linker."))
        if request.inverse_extension is not None:
            st.caption(tr("En extension inverse : − 100 × wlinéarité × défaut du trigger seul."))
        if getattr(candidate, "rbs_calculator", {}):
            st.caption(tr("Indicateur OSTIR : + 100 × wRBS ON/OFF × indice exploratoire."))
