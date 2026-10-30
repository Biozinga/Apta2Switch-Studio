"""Explore marginal motifs and complete sampled structures without rescoring."""
from __future__ import annotations

from functools import partial
from dataclasses import fields
from types import SimpleNamespace
import json

import altair as alt
import pandas as pd
import streamlit as st

from aptaswitch_core.localization import set_language, tr
from aptaswitch_core.thermodynamics.setup import configured_library_path, get_backend
from aptaswitch_core.thermodynamics.types import ModelProvenance
from .scientific_engines_ui import native_backend_identity, result_engine_label, NATIVE_UNAVAILABLE
from aptaswitch_core.preview import candidate_structure_svg
from aptaswitch_core.structure_ensemble import (
    cluster_structures, ensemble_inputs, ensemble_key, ensemble_path, ensemble_sample_variant_keys, load_ensemble, probability_helices, recorded_native_ensemble_inputs,
)
from .design_session import active_job_running
from .ensemble_runtime import EnsembleJob
from .figure_ui import render_diagram
from .saved_switch_design import design_ensemble_sources
from .table_ui import render_readonly_table


REGION_LABELS = {
    "off_stem": "Tiges OFF", "off_toehold": "Toehold OFF",
    "on_binding": "Duplex ON", "on_rbs_linker": "RBS → fin du linker ON",
    "on_rbs_aug": "Préfixe RBS → fin AUG ON",
}


@st.cache_data(show_spinner=False, max_entries=24)
def _families(state, tolerance):
    return cluster_structures(state, tolerance=tolerance)


@st.fragment(run_every=1.0)
def _monitor():
    set_language(st.session_state.get("ui_language", "en"))
    job = st.session_state.get("ensemble_job")
    if job is None:
        return
    snapshot = job.snapshot()
    if snapshot.state != "running":
        st.rerun(scope="app")
    st.progress(snapshot.current / snapshot.total, text=tr(snapshot.message))
    st.button(tr("Annuler l’analyse des familles"), on_click=job.cancel, key="ensemble_cancel", width="stretch",
              help=tr("Annulation effective à la fin du calcul natif en cours."))


def _position(index, state):
    offset = state["switch_offset"]
    if offset and index < offset:
        return f"T{index + 1}"
    return f"S{index - offset + 1}"


def _render_motifs(state, key):
    threshold = st.slider(tr("Probabilité minimale par paire (%)"), 5, 95, 20, 5, key=f"{key}_threshold",
                          help=tr("T = trigger · S = switch · positions à partir de 1. Une tige regroupe des paires consécutives ; leur probabilité moyenne n’est pas la probabilité de la tige entière. Des motifs concurrents peuvent être incompatibles.")) / 100
    helices = probability_helices(state, threshold=threshold)
    rows = []
    for i, helix in enumerate(helices):
        rows.append({tr("Motif"): i + 1,
                     tr("Appariements"): ", ".join(f"{_position(a, state)}–{_position(b, state)}" for a, b in helix["pairs"]),
                     tr("Paires"): helix["length"], tr("P moyenne (%)"): round(100 * helix["mean_pair_probability"], 1),
                     tr("P minimale (%)"): round(100 * helix["min_pair_probability"], 1),
                     tr("Motifs concurrents"): ", ".join(map(str, helix["competing_motifs"])) or "—"})
    if rows:
        render_readonly_table(pd.DataFrame(rows))
    else:
        st.info(tr("Aucune paire n’atteint ce seuil."))
    if st.checkbox(tr("Afficher la matrice de probabilités"), key=f"{key}_matrix",
                   help=tr("La diagonale donne P(non apparié). Hors diagonale : P(i apparié à j). Les valeurs < 0,001 sont masquées ; la matrice complète est conservée dans le JSON.")):
        matrix = state["pairs"]
        # Zero entries carry no information; omit them to keep the plot small.
        records = [{"i": i+1, "j": j+1, "base_i": _position(i, state), "base_j": _position(j, state), "p": p}
                   for i, row in enumerate(matrix) for j, p in enumerate(row) if p >= 0.001]
        frame = pd.DataFrame(records, columns=["i", "j", "base_i", "base_j", "p"])
        chart = alt.Chart(frame).mark_rect().encode(
            x=alt.X("i:O", title=tr("Position dans le complexe"), axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("j:O", title=tr("Position dans le complexe"), axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color("p:Q", scale=alt.Scale(domain=[0, 1], scheme="blues"), title="P"),
            tooltip=[alt.Tooltip("base_i:N", title="i"), alt.Tooltip("base_j:N", title="j"),
                     alt.Tooltip("p:Q", format=".3f", title="P")],
        ).properties(height=420)
        # Inline data avoid Altair's 5000-row transformer limit for long RNAs.
        with alt.data_transformers.enable("default", max_rows=None):
            st.altair_chart(chart, width="stretch", key=f"{key}_heatmap")


def _render_families(result, candidate, analysis, state_name, key):
    state = analysis["states"][state_name]
    tolerance = st.slider(tr("Tolérance du regroupement (%)"), 0, 20, 5, 1, key=f"{key}_tolerance",
                          help=tr("Regroupement indépendant des poids du score. Les différences hors de ces zones ne séparent pas les familles. Les fréquences dépendent de l’échantillonnage et du seuil ; ce ne sont pas des états cinétiques.")) / 100
    grouped = _families(state, tolerance)
    families = grouped["families"]
    st.caption(tr("{samples} tirages · {families} familles · {bases} bases fonctionnelles comparées. Jusqu’à {changes} bases peuvent changer de partenaire par rapport au motif représentatif. À 0 %, les états d’appariement des zones fonctionnelles sont identiques.",
                  samples=grouped["sample_count"], families=len(families), bases=grouped["selected_count"], changes=grouped["max_changes"]))
    shown = families[:20]
    rows = []
    for i, family in enumerate(shown):
        row = {tr("Famille"): i+1, tr("Fréquence estimée (%)"): round(100 * family["fraction"], 2),
               tr("Tirages"): family["count"], tr("Structures distinctes"): family["distinct_structures"]}
        row.update({tr(REGION_LABELS[name]) + " · " + tr("écart (%)"): round(value * 100, 1)
                    for name, value in family["regional_defects"].items()})
        rows.append(row)
    st.markdown("**" + tr("Familles de structures") + "**",
                help=tr("Écart : part moyenne des bases de la zone dont l’état d’appariement diffère de la cible, parmi les tirages de cette famille."))
    render_readonly_table(pd.DataFrame(rows))
    if len(families) > 20:
        st.caption(tr("Les 20 premières familles couvrent {coverage:.1f} % des tirages. Les autres restent accessibles dans le sélecteur et le JSON.", coverage=100 * sum(f["fraction"] for f in shown)))
    labels = {i: tr("Famille {number} · {frequency:.2f} %", number=i+1, frequency=100*family["fraction"])
              for i, family in enumerate(families)}
    index = st.selectbox(tr("Famille à visualiser"), range(len(families)),
                         format_func=labels.__getitem__,
                         key=f"{key}_family_{int(tolerance*100)}",
                         help=tr("Structure réellement échantillonnée, choisie dans le motif fonctionnel le plus fréquent de la famille. Les autres structures de la famille peuvent différer."))
    representative = families[index]["representative"]
    render_diagram(partial(candidate_structure_svg, result.request.architecture, candidate.trigger_seq,
        candidate.switch_seq, representative, state=state_name,
        reporter_sequence=result.request.reporter.sequence_after_start_rna,
        reporter_label=candidate.reporter_label or result.request.reporter.label,
        width=1000, height=500), structure=representative,
        filename=f"candidate_{candidate.trial}_{state_name}_family_{index+1}", key=f"{key}_representative")
    exported = {"analysis": analysis, "state": state_name, **grouped}
    st.download_button(tr("Exporter les familles (JSON)"), json.dumps(exported, ensure_ascii=False, indent=2),
                       file_name=f"candidate_{candidate.trial}_{state_name}_families.json", mime="application/json",
                       key=f"{key}_families_json", on_click="ignore", width="stretch")


def _ensemble_input_options(result, candidate, num_sample, *, allow_calculation):
    """Search current and archived identities separately; never certify a path-only hit."""
    options = []
    if (allow_calculation and getattr(result.request, "engine", "") == "native-hybrid"
            and native_backend_identity(configured_library_path())):
        options.append(ensemble_inputs(result.request, candidate, num_sample, backend=get_backend()))
    stored = getattr(result, "thermodynamics", None)
    if stored:
        archived = SimpleNamespace(provenance=ModelProvenance(**{field.name: stored[field.name] for field in fields(ModelProvenance)}))
        options.append(ensemble_inputs(result.request, candidate, num_sample, backend=archived))
    if not options:
        options.append(ensemble_inputs(result.request, candidate, num_sample))
    return list({ensemble_key(value): value for value in options}.values())


def render_structure_ensemble(result, candidate, *, key_prefix="ensemble", job_factory=None, description=None,
                              allow_calculation=True):
    """Shared family explorer for design fragments and saved full transcripts."""
    allow_calculation = allow_calculation and getattr(result.request, "engine", "") == "native-hybrid"
    with st.expander(tr("Sous-structures et familles OFF/ON"), expanded=not allow_calculation):
        if description:
            st.caption(description)
        sample_options = (1000, 3000, 5000)
        baseline_options = _ensemble_input_options(result, candidate, 1000, allow_calculation=allow_calculation)
        baseline_key = ensemble_key(baseline_options[0])
        saved_times = {}
        recorded_options = {}
        source_directories = design_ensemble_sources(result)
        paths = tuple(path for directory in source_directories
                      for path in ensemble_path(directory, baseline_key).parent.glob("*.json"))
        for value in (baseline_options if paths else ()):
            variants = ensemble_sample_variant_keys(json.dumps(value, sort_keys=True, allow_nan=False))
            for path in paths:
                count = variants.get(path.stem)
                if count is not None and not path.is_symlink() and path.is_file():
                    saved_times[count] = max(saved_times.get(count, 0), path.stat().st_mtime_ns)
        for path in paths:
            try:
                recorded = recorded_native_ensemble_inputs(path, result.request, candidate)
                if recorded is not None:
                    count = recorded["num_sample"]
                    modified = path.stat().st_mtime_ns
                    saved_times[count] = max(saved_times.get(count, 0), modified)
                    recorded_options.setdefault(count, []).append((modified, recorded))
            except (OSError, ValueError, TypeError, KeyError, RecursionError):
                continue
        if not allow_calculation:
            if not saved_times:
                st.info(tr("Aucune sous-famille enregistrée pour ce candidat dans le calcul de design."))
                return
            sample_options = tuple(sorted(saved_times))
        else:
            sample_options = tuple(sorted(set(sample_options) | set(saved_times)))
        preferred = max(saved_times, key=saved_times.get) if saved_times else 1000
        sample_label = "Structures à échantillonner par état" if allow_calculation else "Échantillonnage enregistré (structures par état)"
        sample_count = st.selectbox(tr(sample_label), sample_options,
                                   index=sample_options.index(preferred), key=f"{key_prefix}_{baseline_key}_samples",
                                   help=tr("La matrice révèle les appariements probables, mais pas leur coexistence. Les familles sont estimées à partir de structures complètes échantillonnées par le moteur enregistré, dans les conditions enregistrées de ce calcul."))
        input_options = [value for _, value in sorted(recorded_options.get(sample_count, []), key=lambda item: item[0], reverse=True)]
        input_options.extend(_ensemble_input_options(result, candidate, sample_count, allow_calculation=allow_calculation))
        input_options = list({ensemble_key(value): value for value in input_options}.values())
        inputs = input_options[0]
        key = key_prefix + "_" + ensemble_key(inputs)
        analysis = None
        read_error = None
        for saved_inputs in input_options:
            for directory in source_directories:
                try:
                    analysis = load_ensemble(directory, saved_inputs)
                except (OSError, ValueError, TypeError, KeyError) as exc:
                    read_error = exc
                    continue
                if analysis is not None:
                    key = key_prefix + "_" + ensemble_key(saved_inputs)
                    break
            if analysis is not None:
                break
        if analysis is None and read_error is not None:
            st.warning(tr("Analyse enregistrée illisible : {error}", error=read_error))
        if allow_calculation:
            busy = active_job_running(st.session_state)
            job = st.session_state.get("ensemble_job")
            if job and job.snapshot().state == "running":
                _monitor()
            elif job and st.session_state.get("ensemble_job_key") == key:
                snapshot = job.snapshot()
                if snapshot.state == "failed":
                    st.error(tr("Analyse des familles impossible : {error}", error=snapshot.error))
                elif snapshot.state == "cancelled":
                    st.info(tr("Analyse des familles annulée."))
        if analysis is None:
            if allow_calculation and not native_backend_identity(configured_library_path()):
                st.info(tr(NATIVE_UNAVAILABLE))
            if not allow_calculation:
                return
            if st.button(tr("Analyser les sous-structures et familles"), key=key_prefix + "_start",
                         disabled=busy or not native_backend_identity(configured_library_path()), type="primary", width="stretch"):
                import_path = configured_library_path()
                if not import_path:
                    st.error(tr(NATIVE_UNAVAILABLE))
                else:
                    try:
                        factory = job_factory or EnsembleJob
                        job = factory(result, candidate, import_path, num_sample=sample_count)
                        job.start()
                    except (OSError, ValueError, RuntimeError) as exc:
                        st.error(tr("Analyse des familles impossible : {error}", error=exc))
                    else:
                        st.session_state.ensemble_job = job
                        st.session_state.ensemble_job_key = key
                        st.rerun()
            return
        thermo = result.request.thermo
        st.caption(tr("Analyse enregistrée · {engine} · {version} · {material} · {temperature} °C · Na⁺ {sodium} M · Mg²⁺ {magnesium} M",
                      engine=result_engine_label(analysis), version=analysis["model"].get("version", "?"), material=analysis["model"].get("model", analysis["model"].get("resolved_material", thermo.material)),
                      temperature=thermo.temperature_c, sodium=thermo.sodium_m, magnesium=thermo.magnesium_m))
        state_name = st.radio(tr("Ensemble à explorer"), ("off", "on"), format_func=str.upper,
                              horizontal=True, key=f"{key}_state",
                              help=tr("OFF : switch seul. ON : complexe trigger + switch déjà formé. Une fréquence ON n’est ni le rendement du tube ni une probabilité d’activation."))
        state_key = f"{key}_{state_name}"
        families_tab, motifs_tab = st.tabs([tr("Familles de structures"), tr("Tiges et matrice de probabilités")])
        with families_tab:
            _render_families(result, candidate, analysis, state_name, state_key)
        with motifs_tab:
            _render_motifs(analysis["states"][state_name], state_key)
        st.download_button(tr("Exporter l’ensemble complet (JSON)"), json.dumps(analysis, ensure_ascii=False),
                           file_name=f"candidate_{candidate.trial}_ensembles.json", mime="application/json",
                           key=f"{key}_json", on_click="ignore", width="stretch")
