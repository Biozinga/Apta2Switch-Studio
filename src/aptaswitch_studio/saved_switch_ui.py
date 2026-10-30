"""Retain design candidates and explore their full reporter context."""
from __future__ import annotations

from dataclasses import fields
from functools import partial
import json
from pathlib import Path
import re

import altair as alt
import streamlit as st

from aptaswitch_core.localization import tr
from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.thermodynamics.setup import configured_library_path, get_backend
from .scientific_engines_ui import native_backend_identity, result_engine_label, NATIVE_MODEL_HELP, NATIVE_UNAVAILABLE
from aptaswitch_core.preview import candidate_display_structure, candidate_structure_svg
from aptaswitch_core.saved_switch_diagrams import context_complex_diagram
from aptaswitch_core.saved_switch_analysis import MAX_COMPLEX_SIZE, MAX_SWEEP_POINTS, build_saved_switch_context, validate_analysis_config
from .candidate_summary_ui import render_saved_candidate_summary
from .design_session import active_job_running
from .concentration_ui import render_concentration_control, render_concentration_range
from .ensemble_ui import render_structure_ensemble
from .inverse_extension_ui import render_inverse_extension_candidate
from .figure_ui import MOLECULAR_VIEW_LABELS, render_diagram
from .saved_switch_design import load_saved_design
from .saved_switch_ffs_ui import render_saved_switch_ffs
from .saved_switch_free_md_ui import render_saved_switch_free_md
from .saved_switch_free_md_runtime import active_free_md
from aptaswitch_core.oxdna_free_md import available_free_md
from .saved_switch_ffs_runtime import active_ffs
from .saved_switch_library import candidate_snapshot, list_saved_switches, saved_switch_root, save_switch, SavedSwitchExistsError
from .saved_switch_runtime import active_analysis, configuration_key, load_context, save_context, saved_analyses, start_analysis
from .saved_switch_ensemble_runtime import create_context_ensemble_job, load_context_ensemble_source
from .table_ui import render_readonly_table


SAVED_SWITCH_SECTIONS = ("Candidat", "Thermodynamique", "Cinétique et film")
SAVED_SWITCH_SECTION_LABELS = {
    "Candidat": "Séquences et structures",
    "Thermodynamique": "Contexte du toehold switch avec rapporteur complet",
    "Cinétique et film": "Simulation cinétique et film",
    "MD libre ADN–ADN": "MD libre ADN–ADN",
}


def _root() -> Path:
    return saved_switch_root(st.session_state.get("active_design_path"), export_root=st.session_state.get("output_dir"))


def _open_saved(path=None):
    if path:
        st.session_state.saved_switch_selected = str(path)
    st.session_state.nav = "Toehold switches enregistrés"


def render_save_candidate(result, candidate) -> None:
    st.divider()
    st.subheader(tr("Retenir ce toehold switch"))
    key = f"saved_switch_save_{result.run_id}_{candidate.trial}"
    name = st.text_input(tr("Nom du toehold switch à enregistrer"), value=f"Toehold switch {candidate.trial}", key=key + "_name")
    if candidate.excluded:
        st.warning(tr("Ce candidat est exclu du classement. Son motif d’exclusion sera conservé."))
    pending = st.session_state.get(key + "_pending")
    if pending and pending["name"] != name.strip():
        st.session_state.pop(key + "_pending", None)
        pending = None
    overwrite = False
    if pending:
        st.warning(tr("Ce nom existe déjà. L’écrasement conserve l’ancien dossier dans .backups."))
        overwrite = st.button(tr("Confirmer l’écrasement du toehold switch"), key=key + "_overwrite", disabled=active_job_running(st.session_state))
    save = st.button(tr("Enregistrer ce toehold switch"), type="primary", key=key + "_save")
    if save or overwrite:
        try:
            root = _root()
            existing = list_saved_switches(root)
            # Do not rename a folder while its native worker is still writing.
            for item in existing:
                if active_analysis(item) or active_ffs(item) or active_free_md(item):
                    raise ValueError(tr("Un calcul utilise un toehold switch enregistré. Attendez sa fin avant d’enregistrer ou d’écraser un candidat."))
            original = None
            design_path = st.session_state.get("active_design_path")
            if design_path and Path(result.output_dir).resolve().is_relative_to(Path(design_path).resolve()):
                aptamer = st.session_state.get("design_original_aptamer")
                trigger = st.session_state.get("design_original_trigger")
                if aptamer and trigger:
                    original = {"aptamer": aptamer, "trigger": trigger}
            snapshot = candidate_snapshot(result, candidate, name, original_sequences=original)
            saved = save_switch(root, name, snapshot, overwrite=overwrite, source_output_dir=result.output_dir)
            st.session_state.pop(key + "_pending", None)
            st.session_state[key + "_saved"] = str(saved.path)
            st.success(tr("Toehold switch enregistré : {name}", name=saved.name))
        except SavedSwitchExistsError:
            st.session_state[key + "_pending"] = {"name": name.strip()}
            st.rerun()
        except (ValueError, OSError) as exc:
            st.error(tr(str(exc)))
    if st.session_state.get(key + "_saved"):
        st.button(tr("Ouvrir les toehold switches enregistrés"), key=key + "_open", on_click=_open_saved, args=(st.session_state[key + "_saved"],))


def _bounded_number(label, value, lower, upper, key, *, help=None):
    value = min(upper, max(lower, float(value)))
    return st.number_input(tr(label), value=value, min_value=float(lower), max_value=float(upper), key=key, format="%.3f", help=tr(help) if help else None)


def _sweep(label, axis, config, minimum, maximum, key, *, concentration=False):
    values = config.get(axis, [])
    enabled = st.checkbox(tr(label), value=bool(values), disabled=minimum == maximum, key=key + "_enabled",
                          help=tr("Chaque courbe fait varier un seul paramètre ; les autres restent aux valeurs ci-dessus. Jusqu’à 51 points au total."))
    if not enabled or minimum == maximum:
        return []
    # Default to a modest interval; each point evaluates two equilibrium tubes.
    lo = values[0] if values else minimum
    hi = values[-1] if values else min(maximum, minimum + (.2 if concentration else 40.))
    if concentration:
        bounds, c = st.columns([.78, .22], vertical_alignment="bottom")
        with bounds:
            start, stop = render_concentration_range(label, key=key, start_m=lo, stop_m=hi, min_m=minimum, max_m=maximum)
    else:
        a, b, c = st.columns(3)
        with a:
            start = _bounded_number("Début", lo, minimum, maximum, key + "_min")
        with b:
            stop = _bounded_number("Fin", hi, minimum, maximum, key + "_max")
    with c:
        count = st.number_input(tr("Nombre de points"), min_value=2, max_value=MAX_SWEEP_POINTS, value=max(2, len(values)) if values else 5, key=key + "_count")
    if stop < start:
        raise ValueError(tr("La fin du balayage doit être supérieure ou égale au début."))
    return [start + (stop - start) * index / (int(count) - 1) for index in range(int(count))]


def _context_controls(saved):
    key = "saved_switch_" + saved.snapshot["snapshot_id"]
    drafts = st.session_state.setdefault("saved_switch_context_drafts", {})
    config = dict(drafts.get(saved.snapshot["snapshot_id"], load_context(saved).get("config", {})))
    request = saved.snapshot["request"]
    thermo = request["thermo"]
    backend = get_backend() if native_backend_identity(configured_library_path()) else None
    limits = {"temperature_min_c": 25., "temperature_max_c": 55., "sodium_min_m": .12, "sodium_max_m": 1., "magnesium_max_m": 0.}
    st.caption(tr(NATIVE_MODEL_HELP))
    if backend is None:
        st.info(tr(NATIVE_UNAVAILABLE))
    config = {**config, "max_complex_size": MAX_COMPLEX_SIZE}
    a, b, c = st.columns(3)
    with a:
        config["temperature_c"] = _bounded_number("Température (°C)", config.get("temperature_c", thermo["temperature_c"]), limits["temperature_min_c"], limits["temperature_max_c"], key + "_temperature")
    with b:
        config["sodium_m"] = render_concentration_control(
            "Na⁺", value_key=key + "_sodium", unit_key=key + "_sodium_unit",
            min_m=limits["sodium_min_m"], max_m=limits["sodium_max_m"], default_unit="mM",
            value_m=config.get("sodium_m", thermo["sodium_m"]),
        )
    with c:
        magnesium_disabled = not limits["magnesium_max_m"]
        if magnesium_disabled:
            st.session_state[key + "_magnesium"] = 0.
        config["magnesium_m"] = render_concentration_control(
            "Mg²⁺", value_key=key + "_magnesium", unit_key=key + "_magnesium_unit",
            min_m=0., max_m=limits["magnesium_max_m"] or None, default_unit="mM",
            value_m=config.get("magnesium_m", thermo["magnesium_m"]), disabled=magnesium_disabled,
            help="Mg²⁺ : 0 mM pour ce modèle." if magnesium_disabled else None,
        )
    a, b = st.columns(2)
    for column, field, label in ((a, "switch_concentration_m", "Switch"), (b, "trigger_concentration_m", "Trigger")):
        with column:
            fallback = thermo.get(field) or thermo.get("concentration_m", 5e-6)
            config[field] = render_concentration_control(label, value_key=key + field, unit_key=key + field + "_unit",
                                                        min_m=1e-12, max_m=1., value_m=config.get(field, fallback))
    config["include_aptamer"] = st.checkbox(tr("Ajouter l’aptamère aux deux tubes"), value=config.get("include_aptamer", False), key=key + "_aptamer",
                                          help=tr("OFF : switch sans trigger · ON : switch + trigger."))
    if config["include_aptamer"]:
        config["aptamer_concentration_m"] = render_concentration_control(
            "Aptamère", value_key=key + "_aptamer_concentration", unit_key=key + "_aptamer_concentration_unit",
            min_m=1e-12, max_m=1., value_m=config.get("aptamer_concentration_m", thermo.get("concentration_m", 5e-6)),
        )
    with st.expander(tr("Effet de la température et des sels")):
        config["temperature_values_c"] = _sweep("Balayer la température (°C)", "temperature_values_c", config, limits["temperature_min_c"], limits["temperature_max_c"], key + "_temp_scan")
        config["sodium_values_m"] = _sweep("Balayer Na⁺", "sodium_values_m", config, limits["sodium_min_m"], limits["sodium_max_m"], key + "_na_scan", concentration=True)
        config["magnesium_values_m"] = (_sweep("Balayer Mg²⁺", "magnesium_values_m", config, 0., limits["magnesium_max_m"], key + "_mg_scan", concentration=True)
                                         if limits["magnesium_max_m"] else [])
    context = None
    error = None
    try:
        config = validate_analysis_config(saved.snapshot, config, backend=backend)
        context = build_saved_switch_context(saved.snapshot, config, backend=backend)
    except ValueError as exc:
        error = tr(str(exc))
        st.info(error)
    # Drafts survive section navigation, but other tools use only saved values.
    drafts[saved.snapshot["snapshot_id"]] = dict(config)
    # The detached FFS keeps its own immutable inputs and output directory.
    # Native thermodynamics performs fresh host-memory admission per request.
    busy = active_job_running(st.session_state, ignore_ffs=True)
    a, b = st.columns(2)
    with a:
        if st.button(tr("Enregistrer le contexte"), key=key + "_save_context", disabled=bool(error) or busy or backend is None):
            try:
                save_context(saved, config, backend=backend)
                st.success(tr("Contexte enregistré dans le dossier du toehold switch."))
            except (ValueError, OSError) as exc:
                st.error(tr(str(exc)))
    with b:
        if st.button(tr("Lancer l’analyse thermodynamique"), type="primary", key=key + "_analyze", disabled=bool(error) or busy or backend is None,
                     help=tr("Nouveau calcul : complexes de un à {count} brins, avec leurs ordres cycliques distincts.", count=MAX_COMPLEX_SIZE)):
            _launch_context_analysis(saved, config)
    return config, context


def _launch_context_analysis(saved, config):
    try:
        job, _ = start_analysis(saved, config, configured_library_path() or "")
        st.session_state.saved_switch_analysis_job = job
        st.session_state["saved_switch_analysis_pending_" + saved.snapshot["snapshot_id"]] = job.run_id
        st.rerun()
    except (ValueError, OSError) as exc:
        st.error(tr(str(exc)))


@st.fragment(run_every=2.)
def _analysis_monitor(job):
    snapshot = job.snapshot()
    st.progress(min(1., snapshot.current / max(1, snapshot.total)), text=tr(snapshot.message) or tr("Préparation du calcul…"))
    st.caption(tr("{done}/{total} points · {minutes:.1f} min", done=snapshot.current, total=snapshot.total, minutes=snapshot.elapsed_seconds / 60.))
    if snapshot.state == "running":
        st.button(tr("Annuler l’analyse thermodynamique"), key="saved_switch_cancel_analysis", on_click=job.cancel)
    elif job.consume_once():
        st.rerun()


def _complex_choices(point, state):
    """Keep the OFF/ON reference selected, with every enumerated species accessible."""
    def counts(row):
        return tuple(row.get(key, 0) for key in ("switch_count", "trigger_count", "aptamer_count"))
    target = (1, 0, 0) if state == "off" else (1, 1, 0)
    rows = point[state]["population"]["complexes"]
    primary = next((row for row in rows if counts(row) == target), None)
    if primary:
        first = {**primary, "mfe": primary.get("mfe") or point[state]["mfe"]}
    else:
        first = {"label": "Switch" if state == "off" else "Trigger + Switch", "mfe": point[state]["mfe"]}
    choices = [first, *(row for row in rows if row is not primary and sum(counts(row)) > 0)]
    # Physical identity is the cyclic strand order, never composition alone.
    # Older saved structures still carry their exact strand order in the MFE.
    choices = [{**row, "complex_key": row.get("complex_key") or "|".join(
        row.get("strand_order") or (row.get("mfe") or {}).get("strand_order") or [f"{state}_{index}"])}
        for index, row in enumerate(choices)]
    return choices


def _render_analysis(analysis, saved):
    snapshot = saved.snapshot
    st.caption(result_engine_label(analysis.get("baseline") or analysis))
    if analysis.get("status") == "failed":
        st.error(tr(analysis.get("error", "Le calcul a échoué.")))
    elif analysis.get("status") == "cancelled":
        st.info(tr("Calcul annulé ; les points terminés restent disponibles."))
    baseline = analysis.get("baseline")
    if baseline:
        a, b, c = st.columns(3)
        metrics = baseline["metrics"]
        a.metric(tr("Switch lié au trigger (ON)"), f"{metrics['on_yield_pct']:.2f} %")
        b.metric(tr("Trigger lié au switch (ON)"), f"{metrics['on_trigger_bound_pct']:.2f} %")
        c.metric(tr("Switch monomère (OFF)"), f"{metrics['off_switch_monomer_pct']:.2f} %")
        st.caption(tr("Fractions de brins dans les complexes à l’équilibre. Ces pourcentages ne prédisent pas directement l’expression du rapporteur."))
        for state in ("off", "on"):
            with st.expander(tr("Tube {state} : composition et structure MFE", state=state.upper())):
                rows = [{tr("Complexe"): row["label"], tr("Concentration (µM)"): row["concentration_m"] * 1e6,
                         tr("Fraction des switches (%)"): row.get("switch_fraction_pct"),
                         tr("Fraction des triggers (%)"): row.get("trigger_fraction_pct")}
                        for row in baseline[state]["population"]["complexes"]]
                render_readonly_table(rows, height=420 if len(rows) > 10 else None)
                key = f"saved_context_{snapshot['snapshot_id']}_{analysis['run_id']}_{state}"
                choices = _complex_choices(baseline, state)
                labels = {choice["complex_key"]: choice["label"] for choice in choices}
                selection = st.selectbox(tr("Complexe à visualiser"), list(labels),
                                         format_func=labels.__getitem__, key=key + "_complex")
                choice_index = next(index for index, choice in enumerate(choices) if choice["complex_key"] == selection)
                mfe = choices[choice_index].get("mfe")
                if not mfe:
                    st.info(tr("La structure de ce complexe n’a pas encore été calculée."))
                    if st.button(tr("Calculer les structures de ces complexes"), key=key + "_structures",
                                 disabled=active_job_running(st.session_state, ignore_ffs=True) or not native_backend_identity(configured_library_path())):
                        config = {**analysis["config"], **baseline["conditions"],
                                  "max_complex_size": MAX_COMPLEX_SIZE, "temperature_values_c": [],
                                  "sodium_values_m": [], "magnesium_values_m": []}
                        _launch_context_analysis(saved, config)
                    st.caption(tr("Le nouveau calcul reprend les conditions de référence de ce résultat, avec des complexes jusqu’à {count} brins, sans relancer les balayages.",
                                  count=MAX_COMPLEX_SIZE))
                    continue
                st.caption(tr("Brins de cette MFE : {strands}", strands=" + ".join(mfe["strand_order"])))
                try:
                    diagram = context_complex_diagram(snapshot, analysis["context"], mfe)
                    diagram_key = key + f"_complex_{choice_index}"
                    structure_tab, linear_tab = st.tabs([tr(label) for label in MOLECULAR_VIEW_LABELS])
                    with structure_tab:
                        render_diagram(diagram.structure_svg, height=650, structure=diagram.structure,
                                       filename=f"full_reporter_{state}_complex_{choice_index}", key=diagram_key)
                    with linear_tab:
                        st.caption(tr("Traits bleus : appariements entre brins. Arcs gris : appariements au sein d’un brin. Ils proviennent de la structure MFE de ce complexe."))
                        render_diagram(diagram.linear_svg, structure=diagram.structure,
                                       filename=f"full_reporter_{state}_complex_{choice_index}_linear", key=diagram_key + "_linear",
                                       download_labels=("Séquence linéaire SVG", "Séquence linéaire PNG"))
                except (KeyError, TypeError, ValueError) as exc:
                    st.error(tr('Aperçu indisponible : {exc}', exc=exc))
                st.caption(tr("MFE conditionnelle du complexe, indépendamment de son abondance dans le tube. Énergie : {energy:g} kcal/mol.", energy=mfe["energy_kcal_mol"]))
    for axis, label, scale in (("temperature_c", "Température (°C)", 1.), ("sodium_m", "Na⁺ (mM)", 1000.), ("magnesium_m", "Mg²⁺ (mM)", 1000.)):
        points = analysis.get("sweeps", {}).get(axis, [])
        if not points:
            continue
        rows = [{"x": point["conditions"][axis] * scale, "fraction": point["metrics"][field], "species": tr(name)} for point in points for field, name in (("on_yield_pct", "Switch lié au trigger (ON)"), ("off_switch_monomer_pct", "Switch monomère (OFF)"))]
        chart = alt.Chart(alt.Data(values=rows)).mark_line(point=True).encode(x=alt.X("x:Q", title=tr(label)), y=alt.Y("fraction:Q", title=tr("Fraction (%)"), scale=alt.Scale(domain=[0, 100])), color=alt.Color("species:N", title=None), tooltip=[alt.Tooltip("x:Q", title=tr(label)), alt.Tooltip("fraction:Q", title="%", format=".2f"), alt.Tooltip("species:N", title=tr("Complexe"))]).interactive()
        st.altair_chart(chart, width="stretch")
    with st.expander(tr("Conditions, portée du calcul et provenance")):
        if baseline:
            st.caption(tr("Conditions de ce résultat : {temperature:g} °C · Na⁺ {sodium:g} mM · Mg²⁺ {magnesium:g} mM", temperature=baseline["conditions"]["temperature_c"], sodium=1000 * baseline["conditions"]["sodium_m"], magnesium=1000 * baseline["conditions"]["magnesium_m"]))
            maximum = (analysis.get("config") or {}).get("max_complex_size")
            if type(maximum) is int and maximum > 0:
                st.caption(tr("Complexes explorés dans ce résultat : jusqu’à {count} brins.", count=maximum))
            else:
                st.caption(tr("La taille maximale des complexes n’est pas enregistrée dans ce résultat."))
        for note in analysis.get("notes", []):
            st.write(tr(note))
        st.json(analysis.get("config", {}), expanded=False)
        if baseline:
            st.json(baseline.get("thermodynamics", {}), expanded=False)
        context = analysis.get("context") or {}
        if context.get("switch_sequence") and context.get("trigger_sequence"):
            fasta = f">{saved.name} | full_switch\n{context['switch_sequence']}\n>{saved.name} | trigger\n{context['trigger_sequence']}\n"
            st.download_button(tr("Télécharger le contexte FASTA"), fasta, "full_context.fasta", key="saved_switch_context_fasta_" + analysis.get("run_id", ""))
        st.download_button(tr("Télécharger cette analyse JSON"), json.dumps(analysis, ensure_ascii=False, indent=2), "context_analysis.json", "application/json", key="saved_switch_analysis_export_" + analysis.get("run_id", ""))


def _saved_thermo(saved):
    """Load persisted conditions without rendering or borrowing thermodynamic drafts."""
    config = load_context(saved).get("config", {})
    values = {**saved.snapshot["request"]["thermo"],
              **{name: config[name] for name in ("temperature_c", "sodium_m", "magnesium_m") if name in config}}
    return ThermoConditions(**{field.name: values[field.name] for field in fields(ThermoConditions) if field.name in values})


def _render_candidate(saved):
    snapshot = saved.snapshot
    result, candidate = load_saved_design(saved)
    request = result.request
    st.subheader(tr("Structures 2D du design"), anchor=False)
    has_extension = request.inverse_extension is not None or bool(candidate.inverse_extension_metrics)
    labels = [tr("Structures MFE"), tr("Sous-familles du design")]
    if has_extension:
        labels.append(tr("Comparaison de l’extension"))
    tabs = st.tabs(labels)
    key = "saved_design_" + snapshot["snapshot_id"]
    with tabs[0]:
        structures = {state: structure for state, structure in (("off", candidate.structure_off), ("on", candidate.structure_on)) if structure}
        if not structures:
            st.info(tr("Aucune structure MFE n’a été conservée pour ce candidat."))
        else:
            state = st.radio(tr("État MFE du design"), list(structures), format_func=str.upper,
                             horizontal=True, key=key + "_mfe_state")
            try:
                structure = candidate_display_structure(structures[state])
                render_diagram(
                    partial(candidate_structure_svg, request.architecture, candidate.trigger_seq,
                            candidate.switch_seq, structure, state=state,
                            reporter_sequence=request.reporter.sequence_after_start_rna,
                            reporter_label=candidate.reporter_label or request.reporter.label,
                            width=1000, height=500),
                    structure=structure, filename=f"candidate_{candidate.trial}_{state}_mfe", key=key + "_mfe",
                )
            except (TypeError, ValueError) as exc:
                st.error(tr('Structure calculée indisponible : {exc}', exc=exc))
    with tabs[1]:
        render_structure_ensemble(
            result, candidate, key_prefix=key, allow_calculation=False,
            description=tr("Sous-familles conservées depuis le calcul de design de ce candidat, avec le fragment de rapporteur utilisé à cette étape."),
        )
    if has_extension:
        with tabs[2]:
            render_inverse_extension_candidate(result, candidate, key_prefix=key + "_extension")
    render_saved_candidate_summary(saved, result, candidate)
    with st.expander(tr("Dossier et données du candidat")):
        st.caption(str(saved.path))
        st.caption(tr("Enregistré le {date}", date=saved.saved_at))
        st.download_button(tr("Télécharger les séquences FASTA"), saved.fasta_path.read_bytes(), "sequences.fasta", key="saved_switch_fasta_" + snapshot["snapshot_id"])
        st.download_button(tr("Télécharger le candidat JSON"), saved.manifest_path.read_bytes(), "candidate.json", "application/json", key="saved_switch_json_" + snapshot["snapshot_id"])


def _render_thermodynamics(saved):
    st.subheader(tr("Contexte du toehold switch avec rapporteur complet"))
    job = active_analysis(saved)
    if job:
        st.session_state.saved_switch_analysis_job = job
    try:
        config, _ = _context_controls(saved)
    except (ValueError, OSError) as exc:
        st.error(tr(str(exc)))
        config = {}
    if job:
        _analysis_monitor(job)
    analyses = saved_analyses(saved)
    if analyses:
        selection_key = "saved_switch_analysis_select_" + saved.snapshot["snapshot_id"]
        pending_key = "saved_switch_analysis_pending_" + saved.snapshot["snapshot_id"]
        pending = st.session_state.get(pending_key)
        for index, item in enumerate(analyses):
            if pending is not None and item.get("run_id") == pending:
                st.session_state[selection_key] = index
                st.session_state.pop(pending_key, None)
                break
        selection = st.selectbox(tr("Analyse thermodynamique enregistrée"), list(range(len(analyses))), format_func=lambda index: f"{analyses[index].get('completed_at', '')} · {analyses[index]['status']}", key=selection_key)
        analysis = analyses[selection]
        if config and native_backend_identity(configured_library_path()) and analysis.get("configuration_key") != configuration_key(saved, config):
            st.info(tr("Ce résultat correspond à des conditions précédemment enregistrées, différentes des réglages affichés."))
        _render_analysis(analysis, saved)
        if analysis.get("baseline"):
            try:
                source, candidate = load_context_ensemble_source(saved, analysis)
                description = tr("Familles du contexte complet : {length} nucléotides dans le switch, rapporteur entier inclus. Le calcul utilise les conditions de cette analyse enregistrée.",
                                 length=len(candidate.switch_seq))
                if analysis.get("config", {}).get("include_aptamer"):
                    description += " " + tr("L’aptamère influence les populations des tubes ; ces familles décrivent le switch monomère et le complexe trigger + switch, sans aptamère dans ces deux complexes.")
                render_structure_ensemble(
                    source, candidate,
                    key_prefix=f"context_ensemble_{saved.snapshot['snapshot_id']}_{analysis['run_id']}",
                    job_factory=partial(create_context_ensemble_job, saved, analysis),
                    description=description,
                    allow_calculation=bool(analysis.get("thermodynamics") or analysis.get("baseline", {}).get("thermodynamics")),
                )
            except (OSError, ValueError, TypeError, KeyError) as exc:
                st.error(tr("Analyse des familles impossible : {error}", error=exc))


def _select_section():
    st.session_state.saved_switch_section = st.session_state.saved_switch_section_selector


def render_saved_switches_page():
    st.title(tr("Toehold switches enregistrés"), anchor=False)
    try:
        switches = list_saved_switches(_root())
    except (ValueError, OSError) as exc:
        st.error(tr(str(exc)))
        return
    if not switches:
        st.info(tr("Aucun toehold switch enregistré. Dans Résultats, choisissez un candidat puis cliquez sur « Enregistrer ce toehold switch »."))
        return
    paths = [str(item.path) for item in switches]
    if st.session_state.get("saved_switch_selected") not in paths:
        st.session_state.saved_switch_selected = paths[0]
    names = {}
    for item in switches:
        running = active_analysis(item) or active_ffs(item) or active_free_md(item)
        name = re.sub(r"(?i)^t[- ]?switch(?= \d+$)", "Toehold switch", item.name)
        names[str(item.path)] = name + (" · " + tr("Calcul en cours") if running else "")
    # Keep the selection outside widget state so it survives navigation away.
    st.session_state.saved_switch_selector = st.session_state.saved_switch_selected
    st.selectbox(
        tr("Toehold switch enregistré"), paths, format_func=names.__getitem__,
        key="saved_switch_selector", width="stretch",
        on_change=lambda: _open_saved(st.session_state.saved_switch_selector),
    )
    saved = switches[paths.index(st.session_state.saved_switch_selected)]
    snapshot = saved.snapshot
    if snapshot["candidate"].get("excluded"):
        st.warning(tr(snapshot["candidate"].get("exclusion_reason") or "Candidat exclu du classement."))
    sections = SAVED_SWITCH_SECTIONS + (("MD libre ADN–ADN",) if available_free_md(snapshot) else ())
    section = st.session_state.get("saved_switch_section", SAVED_SWITCH_SECTIONS[0])
    if section not in sections:
        section = SAVED_SWITCH_SECTIONS[0]
    # Widget keys may be discarded when this page is absent. Keep the selection
    # separately, just like the candidate selector; navigation never owns jobs.
    st.session_state.saved_switch_section = section
    st.session_state.saved_switch_section_selector = section
    section_labels = {name: tr(SAVED_SWITCH_SECTION_LABELS[name]) for name in sections}
    st.radio(tr("Section du candidat"), sections, horizontal=True,
             format_func=section_labels.__getitem__, key="saved_switch_section_selector", on_change=_select_section,
             label_visibility="collapsed")
    try:
        if section == "Candidat":
            _render_candidate(saved)
        elif section == "Thermodynamique":
            _render_thermodynamics(saved)
        elif section == "Cinétique et film":
            render_saved_switch_ffs(saved, thermo=_saved_thermo(saved))
        elif section == "MD libre ADN–ADN":
            render_saved_switch_free_md(saved, thermo=_saved_thermo(saved))
    except (ValueError, OSError) as exc:
        st.error(tr(str(exc)))
