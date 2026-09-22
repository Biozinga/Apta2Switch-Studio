"""NUPACK trigger-extension controls and dedicated results view."""

from __future__ import annotations

from functools import partial
from pathlib import Path

import streamlit as st

from aptaswitch_core.localization import set_language, tr

from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.nupack_setup import configured_nupack_import_path
from aptaswitch_core.sequences import normalize_dna
from aptaswitch_core.extension_preview import extension_preview_svg, aptamer_trigger_linear_svg
from aptaswitch_core.trigger_extension import EXTENSION_SIDES, extension_split
from aptaswitch_studio.extension_runtime import ExtensionRunJob, ExtensionRunRequest
from aptaswitch_studio.concentration_ui import render_concentration_control
from aptaswitch_studio.web_runtime import (
    concentration_to_m, ensure_writable_export_dir, exports_zip_bytes, normalized_export_dir,
)

from aptaswitch_studio.figure_ui import (
    ExportOption, MOLECULAR_VIEW_LABELS, render_diagram, render_export_menu,
    render_svg, render_svg_exports,
)


EXTENSION_DEFAULTS = {
    "extension_side": "5prime",
    "extension_mode": "Deux versions : 24 et 30 nt",
    "extension_target_length": 24,
    "extension_allow_long": False,
    "extension_top": 100,
    "extension_ensemble_cap": 50000,
    "extension_temperature": 37.0,
    "extension_sodium": 120.0,
    "extension_sodium_unit": "mM",
    "extension_magnesium": 0.0,
    "extension_magnesium_unit": "mM",
    "extension_material": "ADN · dna04",
    "extension_job": None,
    "extension_result": None,
    "extension_exports": {},
    "extension_error": "",
    "extension_notice": "",
    "extension_ligand_snapshot": "",
    "extension_ligand_positions": [],
    "system_analysis_signature": None,
    "system_analysis_result": None,
    "system_analysis_error": "",
}


def sync_ligand_selection(aptamer: str) -> None:
    if st.session_state.extension_ligand_snapshot != aptamer:
        for key in list(st.session_state):
            if key.startswith("extension_ligand_base_"):
                del st.session_state[key]
        st.session_state.extension_ligand_positions = []
        st.session_state.extension_ligand_snapshot = aptamer


def _clear_ligand_selection() -> None:
    for key in list(st.session_state):
        if key.startswith("extension_ligand_base_"):
            st.session_state[key] = False
    st.session_state.extension_ligand_positions = []


def render_ligand_selection(aptamer: str) -> None:
    if not aptamer:
        st.caption(tr("Saisissez l’aptamère pour cocher les bases impliquées dans la liaison au ligand."))
        return
    rna = st.session_state.extension_material == "ARN · rna06"
    display = aptamer.replace("T", "U") if rna else aptamer
    st.caption(tr("Cochez les bases connues, de 5′ vers 3′. Les positions commencent à 1 ; plusieurs régions séparées peuvent être sélectionnées."))
    selected = []
    for start in range(0, len(display), 12):
        columns = st.columns(12)
        for offset, base in enumerate(display[start:start + 12]):
            position = start + offset + 1
            key = f"extension_ligand_base_{position}"
            if key not in st.session_state:
                st.session_state[key] = position in st.session_state.extension_ligand_positions
            if columns[offset].checkbox(f"{position} · {base}", key=key):
                selected.append(position)
    st.session_state.extension_ligand_positions = selected
    if selected:
        st.caption(tr("Bases sélectionnées : ") + ", ".join(map(str, selected)))
    else:
        st.caption(tr("Aucune base sélectionnée. L’extension reste possible sans annotation du ligand."))
    st.button(tr("Effacer la sélection du ligand"), on_click=_clear_ligand_selection, disabled=not selected)
    st.caption(tr("Ces annotations repèrent la région de liaison dans les schémas et les résultats ; elles ne modélisent pas le ligand et ne modifient pas le score d’extension."))


def render_extension_preview(aptamer: str, trigger: str) -> None:
    if not (aptamer and trigger):
        return
    selected = st.session_state.extension_ligand_positions
    rna = st.session_state.extension_material == "ARN · rna06"
    lengths = extension_targets(len(trigger), st.session_state.extension_mode, st.session_state.extension_target_length)
    if lengths:
        st.markdown(tr("#### Aperçu en direct de l’extension"))
        for length in lengths:
            st.markdown(tr("**Trigger final · {length} nt**", length=length))
            diagram = partial(extension_preview_svg, aptamer, trigger, length, st.session_state.extension_side, selected, rna=rna)
            render_diagram(
                diagram, caption=tr("N indique une base laissée libre pour l’optimisation NUPACK."),
                filename=f"apercu_extension_{length}nt", key=f"extension_preview_{length}",
                download_labels=(tr("Aperçu de l’extension SVG"), tr("Aperçu de l’extension PNG")),
            )


def current_system_thermo() -> ThermoConditions:
    material, chemistry = {"ADN · dna04": ("dna04", "dna"), "ARN · rna06": ("rna06", "rna")}[st.session_state.extension_material]
    return ThermoConditions(
        material=material, switch_material=chemistry, trigger_material=chemistry,
        temperature_c=float(st.session_state.extension_temperature),
        sodium_m=concentration_to_m(st.session_state.extension_sodium, st.session_state.get("extension_sodium_unit", "mM")),
        magnesium_m=concentration_to_m(st.session_state.extension_magnesium, st.session_state.get("extension_magnesium_unit", "mM")),
    )


def render_system_conditions() -> None:
    """Show experimental conditions using the shared value/unit controls."""
    st.markdown(tr("#### Conditions du système"))
    material_labels = {value: tr(value) for value in ("ADN · dna04", "ARN · rna06")}
    st.selectbox(tr("Modèle thermodynamique"), tuple(material_labels), key="extension_material", format_func=material_labels.__getitem__,
                 help=tr("ADN : U est converti en T. ARN : T est converti en U, pour les deux brins."))
    st.number_input(tr("Température du système (°C)"), min_value=0.0, max_value=100.0, key="extension_temperature")
    with st.container(key="system_salts"):
        salts = st.columns(2)
        with salts[0]:
            render_concentration_control("Na⁺", value_key="extension_sodium", unit_key="extension_sodium_unit", min_m=0.0, max_m=None)
        with salts[1]:
            render_concentration_control("Mg²⁺", value_key="extension_magnesium", unit_key="extension_magnesium_unit", min_m=0.0, max_m=None)
    st.caption(tr("Choisissez des conditions aussi proches que possible des conditions expérimentales ayant démontré l’efficacité du système aptamère–trigger pour votre cible. Elles seront utilisées pour la structure initiale et l’extension éventuelle."))


@st.cache_data(show_spinner=False, max_entries=32)
def _analyze_initial_system(import_path: str, aptamer: str, trigger: str, thermo: ThermoConditions) -> dict:
    from aptaswitch_core.system_analysis import analyze_system
    return analyze_system(import_path, aptamer, trigger, thermo)


def _retry_system_analysis() -> None:
    st.session_state.system_analysis_signature = None


def render_initial_structure(aptamer: str, trigger: str) -> None:
    from aptaswitch_core.extension_visuals import render_extension_structure_svg

    st.markdown(tr("#### Aperçus du système initial"))
    import_path = configured_nupack_import_path()
    thermo = current_system_thermo()
    signature = (aptamer, trigger, import_path, thermo.material, thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m)
    if st.session_state.system_analysis_signature != signature:
        st.session_state.system_analysis_result = None
        st.session_state.system_analysis_error = ""
        st.session_state.aptamer_indices = []
        st.session_state.aptamer_status = ""
        if aptamer and trigger and import_path:
            if _is_running("extension_job") or _is_running("design_job"):
                st.info(tr("La structure initiale sera actualisée à la fin du calcul en cours."))
                return
            with st.spinner(tr("Calcul de la structure initiale avec NUPACK…")):
                try:
                    result = _analyze_initial_system(import_path, aptamer, trigger, thermo)
                except Exception as exc:
                    st.session_state.system_analysis_error = str(exc)
                else:
                    st.session_state.system_analysis_result = result
                    st.session_state.aptamer_indices = list(result["paired_trigger_indices"])
        st.session_state.system_analysis_signature = signature
    if not (aptamer and trigger):
        st.caption(tr("Saisissez l’aptamère et le trigger pour afficher leur structure 2D."))
        return
    if not import_path:
        st.info(tr("NUPACK est requis pour la structure 2D. Configurez-le dans Réglages → Dépendances."))
        return
    if st.session_state.system_analysis_error:
        st.error(tr("La structure initiale n’a pas pu être calculée : ") + st.session_state.system_analysis_error)
        st.button(tr("Réessayer l’analyse du système"), on_click=_retry_system_analysis)
        return
    result = st.session_state.system_analysis_result
    if not result:
        return
    # The color selection never changes the cache key: clicking a base only
    # redraws the figure and does not repeat the folding calculation.
    energy = result["mfe_energy_kcal_mol"]
    conditions = (f"{thermo.material} · {thermo.temperature_c:g} °C · "
                  f"Na⁺ {st.session_state.extension_sodium:g} {st.session_state.extension_sodium_unit} · "
                  f"Mg²⁺ {st.session_state.extension_magnesium:g} {st.session_state.extension_magnesium_unit}")
    diagram = partial(render_extension_structure_svg,
        result["aptamer_sequence"], result["trigger_sequence"], result["mfe_structure"],
        ligand_aptamer_positions=st.session_state.extension_ligand_positions,
        title=tr("Système initial · aptamère + trigger"), subtitle=f"{conditions} · ΔG = {energy:.2f} kcal/mol",
        width=1000, height=500,
    )
    structure_tab, linear_tab = st.tabs([tr(label) for label in MOLECULAR_VIEW_LABELS])
    with structure_tab:
        render_diagram(
            diagram, filename="systeme_initial", key="system_initial", structure=result["mfe_structure"],
            download_labels=(tr("Structure initiale SVG"), tr("Structure initiale PNG")),
            caption=tr("Structure MFE calculée par NUPACK, ensemble stacking. Rose : bases de l’aptamère annotées pour le ligand. Le ligand lui-même n’est pas simulé."),
        )
    with linear_tab:
        linear = partial(aptamer_trigger_linear_svg,
            result["aptamer_sequence"], result["trigger_sequence"],
            ligand_aptamer_positions=st.session_state.extension_ligand_positions,
            title=tr("Séquences du système initial"),
        )
        render_diagram(
            linear, filename="systeme_initial_lineaire", key="system_initial_linear", structure=result["mfe_structure"],
            download_labels=(tr("Séquence initiale SVG"), tr("Séquence initiale PNG")),
        )
    st.caption(tr("{paired}/{total} bases du trigger appariées à l’aptamère · {conditions}", paired=len(result["paired_trigger_indices"]), total=len(trigger), conditions=conditions))
    with st.expander(tr("Séquences du système initial")):
        st.code(result["aptamer_sequence"] + " + " + result["trigger_sequence"], language=None, wrap_lines=True)


def extension_targets(trigger_length: int, mode: str, custom_length: int) -> tuple[int, ...]:
    if mode == "Deux versions : 24 et 30 nt":
        return tuple(length for length in (24, 30) if length > trigger_length)
    return (int(custom_length),) if custom_length > trigger_length else ()


def extension_disabled_for_length(trigger_length: int) -> bool:
    return 24 <= trigger_length <= 30


def _is_running(key: str) -> bool:
    job = st.session_state.get(key)
    return job is not None and job.snapshot().state == "running"


def _start_extension() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.extension_error = ""
    st.session_state.extension_notice = ""
    try:
        if _is_running("extension_job") or _is_running("design_job"):
            raise ValueError(tr("Un calcul est déjà en cours. Attendez sa fin ou annulez-le dans Résultats."))
        aptamer = normalize_dna(st.session_state.aptamer_input, field_name=tr("aptamère"))
        trigger = normalize_dna(st.session_state.trigger_input, field_name="trigger")
        if extension_disabled_for_length(len(trigger)):
            raise ValueError(tr("L’extension est désactivée pour les triggers de 24 à 30 nt. Vous pouvez passer à la conception du switch."))
        sync_ligand_selection(aptamer)
        lengths = extension_targets(len(trigger), st.session_state.extension_mode, st.session_state.extension_target_length)
        if not lengths:
            raise ValueError(tr("Choisissez une longueur finale supérieure à celle du trigger original."))
        if max(lengths) > 30 and not st.session_state.extension_allow_long:
            raise ValueError(tr("La longueur cible est limitée à 30 nt. Activez l'option de dépassement pour aller au-delà."))
        import_path = configured_nupack_import_path()
        if not import_path:
            raise ValueError(tr("Configurez NUPACK dans Réglages → Dépendances pour calculer les extensions."))
        output, fallback = ensure_writable_export_dir(normalized_export_dir(st.session_state.output_dir))
        if fallback:
            st.session_state.extension_notice = tr("Les exports seront enregistrés dans {output}.", output=output)
        request = ExtensionRunRequest(
            aptamer_dna=aptamer, trigger_dna=trigger, target_lengths=lengths,
            molecule_name=st.session_state.molecule_name, import_path=import_path,
            output_dir=output, top=int(st.session_state.extension_top),
            max_ensemble_candidates=int(st.session_state.extension_ensemble_cap),
            extension_side=st.session_state.extension_side,
            ligand_aptamer_positions=tuple(st.session_state.extension_ligand_positions),
            thermo=current_system_thermo(),
        )
        job = ExtensionRunJob(request)
    except Exception as exc:
        st.session_state.extension_error = str(exc)
        return
    st.session_state.extension_job = job
    st.session_state.extension_result = None
    st.session_state.extension_exports = {}
    for key in list(st.session_state):
        if key.startswith("extension_selected_"):
            del st.session_state[key]
    st.session_state.results_kind = "Extension de triggers"
    st.session_state.nav = "Résultats"
    job.start()


def render_extension_controls(aptamer: str, trigger: str) -> None:
    length_disabled = extension_disabled_for_length(len(trigger))
    with st.container(border=True, key="trigger_extension_controls"):
        if length_disabled:
            st.html("""
                <style>
                .st-key-trigger_extension_controls {
                    background: #edf1f5 !important;
                    border-color: #cbd5e1 !important;
                }
                .st-key-trigger_extension_controls h3,
                .st-key-trigger_extension_controls [data-testid="stCaptionContainer"] {
                    color: #64748b;
                }
                .st-key-trigger_extension_controls button:disabled {
                    background: #e2e8f0 !important;
                    color: #64748b !important;
                    border-color: #cbd5e1 !important;
                    opacity: 1;
                }
                </style>
            """)
        st.subheader(tr("2 · Extension du trigger (facultatif)"))
        if length_disabled:
            st.caption(tr("Votre trigger mesure {length} nt. L’extension est désactivée entre 24 et 30 nt inclus ; vous pouvez passer à la conception du switch.", length=len(trigger)))
        else:
            st.caption(tr("Conservez votre trigger tel quel ou choisissez une extension. Le calcul utilisera les séquences, les annotations et les conditions de la section 1."))
        side_labels = {side: tr(label) for side, label in EXTENSION_SIDES.items()}
        st.radio(tr("Extrémité à étendre"), options=list(EXTENSION_SIDES),
                 format_func=side_labels.__getitem__, key="extension_side", horizontal=True,
                 disabled=length_disabled)
        if st.session_state.extension_side == "both":
            st.caption(tr("Répartition égale des bases ajoutées ; si le nombre est impair, une base de plus est ajoutée en 5′."))
        mode_labels = {value: tr(value) for value in ("Deux versions : 24 et 30 nt", "Longueur personnalisée")}
        st.radio(tr("Longueurs finales"), tuple(mode_labels),
                 key="extension_mode", format_func=mode_labels.__getitem__, disabled=length_disabled)
        if st.session_state.extension_mode == "Longueur personnalisée":
            st.checkbox(tr("Autoriser une longueur supérieure à 30 nt"), key="extension_allow_long",
                        disabled=length_disabled)
            if not st.session_state.extension_allow_long:
                st.session_state.extension_target_length = min(30, st.session_state.extension_target_length)
            st.number_input(
                tr("Longueur cible (nt)"), min_value=2,
                max_value=None if st.session_state.extension_allow_long else 30,
                step=1, key="extension_target_length", disabled=length_disabled,
            )
            if not length_disabled and st.session_state.extension_target_length > 30:
                st.warning(tr("Attention : au-delà de 30 nt, le calcul peut prendre beaucoup de temps. Chaque base ajoutée multiplie la recherche par 4."))
        lengths = extension_targets(len(trigger), st.session_state.extension_mode, st.session_state.extension_target_length)
        if not length_disabled and trigger and lengths:
            for length in lengths:
                added = length - len(trigger)
                prefix, suffix = extension_split(added, st.session_state.extension_side)
                count_text = f"{4 ** added:,}".replace(",", " ") if added <= 40 else f"4^{added}"
                st.caption(tr("{length} nt : +{prefix} bases en 5′ et +{suffix} bases en 3′ · {count} extensions à tester.", length=length, prefix=prefix, suffix=suffix, count=count_text))
            if max(lengths) - len(trigger) >= 10:
                st.warning(tr("Cette recherche exhaustive comporte des millions de candidats ou davantage ; elle peut durer plusieurs heures ou jours."))
        elif not length_disabled and trigger:
            st.info(tr("Choisissez une longueur personnalisée supérieure à la longueur actuelle."))
        with st.expander(tr("Paramètres du calcul d’extension")):
            st.number_input(tr("Candidats conservés par longueur"), min_value=1, max_value=1000, step=1,
                            key="extension_top", disabled=length_disabled)
            st.number_input(tr("Maximum de candidats en analyse d’ensemble par longueur"), min_value=1, step=1000,
                            key="extension_ensemble_cap", disabled=length_disabled)
            st.caption(tr("Toutes les extensions passent le criblage MFE. Si cette limite est atteinte, seule une partie des meilleurs ex æquo passe l’analyse d’ensemble ; cette limitation sera indiquée dans les résultats."))
        available = bool(configured_nupack_import_path())
        running = _is_running("extension_job") or _is_running("design_job")
        st.button(tr("Étendre le trigger"), type="primary", width="stretch", on_click=_start_extension,
                  disabled=length_disabled or not (aptamer and trigger and lengths and available) or running)
        if not available:
            st.caption(tr("NUPACK est requis : configurez-le dans Réglages → Dépendances."))
        elif running:
            st.caption(tr("Un calcul est en cours. Son suivi et son annulation sont disponibles dans Résultats."))
        if st.session_state.extension_error:
            st.error(tr(st.session_state.extension_error))
        if not length_disabled:
            render_extension_preview(aptamer, trigger)


def _cancel_extension() -> None:
    job = st.session_state.get("extension_job")
    if job is not None:
        job.cancel()


@st.fragment(run_every=1.0)
def _extension_monitor(format_duration) -> None:
    # A timed fragment can run without executing the main page again.
    set_language(st.session_state.get("ui_language", "en"))
    job = st.session_state.get("extension_job")
    if job is None:
        return
    snapshot = job.snapshot()
    if snapshot.state == "done":
        if job.consume_once():
            st.session_state.extension_result = snapshot.result
            st.session_state.extension_exports = snapshot.exports or {}
            st.rerun(scope="app")
        if snapshot.result.get("status") == "cancelled":
            st.caption(tr("Extension annulée · résultats partiels."))
        # Once complete, the preserved log lives in the collapsed run details.
        return
    columns = st.columns(4)
    columns[0].metric(tr("Moteur"), "NUPACK")
    columns[1].metric(tr("Progression de l’étape"), f"{snapshot.current:,} / {snapshot.total:,}".replace(",", " "))
    columns[2].metric(tr("Temps écoulé"), format_duration(snapshot.elapsed_seconds))
    columns[3].metric(tr("Temps restant de l’étape"), format_duration(snapshot.eta_seconds))
    st.progress(min(1.0, snapshot.current / max(1, snapshot.total)), text=tr(snapshot.message))
    st.markdown(tr("#### Console de l’extension"))
    st.code("\n".join(tr(line) for line in snapshot.log) or tr("Préparation…"), language=None, wrap_lines=True, height=240)
    if snapshot.state == "running":
        st.button(tr("Annuler l’extension"), on_click=_cancel_extension)
        st.caption(tr("Suivi du criblage MFE, puis des probabilités d’ensemble pour chaque longueur. L’estimation porte sur l’étape en cours."))
        return
    st.error(tr(snapshot.error) if snapshot.error else tr("Le calcul a échoué."))
    if job.consume_once():
        st.rerun(scope="app")


def _use_extension(candidate: dict, aptamer: str, molecule_name: str, ligand_aptamer_positions=()) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.trigger_input = candidate["extended_trigger"]
    st.session_state.aptamer_input = aptamer
    st.session_state.molecule_name = molecule_name
    st.session_state.trigger_snapshot = None
    st.session_state.aptamer_indices = []
    st.session_state.aptamer_status = ""
    _clear_ligand_selection()
    st.session_state.extension_ligand_snapshot = normalize_dna(aptamer, field_name=tr("aptamère"))
    st.session_state.extension_ligand_positions = list(ligand_aptamer_positions)
    for position in ligand_aptamer_positions:
        st.session_state[f"extension_ligand_base_{position}"] = True
    st.session_state.nav = "Conception"


def render_extension_results(*, format_duration) -> None:
    st.caption(tr("Extension du trigger · conservation du complexe aptamère–trigger · MFE et probabilités d’ensemble NUPACK"))
    if st.session_state.extension_notice:
        st.info(tr(st.session_state.extension_notice))
    _extension_monitor(format_duration)
    result = st.session_state.get("extension_result")
    if result is None:
        if st.session_state.get("extension_job") is None:
            st.info(tr("Lancez une extension depuis les séquences de la page Conception."))
        return
    _render_completed_extension(result)


def _render_completed_extension(result: dict) -> None:
    # Filled from the engine's public, JSON-serializable result contract.
    from aptaswitch_core.extension_visuals import (
        extension_selection_landscape_svg, reference_pair_probability_svg, render_extension_structure_svg,
    )

    inputs = result["input"]
    aptamer, original = inputs["aptamer_dna"], inputs["trigger_dna"]
    reference = result.get("reference")
    candidates = result.get("candidates", [])
    st.subheader(tr("{name} · extensions", name=result.get("request", {}).get("molecule_name", tr("Système personnalisé"))))
    selected_ligand = inputs.get("ligand_aptamer_positions", ())
    if selected_ligand:
        st.caption(tr("Bases de l’aptamère liées au ligand (cercles roses) : ") + ", ".join(map(str, selected_ligand)))
    model = result.get("model", {})
    stats = result.get("mfe_screen", {})
    summary = st.columns(3)
    summary[0].metric(tr("Extensions criblées"), f"{stats.get('total_candidates_screened', 0):,}".replace(",", " "))
    summary[1].metric(tr("Analyses d’ensemble"), stats.get("ensemble_candidates_analyzed", 0))
    summary[2].metric(tr("Candidats conservés"), len(candidates))
    st.caption(tr("Un score plus bas est meilleur. Ces calculs évaluent la conservation de structure du complexe à deux brins ; ils ne prédisent pas le rendement du switch complet."))
    if reference:
        st.markdown(tr("#### Aperçus · originale et versions étendues"))
        lengths = sorted({candidate["trigger_length"] for candidate in candidates})
        selected = []
        for length in lengths:
            matches = [candidate for candidate in candidates if candidate["trigger_length"] == length]
            candidate_labels = [
                tr("Rang {rank} · {extension} · score {score:.4f}",
                   rank=row["rank_within_length"], extension=_extension_label(row), score=row["final_score"])
                for row in matches
            ]
            index = st.selectbox(
                tr("Candidat à comparer · {length} nt", length=length), list(range(len(matches))),
                format_func=candidate_labels.__getitem__,
                key=f"extension_selected_{length}",
            )
            selected.append(matches[index])
        panels = [(tr("Original · {length} nt", length=len(original)), original, reference, 0)] + [
            (tr("Version {length} nt · rang {rank}", length=c["trigger_length"], rank=c["rank_within_length"]), c["extended_trigger"], c, c["extension_length"])
            for c in selected
        ]
        for column, (title, sequence, data, added) in zip(st.columns(len(panels)), panels):
            with column:
                st.markdown(f"**{title}**")
                energy = data["mfe_energy_kcal_mol"]
                diagram = partial(render_extension_structure_svg,
                    aptamer, sequence, data["mfe_structure"], extension_length=added,
                    extension_3prime_length=len(data.get("extension_3prime", "")),
                    title=title, subtitle=f"ΔG MFE = {energy:.2f} kcal/mol", width=640, height=620,
                    nucleotide_probabilities=data.get("nucleotide_probabilities"),
                    ligand_aptamer_positions=selected_ligand,
                )
                structure_tab, linear_tab = st.tabs([tr(label) for label in MOLECULAR_VIEW_LABELS])
                with structure_tab:
                    render_diagram(
                        diagram, height=430, filename=f"trigger_{len(sequence)}nt", key=f"extension_structure_{added}",
                        structure=data["mfe_structure"],
                        download_labels=("Structure SVG", "Structure PNG"),
                    )
                with linear_tab:
                    linear = partial(aptamer_trigger_linear_svg,
                        aptamer, sequence, extension_length=added,
                        extension_3prime_length=len(data.get("extension_3prime", "")),
                        ligand_aptamer_positions=selected_ligand, title=title,
                    )
                    render_diagram(
                        linear, filename=f"trigger_{len(sequence)}nt_lineaire", key=f"extension_linear_{added}",
                        structure=data["mfe_structure"],
                        download_labels=(tr("Séquence linéaire SVG"), tr("Séquence linéaire PNG")),
                    )
                st.code(sequence, language=None, wrap_lines=True)
                st.caption(f"ΔG MFE = {energy:.2f} kcal/mol")
                if added:
                    st.caption(_extension_label(data))
                    can_use_for_switch = model.get("aptamer_material", "dna") == "dna"
                    st.button(tr("Utiliser pour le design du switch"), on_click=_use_extension,
                              args=(data, aptamer, result.get("request", {}).get("molecule_name", tr("Système personnalisé")), selected_ligand),
                              key=f"extension_use_{added}", width="stretch", disabled=not can_use_for_switch)
                    if not can_use_for_switch:
                        st.caption(tr("Le design de switches utilise actuellement un aptamère ADN. Les résultats ARN restent exportables."))
        for candidate in selected:
            length = candidate["trigger_length"]
            st.markdown(tr("#### Analyse des candidats · {length} nt", length=length))
            rows = [c for c in candidates if c["trigger_length"] == length]
            st.caption(tr("Le graphique montre les 10 % de candidats disponibles ayant les meilleurs scores pour cette longueur (effectif arrondi au supérieur). Le tableau conserve tous les candidats."))
            landscape = extension_selection_landscape_svg(rows)
            render_svg(landscape, height=600)
            render_svg_exports(
                landscape, filename=f"selection_{length}nt", key=f"extension_landscape_{length}",
                download_labels=(tr("Paysage de sélection SVG"), tr("Paysage de sélection PNG")),
            )
            pairs = candidate.get("reference_pairs", [])
            if pairs:
                pair_svg = reference_pair_probability_svg(pairs, candidate_label=tr("Candidat {length} nt · {extension}", length=length, extension=_extension_label(candidate)))
                render_svg(pair_svg, height=600)
                render_svg_exports(
                    pair_svg, filename=f"paires_{length}nt", key=f"extension_pairs_{length}",
                    download_labels=(tr("Probabilités des paires SVG"), tr("Probabilités des paires PNG")),
                )
            else:
                st.info(tr("La référence ne contient aucune paire MFE à comparer."))
    st.markdown(tr("#### Tableau des triggers étendus"))
    if candidates:
        st.table([
            {tr("Rang / longueur"): c["rank_within_length"], tr("Longueur (nt)"): c["trigger_length"], "Score ↓": c["final_score"],
             "Extension 5′": c["extension_5prime"] or "—", "Extension 3′": c.get("extension_3prime", "") or "—",
             tr("Trigger étendu 5′ → 3′"): c["extended_trigger"],
             "ΔG MFE (kcal/mol)": c["mfe_energy_kcal_mol"], tr("RMSE cœur"): c["core_pair_probability_rmse"],
             tr("Perte des paires"): c["reference_pair_probability_loss"], tr("Appariement extension"): c["extension_pair_probability_sum"],
             tr("Extension ↔ aptamère"): c["extension_to_aptamer_probability_sum"],
             tr("Paires MFE perdues"): c["missing_reference_mfe_pairs"], tr("Nouvelles paires du cœur"): c["new_core_mfe_pairs"],
             tr("Structure 2D (dot-bracket)"): c["mfe_structure"]}
            for c in candidates
        ], hide_index=True, width="stretch")
    else:
        st.warning(tr("Aucun candidat final disponible pour ce calcul."))
    exports = st.session_state.get("extension_exports", {})
    if exports:
        archive_path = Path(exports["zip"]) if exports.get("zip") else None
        archive = archive_path.read_bytes() if archive_path and archive_path.is_file() else exports_zip_bytes(exports)
        export_options = [ExportOption(
            tr("Tous les exports d’extension (.zip)"), archive,
            f"{result.get('run_id', 'extension')}_exports.zip", "application/zip", key="extension_export_archive",
        )]
        for kind, path_text in exports.items():
            path = Path(path_text)
            if kind != "zip" and path.is_file():
                export_options.append(ExportOption(kind.upper(), path.read_bytes(), path.name, key=f"extension_export_{kind}"))
        render_export_menu(export_options, label=tr("Exporter les résultats d’extension"))
        st.caption(tr("Copie locale : {path}", path=result.get("output_dir", "")))


def _extension_label(candidate: dict) -> str:
    return f"5′ : {candidate.get('extension_5prime') or '—'} · 3′ : {candidate.get('extension_3prime') or '—'}"
