"""Native hybrid trigger-extension controls and dedicated results view."""

from __future__ import annotations

from functools import partial
from dataclasses import replace
import csv
import io
import json
from pathlib import Path

import streamlit as st
from .table_ui import render_readonly_table

from aptaswitch_core.localization import set_language, tr

from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.tube_analysis import TubeAnalysisConditions
from aptaswitch_core.thermodynamics.setup import MODEL_ID, configured_library_path, get_backend
from .scientific_engines_ui import native_backend_identity, result_engine_label, NATIVE_MODEL_HELP, NATIVE_UNAVAILABLE
from aptaswitch_core.sequences import normalize_dna
from aptaswitch_core.extension_preview import extension_preview_svg, aptamer_trigger_linear_svg
from aptaswitch_core.trigger_extension import EXTENSION_SIDES, extension_split
from aptaswitch_studio.extension_runtime import ExtensionRunJob, ExtensionRunRequest, request_from_extension_result
from aptaswitch_studio.concentration_ui import render_concentration_control
from aptaswitch_studio.system_tube_runtime import SystemTubeJob
from aptaswitch_studio.tube_analysis_ui import render_tube_analysis
from aptaswitch_studio.design_session import active_job_running, remember_original_sequences, save_active_design_safely
from aptaswitch_studio.web_runtime import (
    concentration_from_m, concentration_to_m, ensure_writable_export_dir, exports_zip_bytes, normalized_export_dir,
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
    "extension_parasite_top": 5,
    "system_aptamer_concentration": 5.0,
    "system_aptamer_concentration_unit": "µM",
    "system_trigger_concentration": 5.0,
    "system_trigger_concentration_unit": "µM",
    "system_tube_volume_ul": 20.0,
    "system_tube_job": None,
    "system_tube_signature": None,
    "extension_temperature": 37.0,
    "extension_sodium": 120.0,
    "extension_sodium_unit": "mM",
    "extension_magnesium": 0.0,
    "extension_magnesium_unit": "mM",
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
    "system_analysis_requires_refresh": False,
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
    st.markdown("**" + tr("Bases liées au ligand") + "**", help="\n\n".join((
        tr("Cochez les bases connues, de 5′ vers 3′. Les positions commencent à 1 ; plusieurs régions séparées peuvent être sélectionnées."),
        tr("Ces annotations repèrent la région de liaison dans les schémas et les résultats ; elles ne modélisent pas le ligand et ne modifient pas le score d’extension."),
    )))
    selected = []
    for start in range(0, len(aptamer), 12):
        columns = st.columns(12)
        for offset, base in enumerate(aptamer[start:start + 12]):
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
    st.button(tr("Effacer la sélection du ligand"), on_click=_clear_ligand_selection, disabled=not selected, width="stretch")


def render_extension_preview(aptamer: str, trigger: str) -> None:
    if not (aptamer and trigger):
        return
    selected = st.session_state.extension_ligand_positions
    lengths = extension_targets(len(trigger), st.session_state.extension_mode, st.session_state.extension_target_length)
    if lengths:
        st.markdown(tr("#### Aperçu en direct de l’extension"))
        for length in lengths:
            st.markdown(tr("**Trigger final · {length} nt**", length=length))
            diagram = partial(extension_preview_svg, aptamer, trigger, length, st.session_state.extension_side, selected)
            render_diagram(
                diagram, caption=tr("X indique une base laissée libre pour l’optimisation du moteur hybride expérimental."),
                filename=f"apercu_extension_{length}nt", key=f"extension_preview_{length}",
                download_labels=(tr("Aperçu de l’extension SVG"), tr("Aperçu de l’extension PNG")),
            )


def current_system_thermo() -> ThermoConditions:
    return ThermoConditions(
        material=MODEL_ID, switch_material="dna", trigger_material="dna",
        temperature_c=float(st.session_state.extension_temperature),
        sodium_m=concentration_to_m(st.session_state.extension_sodium, st.session_state.get("extension_sodium_unit", "mM")),
        magnesium_m=0.0,
    )


def render_system_conditions() -> None:
    """Show experimental conditions using the shared value/unit controls."""
    st.markdown(tr("#### Conditions du système"), help="\n\n".join((
        tr("Choisissez des conditions aussi proches que possible des conditions expérimentales ayant démontré l’efficacité du système aptamère–trigger pour votre cible. Elles seront utilisées pour la structure initiale et l’extension éventuelle."),
        tr(NATIVE_MODEL_HELP),
    )))
    if not 25. <= float(st.session_state.extension_temperature) <= 55.:
        st.session_state.extension_temperature = max(25., min(55., float(st.session_state.extension_temperature)))
        st.info(tr("Les conditions proposées sont ramenées aux limites du modèle."))
    st.number_input(tr("Température du système (°C)"), min_value=25.0, max_value=55.0, key="extension_temperature")
    st.session_state.extension_magnesium = 0.0
    with st.container(key="system_salts"):
        salts = st.columns(2)
        with salts[0]:
            render_concentration_control("Na⁺", value_key="extension_sodium", unit_key="extension_sodium_unit", min_m=.12, max_m=1.0)
        with salts[1]:
            render_concentration_control("Mg²⁺", value_key="extension_magnesium", unit_key="extension_magnesium_unit", min_m=0.0, max_m=None, disabled=True)


@st.cache_data(show_spinner=False, max_entries=32)
def _analyze_initial_system(import_path: str, aptamer: str, trigger: str, thermo: ThermoConditions, backend_identity: str) -> dict:
    from aptaswitch_core.system_analysis import analyze_system
    return analyze_system(None, aptamer, trigger, thermo, backend=get_backend(import_path))


def _retry_system_analysis() -> None:
    st.session_state.system_analysis_signature = None


def current_tube_conditions() -> TubeAnalysisConditions:
    return TubeAnalysisConditions(
        aptamer_concentration_m=concentration_to_m(
            st.session_state.system_aptamer_concentration, st.session_state.system_aptamer_concentration_unit),
        trigger_concentration_m=concentration_to_m(
            st.session_state.system_trigger_concentration, st.session_state.system_trigger_concentration_unit),
        volume_l=float(st.session_state.system_tube_volume_ul) * 1e-6,
        max_complex_size=3,
    )


def _start_system_tubes() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    if active_job_running(st.session_state):
        return
    try:
        signature = (normalize_dna(st.session_state.aptamer_input, field_name=tr("aptamère")),
                     normalize_dna(st.session_state.trigger_input, field_name="trigger"),
                     configured_library_path(), current_system_thermo(), current_tube_conditions(), native_backend_identity(configured_library_path()))
    except ValueError:
        return
    if not all(signature[:3]) or not signature[5]:
        return
    st.session_state.system_tube_signature = signature
    job = SystemTubeJob(signature)
    st.session_state.system_tube_job = job
    job.start()


def _cancel_system_tubes() -> None:
    job = st.session_state.get("system_tube_job")
    if job is not None:
        job.cancel()


@st.fragment(run_every=1.0)
def _system_tube_monitor() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    job = st.session_state.get("system_tube_job")
    if job is None:
        return
    snapshot = job.snapshot()
    if snapshot.state != "running" and job.consume_once():
        save_active_design_safely(st.session_state)
        st.rerun(scope="app")
    if snapshot.state == "running":
        st.progress(min(1.0, snapshot.current / max(1, snapshot.total)), text=tr(snapshot.message))
        st.caption(tr("Temps écoulé : {seconds:.0f} s · deux tubes, complexes jusqu’à 3 brins.", seconds=snapshot.elapsed_seconds))
        with st.expander(tr("Console de l’analyse des tubes"), expanded=True):
            st.code("\n".join(tr(line) for line in snapshot.log), language=None, wrap_lines=True)
        st.button(tr("Annuler l’analyse des tubes"), on_click=_cancel_system_tubes, key="cancel_system_tubes", width="stretch")
        return
    if job.signature != st.session_state.system_tube_signature:
        current = st.session_state.system_tube_signature
        old = job.signature
        same_inputs = (old[:2] == current[:2]
                       and replace(old[3], material=MODEL_ID,
                                   screening_max_complex_size=current[3].screening_max_complex_size,
                                   validation_max_complex_size=current[3].validation_max_complex_size)
                           == replace(current[3], material=MODEL_ID)
                       and all(getattr(old[4], field) == getattr(current[4], field)
                               for field in ("aptamer_concentration_m", "trigger_concentration_m", "volume_l")))
        if not getattr(job, "restored", False) or not same_inputs:
            return
        st.caption(tr("Résultat archivé : conditions et moteur conservés dans l’export."))
    if snapshot.error:
        st.error(tr("L’analyse des tubes a échoué : {error}", error=tr(snapshot.error)))
        return
    if snapshot.result:
        render_tube_analysis(snapshot.result, key="system_tubes")
        rows = [dict(tube=name, **row) for name, tube in snapshot.result.get("tubes", {}).items()
                for row in tube.get("complexes", [])]
        options = [ExportOption(
            "JSON", json.dumps(snapshot.result, ensure_ascii=False, indent=2).encode("utf-8"),
            "initial_system_tubes.json", "application/json", key="system_tubes_json",
        )]
        if rows:
            buffer = io.StringIO()
            writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
            options.append(ExportOption("CSV", buffer.getvalue().encode("utf-8-sig"),
                                        "initial_system_tubes.csv", "text/csv", key="system_tubes_csv"))
        render_export_menu(options, label=tr("Exporter l’analyse des tubes"))
        with st.expander(tr("Console de l’analyse des tubes")):
            st.code("\n".join(tr(line) for line in snapshot.log), language=None, wrap_lines=True)


def render_initial_tubes(aptamer: str, trigger: str) -> None:
    st.markdown(tr("#### Complexes parasites · système initial"),
                help=tr("Deux tubes à l’équilibre : trigger seul et aptamère + trigger, avec des complexes de 1 à 3 brins. L’analyse est lancée à la demande."))
    if aptamer and trigger:
        with st.expander(tr("Concentrations et volume des tubes")):
            columns = st.columns(2)
            with columns[0]:
                render_concentration_control("Concentration initiale d’aptamère",
                    value_key="system_aptamer_concentration", unit_key="system_aptamer_concentration_unit",
                    min_m=1e-12, max_m=None)
            with columns[1]:
                render_concentration_control("Concentration initiale de trigger",
                    value_key="system_trigger_concentration", unit_key="system_trigger_concentration_unit",
                    min_m=1e-12, max_m=None)
            st.number_input(tr("Volume de chaque tube (µL)"), min_value=0.001, step=1.0,
                            format="%.6g", key="system_tube_volume_ul",
                            help=tr("Ces concentrations et ce volume servent aussi au criblage des extensions. Le volume convertit les concentrations en nombres de molécules ; il ne change pas les pourcentages à l’équilibre."))
    signature = (aptamer, trigger, configured_library_path(), current_system_thermo(), current_tube_conditions(), native_backend_identity(configured_library_path()))
    st.session_state.system_tube_signature = signature
    job = st.session_state.get("system_tube_job")
    if job is not None and job.signature != signature and _is_running("system_tube_job"):
        job.cancel()
    if not (aptamer and trigger):
        st.caption(tr("Saisissez l’aptamère et le trigger pour analyser les complexes présents dans les tubes."))
        _system_tube_monitor()
        return
    st.button(tr("Analyser les complexes du système initial"), key="analyze_system_tubes", on_click=_start_system_tubes,
              disabled=not signature[5] or active_job_running(st.session_state), type="primary", width="stretch")
    _system_tube_monitor()


def render_initial_structure(aptamer: str, trigger: str) -> None:
    from aptaswitch_core.extension_visuals import render_extension_structure_svg

    st.markdown(tr("#### Aperçus du système initial"))
    import_path = configured_library_path()
    thermo = current_system_thermo()
    signature = (aptamer, trigger, import_path, thermo.material, thermo.trigger_material, thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m, native_backend_identity(import_path))
    refresh = bool(import_path and st.session_state.get("system_analysis_requires_refresh"))
    previous = st.session_state.system_analysis_signature
    archived_matches = (not signature[-1] and previous and previous[:2] == signature[:2]
                        and previous[4:8] == signature[4:8] and st.session_state.system_analysis_result is not None)
    if (st.session_state.system_analysis_signature != signature or refresh) and not archived_matches:
        st.session_state.system_analysis_result = None
        st.session_state.system_analysis_error = ""
        st.session_state.aptamer_indices = []
        if aptamer and trigger and signature[-1]:
            if active_job_running(st.session_state):
                st.info(tr("La structure initiale sera actualisée à la fin du calcul en cours."))
                return
            with st.spinner(tr("Calcul de la structure initiale avec le moteur hybride natif…")):
                try:
                    if refresh:
                        # A reopened project can originate from another engine
                        # version even when the installation path is unchanged.
                        _analyze_initial_system.clear()
                    result = _analyze_initial_system(import_path, aptamer, trigger, thermo, signature[-1])
                except Exception as exc:
                    st.session_state.system_analysis_error = str(exc)
                else:
                    st.session_state.system_analysis_result = result
                    st.session_state.aptamer_indices = list(result["paired_trigger_indices"])
        st.session_state.system_analysis_signature = signature
        st.session_state.system_analysis_requires_refresh = False
    if not (aptamer and trigger):
        st.caption(tr("Saisissez l’aptamère et le trigger pour afficher leur structure 2D."))
        return
    if not signature[-1] and not st.session_state.system_analysis_result:
        st.info(tr(NATIVE_UNAVAILABLE))
        return
    if st.session_state.system_analysis_error:
        st.error(tr("La structure initiale n’a pas pu être calculée : ") + st.session_state.system_analysis_error)
        st.button(tr("Réessayer l’analyse du système"), on_click=_retry_system_analysis, width="stretch")
        return
    result = st.session_state.system_analysis_result
    if not result:
        return
    if not import_path or not native_backend_identity(import_path):
        st.caption(tr("Résultat archivé : le moteur hybride natif est indisponible pour le recalculer."))
    # The color selection never changes the cache key: clicking a base only
    # redraws the figure and does not repeat the folding calculation.
    energy = result["mfe_energy_kcal_mol"]
    recorded = result.get("model", {})
    if all(key in recorded for key in ("material", "temperature_c", "sodium_m", "magnesium_m")):
        sodium = concentration_from_m(recorded["sodium_m"], st.session_state.extension_sodium_unit)
        magnesium = concentration_from_m(recorded["magnesium_m"], st.session_state.extension_magnesium_unit)
        conditions = (f"{tr('ADN')} · {recorded['material']} · {recorded['temperature_c']:g} °C · "
                      f"Na⁺ {sodium:g} {st.session_state.extension_sodium_unit} · "
                      f"Mg²⁺ {magnesium:g} {st.session_state.extension_magnesium_unit}")
    else:
        conditions = tr("Conditions non renseignées dans le résultat archivé.")
    st.caption(result_engine_label(result))
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
            caption=tr("Structure MFE calculée par {engine}. Rose : bases de l’aptamère annotées pour le ligand. Le ligand lui-même n’est pas simulé.", engine=result_engine_label(result)),
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
    with st.expander(tr("Entrées, provenance et résultat exacts")):
        st.caption(tr("Ordre des brins : aptamère puis trigger, chacun de 5′ vers 3′. Comparez les paires et la notation ci-dessous : l’orientation du dessin 2D peut être différente. Les conditions du système sont distinctes des conditions du design de switch."))
        fasta = f">aptamer\n{result['aptamer_sequence']}\n>trigger\n{result['trigger_sequence']}\n"
        st.json(result.get("thermodynamics", {}), expanded=False)
        st.code(fasta, language=None)
        st.code(result["mfe_structure"], language=None, wrap_lines=True)
        alternatives = result.get("mfe_structures", [])
        if len(alternatives) > 1:
            st.info(tr("Le moteur a retourné {count} structures MFE équivalentes ; la première est affichée et toutes figurent dans l’export.", count=len(alternatives)))
        render_export_menu([
            ExportOption(tr("Entrées et résultat JSON"), json.dumps(result, ensure_ascii=False, indent=2),
                         "systeme_initial.json", "application/json", key="system_initial_json"),
            ExportOption(tr("Séquences FASTA et structure"), fasta + "\n" + result["mfe_structure"] + "\n",
                         "systeme_initial.txt", "text/plain", key="system_initial_text"),
        ])


def extension_targets(trigger_length: int, mode: str, custom_length: int) -> tuple[int, ...]:
    if mode == "Deux versions : 24 et 30 nt":
        return tuple(length for length in (24, 30) if length > trigger_length)
    return (int(custom_length),) if custom_length > trigger_length else ()


def _is_running(key: str) -> bool:
    job = st.session_state.get(key)
    return job is not None and job.snapshot().state == "running"


def _start_extension() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.extension_error = ""
    st.session_state.extension_notice = ""
    try:
        if _is_running("system_tube_job"):
            raise ValueError(tr("L’analyse des tubes est en cours. Suivez-la ou annulez-la dans la section du système initial."))
        if _is_running("saved_switch_analysis_job"):
            raise ValueError(tr("Un calcul de switch enregistré est en cours. Suivez-le ou annulez-le dans l’onglet Toehold switches enregistrés."))
        if active_job_running(st.session_state):
            raise ValueError(tr("Un calcul est déjà en cours. Attendez sa fin ou annulez-le dans Résultats."))
        aptamer = normalize_dna(st.session_state.aptamer_input, field_name=tr("aptamère"))
        trigger = normalize_dna(st.session_state.trigger_input, field_name="trigger")
        sync_ligand_selection(aptamer)
        lengths = extension_targets(len(trigger), st.session_state.extension_mode, st.session_state.extension_target_length)
        if not lengths:
            raise ValueError(tr("Choisissez une longueur finale supérieure à celle du trigger original."))
        if max(lengths) > 30 and not st.session_state.extension_allow_long:
            raise ValueError(tr("La longueur cible est limitée à 30 nt. Activez l'option de dépassement pour aller au-delà."))
        import_path = configured_library_path()
        if not import_path or not native_backend_identity(import_path):
            raise ValueError(tr(NATIVE_UNAVAILABLE))
        active = st.session_state.get("active_design_path")
        requested = Path(active) / "calculations" if active else normalized_export_dir(st.session_state.output_dir)
        output, fallback = ensure_writable_export_dir(requested, fallback=requested) if active else ensure_writable_export_dir(requested)
        if fallback:
            st.session_state.extension_notice = tr("Les exports seront enregistrés dans {output}.", output=output)
        request = ExtensionRunRequest(
            aptamer_dna=aptamer, trigger_dna=trigger, target_lengths=lengths,
            molecule_name=st.session_state.molecule_name, import_path="", thermodynamics_library_path=import_path,
            output_dir=output, top=int(st.session_state.extension_top),
            max_ensemble_candidates=int(st.session_state.extension_ensemble_cap),
            extension_side=st.session_state.extension_side,
            ligand_aptamer_positions=tuple(st.session_state.extension_ligand_positions),
            thermo=current_system_thermo(),
            parasite_top=int(st.session_state.extension_parasite_top),
            tube_conditions=current_tube_conditions(),
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
    with st.container(border=True, key="trigger_extension_controls"):
        st.subheader(tr("Extension exhaustive"),
                     help=tr("Testez toutes les combinaisons de bases ajoutées, puis choisissez un trigger étendu pour concevoir le switch."))
        side_labels = {side: tr(label) for side, label in EXTENSION_SIDES.items()}
        st.radio(tr("Extrémité à étendre"), options=list(EXTENSION_SIDES),
                 format_func=side_labels.__getitem__, key="extension_side", horizontal=True,
                 help=tr("Répartition égale des bases ajoutées ; si le nombre est impair, une base de plus est ajoutée en 5′."))
        mode_labels = {value: tr(value) for value in ("Deux versions : 24 et 30 nt", "Longueur personnalisée")}
        st.radio(tr("Longueurs finales"), tuple(mode_labels),
                 key="extension_mode", format_func=mode_labels.__getitem__)
        if st.session_state.extension_mode == "Longueur personnalisée":
            st.checkbox(tr("Autoriser une longueur supérieure à 30 nt"), key="extension_allow_long")
            if not st.session_state.extension_allow_long:
                st.session_state.extension_target_length = min(30, st.session_state.extension_target_length)
            st.number_input(
                tr("Longueur cible (nt)"), min_value=2,
                max_value=None if st.session_state.extension_allow_long else 30,
                step=1, key="extension_target_length",
            )
            if st.session_state.extension_target_length > 30:
                st.warning(tr("Attention : au-delà de 30 nt, le calcul peut prendre beaucoup de temps. Chaque base ajoutée multiplie la recherche par 4."))
        lengths = extension_targets(len(trigger), st.session_state.extension_mode, st.session_state.extension_target_length)
        if trigger and lengths:
            for length in lengths:
                added = length - len(trigger)
                prefix, suffix = extension_split(added, st.session_state.extension_side)
                count_text = f"{4 ** added:,}".replace(",", " ") if added <= 40 else f"4^{added}"
                st.caption(tr("{length} nt : +{prefix} bases en 5′ et +{suffix} bases en 3′ · {count} extensions à tester.", length=length, prefix=prefix, suffix=suffix, count=count_text))
                if 15 <= added <= 40:
                    st.caption(tr("Soit {billions:.2f} milliards de possibilités. La limite d’analyse d’ensemble ne limite pas ce criblage.", billions=4 ** added / 1e9))
            if max(lengths) - len(trigger) >= 10:
                st.warning(tr("Cette recherche exhaustive comporte des millions de candidats ou davantage ; elle peut durer plusieurs heures ou jours."))
        elif trigger:
            if len(trigger) >= 30 and not st.session_state.extension_allow_long:
                st.info(tr("Pour étendre ce trigger de {length} nt, choisissez une longueur personnalisée supérieure et activez « Autoriser une longueur supérieure à 30 nt ».", length=len(trigger)))
            else:
                st.info(tr("Choisissez une longueur personnalisée supérieure à la longueur actuelle."))
        with st.expander(tr("Paramètres du calcul d’extension")):
            st.number_input(tr("Candidats conservés par longueur"), min_value=1, max_value=1000, step=1,
                            key="extension_top")
            st.number_input(tr("Maximum de candidats en analyse d’ensemble par longueur"), min_value=1, step=1000,
                            key="extension_ensemble_cap",
                            help=tr("Toutes les extensions passent le criblage MFE. Si cette limite est atteinte, seule une partie des meilleurs ex æquo passe l’analyse d’ensemble ; cette limitation sera indiquée dans les résultats."))
            st.number_input(tr("Candidats analysés pour les complexes parasites par longueur"),
                            min_value=0, max_value=1000, step=1, key="extension_parasite_top",
                            help="\n\n".join((
                                tr("Les meilleurs candidats selon le score d’extension passent ensuite l’analyse des deux tubes, jusqu’à 3 brins par complexe. 0 désactive cette étape. Le classement reste fondé sur la conservation de structure."),
                                tr("Par défaut : les 5 meilleurs par longueur, et le trigger original comme référence. Les concentrations et le volume se règlent dans le système initial. Cette étape peut être longue."),
                            )))
        available = bool(native_backend_identity(configured_library_path()))
        running = active_job_running(st.session_state)
        st.button(tr("Étendre le trigger"), type="primary", width="stretch", on_click=_start_extension,
                  disabled=not (aptamer and trigger and lengths and available) or running)
        if not available:
            st.caption(tr(NATIVE_UNAVAILABLE))
        elif _is_running("system_tube_job"):
            st.caption(tr("L’analyse des tubes est en cours. Suivez-la ou annulez-la dans la section du système initial."))
        elif _is_running("saved_switch_analysis_job"):
            st.caption(tr("Un calcul de switch enregistré est en cours. Suivez-le ou annulez-le dans l’onglet Toehold switches enregistrés."))
        elif running:
            st.caption(tr("Un calcul est en cours. Son suivi et son annulation sont disponibles dans Résultats."))
        if st.session_state.extension_error:
            st.error(tr(st.session_state.extension_error))
        render_extension_preview(aptamer, trigger)


def _cancel_extension() -> None:
    job = st.session_state.get("extension_job")
    if job is not None:
        job.cancel()


def _advance_extension() -> None:
    job = st.session_state.get("extension_job")
    if job is not None:
        job.advance_to_ensemble()


def _resume_extension(result: dict) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    try:
        if active_job_running(st.session_state):
            raise ValueError(tr("Un calcul est déjà en cours. Attendez sa fin ou annulez-le dans Résultats."))
        import_path = configured_library_path()
        if not import_path or not native_backend_identity(import_path):
            raise ValueError(tr(NATIVE_UNAVAILABLE))
        active = st.session_state.get("active_design_path")
        output = Path(active) / "calculations" if active else normalized_export_dir(st.session_state.output_dir)
        request = request_from_extension_result(result, import_path="", thermodynamics_library_path=import_path, output_dir=output)
        job = ExtensionRunJob(request, resume_result=result)
        st.session_state.extension_job = job
        st.session_state.extension_error = ""
        st.session_state.extension_notice = tr("Reprise avec les séquences et conditions enregistrées ; le criblage complet ne sera pas recommencé.")
        job.start()
    except Exception as exc:
        st.session_state.extension_error = str(exc)


def _resume_extension_checkpoint(path: Path) -> None:
    try:
        _resume_extension(json.loads(path.read_text(encoding="utf-8")))
    except Exception as exc:
        st.session_state.extension_error = str(exc)


def _render_extension_checkpoints() -> None:
    active = st.session_state.get("active_design_path")
    if not active:
        return
    checkpoints = sorted((Path(active) / "calculations" / "extensions").glob("*/mfe_checkpoint.json"), reverse=True)
    # The final journal is written after the numerical exports. A JSON file
    # alone may be truncated if the process stopped during export.
    checkpoints = [path for path in checkpoints if not (path.parent / f"{path.parent.name}.log").is_file()]
    if not checkpoints:
        return
    with st.expander(tr("Calculs interrompus · sauvegardes MFE")):
        selected = st.selectbox(tr("Sauvegarde à reprendre"), checkpoints,
                                format_func=lambda path: path.parent.name, key="extension_checkpoint")
        st.button(tr("Analyser les candidats de cette sauvegarde"),
                  on_click=_resume_extension_checkpoint, args=(selected,),
                  disabled=active_job_running(st.session_state), width="stretch")


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
            save_active_design_safely(st.session_state)
            st.rerun(scope="app")
        if snapshot.result.get("status") == "cancelled":
            st.caption(tr("Extension annulée · résultats partiels."))
        # Once complete, the preserved log lives in the collapsed run details.
        return
    columns = st.columns(4)
    columns[0].metric(tr("Moteur"), tr("Apta2Switch hybride expérimental"))
    columns[1].metric(tr("Progression de l’étape"), f"{snapshot.current:,} / {snapshot.total:,}".replace(",", " "))
    columns[2].metric(tr("Temps écoulé"), format_duration(snapshot.elapsed_seconds))
    columns[3].metric(tr("Temps restant de l’étape"), format_duration(snapshot.eta_seconds),
                      help=tr("Suivi du criblage MFE, des probabilités d’ensemble, puis des complexes parasites des meilleurs candidats. L’estimation porte sur l’étape en cours."))
    st.progress(min(1.0, snapshot.current / max(1, snapshot.total)), text=tr(snapshot.message))
    st.markdown(tr("#### Console de l’extension"))
    st.code("\n".join(tr(line) for line in snapshot.log) or tr("Préparation…"), language=None, wrap_lines=True, height=240)
    if snapshot.state == "running":
        actions = st.columns(2)
        actions[0].button(tr("Annuler l’extension"), on_click=_cancel_extension, width="stretch")
        actions[1].button(tr("Passer à l’analyse des meilleurs candidats"), on_click=_advance_extension,
                          disabled=not getattr(job, "can_advance", False), key="advance_extension", width="stretch",
                          help=tr("Passer à la suite arrête le criblage MFE et analyse les candidats déjà retenus, dans la limite choisie. La recherche restera indiquée comme partielle."))
        return
    st.error(tr(snapshot.error) if snapshot.error else tr("Le calcul a échoué."))
    if job.consume_once():
        st.rerun(scope="app")


def _use_extension(candidate: dict, aptamer: str, molecule_name: str, ligand_aptamer_positions=()) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    remember_original_sequences(st.session_state)
    st.session_state.design_using_extension = True
    st.session_state.extension_enabled = False
    st.session_state.inverse_extension_enabled = False
    st.session_state.design_section = "switch"
    st.session_state.trigger_input = candidate["extended_trigger"]
    st.session_state.aptamer_input = aptamer
    st.session_state.molecule_name = molecule_name
    st.session_state.trigger_snapshot = None
    st.session_state.aptamer_indices = []
    _clear_ligand_selection()
    st.session_state.extension_ligand_snapshot = normalize_dna(aptamer, field_name=tr("aptamère"))
    st.session_state.extension_ligand_positions = list(ligand_aptamer_positions)
    for position in ligand_aptamer_positions:
        st.session_state[f"extension_ligand_base_{position}"] = True
    st.session_state.nav = "Conception"
    save_active_design_safely(st.session_state)


def render_extension_results(*, format_duration) -> None:
    st.caption(tr("Extension du trigger · conservation du complexe aptamère–trigger · MFE et probabilités d’ensemble"))
    if st.session_state.extension_notice:
        st.info(tr(st.session_state.extension_notice))
    if st.session_state.extension_error:
        st.error(tr(st.session_state.extension_error))
    _extension_monitor(format_duration)
    _render_extension_checkpoints()
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
    st.subheader(tr("{name} · extensions", name=result.get("request", {}).get("molecule_name", tr("Système personnalisé"))),
                 help=tr("Un score plus bas est meilleur. Ces calculs évaluent la conservation de structure du complexe à deux brins ; ils ne prédisent pas le rendement du switch complet."))
    selected_ligand = inputs.get("ligand_aptamer_positions", ())
    if selected_ligand:
        st.caption(tr("Bases de l’aptamère liées au ligand (cercles roses) : ") + ", ".join(map(str, selected_ligand)))
    stats = result.get("mfe_screen", {})
    if (result.get("request") and stats.get("mfe_survivors_for_ensemble", 0)
            and not stats.get("ensemble_complete", False)):
        st.button(tr("Reprendre l’analyse des candidats retenus"),
                  on_click=_resume_extension, args=(result,), key="resume_extension",
                  disabled=active_job_running(st.session_state), type="primary", width="stretch",
                  help=tr("Les séquences et conditions d’origine seront réutilisées. L’analyse d’ensemble sera relancée sur le groupe MFE sauvegardé ; les anciens calculs sans sauvegarde nécessitent une récupération des premiers ex æquo."))
    summary = st.columns(3)
    summary[0].metric(tr("Extensions criblées"), f"{stats.get('total_candidates_screened', 0):,}".replace(",", " "))
    summary[1].metric(tr("Analyses d’ensemble"), stats.get("ensemble_candidates_analyzed", 0))
    summary[2].metric(tr("Candidats conservés"), len(candidates))
    if stats.get("total_candidates_expected"):
        st.caption(tr("Espace de recherche : {count} possibilités · {percent:.4f} % examinées.",
                      count=f"{stats['total_candidates_expected']:,}".replace(",", " "),
                      percent=100 * stats.get("total_candidates_screened", 0) / stats["total_candidates_expected"]))
    for warning in dict.fromkeys(result.get("warnings", [])):
        st.warning(tr(warning))
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
                    st.button(tr("Utiliser pour le design du switch"), on_click=_use_extension,
                              args=(data, aptamer, result.get("request", {}).get("molecule_name", tr("Système personnalisé")), selected_ligand),
                              key=f"extension_use_{added}", width="stretch")
        screening = result.get("parasite_analysis", {})
        if screening and screening.get("status") != "disabled":
            st.markdown(tr("#### Complexes parasites · originale et extensions"))
            st.caption(tr("{completed}/{selected} candidats analysés · {top} meilleurs par longueur demandés. Ces mesures ne modifient pas le score d’extension.",
                          completed=screening.get("completed_count", 0), selected=screening.get("selected_count", 0),
                          top=screening.get("top_per_length", 5)))
            tube_tabs = st.tabs([panel[0] for panel in panels])
            for tab, (_, _, data, added) in zip(tube_tabs, panels):
                with tab:
                    analysis = data.get("tube_analysis")
                    if analysis:
                        render_tube_analysis(analysis, key=f"extension_tubes_{added}",
                                             reference=reference.get("tube_analysis") if added else None)
                    else:
                        st.caption(tr("Ce candidat n’a pas été analysé dans les tubes. Choisissez un candidat parmi les meilleurs rangs analysés."))
        for candidate in selected:
            length = candidate["trigger_length"]
            st.markdown(tr("#### Analyse des candidats · {length} nt", length=length),
                        help=tr("Le graphique montre les 10 % de candidats disponibles ayant les meilleurs scores pour cette longueur (effectif arrondi au supérieur). Le tableau conserve tous les candidats."))
            rows = [c for c in candidates if c["trigger_length"] == length]
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
        render_readonly_table([
            {tr("Rang / longueur"): c["rank_within_length"], tr("Longueur (nt)"): c["trigger_length"], "Score ↓": c["final_score"],
             "Extension 5′": c["extension_5prime"] or "—", "Extension 3′": c.get("extension_3prime", "") or "—",
             tr("Trigger étendu 5′ → 3′"): c["extended_trigger"],
             "ΔG MFE (kcal/mol)": c["mfe_energy_kcal_mol"], tr("RMSE cœur"): c["core_pair_probability_rmse"],
             tr("Perte des paires"): c["reference_pair_probability_loss"], tr("Appariement extension"): c["extension_pair_probability_sum"],
             tr("Extension ↔ aptamère"): c["extension_to_aptamer_probability_sum"],
             tr("Paires MFE perdues"): c["missing_reference_mfe_pairs"], tr("Nouvelles paires du cœur"): c["new_core_mfe_pairs"],
             tr("Structure 2D (dot-bracket)"): c["mfe_structure"],
             **({tr("Parasites · trigger seul (%)"): c.get("trigger_only_parasite_pct"),
                 tr("Parasites · trigger avec aptamère (%)"): c.get("mixed_trigger_parasite_pct"),
                 tr("Trigger lié 1:1 (%)"): c.get("mixed_trigger_target_pct")}
                if result.get("parasite_analysis", {}).get("status") not in (None, "disabled") else {})}
            for c in candidates
        ])
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
