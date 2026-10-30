"""Local Streamlit interface for Apta2Switch-Studio."""

from __future__ import annotations

import base64
import sys
from dataclasses import replace
from functools import partial
from pathlib import Path

# ``streamlit run src/aptaswitch_studio/web_app.py`` executes this file outside
# the package. Add ``src`` in that development-only case; installed entry points
# already resolve the packages normally.
if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st
from aptaswitch_core.localization import LANGUAGES, set_language, tr

from aptaswitch_core.architecture import DEFAULT_TRIGGER_LENGTH, default_architecture_for_trigger
from aptaswitch_core.models import InverseExtensionSettings, DEFAULT_RBS_SEQUENCE, FunctionalDefectWeights, ReporterContext, RunRequest, ScoringWeights, SequenceInput, SwitchArchitecture, ThermoConditions, STOP_FILTER_DISABLED_WARNING
from aptaswitch_core.settings import load_settings, save_settings
from aptaswitch_studio.licenses_ui import render_licenses
from aptaswitch_core.thermodynamics.setup import ENGINE_ID, MODEL_ID, configured_library_path
from aptaswitch_studio.scientific_engines_ui import native_backend_identity, NATIVE_UNAVAILABLE

from aptaswitch_core.preview import (
    DOMAIN_COLORS,
    TRIGGER_COLOR,
    build_switch_model,
    candidate_display_structure,
    candidate_structure_svg,
    linear_sequence_svg,
    switch_preview_svg,
)
from aptaswitch_core.reading import prepare_reading_request
from aptaswitch_core.reporter_library import add_custom_reporter, all_reporters, is_builtin_reporter, remove_custom_reporter
from aptaswitch_core.sequences import SequenceValidationError, normalize_dna
from aptaswitch_core.scoring import candidate_score_100
from aptaswitch_studio.extension_ui import (
    EXTENSION_DEFAULTS,
    current_system_thermo,
    render_extension_results,
    render_ligand_selection,
    render_system_conditions,
    render_initial_structure,
    render_initial_tubes,
    sync_ligand_selection,
)
from aptaswitch_studio.inverse_extension_ui import (
    INVERSE_EXTENSION_DEFAULTS, render_inverse_extension_results,
)
from aptaswitch_studio.design_workflow_ui import (
    render_design_navigation, render_extension_page, select_design_section,
)
from aptaswitch_studio.design_run_ui import render_design_run_controls, render_design_continuation
from aptaswitch_studio.defect_ui import (
    FUNCTIONAL_DEFECT_DEFAULTS, render_candidate_defects,
    render_defect_method, render_regional_defect_controls,
)

from aptaswitch_studio.figure_ui import (
    MOLECULAR_VIEW_LABELS, ExportOption, render_diagram, render_export_menu,
)
from aptaswitch_studio.architecture_ui import render_architecture_selector
from aptaswitch_studio.concentration_ui import (
    render_concentration_control as _render_concentration_control,
)
from aptaswitch_studio.run_library_ui import render_run_library, render_saved_run_details
from aptaswitch_studio.results_table_ui import candidate_status, render_results_table
from aptaswitch_studio.scoring_help import (
    ASSOCIATION_WITH_APTAMER_LABEL, ASSOCIATION_WITH_APTAMER_WEIGHT_LABEL,
    association_with_aptamer_help,
)
from aptaswitch_studio.saved_switch_ui import render_save_candidate, render_saved_switches_page
from aptaswitch_studio.design_library import default_design_root, design_output_dir
from aptaswitch_studio.design_library_ui import persist_active_design, render_design_library, render_design_toolbar
from aptaswitch_studio.design_session import (
    FORM_FIELDS, project_reporters, attach_saved_switch_jobs, active_job_running,
)
from aptaswitch_studio.blast_ui import render_blast_screen
from aptaswitch_studio.web_runtime import (
    DesignRunJob,
    concentration_from_m,
    concentration_to_m,
    default_export_dir,
    ensure_writable_export_dir,
    exports_zip_bytes,
    legacy_export_dir,
    normalized_export_dir,
    preferred_concentration_unit,
)


ASSETS = Path(__file__).resolve().parent / "assets"
LOGO_PATH = ASSETS / "logo.png"
APP_ICON_PATH = ASSETS / "app_icon.png"
SIDEBAR_LOGO_PATH = ASSETS / "igem_SU_logo.png"
NAV_ITEMS = ("Mes designs", "Conception", "Résultats", "Toehold switches enregistrés", "Réglages")
SETTINGS_TABS = ("Licences et remerciements", "Reporters", "Général")
ADD_REPORTER = "__add_reporter__"
STANDARD_ACTIVATOR = default_architecture_for_trigger(DEFAULT_TRIGGER_LENGTH)
DEFAULT_SCORING = ScoringWeights()


def _default_output_dir() -> str:
    configured = load_settings().get("default_output_dir", "")
    if configured and Path(configured).expanduser() == legacy_export_dir():
        configured = ""
    return str(normalized_export_dir(configured))


def _initialize_session() -> None:
    if "ui_language" not in st.session_state:
        saved_language = load_settings().get("ui_language", "en")
        st.session_state.ui_language = saved_language if saved_language in LANGUAGES else "en"
    set_language(st.session_state.ui_language)
    legacy_concentration_value = st.session_state.get("species_concentration_value", 5.0)
    legacy_concentration_unit = st.session_state.get("species_concentration_unit", "µM")
    defaults = {
        "nav": NAV_ITEMS[0],
        "molecule_name": "custom target",
        "aptamer_input": "",
        "trigger_input": "",
        "design_mode": "design",
        "reading_switch_input": "",
        "reading_use_full_reporter": False,
        "trigger_snapshot": None,
        "t7_leader": "GGG",
        "rbs_loop": DEFAULT_RBS_SEQUENCE,
        "linker_auto": True,
        "frame_linker": "",
        "trials": 50,
        "optimization_max_evaluations": 1000,
        "optimization_population_size": 64,
        "optimization_lineages": 8,
        "optimization_local_tolerance": 0.005,
        "exclude_stop_candidates": True,
        "exclude_non_linear_candidates": True,
        "temperature_c": 29.0,
        "trigger_concentration_value": legacy_concentration_value,
        "trigger_concentration_unit": legacy_concentration_unit,
        "switch_concentration_value": legacy_concentration_value,
        "switch_concentration_unit": legacy_concentration_unit,
        "sodium_concentration_value": 120.0,
        "sodium_concentration_unit": "mM",
        "magnesium_concentration_value": 0.0,
        "magnesium_concentration_unit": "mM",
        "reporter_key": "sfgfp",
        "last_reporter_key": st.session_state.get("reporter_key", "sfgfp"),
        "settings_tab": "Licences et remerciements",
        "reporter_return_to_design": False,
        "weight_on": DEFAULT_SCORING.on_yield,
        "weight_leak": DEFAULT_SCORING.leak,
        "weight_defect": DEFAULT_SCORING.defect,
        "weight_ddg": DEFAULT_SCORING.ddg,
        "weight_rbs": DEFAULT_SCORING.rbs,
        "weight_trigger_linearity": DEFAULT_SCORING.trigger_linearity,
        "weight_rbs_calculator": DEFAULT_SCORING.rbs_calculator,
        **FUNCTIONAL_DEFECT_DEFAULTS,
        "output_dir": _default_output_dir(),
        "aptamer_indices": [],
        "run_error": "",
        "run_notice": "",
        "long_run_confirm": False,
        "preview_state": "OFF",
        "extension_enabled": st.session_state.get("inverse_extension_enabled", False),
        "design_section": "system",
        "last_result": None,
        "last_exports": {},
        "design_job": None,
        "results_kind": "Design de switches",
        **EXTENSION_DEFAULTS,
        **INVERSE_EXTENSION_DEFAULTS,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value
    if st.session_state.nav not in NAV_ITEMS:
        st.session_state.nav = NAV_ITEMS[3] if st.session_state.get("saved_switch_selected") else NAV_ITEMS[0]

    # Bring an already open browser session to the common app-local default.
    # Custom locations are kept, and old runs stay discoverable in Results.
    if st.session_state.get("run_storage_revision", 0) < 1:
        if Path(st.session_state.output_dir).expanduser() == legacy_export_dir():
            st.session_state.output_dir = str(default_export_dir())
        if st.session_state.get("settings_default_output") == str(legacy_export_dir()):
            st.session_state.settings_default_output = str(default_export_dir())
        st.session_state.run_storage_revision = 1

    # Update the previous default once, while preserving user-supplied RBSs.
    if st.session_state.get("rbs_default_revision", 0) < 1:
        if st.session_state.rbs_loop == "UAGAGGAGAAC":
            st.session_state.rbs_loop = DEFAULT_RBS_SEQUENCE
        st.session_state.rbs_default_revision = 1

    # Existing weights are user choices; new defaults only fill missing keys.

    # Retire demonstration runs already present in an open browser session.
    old_job = st.session_state.get("design_job")
    if old_job is not None and old_job.request.engine != ENGINE_ID:
        old_job.cancel()
        st.session_state.design_job = None
    old_result = st.session_state.get("last_result")
    if old_result is not None and old_result.request.engine != ENGINE_ID:
        st.session_state.last_result = None
        st.session_state.last_exports = {}
    old_kinetic_job = st.session_state.pop("ms_job", None)
    if old_kinetic_job is not None:
        old_kinetic_job.cancel()

    # Streamlit normally deletes a widget's key when that widget disappears.
    # Touch the page-scoped values on every rerun so navigating to another page
    # does not erase a completed design form.
    extra_keys = {
        "settings_default_output",
        "selected_trial",
        "candidate_state",
        "confirm_reset_settings",
        "saved_run_path",
        "run_open_path",
        "ui_language",
        "new_reporter_label",
        "new_reporter_sequence",
        "blast_lookup_query",
        "blast_lookup_selected",
        "blast_scopes",
        "blast_max_hits",
        "blast_consent",
    }
    dynamic_prefixes = ("extension_selected_", "extension_ligand_base_", "blast_alignment_")
    for key in list(st.session_state):
        if key in defaults or key in FORM_FIELDS or key in extra_keys or key.startswith(dynamic_prefixes):
            st.session_state[key] = st.session_state[key]


def _inject_styles() -> None:
    st.markdown(
        """
        <style>
        .stApp {
            background:
                radial-gradient(circle at 88% -10%, rgba(16,167,200,.12), transparent 28rem),
                linear-gradient(180deg, #f8fbff 0%, #f3f6fa 100%);
        }
        /* Keep only the sidebar toggle from Streamlit's built-in header. */
        [data-testid="stHeader"] {
            visibility:hidden; height:0; min-height:0;
            background:transparent; pointer-events:none;
        }
        [data-testid="stToolbar"] { pointer-events:none; }
        /* Keep chart data/download/fullscreen actions in the app's own UI.
           Scientific hover tooltips remain available on the plot itself. */
        [data-testid="stElementContainer"]:has([data-testid="stVegaLiteChart"])
        [data-testid="stElementToolbar"] { display:none !important; }
        [data-testid="stExpandSidebarButton"] {
            visibility:visible; pointer-events:auto;
            position:fixed; top:.5rem; left:.5rem;
            background:#f8fbff; border-radius:.5rem;
        }
        [data-testid="stSidebar"] { background:#0d1f3c; }
        [data-testid="stSidebar"] * { color:#edf5ff; }
        [data-testid="stSidebar"] [data-testid="stRadio"] label {
            padding:.35rem .45rem; border-radius:.55rem;
        }
        .block-container, [data-testid="stMainBlockContainer"] {
            max-width:1440px; padding-top:1.4rem; padding-bottom:4rem;
        }
        h1, h2, h3, h4 { color:#0d1f3c; letter-spacing:-.025em; }
        h1 { font-size:1.9rem !important; padding-top:.35rem !important; }
        h2 { font-size:1.5rem !important; }
        h3 { font-size:1.25rem !important; }
        h4 { font-size:1.08rem !important; }
        /* Hide only the generated heading permalinks, including in dialogs.
           Help buttons and real documentation links remain available. */
        :is(h1,h2,h3,h4,h5,h6) a[href^="#"] { display:none !important; }
        .apta-brand { display:flex; align-items:center; gap:.7rem; margin:0 0 .4rem; }
        .apta-brand img { flex-shrink:0; }
        .apta-brand-name { color:#243f61; font-size:1.15rem; font-weight:700; }
        .st-key-active_design_banner {
            background:#eaf2ff; border:1px solid #d0e0f6; border-left:4px solid #1258dc;
            border-radius:.8rem; padding:1rem 1.1rem; margin:.2rem 0 .5rem;
        }
        .active-design-label { color:#45618a; font-size:.75rem; font-weight:700;
            letter-spacing:.05em; text-transform:uppercase; margin-bottom:.2rem; }
        .active-design-name { color:#123a70; font-size:1.3rem; font-weight:700;
            line-height:1.35; overflow-wrap:anywhere; white-space:normal; }
        div[data-testid="stVerticalBlockBorderWrapper"] {
            background:rgba(255,255,255,.88); border:1px solid #dbe5f0;
            border-radius:1rem; box-shadow:0 8px 24px rgba(13,31,60,.045);
        }
        .stButton > button, .stDownloadButton > button, .stFormSubmitButton > button,
        .stLinkButton > a, [data-testid="stPopoverButton"] {
            border-radius:.65rem; font-weight:650;
        }
        .stButton > button[kind="primary"] { background:#1258dc; border-color:#1258dc; }
        .legend-chip { display:inline-flex; align-items:center; gap:.3rem; margin:.15rem .7rem .15rem 0;
            color:#425269; font-size:.82rem; }
        .legend-dot { width:.65rem; height:.65rem; border-radius:50%; display:inline-block; }
        code { white-space:pre-wrap !important; overflow-wrap:anywhere; }
        @media (max-width:760px) {
            .st-key-active_design_banner > [data-testid="stHorizontalBlock"] {
                flex-wrap:wrap; gap:.75rem;
            }
            .st-key-active_design_banner > [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {
                width:auto !important; flex:1 1 calc(50% - .75rem) !important; min-width:0 !important;
            }
            .st-key-active_design_banner > [data-testid="stHorizontalBlock"] > [data-testid="stColumn"]:first-child {
                width:100% !important; flex:1 1 100% !important;
            }
        }
        @media (max-width:1050px) {
            .st-key-system_salts > [data-testid="stHorizontalBlock"] {
                flex-direction:column; gap:.75rem;
            }
            .st-key-system_salts > [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {
                width:100%; flex:1 1 100%; min-width:0;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def _render_brand_image(path: Path, *, width: int, alt: str, backdrop: bool = False) -> None:
    # st.image(width=...) downsamples the actual bitmap to CSS pixels. Embed
    # the original PNG so Retina displays retain the source's full resolution.
    from html import escape

    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    card = "background:#fff;padding:8px;border-radius:16px;" if backdrop else "background:transparent;"
    st.html(
        f'<div style="width:fit-content;max-width:100%;{card}">'
        f'<img src="data:image/png;base64,{encoded}" alt="{escape(alt, quote=True)}" '
        f'width="{width}" style="display:block;width:{width}px;max-width:100%;height:auto;background:transparent;object-fit:contain;" />'
        '</div>'
    )


def _render_header() -> None:
    logo = ""
    if LOGO_PATH.exists():
        encoded = base64.b64encode(LOGO_PATH.read_bytes()).decode("ascii")
        logo = (f'<img src="data:image/png;base64,{encoded}" alt="Logo Apta2Switch-Studio" '
                'width="42" height="42" style="object-fit:contain" />')
    st.html(f'<div class="apta-brand">{logo}<span class="apta-brand-name">Δpta2Switch-Studio</span></div>')


def _render_sidebar() -> str:
    with st.sidebar:
        if SIDEBAR_LOGO_PATH.exists():
            _render_brand_image(SIDEBAR_LOGO_PATH, width=120, alt="Logo iGEM Sorbonne Université", backdrop=True)
        st.markdown("### Δpta2Switch")
        if st.session_state.get("active_design_name"):
            st.caption(st.session_state.active_design_name)
        pages = NAV_ITEMS if st.session_state.get("active_design_path") else ("Mes designs", "Réglages")
        page = st.radio(tr('Navigation'), pages, key="nav", label_visibility="collapsed", format_func={label: tr(label) for label in pages}.__getitem__, on_change=persist_active_design)
        job = st.session_state.get("design_job")
        if job is not None and job.snapshot().state == "running":
            snapshot = job.snapshot()
            st.divider()
            st.caption(tr('Calcul en cours · {current}/{total}', current=snapshot.current, total=snapshot.total))
    return page


def _normalize_optional_sequence(text: str, field_name: str) -> tuple[str, str | None]:
    try:
        return normalize_dna(text, field_name=field_name), None
    except SequenceValidationError as exc:
        return "", str(exc)


def _sync_trigger_state(trigger: str) -> None:
    previous = st.session_state.get("trigger_snapshot")
    if previous == trigger:
        return
    st.session_state.trigger_snapshot = trigger
    st.session_state.aptamer_indices = []


def _inverse_extension_settings() -> InverseExtensionSettings | None:
    if st.session_state.get("design_mode") == "reading":
        return None
    if not (st.session_state.get("extension_enabled", False)
            and st.session_state.get("inverse_extension_enabled", False)):
        return None
    return InverseExtensionSettings(
        target_length=int(st.session_state.inverse_extension_target_length),
        extension_side=st.session_state.inverse_extension_side,
        top_candidates=int(st.session_state.inverse_extension_top),
        thermo=current_system_thermo(),
    )


def _effective_trigger(trigger: str) -> str:
    if st.session_state.get("design_mode") == "reading":
        return trigger
    if not trigger:
        return ""
    inverse = _inverse_extension_settings()
    if inverse is not None:
        # Invalid drafts remain visible without inventing or cropping bases.
        if inverse.target_length <= len(trigger):
            return ""
        prefix, suffix = inverse.extension_lengths(len(trigger))
        return "X" * prefix + trigger + "X" * suffix
    return trigger


def _clean_rna(text: str) -> str:
    # Retain invalid symbols so validation can reject them instead of silently
    # changing the user's architecture (in particular an obsolete N motif).
    return "".join((text or "").split()).upper().replace("T", "U")


def _current_architecture(trigger_length: int) -> SwitchArchitecture:
    minimum_length = 2 if st.session_state.get("design_mode") == "reading" else STANDARD_ACTIVATOR.trigger_region_length
    design_length = trigger_length if trigger_length >= minimum_length else STANDARD_ACTIVATOR.trigger_region_length
    default_architecture = default_architecture_for_trigger(design_length)
    rbs = _clean_rna(st.session_state.rbs_loop) or DEFAULT_RBS_SEQUENCE
    frame_linker = default_architecture.frame_linker
    if not st.session_state.linker_auto:
        frame_linker = _clean_rna(st.session_state.frame_linker)
    return replace(
        default_architecture,
        rbs_loop=rbs,
        t7_leader=_clean_rna(st.session_state.t7_leader),
        frame_linker=frame_linker,
        architecture_name=f"toehold_activator_standard_{design_length}nt",
    )


def _format_duration(seconds: float | None) -> str:
    if seconds is None:
        return tr("Estimation en cours…")
    total_seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, remaining_seconds = divmod(remainder, 60)
    if hours:
        return f"{hours} h {minutes:02d} min"
    if minutes:
        return f"{minutes} min {remaining_seconds:02d} s"
    return f"{remaining_seconds} s"


def _current_thermo() -> ThermoConditions:
    trigger_concentration_m = concentration_to_m(
        st.session_state.trigger_concentration_value,
        st.session_state.trigger_concentration_unit,
    )
    switch_concentration_m = concentration_to_m(
        st.session_state.switch_concentration_value,
        st.session_state.switch_concentration_unit,
    )
    return ThermoConditions(
        temperature_c=float(st.session_state.temperature_c),
        sodium_m=concentration_to_m(
            st.session_state.sodium_concentration_value,
            st.session_state.sodium_concentration_unit,
        ),
        magnesium_m=0.0,
        concentration_m=trigger_concentration_m,
        trigger_concentration_m=trigger_concentration_m,
        switch_concentration_m=switch_concentration_m,
        material=MODEL_ID,
        switch_material="rna",
        trigger_material="dna",
    )


def _current_reporter() -> ReporterContext:
    reporters = project_reporters(st.session_state)
    key = st.session_state.reporter_key
    if key not in reporters:
        key = "sfgfp"
        st.session_state.reporter_key = key
    return reporters[key]


def _build_request() -> RunRequest:
    reading = st.session_state.get("design_mode") == "reading"
    aptamer = normalize_dna(st.session_state.aptamer_input, field_name="aptamer")
    trigger = normalize_dna(st.session_state.trigger_input, field_name="trigger")
    inverse = _inverse_extension_settings()
    design_length = inverse.target_length if inverse is not None else len(trigger)
    if not reading and design_length < STANDARD_ACTIVATOR.trigger_region_length:
        raise ValueError(tr('Le toehold activateur standard nécessite un trigger d’au moins {trigger_region_length} nt. Activez l’extension du trigger.', trigger_region_length=STANDARD_ACTIVATOR.trigger_region_length))
    effective = _effective_trigger(trigger)
    if len(effective) < 2:
        raise ValueError(tr('La région du trigger utilisée pour le switch doit contenir au moins 2 bases.'))
    architecture = (default_architecture_for_trigger(len(effective)) if reading
                    else _current_architecture(len(effective)))
    reporter = _current_reporter()
    library_path = configured_library_path()
    if not library_path or not native_backend_identity(library_path):
        raise ValueError(tr(NATIVE_UNAVAILABLE))
    request = RunRequest(
        sequences=SequenceInput(
            aptamer_dna=aptamer,
            trigger_dna=trigger,
            molecule_name=st.session_state.molecule_name,
            binding_start=0,
            binding_length=len(trigger) if inverse is not None else len(effective),
        ),
        architecture=architecture,
        thermo=_current_thermo(),
        reporter=reporter,
        scoring=ScoringWeights(
            on_yield=float(st.session_state.weight_on),
            leak=float(st.session_state.weight_leak),
            defect=float(st.session_state.weight_defect),
            ddg=float(st.session_state.weight_ddg),
            rbs=float(st.session_state.weight_rbs),
            trigger_linearity=float(st.session_state.weight_trigger_linearity),
            rbs_calculator=float(st.session_state.weight_rbs_calculator),
        ),
        trials=int(st.session_state.trials),
        optimization_options={
            "max_evaluations": int(st.session_state.optimization_max_evaluations),
            "max_proposals": max(5000, 5 * int(st.session_state.optimization_max_evaluations)),
            "population_size": int(st.session_state.optimization_population_size),
            "lineages": int(st.session_state.optimization_lineages),
            "local_tolerance": float(st.session_state.optimization_local_tolerance),
        },
        output_dir=normalized_export_dir(st.session_state.output_dir),
        engine=ENGINE_ID,
        thermodynamics_library_path=library_path,
        exclude_stop_candidates=bool(st.session_state.exclude_stop_candidates),
        exclude_non_linear_candidates=bool(st.session_state.exclude_non_linear_candidates),
        inverse_extension=inverse,
        defect_metric="functional",
        score_scale="bounded_v1",
        functional_defect_weights=FunctionalDefectWeights(
            off_toehold=float(st.session_state.weight_functional_off_toehold),
            off_stem=float(st.session_state.weight_functional_off_stem),
            on_rbs_aug=float(st.session_state.weight_functional_on_rbs_aug),
            on_rbs_linker=float(st.session_state.weight_functional_on_rbs_linker),
        ),
    )
    if reading:
        return prepare_reading_request(
            request, st.session_state.reading_switch_input,
            reporter_is_complete=bool(st.session_state.reading_use_full_reporter),
        )
    return request


def _start_inverse_reanalysis(result, top: int) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.run_error = ""
    try:
        if active_job_running(st.session_state):
            raise ValueError(tr("Un calcul est déjà en cours. Suivez-le dans Résultats."))
        library_path = configured_library_path()
        if not library_path or not native_backend_identity(library_path):
            raise ValueError(tr(NATIVE_UNAVAILABLE))
        comparison = result.request.inverse_extension
        request = replace(result.request,
                          inverse_extension=replace(
                              comparison, top_candidates=int(top),
                              thermo=replace(comparison.thermo, material=MODEL_ID),
                          ),
                          engine=ENGINE_ID, thermo=replace(result.request.thermo, material=MODEL_ID),
                          thermodynamics_library_path=library_path)
        if st.session_state.get("active_design_path"):
            request = replace(request, output_dir=design_output_dir(Path(st.session_state.active_design_path)))
        request = request.normalized()
        job = DesignRunJob(request, analysis_source=result)
        job.start()
        st.session_state.design_job = job
        st.session_state.run_notice = tr("Comparaison des switches enregistrés : leurs séquences ne seront pas redessinées.")
    except Exception as exc:
        st.session_state.run_error = str(exc)


def _start_run() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.run_error = ""
    st.session_state.run_notice = ""
    tube_job = st.session_state.get("system_tube_job")
    if tube_job is not None and tube_job.snapshot().state == "running":
        st.session_state.run_error = tr("L’analyse des tubes est en cours. Suivez-la ou annulez-la dans la section du système initial.")
        return
    if active_job_running(st.session_state):
        st.session_state.run_error = "Un calcul est déjà en cours. Suivez-le dans Résultats."
        return
    try:
        request = _build_request().normalized()
        if st.session_state.get("active_design_path"):
            request = replace(request, output_dir=design_output_dir(Path(st.session_state.active_design_path)))
    except Exception as exc:
        st.session_state.run_error = str(exc)
        return
    requested_output = request.output_dir
    try:
        if st.session_state.get("active_design_path"):
            if not persist_active_design():
                return
            requested_output.mkdir(parents=True, exist_ok=True)
            writable_output, used_fallback = requested_output, False
        else:
            writable_output, used_fallback = ensure_writable_export_dir(requested_output)
    except OSError as exc:
        st.session_state.run_error = f"Impossible d'écrire les exports : {exc}"
        return
    if used_fallback:
        request = replace(request, output_dir=writable_output)
        st.session_state.output_dir = str(writable_output)
        st.session_state.run_notice = (
            f"Le dossier {requested_output} n'était pas inscriptible. "
            f"Les exports seront enregistrés dans {writable_output}."
        )
    try:
        job = DesignRunJob(request)
        job.start()
    except (OSError, ValueError, RuntimeError) as exc:
        st.session_state.run_error = str(exc)
        return
    st.session_state.design_job = job
    st.session_state.last_result = None
    st.session_state.last_exports = {}
    st.session_state.results_kind = "Design de switches"
    st.session_state.nav = "Résultats"


def _reset_all_settings() -> None:
    save_settings({"ui_language": st.session_state.ui_language})
    for key in (
        "settings_default_output",
    ):
        st.session_state.pop(key, None)
    st.session_state.reporter_key = "sfgfp"
    st.session_state.output_dir = str(default_export_dir())


def _switch_legend() -> list[tuple[str, str]]:
    return [
        ("Toehold", DOMAIN_COLORS["Toehold"]),
        ("Stem", DOMAIN_COLORS["Stem L"]),
        ("Bulge", DOMAIN_COLORS["Bulge"]),
        ("RBS", DOMAIN_COLORS["RBS"]),
        ("AUG", DOMAIN_COLORS["AUG"]),
        ("Reporter", DOMAIN_COLORS["Reporter"]),
        ("Trigger", TRIGGER_COLOR),
    ]


def _render_initial_system() -> tuple[str, str]:
    with st.container(border=True):
        st.subheader(tr('Système aptamère–trigger'), help=tr('Renseignez les séquences disponibles, les conditions du système et, si elles sont connues, les bases impliquées dans la liaison au ligand.'))
        if st.session_state.get("design_using_extension"):
            st.info(tr("Vous utilisez un trigger étendu. Les séquences d’origine restent conservées dans ce design."))
            with st.expander(tr("Séquences d’origine")):
                st.caption(tr("Aptamère original"))
                st.code(st.session_state.get("design_original_aptamer", ""), language=None)
                st.caption(tr("Trigger original"))
                st.code(st.session_state.get("design_original_trigger", ""), language=None)
                st.button(tr("Revenir aux séquences d’origine"), key="restore_original_sequences",
                          on_click=_restore_original_sequences, disabled=active_job_running(st.session_state))
        left, right = st.columns(2, gap="large", vertical_alignment="top")
        with left:
            st.markdown(tr('#### Séquences du système'))
            st.text_input(tr('Molécule cible'), key="molecule_name")
            st.text_area(
                tr("Séquence de l'aptamère"), key="aptamer_input", height=110,
                placeholder=tr('Séquence ADN — espaces et retours à la ligne acceptés'),
            )
            st.text_area(tr('Séquence du trigger'), key="trigger_input", height=110,
                         placeholder=tr('Séquence ADN'))
            trigger, trigger_error = _normalize_optional_sequence(st.session_state.trigger_input, "trigger")
            aptamer, aptamer_error = _normalize_optional_sequence(st.session_state.aptamer_input, "aptamer")
            for raw, error in ((st.session_state.aptamer_input, aptamer_error), (st.session_state.trigger_input, trigger_error)):
                if raw.strip() and error:
                    st.warning(tr(error))
            if aptamer or trigger:
                st.caption(tr('Aptamère : {value} nt · Trigger : {value1} nt', value=len(aptamer), value1=len(trigger)))
            if trigger and len(trigger) < 24:
                st.warning(tr('Attention : votre trigger mesure {value} nt. Sa longueur peut être trop faible ; vous pouvez activer l’extension du trigger.', value=len(trigger)))
        with right:
            render_system_conditions()
        _sync_trigger_state(trigger)
        sync_ligand_selection(aptamer)
        with st.expander(tr('Bases de l’aptamère impliquées dans la liaison au ligand (si connues)'), expanded=False):
            render_ligand_selection(aptamer)
        render_initial_structure(aptamer, trigger)
        render_initial_tubes(aptamer, trigger)
    return aptamer, trigger


def _restore_original_sequences() -> None:
    st.session_state.aptamer_input = st.session_state.get("design_original_aptamer", "")
    st.session_state.trigger_input = st.session_state.get("design_original_trigger", "")
    st.session_state.design_using_extension = False
    st.session_state.trigger_snapshot = None
    persist_active_design()


def _use_custom_linker() -> None:
    st.session_state.linker_auto = False


def _restore_architecture_defaults(trigger: str) -> None:
    st.session_state.t7_leader = STANDARD_ACTIVATOR.t7_leader
    st.session_state.rbs_loop = STANDARD_ACTIVATOR.rbs_loop
    st.session_state.linker_auto = True
    effective = _effective_trigger(trigger)
    standard = default_architecture_for_trigger(max(len(effective), STANDARD_ACTIVATOR.trigger_region_length))
    st.session_state.frame_linker = standard.resolved_frame_linker()


def _optimize_linker() -> None:
    st.session_state.frame_linker = "X" * len(_clean_rna(st.session_state.frame_linker))
    st.session_state.linker_auto = False


def _render_architecture(trigger: str) -> tuple[SwitchArchitecture, str]:
    with st.container(border=True):
        st.subheader(tr('Architecture du toehold switch'))
        render_architecture_selector()
        sequence_columns = st.columns(2)
        sequence_columns[0].text_input(tr('Leader 5′'), key="t7_leader")
        sequence_columns[1].text_input(
            tr('RBS'), key="rbs_loop",
            help=tr('La boucle cible contient 3 bases libres (XXX), optimisées par le moteur hybride expérimental et non appariées, suivies de votre RBS.'),
        )
        effective = _effective_trigger(trigger)
        if st.session_state.get("design_mode") != "reading" and trigger and len(effective) < STANDARD_ACTIVATOR.trigger_region_length:
            st.info(tr('Le modèle standard nécessite un trigger d’au moins {trigger_region_length} nt. Activez l’extension du trigger.', trigger_region_length=STANDARD_ACTIVATOR.trigger_region_length))
        architecture = _current_architecture(len(effective))
        if st.session_state.linker_auto:
            st.session_state.frame_linker = architecture.resolved_frame_linker()
        st.text_input(
            tr('Linker'), key="frame_linker", on_change=_use_custom_linker,
            help=tr('Cliquez dans la séquence pour remplacer le linker par défaut. A, C, G, U ou T imposent une base ; X laisse une base à optimiser. Une séquence vide supprime le linker.'),
        )
        linker_pattern = _clean_rna(st.session_state.frame_linker)
        linker_actions = st.columns(2)
        linker_actions[0].button(
            tr('Optimiser le linker (XXX…)'), key="optimize_linker", on_click=_optimize_linker,
            disabled=not linker_pattern or bool(set(linker_pattern) - set("ACGUX")) or set(linker_pattern) == {"X"},
            help=tr('Remplace toutes les bases du linker par X sans changer sa longueur. Sa séquence sera choisie pendant la conception.'),
            width="stretch",
        )
        linker_actions[1].button(
            tr('Rétablir les valeurs par défaut'), key="restore_architecture_defaults",
            on_click=_restore_architecture_defaults, args=(trigger,),
            disabled=(st.session_state.linker_auto
                      and _clean_rna(st.session_state.t7_leader) == STANDARD_ACTIVATOR.t7_leader
                      and _clean_rna(st.session_state.rbs_loop) == STANDARD_ACTIVATOR.rbs_loop),
            help=tr('Rétablit le leader GGG, le RBS, le linker et l’architecture standard adaptée au trigger.'),
            width="stretch",
        )
        try:
            architecture.validate_for_trigger(architecture.trigger_region_length)
        except ValueError as exc:
            st.warning(tr(str(exc)))
            return architecture, ""
        linker = architecture.resolved_frame_linker()
        st.caption(tr('Linker : {length} nt · après les {stem_length} nt de tige qui suivent AUG, avant le rapporteur.',
                      length=len(linker), stem_length=architecture.upper_stem2_length + architecture.stem_length))

    # A short input cannot be drawn as if it filled the standard trigger region.
    preview_trigger = effective if len(effective) == architecture.trigger_region_length else ""
    return architecture, preview_trigger


def _render_live_previews(architecture: SwitchArchitecture, effective_trigger: str) -> None:
    aptamer, aptamer_error = _normalize_optional_sequence(st.session_state.aptamer_input, "aptamer")
    trigger, trigger_error = _normalize_optional_sequence(st.session_state.trigger_input, "trigger")
    if (aptamer_error or trigger_error or not aptamer or not trigger
            or not effective_trigger or len(effective_trigger) != architecture.trigger_region_length):
        return
    reporter = _current_reporter()
    try:
        model = build_switch_model(
            architecture, effective_trigger, reporter.sequence_after_start_rna,
            reporter.label,
        )
    except ValueError as exc:
        st.warning(tr(str(exc)))
        return
    with st.container(border=True):
        st.subheader(tr('Aperçus du toehold switch'), help=tr('Structures cibles construites à partir de vos séquences. X indique les bases à optimiser par le moteur hybride expérimental ; les structures prédites seront disponibles après le calcul.'))
        structure_tab, linear_tab = st.tabs([tr(label) for label in MOLECULAR_VIEW_LABELS])
        with structure_tab:
            state = st.radio(
                tr('État du switch'),
                ("OFF", "ON + trigger"),
                horizontal=True,
                key="preview_state",
            )
            diagram = partial(
                switch_preview_svg,
                architecture,
                effective_trigger,
                reporter.sequence_after_start_rna,
                reporter.label,
                state="on" if state.startswith("ON") else "off",
                width=1000,
                height=500,
            )
            state_name = "on" if state.startswith("ON") else "off"
            render_diagram(
                diagram, legend=_switch_legend(), filename=f"switch_{state_name}", key="switch_preview",
                structure=model.on_structure if state_name == "on" else model.off_structure,
                download_labels=("Télécharger cet aperçu SVG", "Télécharger cet aperçu PNG"),
            )
        with linear_tab:
            diagram = partial(
                linear_sequence_svg,
                architecture,
                effective_trigger,
                reporter.sequence_after_start_rna,
                reporter.label,
            )
            render_diagram(
                diagram, filename="switch_lineaire", key="switch_preview_linear",
                caption="X indique une base laissée libre pour l’optimisation du moteur hybride expérimental.",
                download_labels=("Séquence linéaire SVG", "Séquence linéaire PNG"),
            )


def _scoring_help(thermo: ThermoConditions) -> dict[str, str]:
    """Describe the measured quantities using the actual analysis-tube inputs."""
    def concentration(value: float) -> str:
        unit = preferred_concentration_unit(value)
        return f"{concentration_from_m(value, unit):g} {unit}"

    trigger_m = thermo.resolved_trigger_concentration_m
    switch_m = thermo.resolved_switch_concentration_m
    ratio = f"{trigger_m / switch_m:g}" if switch_m > 0 else tr("indéfini")
    trigger_value, switch_value = concentration(trigger_m), concentration(switch_m)
    return {
        "weight_rbs": (
            tr('**Stabilité de la région RBS–linker.** Énergie MFE (kcal/mol) calculée avec ViennaRNA sur le fragment ARN isolé, du préfixe de boucle avant le RBS à la fin du linker. Un ΔG moins négatif est favorisé. Ce poids guide aussi l’affinage après les défauts locaux : coût = défaut global + w × max(0, −ΔG) / 100, en respectant les limites locales et les gardes de concentration et de duplex. Au classement final, gagner 1 kcal/mol ajoute w points bruts. Un poids nul désactive cette préférence. Le seuil −8 kcal/mol reste indicatif ; aucune exclusion automatique n’en découle. Ce fragment isolé ne mesure pas l’ouverture simultanée du RBS dans le complexe ON.')
        ),
        "weight_rbs_calculator": tr("**RBS ON/OFF · prédiction exploratoire.** OSTIR 1.1.3 avec ViennaRNA 2.7.2 estime l’initiation de la traduction pour E. coli à 37 °C, selon sa référence ARN à 1,021 M de sel. Le calcul conserve la séquence du candidat dans une fenêtre de 35 bases de chaque côté de l’AUG voulu. OFF utilise l’ARN sans contrainte. En ON, seules les bases prévues pour le duplex avec le trigger sont rendues indisponibles aux appariements ARN internes ; aucune énergie de liaison ADN–ARN n’est ajoutée au modèle OSTIR. C’est une approximation de l’état occupé, pas une validation du complexe hybride. Le bonus combine un ON élevé et un contraste ON/OFF : qON = TIR_ON / (TIR_ON + K), avec K = exp(7,279194329), puis qcontraste = max(0, (TIR_ON − TIR_OFF) / (TIR_ON + TIR_OFF)) et mérite = qON × qcontraste. TIR_ON et TIR_OFF désignent les initiations relatives prédites par OSTIR, pas les pourcentages d’association en tube. Le score brut ajoute poids × 100 × mérite ; poids 3 par défaut. L’affinage ajoute un coût poids × (1 − mérite) en préservant les priorités locales. Cette normalisation et ce poids sont des choix de classement réglables, pas des paramètres optimaux établis expérimentalement. Les valeurs OSTIR ne sont ni un taux absolu de protéine ni un ratio d’induction validé. Un poids nul désactive ce critère."),
        "weight_defect": (
            tr('**Compromis régional entre les cibles OFF et ON.** Le défaut structural d’ensemble du moteur hybride expérimental mesure, dans chaque zone, la fraction attendue de nucléotides qui ne respectent pas leur état cible. OFF combine « toehold libre » et « tiges avant/après AUG maintenues » ; ON combine « préfixe RBS → fin AUG libre » et « préfixe RBS → fin linker libre ». Dans le classique 30 nt, les tiges OFF comptent 6 paires avant AUG et 9 après AUG. Les trois bases du préfixe de boucle avant le RBS sont incluses dans les deux zones ON. Ces zones se recouvrent : RBS → AUG contribue aux deux critères. Le duplex avec le trigger ne fait plus partie de ce malus. Le malus utilise la moyenne pondérée des quatre zones selon leurs poids relatifs : par défaut, les trois premiers critères pèsent chacun 26,7 % et RBS → fin du linker pèse 20 %, quelle que soit la longueur des zones. Les moyennes OFF et ON suivent les poids de leurs zones ; un état de poids nul est désactivé. Le RBS et AUG peuvent être libres en OFF. Ce calcul porte sur l’ensemble des conformations, pas seulement sur la structure MFE, et ne contient pas le défaut de concentration. Plus faible est meilleur ; 0 est idéal. Les valeurs sont des fractions de 0 à 1 : 0,05 correspond à 5 %. Le score brut retire w × (100 × compromis OFF/ON) ; réduire le compromis d’un point de pourcentage apporte w points bruts. Le coefficient w règle l’amplitude totale du malus, séparément des poids relatifs des zones. Ces quatre poids relatifs guident aussi la première phase d’optimisation ; la seconde réduit le défaut global et la stabilité du RBS–linker isolé, selon son poids, en préservant les défauts locaux et les gardes de concentration et de duplex. Aucun seuil d’exclusion n’est ajouté. Ce n’est pas une mesure d’expression.')
        ),
        "weight_leak": association_with_aptamer_help(thermo, weight=st.session_state.weight_leak),
        "weight_on": (
            tr('**Association ON avec le trigger libre.** Pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre. Le tube contient seulement le trigger libre ({trigger_value}) et le switch ({switch_value}), soit un ratio trigger:switch de {ratio}:1, sans aptamère ni ligand. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus élevé est meilleur pour favoriser la capture du trigger dans le tube de switch, séparé de celui des aptamères. Par défaut, l’association ON participe au score et le critère d’association en présence d’aptamère reste informatif (poids 0). Le score brut ajoute w × ON ; gagner un point de pourcentage apporte w points au score brut. Ce pourcentage est un indicateur du complexe ON ; il ne mesure pas directement l’activation de la traduction ni la quantité de protéine produite.', trigger_value=trigger_value, switch_value=switch_value, ratio=ratio)
        ),
        "weight_ddg": (
            tr('**Énergie de formation du complexe trigger–switch.** ΔΔG = ΔG MFE du complexe ON − ΔG MFE du switch seul − ΔG MFE du trigger seul, en kcal/mol, avec les conditions choisies. Une valeur plus négative favorise énergétiquement l’association. Le score brut retire w × ΔΔG : une diminution de 1 kcal/mol apporte w points au score brut. Son poids est faible par défaut car ce critère complète le rendement en tube sans décrire à lui seul le fonctionnement. Ce n’est ni une barrière d’activation ni une vitesse de réaction.')
        ),
        "weight_trigger_linearity": tr("**Trigger étendu libre, en mode extension inverse.** Le moteur hybride expérimental calcule les probabilités du trigger seul, avec sa chimie et les conditions du design du switch. Le défaut est la fraction moyenne de bases appariées : 1 − moyenne des probabilités non appariées. Le score brut retire w × 100 × ce défaut ; poids 3 par défaut. Une baisse d’un point de pourcentage de bases appariées apporte donc 3 points bruts avec ce réglage. Ce poids est une priorité de classement réglable, pas un optimum expérimental établi. Un poids nul désactive ce malus. Cette mesure porte sur l’ensemble des conformations, pas seulement sur la MFE ; elle ne garantit pas que toutes les bases soient libres simultanément. Elle ne s’applique pas au design sans extension inverse."),
    }


def _on_reporter_selected() -> None:
    if st.session_state.reporter_key == ADD_REPORTER:
        previous = st.session_state.get("last_reporter_key", "sfgfp")
        st.session_state.reporter_key = previous if previous in project_reporters(st.session_state) else "sfgfp"
        st.session_state.nav = "Réglages"
        st.session_state.settings_tab = "Reporters"
        st.session_state.reporter_return_to_design = True
    else:
        st.session_state.last_reporter_key = st.session_state.reporter_key


def _return_to_design() -> None:
    st.session_state.nav = "Conception"
    st.session_state.design_section = "switch"
    st.session_state.reporter_return_to_design = False


def _render_run_configuration(*, reading: bool = False) -> None:
    with st.container(border=True):
        st.subheader(tr('3 · Calcul des métriques' if reading else 'Calcul et classement du toehold switch'))
        conditions_tab, scoring_tab, output_tab = st.tabs([tr('Conditions'), tr('Scoring'), tr('Sortie')])
        with conditions_tab:
            native_ready = bool(native_backend_identity(configured_library_path()))
            if not native_ready:
                st.info(tr(NATIVE_UNAVAILABLE))
            cols = st.columns(3)
            if reading:
                cols[0].caption(tr('Une analyse des séquences fournies, sans optimisation.'))
            else:
                cols[0].number_input(tr('Essais'), min_value=1, max_value=100000, step=1, key="trials")
            if not 25. <= float(st.session_state.temperature_c) <= 55.:
                st.session_state.temperature_c = max(25., min(55., float(st.session_state.temperature_c)))
                st.info(tr("Les conditions proposées sont ramenées aux limites du modèle."))
            cols[1].number_input(
                tr('Température (°C)'), min_value=25.0, max_value=55.0, step=0.1, key="temperature_c"
            )
            reporters = project_reporters(st.session_state)
            if st.session_state.reporter_key not in reporters:
                st.session_state.reporter_key = "sfgfp"
            cols[2].selectbox(
                tr('Reporter'),
                options=[*reporters, ADD_REPORTER],
                format_func={**{key: reporter.label for key, reporter in reporters.items()},
                             ADD_REPORTER: tr("Ajouter un reporter…")}.__getitem__,
                key="reporter_key",
                on_change=_on_reporter_selected,
                disabled=reading and st.session_state.reading_use_full_reporter,
                help=tr('Conservez le rapporteur complet. Le design utilise automatiquement ses 30 premières bases après AUG ; les analyses en contexte utilisent ensuite la séquence complète.'),
            )
            selected_reporter = reporters.get(st.session_state.reporter_key)
            if reading:
                st.checkbox(tr('Le rapporteur est entièrement fourni dans la séquence du switch'),
                            key="reading_use_full_reporter")
                st.caption(tr('Sinon, choisissez le rapporteur correspondant au fragment fourni. Sa séquence complète sera conservée pour les analyses suivantes.'))
            if not reading and selected_reporter and selected_reporter.full_sequence_after_start_rna:
                cols[2].caption(tr('{count} bases après AUG conservées · {prefix} bases utilisées pour le design.',
                                   count=len(selected_reporter.full_sequence_after_start_rna),
                                   prefix=len(selected_reporter.sequence_after_start_rna)))
            st.markdown(tr('##### Concentrations'), help=tr('Le rendement ON est calculé par rapport à la concentration initiale du toehold switch.'))
            concentration_cols = st.columns(4)
            with concentration_cols[0]:
                _render_concentration_control(
                    tr('Trigger'),
                    value_key="trigger_concentration_value",
                    unit_key="trigger_concentration_unit",
                    min_m=1e-9,
                    max_m=1e-2,
                )
            with concentration_cols[1]:
                _render_concentration_control(
                    tr('Toehold switch'),
                    value_key="switch_concentration_value",
                    unit_key="switch_concentration_unit",
                    min_m=1e-9,
                    max_m=1e-2,
                )
            with concentration_cols[2]:
                _render_concentration_control(
                    tr('Sodium'),
                    value_key="sodium_concentration_value",
                    unit_key="sodium_concentration_unit",
                    min_m=0.12,
                    max_m=1.0,
                )
            with concentration_cols[3]:
                st.session_state.magnesium_concentration_value = 0.0
                _render_concentration_control(
                    tr('Magnésium'),
                    value_key="magnesium_concentration_value",
                    unit_key="magnesium_concentration_unit",
                    min_m=0.0,
                    max_m=1.0,
                    disabled=True,
                    help=tr('Le modèle hybride expérimental impose Mg²⁺ = 0.'),
                )
        with scoring_tab:
            st.markdown(tr('##### Poids du classement'), help=tr('Priorités par défaut : accessibilité RBS–linker et compromis structural OFF/ON, puis association ON et ΔΔG. Le bonus RBS ON/OFF a un poids de 3. En extension inverse, le malus de repliement du trigger seul a un poids de 3. L’association trigger–switch en présence d’aptamère reste informative avec un poids de 0. Les poids restent modifiables ; 0 désactive un critère.'))
            scoring_help = _scoring_help(_current_thermo())
            criteria = (
                ("ΔG RBS–linker", "weight_rbs"), ("Défaut OFF/ON · malus", "weight_defect"),
                ("ON % · bonus", "weight_on"), ("ΔΔG · liaison", "weight_ddg"),
                (ASSOCIATION_WITH_APTAMER_WEIGHT_LABEL, "weight_leak"),
                ("RBS ON/OFF · bonus", "weight_rbs_calculator"),
            )
            for offset in (0, 3):
                for column, (label, key) in zip(st.columns(3), criteria[offset:offset+3]):
                    column.number_input(tr(label), 0.0, 20.0, step=0.05, key=key, help=scoring_help[key])
            if st.session_state.inverse_extension_enabled and not reading:
                st.number_input(tr("Linéarité du trigger seul · malus"), 0.0, 20.0, step=0.05,
                                key="weight_trigger_linearity", help=scoring_help["weight_trigger_linearity"])
            render_regional_defect_controls()
            with st.expander(tr('Voir la formule')):
                st.latex(
                    r"S_{brut} = w_{ON}\,ON - w_{A}\,A - w_{defect}(100\,D_{\mathrm{OFF/ON}}) "
                    r"- w_{\Delta\Delta G}\,\Delta\Delta G + w_{RBS}\,\Delta G_{RBS}"
                    r"+100\,w_{RBS\ ON/OFF}\,M_{RBS}"
                )
                st.caption(tr("A : pourcentage de switches associés au trigger en présence d’aptamère."))
                st.caption(tr("M_RBS : mérite normalisé OSTIR, combinant initiation ON et contraste ON/OFF. Voir l’aide du poids RBS ON/OFF pour le modèle et ses limites."))
                if st.session_state.inverse_extension_enabled and not reading:
                    st.latex(r"S_{brut}\leftarrow S_{brut}-w_{trigger}(100\,D_{trigger\ seul})")
                st.latex(r"Score_{/100}=50+\frac{100}{\pi}\arctan\!\left(\frac{S_{brut}}{100}\right)")
                st.latex(
                    r"D_{\mathrm{OFF/ON}}=\frac{\sum_{i=1}^{4} a_i D_i}{\sum_{i=1}^{4} a_i},"
                    r"\qquad a_i\geq0,\quad\sum_i a_i>0"
                )
                st.caption(tr("Le score affiché est borné de 0 à 100 par une transformation fixe de la somme pondérée. Elle conserve l’ordre des candidats ; la somme brute départage les notes arrondies identiques. Un score brut nul donne 50/100. Ce n’est pas un pourcentage de fonctionnement. Les poids agissent sur le score brut, avant conversion."))
            if reading:
                st.caption(tr('Toutes les métriques sont calculées. Un STOP prématuré ou un RBS–linker replié est signalé dans les résultats.'))
            else:
                st.checkbox(tr('Exclure les candidats avec un STOP prématuré'), key="exclude_stop_candidates",
                            help=tr('La génération et les mutations évitent les STOP avant le rapporteur. Si les bases fixées imposent un STOP, la conception s’arrête avant les analyses thermodynamiques.'))
                if not st.session_state.exclude_stop_candidates:
                    st.warning(tr(STOP_FILTER_DISABLED_WARNING))
                st.checkbox(
                    tr('Exclure les candidats sans linéarité parfaite du RBS–linker en ON'),
                    key="exclude_non_linear_candidates",
                    help=(
                        tr('Toutes les bases, de la première base de la boucle contenant le RBS (y compris les bases optimisées avant le RBS) jusqu’à la dernière base avant le rapporteur, doivent être non appariées dans les structures MFE ON renvoyées par le moteur hybride expérimental pour le complexe trigger–switch. Un seul appariement dans cette région suffit à exclure le candidat. Ce critère porte sur ces structures prédites, pas sur l’ensemble des conformations.')
                    ),
                )
        with output_tab:
            if st.session_state.get("active_design_path"):
                st.caption(tr("Tous les calculs sont enregistrés dans le dossier du design actif."))
                st.code(st.session_state.output_dir, language=None)
            else:
                st.text_input(tr('Dossier local des runs'), key="output_dir")
            st.caption(tr('Exports : JSON, CSV, XLSX, SVG, journal et projet .aptaswitch.json.'))

        if not reading:
            with st.expander(tr("Optimisation par essai")):
                st.caption(tr("Une population variée est explorée, puis plusieurs lignées réduisent les quatre défauts locaux selon vos poids. La dernière phase réduit le défaut global sous des plafonds locaux et des gardes de concentration et de duplex fixes."))
                st.number_input(tr("Évaluations maximales par essai"), min_value=1, max_value=100000,
                                step=100, key="optimization_max_evaluations",
                                help=tr("Chaque essai partage son budget entre la population initiale et les phases locale et globale. Augmenter ce budget autorise davantage de calculs sans garantir une meilleure séquence."))
                population = st.number_input(tr("Population initiale"), min_value=1, max_value=512, step=1,
                                             key="optimization_population_size")
                st.session_state.optimization_lineages = min(int(st.session_state.optimization_lineages), int(population))
                st.number_input(tr("Lignées explorées"), min_value=1, max_value=min(64, int(population)), step=1,
                                key="optimization_lineages")
                st.number_input(tr("Tolérance sur les défauts locaux"), min_value=0.0, max_value=1.0, step=0.001,
                                format="%.3f", key="optimization_local_tolerance",
                                help=tr("Fraction de défaut autorisée au-dessus de la référence locale figée pendant la phase globale. 0,005 représente 0,5 point de pourcentage."))
                st.caption(tr("Les défauts obtenus et le motif d’arrêt restent enregistrés, même si les cibles ne sont pas atteintes."))
            if int(st.session_state.trials) > 50:
                st.warning(tr('Avec le moteur hybride expérimental, plus de 50 essais peuvent prendre plusieurs heures ou plusieurs jours.'))
                st.checkbox(tr('Je confirme ce calcul long'), key="long_run_confirm")
        can_start = reading or int(st.session_state.trials) <= 50 or st.session_state.get("long_run_confirm", False)
        try:
            configured_request = _build_request().normalized()
        except ValueError as exc:
            can_start = False
            if (reading and all(st.session_state.get(key, "").strip()
                                for key in ("aptamer_input", "trigger_input", "reading_switch_input"))):
                st.warning(tr(str(exc)))
            elif (st.session_state.get("inverse_extension_enabled")
                    and st.session_state.aptamer_input and st.session_state.trigger_input):
                st.caption(tr(str(exc)))
        else:
            if reading:
                detected = configured_request.architecture
                st.caption(tr("Découpage reconnu automatiquement · toehold : {toehold} nt · tige inférieure : {stem} nt · linker : {linker} nt.",
                              toehold=detected.toehold_length, stem=detected.stem_length,
                              linker=len(detected.resolved_frame_linker())))
        running = active_job_running(st.session_state)
        st.button(
            tr('Calculer les métriques' if reading else 'Lancer la conception'),
            type="primary",
            on_click=_start_run,
            disabled=not can_start or running,
            width="stretch",
        )
        if st.session_state.run_error:
            st.error(tr(st.session_state.run_error))
        if st.session_state.run_notice:
            st.info(tr(st.session_state.run_notice))


def _render_reading_page() -> None:
    st.caption(tr('Fournissez l’aptamère, le trigger et le switch en 5′ → 3′. Les séquences saisies sont conservées pour le calcul et l’enregistrement.'))
    with st.container(border=True):
        st.subheader(tr('1 · Séquences à analyser'))
        st.text_input(tr('Molécule cible'), key="molecule_name")
        left, right = st.columns(2)
        left.text_area(tr("Séquence de l'aptamère"), key="aptamer_input", height=110,
                       help=tr('Aptamère ADN · espaces et retours à la ligne acceptés.'))
        right.text_area(tr('Séquence du trigger'), key="trigger_input", height=110)
        st.text_area(tr('Séquence du toehold switch'), key="reading_switch_input", height=150,
                     help=tr('Incluez le leader, le switch ARN, le linker et le fragment de rapporteur déjà présent. Le trigger et l’aptamère sont en ADN.'))
        aptamer, aptamer_error = _normalize_optional_sequence(st.session_state.aptamer_input, "aptamer")
        trigger, trigger_error = _normalize_optional_sequence(st.session_state.trigger_input, "trigger")
        for raw, error in ((st.session_state.aptamer_input, aptamer_error),
                           (st.session_state.trigger_input, trigger_error)):
            if raw.strip() and error:
                st.warning(tr(error))
        _sync_trigger_state(trigger)
        st.caption(tr('Aptamère : {aptamer} nt · Trigger : {trigger} nt · Switch : {switch} nt',
                      aptamer=len(aptamer), trigger=len(trigger),
                      switch=len("".join(st.session_state.reading_switch_input.split()))))
    st.caption(tr("Le découpage est reconnu parmi les architectures prises en charge, à partir du switch et du trigger fournis. Le rapporteur sélectionné permet d’identifier leur jonction, sans modifier vos séquences."))
    _render_run_configuration(reading=True)


def _render_design_page() -> None:
    st.title(tr('Conception'), anchor=False)
    modes = {"design": tr("Conception"), "reading": tr("Lecture de séquences")}
    st.radio(tr("Mode du design"), tuple(modes), key="design_mode",
             format_func=modes.__getitem__, horizontal=True,
             disabled=active_job_running(st.session_state))
    if st.session_state.design_mode == "reading":
        _render_reading_page()
        return
    section = render_design_navigation()
    if section == "system":
        _, trigger = _render_initial_system()
        render_blast_screen(trigger)
        next_section = "extension" if st.session_state.extension_enabled else "switch"
        label = "Continuer vers l’extension" if next_section == "extension" else "Continuer vers la conception du switch"
        st.button(tr(label), on_click=select_design_section, args=(next_section,),
                  type="primary", width="stretch")
        return
    aptamer, _ = _normalize_optional_sequence(st.session_state.aptamer_input, "aptamer")
    trigger, _ = _normalize_optional_sequence(st.session_state.trigger_input, "trigger")
    _sync_trigger_state(trigger)
    sync_ligand_selection(aptamer)
    if section == "extension":
        render_extension_page(aptamer, trigger)
        return
    architecture, effective = _render_architecture(trigger)
    _render_live_previews(architecture, effective)
    _render_run_configuration()


@st.fragment(run_every=1.0)
def _design_job_monitor() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    job = st.session_state.get("design_job")
    if job is None:
        return
    snapshot = job.snapshot()
    if snapshot.state == "running":
        total = max(1, snapshot.total)
        metrics = st.columns(4)
        metrics[0].metric(tr('Moteur'), job.request.engine.upper())
        metrics[1].metric(tr('Progression'), tr('{current} / {total}', current=snapshot.current, total=snapshot.total))
        metrics[2].metric(tr('Temps écoulé'), _format_duration(snapshot.elapsed_seconds))
        metrics[3].metric(tr('Temps restant estimé'), _format_duration(snapshot.eta_seconds))
        # Finishing the last trial still leaves final ranking and exports.
        st.progress(min(0.99, snapshot.current / total), text=tr(snapshot.message))
        with st.expander(tr('Console du moteur')):
            st.code("\n".join(tr(line) for line in snapshot.log) if snapshot.log else tr("En attente des premiers messages…"),
                    language=None, wrap_lines=True, height=220)
            st.caption(tr('Le temps restant est recalculé après chaque essai terminé.'))
        render_design_run_controls(job, snapshot)
        if st.session_state.get("run_error"):
            st.error(tr(st.session_state.run_error))
        return
    if snapshot.state == "failed":
        st.error(tr(snapshot.error or "Le calcul a échoué."))
    if job.consume_once():
        if snapshot.result is not None:
            st.session_state.last_result = snapshot.result
            st.session_state.last_exports = snapshot.exports or {}
        if snapshot.state == "failed":
            st.session_state.run_error = snapshot.error or "Le calcul a échoué."
        persist_active_design()
        st.rerun(scope="app")


def _render_results_page() -> None:
    st.title(tr('Résultats'), anchor=False)
    render_run_library()
    kind = st.radio(
        tr('Type de résultats'), ("Design de switches", "Extension de triggers"),
        horizontal=True, key="results_kind",
        format_func={label: tr(label) for label in ("Design de switches", "Extension de triggers")}.__getitem__,
    )
    if kind == "Extension de triggers":
        render_saved_run_details("extension")
        render_extension_results(format_duration=_format_duration)
        return
    _design_job_monitor()
    result = st.session_state.get("last_result")
    if result is None:
        if not st.session_state.get("design_job"):
            st.info(tr('Ouvrez un run enregistré ci-dessus ou lancez un calcul depuis Conception.'))
        return

    excluded_count = sum(getattr(candidate, "excluded", False) for candidate in result.candidates)
    st.caption(tr("{run} · {engine} · {ranked} classés · {excluded} exclus", run=result.run_id,
                  engine=result.request.engine.upper(), ranked=len(result.candidates) - excluded_count,
                  excluded=excluded_count))
    if result.status == "cancelled":
        st.info(tr("Calcul arrêté : les candidats terminés restent consultables."))
    partial = render_design_continuation(result)
    if result.candidates and excluded_count == len(result.candidates):
        st.warning(tr('Tous les candidats sont exclus par le criblage. Aucun candidat n’a été retenu pour le classement.'))
    if st.session_state.get("run_error"):
        st.error(tr(st.session_state.run_error))

    _render_optimization_status(result)

    has_inverse = result.request.inverse_extension is not None or bool(result.inverse_extension_analysis)
    section = "Classement"
    if has_inverse:
        sections = ("Classement", "Comparaison de l’extension")
        section = st.radio(tr("Vue des résultats"), sections,
                           horizontal=True, key="switch_results_section", format_func={name: tr(name) for name in sections}.__getitem__,
                           label_visibility="collapsed")
    if section == "Comparaison de l’extension":
        if partial:
            st.info(tr("La comparaison des extensions sera disponible après la finalisation de ce calcul."))
            return
        render_inverse_extension_results(result, on_reanalyze=_start_inverse_reanalysis)
        return

    render_results_table(result)
    with st.expander(tr("Comprendre le classement")):
        st.caption(tr("Score sur 100 : plus élevé est meilleur. La somme pondérée brute reste consultable dans le détail du candidat."))
        render_defect_method(result)
    if result.warnings:
        with st.expander(tr("Points à vérifier ({count})", count=len(result.warnings))):
            for warning in result.warnings:
                st.warning(tr(warning))
    render_saved_run_details("design")
    exports = st.session_state.get("last_exports", {})
    # Loading a checkpoint can also discover older JSON/CSV files beside it.
    # None of these is a final export package for an unfinished calculation.
    exports = {} if partial else {kind: path for kind, path in exports.items() if kind != "checkpoint"}
    if exports:
        zip_data = exports_zip_bytes(exports)
        export_options = [ExportOption("Tous les exports (.zip)", zip_data,
                                      f"{result.run_id}_exports.zip", "application/zip", "switch_run_zip")]
        for kind, path_text in exports.items():
            path = Path(path_text)
            if path.is_file():
                export_options.append(ExportOption(kind.upper(), path.read_bytes(), path.name,
                                                   "application/octet-stream", f"switch_run_{kind}"))
        render_export_menu(export_options, label="Exporter les résultats")

    if not result.candidates:
        if partial:
            st.info(tr("Aucun candidat terminé n’est encore disponible dans ce point de reprise."))
        else:
            st.warning(tr("Le calcul n'a produit aucun candidat."))
        return

    candidates_by_trial = {candidate.trial: candidate for candidate in result.candidates}
    trial_options = list(candidates_by_trial)
    candidate_labels = {}
    rank = 0
    for candidate in result.candidates:
        if candidate.excluded:
            candidate_labels[candidate.trial] = tr("{status} · essai {trial}", status=candidate_status(candidate), trial=candidate.trial)
        else:
            rank += 1
            score = candidate_score_100(candidate, result.request.score_scale)
            candidate_labels[candidate.trial] = tr("Rang {rank} · essai {trial} · score {score}",
                                                   rank=rank, trial=candidate.trial,
                                                   score="—" if score is None else f"{score:.2f}")
    if st.session_state.get("selected_trial") not in trial_options:
        remembered = st.session_state.get("results_candidate_selection", {})
        same_run = remembered.get("run") == (result.run_id, str(result.output_dir))
        st.session_state.selected_trial = (remembered["trial"] if same_run and remembered.get("trial") in trial_options
                                           else trial_options[0])
    selected_trial = st.selectbox(
        tr('Candidat à explorer'),
        options=trial_options,
        format_func=candidate_labels.__getitem__,
        key="selected_trial",
    )
    st.session_state.results_candidate_selection = {"run": (result.run_id, str(result.output_dir)), "trial": selected_trial}
    candidate = candidates_by_trial[selected_trial]
    _render_candidate_navigation(trial_options, selected_trial)
    _render_switch_candidate(result, candidate)
    render_save_candidate(result, candidate)


def _render_optimization_status(result) -> None:
    trials = getattr(result, "optimization", {}).get("trials", [])
    if not trials:
        return
    reached = sum(item.get("stop_reason") == "target_reached" for item in trials)
    with st.expander(tr("Convergence des essais")):
        st.caption(tr("{reached}/{total} essais ont atteint les cibles locale et globale. Un calcul terminé n’implique pas que ces cibles ont été atteintes.", reached=reached, total=len(trials)))
        if result.optimization.get("rbs_stability"):
            st.caption(tr("L’affinage combine le défaut global, la stabilité RBS–linker, l’indice OSTIR et, en extension inverse, la linéarité du trigger, selon leurs poids. Les limites des défauts locaux restent prioritaires. Ce coût n’est pas une probabilité."))
        reasons = {"target_reached": "Cibles locale et globale atteintes", "budget_exhausted": "Limite de recherche atteinte",
                   "space_exhausted": "Espace de séquences entièrement exploré",
                   "no_mutable_sites": "Aucune base modifiable", "no_admissible_sequence": "Aucune séquence admissible trouvée dans le budget", "cancelled": "Annulé"}
        def stop_label(item):
            if item.get("stop_reason") == "budget_exhausted":
                if item.get("budget_kind") == "evaluations":
                    return tr("Limite d’évaluations atteinte")
                if item.get("budget_kind") == "proposals":
                    return tr("Limite de propositions atteinte")
            return tr(reasons.get(item.get("stop_reason"), item.get("stop_reason", "—")))
        from aptaswitch_studio.table_ui import render_readonly_table
        def refinement_cost(item):
            objective = item.get("best_objective") or {}
            if "rbs_linker_energy" not in objective:
                return None
            return objective["global_defect"] + objective["rbs_weight"] * (max(0., -objective["rbs_linker_energy"]) / 100.)
        render_readonly_table([{
            tr("Essai"): item.get("trial"),
            tr("Évaluations"): item.get("evaluations"),
            tr("Phase"): tr({"initialization": "Population initiale", "local": "Défauts locaux", "global": "Affinage du score"}.get(item.get("phase"), "—")),
            tr("Défaut local (%)"): None if item.get("best_loss") is None else 100 * item["best_loss"],
            tr("Défaut global (%)"): None if item.get("best_global_loss") is None else 100 * item["best_global_loss"],
            tr("Coût d’affinage"): refinement_cost(item),
            tr("Motif d’arrêt"): stop_label(item),
        } for item in trials])


def _render_switch_candidate(result, candidate) -> None:
    excluded = getattr(candidate, "excluded", False)
    if excluded and not candidate.structure_on:
        with st.container(border=True):
            st.error(tr('Candidat exclu · {exclusion_reason}', exclusion_reason=tr(candidate.exclusion_reason)))
            st.caption(tr('Les calculs de structure et d’association ont été ignorés pour ce candidat. Sa séquence reste consultable et exportée.'))
            st.markdown(tr('**Trigger 5′ → 3′**'))
            st.code(candidate.trigger_seq, language=None)
            st.markdown(tr('**Switch 5′ → 3′**'))
            st.code(candidate.switch_seq, language=None)
        return

    with st.container(border=True):
        if not excluded:
            summary = st.columns(3)
            score = candidate_score_100(candidate, result.request.score_scale)
            summary[0].metric(tr("Score /100"), "—" if score is None else f"{score:.2f}")
            summary[1].metric(tr("ON %"), "—" if candidate.on_yield_pct is None else f"{candidate.on_yield_pct:.2f}")
            summary[2].metric(tr(ASSOCIATION_WITH_APTAMER_LABEL), "—" if candidate.leak_pct is None else f"{candidate.leak_pct:.2f}",
                              help=association_with_aptamer_help(result.request.thermo, weight=result.request.scoring.leak))
        if excluded:
            st.error(tr('Candidat exclu · {exclusion_reason}', exclusion_reason=tr(candidate.exclusion_reason)))
            st.caption(tr('Structure ON conservée lors du criblage. Les autres analyses ont été ignorées ; la séquence et cette structure restent consultables et exportables.'))
            state = "ON + trigger"
        else:
            state = st.radio(tr('Structure'), ("OFF", "ON + trigger"), horizontal=True, key="candidate_state")
            if candidate.has_stop_codon:
                st.error(tr('Codon STOP prématuré — codon {first_stop}', first_stop=candidate.first_stop))
        try:
            structure = candidate_display_structure(candidate.structure_on if state.startswith("ON") else candidate.structure_off)
        except ValueError as exc:
            st.error(tr('Structure calculée indisponible : {exc}', exc=exc))
            return
        state_name = "on" if state.startswith("ON") else "off"
        structure_tab, linear_tab = st.tabs([tr(label) for label in MOLECULAR_VIEW_LABELS])
        with structure_tab:
            try:
                render_diagram(
                    partial(
                        candidate_structure_svg,
                        result.request.architecture, candidate.trigger_seq, candidate.switch_seq, structure,
                        state=state_name,
                        reporter_sequence=result.request.reporter.sequence_after_start_rna,
                        reporter_label=candidate.reporter_label or result.request.reporter.label,
                        width=1000, height=500,
                    ),
                    structure=structure, filename=f"candidate_{candidate.trial}_{state_name}", key="switch_candidate",
                )
            except Exception as exc:
                st.error(tr('Aperçu indisponible : {exc}', exc=exc))
        with linear_tab:
            try:
                render_diagram(
                    partial(
                        linear_sequence_svg,
                        result.request.architecture, candidate.trigger_seq,
                        result.request.reporter.sequence_after_start_rna,
                        candidate.reporter_label or result.request.reporter.label,
                        switch_sequence=candidate.switch_seq,
                    ),
                    structure=structure, filename=f"candidate_{candidate.trial}_lineaire", key="switch_candidate_linear",
                    download_labels=("Séquence linéaire SVG", "Séquence linéaire PNG"),
                )
            except ValueError as exc:
                st.error(tr('Séquence linéaire indisponible : {exc}', exc=exc))
        with st.expander(tr('Séquences du candidat')):
            st.markdown(tr('**Trigger 5′ → 3′**'))
            st.code(candidate.trigger_seq, language=None)
            st.markdown(tr('**Switch 5′ → 3′**'))
            st.code(candidate.switch_seq, language=None)

        from aptaswitch_studio.score_ui import render_candidate_score
        render_candidate_score(result, candidate)
        render_candidate_defects(candidate)
        from aptaswitch_studio.accessibility_ui import render_rbs_accessibility
        from aptaswitch_studio.rbs_calculator_ui import render_rbs_calculator
        render_rbs_calculator(candidate)
        render_rbs_accessibility(candidate)
        if not excluded:
            from aptaswitch_studio.ensemble_ui import render_structure_ensemble
            render_structure_ensemble(result, candidate)


def _select_ranked_candidate(trial: int) -> None:
    st.session_state.selected_trial = trial


def _render_candidate_navigation(trial_options: list[int], selected_trial: int) -> None:
    position = trial_options.index(selected_trial)
    previous, following = st.columns(2)
    previous.button(
        tr('← Candidat précédent'), key="previous_candidate", width="stretch",
        help=tr('Afficher le candidat précédent dans le classement.'),
        disabled=position == 0, on_click=_select_ranked_candidate,
        args=(trial_options[max(0, position - 1)],),
    )
    following.button(
        tr('Candidat suivant →'), key="next_candidate", width="stretch",
        help=tr('Afficher le candidat suivant dans le classement.'),
        disabled=position == len(trial_options) - 1, on_click=_select_ranked_candidate,
        args=(trial_options[min(len(trial_options) - 1, position + 1)],),
    )





def _render_reporter_settings() -> None:
    st.subheader(tr('Bibliothèque de reporters'), help=tr('Conservez le rapporteur complet. Le design utilise automatiquement ses 30 premières bases après AUG ; les analyses en contexte utilisent ensuite la séquence complète.'))
    if st.session_state.reporter_return_to_design:
        st.button(tr('← Retour à la conception'), on_click=_return_to_design, key="return_from_reporters", width="stretch")
    reporters = all_reporters()
    for reporter in reporters.values():
        cols = st.columns([0.2, 0.62, 0.18], vertical_alignment="center")
        cols[0].markdown(tr('**{label}**', label=reporter.label))
        full_sequence = reporter.full_sequence_after_start_rna
        cols[1].code("AUG" + full_sequence, language=None, wrap_lines=True)
        cols[1].caption(tr('Séquence codante complète · {count} nt', count=len(full_sequence) + 3))
        if reporter.source_url:
            cols[1].markdown(f"[NCBI {reporter.source_accession}]({reporter.source_url})")
        if is_builtin_reporter(reporter.key):
            cols[2].caption(tr('Intégré'))
        elif cols[2].button(tr('Supprimer'), key=f"remove_reporter_{reporter.key}", width="stretch"):
            remove_custom_reporter(reporter.key)
            if st.session_state.reporter_key == reporter.key:
                st.session_state.reporter_key = "sfgfp"
            st.rerun()

    with st.form("add_reporter_form", clear_on_submit=True):
        st.markdown(tr('#### Ajouter un reporter'))
        label = st.text_input(tr('Nom'), placeholder=tr('ex. mCherry'), key="new_reporter_label")
        sequence = st.text_area(tr('Séquence codante complète ARN/ADN, commençant par AUG/ATG'), height=140, key="new_reporter_sequence",
                                help=tr('Le premier AUG/ATG est retiré automatiquement : le switch possède déjà son codon de départ. Le STOP terminal, s’il est fourni, est conservé.'))
        submitted = st.form_submit_button(tr('Ajouter'), type="primary", width="stretch")
    if submitted:
        try:
            reporter = add_custom_reporter(label, sequence)
        except (ValueError, SequenceValidationError) as exc:
            st.error(tr(str(exc)))
        else:
            st.session_state.reporter_key = reporter.key
            st.session_state.last_reporter_key = reporter.key
            st.success(tr('Reporter {label} ajouté.', label=reporter.label))
            st.rerun()


def _render_general_settings() -> None:
    st.subheader(tr('Préférences locales'))
    st.selectbox(tr("Langue"), options=list(LANGUAGES), format_func=LANGUAGES.get,
                 key="ui_language", on_change=_save_language)
    if "settings_default_output" not in st.session_state:
        st.session_state.settings_default_output = str(default_design_root())
    with st.form("default_output_form"):
        st.text_input(tr("Dossier des designs"), key="settings_default_output",
                      help=tr("Chaque design possède son dossier avec ses séquences, ses calculs et ses structures 3D."))
        save_default = st.form_submit_button(tr('Enregistrer'), type="primary", width="stretch")
    if save_default:
        path = st.session_state.settings_default_output.strip()
        if not path:
            st.error(tr('Le dossier ne peut pas être vide.'))
        else:
            settings = load_settings()
            settings["design_library_dir"] = path
            save_settings(settings)
            st.success(tr("Dossier des designs enregistré. Les designs existants restent à leur emplacement actuel."))

    st.divider()
    st.markdown(tr('#### Réinitialisation'))
    st.caption(tr('Efface les reporters personnalisés et le dossier par défaut. Les exports ne sont pas supprimés.'))
    confirm = st.checkbox(tr('Je confirme la réinitialisation de tous les réglages'), key="confirm_reset_settings")
    st.button(
        tr('Réinitialiser les réglages'),
        disabled=not confirm,
        on_click=_reset_all_settings,
        width="stretch",
    )


def _save_language() -> None:
    set_language(st.session_state.ui_language)
    settings = load_settings()
    settings["ui_language"] = st.session_state.ui_language
    save_settings(settings)


def _remember_settings_tab() -> None:
    for label in SETTINGS_TABS:
        if st.session_state.settings_tab_label in (label, tr(label)):
            st.session_state.settings_tab = label
            break


def _render_settings_page() -> None:
    st.title(tr('Réglages'), anchor=False)
    if st.session_state.settings_tab not in SETTINGS_TABS:
        st.session_state.settings_tab = SETTINGS_TABS[0]
    st.session_state.settings_tab_label = tr(st.session_state.settings_tab)
    licenses, reporters, general = st.tabs(
        [tr(label) for label in SETTINGS_TABS],
        key="settings_tab_label", on_change=_remember_settings_tab,
        default=tr(st.session_state.settings_tab),
    )
    with licenses:
        render_licenses()
    with reporters:
        _render_reporter_settings()
    with general:
        _render_general_settings()


def main() -> None:
    st.set_page_config(
        page_title="Δpta2Switch-Studio",
        page_icon=str(APP_ICON_PATH) if APP_ICON_PATH.exists() else "🧬",
        layout="wide",
        initial_sidebar_state="expanded",

    )
    _initialize_session()
    if st.session_state.get("active_design_path"):
        st.session_state.output_dir = str(design_output_dir(Path(st.session_state.active_design_path)))
    elif st.session_state.nav not in ("Mes designs", "Réglages"):
        st.session_state.nav = "Mes designs"
    attach_saved_switch_jobs(st.session_state)
    _inject_styles()
    _render_header()
    page = _render_sidebar()
    if page in ("Conception", "Résultats", "Toehold switches enregistrés") and not st.session_state.get("active_design_path"):
        render_design_library()
    elif page == "Mes designs":
        render_design_library()
    elif page == "Conception":
        render_design_toolbar()
        _render_design_page()
    elif page == "Résultats":
        render_design_toolbar()
        _render_results_page()
    elif page == "Toehold switches enregistrés":
        render_design_toolbar()
        render_saved_switches_page()
    else:
        _render_settings_page()
    persist_active_design()
    if st.session_state.get("design_save_error"):
        st.error(tr("Le design n’a pas pu être sauvegardé : {error}", error=st.session_state.design_save_error))


if __name__ == "__main__":
    main()
