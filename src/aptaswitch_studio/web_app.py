"""Local Streamlit interface for Apta2Switch-Studio."""

from __future__ import annotations

import base64
import os
import subprocess
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
from aptaswitch_core.localization import LANGUAGES, get_language, set_language, tr

from aptaswitch_core.architecture import DEFAULT_TRIGGER_LENGTH, architecture_summary, default_architecture_for_trigger
from aptaswitch_core.models import DEFAULT_RBS_SEQUENCE, ReporterContext, RunRequest, ScoringWeights, SequenceInput, SwitchArchitecture, ThermoConditions, STOP_FILTER_DISABLED_WARNING
from aptaswitch_core.nupack_setup import (
    configured_nupack_import_path,
    detect_nupack,
    load_settings,
    nupack_drop_directory,
    save_nupack_settings,
    save_settings,
    validate_nupack_path,
)
from aptaswitch_core.preview import (
    DOMAIN_COLORS,
    TRIGGER_COLOR,
    build_switch_model,
    candidate_display_structure,
    candidate_structure_svg,
    linear_sequence_svg,
    switch_animation_frames,
    switch_preview_svg,
)
from aptaswitch_core.reporter_library import add_custom_reporter, all_reporters, is_builtin_reporter, remove_custom_reporter
from aptaswitch_core.sequences import SequenceValidationError, normalize_dna
from aptaswitch_studio.extension_ui import (
    EXTENSION_DEFAULTS,
    render_extension_controls,
    render_extension_results,
    render_ligand_selection,
    render_system_conditions,
    render_initial_structure,
    sync_ligand_selection,
)
from aptaswitch_studio.figure_ui import (
    MOLECULAR_VIEW_LABELS, ExportOption, render_diagram, render_export_menu,
    render_svg as _render_svg,
)
from aptaswitch_studio.architecture_ui import render_architecture_selector
from aptaswitch_studio.desktop import external_process_environment
from aptaswitch_studio.concentration_ui import (
    convert_concentration_display_unit as _convert_concentration_display_unit,
    render_concentration_control as _render_concentration_control,
)
from aptaswitch_studio.run_library_ui import render_run_library, render_saved_run_details
from aptaswitch_studio.web_runtime import (
    DesignRunJob,
    concentration_from_m,
    concentration_to_m,
    default_export_dir,
    ensure_writable_export_dir,
    exports_zip_bytes,
    largest_unbound_window,
    legacy_export_dir,
    normalized_export_dir,
    preferred_concentration_unit,
)


ASSETS = Path(__file__).resolve().parent / "assets"
LOGO_PATH = ASSETS / "logo.png"
APP_ICON_PATH = ASSETS / "app_icon.png"
SIDEBAR_LOGO_PATH = ASSETS / "igem_SU_logo.png"
NAV_ITEMS = ("Conception", "Résultats", "Réglages")
ADD_REPORTER = "__add_reporter__"
MATERIAL_OPTIONS = {
    "Switch ARN + trigger ADN (rna-dna06)": ("rna-dna06", "rna", "dna"),
    "Switch ADN + trigger ARN (rna-dna06)": ("rna-dna06", "dna", "rna"),
    "Switch ADN + trigger ADN (dna04)": ("dna04", "dna", "dna"),
    "Switch ARN + trigger ARN (rna06)": ("rna06", "rna", "rna"),
    "Switch ARN + trigger ARN (rna95)": ("rna95", "rna", "rna"),
    "Switch ARN + trigger ARN (rna99)": ("rna99", "rna", "rna"),
}
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
    if "customize_tswitch" not in st.session_state:
        # Existing browser sessions also start with the new standard activator.
        st.session_state.toehold_length = STANDARD_ACTIVATOR.toehold_length
    legacy_concentration_value = st.session_state.get("species_concentration_value", 5.0)
    legacy_concentration_unit = st.session_state.get("species_concentration_unit", "µM")
    defaults = {
        "nav": NAV_ITEMS[0],
        "molecule_name": "custom target",
        "aptamer_input": "",
        "trigger_input": "",
        "binding_start": 1,
        "binding_end": 1,
        "trigger_snapshot": None,
        "toehold_length": STANDARD_ACTIVATOR.toehold_length,
        "customize_tswitch": False,
        "t7_leader": "GGG",
        "rbs_loop": DEFAULT_RBS_SEQUENCE,
        "linker_auto": True,
        "frame_linker": "",
        "rbs_aug_distance": STANDARD_ACTIVATOR.rbs_to_aug_distance,
        "material_label": next(iter(MATERIAL_OPTIONS)),
        "trials": 50,
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
        "settings_tab": "Dépendances",
        "reporter_return_to_design": False,
        "weight_on": DEFAULT_SCORING.on_yield,
        "weight_leak": DEFAULT_SCORING.leak,
        "weight_defect": DEFAULT_SCORING.defect,
        "weight_ddg": DEFAULT_SCORING.ddg,
        "weight_rbs": DEFAULT_SCORING.rbs,
        "engine": "nupack",
        "output_dir": _default_output_dir(),
        "aptamer_indices": [],
        "aptamer_status": "",
        "run_error": "",
        "run_notice": "",
        "sequence_action_error": "",
        "long_run_confirm": False,
        "preview_state": "OFF",
        "last_result": None,
        "last_exports": {},
        "design_job": None,
        "animation_frame": 0,
        "animation_playing": False,
        "results_kind": "Design de switches",
        **EXTENSION_DEFAULTS,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

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

    scoring_revision = st.session_state.get("scoring_defaults_revision", 0)
    if scoring_revision < 2:
        previous_defaults = [
            {"weight_on": 1.0, "weight_leak": 3.0, "weight_defect": 6.0,
             "weight_ddg": 0.25, "weight_rbs": 6.0},
        ]
        if scoring_revision < 1:
            previous_defaults.append({"weight_on": 1.0, "weight_leak": 3.0, "weight_defect": 3.0,
                                      "weight_ddg": 0.5, "weight_rbs": 4.0})
        if any(all(st.session_state[key] == value for key, value in weights.items())
               for weights in previous_defaults):
            for key in previous_defaults[0]:
                st.session_state[key] = defaults[key]
        st.session_state.scoring_defaults_revision = 2

    # Retire demonstration runs already present in an open browser session.
    st.session_state.engine = "nupack"
    old_job = st.session_state.get("design_job")
    if old_job is not None and old_job.request.engine != "nupack":
        old_job.cancel()
        st.session_state.design_job = None
    old_result = st.session_state.get("last_result")
    if old_result is not None and old_result.request.engine != "nupack":
        st.session_state.last_result = None
        st.session_state.last_exports = {}
    old_kinetic_job = st.session_state.pop("ms_job", None)
    if old_kinetic_job is not None:
        old_kinetic_job.cancel()

    # Streamlit normally deletes a widget's key when that widget disappears.
    # Touch the page-scoped values on every rerun so navigating to another page
    # does not erase a completed design form.
    extra_keys = {
        "settings_nupack_path",
        "settings_default_output",
        "selected_trial",
        "candidate_state",
        "confirm_reset_settings",
        "saved_run_path",
        "run_open_path",
        "ui_language",
        "new_reporter_label",
        "new_reporter_sequence",
    }
    dynamic_prefixes = ("extension_selected_", "extension_ligand_base_")
    for key in list(st.session_state):
        if key in defaults or key in extra_keys or key.startswith(dynamic_prefixes):
            st.session_state[key] = st.session_state[key]


def _auto_detect_nupack_once() -> None:
    if st.session_state.get("dependency_detection_done"):
        return
    st.session_state.dependency_detection_done = True
    if configured_nupack_import_path():
        return
    detected = detect_nupack()
    if detected and detected.ok:
        save_nupack_settings(detected, True)


def _inject_styles() -> None:
    st.markdown(
        """
        <style>
        :root { --apta-blue:#1258dc; --apta-navy:#0d1f3c; --apta-cyan:#10a7c8; }
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
        [data-testid="stExpandSidebarButton"] {
            visibility:visible; pointer-events:auto;
            position:fixed; top:.5rem; left:.5rem;
            background:#f8fbff; border-radius:.5rem;
        }
        [data-testid="stSidebar"] { background:#0d1f3c; }
        [data-testid="stSidebar"] * { color:#edf5ff; }
        [data-testid="stSidebar"] [data-testid="stImage"] {
            width:fit-content; background:#ffffff; padding:.5rem;
            border-radius:1rem; box-shadow:0 6px 18px rgba(0,0,0,.22);
        }
        [data-testid="stSidebar"] [data-testid="stImage"] img {
            display:block; border-radius:.65rem;
        }
        [data-testid="stSidebar"] [data-testid="stRadio"] label {
            padding:.35rem .45rem; border-radius:.55rem;
        }
        .block-container, [data-testid="stMainBlockContainer"] {
            max-width:1440px; padding-top:1.4rem; padding-bottom:4rem;
        }
        h1, h2, h3 { color:#0d1f3c; letter-spacing:-.025em; }
        h1 { font-size:2.25rem !important; }
        div[data-testid="stVerticalBlockBorderWrapper"] {
            background:rgba(255,255,255,.88); border:1px solid #dbe5f0;
            border-radius:1rem; box-shadow:0 8px 24px rgba(13,31,60,.045);
        }
        .stButton > button, .stDownloadButton > button { border-radius:.65rem; font-weight:650; }
        .stButton > button[kind="primary"] { background:#1258dc; border-color:#1258dc; }
        .apta-subtitle { color:#55657a; font-size:1.02rem; margin-top:-.7rem; margin-bottom:1.2rem; }
        .status-pill { display:inline-block; padding:.26rem .62rem; margin:.1rem .2rem .1rem 0;
            border-radius:999px; font-size:.78rem; font-weight:700; }
        .status-ok { background:#dcfce7; color:#166534; }
        .status-muted { background:#e8eef5; color:#475569; }
        .status-warn { background:#fff3cd; color:#805b00; }
        .legend-chip { display:inline-flex; align-items:center; gap:.3rem; margin:.15rem .7rem .15rem 0;
            color:#425269; font-size:.82rem; }
        .legend-dot { width:.65rem; height:.65rem; border-radius:50%; display:inline-block; }
        code { white-space:pre-wrap !important; overflow-wrap:anywhere; }
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
    logo_col, title_col = st.columns([0.12, 0.88], vertical_alignment="center")
    with logo_col:
        if LOGO_PATH.exists():
            _render_brand_image(LOGO_PATH, width=112, alt="Logo Apta2Switch-Studio")
    with title_col:
        st.markdown("# Δpta2Switch-Studio")
        st.markdown(
            tr('<div class="apta-subtitle">Concevoir, classer et explorer des toehold switches activés par aptamère.</div>'),
            unsafe_allow_html=True,
        )


def _render_sidebar() -> str:
    with st.sidebar:
        if SIDEBAR_LOGO_PATH.exists():
            _render_brand_image(SIDEBAR_LOGO_PATH, width=120, alt="Logo iGEM Sorbonne Université", backdrop=True)
        st.markdown("### Δpta2Switch")
        page = st.radio(tr('Navigation'), NAV_ITEMS, key="nav", label_visibility="collapsed", format_func={label: tr(label) for label in NAV_ITEMS}.__getitem__)
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
    previous_length = len(previous) if previous else 0
    current_length = len(trigger)
    st.session_state.trigger_snapshot = trigger
    st.session_state.aptamer_indices = []
    st.session_state.aptamer_status = ""
    if previous_length != current_length:
        custom = bool(st.session_state.customize_tswitch)
        entire_previous_trigger = (st.session_state.binding_start == 1
                                   and st.session_state.binding_end == max(1, previous_length))
        if not custom or not previous_length or entire_previous_trigger:
            st.session_state.binding_start = 1
            st.session_state.binding_end = max(1, current_length)
        else:
            st.session_state.binding_start = min(st.session_state.binding_start, max(1, current_length))
            st.session_state.binding_end = max(st.session_state.binding_start,
                                               min(st.session_state.binding_end, max(1, current_length)))
        if current_length >= 2:
            st.session_state.toehold_length = (
                min(max(1, int(st.session_state.toehold_length)), current_length - 1) if custom
                else default_architecture_for_trigger(current_length).toehold_length
            )


def _effective_trigger(trigger: str) -> str:
    if not trigger:
        return ""
    if not st.session_state.customize_tswitch:
        return trigger
    start = max(1, min(int(st.session_state.binding_start), len(trigger)))
    end = max(start, min(int(st.session_state.binding_end), len(trigger)))
    return trigger[start - 1 : end]


def _clean_rna(text: str, *, allow_free_bases: bool = False) -> str:
    alphabet = "ACGUN" if allow_free_bases else "ACGU"
    return "".join(char for char in (text or "").upper().replace("T", "U") if char in alphabet)


def _current_architecture(trigger_length: int) -> SwitchArchitecture:
    custom = bool(st.session_state.customize_tswitch)
    minimum_length = 2 if custom else STANDARD_ACTIVATOR.trigger_region_length
    design_length = trigger_length if trigger_length >= minimum_length else STANDARD_ACTIVATOR.trigger_region_length
    default_architecture = default_architecture_for_trigger(design_length)
    toehold = (max(1, min(int(st.session_state.toehold_length), design_length - 1)) if custom
               else default_architecture.toehold_length)
    rbs = _clean_rna(st.session_state.rbs_loop) or DEFAULT_RBS_SEQUENCE
    frame_linker = _clean_rna(st.session_state.frame_linker, allow_free_bases=True) if custom and not st.session_state.linker_auto else None
    return SwitchArchitecture(
        toehold_length=toehold,
        stem_length=design_length - toehold,
        upper_stem_length=int(st.session_state.rbs_aug_distance) if custom else STANDARD_ACTIVATOR.rbs_to_aug_distance,
        upper_stem2_length=default_architecture.upper_stem2_length,
        rbs_loop=rbs,
        t7_leader=_clean_rna(st.session_state.t7_leader),
        frame_linker=frame_linker,
        aug_spacer=0,
        rbs_prefix_length=STANDARD_ACTIVATOR.rbs_prefix_length,
        architecture_name=(f"toehold_activator_custom_{design_length}nt_toehold{toehold}_stem{design_length - toehold}"
                           if custom else f"toehold_activator_standard_{design_length}nt"),
    )


def _current_material() -> tuple[str, str, str]:
    return MATERIAL_OPTIONS.get(st.session_state.material_label, next(iter(MATERIAL_OPTIONS.values())))


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
    material, switch_material, trigger_material = _current_material()
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
        magnesium_m=concentration_to_m(
            st.session_state.magnesium_concentration_value,
            st.session_state.magnesium_concentration_unit,
        ),
        concentration_m=trigger_concentration_m,
        trigger_concentration_m=trigger_concentration_m,
        switch_concentration_m=switch_concentration_m,
        material=material,
        switch_material=switch_material,
        trigger_material=trigger_material,
    )


def _current_reporter() -> ReporterContext:
    reporters = all_reporters()
    key = st.session_state.reporter_key
    if key not in reporters:
        key = "sfgfp"
        st.session_state.reporter_key = key
    return reporters[key]


def _build_request() -> RunRequest:
    aptamer = normalize_dna(st.session_state.aptamer_input, field_name="aptamer")
    trigger = normalize_dna(st.session_state.trigger_input, field_name="trigger")
    if not st.session_state.customize_tswitch and len(trigger) < STANDARD_ACTIVATOR.trigger_region_length:
        raise ValueError(tr('Le toehold activateur standard nécessite un trigger d’au moins {trigger_region_length} nt. Étendez le trigger en section 2 ou activez « Je veux personnaliser mon tSwitch ».', trigger_region_length=STANDARD_ACTIVATOR.trigger_region_length))
    effective = _effective_trigger(trigger)
    if len(effective) < 2:
        raise ValueError(tr('La région du trigger utilisée pour le switch doit contenir au moins 2 bases.'))
    architecture = _current_architecture(len(effective))
    reporter = _current_reporter()
    nupack_path = configured_nupack_import_path()
    if not nupack_path:
        raise ValueError(tr('Configurez NUPACK dans Réglages pour lancer la conception.'))
    return RunRequest(
        sequences=SequenceInput(
            aptamer_dna=aptamer,
            trigger_dna=trigger,
            molecule_name=st.session_state.molecule_name,
            binding_start=int(st.session_state.binding_start) - 1 if st.session_state.customize_tswitch else 0,
            binding_length=len(effective),
        ),
        architecture=architecture,
        thermo=_current_thermo(),
        reporter=ReporterContext(
            key=reporter.key,
            label=reporter.label,
            sequence_after_start_rna=reporter.sequence_after_start_rna,
        ),
        scoring=ScoringWeights(
            on_yield=float(st.session_state.weight_on),
            leak=float(st.session_state.weight_leak),
            defect=float(st.session_state.weight_defect),
            ddg=float(st.session_state.weight_ddg),
            rbs=float(st.session_state.weight_rbs),
        ),
        trials=int(st.session_state.trials),
        output_dir=normalized_export_dir(st.session_state.output_dir),
        engine="nupack",
        nupack_path=nupack_path,
        exclude_stop_candidates=bool(st.session_state.exclude_stop_candidates),
        exclude_non_linear_candidates=bool(st.session_state.exclude_non_linear_candidates),
    )


def _exclude_aptamer_region() -> None:
    trigger, error = _normalize_optional_sequence(st.session_state.trigger_input, "trigger")
    if error:
        st.session_state.sequence_action_error = error
        return
    window = largest_unbound_window(len(trigger), st.session_state.aptamer_indices)
    if window is None:
        st.session_state.sequence_action_error = "Aucune région libre d'au moins 2 nt n'a été trouvée."
        return
    start, length = window
    st.session_state.binding_start = start + 1
    st.session_state.binding_end = start + length
    st.session_state.sequence_action_error = ""


def _start_run() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.run_error = ""
    st.session_state.run_notice = ""
    for job_key in ("design_job", "extension_job"):
        active_job = st.session_state.get(job_key)
        if active_job is not None and active_job.snapshot().state == "running":
            st.session_state.run_error = "Un calcul est déjà en cours. Suivez-le ou annulez-le dans Résultats."
            return
    try:
        request = _build_request().normalized()
    except Exception as exc:
        st.session_state.run_error = str(exc)
        return
    requested_output = request.output_dir
    try:
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
    job = DesignRunJob(request)
    st.session_state.design_job = job
    st.session_state.last_result = None
    st.session_state.last_exports = {}
    st.session_state.results_kind = "Design de switches"
    st.session_state.nav = "Résultats"
    job.start()


def _cancel_design_run() -> None:
    job = st.session_state.get("design_job")
    if job is not None:
        job.cancel()


def _reset_all_settings() -> None:
    save_settings({"ui_language": st.session_state.ui_language})
    for key in (
        "settings_nupack_path",
        "settings_multistrand_path",
        "settings_default_output",
    ):
        st.session_state.pop(key, None)
    st.session_state.dependency_detection_done = True
    st.session_state.reporter_key = "sfgfp"
    st.session_state.output_dir = str(default_export_dir())
    st.session_state.engine = "nupack"


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


@st.cache_data(show_spinner=False, max_entries=12)
def _cached_animation_frames(
    architecture: SwitchArchitecture,
    trigger: str,
    reporter_sequence: str,
    reporter_label: str,
    switch_material: str,
    trigger_material: str,
    language: str,
) -> list[str]:
    """Cache horizontal ON frames per language with free linker bases shown as N."""
    set_language(language)
    return switch_animation_frames(
        architecture,
        trigger or None,
        reporter_sequence,
        reporter_label,
        switch_material=switch_material,
        trigger_material=trigger_material,
        width=1000,
        height=500,
        show_annotations=False,
    )


@st.fragment(run_every=0.18)
def _animation_player(frames: list[str]) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    if not frames:
        st.info(tr("Ajoutez un trigger valide pour générer l'animation."))
        return
    if st.session_state.animation_playing:
        st.session_state.animation_frame = (int(st.session_state.animation_frame) + 1) % len(frames)
    controls = st.columns([0.35, 0.65], vertical_alignment="center")
    with controls[0]:
        st.toggle(tr('Lecture automatique'), key="animation_playing")
    with controls[1]:
        st.slider(tr('Étape'), 0, len(frames) - 1, key="animation_frame", label_visibility="collapsed")
    _render_svg(frames[int(st.session_state.animation_frame)], height=510)


def _render_sequence_and_binding() -> tuple[str, str]:
    with st.container(border=True):
        st.subheader(tr('1 · Système initial : données et structure 2D'))
        st.caption(tr('Renseignez les séquences disponibles, les conditions du système et, si elles sont connues, les bases impliquées dans la liaison au ligand.'))
        left, right = st.columns(2, gap="large", vertical_alignment="top")
        with left:
            st.markdown(tr('#### Séquences du système'))
            st.text_input(tr('Molécule cible'), key="molecule_name")
            st.text_area(
                tr("Séquence de l'aptamère"), key="aptamer_input", height=110,
                placeholder=tr('Séquence ADN ou ARN — espaces et retours à la ligne acceptés'),
            )
            st.text_area(tr('Séquence du trigger'), key="trigger_input", height=110,
                         placeholder=tr('Séquence ADN ou ARN'))
            trigger, trigger_error = _normalize_optional_sequence(st.session_state.trigger_input, "trigger")
            aptamer, aptamer_error = _normalize_optional_sequence(st.session_state.aptamer_input, "aptamer")
            for raw, error in ((st.session_state.aptamer_input, aptamer_error), (st.session_state.trigger_input, trigger_error)):
                if raw.strip() and error:
                    st.warning(tr(error))
            if aptamer or trigger:
                st.caption(tr('Aptamère : {value} nt · Trigger : {value1} nt', value=len(aptamer), value1=len(trigger)))
            if trigger and len(trigger) < 24:
                st.warning(tr('Attention : votre trigger mesure {value} nt. Sa longueur peut être trop faible ; vous pouvez l’étendre dans la section 2.', value=len(trigger)))
        with right:
            render_system_conditions()
        _sync_trigger_state(trigger)
        sync_ligand_selection(aptamer)
        with st.expander(tr('Bases de l’aptamère impliquées dans la liaison au ligand (si connues)'), expanded=False):
            render_ligand_selection(aptamer)
        render_initial_structure(aptamer, trigger)
    return aptamer, trigger


def _render_switch_binding_controls(trigger: str) -> None:
    st.markdown(tr('#### Région du trigger utilisée pour le switch'))
    max_index = max(1, len(trigger))
    st.session_state.binding_start = min(max(1, int(st.session_state.binding_start)), max_index)
    st.session_state.binding_end = min(max(st.session_state.binding_start, int(st.session_state.binding_end)), max_index)
    columns = st.columns([1, 1, 1.2], vertical_alignment="bottom")
    columns[0].number_input(tr('Première base'), min_value=1, max_value=max_index, step=1,
                            key="binding_start", disabled=not trigger)
    columns[1].number_input(tr('Dernière base'), min_value=1, max_value=max_index, step=1,
                            key="binding_end", disabled=not trigger)
    columns[2].button(tr('Exclure la zone liée'), on_click=_exclude_aptamer_region,
                      disabled=not st.session_state.aptamer_indices, width="stretch")
    st.caption(tr('La zone liée est identifiée dans la structure du système initial, calculée en section 1.'))
    if st.session_state.get("sequence_action_error"):
        st.error(tr(st.session_state.sequence_action_error))


def _render_architecture(trigger: str) -> tuple[SwitchArchitecture, str]:
    with st.container(border=True):
        st.subheader(tr('3 · Architecture du switch'))
        render_architecture_selector()
        sequence_columns = st.columns(2)
        sequence_columns[0].text_input(tr('Leader 5′'), key="t7_leader")
        sequence_columns[1].text_input(
            tr('RBS'), key="rbs_loop",
            help=tr('La boucle cible contient 3 bases libres (NNN), optimisées par NUPACK et non appariées, suivies de votre RBS.'),
        )
        st.selectbox(tr('Chimie / modèle'), options=list(MATERIAL_OPTIONS), key="material_label",
                     format_func={label: tr(label) for label in MATERIAL_OPTIONS}.__getitem__)
        custom = st.checkbox(tr('Je veux personnaliser mon tSwitch'), key="customize_tswitch")
        if custom:
            st.warning(tr('La personnalisation de l’architecture est déconseillée : elle peut modifier le repliement et le fonctionnement du switch. Conservez les paramètres standards sauf besoin particulier.'))
            _render_switch_binding_controls(trigger)
            effective = _effective_trigger(trigger)
            design_length = len(effective) if len(effective) >= 2 else STANDARD_ACTIVATOR.trigger_region_length
            st.session_state.toehold_length = max(1, min(int(st.session_state.toehold_length), design_length - 1))
            st.slider(
                tr('Longueur du toehold'),
                min_value=1,
                max_value=design_length - 1,
                key="toehold_length",
                help=tr('Le reste de la région du trigger forme le stem inférieur.'),
            )
            architecture = _current_architecture(len(effective))
            st.info(tr(architecture_summary(architecture)))
            toe_pct = 100 * architecture.toehold_length / architecture.trigger_region_length
            st.markdown(
                tr('\n                <div style="height:42px;display:flex;border-radius:8px;overflow:hidden;color:white;font-weight:700;">\n                  <div style="width:{toe_pct:.2f}%;background:#1258dc;padding:.65rem;white-space:nowrap">Toehold · {toehold_length} nt</div>\n                  <div style="flex:1;background:#0e9f76;padding:.65rem;white-space:nowrap">Stem · {stem_length} nt</div>\n                </div>\n                ', toe_pct=toe_pct, toehold_length=architecture.toehold_length, stem_length=architecture.stem_length),
                unsafe_allow_html=True,
            )
            st.number_input(
                tr('Distance entre le RBS et AUG (nt)'), min_value=1, max_value=30, step=1,
                key="rbs_aug_distance",
                help=tr('Nombre total de bases N entre la fin du RBS et AUG. Ces bases, optimisées par NUPACK, forment le bras droit de la tige supérieure dans la cible OFF.'),
            )
            link_cols = st.columns([0.45, 0.55], vertical_alignment="bottom")
            link_cols[0].toggle(
                tr('Linker automatique'), key="linker_auto",
                help=tr('La longueur est calculée pour préserver le cadre de lecture (0, 1 ou 2 bases). NUPACK choisit ces bases pour se rapprocher des structures cibles.'),
            )
            link_cols[1].text_input(
                tr('Linker personnalisé'),
                key="frame_linker",
                disabled=st.session_state.linker_auto,
                help=tr('N laisse une base libre pour NUPACK. A, C, G, U ou T imposent une base si vous souhaitez contraindre le linker.'),
            )
        else:
            effective = _effective_trigger(trigger)
            if trigger and len(trigger) < STANDARD_ACTIVATOR.trigger_region_length:
                st.info(tr('Le modèle standard nécessite un trigger d’au moins {trigger_region_length} nt. Étendez-le en section 2 ou personnalisez le tSwitch.', trigger_region_length=STANDARD_ACTIVATOR.trigger_region_length))
        architecture = _current_architecture(len(effective))

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
    _, switch_material, trigger_material = _current_material()
    model = build_switch_model(
        architecture, effective_trigger, reporter.sequence_after_start_rna,
        reporter.label, switch_material, trigger_material,
    )
    with st.container(border=True):
        st.subheader(tr('4 · Aperçus du switch en direct'))
        st.caption(tr('Structures cibles construites à partir de vos séquences. N indique les bases à optimiser par NUPACK ; les structures prédites seront disponibles après le calcul.'))
        structure_tab, linear_tab, animation_tab = st.tabs(
            [*(tr(label) for label in MOLECULAR_VIEW_LABELS), tr('Animation de liaison')]
        )
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
                switch_material=switch_material,
                trigger_material=trigger_material,
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
                switch_material,
                trigger_material,
            )
            render_diagram(
                diagram, filename="switch_lineaire", key="switch_preview_linear",
                caption="N indique une base laissée libre pour l’optimisation NUPACK.",
                download_labels=("Séquence linéaire SVG", "Séquence linéaire PNG"),
            )
        with animation_tab:
            make_animation = st.checkbox(
                tr("Générer l'animation schématique"),
                value=False,
                help=tr('La génération peut prendre quelques secondes pour les séquences longues.'),
            )
            if make_animation:
                with st.spinner(tr('Préparation des images…')):
                    frames = _cached_animation_frames(
                        architecture,
                        effective_trigger,
                        reporter.sequence_after_start_rna,
                        reporter.label,
                        switch_material,
                        trigger_material,
                        get_language(),
                    )
                if st.session_state.animation_frame >= len(frames):
                    st.session_state.animation_frame = 0
                _animation_player(frames)
            else:
                st.caption(tr("Activez l'animation pour suivre l'accrochage au toehold puis l'ouverture du stem."))


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
            tr('**Accessibilité de la région RBS–linker.** Énergie libre de repliement MFE, en kcal/mol, calculée par NUPACK sur cette région isolée : de la première base de la boucle contenant le RBS jusqu’à la base précédant le premier nucléotide du rapporteur. Les trois bases optimisées avant le RBS (NNN dans la cible) sont incluses ; la tige située avant la boucle et le rapporteur sont exclus. La région inclut le RBS, l’espace RBS–AUG, AUG, le bras aval du stem et le linker éventuel. Un ΔG moins négatif, proche de 0, est favorisé : moins de structure secondaire à défaire pour la traduction. Ce poids ajoute w × ΔG au score brut : gagner 1 kcal/mol ajoute w points. [Définition de Green, tableau S3](https://yin.hms.harvard.edu/publications/2014.tswitch1.table.s3.xlsx).')
        ),
        "weight_defect": (
            tr('**Fidélité aux cibles OFF et ON.** Le défaut d’ensemble normalisé retourné par NUPACK mesure l’écart aux structures et aux concentrations cibles des tubes de design : switch seul en OFF ; trigger + switch en ON. Il tient compte des mauvais appariements et du manque de complexes aux concentrations visées. Plus faible est meilleur ; 0 est idéal. Un défaut de 0,05 correspond à 5 %. Le score retire w × (100 × défaut) : réduire le défaut d’un point de pourcentage apporte w points au score brut. Ce n’est pas une mesure d’expression.')
        ),
        "weight_leak": (
            tr('**Fuite : critère secondaire de départage.** Dans votre protocole, les aptamères sont dans un autre tube. Le calcul ci-dessous réunit les trois espèces : c’est un scénario de comparaison, pas une simulation de la fuite ou du transfert entre vos tubes. Il donne le pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre, malgré la présence de l’aptamère censé retenir le trigger. Le tube contient aptamère ({trigger_value}), trigger ({trigger_value}) et switch ({switch_value}), soit un ratio aptamère:trigger:switch de {ratio}:{ratio}:1. La concentration d’aptamère est actuellement égale à celle du trigger. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus faible est meilleur dans ce scénario secondaire. Le score retire w × fuite ; réduire la fuite d’un point de pourcentage apporte w points au score. Par défaut, le poids de 0,05 limite l’écart dû à la fuite à 5 points sur toute la plage 0–100 %, pour départager surtout des candidats très proches. C’est un indicateur d’association, pas une mesure de traduction ; le ligand n’est pas simulé.', trigger_value=trigger_value, switch_value=switch_value, ratio=ratio)
        ),
        "weight_on": (
            tr('**Association ON avec le trigger libre.** Pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre. Le tube contient seulement le trigger libre ({trigger_value}) et le switch ({switch_value}), soit un ratio trigger:switch de {ratio}:1, sans aptamère ni ligand. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus élevé est meilleur pour favoriser la capture du trigger dans le tube de switch, séparé de celui des aptamères. Ce critère est prioritaire sur la fuite : par défaut, un point de pourcentage de ON pèse 60 fois plus qu’un point de fuite. Le score ajoute w × ON ; gagner un point de pourcentage apporte w points au score brut. Ce pourcentage est un indicateur du complexe ON ; il ne mesure pas directement l’activation de la traduction ni la quantité de protéine produite.', trigger_value=trigger_value, switch_value=switch_value, ratio=ratio)
        ),
        "weight_ddg": (
            tr('**Énergie de formation du complexe trigger–switch.** ΔΔG = ΔG MFE du complexe ON − ΔG MFE du switch seul − ΔG MFE du trigger seul, en kcal/mol, avec les conditions choisies. Une valeur plus négative favorise énergétiquement l’association. Le score retire w × ΔΔG : une diminution de 1 kcal/mol apporte w points au score brut. Son poids est faible par défaut car ce critère complète le rendement en tube sans décrire à lui seul le fonctionnement. Ce n’est ni une barrière d’activation ni une vitesse de réaction.')
        ),
    }


def _on_reporter_selected() -> None:
    if st.session_state.reporter_key == ADD_REPORTER:
        previous = st.session_state.get("last_reporter_key", "sfgfp")
        st.session_state.reporter_key = previous if previous in all_reporters() else "sfgfp"
        st.session_state.nav = "Réglages"
        st.session_state.settings_tab = "Reporters"
        st.session_state.reporter_return_to_design = True
    else:
        st.session_state.last_reporter_key = st.session_state.reporter_key


def _return_to_design() -> None:
    st.session_state.nav = "Conception"
    st.session_state.reporter_return_to_design = False


def _render_run_configuration() -> None:
    with st.container(border=True):
        st.subheader(tr('5 · Calcul et scoring du switch'))
        conditions_tab, scoring_tab, output_tab = st.tabs([tr('Conditions'), tr('Scoring'), tr('Sortie')])
        with conditions_tab:
            nupack_ready = configured_nupack_import_path() is not None
            if nupack_ready:
                st.caption(tr('Moteur de calcul : NUPACK'))
            else:
                st.info(tr('Configurez NUPACK dans Réglages pour lancer la conception.'))
            cols = st.columns(3)
            cols[0].number_input(tr('Essais'), min_value=1, max_value=100000, step=1, key="trials")
            cols[1].number_input(
                tr('Température (°C)'), min_value=-20.0, max_value=100.0, step=0.1, key="temperature_c"
            )
            reporters = all_reporters()
            if st.session_state.reporter_key not in reporters:
                st.session_state.reporter_key = "sfgfp"
            cols[2].selectbox(
                tr('Reporter'),
                options=[*reporters, ADD_REPORTER],
                format_func={**{key: reporter.label for key, reporter in reporters.items()},
                             ADD_REPORTER: tr("Ajouter un reporter…")}.__getitem__,
                key="reporter_key",
                on_change=_on_reporter_selected,
            )
            st.markdown(tr('##### Concentrations'))
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
                    min_m=0.0,
                    max_m=2.0,
                )
            with concentration_cols[3]:
                _render_concentration_control(
                    tr('Magnésium'),
                    value_key="magnesium_concentration_value",
                    unit_key="magnesium_concentration_unit",
                    min_m=0.0,
                    max_m=1.0,
                )
            st.caption(
                tr('Le rendement ON est calculé par rapport à la concentration initiale du toehold switch.')
            )
        with scoring_tab:
            st.caption(tr('Priorités par défaut : accessibilité RBS–linker et fidélité aux cibles, puis association ON et ΔΔG. La fuite a un poids très faible pour départager les candidats proches, les aptamères étant dans un autre tube. Les poids restent modifiables ; 0 désactive un critère.'))
            scoring_help = _scoring_help(_current_thermo())
            cols = st.columns(5)
            for column, label, key in zip(cols, (
                "ΔG RBS–linker", "Defect · malus", "ON % · bonus", "ΔΔG · liaison", "Fuite · départage",
            ), ("weight_rbs", "weight_defect", "weight_on", "weight_ddg", "weight_leak")):
                column.number_input(tr(label), 0.0, 20.0, step=0.05, key=key, help=tr(scoring_help[key]))
            with st.expander(tr('Voir la formule')):
                st.latex(
                    r"score = w_{ON}\,ON - w_{leak}\,leak - w_{defect}(100\,defect) "
                    r"- w_{\Delta\Delta G}\,\Delta\Delta G + w_{RBS}\,\Delta G_{RBS}"
                )
                st.caption(tr('Le score conserve la somme pondérée, sans plafond ni plancher : plus il est élevé, mieux le candidat est classé. Il peut dépasser 100 ou être négatif ; ce n’est pas un pourcentage. Les poids sont des coefficients de classement, pas des pourcentages d’importance ni des coefficients publiés par Green.'))
            st.checkbox(tr('Exclure les candidats avec un STOP prématuré'), key="exclude_stop_candidates")
            if st.session_state.exclude_stop_candidates:
                st.caption(tr('Dès la séquence générée, les candidats avec un STOP avant le reporter sont exclus. Ils restent en fin de tableau, sans score ni analyses thermodynamiques supplémentaires.'))
            else:
                st.warning(tr(STOP_FILTER_DISABLED_WARNING))
            st.checkbox(
                tr('Exclure les candidats sans linéarité parfaite du RBS–linker en ON'),
                key="exclude_non_linear_candidates",
                help=(
                    tr('Toutes les bases, de la première base de la boucle contenant le RBS (y compris les bases optimisées avant le RBS) jusqu’à la dernière base avant le rapporteur, doivent être non appariées dans les structures MFE ON renvoyées par NUPACK pour le complexe trigger–switch. Un seul appariement dans cette région suffit à exclure le candidat. Ce critère porte sur ces structures prédites, pas sur l’ensemble des conformations.')
                ),
            )
            if st.session_state.exclude_non_linear_candidates:
                st.caption(tr('Criblage de la structure ON après le filtre STOP, avant les autres analyses. Les candidats exclus restent en fin de tableau, sans score ni calcul de rendement ou de fuite.'))
        with output_tab:
            st.text_input(
                tr('Dossier local des runs'),
                key="output_dir",
                help=tr('Les designs et extensions sont rangés dans les sous-dossiers designs/ et extensions/. Vous pouvez les rouvrir dans Résultats.'),
            )
            st.caption(tr('Exports : JSON, CSV, XLSX, SVG, journal et projet .aptaswitch.json.'))

        if int(st.session_state.trials) > 50:
            st.warning(tr('Avec NUPACK, plus de 50 essais peuvent prendre plusieurs heures ou plusieurs jours.'))
            st.checkbox(tr('Je confirme ce calcul long'), key="long_run_confirm")
        can_start = int(st.session_state.trials) <= 50 or st.session_state.get("long_run_confirm", False)
        try:
            _build_request().normalized()
        except ValueError:
            can_start = False
        current_job = st.session_state.get("design_job")
        running = current_job is not None and current_job.snapshot().state == "running"
        extension_job = st.session_state.get("extension_job")
        running = running or (extension_job is not None and extension_job.snapshot().state == "running")
        st.button(
            tr('Lancer la conception'),
            type="primary",
            on_click=_start_run,
            disabled=not can_start or running,
            width="stretch",
        )
        if st.session_state.run_error:
            st.error(tr(st.session_state.run_error))
        if st.session_state.run_notice:
            st.info(tr(st.session_state.run_notice))


def _render_design_page() -> None:
    st.markdown(tr('## Conception'))
    st.markdown(
        tr('<div class="apta-subtitle">Définissez les séquences, contrôlez l’architecture puis lancez le classement.</div>'),
        unsafe_allow_html=True,
    )
    aptamer, trigger = _render_sequence_and_binding()
    render_extension_controls(aptamer, trigger)
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
    total = max(1, snapshot.total)
    metrics = st.columns(4)
    metrics[0].metric(tr('Moteur'), job.request.engine.upper())
    metrics[1].metric(tr('Progression'), tr('{current} / {total}', current=snapshot.current, total=snapshot.total))
    metrics[2].metric(tr('Temps écoulé'), _format_duration(snapshot.elapsed_seconds))
    metrics[3].metric(tr('Temps restant estimé'), _format_duration(snapshot.eta_seconds))
    st.progress(min(1.0, snapshot.current / total), text=tr(snapshot.message))
    st.markdown(tr('#### Console du moteur'))
    st.code(
        "\n".join(tr(line) for line in snapshot.log) if snapshot.log else tr("En attente des premiers messages…"),
        language=None,
        wrap_lines=True,
        height=220,
    )
    if snapshot.state == "running":
        st.caption(tr('Le temps restant est recalculé après chaque candidat terminé.'))
    if snapshot.state == "running":
        st.button(tr('Annuler le calcul'), on_click=_cancel_design_run)
        return
    if snapshot.state == "failed":
        st.error(tr(snapshot.error or "Le calcul a échoué."))
    if snapshot.state == "done":
        status = "warning" if snapshot.result.status == "cancelled" else "success"
        getattr(st, status)(tr(snapshot.message))
    if job.consume_once():
        if snapshot.state == "done":
            st.session_state.last_result = snapshot.result
            st.session_state.last_exports = snapshot.exports or {}
        else:
            st.session_state.run_error = snapshot.error or "Le calcul a échoué."
        st.rerun(scope="app")


def _candidate_status(candidate) -> str:
    if not getattr(candidate, "excluded", False):
        return tr("Évalué — STOP prématuré" if candidate.has_stop_codon else "Évalué")
    reason = getattr(candidate, "exclusion_reason", "") or "critère de criblage"
    reason = reason.partition(" — ")[0]
    if reason.startswith("Codon STOP prématuré"):
        reason = "STOP prématuré"
    return tr("Exclu — {reason}", reason=tr(reason))


def _result_rows(result) -> list[dict]:
    def display(value) -> str:
        return "—" if value is None else str(value)

    rows = []
    for rank, candidate in enumerate(result.candidates, start=1):
        excluded = getattr(candidate, "excluded", False)
        rows.append(
            {
                tr('Rang'): "—" if excluded else str(rank),
                tr('Statut'): _candidate_status(candidate),
                tr('Score'): display(candidate.score),
                tr('Defect'): display(candidate.defect),
                tr('ON %'): display(candidate.on_yield_pct),
                tr('Leak %'): display(candidate.leak_pct),
                tr('ΔG OFF'): display(candidate.delta_g_off),
                tr('ΔG ON'): display(candidate.delta_g_on),
                tr('ΔΔG'): display(candidate.ddg_activation),
                tr('ΔG RBS–linker'): display(candidate.delta_g_rbs_linker),
                tr('STOP'): tr("oui" if candidate.has_stop_codon else "non"),
                tr('Premier STOP'): display(candidate.first_stop) if candidate.has_stop_codon else "—",
                tr('Essai'): candidate.trial,
            }
        )
    return rows


def _render_results_page() -> None:
    st.markdown(tr('## Résultats'))
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
    st.markdown(
        tr('<div class="apta-subtitle">Classement, structures calculées et exports.</div>'),
        unsafe_allow_html=True,
    )
    render_saved_run_details("design")
    _design_job_monitor()
    result = st.session_state.get("last_result")
    if result is None:
        if not st.session_state.get("design_job"):
            st.info(tr('Ouvrez un run enregistré ci-dessus ou lancez un calcul depuis Conception.'))
        return

    with st.container(border=True):
        summary_cols = st.columns(4)
        summary_cols[0].metric(tr('Run'), result.run_id)
        summary_cols[1].metric(tr('Candidats'), len(result.candidates))
        summary_cols[2].metric(tr('Moteur'), result.request.engine.upper())
        summary_cols[3].metric(tr('Statut'), tr("Annulé" if result.status == "cancelled" else "Terminé"))
        for warning in result.warnings:
            st.warning(tr(warning))
        excluded_count = sum(getattr(candidate, "excluded", False) for candidate in result.candidates)
        if excluded_count:
            st.info(tr('{excluded_count} candidat(s) exclu(s) par le criblage, affiché(s) en fin de tableau avec leur motif. — = non calculé ; les exclus ne reçoivent ni rang ni score.', excluded_count=excluded_count))
        if result.candidates and excluded_count == len(result.candidates):
            st.warning(tr('Tous les candidats sont exclus par le criblage. Aucun candidat n’a été retenu pour le classement.'))

        st.caption(tr('Tableau de classement en lecture seule.'))
        st.table(_result_rows(result))

        exports = st.session_state.get("last_exports", {})
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
            st.caption(tr('Copie locale : {output_dir}', output_dir=result.output_dir))

    if not result.candidates:
        st.warning(tr("Le calcul n'a produit aucun candidat."))
        return

    candidates_by_trial = {candidate.trial: candidate for candidate in result.candidates}
    trial_options = list(candidates_by_trial)
    if st.session_state.get("selected_trial") not in trial_options:
        st.session_state.selected_trial = trial_options[0]
    selected_trial = st.selectbox(
        tr('Candidat à explorer'),
        options=trial_options,
        format_func={trial: (
            tr("{status} · essai {trial}", status=_candidate_status(candidates_by_trial[trial]), trial=trial)
            if getattr(candidates_by_trial[trial], "excluded", False) else
            tr("Rang {rank} · essai {trial} · score {score:.2f}",
               rank=trial_options.index(trial) + 1, trial=trial, score=candidates_by_trial[trial].score)
        ) for trial in trial_options}.__getitem__,
        key="selected_trial",
    )
    candidate = candidates_by_trial[selected_trial]
    _render_switch_candidate(result, candidate)
    _render_candidate_navigation(trial_options, selected_trial)


def _render_switch_candidate(result, candidate) -> None:
    excluded = getattr(candidate, "excluded", False)
    if excluded and not candidate.structure_on:
        with st.container(border=True):
            st.error(tr('Candidat exclu · {exclusion_reason}', exclusion_reason=tr(candidate.exclusion_reason)))
            st.caption(tr('Les calculs de structure, rendement et fuite ont été ignorés pour ce candidat. Sa séquence reste consultable et exportée.'))
            st.markdown(tr('**Trigger 5′ → 3′**'))
            st.code(candidate.trigger_seq, language=None)
            st.markdown(tr('**Switch 5′ → 3′**'))
            st.code(candidate.switch_seq, language=None)
        return

    with st.container(border=True):
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
                        switch_material=result.request.thermo.switch_material,
                        trigger_material=result.request.thermo.trigger_material,
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
                        result.request.thermo.switch_material, result.request.thermo.trigger_material,
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


def _open_nupack_directory() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.nupack_folder_error = ""
    try:
        folder = nupack_drop_directory(create=True)
        st.session_state.settings_nupack_path = str(folder)
        if sys.platform == "darwin":
            subprocess.run(["open", str(folder)], check=True, capture_output=True,
                           env=external_process_environment())
        elif sys.platform == "win32":
            os.startfile(str(folder))
        else:
            subprocess.run(["xdg-open", str(folder)], check=True, capture_output=True,
                           env=external_process_environment())
    except (OSError, subprocess.SubprocessError) as exc:
        st.session_state.nupack_folder_error = tr("Impossible d’ouvrir le dossier NUPACK : {error}", error=exc)


def _render_dependency_settings() -> None:
    st.subheader(tr('Moteurs externes'))
    st.caption(tr("Ces outils ne sont pas distribués avec Apta2Switch-Studio ; leurs licences propres s'appliquent."))
    settings = load_settings()
    nupack_settings = settings.get("nupack", {})
    if "settings_nupack_path" not in st.session_state:
        st.session_state.settings_nupack_path = nupack_settings.get("path") or str(nupack_drop_directory())
    st.caption(tr("Téléchargez NUPACK, décompressez son archive et placez le dossier obtenu dans le dossier ci-dessous. Puis confirmez votre licence et cliquez sur Valider NUPACK. Aucune commande ni installation de Python n’est nécessaire."))
    st.caption(tr("Dossier NUPACK : {path}", path=nupack_drop_directory()))
    download, folder = st.columns(2)
    download.link_button(tr("Télécharger NUPACK"), "https://www.nupack.org/downloads", width="stretch")
    folder.button(tr("Ouvrir le dossier NUPACK"), on_click=_open_nupack_directory,
                  key="open_nupack_directory", width="stretch")
    if st.session_state.get("nupack_folder_error"):
        st.error(tr(st.session_state.nupack_folder_error))
    with st.form("nupack_settings_form"):
        st.markdown(tr('#### NUPACK · design et thermodynamique'))
        st.text_input(
            tr('Dossier NUPACK ou fichier .whl'),
            key="settings_nupack_path",
            placeholder=tr('/chemin/vers/nupack'),
        )
        license_confirmed = st.checkbox(
            tr("Je confirme être responsable de la licence et de l'autorisation d'utilisation de NUPACK."),
            value=bool(nupack_settings.get("license_confirmed")),
        )
        validate_nupack = st.form_submit_button(tr('Valider NUPACK'), type="primary")
    if validate_nupack:
        if not st.session_state.settings_nupack_path.strip():
            st.error(tr("Indiquez d'abord un chemin NUPACK."))
        elif not license_confirmed:
            st.error(tr('Confirmez la responsabilité de licence avant de continuer.'))
        else:
            with st.spinner(tr("Validation de l'import NUPACK…")):
                validation = validate_nupack_path(st.session_state.settings_nupack_path)
            if validation.ok:
                save_nupack_settings(validation, True)
                st.success(tr(validation.message))
            else:
                st.error(tr(validation.message))
    nupack_path = configured_nupack_import_path()
    if nupack_path:
        st.success(tr('NUPACK activé depuis {nupack_path}', nupack_path=nupack_path))
    else:
        st.info(tr('Configurez NUPACK pour lancer les calculs de structure, d’extension et de conception.'))



def _render_reporter_settings() -> None:
    st.subheader(tr('Bibliothèque de reporters'))
    if st.session_state.reporter_return_to_design:
        st.button(tr('← Retour à la conception'), on_click=_return_to_design, key="return_from_reporters")
    st.caption(tr('La séquence commence juste après le codon de départ AUG du switch.'))
    reporters = all_reporters()
    for reporter in reporters.values():
        cols = st.columns([0.2, 0.62, 0.18], vertical_alignment="center")
        cols[0].markdown(tr('**{label}**', label=reporter.label))
        cols[1].code(reporter.sequence_after_start_rna, language=None)
        if is_builtin_reporter(reporter.key):
            cols[2].caption(tr('Intégré'))
        elif cols[2].button(tr('Supprimer'), key=f"remove_reporter_{reporter.key}"):
            remove_custom_reporter(reporter.key)
            if st.session_state.reporter_key == reporter.key:
                st.session_state.reporter_key = "sfgfp"
            st.rerun()

    with st.form("add_reporter_form", clear_on_submit=True):
        st.markdown(tr('#### Ajouter un reporter'))
        label = st.text_input(tr('Nom'), placeholder=tr('ex. mCherry'), key="new_reporter_label")
        sequence = st.text_area(tr('Séquence ARN/ADN après AUG'), height=90, key="new_reporter_sequence")
        submitted = st.form_submit_button(tr('Ajouter'))
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
    st.caption(tr('Par défaut, les runs sont enregistrés dans le dossier de l’application, sous runs/designs ou runs/extensions.'))
    if "settings_default_output" not in st.session_state:
        st.session_state.settings_default_output = _default_output_dir()
    with st.form("default_output_form"):
        st.text_input(tr("Dossier d'export par défaut"), key="settings_default_output")
        save_default = st.form_submit_button(tr('Enregistrer'))
    if save_default:
        path = st.session_state.settings_default_output.strip()
        if not path:
            st.error(tr('Le dossier ne peut pas être vide.'))
        else:
            settings = load_settings()
            settings["default_output_dir"] = path
            save_settings(settings)
            st.session_state.output_dir = path
            st.success(tr('Dossier par défaut enregistré.'))

    st.divider()
    st.markdown(tr('#### Réinitialisation'))
    st.caption(tr('Efface la configuration des moteurs, les reporters personnalisés et le dossier par défaut. Les exports ne sont pas supprimés.'))
    confirm = st.checkbox(tr('Je confirme la réinitialisation de tous les réglages'), key="confirm_reset_settings")
    st.button(
        tr('Réinitialiser les réglages'),
        disabled=not confirm,
        on_click=_reset_all_settings,
    )


def _save_language() -> None:
    set_language(st.session_state.ui_language)
    settings = load_settings()
    settings["ui_language"] = st.session_state.ui_language
    save_settings(settings)


def _remember_settings_tab() -> None:
    for label in ("Dépendances", "Reporters", "Général"):
        if st.session_state.settings_tab_label in (label, tr(label)):
            st.session_state.settings_tab = label
            break


def _render_settings_page() -> None:
    st.markdown(tr('## Réglages'))
    st.markdown(
        tr('<div class="apta-subtitle">Dépendances scientifiques, reporters et préférences enregistrées sur cette machine.</div>'),
        unsafe_allow_html=True,
    )
    st.selectbox(tr("Langue"), options=list(LANGUAGES), format_func=LANGUAGES.get,
                 key="ui_language", on_change=_save_language)
    st.session_state.settings_tab_label = tr(st.session_state.settings_tab)
    dependencies, reporters, general = st.tabs(
        [tr(label) for label in ("Dépendances", "Reporters", "Général")],
        key="settings_tab_label", on_change=_remember_settings_tab,
        default=tr(st.session_state.settings_tab),
    )
    with dependencies:
        _render_dependency_settings()
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
    _auto_detect_nupack_once()
    _inject_styles()
    _render_header()
    page = _render_sidebar()
    if page == "Conception":
        _render_design_page()
    elif page == "Résultats":
        _render_results_page()
    else:
        _render_settings_page()


if __name__ == "__main__":
    main()
