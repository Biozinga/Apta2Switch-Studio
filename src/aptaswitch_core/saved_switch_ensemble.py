"""Adapt a saved full-reporter reference to the shared OFF/ON ensemble analysis.

The selected analysis owns its conditions and sequences; current form values
and the design's short reporter prefix are never substituted. Sampling remains
conditional on the switch monomer or the formed Trigger+Switch complex. Adding
an aptamer to the equilibrium tubes does not add a strand to those complexes.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from types import SimpleNamespace

from .defect_metrics import _pair_targets
from .localization import tr
from .models import ReporterContext, SwitchArchitecture, ThermoConditions
from .saved_switch_analysis import build_saved_switch_context, validate_analysis_config


@dataclass(frozen=True)
class ContextEnsembleRequest:
    architecture: SwitchArchitecture
    thermo: ThermoConditions
    reporter: ReporterContext
    engine: str = "native-hybrid"
    thermodynamics_library_path: str | None = None


@dataclass(frozen=True)
class ContextEnsembleSource:
    request: ContextEnsembleRequest
    output_dir: Path
    thermodynamics: dict | None = None


@dataclass(frozen=True)
class ContextEnsembleCandidate:
    trial: int
    switch_seq: str
    trigger_seq: str
    reporter_label: str


def context_ensemble_source(
    snapshot: dict, analysis: dict, output_dir: Path, *, backend=None,
) -> tuple[ContextEnsembleSource, ContextEnsembleCandidate]:
    """Validate a saved reference and detach inputs for the shared sampler.

    A completed baseline remains usable if a later temperature/salt sweep was
    cancelled or failed. The reference's conditions must agree with its saved
    configuration, and its strand identities/order must match the immutable
    candidate including the entire reporter. No filesystem or native call is
    performed here. ``output_dir`` is supplied by the trusted caller, never read
    from the imported analysis.
    """
    if (not isinstance(analysis, dict) or not isinstance(analysis.get("config"), dict)
            or not isinstance(analysis.get("baseline"), dict)):
        raise ValueError(tr("Le calcul du contexte complet ne contient pas de référence OFF/ON exploitable."))
    # Reopening a saved artifact requires its recorded model, not an installed
    # runtime. An explicit backend, used for new sampling, must match exactly.
    recorded = analysis.get("thermodynamics")
    if backend is None and recorded is not None:
        from .thermodynamics.types import ModelProvenance
        if not isinstance(recorded, dict) or not isinstance(analysis.get("model_limits"), dict):
            raise ValueError("The saved native context is missing its model identity or domain")
        backend = SimpleNamespace(provenance=ModelProvenance(**recorded),
                                  condition_limits=analysis["model_limits"])
    config = validate_analysis_config(snapshot, analysis["config"], backend=backend)
    expected = build_saved_switch_context(snapshot, config, backend=backend)
    context = analysis.get("context")
    baseline = analysis["baseline"]
    if backend is not None and (analysis.get("thermodynamics") != asdict(backend.provenance)
                               or baseline.get("thermodynamics") != asdict(backend.provenance)):
        raise ValueError("Les sous-familles doivent utiliser le modèle exact de la référence enregistrée.")
    conditions = {key: config[key] for key in ("temperature_c", "sodium_m", "magnesium_m")}
    # Acquisition provenance and display names do not change the molecular
    # system: an earlier analysis may have obtained the same complete reporter
    # from user input instead of the saved reporter library. Validate every
    # scientific field while leaving that original provenance untouched.
    descriptive_fields = {"candidate_name", "reporter_source"}
    if (not isinstance(context, dict)
            or any(context.get(key) != value for key, value in expected.items()
                   if key not in descriptive_fields)
            or baseline.get("conditions") != conditions):
        raise ValueError(tr("Les données du contexte complet ne correspondent pas au candidat et aux conditions enregistrés."))

    for state, names in (("off", ["Switch"]), ("on", ["Trigger", "Switch"])):
        sequences = ([context["trigger_sequence"]] if state == "on" else []) + [context["switch_sequence"]]
        try:
            mfe = baseline[state]["mfe"]
            if mfe["strand_order"] != names or mfe["sequences"] != sequences:
                raise ValueError("Different strand identities or order")
            _pair_targets(mfe["structure"], list(map(len, sequences)))
        except (KeyError, TypeError, AttributeError, ValueError) as error:
            raise ValueError(tr("La structure enregistrée ne correspond pas aux séquences du contexte complet.")) from error

    saved_request, saved_candidate = snapshot["request"], snapshot["candidate"]
    architecture = SwitchArchitecture(**saved_request["architecture"])
    architecture.validate_for_trigger(len(context["trigger_sequence"]))
    thermo = ThermoConditions(
        **conditions, material=context["material"],
        concentration_m=config["switch_concentration_m"],
        switch_concentration_m=config["switch_concentration_m"],
        trigger_concentration_m=config["trigger_concentration_m"],
        screening_max_complex_size=config["max_complex_size"],
        validation_max_complex_size=config["max_complex_size"],
    )
    if backend is None:
        if saved_request.get("engine") != "native-hybrid":
            raise ValueError("The archived calculation cannot be sampled with the supported engine.")
        thermo.validate_model()
    reporter = replace(
        ReporterContext(**saved_request["reporter"]),
        sequence_after_start_rna=context["reporter_sequence_after_start"].replace("T", "U"),
    )
    reporter.validate()
    trial = saved_candidate.get("trial")
    if type(trial) is not int or trial < 1:
        raise ValueError(tr("Les données du contexte complet ne correspondent pas au candidat et aux conditions enregistrés."))
    engine = saved_request["engine"] if backend is None else backend.provenance.backend
    if backend is not None:
        from .thermodynamics.setup import MODEL_ID
        if (backend.provenance.backend == "AptaSwitch-native-hybrid"
                and backend.provenance.model == MODEL_ID):
            engine = "native-hybrid"
    request = ContextEnsembleRequest(
        architecture=architecture, thermo=thermo, reporter=reporter,
        engine=engine,
    )
    candidate = ContextEnsembleCandidate(
        trial, context["switch_sequence"], context["trigger_sequence"],
        saved_candidate.get("reporter_label") or reporter.label,
    )
    return ContextEnsembleSource(request, Path(output_dir),
                                 asdict(backend.provenance) if backend is not None else None), candidate
