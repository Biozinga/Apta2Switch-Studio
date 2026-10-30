"""Read saved scientific runs without recomputing or changing their ranking."""

from __future__ import annotations

from aptaswitch_core.localization import tr

import json
import math
import os
from dataclasses import asdict, dataclass, fields, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

from aptaswitch_core.models import (
    DesignCandidate, ReporterContext, RunRequest, RunResult, ScoringWeights,
    SequenceInput, SwitchArchitecture, ThermoConditions, InverseExtensionSettings,
    FunctionalDefectWeights,
)
from aptaswitch_core.scoring import selection_score_100


@dataclass(frozen=True)
class SavedRun:
    path: Path
    kind: str
    run_id: str
    started_at: str
    molecule_name: str
    candidate_count: int


@dataclass(frozen=True)
class LoadedRun:
    kind: str
    result: RunResult | dict[str, Any]
    exports: dict[str, str]
    source: Path


def _object(value: Any, label: str) -> dict:
    if not isinstance(value, dict):
        raise ValueError(tr('{v0} : un objet JSON est attendu.', v0=label))
    return value


def _array(value: Any, label: str) -> list:
    if not isinstance(value, list):
        raise ValueError(tr('{v0} : une liste est attendue.', v0=label))
    return value


def _text(value: Any, label: str, *, empty: bool = False) -> str:
    if not isinstance(value, str) or (not empty and not value.strip()):
        raise ValueError(tr('{v0} : texte absent ou invalide.', v0=label))
    return value


def _number(value: Any, label: str, *, optional: bool = False) -> None:
    if value is None and optional:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(tr('{v0} : valeur numérique absente ou invalide.', v0=label))


def _finite_nested_numbers(value: Any, label: str) -> None:
    """Keep extensible diagnostic metadata, but never accept NaN/Infinity."""
    if isinstance(value, dict):
        for key, item in value.items():
            _finite_nested_numbers(item, f"{label}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _finite_nested_numbers(item, f"{label}[{index}]")
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        _number(value, label)


def _validate_selection_score(candidate: DesignCandidate) -> None:
    """Validate the bounded score against its recorded comparison sum."""
    _number(candidate.score_raw, "Score brut", optional=True)
    if candidate.score is not None:
        if not 0 <= candidate.score <= 100:
            raise ValueError("The recorded bounded selection score must be between 0 and 100")
        if candidate.score_raw is not None and not math.isclose(
            candidate.score, selection_score_100(candidate.score_raw), rel_tol=0, abs_tol=1e-4,
        ):
            raise ValueError("The bounded selection score does not match its recorded raw score")


def _validate_trigger_linearity(candidate: DesignCandidate) -> None:
    """Missing measurements stay missing; recorded probabilities remain auditable."""
    diagnostic = _object(candidate.trigger_linearity, "Linéarité du trigger seul")
    if not diagnostic:
        return
    _finite_nested_numbers(diagnostic, "Linéarité du trigger seul")
    for name in ("normalized_defect", "unpaired_fraction"):
        _probability(diagnostic.get(name), f"Linéarité du trigger · {name}")
    probabilities = _array(diagnostic.get("per_base_unpaired"), "Probabilités du trigger libre")
    length = len(candidate.trigger_seq)
    if len(probabilities) != length:
        raise ValueError("Trigger unpaired probabilities must match the trigger sequence length")
    for probability in probabilities:
        _probability(probability, "Probabilité du trigger libre")
    mean_unpaired = math.fsum(probabilities) / length
    if not math.isclose(diagnostic["unpaired_fraction"], mean_unpaired, abs_tol=1e-9):
        raise ValueError("The trigger unpaired fraction must match its per-base probabilities")
    if not math.isclose(diagnostic["normalized_defect"], 1 - mean_unpaired, abs_tol=1e-9):
        raise ValueError("The trigger linearity defect must be the expected paired fraction")
    structure = _text(diagnostic.get("structure_mfe"), "Structure du trigger seul")
    if len(structure) != length or set(structure) - set(".()"):
        raise ValueError("The isolated-trigger structure must match the trigger sequence")
    depth = 0
    for base in structure:
        depth += (base == "(") - (base == ")")
        if depth < 0:
            raise ValueError("The isolated-trigger structure must be balanced")
    if depth:
        raise ValueError("The isolated-trigger structure must be balanced")
    if "nucleotide_count" in diagnostic:
        _integer(diagnostic["nucleotide_count"], "Nombre de bases du trigger", minimum=1)
        if diagnostic["nucleotide_count"] != length:
            raise ValueError("The trigger linearity nucleotide count must match the sequence length")
    if "expected_paired_nt" in diagnostic:
        _number(diagnostic["expected_paired_nt"], "Bases appariées attendues du trigger")
        if not math.isclose(diagnostic["expected_paired_nt"], length * (1 - mean_unpaired), abs_tol=1e-8):
            raise ValueError("The expected paired trigger nucleotide count is inconsistent")
    if "model" in diagnostic:
        _object(diagnostic["model"], "Modèle du trigger seul")


def _validate_rbs_accessibility(candidate, request, thermodynamics):
    from aptaswitch_core.architecture import architecture_layout
    from aptaswitch_core.thermodynamics.accessibility import UnpairedWindowAccessibility
    from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, Strand
    diagnostic = _object(candidate.rbs_accessibility, "Ouverture RBS → AUG")
    if not diagnostic:
        return
    measured = UnpairedWindowAccessibility.from_dict(diagnostic)
    spec = ComplexSpec((Strand(candidate.trigger_seq, "dna", "Trigger"),
                        Strand(candidate.switch_seq, "rna", "Switch")))
    thermo = request.thermo
    conditions = Conditions(thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m)
    layout = architecture_layout(request.architecture)
    offset = len(candidate.trigger_seq)
    if (candidate.excluded or measured.complex_spec != spec or measured.conditions != conditions
            or measured.start != offset + layout["rbs_prefix_start"]
            or measured.end != offset + layout["aug_start"] + 3
            or any(thermodynamics.get(key) != value for key, value in asdict(measured.provenance).items())):
        raise ValueError("Recorded RBS accessibility differs from the candidate's ON context")


def _validate_defect_metrics(candidate: DesignCandidate, request: RunRequest) -> None:
    """Validate optional new measurements without inventing absent values."""
    for name in ("defect_structural", "defect_concentration"):
        value = getattr(candidate, name)
        if value is not None:
            _probability(value, f"Candidat · {name}")
    breakdown = _object(candidate.defect_breakdown, "Décomposition du défaut")
    regional = _object(candidate.regional_defects, "Défauts régionaux")
    for mapping, label in ((breakdown, "Décomposition du défaut"), (regional, "Défauts régionaux")):
        _finite_nested_numbers(mapping, label)
    for key in ("global", "structural", "concentration"):
        if key in breakdown:
            _probability(breakdown[key], f"Décomposition du défaut · {key}")
    if "tubes" in breakdown:
        for tube in _array(breakdown["tubes"], "Décomposition du défaut · tubes"):
            tube = _object(tube, "Défaut du tube")
            if "name" in tube:
                _text(tube["name"], "Nom du tube")
            for key in ("global", "structural", "concentration"):
                if key in tube:
                    _probability(tube[key], f"Défaut du tube · {key}")
            for key in ("dimensional_defect", "dimensional_structural", "dimensional_concentration", "target_nucleotide_concentration_m"):
                if key in tube:
                    _number(tube[key], f"Défaut du tube · {key}", optional=key == "target_nucleotide_concentration_m")
                    if tube[key] is not None and tube[key] < 0:
                        raise ValueError(f"Défaut du tube · {key} must not be negative")
            if "complexes" in tube:
                for complex_ in _array(tube["complexes"], "Défaut du tube · complexes"):
                    complex_ = _object(complex_, "Contribution du complexe")
                    if "name" in complex_:
                        _text(complex_["name"], "Nom du complexe")
                    for key in ("structural", "concentration", "total"):
                        if key in complex_:
                            _number(complex_[key], f"Contribution du complexe · {key}")
                            if complex_[key] < 0:
                                raise ValueError(f"Contribution du complexe · {key} must not be negative")
    if not regional:
        return
    if regional.get("scheme") != "functional_v2":
        raise ValueError("Functional defect measurements must identify the functional_v2 scheme")
    weights_type = FunctionalDefectWeights
    region_keys = tuple(field.name for field in fields(weights_type))
    saved_weights = _object(regional.get("weights", {}), "Poids régionaux enregistrés")
    if set(saved_weights) != set(region_keys):
        raise ValueError("Functional defect weights must name exactly the functional regions")
    region_weights = _model(weights_type, saved_weights, "Poids régionaux enregistrés")
    region_weights.validate()
    if region_weights != request.functional_defect_weights:
        raise ValueError("Functional defect weights do not match the recorded request")
    if "normalized_weights" in regional:
        normalized = _object(regional["normalized_weights"], "Parts régionales enregistrées")
        if set(normalized) != set(region_keys):
            raise ValueError("Functional weight fractions must name exactly the functional regions")
        expected = region_weights.fractions()
        for key in region_keys:
            _probability(normalized.get(key), f"Part régionale · {key}")
            if not math.isclose(normalized[key], expected[key], rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError("Regional weight fractions do not match the saved weights")
    for key in (*region_keys, "off", "on", "balanced"):
        if key in regional:
            if key in {"off", "on"} and regional[key] is None:
                if any(getattr(region_weights, zone) > 0 for zone in region_keys if zone.startswith(key + "_")):
                    raise ValueError("An enabled state's regional defect cannot be null")
                continue
            _probability(regional[key], f"Défaut régional · {key}")
    region_diagnostics = dict(_object(regional.get("regions", {}), "Régions du défaut"))
    if any(key not in regional for key in (*region_keys, "off", "on", "balanced")):
        raise ValueError("Functional defect measurements are incomplete")
    if "on_binding" in regional or set(region_diagnostics) != set(region_keys):
        raise ValueError("Functional defect measurements contain incompatible regions")
    if "off_stem_parts" in regional:
        parts = _object(regional["off_stem_parts"], "Parties de la tige OFF")
        if set(parts) != {"before_aug", "after_aug"}:
            raise ValueError("OFF stem diagnostics must contain the before- and after-AUG regions")
        region_diagnostics.update({f"off_stem_{key}": value for key, value in parts.items()})
    if region_diagnostics:
        for key, value in region_diagnostics.items():
            region = _object(value, f"Région · {key}")
            if "indices_switch" in region:
                indices = _array(region["indices_switch"], f"Région · {key} · indices")
                for index in indices:
                    _integer(index, f"Région · {key} · indice")
                    if index >= len(candidate.switch_seq):
                        raise ValueError("A regional defect index exceeds the switch length")
                if len(set(indices)) != len(indices):
                    raise ValueError("Regional defect indices must be unique")
            if "nucleotide_count" in region:
                _integer(region["nucleotide_count"], f"Région · {key} · nombre de bases")
                if region["nucleotide_count"] > len(candidate.switch_seq):
                    raise ValueError("Regional nucleotide count exceeds the switch length")
                if "indices_switch" in region and region["nucleotide_count"] != len(region["indices_switch"]):
                    raise ValueError("Regional defect nucleotide count does not match its indices")
            if "pair_count" in region:
                _integer(region["pair_count"], f"Région · {key} · nombre de paires")
                if "nucleotide_count" in region and 2 * region["pair_count"] != region["nucleotide_count"]:
                    raise ValueError("Stem pair count must describe both arms")
            if "expected_incorrect_nt" in region:
                _number(region["expected_incorrect_nt"], f"Région · {key} · bases incorrectes attendues")
                if not 0 <= region["expected_incorrect_nt"] <= region.get("nucleotide_count", len(candidate.switch_seq)) + 1e-9:
                    raise ValueError("Expected incorrect nucleotide count is outside the regional length")
            for name in ("normalized_defect", "share_of_switch_defect"):
                if name in region and (name != "share_of_switch_defect" or region[name] is not None):
                    _probability(region[name], f"Région · {key} · {name}")
    for key in ("switch_structural", "switch_expected_incorrect_nt", "per_base"):
        if key not in regional:
            continue
        for state, value in _object(regional[key], f"Défaut régional · {key}").items():
            if key == "per_base":
                values = _array(value, f"Défaut régional · {key} · {state}")
                if len(values) != len(candidate.switch_seq):
                    raise ValueError("Per-base defect count does not match the switch length")
                for probability in values:
                    _probability(probability, f"Défaut régional · {key} · {state}")
            elif key == "switch_structural":
                _probability(value, f"Défaut régional · {key} · {state}")
            else:
                _number(value, f"Défaut régional · {key} · {state}")
                if not 0 <= value <= len(candidate.switch_seq) + 1e-9:
                    raise ValueError("Expected incorrect switch nucleotide count is outside its length")


def _integer(value: Any, label: str, *, minimum: int = 0) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(tr('{v0} : entier attendu (minimum {v1}).', v0=label, v1=minimum))


def _sequence(value: Any, label: str, *, empty: bool = False) -> str:
    sequence = _text(value, label, empty=empty)
    if set(sequence.upper()) - set("ACGTU"):
        raise ValueError(tr('{v0} : séquence de bases invalide.', v0=label))
    return sequence


def _read_json(path: Path) -> dict:
    try:
        return _object(json.loads(path.read_text(encoding="utf-8-sig")), tr('Fichier de run'))
    except (OSError, UnicodeError) as exc:
        raise ValueError(tr('Impossible de lire le run « {v0} » : {v1}', v0=path.name, v1=exc)) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(tr('Le fichier « {v0} » n’est pas un JSON valide (ligne {v1}).', v0=path.name, v1=exc.lineno)) from exc
    except RecursionError as exc:
        raise ValueError(tr('Le fichier « {v0} » contient une imbrication JSON invalide.', v0=path.name)) from exc


def _resolve_results(path: Path) -> tuple[Path, dict]:
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        checkpoint = path / "design_checkpoint.json"
        if checkpoint.is_file() and not checkpoint.is_symlink():
            try:
                data = _read_json(checkpoint)
                optimization = _object(data.get("optimization", {}), "Design checkpoint")
                phase = _object(optimization.get("checkpoint", {}), "Design checkpoint").get("phase")
                if phase in {"design", "inverse_extension", "exports"} or data.get("status") in {"running", "paused"}:
                    return checkpoint, data
            except ValueError:
                # A damaged checkpoint must not hide an intact result archive.
                pass
        manifest = path / ".aptaswitch.json"
        if manifest.is_file():
            return _resolve_results(manifest)
        matches = sorted(path.glob("*_results.json"))
        if not matches and checkpoint.is_file() and not checkpoint.is_symlink():
            return checkpoint, _read_json(checkpoint)
        if not matches:
            raise ValueError(tr('Ce dossier ne contient pas de run. Sélectionnez le dossier d’un run ou son fichier de résultats JSON.'))
        if len(matches) > 1:
            raise ValueError(tr('Ce dossier contient plusieurs runs. Sélectionnez le fichier de résultats JSON du run à ouvrir.'))
        path = matches[0].resolve()
    payload = _read_json(path)
    if "results_file" in payload:
        filename = _text(payload["results_file"], tr('Fichier de résultats du projet'))
        candidate = path.parent / filename
        if Path(filename).name != filename or candidate.resolve().parent != path.parent or candidate.resolve() == path:
            raise ValueError(tr('Le manifeste doit désigner un fichier de résultats dans le même dossier.'))
        path = candidate.resolve()
        payload = _read_json(path)
        if "results_file" in payload:
            raise ValueError(tr('Le manifeste ne désigne pas un fichier de résultats.'))
    return path, payload


def _model(model_type, value: Any, label: str, **defaults):
    data = _object(value, label)
    known = {field.name for field in fields(model_type)}
    try:
        return model_type(**{**defaults, **{key: value for key, value in data.items() if key in known}})
    except (TypeError, ValueError) as exc:
        raise ValueError(tr('{v0} : données manquantes ou invalides.', v0=label)) from exc


def _warnings(payload: dict) -> list[str]:
    warnings = _array(payload.get("warnings", []), tr('Avertissements'))
    for warning in warnings:
        _text(warning, tr('Avertissement'), empty=True)
    return warnings


def _recorded_thermo(value: Any, label: str, *, system: bool = False) -> ThermoConditions:
    """Read recorded conditions without imposing today's calculation model.

    Archives keep the model and conditions that produced their measurements.
    New calculations still pass through ThermoConditions.validate_model().
    """
    thermo = _model(ThermoConditions, value, label, material="")
    _text(thermo.material, f"{label} · material")
    if system:
        thermo.validate_system_chemistry()
    else:
        thermo.validate_switch_chemistry()
    for name in ("temperature_c", "sodium_m", "magnesium_m", "concentration_m",
                 "trigger_concentration_m", "switch_concentration_m"):
        _number(getattr(thermo, name), f"{label} · {name}", optional=name.endswith("_concentration_m"))
    return thermo


def _optimization_fractions(value: Any, count: int, label: str) -> list[float]:
    values = _array(value, label)
    if len(values) != count:
        raise ValueError(f"{label}: incorrect number of defect fractions")
    for number in values:
        _number(number, label)
        if not 0 <= number <= 1:
            raise ValueError(f"{label}: defects must be fractions from zero to one")
    return values


def _optimization_objective(value: Any, weights: FunctionalDefectWeights, *, rbs_weight=None, translation_search=None) -> tuple:
    objective = _object(value, "Optimization objective")
    local = _optimization_fractions(objective.get("local_defects"), 4, "Local defects")
    guards = _optimization_fractions(objective.get("guard_defects"), 3, "Guard defects")
    global_defect, = _optimization_fractions([objective.get("global_defect")], 1, "Global defect")
    _optimization_rbs_penalty(objective)
    if rbs_weight is not None:
        _number(objective.get("rbs_linker_energy"), "Optimization RBS energy")
        if objective.get("rbs_weight") != rbs_weight:
            raise ValueError("Recorded RBS objective weight differs from the design request")
    if translation_search is not None:
        for name in ("rbs_calculator_weight", "trigger_linearity_weight"):
            if objective.get(name) != translation_search[name]:
                raise ValueError("Recorded translation objective weights differ from the request")
    recorded_weights = _array(objective.get("local_weights"), "Local objective weights")
    if len(recorded_weights) != 4:
        raise ValueError("The optimization objective requires four regional weights")
    for weight in recorded_weights:
        _number(weight, "Local objective weight")
        if weight < 0:
            raise ValueError("Negative local objective weight")
    scale = max(recorded_weights)
    if scale <= 0:
        raise ValueError("The optimization objective requires an active region")
    scaled = [weight / scale for weight in recorded_weights]
    fractions = [weight / math.fsum(scaled) for weight in scaled]
    if any(not math.isclose(actual, expected, abs_tol=1e-12)
           for actual, expected in zip(fractions, weights.fractions().values())):
        raise ValueError("Recorded objective weights differ from the design request")
    for field in ("site_weights", "global_site_weights"):
        for entry in _array(objective.get(field), "Mutation site weights"):
            if not isinstance(entry, list) or len(entry) != 2:
                raise ValueError("Invalid recorded mutation site weight")
            site = _object(entry[0], "Mutation site")
            _text(site.get("strand"), "Mutation strand")
            _integer(site.get("index"), "Mutation site index")
            _number(entry[1], "Mutation site weight")
            if entry[1] < 0:
                raise ValueError("Negative mutation site weight")
    return math.fsum(weight * defect for weight, defect in zip(fractions, local)), global_defect, local, guards


def _optimization_rbs_penalty(objective):
    energy = objective.get("rbs_linker_energy", 0.)
    weight = objective.get("rbs_weight", 0.)
    _number(energy, "Optimization RBS energy")
    _number(weight, "Optimization RBS weight")
    penalty = weight * (max(0., -energy) / 100.)
    if weight < 0 or not math.isfinite(penalty):
        raise ValueError("Invalid RBS stability optimization penalty")
    for weight_name, value_name, invert in (("rbs_calculator_weight", "rbs_calculator_merit", True),
                                             ("trigger_linearity_weight", "trigger_linearity_defect", False)):
        extra_weight, value = objective.get(weight_name, 0.), objective.get(value_name, 0.)
        _number(extra_weight, weight_name)
        _number(value, value_name)
        if extra_weight < 0 or not 0 <= value <= 1:
            raise ValueError("Invalid translation refinement measurement")
        penalty += extra_weight * (1. - value if invert else value)
    if not math.isfinite(penalty):
        raise ValueError("Invalid translation refinement penalty")
    return penalty


def _fixed_reading_outcomes(request, payload, optimization, *, retained_candidate):
    """Validate actual fixed evaluation metadata without an optimizer trace."""
    evaluations = optimization.get("evaluations")
    rows = _array(payload.get("results"), "Fixed reading candidates")
    if (optimization.get("algorithm") != "fixed_sequence_evaluation"
            or optimization.get("objective") != "supplied_sequences"
            or optimization.get("sequence_optimization") is not False
            or optimization.get("trials") != []
            or type(evaluations) is not int or evaluations not in (0, 1)
            or len(rows) != evaluations):
        raise ValueError("The recorded fixed-sequence evaluation is inconsistent")
    if payload.get("status", "completed") == "completed" and evaluations != 1:
        raise ValueError("Completed fixed reading is missing its evaluation")
    scope = optimization.get("scope", "full_run")
    if scope not in ("full_run", "retained_candidate"):
        raise ValueError("Unknown recorded fixed reading scope")
    if scope == "retained_candidate" and (
        not retained_candidate or evaluations != 1
        or optimization.get("retained_trial") != 1
        or optimization.get("source_trial_count") != 1
        or optimization.get("source_run_id") != payload.get("run_id")
    ):
        raise ValueError("Retained fixed reading identity differs from the selected candidate")
    return {1} if evaluations else set()


def _load_design(payload: dict, source: Path, *, retained_candidate: bool = False) -> RunResult:
    data = _object(payload.get("request"), tr('Paramètres du design'))
    if data.get("engine") != "native-hybrid":
        raise ValueError(tr('Ce run ne provient pas d’un moteur scientifique pris en charge. Les anciens runs de démonstration ne peuvent pas être ouverts comme résultats scientifiques.'))
    sequences = _model(SequenceInput, data.get("sequences"), tr('Séquences'))
    _sequence(sequences.aptamer_dna, tr('Aptamère'))
    _sequence(sequences.trigger_dna, "Trigger")
    _text(sequences.molecule_name, tr('Nom du système'), empty=True)
    _integer(sequences.binding_start, tr('Début de la région reconnue'))
    if sequences.binding_length is not None:
        _integer(sequences.binding_length, tr('Longueur de la région reconnue'), minimum=1)
    # Historical loops included their entire sequence in rbs_loop. Adding the
    # current three-base prefix would silently redraw those saved structures.
    architecture = _model(SwitchArchitecture, data.get("architecture"), "Architecture", rbs_prefix_length=0)
    if architecture.binding_layout not in ("legacy_lower_stem", "full_left_arm"):
        raise ValueError(tr('La disposition de reconnaissance du trigger est inconnue.'))
    for name in ("toehold_length", "stem_length", "upper_stem_length", "upper_stem2_length",
                 "bulge_length", "aug_spacer", "rbs_prefix_length"):
        _integer(getattr(architecture, name), f"Architecture · {name}")
    for name in ("rbs_loop", "t7_leader"):
        _sequence(getattr(architecture, name), f"Architecture · {name}", empty=True)
    if architecture.frame_linker is not None:
        _text(architecture.frame_linker, "Linker", empty=True)
        if set(architecture.frame_linker.upper()) - set("ACGTUX"):
            raise ValueError(tr('Linker : séquence invalide.'))
    thermo = _recorded_thermo(data.get("thermo", {}), "Conditions")
    reporter = _model(ReporterContext, data.get("reporter", {}), tr('Rapporteur'))
    reporter.validate()
    _sequence(reporter.sequence_after_start_rna, tr('Séquence du rapporteur'), empty=True)
    _text(reporter.label, tr('Rapporteur'), empty=True)
    scoring = _model(ScoringWeights, data.get("scoring", {}), tr('Pondérations'), rbs_calculator=0.)
    for field in fields(scoring):
        _number(getattr(scoring, field.name), tr('Pondération · {v0}', v0=field.name))
    if scoring.rbs_calculator < 0:
        raise ValueError("The RBS Calculator weight must not be negative")
    if scoring.trigger_linearity < 0:
        raise ValueError("The trigger linearity weight must not be negative")
    functional_weights_data = _object(data.get("functional_defect_weights", {}), "Poids des défauts fonctionnels")
    if set(functional_weights_data) != {field.name for field in fields(FunctionalDefectWeights)}:
        raise ValueError("Functional defect weights must name exactly the functional regions")
    functional_weights = _model(FunctionalDefectWeights, functional_weights_data, "Poids des défauts fonctionnels")
    functional_weights.validate()
    request = _model(RunRequest, data, tr('Paramètres du design'))
    _text(data.get("defect_metric"), "Defect scoring metric")
    if request.defect_metric != "functional":
        raise ValueError("The recorded defect scoring metric is unknown")
    _text(data.get("score_scale"), "Échelle du score")
    if request.score_scale != "bounded_v1":
        raise ValueError("The recorded selection score scale is unknown")
    inverse = None
    if data.get("inverse_extension") is not None:
        inverse_data = _object(data["inverse_extension"], "Extension inverse")
        inverse = _model(InverseExtensionSettings, inverse_data, "Extension inverse")
        inverse = InverseExtensionSettings(inverse.target_length, inverse.extension_side, inverse.top_candidates,
                                          _recorded_thermo(inverse_data.get("thermo", {}), "Conditions extension inverse", system=True))
        inverse.extension_lengths(len(sequences.trigger_dna))
    request = RunRequest(
        sequences=sequences, architecture=architecture, thermo=thermo,
        reporter=reporter, scoring=scoring, trials=request.trials,
        output_dir=source.parent.parent, engine=request.engine,
        exclude_stop_candidates=request.exclude_stop_candidates,
        exclude_non_linear_candidates=request.exclude_non_linear_candidates,
        inverse_extension=inverse,
        defect_metric=request.defect_metric,
        functional_defect_weights=functional_weights,
        score_scale=request.score_scale,
        thermodynamics_library_path=request.thermodynamics_library_path,
        optimization_options=_object(data.get("optimization_options", {}), "Options d’optimisation"),
        fixed_switch_sequence=request.fixed_switch_sequence,
    )
    if request.fixed_switch_sequence is not None:
        if request.exclude_stop_candidates or request.exclude_non_linear_candidates:
            raise ValueError("Fixed reading archives cannot exclude supplied sequences before measuring them")
        request = request.normalized()
    _integer(request.trials, tr('Nombre d’essais'), minimum=1)
    from aptaswitch_core.thermodynamics.optimization import TrialOptions
    from aptaswitch_core.thermodynamics.types import Conditions, ModelProvenance
    try:
        TrialOptions(**request.optimization_options)
    except TypeError as exc:
        raise ValueError("Unknown recorded optimization options") from exc
    thermodynamics = _object(payload.get("thermodynamics", {}), "Thermodynamic provenance")
    optimization = _object(payload.get("optimization", {}), "Optimization outcomes")
    _finite_nested_numbers(thermodynamics, "Thermodynamic provenance")
    _finite_nested_numbers(optimization, "Optimization outcomes")
    if thermodynamics:
        model_fields = {field.name for field in fields(ModelProvenance)}
        try:
            recorded_model = ModelProvenance(**{key: thermodynamics[key] for key in model_fields})
        except (KeyError, TypeError) as exc:
            raise ValueError("Incomplete recorded thermodynamic provenance") from exc
        if recorded_model.model != request.thermo.material:
            raise ValueError("Recorded model provenance differs from the design request")
        if thermodynamics.get("engine", request.engine) != request.engine:
            raise ValueError("Recorded engine differs from the design request")
    if request.engine == "native-hybrid":
        if not thermodynamics or not optimization:
            raise ValueError("Native design archives require thermodynamic provenance and optimization outcomes")
        if thermodynamics.get("engine") != request.engine:
            raise ValueError("Recorded engine differs from the design request")
        expected_conditions = {"celsius": thermo.temperature_c, "sodium": thermo.sodium_m,
                               "magnesium": thermo.magnesium_m}
        try:
            recorded_conditions = Conditions(**_object(thermodynamics.get("conditions"), "Thermodynamic conditions"))
        except TypeError as exc:
            raise ValueError("Invalid recorded thermodynamic conditions") from exc
        if asdict(recorded_conditions) != expected_conditions:
            raise ValueError("Recorded thermodynamic conditions differ from the design request")
        if request.fixed_switch_sequence is not None:
            report_trials = _fixed_reading_outcomes(request, payload, optimization, retained_candidate=retained_candidate)
        else:
            if (optimization.get("algorithm") != "diverse_local_first"
                    or optimization.get("objective") != "functional_regions_then_tube_defect"
                    or optimization.get("seed_rule") != "seed + trial - 1"):
                raise ValueError("Unknown recorded native optimization method")
            translation_search = optimization.get("translation_search")
            if translation_search is not None:
                from aptaswitch_core.thermodynamics.optimization import translation_search_settings
                if translation_search != translation_search_settings(request.scoring, request.inverse_extension is not None):
                    raise ValueError("Recorded translation search method differs from the design request")
            expected_rbs_weight = None
            if "rbs_stability" in optimization:
                from aptaswitch_core.thermodynamics.optimization import rbs_stability_settings
                if optimization["rbs_stability"] != rbs_stability_settings(request.scoring.rbs):
                    raise ValueError("Recorded RBS refinement method differs from the design request")
                expected_rbs_weight = request.scoring.rbs
            try:
                recorded_weights = FunctionalDefectWeights(**_object(optimization.get("functional_weights"), "Optimization weights"))
                recorded_weights.validate()
            except TypeError as exc:
                raise ValueError("Invalid recorded optimization weights") from exc
            if recorded_weights != request.functional_defect_weights:
                raise ValueError("Recorded optimization weights differ from the design request")
            try:
                recorded_options = TrialOptions(**_object(optimization.get("options"), "Optimization options"))
            except TypeError as exc:
                raise ValueError("Unknown recorded optimization options") from exc
            if asdict(recorded_options) != asdict(TrialOptions(**request.optimization_options)):
                raise ValueError("Recorded optimization options differ from the design request")
            reports = _array(optimization.get("trials"), "Optimization trials")
            scope = optimization.get("scope", "full_run")
            if scope not in ("full_run", "retained_candidate"):
                raise ValueError("Unknown recorded optimization scope")
            selected_trial = None
            if scope == "retained_candidate":
                if not retained_candidate:
                    raise ValueError("A retained candidate is not a complete design run")
                selected_trial = optimization.get("retained_trial")
                _integer(selected_trial, "Retained optimization trial", minimum=1)
                source_trial_count = optimization.get("source_trial_count")
                _integer(source_trial_count, "Source design trial count", minimum=1)
                rows = _array(payload.get("results"), "Retained candidate")
                if (source_trial_count != request.trials or selected_trial > source_trial_count
                        or optimization.get("source_run_id") != payload.get("run_id")
                        or len(rows) != 1 or _object(rows[0], "Retained candidate").get("trial") != selected_trial
                        or len(reports) != 1):
                    raise ValueError("Retained optimization identity differs from the selected candidate")
            report_trials = set()
            for report in reports:
                report = _object(report, "Optimization trial")
                number = report.get("trial")
                _integer(number, "Optimization trial", minimum=1)
                if number in report_trials or number > request.trials:
                    raise ValueError("Recorded optimization trial identities are inconsistent")
                report_trials.add(number)
                if selected_trial is not None and number != selected_trial:
                    raise ValueError("Retained candidate has a different optimization trial")
                try:
                    trial_options = TrialOptions(**_object(report.get("options"), "Trial options"))
                except TypeError as exc:
                    raise ValueError("Unknown recorded trial options") from exc
                if trial_options != replace(recorded_options, seed=recorded_options.seed + number - 1):
                    raise ValueError("Recorded trial options differ from the optimization settings")
                stop_reason = _text(report.get("stop_reason"), "Optimization stopping condition")
                if stop_reason not in {"target_reached", "budget_exhausted", "space_exhausted", "no_mutable_sites", "no_admissible_sequence", "cancelled"}:
                    raise ValueError("Unknown recorded optimization stopping condition")
                for field, maximum in (("evaluations", trial_options.max_evaluations),
                                       ("proposals", trial_options.max_proposals),
                                       ("population_evaluated", trial_options.population_size),
                                       ("lineages_started", trial_options.lineages)):
                    _integer(report.get(field), "Optimization " + field)
                    if report[field] > maximum:
                        raise ValueError("Recorded optimization usage exceeds its budget")
                _integer(report.get("cache_hits"), "Optimization cache hits")
                if report.get("phase") not in {"initialization", "local", "global"}:
                    raise ValueError("Unknown recorded optimization phase")
                if report["population_evaluated"] > report["evaluations"]:
                    raise ValueError("Recorded population exceeds the evaluation count")
                _number(report.get("best_loss"), "Optimization best loss", optional=True)
                loss = report.get("best_loss")
                if (loss is not None and not 0 <= loss <= 1) or (loss is None) != (report["evaluations"] == 0):
                    raise ValueError("Recorded optimization loss is inconsistent with its evaluations")
                global_loss = report.get("best_global_loss")
                _number(global_loss, "Optimization global loss", optional=True)
                if (global_loss is None) != (loss is None) or (global_loss is not None and not 0 <= global_loss <= 1):
                    raise ValueError("Recorded global defect is inconsistent with the local defect")
                objective = report.get("best_objective")
                if (objective is None) != (loss is None):
                    raise ValueError("Recorded optimization objective is missing")
                if objective is not None:
                    measured_local, measured_global, local_defects, guards = _optimization_objective(objective, functional_weights, rbs_weight=expected_rbs_weight, translation_search=translation_search)
                    if not math.isclose(loss, measured_local, abs_tol=1e-12) or not math.isclose(global_loss, measured_global, abs_tol=1e-12):
                        raise ValueError("Recorded optimization losses differ from the measured objective")
                limits = {}
                for field, count in (("local_limits", 4), ("guard_limits", 3)):
                    limits[field] = None if report.get(field) is None else _optimization_fractions(report[field], count, field)
                if limits["local_limits"] is not None and limits["guard_limits"] is None:
                    raise ValueError("Local optimization limits require fixed concentration and duplex guards")
                if report["phase"] == "global":
                    if any(value is None for value in limits.values()):
                        raise ValueError("The global phase requires fixed local and guard limits")
                    if objective is not None and (any(value > limit + 1e-12 for value, limit in zip(local_defects, limits["local_limits"]))
                                                  or any(value > limit + 1e-12 for value, limit in zip(guards, limits["guard_limits"]))):
                        raise ValueError("The retained objective violates the fixed optimization limits")
                if report["stop_reason"] == "target_reached" and (loss is None or loss > trial_options.target_local_loss
                        or global_loss + _optimization_rbs_penalty(objective) > trial_options.target_global_loss):
                    raise ValueError("Recorded optimization target was not reached")
                elites = _array(report.get("elites"), "Optimization elites")
                if len(elites) > trial_options.lineages:
                    raise ValueError("Too many recorded optimization elites")
                for elite in elites:
                    _, _, elite_local, elite_guards = _optimization_objective(_object(elite, "Optimization elite").get("objective"), functional_weights, rbs_weight=expected_rbs_weight, translation_search=translation_search)
                    if limits["local_limits"] is not None and (any(value > limit + 1e-12 for value, limit in zip(elite_local, limits["local_limits"]))
                                                               or any(value > limit + 1e-12 for value, limit in zip(elite_guards, limits["guard_limits"]))):
                        raise ValueError("A recorded elite violates the fixed optimization limits")
                _integer(report.get("candidate_screened"), "Screened optimization candidates")
                selected_objective = report.get("selected_objective")
                if selected_objective is not None:
                    _, _, selected_local, selected_guards = _optimization_objective(selected_objective, functional_weights, rbs_weight=expected_rbs_weight, translation_search=translation_search)
                    if report["candidate_screened"] < 1:
                        raise ValueError("The selected candidate has no screening record")
                    if limits["local_limits"] is not None and (any(value > limit + 1e-12 for value, limit in zip(selected_local, limits["local_limits"]))
                                                               or any(value > limit + 1e-12 for value, limit in zip(selected_guards, limits["guard_limits"]))):
                        raise ValueError("The selected candidate violates the fixed optimization limits")
                budget = report.get("budget_kind")
                if budget is not None:
                    _text(budget, "Optimization budget kind")
                if report["stop_reason"] == "budget_exhausted":
                    if (budget not in {"evaluations", "proposals"}
                            or report[budget] != getattr(trial_options, "max_" + budget)):
                        raise ValueError("Recorded exhausted optimization budget is inconsistent")
                elif budget is not None:
                    raise ValueError("Recorded optimization budget does not match its stopping condition")
            if (scope == "full_run" and payload.get("status", "completed") == "completed"
                    and len(report_trials) != request.trials):
                raise ValueError("Completed native design is missing optimization trials")
    candidates = []
    trials = set()
    for row in _array(payload.get("results"), tr('Candidats du design')):
        candidate = _model(DesignCandidate, row, tr('Candidat du design'))
        _integer(candidate.trial, tr('Numéro d’essai'), minimum=1)
        if candidate.trial in trials:
            raise ValueError(tr('Plusieurs candidats portent le même numéro d’essai.'))
        trials.add(candidate.trial)
        if request.engine == "native-hybrid" and candidate.trial not in report_trials:
            raise ValueError("Native candidate has no matching optimization outcome")
        _sequence(candidate.switch_seq, tr('Séquence du switch'))
        _sequence(candidate.trigger_seq, tr('Séquence du trigger candidat'))
        if request.fixed_switch_sequence is not None and (
            candidate.switch_seq != request.fixed_switch_sequence
            or candidate.trigger_seq != request.sequences.trigger_dna
            or candidate.excluded or candidate.trial != 1
        ):
            raise ValueError("The fixed reading candidate changed the supplied sequences")
        for name in ("excluded", "has_stop_codon"):
            if not isinstance(getattr(candidate, name), bool):
                raise ValueError(tr('Candidat · {v0} : booléen attendu.', v0=name))
        for name in ("score", "defect", "on_yield_pct", "leak_pct", "delta_g_off",
                     "delta_g_on", "ddg_activation", "delta_g_rbs_linker"):
            _number(getattr(candidate, name), tr('Candidat · {v0}', v0=name), optional=candidate.excluded)
        for name in ("structure_off", "structure_on", "reporter_label", "exclusion_reason"):
            _text(getattr(candidate, name), tr('Candidat · {v0}', v0=name), empty=True)
        candidates.append(candidate)
        _object(candidate.inverse_extension_metrics, "Comparaison extension inverse")
        _validate_defect_metrics(candidate, request)
        _validate_selection_score(candidate)
        _validate_trigger_linearity(candidate)
        _validate_rbs_accessibility(candidate, request, thermodynamics)
        prediction = _object(candidate.rbs_calculator, "RBS Calculator")
        if prediction:
            from aptaswitch_core.architecture import architecture_layout
            from aptaswitch_core.rbs_calculator import validate_rbs_prediction
            layout = architecture_layout(request.architecture)
            validate_rbs_prediction(prediction, candidate.switch_seq, layout["aug_start"],
                tuple(range(layout["len_t7"], layout["len_t7"] + layout["len_trigger"])))
        elif scoring.rbs_calculator > 0 and not candidate.excluded:
            raise ValueError("A weighted RBS Calculator measurement is missing")
    return RunResult(
        run_id=_run_id(payload, source), request=request, candidates=candidates,
        output_dir=source.parent, started_at=_text(payload.get("started_at", ""), tr('Date de début'), empty=True),
        completed_at=_text(payload.get("completed_at", ""), tr('Date de fin'), empty=True),
        status=_text(payload.get("status", "completed"), tr('Statut')), warnings=_warnings(payload),
        inverse_extension_analysis=_object(payload.get("inverse_extension_analysis", {}), "Analyse extension inverse"),
        thermodynamics=thermodynamics,
        optimization=optimization,
    )


def _structure(data: dict, lengths: tuple[int, int], label: str) -> None:
    structure = _text(data.get("mfe_structure"), f"{label} · structure")
    if set(structure) - set(".()+") or list(map(len, structure.split("+"))) != list(lengths):
        raise ValueError(tr('{v0} : la structure ne correspond pas aux deux séquences.', v0=label))
    depth = 0
    for base in structure:
        depth += (base == "(") - (base == ")")
        if depth < 0:
            raise ValueError(tr('{v0} : structure non équilibrée.', v0=label))
    if depth:
        raise ValueError(tr('{v0} : structure non équilibrée.', v0=label))
    _number(data.get("mfe_energy_kcal_mol"), tr('{v0} · énergie MFE', v0=label))
    probabilities = data.get("nucleotide_probabilities")
    if probabilities is not None:
        if len(_array(probabilities, tr('{v0} · probabilités', v0=label))) != sum(lengths):
            raise ValueError(tr('{v0} : il faut une probabilité par base.', v0=label))
        for probability in probabilities:
            _probability(probability, label)


def _probability(value: Any, label: str) -> None:
    _number(value, tr('{v0} · probabilité', v0=label))
    if not -1e-12 <= value <= 1 + 1e-12:
        raise ValueError(tr('{v0} : probabilité en dehors de [0, 1].', v0=label))


def _load_extension(payload: dict, source: Path) -> dict:
    model = _object(payload.get("model"), tr('Modèle de l’extension'))
    if model.get("aptamer_material") != "dna" or model.get("trigger_material") != "dna":
        raise ValueError(tr("L’aptamère et le trigger doivent être en ADN."))
    inputs = _object(payload.get("input"), tr('Entrées de l’extension'))
    aptamer = _sequence(inputs.get("aptamer_dna"), tr('Aptamère'))
    original = _sequence(inputs.get("trigger_dna"), tr('Trigger original'))
    if inputs.get("extension_side", "5prime") not in ("5prime", "3prime", "both"):
        raise ValueError(tr('L’extrémité choisie pour l’extension est inconnue.'))
    for position in _array(inputs.get("ligand_aptamer_positions", []), tr('Bases liées au ligand')):
        _integer(position, tr('Position liée au ligand'), minimum=1)
        if position > len(aptamer):
            raise ValueError(tr('Une base liée au ligand dépasse la longueur de l’aptamère.'))
    for key in ("request", "model", "mfe_screen", "curves"):
        if key in payload:
            _object(payload[key], key)
    for key in ("total_candidates_screened", "ensemble_candidates_analyzed"):
        stats = payload.get("mfe_screen", {})
        if key in stats:
            _integer(stats[key], tr('Criblage · {v0}', v0=key))
    if "request" in payload and "molecule_name" in payload["request"]:
        _text(payload["request"]["molecule_name"], tr('Nom du système'), empty=True)
    _warnings(payload)
    reference = payload.get("reference")
    if reference is not None:
        _structure(_object(reference, tr('Référence')), (len(aptamer), len(original)), tr('Référence'))
    candidates = _array(payload.get("candidates"), tr('Candidats étendus'))
    if candidates and not reference:
        raise ValueError(tr('La structure de référence du run d’extension est absente.'))
    for row in candidates:
        candidate = _object(row, tr('Candidat étendu'))
        sequence = _sequence(candidate.get("extended_trigger"), tr('Trigger étendu'))
        prefix = _sequence(candidate.get("extension_5prime"), "Extension 5′", empty=True)
        suffix = _sequence(candidate.get("extension_3prime", ""), "Extension 3′", empty=True)
        if (prefix + original + suffix).upper() != sequence.upper():
            raise ValueError(tr('Le trigger étendu ne conserve pas le trigger original entre ses extensions.'))
        for key in ("rank_within_length", "trigger_length", "extension_length"):
            _integer(candidate.get(key), tr('Candidat · {v0}', v0=key), minimum=1)
        if candidate["trigger_length"] != len(sequence) or candidate["extension_length"] != len(prefix) + len(suffix):
            raise ValueError(tr('Les longueurs enregistrées ne correspondent pas au trigger étendu.'))
        for key in ("final_score", "core_pair_probability_rmse", "reference_pair_probability_loss",
                    "extension_pair_probability_sum", "extension_to_aptamer_probability_sum",
                    "missing_reference_mfe_pairs", "new_core_mfe_pairs"):
            _number(candidate.get(key), tr('Candidat · {v0}', v0=key))
        _structure(candidate, (len(aptamer), len(sequence)), tr('Candidat étendu'))
        for pair in _array(candidate.get("reference_pairs", []), tr('Paires de référence')):
            pair = _object(pair, tr('Paire de référence'))
            for key in ("reference_probability", "candidate_probability"):
                _probability(pair.get(key), tr('Paire de référence'))
    # Only location metadata changes when a directory has been moved. All
    # measured values, recorded conditions, candidate order and scores survive.
    payload["run_id"] = _run_id(payload, source)
    payload["output_dir"] = str(source.parent)
    if "request" in payload and "output_dir" in payload["request"]:
        payload["request"]["output_dir"] = str(source.parent.parent)
    for key in ("started_at", "completed_at"):
        _text(payload.get(key, ""), key, empty=True)
    _text(payload.get("status", "completed"), tr('Statut'))
    return payload


def _run_id(payload: dict, source: Path) -> str:
    run_id = _text(payload.get("run_id", source.stem.removesuffix("_results")), tr('Identifiant du run'))
    if run_id in (".", "..") or "/" in run_id or "\\" in run_id:
        raise ValueError(tr('L’identifiant du run contient un séparateur de dossier.'))
    return run_id


def _existing_exports(source: Path, run_id: str, kind: str) -> dict[str, str]:
    exports = {"json": str(source)}
    names = {
        "csv": f"{run_id}_ranked.csv", "xlsx": f"{run_id}_ranked.xlsx",
        "log": f"{run_id}.log", "project": ".aptaswitch.json",
        "zip": f"{run_id}_complete.zip",
    }
    if kind == "design":
        names["svg"] = f"{run_id}_top1.svg"
        names["inverse_csv"] = f"{run_id}_inverse_extension.csv"
    else:
        names["tube_csv"] = f"{run_id}_tube_complexes.csv"
        for path in sorted(source.parent.iterdir()):
            if path.name.startswith(("svg_", "png_")) and path.suffix.lower() in (".svg", ".png"):
                key = path.stem if path.suffix.lower() == ".svg" else path.stem.replace("svg_", "png_", 1)
                names[key] = path.name
    for key, name in names.items():
        path = source.parent / name
        if path.is_file() and path.resolve().parent == source.parent:
            exports[key] = str(path)
    return exports


def load_run(path: Path) -> LoadedRun:
    """Open a run directory, results JSON, or sibling project manifest.

    Files are read only. Reopening does not load a scientific engine, rerank candidates,
    regenerate exports, or mutate stored parameters and results.
    """
    try:
        source, payload = _resolve_results(path)
    except (OSError, RuntimeError) as exc:
        raise ValueError(tr('Impossible d’ouvrir ce chemin de run : {v0}', v0=exc)) from exc
    if "results" in payload and "request" in payload:
        result = _load_design(payload, source)
        kind, run_id = "design", result.run_id
    elif "input" in payload and "candidates" in payload:
        result = _load_extension(payload, source)
        kind, run_id = "extension", result["run_id"]
    else:
        raise ValueError(tr('Ce JSON n’est pas un résultat de design de switches ou d’extension de triggers. Ouvrez le fichier *_results.json du run.'))
    return LoadedRun(kind, result, _existing_exports(source, run_id, kind), source)


def discover_runs(roots: Sequence[Path]) -> list[SavedRun]:
    """Find valid saved runs under output roots, categories or run folders."""
    sources: set[Path] = set()
    for root in roots:
        root = Path(root).expanduser()
        if root.is_file():
            sources.add(root)
            continue
        for directory, folders, filenames in os.walk(root, followlinks=False):
            folders[:] = [name for name in folders if not name.startswith(".") and not (Path(directory) / name).is_symlink()]
            sources.update(Path(directory) / name for name in filenames
                           if (name.endswith("_results.json") or name in {".aptaswitch.json", "design_checkpoint.json"})
                           and not (Path(directory) / name).is_symlink())
    saved = {}
    seen: set[Path] = set()
    for path in sorted(sources):
        try:
            loaded = load_run(path)
            if loaded.source in seen:
                continue
            seen.add(loaded.source)
            result = loaded.result
            if isinstance(result, RunResult):
                row = SavedRun(loaded.source, loaded.kind, result.run_id, result.started_at,
                               result.request.sequences.molecule_name, len(result.candidates))
            else:
                row = SavedRun(loaded.source, loaded.kind, result["run_id"], result.get("started_at", ""),
                               result.get("request", {}).get("molecule_name", "Système personnalisé"), len(result["candidates"]))
            # One run can have a checkpoint and final exports. A newer pending
            # checkpoint takes precedence over exports from an earlier pause;
            # once finalized, use the ordinary results archive again.
            is_checkpoint = loaded.source.name == "design_checkpoint.json"
            pending = isinstance(result, RunResult) and (
                result.optimization.get("checkpoint", {}).get("phase") in {"design", "inverse_extension", "exports"}
                or result.status in {"running", "paused"})
            priority = 2 if is_checkpoint and pending else 0 if is_checkpoint else 1
            identity = (loaded.kind, loaded.source.parent, row.run_id)
            if identity not in saved or priority > saved[identity][0]:
                saved[identity] = (priority, row)
        except (ValueError, OSError, RecursionError):
            # One damaged/incomplete file must not hide the remaining runs.
            continue

    def date_key(row: SavedRun) -> tuple[float, str]:
        try:
            stamp = datetime.fromisoformat(row.started_at.replace("Z", "+00:00")).timestamp()
        except (ValueError, OverflowError, OSError):
            try:
                stamp = row.path.stat().st_mtime
            except OSError:
                stamp = 0.0
        return stamp, row.run_id

    return sorted((row for _, row in saved.values()), key=date_key, reverse=True)
