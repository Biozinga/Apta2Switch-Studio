"""Read saved scientific runs without recomputing or changing their ranking."""

from __future__ import annotations

from aptaswitch_core.localization import tr

import json
import math
import os
from dataclasses import dataclass, fields
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

from aptaswitch_core.models import (
    DesignCandidate, ReporterContext, RunRequest, RunResult, ScoringWeights,
    SequenceInput, SwitchArchitecture, ThermoConditions,
)


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
        manifest = path / ".aptaswitch.json"
        if manifest.is_file():
            return _resolve_results(manifest)
        matches = sorted(path.glob("*_results.json"))
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


def _load_design(payload: dict, source: Path) -> RunResult:
    data = _object(payload.get("request"), tr('Paramètres du design'))
    if data.get("engine") != "nupack":
        raise ValueError(tr('Ce run n’est pas un calcul NUPACK. Les anciens runs de démonstration ne peuvent pas être ouverts comme résultats scientifiques.'))
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
    for name in ("toehold_length", "stem_length", "upper_stem_length", "upper_stem2_length",
                 "bulge_length", "aug_spacer", "rbs_prefix_length"):
        _integer(getattr(architecture, name), f"Architecture · {name}")
    for name in ("rbs_loop", "t7_leader"):
        _sequence(getattr(architecture, name), f"Architecture · {name}", empty=True)
    if architecture.frame_linker is not None:
        _text(architecture.frame_linker, "Linker", empty=True)
        if set(architecture.frame_linker.upper()) - set("ACGTUN"):
            raise ValueError(tr('Linker : séquence invalide.'))
    thermo = _model(ThermoConditions, data.get("thermo", {}), "Conditions")
    for name in ("temperature_c", "sodium_m", "magnesium_m", "concentration_m",
                 "trigger_concentration_m", "switch_concentration_m"):
        _number(getattr(thermo, name), f"Conditions · {name}", optional=name.endswith("_concentration_m"))
    reporter = _model(ReporterContext, data.get("reporter", {}), tr('Rapporteur'))
    _sequence(reporter.sequence_after_start_rna, tr('Séquence du rapporteur'), empty=True)
    _text(reporter.label, tr('Rapporteur'), empty=True)
    scoring = _model(ScoringWeights, data.get("scoring", {}), tr('Pondérations'))
    for field in fields(scoring):
        _number(getattr(scoring, field.name), tr('Pondération · {v0}', v0=field.name))
    # Historical runs did not apply filters introduced after they were saved.
    # Preserve recorded policies without applying today's defaults or reranking.
    request = _model(
        RunRequest, data, tr('Paramètres du design'),
        exclude_stop_candidates=False, exclude_non_linear_candidates=False,
    )
    request = RunRequest(
        sequences=sequences, architecture=architecture, thermo=thermo,
        reporter=reporter, scoring=scoring, trials=request.trials,
        output_dir=source.parent.parent, engine="nupack", nupack_path=request.nupack_path,
        exclude_stop_candidates=request.exclude_stop_candidates,
        exclude_non_linear_candidates=request.exclude_non_linear_candidates,
    )
    _integer(request.trials, tr('Nombre d’essais'), minimum=1)
    candidates = []
    trials = set()
    for row in _array(payload.get("results"), tr('Candidats du design')):
        candidate = _model(DesignCandidate, row, tr('Candidat du design'))
        _integer(candidate.trial, tr('Numéro d’essai'), minimum=1)
        if candidate.trial in trials:
            raise ValueError(tr('Plusieurs candidats portent le même numéro d’essai.'))
        trials.add(candidate.trial)
        _sequence(candidate.switch_seq, tr('Séquence du switch'))
        _sequence(candidate.trigger_seq, tr('Séquence du trigger candidat'))
        for name in ("excluded", "has_stop_codon"):
            if not isinstance(getattr(candidate, name), bool):
                raise ValueError(tr('Candidat · {v0} : booléen attendu.', v0=name))
        for name in ("score", "defect", "on_yield_pct", "leak_pct", "delta_g_off",
                     "delta_g_on", "ddg_activation", "delta_g_rbs_linker"):
            _number(getattr(candidate, name), tr('Candidat · {v0}', v0=name), optional=candidate.excluded)
        for name in ("structure_off", "structure_on", "reporter_label", "exclusion_reason"):
            _text(getattr(candidate, name), tr('Candidat · {v0}', v0=name), empty=True)
        candidates.append(candidate)
    return RunResult(
        run_id=_run_id(payload, source), request=request, candidates=candidates,
        output_dir=source.parent, started_at=_text(payload.get("started_at", ""), tr('Date de début'), empty=True),
        completed_at=_text(payload.get("completed_at", ""), tr('Date de fin'), empty=True),
        status=_text(payload.get("status", "completed"), tr('Statut')), warnings=_warnings(payload),
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
    else:
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

    Files are read only. Reopening does not import NUPACK, rerank candidates,
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
                           if (name.endswith("_results.json") or name == ".aptaswitch.json")
                           and not (Path(directory) / name).is_symlink())
    saved = []
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
            saved.append(row)
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

    return sorted(saved, key=date_key, reverse=True)
