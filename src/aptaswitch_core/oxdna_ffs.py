"""Guided preparation and resumable native direct-FFS orchestration."""
from __future__ import annotations

from dataclasses import fields
import hashlib
import copy
import json
import math
import os
from pathlib import Path
import secrets
import subprocess
import time

from .architecture import architecture_layout, target_switch_structures
from .localization import tr
from .models import SwitchArchitecture, ThermoConditions
from .oxdna_box import estimate_ffs_box
from .oxdna_kinetics import (
    FFS_SOURCE_SAMPLES_DEFAULT, FFS_SOURCE_SAMPLES_MIN, FFS_SOURCE_SAMPLES_MAX,
    FFS_TRIALS_DEFAULT, FFS_TRIALS_MIN, FFS_TRIALS_MAX,
    OPENING_PROTOCOL_VERSION, PROXIMITY_PROTOCOL_VERSION,
    PROXIMITY_DISTANCE_MIN, PROXIMITY_DISTANCE_MAX, validate_kinetic_protocol,
)
from .oxdna_setup import PARAMETERS
from .preview import parse_structure
from .scientific_process import external_process_environment, process_group_options, stop_process


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity(snapshot: dict) -> str:
    # Same immutable identity as candidate-scoped 3D artifacts.
    data = {key: snapshot.get(key) for key in ("snapshot_id", "saved_at", "candidate", "request", "sequences")}
    return hashlib.sha256(json.dumps(data, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def _sequence(value: str, material: str) -> str:
    seq = "".join(str(value).upper().split())
    if not seq or set(seq) - set("ACGTU"):
        raise ValueError(tr("La simulation exige des séquences ADN/ARN entièrement définies."))
    return seq.replace("U", "T") if material == "dna" else seq.replace("T", "U")


def _guide(value: str, sequence: str) -> str:
    if not isinstance(value, str) or len(value) != len(sequence) or set(value) - set("()."):
        raise ValueError(tr("Le guide oxDNA doit contenir une position par nucléotide, en notation point-parenthèse."))
    depth = 0
    for symbol in value:
        depth += (symbol == "(") - (symbol == ")")
        if depth < 0:
            raise ValueError(tr("Le guide oxDNA contient des parenthèses non équilibrées."))
    if depth:
        raise ValueError(tr("Le guide oxDNA contient des parenthèses non équilibrées."))
    return value


def build_ffs_plan(snapshot: dict) -> dict:
    request, candidate = snapshot["request"], snapshot["candidate"]
    arch = SwitchArchitecture(**{key: value for key, value in request["architecture"].items()
                                 if key in {field.name for field in fields(SwitchArchitecture)}})
    layout = architecture_layout(arch)
    sw, trigger = candidate["switch_seq"], candidate["trigger_seq"]
    n, nt = len(sw), len(trigger)
    reporter_prefix = request["reporter"]["sequence_after_start_rna"]
    if (nt != layout["len_trigger"] or not reporter_prefix
            or n != layout["reporter_start"] + len(reporter_prefix)
            or sw[layout["reporter_start"]:].upper().replace("T", "U") != reporter_prefix.upper().replace("T", "U")):
        raise ValueError(tr("L’architecture du candidat ne correspond pas aux séquences à simuler."))
    off, on = target_switch_structures(arch, n - layout["reporter_start"])
    stem_pairs = [list(pair) for pair in parse_structure(off)[0]]
    binding = []
    for i, j in parse_structure(on)[0]:
        if i < nt <= j:
            binding.append([j - nt, n + i])
    bases = (sw + trigger).upper().replace("U", "T")
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    if not binding or any(complement.get(bases[i]) != bases[j] for i, j in binding):
        raise ValueError(tr("Les bases du trigger ne correspondent pas au duplex cible enregistré."))
    # Initiation region: RBS prefix through AUG, excluding the reporter.
    region = list(range(layout["rbs_prefix_start"], layout["aug_start"] + 3))
    opening = sorted({(min(i, j), max(i, j)) for i in region for j in range(n + nt)
                      if i != j and not (j < n and abs(i - j) == 1)})
    # Biological endpoint; the native protocol supplies its fixed nested interfaces.
    interfaces = [{"name": "opening", "min_binding": len(binding),
                   "max_stem": 0, "max_opening": 0}]
    toehold = [pair for pair in binding
               if layout["len_t7"] <= pair[0] < layout["len_t7"] + layout["len_toehold"]]
    if not toehold:
        raise ValueError(tr("La simulation nécessite une région toehold accessible non vide."))
    return {"binding_pairs": binding, "toehold_pairs": toehold,
            "stem_pairs": stem_pairs, "opening_pairs": [list(pair) for pair in opening],
            "switch_length": n, "nucleotide_count": n + nt, "distance_metric": "nucleotide_centres",
            # Both strands are indexed 5′→3′. The antiparallel entry pair is
            # the first base after the leader and the trigger's final base.
            "distance_scope": "trigger_3prime_to_toehold_start",
            "distance_pairs": [[layout["len_t7"], n + nt - 1]],
            "separation_pairs": [[i, j] for i in range(n) for j in range(n, n + nt)],
            "opening_region": region, "interfaces": interfaces, "basin_distance": 9.,
            "basin_min_stem": max(1, math.ceil(len(stem_pairs) * .8)),
            "hb_threshold": -.1, "indexing": "0-based switch then trigger, each 5-to-3",
            "target_source": "saved_design_architecture", "rate_estimation": True}


def _number(options: dict, name: str, default, lower, upper, *, integer=False):
    value = options.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not lower <= value <= upper:
        raise ValueError(tr("Paramètre FFS invalide : {name}.", name=name))
    if integer and int(value) != value:
        raise ValueError(tr("Paramètre FFS entier requis : {name}.", name=name))
    return int(value) if integer else float(value)


def ffs_visualization(snapshot: dict) -> dict:
    """Use the 2D domain palette, indexed switch then trigger as in oxDNA."""
    from .preview import build_switch_model, TRIGGER_COLOR, EXTENSION_COLOR, EXTENSION_3PRIME_COLOR
    request, candidate = snapshot["request"], snapshot["candidate"]
    reporter = request["reporter"]
    model = build_switch_model(SwitchArchitecture(**request["architecture"]),
                               candidate["trigger_seq"], reporter["sequence_after_start_rna"],
                               reporter["label"])
    if len(model.switch_colors) != len(candidate["switch_seq"]):
        raise ValueError("The oxDNA domain colors do not match the saved switch length")
    labels = [domain.label for domain in model.domains for _ in domain.sequence]
    colors = list(model.switch_colors)
    trigger_colors = [TRIGGER_COLOR] * len(candidate["trigger_seq"])
    trigger_labels = ["Trigger"] * len(trigger_colors)
    inverse = request.get("inverse_extension")
    if inverse:
        from .trigger_extension import extension_split
        original = request["sequences"]["trigger_dna"]
        prefix, suffix = extension_split(len(trigger_colors) - len(original), inverse["extension_side"])
        if candidate["trigger_seq"][prefix:len(trigger_colors)-suffix if suffix else None].replace("U", "T") != original.replace("U", "T"):
            raise ValueError("The trigger extension does not match its saved original sequence")
        for indices, color, label in ((range(prefix), EXTENSION_COLOR, "Extension 5′"),
                                      (range(len(trigger_colors)-suffix, len(trigger_colors)), EXTENSION_3PRIME_COLOR, "Extension 3′")):
            for index in indices:
                trigger_colors[index], trigger_labels[index] = color, label
    colors.extend(trigger_colors)
    labels.extend(trigger_labels)
    legend = list(dict.fromkeys(zip(labels, colors)))
    return {"sequence": list(candidate["switch_seq"] + candidate["trigger_seq"]),
            "colors": colors, "domains": labels, "legend": legend}


def prepare_ffs_request(snapshot: dict, runtime: dict, options: dict, *, guides: dict | None = None, backend=None) -> dict:
    """Freeze a candidate and options without folding on the UI thread.

    Guides are normally resolved at job execution. Runs have no time limit by default.
    Explicit guides are for validated internal callers and reproducible tests.
    """
    if not snapshot.get("snapshot_id"):
        raise ValueError(tr("La FFS est réservée aux toehold switches enregistrés."))
    thermo = snapshot["request"]["thermo"]
    conditions = ThermoConditions(**{field.name: thermo[field.name] for field in fields(ThermoConditions) if field.name in thermo})
    conditions.validate_model()
    conditions.validate_switch_chemistry()
    sw = _sequence(snapshot["candidate"]["switch_seq"], "rna")
    trigger = _sequence(snapshot["candidate"]["trigger_seq"], "dna")
    if len(trigger) != 30:
        raise ValueError(tr("La cinétique exige un switch ARN et un trigger ADN de 30 nucléotides."))
    allowed = {"kinetics_mode", "kinetics_source_samples", "kinetics_trials", "kinetics_replicas", "kinetics_workers",
               "seed", "temperature_c", "salt_m", "import_path", "thermodynamics_library_path"}
    unknown = set(options) - allowed
    if unknown:
        raise ValueError(tr("Paramètre cinétique inconnu : {name}.", name=sorted(unknown)[0]))
    if max(len(sw), len(trigger)) > 2000:
        raise ValueError(tr("Chaque brin oxDNA est limité à 2 000 nucléotides."))
    mode = options.get("kinetics_mode", "opening")
    if mode not in ("opening", "proximity"):
        raise ValueError(tr("Mode de simulation FFS inconnu."))
    kinetics = {"version": OPENING_PROTOCOL_VERSION, "mode": mode,
                "source_samples": (_number(options, "kinetics_source_samples", FFS_SOURCE_SAMPLES_DEFAULT,
                                          FFS_SOURCE_SAMPLES_MIN, FFS_SOURCE_SAMPLES_MAX, integer=True)
                                   if mode == "opening" else 0),
                "trials": _number(options, "kinetics_trials", FFS_TRIALS_DEFAULT,
                                  FFS_TRIALS_MIN, FFS_TRIALS_MAX, integer=True),
                "replicas": _number(options, "kinetics_replicas", 10, 4, 100, integer=True),
                "workers": _number(options, "kinetics_workers", min(10, os.cpu_count() or 1), 1, 64, integer=True),
                "equilibration_steps": 200000, "chunk_steps": 1000}
    if mode == "proximity":
        kinetics.update(version=PROXIMITY_PROTOCOL_VERSION, source_samples=0,
                        equilibration_method="isolated_monomers_before_close_placement_v1",
                        initial_distance_min=PROXIMITY_DISTANCE_MIN,
                        initial_distance_max=PROXIMITY_DISTANCE_MAX)
    kinetics["workers"] = min(kinetics["workers"], kinetics["replicas"])
    sampling = {"seed": _number(options, "seed", secrets.randbelow(2**31 - 1) + 1, 1, 2**31 - 1, integer=True),
                "dt": .003, "frame_interval": 1000, "max_movie_frames": 600,
                "temperature_c": _number(options, "temperature_c", thermo["temperature_c"], 0., 90.),
                "salt_m": _number(options, "salt_m", thermo["sodium_m"], .01, 2.)}
    from .thermodynamics.setup import get_backend, configured_library_path
    from dataclasses import asdict
    library_path = (options.get("thermodynamics_library_path")
                    or snapshot["request"].get("thermodynamics_library_path") or configured_library_path())
    backend = backend or get_backend(library_path)
    limits = backend.condition_limits
    if not limits["sodium_min_m"] <= sampling["salt_m"] <= limits["sodium_max_m"]:
        raise ValueError(tr("Le sel choisi dépasse les limites du modèle thermodynamique utilisé pour préparer la 2D."))
    if not limits["temperature_min_c"] <= sampling["temperature_c"] <= limits["temperature_max_c"]:
        raise ValueError(tr("La température choisie dépasse les limites du modèle thermodynamique utilisé pour préparer la 2D."))
    # The finite periodic volume remains identical in source flux and shooting.
    box_policy = estimate_ffs_box(len(sw), len(trigger), "rna", "dna",
                                  sampling["temperature_c"], sampling["salt_m"])
    box = box_policy["box_size"]
    box_policy["estimated_box_size"] = box
    result = {"format": "aptaswitch-oxdna-ffs", "version": 2, "snapshot_id": snapshot["snapshot_id"],
            "saved_switch_identity": _identity(snapshot), "runtime": dict(runtime),
            "preparation_mode": "native_hybrid_guided_oxna", "kinetics": kinetics,
            "switch_sequence": sw, "trigger_sequence": trigger, "switch_material": "rna", "trigger_material": "dna",
            "temperature_c": sampling["temperature_c"], "salt_m": sampling["salt_m"],
            "seed": sampling["seed"], "use_average_seq": False, "model": "NA", "box_size": box,
            "box_policy": box_policy,
            "initial_separation": 10., "sampling": sampling, "plan": build_ffs_plan(snapshot),
            "visualization": ffs_visualization(snapshot),
            "guide_model": backend.provenance.model,
            "guide_thermodynamics": asdict(backend.provenance),
            "thermodynamics_library_path": str(library_path or ""),
            "guide_conditions": {"temperature_c": sampling["temperature_c"], "sodium_m": sampling["salt_m"],
                                 "magnesium_m": 0., "ensemble": backend.provenance.ensemble},
            "concentration_each_m": box_policy["concentration_each_m"],
            "magnesium_modeled": False, "experimental_magnesium_m": thermo.get("magnesium_m", 0.),
            "purpose": "native_direct_ffs_flux_and_shooting"}
    if mode == "proximity":
        result["purpose"] = "native_ffs_conditional_opening_from_unbound_proximity"
        result["plan"].update(rate_estimation=False, distance_scope="trigger_to_toehold",
            distance_pairs=[[i, j] for i in sorted({pair[0] for pair in result["plan"]["toehold_pairs"]})
                            for j in range(len(sw), len(sw) + len(trigger))])
    if guides is not None:
        for role in ("switch", "trigger"):
            result[role + "_guide"] = _guide(guides[role], result[role + "_sequence"])
        result["guide_source"] = "provided_internal_guides"
    return result


def resolve_ffs_guides(request: dict, *, cancel_check=None, progress_callback=None, backend=None) -> dict:
    """Resolve new isolated-strand guides; never recompute saved run guides."""
    from .scientific_process import CalculationCancelled
    from .thermodynamics.setup import get_backend
    from .thermodynamics.application_analysis import analyze_strand
    from dataclasses import asdict
    result = copy.deepcopy(request)

    def check():
        if cancel_check and cancel_check():
            raise CalculationCancelled(tr("Préparation oxDNA annulée."))

    check()
    if all(role + "_guide" in result for role in ("switch", "trigger")):
        for role in ("switch", "trigger"):
            _guide(result[role + "_guide"], result[role + "_sequence"])
        return result
    backend = backend or get_backend(result.get("thermodynamics_library_path"))
    provenance = asdict(backend.provenance)
    if (result.get("guide_thermodynamics") != provenance
            or result.get("guide_model") != backend.provenance.model):
        raise ValueError(tr("Le modèle des guides 2D diffère du moteur sélectionné. Créez une nouvelle simulation ; les anciens checkpoints restent inchangés."))
    values = result["guide_conditions"]
    conditions = ThermoConditions(material=backend.provenance.model,
        temperature_c=values["temperature_c"], sodium_m=values["sodium_m"], magnesium_m=values["magnesium_m"])
    guides = {}
    for role in ("switch", "trigger"):
        check()
        if progress_callback:
            progress_callback(tr("Guide 2D du {role} isolé…", role=role))
        sequence = result[role + "_sequence"]
        record = analyze_strand(backend, sequence, result[role + "_material"], conditions)
        check()
        result[role + "_guide"] = _guide(record["structure"], sequence)
        guides[role] = {**record, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}
    result.update(guide_source="native_hybrid_mfe_isolated_strands", guide_provenance=guides)
    return result


def _write(path: Path, data: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    os.replace(temporary, path)


def load_ffs_resume_request(output_dir: Path) -> dict:
    """Reuse the frozen experiment and its prepared coordinates, never redesign it."""
    output = Path(output_dir).resolve()
    path = output / "worker_request.json"
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 32 * 1024 * 1024:
        raise ValueError(tr("Aucune préparation FFS sauvegardée ne permet de reprendre ce calcul."))
    request = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(request, dict) or not isinstance(request.get("prepared"), dict):
        raise ValueError(tr("La préparation FFS sauvegardée est invalide."))
    if not isinstance(request.get("kinetics"), dict):
        raise ValueError(tr("Ce dossier ne contient pas de simulation FFS native reprenable."))
    validate_kinetic_protocol(request["kinetics"])
    prepared = request["prepared"]
    # Candidate folders can be moved. Only rebase known preparation artifacts;
    # no file supplied by a saved run may escape its own calculation folder.
    old_preparation = Path(prepared["topology_file"]).parent
    for field in ("configuration_file", "topology_file"):
        value = Path(prepared[field])
        relative = value.relative_to(old_preparation) if value.is_absolute() else value
        if ".." in relative.parts:
            raise ValueError("Unsafe FFS preparation path")
        target = output / "preparation" / relative
        if not target.is_file() or target.is_symlink() or not target.resolve().is_relative_to(output):
            raise ValueError(tr("Un fichier de préparation nécessaire à la reprise est manquant."))
        prepared[field] = str(target)
    request["resume"] = True
    return request


def ffs_worker_running(output_dir: Path) -> bool:
    """The native worker owns this lock, independently of any browser session."""
    import fcntl
    path = Path(output_dir) / "worker.lock"
    if not path.is_file() or path.is_symlink():
        return False
    try:
        with path.open("r") as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return True
            fcntl.flock(handle, fcntl.LOCK_UN)
    except OSError:
        return False
    return False


def run_oxdna_ffs(request: dict, output_dir: Path, progress_callback=None, cancel_check=None) -> dict:
    """Run until sampling finishes or cancellation, preserving partial artifacts."""
    from .oxdna_guided import prepare_guided_system
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if ffs_worker_running(output):
        raise ValueError(tr("Cette simulation est déjà en cours."))
    resuming = bool(request.get("resume"))
    prior_elapsed = 0.
    previous = {}
    if resuming:
        request = load_ffs_resume_request(output)
        try:
            previous = json.loads((output / "result.json").read_text())
            prior_elapsed = float(previous.get("elapsed_seconds", 0.))
        except (OSError, ValueError, TypeError):
            pass
        (output / "CANCEL").unlink(missing_ok=True)
        (output / "resume_error.json").unlink(missing_ok=True)
    started = time.monotonic()
    process, result, prepared = None, {}, None
    total = request["kinetics"]["replicas"]

    def progress(message):
        if progress_callback:
            progress_callback(0, total, message)

    try:
        if resuming:
            prepared = request["prepared"]
            progress(tr("Reprise de la simulation depuis les configurations sauvegardées…"))
        else:
            request = resolve_ffs_guides(request, cancel_check=cancel_check, progress_callback=progress)
            _write(output / "resolved_request.json", request)
            prepared = prepare_guided_system(request, output / "preparation", progress_callback=progress,
                                            cancel_check=cancel_check, deadline=None)
        if not resuming and prepared.get("box_size") is not None:
            box = float(prepared["box_size"])
            if not math.isfinite(box) or box < request["box_size"]:
                raise ValueError("The prepared box is smaller than the automatic contour requirement.")
            request["box_size"] = box
            request["concentration_each_m"] = 1 / (6.02214076e23 * (box * .8518e-9)**3 * 1000)
            request["box_policy"].update(box_size=box, box_size_nm=box*.8518,
                concentration_each_m=request["concentration_each_m"],
                placement_adjusted=box > request["box_policy"]["estimated_box_size"])
            _write(output / "resolved_request.json", request)
        if cancel_check and cancel_check():
            raise RuntimeError(tr("Simulation annulée après la préparation."))
        if prepared.get("ffs_eligible") is False:
            if prepared.get("guide_accepted") is False:
                raise ValueError(tr("Après {count} essais de préparation, aucune conformation ne respecte exactement le guide 2D sans contraintes. La meilleure reste consultable avec ses écarts ; la FFS n’a pas été lancée.", count=len(prepared.get("attempts", []))))
            raise ValueError(tr("La conformation préparée ne satisfait pas l’état OFF requis pour la FFS. Consultez le contrôle des appariements."))
        runtime = dict(request["runtime"])
        for kind in PARAMETERS:
            parameter = output / "preparation" / f"seq_dep_file_{kind}.txt"
            expected = runtime.get("parameter_sha256", {}).get(kind)
            if expected and _sha(parameter) != expected:
                raise ValueError(tr("Les paramètres du moteur oxNA ont changé depuis la préparation du calcul."))
            runtime["seq_dep_file_" + kind] = str(parameter)
        worker_request = {**request, "runtime": runtime, "prepared": prepared}
        _write(output / "worker_request.json", worker_request)
        worker = Path(__file__).with_name("oxdna_ffs_worker.py")
        command = [request["runtime"]["python"], str(worker), str(output / "worker_request.json")]
        result_path = output / "result.json"
        old_result_stamp = result_path.stat().st_mtime_ns if result_path.is_file() else None
        with (output / "engine.log").open("ab") as log:
            process = subprocess.Popen(command, cwd=output, stdout=log, stderr=log, stdin=subprocess.DEVNULL,
                                       env=external_process_environment(), **process_group_options())
            last_progress, stopping = "", None
            while process.poll() is None:
                now = time.monotonic()
                if cancel_check and cancel_check():
                    if stopping is None:
                        (output / "CANCEL").touch()
                        stopping = now
                    if now - stopping > 30:
                        stop_process(process)
                        break
                path = output / "progress.json"
                if path.is_file() and path.stat().st_size < 32 * 1024 * 1024:
                    try:
                        text = path.read_text()
                        if text != last_progress:
                            data = json.loads(text)
                            if progress_callback:
                                progress_callback(min(total, int(data.get("current", 1))), total, data.get("message", "FFS"))
                            last_progress = text
                    except (OSError, ValueError, TypeError):
                        pass
                time.sleep(.1)
        path = output / "result.json"
        if path.is_file() and path.stat().st_mtime_ns != old_result_stamp:
            result = json.loads(path.read_text())
            if process.returncode and result.get("status") == "completed":
                result.update(status="failed", error=tr("Le moteur oxDNA s’est arrêté ; consultez engine.log."))
        else:
            result = {**previous, "status": "failed",
                      "error": tr("Le moteur oxDNA s’est arrêté ; consultez engine.log.")}
            diagnostic = output / "resume_error.json"
            if diagnostic.is_file():
                result["error"] = json.loads(diagnostic.read_text()).get("error") or result["error"]
        if cancel_check and cancel_check():
            result["status"] = "cancelled"
    except Exception as exc:
        cancelled = bool(cancel_check and cancel_check())
        result = {"status": "cancelled" if cancelled else "failed",
                  "error": str(exc)}
    finally:
        if process is not None and process.poll() is None:
            stop_process(process)
    return _finish_result(result, request, output, started - prior_elapsed, prepared)


def _finish_result(result, request, output, started, prepared):
    if prepared is not None:
        result["preparation"] = prepared
        config = Path(prepared["configuration_file"])
        if config.is_file() and config.resolve().is_relative_to(output.resolve()):
            result["prepared_movie_file"] = config.resolve().relative_to(output.resolve()).as_posix()
    movie = result.get("trajectory_file") or result.get("prepared_movie_file", "")
    result.update(model="NA",
                  elapsed_seconds=max(float(result.get("elapsed_seconds") or 0.), time.monotonic() - started),
                  visualization=request.get("visualization"),
                  guide_source=request.get("guide_source"), guide_conditions=request.get("guide_conditions"),
                  guide_provenance=request.get("guide_provenance", {}),
                  switch_guide=request.get("switch_guide"), trigger_guide=request.get("trigger_guide"),
                  topology_file="preparation/system.top", movie_file=movie,
                  movie_kind="ffs" if result.get("trajectory_file") else "prepared")
    frozen_defaults = {"sampling": request["sampling"], "plan": request["plan"],
                       "box_size": request["box_size"], "box_policy": request.get("box_policy"),
                       "concentration_each_m": request["concentration_each_m"], "rate_estimate": None}
    # Retain the native controller's interfaces, finite volume and rate report.
    for key, value in frozen_defaults.items():
        result.setdefault(key, value)
    result.setdefault("kinetics", request["kinetics"])
    result["method"] = "NATIVE_FFS"
    result["simulation_mode"] = request["kinetics"]["mode"]
    _write(output / "result.json", result)
    return result
