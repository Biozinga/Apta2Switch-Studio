"""Exploratory OFF/ON translation indicator using unmodified OSTIR 1.1.3.

OFF is the ordinary RNA-only prediction. ON is an ideal occupancy projection:
the switch bases bound by the trigger cannot pair with the RNA or ribosome.
It neither calculates a DNA duplex nor predicts induction in an experiment.
The original RNA sequence is retained; OSTIR's own +/-35 nt AUG window is
made explicit before constraints are applied, avoiding upstream offset bugs.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import importlib
from importlib.metadata import version
import json
import math
from numbers import Real
from pathlib import Path
import warnings


OSTIR_VERSION = "1.1.3"
VIENNA_VERSION = "2.7.2"
MODEL_METHOD = "ostir_ideal_occupancy_v1"
ANTI_SHINE_DALGARNO = "ACCUCCUUA"
BETA = 0.40002512
LOG_K = 7.279194329
PARAMETER_SHA256 = "0c671a7455f5686bc160fc655194b11888594ab7d506626cdad1470bcadd17d9"
FORMULA = "on/(on+K) * max(0,(on-off)/(on+off))"
SCOPE = "RNA-only OFF versus ideal trigger-occupied RNA projection; no DNA binding energy or concentration"
NORMALIZATION = "K is a mathematical pivot, not an experimentally validated ON target; the merit is a ranking heuristic"
ENERGY_FIELDS = ("dG_total", "dG_rRNA:mRNA", "dG_mRNA", "dG_spacing", "dG_standby",
                 "dG_start_codon", "RBS_distance_bp")


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return value


def _sigmoid(value):
    if value >= 0:
        return 1. / (1. + math.exp(-value))
    exp = math.exp(value)
    return exp / (1. + exp)


def compute_rbs_merit(off_dg: float, on_dg: float) -> dict:
    """Return a bounded ranking heuristic, not an experimental calibration.

    qON = on/(on+K); qContrast = max(0,(on-off)/(on+off)); merit = qON*qContrast.
    K is OSTIR's fitted rate at dG_total=0, used only as a normalization pivot.
    Log arithmetic retains contrast even when reported rounded rates are zero.
    An underflow/overflow cannot silently become a valid zero/infinite rate.
    """
    off_dg, on_dg = _finite(off_dg, "OFF energy"), _finite(on_dg, "ON energy")
    off_log, on_log = LOG_K - BETA * off_dg, LOG_K - BETA * on_dg
    difference = on_log - off_log
    if not all(math.isfinite(v) for v in (off_log, on_log, difference)):
        raise ValueError("OSTIR energies exceed the supported numerical range")
    try:
        off_rate, on_rate = math.exp(off_log), math.exp(on_log)
    except OverflowError as error:
        raise ValueError("OSTIR rates exceed the supported numerical range") from error
    if not all(math.isfinite(v) and v > 0 for v in (off_rate, on_rate)):
        raise ValueError("OSTIR rates must be positive finite numbers")
    on_fraction = _sigmoid(-BETA * on_dg)
    contrast = max(0., math.tanh(difference / 2.))
    return {
        "off_rate": off_rate, "on_rate": on_rate,
        "log10_on_off": difference / math.log(10.),
        "on_fraction": on_fraction, "contrast": contrast,
        "merit": on_fraction * contrast, "pivot_rate": math.exp(LOG_K),
        "formula": FORMULA,
    }


@lru_cache(maxsize=1)
def _runtime():
    try:
        actual = version("ostir")
        RNA = importlib.import_module("RNA")
        with warnings.catch_warnings():
            # OSTIR checks obsolete CLI executables at import, although all
            # calculations here use the verified Python RNA binding.
            warnings.filterwarnings("ignore", message="RBS Calculator Vienna is missing dependency ViennaRNA!",
                                    category=UserWarning, module=r"ostir\.ViennaRNA")
            factory = importlib.import_module("ostir.ostir_factory")
            calculations = importlib.import_module("ostir.ostir_calculations")
            folding = importlib.import_module("ostir.ViennaRNA")
    except (ImportError, ModuleNotFoundError) as error:
        raise RuntimeError("OSTIR 1.1.3 and ViennaRNA 2.7.2 are required for the RBS indicator") from error
    if actual != OSTIR_VERSION or RNA.__version__ != VIENNA_VERSION:
        raise RuntimeError("The RBS indicator requires OSTIR 1.1.3 and ViennaRNA 2.7.2")
    constants = calculations.ostir_constants
    if (constants.Beta != BETA or constants.logK != LOG_K or constants.temp != 37.
            or constants.cutoff != 35 or folding.vienna_constants.material != "rna2004"):
        raise RuntimeError("OSTIR's reference model has been modified")
    parameter_file = Path(folding.get_paramater_file("rna2004"))
    parameter_digest = hashlib.sha256(parameter_file.read_bytes()).hexdigest()
    if parameter_digest != PARAMETER_SHA256:
        raise RuntimeError("OSTIR's reference RNA parameter file has been modified")
    return RNA, factory.OSTIRFactory, folding, parameter_file, parameter_digest


def _model_result(factory_type, sequence, aug_start, constraints, dangles):
    factory = factory_type(sequence, [aug_start + 1, aug_start + 1],
                           ANTI_SHINE_DALGARNO, constraints, verbosity=0,
                           decimal_places=8, circular=False, name="toehold_switch")
    # Preserve the full-sequence upstream rule, even though our local AUG is
    # always <=35 and would otherwise change dangles from none to all.
    factory.auto_dangles = False
    factory.dangles_default = dangles
    factory.threads = 1
    factory.calc_dG()
    if len(factory.results) != 1:
        from .localization import tr
        raise ValueError(tr("OSTIR ne trouve aucun site d’initiation calculable autour de l’AUG choisi ({state}). Aucun score nul n’est inventé. Vérifiez le RBS personnalisé ou mettez son poids RBS ON/OFF à 0 pour concevoir sans ce critère.",
                            state="ON" if constraints else "OFF"))
    result = factory.results[0].results()
    if result.get("start_position") != aug_start + 1 or result.get("start_codon") != "AUG":
        raise ValueError("OSTIR returned a different initiation site")
    return {key: _finite(result[key], key) for key in ENERGY_FIELDS}


def _reference(dangles):
    return {
        "ostir_version": OSTIR_VERSION, "vienna_version": VIENNA_VERSION,
        "parameter_set": "RNA-Turner2004", "parameter_sha256": PARAMETER_SHA256,
        "celsius": 37., "sodium_m": 1.021, "magnesium_m": 0.,
        "anti_shine_dalgarno": ANTI_SHINE_DALGARNO,
        "dangles": dangles, "window_radius_nt": 35,
        "beta": BETA, "log_k": LOG_K, "rate_units": "OSTIR arbitrary units",
    }


def _input_context(sequence, aug_start, occupied_indices):
    if not isinstance(sequence, str) or not sequence or set(sequence) - set("ACGU"):
        raise ValueError("The RBS indicator requires the actual uppercase RNA sequence")
    if type(aug_start) is not int or not 0 < aug_start <= len(sequence) - 3:
        raise ValueError("A zero-based AUG coordinate with an upstream RNA region is required")
    if sequence[aug_start:aug_start + 3] != "AUG":
        raise ValueError("The selected initiation site must be AUG")
    if (not isinstance(occupied_indices, tuple)
            or any(type(i) is not int or not 0 <= i < aug_start for i in occupied_indices)
            or len(set(occupied_indices)) != len(occupied_indices)):
        raise ValueError("Occupied sites must be unique zero-based RNA indices upstream of AUG")
    occupied_indices = tuple(sorted(occupied_indices))
    begin, end = max(0, aug_start - 35), min(len(sequence), aug_start + 35)
    window = sequence[begin:end]
    local_occupied = [i - begin for i in occupied_indices if begin <= i < end]
    return {
        "sequence": sequence, "aug_start": aug_start, "occupied_indices": list(occupied_indices),
        "window_start": begin, "window_end": end, "window_sequence": window,
        "window_aug_start": aug_start - begin, "window_occupied_indices": local_occupied,
        "on_constraints": "".join("x" if i in local_occupied else "." for i in range(len(window))),
    }


@lru_cache(maxsize=256)
def _prediction(sequence, aug_start, occupied_indices):
    RNA, factory_type, folding, parameter_file, _ = _runtime()
    context = _input_context(sequence, aug_start, occupied_indices)
    window, local_aug = context["window_sequence"], context["window_aug_start"]
    mask = context["on_constraints"]
    dangles = "none" if aug_start > 35 else "all"
    # OSTIR changes ViennaRNA's process-global energy tables. Cooperate with
    # both homogeneous adapters, and reload before every paired prediction:
    # OSTIR caches model details, which otherwise leaves prior DNA tables live.
    from .thermodynamics.vienna import _NATIVE_LOCK
    with _NATIVE_LOCK:
        RNA.params_load(str(parameter_file))
        details = folding.get_paramater_object("rna2004", 37., dangles)
        if details.salt != 1.021 or details.temperature != 37.:
            raise RuntimeError("OSTIR's RNA reference conditions have been modified")
        off = _model_result(factory_type, window, local_aug, None, dangles)
        on = (_model_result(factory_type, window, local_aug, mask, dangles)
              if context["window_occupied_indices"] else dict(off))
    metrics = compute_rbs_merit(off["dG_total"], on["dG_total"])
    result = {
        "schema_version": 1, "method": MODEL_METHOD,
        **context, "off": off, "on": on, **metrics,
        "reference": _reference(dangles), "scope": SCOPE, "normalization": NORMALIZATION,
    }
    # Cache immutable serialized data: UI callers cannot corrupt future scores.
    return json.dumps(result, allow_nan=False, separators=(",", ":"))


def predict_rbs_switch(sequence: str, aug_start: int,
                       occupied_indices: tuple[int, ...]) -> dict:
    """Predict a reference OFF/ideal-ON indicator without changing the candidate.

    All indices are zero-based in the complete RNA switch. Occupied sites must
    be upstream of the chosen AUG; the caller supplies only the architecture's
    trigger recognition sites, never forcibly opens RBS/AUG to improve a score.
    Sites outside OSTIR's native window are retained in provenance but cannot
    influence its local prediction. No experiment or hybrid rate is inferred.
    """
    context = _input_context(sequence, aug_start, occupied_indices)
    return json.loads(_prediction(sequence, aug_start, tuple(context["occupied_indices"])))


def validate_rbs_prediction(data: dict, sequence: str, aug_start: int,
                            occupied_indices: tuple[int, ...]) -> None:
    """Validate recorded context and derived values without importing a runtime.

    This checks internal consistency and provenance, not the truth of a saved
    energy calculation; it never recalculates or repairs historical results.
    """
    context = _input_context(sequence, aug_start, occupied_indices)
    expected = {
        "schema_version": 1, "method": MODEL_METHOD, **context,
        "reference": _reference("none" if aug_start > 35 else "all"),
        "scope": SCOPE, "normalization": NORMALIZATION,
    }
    if not isinstance(data, dict):
        raise ValueError("Recorded RBS prediction must be an object")
    try:
        json.dumps(data, allow_nan=False)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("Recorded RBS prediction is not strict JSON") from error
    for key, value in expected.items():
        # Type-sensitive JSON comparison excludes booleans masquerading as
        # integer coordinates while allowing JSON numeric 37 versus 37.0.
        actual = data.get(key)
        if actual != value or (type(value) is int and type(actual) is not int):
            raise ValueError(f"Recorded RBS prediction has a different {key}")
    for name in ("occupied_indices", "window_occupied_indices"):
        if not isinstance(data[name], list) or any(type(i) is not int for i in data[name]):
            raise ValueError("Recorded RBS indices must be integer arrays")
    for key in ("celsius", "sodium_m", "magnesium_m", "window_radius_nt", "beta", "log_k"):
        _finite(data["reference"][key], key)
    for state in ("off", "on"):
        row = data.get(state)
        if not isinstance(row, dict) or set(row) != set(ENERGY_FIELDS):
            raise ValueError("Recorded RBS energy components are incomplete or unrecognized")
        for key in ENERGY_FIELDS:
            _finite(row[key], key)
        total = row["dG_rRNA:mRNA"] + row["dG_spacing"] + row["dG_start_codon"] - row["dG_mRNA"] - row["dG_standby"]
        if not math.isclose(row["dG_total"], total, rel_tol=1e-12, abs_tol=1e-6):
            raise ValueError("Recorded OSTIR energy components do not sum to its total")
    metrics = compute_rbs_merit(data["off"]["dG_total"], data["on"]["dG_total"])
    if set(data) != set(expected) | set(metrics) | {"off", "on"}:
        raise ValueError("Recorded RBS fields are incomplete or unrecognized")
    if data["formula"] != FORMULA:
        raise ValueError("Recorded RBS normalization formula differs")
    for key, value in metrics.items():
        if key != "formula" and not math.isclose(_finite(data[key], key), value, rel_tol=1e-12, abs_tol=0.):
            raise ValueError(f"Recorded RBS {key} differs from its energies")
    if not context["window_occupied_indices"] and data["on"] != data["off"]:
        raise ValueError("An unoccupied ON projection must equal OFF")
