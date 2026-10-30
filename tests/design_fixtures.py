"""Synthetic data for rendering/export tests only, never a design engine.

The fixed metrics below are arbitrary test values. The structures are target
schematics, not calculated predictions; this module is outside ``src`` and is
excluded from the packaged application.
"""

from dataclasses import asdict, replace

from aptaswitch_core.models import DesignCandidate, RunRequest, RunResult
from aptaswitch_core.preview import build_switch_model
from aptaswitch_core.thermodynamics.optimization import TrialOptions


def synthetic_native_metadata(request: RunRequest) -> dict:
    """Explicitly identified fixture metadata, never a scientific measurement."""
    if request.engine != "native-hybrid":
        return {"thermodynamics": {}, "optimization": {}}
    options = TrialOptions(**request.optimization_options)
    return {
        "thermodynamics": {
            "engine": request.engine,
            "backend": "synthetic-rendering-fixture",
            "version": "fixture-1",
            "model": request.thermo.material,
            "parameter_digest": "fixture:synthetic-no-physical-parameters",
            "ensemble": "synthetic target schematics; no calculated ensemble",
            "standard_state": "synthetic fixture; no equilibrium prediction",
            "conditions": {"celsius": request.thermo.temperature_c,
                           "sodium": request.thermo.sodium_m, "magnesium": request.thermo.magnesium_m},
        },
        "optimization": {
            "algorithm": "diverse_local_first",
            "objective": "functional_regions_then_tube_defect",
            "seed_rule": "seed + trial - 1",
            "functional_weights": asdict(request.functional_defect_weights),
            "analysis_cache": {},
            "options": asdict(options),
            "trials": [{"trial": trial, "options": asdict(replace(options, seed=options.seed + trial - 1)),
                        "evaluations": 1, "proposals": 0, "cache_hits": 0,
                        "phase": "initialization", "population_evaluated": 1, "lineages_started": 0,
                        "best_global_loss": 0., "best_objective": {
                            "local_defects": [0.] * 4, "global_defect": 0., "guard_defects": [0.] * 3,
                            "local_weights": list(asdict(request.functional_defect_weights).values()), "site_weights": [], "global_site_weights": [],
                        },
                        "selected_assignment": None, "candidate_screened": 1,
                        "selected_objective": {
                            "local_defects": [0.] * 4, "global_defect": 0., "guard_defects": [0.] * 3,
                            "local_weights": list(asdict(request.functional_defect_weights).values()), "site_weights": [], "global_site_weights": [],
                        },
                        "local_limits": None, "guard_limits": None, "elites": [],
                        "stop_reason": "target_reached", "budget_kind": None,
                        "best_loss": 0., "fixture": "synthetic, not calculated"}
                       for trial in range(1, request.trials + 1)],
        },
    }


def design_result_fixture(request: RunRequest) -> RunResult:
    # This fixture contains synthetic tube metrics, no OSTIR prediction.
    request = replace(request, scoring=replace(request.scoring, rbs_calculator=0.))
    request = request.normalized()
    model = build_switch_model(
        request.architecture, request.sequences.effective_trigger(),
        reporter_sequence=request.reporter.sequence_after_start_rna,
        reporter_label=request.reporter.label,
    )
    candidates = [
        DesignCandidate(
            trial=trial, switch_seq=model.switch_sequence.replace("X", "A"),
            trigger_seq=model.trigger_sequence, score=50.0, defect=0.1,
            on_yield_pct=60.0, leak_pct=10.0, delta_g_off=-15.0,
            delta_g_on=-25.0, ddg_activation=-8.0, delta_g_rbs_linker=-2.0,
            has_bad_rbs_linker=False, has_stop_codon=False,
            structure_off=model.off_structure, structure_on=model.on_structure,
            reporter_label=request.reporter.label,
        )
        for trial in range(1, request.trials + 1)
    ]
    return RunResult(
        run_id="test_fixture", request=request, candidates=candidates,
        output_dir=request.output_dir / "test_fixture",
        started_at="2026-01-01T00:00:00+00:00",
        completed_at="2026-01-01T00:00:01+00:00",
        warnings=["Synthetic rendering/export test fixture; no predictions."],
        **synthetic_native_metadata(request),
    )
