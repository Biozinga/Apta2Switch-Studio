"""Synthetic data for rendering/export tests only, never a design engine.

The fixed metrics below are arbitrary test values. The structures are target
schematics, not calculated predictions; this module is outside ``src`` and is
excluded from the packaged application.
"""

from aptaswitch_core.models import DesignCandidate, RunRequest, RunResult
from aptaswitch_core.preview import build_switch_model


def design_result_fixture(request: RunRequest) -> RunResult:
    request = request.normalized()
    model = build_switch_model(
        request.architecture, request.sequences.effective_trigger(),
        reporter_sequence=request.reporter.sequence_after_start_rna,
        reporter_label=request.reporter.label,
        switch_material=request.thermo.switch_material,
        trigger_material=request.thermo.trigger_material,
    )
    candidates = [
        DesignCandidate(
            trial=trial, switch_seq=model.switch_sequence.replace("N", "A"),
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
    )
