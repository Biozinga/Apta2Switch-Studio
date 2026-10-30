"""Application design dispatch and persistent native-model provenance."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime
import json

from .localization import tr
from .models import RunResult, STOP_FILTER_DISABLED_WARNING
from .run_storage import new_run_directory
from .design_checkpoint import control_action, save_design_checkpoint, validate_resume_result
from .scoring import sort_candidates
from .thermodynamics.design_engine import design_candidates
from .thermodynamics.optimization import TrialOptions, rbs_stability_settings, translation_search_settings
from .thermodynamics.types import Conditions


def _run_fixed_reading(request, backend, thermodynamics, *, progress_callback=None,
                       cancel_check=None, control_check=None, resume_result=None,
                       checkpoint_callback=None):
    """Evaluate one supplied assignment with the existing candidate calculator."""
    from .thermodynamics.design_constraints import DesignAssignment, build_design_problem
    from .thermodynamics.design_engine import evaluate_candidate
    from .thermodynamics.types import Strand

    if resume_result is None:
        started = datetime.now().astimezone()
        run_id, output_dir = new_run_directory(request.output_dir, "design", started)
        result = RunResult(run_id=run_id, request=request, candidates=[],
            output_dir=output_dir.resolve(), started_at=started.isoformat(), completed_at='', status='running',
            thermodynamics=thermodynamics,
            optimization={"algorithm": "fixed_sequence_evaluation", "objective": "supplied_sequences",
                          "sequence_optimization": False, "evaluations": 0, "trials": []})
    else:
        result = validate_resume_result(resume_result, request, thermodynamics)
        result.status, result.completed_at = 'running', ''

    result.warnings = [tr("Modèle natif ADN/ARN expérimental : les paramètres des boucles hybrides sont approximatifs ; le classement calculé ne valide pas les performances expérimentales.")]
    phase = 'exports' if result.optimization['evaluations'] else 'design'
    save_design_checkpoint(result, phase=phase, callback=checkpoint_callback)

    def stopped():
        if cancel_check and cancel_check():
            return 'cancelled'
        action = control_action(control_check)
        return 'paused' if action == 'pause' or (action == 'finish' and not result.optimization['evaluations']) else None

    status = stopped()
    if status is None and not result.optimization['evaluations']:
        if progress_callback:
            progress_callback(0, 1, tr("Moteur natif ADN/ARN · lecture des séquences fournies, sans optimisation…"))
        problem = build_design_problem(request)
        assignment = DesignAssignment((
            Strand(request.sequences.trigger_dna, "dna", "Trigger"),
            Strand(request.fixed_switch_sequence, "rna", "Switch"),
        ))
        candidate = evaluate_candidate(request, problem, assignment, backend, trial=1)
        result.candidates = [candidate]
        result.optimization['evaluations'] = 1
        phase = 'exports'
        status = stopped()
    if result.candidates:
        candidate = result.candidates[0]
        if candidate.has_stop_codon:
            result.warnings.append(tr("Mode lecture : le switch fourni contient un STOP prématuré ; ses métriques ont été calculées sans modifier sa séquence."))
        from .architecture import architecture_layout
        layout = architecture_layout(request.architecture)
        on_switch = candidate.structure_on.split('+')[-1]
        if any(symbol != '.' for symbol in on_switch[layout['rbs_linker_start']:layout['rbs_linker_end']]):
            result.warnings.append(tr("Mode lecture : la région RBS–linker est appariée dans la structure MFE ON enregistrée ; ses métriques ont été conservées."))
    result.status = status or 'completed'
    result.completed_at = datetime.now().astimezone().isoformat()
    save_design_checkpoint(result, phase=phase, callback=checkpoint_callback)
    if progress_callback:
        progress_callback(result.optimization['evaluations'], 1,
                          tr("Moteur natif ADN/ARN · lecture terminée") if result.status == 'completed' else
                          tr("Moteur natif ADN/ARN · lecture interrompue, état sauvegardé"))
    return result


def run_design(request, progress_callback=None, cancel_check=None, *, control_check=None,
               resume_result=None, checkpoint_callback=None):
    """Compute with the selected native library; never relabel historical runs."""
    from .thermodynamics.setup import ENGINE_ID, MODEL_ID, backend_for_request

    request = request.normalized()
    if request.engine != ENGINE_ID or request.thermo.material != MODEL_ID:
        raise ValueError("New designs require the native-hybrid engine and its declared model; historical runs can be opened without recalculation.")
    backend = backend_for_request(request)
    options = TrialOptions(**request.optimization_options)
    thermodynamics = {"engine": ENGINE_ID, **asdict(backend.provenance),
                      "conditions": asdict(Conditions(request.thermo.temperature_c,
                          request.thermo.sodium_m, request.thermo.magnesium_m))}
    if request.fixed_switch_sequence is not None:
        return _run_fixed_reading(request, backend, thermodynamics,
            progress_callback=progress_callback, cancel_check=cancel_check, control_check=control_check,
            resume_result=resume_result, checkpoint_callback=checkpoint_callback)
    if resume_result is None:
        started = datetime.now().astimezone()
        run_id, output_dir = new_run_directory(request.output_dir, "design", started)
        result = RunResult(run_id=run_id, request=request, candidates=[],
            output_dir=output_dir.resolve(), started_at=started.isoformat(), completed_at='', status='running',
            thermodynamics=thermodynamics,
            optimization={"algorithm": "diverse_local_first", "objective": "functional_regions_then_tube_defect",
                          "options": asdict(options), "seed_rule": "seed + trial - 1",
                          "functional_weights": asdict(request.functional_defect_weights),
                          "rbs_stability": rbs_stability_settings(request.scoring.rbs),
                          "translation_search": translation_search_settings(request.scoring, request.inverse_extension is not None),
                          "analysis_cache": {}, "trials": []})
    else:
        result = validate_resume_result(resume_result, request, thermodynamics)
        result.status = 'running'
        result.completed_at = ''
        if result.optimization['checkpoint']['phase'] != 'design':
            save_design_checkpoint(result, callback=checkpoint_callback)
            return result
    prior_count = len(result.optimization['trials'])
    prior_cache = dict(result.optimization.get('analysis_cache', {}))

    def update_warnings():
        result.warnings = [] if request.exclude_stop_candidates else [tr(STOP_FILTER_DISABLED_WARNING)]
        unmet = sum(report['stop_reason'] not in {"target_reached", "cancelled"}
                    for report in result.optimization['trials'])
        if unmet:
            result.warnings.append(tr("{count} essais se sont terminés sans atteindre les cibles locale et globale ; leurs conditions d’arrêt sont enregistrées.", count=unmet))
        result.warnings.append(tr("Modèle natif ADN/ARN expérimental : les paramètres des boucles hybrides sont approximatifs ; le classement calculé ne valide pas les performances expérimentales."))

    def update_cache(statistics):
        result.optimization['analysis_cache'] = dict(statistics)
        for key in ('backend_calls', 'cache_hits', 'analysis_calls', 'minimum_calls', 'screen_calls', 'partition_calls'):
            result.optimization['analysis_cache'][key] = prior_cache.get(key, 0)+statistics.get(key, 0)

    def completed_trial(trial, candidate, report, statistics):
        if trial != len(result.optimization['trials'])+1 or report.stop_reason == 'cancelled':
            raise ValueError('Completed design trials must form one uninterrupted sequence.')
        result.optimization['trials'].append(json.loads(json.dumps(
            {'trial': trial, **asdict(report)}, allow_nan=False)))
        result.optimization.pop('interrupted_trial', None)
        if candidate is not None:
            result.candidates.append(candidate)
            result.candidates = sort_candidates(result.candidates, defect_weight=request.scoring.defect)
        update_cache(statistics)
        update_warnings()
        save_design_checkpoint(result, phase='design', callback=checkpoint_callback)

    update_warnings()
    # Identity, settings and the zero-trial boundary exist before any search.
    save_design_checkpoint(result, phase='design', callback=checkpoint_callback)
    if progress_callback:
        progress_callback(prior_count, request.trials, tr("Moteur natif ADN/ARN · préparation des contraintes et des tubes…"))

    def progress(trial, total, state):
        if progress_callback:
            local = "—" if state.best_loss is None else f"{100 * state.best_loss:.2f} %"
            global_loss = "—" if state.best_global_loss is None else f"{100 * state.best_global_loss:.2f} %"
            phase = tr({"initialization": "Population initiale", "local": "Défauts locaux", "global": "Défaut global et stabilité RBS–linker"}[state.phase])
            progress_callback(trial - 1, total,
                              tr("Essai {trial}/{total} · {phase} · {evaluations}/{maximum} évaluations · défaut local {local} · global {global_loss}",
                                 trial=trial, total=total, phase=phase, evaluations=state.evaluations,
                                 maximum=options.max_evaluations, local=local, global_loss=global_loss)
                              + (tr(" · coût d’affinage {cost:.4f}", cost=state.best_objective.refinement_loss)
                                 if state.best_objective is not None else ""))

    batch = design_candidates(request, backend, options=options,
                              progress_callback=progress, cancel_check=cancel_check,
                              start_trial=prior_count+1, trial_callback=completed_trial,
                              control_check=control_check)
    # The callback persisted every finished trial before any next trial starts.
    # Interrupted optimization or finalist screening must repeat at resume.
    finished_here = len(result.optimization['trials']) - prior_count
    if len(batch.optimization_trials) > finished_here:
        interrupted = batch.optimization_trials[finished_here]
        result.optimization['interrupted_trial'] = {
            'trial': len(result.optimization['trials'])+1,
            'stage': 'optimization' if interrupted.stop_reason == 'cancelled' else 'screening',
            **json.loads(json.dumps(asdict(interrupted), allow_nan=False))}
    else:
        result.optimization.pop('interrupted_trial', None)
    update_cache(batch.analysis_cache)
    result.status = batch.status
    result.completed_at = batch.completed_at
    update_warnings()
    phase = ('design' if batch.status in ('paused', 'cancelled') else
             'inverse_extension' if request.inverse_extension is not None else 'exports')
    save_design_checkpoint(result, phase=phase,
                           final_status='stopped' if batch.status == 'stopped' else 'completed',
                           callback=checkpoint_callback)
    if progress_callback:
        message = {"cancelled": "Moteur natif ADN/ARN · calcul annulé",
                   "paused": "Conception en pause · les essais terminés sont sauvegardés.",
                   "stopped": "Conception arrêtée à la fin de l’essai · classement des résultats…"}.get(
                       result.status, "Moteur natif ADN/ARN · conception terminée")
        progress_callback(len(result.optimization['trials']), request.trials, tr(message))
    return result
