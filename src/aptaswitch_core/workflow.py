"""High-level run workflow."""

from __future__ import annotations

from .localization import tr

from .engines import run_design
from .exports import export_run
from .design_checkpoint import (CHECKPOINT_NAME, control_action, save_design_checkpoint,
                                seal_design_checkpoint)


def reanalyze_inverse_and_export(source, request, progress_callback=None, cancel_check=None):
    """Recompare saved switches without running sequence design again."""
    from copy import deepcopy
    from dataclasses import replace
    from datetime import datetime
    from .inverse_extension import analyze_inverse_extensions
    from .thermodynamics.setup import backend_for_request
    from .run_storage import new_run_directory
    if source.request.inverse_extension is None or request.inverse_extension is None:
        raise ValueError(tr("Ce calcul ne contient pas d’extension conçue avec le switch."))
    result = deepcopy(source)
    # This is a new comparison-only run, not continuation of the source job.
    result.optimization.pop('checkpoint', None)
    # The new comparison must not relabel the retained switch-design metrics.
    request.inverse_extension.validate(len(source.request.sequences.trigger_dna))
    result.request = replace(source.request, inverse_extension=request.inverse_extension,
                             output_dir=request.output_dir,
                             thermodynamics_library_path=request.thermodynamics_library_path)
    started = datetime.now().astimezone()
    result.run_id, result.output_dir = new_run_directory(request.output_dir, "design", started)
    result.started_at = started.isoformat()
    analyze_inverse_extensions(result, backend=backend_for_request(request),
                               progress_callback=progress_callback, cancel_check=cancel_check)
    result.inverse_extension_analysis["source_run_id"] = source.run_id
    result.inverse_extension_analysis["switch_design_repeated"] = False
    result.completed_at = datetime.now().astimezone().isoformat()
    return result, export_run(result)


def run_and_export(request, progress_callback=None, cancel_check=None, *, control_check=None,
                   resume_result=None, checkpoint_callback=None):
    """Resume completed trial boundaries, then finish analysis and exports once."""
    from datetime import datetime

    result = run_design(request, progress_callback=progress_callback, cancel_check=cancel_check,
                        control_check=control_check, resume_result=resume_result,
                        checkpoint_callback=checkpoint_callback)
    metadata = result.optimization['checkpoint']

    def paused():
        if result.status == 'paused' or control_action(control_check) == 'pause':
            result.status = 'paused'
            result.completed_at = datetime.now().astimezone().isoformat()
            save_design_checkpoint(result, callback=checkpoint_callback)
            return True
        return False

    if result.status != 'cancelled' and paused():
        return result, {'checkpoint': str(result.output_dir/CHECKPOINT_NAME)}
    if metadata['phase'] == 'inverse_extension' and result.status != 'cancelled':
        from .inverse_extension import analyze_inverse_extensions
        from .thermodynamics.setup import backend_for_request
        result.status = 'running'
        save_design_checkpoint(result, phase='inverse_extension', callback=checkpoint_callback)
        # This stage is resumed independently of sequence design. A pause waits
        # for the comparison stage; the existing hard cancellation stays active.
        analyze_inverse_extensions(result,
                                   backend=backend_for_request(result.request),
                                   progress_callback=progress_callback, cancel_check=cancel_check)
        result.completed_at = datetime.now().astimezone().isoformat()
        if result.inverse_extension_analysis.get('status') == 'cancelled':
            result.status = 'cancelled'
            save_design_checkpoint(result, phase='inverse_extension', callback=checkpoint_callback)
        elif result.inverse_extension_analysis.get('status') != 'completed':
            # The comparison records its own failures rather than raising.
            # Preserve them at this phase and explicitly fail the outer job;
            # otherwise a failed comparison would be sealed as a finished run.
            result.status = 'failed'
            result.inverse_extension_analysis['status'] = 'failed'
            save_design_checkpoint(result, phase='inverse_extension', callback=checkpoint_callback)
            message = tr('La comparaison des extensions inverses a échoué. Les essais de conception sont sauvegardés ; reprenez le calcul pour réessayer la comparaison.')
            detail = result.inverse_extension_analysis.get('error')
            raise RuntimeError(message + (' '+str(detail) if detail else ''))
        else:
            save_design_checkpoint(result, phase='exports', callback=checkpoint_callback)
            if paused():
                return result, {'checkpoint': str(result.output_dir/CHECKPOINT_NAME)}
    if result.status != 'cancelled':
        result.status = metadata['final_status']
        result.completed_at = datetime.now().astimezone().isoformat()
        # On export failure the persisted checkpoint remains at 'exports'.
        # Only publish 'done' after every requested export has been written.
        save_design_checkpoint(result, phase='exports', callback=checkpoint_callback)
        seal_design_checkpoint(result, phase='done')
    if progress_callback:
        progress_callback(len(result.candidates) if request.fixed_switch_sequence is not None else len(result.optimization['trials']), request.trials, tr('Écriture des exports…'))
    exported = export_run(result)
    if result.status != 'cancelled':
        save_design_checkpoint(result, phase='done', callback=checkpoint_callback)
    if progress_callback:
        progress_callback(len(result.candidates) if request.fixed_switch_sequence is not None else len(result.optimization['trials']), request.trials, tr('Exports terminés.'))
    return result, exported
