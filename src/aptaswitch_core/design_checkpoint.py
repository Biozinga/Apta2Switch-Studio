"""Atomic, verifiable checkpoints at completed design-trial boundaries."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, replace
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

from .models import RunResult
from .localization import tr
from .run_storage import run_output_directory
from .scoring import sort_candidates
from .thermodynamics.optimization import TrialOptions, rbs_stability_settings, translation_search_settings


CHECKPOINT_NAME = 'design_checkpoint.json'
CHECKPOINT_VERSION = 1
PHASES = {'design', 'inverse_extension', 'exports', 'done'}


def control_action(control_check):
    action = control_check() if control_check else None
    if action not in (None, 'pause', 'finish'):
        raise ValueError('Design control must return None, pause or finish.')
    return action


def _digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def request_digest(request):
    payload = request.normalized().to_dict()
    # A rebuilt app may load identical scientific parameters at another path.
    # The independently checked backend provenance must remain identical.
    payload.pop('thermodynamics_library_path', None)
    payload['output_dir'] = str(run_output_directory(payload['output_dir'], 'design').resolve())
    return _digest(payload)


def _payload_digest(payload):
    payload = deepcopy(payload)
    payload['optimization']['checkpoint'].pop('payload_digest', None)
    # The archive reader resolves storage and runtime locations for this
    # machine. Scientific identity and run location are checked independently.
    payload['request'].pop('thermodynamics_library_path', None)
    payload['request']['output_dir'] = str(run_output_directory(
        payload['request']['output_dir'], 'design').resolve())
    payload['output_dir'] = str(Path(payload['output_dir']).resolve())
    return _digest(payload)


def seal_design_checkpoint(result: RunResult, *, phase=None, final_status=None) -> None:
    """Update identity and integrity metadata without publishing unfinished exports."""
    metadata = result.optimization.setdefault('checkpoint', {})
    metadata.update(version=CHECKPOINT_VERSION,
                    completed_trials=len(result.optimization['trials']),
                    request_digest=request_digest(result.request),
                    thermodynamics_digest=_digest(result.thermodynamics))
    if phase is not None:
        if phase not in PHASES:
            raise ValueError('Unknown design checkpoint phase.')
        metadata['phase'] = phase
    if final_status is not None:
        if final_status not in ('completed', 'stopped'):
            raise ValueError('Unknown final design status.')
        metadata['final_status'] = final_status
    metadata.setdefault('final_status', 'completed')
    metadata['payload_digest'] = _payload_digest(result.to_dict())


def save_design_checkpoint(result: RunResult, *, phase=None, final_status=None, callback=None) -> Path:
    """Persist one complete boundary before publishing it to the UI callback."""
    seal_design_checkpoint(result, phase=phase, final_status=final_status)
    payload = json.dumps(result.to_dict(), indent=2, allow_nan=False)
    directory = Path(result.output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory/CHECKPOINT_NAME
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=directory,
                                         prefix='.design_checkpoint.', suffix='.tmp', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(payload+'\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
        # Keep the rename durable where directory fsync is supported.
        try:
            descriptor = os.open(directory, os.O_RDONLY)
            try:
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        except OSError:
            pass
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    if callback:
        callback(deepcopy(result))
    return target


def validate_resume_result(source: RunResult, request, thermodynamics: dict) -> RunResult:
    """Accept only the latest intact checkpoint for this request and model."""
    if not isinstance(source, RunResult):
        raise ValueError('Resuming design requires a loaded RunResult checkpoint.')
    result = deepcopy(source)
    optimization = result.optimization
    if not isinstance(optimization, dict):
        raise ValueError('The checkpoint optimization metadata is invalid.')
    metadata = optimization.get('checkpoint', {})
    if (not isinstance(metadata, dict) or type(metadata.get('version')) is not int
            or metadata.get('version') != CHECKPOINT_VERSION or metadata.get('phase') not in PHASES
            or metadata.get('final_status') not in ('completed', 'stopped')):
        raise ValueError('This result is not a supported design checkpoint.')
    if metadata['phase'] == 'done' or result.status not in ('running', 'paused', 'cancelled', 'stopped', 'completed', 'failed'):
        raise ValueError('This design has no unfinished checkpoint to resume.')
    if (request_digest(request) != metadata.get('request_digest')
            or request_digest(result.request) != metadata['request_digest']):
        raise ValueError('The design settings differ from the saved checkpoint.')
    if (_digest(thermodynamics) != metadata.get('thermodynamics_digest')
            or result.thermodynamics != thermodynamics):
        raise ValueError('The thermodynamic engine or parameters differ from the saved checkpoint.')
    output = Path(result.output_dir).resolve()
    expected = (run_output_directory(request.output_dir, 'design')/result.run_id).resolve()
    if output != expected or output.name != result.run_id:
        raise ValueError('The design checkpoint belongs to another output directory.')
    payload = result.to_dict()
    if _payload_digest(payload) != metadata.get('payload_digest'):
        raise ValueError('The design checkpoint contents failed their integrity check.')
    disk = json.loads((output/CHECKPOINT_NAME).read_text(encoding='utf-8'))
    disk_metadata = disk.get('optimization', {}).get('checkpoint', {})
    if (_payload_digest(disk) != disk_metadata.get('payload_digest')
            or disk_metadata.get('payload_digest') != metadata['payload_digest']):
        raise ValueError('A newer or different design checkpoint is already saved.')
    if request.fixed_switch_sequence is not None:
        evaluations = optimization.get('evaluations')
        if (optimization.get('algorithm') != 'fixed_sequence_evaluation'
                or optimization.get('objective') != 'supplied_sequences'
                or optimization.get('sequence_optimization') is not False
                or optimization.get('trials') != []
                or metadata.get('completed_trials') != 0
                or type(evaluations) is not int or evaluations not in (0, 1)
                or len(result.candidates) != evaluations):
            raise ValueError('The fixed-sequence evaluation checkpoint is inconsistent.')
        if result.candidates:
            candidate = result.candidates[0]
            if (candidate.trial != 1 or candidate.excluded
                    or candidate.switch_seq != request.fixed_switch_sequence
                    or candidate.trigger_seq != request.sequences.trigger_dna):
                raise ValueError('The fixed-sequence checkpoint changed the supplied sequences.')
        if metadata['phase'] != 'design' and evaluations != 1:
            raise ValueError('The fixed-sequence checkpoint omits its evaluation.')
        result.request, result.output_dir = request.normalized(), output
        return result
    options = TrialOptions(**request.optimization_options)
    if (optimization.get('rbs_stability') != rbs_stability_settings(request.scoring.rbs)
            or optimization.get('translation_search') != translation_search_settings(request.scoring, request.inverse_extension is not None)):
        raise ValueError(tr('Ce checkpoint utilise une autre méthode de conception. Ses candidats restent consultables ; lancez un nouveau design pour utiliser la méthode actuelle.'))
    if (optimization.get('algorithm') != 'diverse_local_first'
            or optimization.get('objective') != 'functional_regions_then_tube_defect'
            or optimization.get('seed_rule') != 'seed + trial - 1'
            or optimization.get('functional_weights') != asdict(request.functional_defect_weights)
            or optimization.get('options') != asdict(options)):
        raise ValueError('The checkpoint optimization settings are inconsistent.')
    reports = optimization.get('trials')
    completed = metadata.get('completed_trials')
    if (not isinstance(reports, list) or type(completed) is not int
            or completed != len(reports) or not 0 <= completed <= request.trials):
        raise ValueError('The checkpoint completed-trial count is invalid.')
    for trial, report in enumerate(reports, 1):
        if (not isinstance(report, dict) or report.get('trial') != trial or report.get('stop_reason') == 'cancelled'
                or report.get('options') != asdict(replace(options, seed=options.seed+trial-1))):
            raise ValueError('The checkpoint trial sequence or random seeds are inconsistent.')
    trials = [candidate.trial for candidate in result.candidates]
    if (len(set(trials)) != len(trials)
            or any(type(trial) is not int or not 1 <= trial <= completed for trial in trials)
            or trials != [candidate.trial for candidate in sort_candidates(result.candidates,
                                  defect_weight=request.scoring.defect)]):
        raise ValueError('The checkpoint candidates are inconsistent or duplicated.')
    if (metadata['phase'] != 'design' and metadata['final_status'] == 'completed'
            and completed != request.trials):
        raise ValueError('This checkpoint omits unfinished design trials.')
    result.request = request.normalized()
    result.output_dir = output
    return result
