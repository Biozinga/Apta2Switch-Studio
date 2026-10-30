"""Frozen saved-switch DNA–DNA experiments and external oxDNA orchestration."""
from __future__ import annotations

import copy
from dataclasses import asdict
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import secrets
import shutil
import subprocess
import time
import zipfile

from .oxdna_box import estimate_ffs_box
from .oxdna_free_md_metrics import TIME_US
from .scientific_process import external_process_environment, process_group_options, stop_process


def free_md_sequences(snapshot: dict) -> dict:
    """Use the immutable candidate's original and genuinely extended DNA trigger."""
    sequences = snapshot.get('sequences') or {}
    thermo = (snapshot.get('request') or {}).get('thermo') or {}
    if thermo.get('trigger_material') != 'dna':
        raise ValueError('La MD libre ADN–ADN exige un trigger ADN.')
    values = {name: ''.join(str(sequences.get(field) or '').split()).upper()
              for name, field in (('aptamer', 'aptamer'), ('original', 'original_trigger'), ('extended', 'trigger'))}
    if any(not value or set(value) - set('ACGT') or len(value) > 2000 for value in values.values()):
        raise ValueError('Un aptamère ADN et deux triggers ADN valides sont requis (2 000 bases maximum par brin).')
    if len(values['extended']) <= len(values['original']) or values['original'] not in values['extended']:
        raise ValueError('Ce switch ne possède pas de trigger ADN étendu contenant le trigger original.')
    return values


def available_free_md(snapshot: dict) -> bool:
    try:
        free_md_sequences(snapshot)
        return True
    except (ValueError, TypeError, AttributeError):
        return False


def parse_monitored_bases(text: str, length: int) -> list[int]:
    """Accept 1-based comma/space/semicolon separated indices and ranges."""
    if not isinstance(text, str) or type(length) is not int or length < 1:
        raise ValueError('Indiquez les positions des bases de l’aptamère, numérotées à partir de 1.')
    normalized = re.sub(r'\s*[-–]\s*', '-', text.strip())
    tokens = [token for token in re.split(r'[\s,;]+', normalized) if token]
    values = set()
    for token in tokens:
        if not re.fullmatch(r'[0-9]+(?:-[0-9]+)?', token):
            raise ValueError('Positions invalides : utilisez par exemple 3, 7-12, 20.')
        bounds = [int(part) for part in token.split('-')]
        first, last = bounds[0], bounds[-1]
        if not 1 <= first <= last <= length:
            raise ValueError(f'Les positions doivent être comprises entre 1 et {length}, dans l’ordre croissant.')
        values.update(range(first, last + 1))
    if not values:
        raise ValueError('Sélectionnez au moins une base de l’aptamère à surveiller.')
    return sorted(values)


def _number(options, key, default, minimum, maximum, *, integer=False):
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'Paramètre invalide : {key}.')
    if not minimum <= value <= maximum or (integer and int(value) != value):
        raise ValueError(f'{key} doit être compris entre {minimum} et {maximum}.')
    return int(value) if integer else float(value)


def prepare_free_md_request(snapshot: dict, runtime: dict, options: dict, *, backend=None) -> dict:
    """Freeze the native model; folding and oxDNA preparation run off the UI thread."""
    sequences = free_md_sequences(snapshot)
    if not snapshot.get('snapshot_id'):
        raise ValueError('Enregistrez le switch avant de lancer la MD libre.')
    allowed = {'monitored_bases', 'duration_us', 'replicas', 'workers', 'temperature_c', 'salt_m',
               'contact_cutoff_nm', 'seed', 'thermodynamics_library_path'}
    if set(options) - allowed:
        raise ValueError('Paramètre MD inconnu : ' + sorted(set(options) - allowed)[0])
    bases = options.get('monitored_bases', [])
    if isinstance(bases, str):
        bases = parse_monitored_bases(bases, len(sequences['aptamer']))
    if not isinstance(bases, (list, tuple)) or not bases or any(type(i) is not int or
            not 1 <= i <= len(sequences['aptamer']) for i in bases):
        raise ValueError('Les bases surveillées doivent être des positions de l’aptamère, à partir de 1.')
    bases = sorted(set(bases))
    replicas = _number(options, 'replicas', 10, 2, 100, integer=True)
    if replicas % 2:
        raise ValueError('Choisissez un nombre pair de répliques pour une répartition 50 / 50.')
    workers = _number(options, 'workers', min(4, replicas, os.cpu_count() or 1), 1, replicas, integer=True)
    thermo = snapshot['request']['thermo']
    temperature = _number(options, 'temperature_c', thermo.get('temperature_c', 37.), 25., 55.)
    salt = _number(options, 'salt_m', thermo.get('sodium_m', .15), .12, 1.)
    duration = _number(options, 'duration_us', .1, .001, 100.)
    cutoff = _number(options, 'contact_cutoff_nm', 1.5, .1, 10.)
    seed = _number(options, 'seed', secrets.randbelow(2**31), 0, 2**31 - 1, integer=True)
    if not isinstance(runtime, dict) or not runtime.get('python') or not runtime.get('seq_dep_file_DNA'):
        raise ValueError('Configurez le moteur oxDNA et ses paramètres ADN avant de lancer le calcul.')
    from .thermodynamics.setup import get_backend, configured_library_path
    library_path = (options.get('thermodynamics_library_path')
                    or snapshot['request'].get('thermodynamics_library_path') or configured_library_path())
    backend = backend or get_backend(library_path)
    dt = .003
    steps = math.ceil(duration / (dt * TIME_US))
    identity_data = {key: snapshot.get(key) for key in ('snapshot_id', 'saved_at', 'candidate', 'request', 'sequences')}
    identity = hashlib.sha256(json.dumps(identity_data, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
    box = estimate_ffs_box(len(sequences['aptamer']), len(sequences['extended']), 'dna', 'dna', temperature, salt)
    return {'format': 'aptaswitch-oxdna-free-md', 'version': 1, 'method': 'FREE_MD_DNA_DNA',
            'snapshot_id': snapshot['snapshot_id'], 'saved_switch_identity': identity,
            'runtime': copy.deepcopy(runtime), 'sequences': sequences,
            'variants': [{'name': name, 'aptamer_sequence': sequences['aptamer'],
                          'trigger_sequence': sequences[name]} for name in ('original', 'extended')],
            'monitored_bases': bases, 'replicas': replicas, 'workers': workers,
            'sampling': {'dt': dt, 'production_steps': steps, 'sample_every': max(1, math.ceil(steps / 1000)),
                         'equilibration_steps': 50000, 'temperature_c': temperature, 'salt_m': salt,
                         'seed': seed, 'duration_us': duration, 'actual_duration_us': steps * dt * TIME_US},
            'contact_cutoff_nm': cutoff, 'box_policy': box,
            'thermodynamics_library_path': str(library_path or ''),
            'guide_model': backend.provenance.model, 'guide_thermodynamics': asdict(backend.provenance),
            'guide_conditions': {'temperature_c': temperature, 'sodium_m': salt, 'magnesium_m': 0.,
                                 'ensemble': backend.provenance.ensemble},
            'model': 'DNA2', 'time_unit_ps': 3.031, 'length_unit_nm': .8518,
            'ligand_simulated': False, 'production_external_forces': False,
            'time_interpretation': 'nominal_model_time_not_calibrated_experimental_kinetics',
            'initial_state': 'native_hybrid_complex_guided_then_unrestrained_equilibration',
            'equilibrium_validated': False}


def _write(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _validate_guide_request(request):
    """Never relabel an older folding model or alter frozen preparation conditions."""
    from .thermodynamics.setup import MODEL_ID
    expected_conditions = {'temperature_c': request['sampling']['temperature_c'],
                           'sodium_m': request['sampling']['salt_m'], 'magnesium_m': 0.,
                           'ensemble': request.get('guide_thermodynamics', {}).get('ensemble')}
    provenance = request.get('guide_thermodynamics', {})
    if (request.get('guide_model') != MODEL_ID or provenance.get('model') != MODEL_ID
            or provenance.get('backend') != 'AptaSwitch-native-hybrid'
            or request.get('guide_conditions') != expected_conditions
            or request.get('initial_state') != 'native_hybrid_complex_guided_then_unrestrained_equilibration'):
        raise ValueError('Le modèle des guides MD est incompatible ; créez un nouveau calcul avec le moteur intégré.')
    if request.get('guide_source') not in (None, 'native_hybrid_mfe_two_strand_dna_complex'):
        raise ValueError('Un ancien calcul MD ne peut pas être converti vers un autre moteur thermodynamique.')
    if request.get('guide_provenance', provenance) != provenance:
        raise ValueError('La provenance des guides MD diffère du modèle enregistré.')


def prepare_free_md_system(request, output, *, progress_callback=None, cancel_check=None, backend=None):
    """Build connected-complex seeds; native replicas release all guide traps."""
    import numpy as np
    from .thermodynamics.setup import get_backend
    from .thermodynamics.types import ComplexSpec, Conditions, Strand
    from .oxdna_guided import seed_strand, _configuration
    from .preview import parse_structure
    result = copy.deepcopy(request)
    output = Path(output).resolve()
    preparation = output / 'preparation'
    # A paused folding/geometry preparation has no native state yet. Preserve
    # its artifacts and retry with the same frozen request in a fresh folder.
    attempt = 1
    while preparation.exists():
        attempt += 1
        preparation = output / f'preparation_{attempt:03d}'
    preparation.mkdir(parents=True)
    def check():
        if cancel_check and cancel_check():
            raise RuntimeError('Préparation MD mise en pause avant le lancement des répliques.')
    check()
    _validate_guide_request(result)
    backend = backend or get_backend(result.get('thermodynamics_library_path'))
    if asdict(backend.provenance) != result['guide_thermodynamics']:
        raise ValueError('Le moteur thermodynamique a changé depuis la préparation MD ; créez un nouveau calcul.')
    conditions = Conditions(result['sampling']['temperature_c'], result['sampling']['salt_m'], 0.)
    parameter = Path(result['runtime']['seq_dep_file_DNA']).expanduser().resolve()
    expected = result['runtime'].get('parameter_sha256', {}).get('DNA')
    if expected and _sha(parameter) != expected:
        raise ValueError('Les paramètres oxDNA ont changé depuis leur détection.')
    frozen_parameter = preparation / 'seq_dep_file_DNA.txt'
    shutil.copyfile(parameter, frozen_parameter)
    result['runtime']['seq_dep_file_DNA'] = str(frozen_parameter)
    seeds = []
    result['guide_provenance'] = asdict(backend.provenance)
    for variant in result['variants']:
        check()
        if progress_callback:
            progress_callback('Moteur thermodynamique intégré · préparation du complexe ' + variant['name'] + '…')
        aptamer, trigger = variant['aptamer_sequence'], variant['trigger_sequence']
        spec = ComplexSpec((Strand(aptamer, 'dna', 'aptamer'), Strand(trigger, 'dna', 'trigger')))
        minimum = backend.minimum(spec, conditions)
        check()
        if (minimum.context.complex_spec != spec or minimum.context.conditions != conditions
                or minimum.context.provenance != backend.provenance or not math.isfinite(minimum.energy_kcal_mol)):
            raise ValueError('Le moteur thermodynamique n’a pas renvoyé de complexe ADN cohérent.')
        guide = minimum.structure
        segments = guide.split('+')
        if len(segments) != 2 or [len(part) for part in segments] != [len(aptamer), len(trigger)]:
            raise ValueError('Le guide thermodynamique ne correspond pas aux deux brins ADN.')
        # A temporary unpaired linker makes a nested complex representable by
        # the geometry builder. Its beads are removed BEFORE creating topology;
        # no linker, covalent bridge or ligand is present in any simulation.
        combined = seed_strand(aptamer + 'T' * 8 + trigger, segments[0] + '.' * 8 + segments[1], 'dna')
        seed = np.concatenate([combined[:len(aptamer)], combined[len(aptamer) + 8:]])
        pairs, count, cuts = parse_structure(guide)
        if count != len(aptamer + trigger) or cuts != {len(aptamer) - 1}:
            raise ValueError('Le guide thermodynamique ne correspond pas aux deux brins ADN.')
        pairs = [[int(i), int(j)] for i, j in pairs]
        if not any(i < len(aptamer) <= j for i, j in pairs):
            raise ValueError('Aucun appariement aptamère–trigger dans le guide du groupe ' + variant['name'] + '.')
        variant.update(guide=guide, guide_pairs=pairs, mfe_energy_kcal_mol=minimum.energy_kcal_mol,
                       guide_solver=minimum.context.solver)
        seeds.append(seed)
    check()
    box = max(result['box_policy']['box_size'],
              math.ceil(max(float(np.ptp(seed[:, :3], axis=0).max()) for seed in seeds)
                        + 2 * result['box_policy']['minimum_periodic_clearance']))
    result['box_policy'].update(box_size=box, box_size_nm=box * .8518, common_to_both_groups=True)
    # A descriptive number density only; no association rate is inferred here.
    result['box_policy']['concentration_each_m'] = 1 / (6.02214076e23 * (box * .8518e-9) ** 3 * 1000)
    hashes = {}
    for variant, seed in zip(result['variants'], seeds):
        directory = preparation / variant['name']
        directory.mkdir()
        seed[:, :3] += box / 2 - (seed[:, :3].min(axis=0) + seed[:, :3].max(axis=0)) / 2
        config, topology = directory / 'seed.dat', directory / 'system.top'
        _configuration(config, seed, box)
        aptamer, trigger = variant['aptamer_sequence'], variant['trigger_sequence']
        topology.write_text(f'{len(aptamer) + len(trigger)} 2 5->3\n{aptamer} type=DNA\n{trigger} type=DNA\n', encoding='utf-8')
        variant.update(configuration_file=str(config), topology_file=str(topology), box_size=box)
        for path in (config, topology):
            hashes[path.relative_to(output).as_posix()] = _sha(path)
    hashes[frozen_parameter.relative_to(output).as_posix()] = _sha(frozen_parameter)
    result['prepared_file_sha256'] = hashes
    result['guide_source'] = 'native_hybrid_mfe_two_strand_dna_complex'
    result['preparation_geometry'] = 'temporary_unpaired_geometry_linker_removed_before_native_simulation'
    _write(output / 'resolved_request.json', result)
    return result


def free_md_worker_running(output_dir):
    path = Path(output_dir) / 'worker.lock'
    if not path.is_file():
        return False
    try:
        import fcntl
        with path.open('r') as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return True
            fcntl.flock(handle, fcntl.LOCK_UN)
        return False
    except ImportError:
        import msvcrt
        try:
            with path.open('r+b') as handle:
                handle.seek(0)
                try:
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                except OSError:
                    return True
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            return False
        except OSError:
            return False
    except OSError:
        return False


def load_free_md_resume_request(output_dir):
    output = Path(output_dir).resolve()
    if free_md_worker_running(output):
        raise ValueError('Cette simulation est déjà en cours.')
    if not (output / 'resume_state.json').is_file():
        raise ValueError('Aucune sauvegarde de MD libre ne permet de reprendre ce calcul.')
    resume = json.loads((output / 'resume_state.json').read_text())
    if resume.get('stage') == 'preparation':
        path = output / 'preparation_request.json'
        if path.is_symlink() or _sha(path) != resume.get('request_sha256'):
            raise ValueError('Les paramètres de préparation ont changé depuis leur sauvegarde.')
        if any(output.glob('replica_*/state.json')):
            raise ValueError('Une préparation incomplète ne peut pas remplacer une trajectoire native.')
        request = json.loads(path.read_text())
        _validate_guide_request(request)
        return {**request, 'resume': True, '_resume_preparation': True}
    request = json.loads((output / 'worker_request.json').read_text())
    if request.get('format') != 'aptaswitch-oxdna-free-md' or request.get('version') != 1:
        raise ValueError('Format de reprise MD incompatible.')
    _validate_guide_request(request)
    hashes = request.get('prepared_file_sha256')
    if not hashes:
        raise ValueError('Les empreintes de la préparation MD sont absentes.')
    for filename, expected in hashes.items():
        path = output / filename
        if path.is_symlink() or not path.resolve().is_relative_to(output) or _sha(path) != expected:
            raise ValueError('Un fichier de préparation MD a changé ; reprise refusée.')
    paths = [Path(v[field]) for v in request['variants'] for field in ('configuration_file', 'topology_file')]
    paths.append(Path(request['runtime']['seq_dep_file_DNA']))
    if any(not path.resolve().is_relative_to(output) or path.is_symlink() or
           path.resolve().relative_to(output).as_posix() not in hashes for path in paths):
        raise ValueError('Un fichier de préparation est extérieur à ce calcul.')
    from .oxdna_free_md_worker import _identity, _validate
    _validate(request)
    if resume.get('identity') != _identity(request):
        raise ValueError('Les conditions ont changé depuis la sauvegarde de la MD libre.')
    return {**request, 'resume': True}


def run_free_md(request, output_dir, progress_callback=None, cancel_check=None):
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if free_md_worker_running(output):
        raise ValueError('Cette simulation est déjà en cours.')
    request, process = copy.deepcopy(request), None
    started, result, previous = time.monotonic(), {}, {}
    def progress(message):
        if progress_callback:
            progress_callback(0, request['replicas'], message)
    try:
        resuming = request.pop('resume', False)
        pending_preparation = not resuming
        if resuming:
            previous_path = output / 'result.json'
            if previous_path.is_file():
                previous = json.loads(previous_path.read_text())
            request = load_free_md_resume_request(output)
            request.pop('resume')
            pending_preparation = request.pop('_resume_preparation', False)
            progress('Reprise de la préparation avec les paramètres enregistrés…' if pending_preparation
                     else 'Reprise des configurations et états aléatoires sauvegardés…')
        else:
            if (output / 'worker_request.json').exists() or (output / 'resume_state.json').exists():
                raise ValueError('Ce calcul existe déjà ; utilisez la reprise.')
            _write(output / 'preparation_request.json', request)
            _write(output / 'resume_state.json', {'format': 'aptaswitch-oxdna-free-md-resume', 'version': 1,
                   'stage': 'preparation', 'request_sha256': _sha(output / 'preparation_request.json')})
        if pending_preparation:
            request = prepare_free_md_system(request, output, progress_callback=progress, cancel_check=cancel_check)
            _write(output / 'worker_request.json', request)
            from .oxdna_free_md_worker import _identity
            _write(output / 'resume_state.json', {'format': 'aptaswitch-oxdna-free-md-resume', 'version': 1,
                   'identity': _identity(request), 'snapshot_id': request['snapshot_id']})
        if cancel_check and cancel_check():
            raise RuntimeError('MD libre mise en pause avant le lancement des répliques.')
        (output / 'CANCEL').unlink(missing_ok=True)
        _write(output / 'worker_request.json', request)
        result_path = output / 'result.json'
        old_stamp = result_path.stat().st_mtime_ns if result_path.exists() else None
        worker = Path(__file__).with_name('oxdna_free_md_worker.py')
        with (output / 'engine.log').open('ab') as log:
            process = subprocess.Popen([request['runtime']['python'], str(worker), str(output / 'worker_request.json')],
                cwd=output, stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                env=external_process_environment(), **process_group_options())
            last_progress, stopping = '', None
            while process.poll() is None:
                if cancel_check and cancel_check():
                    if stopping is None:
                        (output / 'CANCEL').touch()
                        stopping = time.monotonic()
                    if time.monotonic() - stopping > 30:
                        stop_process(process)
                        break
                path = output / 'progress.json'
                try:
                    if path.is_file() and path.stat().st_size < 32 * 1024 * 1024:
                        text = path.read_text()
                        if text != last_progress:
                            data = json.loads(text)
                            if progress_callback:
                                progress_callback(data.get('production_steps', 0), max(1, data.get('target_steps', 1)),
                                    data.get('message') or f"MD libre · {data.get('completed_replicas', 0)}/{request['replicas']} répliques terminées")
                            last_progress = text
                except (OSError, ValueError, TypeError):
                    pass
                time.sleep(.1)
        if result_path.is_file() and result_path.stat().st_mtime_ns != old_stamp:
            result = json.loads(result_path.read_text())
        else:
            result = {'status': 'failed', 'error': 'Le moteur oxDNA s’est arrêté ; consultez engine.log.'}
        if process.returncode and result.get('status') == 'completed':
            result.update(status='failed', sampling_complete=False, error='Le moteur oxDNA s’est arrêté ; consultez engine.log.')
        if cancel_check and cancel_check():
            result.update(status='cancelled', sampling_complete=False)
    except Exception as exc:
        cancelled = bool(cancel_check and cancel_check())
        result = {**previous, 'status': 'cancelled' if cancelled else 'failed', 'sampling_complete': False}
        if cancelled:
            result.pop('error', None)
            result['interruption_reason'] = str(exc)
        else:
            result['error'] = str(exc)
    finally:
        if process is not None and process.poll() is None:
            stop_process(process)
    # A hard interruption may stop the supervisor before its final summary.
    # Committed per-replica checkpoints still define a useful partial result.
    if 'replicas' not in result and any(output.glob('replica_*/state.json')):
        try:
            from .oxdna_free_md_metrics import summarize_run
            result = {**summarize_run(output, request), **result}
        except (OSError, ValueError, KeyError):
            pass
    if result.get('status') != 'completed':
        result['sampling_complete'] = False
        if result.get('status') == 'failed' and result.get('failures') and not result.get('error'):
            details = ['Réplique ' + str(row['replica_index'] + 1) + ' : ' + str(row.get('error', 'arrêt du moteur'))
                       for row in result['failures'][:3]]
            result['error'] = 'Certaines répliques ont échoué. ' + ' ; '.join(details)
    elif not result.get('sampling_complete'):
        result.update(status='failed', error='Le moteur a terminé sans produire toutes les répliques demandées.')
    result.update(format='aptaswitch-oxdna-free-md-result', version=1, method='FREE_MD_DNA_DNA',
                  snapshot_id=request['snapshot_id'], candidate_identity=request['saved_switch_identity'],
                  request=request, model='DNA2', elapsed_seconds=time.monotonic() - started,
                  output_dir=str(output), ligand_simulated=False)
    _write(output / 'result.json', result)
    return result


def make_free_md_archive(run_dir):
    """Export only regular files physically contained in this calculation."""
    root = Path(run_dir).resolve()
    if not root.is_dir():
        raise ValueError('Dossier de calcul MD introuvable.')
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(root.rglob('*')):
            if (path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root)
                    or path.name in {'worker.lock', 'CANCEL'} or path.name.endswith('.tmp')):
                continue
            archive.write(path, path.relative_to(root).as_posix())
    return buffer.getvalue()
