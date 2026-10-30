"""External, resumable DNA2 trajectories of saved aptamer–trigger complexes.

Every replica owns an oxpy process. Preparation restraints are removed before
equilibration and production. Native binary checkpoints preserve the oxDNA RNG
and phase point; the human-readable .dat files retain positions and momenta.
"""
from __future__ import annotations

from contextlib import contextmanager
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np


class WorkerBusyError(RuntimeError):
    """Another process owns the output; never overwrite its state or result."""


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', encoding='utf-8') as handle:
        json.dump(value, handle, ensure_ascii=False, allow_nan=False, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except FileNotFoundError:
        return default


def _identity(request: dict) -> str:
    scientific = {key: value for key, value in request.items() if key not in {'workers'}}
    content = json.dumps(scientific, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    digest = hashlib.sha256(content)
    for variant in request['variants']:
        for field in ('configuration_file', 'topology_file'):
            digest.update(Path(variant[field]).read_bytes())
    digest.update(Path(request['runtime']['seq_dep_file_DNA']).read_bytes())
    return digest.hexdigest()


def _validate(request: dict) -> None:
    if request.get('format') != 'aptaswitch-oxdna-free-md' or request.get('version') != 1:
        raise ValueError('Unsupported free oxDNA MD request.')
    replicas = request.get('replicas')
    if type(replicas) is not int or replicas < 2 or replicas % 2:
        raise ValueError('An even number of replicas is required for a balanced comparison.')
    if type(request.get('workers')) is not int or not 1 <= request['workers'] <= replicas:
        raise ValueError('Invalid number of free MD workers.')
    variants = request.get('variants', [])
    if len(variants) != 2 or {v.get('name') for v in variants} != {'original', 'extended'}:
        raise ValueError('The original and extended trigger are both required.')
    sampling = request['sampling']
    for key in ('production_steps', 'sample_every'):
        if type(sampling.get(key)) is not int or sampling[key] < 1:
            raise ValueError(f'Invalid free MD {key}.')
    if type(sampling.get('equilibration_steps')) is not int or sampling['equilibration_steps'] < 0:
        raise ValueError('Invalid free MD equilibration_steps.')
    for key in ('dt', 'salt_m'):
        if not math.isfinite(float(sampling[key])) or float(sampling[key]) <= 0:
            raise ValueError(f'Invalid free MD {key}.')
    if not math.isfinite(float(sampling['temperature_c'])) or sampling['temperature_c'] <= -273.15:
        raise ValueError('Invalid free MD temperature.')
    if type(sampling.get('seed')) is not int or not 0 <= sampling['seed'] < 2 ** 31:
        raise ValueError('Invalid free MD seed.')
    for variant in variants:
        aptamer, trigger = variant['aptamer_sequence'], variant['trigger_sequence']
        if not aptamer or not trigger or set(aptamer + trigger) - set('ACGT'):
            raise ValueError('Free DNA2 MD requires two nonempty DNA sequences.')
        if any(type(i) is not int or not 1 <= i <= len(aptamer)
               for i in request.get('monitored_bases', [])) or not request.get('monitored_bases'):
            raise ValueError('Monitored bases must be 1-based aptamer positions.')
        for pair in variant['guide_pairs']:
            if len(pair) != 2 or any(type(i) is not int or not 0 <= i < len(aptamer + trigger) for i in pair) or pair[0] == pair[1]:
                raise ValueError('Invalid complex guide pair.')


def parse_frame(frame: str) -> tuple[np.ndarray, np.ndarray]:
    lines = [line for line in frame.splitlines() if line.strip()]
    box = np.asarray([float(value) for value in lines[1].split('=', 1)[1].split()])
    rows = np.asarray([[float(value) for value in line.split()] for line in lines[3:]])
    if rows.ndim != 2 or rows.shape[1] != 15 or box.shape != (3,) or not np.isfinite(rows).all() or not np.isfinite(box).all() or np.any(box <= 0):
        raise ValueError('Invalid native oxDNA checkpoint.')
    return rows, box


def guide_evidence(variant: dict, bonds) -> dict:
    expected = {tuple(sorted(pair)) for pair in variant['guide_pairs']}
    actual = {tuple(sorted(pair)) for pair in bonds}
    boundary = len(variant['aptamer_sequence'])
    return {'expected_pairs': len(expected), 'native_guide_pairs': len(expected & actual),
            'missing_pairs': [list(pair) for pair in sorted(expected - actual)],
            'additional_pairs': [list(pair) for pair in sorted(actual - expected)],
            'hydrogen_bonds': [list(pair) for pair in sorted(actual)],
            'interstrand_pairs': sum(i < boundary <= j for i, j in actual),
            'equilibrium_validated': False}


class DNAAdapter:
    """One native oxDNA manager; binary reload restores the native RNG."""

    def __init__(self, request, variant, configuration, *, directory, seed, phase,
                 initial=False, binary=False, forces=None):
        import oxpy
        self.context = oxpy.Context()
        self.context.__enter__()
        self.manager = self.info = self.conf = self.hb = None
        self.directory = Path(directory)
        self.binary_file = self.directory / 'native_current.bin'
        sampling = request['sampling']
        relaxed = phase in {'minimize_guided', 'relax_guided'}
        values = {
            'backend': 'CPU', 'sim_type': 'min' if phase == 'minimize_guided' else 'MD',
            'interaction_type': 'DNA_relax' if relaxed else 'DNA2',
            'topology': str(Path(variant['topology_file']).resolve()),
            'T': f"{sampling['temperature_c']}C", 'salt_concentration': str(sampling['salt_m']),
            'dt': str(sampling['dt'] if phase == 'production' else .001),
            'thermostat': 'brownian', 'newtonian_steps': '103', 'diff_coeff': '2.5',
            'refresh_vel': 'true' if initial else 'false', 'restart_step_counter': 'false',
            'time_scale': 'linear', 'verlet_skin': '.3', 'steps': '0',
            'print_conf_interval': '1000000000', 'print_energy_every': '1000000000',
            'trajectory_file': os.devnull, 'energy_file': os.devnull,
            'lastconf_file': os.devnull, 'lastconf_file_bin': str(self.binary_file),
            'log_file': str(self.directory / ('native_' + phase + '.log')),
            'use_average_seq': 'false', 'seq_dep_file': str(Path(request['runtime']['seq_dep_file_DNA']).resolve()),
            'external_forces': 'true' if forces else 'false',
            # Avoid manager-local periodic centre-of-mass shifts on restart.
            'fix_diffusion': 'false', 'reset_initial_com_momentum': 'false',
            'data_output_1': '{\nname = ' + os.devnull + '\nprint_every = 1000000000\nonly_last = true\ncol_1 = {\ntype = Configuration\nid = free_md_conf\n}\ncol_2 = {\ntype = hb_list\nid = free_md_hb\n}\n}',
        }
        if relaxed:
            # DNA_relax inherits the DNA1 default (equal grooves). Keep the
            # DNA2 backbone interaction-site geometry throughout preparation,
            # otherwise restoring DNA2 can turn relaxed tails into broken FENE
            # bonds despite unchanged nucleotide positions and orientations.
            values.update(relax_type='harmonic_force', relax_strength='32',
                          major_minor_grooving='true')
        if forces:
            values['external_forces_file'] = str(forces)
        if binary:
            values['reload_from'] = str(Path(configuration).resolve())
        else:
            values.update(conf_file=str(Path(configuration).resolve()), seed=str(seed))
        (self.directory / ('native_' + phase + '_input.txt')).write_text(
            '\n'.join(f'{key} = {value}' for key, value in values.items()) + '\n', encoding='utf-8')
        try:
            inp = oxpy.InputFile()
            for key, value in values.items():
                inp[key] = value
            self.manager = oxpy.OxpyManager(inp)
            self.info = self.manager.config_info()
            self.conf = self.info.get_observable_by_id('free_md_conf')
            self.hb = self.info.get_observable_by_id('free_md_hb')
        except BaseException:
            self.close()
            raise

    def run(self, steps):
        before = self.manager.current_step
        self.manager.run(steps, print_output=False)
        if self.manager.current_step - before != steps:
            raise RuntimeError('Native free MD stopped before completing the requested chunk.')

    def frame(self):
        self.manager.update_observable_data()
        return self.conf.get_output_string(self.manager.current_step).rstrip() + '\n'

    def bonds(self):
        self.info.interaction.begin_energy_computation()
        return {tuple(sorted(pair)) for pair in self.hb.hb_list()}

    def energy(self):
        result = float(self.manager.system_energy())
        if not math.isfinite(result):
            raise ValueError('Nonfinite native free MD energy.')
        return result

    def checkpoint(self, target):
        self.manager.update_observable_data()
        self.manager.print_configuration(also_last=True)
        if not self.binary_file.is_file() or not self.binary_file.stat().st_size:
            raise RuntimeError('Native oxDNA did not write the binary RNG checkpoint.')
        shutil.copyfile(self.binary_file, target)

    def close(self):
        self.conf = self.hb = self.info = self.manager = None
        gc.collect()
        context, self.context = self.context, None
        if context:
            context.__exit__(None, None, None)


class Replica:
    """Transactional trajectory/sample appends and a resumable native phase."""

    def __init__(self, request, output, index, *, adapter_factory=DNAAdapter):
        _validate(request)
        if type(index) is not int or not 0 <= index < request['replicas']:
            raise ValueError('Invalid replica index.')
        self.request, self.root, self.index = request, Path(output).resolve(), index
        self.output = self.root / f'replica_{index:03d}'
        self.output.mkdir(parents=True, exist_ok=True)
        # Interleaving keeps a paused partial comparison balanced in launch order.
        name = 'original' if index % 2 == 0 else 'extended'
        self.variant = next(item for item in request['variants'] if item['name'] == name)
        self.adapter_factory = adapter_factory
        self.identity = _identity(request)
        self.state_path = self.output / 'state.json'
        self.state = read_json(self.state_path)
        if self.state:
            if self.state.get('identity') != self.identity or self.state.get('replica_index') != index or self.state.get('variant') != name:
                raise ValueError('The saved free MD checkpoint belongs to different conditions.')
            self._rollback_appends()
        else:
            seed = (request['sampling']['seed'] + 104729 * (index + 1)) % (2 ** 31 - 1)
            self.state = {'version': 1, 'identity': self.identity, 'replica_index': index,
                          'variant': name, 'phase': 'minimize_guided', 'phase_steps': 0,
                          'production_steps': 0, 'samples_count': 0, 'status': 'running',
                          'seed': seed, 'generation': 0, 'trajectory_bytes': 0, 'samples_bytes': 0,
                          'topology_file': self.variant['topology_file'], 'preparation': [],
                          'rng_provenance': 'independent replica seed; native binary checkpoint/reload_from',
                          'equilibrium_validated': False, 'production_external_forces': False}
            if any((self.output / name).exists() for name in ('trajectory.dat', 'samples.jsonl')):
                raise ValueError('Free MD output exists without a committed state; refusing to overwrite it.')
            write_json(self.state_path, self.state)
        self.forces = self.output / 'guide_forces.txt'
        self.forces.write_text('\n'.join(
            '{\ntype = mutual_trap\nparticle = %d\nref_particle = %d\nr0 = 1.2\nstiff = 0.9\nPBC = true\n}\n' % (i, j)
            for pair in self.variant['guide_pairs'] for i, j in (pair, pair[::-1])), encoding='utf-8')
        self.last_save = time.monotonic()

    def _rollback_appends(self):
        # A process may die after appending data but before committing its state.
        # Restore both streams to the same committed sample, including partial lines.
        for name, key in (('trajectory.dat', 'trajectory_bytes'), ('samples.jsonl', 'samples_bytes')):
            path = self.output / name
            expected = self.state.get(key, 0)
            if not path.exists():
                if expected:
                    raise ValueError(f'Missing committed free MD data: {name}.')
                continue
            if path.stat().st_size < expected:
                raise ValueError(f'Truncated committed free MD data: {name}.')
            with path.open('r+b') as handle:
                handle.truncate(expected)

    def _commit(self, adapter):
        generation = self.state['generation'] + 1
        slot = generation % 2
        dat = self.output / f'checkpoint_{slot}.dat'
        binary = self.output / f'checkpoint_{slot}.bin'
        frame = adapter.frame()
        rows, box = parse_frame(frame)
        self._geometry_check(rows, box)
        dat.write_text(frame, encoding='utf-8')
        adapter.checkpoint(binary)
        self.state.update(generation=generation, checkpoint_file=dat.name,
                          checkpoint_binary=binary.name, native_step=int(adapter.manager.current_step))
        write_json(self.state_path, self.state)
        # Convenience files are never authoritative for recovery.
        shutil.copyfile(dat, self.output / 'checkpoint.dat')
        self.last_save = time.monotonic()

    def _geometry_check(self, rows, box):
        n = len(self.variant['aptamer_sequence'])
        if len(rows) != n + len(self.variant['trigger_sequence']):
            raise ValueError('The native free MD configuration has an unexpected nucleotide count.')
        if not np.allclose(box, self.variant['box_size'], rtol=0, atol=1e-8):
            raise ValueError('The native free MD box changed from the fixed common simulation volume.')
        policy = self.request.get('box_policy')
        if policy and policy.get('geometry_monitor_required'):
            try:
                from .oxdna_box import check_periodic_self_clearance
            except ImportError:
                from oxdna_box import check_periodic_self_clearance
            try:
                evidence = check_periodic_self_clearance(rows[:, :3].tolist(), [n, len(rows) - n],
                    box.tolist(), policy['interaction_cutoff'], policy.get('periodic_margin', 1.))
            except ValueError as exc:
                raise ValueError(str(exc).replace('This configuration must not be retained for FFS.',
                    'Free MD stopped before retaining this configuration.')) from exc
            self.state['periodic_self_check'] = evidence

    def _sample(self, adapter):
        try:
            from .oxdna_free_md_metrics import frame_metrics
        except ImportError:
            from oxdna_free_md_metrics import frame_metrics
        frame = adapter.frame()
        rows, box = parse_frame(frame)
        self._geometry_check(rows, box)
        reference_path = self.output / 'reference.dat'
        if not reference_path.exists():
            reference_path.write_text(frame, encoding='utf-8')
        reference, _ = parse_frame(reference_path.read_text(encoding='utf-8'))
        n = len(self.variant['aptamer_sequence'])
        bonds = adapter.bonds()
        data = frame_metrics(rows, box, bonds, n, self.request['monitored_bases'],
                             reference[:n, :3], float(self.request.get('contact_cutoff_nm', 1.5)))
        step, dt = self.state['production_steps'], self.request['sampling']['dt']
        data.update(variant=self.variant['name'], replica_index=self.index, step=step,
                    native_step=int(adapter.manager.current_step), time_model=step * dt,
                    time_us=step * dt * 3.031e-6, energy=adapter.energy(),
                    hydrogen_bonds=[list(pair) for pair in sorted(bonds)])
        encoded = (json.dumps(data, allow_nan=False, separators=(',', ':')) + '\n').encode()
        for filename, payload, key in (('trajectory.dat', frame.encode(), 'trajectory_bytes'),
                                       ('samples.jsonl', encoded, 'samples_bytes')):
            with (self.output / filename).open('ab') as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
                self.state[key] = handle.tell()
        self.state['samples_count'] += 1
        self._commit(adapter)

    def run(self):
        if self.state['status'] == 'completed':
            return self.state
        phases = [('minimize_guided', 2000, True), ('relax_guided', 20000, True),
                  ('guide_physical', 30000, True), ('release', 10000, False),
                  ('equilibration', self.request['sampling']['equilibration_steps'], False),
                  ('production', self.request['sampling']['production_steps'], False)]
        while self.state['phase'] != 'completed':
            if (self.root / 'CANCEL').exists():
                self.state['status'] = 'cancelled'
                write_json(self.state_path, self.state)
                return self.state
            phase_index = next(i for i, item in enumerate(phases) if item[0] == self.state['phase'])
            phase, total, restrained = phases[phase_index]
            binary = self.state.get('checkpoint_binary')
            configuration = self.output / binary if binary else Path(self.variant['configuration_file'])
            # The minimizer does not thermalize momenta. Refresh once when MD
            # starts; subsequent phase transitions and resumes preserve them.
            initial = (not binary or
                       (phase == 'relax_guided' and self.state['phase_steps'] == 0 and not self.state.get('md_initialized')))
            adapter = self.adapter_factory(self.request, self.variant, configuration,
                directory=self.output, seed=self.state['seed'], phase=phase, initial=initial,
                binary=bool(binary), forces=self.forces if restrained and self.variant['guide_pairs'] else None)
            self.state['status'] = 'running'
            try:
                if initial and phase == 'relax_guided':
                    self.state['md_initialized'] = True
                if phase == 'production' and self.state['samples_count'] == 0:
                    self._sample(adapter)
                while self.state['phase_steps'] < total:
                    if (self.root / 'CANCEL').exists():
                        self.state['status'] = 'cancelled'
                        self._commit(adapter)
                        return self.state
                    chunk = min(1000, total - self.state['phase_steps'])
                    if phase == 'production':
                        interval = self.request['sampling']['sample_every']
                        chunk = min(chunk, interval - self.state['production_steps'] % interval)
                    adapter.run(chunk)
                    self.state['phase_steps'] += chunk
                    if phase == 'production':
                        self.state['production_steps'] += chunk
                        if self.state['production_steps'] % interval == 0 or self.state['production_steps'] == total:
                            self._sample(adapter)
                    if time.monotonic() - self.last_save >= 2:
                        self._commit(adapter)
                if phase != 'production':
                    self.state['preparation'].append({'phase': phase, 'steps': total,
                        'model': 'DNA_relax' if phase in {'minimize_guided', 'relax_guided'} else 'DNA2',
                        'external_forces': bool(restrained and self.variant['guide_pairs']),
                        'validation': guide_evidence(self.variant, adapter.bonds())})
                self.state['phase'] = phases[phase_index + 1][0] if phase_index + 1 < len(phases) else 'completed'
                self.state['phase_steps'] = 0
                if self.state['phase'] == 'completed':
                    self.state['status'] = 'completed'
                self._commit(adapter)
            finally:
                adapter.close()
        return self.state


@contextmanager
def worker_lock(output):
    """One supervisor per run; an abandoned Unix lock releases with the PID."""
    path = Path(output) / 'worker.lock'
    if os.name != 'nt':
        import fcntl
        with path.open('a+') as handle:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise WorkerBusyError('A free MD worker already owns this calculation.') from exc
            handle.seek(0)
            handle.truncate()
            handle.write(str(os.getpid()))
            handle.flush()
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    else:
        # msvcrt locks are released by Windows on process death too.
        import msvcrt
        with path.open('a+b') as handle:
            if path.stat().st_size == 0:
                handle.write(b'0')
                handle.flush()
            handle.seek(0)
            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as exc:
                raise WorkerBusyError('A free MD worker already owns this calculation.') from exc
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)


def run_worker(request, output):
    """Supervise isolated replicas, with a fixed cap on concurrent native MD."""
    _validate(request)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with worker_lock(output):
        started = time.monotonic()
        identity = _identity(request)
        manifest_path = output / 'resume_state.json'
        manifest = read_json(manifest_path)
        if manifest and manifest.get('identity') != identity:
            raise ValueError('The saved free MD calculation belongs to different conditions.')
        request_path = output / 'worker_request.json'
        if request_path.exists() and _identity(read_json(request_path)) != identity:
            raise ValueError('A different free MD request already occupies this output directory.')
        write_json(request_path, request)
        write_json(manifest_path, {'format': 'aptaswitch-oxdna-free-md-resume', 'version': 1,
                                  'identity': identity, 'snapshot_id': request.get('snapshot_id')})
        jobs, pending = {}, list(range(request['replicas']))
        failures = []
        interrupted = False
        try:
            while pending or jobs:
                cancelled = (output / 'CANCEL').exists()
                if cancelled:
                    pending.clear()
                while pending and len(jobs) < request['workers'] and not cancelled:
                    index = pending.pop(0)
                    previous = read_json(output / f'replica_{index:03d}' / 'state.json', {})
                    if previous.get('status') == 'completed':
                        if previous.get('identity') != _identity(request):
                            raise ValueError('A completed replica belongs to another free MD request.')
                        continue
                    log = (output / f'replica_{index:03d}.log').open('ab')
                    environment = dict(os.environ)
                    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
                        environment[key] = '1'
                    process = subprocess.Popen([request['runtime']['python'], str(Path(__file__).resolve()),
                        str(request_path), '--replica', str(index)], stdout=log, stderr=subprocess.STDOUT,
                        cwd=output, env=environment)
                    jobs[index] = (process, log)
                for index, (process, log) in list(jobs.items()):
                    code = process.poll()
                    if code is not None:
                        log.close()
                        del jobs[index]
                        if code:
                            state = read_json(output / f'replica_{index:03d}' / 'state.json', {})
                            failures.append({'replica_index': index, 'returncode': code,
                                             'error': state.get('error', f'Native worker exited with code {code}.')})
                states = [read_json(output / f'replica_{index:03d}' / 'state.json', {}) for index in range(request['replicas'])]
                completed = sum(state.get('status') == 'completed' for state in states)
                total_steps = sum(state.get('production_steps', 0) for state in states)
                progress = {'status': 'cancelling' if cancelled and jobs else ('cancelled' if cancelled else 'running'),
                    'completed_replicas': completed, 'replicas': request['replicas'], 'active_replicas': sorted(jobs),
                    'current': completed, 'total': request['replicas'],
                    'message': f'oxDNA DNA2 · {completed}/{request["replicas"]} réplicats terminés',
                    'elapsed_seconds': time.monotonic() - started,
                    'production_steps': total_steps, 'target_steps': request['sampling']['production_steps'] * request['replicas'],
                    'fraction': total_steps / (request['sampling']['production_steps'] * request['replicas']),
                    'replica_states': [{key: state.get(key) for key in ('replica_index', 'variant', 'phase', 'status', 'production_steps', 'samples_count')} for state in states],
                    'failures': failures}
                write_json(output / 'progress.json', progress)
                if jobs:
                    time.sleep(.25)
        except BaseException:
            interrupted = True
            (output / 'CANCEL').touch()
            raise
        finally:
            if interrupted:
                for process, log in jobs.values():
                    try:
                        process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        process.terminate()
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait()
                    log.close()
        try:
            from .oxdna_free_md_metrics import summarize_run
        except ImportError:
            from oxdna_free_md_metrics import summarize_run
        report = summarize_run(output, request)
        status = 'failed' if failures else ('cancelled' if (output / 'CANCEL').exists() else 'completed')
        report.update(status=status, failures=failures)
        write_json(output / 'result.json', report)
        progress['status'] = status
        write_json(output / 'progress.json', progress)
        return report


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('request')
    parser.add_argument('--replica', type=int)
    args = parser.parse_args()
    path = Path(args.request).resolve()
    request = read_json(path)
    try:
        if args.replica is not None:
            directory = path.parent / f'replica_{args.replica:03d}'
            directory.mkdir(parents=True, exist_ok=True)
            with worker_lock(directory):
                result = Replica(request, path.parent, args.replica).run()
        else:
            result = run_worker(request, path.parent)
        return 1 if result.get('status') == 'failed' else 0
    except WorkerBusyError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except Exception as exc:
        if args.replica is not None:
            state_path = path.parent / f'replica_{args.replica:03d}' / 'state.json'
            state_path.parent.mkdir(parents=True, exist_ok=True)
            state = read_json(state_path, {})
            state.update(status='failed', error=str(exc))
            write_json(state_path, state)
        else:
            write_json(path.parent / 'result.json', {'status': 'failed', 'error': str(exc)})
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
