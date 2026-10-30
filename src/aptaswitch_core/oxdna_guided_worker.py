"""External oxpy preparation. Run only with the native oxDNA Python."""
from __future__ import annotations
import gc
import json
import math
import os
from pathlib import Path
import sys
import time

try:
    from .oxdna_geometry import _place_strands, configuration_geometry
    from .oxdna_loop_geometry import configuration_loop_threading
    from .oxdna_box import check_periodic_self_clearance
    from .oxdna_ffs_worker import distance_pair_groups
except ImportError:  # Direct invocation in the external oxpy interpreter.
    from oxdna_geometry import _place_strands, configuration_geometry
    from oxdna_loop_geometry import configuration_loop_threading
    from oxdna_box import check_periodic_self_clearance
    from oxdna_ffs_worker import distance_pair_groups


def write(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def validation_from_bonds(request, bonds, distance):
    if distance is not None and not math.isfinite(distance):
        raise ValueError('Native oxDNA preparation produced a non-finite strand distance; this checkpoint is rejected.')
    guides = {}
    expected = set()
    for role, guide in request['guides'].items():
        pairs = {tuple(pair) for pair in guide['pairs']}
        expected |= pairs
        guides[role] = {'expected_pairs':len(pairs),'native_pairs':len(pairs & bonds),
                        'missing_pairs':[list(pair) for pair in sorted(pairs-bonds)],
                        'fraction': len(pairs & bonds)/len(pairs) if pairs else 1.}
    boundary = len(request['switch_sequence'])
    stem = {tuple(pair) for pair in request.get('plan',{}).get('stem_pairs',[])}
    native_stem = len(stem & bonds)
    minimum = request.get('plan',{}).get('basin_min_stem',0)
    associations = [list(pair) for pair in sorted(bonds) if pair[0]<boundary<=pair[1]]
    separated = distance is None or distance >= request.get('plan', {}).get('basin_distance', 4.)
    guide_accepted = bonds == expected
    off_basin_eligible = native_stem>=minimum and not associations and separated
    return {'guides':guides,'hydrogen_bonds':[list(pair) for pair in sorted(bonds)],
            'unexpected_pairs':[list(pair) for pair in sorted(bonds-expected)],
            'interstrand_bonds':associations,'off_stem_pairs':native_stem,'off_stem_required':minimum,
            'guide_exact_match': bonds == expected,
            'basin_distance':distance,
            'off_basin_eligible':off_basin_eligible, 'guide_accepted':guide_accepted,
            'ffs_eligible':off_basin_eligible and guide_accepted}


class Stage:
    def __init__(self, request, configuration, *, kind='MD', relaxed=False, forces=None, dt=.001, seed=1):
        import oxpy
        self.context = oxpy.Context()
        self.context.__enter__()
        self.manager = self.info = self.conf = self.hb = None
        values = {
            'backend':'CPU', 'interaction_type':'NA_relax' if relaxed else 'NA',
            'sim_type':kind, 'topology':request['topology_file'], 'conf_file':str(configuration),
            'T':f"{request['temperature_c']}C", 'salt_concentration':str(request['salt_m']),
            'seed':str(seed), 'steps':'0', 'dt':str(dt), 'thermostat':'brownian',
            'newtonian_steps':'103', 'diff_coeff':'2.5', 'refresh_vel':'true',
            'restart_step_counter':'true', 'time_scale':'linear', 'verlet_skin':'.3',
            'print_conf_interval':'1000000000','print_energy_every':'1000000000',
            'trajectory_file':os.devnull,'energy_file':os.devnull,'lastconf_file':os.devnull,'log_file':os.devnull,
            'use_average_seq':'false', 'external_forces':'true' if forces else 'false',
            'data_output_1':'{\nname = '+os.devnull+'\nprint_every = 1000000000\nonly_last = true\ncol_1 = {\ntype = Configuration\nid = guided_conf\n}\ncol_2 = {\ntype = hb_list\nid = guided_hb\n}\n}',
        }
        if relaxed:
            values.update(relax_type='harmonic_force',relax_strength='32')
        if forces:
            values['external_forces_file'] = str(forces)
        for chemistry in ('DNA','RNA','NA'):
            key = 'seq_dep_file_'+chemistry
            values[key] = request['runtime'][key]
        try:
            settings = oxpy.InputFile()
            for key,value in values.items():
                settings[key]=value
            self.manager = oxpy.OxpyManager(settings)
            self.info = self.manager.config_info()
            self.conf = self.info.get_observable_by_id('guided_conf')
            self.hb = self.info.get_observable_by_id('guided_hb')
        except BaseException:
            self.close()
            raise
    def run(self, steps, deadline):
        remaining = steps
        while remaining:
            if deadline is not None and time.monotonic() >= deadline:
                raise TimeoutError('oxDNA guided preparation reached its time budget.')
            chunk = min(100, remaining)
            self.manager.run(chunk, print_output=False)
            remaining -= chunk
    def frame(self):
        self.manager.update_observable_data()
        frame = self.conf.get_output_string(self.manager.current_step).rstrip()+'\n'
        if any(not math.isfinite(float(value)) for line in frame.splitlines()[3:] for value in line.split()):
            raise ValueError('Native oxDNA preparation produced non-finite coordinates; this checkpoint is rejected.')
        return frame
    def observe(self,request):
        self.info.interaction.begin_energy_computation()
        bonds = {tuple(sorted(pair)) for pair in self.hb.hb_list()}
        groups = distance_pair_groups(request.get('plan', {}))
        separation = groups.get('separation', groups['distance'])
        particles = self.info.particles()
        distance = min((math.sqrt(self.info.box.sqr_min_image_distance(particles[i], particles[j]))
                        for i, j in separation), default=None)
        observation = validation_from_bonds(request, bonds, distance)
        frame = self.frame()
        lengths = [len(request[role+'_sequence']) for role in ('switch', 'trigger')]
        materials = [request[role+'_material'] for role in ('switch', 'trigger')]
        observation['geometry'] = configuration_geometry(frame, lengths, materials)
        observation['loop_threading'] = configuration_loop_threading(frame, lengths, materials,
            [request[role+'_guide'] for role in ('switch', 'trigger')])
        observation['ffs_eligible'] &= observation['geometry']['passed'] and observation['loop_threading']['passed']
        return observation
    def close(self):
        self.conf=self.hb=self.info=self.manager=None
        gc.collect()
        context,self.context=self.context,None
        if context:
            context.__exit__(None,None,None)


def separate_monomers(source, target, switch_length, separation, *, box_policy=None):
    """Translate whole relaxed monomers into a separated starting basin.

    Unwrap each covalent chain first. This is preparation, not a trajectory
    segment: rotations preserve chirality and internal geometry.
    The subsequent native NA initialization remains authoritative.
    """
    lines = source.read_text().splitlines()
    box = [float(value) for value in lines[1].split('=', 1)[1].split()]
    rows = [[float(value) for value in line.split()] for line in lines[3:] if line.strip()]
    if (len(box) != 3 or any(not math.isfinite(value) or value <= 0 for value in box)
            or not 0 < switch_length < len(rows)
            or any(len(row) != 15 or not all(math.isfinite(value) for value in row) for row in rows)
            or not math.isfinite(separation) or separation <= 0):
        raise ValueError('Malformed relaxed oxDNA configuration or separation.')
    strands = [rows[:switch_length], rows[switch_length:]]
    for strand in strands:
        for previous, row in zip(strand, strand[1:]):
            for axis in range(3):
                delta = row[axis] - previous[axis]
                row[axis] = previous[axis] + (delta+box[axis]/2) % box[axis] - box[axis]/2
    import numpy as np
    if any(abs(side-box[0]) > 1e-8 for side in box):
        raise ValueError('Guided preparation requires a cubic box.')
    options = {'reset_velocities': False}
    if box_policy:
        options.update(periodic_margin=box_policy['minimum_periodic_clearance'],
                       auto_size=True, reorient=True)
    placed, side = _place_strands(np.asarray(strands[0]), np.asarray(strands[1]),
                                  separation, box[0], **options)
    header = [lines[0], f'b = {side:.17g} {side:.17g} {side:.17g}', lines[2]]
    target.write_text('\n'.join(header)+'\n'+
                      '\n'.join(' '.join(f'{value:.17g}' for value in row) for row in placed)+'\n')
    guard = (check_periodic_self_clearance(placed[:, :3].tolist(),
                 [switch_length, len(placed)-switch_length], [side]*3,
                 box_policy['interaction_cutoff']) if box_policy else None)
    return {'method':('rigid monomer rotation and translation after periodic unwrapping' if box_policy
                      else 'rigid monomer translation after periodic unwrapping'),
            'box_size': side, 'periodic_self_check': guard,
            'gap_between_com_bounding_boxes':separation,
            'minimum_periodic_x_gap':side-float(np.ptp(placed[:, 0])),
            'intrachain_geometry_changed':False,'production_external_forces':False}


def run(request, output):
    """Find a guide-exact *prepared* checkpoint, never claim equilibration.

    A restrained or relaxed observation can nominate a candidate, but only a
    separately initialized, strict NA trajectory without forces can accept it.
    Failed probes resume the saved preparation checkpoint with a fresh seed.
    """
    budget = request.get('budget_seconds')
    deadline = None if budget is None else time.monotonic()+float(budget)
    attempt = request.get('preparation_attempt', {})
    stiffness = float(attempt.get('restraint_stiffness', .9))
    if not math.isfinite(stiffness) or stiffness <= 0:
        raise ValueError('Preparation restraint stiffness must be finite and positive.')
    sequence_number = 0

    def next_seed():
        nonlocal sequence_number
        sequence_number += 1
        # Distinct streams for every segment and probe, reproducible per attempt.
        return (int(request['seed'])+104729*sequence_number) % (2**31-1) or 1

    def progress(message, phase):
        write(output/'guided_progress.json', {'message': message, 'phase': phase})

    pairs = [pair for guide in request['guides'].values() for pair in guide['pairs']]
    forces = output/'guide_forces.txt'
    forces.write_text('\n'.join(
        '{\ntype = mutual_trap\nparticle = %d\nref_particle = %d\nr0 = 1.2\nstiff = %.17g\nPBC = true\n}\n' % (i, j, stiffness)
        for pair in pairs for i, j in (pair, pair[::-1])))
    configuration = Path(request['configuration_file'])
    stages, probes = [], []
    best_probe = None
    best_probe_rank = None

    def record():
        write(output/'guided_stages.json', stages)
        write(output/'guided_probes.json', probes)

    def probe(source, name, step, *, purpose='candidate_validation'):
        """Called only with no other oxpy context alive."""
        nonlocal best_probe, best_probe_rank
        number = len(probes)+1
        entry = {'index': number, 'source_stage': name, 'source_step': step,
                 'seed': next_seed(), 'model': 'NA', 'restraints': False,
                 'planned_steps': 1000, 'sampling_steps': 100, 'steps': 0, 'accepted': False,
                 'geometry_checks': [], 'purpose': purpose}
        probes.append(entry)
        record()
        progress(f'oxDNA · strict NA check {number}', 'probe')
        stage = None
        try:
            separated = output/f'probe_{number:03d}_separated.dat'
            placement = separate_monomers(source, separated, len(request['switch_sequence']),
                float(request.get('initial_separation', 6.)), box_policy=request.get('box_policy'))
            entry['placement'] = placement
            stage = Stage(request, separated, seed=entry['seed'])
            entry['initial'] = stage.observe(request)
            if not entry['initial']['geometry']['passed']:
                raise ValueError('The prepared conformation has repulsive-core contacts or invalid backbone bonds.')
            if not entry['initial']['loop_threading']['passed']:
                raise ValueError('The prepared conformation threads a guide loop or its loop surface cannot be checked.')
            # A short unforced trajectory checks relaxation beyond one snapshot.
            # These checks apply to preparation only, never to FFS acceptance.
            while entry['steps'] < entry['planned_steps']:
                stage.run(entry['sampling_steps'], deadline)
                entry['steps'] += entry['sampling_steps']
                observation = stage.observe(request)
                entry['geometry_checks'].append({'step': entry['steps'], 'geometry': observation['geometry'],
                                                'loop_threading': observation['loop_threading']})
                if not observation['geometry']['passed']:
                    raise ValueError('The unforced preparation probe developed a severe geometric clash.')
                if not observation['loop_threading']['passed']:
                    raise ValueError('The unforced preparation probe threads a guide loop or its loop surface cannot be checked.')
            entry['final'] = observation
            path = output/f'probe_{number:03d}_validated.dat'
            path.write_text(stage.frame())
            entry['configuration_file'] = str(path)
            entry['accepted'] = observation['ffs_eligible']
            entry['strict_force_field_initialization_passed'] = True
            result = (path, observation, placement, {'stage': name, 'step': step, 'probe_index': number})
            missing = sum(value['expected_pairs']-value['native_pairs'] for value in observation['guides'].values())
            rank = (not observation['ffs_eligible'], missing+len(observation['unexpected_pairs']),
                    len(observation['interstrand_bonds']), missing, not observation['off_basin_eligible'])
            if best_probe_rank is None or rank < best_probe_rank:
                best_probe, best_probe_rank = result, rank
            return result if entry['accepted'] else None
        except (TimeoutError, ImportError, FileNotFoundError, PermissionError):
            raise
        except Exception as exc:
            # A stretched, clashing, or fraying checkpoint is a failed attempt
            # to validate, not permission to silently weaken the model.
            entry['error'] = str(exc)
            return None
        finally:
            if stage is not None:
                stage.close()
            record()

    selected = None
    settings = [('minimize_guided', 'min', True, 2000, True, 250),
                ('relax_guided', 'MD', True, 20000, True, 1000),
                ('guide_physical', 'MD', False, 30000, True, 1000),
                ('release', 'MD', False, 10000, False, 1000)]
    for name, kind, relaxed, planned_steps, use_forces, sampling in settings:
        progress('oxDNA · '+name, name)
        entry = {'name': name, 'planned_steps': planned_steps, 'steps': 0,
                 'sampling_steps': sampling, 'model': 'NA_relax' if relaxed else 'NA',
                 'restraints': use_forces and bool(pairs),
                 'restraint_stiffness': stiffness if use_forces and pairs else None,
                 'segments': [], 'checkpoints': [], 'selected_checkpoint_step': None}
        stages.append(entry)
        stage = None
        try:
            while entry['steps'] <= planned_steps:
                if stage is None:
                    seed = next_seed()
                    stage = Stage(request, configuration, kind=kind, relaxed=relaxed,
                        forces=forces if use_forces and pairs else None, seed=seed)
                    entry['segments'].append({'seed': seed, 'start_step': entry['steps'], 'steps': 0})
                observation = stage.observe(request)
                if 'initial' not in entry:
                    entry['initial'] = observation
                frame = stage.frame()
                entry['final'] = observation
                entry['checkpoints'].append({'step': entry['steps'], 'validation': observation})
                configuration = output/(name+'.dat')
                configuration.write_text(frame)
                record()
                # The OFF distance is only meaningful after rigidly separating
                # the monomers; do not reject useful internal conformations here.
                if (observation['guide_exact_match'] and observation['geometry']['passed']
                        and observation['loop_threading']['passed']):
                    stage.close()
                    stage = None
                    selected = probe(configuration, name, entry['steps'])
                    if selected is not None:
                        entry['selected_checkpoint_step'] = entry['steps']
                        record()
                        break
                if entry['steps'] == planned_steps:
                    break
                if stage is None:
                    seed = next_seed()
                    stage = Stage(request, configuration, kind=kind, relaxed=relaxed,
                        forces=forces if use_forces and pairs else None, seed=seed)
                    entry['segments'].append({'seed': seed, 'start_step': entry['steps'], 'steps': 0})
                count = min(sampling, planned_steps-entry['steps'])
                stage.run(count, deadline)
                entry['steps'] += count
                entry['segments'][-1]['steps'] += count
        finally:
            if stage is not None:
                stage.close()
            record()
        if selected is not None:
            break
    if selected is None:
        # Even a failed guide must remain inspectable when strict NA can run.
        # It is never labelled eligible merely because coordinates exist.
        selected = probe(configuration, stages[-1]['name'], stages[-1]['steps'], purpose='final_diagnostic')
    if selected is None:
        selected = best_probe
    if selected is None:
        raise ValueError('No prepared checkpoint could run under strict NA without forces.')
    configuration, observation, placement, checkpoint = selected
    validated = output/'validated.dat'
    validated.write_text(configuration.read_text())
    missing = sum(value['expected_pairs']-value['native_pairs'] for value in observation['guides'].values())
    warnings = ['Guided preparation is a nonequilibrium starting model; no kinetic rate is inferred.']
    if missing:
        warnings.append(f'{missing} expected 2D guide pairs are absent. This conformation is available for inspection, but FFS is disabled.')
    if observation['unexpected_pairs']:
        warnings.append('Additional native base pairs formed outside the 2D guide; FFS is disabled for this starting conformation.')
    if observation['interstrand_bonds']:
        warnings.append('The strands associated during preparation; FFS is disabled for this starting conformation.')
    return {'status': 'completed', 'configuration_file': str(validated), 'validation': observation,
            'box_size': placement['box_size'], 'final_placement': placement,
            'guide_minimum_fraction': 1., 'guide_accepted': observation['guide_accepted'],
            'off_basin_eligible': observation['off_basin_eligible'],
            'strict_force_field_initialization_passed': True,
            'production_external_forces': False, 'ffs_eligible': observation['ffs_eligible'],
            'preparation_attempt': attempt, 'selected_checkpoint': checkpoint,
            'relaxation': {'stages': stages, 'probes': probes, 'equilibrium_validated': False},
            'warnings': warnings}


def main():
    request_path=Path(sys.argv[1]).resolve()
    request=json.loads(request_path.read_text())
    output=request_path.parent
    try:
        result=run(request,output)
    except Exception as exc:
        result={'status':'failed','error':str(exc),
                'failure_kind': ('timeout' if isinstance(exc, TimeoutError) else
                                'environment' if isinstance(exc, (ImportError, FileNotFoundError, PermissionError))
                                else 'validation')}
    write(output/'guided_result.json',result)
    if result['status']!='completed':
        print(result['error'],file=sys.stderr)
        return 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
