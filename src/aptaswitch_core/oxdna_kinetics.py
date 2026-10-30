"""Direct FFS protocol and finite-box rate estimates in oxDNA time units.

The source flux counts first crossings coming from A, with B-history exposure
excluded. Native FFS_MD stops all trajectories. No NAM approximation or prepared
bound seed participates in this protocol.
"""
from __future__ import annotations
import hashlib
import json
import math
import statistics
from collections.abc import Mapping, Sequence

OPENING_PROTOCOL_VERSION = 3
PROXIMITY_PROTOCOL_VERSION = 4
PROXIMITY_EQUILIBRATION_METHOD = 'isolated_monomers_before_close_placement_v1'
PROXIMITY_DISTANCE_MIN = 2.0
PROXIMITY_DISTANCE_MAX = 2.5

# Public production controls. Checkpoint loading and the native worker retain
# their frozen allocation instead of applying new-form defaults retroactively.
FFS_SOURCE_SAMPLES_MIN = 75
FFS_SOURCE_SAMPLES_MAX = 200
FFS_SOURCE_SAMPLES_DEFAULT = FFS_SOURCE_SAMPLES_MIN
FFS_TRIALS_MIN = 250
FFS_TRIALS_MAX = 1000
FFS_TRIALS_DEFAULT = FFS_TRIALS_MIN


def validate_kinetic_protocol(options: Mapping) -> None:
    mode, version = options.get('mode'), options.get('version')
    if (mode not in ('opening', 'proximity') or type(version) is not int
            or version not in ((OPENING_PROTOCOL_VERSION,) if mode == 'opening' else
                               (OPENING_PROTOCOL_VERSION, PROXIMITY_PROTOCOL_VERSION))):
        raise ValueError('This run does not use native direct FFS. Start a new simulation.')
    if mode == 'proximity':
        expected = {'source_samples': 0,
                    'initial_distance_min': PROXIMITY_DISTANCE_MIN,
                    'initial_distance_max': PROXIMITY_DISTANCE_MAX}
        if version == OPENING_PROTOCOL_VERSION:
            # Historical sources retain their original no-MD preparation.
            expected['equilibration_steps'] = 0
        elif (type(options.get('equilibration_steps')) is not int or options['equilibration_steps'] <= 0
              or options.get('equilibration_method') != PROXIMITY_EQUILIBRATION_METHOD):
            raise ValueError('New proximity sources require isolated native monomer equilibration before close placement.')
        if any(isinstance(options.get(key), bool) or options.get(key) != value
               for key, value in expected.items()):
            raise ValueError('The proximity protocol requires its frozen unbound initial geometry and no source flux.')


def _count(value, name: str, *, positive=False) -> int:
    if type(value) is not int or value < int(positive):
        raise ValueError(f'{name} must be a nonnegative integer (positive when required).')
    return value


def candidate_interfaces(plan, *, mode='opening', initial_distance_max=PROXIMITY_DISTANCE_MAX):
    # Nested unions are used by the native condition writer. A real intended
    # bond therefore qualifies even in an unusual distance configuration.
    items = [{'name': f'contact_{d}', 'distance_max': float(d), 'min_binding': 0}
             for d in (8, 7, 6, 5, 4, 3, 2)]
    items.append({'name': 'duplex_1', 'min_binding': 1})
    final = dict(plan['interfaces'][-1])
    items += [{'name': f'duplex_{i}', 'min_binding': i}
              for i in range(2, final['min_binding'] + 1)]
    items += [{'name': f'stem_{i}', 'min_binding': final['min_binding'], 'max_stem': i}
              for i in range(len(plan['stem_pairs'])-1, final['max_stem'], -1)]
    items += [{**final, 'name': f'opening_{i}', 'max_opening': i} for i in (6, 3, 1)]
    items.append(final)
    if mode == 'proximity':
        return [{'name': 'near_toehold', 'distance_max': initial_distance_max, 'min_binding': 0},
                *[item for item in items if 'distance_max' not in item or item['distance_max'] < initial_distance_max]]
    if mode != 'opening':
        raise ValueError('Unknown FFS simulation mode.')
    return items


def progress_index(observation, interfaces):
    def get(key, default):
        return observation.get(key, default) if isinstance(observation, dict) else getattr(observation, key, default)
    for index in range(len(interfaces)-1, -1, -1):
        target = interfaces[index]
        if (get('binding', 0) >= target.get('min_binding', 0)
                and all(key not in target or get(field, math.inf) <= target[key]
                        for key, field in [('distance_max', 'distance'), ('max_stem', 'stem'), ('max_opening', 'opening')])):
            return index
    return -1


def condition(target):
    expressions = []
    for key, name, operator in [('min_binding', 'binding', '>='), ('distance_max', 'distance', '<='),
                                 ('max_stem', 'stem', '<='), ('max_opening', 'opening', '<=')]:
        if key in target and (key != 'min_binding' or target[key] > 0):
            expressions.append(f'{name} {operator} {target[key]}')
    return expressions


def source_fingerprint(request):
    # Random seed and allocation are absent from this physical provenance.
    # The full checkpoint identity additionally freezes those choices.
    values = {key: request.get(key) for key in ('saved_switch_identity', 'switch_sequence', 'trigger_sequence',
        'switch_material', 'trigger_material', 'model', 'box_size', 'switch_guide', 'trigger_guide', 'plan')}
    values['sampling'] = {key: request['sampling'].get(key, default) for key, default in (
        ('temperature_c', None), ('salt_m', None), ('dt', .003), ('newtonian_steps', 103), ('diff_coeff', 2.5))}
    values['parameters'] = request['runtime'].get('parameter_sha256')
    values['protocol'] = request.get('kinetics', {}).get('version', OPENING_PROTOCOL_VERSION)
    if (request.get('kinetics') or {}).get('mode') == 'proximity':
        values['initial_condition'] = {key: request['kinetics'][key]
                                       for key in ('mode', 'initial_distance_min', 'initial_distance_max')}
        values['initial_condition']['ensemble'] = 'prepared_unbound_proximity_v1'
        if request['kinetics'].get('version') == PROXIMITY_PROTOCOL_VERSION:
            values['initial_condition'].update(
                ensemble='isolated_monomer_md_then_unbound_proximity_v2',
                equilibration_method=request['kinetics']['equilibration_method'],
                equilibration_steps_per_monomer=request['kinetics']['equilibration_steps'])
    if request.get('protocol_amendments'):
        # An exploratory continuation cannot seed a fresh kinetic estimate.
        values['protocol_amendments'] = request['protocol_amendments']
    return hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()


def fixed_trial_probability(stages: Sequence[Mapping]) -> dict:
    """Product of conditional probabilities from frozen production batches.

    Each stage supplies target_trials, attempted, successes, failures,
    censored and pending (the latter two default to zero). Only explicit
    failures count as failures. Unknown or unfinished shots never do.
    A completed zero-success prefix yields a zero point estimate without
    inventing downstream trials, and is labelled ``no_success``.
    """
    if not isinstance(stages, Sequence) or isinstance(stages, (str, bytes)) or not stages:
        raise ValueError("At least one fixed-trial stage is required.")
    rows = []
    for index, stage in enumerate(stages):
        if not isinstance(stage, Mapping):
            raise ValueError("Every stage must contain explicit trial counts.")
        try:
            counts = {key: _count(stage[key], key, positive=key == "target_trials")
                      for key in ("target_trials", "attempted", "successes", "failures")}
        except KeyError as exc:
            raise ValueError(f"Missing fixed-trial count: {exc.args[0]}.") from exc
        counts.update({key: _count(stage.get(key, 0), key) for key in ("censored", "pending")})
        terminal = counts["successes"] + counts["failures"]
        if terminal + counts["censored"] + counts["pending"] != counts["attempted"]:
            raise ValueError("Every attempted trial needs an explicit outcome or pending/censored status.")
        if counts["attempted"] > counts["target_trials"]:
            raise ValueError("Production attempts exceed the frozen trial allocation.")
        complete = terminal == counts["target_trials"]
        rows.append({"name": str(stage.get("name", index)), **counts,
                     "complete": complete,
                     "probability": counts["successes"] / terminal if complete else None})

    log_beta = 0.
    status, beta, zero_index = "complete", None, None
    for index, row in enumerate(rows):
        if zero_index is not None:
            if row["attempted"]:
                raise ValueError("Downstream trials cannot follow an empty parent bank.")
            continue
        if not row["complete"]:
            status = "incomplete"
            break
        if row["successes"] == 0:
            zero_index = index
            status, beta, log_beta = "no_success", 0., None
        else:
            log_beta += math.log(row["probability"])
    if status == "complete":
        beta = math.exp(log_beta)
        if beta == 0:
            status, beta = "underflow", None
    return {"status": status, "beta": beta, "log_beta": log_beta if status != "incomplete" else None,
            "stages": rows, "zero_success_stage": zero_index,
            "censored": sum(row["censored"] for row in rows),
            "pending": sum(row["pending"] for row in rows),
            "uncertainty_unit": "independent_replica_not_descendant_trial"}


# Two-sided 95% Student-t critical values, df 1..30. Above this, the
# Cornish-Fisher expansion converges rapidly; the interval remains approximate.
_T95 = (12.706205, 4.302653, 3.182446, 2.776445, 2.570582, 2.446912,
        2.364624, 2.306004, 2.262157, 2.228139, 2.200985, 2.178813,
        2.160369, 2.144787, 2.131450, 2.119905, 2.109816, 2.100922,
        2.093024, 2.085963, 2.079614, 2.073873, 2.068658, 2.063899,
        2.059539, 2.055529, 2.051831, 2.048407, 2.045230, 2.042272)


def _t95(df: int) -> float:
    if df <= len(_T95):
        return _T95[df - 1]
    z = statistics.NormalDist().inv_cdf(.975)
    return (z + (z**3 + z) / (4 * df)
            + (5*z**5 + 16*z**3 + 3*z) / (96 * df**2)
            + (3*z**7 + 19*z**5 + 17*z**3 - 15*z) / (384 * df**3))


def direct_ffs_report(replicas, *, dt, concentration_m, protocol_amendments=None):
    for value, label in ((dt, 'time step'), (concentration_m, 'box concentration')):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'The {label} must be finite and positive.')
    rows, identities = [], set()
    for replica in replicas:
        identity = replica.get('replica_id')
        if not isinstance(identity, str) or not identity or identity in identities:
            raise ValueError('Every independent replica needs a distinct identity.')
        identities.add(identity)
        flux = replica.get('flux', {})
        steps, crossings = _count(flux.get('exposure_steps', 0), 'source exposure'), len(flux.get('bank', []))
        probability = fixed_trial_probability(replica['stages']) if replica.get('stages') else {'beta': None}
        completed = bool(replica.get('completed', False) and flux.get('complete') and steps > 0 and crossings > 0)
        phi = crossings / (steps * dt) if steps > 0 else None
        beta = probability['beta'] if completed else None
        value = phi * beta if phi is not None and beta is not None else None
        rows.append({'replica_id': replica['replica_id'], 'flux_crossings': crossings,
                     'a_history_time': steps * dt, 'flux': phi, 'probability': beta,
                     'value_model': value, 'complete': completed, 'stages': probability.get('stages', [])})
    complete = bool(rows) and all(r['value_model'] is not None for r in rows)
    values = [r['value_model'] for r in rows] if complete else []
    mean = statistics.fmean(values) if values else None
    interval = None
    # Across whole repetitions: descendants are correlated through parents.
    # An all-zero/identical sample must not manufacture zero-width certainty.
    if len(values) >= 4 and any(values) and statistics.stdev(values) > 0:
        half = _t95(len(values)-1) * statistics.stdev(values) / math.sqrt(len(values))
        interval = [max(0., mean-half), mean+half]
    report = {'status': ('model_estimate' if any(values) else 'no_success') if complete else 'incomplete',
            'value_model': mean, 'inverse_model_time': mean, 'interval95_model': interval,
            'unit': '1/oxDNA_time', 'time_basis': 'model_time_not_experimental_seconds',
            'value_m_inv_model_time': mean/concentration_m if mean is not None else None,
            'concentration_each_m': concentration_m, 'finite_box_estimate': True,
            'replicas': rows, 'experimental_rate': None,
            'uncertainty_method': 'approximate_student_t_across_whole_replicas',
            'formula': 'mean_r(flux_r * product_i(P_r_i))',
            'limitations': ['Finite-box model rate; physical-time calibration is not established.',
                            'Equilibration and source decorrelation require convergence assessment.',
                            'Between-replica uncertainty does not measure force-field accuracy.']}
    if protocol_amendments:
        report.update(status='exploratory_modified_protocol', value_model=None, inverse_model_time=None,
                      interval95_model=None, value_m_inv_model_time=None, uncertainty_method=None,
                      protocol_amendments=protocol_amendments)
        for row in rows:
            row.update(flux=None, probability=None, value_model=None)
            for stage in row['stages']:
                stage['probability'] = None
        report['limitations'].insert(0, 'Protocol changed during this test run; retained configurations do not support a kinetic estimate.')
    return report


def proximity_ffs_report(replicas, *, protocol_amendments=None):
    """Opening-before-separation probability for the prepared, unbound ensemble.

    Each independently placed source contributes one whole branching estimate.
    There is no encounter flux, waiting-time estimate or concentration scaling.
    """
    replicas = list(replicas)
    equilibrated = [bool(replica.get('monomer_equilibration')) for replica in replicas]
    if any(equilibrated) and not all(equilibrated):
        raise ValueError('Legacy and equilibrated proximity sources cannot be averaged together.')
    rows, identities = [], set()
    for replica in replicas:
        identity = replica.get('replica_id')
        if not isinstance(identity, str) or not identity or identity in identities:
            raise ValueError('Every independent replica needs a distinct identity.')
        identities.add(identity)
        probability = fixed_trial_probability(replica['stages']) if replica.get('stages') else {'beta': None}
        complete = bool(replica.get('completed') and replica.get('initial_source'))
        if replica.get('monomer_equilibration'):
            complete &= bool(replica['monomer_equilibration'].get('completed'))
        value = probability['beta'] if complete else None
        rows.append({'replica_id': identity, 'probability': value, 'complete': complete and value is not None,
                     'stages': probability.get('stages', [])})
    complete = bool(rows) and all(row['complete'] for row in rows)
    values = [row['probability'] for row in rows] if complete else []
    mean, interval = (statistics.fmean(values) if values else None), None
    if len(values) >= 4 and statistics.stdev(values) > 0:
        half = _t95(len(values)-1) * statistics.stdev(values) / math.sqrt(len(values))
        interval = [max(0., mean-half), min(1., mean+half)]
    report = {
        'status': ('conditional_estimate' if any(values) else 'no_success') if complete else 'incomplete',
        'value': mean, 'interval95': interval, 'unit': 'probability', 'endpoint': 'opening',
        'initial_condition': {'ensemble': 'prepared_unbound_proximity_v1',
                              'distance_min': PROXIMITY_DISTANCE_MIN, 'distance_max': PROXIMITY_DISTANCE_MAX,
                              'distance_unit': 'oxDNA_length', 'interstrand_bonds': 0,
                              'equilibrium_sample': False},
        'replicas': rows, 'formula': 'mean_r(product_i(P_r_i))',
        'uncertainty_method': 'approximate_student_t_across_whole_replicas',
        'limitations': ['Conditional on the artificial unbound proximity preparation, not a natural encounter ensemble.',
                        'Opening is counted before return to separated strands; failures and reversals are retained.',
                        'No encounter flux, full kinetic rate or experimental waiting time is estimated.',
                        'Source diversity and between-replica convergence require assessment.'],
    }
    if any(equilibrated):
        targets = {replica['monomer_equilibration'].get('target_steps_per_monomer',
                    replica['monomer_equilibration'].get('target_steps_each')) for replica in replicas}
        if len(targets) != 1:
            raise ValueError('Proximity sources must use the same isolated-monomer burn-in allocation.')
        report['initial_condition'].update(
            ensemble='isolated_monomer_md_then_unbound_proximity_v2',
            equilibration_method=PROXIMITY_EQUILIBRATION_METHOD,
            equilibration_steps_per_monomer=targets.pop(), monomer_equilibrium_validated=False)
        report['limitations'].append('Finite isolated-monomer burn-in is recorded; convergence is not established by its step count.')
    if protocol_amendments:
        report.update(status='exploratory_modified_protocol', value=None, interval95=None,
                      uncertainty_method=None, protocol_amendments=protocol_amendments)
        for row in rows:
            row['probability'] = None
            for stage in row['stages']:
                stage['probability'] = None
    return report
