"""Descriptive, PBC-aware observables for free DNA–DNA trajectories.

Distances use nucleotide centres, not a ligand surface or binding affinity.
Time-series dispersion is between independent replicas, never between frames.
"""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import numpy as np

LENGTH_NM = .8518
TIME_US = 3.031e-6
METRICS = {
    'minimum_trigger_site_distance_nm': {'label': 'Distance minimale trigger–zone', 'unit': 'nm'},
    'site_contact_fraction': {'label': 'Fraction de bases proches du trigger', 'unit': 'fraction'},
    'site_paired_fraction': {'label': 'Fraction de bases appariées dans la zone', 'unit': 'fraction'},
    'site_trigger_paired_fraction': {'label': 'Fraction de bases appariées au trigger', 'unit': 'fraction'},
    'interstrand_pairs': {'label': 'Paires aptamère–trigger', 'unit': 'pairs'},
    'rg_aptamer_nm': {'label': 'Rayon de giration de l’aptamère', 'unit': 'nm'},
    'rmsd_aptamer_nm': {'label': 'RMSD de l’aptamère', 'unit': 'nm'},
    'rmsd_site_nm': {'label': 'RMSD de la zone surveillée', 'unit': 'nm'},
}


def _unwrap(positions, box):
    positions = np.asarray(positions, dtype=float)
    increments = np.diff(positions, axis=0)
    increments -= box * np.rint(increments / box)
    return np.vstack([positions[0], positions[0] + np.cumsum(increments, axis=0)])


def frame_metrics(rows, box, hb_pairs, aptamer_length, monitored_bases,
                  reference_positions, contact_cutoff_nm=1.5):
    """Measure a frame; aptamer is strand 0, monitored indices are 1-based."""
    rows, box = np.asarray(rows, dtype=float), np.asarray(box, dtype=float)
    n = aptamer_length
    reference = np.asarray(reference_positions, dtype=float)
    if (rows.ndim != 2 or rows.shape[1] < 3 or type(n) is not int or not 1 <= n < len(rows)
            or box.shape != (3,) or not np.isfinite(rows).all() or not np.isfinite(box).all()
            or np.any(box <= 0) or reference.shape != (n, 3) or not np.isfinite(reference).all()):
        raise ValueError('Invalid DNA–DNA frame or reference.')
    if (not monitored_bases or len(set(monitored_bases)) != len(monitored_bases)
            or any(type(i) is not int or not 1 <= i <= n for i in monitored_bases)
            or not np.isfinite(contact_cutoff_nm) or contact_cutoff_nm <= 0):
        raise ValueError('Invalid monitored bases or proximity threshold.')
    indices = np.array(monitored_bases) - 1
    aptamer, reference = _unwrap(rows[:n, :3], box), _unwrap(reference, box)
    centre, reference_centre = aptamer.mean(axis=0), reference.mean(axis=0)
    centred, ref_centred = aptamer - centre, reference - reference_centre
    u, _, vt = np.linalg.svd(centred.T @ ref_centred)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    aligned = centred @ u @ correction @ vt + reference_centre
    squared = np.sum((aligned - reference) ** 2, axis=1)
    differences = rows[indices, None, :3] - rows[None, n:, :3]
    differences -= box * np.rint(differences / box)
    distances = np.sqrt(np.min(np.sum(differences ** 2, axis=2), axis=1)) * LENGTH_NM
    pairs = {tuple(sorted(pair)) for pair in hb_pairs}
    if any(len(pair) != 2 or any(type(i) not in (int, np.int64, np.int32)
            or not 0 <= i < len(rows) for i in pair) or pair[0] == pair[1] for pair in pairs):
        raise ValueError('Invalid native hydrogen-bond indices.')
    paired = {i for pair in pairs for i in pair}
    interstrand = {(i, j) for i, j in pairs if i < n <= j}
    paired_trigger = {i for i, _ in interstrand}
    site = [{'base': int(i + 1), 'min_trigger_distance_nm': float(distance),
             'near_trigger': bool(distance < contact_cutoff_nm), 'paired': int(i) in paired,
             'paired_to_trigger': int(i) in paired_trigger,
             'aligned_position_nm': (aligned[i] * LENGTH_NM).tolist()}
            for i, distance in zip(indices, distances)]
    return {
        'minimum_trigger_site_distance_nm': float(distances.min()),
        'site_contact_fraction': float(np.mean([row['near_trigger'] for row in site])),
        'site_paired_fraction': float(np.mean([row['paired'] for row in site])),
        'site_trigger_paired_fraction': float(np.mean([row['paired_to_trigger'] for row in site])),
        'interstrand_pairs': len(interstrand),
        'rg_aptamer_nm': float(np.sqrt(np.mean(np.sum(centred ** 2, axis=1))) * LENGTH_NM),
        'rmsd_aptamer_nm': float(np.sqrt(squared.mean()) * LENGTH_NM),
        'rmsd_site_nm': float(np.sqrt(squared[indices].mean()) * LENGTH_NM),
        'per_base': site,
    }


def _stats(values):
    values = np.asarray(values, dtype=float)
    return {'mean': float(values.mean()) if len(values) else None,
            'sd': float(values.std(ddof=1)) if len(values) > 1 else None, 'n': len(values)}


def _write_csv(path, rows, fieldnames=None):
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(stream.getvalue(), encoding='utf-8')
    temporary.replace(path)


def _read_samples(directory, state):
    """Only read transactionally committed observations, including on pause."""
    path = directory / 'samples.jsonl'
    count = int(state.get('samples_count', 0))
    if not count:
        return []
    previous_step = None
    with path.open(encoding='utf-8') as handle:
        for _ in range(count):
            row = json.loads(handle.readline())
            if any(key not in row or not np.isfinite(row[key]) for key in METRICS):
                raise ValueError('Invalid committed free MD observations.')
            if previous_step is not None and row['step'] <= previous_step:
                raise ValueError('Duplicate or out-of-order free MD observations.')
            previous_step = row['step']
            yield row


def _summarize_replica(directory, state, bases):
    """Stream site coordinates; memory does not grow with trajectory duration."""
    samples = []
    sums = {base: np.zeros(4) for base in bases}
    means = {base: np.zeros(3) for base in bases}
    m2 = {base: np.zeros(3) for base in bases}
    temporary = directory / 'metrics.csv.tmp'
    writer = None
    if not state.get('samples_count'):
        return [], []
    with temporary.open('w', encoding='utf-8', newline='') as handle:
        for count, sample in enumerate(_read_samples(directory, state), 1):
            row = {key: sample[key] for key in ('replica_index', 'variant', 'step', 'time_us', *METRICS)
                   if key in sample}
            row.update({key: sample[key] for key in ('native_step', 'time_model', 'energy') if key in sample})
            site_rows = {site['base']: site for site in sample['per_base']}
            if set(site_rows) != set(bases):
                raise ValueError('Committed site observations do not match monitored aptamer bases.')
            for base in bases:
                site = site_rows[base]
                values = [site[key] for key in ('min_trigger_distance_nm', 'near_trigger', 'paired', 'paired_to_trigger')]
                position = np.asarray(site['aligned_position_nm'])
                if position.shape != (3,) or not np.isfinite(position).all() or not np.isfinite(values).all():
                    raise ValueError('Invalid committed site coordinates.')
                sums[base] += values
                delta = position - means[base]
                means[base] += delta / count
                m2[base] += delta * (position - means[base])
                for key in ('min_trigger_distance_nm', 'near_trigger', 'paired', 'paired_to_trigger'):
                    row[f'base_{base}_{key}'] = site[key]
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
            writer.writerow(row)
            samples.append({key: sample[key] for key in ('step', *METRICS)})
    temporary.replace(directory / 'metrics.csv')
    count = len(samples)
    site = [{'base': base, **dict(zip(('mean_distance_nm', 'contact_fraction', 'paired_fraction',
             'paired_to_trigger_fraction'), (sums[base] / count).tolist())),
             'rmsf_nm': float(np.sqrt(max(0., m2[base].sum() / count))) if count > 1 else None}
            for base in bases]
    _write_csv(directory / 'site.csv', site)
    return samples, site


def _visualization(variant, bases):
    aptamer, trigger = variant['aptamer_sequence'], variant['trigger_sequence']
    return {'sequence': list(aptamer + trigger),
            'colors': ['#C44788' if i + 1 in bases else '#167D9A' for i in range(len(aptamer))]
                      + ['#D47E1E'] * len(trigger),
            'domains': ['Zone surveillée' if i + 1 in bases else 'Aptamère' for i in range(len(aptamer))]
                       + ['Trigger original' if variant['name'] == 'original' else 'Trigger étendu'] * len(trigger)}


def summarize_run(output, request):
    """Create per-replica CSVs, between-replica curves and portable movie paths."""
    output = Path(output).resolve()
    variants = {variant['name']: variant for variant in request['variants']}
    replicas, all_samples, curves = [], {}, []
    for index in range(request['replicas']):
        directory = output / f'replica_{index:03d}'
        state_path = directory / 'state.json'
        state = json.loads(state_path.read_text()) if state_path.is_file() else {}
        # The native controller uses alternating arms to balance partial runs.
        name = state.get('variant', 'original' if index % 2 == 0 else 'extended')
        variant = variants[name]
        samples, site = _summarize_replica(directory, state, request['monitored_bases'])
        all_samples[index] = samples
        relative = directory.relative_to(output).as_posix()
        topology = Path(variant['topology_file']).resolve()
        if not topology.is_relative_to(output):
            raise ValueError('Free MD topology is outside the calculation directory.')
        replica = {'replica_index': index, 'variant': name, 'status': state.get('status', 'pending'),
                   'production_steps': state.get('production_steps', 0), 'sample_count': len(samples),
                   'summary': {}, 'site': site, 'topology_file': topology.relative_to(output).as_posix(),
                   'visualization': _visualization(variant, request['monitored_bases']),
                   'preparation': state.get('preparation', [])}
        if samples:
            replica.update(trajectory_file=relative + '/trajectory.dat', metrics_file=relative + '/metrics.csv',
                           movie_kind='free_md')
            replica['summary'] = {key: {'mean': float(np.mean([s[key] for s in samples])),
                'min': float(min(s[key] for s in samples)), 'max': float(max(s[key] for s in samples))}
                for key in METRICS}
        replicas.append(replica)
    groups = {}
    for name in ('original', 'extended'):
        members = [replica for replica in replicas if replica['variant'] == name]
        observed = [replica for replica in members if replica['sample_count']]
        group = {'completed': sum(row['status'] == 'completed' for row in members),
                 'available': len(observed), 'summary': {}, 'site': []}
        for key in METRICS:
            group['summary'][key] = _stats([row['summary'][key]['mean'] for row in observed])
        for base in request['monitored_bases']:
            rows = [next(row for row in replica['site'] if row['base'] == base) for replica in observed]
            row = {'base': base, 'n': len(rows)}
            for key in ('mean_distance_nm', 'contact_fraction', 'paired_fraction', 'paired_to_trigger_fraction', 'rmsf_nm'):
                stats = _stats([r[key] for r in rows if r[key] is not None])
                row.update({key: stats['mean'], key + '_sd': stats['sd'], key + '_n': stats['n']})
            group['site'].append(row)
        groups[name] = group
        by_step = {}
        for replica in observed:
            for sample in all_samples[replica['replica_index']]:
                by_step.setdefault(sample['step'], []).append(sample)
        for step, observations in sorted(by_step.items()):
            for key in METRICS:
                curves.append({'variant': name, 'step': step,
                               'time_us': step * request['sampling']['dt'] * TIME_US,
                               'metric': key, **_stats([sample[key] for sample in observations])})
    complete = all(row['status'] == 'completed' and row['production_steps'] == request['sampling']['production_steps']
                   for row in replicas)
    comparison = {'status': 'complete' if complete else 'partial', 'expected_per_group': request['replicas'] // 2,
                  'groups': groups, 'curves': curves, 'metrics': METRICS,
                  'dispersion': 'sample_standard_deviation_between_independent_replicas',
                  'summary_window': 'observed_production_frames_per_replica',
                  'rmsf_reference': 'aptamer_alignment_to_first_production_frame',
                  'distance_definition': 'minimum_periodic_distance_between_nucleotide_centres',
                  'contact_cutoff_nm': request.get('contact_cutoff_nm', 1.5),
                  'ligand_simulated': False}
    _write_csv(output / 'comparison_curves.csv', curves,
               ['variant', 'step', 'time_us', 'metric', 'mean', 'sd', 'n'])
    _write_csv(output / 'comparison_summary.csv', [dict(variant=name, metric=key, **stats)
               for name, group in groups.items() for key, stats in group['summary'].items()])
    _write_csv(output / 'comparison_site.csv', [dict(variant=name, **row)
               for name, group in groups.items() for row in group['site']])
    return {'sampling_complete': complete, 'request': request, 'replicas': replicas, 'comparison': comparison}
