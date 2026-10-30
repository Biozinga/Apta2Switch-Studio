"""Small, read-only summaries and one interface chart for a selected FFS replica."""
from __future__ import annotations

from collections.abc import Mapping
import math

import altair as alt
import pandas as pd

from aptaswitch_core.localization import tr
from aptaswitch_core.oxdna_box import OXDNA_LENGTH_NM
from aptaswitch_core.oxdna_kinetics import candidate_interfaces


_CHART_WINDOW = 10


def distance_label(plan, *, column=False):
    """Describe the frozen distance observable without relabelling older runs."""
    labels = {
        'trigger_3prime_to_toehold_start': (
            'Distance trigger 3′–début du toehold', 'Distance 3′–toehold (nm)'),
        'trigger_to_toehold': (
            'Distance minimale trigger–toehold', 'Distance toehold (nm)'),
    }
    full, compact = labels.get((plan or {}).get('distance_scope'), (
        'Distance minimale entre partenaires cibles', 'Distance (nm)'))
    return tr(compact if column else full)


def distance_description(plan):
    scope = (plan or {}).get('distance_scope')
    if scope == 'trigger_3prime_to_toehold_start':
        return tr('La distance mesure uniquement l’écart entre la dernière base du trigger côté 3′ et la première base du toehold après le leader du switch. Elle ne contraint pas leur orientation. Un essai s’arrête sur un retour à l’état séparé uniquement quand aucune liaison interbrin ne subsiste et que tous les nucléotides des deux brins sont éloignés.')
    if scope == 'trigger_to_toehold':
        return tr('La distance suit le contact le plus proche entre le trigger et le toehold. Un essai s’arrête sur un retour à l’état séparé uniquement quand aucune liaison interbrin ne subsiste et que tous les nucléotides des deux brins sont éloignés.')
    return ''


def _count(value):
    return value if type(value) is int and value >= 0 else None


def _number(value):
    return (value if isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value >= 0 else None)


def _interfaces(request):
    plan = request.get('plan') or {}
    if not plan.get('interfaces'):
        return []
    if plan['interfaces'][0].get('name') in {'contact_8', 'near_toehold'}:
        return plan['interfaces']
    return candidate_interfaces(plan, mode=(request.get('kinetics') or {}).get('mode', 'opening'))


def _interface_label(interface, plan=None):
    """Keep physical units and predicates visible instead of internal IDs."""
    name = interface['name']
    if name == 'near_toehold':
        return tr('Départ proche du toehold')
    if name == 'opening':
        return tr('Ouverture RBS')
    if 'max_opening' in interface:
        return tr('RBS : ≤ {count} liaison' if interface['max_opening'] == 1 else 'RBS : ≤ {count} liaisons',
                  count=interface['max_opening'])
    if 'max_stem' in interface:
        return tr('Tige : ≤ {count} pb', count=interface['max_stem'])
    if 'distance_max' in interface:
        if (plan or {}).get('distance_scope') == 'trigger_3prime_to_toehold_start':
            return tr('3′–toehold : ≤ {distance:.2f} nm', distance=interface['distance_max'] * OXDNA_LENGTH_NM)
        return tr('Distance cible : ≤ {distance:.2f} nm', distance=interface['distance_max'] * OXDNA_LENGTH_NM)
    if 'min_binding' in interface:
        return tr('Duplex : {count} pb', count=interface['min_binding'])
    return name


def _zero_success(state):
    return any(stage.get('completed') and _count(stage.get('successes')) == 0
               and (_count(stage.get('target_trials')) or 0) > 0
               and _count(stage.get('failures')) == stage['target_trials']
               for stage in state.get('stages') or [])


def _replica_status(state, status):
    own_status = state.get('status')
    if own_status == 'failed':
        return tr('Échec')
    if state.get('completed'):
        return tr('Terminée sans succès' if _zero_success(state) else 'Terminée')
    if not state:
        return tr('En attente' if status == 'running' else 'Non démarrée')
    if own_status in {'cancelled', 'paused'} or status in {'cancelled', 'paused'}:
        return tr('En pause')
    if own_status == 'interrupted' or status in {'failed', 'interrupted'}:
        return tr('Interrompue')
    return tr({'initialization': 'Préparation du départ proche',
               'equilibration_switch': 'Équilibrage du switch isolé',
               'equilibration_trigger': 'Équilibrage du trigger isolé',
               'equilibration': 'Équilibrage', 'flux': 'Collecte initiale',
               'production': 'Essais FFS'}.get(state.get('phase'), 'En attente'))


def _fraction(done, target, template):
    if _count(done) is None or not _count(target):
        return '—'
    return tr(template, done=done, target=target)


def _current_stage(state):
    stages = state.get('stages') or []
    return next((stage for stage in stages if not stage.get('completed')), stages[-1] if stages else {})


def _monomer_counts(state):
    values = state.get('monomer_equilibration') or {}
    if not isinstance(values, Mapping):
        return None
    target = _count(values.get('target_steps_each'))
    switch, trigger = _count(values.get('switch_steps')), _count(values.get('trigger_steps'))
    if (not target or switch is None or trigger is None or switch > target or trigger > target
            or _count(values.get('total_steps')) != 2 * target
            or _count(values.get('completed_steps')) != switch + trigger):
        return None
    return switch, trigger, target


def replica_rows(states, request, *, status='running'):
    """One row per configured replica, including queued and paused replicas.

    ``states`` contains the small native ``progress.json`` payloads, either a
    mapping by zero-based replica index or a list with each record's ``index``.
    No trajectory, checkpoint or file is read, and no global percentage is made.
    """
    options = request.get('kinetics') or {}
    count = _count(options.get('replicas')) or 0
    records = {}
    entries = states.items() if isinstance(states, Mapping) else enumerate(states or [])
    for key, state in entries:
        if not isinstance(state, Mapping):
            continue
        index = state.get('index', key)
        if type(index) is int and 0 <= index < count:
            records[index] = state
    interfaces = _interfaces(request)
    labels = {interface['name']: _interface_label(interface, request.get('plan')) for interface in interfaces}
    monomer_columns = options.get('mode') == 'proximity' and (
        options.get('equilibration_method') == 'isolated_monomers_before_close_placement_v1'
        or any(state.get('monomer_equilibration') for state in records.values()))
    rows = []
    for index in range(count):
        state = records.get(index, {})
        phase = state.get('phase')
        step, advancement, successes = '—', '—', '—'
        if phase == 'initialization':
            step = tr('Départ proche du toehold')
        elif phase in {'equilibration_switch', 'equilibration_trigger'}:
            step = tr('Équilibrage du switch isolé' if phase == 'equilibration_switch'
                      else 'Équilibrage du trigger isolé')
            if counts := _monomer_counts(state):
                switch, trigger, target = counts
                advancement = _fraction(switch if phase == 'equilibration_switch' else trigger,
                                        target, '{done}/{target} pas')
        elif phase == 'equilibration':
            step = tr('Équilibrage')
            advancement = _fraction(state.get('burn_steps'), options.get('equilibration_steps'),
                                    '{done}/{target} pas')
        elif phase == 'flux':
            flux = state.get('flux') or {}
            step = tr({'initial_A': 'Attente de l’état séparé',
                       'forward': 'Collecte initiale', 'return': 'Retour vers l’état séparé',
                       'B_history': 'Retour naturel depuis ON'}.get(flux.get('history'), 'Collecte initiale'))
            advancement = _fraction(flux.get('crossings'), options.get('source_samples'),
                                    '{done}/{target} franchissements')
        elif phase in {'production', 'completed', 'report'} or state.get('completed'):
            stage = _current_stage(state)
            step = labels.get(stage.get('name'), stage.get('name', '—'))
            done = [_count(stage.get(field)) for field in ('successes', 'failures')]
            if all(value is not None for value in done):
                advancement = _fraction(sum(done), stage.get('target_trials', options.get('trials')),
                                        '{done}/{target} essais')
                successes = str(done[0])
        observation = state.get('observation') or {}
        binding, distance = _count(observation.get('binding')), _number(observation.get('distance'))
        row = {tr('Réplique'): index + 1, tr('État'): _replica_status(state, status),
                     tr('Étape'): step, tr('Avancement étape'): advancement,
                     tr('Succès interface'): successes, tr('Paires cibles'): str(binding) if binding is not None else '—',
                     distance_label(request.get('plan'), column=True): f'{distance * OXDNA_LENGTH_NM:.2f}' if distance is not None else '—'}
        if monomer_columns:
            counts = _monomer_counts(state)
            for label, done, target in (
                ('Équilibrage switch', counts[0] if counts else None, counts[2] if counts else None),
                ('Équilibrage trigger', counts[1] if counts else None, counts[2] if counts else None),
                ('Équilibrage total', sum(counts[:2]) if counts else None, 2 * counts[2] if counts else None),
            ):
                row[tr(label)] = _fraction(done, target, '{done}/{target} pas')
        rows.append(row)
    return rows


def selected_replica_chart(state, request, *, status='running'):
    """Recent production-interface successes; None until production starts.

    Show at most ten adjacent interfaces, including the next unstarted one.
    The count axis stays fixed at this replica's frozen allocation. Flux λ₀
    crossings are deliberately excluded: those are source events, not trials.
    """
    if state.get('phase') in {'equilibration_switch', 'equilibration_trigger'}:
        return None
    interfaces = _interfaces(request)[1:]
    stages = {stage['name']: stage for stage in state.get('stages') or []}
    begun = [i for i, interface in enumerate(interfaces) if interface['name'] in stages]
    if not begun:
        return None
    last = max(begun)
    end = min(len(interfaces), last + 2)
    start = max(0, end - _CHART_WINDOW)
    visible = interfaces[start:end]
    allocation = _count((request.get('kinetics') or {}).get('trials')) or 0
    maximum = max([allocation, *[(_count(stage.get('target_trials')) or 0) for stage in stages.values()], 1])
    rows = []
    for interface in visible:
        stage = stages.get(interface['name'], {})
        successes = _count(stage.get('successes')) or 0
        failures = _count(stage.get('failures')) or 0
        target = _count(stage.get('target_trials')) or allocation
        if not stage:
            state_label = tr('Non atteinte' if state.get('completed') and _zero_success(state) else 'À venir')
        elif stage.get('completed'):
            state_label = tr('Sans succès' if not successes else 'Terminée')
        elif _replica_status(state, status) in {tr('En pause'), tr('Interrompue'), tr('Échec')}:
            state_label = _replica_status(state, status)
        else:
            state_label = tr('En cours')
        rows.append({'interface': interface['name'], 'label': _interface_label(interface, request.get('plan')),
                     'successes': successes, 'failures': failures, 'finished': successes + failures,
                     'target': target, 'pending': _count(stage.get('pending')) or 0, 'state': state_label,
                     'caption': tr('{successes} succès · {done}/{target} essais',
                                   successes=successes, done=successes + failures, target=target),
                     'caption_position': maximum * .025})
    palette = {tr(label): color for label, color in [
        ('Terminée', '#16845b'), ('En cours', '#1258dc'), ('Sans succès', '#b45309'),
        ('En pause', '#64748b'), ('Interrompue', '#64748b'), ('Échec', '#b91c1c'),
        ('À venir', '#cbd5e1'), ('Non atteinte', '#cbd5e1')]}
    present = [label for label in palette if any(row['state'] == label for row in rows)]
    colors = alt.Color('state:N', title=None, legend=alt.Legend(orient='bottom'), scale=alt.Scale(
        domain=present, range=[palette[label] for label in present]))
    return alt.layer(
        alt.Chart().mark_bar(color='#edf2f7').encode(x='target:Q'),
        alt.Chart().mark_bar(opacity=.7).encode(x='successes:Q', color=colors),
        alt.Chart().mark_point(filled=True, size=32).encode(x='successes:Q', color=colors),
        alt.Chart().mark_text(align='left', color='#14243b', fontSize=11).encode(x='caption_position:Q', text='caption:N'),
        data=pd.DataFrame(rows),
    ).encode(
        y=alt.Y('label:N', title=None, sort=[row['label'] for row in rows],
                axis=alt.Axis(labelOverlap=False, labelLimit=250)),
        x=alt.X(type='quantitative', title=tr('Succès conservés'), scale=alt.Scale(domain=[0, maximum], nice=False),
                axis=alt.Axis(tickMinStep=1)),
        tooltip=[alt.Tooltip('label:N', title=tr('Interface')), alt.Tooltip('state:N', title=tr('État')),
                 alt.Tooltip('successes:Q', title=tr('Succès'), format='d'),
                 alt.Tooltip('finished:Q', title=tr('Essais terminés'), format='d'),
                 alt.Tooltip('target:Q', title=tr('Essais prévus'), format='d'),
                 alt.Tooltip('failures:Q', title=tr('Retours à l’état séparé'), format='d'),
                 alt.Tooltip('pending:Q', title=tr('Essais en cours'), format='d')],
    ).properties(height=max(105, 28 * len(rows)),
                 autosize=alt.AutoSizeParams(type='fit-x', contains='padding'))


def trajectory_observables_chart(observables, sampling):
    """Bond counts along the selected successful FFS lineage in model time."""
    dt = _number((sampling or {}).get('dt'))
    if not dt or not isinstance(observables, list) or len(observables) < 2:
        return None
    steps = []
    for point in observables:
        if not isinstance(point, Mapping):
            return None
        step = _count(point.get('step'))
        if step is None or (steps and step < steps[-1]):
            return None
        if any(_count(point.get(field)) is None for field in ('binding', 'stem', 'opening')):
            return None
        steps.append(step)
    if steps[-1] == steps[0]:
        return None
    series = [('binding', tr('Paires trigger–switch'), '#2563eb'),
              ('stem', tr('Paires de la tige OFF'), '#e07819'),
              ('opening', tr('Liaisons RBS–AUG'), '#7c3aed')]
    rows = [{'time': (point['step'] - steps[0]) * dt, 'series': label, 'bonds': point[field]}
            for point in observables for field, label, _ in series]
    return alt.Chart(pd.DataFrame(rows)).mark_line(strokeWidth=2).encode(
        x=alt.X('time:Q', title=tr('Temps oxDNA le long du trajet FFS')),
        y=alt.Y('bonds:Q', title=tr('Nombre de paires / liaisons'), scale=alt.Scale(zero=True)),
        color=alt.Color('series:N', title=None,
                        scale=alt.Scale(domain=[label for _, label, _ in series],
                                        range=[color for _, _, color in series]),
                        legend=alt.Legend(orient='bottom')),
        tooltip=[alt.Tooltip('time:Q', title=tr('Temps oxDNA'), format='.5g'),
                 alt.Tooltip('series:N', title=tr('Mesure')),
                 alt.Tooltip('bonds:Q', title=tr('Nombre'), format='d')],
    ).properties(height=260)
