"""Native oxDNA FFS: source flux, first-contact bank and opening trajectories."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import streamlit as st

from aptaswitch_core.localization import set_language, tr
from aptaswitch_core.thermodynamics.setup import configured_library_path
from aptaswitch_core.oxdna_ffs import ffs_visualization
from aptaswitch_core.oxdna_kinetics import (
    FFS_SOURCE_SAMPLES_DEFAULT, FFS_SOURCE_SAMPLES_MIN, FFS_SOURCE_SAMPLES_MAX,
    FFS_TRIALS_DEFAULT, FFS_TRIALS_MIN, FFS_TRIALS_MAX,
)
from aptaswitch_core.oxdna_setup import detect_oxdna_runtime
from .scientific_engines_ui import native_backend_identity
from .design_session import active_job_running
from .ffs_progress_charts import distance_description, distance_label
from .concentration_ui import render_concentration_control
from .saved_switch_ffs_runtime import active_ffs, resume_ffs, saved_ffs_results, start_ffs
from .table_ui import render_readonly_table


def render_oxdna_movie(result: dict, run_dir: Path, key: str):
    from .oxdna_movie import render_oxdna_movie as render

    return render(result, run_dir, key)


def _key(saved):
    return "saved_switch_ffs_" + hashlib.sha256(str(saved.path.resolve()).encode()).hexdigest()[:20]


def _start(saved, runtime, options):
    try:
        if active_job_running(st.session_state):
            return
        job, _ = start_ffs(saved, runtime, {**options, "thermodynamics_library_path": configured_library_path() or ""})
        st.session_state.saved_switch_ffs_job = job
        st.session_state.pop(_key(saved) + "_error", None)
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        st.session_state[_key(saved) + "_error"] = str(exc)


def _resume(saved, run_dir):
    """Resume the frozen run, independently of the current form values."""
    try:
        if active_job_running(st.session_state):
            return
        job, _ = resume_ffs(saved, Path(run_dir))
        st.session_state.saved_switch_ffs_job = job
        st.session_state.pop(_key(saved) + "_error", None)
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        st.session_state[_key(saved) + "_error"] = str(exc)


def _number(value, *, minimum=0.):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value >= minimum


def _read_json(path):
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _contained_movie_file(folder, value):
    """Accept only regular movie artifacts contained by this calculation."""
    from .oxdna_movie import _contained
    if not isinstance(value, str) or not value:
        return None
    root = Path(folder).resolve()
    path = Path(value)
    try:
        filename = str(path.relative_to(root)) if path.is_absolute() else value
        _contained(root, filename)
        return filename
    except (OSError, ValueError):
        return None


def _render_prepared_view(folder, request, data, key, *, visualization=None):
    """Show a validated preparation while a replica has no 3D artifact yet."""
    report_file = _contained_movie_file(folder, 'preparation/preparation.json')
    report = _read_json(Path(folder) / report_file) if report_file else {}
    prepared = data.get('preparation') or request.get('prepared') or report
    topology = _contained_movie_file(folder, prepared.get('topology_file')) or _contained_movie_file(folder, 'preparation/system.top')
    if topology is None:
        return False
    candidates = [data.get('prepared_movie_file'),
                  data.get('movie_file') if data.get('movie_kind') == 'prepared' else None,
                  prepared.get('configuration_file'), report.get('configuration_file')]
    for candidate in candidates:
        if filename := _contained_movie_file(folder, candidate):
            st.caption(tr('Configuration 3D préparée · les essais FFS de cette réplique ne sont pas encore disponibles.'))
            render_oxdna_movie({'status': 'prepared', 'movie_kind': 'prepared', 'movie_file': filename,
                               'topology_file': topology, 'visualization': visualization or data.get('visualization')},
                              folder, key + '_prepared')
            return True
    return False


def _remember_option(key, field):
    draft = dict(st.session_state.get(key + "_draft", {}))
    draft[field] = st.session_state[key + "_" + field]
    st.session_state[key + "_draft"] = draft


def _input(container, key, field, label, *, value, minimum, maximum, step, disabled, help=None):
    draft = st.session_state.get(key + "_draft", {})
    value = max(minimum, min(maximum, draft.get(field, value)))
    return container.number_input(tr(label), min_value=minimum, max_value=maximum, value=value,
                                  step=step, key=key + "_" + field, disabled=disabled,
                                  on_change=_remember_option, args=(key, field), help=help)


def _render_kinetics(data):
    kinetics = data.get("kinetics") or {}
    proximity = _simulation_mode(data) == "proximity"
    phase = kinetics.get("phase")
    phases = {"flux": "Franchissements naturels depuis les brins séparés",
              "equilibration_switch": "Équilibrage du switch isolé",
              "equilibration_trigger": "Équilibrage du trigger isolé",
              "initialization": "Préparation du départ proche sans liaison interbrin",
              "production": "Franchissements des interfaces FFS", "report": "Analyse des résultats"}
    details = [tr(phases[phase])] if phase in phases and not (proximity and phase == "flux") else []
    completed, replicas = kinetics.get("completed_replicas"), kinetics.get("replicas")
    if _number(completed) and _number(replicas, minimum=1):
        details.append(tr("{completed}/{replicas} répétitions terminées",
                          completed=int(completed), replicas=int(replicas)))
    workers = kinetics.get("workers")
    if _number(workers, minimum=1):
        details.append(tr("{count} simulations en parallèle", count=int(workers)))
    crossings, target = kinetics.get("flux_crossings"), kinetics.get("flux_target")
    if not proximity and _number(crossings) and _number(target, minimum=1):
        details.append(tr("Flux initial : {count}/{target} franchissements",
                          count=int(crossings), target=int(target)))
    if details:
        st.caption(" · ".join(details))
    if proximity:
        _render_conditional_probability(data)
        return
    rate = kinetics.get("rate") or {}
    status, value = rate.get("status"), rate.get("value_model")
    if status == "model_estimate" and _number(value) and value > 0 and rate.get("unit") == "1/oxDNA_time":
        st.metric(tr("Taux d’ouverture dans la boîte"), f"{value:.4g} / τ",
                  help=tr("Unité : 1 / temps oxDNA. Il s’agit d’un taux du modèle, pas d’une cinétique expérimentale en secondes."))
        interval = rate.get("interval95_model")
        if (isinstance(interval, (list, tuple)) and len(interval) == 2
                and all(_number(bound) for bound in interval) and interval[0] <= value <= interval[1]):
            st.caption(tr("Intervalle à 95 % approximatif entre répétitions : [{lower:.4g}, {upper:.4g}] / temps oxDNA.",
                          lower=interval[0], upper=interval[1]))
        else:
            st.caption(tr("L’incertitude du taux n’est pas encore estimable."))
        molar_rate = rate.get("value_m_inv_model_time")
        if _number(molar_rate) and molar_rate > 0:
            st.caption(tr("Équivalent pour la concentration de la boîte : {value:.4g} M⁻¹ / temps oxDNA.", value=molar_rate))
    elif status == "no_success":
        st.info(tr("Aucun succès final n’a été observé dans les essais prévus. Cela ne démontre pas que le taux est nul."))
    elif status == "incomplete":
        st.info(tr("Le calcul est incomplet ; le taux sera disponible une fois les étapes et les répétitions terminées."))


def _stage_rows(stages, replica):
    return [{tr("Répétition"): replica, tr("Étape"): stage.get("name", ""),
             tr("Essais lancés"): stage.get("attempted", 0),
             tr("Essais terminés"): stage.get("successes", 0) + stage.get("failures", 0),
             tr("Essais prévus"): stage.get("target_trials", 0), tr("Succès"): stage.get("successes", 0),
             tr("Retours à l’état séparé"): stage.get("failures", 0),
             tr("Essais en cours"): stage.get("pending", 0),
             tr("État"): tr("Terminée" if stage.get("completed") else "En cours")}
            for stage in stages or [] if isinstance(stage, dict)]


def _diagnostics(data, *, log=()):
    """Keep numerical details available without crowding the film and rate."""
    with st.expander(tr("Détails du calcul")):
        kinetics = data.get("kinetics") or {}
        repetitions = data.get("repetitions") or []
        rows = ([row for index, replica in enumerate(repetitions, 1)
                 for row in _stage_rows(replica.get("stages"), replica.get("index", index - 1) + 1)] if repetitions
                else _stage_rows(data.get("stages"), tr("Toutes les répétitions")))
        if rows:
            render_readonly_table(rows)
        observation = data.get("observation") or {}
        distance = observation.get("distance")
        if _number(distance):
            st.caption(tr("{label} : {distance:.2f} unités oxDNA (≈ {nm:.2f} nm).",
                          label=distance_label(data.get("plan")), distance=distance, nm=distance * .8518))
        separation = observation.get("separation")
        if _number(separation):
            st.caption(tr("{label} : {distance:.2f} unités oxDNA (≈ {nm:.2f} nm).",
                          label=tr("Distance minimale entre les deux brins"), distance=separation, nm=separation * .8518))
        report = ((kinetics.get("conditional_probability") or data.get("conditional_probability") or {})
                  if _simulation_mode(data) == "proximity" else kinetics.get("rate") or {})
        for limitation in report.get("limitations") or []:
            st.write(tr(limitation))
        st.json(data, expanded=False)
        if log:
            st.code("\n".join(log), language=None)


def _replica_progress(folder, request, data):
    """Read the small native progress files, never the trajectory checkpoints."""
    count = (request.get("kinetics") or {}).get("replicas") or (data.get("kinetics") or {}).get("replicas", 0)
    if not _number(count, minimum=1):
        return {}
    states = {index: state for index in range(min(100, int(count)))
              if (state := _read_json(folder / f"replica_{index + 1:03d}" / "progress.json"))}
    for position, state in enumerate(data.get("repetitions") or []):
        suffix = str(state.get("replica_id", "")).removeprefix("replica_")
        index = state.get("index", int(suffix) - 1 if suffix.isdigit() else position)
        if index not in states and 0 <= index < int(count):
            states[index] = state
    if data.get("status") in {"completed", "cancelled", "interrupted", "failed"}:
        for index, state in states.items():
            result = _read_json(folder / f"replica_{index + 1:03d}" / "result.json")
            # A resumed live run may retain an older result; use it only once
            # the master is terminal. Per-replica results are small summaries.
            for field in ("status", "error"):
                if field in result:
                    state[field] = result[field]
    return states


def _replica_dashboard(folder, request, data, key, *, visualization=None, live=False, prepared_fallback=True):
    from .ffs_progress_charts import replica_rows, selected_replica_chart

    states = _replica_progress(folder, request, data)
    # The master publishes the actual native interfaces, unlike the input plan.
    request = {**request, "plan": data.get("plan") or request.get("plan") or {},
               "kinetics": {**(data.get("kinetics") or {}), **(request.get("kinetics") or {})}}
    rows = replica_rows(states, request, status=data.get("status", "running" if live else "completed"))
    if not rows:
        if prepared_fallback:
            _render_prepared_view(folder, request, data, key, visualization=visualization)
        return states
    explanations = [tr("Une ligne = une simulation indépendante. L’avancement porte sur l’étape indiquée, pas sur la durée restante. Les paires et la distance décrivent la configuration actuelle.")]
    if description := distance_description(request["plan"]):
        explanations.append(description)
    if any(state.get("phase") == "flux" for state in states.values()):
        explanations.append(tr("Collecte initiale : configurations franchissant le premier seuil de rapprochement après un retour à l’état séparé. Ce compteur ne compte pas des paires de bases."))
    st.markdown("**" + tr("Suivi par réplique") + "**", help="\n\n".join(explanations))
    if request.get("protocol_amendments") or data.get("protocol_amendments"):
        st.info(tr("Les critères ont changé pendant ce calcul ; aucune probabilité conditionnelle n’est estimée."
                   if _simulation_mode(request) == "proximity" else
                   "Run de test : les critères ont été modifiés en cours de simulation. Les configurations conservées restent explorables ; aucun taux cinétique ne sera estimé pour ce run."))
    render_readonly_table(rows)
    labels = {index: tr("Réplique {number}", number=index + 1) for index in range(len(rows))}
    selected = st.selectbox(tr("Réplique à visualiser"), list(labels), key=key + "_replica",
                            format_func=labels.__getitem__)
    state = states.get(selected, {})
    has_view = False
    if live:
        filename = f"replica_{selected + 1:03d}/live.dat"
        if _contained_movie_file(folder, filename):
            render_oxdna_movie({"status": "running", "movie_file": filename, "movie_kind": "live",
                               "replica_number": selected + 1,
                               "topology_file": "preparation/system.top", "visualization": visualization},
                              folder, key + f"_live_{selected}")
            has_view = True
        else:
            st.caption(tr("La vue live apparaîtra lorsque cette réplique aura démarré."))
    if state:
        filename = f"replica_{selected + 1:03d}/movie.dat"
        if _contained_movie_file(folder, filename):
            metadata = _read_json(folder / f"replica_{selected + 1:03d}" / "movie_metadata.json")
            render_oxdna_movie({**metadata, "movie_file": filename,
                               "replica_number": selected + 1,
                               "topology_file": "preparation/system.top", "visualization": visualization},
                              folder, key + f"_movie_{selected}")
            _render_trajectory_observables(metadata, request, key + f"_trajectory_{selected}")
            has_view = True
    if not live and not has_view:
        filename = f"replica_{selected + 1:03d}/live.dat"
        if _contained_movie_file(folder, filename):
            st.caption(tr('Dernière configuration 3D de la réplique · cet instantané ne représente pas une trajectoire FFS réussie.'))
            render_oxdna_movie({'status': 'prepared', 'movie_kind': 'prepared', 'movie_file': filename,
                               'replica_number': selected + 1, 'final_interface_reached': False,
                               'topology_file': 'preparation/system.top', 'visualization': visualization},
                              folder, key + f'_snapshot_{selected}')
            has_view = True
    if not has_view and prepared_fallback:
        has_view = _render_prepared_view(folder, request, data, key, visualization=visualization)
    if not live and state and not has_view:
        st.caption(tr("Aucun film de trajectoire n’est disponible pour cette réplique."))
    chart = selected_replica_chart(state, request, status=data.get("status", "running" if live else "completed"))
    if chart is not None:
        st.markdown("**" + tr("Succès aux interfaces") + "**",
                    help=tr("Chaque barre compte les succès à une interface de cette réplique. Les autres essais peuvent échouer ou être encore en cours ; ce n’est pas un pourcentage d’ouverture."))
        st.altair_chart(chart, width="stretch", key=key + f"_interfaces_{selected}")
    return states


@st.fragment(run_every=1.)
def _monitor(job):
    # Fragment reruns have their own execution context; restore the UI language.
    set_language(st.session_state.get("ui_language", "en"))
    snapshot = job.snapshot()
    if snapshot.state != "running":
        if job.consume_once():
            st.rerun()
        return
    folder = job.folder
    progress = _read_json(folder / "progress.json")
    request = getattr(job, "request", None) or _read_json(folder / "worker_request.json")
    st.markdown("**" + _mode_title(_simulation_mode(request or progress)) + "**")
    kinetics = progress.get("kinetics") or {}
    completed, replicas = kinetics.get("completed_replicas"), kinetics.get("replicas")
    if _number(completed) and _number(replicas, minimum=1):
        st.caption(tr("{completed}/{replicas} répliques terminées · {workers} calculs simultanés",
                      completed=int(completed), replicas=int(replicas), workers=kinetics.get("workers", 0)))
    else:
        st.caption(tr(snapshot.message))
    st.caption(tr("Simulation en cours · {minutes:.1f} min", minutes=snapshot.elapsed_seconds / 60.))
    st.button(tr("Mettre en pause"), key=_key(job.saved) + "_pause", on_click=job.cancel,
              type="secondary", width="stretch")
    states = _replica_dashboard(folder, request, progress, _key(job.saved) + "_" + job.run_id,
                                visualization=ffs_visualization(job.saved.snapshot), live=True)
    _diagnostics({**progress, "plan": progress.get("plan") or request.get("plan") or {},
                  "repetitions": list((states or {}).values())}, log=snapshot.log)


def _result_label(result):
    labels = {"completed": "Terminée", "cancelled": "En pause", "interrupted": "Interrompue", "failed": "Échec"}
    return f"{_mode_title(_simulation_mode(result))} · {result.get('completed_at') or result.get('started_at', '')} · {tr(labels.get(result.get('status'), 'Interrompue'))}"


def _render_result(result, key):
    request = result.get("requested_parameters") or _read_json(Path(result["run_dir"]) / "worker_request.json")
    status = result.get("status")
    if result.get("error") or status == "failed":
        st.error(result.get("error") or tr("La simulation a échoué."))
    elif status in {"cancelled", "interrupted"}:
        st.info(tr("Simulation en pause ; vous pouvez reprendre le calcul enregistré.") if result.get("resumable")
                else tr("Simulation interrompue ; les résultats déjà produits sont conservés."))
    folder = Path(result["run_dir"])
    global_view = False
    if not _replica_progress(folder, request, result):
        if result.get("movie_file"):
            render_oxdna_movie(result, folder, key + "_movie")
            _render_trajectory_observables(result, request, key + "_trajectory")
            global_view = True
        else:
            st.caption(tr("Aucun film n’est encore disponible pour ce calcul."))
    result = {**result, "requested_parameters": request}
    _render_kinetics(result)
    _replica_dashboard(folder, request, result, key, visualization=result.get("visualization"),
                       prepared_fallback=not global_view)
    _diagnostics({**result, "plan": result.get("plan") or request.get("plan") or {}})
    st.download_button(tr("Télécharger le résultat de simulation"),
                       json.dumps(result, ensure_ascii=False, indent=2).encode("utf-8"),
                       "kinetics_result.json", "application/json", key=key + "_json", width="stretch")


def render_saved_switch_ffs(saved, *, thermo=None) -> None:
    key = _key(saved)
    st.subheader(tr("Cinétique · rapprochement, hybridation et ouverture"),
                 help=tr("Choisissez le départ de la simulation. Le calcul continue lorsque vous changez de page et peut être mis en pause puis repris."))
    job = active_ffs(saved)
    if job:
        st.session_state.saved_switch_ffs_job = job
    busy = active_job_running(st.session_state)
    modes = {"opening": tr("Cinétique complète"), "proximity": tr("Départ proche du toehold")}
    default_mode = st.session_state.get(key + "_draft", {}).get("kinetics_mode", "opening")
    mode = st.radio(tr("Mode de simulation"), list(modes),
                    index=list(modes).index(default_mode) if default_mode in modes else 0,
                    format_func=modes.__getitem__, key=key + "_kinetics_mode", disabled=busy,
                    on_change=_remember_option, args=(key, "kinetics_mode"))
    if mode == "proximity":
        st.caption(tr('Les deux brins passent par la même préparation 3D guidée. Chaque réplique équilibre ensuite le switch et le trigger séparément pendant 200 000 pas chacun, sans contraintes, avant leur placement proche du toehold. Aucun flux initial n’est mesuré.'))
        st.caption(tr("Le trigger ADN démarre près du toehold, à 1,7–2,1 nm, sans aucune liaison interbrin. Chaque réplique utilise une position et une orientation différentes. Les reculs et les dissociations restent possibles pendant la simulation."))
        st.caption(tr("Ce mode estime la probabilité d’ouverture depuis ce départ proche, sans attendre la collecte des franchissements initiaux. Il ne mesure pas le taux global d’association depuis des brins libres."))
    else:
        st.caption(tr("Chaque réplique démarre avec les deux brins séparés, à des positions et orientations aléatoires dans le solvant. Après équilibrage, la FFS mesure le flux initial, conserve les franchissements de seuils et poursuit les essais jusqu’à l’hybridation complète et l’ouverture."))
    if busy and not job:
        st.info(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
    conditions = saved.snapshot["request"]["thermo"]
    limits = {"temperature_min_c": 25., "temperature_max_c": 55., "sodium_min_m": .12, "sodium_max_m": 1.}
    temperature_bounds = max(0., float(limits["temperature_min_c"])), min(90., float(limits["temperature_max_c"]))
    salt_bounds = max(.01, float(limits["sodium_min_m"])), min(2., float(limits["sodium_max_m"]))
    temperature_value = float(getattr(thermo, "temperature_c", conditions["temperature_c"]))
    salt_value = float(getattr(thermo, "sodium_m", conditions["sodium_m"]))
    temperature_default = max(temperature_bounds[0], min(temperature_bounds[1], temperature_value))
    salt_default = max(salt_bounds[0], min(salt_bounds[1], salt_value))
    columns = st.columns(2)
    temperature = _input(columns[0], key, "temperature_c", "Température de simulation (°C)",
                         minimum=temperature_bounds[0], maximum=temperature_bounds[1],
                         value=temperature_default, step=1., disabled=busy)
    draft = dict(st.session_state.get(key + "_draft", {}))
    salt_initial = max(salt_bounds[0], min(salt_bounds[1], draft.get("salt_m", salt_default)))
    with columns[1]:
        salt = render_concentration_control(
            "Sel monovalent Na⁺", value_key=key + "_salt_value", unit_key=key + "_salt_unit",
            min_m=salt_bounds[0], max_m=salt_bounds[1], value_m=salt_initial,
            default_unit="mM", disabled=busy,
        )
    # Widget values are display units; the persisted draft and oxDNA options
    # always retain the physical mol/L value under salt_m.
    draft["salt_m"] = salt
    st.session_state[key + "_draft"] = draft
    if temperature_default != temperature_value or salt_default != salt_value:
        st.info(tr("Les conditions proposées sont ramenées aux limites du modèle."))
    with st.expander(tr("Paramètres avancés")):
        columns = st.columns(2)
        replicas = _input(columns[0], key, "kinetics_replicas", "Répliques indépendantes",
                          minimum=4, maximum=100, value=10, step=1, disabled=busy)
        capacity = min(64, max(1, os.cpu_count() or 1))
        workers = _input(columns[1], key, "kinetics_workers", "Calculs simultanés",
                         minimum=1, maximum=min(capacity, int(replicas)),
                         value=min(10, capacity), step=1, disabled=busy,
                         help=tr("{cores} processeurs logiques détectés. {replicas} répliques seront calculées, dont au maximum {workers} en même temps. Plus de répliques ajoute de l’échantillonnage et du travail ; plus de calculs simultanés partage les ressources du processeur.",
                                 cores=capacity, replicas=int(replicas), workers=min(capacity, int(replicas))))
        if mode == "opening":
            source_samples = _input(st, key, "kinetics_source_samples", "Franchissements initiaux par répétition",
                                    minimum=FFS_SOURCE_SAMPLES_MIN, maximum=FFS_SOURCE_SAMPLES_MAX,
                                    value=FFS_SOURCE_SAMPLES_DEFAULT, step=1, disabled=busy,
                                    help=tr("Les premiers contacts sont collectés pendant la dynamique libre. Chaque nouveau franchissement est compté après un retour à l’état séparé ; la durée simulée permet de mesurer le flux initial."))
        trials = _input(st, key, "kinetics_trials", "Essais de production par interface",
                        minimum=FFS_TRIALS_MIN, maximum=FFS_TRIALS_MAX,
                        value=FFS_TRIALS_DEFAULT, step=1, disabled=busy)
    with st.expander(tr("Modèle et critères")):
        if mode == "opening":
            st.write(tr("Chaque réplique démarre avec les deux brins séparés, à des positions et orientations aléatoires dans le solvant. Après équilibrage, la FFS mesure le flux initial, conserve les franchissements de seuils et poursuit les essais jusqu’à l’hybridation complète et l’ouverture."))
        st.write(distance_description({"distance_scope": "trigger_to_toehold" if mode == "proximity" else "trigger_3prime_to_toehold_start"}))
        st.write(tr("Le modèle oxDNA/NA simule le fragment de switch ARN et son trigger ADN de 30 nucléotides, sans le rapporteur complet. Le sel monovalent est pris en compte ; le magnésium ne l’est pas."))
        st.write(tr("L’état final exige les 30 appariements cibles du trigger, la tige OFF entièrement libérée et aucune liaison hydrogène dans la région RBS → AUG."))
        if mode == "proximity":
            st.write(tr("Le moteur vérifie l’absence de liaison interbrin au départ. La probabilité dépend de cet ensemble de configurations proches préparées ; il ne représente pas les rencontres naturelles depuis les brins libres."))
        else:
            st.write(tr("Le taux combine le flux initial et les probabilités de franchissement jusqu’à l’ouverture, avec les mêmes conditions et configurations de départ."))
        st.write(tr("La FFS conserve les franchissements des seuils pour relancer les essais suivants. Dans chaque essai, les brins peuvent reculer ; les dissociations sont comptées comme des échecs."))
        st.write(tr("Les essais prévus sont comptés jusqu’au bout ; aucun plafond de temps n’interrompt le calcul."))
        st.write(tr("Les interfaces et les conditions d’arrêt sont évaluées par la fonction FFS_MD d’oxDNA." if mode == "proximity"
                    else "Les interfaces et les conditions d’arrêt sont évaluées par la fonction FFS_MD d’oxDNA. Les taux restent exprimés dans le temps du modèle."))
    try:
        runtime = detect_oxdna_runtime()
    except (OSError, ValueError, RuntimeError):
        runtime = None
    native_ready = bool(native_backend_identity(configured_library_path()))
    if not runtime or not native_ready:
        st.info(tr("Les moteurs de simulation intégrés sont indisponibles. Réinstallez le paquet complet de l’application."))
    supported = (conditions["switch_material"] == "rna" and conditions["trigger_material"] == "dna"
                 and len(saved.snapshot["candidate"]["trigger_seq"]) == 30)
    if not supported:
        st.info(tr("Ce calcul nécessite un switch ARN et un trigger ADN de 30 nucléotides."))
    options = {"kinetics_mode": mode, "kinetics_trials": int(trials),
               "kinetics_replicas": int(replicas), "kinetics_workers": int(workers),
               "temperature_c": float(temperature), "salt_m": float(salt)}
    if mode == "opening":
        options["kinetics_source_samples"] = int(source_samples)
    st.button(tr("Lancer depuis le toehold" if mode == "proximity" else "Lancer la simulation cinétique"), key=key + "_start", type="primary",
              disabled=not runtime or not native_ready or busy or not supported,
              on_click=_start, args=(saved, runtime, options), width="stretch")
    if st.session_state.get(key + "_error"):
        st.error(st.session_state[key + "_error"])
    if job:
        _monitor(job)
    results = saved_ffs_results(saved)
    if not results:
        return
    if job and not st.toggle(tr("Afficher les simulations précédentes"), key=key + "_history"):
        return
    if job:
        st.caption(tr("Historique : ces résultats ne décrivent pas le calcul en cours."))
    choices = list(range(len(results)))
    if st.session_state.get(key + "_result") not in choices:
        st.session_state[key + "_result"] = 0
    labels = [_result_label(result) for result in results]
    selection = st.selectbox(tr("Simulation enregistrée"), choices, key=key + "_result", format_func=labels.__getitem__)
    result = {**results[selection], "visualization": ffs_visualization(saved.snapshot)}
    if result.get("resumable"):
        st.button(tr("Reprendre la simulation"), key=key + "_resume_" + result["run_id"],
                  disabled=busy, on_click=_resume, args=(saved, result["run_dir"]), type="primary", width="stretch")
    _render_result(result, key + "_" + result["run_id"])


def _simulation_mode(data):
    request = data.get("requested_parameters") or data
    return ((request.get("kinetics") or {}).get("mode")
            or (request.get("options") or {}).get("kinetics_mode")
            or data.get("simulation_mode") or "opening")


def _mode_title(mode):
    return tr("Départ proche du toehold" if mode == "proximity" else "Simulation cinétique complète")


def _render_conditional_probability(data):
    probability = (data.get("kinetics") or {}).get("conditional_probability") or data.get("conditional_probability") or {}
    status, value = probability.get("status"), probability.get("value")
    st.caption(tr("Probabilité conditionnelle d’atteindre l’ouverture depuis les configurations proches préparées. Le temps de rencontre des brins n’est pas mesuré ; aucun taux global d’association n’est calculé."))
    if status == "conditional_estimate" and _number(value) and value <= 1 and probability.get("unit") == "probability":
        st.metric(tr("Probabilité d’ouverture depuis le départ proche"), f"{value:.3%}")
        interval = probability.get("interval95")
        if (isinstance(interval, (list, tuple)) and len(interval) == 2
                and all(_number(bound) and bound <= 1 for bound in interval)
                and interval[0] <= value <= interval[1]):
            st.caption(tr("Intervalle à 95 % approximatif entre répétitions : [{lower:.3%}, {upper:.3%}].",
                          lower=interval[0], upper=interval[1]))
        else:
            st.caption(tr("L’incertitude de la probabilité n’est pas encore estimable."))
    elif status == "no_success":
        st.info(tr("Aucune ouverture n’a été observée dans les essais prévus depuis le départ proche. Cela ne démontre pas que sa probabilité est nulle."))
    elif status == "incomplete":
        st.info(tr("Le calcul est incomplet ; la probabilité conditionnelle sera disponible une fois les étapes et les répétitions terminées."))
    elif status == "exploratory_modified_protocol":
        st.info(tr("Les critères ont changé pendant ce calcul ; aucune probabilité conditionnelle n’est estimée."))


def _render_trajectory_observables(data, request, key):
    if _simulation_mode(request or data) != "proximity" or not data.get("final_interface_reached"):
        return
    from .ffs_progress_charts import trajectory_observables_chart

    chart = trajectory_observables_chart(data.get("trajectory_observables"),
                                         request.get("sampling") or data.get("sampling"))
    if chart is not None:
        st.markdown("**" + tr("Appariements le long du film") + "**")
        st.altair_chart(chart, width="stretch", key=key)
        st.caption(tr("Trajectoire FFS réussie assemblée à partir de ses segments successifs. L’axe suit le temps du modèle le long de ce trajet, sans le temps de rencontre des brins ni les essais échoués ; il ne représente pas une cinétique expérimentale en microsecondes."))
