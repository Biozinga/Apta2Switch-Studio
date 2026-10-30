"""Saved DNA-DNA free-MD comparisons and their measured steric proxies."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import altair as alt
import streamlit as st

from aptaswitch_core.localization import set_language, tr
from aptaswitch_core.thermodynamics.setup import configured_library_path
from .scientific_engines_ui import native_backend_identity, NATIVE_UNAVAILABLE
from aptaswitch_core.oxdna_setup import detect_oxdna_runtime
from aptaswitch_core.oxdna_free_md import (
    available_free_md, free_md_sequences, make_free_md_archive, parse_monitored_bases,
)
from .design_session import active_job_running
from .saved_switch_free_md_runtime import (
    active_free_md, resume_free_md, saved_free_md_results, start_free_md,
)
from .oxdna_movie import render_oxdna_movie

_TIME_HELP = ("oxDNA utilise sa propre échelle de temps. La conversion nominale τ ≈ 3,031 ps "
              "donne une durée indicative en µs ; elle n’est pas une calibration de la cinétique expérimentale. "
              "La durée réelle du calcul dépend du nombre de nucléotides et de votre processeur.")
_VARIANTS = {"original": "Trigger original", "extended": "Trigger étendu"}
_METRICS = {
    "minimum_trigger_site_distance_nm": ("Distance minimale trigger–zone", "nm"),
    "site_contact_fraction": ("Fraction de bases proches du trigger", "fraction"),
    "site_paired_fraction": ("Fraction de bases appariées dans la zone", "fraction"),
    "site_trigger_paired_fraction": ("Fraction de bases appariées au trigger", "fraction"),
    "interstrand_pairs": ("Paires aptamère–trigger", "paires"),
    "rg_aptamer_nm": ("Rayon de giration de l’aptamère", "nm"),
    "rmsd_aptamer_nm": ("RMSD de l’aptamère", "nm"),
    "rmsd_site_nm": ("RMSD de la zone surveillée", "nm"),
}


def _key(saved):
    return "saved_switch_free_md_" + hashlib.sha256(str(saved.path.resolve()).encode()).hexdigest()[:20]


def _read_json(path):
    try:
        if path.is_symlink() or path.stat().st_size > 32 * 1024 * 1024:
            return {}
        data = json.loads(path.read_text())
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _finite(value):
    try:
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
    except OverflowError:
        return False


def _run_file(folder, filename, *, maximum_bytes=32 * 1024 * 1024):
    """Accept only regular files inside this run, without following symlinks."""
    if not isinstance(filename, (str, Path)) or not str(filename):
        return None
    try:
        folder = Path(folder)
        if folder.is_symlink():
            return None
        root = folder.resolve()
        candidate = Path(filename)
        if ".." in candidate.parts:
            return None
        candidate = candidate if candidate.is_absolute() else root / candidate
        relative = candidate.relative_to(root)
        current = root
        for part in relative.parts:
            current = current / part
            if current.is_symlink():
                return None
        if (not current.resolve().is_relative_to(root) or not current.is_file()
                or current.stat().st_size > maximum_bytes):
            return None
        return current
    except (OSError, ValueError, TypeError):
        return None


def _run_json(folder, filename):
    path = _run_file(folder, filename)
    return _read_json(path) if path else {}


def _nonnegative_int(value, default=0):
    return max(0, int(value)) if _finite(value) else default


def _production_progress(progress, request):
    """Production steps are comparable across live and reattached jobs."""
    sampling = request.get("sampling")
    sampling = sampling if isinstance(sampling, dict) else {}
    replicas = min(100, _nonnegative_int(request.get("replicas", progress.get("replicas", 0))))
    per_replica = _nonnegative_int(sampling.get("production_steps"))
    total = _nonnegative_int(progress.get("target_steps")) or per_replica * replicas
    steps = min(total, _nonnegative_int(progress.get("production_steps")))
    return {"steps": steps, "total": total, "fraction": steps / total if total else 0.,
            "completed": min(replicas, _nonnegative_int(progress.get("completed_replicas"))),
            "replicas": replicas, "per_replica": per_replica}


def _replica_progress(folder, request, data):
    """Read bounded state files, including replicas that have not started yet."""
    count = _production_progress(data, request)["replicas"]
    compact = data.get("replica_states") or []
    states = {row["replica_index"]: row for row in compact if isinstance(row, dict)
              and type(row.get("replica_index")) is int and 0 <= row["replica_index"] < count} if isinstance(compact, list) else {}
    rows = []
    for index in range(count):
        variant = "original" if index % 2 == 0 else "extended"
        state = states.get(index, {})
        if state.get("variant") not in (None, variant):
            state = {}
        committed = _run_json(folder, f"replica_{index:03d}/state.json")
        if (type(committed.get("replica_index")) is int and committed["replica_index"] == index
                and committed.get("variant") == variant):
            state = {**state, **committed}
        status = state.get("status")
        rows.append({**state, "replica_index": index, "variant": variant,
                     "status": status if isinstance(status, str) and status else "pending"})
    return rows


_PHASES = {
    "minimize_guided": "Préparation : minimisation guidée",
    "relax_guided": "Préparation : relaxation guidée",
    "guide_physical": "Préparation : stabilisation du complexe",
    "release": "Préparation : retrait des contraintes",
    "equilibration": "Préparation : équilibration libre",
    "production": "Production MD libre",
    "completed": "Terminée",
}
_STATUSES = {"pending": "En attente", "running": "En cours", "completed": "Terminée",
             "cancelled": "En pause", "interrupted": "Interrompue", "failed": "Échec"}


def _phase_label(row):
    phase = row.get("phase")
    phase = phase if isinstance(phase, str) else ""
    if row.get("status") == "completed":
        return tr("Terminée")
    return tr(_PHASES.get(phase, "Préparation" if row.get("status") == "running" else "En attente"))


def _render_live_replicas(data, folder, request, key):
    replicas = _replica_progress(folder, request, data)
    if not replicas:
        return
    target = _production_progress(data, request)["per_replica"]
    rows = []
    for replica in replicas:
        steps = min(target, _nonnegative_int(replica.get("production_steps")))
        rows.append({tr("Réplique"): replica["replica_index"] + 1,
                     tr("Groupe"): tr(_VARIANTS[replica["variant"]]),
                     tr("État"): tr(_STATUSES.get(replica["status"], "Interrompue")),
                     tr("Phase"): _phase_label(replica),
                     tr("Production (%)"): round(100 * steps / target, 2) if target else 0.,
                     tr("Pas MD de production"): f"{steps} / {target}" if target else "—",
                     tr("Pas de la phase"): _nonnegative_int(replica.get("phase_steps")),
                     tr("Observations"): _nonnegative_int(replica.get("samples_count", replica.get("sample_count")))})
    st.markdown("**" + tr("Suivi par réplique") + "**")
    st.caption(tr("Le pourcentage correspond aux pas de production réalisés. La préparation et l’équilibration sont affichées séparément et ne comptent pas dans ce pourcentage."))
    st.dataframe(rows, width="stretch", hide_index=True)
    labels = {row["replica_index"]: tr("Réplique {number}", number=row["replica_index"] + 1)
              + " · " + tr(_VARIANTS[row["variant"]]) for row in replicas}
    selection_key = key + "_replica"
    if st.session_state.get(selection_key) not in labels:
        st.session_state[selection_key] = next(iter(labels))
    selected = st.selectbox(tr("Réplique à visualiser"), list(labels), key=selection_key,
                            format_func=labels.__getitem__)
    state = replicas[selected]
    steps = min(target, _nonnegative_int(state.get("production_steps")))
    fraction = steps / target if target else 0.
    st.progress(fraction, text=tr("{percent:.2f} % de production · {steps} / {total} pas MD",
                                percent=fraction * 100, steps=steps, total=target))
    st.caption(_phase_label(state))
    checkpoint = state.get("checkpoint_file")
    path = (_run_file(folder, f"replica_{selected:03d}/{checkpoint}")
            if checkpoint in ("checkpoint_0.dat", "checkpoint_1.dat") else None)
    variants = request.get("variants")
    variant = next((row for row in variants if isinstance(row, dict)
                    and row.get("name") == state["variant"]), {}) if isinstance(variants, list) else {}
    topology = _run_file(folder, variant.get("topology_file") or state.get("topology_file"))
    if not path or not topology:
        st.caption(tr("La vue live apparaîtra après la première sauvegarde de cette réplique. Vous pouvez déjà sélectionner une autre réplique."))
        return replicas
    visualization = variant.get("visualization")
    if not isinstance(visualization, dict):
        from aptaswitch_core.oxdna_free_md_metrics import _visualization
        sequences = (variant.get("aptamer_sequence"), variant.get("trigger_sequence"))
        bases = request.get("monitored_bases")
        if (all(isinstance(sequence, str) and sequence for sequence in sequences)
                and sum(map(len, sequences)) <= 5000 and isinstance(bases, list)):
            visualization = _visualization(variant, bases)
        else:
            visualization = None
    root = Path(folder).resolve()
    render_oxdna_movie({"status": state["status"], "movie_kind": "live", "trajectory_endpoint": "free_md",
                       "movie_file": path.relative_to(root).as_posix(),
                       "topology_file": topology.relative_to(root).as_posix(),
                       "replica_number": selected + 1, "visualization": visualization},
                      root, key + f"_live_{selected}")
    return replicas


def _remember(key, field):
    st.session_state[key + "_draft"] = {
        **st.session_state.get(key + "_draft", {}), field: st.session_state[key + "_" + field]}


def _number(container, key, field, label, *, value, minimum, maximum, step, disabled, **extra):
    value = max(minimum, min(maximum, st.session_state.get(key + "_draft", {}).get(field, value)))
    return container.number_input(tr(label), value=value, min_value=minimum, max_value=maximum,
                                  step=step, disabled=disabled, key=key + "_" + field,
                                  on_change=_remember, args=(key, field), **extra)


def _start(saved, runtime, options):
    try:
        if active_job_running(st.session_state):
            return
        job, _ = start_free_md(saved, runtime, {**options, "thermodynamics_library_path": configured_library_path() or ""})
        st.session_state.saved_switch_free_md_job = job
        st.session_state.pop(_key(saved) + "_error", None)
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        st.session_state[_key(saved) + "_error"] = str(exc)


def _resume(saved, folder):
    try:
        if active_job_running(st.session_state):
            return
        job, _ = resume_free_md(saved, Path(folder))
        st.session_state.saved_switch_free_md_job = job
        st.session_state.pop(_key(saved) + "_error", None)
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        st.session_state[_key(saved) + "_error"] = str(exc)


def comparison_chart(comparison, metric):
    """Between-replica spread, never a confidence interval from correlated frames."""
    rows = []
    for row in comparison.get("curves", []):
        if (row.get("metric") != metric or row.get("variant") not in _VARIANTS
                or not all(_finite(row.get(name)) for name in ("time_us", "mean", "n"))):
            continue
        sd = row.get("sd") if _finite(row.get("sd")) and row["sd"] >= 0 else 0.
        rows.append({"time_us": row["time_us"], "mean": row["mean"], "sd": sd,
                     "lower": row["mean"] - sd, "upper": row["mean"] + sd,
                     "n": row["n"], "variant": tr(_VARIANTS[row["variant"]])})
    if not rows:
        return None
    label, unit = _METRICS[metric]
    ytitle = tr(label) + (" (" + tr(unit) + ")" if unit else "")
    color = alt.Color("variant:N", title=tr("Groupe"),
                      scale=alt.Scale(domain=[tr(_VARIANTS[name]) for name in _VARIANTS],
                                      range=["#167D9A", "#D47E1E"]))
    base = alt.Chart(alt.Data(values=rows)).encode(
        x=alt.X("time_us:Q", title=tr("Temps approximatif du modèle (µs)")), color=color)
    band = base.mark_area(opacity=.14).encode(y=alt.Y("lower:Q", title=ytitle), y2="upper:Q")
    line = base.mark_line(strokeWidth=2.5).encode(y=alt.Y("mean:Q", title=ytitle), tooltip=[
        alt.Tooltip("variant:N", title=tr("Groupe")),
        alt.Tooltip("time_us:Q", title=tr("Temps approximatif (µs)"), format=".4g"),
        alt.Tooltip("mean:Q", title=tr("Moyenne"), format=".4g"),
        alt.Tooltip("sd:Q", title=tr("Écart-type entre répliques"), format=".4g"),
        alt.Tooltip("n:Q", title=tr("Répliques observées")),
    ])
    return (band + line).properties(height=300).interactive()


def site_chart(comparison, metric):
    rows = [{"base": row["base"], "value": row[metric], "variant": tr(_VARIANTS[variant])}
            for variant, group in comparison.get("groups", {}).items() if variant in _VARIANTS
            for row in group.get("site", []) if _finite(row.get("base")) and _finite(row.get(metric))]
    if not rows:
        return None
    title = tr("Fluctuation RMSF (nm)" if metric == "rmsf_nm" else "Fréquence de proximité du trigger")
    return alt.Chart(alt.Data(values=rows)).mark_line(point=True).encode(
        x=alt.X("base:O", title=tr("Base de l’aptamère (1 = extrémité 5′)")),
        y=alt.Y("value:Q", title=title, scale=alt.Scale(zero=True)),
        color=alt.Color("variant:N", title=tr("Groupe"),
                        scale=alt.Scale(domain=[tr(value) for value in _VARIANTS.values()],
                                        range=["#167D9A", "#D47E1E"])),
        tooltip=[alt.Tooltip("variant:N", title=tr("Groupe")),
                 alt.Tooltip("base:O", title=tr("Base")),
                 alt.Tooltip("value:Q", title=title, format=".4g")],
    ).properties(height=230)


def _render_comparison(comparison, key):
    groups = comparison.get("groups") or {}
    if not groups:
        return
    st.markdown("**" + tr("Comparaison original / étendu") + "**")
    st.dataframe([{tr("Groupe"): tr(_VARIANTS[variant]),
                   tr("Répliques terminées"): group.get("completed", 0),
                   tr("Répliques avec données"): group.get("available", 0),
                   tr("Répliques prévues"): comparison.get("expected_per_group", 0)}
                  for variant, group in groups.items() if variant in _VARIANTS],
                 hide_index=True, width="stretch")
    if comparison.get("status") != "complete":
        st.info(tr("Comparaison partielle : les groupes peuvent couvrir des durées et des nombres de répliques différents."))
    metrics = [metric for metric in _METRICS if any(row.get("metric") == metric for row in comparison.get("curves", []))]
    if metrics:
        metric_labels = {name: tr(_METRICS[name][0]) for name in metrics}
        metric = st.selectbox(tr("Mesure à comparer"), metrics, key=key + "_metric",
                              format_func=metric_labels.__getitem__)
        chart = comparison_chart(comparison, metric)
        if chart is not None:
            st.altair_chart(chart, width="stretch", key=key + "_comparison")
        st.caption(tr("Courbes : moyenne des répliques à chaque temps ; bande : ± un écart-type entre répliques. Ce n’est pas un intervalle de confiance. Le nombre de répliques observées est indiqué au survol."))
    for metric in ("contact_fraction", "rmsf_nm"):
        chart = site_chart(comparison, metric)
        if chart is not None:
            st.altair_chart(chart, width="stretch", key=key + "_" + metric)
    st.caption(tr("Les fréquences portent sur les trames enregistrées. Les RMSD et RMSF sont calculés après alignement de l’aptamère sur sa première trame de production ; ils excluent sa translation et sa rotation globales."))
    with st.expander(tr("Valeurs résumées par groupe")):
        rows = []
        for variant, group in groups.items():
            for metric, stats in group.get("summary", {}).items():
                if variant not in _VARIANTS or metric not in _METRICS:
                    continue
                rows.append({tr("Groupe"): tr(_VARIANTS[variant]), tr("Mesure"): tr(_METRICS[metric][0]),
                             tr("Unité"): tr(_METRICS[metric][1]), tr("Moyenne"): stats.get("mean"),
                             tr("Écart-type entre répliques"): stats.get("sd"), tr("Répliques observées"): stats.get("n")})
        if rows:
            st.dataframe(rows, width="stretch", hide_index=True)


def _render_replicas(data, folder, key):
    replicas = data.get("replica_states") or data.get("replicas") or []
    if not isinstance(replicas, list) or not replicas:
        return
    replicas = [row for row in replicas if isinstance(row, dict) and row.get("replica_index") is not None]
    rows = [{tr("Réplique"): row.get("replica_index", index) + 1,
             tr("Groupe"): tr(_VARIANTS.get(row.get("variant"), row.get("variant", ""))),
             tr("État"): tr(row.get("status", "")), tr("Pas MD"): row.get("production_steps", 0),
             tr("Observations"): row.get("sample_count", row.get("samples_count", 0))} for index, row in enumerate(replicas)]
    with st.expander(tr("Suivi des répliques")):
        st.dataframe(rows, width="stretch", hide_index=True)
    movies = [row for row in replicas if row.get("trajectory_file") and row.get("topology_file")]
    if not movies:
        return
    labels = [tr("Réplique {number}", number=row.get("replica_index", index) + 1)
              + " · " + tr(_VARIANTS.get(row.get("variant"), "")) for index, row in enumerate(movies)]
    selected = st.selectbox(tr("Réplique à visualiser"), list(range(len(movies))), key=key + "_replica",
                           format_func=labels.__getitem__)
    row = movies[selected]
    render_oxdna_movie({**row, "trajectory_endpoint": "free_md", "movie_file": row["trajectory_file"]},
                      folder, key + "_movie_" + str(selected))


def _archive_revision(result):
    """A resumed run keeps its ID but produces a new set of raw artifacts."""
    revision = {name: result.get(name) for name in
                ("status", "completed_at", "resumed_at", "elapsed_seconds", "sampling_complete")}
    replicas = result.get("replica_states") or result.get("replicas") or []
    if not isinstance(replicas, (list, tuple)):
        replicas = []
    revision["replicas"] = [
        {name: row.get(name) for name in ("replica_index", "status", "production_steps", "sample_count", "samples_count")}
        for row in replicas if isinstance(row, dict)
    ]
    return hashlib.sha256(json.dumps(revision, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _archive(folder, key, revision):
    try:
        st.session_state[key + "_archive"] = make_free_md_archive(Path(folder))
        st.session_state[key + "_archive_revision"] = revision
        st.session_state.pop(key + "_archive_error", None)
    except (OSError, ValueError) as exc:
        st.session_state[key + "_archive_error"] = str(exc)


def _render_result(result, key):
    revision = _archive_revision(result)
    if st.session_state.get(key + "_archive_revision") != revision:
        st.session_state.pop(key + "_archive", None)
        st.session_state.pop(key + "_archive_error", None)
        st.session_state[key + "_archive_revision"] = revision
    if result.get("error"):
        st.error(tr(result["error"]))
    if result.get("status") in {"cancelled", "interrupted"}:
        st.info(tr("Simulation en pause ; vous pouvez reprendre le calcul enregistré.") if result.get("resumable")
                else tr("Simulation interrompue ; les résultats déjà produits sont conservés."))
    st.caption(tr(_TIME_HELP))
    request = result.get("request") or result.get("requested_parameters") or {}
    sampling = request.get("sampling") or {}
    if _finite(sampling.get("duration_us")) and request.get("monitored_bases"):
        st.caption(tr("Paramètres enregistrés : {duration:g} µs approximatifs par réplique · bases {bases} · seuil de proximité {cutoff:g} nm.",
                      duration=sampling["duration_us"], bases=", ".join(map(str, request["monitored_bases"])),
                      cutoff=request.get("contact_cutoff_nm", 1.5)))
    _render_comparison(result.get("comparison") or {}, key)
    _render_replicas(result, Path(result["run_dir"]), key)
    st.download_button(tr("Télécharger le résultat de simulation"), json.dumps(result, ensure_ascii=False, indent=2),
                       "oxDNA_free_MD_result.json", "application/json", key=key + "_json")
    st.button(tr("Préparer l’archive des données brutes"), key=key + "_make_archive",
              on_click=_archive, args=(result["run_dir"], key, revision))
    if key + "_archive" in st.session_state:
        st.download_button(tr("Télécharger toutes les données brutes (ZIP)"), st.session_state[key + "_archive"],
                           "oxDNA_free_MD_raw.zip", "application/zip", key=key + "_zip")
    if st.session_state.get(key + "_archive_error"):
        st.error(st.session_state[key + "_archive_error"])


@st.fragment(run_every=2)
def _monitor(job):
    # Fragment reruns have their own context, including the selected language.
    set_language(st.session_state.get("ui_language", "en"))
    snapshot = job.snapshot()
    progress = _run_json(job.folder, "progress.json")
    request = _run_json(job.folder, "worker_request.json") or getattr(job, "request", None) or {}
    if not isinstance(request, dict):
        request = {}
    values = _production_progress(progress, request)
    st.markdown("**" + tr("MD libre ADN–ADN en cours") + "**")
    st.progress(values["fraction"], text=tr("oxDNA DNA2 · {percent:.2f} % · {completed}/{replicas} réplicats terminés",
                percent=values["fraction"] * 100, completed=values["completed"], replicas=values["replicas"]))
    st.caption(tr("{steps} / {total} pas MD de production", steps=values["steps"], total=values["total"]))
    if snapshot.state == "running":
        st.button(tr("Mettre en pause"), key="free_md_pause_" + job.run_id, on_click=job.cancel)
        if not progress or progress.get("status") == "cancelling":
            st.caption(tr(snapshot.message))
        _render_live_replicas(progress, job.folder, request, "free_md_live_" + job.run_id)
    elif job.consume_once():
        st.rerun()


def render_saved_switch_free_md(saved, *, thermo=None):
    if not available_free_md(saved.snapshot):
        return
    sequences = free_md_sequences(saved.snapshot)
    key = _key(saved)
    st.subheader(tr("MD libre ADN–ADN · zone de liaison du ligand"))
    st.caption(tr("Comparez l’aptamère avec son trigger original et avec le trigger étendu enregistré. Chaque groupe reçoit la moitié des répliques, avec la même durée et les mêmes conditions."))
    st.caption(tr("Les complexes ADN–ADN sont préparés, puis évoluent librement avec oxDNA2. Le ligand n’est pas simulé : les distances, contacts et fluctuations mesurent des indices d’encombrement, pas un docking ni une affinité."))
    job = active_free_md(saved)
    if job:
        st.session_state.saved_switch_free_md_job = job
    busy = active_job_running(st.session_state)
    if busy:
        st.info(tr("Un calcul est en cours. Mettez la simulation en pause avant d’en lancer une autre."))
    with st.expander(tr("Séquences comparées")):
        for field, label in (("aptamer", "Aptamère"), ("original", "Trigger original"), ("extended", "Trigger étendu")):
            st.write(tr(label))
            st.code(sequences[field], language=None, wrap_lines=True)
    draft = st.session_state.get(key + "_draft", {})
    bases = st.text_input(tr("Bases de l’aptamère à surveiller"), value=draft.get("monitored_text", ""),
                         help=tr("Positions à partir de 1, de 5′ vers 3′. Exemple : 4-8, 12, 18-20."),
                         key=key + "_monitored_text", disabled=busy, on_change=_remember, args=(key, "monitored_text"))
    selected, valid = [], False
    if bases.strip():
        try:
            selected = parse_monitored_bases(bases, len(sequences["aptamer"]))
            valid = True
            st.caption(tr("{count} bases surveillées sur {length} nucléotides.", count=len(selected), length=len(sequences["aptamer"])))
        except ValueError as exc:
            st.error(tr(str(exc)))
    columns = st.columns(2)
    duration = _number(columns[0], key, "duration_us", "Durée approximative par réplique (µs)", value=.1,
                       minimum=.001, maximum=100., step=.01, format="%.3f", disabled=busy, help=tr(_TIME_HELP))
    replicas = _number(columns[1], key, "replicas", "Nombre total de répliques (pair)", value=10,
                       minimum=2, maximum=100, step=2, disabled=busy)
    even = int(replicas) % 2 == 0
    if not even:
        st.error(tr("Choisissez un nombre pair de répliques pour répartir les deux groupes à 50 / 50."))
    st.caption(tr("{count} répliques avec le trigger original et {count} avec le trigger étendu.", count=int(replicas) // 2))
    conditions = saved.snapshot["request"]["thermo"]
    with st.expander(tr("Paramètres avancés")):
        columns = st.columns(2)
        temperature = _number(columns[0], key, "temperature_c", "Température de simulation (°C)", value=float(getattr(thermo, "temperature_c", conditions.get("temperature_c", 37.))), minimum=25., maximum=55., step=1., disabled=busy)
        salt = _number(columns[1], key, "salt_m", "Sel monovalent Na⁺ (M)", value=float(getattr(thermo, "sodium_m", conditions.get("sodium_m", .15))), minimum=.12, maximum=1., step=.01, format="%.3f", disabled=busy)
        cutoff = _number(st, key, "contact_cutoff_nm", "Seuil de proximité exploratoire (nm)", value=1.5, minimum=.1, maximum=10., step=.1, disabled=busy,
                         help=tr("Ce seuil géométrique compare les deux groupes ; ce n’est pas une limite de contact du ligand validée expérimentalement."))
        workers = _number(st, key, "workers", "Calculs simultanés", value=min(4, os.cpu_count() or 1), minimum=1, maximum=min(int(replicas), os.cpu_count() or 1), step=1, disabled=busy)
    runtime = detect_oxdna_runtime()
    native_ready = bool(native_backend_identity(configured_library_path()))
    ready = bool(runtime and native_ready)
    if not runtime:
        st.info(tr("Configurez oxDNA dans Réglages → Dépendances pour lancer le calcul."))
    if not native_ready:
        st.info(tr(NATIVE_UNAVAILABLE))
    options = {"monitored_bases": selected, "duration_us": float(duration), "replicas": int(replicas),
               "workers": int(workers), "temperature_c": float(temperature), "salt_m": float(salt),
               "contact_cutoff_nm": float(cutoff)}
    st.button(tr("Lancer la MD libre ADN–ADN"), key=key + "_start", type="primary",
              disabled=busy or not ready or not valid or not even, on_click=_start, args=(saved, runtime, options))
    if st.session_state.get(key + "_error"):
        st.error(st.session_state[key + "_error"])
    if job:
        _monitor(job)
    results = saved_free_md_results(saved)
    if not results or (job and not st.toggle(tr("Afficher les simulations précédentes"), key=key + "_history")):
        return
    labels = {"completed": "Terminée", "cancelled": "En pause", "interrupted": "Interrompue", "failed": "Échec"}
    result_labels = [str(row.get("started_at", "")) + " · " + tr(labels.get(row.get("status"), "Interrompue"))
                     for row in results]
    selection = st.selectbox(tr("MD libre enregistrée"), list(range(len(results))), key=key + "_result",
                            format_func=result_labels.__getitem__)
    result = results[selection]
    if result.get("resumable"):
        st.button(tr("Reprendre la simulation"), key=key + "_resume_" + result["run_id"], disabled=busy,
                  on_click=_resume, args=(saved, result["run_dir"]))
    _render_result(result, key + "_" + result["run_id"])
