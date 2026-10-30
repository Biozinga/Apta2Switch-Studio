"""Explicit, consented NCBI screening of a design's current trigger."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

import streamlit as st
import pandas as pd

from aptaswitch_core.blast_specificity import BlastScreenRequest, Taxon, search_taxa
from aptaswitch_core.localization import set_language, tr
from .blast_runtime import BlastScreenJob
from .design_session import active_job_running, save_active_design_safely
from .table_ui import render_readonly_table


SCOPE_LABELS = {"genomic": "Génomes de référence", "transcript": "Transcrits ARN de référence"}


def _request() -> BlastScreenRequest:
    return BlastScreenRequest(
        sequence=st.session_state.get("trigger_input", ""),
        taxa=tuple(Taxon(int(taxon["taxid"]), taxon["name"])
                   for taxon in st.session_state.get("blast_taxa", [])),
        scopes=tuple(st.session_state.get("blast_scopes", ["genomic", "transcript"])),
        max_hits=int(st.session_state.get("blast_max_hits", 100)),
    ).normalized()


def _find_taxa() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.blast_lookup_results = []
    st.session_state.pop("blast_lookup_selected", None)
    st.session_state.pop("blast_lookup_error", None)
    try:
        query = st.session_state.get("blast_lookup_query", "").strip()
        if not query:
            raise ValueError(tr("Saisissez un organisme à rechercher."))
        st.session_state.blast_lookup_results = [asdict(taxon) for taxon in search_taxa(query)]
        if not st.session_state.blast_lookup_results:
            st.session_state.blast_lookup_error = tr("Aucun organisme trouvé. Essayez son nom scientifique ou son TaxID NCBI.")
    except (OSError, ValueError, RuntimeError) as exc:
        st.session_state.blast_lookup_error = str(exc)


def _add_taxon() -> None:
    selected = st.session_state.get("blast_lookup_selected")
    candidates = st.session_state.get("blast_lookup_results", [])
    chosen = next((taxon for taxon in candidates if taxon["taxid"] == selected), None)
    if chosen and all(taxon["taxid"] != selected for taxon in st.session_state.get("blast_taxa", [])):
        if len(st.session_state.get("blast_taxa", [])) >= 20:
            st.session_state.blast_lookup_error = tr("Sélectionnez au maximum 20 organismes par contrôle.")
            return
        st.session_state.blast_taxa = [*st.session_state.get("blast_taxa", []), chosen]
        st.session_state.blast_consent = False
        save_active_design_safely(st.session_state)


def _remove_taxon(taxid: int) -> None:
    st.session_state.blast_taxa = [taxon for taxon in st.session_state.get("blast_taxa", [])
                                 if taxon["taxid"] != taxid]
    st.session_state.blast_consent = False
    save_active_design_safely(st.session_state)


def _start_blast() -> None:
    set_language(st.session_state.get("ui_language", "en"))
    st.session_state.pop("blast_start_error", None)
    try:
        request = _request()
        if (not st.session_state.get("blast_consent")
                or st.session_state.get("blast_consent_signature") != request.fingerprint()):
            raise ValueError(tr("Confirmez l’envoi de ce trigger à NCBI avant de lancer le contrôle."))
        if active_job_running(st.session_state):
            raise ValueError(tr("Attendez la fin du calcul en cours avant de lancer BLASTn."))
        active = st.session_state.get("active_design_path")
        if not active:
            raise ValueError(tr("Ouvrez un design avant de lancer BLASTn."))
        # Consume permission for this specific submission. Changing a sequence,
        # organisms or databases always needs a new explicit choice.
        st.session_state.blast_consent = False
        save_active_design_safely(st.session_state)
        if st.session_state.get("design_save_error"):
            raise ValueError(st.session_state.design_save_error)
        job = BlastScreenJob(request, Path(active) / "calculations")
        st.session_state.blast_job = job
        job.start()
    except (OSError, ValueError, RuntimeError) as exc:
        st.session_state.blast_start_error = str(exc)


def _cancel_blast() -> None:
    job = st.session_state.get("blast_job")
    if job is not None:
        job.cancel()


def _same_query(result: dict, fingerprint: str | None) -> bool:
    return bool(fingerprint and result.get("fingerprint") == fingerprint)


def _show_result(result: dict, fingerprint: str | None) -> None:
    if not _same_query(result, fingerprint):
        st.warning(tr("Ce contrôle concerne un ancien trigger ou une autre sélection d’organismes. Relancez-le pour les paramètres actuels."))
    status = result.get("status")
    if status in {"failed", "partial", "cancelled"}:
        st.warning(tr("Le contrôle est incomplet. Les recherches en erreur ou annulées ne permettent pas de conclure à une absence de correspondance."))
    query = result.get("request", {}).get("sequence", "")
    st.caption(tr("Trigger contrôlé : {sequence}", sequence=query))
    if result.get("completed_at"):
        try:
            date = datetime.fromisoformat(result["completed_at"]).astimezone().strftime("%Y-%m-%d %H:%M")
        except (TypeError, ValueError):
            date = str(result["completed_at"])
        st.caption(tr("Contrôle enregistré : {date}", date=date))
    tabs = st.tabs([tr(SCOPE_LABELS[scope]) for scope in ("genomic", "transcript")])
    for scope, tab in zip(("genomic", "transcript"), tabs):
        with tab:
            runs = [run for run in result.get("runs", []) if run.get("scope") == scope]
            if not runs:
                st.caption(tr("Cette catégorie n’a pas été recherchée."))
                continue
            for index, run in enumerate(runs):
                taxon = run.get("taxon", {})
                st.markdown(f"**{taxon.get('name', '')}** · TaxID {taxon.get('taxid', '')}")
                st.caption(tr("Base : {database} · Requête NCBI : {rid}", database=run.get("database", ""), rid=run.get("rid") or "—"))
                if run.get("scope_warning"):
                    st.warning(tr(run["scope_warning"]))
                if run.get("limit_warning"):
                    st.warning(tr(run["limit_warning"]))
                if run.get("status") != "completed":
                    st.warning(tr("Recherche non terminée : {error}", error=run.get("error") or run.get("status", "")))
                    continue
                hits = run.get("hits", [])
                if not hits:
                    st.info(tr("Aucune correspondance rapportée par cette recherche. Cela ne valide pas l’absence de faux positifs."))
                    continue
                rows = []
                for hit in hits:
                    rows.append({
                        tr("Accession"): hit.get("accession"),
                        tr("Organisme"): hit.get("scientific_name") or tr("Non renseigné par NCBI"),
                        "TaxID": hit.get("taxid"),
                        tr("Identité (%)"): hit.get("identity_pct"),
                        tr("Trigger couvert (%)"): hit.get("query_coverage_pct"),
                        tr("Longueur alignée"): hit.get("alignment_length"),
                        tr("Mésappariements"): hit.get("mismatches"),
                        tr("Gaps"): hit.get("gaps"),
                        tr("Brin aligné"): hit.get("hit_strand"),
                        tr("Relation au trigger"): tr({"same_strand": "Même sens", "reverse_complement": "Complément inversé"}.get(hit.get("match_orientation"), "Non précisée")),
                        "E-value": hit.get("evalue"),
                        tr("Description"): hit.get("title"),
                    })
                render_readonly_table(pd.DataFrame(rows).style.format({
                    tr("Identité (%)"): "{:.1f}",
                    tr("Trigger couvert (%)"): "{:.1f}",
                    "E-value": "{:.2e}",
                }, na_rep="—"), height=420 if len(rows) > 10 else None)
                st.caption(tr("Identité et couverture décrivent l’alignement, pas une probabilité d’activation du switch. La liste retournée peut être limitée."))
                with st.expander(tr("Voir les alignements et les séquences candidates")):
                    key = f"blast_alignment_{scope}_{index}"
                    if st.session_state.get(key, 0) not in range(len(hits)):
                        st.session_state[key] = 0
                    selected = st.selectbox(tr("Alignement"), range(len(hits)), key=key,
                                            format_func=lambda i: f"{hits[i].get('accession', '')} · {i + 1}")
                    hit = hits[selected]
                    accession = str(hit.get("accession", ""))
                    if accession:
                        st.link_button(tr("Ouvrir la séquence dans NCBI"), "https://www.ncbi.nlm.nih.gov/nuccore/" + quote(accession, safe=""))
                    st.code(f"Trigger   {hit.get('qseq', '')}\n          {hit.get('midline', '')}\n{tr('Séquence'):<9} {hit.get('hseq', '')}", language=None)
                    st.caption(tr("Le brin complémentaire peut présenter un risque d’hybridation ; il n’est pas automatiquement équivalent à un trigger activateur dans le sens direct."))
                    if scope == "transcript":
                        predicted = accession.startswith(("XM_", "XR_"))
                        st.caption(tr("Transcrit prédit par annotation : sa présence dans votre échantillon n’est pas démontrée.") if predicted
                                   else tr("Référence de transcrit : l’expression dans votre tissu et sa libération en solution restent à vérifier."))
    exports = st.session_state.get("blast_exports", {})
    controls = st.columns(2)
    for column, kind, label, mime in ((controls[0], "json", "Télécharger le contrôle (JSON)", "application/json"),
                                     (controls[1], "csv", "Télécharger les correspondances (CSV)", "text/csv")):
        path = exports.get(kind)
        if path and Path(path).is_file():
            column.download_button(tr(label), Path(path).read_bytes(), file_name=Path(path).name,
                                   mime=mime, key=f"blast_download_{kind}")


@st.fragment(run_every=2.0)
def _blast_monitor(fingerprint: str | None) -> None:
    set_language(st.session_state.get("ui_language", "en"))
    job = st.session_state.get("blast_job")
    if job is not None:
        snapshot = job.snapshot()
        if snapshot.state == "running":
            st.progress(min(1.0, snapshot.current / max(1, snapshot.total)), text=tr(snapshot.message))
            st.caption(tr("NCBI traite les recherches à distance. Cela peut prendre plusieurs minutes ; l’application respecte les délais de consultation du service."))
            st.button(tr("Annuler le contrôle BLASTn"), key="blast_cancel", on_click=_cancel_blast)
        else:
            if job.consume_once():
                if snapshot.result is not None:
                    st.session_state.blast_result = snapshot.result
                    st.session_state.blast_exports = snapshot.exports or {}
                if snapshot.error:
                    st.session_state.blast_start_error = snapshot.error
                save_active_design_safely(st.session_state)
                st.rerun(scope="app")
    result = st.session_state.get("blast_result")
    if result:
        if job is not None and job.snapshot().state == "running":
            st.caption(tr("Dernier contrôle enregistré — une nouvelle recherche est en cours."))
        _show_result(result, fingerprint)


def render_blast_screen(trigger: str) -> None:
    defaults = {"blast_taxa": [], "blast_scopes": ["genomic", "transcript"],
                "blast_max_hits": 100, "blast_consent": False}
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value
    with st.container(border=True):
        st.subheader(tr("Contrôle du trigger dans l’échantillon · BLASTn"),
                     help=tr("Avant l’extension ou la conception du toehold switch, recherchez les séquences ressemblantes dans les organismes présents dans l’échantillon."))
        with st.expander(tr("Ce que ce contrôle peut établir")):
            st.markdown(tr("**Génomes** : recherche dans les génomes de référence, donc l’ADN génomique pour les plantes et champignons. Certains virus possèdent un génome ARN.\n\n**Transcrits ARN** : recherche dans les transcrits de référence, parfois prédits.\n\n**En solution** : BLAST ne mesure ni l’expression, ni la quantité, ni la libération de ce brin. Cela dépend notamment du tissu, de la lyse, de l’extraction et de la dégradation. Un alignement seul ne démontre pas un faux positif."))
            st.markdown(tr("Documentation : [BLAST NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) · [Références génomiques et transcrits RefSeq](https://www.ncbi.nlm.nih.gov/refseq/about/)"))
        running = active_job_running(st.session_state)
        st.text_input(tr("Rechercher un organisme (nom ou TaxID)"), key="blast_lookup_query",
                      placeholder=tr("Triticum aestivum (blé), Hordeum vulgare (orge), Fusarium…"))
        st.button(tr("Rechercher dans NCBI Taxonomy"), key="blast_find_taxa", on_click=_find_taxa, disabled=running, width="stretch")
        choices = {item["taxid"]: item for item in st.session_state.get("blast_lookup_results", [])}
        if choices:
            if st.session_state.get("blast_lookup_selected") not in choices:
                st.session_state.blast_lookup_selected = next(iter(choices))
            st.selectbox(tr("Organisme à ajouter"), list(choices), key="blast_lookup_selected",
                         format_func=lambda taxid: f"{choices[taxid]['name']} · TaxID {taxid}", disabled=running)
            st.button(tr("Ajouter à l’échantillon"), key="blast_add_taxon", on_click=_add_taxon, disabled=running, width="stretch")
        if st.session_state.get("blast_lookup_error"):
            st.warning(st.session_state.blast_lookup_error)
        for taxon in st.session_state.blast_taxa:
            name, remove = st.columns([5, 1])
            name.write(f"{taxon['name']} · TaxID {taxon['taxid']}")
            remove.button(tr("Retirer"), key=f"blast_remove_{taxon['taxid']}",
                          on_click=_remove_taxon, args=(taxon["taxid"],), disabled=running, width="stretch")
        if not st.session_state.blast_taxa:
            st.caption(tr("Ajoutez au moins un organisme. Choisissez une espèce précise lorsque vous la connaissez."))
        scope_labels = {scope: tr(label) for scope, label in SCOPE_LABELS.items()}
        st.multiselect(tr("Séquences à rechercher"), ["genomic", "transcript"], key="blast_scopes",
                       format_func=scope_labels.__getitem__, disabled=running)
        st.number_input(tr("Maximum de séquences retournées par organisme et catégorie"),
                        min_value=1, max_value=500, step=1, key="blast_max_hits", disabled=running)
        fingerprint = None
        request_error = ""
        try:
            request = _request()
            fingerprint = request.fingerprint()
        except ValueError as exc:
            request = None
            request_error = str(exc)
        if st.session_state.get("blast_consent_signature") != fingerprint:
            st.session_state.blast_consent = False
            st.session_state.blast_consent_signature = fingerprint
        if trigger:
            st.caption(tr("Trigger à contrôler : {length} nt", length=len(trigger)))
            st.code(trigger, language=None)
        if trigger and len(trigger) < 18:
            st.caption(tr("Un trigger très court peut correspondre à de nombreuses séquences. L’absence de résultat BLAST ne suffit pas à établir sa spécificité."))
        if request_error and trigger and st.session_state.blast_taxa:
            st.caption(tr(request_error))
        st.checkbox(tr("J’autorise l’envoi de ce trigger et de cette sélection d’organismes à NCBI pour ce contrôle."),
                    key="blast_consent", disabled=running or request is None)
        st.button(tr("Contrôler le trigger avec BLASTn"), key="blast_start", on_click=_start_blast,
                  disabled=running or request is None or not st.session_state.blast_consent, type="primary")
        if st.session_state.get("blast_start_error"):
            st.error(st.session_state.blast_start_error)
        _blast_monitor(fingerprint)
