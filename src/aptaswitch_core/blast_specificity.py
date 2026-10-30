"""Explicit remote BLASTn screening against independently selected references.

No request is made at import, normalization, fingerprinting or result parsing.
Only ``search_taxa`` and an explicitly invoked ``BlastClient.run`` use NCBI.
Protocol: https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
Rate limits: https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import threading
import time
from typing import Callable
import unicodedata
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from .localization import tr


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATABASES = {"genomic": "refseq_genomes", "transcript": "refseq_rna"}
MIN_QUERY_LENGTH = 7
MAX_HITS = 500
MAX_TAXA = 20
TOOL_NAME = "Apta2Switch-Studio"
MAX_RESPONSE_BYTES = 32 * 1024 * 1024
SEARCH_PARAMETERS_VERSION = 1


class BlastRemoteError(RuntimeError):
    def __init__(self, message: str, *, raw_report: str = ""):
        super().__init__(message)
        self.raw_report = raw_report


class BlastCancelled(RuntimeError):
    """Raised internally when cancellation is observed before a network call."""


@dataclass(frozen=True)
class Taxon:
    taxid: int
    name: str

    def normalized(self) -> "Taxon":
        if type(self.taxid) is not int or not 1 <= self.taxid <= 2**31 - 1:
            raise ValueError(tr("L’identifiant taxonomique NCBI doit être un entier positif valide."))
        if not isinstance(self.name, str) or not self.name.strip() or len(self.name) > 300:
            raise ValueError(tr("Sélectionnez un organisme identifié dans NCBI Taxonomy."))
        return Taxon(self.taxid, " ".join(self.name.split()))

    def to_dict(self) -> dict:
        return asdict(self.normalized())


@dataclass(frozen=True)
class BlastScreenRequest:
    sequence: str
    taxa: tuple[Taxon, ...]
    scopes: tuple[str, ...] = ("genomic", "transcript")
    max_hits: int = 100

    def normalized(self) -> "BlastScreenRequest":
        if not isinstance(self.sequence, str):
            raise ValueError(tr("La séquence du trigger doit contenir des bases A, C, G, T ou U."))
        sequence = "".join(self.sequence.split()).upper().replace("U", "T")
        if len(sequence) < MIN_QUERY_LENGTH or set(sequence) - set("ACGT"):
            raise ValueError(tr("BLASTn nécessite au moins sept bases non ambiguës pour le trigger."))
        if not isinstance(self.taxa, (tuple, list)) or not 1 <= len(self.taxa) <= MAX_TAXA:
            raise ValueError(tr("Sélectionnez entre un et vingt organismes pour BLASTn."))
        taxa = []
        for taxon in self.taxa:
            if not isinstance(taxon, Taxon):
                raise ValueError(tr("Sélectionnez un organisme identifié dans NCBI Taxonomy."))
            taxon = taxon.normalized()
            if taxon.taxid not in {entry.taxid for entry in taxa}:
                taxa.append(taxon)
        if (not isinstance(self.scopes, (tuple, list)) or not self.scopes
                or any(not isinstance(scope, str) or scope not in DATABASES for scope in self.scopes)):
            raise ValueError(tr("Sélectionnez les génomes de référence, les transcrits ARN, ou les deux."))
        if type(self.max_hits) is not int or not 1 <= self.max_hits <= MAX_HITS:
            raise ValueError(tr("Le nombre maximal de séquences BLAST doit être compris entre 1 et 500."))
        return BlastScreenRequest(sequence, tuple(taxa), tuple(dict.fromkeys(self.scopes)), self.max_hits)

    def to_dict(self) -> dict:
        request = self.normalized()
        return {"sequence": request.sequence, "taxa": [taxon.to_dict() for taxon in request.taxa],
                "scopes": list(request.scopes), "max_hits": request.max_hits}

    def fingerprint(self) -> str:
        request = self.normalized()
        data = {"sequence": request.sequence, "taxids": sorted(taxon.taxid for taxon in request.taxa),
                "scopes": sorted(request.scopes), "max_hits": request.max_hits,
                "parameters_version": SEARCH_PARAMETERS_VERSION}
        return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def build_blast_parameters(request: BlastScreenRequest, taxon: Taxon, scope: str) -> dict[str, str]:
    """Return the exact short-query settings and taxon restriction to submit."""
    request, taxon = request.normalized(), taxon.normalized()
    if scope not in request.scopes or taxon.taxid not in {item.taxid for item in request.taxa}:
        raise ValueError(tr("La recherche BLAST ne correspond pas à la sélection autorisée."))
    restriction = f"txid{taxon.taxid}[Organism:exp]"
    if scope == "genomic":
        restriction += " AND biomol_genomic[PROP]"
    return {
        "CMD": "Put", "PROGRAM": "blastn", "DATABASE": DATABASES[scope],
        "QUERY": f">trigger\n{request.sequence}\n", "ENTREZ_QUERY": restriction,
        "WORD_SIZE": "7", "EXPECT": "1000", "FILTER": "F",
        "NUCL_REWARD": "1", "NUCL_PENALTY": "-3", "GAPCOSTS": "5 2",
        "SHORT_QUERY_ADJUST": "false", "HITLIST_SIZE": str(request.max_hits),
        "FORMAT_TYPE": "JSON2_S", "tool": TOOL_NAME,
    }


def _read_response(opener: Callable, request: Request, timeout: float) -> str:
    try:
        response = opener(request, timeout=timeout)
        try:
            data = response.read(MAX_RESPONSE_BYTES + 1)
        finally:
            response.close()
    except HTTPError as exc:
        raw = exc.read(MAX_RESPONSE_BYTES).decode("utf-8", errors="replace")
        raise BlastRemoteError(tr("NCBI a renvoyé une erreur HTTP {code}.", code=exc.code), raw_report=raw) from exc
    except (URLError, OSError, TimeoutError) as exc:
        # URLError text can include URLs. Our BLAST URL never contains a query
        # sequence: all submissions and retrievals use POST form data.
        raise BlastRemoteError(tr("La connexion à NCBI a échoué : {error}", error=exc)) from exc
    if len(data) > MAX_RESPONSE_BYTES:
        raise BlastRemoteError(tr("La réponse NCBI dépasse la taille maximale autorisée."))
    return data.decode("utf-8", errors="replace")


def _request(url: str, parameters: dict, *, post: bool = False) -> Request:
    encoded = urlencode(parameters).encode("ascii")
    headers = {"User-Agent": TOOL_NAME + "/BLASTn", "Accept": "application/json, application/xml, text/plain"}
    if post:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
        return Request(url, data=encoded, headers=headers, method="POST")
    return Request(url + "?" + encoded.decode("ascii"), headers=headers)


def search_taxa(query: str, *, opener: Callable | None = None, timeout: float = 30,
                max_results: int = 20, sleep: Callable[[float], None] = time.sleep) -> list[Taxon]:
    """Return actual NCBI Taxonomy choices; never choose an organism implicitly."""
    if not isinstance(query, str) or not query.strip() or len(query) > 200:
        raise ValueError(tr("Indiquez un nom d’organisme ou un identifiant NCBI Taxonomy."))
    if type(max_results) is not int or not 1 <= max_results <= 50:
        raise ValueError(tr("Le nombre de résultats taxonomiques doit être compris entre 1 et 50."))
    term = " ".join(query.split())
    alias = "".join(char for char in unicodedata.normalize("NFKD", term.casefold())
                    if not unicodedata.combining(char))
    term = {"ble": "Triticum aestivum", "bles": "Triticum aestivum", "ble tendre": "Triticum aestivum",
            "wheat": "Triticum aestivum", "orge": "Hordeum vulgare", "barley": "Hordeum vulgare"}.get(alias, term)
    opener = opener or urlopen
    parameters = {"db": "taxonomy", "term": term if term.isdecimal() else f"{term}[All Names]",
                  "retmode": "json", "retmax": max_results, "tool": TOOL_NAME}
    raw = _read_response(opener, _request(EUTILS_URL + "esearch.fcgi", parameters), timeout)
    try:
        data = json.loads(raw)
        if data.get("error") or "esearchresult" not in data:
            raise ValueError("missing search result")
        search = data["esearchresult"]
        if search.get("ERROR") or search.get("error"):
            raise ValueError("taxonomy error")
        ids = search["idlist"]
        if not isinstance(ids, list) or any(not str(item).isdecimal() for item in ids):
            raise ValueError("invalid taxonomy identifiers")
    except (ValueError, TypeError, KeyError, AttributeError) as exc:
        raise BlastRemoteError(tr("NCBI Taxonomy a renvoyé une réponse non valide."), raw_report=raw) from exc
    if not ids:
        return []
    sleep(0.35)  # E-utilities permits three requests per second without an API key.
    raw = _read_response(opener, _request(EUTILS_URL + "efetch.fcgi", {
        "db": "taxonomy", "id": ",".join(map(str, ids)), "retmode": "xml", "tool": TOOL_NAME,
    }), timeout)
    try:
        root = ET.fromstring(raw)
        if root.tag != "TaxaSet" or root.findall(".//ERROR"):
            raise ValueError("missing TaxaSet")
        result = [Taxon(int(node.findtext("TaxId", "")), node.findtext("ScientificName", "")).normalized()
                  for node in root.findall("./Taxon")]
        if not result:
            raise ValueError("missing taxon records")
        return result
    except (ET.ParseError, ValueError, TypeError) as exc:
        raise BlastRemoteError(tr("NCBI Taxonomy a renvoyé une réponse non valide."), raw_report=raw) from exc


def _number(value, *, integer: bool = False):
    if isinstance(value, bool):
        raise ValueError("invalid numeric BLAST field")
    result = int(value) if integer else float(value)
    if not math.isfinite(result) or (integer and float(value) != result):
        raise ValueError("invalid numeric BLAST field")
    return result


def _strand(value, start: int, end: int) -> str:
    if str(value).casefold() in {"plus", "+", "1"}:
        return "+"
    if str(value).casefold() in {"minus", "-", "-1"}:
        return "-"
    return "+" if end >= start else "-"


def _longest_match(query: str, subject: str) -> int:
    longest = current = 0
    for left, right in zip(query.upper(), subject.upper()):
        current = current + 1 if left == right and left in "ACGTU" else 0
        longest = max(longest, current)
    return longest


def _canonical_filter(value: str) -> str:
    return re.sub(r"\s+", "", value.casefold()).replace("[organism:exp]", "[orgn]")


def parse_blast_json(raw: str, *, query_length: int, expected_database: str = "",
                     expected_filter: str = "", requested_taxid: int | None = None) -> dict:
    """Read JSON2/JSON2_S without treating errors or missing records as no hits."""
    try:
        data = json.loads(raw)
        records = data.get("BlastOutput2", data) if isinstance(data, dict) else data
        if isinstance(records, dict):
            records = [records]
        if not isinstance(records, list) or len(records) != 1:
            raise ValueError("expected exactly one query report")
        record = records[0]
        if record.get("error"):
            raise ValueError(str(record["error"]))
        report = record["report"]
        if report.get("error"):
            raise ValueError(str(report["error"]))
        if report.get("program") != "blastn":
            raise ValueError("missing or incompatible BLAST program")
        search = report["results"]["search"]
        if search.get("error"):
            raise ValueError(str(search["error"]))
        reported_length = _number(search["query_len"], integer=True)
        if reported_length != query_length:
            raise ValueError("query length differs from submitted trigger")
        hits = search["hits"]
        if not isinstance(hits, list):
            raise ValueError("missing hit list")
        if re.search(r"\b(error|failed|failure|invalid)\b", str(search.get("message", "")), re.IGNORECASE):
            raise ValueError(str(search["message"]))
        database = str(report.get("search_target", {}).get("db", ""))
        if expected_database and database != expected_database:
            raise ValueError("reported database differs from submitted database")
        params = report.get("params", {})
        echoed_filter = params.get("entrez_query", params.get("entrez-query"))
        filter_verified = bool(expected_filter and echoed_filter is not None
                               and _canonical_filter(str(echoed_filter)) == _canonical_filter(expected_filter))
        if expected_filter and echoed_filter is not None and not filter_verified:
            raise ValueError("reported taxonomic restriction differs from submitted restriction")
        rows = []
        for hit_index, hit in enumerate(hits, 1):
            descriptions = hit.get("description", [])
            if not isinstance(descriptions, list) or not descriptions or any(not isinstance(item, dict) for item in descriptions):
                raise ValueError("missing subject descriptions")
            primary = descriptions[0]
            taxid = primary.get("taxid")
            taxid = None if taxid is None else _number(taxid, integer=True)
            if taxid is not None and taxid <= 0:
                raise ValueError("invalid subject taxonomy identifier")
            subject_length = _number(hit["len"], integer=True)
            hsps = hit["hsps"]
            if not isinstance(hsps, list) or not hsps or subject_length <= 0:
                raise ValueError("missing alignments")
            for hsp_index, hsp in enumerate(hsps, 1):
                qseq, hseq = str(hsp["qseq"]), str(hsp["hseq"])
                length = _number(hsp["align_len"], integer=True)
                aligned_identities = sum(left.upper() == right.upper() and left != "-"
                                         for left, right in zip(qseq, hseq))
                aligned_gaps = sum(left == "-" or right == "-" for left, right in zip(qseq, hseq))
                identities = _number(hsp.get("identity", aligned_identities), integer=True)
                gaps = _number(hsp.get("gaps", aligned_gaps), integer=True)
                qfrom, qto, hfrom, hto = [_number(hsp[field], integer=True)
                                         for field in ("query_from", "query_to", "hit_from", "hit_to")]
                if (length <= 0 or len(qseq) != length or len(hseq) != length
                        or not 0 <= identities <= length or not 0 <= gaps <= length - identities
                        or identities != aligned_identities or gaps != aligned_gaps
                        or any(left == right == "-" for left, right in zip(qseq, hseq))
                        or not 1 <= min(qfrom, qto) <= max(qfrom, qto) <= query_length
                        or not 1 <= min(hfrom, hto) <= max(hfrom, hto) <= subject_length
                        or len(qseq.replace("-", "")) != abs(qto - qfrom) + 1
                        or len(hseq.replace("-", "")) != abs(hto - hfrom) + 1):
                    raise ValueError("inconsistent BLAST alignment coordinates")
                evalue = _number(hsp["evalue"])
                bit_score = _number(hsp["bit_score"])
                if evalue < 0:
                    raise ValueError("negative e-value")
                row = {
                    "accession": str(primary.get("accession", "")),
                    "accessions": list(dict.fromkeys(str(item["accession"]) for item in descriptions if item.get("accession"))),
                    "subject_id": str(primary.get("id", "")), "title": str(primary.get("title", "")),
                    "taxid": taxid, "scientific_name": str(primary.get("sciname", "")),
                    "descriptions": descriptions, "subject_length": subject_length,
                    "hit_number": hit_index, "hsp_number": hsp_index,
                    "query_from": qfrom, "query_to": qto, "hit_from": hfrom, "hit_to": hto,
                    "query_strand": _strand(hsp.get("query_strand", hsp.get("query_frame")), qfrom, qto),
                    "hit_strand": _strand(hsp.get("hit_strand", hsp.get("hit_frame")), hfrom, hto),
                    "alignment_length": length, "identities": identities,
                    "identity_pct": 100.0 * identities / length,
                    "query_coverage_pct": 100.0 * (abs(qto - qfrom) + 1) / query_length,
                    "mismatches": length - identities - gaps, "gaps": gaps,
                    "evalue": evalue, "bit_score": bit_score, "score": _number(hsp["score"]),
                    "longest_contiguous_match": _longest_match(qseq, hseq),
                    "full_length_exact": (min(qfrom, qto) == 1 and max(qfrom, qto) == query_length
                                          and length == query_length and identities == query_length and gaps == 0),
                    "qseq": qseq, "hseq": hseq, "midline": str(hsp.get("midline", "")),
                    "taxon_match": "exact" if requested_taxid is not None and taxid == requested_taxid else "descendant_or_unverified",
                }
                row["match_orientation"] = ("same_strand" if row["query_strand"] == row["hit_strand"]
                                            else "reverse_complement")
                rows.append(row)
        return {
            "hits": rows, "subject_hit_count": len(hits), "query_length": reported_length,
            "reported_database": database, "reported_filter": echoed_filter,
            "database_verified": bool(expected_database and database == expected_database),
            "filter_verified": filter_verified,
            "scope_warning": "" if filter_verified else tr("NCBI n’a pas confirmé le filtre taxonomique dans le rapport ; la portée des résultats reste à vérifier."),
            "database_statistics": search.get("stat", {}), "reported_parameters": params,
            "blast_version": report.get("version", ""), "reported_message": search.get("message", ""),
        }
    except (ValueError, KeyError, TypeError, AttributeError, OverflowError) as exc:
        raise BlastRemoteError(tr("Le rapport BLAST reçu est invalide ou ne correspond pas à la recherche : {error}", error=exc), raw_report=raw) from exc


class _RateGate:
    def __init__(self, clock: Callable = time.monotonic, sleep: Callable = time.sleep):
        self.clock, self.sleep = clock, sleep
        self.last_contact: float | None = None
        self.lock = threading.Lock()

    def wait(self, cancel_check: Callable[[], bool], *, earliest: float = 0) -> None:
        while not self.lock.acquire(timeout=0.1):
            if cancel_check():
                raise BlastCancelled()
        try:
            target = max(earliest, self.last_contact + 10 if self.last_contact is not None else 0)
            while self.clock() < target:
                if cancel_check():
                    raise BlastCancelled()
                self.sleep(min(0.25, target - self.clock()))
            if cancel_check():
                raise BlastCancelled()
            self.last_contact = self.clock()
        finally:
            self.lock.release()


_BLAST_GATE = _RateGate()


class BlastClient:
    """One explicit execution, with cancellable waits and no automatic resubmit."""

    def __init__(self, *, opener: Callable | None = None, clock: Callable = time.monotonic,
                 sleep: Callable = time.sleep, timeout: float = 45, max_wait_seconds: float = 1800,
                 email: str = "", rate_gate: _RateGate | None = None):
        if any(not math.isfinite(value) or value <= 0 for value in (timeout, max_wait_seconds)):
            raise ValueError(tr("Les délais de connexion et d’attente BLAST doivent être strictement positifs."))
        self.opener, self.clock, self.sleep = opener or urlopen, clock, sleep
        self.timeout, self.max_wait_seconds, self.email = timeout, max_wait_seconds, email
        self.gate = rate_gate or (_BLAST_GATE if clock is time.monotonic and sleep is time.sleep else _RateGate(clock, sleep))

    def _contact(self, params: dict, cancel_check: Callable[[], bool], *, earliest: float = 0) -> str:
        self.gate.wait(cancel_check, earliest=earliest)
        params = dict(params)
        if self.email:
            params["email"] = self.email
        return _read_response(self.opener, _request(BLAST_URL, params, post=True), self.timeout)

    def run(self, request: BlastScreenRequest, *, progress_callback: Callable[[int, int, str], None] | None = None,
            cancel_check: Callable[[], bool] | None = None) -> dict:
        request = request.normalized()
        cancelled = cancel_check or (lambda: False)
        combinations = [(taxon, scope) for taxon in request.taxa for scope in request.scopes]
        result = {"status": "completed", "request": request.to_dict(),
                  "fingerprint": request.fingerprint(), "runs": []}
        for index, (taxon, scope) in enumerate(combinations):
            parameters = build_blast_parameters(request, taxon, scope)
            if self.email:
                parameters["email"] = self.email
            run = {"taxon": taxon.to_dict(), "scope": scope, "evidence_type": scope,
                   "database": DATABASES[scope], "params": parameters, "status": "failed",
                   "rid": "", "hits": [], "raw_report": "", "raw_format": "json",
                   "queried_at": datetime.now(timezone.utc).isoformat()}
            result["runs"].append(run)
            try:
                if cancelled():
                    raise BlastCancelled()
                if progress_callback:
                    progress_callback(index, len(combinations), tr("BLASTn · envoi pour {organism} · {scope}", organism=taxon.name, scope=scope))
                submitted = self._contact(parameters, cancelled)
                rid = re.search(r"\bRID\s*=\s*([A-Za-z0-9-]+)", submitted)
                rtoe = re.search(r"\bRTOE\s*=\s*(\d+)", submitted)
                if not rid:
                    raise BlastRemoteError(tr("NCBI n’a pas fourni d’identifiant pour la recherche BLAST."), raw_report=submitted)
                run["rid"] = rid.group(1)
                run["estimated_seconds"] = int(rtoe.group(1)) if rtoe else 60
                last_poll = self.clock()
                deadline = last_poll + self.max_wait_seconds
                next_poll = last_poll + max(60, run["estimated_seconds"])
                while True:
                    if next_poll > deadline:
                        raise BlastRemoteError(tr("Le délai d’attente BLAST est dépassé. La recherche distante peut encore être en cours."))
                    if progress_callback:
                        progress_callback(index, len(combinations), tr("BLASTn · attente NCBI · {organism} · RID {rid}", organism=taxon.name, rid=run["rid"]))
                    raw = self._contact({"CMD": "Get", "RID": run["rid"], "FORMAT_TYPE": "JSON2_S", "tool": TOOL_NAME},
                                        cancelled, earliest=next_poll)
                    last_poll = self.clock()
                    if cancelled():
                        raise BlastCancelled()
                    status = (None if raw.lstrip().startswith(("{", "[")) else
                              re.search(r"\bStatus\s*=\s*(WAITING|READY|FAILED|UNKNOWN)\b", raw, re.IGNORECASE))
                    if status and status.group(1).upper() in {"WAITING", "READY"}:
                        next_poll = last_poll + 60
                        continue
                    run["raw_report"] = raw
                    if status:
                        raise BlastRemoteError(tr("NCBI indique que la recherche BLAST a échoué ou a expiré : {status}.", status=status.group(1).upper()), raw_report=raw)
                    parsed = parse_blast_json(raw, query_length=len(request.sequence), expected_database=run["database"],
                                              expected_filter=parameters["ENTREZ_QUERY"], requested_taxid=taxon.taxid)
                    run.update(parsed, status="completed")
                    run["hit_limit_reached"] = parsed["subject_hit_count"] >= request.max_hits
                    if run["hit_limit_reached"]:
                        run["limit_warning"] = tr("Le nombre maximal de séquences retournées est atteint ; d’autres correspondances peuvent ne pas être affichées.")
                    break
            except BlastCancelled:
                run["status"] = "cancelled"
                result["status"] = "cancelled"
                break
            except (BlastRemoteError, ValueError) as exc:
                run["error"] = str(exc)
                if getattr(exc, "raw_report", ""):
                    run["raw_report"] = exc.raw_report
                    if not exc.raw_report.lstrip().startswith(("{", "[")):
                        run["raw_format"] = "text"
            if progress_callback:
                progress_callback(index + 1, len(combinations), tr("BLASTn · recherche terminée pour {organism} · {scope}", organism=taxon.name, scope=scope))
        if result["status"] != "cancelled":
            completed = sum(run["status"] == "completed" for run in result["runs"])
            result["status"] = "completed" if completed == len(combinations) else "partial" if completed else "failed"
        return result
