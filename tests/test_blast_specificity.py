from __future__ import annotations

import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlsplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.blast_specificity import (
    BLAST_URL, BlastClient, BlastRemoteError, BlastScreenRequest, Taxon, _RateGate,
    build_blast_parameters, parse_blast_json, search_taxa,
)


WHEAT = Taxon(4565, "Triticum aestivum")
BARLEY = Taxon(4513, "Hordeum vulgare")
SEQUENCE = "ACGTACGTACGT"  # Synthetic fixture only; no test performs network I/O.


def request(**changes):
    return BlastScreenRequest(**{"sequence": SEQUENCE, "taxa": (WHEAT,), "scopes": ("genomic",), **changes})


def hsp(**changes):
    return {
        "num": 1, "bit_score": 24.0, "score": 12, "evalue": 0.025,
        "query_from": 1, "query_to": 12, "hit_from": 101, "hit_to": 112,
        "query_strand": "Plus", "hit_strand": "Plus", "identity": 12,
        "align_len": 12, "gaps": 0, "qseq": SEQUENCE, "hseq": SEQUENCE,
        "midline": "|" * 12, **changes,
    }


def hit(*hsps, **changes):
    return {
        "num": 1, "description": [
            {"id": "ref|NC_000001.2|", "accession": "NC_000001.2", "title": "Synthetic reference chromosome",
             "taxid": 4565, "sciname": "Triticum aestivum"},
            {"id": "gb|AB000001.3|", "accession": "AB000001.3", "title": "Alternate accession",
             "taxid": 4565, "sciname": "Triticum aestivum"},
        ], "len": 1000, "hsps": list(hsps) or [hsp()], **changes,
    }


def report(*, hits=None, database="refseq_genomes", taxid=4565, genomic=True, echo=True):
    params = {"expect": 1000, "sc_match": 1, "sc_mismatch": -3, "gap_open": 5, "gap_extend": 2}
    if echo:
        params["entrez_query"] = f"txid{taxid}[Organism:exp]" + (" AND biomol_genomic[PROP]" if genomic else "")
    return {"BlastOutput2": [{"report": {
        "program": "blastn", "version": "BLASTN 2.16.0+", "search_target": {"db": database},
        "params": params, "results": {"search": {
            "query_id": "Query_1", "query_title": "trigger", "query_len": len(SEQUENCE),
            "hits": [hit()] if hits is None else hits,
            "stat": {"db_num": 10, "db_len": 2000},
        }},
    }}]}


class FakeClock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


class FakeOpener:
    def __init__(self, responses, clock=None):
        self.responses = list(responses)
        self.clock = clock or (lambda: 0)
        self.calls = []

    def __call__(self, request, *, timeout):
        self.calls.append((self.clock(), request, timeout))
        if not self.responses:
            raise AssertionError("Unexpected network request beyond the explicit fixture")
        result = self.responses.pop(0)
        if isinstance(result, BaseException):
            raise result
        return io.BytesIO(result.encode() if isinstance(result, str) else result)

    def params(self, index):
        req = self.calls[index][1]
        return {key: values[0] for key, values in parse_qs(
            req.data.decode() if req.data else urlsplit(req.full_url).query,
        ).items()}


class BlastRequestTests(unittest.TestCase):
    def test_normalization_and_fingerprint_are_offline_and_ignore_names_and_selection_order(self):
        with patch("aptaswitch_core.blast_specificity.urlopen", side_effect=AssertionError("network forbidden")):
            first = request(sequence=" acgu acgu\nacgu ", taxa=(WHEAT, BARLEY, WHEAT),
                            scopes=("transcript", "genomic", "genomic"))
            second = request(taxa=(Taxon(4513, "Barley"), Taxon(4565, "Wheat")),
                             scopes=("genomic", "transcript"))
            self.assertEqual(first.fingerprint(), second.fingerprint())
            self.assertEqual(first.normalized().sequence, SEQUENCE)
            self.assertEqual(first.to_dict()["taxa"], [WHEAT.to_dict(), BARLEY.to_dict()])
            self.assertEqual(first.to_dict()["scopes"], ["transcript", "genomic"])
            for changed in (request(sequence="ACGTACGTACGA"), request(max_hits=99),
                            request(taxa=(BARLEY,)), request(scopes=("transcript",))):
                self.assertNotEqual(request().fingerprint(), changed.fingerprint())

    def test_invalid_requests_do_not_contact_ncbi(self):
        invalid = [
            {"sequence": "ACGTAC"}, {"sequence": "ACGTNCGT"}, {"sequence": ">trigger\nACGTACGT"},
            {"sequence": None}, {"taxa": ()}, {"taxa": (WHEAT,) * 21}, {"taxa": ({"taxid": 4565},)},
            {"taxa": (Taxon(True, "Invalid"),)}, {"taxa": (Taxon(0, "Invalid"),)},
            {"taxa": (Taxon(2**31, "Invalid"),)}, {"taxa": (Taxon(1, ""),)},
            {"scopes": ()}, {"scopes": ("nt",)}, {"scopes": ({},)},
            {"max_hits": True}, {"max_hits": 0}, {"max_hits": 501}, {"max_hits": 1.2},
        ]
        opener = FakeOpener([])
        for changes in invalid:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                BlastClient(opener=opener).run(request(**changes))
        self.assertEqual(opener.calls, [])

    def test_short_query_parameters_and_independent_reference_scopes(self):
        query = request(scopes=("genomic", "transcript"), taxa=(WHEAT, BARLEY), max_hits=250)
        genomic = build_blast_parameters(query, WHEAT, "genomic")
        self.assertEqual(genomic["QUERY"], ">trigger\nACGTACGTACGT\n")
        self.assertEqual(genomic["DATABASE"], "refseq_genomes")
        self.assertEqual(genomic["ENTREZ_QUERY"], "txid4565[Organism:exp] AND biomol_genomic[PROP]")
        for key, value in {"PROGRAM": "blastn", "WORD_SIZE": "7", "EXPECT": "1000", "FILTER": "F",
                           "NUCL_REWARD": "1", "NUCL_PENALTY": "-3", "GAPCOSTS": "5 2",
                           "SHORT_QUERY_ADJUST": "false", "HITLIST_SIZE": "250"}.items():
            self.assertEqual(genomic[key], value)
        transcript = build_blast_parameters(query, BARLEY, "transcript")
        self.assertEqual(transcript["DATABASE"], "refseq_rna")
        self.assertEqual(transcript["ENTREZ_QUERY"], "txid4513[Organism:exp]")
        with self.assertRaises(ValueError):
            build_blast_parameters(request(), BARLEY, "genomic")
        with self.assertRaises(ValueError):
            build_blast_parameters(request(), WHEAT, "transcript")


class TaxonomyTests(unittest.TestCase):
    def test_search_returns_confirmed_taxa_and_ignores_nested_lineage_records(self):
        xml = """<TaxaSet><Taxon><TaxId>4565</TaxId><ScientificName>Triticum aestivum</ScientificName>
        <LineageEx><Taxon><TaxId>4564</TaxId><ScientificName>Triticum</ScientificName></Taxon></LineageEx>
        </Taxon><Taxon><TaxId>4513</TaxId><ScientificName>Hordeum vulgare</ScientificName></Taxon></TaxaSet>"""
        for name, term in ((" blés ", "Triticum aestivum"), ("ble tendre", "Triticum aestivum"),
                           ("barley", "Hordeum vulgare"), ("champignon", "champignon"), ("4565", "4565")):
            opener = FakeOpener([json.dumps({"esearchresult": {"idlist": ["4565", "4513"]}}), xml])
            pauses = []
            self.assertEqual(search_taxa(name, opener=opener, sleep=pauses.append), [WHEAT, BARLEY])
            self.assertEqual(opener.params(0)["term"], term if term.isdecimal() else term + "[All Names]")
            self.assertEqual(opener.params(1)["id"], "4565,4513")
            self.assertGreaterEqual(pauses[0], 1 / 3)

    def test_empty_taxonomy_result_is_distinct_from_invalid_or_remote_error(self):
        opener = FakeOpener(['{"esearchresult":{"idlist":[]}}'])
        self.assertEqual(search_taxa("no such organism", opener=opener), [])
        self.assertEqual(len(opener.calls), 1)
        for raw in ("not json", "[]", '{"error":"busy"}', '{"esearchresult":{"idlist":["not an id"]}}'):
            with self.subTest(raw=raw), self.assertRaises(BlastRemoteError):
                search_taxa("wheat", opener=FakeOpener([raw]))
        with self.assertRaises(BlastRemoteError):
            search_taxa("wheat", opener=FakeOpener(['{"esearchresult":{"idlist":["4565"]}}',
                                                  "<TaxaSet><ERROR>Busy</ERROR></TaxaSet>"]), sleep=lambda _: None)


class BlastParsingTests(unittest.TestCase):
    def parse(self, data, **options):
        return parse_blast_json(json.dumps(data), query_length=12, expected_database="refseq_genomes",
                                expected_filter="txid4565[Organism:exp] AND biomol_genomic[PROP]",
                                requested_taxid=4565, **options)

    def test_complete_hsp_preserves_accession_versions_alignment_and_provenance(self):
        parsed = self.parse(report())
        row = parsed["hits"][0]
        self.assertTrue(parsed["filter_verified"])
        self.assertTrue(parsed["database_verified"])
        self.assertEqual(parsed["blast_version"], "BLASTN 2.16.0+")
        self.assertEqual(row["accessions"], ["NC_000001.2", "AB000001.3"])
        self.assertEqual(row["taxid"], 4565)
        self.assertEqual(row["scientific_name"], "Triticum aestivum")
        self.assertEqual((row["query_coverage_pct"], row["identity_pct"]), (100, 100))
        self.assertEqual(row["longest_contiguous_match"], 12)
        self.assertEqual((row["qseq"], row["hseq"], row["midline"]), (SEQUENCE, SEQUENCE, "|" * 12))
        self.assertTrue(row["full_length_exact"])
        self.assertEqual(row["match_orientation"], "same_strand")

    def test_reverse_complement_is_explicit_and_not_presented_as_same_strand(self):
        for strand_fields in ({"hit_strand": "Minus"}, {"hit_strand": None, "hit_frame": -1}):
            data = report(hits=[hit(hsp(hit_from=112, hit_to=101, **strand_fields))])
            row = self.parse(data)["hits"][0]
            self.assertEqual(row["hit_strand"], "-")
            self.assertEqual(row["match_orientation"], "reverse_complement")
            self.assertTrue(row["full_length_exact"])

    def test_gaps_partial_coverage_and_longest_contiguous_match(self):
        data = report(hits=[hit(hsp(query_from=2, query_to=10, hit_from=101, hit_to=109,
                                    qseq="CGTA-CGTAC", hseq="CGTATCG-AT", align_len=10,
                                    identity=7, gaps=2, midline="|||| || | "))])
        row = self.parse(data)["hits"][0]
        self.assertEqual((row["mismatches"], row["gaps"], row["identity_pct"]), (1, 2, 70))
        self.assertEqual(row["query_coverage_pct"], 75)
        self.assertEqual(row["longest_contiguous_match"], 4)
        self.assertFalse(row["full_length_exact"])
        del data["BlastOutput2"][0]["report"]["results"]["search"]["hits"][0]["hsps"][0]["gaps"]
        self.assertEqual(self.parse(data)["hits"][0]["gaps"], 2)

    def test_disjoint_hsps_do_not_become_a_single_full_exact_match(self):
        first = hsp(query_from=1, query_to=6, hit_from=101, hit_to=106, align_len=6, identity=6,
                    qseq="ACGTAC", hseq="ACGTAC", midline="||||||")
        second = hsp(query_from=7, query_to=12, hit_from=501, hit_to=506, align_len=6, identity=6,
                     qseq="GTACGT", hseq="GTACGT", midline="||||||")
        parsed = self.parse(report(hits=[hit(first, second)]))
        self.assertEqual(parsed["subject_hit_count"], 1)
        self.assertEqual(len(parsed["hits"]), 2)
        self.assertTrue(all(row["query_coverage_pct"] == 50 and not row["full_length_exact"] for row in parsed["hits"]))

    def test_no_hits_and_unconfirmed_scope_remain_explicit(self):
        data = report(hits=[], echo=False)
        parsed = self.parse(data)
        self.assertEqual(parsed["hits"], [])
        self.assertEqual(parsed["subject_hit_count"], 0)
        self.assertFalse(parsed["filter_verified"])
        self.assertIsNone(parsed["reported_filter"])
        self.assertTrue(parsed["scope_warning"])

    def test_missing_or_different_subject_taxon_is_never_replaced_with_selected_taxon(self):
        for description in ({"accession": "NC_000001.2"}, {"accession": "NC_000001.2", "taxid": 1234, "sciname": "Other taxon"}):
            row = self.parse(report(hits=[hit(description=[description])]))["hits"][0]
            self.assertEqual(row["taxid"], description.get("taxid"))
            self.assertEqual(row["scientific_name"], description.get("sciname", ""))
            self.assertEqual(row["taxon_match"], "descendant_or_unverified")

    def test_supported_json_envelopes_preserve_one_query_only(self):
        data = report()
        singleton = {"BlastOutput2": data["BlastOutput2"][0]}
        for envelope in (singleton, data["BlastOutput2"], data["BlastOutput2"][0]):
            self.assertEqual(len(self.parse(envelope)["hits"]), 1)
        with self.assertRaises(BlastRemoteError):
            self.parse({"BlastOutput2": data["BlastOutput2"] * 2})

    def test_error_malformed_report_wrong_scope_and_inconsistent_metrics_are_failures(self):
        invalid = [report(database="nt"), report(database=""), report(taxid=4513),
                   {"BlastOutput2": [{"error": "timeout"}]}, {}, []]
        for field in ("search_target", "program"):
            data = report(hits=[])
            del data["BlastOutput2"][0]["report"][field]
            invalid.append(data)
        for change in ({"query_len": 13}, {"hits": None}, {"message": "Search failed"}, {"error": "NCBI error"}):
            data = report(hits=[])
            data["BlastOutput2"][0]["report"]["results"]["search"].update(change)
            invalid.append(data)
        for changes in ({"identity": 13}, {"identity": 11}, {"gaps": 1}, {"query_to": 11},
                        {"evalue": -1}, {"bit_score": float("nan")}, {"align_len": 11}, {"hit_to": 1001}):
            invalid.append(report(hits=[hit(hsp(**changes))]))
        for data in invalid:
            raw = json.dumps(data)
            with self.subTest(data=data), self.assertRaises(BlastRemoteError) as caught:
                self.parse(data)
            self.assertEqual(caught.exception.raw_report, raw)


class BlastClientTests(unittest.TestCase):
    def client(self, responses, **options):
        clock = FakeClock()
        opener = FakeOpener(responses, clock)
        return BlastClient(opener=opener, clock=clock, sleep=clock.sleep, **options), clock, opener

    def test_submit_poll_get_are_post_and_respect_global_and_rid_limits(self):
        client, clock, opener = self.client([
            "RID = FIRST\nRTOE = 90", "Status=WAITING", json.dumps(report()),
            "RID = SECOND\nRTOE = 1", json.dumps(report(database="refseq_rna", genomic=False)),
        ])
        events = []
        query = request(scopes=("genomic", "transcript"))
        result = client.run(query, progress_callback=lambda *args: events.append(args))
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["request"], query.to_dict())
        self.assertEqual(result["fingerprint"], query.fingerprint())
        self.assertEqual([call[0] for call in opener.calls], [0, 90, 150, 160, 220])
        self.assertTrue(all(req.get_method() == "POST" and req.full_url == BLAST_URL for _, req, _ in opener.calls))
        self.assertEqual(opener.params(0)["QUERY"], ">trigger\n" + SEQUENCE + "\n")
        self.assertNotIn("QUERY", opener.params(1))
        self.assertEqual(result["runs"][0]["raw_report"], json.dumps(report()))
        self.assertEqual(result["runs"][0]["params"], opener.params(0))
        self.assertEqual(events[-1][:2], (2, 2))

    def test_each_taxon_and_scope_has_an_independent_report_and_partial_error(self):
        responses = []
        for taxon in (WHEAT, BARLEY):
            for scope in ("genomic", "transcript"):
                responses.append(f"RID = TAX{taxon.taxid}{scope}\nRTOE = 0")
                responses.append("Status=FAILED" if taxon == BARLEY and scope == "genomic" else json.dumps(report(
                    hits=[], taxid=taxon.taxid, genomic=scope == "genomic",
                    database="refseq_genomes" if scope == "genomic" else "refseq_rna")))
        client, _, opener = self.client(responses)
        result = client.run(request(taxa=(WHEAT, BARLEY), scopes=("genomic", "transcript")))
        self.assertEqual(result["status"], "partial")
        self.assertEqual([run["status"] for run in result["runs"]], ["completed", "completed", "failed", "completed"])
        self.assertEqual(result["runs"][2]["raw_format"], "text")
        self.assertEqual(result["runs"][2]["raw_report"], "Status=FAILED")
        self.assertEqual(len(opener.calls), 8)

    def test_cancel_before_submission_and_during_wait_never_sends_an_extra_request(self):
        client, _, opener = self.client([])
        self.assertEqual(client.run(request(), cancel_check=lambda: True)["status"], "cancelled")
        self.assertEqual(opener.calls, [])
        client, clock, opener = self.client(["RID = CANCEL\nRTOE = 1"])
        result = client.run(request(), cancel_check=lambda: clock.now >= 1)
        self.assertEqual(result["status"], "cancelled")
        self.assertEqual(result["runs"][0]["rid"], "CANCEL")
        self.assertEqual(len(opener.calls), 1)
        self.assertLess(clock.now, 2)

    def test_timeout_network_failure_malformed_and_expired_reports_are_not_no_hits(self):
        errors = ["no RID returned", URLError("offline"), HTTPError(BLAST_URL, 503, "busy", {}, io.BytesIO(b"busy"))]
        for response in errors:
            client, _, opener = self.client([response])
            result = client.run(request())
            self.assertEqual(result["status"], "failed")
            self.assertIn("error", result["runs"][0])
            self.assertEqual(len(opener.calls), 1)
        for response in ("Status=UNKNOWN", "not JSON", '{"BlastOutput2":[{"error":"Service failure"}]}'):
            client, _, _ = self.client(["RID = ERROR", response])
            result = client.run(request())
            self.assertEqual(result["status"], "failed")
            self.assertEqual(result["runs"][0]["raw_report"], response)
        client, _, opener = self.client(["RID = TIMEOUT\nRTOE = 120"], max_wait_seconds=100)
        self.assertEqual(client.run(request())["status"], "failed")
        self.assertEqual(len(opener.calls), 1)

    def test_hit_limit_flag_and_status_word_inside_title(self):
        data = report()
        data["BlastOutput2"][0]["report"]["results"]["search"]["hits"][0]["description"][0]["title"] = "Status=READY gene"
        client, _, _ = self.client(["RID = LIMIT", json.dumps(data)])
        result = client.run(request(max_hits=1))
        self.assertEqual(result["status"], "completed")
        self.assertTrue(result["runs"][0]["hit_limit_reached"])
        self.assertTrue(result["runs"][0]["limit_warning"])

    def test_default_clients_share_one_process_gate_and_injected_gate_spaces_clients(self):
        self.assertIs(BlastClient().gate, BlastClient().gate)
        clock = FakeClock()
        gate = _RateGate(clock, clock.sleep)
        opener = FakeOpener(["invalid", "invalid"], clock)
        clients = [BlastClient(opener=opener, clock=clock, sleep=clock.sleep, rate_gate=gate) for _ in range(2)]
        for client in clients:
            client.run(request())
        self.assertEqual([call[0] for call in opener.calls], [0, 10])

    def test_email_is_explicitly_recorded_and_timeouts_must_be_positive(self):
        client, _, opener = self.client(["invalid"], email="maintainer@example.org")
        run = client.run(request())["runs"][0]
        self.assertEqual(run["params"], opener.params(0))
        self.assertEqual(run["params"]["email"], "maintainer@example.org")
        for options in ({"timeout": 0}, {"max_wait_seconds": -1}, {"timeout": float("nan")}):
            with self.assertRaises(ValueError):
                BlastClient(**options)


if __name__ == "__main__":
    unittest.main()
