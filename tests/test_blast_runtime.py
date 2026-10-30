"""BLAST execution is explicit; archived failures, provenance and rows survive moves."""

from __future__ import annotations

import copy
import csv
import io
import json
import shutil
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch
from urllib.parse import parse_qs

from aptaswitch_core.blast_specificity import BlastClient, BlastScreenRequest, Taxon
from aptaswitch_studio.blast_runtime import BlastScreenJob, load_blast_result, run_blast_and_export


def screen_request():
    return BlastScreenRequest("ACGTACGTACGTACGTACGTACGT", (Taxon(9606, "Homo sapiens"),),
                              scopes=("genomic", "transcript"), max_hits=100)


def client_result(*, status="completed"):
    hit = {
        "accession": "NC_TEST.1", "accessions": ["NC_TEST.1"], "title": "Recorded reference",
        "taxid": 9606, "scientific_name": "Homo sapiens", "subject_id": "ref|NC_TEST.1|",
        "hit_number": 1, "hsp_number": 1, "query_from": 1, "query_to": 24,
        "hit_from": 200, "hit_to": 177, "query_strand": "Plus", "hit_strand": "Minus",
        "alignment_length": 24, "identities": 23, "identity_pct": 100 * 23 / 24,
        "query_coverage_pct": 100., "mismatches": 1, "gaps": 0,
        "evalue": 0.002, "bit_score": 40.5, "score": 22,
        "longest_contiguous_match": 12, "qseq": screen_request().sequence,
        "hseq": "ACGTACGTACGTTCGTACGTACGT", "midline": "|||||||||||| |||||||||||",
    }
    runs = []
    for scope in screen_request().scopes:
        runs.append({
            "taxon": {"taxid": 9606, "name": "Homo sapiens"}, "scope": scope,
            "status": "completed", "rid": f"RID_{scope.upper()}", "database": "refseq_rna" if scope == "transcript" else "refseq_genomes",
            "query_length": 24, "queried_at": "2026-09-23T12:00:00+00:00",
            "reported_database": "recorded database", "version": "BLASTN test",
            "params": {"PROGRAM": "blastn", "WORD_SIZE": "7", "HITLIST_SIZE": "100", "ENTREZ_QUERY": "txid9606[ORGN]"},
            "subject_hit_count": 1 if scope == "genomic" else 0,
            "hits": [hit] if scope == "genomic" else [],
            "raw_format": "json", "raw_report": json.dumps({"fixture": scope, "hits": [hit] if scope == "genomic" else []}),
        })
    if status == "partial":
        runs[1].update(status="failed", error="Remote result expired", raw_report="{invalid remote response")
    return {"status": status, "runs": runs}


class BlastRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.output = self.root / "design" / "calculations"
        self.request = screen_request()

    def client(self, result=None):
        return Mock(run=Mock(return_value=client_result() if result is None else result))

    def test_constructor_never_creates_client_sends_query_or_writes_files(self):
        factory = Mock()
        job = BlastScreenJob(self.request, self.output, client_factory=factory)
        factory.assert_not_called()
        self.assertFalse(self.output.exists())
        self.assertFalse(job.consume_once())

    def test_explicit_execution_saves_exact_provenance_raw_reports_and_hsp_csv(self):
        source = client_result()
        client = self.client(source)
        progress = Mock()
        result, exports = run_blast_and_export(self.request, self.output, client=client, progress_callback=progress)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["request"], self.request.normalized().to_dict())
        self.assertEqual(result["fingerprint"], self.request.fingerprint())
        self.assertEqual(result["runs"][0]["params"], source["runs"][0]["params"])
        self.assertEqual(result["runs"][0]["version"], "BLASTN test")
        self.assertTrue(Path(result["output_dir"]).is_relative_to(self.output / "blast"))
        self.assertIn("started_at", result)
        self.assertIn("completed_at", result)
        for index, run in enumerate(result["runs"], 1):
            self.assertNotIn("raw_report", run)
            self.assertEqual(run["raw_report_file"], f"search_{index:03d}.json")
            self.assertEqual(Path(exports[f"raw_{index}"]).read_text(), source["runs"][index - 1]["raw_report"])
        self.assertIn("raw_report", source["runs"][0])
        with Path(exports["csv"]).open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["accession"], "NC_TEST.1")
        self.assertEqual(rows[0]["hit_strand"], "Minus")
        self.assertEqual(rows[1]["reported_hsp_count"], "0")
        self.assertEqual(rows[1]["search_status"], "completed")
        self.assertEqual(rows[1]["search_error"], "")
        with zipfile.ZipFile(exports["zip"]) as archive:
            self.assertTrue({"results.json", "request.json", "hits.csv", "blast.log", "search_001.json", "search_002.json"} <= set(archive.namelist()))
        client.run.assert_called_once()
        self.assertGreaterEqual(progress.call_count, 2)

    def test_network_failure_is_archived_as_failure_never_zero_hit_success(self):
        client = self.client()
        client.run.side_effect = TimeoutError("Network timeout")
        result, exports = run_blast_and_export(self.request, self.output, client=client)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["error"], "Network timeout")
        self.assertEqual([run["status"] for run in result["runs"]], ["failed", "failed"])
        with Path(exports["csv"]).open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertTrue(all(row["search_status"] == "failed" and row["search_error"] == "Network timeout" for row in rows))
        self.assertEqual(load_blast_result(Path(exports["json"]))[0]["status"], "failed")

    def test_real_client_protocol_replay_integrates_with_archives_without_network(self):
        now = [0.]
        submitted = {}
        calls = []

        def sleep(seconds):
            now[0] += seconds

        def opener(request, *, timeout):
            params = {key: values[0] for key, values in parse_qs(request.data.decode()).items()}
            calls.append(params)
            if params["CMD"] == "Put":
                rid = f"TEST{len(submitted) + 1}"
                submitted[rid] = params
                return io.BytesIO(f"RID = {rid}\nRTOE = 1\n".encode())
            submitted_params = submitted[params["RID"]]
            raw = {"BlastOutput2": [{"report": {
                "program": "blastn", "version": "BLASTN fixture", "search_target": {"db": submitted_params["DATABASE"]},
                "params": {"entrez_query": submitted_params["ENTREZ_QUERY"]},
                "results": {"search": {"query_len": 24, "hits": [], "stat": {"db_num": 100, "db_len": 50000}}},
            }}]}
            return io.BytesIO(json.dumps(raw).encode())

        client = BlastClient(opener=opener, clock=lambda: now[0], sleep=sleep)
        result, exports = run_blast_and_export(self.request, self.output, client=client)
        self.assertEqual(result["status"], "completed")
        self.assertEqual([call["CMD"] for call in calls], ["Put", "Get", "Put", "Get"])
        self.assertEqual([row["database"] for row in result["runs"]], ["refseq_genomes", "refseq_rna"])
        self.assertTrue(all(row["filter_verified"] and row["subject_hit_count"] == 0 for row in result["runs"]))
        self.assertTrue(all(row["blast_version"] == "BLASTN fixture" for row in result["runs"]))
        self.assertEqual(load_blast_result(Path(exports["json"]))[0]["runs"], result["runs"])

    def test_partial_report_preserves_failed_search_raw_text_and_successful_hits(self):
        result, exports = run_blast_and_export(self.request, self.output, client=self.client(client_result(status="partial")))
        self.assertEqual(result["status"], "partial")
        self.assertEqual(len(result["runs"][0]["hits"]), 1)
        self.assertEqual(result["runs"][1]["status"], "failed")
        self.assertEqual(Path(exports["raw_2"]).read_text(), "{invalid remote response")
        self.assertEqual(load_blast_result(Path(exports["json"]))[0]["runs"][1]["error"], "Remote result expired")

    def test_incomplete_client_success_is_rejected_as_failure(self):
        for bad in ({"status": "completed"}, {"status": "completed", "runs": []},
                    {"status": "completed", "runs": [dict(client_result()["runs"][0], status="failed")]}):
            with self.subTest(result=bad):
                result, _ = run_blast_and_export(self.request, self.output, client=self.client(bad))
                self.assertEqual(result["status"], "failed")

    def test_cancel_before_start_never_constructs_or_calls_client(self):
        factory = Mock()
        job = BlastScreenJob(self.request, self.output, client_factory=factory)
        job.cancel()
        job.start()
        job._thread.join(3)
        self.assertFalse(job._thread.is_alive())
        factory.assert_not_called()
        saved = job.snapshot()
        self.assertEqual(saved.result["status"], "cancelled")
        self.assertEqual(saved.state, "done")
        self.assertIsNotNone(saved.exports["json"])
        self.assertTrue(job.consume_once())
        self.assertFalse(job.consume_once())

    def test_running_job_cancellation_preserves_finished_search_and_stops_client_wait(self):
        entered = threading.Event()

        class WaitingClient:
            def run(inner_self, request, *, progress_callback, cancel_check):
                progress_callback(1, 2, "Waiting for remote report")
                entered.set()
                for _ in range(200):
                    if cancel_check():
                        result = client_result()
                        result["status"] = "cancelled"
                        result["runs"][1].update(status="cancelled", hits=[])
                        return result
                    threading.Event().wait(.01)
                raise AssertionError("Cancellation did not reach the client")

        job = BlastScreenJob(self.request, self.output, client_factory=WaitingClient)
        job.start()
        self.assertTrue(entered.wait(2))
        self.assertEqual(job.snapshot().state, "running")
        job.cancel()
        job._thread.join(3)
        self.assertFalse(job._thread.is_alive())
        snapshot = job.snapshot()
        self.assertEqual(snapshot.result["status"], "cancelled")
        self.assertEqual(snapshot.result["runs"][0]["status"], "completed")
        self.assertEqual(snapshot.result["runs"][1]["status"], "cancelled")
        self.assertEqual(snapshot.current, 1)

    def test_loader_relocates_files_and_keeps_rids_without_network(self):
        result, exports = run_blast_and_export(self.request, self.output, client=self.client())
        old = Path(result["output_dir"])
        destination = self.root / "moved" / old.name
        destination.parent.mkdir()
        shutil.move(str(old), destination)
        with patch("aptaswitch_studio.blast_runtime.BlastClient", side_effect=AssertionError("Offline restore only")):
            restored, relocated = load_blast_result(destination)
        self.assertEqual(restored["output_dir"], str(destination))
        self.assertEqual(restored["runs"], result["runs"])
        self.assertTrue(all(Path(value).parent == destination for value in relocated.values()))
        self.assertEqual(set(relocated), set(exports))

    def test_loader_rejects_changed_fingerprint_and_external_raw_files(self):
        result, exports = run_blast_and_export(self.request, self.output, client=self.client())
        path = Path(exports["json"])
        tampered = copy.deepcopy(result)
        tampered["request"]["sequence"] = "T" * 24
        path.write_text(json.dumps(tampered))
        with self.assertRaises(ValueError):
            load_blast_result(path)
        external = self.root / "external.json"
        external.write_text("external")
        result["runs"][0]["raw_report_file"] = "../../../../native.json"
        result["runs"][1]["raw_report_file"] = str(external)
        path.write_text(json.dumps(result))
        _, restored = load_blast_result(path)
        self.assertFalse(any(key.startswith("raw_") for key in restored))

    def test_successive_screens_keep_previous_report_unchanged(self):
        _, first = run_blast_and_export(self.request, self.output, client=self.client())
        previous = Path(first["json"]).read_bytes()
        _, second = run_blast_and_export(self.request, self.output, client=self.client())
        self.assertNotEqual(first["json"], second["json"])
        self.assertEqual(Path(first["json"]).read_bytes(), previous)


if __name__ == "__main__":
    unittest.main()
