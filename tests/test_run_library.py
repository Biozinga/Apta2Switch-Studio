"""Opening run archives must preserve recorded results without new calculations."""

from __future__ import annotations

import csv
import json
import shutil
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.exports import export_run
from aptaswitch_core.models import DesignCandidate, RunRequest, ScoringWeights, SequenceInput
from aptaswitch_studio.run_library import discover_runs, load_run
from tests.design_fixtures import design_result_fixture


FIXTURE = Path(__file__).parent / "fixtures" / "extension_nupack_smoke.json"


class RunLibraryTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()

    def design(self, run_id="saved_design"):
        request = RunRequest(
            sequences=SequenceInput("ACGTACGT", "ACGT" * 6, "Mon système"),
            architecture=default_architecture_for_trigger(24),
            scoring=ScoringWeights(on_yield=1.2, leak=0.1, defect=7, rbs=4),
            output_dir=self.root / "designs", trials=2,
        )
        result = design_result_fixture(request)
        result.run_id = run_id
        result.output_dir = request.output_dir / run_id
        result.candidates[0].score = -103.759
        result.candidates.append(DesignCandidate.excluded_for_stop(
            trial=3, switch_seq=result.candidates[0].switch_seq,
            trigger_seq=result.candidates[0].trigger_seq, first_stop=2, reporter_label="sfGFP",
        ))
        exports = export_run(result)
        return result, Path(exports["json"])

    def extension(self, run_id="saved_extension"):
        result = json.loads(FIXTURE.read_text())
        result["run_id"] = run_id
        result["output_dir"] = "/old/unavailable/location"
        directory = self.root / "extensions" / run_id
        directory.mkdir(parents=True)
        path = directory / f"{run_id}_results.json"
        path.write_text(json.dumps(result))
        return result, path

    def write(self, path, data):
        path.write_text(json.dumps(data))

    def test_design_roundtrip_preserves_order_scores_weights_and_exclusions(self):
        original, path = self.design()
        with patch("aptaswitch_core.nupack_engine.run_nupack_engine", side_effect=AssertionError("Must not compute")):
            loaded = load_run(path)
        self.assertEqual(loaded.kind, "design")
        self.assertEqual(loaded.source, path)
        self.assertEqual(loaded.result.to_dict(), original.to_dict())
        self.assertLess(loaded.result.candidates[0].score, 0)
        self.assertIsNone(loaded.result.candidates[-1].score)
        self.assertTrue(loaded.result.candidates[-1].excluded)
        self.assertEqual(loaded.result.request.scoring.on_yield, 1.2)
        self.assertEqual(loaded.exports["json"], str(path))

    def test_open_directory_and_manifest_and_extension_directory(self):
        original, path = self.design()
        for location in (path.parent, path.parent / ".aptaswitch.json"):
            self.assertEqual(load_run(location).result.to_dict(), original.to_dict())
        _, extension_path = self.extension()
        self.assertEqual(load_run(extension_path.parent).kind, "extension")

    def test_moved_design_relocates_all_export_paths_and_ignores_missing_exports(self):
        original, path = self.design()
        destination = self.root / "moved" / original.run_id
        shutil.move(str(path.parent), destination)
        moved_json = destination / path.name
        (destination / f"{original.run_id}_ranked.csv").unlink()
        loaded = load_run(destination)
        self.assertEqual(loaded.source, moved_json)
        self.assertEqual(loaded.result.output_dir, destination)
        self.assertEqual(loaded.result.request.output_dir, destination.parent)
        self.assertNotIn("csv", loaded.exports)
        self.assertIn("project", loaded.exports)
        for exported in loaded.exports.values():
            self.assertEqual(Path(exported).parent, destination)
            self.assertTrue(Path(exported).is_file())
        self.assertEqual(json.loads(moved_json.read_text())["output_dir"], str(original.output_dir))

    def test_old_architecture_does_not_gain_a_new_rbs_prefix(self):
        _, path = self.design()
        data = json.loads(path.read_text())
        architecture = data["request"]["architecture"]
        del architecture["rbs_prefix_length"]
        architecture["rbs_loop"] = "UAGAGGAGAAC"
        architecture["frame_linker"] = "CC"
        for row in data["results"]:
            row.pop("excluded")
            row.pop("exclusion_reason")
        data["results"] = data["results"][:2]
        self.write(path, data)
        loaded = load_run(path)
        self.assertEqual(loaded.result.request.architecture.rbs_prefix_length, 0)
        self.assertEqual(loaded.result.request.architecture.rbs_loop, "UAGAGGAGAAC")
        self.assertEqual(loaded.result.request.architecture.frame_linker, "CC")
        self.assertEqual(loaded.result.candidates[0].score, data["results"][0]["score"])
        self.assertFalse(loaded.result.candidates[0].excluded)

    def test_extension_preserves_exact_data_and_old_outputs_without_request(self):
        original, path = self.extension()
        del original["request"]
        self.write(path, original)
        loaded = load_run(path)
        self.assertEqual(loaded.kind, "extension")
        expected = dict(original, output_dir=str(path.parent))
        self.assertEqual(loaded.result, expected)
        self.assertEqual(loaded.exports, {"json": str(path)})
        self.assertEqual(discover_runs([self.root])[0].molecule_name, "Système personnalisé")

    def test_legacy_stop_policy_is_disabled_without_reranking_evaluated_candidates(self):
        _, path = self.design()
        data = json.loads(path.read_text())
        data["results"] = data["results"][:2]
        data["results"][0]["has_stop_codon"] = True
        data["results"][0]["first_stop"] = 2
        for policy in (None, False, True):
            with self.subTest(recorded_policy=policy):
                if policy is None:
                    data["request"].pop("exclude_stop_candidates", None)
                else:
                    data["request"]["exclude_stop_candidates"] = policy
                self.write(path, data)
                loaded = load_run(path)
                self.assertEqual(loaded.result.request.exclude_stop_candidates,
                                 False if policy is None else policy)
                self.assertEqual([candidate.to_dict() for candidate in loaded.result.candidates], data["results"])
                self.assertTrue(loaded.result.candidates[0].has_stop_codon)
                self.assertFalse(loaded.result.candidates[0].excluded)

    def test_linearity_policy_preserves_archives_without_recomputing_or_reranking(self):
        _, path = self.design()
        data = json.loads(path.read_text())
        data["results"][0]["has_bad_rbs_linker"] = True
        for policy in (None, False, True):
            with self.subTest(recorded_policy=policy):
                if policy is None:
                    data["request"].pop("exclude_non_linear_candidates", None)
                else:
                    data["request"]["exclude_non_linear_candidates"] = policy
                self.write(path, data)
                before = path.read_bytes()
                with patch("aptaswitch_core.nupack_engine.run_nupack_engine", side_effect=AssertionError("Must not compute")):
                    loaded = load_run(path)
                self.assertEqual(loaded.result.request.exclude_non_linear_candidates,
                                 False if policy is None else policy)
                self.assertEqual([candidate.to_dict() for candidate in loaded.result.candidates], data["results"])
                self.assertFalse(loaded.result.candidates[0].excluded)
                self.assertEqual(path.read_bytes(), before)

    def test_linearity_policy_and_partially_screened_exclusion_survive_exports(self):
        result, _ = self.design()
        candidate = result.candidates[0]
        result.candidates = [DesignCandidate.excluded_for_non_linear(
            trial=candidate.trial, switch_seq=candidate.switch_seq,
            trigger_seq=candidate.trigger_seq, structure_on=candidate.structure_on,
            delta_g_on=-1.75, has_stop_codon=False, first_stop=None,
            reporter_label=candidate.reporter_label,
        )]
        reason = result.candidates[0].exclusion_reason
        for policy in (True, False):
            with self.subTest(recorded_policy=policy):
                result.request = replace(result.request, exclude_non_linear_candidates=policy)
                exports = export_run(result)
                self.assertNotIn("svg", exports)
                for key in ("json", "project"):
                    saved = json.loads(Path(exports[key]).read_text())
                    self.assertEqual(saved["request"]["exclude_non_linear_candidates"], policy)
                journal = Path(exports["log"]).read_text()
                self.assertIn(f"Exclure les candidats au RBS–linker non linéaire : {policy}", journal)
                self.assertIn("Candidats exclus (analyses restantes ignorées) : 1", journal)
                with Path(exports["csv"]).open(newline="") as handle:
                    row, = list(csv.DictReader(handle))
                self.assertEqual(row["Rank"], "")
                self.assertEqual(row["Score_selection"], "")
                self.assertEqual(row["Status"], "excluded")
                self.assertEqual(row["Exclusion_reason"], reason)
                self.assertEqual(row["STOP_codon"], "no")
                self.assertEqual(row["dG_ON_kcal_mol"], "-1.75")
                self.assertEqual(row["dG_RBS_linker_kcal_mol"], "")
                loaded = load_run(exports["project"])
                self.assertEqual(loaded.result.to_dict(), result.to_dict())

    def test_extension_exports_are_recovered_from_actual_siblings(self):
        original, path = self.extension()
        names = ["svg_reference.svg", "svg_candidate_24.svg", "svg_candidate_24.png", "png_comparison.png",
                 f"{original['run_id']}_ranked.csv", f"{original['run_id']}_complete.zip",
                 f"{original['run_id']}.log"]
        for name in names:
            (path.parent / name).write_text("existing export")
        (path.parent / "other_data.csv").write_text("unrelated")
        loaded = load_run(path)
        self.assertEqual(set(loaded.exports), {"json", "csv", "zip", "log", "svg_reference", "svg_candidate_24", "png_candidate_24", "png_comparison"})

    def test_discovery_deduplicates_roots_and_manifest_sorts_dates_with_timezones(self):
        design, design_path = self.design()
        extension, extension_path = self.extension()
        design_data = design.to_dict()
        design_data["started_at"] = "2026-09-22T15:00:00+00:00"
        extension["started_at"] = "2026-09-22T16:00:00+02:00"
        self.write(design_path, design_data)
        self.write(extension_path, extension)
        (self.root / "incomplete_results.json").write_text("{")
        (self.root / "irrelevant.json").write_text("{}")
        (self.root / "recursive_link").symlink_to(self.root, target_is_directory=True)
        rows = discover_runs([self.root, design_path.parent, design_path])
        self.assertEqual([row.kind for row in rows], ["design", "extension"])
        self.assertEqual(rows[0].candidate_count, 3)
        self.assertEqual(rows[0].molecule_name, "Mon système")
        self.assertEqual(len({row.path for row in rows}), 2)

    def test_discovery_accepts_manifest_with_nonstandard_result_filename(self):
        _, path = self.design()
        renamed = path.with_name("my_results_backup.json")
        path.rename(renamed)
        manifest = path.parent / ".aptaswitch.json"
        data = json.loads(manifest.read_text())
        data["results_file"] = renamed.name
        self.write(manifest, data)
        self.assertEqual([row.path for row in discover_runs([self.root])], [renamed])

    def test_mock_runs_are_rejected_and_hidden_from_library(self):
        _, path = self.design()
        data = json.loads(path.read_text())
        data["request"]["engine"] = "mock"
        self.write(path, data)
        with self.assertRaisesRegex(ValueError, "démonstration"):
            load_run(path)
        self.assertEqual(discover_runs([self.root]), [])

    def test_bad_files_have_readable_errors(self):
        path = self.root / "broken_results.json"
        for content, message in (("{", "JSON valide"), ("[]", "objet JSON"), ("{}", "n’est pas un résultat")):
            path.write_text(content)
            with self.subTest(content=content), self.assertRaisesRegex(ValueError, message):
                load_run(path)
        with self.assertRaisesRegex(ValueError, "Impossible de lire"):
            load_run(self.root / "missing.json")
        empty = self.root / "empty"
        empty.mkdir()
        with self.assertRaisesRegex(ValueError, "ne contient pas de run"):
            load_run(empty)

    def test_manifest_cannot_escape_its_run_folder(self):
        _, path = self.design()
        manifest = path.parent / ".aptaswitch.json"
        for target in ("../outside.json", "/absolute.json", ".aptaswitch.json"):
            with self.subTest(target=target):
                self.write(manifest, {"results_file": target})
                with self.assertRaisesRegex(ValueError, "même dossier"):
                    load_run(manifest)
        alias = path.parent / "external.json"
        outside = self.root / "outside.json"
        outside.write_text(path.read_text())
        alias.symlink_to(outside)
        self.write(manifest, {"results_file": alias.name})
        with self.assertRaisesRegex(ValueError, "même dossier"):
            load_run(manifest)

    def test_directory_with_multiple_runs_requests_explicit_file(self):
        _, path = self.extension()
        path.with_name("second_results.json").write_text(path.read_text())
        with self.assertRaisesRegex(ValueError, "plusieurs runs"):
            load_run(path.parent)

    def test_malformed_candidate_data_rejected_before_ui_render(self):
        _, design_path = self.design()
        original_design = json.loads(design_path.read_text())
        for field, value in (("score", None), ("trial", "one"), ("switch_seq", "NNNN"), ("excluded", "false")):
            data = json.loads(json.dumps(original_design))
            data["results"][0][field] = value
            self.write(design_path, data)
            with self.subTest(field=field), self.assertRaises(ValueError):
                load_run(design_path)
        original_extension, path = self.extension()
        for field, value in (("mfe_structure", "....+...."), ("final_score", None),
                             ("extended_trigger", "AACGT"), ("nucleotide_probabilities", [0.5])):
            data = json.loads(json.dumps(original_extension))
            data["candidates"][0][field] = value
            self.write(path, data)
            with self.subTest(field=field), self.assertRaises(ValueError):
                load_run(path)
        self.assertEqual(discover_runs([self.root]), [])

    def test_cancelled_extension_without_reference_or_candidates_can_be_reopened(self):
        data, path = self.extension()
        data.update(status="cancelled", reference=None, candidates=[])
        self.write(path, data)
        self.assertEqual(load_run(path).result["status"], "cancelled")


if __name__ == "__main__":
    unittest.main()
