"""Comparison-only runs retain exact, independently stored design families."""
from copy import deepcopy
from dataclasses import fields, replace
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.structure_ensemble import load_ensemble, save_ensemble
from aptaswitch_core.thermodynamics.types import ModelProvenance
from aptaswitch_studio import saved_switch_design as design
from aptaswitch_studio import saved_switch_library as library
from tests.test_structure_ensemble import fake_analysis, fixture
from tests.ui_helpers import FrenchAppTest as AppTest


def comparison_families_app():
    import streamlit as st
    from aptaswitch_studio.ensemble_ui import render_structure_ensemble
    result = st.session_state.result
    render_structure_ensemble(result, result.candidates[0], allow_calculation=False)


class DesignEnsembleLineageTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.project = self.folder / "Project"
        self.runs = self.project / "calculations" / "designs"
        self.template = fixture(self.project / "calculations")
        self.root = self.project / library.LIBRARY_FOLDER
        self.enterContext(patch.object(design, "default_export_dir", return_value=self.folder / "unused"))
        self.enterContext(patch.object(design, "legacy_export_dir", return_value=self.folder / "legacy"))
        self.enterContext(patch("aptaswitch_studio.ensemble_ui.get_backend",
                                side_effect=AssertionError("No runtime needed")))
        self.enterContext(patch("aptaswitch_studio.ensemble_ui.EnsembleJob",
                                side_effect=AssertionError("No recalculation")))

    def run_result(self, run_id, parent=None):
        result = deepcopy(self.template)
        result.run_id = run_id
        result.output_dir = self.runs / run_id
        if parent is not None:
            result.inverse_extension_analysis = {"source_run_id": parent, "switch_design_repeated": False}
        result.output_dir.mkdir(parents=True)
        (result.output_dir / f"{run_id}_results.json").write_text(json.dumps(result.to_dict()))
        return result

    def native_families(self, result, count=137, seed=23):
        backend = SimpleNamespace(provenance=ModelProvenance(**{
            field.name: result.thermodynamics[field.name] for field in fields(ModelProvenance)}))
        analysis = fake_analysis(result, count, backend=backend, seed=seed)
        save_ensemble(result.output_dir, analysis)
        return analysis

    def save(self, result, *, copy_source=False):
        snapshot = library.candidate_snapshot(result, result.candidates[0], "Retained")
        return library.save_switch(self.root, "Retained", snapshot,
                                   **({"source_output_dir": result.output_dir} if copy_source else {}))

    def test_save_detaches_native_counts_and_seeds_across_comparison_chain(self):
        original = self.run_result("original")
        first = self.run_result("comparison1", "original")
        current = self.run_result("comparison2", "comparison1")
        analyses = [self.native_families(original), self.native_families(first, 3000, 47)]
        before = current.to_dict()
        saved = self.save(current, copy_source=True)
        self.assertEqual(current.to_dict(), before)
        shutil.rmtree(self.runs)
        restored, candidate = design.load_saved_design(saved)
        self.assertEqual(candidate.to_dict(), current.candidates[0].to_dict())
        for analysis in analyses:
            self.assertEqual(load_ensemble(restored.output_dir, analysis["inputs"]), analysis)
        self.assertEqual(saved.snapshot["source"]["inverse_extension_analysis"],
                         current.inverse_extension_analysis)

    def test_existing_saved_candidate_recovers_its_original_run_families(self):
        original = self.run_result("original")
        current = self.run_result("comparison", "original")
        analysis = self.native_families(original)
        saved = self.save(current)
        manifest = saved.manifest_path.read_bytes()
        restored, _ = design.load_saved_design(saved)
        self.assertEqual(load_ensemble(restored.output_dir, analysis["inputs"]), analysis)
        self.assertEqual(saved.manifest_path.read_bytes(), manifest)

    def test_readonly_comparison_explorer_finds_native_parent_without_copying(self):
        original = self.run_result("original")
        current = self.run_result("comparison", "original")
        analysis = self.native_families(original)
        before = {p: p.read_bytes() for p in self.runs.rglob("*") if p.is_file()}
        app = AppTest.from_function(comparison_families_app, default_timeout=30)
        app.session_state.result = current
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.warning)
        count = next(item for item in app.selectbox if item.label == "Échantillonnage enregistré (structures par état)")
        self.assertEqual(count.options, ["137"])
        self.assertEqual(app.table[0].value["Tirages"].tolist(), [137])
        self.assertTrue(any(analysis["model"]["version"] in item.value for item in app.caption))
        self.assertEqual(before, {p: p.read_bytes() for p in self.runs.rglob("*") if p.is_file()})

    def test_parent_cache_still_requires_exact_sequences_and_conditions(self):
        original = self.run_result("original")
        current = self.run_result("comparison", "original")
        original.candidates[0].switch_seq = "C" + original.candidates[0].switch_seq[1:]
        self.native_families(original)
        original.candidates[0] = deepcopy(current.candidates[0])
        original.request = replace(original.request, thermo=replace(original.request.thermo, temperature_c=42))
        self.native_families(original, 3000)
        saved = self.save(current, copy_source=True)
        self.assertFalse(list(saved.path.rglob("structure_ensembles/*.json")))

    def test_corrupt_local_copy_does_not_hide_valid_parent_families(self):
        original = self.run_result("original")
        current = self.run_result("comparison", "original")
        analysis = self.native_families(original)
        broken = save_ensemble(current.output_dir, analysis)
        broken.write_text("{invalid JSON")
        app = AppTest.from_function(comparison_families_app, default_timeout=30)
        app.session_state.result = current
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.warning)
        self.assertEqual(app.table[0].value["Tirages"].tolist(), [137])
        self.assertEqual(broken.read_text(), "{invalid JSON")

    def test_only_explicit_safe_comparison_links_are_followed(self):
        original = self.run_result("original")
        current = self.run_result("comparison")
        for analysis in ({"source_run_id": "original"},
                         {"source_run_id": "original", "switch_design_repeated": True},
                         {"source_run_id": "../original", "switch_design_repeated": False},
                         {"source_run_id": str(original.output_dir), "switch_design_repeated": False}):
            current.inverse_extension_analysis = analysis
            with self.subTest(analysis=analysis):
                self.assertEqual(design.design_ensemble_sources(current), (current.output_dir,))
        current.inverse_extension_analysis = {"source_run_id": "linked", "switch_design_repeated": False}
        (self.runs / "linked").symlink_to(original.output_dir, target_is_directory=True)
        self.assertEqual(design.design_ensemble_sources(current), (current.output_dir,))

    def test_cycles_stop_and_missing_parent_record_does_not_hide_its_cache(self):
        original = self.run_result("original", "comparison")
        current = self.run_result("comparison", "original")
        self.assertEqual(design.design_ensemble_sources(current), (current.output_dir, original.output_dir))
        (original.output_dir / "original_results.json").unlink()
        self.assertEqual(design.design_ensemble_sources(current), (current.output_dir, original.output_dir))


if __name__ == "__main__":
    unittest.main()
