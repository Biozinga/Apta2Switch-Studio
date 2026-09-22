"""Position-aware screening, ensemble metrics and UI for both trigger ends."""

from __future__ import annotations

import itertools
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from test_trigger_extension import FakeNupack, pair_matrix
from aptaswitch_core import trigger_extension as extension
from aptaswitch_core.extension_visuals import (
    APTAMER_COLOR, TRIGGER_COLOR, EXTENSION_COLOR, EXTENSION_3PRIME_COLOR,
    render_extension_structure_svg, export_extension_visuals,
)
from aptaswitch_core.models import ThermoConditions
from aptaswitch_studio.extension_runtime import ExtensionRunRequest, run_extension_and_export


class ExtensionSideTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_split_even_odd_and_single_base(self):
        self.assertEqual(extension.extension_split(6, "5prime"), (6, 0))
        self.assertEqual(extension.extension_split(6, "3prime"), (0, 6))
        self.assertEqual(extension.extension_split(6, "both"), (3, 3))
        self.assertEqual(extension.extension_split(5, "both"), (3, 2))
        self.assertEqual(extension.extension_split(1, "both"), (1, 0))
        with self.assertRaises(ValueError):
            extension.extension_split(3, "invalid")

    def test_all_search_modes_preserve_original_core_at_correct_offset(self):
        for side, prefix, suffix in (("5prime", 3, 0), ("3prime", 0, 3), ("both", 2, 1)):
            with self.subTest(side=side):
                fake = FakeNupack()
                with patch.object(extension, "nupack_import_context") as context:
                    context.return_value.__enter__.return_value = fake
                    result = extension.design_trigger_extensions(
                        aptamer_dna="AC", trigger_dna="GT", target_lengths=[5],
                        thermo=ThermoConditions(material="dna04"), import_path="/test", top=100,
                        extension_side=side,
                    )
                expected = {"d" + "".join(bases[:prefix]) + "GT" + "".join(bases[prefix:])
                            for bases in itertools.product("ACGT", repeat=3)}
                self.assertEqual({seq for _, seq in fake.mfe_sequences}, expected)
                self.assertEqual(len(fake.mfe_sequences), 64)
                self.assertEqual(len(result["candidates"]), 64)
                self.assertTrue(result["mfe_screen"]["exhaustive"])
                self.assertEqual(result["input"]["extension_side"], side)
                for candidate in result["candidates"]:
                    self.assertEqual(len(candidate["extension_5prime"]), prefix)
                    self.assertEqual(len(candidate["extension_3prime"]), suffix)
                    self.assertEqual(candidate["extended_trigger"], candidate["extension_5prime"] + "GT" + candidate["extension_3prime"])
                    positions = candidate["extension_base_probabilities"]
                    self.assertEqual([row["trigger_index"] for row in positions], list(range(1, prefix + 1)) + list(range(prefix + 3, 6)))
                    self.assertEqual([row["side"] for row in positions], ["5prime"] * prefix + ["3prime"] * suffix)

    def test_mfe_projection_and_ensemble_metrics_use_both_extension_segments(self):
        reference = extension.ComplexAnalysisResult("((+))", -3.2, -3.5, .6, .6,
                                                    pair_matrix(4, {(0, 3): .8, (1, 2): .7}))
        for side, prefix, suffix in (("5prime", 2, 0), ("3prime", 0, 2), ("both", 1, 1)):
            with self.subTest(side=side):
                core = [2 + prefix, 3 + prefix]
                added = list(range(2, 2 + prefix)) + list(range(4 + prefix, 6))
                target = "((+" + "." * prefix + "))" + "." * suffix
                self.assertEqual(extension.padded_reference_structure("((+))", 2, side), target)
                matrix = pair_matrix(6, {
                    (0, core[1]): .75, (1, core[0]): .65,
                    (0, added[0]): .03, (1, added[1]): .04,
                    (added[0], core[0]): .02, (added[1], core[1]): .01,
                    (added[0], added[1]): .05,
                })
                candidate = extension.score_mfe_candidate(
                    "AC", extension.MfeResult(target, -3.1), extension.MfeResult("((+))", -3.2),
                    {(0, core[1]), (1, core[0])}, 2, "GT", extension_side=side,
                )
                self.assertEqual(extension.mfe_key(candidate), (0, 0, 0, 0, 0))
                analysis = extension.ComplexAnalysisResult(target, -3.1, -3.45, .55, .55, matrix)
                metrics = extension.ensemble_metrics(candidate, reference, analysis, 2, 2)
                self.assertAlmostEqual(metrics["reference_pair_probability_loss"], .1)
                self.assertAlmostEqual(metrics["extension_pair_probability_sum"], .20)
                self.assertAlmostEqual(metrics["extension_to_aptamer_probability_sum"], .07)
                self.assertAlmostEqual(metrics["extension_to_trigger_core_probability_sum"], .03)
                mapped_matrix = matrix[np.ix_([0, 1] + core, [0, 1] + core)]
                expected_rmse = np.sqrt(np.mean((mapped_matrix - reference.pair_matrix) ** 2))
                self.assertAlmostEqual(metrics["core_pair_probability_rmse"], expected_rmse, places=11)
                pairs = extension._reference_pair_rows(reference, 2, analysis, prefix)
                self.assertEqual([row["candidate_probability"] for row in pairs], [.75, .65])

    def test_mfe_filter_counts_new_pairs_in_suffix_and_between_ends(self):
        for side, structure in (("3prime", "..+(.)."), ("both", "..+(..)")):
            with self.subTest(side=side):
                candidate = extension.score_mfe_candidate(
                    "AC", extension.MfeResult(structure, -1), extension.MfeResult("..+..", 0),
                    set(), 2, "GT", extension_side=side,
                )
                self.assertEqual(candidate.extension_mfe_pairs, 1)
                self.assertEqual(candidate.new_core_mfe_pairs, 0)

    def test_structures_color_prefix_and_suffix_at_their_actual_positions(self):
        from aptaswitch_core.preview import render_structure_svg
        with patch("aptaswitch_core.preview.render_structure_svg", wraps=render_structure_svg) as draw:
            svg = render_extension_structure_svg("AC", "AGTC", "((+. )).".replace(" ", ""), extension_length=2, extension_3prime_length=1)
        self.assertEqual(draw.call_args.args[2], [APTAMER_COLOR] * 2 + [EXTENSION_COLOR] + [TRIGGER_COLOR] * 2 + [EXTENSION_3PRIME_COLOR])
        self.assertIn("5′ 1 nt / 3′ 1 nt", svg)
        with self.assertRaises(ValueError):
            render_extension_structure_svg("AC", "AGTC", "..+....", extension_length=1, extension_3prime_length=2)

    def test_runtime_and_exports_keep_chosen_side_and_actual_sequences(self):
        fake = FakeNupack()
        with tempfile.TemporaryDirectory() as tmp, patch.object(extension, "nupack_import_context") as context:
            context.return_value.__enter__.return_value = fake
            request = ExtensionRunRequest("AC", "GT", (4,), ThermoConditions(material="dna04"), "/test", Path(tmp), top=2, extension_side="both")
            result, paths = run_extension_and_export(request)
            self.assertEqual(result["request"]["extension_side"], "both")
            self.assertEqual(result["input"]["extension_side"], "both")
            self.assertIn("extension_3prime", Path(paths["csv"]).read_text())
            self.assertIn("Moitié 5′ / moitié 3′", Path(paths["log"]).read_text())
            for key in ("svg_candidate_4", "svg_comparison"):
                self.assertIn("5′ 1 nt / 3′ 1 nt", Path(paths[key]).read_text())


@unittest.skipUnless(os.environ.get("APTASWITCH_TEST_NUPACK"), "Set APTASWITCH_TEST_NUPACK for real NUPACK checks")
class RealNupackExtensionSidesTests(unittest.TestCase):
    def test_real_dna_and_rna_searches_at_3prime_and_both_ends(self):
        for side in ("3prime", "both"):
            for material, trigger in (("dna04", "ACGT"), ("rna06", "ACGU")):
                with self.subTest(side=side, material=material):
                    result = extension.design_trigger_extensions(
                        aptamer_dna="ACGTACGT", trigger_dna=trigger, target_lengths=[6],
                        thermo=ThermoConditions(material=material, temperature_c=37),
                        import_path=os.environ["APTASWITCH_TEST_NUPACK"], top=3, extension_side=side,
                    )
                    self.assertEqual(result["status"], "completed")
                    self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 16)
                    self.assertTrue(result["candidates"])
                    for candidate in result["candidates"]:
                        self.assertEqual(candidate["extended_trigger"], candidate["extension_5prime"] + trigger + candidate["extension_3prime"])
                        self.assertEqual(len(candidate["extension_3prime"]), 2 if side == "3prime" else 1)
                        self.assertEqual(len(candidate["reference_pairs"]), len(result["reference"]["mfe_pairs"]))
                    with tempfile.TemporaryDirectory() as tmp:
                        paths = export_extension_visuals(result, tmp)
                        self.assertTrue(Path(paths["svg_candidate_6"]).is_file())


if __name__ == "__main__":
    unittest.main()
