from __future__ import annotations

import re
import sys
import unittest
import xml.etree.ElementTree as ET
from math import hypot
from pathlib import Path
from dataclasses import replace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.architecture import default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.sequences import reverse_complement_rna
from aptaswitch_core.preview import (
    PLACEHOLDER_BASE,
    build_switch_model,
    linear_sequence_svg,
    parse_structure,
    premature_stop_positions,
    radial_layout,
    render_structure_svg,
    switch_preview_svg,
)

TRIGGER = "CTCTCCTTACGCCACCCACACCCGACGTAC"  # Standard 30 nt trigger region
CLASSIC_LINKER = "AACCUGGCGGCAGCGCAAAAG"


class SwitchModelTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))
        self.arch = default_architecture_for_trigger(len(TRIGGER))
        self.model = build_switch_model(self.arch, TRIGGER)

    def test_off_structure_matches_sequence_length(self):
        self.assertEqual(len(self.model.off_structure), len(self.model.switch_sequence))

    def test_off_structure_is_balanced(self):
        pairs, n, cuts = parse_structure(self.model.off_structure)
        self.assertEqual(n, len(self.model.switch_sequence))
        self.assertFalse(cuts)
        self.assertTrue(pairs)

    def test_on_structure_has_strand_break(self):
        pairs, n, cuts = parse_structure(self.model.on_structure)
        self.assertEqual(len(cuts), 1)
        self.assertEqual(n, len(self.model.on_sequence))
        self.assertEqual(len(pairs), self.arch.trigger_region_length)

    def test_optimizable_domains_are_placeholders(self):
        optimizable = [d for d in self.model.domains if d.optimizable]
        self.assertTrue(optimizable)
        for domain in optimizable:
            self.assertEqual(set(domain.sequence), {PLACEHOLDER_BASE})

    def test_frame_linker_is_optimized_and_unpaired_in_both_target_previews(self):
        architecture = replace(self.arch, binding_layout="legacy_lower_stem", stem_length=18,
                               upper_stem2_length=4, frame_linker=None)
        model = build_switch_model(architecture, TRIGGER)
        linker_index = next(i for i, domain in enumerate(model.domains) if domain.label == "Linker")
        linker = model.domains[linker_index]
        self.assertEqual(linker.sequence, "XX")
        self.assertTrue(linker.optimizable)
        linker_start = sum(len(domain.sequence) for domain in model.domains[:linker_index])
        for structure in (model.off_structure, model.on_structure.split("+")[-1]):
            self.assertEqual(structure[linker_start:linker_start + 2], "..")

        fixed = build_switch_model(replace(self.arch, frame_linker="CC"), TRIGGER)
        fixed_linker = next(domain for domain in fixed.domains if domain.label == "Linker")
        self.assertEqual(fixed_linker.sequence, "CC")
        self.assertFalse(fixed_linker.optimizable)

    def test_linear_preview_shows_free_linker_then_actual_candidate_bases(self):
        architecture = replace(self.arch, binding_layout="legacy_lower_stem", stem_length=18,
                               upper_stem2_length=4, frame_linker=None)
        model = build_switch_model(architecture, TRIGGER)
        linker_index = next(i for i, domain in enumerate(model.domains) if domain.label == "Linker")
        linker_start = sum(len(domain.sequence) for domain in model.domains[:linker_index])
        sequence = model.switch_sequence.replace("X", "A")
        sequence = sequence[:linker_start] + "GU" + sequence[linker_start + 2:]
        for actual, expected_linker in ((None, "XX"), (sequence, "GU")):
            with self.subTest(candidate=actual is not None):
                svg, _, _ = linear_sequence_svg(architecture, TRIGGER, switch_sequence=actual)
                root = ET.fromstring(svg)
                bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
                switch = "".join(bases)[len(model.trigger_sequence):]
                self.assertEqual(switch[linker_start:linker_start + 2], expected_linker)
                self.assertEqual(switch, actual or model.switch_sequence)

    def test_standard_30_has_six_upper_pairs_aug_bulge_and_nine_lower_pairs(self):
        domains = self.model.domains
        labels = [domain.label for domain in domains]
        aug_index = labels.index("AUG")
        self.assertEqual(labels[aug_index + 1:], ["Stem R", "Linker", "sfGFP"])
        self.assertNotIn("U2 L", labels)
        self.assertNotIn("U2 R", labels)
        linker = next(domain for domain in domains if domain.label == "Linker")
        self.assertEqual(linker.sequence, CLASSIC_LINKER)
        self.assertFalse(linker.optimizable)
        self.assertEqual(next(domain.sequence for domain in domains if domain.label == "RBS"), "AGAGGAGA")
        aug_end = sum(len(domain.sequence) for domain in domains[:aug_index + 1])
        reporter_start = len(self.model.switch_sequence) - len(domains[-1].sequence)
        self.assertEqual(reporter_start - aug_end, 30)
        self.assertEqual(self.model.switch_sequence[aug_end:reporter_start], "X" * 9 + CLASSIC_LINKER)
        self.assertEqual(self.model.off_structure[aug_end:reporter_start], ")" * 9 + "." * 21)
        self.assertEqual(self.model.on_structure.split("+")[-1][aug_end:reporter_start], "." * 30)
        self.assertEqual(self.model.switch_sequence[aug_end - 3:aug_end], "AUG")
        self.assertEqual(self.model.off_structure[aug_end - 3:aug_end], "...")
        self.assertEqual(self.model.off_structure[aug_end - 9:aug_end - 3], ")" * 6)
        self.assertEqual(len(parse_structure(self.model.off_structure)[0]), 15)
        svg, _, _ = linear_sequence_svg(self.arch, TRIGGER)
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
        self.assertEqual("".join(bases)[len(TRIGGER):], self.model.switch_sequence)
        self.assertNotIn("U2 L", svg)
        self.assertNotIn("U2 R", svg)
        self.assertIn("Linker", svg)

    def test_explicit_empty_linker_and_saved_none_linker_do_not_inherit_classic_sequence(self):
        for linker in ("", None):
            with self.subTest(frame_linker=linker):
                architecture = replace(self.arch, frame_linker=linker)
                model = build_switch_model(architecture, TRIGGER)
                self.assertNotIn("Linker", [domain.label for domain in model.domains])
                self.assertNotIn(CLASSIC_LINKER, model.switch_sequence)
                self.assertEqual(model.switch_sequence[62:], model.domains[-1].sequence)
                self.assertEqual(len(model.switch_sequence), len(self.model.switch_sequence) - 21)

    def test_manual_linker_replaces_the_classic_sequence_in_both_preview_forms(self):
        architecture = replace(self.arch, frame_linker="CCX")
        model = build_switch_model(architecture, TRIGGER)
        self.assertEqual(model.switch_sequence[62:65], "CCX")
        self.assertNotIn(CLASSIC_LINKER, model.switch_sequence)
        self.assertEqual([domain.sequence for domain in model.domains if domain.label == "Linker"], ["CCX"])
        for structure in (model.off_structure, model.on_structure.split("+")[-1]):
            self.assertEqual(structure[62:65], "...")
        svg, _, _ = linear_sequence_svg(architecture, TRIGGER)
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
        self.assertEqual("".join(bases)[len(TRIGGER):], model.switch_sequence)
        self.assertIn("Linker", svg)

    def test_toehold_derived_from_trigger(self):
        toehold = next(d for d in self.model.domains if d.label == "Toehold")
        self.assertNotIn(PLACEHOLDER_BASE, toehold.sequence)
        self.assertEqual(len(toehold.sequence), self.arch.toehold_length)

    def test_full_left_arm_is_fixed_to_trigger_in_all_recognized_domains(self):
        recognized = [domain for domain in self.model.domains
                      if domain.label in {"Toehold", "Stem L", "Bulge", "U L"}]
        self.assertEqual([len(domain.sequence) for domain in recognized], [12, 9, 3, 6])
        self.assertEqual("".join(domain.sequence for domain in recognized), reverse_complement_rna(TRIGGER))
        self.assertTrue(all(not domain.optimizable for domain in recognized))
        self.assertEqual((self.model.off_structure, self.model.on_structure),
                         target_switch_structures(self.arch, len(self.model.domains[-1].sequence)))

    def test_full_left_arm_rejects_trigger_length_mismatch_instead_of_showing_fake_placeholders(self):
        for trigger in (TRIGGER[:-1], TRIGGER + "A"):
            with self.subTest(trigger_length=len(trigger)), self.assertRaisesRegex(ValueError, "région reconnue"):
                build_switch_model(self.arch, trigger)

    def test_missing_trigger_uses_placeholders(self):
        model = build_switch_model(self.arch, None)
        toehold = next(d for d in model.domains if d.label == "Toehold")
        self.assertEqual(set(toehold.sequence), {PLACEHOLDER_BASE})

    def test_old_n_pattern_is_rejected_and_x_is_only_for_design(self):
        for changes in ({"frame_linker": "CN"}, {"rbs_loop": "AGNG"}, {"t7_leader": "GGN"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                build_switch_model(replace(self.arch, **changes), TRIGGER)
        with self.assertRaises(ValueError):
            build_switch_model(self.arch, "N" + TRIGGER[1:])
        model = build_switch_model(replace(self.arch, frame_linker="CX"), "X" + TRIGGER[1:])
        self.assertIn("X", model.switch_sequence)
        self.assertEqual(model.trigger_sequence, "X" + TRIGGER[1:])

    def test_standard_rbs_loop_has_three_optimized_unpaired_bases_before_user_rbs(self):
        architecture = replace(self.arch, rbs_loop="AGGAGG", t7_leader="GGA")
        model = build_switch_model(architecture, TRIGGER)
        domains = model.domains
        rbs_index = next(i for i, domain in enumerate(domains) if domain.label == "RBS")
        prefix = domains[rbs_index - 1]
        self.assertEqual(prefix.sequence, "XXX")
        self.assertTrue(prefix.optimizable)
        self.assertEqual(domains[rbs_index].sequence, "AGGAGG")
        self.assertFalse(domains[rbs_index].optimizable)
        loop_start = sum(len(domain.sequence) for domain in domains[:rbs_index - 1])
        loop_length = len(prefix.sequence) + len(domains[rbs_index].sequence)
        self.assertEqual(model.off_structure[loop_start:loop_start + loop_length], "." * loop_length)
        self.assertEqual(model.on_structure.split("+")[1][loop_start:loop_start + loop_length], "." * loop_length)
        rbs_end = loop_start + loop_length
        aug_start = sum(len(domain.sequence) for domain in domains[:next(i for i, domain in enumerate(domains) if domain.label == "AUG")])
        self.assertEqual(aug_start - rbs_end, 6)
        self.assertEqual(model.switch_sequence[rbs_end:aug_start], "XXXXXX")

    def test_custom_loop_and_spacing_keep_target_structures_consistent(self):
        architecture = replace(self.arch, rbs_prefix_length=5, stem_length=8,
                               upper_stem_length=8, aug_spacer=2, bulge_length=2)
        model = build_switch_model(architecture, TRIGGER)
        self.assertEqual(parse_structure(model.off_structure)[1], len(model.switch_sequence))
        self.assertEqual(parse_structure(model.on_structure)[1], len(model.on_sequence))
        self.assertEqual(
            (model.off_structure, model.on_structure),
            target_switch_structures(architecture, len(model.domains[-1].sequence)),
        )

    def test_stop_highlights_follow_aug_after_added_rbs_loop_bases(self):
        model = self.model
        aug_start = 0
        for domain in model.domains:
            if domain.label == "AUG":
                break
            aug_start += len(domain.sequence)
        resolved = model.switch_sequence.replace("X", "A")
        stop_start = aug_start + 3
        sequence = resolved[:stop_start] + "UAA" + resolved[stop_start + 3:]
        self.assertEqual(premature_stop_positions(self.arch, sequence), list(range(stop_start, stop_start + 3)))

    def test_linear_view_contains_exactly_the_same_loop_bases_as_2d(self):
        svg, _, _ = linear_sequence_svg(self.arch, TRIGGER)
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
        self.assertEqual("".join(bases), self.model.trigger_sequence[::-1] + self.model.switch_sequence)
        self.assertIn("Boucle RBS", svg)

    def test_standard_24_and_32_keep_legacy_target_structures(self):
        for trigger_length in (24, 32):
            with self.subTest(trigger_length=trigger_length):
                architecture = default_architecture_for_trigger(trigger_length)
                trigger = ("ACGT" * 8)[:trigger_length]
                model = build_switch_model(architecture, trigger)
                self.assertEqual(architecture.binding_layout, "legacy_lower_stem")
                self.assertIsNone(architecture.frame_linker)
                self.assertEqual(architecture.resolved_frame_linker(), "XX" if trigger_length == 24 else "")
                self.assertEqual(next(domain.sequence for domain in model.domains if domain.label == "Bulge"), "XXX")
                self.assertEqual(next(domain.sequence for domain in model.domains if domain.label == "U L"), "X" * 6)
                self.assertEqual(
                    (model.off_structure, model.on_structure),
                    target_switch_structures(architecture, len(model.domains[-1].sequence)),
                )


class LayoutTests(unittest.TestCase):
    def test_layout_returns_all_coordinates(self):
        structure = "..((((....))))..."
        points = radial_layout(structure)
        self.assertEqual(len(points), len(structure))
        for x, y in points:
            self.assertTrue(all(map(lambda v: v == v, (x, y))))  # no NaN

    def test_backbone_spacing_is_uniform(self):
        structure = "((((....))))"
        points = radial_layout(structure, spacing=1.0)
        gaps = [hypot(points[i + 1][0] - points[i][0], points[i + 1][1] - points[i][1]) for i in range(len(points) - 1)]
        for gap in gaps:
            self.assertAlmostEqual(gap, 1.0, delta=0.15)


class SvgTests(unittest.TestCase):
    def test_render_is_well_formed_xml(self):
        model = build_switch_model(default_architecture_for_trigger(len(TRIGGER)), TRIGGER)
        svg = render_structure_svg(model.switch_sequence, model.off_structure, model.switch_colors)
        root = ET.fromstring(svg)
        self.assertTrue(root.tag.endswith("svg"))

    def test_switch_preview_svg_off_and_on(self):
        arch = default_architecture_for_trigger(len(TRIGGER))
        for state in ("off", "on"):
            svg = switch_preview_svg(arch, TRIGGER, state=state)
            ET.fromstring(svg)  # raises on malformed XML
            self.assertIn("<circle", svg)

    def test_off_preview_uses_the_horizontal_axis_and_stays_inside(self):
        arch = default_architecture_for_trigger(len(TRIGGER))
        svg = switch_preview_svg(arch, TRIGGER, state="off", width=1000, height=500)
        circles = [
            (float(x), float(y))
            for x, y in re.findall(r'<circle cx="([\d.]+)" cy="([\d.]+)"', svg)
        ]
        xs = [x for x, _ in circles]
        ys = [y for _, y in circles]
        self.assertGreater(max(xs) - min(xs), 1.5 * (max(ys) - min(ys)))
        self.assertGreaterEqual(min(xs), 30)
        self.assertLessEqual(max(xs), 970)
        self.assertGreaterEqual(min(ys), 80)
        self.assertLessEqual(max(ys), 425)

    def test_on_preview_stays_inside_wide_canvas(self):
        arch = default_architecture_for_trigger(len(TRIGGER))
        svg = switch_preview_svg(arch, TRIGGER, state="on", width=1000, height=500)
        circles = [
            (float(x), float(y))
            for x, y in re.findall(r'<circle cx="([\d.]+)" cy="([\d.]+)"', svg)
        ]
        self.assertTrue(circles)
        self.assertGreaterEqual(min(x for x, _ in circles), 30)
        self.assertLessEqual(max(x for x, _ in circles), 970)
        self.assertGreaterEqual(min(y for _, y in circles), 80)
        self.assertLessEqual(max(y for _, y in circles), 425)


if __name__ == "__main__":
    unittest.main()
