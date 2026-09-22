from __future__ import annotations

import re
import sys
import unittest
import xml.etree.ElementTree as ET
from math import hypot
from pathlib import Path
from dataclasses import asdict, replace
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.preview import (
    PLACEHOLDER_BASE,
    TRIGGER_COLOR,
    build_switch_model,
    _cap_structure,
    _invasion_structures,
    linear_sequence_svg,
    parse_structure,
    premature_stop_positions,
    radial_layout,
    render_structure_svg,
    switch_animation_frames,
    switch_preview_svg,
)

TRIGGER = "CTCTCCTTACGCCACCCACACCCGACGTAC"  # Standard 30 nt trigger region


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
        self.assertEqual(len(pairs), self.arch.toehold_length + self.arch.stem_length)

    def test_optimizable_domains_are_placeholders(self):
        optimizable = [d for d in self.model.domains if d.optimizable]
        self.assertTrue(optimizable)
        for domain in optimizable:
            self.assertEqual(set(domain.sequence), {PLACEHOLDER_BASE})

    def test_frame_linker_is_optimized_and_unpaired_in_both_target_previews(self):
        architecture = replace(self.arch, upper_stem2_length=4)
        model = build_switch_model(architecture, TRIGGER)
        linker_index = next(i for i, domain in enumerate(model.domains) if domain.label == "Linker")
        linker = model.domains[linker_index]
        self.assertEqual(linker.sequence, "NN")
        self.assertTrue(linker.optimizable)
        linker_start = sum(len(domain.sequence) for domain in model.domains[:linker_index])
        for structure in (model.off_structure, model.on_structure.split("+")[-1]):
            self.assertEqual(structure[linker_start:linker_start + 2], "..")

        fixed = build_switch_model(replace(self.arch, frame_linker="CC"), TRIGGER)
        fixed_linker = next(domain for domain in fixed.domains if domain.label == "Linker")
        self.assertEqual(fixed_linker.sequence, "CC")
        self.assertFalse(fixed_linker.optimizable)

    def test_linear_preview_shows_free_linker_then_actual_candidate_bases(self):
        architecture = replace(self.arch, upper_stem2_length=4)
        model = build_switch_model(architecture, TRIGGER)
        linker_index = next(i for i, domain in enumerate(model.domains) if domain.label == "Linker")
        linker_start = sum(len(domain.sequence) for domain in model.domains[:linker_index])
        sequence = model.switch_sequence.replace("N", "A")
        sequence = sequence[:linker_start] + "GU" + sequence[linker_start + 2:]
        for actual, expected_linker in ((None, "NN"), (sequence, "GU")):
            with self.subTest(candidate=actual is not None):
                svg, _, _ = linear_sequence_svg(architecture, TRIGGER, switch_sequence=actual)
                root = ET.fromstring(svg)
                bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
                switch = "".join(bases)[len(model.trigger_sequence):]
                self.assertEqual(switch[linker_start:linker_start + 2], expected_linker)
                self.assertEqual(switch, actual or model.switch_sequence)

    def test_standard_30_has_exactly_18_stem_bases_after_aug_before_reporter(self):
        domains = self.model.domains
        labels = [domain.label for domain in domains]
        aug_index = labels.index("AUG")
        self.assertEqual(labels[aug_index + 1:], ["Stem R", "sfGFP"])
        self.assertNotIn("U2 L", labels)
        self.assertNotIn("U2 R", labels)
        self.assertNotIn("Linker", labels)
        self.assertEqual(next(domain.sequence for domain in domains if domain.label == "RBS"), "AGAGGAGA")
        aug_end = sum(len(domain.sequence) for domain in domains[:aug_index + 1])
        reporter_start = len(self.model.switch_sequence) - len(domains[-1].sequence)
        self.assertEqual(reporter_start - aug_end, 18)
        self.assertEqual(self.model.switch_sequence[aug_end:reporter_start], "N" * 18)
        self.assertEqual(self.model.off_structure[aug_end:reporter_start], ")" * 18)
        self.assertEqual(self.model.on_structure.split("+")[-1][aug_end:reporter_start], "." * 18)
        svg, _, _ = linear_sequence_svg(self.arch, TRIGGER)
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
        self.assertEqual("".join(bases)[len(TRIGGER):], self.model.switch_sequence)
        self.assertNotIn("U2 L", svg)
        self.assertNotIn("U2 R", svg)
        self.assertNotIn("Linker", svg)

    def test_toehold_derived_from_trigger(self):
        toehold = next(d for d in self.model.domains if d.label == "Toehold")
        self.assertNotIn(PLACEHOLDER_BASE, toehold.sequence)
        self.assertEqual(len(toehold.sequence), self.arch.toehold_length)

    def test_missing_trigger_uses_placeholders(self):
        model = build_switch_model(self.arch, None)
        toehold = next(d for d in model.domains if d.label == "Toehold")
        self.assertEqual(set(toehold.sequence), {PLACEHOLDER_BASE})

    def test_standard_rbs_loop_has_three_optimized_unpaired_bases_before_user_rbs(self):
        architecture = replace(self.arch, rbs_loop="AGGAGG", t7_leader="GGA")
        model = build_switch_model(architecture, TRIGGER)
        domains = model.domains
        rbs_index = next(i for i, domain in enumerate(domains) if domain.label == "RBS")
        prefix = domains[rbs_index - 1]
        self.assertEqual(prefix.sequence, "NNN")
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
        self.assertEqual(model.switch_sequence[rbs_end:aug_start], "NNNNNN")

    def test_custom_loop_and_spacing_stay_consistent_through_binding_animation(self):
        architecture = replace(self.arch, rbs_prefix_length=5, upper_stem_length=8, aug_spacer=2, bulge_length=2)
        model = build_switch_model(architecture, TRIGGER)
        self.assertEqual(parse_structure(model.off_structure)[1], len(model.switch_sequence))
        steps = _invasion_structures(architecture)
        for trigger, switch in steps:
            tail = "." * (len(model.switch_sequence) - len(switch))
            structure = trigger + "+" + switch + tail
            self.assertEqual(parse_structure(structure)[1], len(model.on_sequence))
        trigger, switch = steps[-1]
        self.assertEqual(trigger + "+" + switch + "." * (len(model.switch_sequence) - len(switch)), model.on_structure)

    def test_stop_highlights_follow_aug_after_added_rbs_loop_bases(self):
        model = self.model
        aug_start = 0
        for domain in model.domains:
            if domain.label == "AUG":
                break
            aug_start += len(domain.sequence)
        resolved = model.switch_sequence.replace("N", "A")
        stop_start = aug_start + 3
        sequence = resolved[:stop_start] + "UAA" + resolved[stop_start + 3:]
        self.assertEqual(premature_stop_positions(self.arch, sequence), list(range(stop_start, stop_start + 3)))

    def test_linear_view_contains_exactly_the_same_loop_bases_as_2d(self):
        svg, _, _ = linear_sequence_svg(self.arch, TRIGGER)
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall('.//{http://www.w3.org/2000/svg}text') if node.attrib.get("font-size") == "13"]
        self.assertEqual("".join(bases), self.model.trigger_sequence[::-1] + self.model.switch_sequence)
        self.assertIn("Boucle RBS", svg)

    def test_pre_reload_architecture_keeps_its_original_sequence_and_cap(self):
        historical = replace(self.arch, upper_stem_length=5, upper_stem2_length=4,
                             rbs_loop="UAGAGGAGAAC", rbs_prefix_length=0, frame_linker="CC")
        values = asdict(historical)
        values.pop("rbs_prefix_length")
        legacy = SimpleNamespace(**values, resolved_frame_linker=historical.resolved_frame_linker)
        expected = build_switch_model(historical, TRIGGER)
        self.assertEqual(build_switch_model(legacy, TRIGGER), expected)
        self.assertEqual(_cap_structure(legacy), _cap_structure(historical))
        self.assertEqual(
            linear_sequence_svg(legacy, TRIGGER, switch_sequence=expected.switch_sequence.replace("N", "A")),
            linear_sequence_svg(historical, TRIGGER, switch_sequence=expected.switch_sequence.replace("N", "A")),
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


class AnimationTests(unittest.TestCase):
    def test_frames_are_well_formed_and_ordered(self):
        arch = default_architecture_for_trigger(len(TRIGGER))
        frames = switch_animation_frames(arch, TRIGGER)
        self.assertGreater(len(frames), 20)
        for svg in (frames[0], frames[len(frames) // 2], frames[-1]):
            root = ET.fromstring(svg)
            self.assertTrue(root.tag.endswith("svg"))
            self.assertIn("<circle", svg)

    def test_trigger_stays_linear_during_approach_and_landing(self):
        # The still-unbound trigger must be a straight line (never coiled by
        # the generic radial loop placement) while it approaches and right
        # after it lands on the toehold.
        arch = default_architecture_for_trigger(len(TRIGGER))
        frames = switch_animation_frames(arch, TRIGGER, approach_steps=6, hold_start=2, tween_steps=1)
        pattern = re.compile(
            r'<circle cx="([\d.]+)" cy="([\d.]+)" r="[\d.]+" fill="' + re.escape(TRIGGER_COLOR) + r'"'
        )
        for svg in frames[:8]:  # approach frames + landed beat
            pts = [(float(x), float(y)) for x, y in pattern.findall(svg)]
            self.assertGreater(len(pts), 2)
            # Use the two most distant points (first/last) as the reference
            # line: with adjacent points only ~1-3px apart, per-coordinate
            # SVG rounding (.1f) would otherwise dominate the measurement.
            (x0, y0), (x1, y1) = pts[0], pts[-1]
            dx, dy = x1 - x0, y1 - y0
            norm = hypot(dx, dy) or 1e-6
            for x, y in pts[1:-1]:
                distance_from_line = abs((x - x0) * dy - (y - y0) * dx) / norm
                self.assertLess(distance_from_line, 1.0)

    def test_last_frame_matches_static_on_view(self):
        arch = default_architecture_for_trigger(len(TRIGGER))
        frames = switch_animation_frames(arch, TRIGGER, hold_end=1)
        on_svg = switch_preview_svg(arch, TRIGGER, state="on")
        # Same number of bases drawn in both the final animation frame and the
        # static ON view (same underlying sequence/structure).
        self.assertEqual(frames[-1].count("<circle"), on_svg.count("<circle"))

    def test_toehold_contact_point_never_moves(self):
        # The point where the trigger's tip first touches the switch (switch
        # index t7_len, i.e. the toehold's outer edge) must stay pinned at the
        # exact same pixel in every single frame: only the rest of the
        # molecule may reconfigure as the stem opens around it, never the
        # contact point itself (otherwise it reads as "the switch sliding").
        arch = default_architecture_for_trigger(len(TRIGGER))
        model = build_switch_model(arch, TRIGGER)
        anchor_index = len(model.trigger_sequence) + len(arch.t7_leader)
        frames = switch_animation_frames(
            arch, TRIGGER, approach_steps=3, hold_start=2, tween_steps=2, hold_end=2
        )
        circle_pattern = re.compile(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="[\d.]+"')
        ref_x = ref_y = None
        for svg in frames:
            pts = circle_pattern.findall(svg)
            x, y = float(pts[anchor_index][0]), float(pts[anchor_index][1])
            if ref_x is None:
                ref_x, ref_y = x, y
            self.assertAlmostEqual(x, ref_x, delta=0.2)
            self.assertAlmostEqual(y, ref_y, delta=0.2)

    def test_frame_count_scales_with_stem_length(self):
        short_arch = default_architecture_for_trigger(len(TRIGGER))
        frames = switch_animation_frames(
            short_arch, TRIGGER, approach_steps=3, tween_steps=2, hold_start=1, hold_end=1
        )
        stem_len = max(1, short_arch.stem_length)
        expected = 3 + 1 + stem_len * 2 + 1
        self.assertEqual(len(frames), expected)


if __name__ == "__main__":
    unittest.main()
