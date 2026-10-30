"""Full reporter drawings retain measured bases/pairs and architecture colors."""
from copy import deepcopy
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.localization import language_context
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.preview import APTAMER_COLOR, DOMAIN_COLORS, TRIGGER_COLOR, render_linear_svg, render_structure_svg
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context
from aptaswitch_core.saved_switch_diagrams import context_complex_diagram
from aptaswitch_studio.saved_switch_library import candidate_snapshot
from design_fixtures import design_result_fixture

NS = {"s": "http://www.w3.org/2000/svg"}


class SavedSwitchDiagramTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context("fr"))
        trigger = "ACGTACGTACGTACGTACGTACGTACGTAC"
        self.request = RunRequest(sequences=SequenceInput("AACCGGTTAACC", trigger), architecture=default_architecture_for_trigger(len(trigger)), trials=1)

    def data(self):
        result = design_result_fixture(self.request)
        snapshot = candidate_snapshot(result, result.candidates[0], "Saved")
        context = build_saved_switch_context(snapshot, {})
        return snapshot, context

    def mfe(self, context, state):
        sequences = ([context["trigger_sequence"]] if state == "on" else []) + [context["switch_sequence"]]
        # Deliberately different from the target, including a hairpin in the full reporter.
        structures = ["." * len(seq) for seq in sequences]
        structures[-1] = structures[-1][:-12] + "(((......)))"
        return {"structure": "+".join(structures), "sequences": sequences,
                "strand_order": ["Trigger", "Switch"] if state == "on" else ["Switch"]}

    def test_all_domains_cover_full_reporter_in_off_and_on_without_changing_mfe(self):
        snapshot, context = self.data()
        layout = architecture_layout(self.request.architecture)
        for state in ("off", "on"):
            mfe = self.mfe(context, state)
            original = deepcopy(mfe)
            diagram = context_complex_diagram(snapshot, context, mfe)
            offset = len(context["trigger_sequence"]) if state == "on" else 0
            self.assertEqual(diagram.sequence, "".join(mfe["sequences"]))
            self.assertEqual(diagram.structure, mfe["structure"])
            self.assertEqual(diagram.colors[offset + layout["reporter_start"]:],
                             [DOMAIN_COLORS["Reporter"]] * context["reporter_length"])
            self.assertEqual(diagram.colors[offset + layout["aug_start"]:offset + layout["aug_start"] + 3], [DOMAIN_COLORS["AUG"]] * 3)
            self.assertEqual(diagram.colors[:offset], [TRIGGER_COLOR] * offset)
            self.assertIn(("Linker", DOMAIN_COLORS["Linker"]), diagram.legend)
            self.assertIn(("RBS", DOMAIN_COLORS["RBS"]), diagram.legend)
            self.assertEqual(any(label == "Trigger" for label, _ in diagram.legend), state == "on")
            with patch("aptaswitch_core.saved_switch_diagrams.render_structure_svg", wraps=render_structure_svg) as renderer:
                diagram.structure_svg(show_annotations=False)
            self.assertEqual(renderer.call_args.args[:2], (diagram.sequence, mfe["structure"]))
            self.assertEqual(mfe, original)

    def test_linear_view_uses_exact_bases_and_measured_hairpins_without_target_pairing(self):
        snapshot, context = self.data()
        for state in ("off", "on"):
            with self.subTest(state=state):
                diagram = context_complex_diagram(snapshot, context, self.mfe(context, state))
                svg, width, _ = diagram.linear_svg(show_annotations=False)
                root = ET.fromstring(svg)
                groups = root.findall(".//s:g", NS)
                self.assertEqual("".join(g.find("s:text", NS).text for g in groups), diagram.sequence)
                self.assertEqual(len(root.findall('.//s:path[@stroke="#94a3b8"]', NS)), 3)
                self.assertFalse(root.findall('.//s:line[@stroke-width="1"]', NS))
                self.assertEqual([row.label for row in diagram.rows], ["Trigger", "Switch"] if state == "on" else ["Switch"])
                self.assertGreater(width, 15000)
                self.assertNotIn("N", diagram.sequence)
                self.assertNotIn("T", context["switch_sequence"])
                self.assertNotIn("U", context["trigger_sequence"])

    def test_diagrams_reject_mismatched_strand_chemistries(self):
        snapshot, context = self.data()
        for role, material in (("switch", "dna"), ("trigger", "rna")):
            with self.subTest(role=role), self.assertRaises(ValueError):
                damaged = {**context, f"{role}_material": material}
                context_complex_diagram(snapshot, damaged, self.mfe(context, "on"))

    def test_pentamer_linear_pairs_join_exact_copies_and_nucleotides(self):
        snapshot, context = self.data()
        names = ["Switch", "Trigger", "Switch", "Aptamer", "Trigger"]
        sequences = [context[f"{name.lower()}_sequence"] for name in names]
        strands = [list("." * len(sequence)) for sequence in sequences]
        expected_pairs = []
        for index, strand in enumerate(strands):
            strand[3], strand[8] = "(", ")"
            expected_pairs.append((index, 3, index, 8))
        for index in range(4):
            strands[index][-1], strands[index + 1][0] = "(", ")"
            expected_pairs.append((index, len(strands[index]) - 1, index + 1, 0))
        # A connection across non-adjacent strands must not be attached to an
        # intermediate row or to the other copy of the same species.
        strands[0][1], strands[-1][-2] = "(", ")"
        expected_pairs.append((0, 1, 4, len(strands[4]) - 2))
        mfe = {"strand_order": names, "sequences": sequences, "structure": "+".join("".join(s) for s in strands)}
        diagram = context_complex_diagram(snapshot, context, mfe)
        self.assertEqual([row.label for row in diagram.rows], ["Switch 1", "Trigger 1", "Switch 2", "Aptamer", "Trigger 2"])
        for show in (False, True):
            with patch("aptaswitch_core.saved_switch_diagrams.render_linear_svg", wraps=render_linear_svg) as render:
                svg, width, height = diagram.linear_svg(show_annotations=show)
            self.assertCountEqual(render.call_args.kwargs["pairs"], expected_pairs)
            root = ET.fromstring(svg)
            boxes = [group.find("s:rect", NS) for group in root.findall(".//s:g", NS)]
            row_boxes, start = [], 0
            for seq in sequences:
                row_boxes.append(boxes[start:start + len(seq)])
                start += len(seq)
            lines = root.findall('.//s:line[@stroke="#2563eb"]', NS)
            self.assertEqual(len(lines), 5)
            self.assertEqual(len(root.findall('.//s:path[@stroke="#94a3b8"]', NS)), 5)
            endpoints = {(float(line.attrib["x1"]), float(line.attrib["y1"]),
                          float(line.attrib["x2"]), float(line.attrib["y2"])) for line in lines}
            for ra, ia, rb, ib in expected_pairs:
                if ra == rb:
                    continue
                a, b = row_boxes[ra][ia].attrib, row_boxes[rb][ib].attrib
                self.assertIn((float(a["x"]) + 10, float(a["y"]) + 26, float(b["x"]) + 10, float(b["y"])), endpoints)
            self.assertTrue(all(float(box.attrib["y"]) + 26 < height for box in boxes))
            self.assertGreater(width, 15000)
            numbering = root.findall('.//s:text[@font-size="9"]', NS)
            self.assertEqual(sum(n.text == "10" for n in numbering), 5)

    def test_unpaired_complex_has_no_invented_hybridization_lines(self):
        snapshot, context = self.data()
        names = ["Trigger", "Trigger", "Trigger"]
        diagram = context_complex_diagram(snapshot, context, self.complex_mfe(context, names))
        root = ET.fromstring(diagram.linear_svg()[0])
        self.assertFalse(root.findall('.//s:line[@stroke="#2563eb"]', NS))
        self.assertFalse(root.findall('.//s:path[@stroke="#94a3b8"]', NS))

    def test_inconsistent_record_is_rejected_before_placeholder_rendering(self):
        snapshot, context = self.data()
        for change in ("length", "sequence", "order", "missing"):
            mfe = self.mfe(context, "on")
            if change == "length":
                mfe["structure"] = "." + mfe["structure"]
            elif change == "sequence":
                mfe["sequences"][0] = "T" * len(mfe["sequences"][0])
            elif change == "order":
                mfe["strand_order"].reverse()
            else:
                mfe["structure"] = ""
            with self.assertRaises(ValueError):
                context_complex_diagram(snapshot, context, mfe)

    def test_live_and_export_legends_include_all_domains_and_wrap(self):
        snapshot, context = self.data()
        snapshot["request"]["reporter"]["label"] = "Rapporteur complet personnalisé pour le candidat"
        diagram = context_complex_diagram(snapshot, context, self.mfe(context, "on"))
        for show in (False, True):
            root = ET.fromstring(diagram.structure_svg(show_annotations=show))
            legend_text = root.findall('.//s:text[@fill="#334155"]', NS)
            self.assertEqual([node.text for node in legend_text], [label for label, _ in diagram.legend])
            self.assertGreater(len({node.attrib["y"] for node in legend_text}), 1)
            self.assertTrue(all(float(node.attrib["y"]) < 650 for node in legend_text))

    def complex_mfe(self, context, names):
        sequences = [context[f"{name.lower()}_sequence"] for name in names]
        structures = ["." * len(sequence) for sequence in sequences]
        if len(names) == 2:
            # Preserve an intermolecular duplex, including between identical copies.
            structures[0] = "(((" + structures[0][3:]
            structures[1] = structures[1][:-3] + ")))"
        return {"strand_order": names, "sequences": sequences, "structure": "+".join(structures)}

    def test_all_monomers_and_dimers_preserve_exact_measured_order_and_pairs(self):
        snapshot, context = self.data()
        names = ("Switch", "Trigger", "Aptamer")
        compositions = [[name] for name in names] + [[left, right] for left in names for right in names]
        for composition in compositions:
            with self.subTest(composition=composition):
                mfe = self.complex_mfe(context, composition)
                original = deepcopy(mfe)
                diagram = context_complex_diagram(snapshot, context, mfe)
                self.assertEqual(diagram.structure, mfe["structure"])
                self.assertEqual(diagram.sequence, "".join(mfe["sequences"]))
                self.assertEqual(len(diagram.colors), len(diagram.sequence))
                labels = ([f"{composition[0]} 1", f"{composition[0]} 2"]
                          if len(composition) == 2 and composition[0] == composition[1] else composition)
                self.assertEqual([row.label for row in diagram.rows], labels)
                for row, sequence in zip(diagram.rows, mfe["sequences"]):
                    self.assertEqual("".join(base.base for base in row.bases), sequence)
                    self.assertEqual([base.position for base in row.bases], list(range(1, len(sequence) + 1)))
                self.assertEqual(mfe, original)

    def test_every_switch_copy_has_full_reporter_domains_and_stop_annotations(self):
        snapshot, context = self.data()
        layout = architecture_layout(self.request.architecture)
        start = layout["aug_start"] + 3
        context["switch_sequence"] = context["switch_sequence"][:start] + "UAG" + context["switch_sequence"][start + 3:]
        single = context_complex_diagram(snapshot, context, self.complex_mfe(context, ["Switch"]))
        self.assertTrue(set(range(start, start + 3)).issubset(single.highlight))
        for composition in (["Switch", "Switch"], ["Trigger", "Switch"], ["Aptamer", "Switch"]):
            diagram = context_complex_diagram(snapshot, context, self.complex_mfe(context, composition))
            expected_highlight, offset = set(), 0
            for row, name in zip(diagram.rows, composition):
                if name == "Switch":
                    self.assertEqual(row.domains, single.rows[0].domains)
                    self.assertEqual([base.color for base in row.bases], single.colors)
                    self.assertEqual([base.color for base in row.bases[layout["reporter_start"]:]],
                                     [DOMAIN_COLORS["Reporter"]] * context["reporter_length"])
                    expected_highlight.update(offset + index for index in single.highlight)
                offset += len(row.bases)
            self.assertEqual(diagram.highlight, expected_highlight)
            self.assertEqual(sum(label == "RBS" for label, _ in diagram.legend), 1)
            self.assertEqual(sum(label == "Linker" for label, _ in diagram.legend), 1)

    def test_aptamer_and_trigger_colors_do_not_invent_switch_or_ligand_annotations(self):
        snapshot, context = self.data()
        diagram = context_complex_diagram(snapshot, context, self.complex_mfe(context, ["Trigger", "Aptamer"]))
        self.assertEqual(diagram.colors, [TRIGGER_COLOR] * len(context["trigger_sequence"])
                         + [APTAMER_COLOR] * len(context["aptamer_sequence"]))
        self.assertEqual(diagram.legend, [("Trigger", TRIGGER_COLOR), ("Aptamer", APTAMER_COLOR)])
        self.assertFalse(diagram.highlight)
        with patch("aptaswitch_core.saved_switch_diagrams.render_structure_svg", wraps=render_structure_svg) as renderer:
            diagram.structure_svg(show_annotations=True)
        self.assertEqual(renderer.call_args.args[:2], (diagram.sequence, diagram.structure))
        self.assertEqual(renderer.call_args.kwargs["legend"], diagram.legend)

    def test_unexpected_complex_rejects_invalid_structure_sequences_and_chemistry(self):
        snapshot, context = self.data()
        for change in ("order", "sequence", "unbalanced", "empty_strand", "short", "chemistry", "size", "unknown", "invalid_name"):
            with self.subTest(change=change):
                current = deepcopy(context)
                mfe = self.complex_mfe(current, ["Aptamer", "Switch"])
                if change == "order":
                    mfe["strand_order"].reverse()
                elif change == "sequence":
                    mfe["sequences"][0] = "T" * len(mfe["sequences"][0])
                elif change == "unbalanced":
                    mfe["structure"] = mfe["structure"][:-1] + "."
                elif change == "empty_strand":
                    mfe["structure"] += "+"
                elif change == "short":
                    mfe["structure"] = mfe["structure"][1:]
                elif change == "chemistry":
                    current["switch_sequence"] = current["switch_sequence"].replace("U", "T")
                    mfe = self.complex_mfe(current, ["Aptamer", "Switch"])
                elif change == "size":
                    mfe = self.complex_mfe(current, ["Switch"] * 6)
                elif change == "unknown":
                    mfe["strand_order"][0] = "Unknown"
                else:
                    mfe["strand_order"][0] = []
                with self.assertRaises(ValueError):
                    context_complex_diagram(snapshot, current, mfe)

if __name__ == "__main__":
    unittest.main()
