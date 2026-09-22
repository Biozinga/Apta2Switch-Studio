"""Initial-system analysis must respect inputs without launching a search."""

from __future__ import annotations

import json
import os
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

from aptaswitch_core import system_analysis
from aptaswitch_core.models import ThermoConditions


class SystemAnalysisTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def analyze(self, *, aptamer="ac gu\n", trigger="tg ca", thermo=None, structure="((..+..))", energy=-3.25):
        fake = SimpleNamespace(Model=Mock(return_value=object()), mfe=Mock(return_value=[SimpleNamespace(structure=structure, energy=energy)]))
        with patch.object(system_analysis, "nupack_import_context") as context:
            context.return_value.__enter__.return_value = fake
            result = system_analysis.analyze_system("/test/nupack", aptamer, trigger, thermo or ThermoConditions(material="dna04"))
        return result, fake

    def test_dna_normalization_exact_thermo_and_single_mfe_call(self):
        conditions = ThermoConditions(material="dna04", temperature_c=26.75, sodium_m=.042, magnesium_m=.007)
        result, fake = self.analyze(thermo=conditions)
        fake.Model.assert_called_once_with(material="dna04", celsius=26.75, sodium=.042, magnesium=.007, ensemble="stacking")
        fake.mfe.assert_called_once_with(strands=["dACGT", "dTGCA"], model=fake.Model.return_value)
        self.assertEqual(result["aptamer_sequence"], "ACGT")
        self.assertEqual(result["trigger_sequence"], "TGCA")
        self.assertEqual(result["paired_trigger_indices"], [2, 3])
        self.assertEqual(result["mfe_structure"], "((..+..))")
        self.assertEqual(result["mfe_energy_kcal_mol"], -3.25)
        self.assertEqual(result["model"]["magnesium_m"], .007)
        self.assertEqual(json.loads(json.dumps(result)), result)

    def test_rna_normalizes_both_strands_and_uses_requested_model(self):
        result, fake = self.analyze(thermo=ThermoConditions(material="rna06", temperature_c=30, sodium_m=.5, magnesium_m=0))
        fake.Model.assert_called_once_with(material="rna06", celsius=30, sodium=.5, magnesium=0, ensemble="stacking")
        fake.mfe.assert_called_once_with(strands=["rACGU", "rUGCA"], model=fake.Model.return_value)
        self.assertEqual(result["model"]["aptamer_material"], "rna")
        self.assertEqual(result["model"]["trigger_material"], "rna")
        self.assertEqual(result["trigger_sequence"], "UGCA")

    def test_trigger_self_pairs_are_not_marked_as_aptamer_binding(self):
        result, _ = self.analyze(aptamer="ACGTACGT", trigger="ACGTACGT", structure="((....((+))....))")
        self.assertEqual(result["paired_trigger_indices"], [0, 1, 6, 7])
        result, _ = self.analyze(aptamer="ACGTACGT", trigger="ACGTACGT", structure="((....))+((....))")
        self.assertEqual(result["paired_trigger_indices"], [])

    def test_invalid_input_never_calls_nupack(self):
        requests = [
            ("", "ACGT", "TGCA", ThermoConditions(material="dna04")),
            ("/test", "ACGN", "TGCA", ThermoConditions(material="dna04")),
            ("/test", "ACGT", "", ThermoConditions(material="dna04")),
            ("/test", "ACGT", "TGCA", ThermoConditions(material="rna-dna06")),
            ("/test", "ACGT", "TGCA", ThermoConditions(material="dna04", sodium_m=-.01)),
            ("/test", "ACGT", "TGCA", ThermoConditions(material="dna04", magnesium_m=-.01)),
            ("/test", "ACGT", "TGCA", ThermoConditions(material="dna04", temperature_c=float("nan"))),
            ("/test", "ACGT", "TGCA", ThermoConditions(material="dna04", sodium_m=float("inf"))),
        ]
        with patch.object(system_analysis, "nupack_import_context") as context:
            for request in requests:
                with self.subTest(request=request), self.assertRaises((ValueError, RuntimeError)):
                    system_analysis.analyze_system(*request)
            context.assert_not_called()

    def test_missing_or_invalid_nupack_output_has_no_schematic_fallback(self):
        fake = SimpleNamespace(Model=Mock(), mfe=Mock(return_value=[]))
        with patch.object(system_analysis, "nupack_import_context") as context:
            context.return_value.__enter__.return_value = fake
            with self.assertRaisesRegex(RuntimeError, "aucune structure"):
                system_analysis.analyze_system("/test", "ACGT", "TGCA", ThermoConditions(material="dna04"))
        for structure, energy in [("....+...", -1), ("........", -1), ("((..+....", -1), ("....+....", float("nan"))]:
            with self.subTest(structure=structure, energy=energy), self.assertRaises((ValueError, RuntimeError)):
                self.analyze(structure=structure, energy=energy)


@unittest.skipUnless(os.environ.get("APTASWITCH_TEST_NUPACK"), "Set APTASWITCH_TEST_NUPACK for real NUPACK checks")
class RealNupackSystemAnalysisTests(unittest.TestCase):
    def test_real_dna_and_rna_structures(self):
        for material in ("dna04", "rna06"):
            with self.subTest(material=material):
                conditions = ThermoConditions(material=material, temperature_c=29.5, sodium_m=.15, magnesium_m=0)
                result = system_analysis.analyze_system(os.environ["APTASWITCH_TEST_NUPACK"], "ACGTACGT", "ACGT", conditions)
                self.assertEqual([len(part) for part in result["mfe_structure"].split("+")], [8, 4])
                self.assertTrue(all(0 <= index < 4 for index in result["paired_trigger_indices"]))
                self.assertEqual(result["model"]["material"], material)
                self.assertEqual(result["model"]["temperature_c"], 29.5)
                self.assertEqual(result["model"]["sodium_m"], .15)
                self.assertNotIn("T" if material == "rna06" else "U", result["aptamer_sequence"])


if __name__ == "__main__":
    unittest.main()
