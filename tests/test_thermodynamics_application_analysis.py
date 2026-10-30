"""Application contracts exercised without a selected production engine."""
from dataclasses import asdict, replace
import json
import math
import os
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.thermodynamics import application_analysis as analysis
from aptaswitch_core.thermodynamics.structures import PairingConstraints
from aptaswitch_core.thermodynamics.tube_reporting import TubeAnalysisConditions
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention, solve_tube
from aptaswitch_core.thermodynamics.types import (
    ComplexSpec, Conditions, EmptyEnsembleError, ModelProvenance, PartitionResult,
    Strand, enumerate_complexes,
)


PROVENANCE = ModelProvenance(
    "synthetic", "1", "application-fixture", "sha256:" + "c" * 64,
    "connected-fixture", "audited-1M-fixture",
)
STANDARD = StandardStateConvention(PROVENANCE.standard_state, 1., True, True, True)


class FixtureBackend:
    """Analytic association constants, not a thermodynamic parameter model."""
    provenance = PROVENANCE

    def __init__(self, association=None):
        self.association = association or {}
        self.calls = []
        self.minimum_calls = []

    def partition(self, spec, conditions):
        self.calls.append((spec, conditions))
        order = "".join(s.name for s in spec.strands)
        key = min(order[i:] + order[:i] for i in range(len(order)))
        constant = 1. if len(spec.strands) == 1 else self.association.get(key, 0.)
        log_z = math.log(constant) if constant else -math.inf
        return PartitionResult(spec, conditions, self.provenance, log_z, -.6 * log_z, .6)

    def minimum(self, spec, conditions):
        self.minimum_calls.append((spec, conditions))
        context = SimpleNamespace(complex_spec=spec, conditions=conditions,
                                  provenance=self.provenance, constraints=None,
                                  thermal_energy_kcal_mol=.6)
        return SimpleNamespace(context=context, structure="(...+...)", energy_kcal_mol=-2.4)


class FixedSystemTests(unittest.TestCase):
    def setUp(self):
        self.backend = FixtureBackend()
        self.thermo = ThermoConditions(temperature_c=31, sodium_m=.15, magnesium_m=0, switch_material="dna")

    def system(self, **kwargs):
        return analysis.analyze_system(self.backend, "ac gu\n", "a cgu", self.thermo, **kwargs)

    def test_original_material_order_conditions_and_json_shape(self):
        result = self.system()
        spec, conditions = self.backend.minimum_calls[0]
        self.assertEqual([s.name for s in spec.strands], ["A", "T"])
        self.assertEqual([s.material for s in spec.strands], ["dna", "dna"])
        self.assertEqual([s.sequence for s in spec.strands], ["ACGT", "ACGT"])
        self.assertEqual(conditions, Conditions(31, .15, 0))
        self.assertEqual(result["mfe_structure"], "(...+...)")
        self.assertEqual(result["mfe_energy_kcal_mol"], -2.4)
        self.assertEqual(result["paired_trigger_indices"], [3])
        self.assertEqual(result["strand_order"], ["aptamer", "trigger"])
        self.assertEqual(result["model"]["material"], PROVENANCE.model)
        self.assertEqual(result["thermodynamics"], asdict(PROVENANCE))

        self.assertEqual(json.loads(json.dumps(result, allow_nan=False)), result)
        self.assertFalse(result["mfe_enumeration"]["exhaustive"])
        self.assertIsNone(result["mfe_enumeration"]["truncated"])

    def test_other_system_chemistries_are_rejected_before_folding(self):
        for aptamer, trigger in (("rna", "rna"), ("rna", "dna"), ("dna", "rna")):
            self.thermo = replace(self.thermo, switch_material=aptamer, trigger_material=trigger)
            with self.subTest(aptamer=aptamer, trigger=trigger), self.assertRaises(ValueError):
                self.system()
        self.assertEqual(self.backend.minimum_calls, [])

    def test_single_minimum_cannot_claim_complete_alternatives(self):
        with self.assertRaisesRegex(ValueError, "cannot enumerate"):
            self.system(mfe_policy="complete")
        self.assertEqual(self.backend.minimum_calls, [])

    def test_bounded_alternatives_are_preserved_and_truncation_is_explicit(self):
        def minima(spec, conditions, *, maximum_count, absolute_tolerance):
            self.assertEqual(maximum_count, 2)
            record = self.backend.minimum(spec, conditions)
            return SimpleNamespace(context=record.context, minimum_energy_kcal_mol=record.energy_kcal_mol,
                                   energies_kcal_mol=(record.energy_kcal_mol, record.energy_kcal_mol + 1e-9),
                                   structures=("(...+...)", ".(..+..)."), truncated=False,
                                   maximum_count=maximum_count, tolerance_kcal_mol=absolute_tolerance)
        self.backend.minimum_structures = minima
        result = self.system(mfe_policy="complete", max_mfe_structures=2)
        self.assertEqual(len(result["mfe_structures"]), 2)
        self.assertEqual(result["mfe_structures"][1]["energy_kcal_mol"], -2.4 + 1e-9)
        self.assertTrue(result["mfe_enumeration"]["exhaustive"])
        self.backend.minimum_structures = lambda *args, **kw: SimpleNamespace(
            **{**minima(*args, **kw).__dict__, "truncated": True})
        result = self.system(max_mfe_structures=2)
        self.assertTrue(result["mfe_enumeration"]["truncated"])
        with self.assertRaisesRegex(ValueError, "completeness"):
            self.system(mfe_policy="complete", max_mfe_structures=2)

    def test_intramolecular_pairs_are_excluded_from_trigger_contacts(self):
        original = self.backend.minimum
        def minimum(spec, conditions):
            record = original(spec, conditions)
            record.structure = "(.....)(+..)."
            return record
        self.backend.minimum = minimum
        result = analysis.analyze_system(self.backend, "ACGTACGT", "ACGT", self.thermo)
        self.assertEqual(result["paired_trigger_indices"], [2])

    def test_invalid_or_disconnected_results_are_not_fabricated_into_a_guide(self):
        original = self.backend.minimum
        for structure, energy in (("....+....", 0), ("((..+...)", -2),
                                  ("....+...", -2), ("(...+...)", math.nan)):
            def minimum(spec, conditions):
                record = original(spec, conditions)
                record.structure, record.energy_kcal_mol = structure, energy
                return record
            self.backend.minimum = minimum
            with self.subTest(structure=structure, energy=energy), self.assertRaises(ValueError):
                self.system()
        self.backend.minimum = Mock(side_effect=EmptyEnsembleError("empty connected ensemble"))
        with self.assertRaises(EmptyEnsembleError):
            self.system()

    def test_changed_model_conditions_coordinates_or_constraints_are_rejected(self):
        original = self.backend.minimum
        changes = {
            "conditions": Conditions(37, .15),
            "provenance": replace(PROVENANCE, version="2"),
            "constraints": PairingConstraints.unconstrained(8),
            "complex_spec": ComplexSpec((Strand("ACGT", "dna", "T"), Strand("ACGT", "dna", "A"))),
        }
        for field, value in changes.items():
            def minimum(spec, conditions):
                record = original(spec, conditions)
                setattr(record.context, field, value)
                return record
            self.backend.minimum = minimum
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "Backend changed"):
                self.system()

    def test_bad_inputs_do_not_call_backend(self):
        for thermo in (replace(self.thermo, temperature_c=math.nan),
                       replace(self.thermo, sodium_m=-1), replace(self.thermo, trigger_material="auto")):
            with self.assertRaises(ValueError):
                analysis.analyze_system(self.backend, "ACGT", "ACGT", thermo)
        for options in ({"mfe_policy": "guess"}, {"max_mfe_structures": 0}, {"max_mfe_structures": True}):
            with self.assertRaises(ValueError):
                self.system(**options)
        self.assertEqual(self.backend.minimum_calls, [])


class FixedTubeTests(unittest.TestCase):
    def setUp(self):
        self.backend = FixtureBackend()
        self.thermo = ThermoConditions(temperature_c=31, sodium_m=.15, switch_material="dna")
        self.conditions = TubeAnalysisConditions()

    def tubes(self, *, sequences=("ac gu", "tg ca"), **kwargs):
        return analysis.analyze_parasite_tubes(self.backend, *sequences, self.thermo,
                                              self.conditions, standard_state=STANDARD, **kwargs)

    def test_complete_three_strand_catalog_shared_across_tubes_and_original_shape(self):
        progress = Mock()
        result = self.tubes(progress_callback=progress)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(len(self.backend.calls), 9)  # binary necklaces of lengths 1..3
        self.assertEqual(len({spec.cyclic_key for spec, _ in self.backend.calls}), 9)
        self.assertTrue(all(len(spec.strands) <= 3 for spec, _ in self.backend.calls))
        self.assertTrue(all(c == Conditions(31, .15, 0) for _, c in self.backend.calls))
        single, mixed = (result["tubes"][name] for name in ("trigger_only", "aptamer_trigger"))
        self.assertEqual(len(single["complexes"]), 3)
        self.assertEqual(len(mixed["complexes"]), 9)
        self.assertEqual(single["equilibrium"]["ordered_complex_count"], 3)
        self.assertEqual(mixed["equilibrium"]["ordered_complex_count"], 9)
        self.assertAlmostEqual(mixed["summary"]["aptamer_free_pct"], 100)
        self.assertAlmostEqual(mixed["summary"]["trigger_free_pct"], 100)
        self.assertEqual(result["conditions"], asdict(self.conditions))
        self.assertEqual(progress.call_args_list[0].args[:2], (0, 11))
        self.assertEqual(progress.call_args_list[-1].args[:2], (11, 11))
        self.assertEqual(json.loads(json.dumps(result, allow_nan=False)), result)

    def test_analytic_heterodimer_and_homodimer_populations_and_mass_balance(self):
        k, a_total, t_total = 8e6, 2e-6, 7e-6
        self.backend.association = {"AT": k}
        self.conditions = replace(self.conditions, aptamer_concentration_m=a_total,
                                  trigger_concentration_m=t_total, max_complex_size=2)
        result = self.tubes()
        scale = a_total + t_total + 1 / k
        bound = 2 * a_total * t_total / (scale + math.sqrt(scale**2 - 4 * a_total * t_total))
        mixed = result["tubes"]["aptamer_trigger"]
        self.assertAlmostEqual(mixed["summary"]["trigger_target_pct"], 100 * bound / t_total)
        self.assertAlmostEqual(mixed["summary"]["aptamer_target_pct"], 100 * bound / a_total)
        row = next(r for r in mixed["complexes"] if r["label"] == "A + T")
        self.assertAlmostEqual(row["concentration_m"], bound, delta=bound * 1e-10)
        self.assertAlmostEqual(row["complex_molecule_count"], bound * self.conditions.volume_l * 6.02214076e23,
                               delta=1.)
        self.backend.association = {"TT": 3e7}
        single = self.tubes()["tubes"]["trigger_only"]
        free = 2 * t_total / (1 + math.sqrt(1 + 8 * 3e7 * t_total))
        self.assertAlmostEqual(single["summary"]["trigger_free_pct"], 100 * free / t_total)
        dimer = next(r for r in single["complexes"] if r["strand_count"] == 2)
        self.assertEqual(dimer["trigger_molecule_count"], 2 * dimer["complex_molecule_count"])

    def test_both_trimer_stoichiometries_are_solved_and_reported(self):
        self.backend.association = {"AAT": 2e12, "ATT": 3e12}
        result = self.tubes()["tubes"]["aptamer_trigger"]
        orders = {"".join(r["strand_order"]): r["concentration_m"]
                  for r in result["equilibrium"]["ordered_complexes"]}
        self.assertAlmostEqual(orders["AAT"] / (orders["A"]**2 * orders["T"]), 2e12, delta=2e4)
        self.assertAlmostEqual(orders["ATT"] / (orders["A"] * orders["T"]**2), 3e12, delta=3e4)
        for label, order in (("2 A + T", "AAT"), ("A + 2 T", "ATT")):
            row = next(r for r in result["complexes"] if r["label"] == label)
            self.assertAlmostEqual(row["concentration_m"], orders[order], delta=1e-18)

    def test_identical_chemical_inputs_merge_before_equilibrium_and_split_origin_labels(self):
        self.backend.association = {"TT": 1e6}
        self.conditions = replace(self.conditions, aptamer_concentration_m=3e-6,
                                  trigger_concentration_m=7e-6, max_complex_size=2)
        result = self.tubes(sequences=("ACGT", "ACGT"))
        self.assertEqual(len(self.backend.calls), 2)
        mixed = result["tubes"]["aptamer_trigger"]
        self.assertTrue(mixed["equilibrium"]["chemical_aliases_merged"])
        physical_dimer = mixed["equilibrium"]["ordered_complexes"][1]["concentration_m"]
        fractions = {(2, 0): .3**2, (1, 1): 2 * .3 * .7, (0, 2): .7**2}
        for row in mixed["complexes"]:
            if row["strand_count"] == 2:
                self.assertAlmostEqual(row["concentration_m"], physical_dimer * fractions[
                    row["aptamer_count"], row["trigger_count"]], delta=1e-18)
        for species in ("aptamer", "trigger"):
            self.assertAlmostEqual(sum(r[f"{species}_fraction_pct"] for r in mixed["complexes"]), 100.)
        self.assertNotEqual(result["tubes"]["trigger_only"]["summary"]["trigger_free_pct"],
                            mixed["summary"]["trigger_free_pct"])

    def test_volume_affects_only_molecule_counts_and_dna_inputs_remain_dna(self):
        first = self.tubes()
        self.conditions = replace(self.conditions, volume_l=2 * self.conditions.volume_l)
        second = self.tubes()
        self.assertEqual(second["aptamer_sequence"], "ACGT")
        self.assertEqual(second["trigger_sequence"], "TGCA")
        self.assertTrue(all(s.material == "dna" for spec, _ in self.backend.calls for s in spec.strands))
        for name in first["tubes"]:
            for a, b in zip(first["tubes"][name]["complexes"], second["tubes"][name]["complexes"]):
                self.assertEqual(a["concentration_m"], b["concentration_m"])
                self.assertEqual(2 * a["complex_molecule_count"], b["complex_molecule_count"])

    def test_legacy_population_fixture_reports_identically(self):
        # Existing application's 60% free / 40% sequestered trigger fixture,
        # now obtained from a shared PF model rather than mocked concentrations.
        self.backend.association = {"TT": 1e-6 / (3e-6)**2}
        single = self.tubes()["tubes"]["trigger_only"]
        self.assertAlmostEqual(single["summary"]["trigger_free_pct"], 60.)
        self.assertAlmostEqual(single["summary"]["trigger_parasite_pct"], 40.)
        self.assertAlmostEqual(single["summary"]["parasite_particle_fraction_pct"], 25.)

    def test_cancellation_before_and_between_partitions_has_no_partial_populations(self):
        result = self.tubes(cancel_check=lambda: True)
        self.assertEqual(self.backend.calls, [])
        self.assertEqual(result["status"], "cancelled")
        result = self.tubes(cancel_check=lambda: len(self.backend.calls) >= 3)
        self.assertEqual(len(self.backend.calls), 3)
        self.assertEqual(result["tubes"], {})
        self.assertEqual(result["status"], "cancelled")

    def test_cancellation_after_first_tube_discards_partial_result(self):
        with patch.object(analysis, "solve_tube", wraps=solve_tube) as solve:
            result = self.tubes(cancel_check=lambda: solve.call_count > 0)
        self.assertEqual(solve.call_count, 1)
        self.assertEqual(result["status"], "cancelled")
        self.assertEqual(result["tubes"], {})

    def test_changed_partition_context_is_rejected_before_equilibrium(self):
        original = self.backend.partition
        self.backend.partition = lambda s, c: replace(original(s, c), conditions=Conditions(55, .15))
        with self.assertRaisesRegex(ValueError, "Backend changed"):
            self.tubes()

    def test_standard_state_must_be_explicit_and_match_the_backend(self):
        with self.assertRaises(TypeError):
            analysis.analyze_parasite_tubes(self.backend, "ACGT", "TGCA", self.thermo)
        with self.assertRaisesRegex(ValueError, "concentration standard"):
            analysis.analyze_parasite_tubes(self.backend, "ACGT", "TGCA", self.thermo,
                standard_state=replace(STANDARD, provenance_label="other"))
        self.assertEqual(self.backend.calls, [])




@unittest.skipUnless(os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY"), "Opt-in compiled native backend")
class RealNativeApplicationTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
        self.backend = NativeHybridBackend(os.environ["APTASWITCH_HYBRID_NATIVE_LIBRARY"])
        self.thermo = ThermoConditions(temperature_c=37, sodium_m=1., switch_material="dna", trigger_material="dna")
        self.standard = self.backend.standard_state

    def test_real_minimum_uses_original_coordinates_and_typed_record(self):
        result = analysis.analyze_system(self.backend, "ACGT", "ACGT", self.thermo)
        spec = ComplexSpec((Strand("ACGT", "dna", "A"), Strand("ACGT", "dna", "T")))
        native = self.backend.minimum(spec, Conditions(37, 1.))
        self.assertEqual(result["mfe_structure"], native.structure)
        self.assertEqual(result["mfe_energy_kcal_mol"], native.energy_kcal_mol)
        self.assertTrue(result["paired_trigger_indices"])

    def test_real_three_strand_catalog_preserves_totals_and_matches_direct_solution(self):
        conditions = TubeAnalysisConditions(max_complex_size=3)
        for temperature, sodium in ((37., 1.), (31., .15)):
            with self.subTest(temperature=temperature, sodium=sodium):
                thermo = replace(self.thermo, temperature_c=temperature, sodium_m=sodium)
                result = analysis.analyze_parasite_tubes(self.backend, "AC", "GT", thermo,
                                                        conditions, standard_state=self.standard)
                self.assertEqual(result["status"], "completed")
                for name, tube in result["tubes"].items():
                    for label in (("trigger",) if name == "trigger_only" else ("aptamer", "trigger")):
                        self.assertAlmostEqual(sum(r[f"{label}_fraction_pct"] for r in tube["complexes"]),
                                               100., places=7)
                mixed = result["tubes"]["aptamer_trigger"]
                self.assertEqual(mixed["equilibrium"]["ordered_complex_count"], 9)
                species = (Strand("AC", "dna", "A"),
                           Strand("GT", "dna", "T"))
                records = tuple(self.backend.partition(spec, Conditions(temperature, sodium))
                                for spec in enumerate_complexes(species, 3))
                direct = solve_tube(species, (5e-6, 5e-6), records, standard_state=self.standard)
                actual = mixed["equilibrium"]["ordered_complexes"]
                for row, record, concentration in zip(actual, records, direct.concentrations_molar):
                    self.assertEqual(row["strand_order"], [s.name for s in record.complex_spec.strands])
                    self.assertEqual(row["concentration_m"], float(concentration))
                self.assertEqual(json.loads(json.dumps(result, allow_nan=False)), result)


if __name__ == "__main__":
    unittest.main()
