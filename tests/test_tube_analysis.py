"""Engine-independent tube reporting and native public entry-point contracts."""
from dataclasses import dataclass, replace
import math
import unittest
from unittest.mock import Mock, patch
from aptaswitch_core import tube_analysis
from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.thermodynamics.tube_reporting import TubeAnalysisConditions, _population
from tests.test_thermodynamics_application_analysis import FixtureBackend, STANDARD

@dataclass(frozen=True)
class Strand:
    name: str

@dataclass(frozen=True)
class Complex:
    strands: tuple[Strand, ...]

class TubeAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.conditions = TubeAnalysisConditions()

    def analyze(self, *, populations=None, conditions=None):
        conditions = conditions or self.conditions
        values = populations if populations is not None else {
            "trigger_only": [("T", conditions.trigger_concentration_m)],
            "aptamer_trigger": [("A", conditions.aptamer_concentration_m), ("T", conditions.trigger_concentration_m)],
        }
        tubes = {name: _population({Complex(tuple(Strand(s) for s in names)): amount
                                   for names, amount in rows}, conditions, mixed=name=="aptamer_trigger")
                 for name, rows in values.items()}
        return {"tubes": tubes}, None, None

    def test_dimers_count_twice_as_trigger_strands_but_once_as_particles(self):
        result, _, _ = self.analyze(populations={
            "trigger_only": [("T", 3e-6), ("TT", 1e-6)],
            "aptamer_trigger": [("A", 2e-6), ("T", 2e-6), ("AT", 1e-6), ("AA", .5e-6), ("AAT", .5e-6), ("TT", .75e-6)],
        })
        single = result["tubes"]["trigger_only"]
        self.assertAlmostEqual(single["summary"]["trigger_free_pct"], 60.)
        self.assertAlmostEqual(single["summary"]["trigger_parasite_pct"], 40.)
        self.assertAlmostEqual(single["summary"]["parasite_particle_fraction_pct"], 25.)
        dimer = next(row for row in single["complexes"] if row["trigger_count"] == 2)
        self.assertEqual(dimer["label"], "2 T")
        self.assertAlmostEqual(dimer["complex_molecule_count"], 1e-6 * 20e-6 * 6.02214076e23)
        self.assertEqual(dimer["trigger_molecule_count"], 2 * dimer["complex_molecule_count"])
        self.assertAlmostEqual(sum(row["particle_fraction_pct"] for row in single["complexes"]), 100.)
        mixed = result["tubes"]["aptamer_trigger"]
        self.assertEqual(next(row for row in mixed["complexes"] if row["label"] == "A + T")["category"], "target")
        for species in ("trigger", "aptamer"):
            self.assertAlmostEqual(mixed["summary"][f"{species}_free_pct"], 40.)
            self.assertAlmostEqual(mixed["summary"][f"{species}_target_pct"], 20.)
            self.assertAlmostEqual(mixed["summary"][f"{species}_parasite_pct"], 40.)


    def test_aptamer_oligomers_are_parasite_particles_but_do_not_capture_trigger(self):
        result, _, _ = self.analyze(populations={
            "trigger_only": [("T", 5e-6)],
            "aptamer_trigger": [("AA", 2.5e-6), ("T", 5e-6)],
        })
        summary = result["tubes"]["aptamer_trigger"]["summary"]
        self.assertEqual(summary["trigger_parasite_pct"], 0.)
        self.assertEqual(summary["trigger_free_pct"], 100.)
        self.assertEqual(summary["aptamer_parasite_pct"], 100.)
        self.assertAlmostEqual(summary["parasite_particle_fraction_pct"], 100. / 3.)


    def test_distinct_rotations_with_same_stoichiometry_are_aggregated(self):
        result, _, _ = self.analyze(populations={
            "trigger_only": [("T", 5e-6)],
            "aptamer_trigger": [("A", 4e-6), ("T", 4.5e-6), ("AAT", .2e-6), ("ATA", .3e-6)],
        })
        rows = result["tubes"]["aptamer_trigger"]["complexes"]
        self.assertEqual(len(rows), 3)
        row = next(row for row in rows if row["strand_count"] == 3)
        self.assertEqual(row["label"], "2 A + T")
        self.assertAlmostEqual(row["concentration_m"], .5e-6)
        self.assertAlmostEqual(row["trigger_fraction_pct"], 10.)
        self.assertAlmostEqual(row["aptamer_fraction_pct"], 20.)


    def test_asymmetric_concentrations_use_each_species_own_denominator(self):
        conditions = replace(self.conditions, aptamer_concentration_m=1e-6, trigger_concentration_m=10e-6)
        result, _, _ = self.analyze(conditions=conditions, populations={
            "trigger_only": [("T", 10e-6)],
            "aptamer_trigger": [("AT", 1e-6), ("T", 9e-6)],
        })
        summary = result["tubes"]["aptamer_trigger"]["summary"]
        self.assertAlmostEqual(summary["aptamer_target_pct"], 100.)
        self.assertAlmostEqual(summary["trigger_target_pct"], 10.)
        self.assertAlmostEqual(summary["trigger_free_pct"], 90.)


    def test_volume_changes_counts_not_concentrations_or_fractions(self):
        first, _, _ = self.analyze()
        second, _, _ = self.analyze(conditions=replace(self.conditions, volume_l=2 * self.conditions.volume_l))
        row1 = first["tubes"]["trigger_only"]["complexes"][0]
        row2 = second["tubes"]["trigger_only"]["complexes"][0]
        self.assertEqual(row1["concentration_m"], row2["concentration_m"])
        self.assertEqual(row1["trigger_fraction_pct"], row2["trigger_fraction_pct"])
        self.assertEqual(2 * row1["complex_molecule_count"], row2["complex_molecule_count"])


    def test_invalid_concentrations_volume_and_complex_size_are_rejected(self):
        for name in ("aptamer_concentration_m", "trigger_concentration_m", "volume_l"):
            for value in (0, -1, math.nan, math.inf, -math.inf, None, True, "bad"):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    replace(self.conditions, **{name: value})
        for value in (0, 4, 5, 6, 2.5, 3., True, None, "3"):
            with self.subTest(max_complex_size=value), self.assertRaises(ValueError):
                replace(self.conditions, max_complex_size=value)


    def test_nonfinite_or_significantly_negative_engine_concentrations_are_errors(self):
        for value in (math.nan, math.inf, -math.inf, -1e-12, "bad"):
            with self.subTest(value=value), self.assertRaises(RuntimeError):
                self.analyze(populations={
                    "trigger_only": [("T", 5e-6), ("TT", value)],
                    "aptamer_trigger": [("A", 5e-6), ("T", 5e-6)],
                })


    def test_only_tiny_negative_roundoff_is_clipped(self):
        result, _, _ = self.analyze(populations={
            "trigger_only": [("T", 5e-6), ("TT", -1e-20)],
            "aptamer_trigger": [("A", 5e-6), ("T", 5e-6)],
        })
        row = next(row for row in result["tubes"]["trigger_only"]["complexes"] if row["trigger_count"] == 2)
        self.assertEqual(row["concentration_m"], 0.)


    def test_mass_imbalance_and_invalid_complexes_are_errors(self):
        for values in ([], [("T", 4e-6)], [("A", 5e-6)], [("Q", 5e-6)], [("TTTTTT", 5e-6 / 6)], [("", 5e-6)]):
            with self.subTest(values=values), self.assertRaises(RuntimeError):
                self.analyze(populations={
                    "trigger_only": values,
                    "aptamer_trigger": [("A", 5e-6), ("T", 5e-6)],
                })


class NativeTubeRouteTests(unittest.TestCase):
    def test_default_runtime_is_selected_once_and_catalogs_conserve_material(self):
        backend = FixtureBackend({'AT': 1e6, 'TT': 2e5})
        backend.standard_state = STANDARD
        with patch('aptaswitch_core.thermodynamics.setup.get_backend', return_value=backend) as load:
            result = tube_analysis.analyze_parasite_tubes('/old/native', 'ACGT', 'TGCA', ThermoConditions(switch_material='dna'))
        load.assert_called_once_with(None)
        self.assertEqual(len(backend.calls), 9)  # all cyclic orders through three strands
        self.assertEqual(result['conditions']['max_complex_size'], 3)
        self.assertTrue(all(len(spec.strands) <= 3 for spec, _ in backend.calls))
        self.assertEqual(result['status'], 'completed')
        for tube in result['tubes'].values():
            for species in ('aptamer', 'trigger'):
                expected = tube['total_strand_concentrations_m'][species]
                recovered = sum(row[species+'_count'] * row['concentration_m'] for row in tube['complexes'])
                self.assertAlmostEqual(recovered, expected, delta=1e-12)

    def test_cancelled_catalog_does_not_publish_partial_tubes(self):
        backend = FixtureBackend()
        backend.standard_state = STANDARD
        result = tube_analysis.analyze_parasite_tubes(None, 'ACGT', 'TGCA', ThermoConditions(switch_material='dna'),
                  backend=backend, cancel_check=Mock(side_effect=[False, False, True]))
        self.assertEqual(result['status'], 'cancelled')
        self.assertEqual(result['tubes'], {})
        self.assertEqual(len(backend.calls), 1)
