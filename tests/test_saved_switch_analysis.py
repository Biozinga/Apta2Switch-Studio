"""Full reporter validation must preserve sequence identity and tube mass."""

from dataclasses import asdict, dataclass
from copy import deepcopy
import json
import os
import unittest
from unittest.mock import patch

from aptaswitch_core import saved_switch_analysis as analysis
from aptaswitch_core.architecture import architecture_layout
from aptaswitch_core.models import SwitchArchitecture


def snapshot(material="experimental-Sugimoto95-Vienna04-loopmean", switch_chem="rna", trigger_chem="dna"):
    arch = SwitchArchitecture(toehold_length=2, stem_length=2, upper_stem_length=2,
                              upper_stem2_length=0, rbs_prefix_length=1, rbs_loop="AGG", frame_linker="A")
    layout = architecture_layout(arch)
    switch = list("G" * layout["reporter_start"])
    switch[layout["aug_start"]:layout["aug_start"] + 3] = "AUG"
    return {
        "name": "Test retained candidate", "candidate": {"switch_seq": "".join(switch) + "CCCGGG", "trigger_seq": "AACCGGTT"},
        "request": {
            "architecture": asdict(arch), "reporter": {"sequence_after_start_rna": "CCCGGG", "full_sequence_after_start_rna": "CCCGGGAUCUAA"},
            "sequences": {"aptamer_dna": "GTCAGCT", "trigger_dna": "CCGG"},
            "thermo": {"material": material, "switch_material": switch_chem, "trigger_material": trigger_chem,
                       "sodium_m": .15, "magnesium_m": 0., "temperature_c": 29., "concentration_m": 5e-6},
        },
    }


def config(**kwargs):
    return kwargs


@dataclass(frozen=True)
class Strand:
    sequence: str
    name: str


@dataclass(frozen=True)
class Complex:
    strands: tuple


class SavedSwitchAnalysisTests(unittest.TestCase):
    def run_fake(self, *, candidate=None, settings=None, cancelled=None, progress=None, backend=None):
        from test_thermodynamics_context_analysis import FixedWeightBackend
        backend = backend or FixedWeightBackend()
        result = analysis.analyze_saved_switch(candidate or snapshot(), None, settings or config(),
                                               cancel_check=cancelled, progress_callback=progress,
                                               backend=backend)
        return result, backend

    def test_full_reporter_extends_exact_saved_candidate_without_duplicate_start(self):
        data = snapshot()
        context = analysis.build_saved_switch_context(data, config())
        self.assertEqual(context["switch_sequence"], data["candidate"]["switch_seq"] + "AUCUAA")
        self.assertEqual(context["switch_sequence"].count("AUG"), 1)
        self.assertEqual(context["trigger_sequence"], "AACCGGTT")
        self.assertNotEqual(context["trigger_sequence"], data["request"]["sequences"]["trigger_dna"])

    def test_other_switch_trigger_chemistries_are_rejected(self):
        for switch, trigger in (("dna", "dna"), ("dna", "rna"), ("rna", "rna")):
            with self.subTest(switch=switch, trigger=trigger), self.assertRaises(ValueError):
                analysis.build_saved_switch_context(snapshot(switch_chem=switch, trigger_chem=trigger), config())

    def test_reporter_from_saved_library_requires_no_second_sequence_input(self):
        data = snapshot()
        data["request"]["reporter"]["full_sequence_after_start_rna"] = "CCCGGGAUCUAA"
        result = analysis.build_saved_switch_context(data, {})
        self.assertEqual(result["reporter_source"], "saved_full_reporter")
        self.assertEqual(result["switch_sequence"], data["candidate"]["switch_seq"] + "AUCUAA")
        result, _ = self.run_fake(candidate=data, settings={"temperature_values_c": []})
        self.assertEqual(result["context"]["reporter_source"], "saved_full_reporter")

    def test_reporter_always_comes_from_saved_candidate_not_context_override(self):
        data = snapshot()
        context = analysis.build_saved_switch_context(data, {'full_reporter_sequence_after_start': 'CCCGGGAAAAAA'})
        self.assertEqual(context['reporter_sequence_after_start'], 'CCCGGGAUCUAA')
        self.assertEqual(context['reporter_source'], 'saved_full_reporter')
        normalized = analysis.validate_analysis_config(data, {'full_reporter_sequence_after_start': 'CCCGGGAAAAAA'})
        self.assertNotIn('full_reporter_sequence_after_start', normalized)

    def test_total_scan_and_full_transcript_limits_are_explicit(self):
        data = snapshot()
        valid = analysis.validate_analysis_config(data, config(temperature_values_c=list(range(25,50)), sodium_values_m=[.12 + i / 100 for i in range(25)]))
        self.assertEqual(1 + len(valid['temperature_values_c']) + len(valid['sodium_values_m']), 51)
        with self.assertRaisesRegex(ValueError, "25 points"):
            analysis.validate_analysis_config(data, config(temperature_values_c=list(range(26))))
        with self.assertRaisesRegex(ValueError, "2000 nucléotides"):
            data["request"]["reporter"]["full_sequence_after_start_rna"] = "CCCGGG" + "A" * 2000
            analysis.build_saved_switch_context(data, config())

    def test_missing_mismatching_or_duplicate_aug_reporter_rejected(self):
        for value in (None, "", "AUGCCCGGGAUCUAA", "CCGGGAUCUAA", "NNNN"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                data = snapshot()
                data["request"]["reporter"]["full_sequence_after_start_rna"] = value
                analysis.build_saved_switch_context(data, config())
        data = snapshot()
        data["candidate"]["switch_seq"] += "A"
        with self.assertRaisesRegex(ValueError, "architecture"):
            analysis.build_saved_switch_context(data, config())

    def test_all_conditions_validated_before_native_import(self):
        invalid = [config(sodium_m=.01), config(magnesium_m=.001), config(temperature_c=float("nan")),
                   config(switch_concentration_m=0), config(trigger_concentration_m=True),
                   config(temperature_values_c=list(range(26))),
                   *(config(max_complex_size=size) for size in (1, 4, 5, 6, 3., True)),
                   config(magnesium_values_m=[0]), config(include_aptamer="yes"),
                   config(sodium_values_m=[.2] * 26), config(temperature_values_c=list(range(25,50)), sodium_values_m=[.12 + i / 100 for i in range(25)], magnesium_values_m=[0.])]
        for settings in invalid:
            from test_thermodynamics_context_analysis import FixedWeightBackend
            backend = FixedWeightBackend()
            with self.subTest(config=settings), self.assertRaises(ValueError):
                analysis.analyze_saved_switch(snapshot(), None, settings, backend=backend)
            self.assertEqual(backend.calls, [])

    def test_only_mixed_model_is_accepted_and_magnesium_is_always_excluded(self):
        for material in ("dna04", "rna06", "rna95", "rna99"):
            with self.subTest(material=material), self.assertRaises(ValueError):
                analysis.validate_analysis_config(snapshot(material), config())
        for settings in (config(magnesium_m=.001), config(magnesium_values_m=[0.])):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                analysis.validate_analysis_config(snapshot(), settings)
        self.assertEqual(analysis.model_limits("experimental-Sugimoto95-Vienna04-loopmean")["sodium_min_m"], .12)

    def test_declared_native_domain_is_checked_before_calculation(self):
        from test_thermodynamics_context_analysis import FixedWeightBackend
        backend = FixedWeightBackend()
        backend.condition_limits = {**backend.condition_limits, "sodium_min_m": 1.}
        with self.assertRaisesRegex(ValueError, "sodium_m"):
            self.run_fake(backend=backend)
        self.assertEqual(backend.calls, [])

    def test_default_tubes_and_actual_native_identity(self):
        result, backend = self.run_fake()
        self.assertEqual(result["status"], "completed", result.get("error"))
        self.assertEqual(result["thermodynamics"], asdict(backend.provenance))
        self.assertEqual(result["context"]["material"], backend.provenance.model)

        for state, names in (("off", {"Switch"}), ("on", {"Switch", "Trigger"})):
            population = result["baseline"][state]["population"]
            self.assertEqual(set(population["total_strand_concentrations_m"]), names)
            self.assertTrue(population["mass_balance_checked"])
        self.assertEqual(result["baseline"]["off"]["mfe"]["strand_order"], ["Switch"])
        self.assertEqual(result["baseline"]["on"]["mfe"]["strand_order"], ["Trigger", "Switch"])
        json.dumps(result, allow_nan=False)

    def test_optional_aptamer_preserves_dna_chemistry_and_equal_concentration(self):
        from test_thermodynamics_context_analysis import FixedWeightBackend
        backend = FixedWeightBackend()
        with patch.object(backend, "partition", wraps=backend.partition) as partition:
            result, _ = self.run_fake(candidate=snapshot(),
                settings=config(include_aptamer=True, aptamer_concentration_m=7e-6), backend=backend)
        self.assertEqual(result["status"], "completed", result.get("error"))
        for state in ("off", "on"):
            self.assertEqual(result["baseline"][state]["population"]["total_strand_concentrations_m"]["Aptamer"], 7e-6)
        strands = [strand for call in partition.call_args_list for strand in call.args[0].strands]
        self.assertTrue(any(s.name == "Aptamer" for s in strands))
        self.assertTrue(all(s.material == ("rna" if s.name == "Switch" else "dna") for s in strands))

    def test_complete_catalog_size_and_baseline_minima_preserve_full_sequences(self):
        for size, count in ((2, 5), (3, 9)):
            with self.subTest(size=size):
                candidate = snapshot()
                before = deepcopy(candidate)
                result, _ = self.run_fake(candidate=candidate, settings=config(max_complex_size=size))
                self.assertEqual(candidate, before)
                self.assertEqual(result["status"], "completed", result.get("error"))
                rows = result["baseline"]["on"]["population"]["complexes"]
                self.assertEqual(len(rows), count)
                self.assertEqual(len({row["complex_key"] for row in rows}), count)
                for row in rows:
                    mfe = row["mfe"]
                    if len(row["strand_order"]) > 2:
                        self.assertIsNone(mfe)
                        self.assertEqual(row["structure_status"], "empty_ensemble")
                    else:
                        self.assertEqual(mfe["strand_order"], row["strand_order"])
                        self.assertEqual(mfe["sequences"], [result["context"][name.lower() + "_sequence"] for name in row["strand_order"]])
                        self.assertEqual(list(map(len, mfe["structure"].split("+"))), list(map(len, mfe["sequences"])))

    def test_new_analysis_defaults_to_three_strands_without_changing_the_saved_request(self):
        candidate = snapshot()
        candidate['request']['thermo']['validation_max_complex_size'] = 5
        before = deepcopy(candidate)
        result, _ = self.run_fake(candidate=candidate)
        self.assertEqual(result['status'], 'completed', result.get('error'))
        self.assertEqual(result['config']['max_complex_size'], 3)
        self.assertEqual(len(result['baseline']['on']['population']['complexes']), 9)
        self.assertTrue(all(len(row['strand_order']) <= 3
                            for row in result['baseline']['on']['population']['complexes']))
        self.assertEqual(candidate, before)

    def test_population_uses_species_denominators_and_distinct_cyclic_orders(self):
        strands = {name: Strand("AC", name) for name in ("Switch", "Trigger", "Aptamer")}
        def population(rows, totals, size=5):
            return analysis._population({Complex(tuple(strands[n] for n in names)): value for names, value in rows}, totals, size)
        p = population([(("Switch",), 4e-6), (("Trigger",), 1e-6), (("Trigger", "Switch"), 1e-6)],
                       {"Switch": 5e-6, "Trigger": 2e-6})
        self.assertAlmostEqual(p["target_switch_fraction_pct"], 20.)
        self.assertAlmostEqual(p["target_trigger_fraction_pct"], 50.)
        for first, second in (
            (("Switch", "Trigger", "Aptamer"), ("Switch", "Aptamer", "Trigger")),
            (("Switch", "Switch", "Trigger", "Trigger"), ("Switch", "Trigger", "Switch", "Trigger")),
            (("Switch", "Switch", "Switch", "Trigger", "Trigger"), ("Switch", "Switch", "Trigger", "Switch", "Trigger")),
        ):
            species = set(first)
            rows = [(first, .4e-6), (second, .7e-6)] + [((n,), 5e-6 - first.count(n) * 1.1e-6) for n in species]
            report = population(rows, {n: 5e-6 for n in species})
            complexes = [row for row in report["complexes"] if len(row["strand_order"]) > 1]
            self.assertEqual([row["concentration_m"] for row in complexes], [.7e-6, .4e-6])
            self.assertEqual(complexes[0]["composition_label"], complexes[1]["composition_label"])
            self.assertNotEqual(complexes[0]["complex_key"], complexes[1]["complex_key"])
            self.assertNotEqual(complexes[0]["label"], complexes[1]["label"])

    def test_rotations_merge_but_invalid_complexes_and_mass_fail(self):
        names = ("Switch", "Trigger", "Aptamer")
        strands = {name: Strand("AC", name) for name in names}
        total = {name: 5e-6 for name in names}
        reports = [analysis._population({Complex(tuple(strands[name] for name in names[i:] + names[:i])): 5e-6}, total, 5)
                   for i in range(len(names))]
        self.assertEqual(reports, [reports[0]] * 3)
        for values, totals in (
            ({Complex((strands["Switch"],) * 6): 5e-6 / 6}, {"Switch": 5e-6}),
            ({Complex((strands["Switch"],)): 1e-6}, {"Switch": 5e-6}),
            ({Complex((strands["Switch"],)): float("nan")}, {"Switch": 5e-6}),
            ({Complex((strands["Switch"],)): -1e-6}, {"Switch": 5e-6}),
        ):
            with self.assertRaises(RuntimeError):
                analysis._population(values, totals, 5)

    def test_reference_shared_between_scans_and_only_one_axis_changes(self):
        result, backend = self.run_fake(settings=config(temperature_values_c=[40., 29., 25.], sodium_values_m=[.15, .3]))
        self.assertEqual(result["status"], "completed", result.get("error"))
        self.assertEqual(len(backend.calls), 9 * 4)
        self.assertEqual([p["conditions"]["temperature_c"] for p in result["sweeps"]["temperature_c"]], [25., 29., 40.])
        self.assertIs(result["sweeps"]["temperature_c"][1], result["baseline"])
        self.assertIs(result["sweeps"]["sodium_m"][0], result["baseline"])
        for point in result["sweeps"]["sodium_m"]:
            self.assertEqual(point["conditions"]["temperature_c"], 29.)
            self.assertEqual(point["conditions"]["magnesium_m"], 0.)

    def test_failed_or_cancelled_competing_minimum_never_publishes_half_reference(self):
        from test_thermodynamics_context_analysis import FixedWeightBackend
        backend = FixedWeightBackend()
        minimum = backend.minimum_structures
        calls = []
        def fail(spec, conditions, **kwargs):
            calls.append(spec)
            if len(calls) == 3:
                raise RuntimeError("competing complex MFE failure")
            return minimum(spec, conditions, **kwargs)
        with patch.object(backend, "minimum_structures", side_effect=fail):
            result, _ = self.run_fake(backend=backend)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["error"], "competing complex MFE failure")
        self.assertIsNone(result["baseline"])
        with patch.object(backend, "minimum_structures", wraps=minimum) as compute:
            result, _ = self.run_fake(backend=backend, cancelled=lambda: compute.call_count >= 3)
        self.assertEqual(result["status"], "cancelled")
        self.assertIsNone(result["baseline"])

    def test_failed_later_point_and_cancel_keep_completed_baseline(self):
        from test_thermodynamics_context_analysis import FixedWeightBackend
        backend = FixedWeightBackend()
        backend.fail_temperature = 40.
        result, _ = self.run_fake(backend=backend, settings=config(temperature_values_c=[40.]))
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["baseline"]["conditions"]["temperature_c"], 29.)
        self.assertEqual(result["sweeps"]["temperature_c"], [])
        cancelled = {"value": False}
        def progress(done, total, message):
            if done == 1:
                cancelled["value"] = True
        result, _ = self.run_fake(settings=config(temperature_values_c=[40.]), progress=progress,
                                  cancelled=lambda: cancelled["value"])
        self.assertEqual(result["status"], "cancelled")
        self.assertIsNotNone(result["baseline"])
        self.assertEqual(result["sweeps"]["temperature_c"], [])
        result, backend = self.run_fake(cancelled=lambda: True)
        self.assertEqual(backend.calls, [])
        self.assertIsNone(result["baseline"])


@unittest.skipUnless(os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY"), "Native hybrid runtime required")
class SavedSwitchNativeTests(unittest.TestCase):
    def test_native_short_fixed_transcripts_and_independent_salt_temperature_scans(self):
        for switch_chem, trigger_chem in (("rna", "dna"),):
            with self.subTest(switch=switch_chem, trigger=trigger_chem):
                result = analysis.analyze_saved_switch(snapshot("experimental-Sugimoto95-Vienna04-loopmean", switch_chem, trigger_chem), None,
                    config(max_complex_size=2, include_aptamer=True, temperature_values_c=[37.], sodium_values_m=[.3]))
                self.assertEqual(result["status"], "completed", result.get("error"))
                self.assertTrue(result["thermodynamics"])

                for point in [result["baseline"], *result["sweeps"]["temperature_c"], *result["sweeps"]["sodium_m"]]:
                    for state in ("off", "on"):
                        self.assertTrue(point[state]["population"]["mass_balance_checked"])
                        self.assertEqual(point[state]["population"]["total_strand_concentrations_m"]["Aptamer"], 5e-6)
                    self.assertGreaterEqual(point["metrics"]["on_yield_pct"], 0.)
                    self.assertLessEqual(point["metrics"]["on_yield_pct"], 100.)
                self.assertEqual(result["sweeps"]["temperature_c"][0]["conditions"]["sodium_m"], .15)
                self.assertEqual(result["sweeps"]["sodium_m"][0]["conditions"]["temperature_c"], 29.)
                json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
