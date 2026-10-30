"""Exact joint structure probabilities and optional native/oracle comparisons."""
from dataclasses import replace
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, EmptyEnsembleError, PartitionResult, Strand


LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")
HAS_ORACLE = importlib.util.find_spec("RNA") is not None


def spec(*items):
    return ComplexSpec(tuple(Strand(sequence, material, f"s{i}")
                             for i, (sequence, material) in enumerate(items)))


class StructureProbabilityGuardsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        path = Path(self.temp.name) / "test.bin"
        path.write_bytes(b"not an executable; manifest checks only")
        prefix = "src/aptaswitch_core/thermodynamics/"
        sources = [prefix + "native/" + name for name in (
            "connected_dp.cpp", "connected_dp.hpp", "hybrid_provider.cpp", "hybrid_provider.hpp",
            "hybrid_api.cpp", "hybrid_api.hpp")] + [prefix + "hybrid_motifs.py"]
        manifest = dict(abi=1, status="built", model="experimental-Sugimoto95-Vienna04-loopmean",
            machine=platform.machine(), sources={name: "a" * 64 for name in sources},
            vienna_version="2.7.2", single_rna_solver="rnalib2-pf-bpp-sample-mfe-v1",
            single_rna_minimum_convention="turner2004-centikcal-salt-rounded-v1",
            native_library_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            generated_parameters_sha256="b" * 64, vienna_params_c_sha256="c" * 64,
            vienna_static_library_sha256="d" * 64, vienna_config_h_sha256="e" * 64)
        path.with_suffix(".bin.json").write_text(json.dumps(manifest))
        self.backend = NativeHybridBackend(path)
        self.current, self.conditions = spec(("AA", "rna"), ("T", "dna")), Conditions()
        self.partition = PartitionResult(self.current, self.conditions, self.backend.provenance,
                                         2., -1.2, .6)

    def test_rejects_invalid_requests_before_loading_runtime(self):
        calls = [
            (None, self.conditions, "(.+)"),
            (self.current, {}, "(.+)"),
            (self.current, Conditions(celsius=24), "(.+)"),
            (self.current, self.conditions, True),
            (self.current, self.conditions, "...+"),
            (self.current, self.conditions, "(.&+)"),
            (self.current, self.conditions, "((("),
            (self.current, self.conditions, PairingConstraints.unconstrained(3)),
            (self.current, self.conditions, parse_structure("..", (2,))),
        ]
        with patch.object(self.backend, "_get_runtime") as runtime:
            for arguments in calls:
                with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                    self.backend.structure_probability(*arguments)
            runtime.assert_not_called()

    def test_reused_partition_requires_exact_unrestricted_context_before_load(self):
        originals = self.current.strands
        alternatives = [
            True,
            replace(self.partition, complex_spec=spec(("AU", "rna"), ("T", "dna"))),
            replace(self.partition, complex_spec=spec(("AA", "dna"), ("T", "dna"))),
            replace(self.partition, complex_spec=ComplexSpec(tuple(reversed(originals)))),
            replace(self.partition, complex_spec=ComplexSpec((replace(originals[0], name="alias"), originals[1]))),
            replace(self.partition, conditions=Conditions(celsius=38)),
            replace(self.partition, provenance=replace(self.backend.provenance, parameter_digest="sha256:other")),
            replace(self.partition, constraints=parse_structure("(.+)", (2, 1))),
        ]
        with patch.object(self.backend, "_get_runtime") as runtime:
            for denominator in alternatives:
                with self.subTest(partition=denominator), self.assertRaises(ValueError):
                    self.backend.structure_probability(self.current, self.conditions, "(.+)",
                                                       partition=denominator)
            runtime.assert_not_called()

    def test_denominator_reuse_does_not_request_an_unrestricted_partition(self):
        state = parse_structure("(.+)", (2, 1))
        target = replace(self.partition, log_partition=1., free_energy_kcal_mol=-.6, constraints=state)
        with patch.object(self.backend, "partition", return_value=target) as partition:
            result = self.backend.structure_probability(self.current, self.conditions, state,
                                                       partition=self.partition)
        partition.assert_called_once_with(self.current, self.conditions, state)
        self.assertAlmostEqual(result, math.exp(-1.), places=15)

    def test_reused_partition_rejects_a_different_thermal_energy_convention(self):
        incompatible = replace(self.partition, thermal_energy_kcal_mol=.7, free_energy_kcal_mol=-1.4)
        with patch.object(self.backend, "partition", return_value=self.partition):
            with self.assertRaisesRegex(ValueError, "thermal energy"):
                self.backend.structure_probability(self.current, self.conditions, "(.+)",
                                                   partition=incompatible)

    def test_empty_total_is_undefined_even_for_a_disconnected_target(self):
        empty = replace(self.partition, log_partition=-math.inf, free_energy_kcal_mol=math.inf)
        with patch.object(self.backend, "_get_runtime") as runtime:
            with self.assertRaises(EmptyEnsembleError):
                self.backend.structure_probability(self.current, self.conditions, "..+.", partition=empty)
            runtime.assert_not_called()

    def test_inconsistent_partition_ratio_is_not_silently_clamped(self):
        target = replace(self.partition, log_partition=3., free_energy_kcal_mol=-1.8)
        with patch.object(self.backend, "partition", return_value=target):
            with self.assertRaisesRegex(ArithmeticError, "exceeds"):
                self.backend.structure_probability(self.current, self.conditions, "(.+)",
                                                   partition=self.partition)


@unittest.skipUnless(LIBRARY, "explicit native hybrid library required")
class NativeStructureProbabilityTests(unittest.TestCase):
    def setUp(self):
        self.backend = NativeHybridBackend(LIBRARY)

    def test_joint_probability_is_not_product_of_marginals_and_reuses_partition(self):
        current, conditions = spec(("AA", "rna"), ("T", "dna")), Conditions()
        analysis = self.backend.analyze(current, conditions, need_pairs=True)
        with patch.object(self.backend, "partition", wraps=self.backend.partition) as partition:
            probability = self.backend.structure_probability(current, conditions, "(.+)",
                                                             partition=analysis.partition)
        self.assertEqual(partition.call_count, 1)
        self.assertTrue(partition.call_args.args[2].complete)
        self.assertAlmostEqual(probability, .5, places=12)
        self.assertNotAlmostEqual(probability, analysis.pairs.matrix[0][2] * analysis.pairs.matrix[1][1])
        self.assertAlmostEqual(self.backend.structure_probability(current, conditions, ".(&)"), .5, places=12)
        self.assertAlmostEqual(self.backend.structure_probability(current, conditions, "(.)"), .5, places=12)

    def test_impossible_and_disconnected_targets_are_zero_but_empty_total_raises(self):
        current, conditions = spec(("GCAAA", "rna"), ("TTTGC", "dna")), Conditions()
        denominator = self.backend.partition(current, conditions)
        for structure in (".....+.....", "(....+)...."):
            with self.subTest(structure=structure):
                self.assertEqual(self.backend.structure_probability(current, conditions, structure,
                                                                   partition=denominator), 0.)
        current = spec(("GCGAAAGC", "rna"))
        self.assertEqual(self.backend.structure_probability(current, conditions, "()......"), 0.)
        self.assertEqual(self.backend.structure_probability(spec(("AAA", "rna")), conditions, "..."), 1.)
        with self.assertRaises(EmptyEnsembleError):
            self.backend.structure_probability(spec(("A", "rna"), ("A", "dna")), conditions, ".+.")

    @unittest.skipUnless(HAS_ORACLE, "optional Python RNA oracle unavailable")
    def test_all_small_ensemble_probabilities_match_exhaustive_energy_oracle(self):
        from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend
        from aptaswitch_core.thermodynamics.hybrid_motifs import HybridMotifModel
        model = HybridMotifModel(sodium_correction="published-pair")
        oracle = ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
            thermal_energy=model.thermal_energy, provenance=model.provenance, min_hairpin_size=3)
        cases = [spec(("GCGAAAGC", "rna")), spec(("GCGAAAGC", "dna")),
                 spec(("GCAAA", "rna"), ("TTTGC", "dna")),
                 spec(("TTTGC", "dna"), ("GCAAA", "rna")),
                 spec(("GC", "rna"), ("GC", "dna"), ("GC", "rna"))]
        for current in cases:
            for conditions in (Conditions(celsius=25), Conditions(celsius=55, sodium=1)):
                with self.subTest(order=current.ordered_key, conditions=conditions):
                    expected = oracle.evaluate(current, conditions)
                    denominator = self.backend.partition(current, conditions)
                    actual = [self.backend.structure_probability(current, conditions, state,
                                                                 partition=denominator)
                              for state in expected.structures]
                    self.assertAlmostEqual(math.fsum(actual), 1., places=10)
                    lengths = tuple(len(strand.sequence) for strand in current.strands)
                    for state, probability, wanted in zip(expected.structures, actual, expected.probabilities):
                        with self.subTest(structure=format_structure(state, lengths)):
                            self.assertAlmostEqual(probability, wanted, places=10)


if __name__ == "__main__":
    unittest.main()
