"""Production boundary contracts and optional real native/oracle comparisons.

Set APTASWITCH_HYBRID_NATIVE_LIBRARY to a complete explicit build. RNA Python
is only an optional test oracle; the independent subprocess forbids it.
"""
from dataclasses import FrozenInstanceError
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.structures import PairingConstraints, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, EmptyEnsembleError, Strand


LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")
HAS_ORACLE = importlib.util.find_spec("RNA") is not None


def spec(*items):
    return ComplexSpec(tuple(Strand(sequence, material, f"s{i}")
                             for i, (sequence, material) in enumerate(items)))


class NativeBackendGuardsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "test.bin"
        self.path.write_bytes(b"not an executable; manifest checks only")
        prefix = "src/aptaswitch_core/thermodynamics/"
        names = [prefix + "native/" + name for name in (
            "connected_dp.cpp", "connected_dp.hpp", "hybrid_provider.cpp", "hybrid_provider.hpp",
            "hybrid_api.cpp", "hybrid_api.hpp")] + [prefix + "hybrid_motifs.py"]
        self.manifest = dict(abi=1, status="built", model="experimental-Sugimoto95-Vienna04-loopmean",
            machine=platform.machine(), sources={name: "a" * 64 for name in names},
            vienna_version="2.7.2", single_rna_solver="rnalib2-pf-bpp-sample-mfe-v1",
            single_rna_minimum_convention="turner2004-centikcal-salt-rounded-v1",
            native_library_sha256=hashlib.sha256(self.path.read_bytes()).hexdigest(),
            generated_parameters_sha256="b" * 64, vienna_params_c_sha256="c" * 64,
            vienna_static_library_sha256="d" * 64, vienna_config_h_sha256="e" * 64)
        self.save_manifest()

    def save_manifest(self):
        self.path.with_suffix(".bin.json").write_text(json.dumps(self.manifest))

    def test_no_load_until_valid_request_and_standard_state_is_explicit(self):
        with patch("aptaswitch_core.thermodynamics.native_backend.NativeRuntime") as runtime:
            backend = NativeHybridBackend(self.path)
            runtime.assert_not_called()
        self.assertIn("experimental", backend.provenance.model)
        self.assertIn(self.manifest["native_library_sha256"], backend.provenance.version)
        self.assertIn("RNA-monomer:ViennaRNA-2.7.2", backend.provenance.ensemble)
        self.assertIn("minimum-RNA:centikcal", backend.provenance.ensemble)
        self.assertEqual(backend.standard_state.concentration_molar, 1.0)
        self.assertEqual(backend.standard_state.provenance_label, backend.provenance.standard_state)
        self.assertTrue(backend.standard_state.association_included)
        self.assertTrue(backend.standard_state.cyclic_symmetry_included)
        self.assertTrue(backend.standard_state.connected_ensembles)

    def test_manifest_rejects_missing_wrong_or_malformed_provenance_before_load(self):
        original = json.loads(json.dumps(self.manifest))
        cases = {"status": "compiling", "model": "external", "abi": True,
                 "machine": "another-machine", "generated_parameters_sha256": "bad",
                 "sources": {}, "vienna_static_library_sha256": None,
                 "vienna_config_h_sha256": "bad", "vienna_version": "2.7.3",
                 "single_rna_solver": "connected-dp", "single_rna_minimum_convention": "pf-weights"}
        for key, bad in cases.items():
            self.manifest = {**original, key: bad}
            self.save_manifest()
            with self.subTest(field=key), patch("aptaswitch_core.thermodynamics.native_backend.NativeRuntime") as load:
                with self.assertRaises(ValueError):
                    NativeHybridBackend(self.path)
                load.assert_not_called()
        self.manifest = original
        self.save_manifest()
        self.path.write_bytes(b"changed binary")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            NativeHybridBackend(self.path)

    def test_manifest_requires_explicit_rna_solver_and_header_metadata(self):
        original = dict(self.manifest)
        for key in ("vienna_version", "single_rna_solver", "single_rna_minimum_convention",
                    "vienna_config_h_sha256"):
            self.manifest = {name: value for name, value in original.items() if name != key}
            self.save_manifest()
            with self.subTest(missing=key), patch("aptaswitch_core.thermodynamics.native_backend.NativeRuntime") as load:
                with self.assertRaises(ValueError):
                    NativeHybridBackend(self.path)
                load.assert_not_called()

    def test_binary_changed_between_construction_and_first_call_is_rejected(self):
        backend = NativeHybridBackend(self.path)
        self.path.write_bytes(b"changed")
        with patch("aptaswitch_core.thermodynamics.native_backend.NativeRuntime") as load:
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                backend.partition(spec(("A", "rna")), Conditions())
            load.assert_not_called()

    def test_provenance_changes_with_salt_policy_parameters_and_source(self):
        initial = NativeHybridBackend(self.path).provenance
        self.assertNotEqual(initial, NativeHybridBackend(self.path, sodium_correction="reference").provenance)
        for key in ("generated_parameters_sha256", "vienna_static_library_sha256", "vienna_params_c_sha256",
                    "vienna_config_h_sha256"):
            old = self.manifest[key]
            self.manifest[key] = "9" * 64
            self.save_manifest()
            self.assertNotEqual(initial.parameter_digest, NativeHybridBackend(self.path).provenance.parameter_digest)
            self.manifest[key] = old
        self.manifest["sources"][next(iter(self.manifest["sources"]))] = "f" * 64
        self.save_manifest()
        self.assertNotEqual(initial.parameter_digest, NativeHybridBackend(self.path).provenance.parameter_digest)

    def test_invalid_inputs_fail_before_loading(self):
        backend = NativeHybridBackend(self.path)
        current = spec(("GCA", "rna"), ("TGC", "dna"))
        for conditions in (Conditions(celsius=24), Conditions(celsius=56), Conditions(sodium=.119),
                           Conditions(sodium=1.01), Conditions(magnesium=.001)):
            with self.subTest(conditions=conditions), self.assertRaises(ValueError):
                backend.partition(current, conditions)
        for constraints in (True, (), PairingConstraints.unconstrained(3)):
            with self.assertRaises(ValueError):
                backend.partition(current, Conditions(), constraints)
        for options in ({"count": True, "seed": 0}, {"count": 1, "seed": -1},
                        {"count": 1, "seed": 2**64}):
            with self.assertRaises(ValueError):
                backend.sample(current, Conditions(), **options)
        for bounds in ((True, 2), (-1, 3), (0, 0), (0, 7)):
            with self.assertRaises(ValueError):
                backend.screen_unpaired_region(current, Conditions(), start=bounds[0], end=bounds[1])
        for field in ("absolute_tolerance", "relative_tolerance"):
            for value in (True, -1, math.nan, math.inf):
                with self.assertRaises(ValueError):
                    backend.screen_unpaired_region(current, Conditions(), start=0, end=2, **{field: value})
        self.assertIsNone(backend._runtime)
        with self.assertRaises(MemoryError):
            backend.sample(current, Conditions(), 2**63, seed=0)
        for value in (0, -1, True, 2**64, 1.5):
            with self.assertRaises(ValueError):
                NativeHybridBackend(self.path, memory_limit_bytes=value)
        with self.assertRaises(ValueError):
            NativeHybridBackend(self.path, sodium_correction="guess")
        for value in (0, True, -1, 2**64):
            with self.assertRaises(ValueError):
                backend.minimum_structures(current, Conditions(), maximum_count=value)
        with self.assertRaises(MemoryError):
            backend.minimum_structures(current, Conditions(), maximum_count=2**63)
        for value in (True, -1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                backend.minimum_structures(current, Conditions(), absolute_tolerance=value)

    def test_runtime_rejects_wrong_actual_abi_and_missing_entrypoints(self):
        from aptaswitch_core.thermodynamics.native_runtime import NativeRuntime
        for version, message in ((2, "ABI"), (1, "entrypoint")):
            library = SimpleNamespace(aptadp_abi_version=lambda: version)
            with patch("aptaswitch_core.thermodynamics.native_runtime.C.CDLL", return_value=library):
                with self.assertRaisesRegex(RuntimeError, message):
                    NativeRuntime(self.path, expected_sha256=self.manifest["native_library_sha256"])
        with patch("aptaswitch_core.thermodynamics.native_runtime.C.CDLL") as loader:
            with self.assertRaisesRegex(RuntimeError, "already loaded"):
                NativeRuntime(self.path, expected_sha256="0" * 64)
            loader.assert_not_called()


@unittest.skipUnless(LIBRARY, "explicit native hybrid library required")
class NativeBackendTests(unittest.TestCase):
    def setUp(self):
        self.backend = NativeHybridBackend(LIBRARY)

    def test_real_analysis_preserves_context_chemistry_and_immutability(self):
        current, c = spec(("GCAAA", "rna"), ("TTTGC", "dna")), Conditions()
        actual = self.backend.analyze(current, c, need_pairs=True)
        self.assertEqual(actual.partition.complex_spec, current)
        self.assertEqual(actual.partition.conditions, c)
        self.assertEqual(actual.partition, actual.pairs.partition)
        self.assertEqual(self.backend.partition(current, c), actual.partition)
        self.assertIsNone(self.backend.analyze(current, c, need_pairs=False).pairs)
        minimum = self.backend.minimum(current, c)
        self.assertEqual(minimum.context.complex_spec.ordered_key, (("rna", "GCAAA"), ("dna", "TTTGC")))
        self.assertEqual(minimum.context.provenance, actual.partition.provenance)
        self.assertEqual(minimum.context.solver, "AptaSwitch-connected-DP")
        self.assertEqual(minimum.context.thermal_energy_kcal_mol, actual.partition.thermal_energy_kcal_mol)
        with self.assertRaises(FrozenInstanceError):
            minimum.structure = "bogus"
        constrained = PairingConstraints.unconstrained(10).assign(0, 9)
        result = self.backend.analyze(current, c, need_pairs=True, constraints=constrained)
        self.assertEqual(result.partition.constraints, constrained)
        self.assertAlmostEqual(result.pairs.matrix[0][9], 1., places=12)
        self.assertIsNone(self.backend.partition(current, c, PairingConstraints.unconstrained(10)).constraints)

    def test_empty_ensemble_has_no_minimum_pairs_samples_or_region_screen(self):
        current, c = spec(("A", "rna"), ("A", "dna")), Conditions()
        actual = self.backend.analyze(current, c, need_pairs=True)
        self.assertTrue(actual.partition.is_empty)
        self.assertIsNone(actual.pairs)
        operations = (lambda: self.backend.pairs(current, c), lambda: self.backend.minimum(current, c),
                      lambda: self.backend.sample(current, c, 1, seed=0),
                      lambda: self.backend.minimum_pair_energies(current, c),
                      lambda: self.backend.screen_unpaired_region(current, c, start=0, end=1))
        for operation in operations:
            with self.assertRaises(EmptyEnsembleError):
                operation()

    def test_wobble_is_rna_only_and_memory_failures_are_explicit(self):
        c = Conditions()
        self.assertFalse(self.backend.partition(spec(("GG", "rna"), ("UU", "rna")), c).is_empty)
        self.assertTrue(self.backend.partition(spec(("GG", "dna"), ("TT", "dna")), c).is_empty)
        self.assertTrue(self.backend.partition(spec(("GG", "rna"), ("TT", "dna")), c).is_empty)
        small = NativeHybridBackend(LIBRARY, memory_limit_bytes=100)
        with self.assertRaises(MemoryError):
            small.partition(spec(("GCGAAAGC", "rna")), c)

    def test_bounded_minimum_alternatives_retain_actual_energies_and_truncation(self):
        current, c = spec(("AA", "rna"), ("T", "dna")), Conditions()
        first = self.backend.minimum(current, c)
        one = self.backend.minimum_structures(current, c, maximum_count=1)
        full = self.backend.minimum_structures(current, c, maximum_count=8)
        self.assertEqual(one.structures, (first.structure,))
        self.assertTrue(one.truncated)
        self.assertEqual(len(full.structures), 2)
        self.assertFalse(full.truncated)
        self.assertEqual(full.structures[0], first.structure)
        self.assertEqual(full.minimum_energy_kcal_mol, first.energy_kcal_mol)
        for energy, structure in zip(full.energies_kcal_mol, full.structures):
            constrained = parse_structure(structure, (2, 1))
            self.assertAlmostEqual(self.backend.minimum(current, c, constrained).energy_kcal_mol, energy, places=12)
        restricted = parse_structure(first.structure, (2, 1))
        only = self.backend.minimum_structures(current, c, restricted, maximum_count=1)
        self.assertEqual(only.structures, (first.structure,))
        self.assertFalse(only.truncated)
        empty = spec(("A", "rna"), ("A", "dna"))
        with self.assertRaises(EmptyEnsembleError):
            self.backend.minimum_structures(empty, c)

    def test_joint_samples_are_reproducible_and_match_pair_marginals(self):
        current, c = spec(("GCAAA", "rna"), ("TTTGC", "dna")), Conditions()
        first = self.backend.sample(current, c, 1000, seed=20261005)
        second = self.backend.sample(current, c, 1000, seed=20261005)
        self.assertEqual(first, second)
        self.assertIsInstance(first.structures, tuple)
        expected = self.backend.pairs(current, c).matrix
        counts = [[0] * 10 for _ in range(10)]
        for structure in first.structures:
            state = parse_structure(structure, (5, 5))
            for i, j in enumerate(state.partners):
                counts[i][j] += 1
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(counts[i][j] / 1000, expected[i][j], delta=.075)

    def test_region_screen_finds_degenerate_failure_missing_from_one_traceback(self):
        current, c = spec(("AA", "rna"), ("T", "dna")), Conditions()
        chosen = self.backend.minimum(current, c)
        state = parse_structure(chosen.structure, (2, 1))
        free = next(i for i in range(2) if state.partners[i] == i)
        result = self.backend.screen_unpaired_region(current, c, start=free, end=free + 1)
        self.assertFalse(result.all_minima_unpaired)
        self.assertAlmostEqual(result.gap_kcal_mol, 0, places=12)
        self.assertIsNotNone(result.witness)
        alternative = parse_structure(result.witness.structure, (2, 1))
        self.assertEqual(alternative.partners[free], 2)
        self.assertNotEqual(result.witness.structure, chosen.structure)
        self.assertEqual(result.witness.context.constraints.partners[free], 2)
        linear = self.backend.screen_unpaired_region(spec(("AAA", "rna")), c, start=0, end=3)
        self.assertTrue(linear.all_minima_unpaired)
        self.assertEqual(linear.gap_kcal_mol, math.inf)
        self.assertIsNone(linear.witness)

    @unittest.skipUnless(HAS_ORACLE, "optional Python RNA oracle unavailable")
    def test_real_native_matches_independent_exhaustive_energy_oracle(self):
        from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend
        from aptaswitch_core.thermodynamics.hybrid_motifs import HybridMotifModel
        model = HybridMotifModel(sodium_correction="published-pair")
        oracle = ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
            thermal_energy=model.thermal_energy, provenance=model.provenance, min_hairpin_size=3)
        cases = [spec(("GCGAAAGC", "rna")), spec(("GCGAAAGC", "dna")),
                 spec(("GCAAA", "rna"), ("TTTGC", "dna")),
                 spec(("TTTGC", "dna"), ("GCAAA", "rna"))]
        cases.extend(spec(*[("GC", "rna" if i % 2 else "dna") for i in range(n)]) for n in (3, 4, 5))
        for current in cases:
            for c in (Conditions(celsius=25), Conditions(celsius=55, sodium=1)):
                with self.subTest(order=current.ordered_key, conditions=c):
                    expected = oracle.evaluate(current, c)
                    actual = self.backend.analyze(current, c, need_pairs=True)
                    self.assertAlmostEqual(actual.partition.log_partition, expected.partition.log_partition, places=10)
                    for row, wanted in zip(actual.pairs.matrix, expected.pairs().matrix):
                        for value, target in zip(row, wanted):
                            self.assertAlmostEqual(value, target, places=10)
                    # RNAlib's single-RNA MFE uses the documented centikcal
                    # convention, unlike this oracle's PF-derived energies.
                    # PF and BPP above still share the same unrounded weights.
                    if len(current.strands) == 1 and current.strands[0].material == "rna":
                        continue
                    minima = self.backend.minimum_pair_energies(current, c)
                    self.assertAlmostEqual(minima.minimum_energy_kcal_mol, expected.mfe()[0], places=10)
                    for i in range(current.n_nucleotides):
                        for j in range(i + 1, current.n_nucleotides):
                            target = min((energy for state, energy in zip(expected.structures, expected.energies_kcal_mol)
                                          if state.partners[i] == j), default=math.inf)
                            if math.isfinite(target):
                                self.assertAlmostEqual(minima.matrix[i][j], target, places=10)
                            else:
                                self.assertEqual(minima.matrix[i][j], math.inf)

    def test_clean_process_runs_without_python_rna_or_validation_tools(self):
        root = Path(__file__).resolve().parents[1]
        program = r'''
import importlib.abc, json, sys
sys.path.insert(0, sys.argv[1])
class RejectScientificImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('RNA', 'tests'):
            raise AssertionError('Forbidden runtime import: ' + fullname)
sys.meta_path.insert(0, RejectScientificImports())
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.types import Strand, ComplexSpec, Conditions
b = NativeHybridBackend(sys.argv[2])
s = ComplexSpec((Strand('GCAAA','rna','r'), Strand('TTTGC','dna','d')))
a = b.analyze(s, Conditions(), need_pairs=True)
m = b.minimum(s, Conditions())
samples = b.sample(s, Conditions(), 3, seed=42)
screen = b.screen_unpaired_region(s, Conditions(), start=0, end=1)
alternatives = b.minimum_structures(s, Conditions(), maximum_count=4)
off = ComplexSpec((Strand('GCGAAAGC','rna','off'),))
off_ensemble = b.ensemble(off, Conditions(), 3, seed=42)
off_minimum = b.minimum(off, Conditions())
assert off_ensemble.analysis.partition.provenance == a.partition.provenance
print(json.dumps({'log_z':a.partition.log_partition, 'minimum':m.structure,
                  'samples':samples.structures, 'screen':screen.all_minima_unpaired,
                  'alternatives':alternatives.structures,
                  'off_log_z':off_ensemble.analysis.partition.log_partition,
                  'off_solver':off_minimum.context.solver,
                  'off_samples':off_ensemble.samples.structures}))
'''
        completed = subprocess.run([sys.executable, "-I", "-c", program, str(root / "src"), str(Path(LIBRARY).resolve())],
                                   capture_output=True, text=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(completed.stdout)
        self.assertTrue(math.isfinite(result["log_z"]))
        self.assertEqual(len(result["samples"]), 3)
        self.assertTrue(math.isfinite(result["off_log_z"]))
        self.assertEqual(result["off_solver"], "ViennaRNA-RNAlib-2.7.2")
        self.assertEqual(len(result["off_samples"]), 3)


if __name__ == "__main__":
    unittest.main()
