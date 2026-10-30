"""The compiled RNAlib route is explicit and restricted to one RNA strand."""
import ctypes as C
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
from threading import Barrier
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.thermodynamics import native_runtime
from aptaswitch_core.thermodynamics.evaluation import ComplexAnalysisCache
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.structures import PairingConstraints, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, Strand


LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")
RNA_SOLVER = "ViennaRNA-RNAlib-2.7.2"
HYBRID_SOLVER = "AptaSwitch-connected-DP"


def spec(*items):
    return ComplexSpec(tuple(Strand(sequence, material, f"s{i}")
                             for i, (sequence, material) in enumerate(items)))


class ViennaOffRoutingTests(unittest.TestCase):
    def setUp(self):
        self.path = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve() / "mock-native.bin"
        self.calls = []
        self.library = SimpleNamespace(aptadp_abi_version=Mock(return_value=1))
        for name in ("aptadp_log_partition", "aptadp_pair_probabilities", "aptadp_minimum",
                     "aptadp_minimum_pair_energies", "aptadp_sample", "aptadp_destroy"):
            setattr(self.library, name, Mock(return_value=0))
        self.library.aptadp_create_rna_v1 = Mock(side_effect=self.factory("rna"))
        self.library.aptadp_create_hybrid_v1 = Mock(side_effect=self.factory("hybrid"))
        self.enterContext(patch.object(native_runtime.C, "CDLL", return_value=self.library))
        self.enterContext(patch.dict(native_runtime._LOADED_IDENTITIES, clear=True))

    def factory(self, name):
        def create(input_ptr, model_ptr, output_ptr, thermal_ptr, error, capacity):
            request = C.cast(input_ptr, C.POINTER(native_runtime._Input)).contents
            model = C.cast(model_ptr, C.POINTER(native_runtime._HybridModel)).contents
            n = request.length
            self.calls.append(dict(factory=name, sequence=request.sequence.decode("ascii"),
                materials=tuple(request.materials[:n]), nicks=tuple(request.nicks_before[:n + 1]),
                pairs=tuple(request.pair_allowed[:n * n]), constraints=tuple(request.constraints[:n]),
                salt_policy=model.salt_policy, celsius=model.celsius, sodium=model.sodium_molar,
                magnesium=model.magnesium_molar, limit=request.memory_limit_bytes))
            C.cast(output_ptr, C.POINTER(C.c_void_p)).contents.value = 17
            C.cast(thermal_ptr, C.POINTER(C.c_double)).contents.value = .6163207755
            return 0
        return create

    def runtime(self):
        return native_runtime.NativeRuntime(self.path, expected_sha256="a" * 64)

    def test_only_single_rna_uses_official_factory_and_preserves_input_contract(self):
        runtime = self.runtime()
        conditions = Conditions(celsius=31, sodium=.27)
        cases = ((spec(("GU", "rna")), "rna", RNA_SOLVER),
                 (spec(("GT", "dna")), "hybrid", HYBRID_SOLVER),
                 (spec(("G", "rna"), ("U", "rna")), "hybrid", HYBRID_SOLVER),
                 (spec(("G", "rna"), ("T", "dna")), "hybrid", HYBRID_SOLVER))
        for current, factory, solver in cases:
            with self.subTest(current=current.ordered_key):
                with runtime.request(current, conditions, None, 2, 123456) as request:
                    self.assertEqual(request.solver_name, solver)
                    self.assertEqual(request.thermal_energy_kcal_mol, .6163207755)
                    self.assertEqual(request.handle.value, 17)
                actual = self.calls[-1]
                self.assertEqual(actual["factory"], factory)
                self.assertEqual(actual["sequence"], current.sequence)
                self.assertEqual(actual["materials"], tuple(2 if s.material == "rna" else 1
                    for s in current.strands for _ in s.sequence))
                self.assertEqual(actual["nicks"], tuple(int(i in current.cutpoints) for i in range(3)))
                self.assertEqual(actual["constraints"], (-1, -1))
                self.assertEqual((actual["salt_policy"], actual["celsius"], actual["sodium"],
                                  actual["magnesium"], actual["limit"]), (2, 31, .27, 0, 123456))
                expected_wobble = int(all(s.material == "rna" for s in current.strands))
                self.assertEqual(actual["pairs"], (0, expected_wobble, expected_wobble, 0))
        self.assertEqual(self.library.aptadp_destroy.call_count, len(cases))

    def test_constrained_rna_keeps_official_route_and_reference_salt_policy(self):
        current = spec(("GCGAAAGC", "rna"))
        constraints = PairingConstraints.unconstrained(8).assign(0, 7)
        with self.runtime().request(current, Conditions(sodium=1), constraints, 1, 123456):
            self.assertEqual(self.calls[-1]["factory"], "rna")
            self.assertEqual(self.calls[-1]["constraints"], constraints.partners)
            self.assertEqual(self.calls[-1]["salt_policy"], 1)
        self.library.aptadp_create_hybrid_v1.assert_not_called()

    def test_missing_rna_factory_is_rejected_before_any_request(self):
        del self.library.aptadp_create_rna_v1
        with self.assertRaisesRegex(RuntimeError, "aptadp_create_rna_v1"):
            self.runtime()
        self.library.aptadp_create_hybrid_v1.assert_not_called()

    def test_rna_factory_failure_does_not_fall_back_and_closes_partial_handle(self):
        def fail(input_ptr, model_ptr, output_ptr, thermal_ptr, error, capacity):
            C.cast(output_ptr, C.POINTER(C.c_void_p)).contents.value = 23
            error.value = b"RNA constructor rejected request"
            return 1
        self.library.aptadp_create_rna_v1.side_effect = fail
        with self.assertRaisesRegex(ValueError, "RNA constructor rejected"):
            self.runtime().request(spec(("AAA", "rna")), Conditions(), None, 2, 123456)
        self.library.aptadp_create_hybrid_v1.assert_not_called()
        self.library.aptadp_destroy.assert_called_once()
        self.assertEqual(self.library.aptadp_destroy.call_args.args[0].value, 23)


@unittest.skipUnless(LIBRARY, "explicit native hybrid library required")
class ViennaOffNativeTests(unittest.TestCase):
    def setUp(self):
        self.backend = NativeHybridBackend(LIBRARY)
        self.off = spec(("GCGAAAGC", "rna"))
        self.on = spec(("GCAAA", "rna"), ("TTTGC", "dna"))
        self.conditions = Conditions()

    def test_solver_contexts_share_the_combined_model_identity_and_cache(self):
        cache = ComplexAnalysisCache(self.backend)
        for current, solver in ((self.off, RNA_SOLVER), (self.on, HYBRID_SOLVER)):
            with self.subTest(current=current.ordered_key):
                analysis = cache(current, self.conditions, need_pairs=True)
                self.assertIs(cache(current, self.conditions, need_pairs=True), analysis)
                minimum = self.backend.minimum(current, self.conditions)
                samples = self.backend.sample(current, self.conditions, 3, seed=17)
                for context in (minimum.context, samples.context):
                    self.assertEqual(context.solver, solver)
                    self.assertEqual(context.provenance, analysis.partition.provenance)
                    self.assertEqual(context.provenance, self.backend.provenance)
        self.assertEqual(cache.calls, 2)
        self.assertEqual(cache.hits, 2)

    def test_off_minimum_uses_centikcal_convention_and_samples_preserve_constraints(self):
        constraints = PairingConstraints.unconstrained(8).assign(0, 7)
        minimum = self.backend.minimum(self.off, self.conditions, constraints)
        self.assertAlmostEqual(minimum.energy_kcal_mol * 100,
                               round(minimum.energy_kcal_mol * 100), delta=1e-4)
        self.assertEqual(minimum.context.solver, RNA_SOLVER)
        first = self.backend.sample(self.off, self.conditions, 50, seed=17, constraints=constraints)
        self.assertEqual(first, self.backend.sample(self.off, self.conditions, 50, seed=17, constraints=constraints))
        for structure in (minimum.structure, *first.structures):
            self.assertEqual(parse_structure(structure, (8,)).partners[0], 7)

    def test_concurrent_mixed_factories_preserve_parameters_and_seeded_rna_streams(self):
        rna = spec(("GGGAAACCCGGGAAACCC", "rna"))
        dna = spec(("GGGTTTCCCGGGTTTCCC", "dna"))
        hybrid = spec(("GCAAA", "rna"), ("TTTGC", "dna"))
        cases = (rna, rna, dna, hybrid)
        conditions = (Conditions(celsius=25, sodium=.12), Conditions(celsius=55, sodium=1))
        seeds = ((1 << 40) + 17, (1 << 32) + 29, 37, 41)

        def calculate(worker, phase):
            current, condition = cases[worker], conditions[phase]
            ensemble = self.backend.ensemble(current, condition, 64, seed=seeds[worker])
            minimum = self.backend.minimum(current, condition)
            return (ensemble.analysis.partition.log_partition, ensemble.analysis.pairs.matrix,
                    ensemble.samples.structures, minimum.energy_kcal_mol, minimum.structure)

        # Separate sequential references catch parameter-table contamination;
        # two simultaneous RNA streams catch a seed/reset race independently.
        expected = {(worker, phase): calculate(worker, phase)
                    for worker in range(4) for phase in range(2)}
        self.assertTrue(any(len(set(expected[worker, phase][2])) > 1
                            for worker in (0, 1) for phase in range(2)))
        start = Barrier(4, timeout=30)

        def concurrent(worker):
            for iteration in range(12):
                phase = (worker + iteration) % 2
                start.wait()
                self.assertEqual(calculate(worker, phase), expected[worker, phase])

        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(concurrent, worker) for worker in range(4)]
            for future in futures:
                future.result(timeout=60)


if __name__ == "__main__":
    unittest.main()
