"""Full-reporter admission without allocating or folding a giant test complex."""
from contextlib import nullcontext
import hashlib
import json
import os
from pathlib import Path
import platform
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.thermodynamics import memory_policy as memory
from aptaswitch_core.thermodynamics import native_backend
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.types import (
    ComplexSpec, Conditions, ModelProvenance, Strand, enumerate_complexes,
)

GiB = 1024**3
MiB = 1024**2


class MemoryPolicyTests(unittest.TestCase):
    def setUp(self):
        self.assertEqual(memory._RESERVED, 0)
        self.addCleanup(lambda: self.assertEqual(memory._RESERVED, 0))

    def test_exact_native_formula_admits_full_reporter_catalog_through_five_strands(self):
        self.assertEqual(memory.workspace_bytes(1594), 664772084)
        strands = (Strand('A' * 797, 'rna', 'Switch'), Strand('T' * 30, 'dna', 'Trigger'),
                   Strand('C' * 80, 'dna', 'Aptamer'))
        catalog = enumerate_complexes(strands, 5)
        self.assertEqual(len(catalog), 95)
        self.assertEqual(max(s.n_nucleotides for s in catalog), 3985)
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            for spec in catalog:
                # This enters only the arithmetic admission layer: no runtime,
                # quadratic allocation, or cubic folding is performed.
                with memory.reserve_workspace(spec.n_nucleotides, explicit_limit=None) as limit:
                    self.assertGreaterEqual(limit, memory.workspace_bytes(spec.n_nucleotides))
                    self.assertLessEqual(limit, 6 * GiB)

    def test_insufficient_host_memory_does_not_shrink_the_requested_catalog(self):
        with patch.object(memory, 'available_memory_bytes', return_value=2 * GiB):
            with self.assertRaisesRegex(MemoryError, '3985.*approximately.*catalog is unchanged'):
                with memory.reserve_workspace(3985, explicit_limit=None):
                    self.fail('The oversized request must not be admitted')

    def test_observable_workspace_counts_only_simultaneously_live_tables(self):
        n = 3955
        one_table = ((n + 2) * (n + 1) // 2) * 21 * 8
        ancillary = n * n * 9 + (n + 1) * 256
        for tables in (1, 2, 3):
            with self.subTest(tables=tables):
                self.assertEqual(memory.workspace_bytes(n, tables=tables), tables * one_table + ancillary)
        self.assertEqual(memory.workspace_bytes(n), memory.workspace_bytes(n, tables=3))
        for tables in (0, 4, -1, True, 1.0, None):
            with self.subTest(tables=tables), self.assertRaises(ValueError):
                memory.workspace_bytes(n, tables=tables)

    def test_3955_nt_catalog_partition_and_minima_fit_reported_headroom_without_changing_catalog(self):
        # Reproduce the user's 3.81 GiB worst-case estimate versus 3.73 GiB
        # available budget, without constructing a giant native DP table.
        available = int(3.73 * GiB * 4 / 3)
        budget = memory.safe_budget_bytes(available)
        self.assertAlmostEqual(budget / GiB, 3.73, places=8)
        strands = (Strand('A' * 791, 'rna', 'Switch'), Strand('T' * 30, 'dna', 'Trigger'),
                   Strand('C' * 80, 'dna', 'Aptamer'))
        catalog = enumerate_complexes(strands, 5)
        original_keys = tuple(spec.ordered_key for spec in catalog)
        self.assertEqual(len(catalog), 95)
        self.assertEqual(max(spec.n_nucleotides for spec in catalog), 3955)
        with patch.object(memory, 'available_memory_bytes', return_value=available):
            with self.assertRaisesRegex(MemoryError, '3955.*3.81.*3.73'):
                with memory.reserve_workspace(3955, explicit_limit=None):
                    self.fail('The all-observable worst-case estimate should exceed this budget')
            for spec in catalog:
                n = spec.n_nucleotides
                with memory.reserve_workspace(n, explicit_limit=None, tables=1) as native_limit:
                    self.assertEqual(native_limit, memory.workspace_bytes(n) + 8 * MiB)
                # The baseline screen also requests up to 50 alternatives.
                # Cover its native frontier as well as its two live DP tables.
                single_rna = len(spec.strands) == 1 and spec.strands[0].material == 'rna'
                output = 50 * (n + 1 + 8)
                with memory.reserve_workspace(n, explicit_limit=None, tables=3 if single_rna else 2,
                        enumerate_minima=not single_rna, output_bytes=output,
                        python_bytes=3 * output + 128 * 50) as native_limit:
                    self.assertEqual(native_limit, memory.workspace_bytes(n) + output + 8 * MiB)
                    self.assertLessEqual(memory._RESERVED, budget)
            self.assertEqual(tuple(spec.ordered_key for spec in catalog), original_keys)
            self.assertEqual(enumerate_complexes(strands, 5), catalog)

    def test_operation_sizing_preserves_native_cap_outputs_and_low_memory_rejection(self):
        n, output, python = 3955, 32 * MiB, 64 * MiB
        reservation = memory.workspace_bytes(n, tables=1) + output + 8 * MiB + python
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            with memory.reserve_workspace(n, explicit_limit=None, tables=1,
                                          output_bytes=output, python_bytes=python) as native_limit:
                self.assertEqual(memory._RESERVED, reservation)
                self.assertEqual(native_limit, memory.workspace_bytes(n) + output + 8 * MiB)
            # A smaller live workspace must not override an explicit native cap
            # that is insufficient for the constructor's unchanged ABI estimate.
            with self.assertRaisesRegex(MemoryError, 'explicit limit'):
                with memory.reserve_workspace(n, explicit_limit=2 * GiB, tables=1):
                    self.fail('The unchanged native constructor cap was bypassed')
        with patch.object(memory, 'available_memory_bytes', return_value=GiB):
            with self.assertRaisesRegex(MemoryError, '3955.*catalog is unchanged'):
                with memory.reserve_workspace(n, explicit_limit=None, tables=1):
                    self.fail('A one-table calculation still needs sufficient host memory')

    def test_minimum_enumeration_reserves_the_native_frontier_and_outputs(self):
        n, count = 3955, 50
        output = count * (n + 1 + 8)
        python = 3 * output + 128 * count
        base = memory.workspace_bytes(n, tables=2) + output + 8 * MiB + python
        with patch.object(memory, 'available_memory_bytes', return_value=int(3.73 * GiB * 4 / 3)):
            with memory.reserve_workspace(n, explicit_limit=None, tables=2, enumerate_minima=True,
                                          output_bytes=output, python_bytes=python) as native_limit:
                # The native frontier can spend the unused third table up to
                # its cap. Only its absent dense pair-output matrix is removed.
                expected = native_limit - 8 * n * n + output + 8 * MiB + python
                self.assertEqual(memory._RESERVED, max(base, expected))
                self.assertGreater(memory._RESERVED, base)
                self.assertLess(memory._RESERVED, 3.73 * GiB)
        for flag in (0, 1, None, 'yes'):
            with self.subTest(flag=flag), self.assertRaises(ValueError):
                with memory.reserve_workspace(n, explicit_limit=None, enumerate_minima=flag):
                    self.fail('Non-boolean enumeration mode was accepted')

    def test_larger_explicit_enumeration_frontier_cannot_overcommit_the_host(self):
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            with self.assertRaisesRegex(MemoryError, '3955.*catalog is unchanged'):
                with memory.reserve_workspace(3955, explicit_limit=12 * GiB,
                                              tables=2, enumerate_minima=True):
                    self.fail('Explicit native cap permits a frontier larger than the host budget')
            # The same cap cannot grow a fixed-size observable's workspace.
            with memory.reserve_workspace(3955, explicit_limit=12 * GiB, tables=2) as native_limit:
                self.assertEqual(native_limit, 12 * GiB)

    def test_lazy_minimum_enumeration_releases_a_table_without_reducing_the_native_cap(self):
        n, count = 3955, 50
        one_table = ((n + 2) * (n + 1) // 2) * 21 * 8
        output = count * (n + 1 + 8)
        python = 3 * output + 128 * count
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            with memory.reserve_workspace(n, explicit_limit=None, tables=2, enumerate_minima=True,
                                          output_bytes=output, python_bytes=python) as old_cap:
                old_reservation = memory._RESERVED
            with memory.reserve_workspace(n, explicit_limit=None, tables=1, enumerate_minima=True,
                                          lazy_minimum=True, output_bytes=output,
                                          python_bytes=python) as lazy_cap:
                expected = lazy_cap - 8 * n * n - one_table + output + 8 * MiB + python
                self.assertEqual(lazy_cap, old_cap)
                self.assertEqual(memory._RESERVED, expected)
                self.assertEqual(old_reservation - memory._RESERVED, one_table)
                self.assertAlmostEqual(memory._RESERVED / GiB, 2.481, places=3)
        # The optimized calculation also fits on a host where the old native
        # frontier estimate is too large. The full catalog still applies.
        with patch.object(memory, 'available_memory_bytes', return_value=4 * GiB):
            with self.assertRaises(MemoryError):
                with memory.reserve_workspace(n, explicit_limit=None, tables=2, enumerate_minima=True,
                                              output_bytes=output, python_bytes=python):
                    self.fail('The old eager solver was admitted with the optimized bound')
            with memory.reserve_workspace(n, explicit_limit=None, tables=1, enumerate_minima=True,
                                          lazy_minimum=True, output_bytes=output, python_bytes=python):
                self.assertLess(memory._RESERVED, 3 * GiB)

    def test_lazy_minimum_flag_is_boolean_and_does_not_discount_other_operations(self):
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            for flag in (0, 1, None, 'yes'):
                with self.subTest(flag=flag), self.assertRaises(ValueError):
                    with memory.reserve_workspace(100, explicit_limit=None, lazy_minimum=flag):
                        self.fail('Non-boolean lazy minimum mode was accepted')
            with memory.reserve_workspace(3955, explicit_limit=None, tables=1, lazy_minimum=True):
                self.assertEqual(memory._RESERVED, memory.workspace_bytes(3955, tables=1) + 8 * MiB)

    def test_lazy_minimum_explicit_frontier_still_cannot_overcommit_the_host(self):
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            with self.assertRaisesRegex(MemoryError, '3955.*catalog is unchanged'):
                with memory.reserve_workspace(3955, explicit_limit=12 * GiB, tables=1,
                                              enumerate_minima=True, lazy_minimum=True):
                    self.fail('An optimized frontier exceeded the host memory budget')

    def test_explicit_native_limit_is_respected_exactly(self):
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            with memory.reserve_workspace(1594, explicit_limit=GiB) as limit:
                self.assertEqual(limit, GiB)
            with self.assertRaisesRegex(MemoryError, 'explicit limit'):
                with memory.reserve_workspace(1594, explicit_limit=512 * MiB):
                    self.fail('Explicit limit was ignored')

    def test_concurrent_reservations_and_exceptions_release_capacity(self):
        with patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            with self.assertRaisesRegex(RuntimeError, 'deliberate'):
                with memory.reserve_workspace(3985, explicit_limit=None):
                    with self.assertRaises(MemoryError):
                        with memory.reserve_workspace(3985, explicit_limit=None):
                            self.fail('Concurrent large requests overcommitted the host')
                    raise RuntimeError('deliberate')
            with memory.reserve_workspace(3985, explicit_limit=None):
                pass

    def test_python_and_sample_outputs_are_included_before_allocating(self):
        with patch.object(memory, 'available_memory_bytes', return_value=GiB):
            for options in ({'python_bytes': GiB}, {'output_bytes': 2**63}):
                with self.subTest(options=options), self.assertRaises(MemoryError):
                    with memory.reserve_workspace(10, explicit_limit=None, **options):
                        self.fail('Output allocation was admitted')

    def test_unknown_and_low_memory_statistics_are_conservative(self):
        self.assertEqual(memory.safe_budget_bytes(None), 512 * MiB)
        self.assertEqual(memory.safe_budget_bytes(8 * GiB), 6 * GiB)
        self.assertEqual(memory.safe_budget_bytes(128 * MiB), 0)
        with patch.object(memory, 'available_memory_bytes', return_value=None):
            with memory.reserve_workspace(100, explicit_limit=None):
                pass
            with self.assertRaisesRegex(MemoryError, 'statistics'):
                with memory.reserve_workspace(1594, explicit_limit=None):
                    pass

    def test_linux_memavailable_and_cgroup_remaining_both_apply(self):
        files = {'/proc/meminfo': 'MemFree: 131072 kB\nMemAvailable: 8388608 kB\n',
                 '/proc/self/cgroup': '0::/group/leaf\n',
                 '/sys/fs/cgroup/group/leaf/memory.max': 'max',
                 '/sys/fs/cgroup/memory.max': 'max',
                 '/sys/fs/cgroup/group/memory.max': str(4 * GiB),
                 '/sys/fs/cgroup/group/memory.current': str(3 * GiB)}
        def read(path):
            if str(path) not in files:
                raise FileNotFoundError(str(path))
            return files[str(path)]
        with patch.object(Path, 'read_text', read):
            self.assertEqual(memory._linux_available(), GiB)
            files['/sys/fs/cgroup/group/memory.current'] = str(5 * GiB)
            self.assertEqual(memory._linux_available(), 0)

    def test_linux_v1_and_memfree_fallback(self):
        files = {'/proc/meminfo': 'MemFree: 8388608 kB\n', '/proc/self/cgroup': '7:memory:/group\n',
                 '/sys/fs/cgroup/memory/group/memory.limit_in_bytes': str(6 * GiB),
                 '/sys/fs/cgroup/memory/group/memory.usage_in_bytes': str(GiB)}
        def read(path):
            if str(path) not in files:
                raise FileNotFoundError(str(path))
            return files[str(path)]
        with patch.object(Path, 'read_text', read):
            self.assertEqual(memory._linux_available(), 5 * GiB)

    def test_mac_counts_exclude_speculative_and_purgeable_double_count(self):
        def statistics(host, flavor, values, count):
            self.assertEqual((host, flavor), (17, 2))
            values._obj[0], values._obj[2] = 10, 20
            values._obj[12], values._obj[14] = 9, 7
            return 0
        fake = SimpleNamespace(host_statistics=statistics)
        with patch.object(memory, '_mach_api', return_value=(fake, 17)), \
             patch.object(memory.os, 'sysconf', return_value=16384) as page_size:
            self.assertEqual(memory._mac_available(), 30 * 16384)
            page_size.assert_called_once_with('SC_PAGESIZE')

    def test_large_native_request_receives_the_admitted_limit_before_loading(self):
        backend = NativeHybridBackend.__new__(NativeHybridBackend)
        backend._memory_limit, backend._sodium_correction = None, 'published-pair'
        backend._lazy_inside = False
        request = object()
        runtime = SimpleNamespace(request=Mock(return_value=nullcontext(request)))
        with patch.object(backend, '_get_runtime', return_value=runtime) as load, \
             patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            spec = SimpleNamespace(n_nucleotides=3985, strands=(object(),) * 5)
            with backend._request(spec, 'conditions', None) as created:
                self.assertIs(created, request)
            self.assertGreater(runtime.request.call_args.args[-1], memory.workspace_bytes(3985))
            load.reset_mock()
            with patch.object(memory, 'available_memory_bytes', return_value=GiB):
                with self.assertRaises(MemoryError):
                    with backend._request(spec, 'conditions', None):
                        pass
            load.assert_not_called()

    def test_backend_routes_each_observable_to_its_simultaneous_table_count(self):
        backend = NativeHybridBackend.__new__(NativeHybridBackend)
        backend._memory_limit, backend._sodium_correction = None, 'published-pair'
        backend._lazy_inside = False
        backend._provenance = ModelProvenance('test-native', '1', 'test-model', 'sha256:test',
                                             'test-ensemble', '1 M')
        spec = ComplexSpec((Strand('GCAAA', 'rna', 'switch'), Strand('TTTGC', 'dna', 'trigger')))
        conditions, n = Conditions(37, .12, 0), spec.n_nucleotides
        structure = '.' * n
        request = SimpleNamespace(
            thermal_energy_kcal_mol=1., solver_name='test-solver',
            log_partition=lambda: 0.,
            pairs=lambda: tuple(tuple(float(i == j) for j in range(n)) for i in range(n)),
            minimum=lambda: (0., structure),
            minimum_pair_energies=lambda: (0., tuple((float('inf'),) * n for _ in range(n))),
            minimum_structures=lambda tolerance, count: (0., (structure,), (0.,), False),
            sample=lambda count, seed: (structure,) * count,
        )
        operations = (
            ('partition', lambda: backend.partition(spec, conditions), 1),
            ('analyze without pairs', lambda: backend.analyze(spec, conditions, need_pairs=False), 1),
            ('sample', lambda: backend.sample(spec, conditions, 2, seed=17), 1),
            ('analyze with pairs', lambda: backend.analyze(spec, conditions, need_pairs=True), 2),
            ('pairs', lambda: backend.pairs(spec, conditions), 2),
            ('ensemble', lambda: backend.ensemble(spec, conditions, 2, seed=17), 2),
            ('minimum', lambda: backend.minimum(spec, conditions), 2),
            ('minimum structures', lambda: backend.minimum_structures(spec, conditions), 2),
            ('minimum pair energies', lambda: backend.minimum_pair_energies(spec, conditions), 3),
        )
        for name, operation, tables in operations:
            with self.subTest(observable=name), \
                 patch.object(backend, '_request', return_value=nullcontext(request)) as reserve:
                operation()
                reserve.assert_called_once()
                self.assertEqual(reserve.call_args.kwargs.get('tables', 3), tables)
                self.assertEqual(reserve.call_args.kwargs.get('enumerate_minima', False),
                                 name == 'minimum structures')
                self.assertEqual(reserve.call_args.kwargs.get('minimum_only', False),
                                 name in ('minimum', 'minimum structures', 'minimum pair energies'))

    def test_only_capable_connected_minimum_requests_drop_the_inside_table(self):
        backend = NativeHybridBackend.__new__(NativeHybridBackend)
        backend._memory_limit, backend._sodium_correction = None, 'published-pair'
        spec = ComplexSpec((Strand('GCAAA', 'rna', 'switch'), Strand('TTTGC', 'dna', 'trigger')))
        runtime = SimpleNamespace(request=Mock(return_value=nullcontext(object())))
        with patch.object(backend, '_get_runtime', return_value=runtime), \
             patch.object(native_backend, 'reserve_workspace', return_value=nullcontext(2**30)) as reserve:
            for lazy_inside in (False, True):
                backend._lazy_inside = lazy_inside
                for minimum_only in (False, True):
                    for tables in (2, 3):
                        for enumerate_minima in (False, True):
                            with self.subTest(lazy_inside=lazy_inside, minimum_only=minimum_only,
                                              tables=tables, enumerate_minima=enumerate_minima):
                                with backend._request(spec, Conditions(37, .12, 0), None, tables=tables,
                                                      minimum_only=minimum_only,
                                                      enumerate_minima=enumerate_minima):
                                    pass
                                optimized = lazy_inside and minimum_only
                                self.assertEqual(reserve.call_args.kwargs['tables'], tables - int(optimized))
                                self.assertEqual(reserve.call_args.kwargs.get('lazy_minimum', False), optimized)
                                self.assertEqual(reserve.call_args.kwargs['enumerate_minima'], enumerate_minima)

    def test_manifest_only_enables_the_declared_lazy_inside_capability(self):
        # Older binaries remain valid, but must retain their eager allocation
        # bound. Unknown storage formats cannot authorize a lower bound.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.bin'
            path.write_bytes(b'manifest-only native fixture')
            manifest = dict(abi=1, status='built', model='experimental-Sugimoto95-Vienna04-loopmean',
                machine=platform.machine(), sources={name: 'a' * 64 for name in native_backend._REQUIRED_SOURCES},
                vienna_version='2.7.2', single_rna_solver='rnalib2-pf-bpp-sample-mfe-v1',
                single_rna_minimum_convention='turner2004-centikcal-salt-rounded-v1',
                native_library_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                generated_parameters_sha256='b' * 64, vienna_params_c_sha256='c' * 64,
                vienna_static_library_sha256='d' * 64, vienna_config_h_sha256='e' * 64)
            for capability in (None, 'eager-inside-v1', 'lazy-inside-v1', 'future-unknown'):
                candidate = dict(manifest)
                if capability is not None:
                    candidate['connected_dp_storage'] = capability
                path.with_suffix('.bin.json').write_text(json.dumps(candidate))
                with self.subTest(capability=capability), \
                     patch.object(native_backend, 'NativeRuntime') as load:
                    if capability == 'future-unknown':
                        with self.assertRaises(ValueError):
                            NativeHybridBackend(path)
                    else:
                        backend = NativeHybridBackend(path)
                        self.assertEqual(backend._lazy_inside, capability == 'lazy-inside-v1')
                    load.assert_not_called()

    def test_single_rna_keeps_conservative_admission_for_its_distinct_solver(self):
        backend = NativeHybridBackend.__new__(NativeHybridBackend)
        backend._memory_limit, backend._sodium_correction = None, 'published-pair'
        backend._lazy_inside = True
        spec = ComplexSpec((Strand('A' * 791, 'rna', 'Switch'),))
        runtime = SimpleNamespace(request=Mock(return_value=nullcontext(object())))
        with patch.object(backend, '_get_runtime', return_value=runtime), \
             patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            for requested_tables in (1, 2, 3):
                for enumerate_minima in (False, True):
                    for minimum_only in (False, True):
                        with self.subTest(requested_tables=requested_tables, enumerate_minima=enumerate_minima,
                                          minimum_only=minimum_only):
                            with backend._request(spec, Conditions(37, .12, 0), None, tables=requested_tables,
                                                  enumerate_minima=enumerate_minima, minimum_only=minimum_only):
                                self.assertEqual(memory._RESERVED, memory.workspace_bytes(791) + 8 * MiB)


@unittest.skipUnless(os.environ.get('APTASWITCH_HYBRID_NATIVE_LIBRARY'), 'Compiled native library required')
class NativeMemoryPolicyTests(unittest.TestCase):
    def test_real_full_reporter_dimer_constructor_passes_old_512_mib_ceiling(self):
        # Exercise the C constructor without running the cubic recurrence.
        # Its conservative workspace is 664,772,084 bytes. Eager builds allocate
        # the inside table here; lazy builds defer it to the first PF request.
        required = memory.workspace_bytes(1594) + 8 * MiB
        if memory.safe_budget_bytes(memory.available_memory_bytes()) < required:
            self.skipTest('Insufficient available RAM for the constructor-only check')
        backend = NativeHybridBackend(os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY'])
        spec = ComplexSpec((Strand('A' * 797, 'rna', 'Switch'),) * 2)
        with patch.object(backend, '_partition', side_effect=AssertionError('No folding in this admission test')):
            with backend._request(spec, Conditions(29, .12, 0), None) as created:
                self.assertEqual(created.n, 1594)
                self.assertTrue(created.handle.value)

    def test_automatic_and_explicit_budgets_do_not_change_scientific_outputs(self):
        path = os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY']
        automatic = NativeHybridBackend(path)
        explicit = NativeHybridBackend(path, memory_limit_bytes=512 * MiB)
        spec = ComplexSpec((Strand('GCAAA', 'rna', 'switch'), Strand('TTTGC', 'dna', 'trigger')))
        c = Conditions(37, .12, 0)
        self.assertEqual(automatic.provenance, explicit.provenance)
        self.assertEqual(automatic.partition(spec, c), explicit.partition(spec, c))
        self.assertEqual(automatic.minimum_structures(spec, c), explicit.minimum_structures(spec, c))
        self.assertEqual(automatic.ensemble(spec, c, 20, seed=17), explicit.ensemble(spec, c, 20, seed=17))

    def test_observable_reservations_preserve_all_outputs_against_worst_case_reservations(self):
        path = os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY']
        operation_sized = NativeHybridBackend(path)
        worst_case = NativeHybridBackend(path)
        original_request = worst_case._request
        def reserve_worst_case(*args, **options):
            return original_request(*args, **dict(options, tables=3,
                enumerate_minima=False, minimum_only=False))
        c = Conditions(37, .12, 0)
        specs = (
            ComplexSpec((Strand('GCAAA', 'rna', 'switch'), Strand('TTTGC', 'dna', 'trigger'))),
            ComplexSpec((Strand('GCAAAAGC', 'rna', 'switch'),)),
        )
        operations = (
            ('partition', lambda backend, spec: backend.partition(spec, c)),
            ('analyze without pairs', lambda backend, spec: backend.analyze(spec, c, need_pairs=False)),
            ('analyze with pairs', lambda backend, spec: backend.analyze(spec, c, need_pairs=True)),
            ('pairs', lambda backend, spec: backend.pairs(spec, c)),
            ('minimum', lambda backend, spec: backend.minimum(spec, c)),
            ('minimum structures', lambda backend, spec: backend.minimum_structures(spec, c)),
            ('minimum pair energies', lambda backend, spec: backend.minimum_pair_energies(spec, c)),
            ('sample', lambda backend, spec: backend.sample(spec, c, 20, seed=17)),
            ('ensemble', lambda backend, spec: backend.ensemble(spec, c, 20, seed=17)),
        )
        with patch.object(worst_case, '_request', side_effect=reserve_worst_case), \
             patch.object(memory, 'available_memory_bytes', return_value=8 * GiB):
            for spec in specs:
                for name, operation in operations:
                    with self.subTest(spec=spec.ordered_key, observable=name):
                        self.assertEqual(operation(operation_sized, spec), operation(worst_case, spec))


if __name__ == '__main__':
    unittest.main()
