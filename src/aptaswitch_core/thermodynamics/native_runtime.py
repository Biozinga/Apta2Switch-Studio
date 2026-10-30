"""Private ctypes boundary for the explicit, separately built hybrid library.

No discovery, installation, Python energy callback or folding-library import.
Each request owns its native handle and parameter snapshots. The higher-level
backend checks build provenance before this module loads the supplied binary.
"""
from __future__ import annotations

import ctypes as C
import math
from pathlib import Path
import threading

from .structures import format_structure, is_connected, parse_structure
from .types import EmptyEnsembleError


# Loader metadata only, never model parameters or scientific result caches.
# dlopen can reuse the old image when its file has been replaced in place.
_LOADED_IDENTITIES: dict[Path, str] = {}
_LOAD_LOCK = threading.Lock()


class _Input(C.Structure):
    _fields_ = [("struct_size", C.c_uint32), ("length", C.c_uint32),
                ("sequence", C.c_char_p), ("materials", C.POINTER(C.c_uint8)),
                ("nicks_before", C.POINTER(C.c_uint8)), ("pair_allowed", C.POINTER(C.c_uint8)),
                ("constraints", C.POINTER(C.c_int32)), ("min_hairpin_size", C.c_int32),
                ("max_interior_size", C.c_int32), ("memory_limit_bytes", C.c_uint64)]


class _HybridModel(C.Structure):
    _fields_ = [("struct_size", C.c_uint32), ("salt_policy", C.c_uint32),
                ("celsius", C.c_double), ("sodium_molar", C.c_double), ("magnesium_molar", C.c_double)]


def _call(function, *arguments):
    error = C.create_string_buffer(2048)
    status = function(*arguments, error, len(error))
    if status:
        exception = {1: ValueError, 2: ArithmeticError, 3: EmptyEnsembleError,
                     5: MemoryError}.get(status, RuntimeError)
        raise exception(error.value.decode("utf-8", errors="replace") or f"Native error {status}")


class NativeRuntime:
    """Configured ABI entrypoints; request handles are never shared or cached."""

    def __init__(self, path: Path, *, expected_sha256: str):
        with _LOAD_LOCK:
            previous = _LOADED_IDENTITIES.get(path)
            if previous is not None and previous != expected_sha256:
                raise RuntimeError("Native binary changed at an already loaded path; restart the process or use a distinct build path")
            self.library = C.CDLL(str(path))
            _LOADED_IDENTITIES[path] = expected_sha256
        lib = self.library
        lib.aptadp_abi_version.argtypes, lib.aptadp_abi_version.restype = [], C.c_uint32
        if lib.aptadp_abi_version() != 1:
            raise RuntimeError("Unsupported connected-DP ABI; expected version 1")
        signatures = {
            "aptadp_create_hybrid_v1": [C.POINTER(_Input), C.POINTER(_HybridModel),
                                       C.POINTER(C.c_void_p), C.POINTER(C.c_double)],
            "aptadp_create_rna_v1": [C.POINTER(_Input), C.POINTER(_HybridModel),
                                    C.POINTER(C.c_void_p), C.POINTER(C.c_double)],
            "aptadp_log_partition": [C.c_void_p, C.POINTER(C.c_double)],
            "aptadp_pair_probabilities": [C.c_void_p, C.POINTER(C.c_double), C.c_size_t,
                                           C.POINTER(C.c_double)],
            "aptadp_minimum": [C.c_void_p, C.POINTER(C.c_double), C.c_char_p, C.c_size_t],
            "aptadp_minimum_pair_energies": [C.c_void_p, C.POINTER(C.c_double), C.c_size_t,
                                             C.POINTER(C.c_double)],
            "aptadp_sample": [C.c_void_p, C.c_uint64, C.c_uint64, C.c_char_p, C.c_size_t],
        }
        for name, arguments in signatures.items():
            try:
                function = getattr(lib, name)
            except AttributeError as error:
                raise RuntimeError(f"Native hybrid library lacks required entrypoint {name}") from error
            function.argtypes, function.restype = arguments + [C.c_char_p, C.c_size_t], C.c_int
        lib.aptadp_destroy.argtypes, lib.aptadp_destroy.restype = [C.c_void_p], None
        # Additive ABI capability. Older builds can still perform their other
        # observables; requesting alternatives fails explicitly if absent.
        self._minimum_structures = getattr(lib, "aptadp_minimum_structures", None)
        if self._minimum_structures is not None:
            self._minimum_structures.argtypes = [C.c_void_p, C.c_double, C.c_uint64,
                C.c_char_p, C.c_size_t, C.POINTER(C.c_double), C.c_size_t,
                C.POINTER(C.c_uint64), C.POINTER(C.c_uint8), C.POINTER(C.c_double),
                C.c_char_p, C.c_size_t]
            self._minimum_structures.restype = C.c_int

    def request(self, spec, conditions, constraints, salt_policy, memory_limit_bytes):
        return _Request(self, spec, conditions, constraints, salt_policy, memory_limit_bytes)


class _Request:
    def __init__(self, runtime, spec, conditions, constraints, salt_policy, memory_limit_bytes):
        self.runtime, self.spec, self.constraints = runtime, spec, constraints
        self.n, self.lengths = spec.n_nucleotides, tuple(len(s.sequence) for s in spec.strands)
        single_rna = len(spec.strands) == 1 and spec.strands[0].material == "rna"
        self.solver_name = "ViennaRNA-RNAlib-2.7.2" if single_rna else "AptaSwitch-connected-DP"
        factory = (runtime.library.aptadp_create_rna_v1 if single_rna
                   else runtime.library.aptadp_create_hybrid_v1)
        self.handle = C.c_void_p()
        chemistry = [1 if strand.material == "dna" else 2 for strand in spec.strands for _ in strand.sequence]
        materials = (C.c_uint8 * self.n)(*chemistry)
        cutpoints = set(spec.cutpoints)
        cuts = (C.c_uint8 * (self.n + 1))(*(int(i in cutpoints) for i in range(self.n + 1)))
        mask = (C.c_uint8 * (self.n * self.n))()
        # T/U share a lookup code, but input sequence and chemistry are retained
        # verbatim in the C request. This never transcribes a requested strand.
        encoded = spec.sequence.replace("U", "T")
        for i in range(self.n):
            for j in range(i + 1, self.n):
                pair = encoded[i] + encoded[j]
                allowed = pair in ("AT", "TA", "CG", "GC") or (
                    chemistry[i] == chemistry[j] == 2 and pair in ("GT", "TG"))
                mask[i * self.n + j] = mask[j * self.n + i] = allowed
        required = (C.c_int32 * self.n)(*(constraints.partners if constraints else (-1,) * self.n))
        request = _Input(C.sizeof(_Input), self.n, spec.sequence.encode("ascii"), materials,
                         cuts, mask, required, 3, 30, memory_limit_bytes)
        model = _HybridModel(C.sizeof(_HybridModel), salt_policy, conditions.celsius,
                             conditions.sodium, conditions.magnesium)
        thermal = C.c_double()
        try:
            _call(factory, C.byref(request), C.byref(model),
                  C.byref(self.handle), C.byref(thermal))
            if not self.handle.value or not math.isfinite(thermal.value) or thermal.value <= 0:
                raise ArithmeticError("Native constructor returned an invalid handle or thermal energy")
            self.thermal_energy_kcal_mol = thermal.value
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.handle.value:
            self.runtime.library.aptadp_destroy(self.handle)
            self.handle = C.c_void_p()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def _structure(self, raw):
        state = parse_structure(raw.decode("ascii"), self.lengths)
        if not is_connected(state, self.lengths) or (self.constraints is not None and any(
                wanted >= 0 and state.partners[i] != wanted
                for i, wanted in enumerate(self.constraints.partners))):
            raise ArithmeticError("Native traceback violates connectivity or the supplied constraints")
        return format_structure(state, self.lengths)

    def log_partition(self):
        result = C.c_double()
        _call(self.runtime.library.aptadp_log_partition, self.handle, C.byref(result))
        return result.value

    def pairs(self):
        matrix, adjustment = (C.c_double * (self.n * self.n))(), C.c_double()
        _call(self.runtime.library.aptadp_pair_probabilities, self.handle,
              matrix, len(matrix), C.byref(adjustment))
        if not math.isfinite(adjustment.value) or not 0 <= adjustment.value <= 1e-12:
            raise ArithmeticError("Native probability endpoint adjustment exceeds its numerical contract")
        return tuple(tuple(matrix[i * self.n + j] for j in range(self.n)) for i in range(self.n))

    def minimum(self):
        energy, structure = C.c_double(), C.create_string_buffer(self.n + 1)
        _call(self.runtime.library.aptadp_minimum, self.handle, C.byref(energy), structure, len(structure))
        return energy.value, self._structure(structure.value)

    def minimum_pair_energies(self):
        matrix, minimum = (C.c_double * (self.n * self.n))(), C.c_double()
        _call(self.runtime.library.aptadp_minimum_pair_energies, self.handle,
              matrix, len(matrix), C.byref(minimum))
        return minimum.value, tuple(tuple(matrix[i * self.n + j] for j in range(self.n))
                                    for i in range(self.n))

    def minimum_structures(self, tolerance, maximum_count):
        function = self.runtime._minimum_structures
        if function is None:
            raise RuntimeError("Native hybrid build lacks aptadp_minimum_structures; an updated explicit build is required")
        output = C.create_string_buffer(maximum_count * (self.n + 1))
        energies = (C.c_double * maximum_count)()
        count, truncated, minimum = C.c_uint64(), C.c_uint8(), C.c_double()
        _call(function, self.handle, tolerance, maximum_count, output, len(output), energies, len(energies),
              C.byref(count), C.byref(truncated), C.byref(minimum))
        if not 1 <= count.value <= maximum_count or truncated.value not in (0, 1):
            raise ArithmeticError("Native minimum-structure enumeration returned invalid counts")
        raw = output.raw
        structures = tuple(self._structure(raw[k * (self.n + 1):k * (self.n + 1) + self.n])
                           for k in range(count.value))
        return minimum.value, structures, tuple(energies[:count.value]), bool(truncated.value)

    def sample(self, count, seed):
        output = C.create_string_buffer(count * (self.n + 1))
        _call(self.runtime.library.aptadp_sample, self.handle, seed, count, output, len(output))
        raw = output.raw
        return tuple(self._structure(raw[k * (self.n + 1):k * (self.n + 1) + self.n])
                     for k in range(count))
