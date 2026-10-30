"""Explicit optional ctypes access to synthetic native grammar fixtures only.

This is a validation tool, not a selectable application/physical backend.
Build with build_connected_native.py and pass the resulting library path.
"""
import ctypes as C
from numbers import Real
from pathlib import Path

from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, EmptyEnsembleError


class NativeInput(C.Structure):
    _fields_ = [("struct_size", C.c_uint32), ("length", C.c_uint32),
                ("sequence", C.c_char_p), ("materials", C.POINTER(C.c_uint8)),
                ("nicks_before", C.POINTER(C.c_uint8)), ("pair_allowed", C.POINTER(C.c_uint8)),
                ("constraints", C.POINTER(C.c_int32)), ("min_hairpin_size", C.c_int32),
                ("max_interior_size", C.c_int32), ("memory_limit_bytes", C.c_uint64)]


class NativeModel(C.Structure):
    _fields_ = [("struct_size", C.c_uint32), ("fixture_provider", C.c_uint32),
                ("thermal_energy_kcal_mol", C.c_double), ("log_global_factor", C.c_double)]


class FixtureNativeGrammar:
    """Own one native handle; call close(), or use a context manager.

    Fixture names must be explicit. No library search, download, implicit
    physical provider selection or Python callback occurs during native DP.
    """
    def __init__(self, library, spec: ComplexSpec, *, fixture: str,
                 min_hairpin_size: int, max_interior_size: int = 30,
                 constraints: PairingConstraints | None = None, pair_allowed=None,
                 thermal_energy_kcal_mol: float = 1, log_global_factor: float = 0,
                 memory_limit_bytes: int = 0):
        if not isinstance(spec, ComplexSpec):
            raise ValueError("spec must be a ComplexSpec")
        if fixture not in ("unit", "asymmetric"):
            raise ValueError("An explicit synthetic fixture is required")
        for value, label in ((min_hairpin_size, "min_hairpin_size"),
                             (max_interior_size, "max_interior_size"), (memory_limit_bytes, "memory_limit_bytes")):
            if type(value) is not int or value < 0:
                raise ValueError(f"{label} must be a nonnegative integer")
        if min_hairpin_size >= 2**31 or max_interior_size > 30 or memory_limit_bytes >= 2**64:
            raise ValueError("Integer option exceeds its explicit ABI range")
        if any(isinstance(value, bool) or not isinstance(value, Real)
               for value in (thermal_energy_kcal_mol, log_global_factor)):
            raise ValueError("Thermal energy and global log factor must be real numbers")
        self.spec, self.n = spec, spec.n_nucleotides
        self.lengths = tuple(len(s.sequence) for s in spec.strands)
        if constraints is not None and (not isinstance(constraints, PairingConstraints)
                                       or len(constraints.partners) != self.n):
            raise ValueError("Constraints must cover the complete complex")
        self.library = C.CDLL(str(Path(library).resolve()))
        self.handle = C.c_void_p()
        self._configure()
        if self.library.aptadp_abi_version() != 1:
            raise RuntimeError("Unsupported connected-DP ABI")
        materials = (C.c_uint8 * self.n)(*(1 if s.material == "dna" else 2 for s in spec.strands for _ in s.sequence))
        cuts = (C.c_uint8 * (self.n + 1))(*(int(k in spec.cutpoints) for k in range(self.n + 1)))
        matrix = (C.c_uint8 * (self.n * self.n))()
        for i in range(self.n):
            for j in range(i + 1, self.n):
                allowed = True if pair_allowed is None else pair_allowed(spec, i, j)
                if type(allowed) is not bool:
                    raise ValueError("pair_allowed must return bool")
                matrix[i * self.n + j] = matrix[j * self.n + i] = allowed
        required = (C.c_int32 * self.n)(*(constraints.partners if constraints is not None else (-1,) * self.n))
        request = NativeInput(C.sizeof(NativeInput), self.n, spec.sequence.encode("ascii"), materials,
                              cuts, matrix, required, min_hairpin_size, max_interior_size, memory_limit_bytes)
        model = NativeModel(C.sizeof(NativeModel), 1 if fixture == "unit" else 2,
                            thermal_energy_kcal_mol, log_global_factor)
        self._call(self.library.aptadp_create_v1, C.byref(request), C.byref(model), C.byref(self.handle))

    def _configure(self):
        lib, error = self.library, [C.c_char_p, C.c_size_t]
        lib.aptadp_abi_version.argtypes, lib.aptadp_abi_version.restype = [], C.c_uint32
        lib.aptadp_create_v1.argtypes = [C.POINTER(NativeInput), C.POINTER(NativeModel), C.POINTER(C.c_void_p)] + error
        lib.aptadp_destroy.argtypes, lib.aptadp_destroy.restype = [C.c_void_p], None
        lib.aptadp_log_partition.argtypes = [C.c_void_p, C.POINTER(C.c_double)] + error
        lib.aptadp_pair_probabilities.argtypes = [C.c_void_p, C.POINTER(C.c_double), C.c_size_t, C.POINTER(C.c_double)] + error
        lib.aptadp_minimum.argtypes = [C.c_void_p, C.POINTER(C.c_double), C.c_char_p, C.c_size_t] + error
        lib.aptadp_minimum_pair_energies.argtypes = [C.c_void_p, C.POINTER(C.c_double), C.c_size_t,
                                                   C.POINTER(C.c_double)] + error
        lib.aptadp_minimum_structures.argtypes = [C.c_void_p, C.c_double, C.c_uint64, C.c_char_p,
                                                 C.c_size_t, C.POINTER(C.c_double), C.c_size_t,
                                                 C.POINTER(C.c_uint64), C.POINTER(C.c_uint8),
                                                 C.POINTER(C.c_double)] + error
        lib.aptadp_sample.argtypes = [C.c_void_p, C.c_uint64, C.c_uint64, C.c_char_p, C.c_size_t] + error
        lib.aptadp_memory_estimate.argtypes = [C.c_void_p, C.POINTER(C.c_uint64)] + error
        for name in ("aptadp_create_v1", "aptadp_log_partition", "aptadp_pair_probabilities",
                     "aptadp_minimum", "aptadp_minimum_pair_energies", "aptadp_minimum_structures",
                     "aptadp_sample", "aptadp_memory_estimate"):
            getattr(lib, name).restype = C.c_int

    @staticmethod
    def _call(function, *arguments):
        error = C.create_string_buffer(1024)
        status = function(*arguments, error, len(error))
        if status:
            exception = {1: ValueError, 2: ArithmeticError, 3: EmptyEnsembleError, 5: MemoryError}.get(status, RuntimeError)
            raise exception(error.value.decode("utf-8", errors="replace"))

    def close(self):
        if getattr(self, "handle", None) and self.handle.value:
            self.library.aptadp_destroy(self.handle)
            self.handle = C.c_void_p()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __del__(self):
        self.close()

    def log_partition(self):
        result = C.c_double()
        self._call(self.library.aptadp_log_partition, self.handle, C.byref(result))
        return result.value

    def pairs(self):
        matrix, adjustment = (C.c_double * (self.n * self.n))(), C.c_double()
        self._call(self.library.aptadp_pair_probabilities, self.handle, matrix, len(matrix), C.byref(adjustment))
        return tuple(tuple(matrix[i * self.n + j] for j in range(self.n)) for i in range(self.n)), adjustment.value

    def minimum(self):
        energy, structure = C.c_double(), C.create_string_buffer(self.n + 1)
        self._call(self.library.aptadp_minimum, self.handle, C.byref(energy), structure, len(structure))
        return energy.value, format_structure(parse_structure(structure.value.decode("ascii"), self.lengths), self.lengths)

    def minimum_pair_energies(self):
        """Return (minimum, total-energy matrix) for the labelled ensemble.

        Diagonal/impossible pairs are +inf. All finite values contain the same
        global association/uniform symmetry term. Values are not energy gaps;
        no rounding or clipping is applied to their difference from minimum.
        """
        matrix, minimum = (C.c_double * (self.n * self.n))(), C.c_double()
        self._call(self.library.aptadp_minimum_pair_energies, self.handle, matrix, len(matrix), C.byref(minimum))
        return minimum.value, tuple(tuple(matrix[i * self.n + j] for j in range(self.n)) for i in range(self.n))

    def minimum_structures(self, *, atol_kcal=0.0, count_limit=100):
        """Return (minimum, ((actual_energy, structure), ...), truncated).

        First structure matches minimum(); remaining structures are unsorted.
        The window is global and truncated proves at least one omitted result.
        """
        if type(count_limit) is not int or not 0 <= count_limit < 2**64:
            raise ValueError("count_limit must be a uint64 integer")
        maximum_bytes = C.c_size_t(-1).value
        if count_limit > min(maximum_bytes // (self.n + 1), maximum_bytes // C.sizeof(C.c_double)):
            raise ValueError("count_limit overflows the native output buffers")
        if isinstance(atol_kcal, bool) or not isinstance(atol_kcal, Real):
            raise ValueError("atol_kcal must be a finite nonnegative real")
        output = C.create_string_buffer(count_limit * (self.n + 1))
        energies = (C.c_double * count_limit)()
        count, truncated, minimum = C.c_uint64(), C.c_uint8(), C.c_double()
        self._call(self.library.aptadp_minimum_structures, self.handle, atol_kcal, count_limit,
                   output, len(output), energies, len(energies), C.byref(count), C.byref(truncated), C.byref(minimum))
        raw = output.raw
        structures = tuple((energies[k], format_structure(parse_structure(
            raw[k * (self.n + 1):k * (self.n + 1) + self.n].decode("ascii"), self.lengths), self.lengths))
                           for k in range(count.value))
        return minimum.value, structures, bool(truncated.value)

    def sample(self, count, *, seed):
        if type(count) is not int or not 0 <= count < 2**64 or type(seed) is not int or not 0 <= seed < 2**64:
            raise ValueError("count and seed must be uint64 integers")
        output = C.create_string_buffer(count * (self.n + 1))
        self._call(self.library.aptadp_sample, self.handle, seed, count, output, len(output))
        raw = output.raw
        return [format_structure(parse_structure(raw[k * (self.n + 1):k * (self.n + 1) + self.n].decode("ascii"),
                                                 self.lengths), self.lengths) for k in range(count)]

    def memory_estimate(self):
        result = C.c_uint64()
        self._call(self.library.aptadp_memory_estimate, self.handle, C.byref(result))
        return result.value
