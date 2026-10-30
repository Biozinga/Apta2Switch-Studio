"""Explicit native hybrid handle for isolated verification, not app integration."""
import ctypes as C
from pathlib import Path

from aptaswitch_core.thermodynamics.structures import PairingConstraints
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions
from tests.native_reference.connected_native import FixtureNativeGrammar, NativeInput


class NativeHybridModel(C.Structure):
    _fields_ = [("struct_size", C.c_uint32), ("salt_policy", C.c_uint32),
                ("celsius", C.c_double), ("sodium_molar", C.c_double), ("magnesium_molar", C.c_double)]


class HybridNativeGrammar(FixtureNativeGrammar):
    """The same handle observables as the fixture, with the physical provider.

    The model remains experimental. Library selection is explicit, and this
    diagnostic wrapper neither discovers binaries nor supplies app provenance.
    """
    def __init__(self, library, spec: ComplexSpec, conditions: Conditions, *,
                 sodium_correction="reference", constraints: PairingConstraints | None = None,
                 memory_limit_bytes=0, factory="hybrid"):
        if not isinstance(spec, ComplexSpec) or not isinstance(conditions, Conditions):
            raise ValueError("Explicit spec and conditions are required")
        if sodium_correction not in ("reference", "published-pair"):
            raise ValueError("Unsupported salt policy")
        if factory not in ("hybrid", "rna"):
            raise ValueError("Unsupported native factory")
        if type(memory_limit_bytes) is not int or not 0 <= memory_limit_bytes < 2**64:
            raise ValueError("memory_limit_bytes must be uint64")
        self.spec, self.n, self.conditions = spec, spec.n_nucleotides, conditions
        self.lengths = tuple(len(s.sequence) for s in spec.strands)
        if constraints is not None and (not isinstance(constraints, PairingConstraints)
                                       or len(constraints.partners) != self.n):
            raise ValueError("Constraints must cover the complex")
        self.library = C.CDLL(str(Path(library).resolve()))
        self.handle = C.c_void_p()
        self._configure()
        if self.library.aptadp_abi_version() != 1:
            raise RuntimeError("Unsupported native ABI")
        create = getattr(self.library, f"aptadp_create_{factory}_v1")
        create.argtypes = [C.POINTER(NativeInput), C.POINTER(NativeHybridModel), C.POINTER(C.c_void_p),
                           C.POINTER(C.c_double), C.c_char_p, C.c_size_t]
        create.restype = C.c_int
        chemistry = [1 if s.material == "dna" else 2 for s in spec.strands for _ in s.sequence]
        materials = (C.c_uint8 * self.n)(*chemistry)
        cuts = (C.c_uint8 * (self.n+1))(*(int(i in spec.cutpoints) for i in range(self.n+1)))
        allowed = (C.c_uint8 * (self.n*self.n))()
        seq = spec.sequence.replace("U", "T")
        for i in range(self.n):
            for j in range(i+1, self.n):
                pair = seq[i]+seq[j]
                ok = pair in ("AT", "TA", "CG", "GC") or (
                    chemistry[i] == chemistry[j] == 2 and pair in ("GT", "TG"))
                allowed[i*self.n+j] = allowed[j*self.n+i] = ok
        required = (C.c_int32 * self.n)(*(constraints.partners if constraints else (-1,)*self.n))
        request = NativeInput(C.sizeof(NativeInput), self.n, spec.sequence.encode("ascii"),
                              materials, cuts, allowed, required, 3, 30, memory_limit_bytes)
        model = NativeHybridModel(C.sizeof(NativeHybridModel),
                                  1 if sodium_correction == "reference" else 2,
                                  conditions.celsius, conditions.sodium, conditions.magnesium)
        thermal = C.c_double()
        self._call(create, C.byref(request), C.byref(model), C.byref(self.handle), C.byref(thermal))
        self.thermal_energy_kcal_mol = thermal.value
