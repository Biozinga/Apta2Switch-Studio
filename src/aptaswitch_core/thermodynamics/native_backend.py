"""Explicit native experimental RNA/DNA backend, independent of Python RNA.

The selected library implements the Sugimoto95/Vienna04 loop-mean model,
including its declared homogeneous-loop approximations for mixed chemistry.
Single RNA strands use RNAlib's official MFE, partition and sampling algorithms
with the same Turner parameters and pair-salt model as the hybrid calculations.
RNAlib MFE uses centikcal energies; ensemble observables retain continuous PF
weights. These numerical conventions are recorded rather than conflated.
Build checks identify the exact binary and parameter inputs, not experimental
validation or publisher authenticity. No binary is discovered or downloaded.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Real
from pathlib import Path
import platform
import re
import threading

from .evaluation import ComplexAnalysis
from .native_runtime import NativeRuntime
from .memory_policy import reserve_workspace
from .structures import PairingConstraints, parse_structure
from .tubes import StandardStateConvention
from .types import ComplexSpec, Conditions, ModelProvenance, PairProbabilities, PartitionResult


_MODEL = "experimental-Sugimoto95-Vienna04-loopmean"
_RNA_SOLVER = "rnalib2-pf-bpp-sample-mfe-v1"
_RNA_MINIMUM = "turner2004-centikcal-salt-rounded-v1"
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_SOURCE_PREFIX = "src/aptaswitch_core/thermodynamics/"
_REQUIRED_SOURCES = tuple(_SOURCE_PREFIX + "native/" + name for name in (
    "connected_dp.cpp", "connected_dp.hpp", "hybrid_provider.cpp", "hybrid_provider.hpp",
    "hybrid_api.cpp", "hybrid_api.hpp")) + (_SOURCE_PREFIX + "hybrid_motifs.py",)


def _digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _uint64(value, name):
    if type(value) is not int or not 0 <= value < 2**64:
        raise ValueError(f"{name} must be an unsigned 64-bit integer")
    return value


def _tolerance(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be finite and nonnegative")
    try:
        result = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be representable as a finite float") from error
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


@dataclass(frozen=True, slots=True)
class NativeContext:
    complex_spec: ComplexSpec
    conditions: Conditions
    provenance: ModelProvenance
    constraints: PairingConstraints | None
    thermal_energy_kcal_mol: float
    solver: str = ""


@dataclass(frozen=True, slots=True)
class NativeMinimum:
    """A minimum under the recorded solver's energy convention."""
    context: NativeContext
    energy_kcal_mol: float
    structure: str


@dataclass(frozen=True, slots=True)
class NativeSamples:
    """Joint connected structures, not independent draws from pair marginals."""
    context: NativeContext
    structures: tuple[str, ...]
    seed: int


@dataclass(frozen=True, slots=True)
class NativeEnsemble:
    """Probabilities and joint samples computed from one native partition."""
    analysis: ComplexAnalysis
    samples: NativeSamples


@dataclass(frozen=True, slots=True)
class NativeMinimumStructures:
    """Bounded alternatives within tolerance; first is an actual model minimum."""
    context: NativeContext
    minimum_energy_kcal_mol: float
    structures: tuple[str, ...]
    energies_kcal_mol: tuple[float, ...]
    truncated: bool
    maximum_count: int
    tolerance_kcal_mol: float


@dataclass(frozen=True, slots=True)
class NativePairMinima:
    """Minimum total energy subject to pair (i,j); diagonal/impossible = +inf."""
    context: NativeContext
    minimum_energy_kcal_mol: float
    matrix: tuple[tuple[float, ...], ...]


@dataclass(frozen=True, slots=True)
class NativeRegionScreen:
    """All model minima leave [start,end) free, within the recorded tolerance.

    A failure may include a numerically indistinguishable near-minimum. This
    screen says nothing about accessibility throughout the thermal ensemble.
    """
    context: NativeContext
    start: int
    end: int
    all_minima_unpaired: bool
    minimum_energy_kcal_mol: float
    gap_kcal_mol: float
    tolerance_kcal_mol: float
    offending_pair: tuple[int, int] | None
    witness: NativeMinimum | None


class NativeHybridBackend:
    """Connected PF, BPP, labelled minimum and joint sampling on real chemistries.

    Constructor verifies an explicit build manifest; the library is loaded
    lazily at the first valid request. Every calculation owns a fresh native
    handle. By default each request is admitted using fresh available host memory
    and its workspace/output estimate; no fixed cap excludes full reporters.
    An explicit memory_limit_bytes remains the native per-request cap. Admission
    reserves additional headroom but cannot bound unrelated process allocations.
    """

    def __init__(self, library_path, manifest_path=None, *,
                 sodium_correction="published-pair", memory_limit_bytes=None):
        if sodium_correction not in ("reference", "published-pair"):
            raise ValueError("sodium_correction must be 'reference' or 'published-pair'")
        if memory_limit_bytes is not None and _uint64(memory_limit_bytes, "memory_limit_bytes") == 0:
            raise ValueError("memory_limit_bytes must be positive")
        self._library_path = Path(library_path).expanduser().resolve(strict=True)
        manifest_path = (Path(manifest_path).expanduser() if manifest_path is not None else
                         self._library_path.with_suffix(self._library_path.suffix + ".json"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict) or manifest.get("status") != "built" or manifest.get("model") != _MODEL:
            raise ValueError("A completed experimental hybrid build manifest is required")
        if type(manifest.get("abi")) is not int or manifest["abi"] != 1:
            raise ValueError("The hybrid build manifest must declare ABI 1")
        if manifest.get("machine") != platform.machine():
            raise ValueError("The native hybrid build targets another machine architecture")
        if (manifest.get("vienna_version") != "2.7.2"
                or manifest.get("single_rna_solver") != _RNA_SOLVER
                or manifest.get("single_rna_minimum_convention") != _RNA_MINIMUM):
            raise ValueError("An updated hybrid build with the ViennaRNA 2.7.2 single-RNA solver is required")
        storage = manifest.get("connected_dp_storage", "eager-inside-v1")
        if storage not in ("eager-inside-v1", "lazy-inside-v1"):
            raise ValueError("Unsupported connected-DP storage convention in the build manifest")
        self._lazy_inside = storage == "lazy-inside-v1"
        hashes = ("native_library_sha256", "generated_parameters_sha256",
                  "vienna_params_c_sha256", "vienna_static_library_sha256", "vienna_config_h_sha256")
        if any(not isinstance(manifest.get(key), str) or not _DIGEST.fullmatch(manifest[key]) for key in hashes):
            raise ValueError("Native build and parameter SHA-256 provenance is incomplete")
        sources = manifest.get("sources")
        if (not isinstance(sources, dict) or any(name not in sources for name in _REQUIRED_SOURCES)
                or any(not isinstance(name, str) or not isinstance(value, str)
                       or not _DIGEST.fullmatch(value) for name, value in sources.items())):
            raise ValueError("Native build source SHA-256 provenance is incomplete")
        self._binary_digest = manifest["native_library_sha256"]
        self._verify_binary()
        identity = {key: manifest[key] for key in hashes if key != "native_library_sha256"}
        identity.update({key: manifest[key] for key in (
            "vienna_version", "single_rna_solver", "single_rna_minimum_convention")})
        if "connected_dp_storage" in manifest:
            identity["connected_dp_storage"] = storage
        identity.update(sources=sources, model=_MODEL, abi=1, sodium_correction=sodium_correction)
        self._provenance = ModelProvenance(
            backend="AptaSwitch-native-hybrid", version="abi=1;sha256:" + self._binary_digest,
            model=_MODEL,
            parameter_digest="sha256:" + hashlib.sha256(json.dumps(identity, sort_keys=True,
                separators=(",", ":")).encode()).hexdigest(),
            ensemble="connected-pseudoknot-free;dangles=0;coaxial=0;RNA-WC-GU;DNA-WC;hybrid-WC;"
                     "min_loop=3;max_interior=30;pf_smooth=1;betaScale=1;"
                     "RNA-monomer:ViennaRNA-2.7.2;other:connected-DP;"
                     "minimum-RNA:centikcal-salt-rounded;minimum-other:continuous-PF-weights",
            standard_state="1 M;" + ("Na=1M-reference" if sodium_correction == "reference" else
                "Na=0.12..1M;published-pair-S69-S75-S90")
                + ";Mg=0;T=25..55C;native-salt=1.021M;association-S93;uniform-cyclic-labelled-state",
        )
        self._standard_state = StandardStateConvention(self.provenance.standard_state, 1.0, True, True, True)
        self._sodium_correction, self._memory_limit = sodium_correction, memory_limit_bytes
        self._runtime = None
        self._load_lock = threading.Lock()

    @property
    def provenance(self):
        return self._provenance

    @property
    def standard_state(self):
        """The native model is defined at 1 M, with association/symmetry included."""
        return self._standard_state

    @property
    def sodium_correction(self):
        return self._sodium_correction

    @property
    def condition_limits(self):
        """Declared model domain; no extrapolation or silent input clamping."""
        return dict(temperature_min_c=25., temperature_max_c=55.,
                    sodium_min_m=1. if self.sodium_correction == "reference" else .12,
                    sodium_max_m=1., magnesium_min_m=0., magnesium_max_m=0.)

    def _verify_binary(self):
        if _digest_file(self._library_path) != self._binary_digest:
            raise ValueError("Native hybrid binary SHA-256 differs from its build manifest")

    def _get_runtime(self):
        with self._load_lock:
            if self._runtime is None:
                self._verify_binary()
                runtime = NativeRuntime(self._library_path, expected_sha256=self._binary_digest)
                self._verify_binary()
                self._runtime = runtime
            return self._runtime

    def _validate(self, spec, conditions, constraints):
        if not isinstance(spec, ComplexSpec) or not isinstance(conditions, Conditions):
            raise ValueError("Explicit ComplexSpec and Conditions are required")
        conditions.validate_for(spec)
        sodium_ok = (conditions.sodium == 1 if self.sodium_correction == "reference" else
                     .12 <= conditions.sodium <= 1)
        if not sodium_ok or conditions.magnesium != 0 or not 25 <= conditions.celsius <= 55:
            raise ValueError("Experimental hybrid model requires T=25..55 C, Mg=0 and "
                             + ("Na=1 M" if self.sodium_correction == "reference" else "Na=0.12..1 M"))
        n = spec.n_nucleotides
        if not 0 < n <= 100000:
            raise ValueError("Native nucleotide count must be between 1 and 100000")
        # Fail before allocating the input pair mask. The C implementation
        # independently checks its larger DP workspace before allocation.
        if self._memory_limit is not None and n * n * 8 > self._memory_limit:
            raise MemoryError("Native input/output matrices exceed the explicit memory limit")
        if constraints is not None:
            if not isinstance(constraints, PairingConstraints) or len(constraints.partners) != n:
                raise ValueError("Constraints must describe every nucleotide of the complex")
            if all(partner == -1 for partner in constraints.partners):
                constraints = None
        return constraints

    @contextmanager
    def _request(self, spec, conditions, constraints, *, output_bytes=0, python_bytes=0,
                 tables=3, enumerate_minima=False, minimum_only=False):
        lazy_minimum = self._lazy_inside and minimum_only
        # RNAlib constructs MFE/PF arrays together and can create a conditioned
        # calculation for pair minima. Retain its existing conservative bound;
        # the connected-DP table count does not describe that solver.
        if len(spec.strands) == 1 and spec.strands[0].material == "rna":
            tables = 3
            enumerate_minima = False
            lazy_minimum = False
        elif lazy_minimum:
            # These entrypoints own a fresh handle and never fill the PF table.
            # Older libraries keep the eager allocation and its full allowance.
            tables -= 1
        with reserve_workspace(spec.n_nucleotides, explicit_limit=self._memory_limit,
                               output_bytes=output_bytes, python_bytes=python_bytes,
                               tables=tables, enumerate_minima=enumerate_minima,
                               lazy_minimum=lazy_minimum) as limit:
            with self._get_runtime().request(spec, conditions, constraints,
                    1 if self.sodium_correction == "reference" else 2, limit) as request:
                yield request

    def _context(self, spec, conditions, constraints, request):
        return NativeContext(spec, conditions, self.provenance, constraints,
                             request.thermal_energy_kcal_mol, request.solver_name)

    def _partition(self, spec, conditions, constraints, request):
        log_z, thermal = request.log_partition(), request.thermal_energy_kcal_mol
        return PartitionResult(spec, conditions, self.provenance, log_z, -thermal * log_z, thermal, constraints)

    def partition(self, spec, conditions, constraints=None):
        constraints = self._validate(spec, conditions, constraints)
        with self._request(spec, conditions, constraints, tables=1) as request:
            return self._partition(spec, conditions, constraints, request)

    def structure_probability(
        self, spec: ComplexSpec, conditions: Conditions,
        structure: str | PairingConstraints, *, partition: PartitionResult | None = None,
    ) -> float:
        """Exact probability of one complete structure in the full ensemble.

        Evaluate the fully constrained partition and divide by the unrestricted
        partition in log space. This is a joint probability, not a product of
        pair marginals. Impossible or disconnected structures have probability
        zero; an empty unrestricted ensemble raises ``EmptyEnsembleError``.

        An existing unrestricted ``partition`` avoids recalculating its PF.
        Its exact supplied complex (including order and labels), conditions,
        provenance and native thermal-energy convention must match. Dot-bracket
        strings may be flat or include the supplied strand breaks. Partial
        matchings are rejected. Extremely small probabilities may underflow.
        """
        self._validate(spec, conditions, None)
        if isinstance(structure, str):
            structure = parse_structure(structure, tuple(len(s.sequence) for s in spec.strands))
        if (not isinstance(structure, PairingConstraints) or not structure.complete
                or len(structure.partners) != spec.n_nucleotides):
            raise ValueError("Structure must completely specify every nucleotide of the complex")
        if partition is not None and (
                not isinstance(partition, PartitionResult) or partition.complex_spec != spec
                or partition.conditions != conditions or partition.provenance != self.provenance
                or partition.constraints is not None):
            raise ValueError("Partition must be unrestricted and match the exact complex, conditions and model")
        if partition is None:
            partition = self.partition(spec, conditions)
        partition.require_nonempty()
        target = self.partition(spec, conditions, structure)
        if target.thermal_energy_kcal_mol != partition.thermal_energy_kcal_mol:
            raise ValueError("Partition thermal energy differs from the native model convention")
        if target.is_empty:
            return 0.0
        log_probability = target.log_partition - partition.log_partition
        tolerance = 1e-10 + 1e-12 * max(abs(target.log_partition), abs(partition.log_partition))
        if log_probability > tolerance:
            raise ArithmeticError("Structure partition exceeds the unrestricted partition")
        # The exact ratio cannot exceed one; only absorb endpoint roundoff.
        return math.exp(min(0.0, log_probability))

    def analyze(self, spec, conditions, *, need_pairs, constraints=None):
        if type(need_pairs) is not bool:
            raise ValueError("need_pairs must be bool")
        constraints = self._validate(spec, conditions, constraints)
        matrix_cells = spec.n_nucleotides**2 if need_pairs else 0
        with self._request(spec, conditions, constraints, output_bytes=8 * matrix_cells,
                           python_bytes=48 * matrix_cells, tables=2 if need_pairs else 1) as request:
            partition = self._partition(spec, conditions, constraints, request)
            pairs = PairProbabilities(partition, request.pairs()) if need_pairs and not partition.is_empty else None
            return ComplexAnalysis(partition, pairs)

    def pairs(self, spec, conditions, constraints=None):
        result = self.analyze(spec, conditions, need_pairs=True, constraints=constraints)
        result.partition.require_nonempty()
        return result.pairs

    def minimum(self, spec, conditions, constraints=None):
        constraints = self._validate(spec, conditions, constraints)
        with self._request(spec, conditions, constraints, tables=2, minimum_only=True) as request:
            energy, structure = request.minimum()
            if not math.isfinite(energy):
                raise ArithmeticError("Native minimum energy is nonfinite")
            return NativeMinimum(self._context(spec, conditions, constraints, request), energy, structure)

    def sample(self, spec, conditions, count, *, seed, constraints=None):
        _uint64(count, "count")
        _uint64(seed, "seed")
        constraints = self._validate(spec, conditions, constraints)
        if self._memory_limit is not None and count * (spec.n_nucleotides + 1) > self._memory_limit:
            raise MemoryError("Sample output exceeds the explicit memory limit")
        output = count * (spec.n_nucleotides + 1)
        with self._request(spec, conditions, constraints, output_bytes=output,
                           python_bytes=3 * output + 64 * count, tables=1) as request:
            self._partition(spec, conditions, constraints, request).require_nonempty()
            structures = request.sample(count, seed)
            return NativeSamples(self._context(spec, conditions, constraints, request), structures, seed)

    def ensemble(self, spec, conditions, count, *, seed, constraints=None):
        """Reuse the inside table for pair probabilities and joint traceback."""
        _uint64(count, "count")
        _uint64(seed, "seed")
        constraints = self._validate(spec, conditions, constraints)
        if self._memory_limit is not None and count * (spec.n_nucleotides + 1) + 8 * spec.n_nucleotides**2 > self._memory_limit:
            raise MemoryError("Ensemble output exceeds the explicit memory limit")
        output = count * (spec.n_nucleotides + 1)
        cells = spec.n_nucleotides**2
        with self._request(spec, conditions, constraints, output_bytes=output + 8 * cells,
                           python_bytes=48 * cells + 3 * output + 64 * count, tables=2) as request:
            partition = self._partition(spec, conditions, constraints, request)
            partition.require_nonempty()
            pairs = PairProbabilities(partition, request.pairs())
            samples = NativeSamples(self._context(spec, conditions, constraints, request),
                                    request.sample(count, seed), seed)
            return NativeEnsemble(ComplexAnalysis(partition, pairs), samples)

    def minimum_structures(self, spec, conditions, constraints=None, *,
                           absolute_tolerance=1e-8, maximum_count=32):
        tolerance = _tolerance(absolute_tolerance, "absolute_tolerance")
        if _uint64(maximum_count, "maximum_count") == 0:
            raise ValueError("maximum_count must be positive")
        constraints = self._validate(spec, conditions, constraints)
        if self._memory_limit is not None and maximum_count * (spec.n_nucleotides + 1 + 8) > self._memory_limit:
            raise MemoryError("Minimum-structure output exceeds the explicit memory limit")
        output = maximum_count * (spec.n_nucleotides + 1 + 8)
        with self._request(spec, conditions, constraints, output_bytes=output,
                           python_bytes=3 * output + 128 * maximum_count, tables=2,
                           enumerate_minima=True, minimum_only=True) as request:
            minimum, structures, energies, truncated = request.minimum_structures(tolerance, maximum_count)
            if not math.isfinite(minimum) or any(not math.isfinite(energy) for energy in energies):
                raise ArithmeticError("Native minimum-structure energies are nonfinite")
            numerical = 1e-9 + 1e-12 * max(abs(minimum), *(abs(e) for e in energies))
            if (abs(energies[0] - minimum) > numerical or
                    any(e < minimum - numerical or e > minimum + tolerance + numerical for e in energies)):
                raise ArithmeticError("Native alternatives violate their global energy window")
            if len(set(structures)) != len(structures) or (truncated and len(structures) != maximum_count):
                raise ArithmeticError("Native alternatives contain duplicates or an inconsistent truncation flag")
            return NativeMinimumStructures(self._context(spec, conditions, constraints, request), minimum,
                structures, energies, truncated, maximum_count, tolerance)

    def minimum_pair_energies(self, spec, conditions, constraints=None):
        constraints = self._validate(spec, conditions, constraints)
        cells = spec.n_nucleotides**2
        with self._request(spec, conditions, constraints, output_bytes=8 * cells,
                           python_bytes=48 * cells, tables=3, minimum_only=True) as request:
            minimum, matrix = request.minimum_pair_energies()
            if not math.isfinite(minimum):
                raise ArithmeticError("Native minimum energy is nonfinite")
            for i, row in enumerate(matrix):
                if row[i] != math.inf or any(math.isnan(e) or e == -math.inf for e in row):
                    raise ArithmeticError("Native pair minima contain invalid energies")
                if any(row[j] != matrix[j][i] for j in range(i)):
                    raise ArithmeticError("Native pair minima are not symmetric")
                if any(e < minimum - (1e-8 + 1e-12 * max(abs(e), abs(minimum))) for e in row):
                    raise ArithmeticError("Pair-conditioned minimum is below the unconstrained minimum")
            return NativePairMinima(self._context(spec, conditions, constraints, request), minimum, matrix)

    def screen_unpaired_region(self, spec, conditions, *, start, end, constraints=None,
                               absolute_tolerance=1e-8, relative_tolerance=1e-12):
        """Check all labelled minima, with global zero-based [start,end) indices.

        This is a minimum-energy screen, not a probability cutoff. A failure
        includes a genuine conditional traceback; no minimum is fabricated by
        concatenating marginal pair choices. Empty ensembles raise explicitly.
        """
        atol = _tolerance(absolute_tolerance, "absolute_tolerance")
        rtol = _tolerance(relative_tolerance, "relative_tolerance")
        constraints = self._validate(spec, conditions, constraints)
        if type(start) is not int or type(end) is not int or not 0 <= start < end <= spec.n_nucleotides:
            raise ValueError("Region must be a nonempty half-open interval in complex coordinates")
        minima = self.minimum_pair_energies(spec, conditions, constraints)
        candidates = ((minima.matrix[i][j], min(i, j), max(i, j))
                      for i in range(start, end) for j in range(spec.n_nucleotides) if i != j)
        energy, i, j = min(candidates, default=(math.inf, -1, -1))
        scale = max(abs(minima.minimum_energy_kcal_mol), abs(energy) if math.isfinite(energy) else 0.)
        tolerance = atol + rtol * scale
        if not math.isfinite(tolerance):
            raise ValueError("Combined energy tolerance overflows")
        gap = energy - minima.minimum_energy_kcal_mol
        if gap < -tolerance:
            raise ArithmeticError("Pair-conditioned minimum is below the global minimum")
        passed, witness, pair = gap > tolerance, None, None
        if not passed:
            pair = (i, j)
            conditioned = constraints or PairingConstraints.unconstrained(spec.n_nucleotides)
            if conditioned.partners[i] != j:
                conditioned = conditioned.assign(i, j)
            witness = self.minimum(spec, conditions, conditioned)
            state = parse_structure(witness.structure, tuple(len(s.sequence) for s in spec.strands))
            if state.partners[i] != j or abs(witness.energy_kcal_mol - energy) > tolerance:
                raise ArithmeticError("Conditional traceback does not reproduce its pair minimum")
        return NativeRegionScreen(minima.context, start, end, passed, minima.minimum_energy_kcal_mol,
                                  gap, tolerance, pair, witness)
