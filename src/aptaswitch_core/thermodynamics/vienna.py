"""Homogeneous connected-ensemble adapter for the audited ViennaRNA 2.7.2.

Mixed complexes must be evaluated by the forthcoming material-aware engine;
this adapter never transcribes DNA into RNA. Native stochastic traceback is
deliberately avoided: connected-complex sampling uses conditional partitions.
"""
from __future__ import annotations

import hashlib
import importlib
import math
import threading

from .structures import PairingConstraints, is_connected, parse_structure
from .sampling import ConditionalSampler
from .runtime_identity import identify_vienna_runtime
from .types import ComplexSpec, Conditions, ModelProvenance, PartitionResult, PairProbabilities


_NATIVE_LOCK = threading.RLock()
_ENSEMBLE = (
    "connected-pseudoknot-free;dangles=0;noLP=0;noGU=0;min_loop=3;"
    "pf_smooth=1;betaScale=1;special_hp=1;noGUclosure=0;logML=0;"
    "circ=0;gquad=0;energy_set=0;max_bp_span=-1"
)


def _model_details(RNA, conditions: Conditions, *, with_pairs: bool):
    """Create fixed model details *after* loading the requested material tables.

    RNA.md() otherwise inherits mutable process defaults, including PF energy
    rounding, pair rules and Boltzmann temperature. Constructor arguments also
    rebuild the internal pairing alphabet; assigning energy_set afterwards
    would leave its derived tables stale. Material loaders set the geometry
    and duplex salt factor, which must remain material-specific: native
    parameter caching compares model details, not the loaded energy tables.
    """
    return RNA.md(
        temperature=conditions.celsius, salt=conditions.sodium,
        betaScale=1.0, pf_smooth=1, dangles=0, special_hp=1,
        noLP=0, noGU=0, noGUclosure=0, logML=0, circ=0,
        circ_penalty=1, gquad=0, uniq_ML=1, energy_set=0,
        backtrack=1, backtrack_type="F", compute_bpp=int(with_pairs),
        max_bp_span=-1, min_loop_size=3, window_size=-1,
        oldAliEn=0, ribo=0, cv_fact=1.0, nc_fact=1.0, sfact=1.07,
        saltMLLower=6, saltMLUpper=24, saltDPXInit=99999,
    )


class UnsupportedChemistryError(ValueError):
    """A backend cannot represent the requested chemical ensemble."""


class ViennaBackend:
    """Partition and marginal calculations with explicit parameter provenance.

    dangles=0 keeps a secondary structure as one ensemble state, without
    additional stacking substates. Other ensemble conventions require a
    separate validation and are not silently enabled here.
    """

    def __init__(self, *, build_manifest=None):
        try:
            self.RNA = importlib.import_module("RNA")
        except ImportError as error:
            raise RuntimeError("ViennaRNA 2.7.2 Python bindings are required") from error
        if self.RNA.__version__ != "2.7.2":
            raise RuntimeError("This adapter is validated for ViennaRNA 2.7.2 only")
        self.runtime_identity = identify_vienna_runtime(self.RNA, build_manifest)

    def _calculate(self, complex_spec: ComplexSpec, conditions: Conditions,
                   constraints: PairingConstraints | None, *, with_pairs: bool):
        conditions.validate_for(complex_spec)
        materials = {s.material for s in complex_spec.strands}
        if len(materials) != 1:
            raise UnsupportedChemistryError("ViennaRNA's homogeneous backend cannot evaluate a DNA–RNA complex")
        if conditions.magnesium != 0 or conditions.sodium <= 0:
            raise ValueError("This backend requires positive sodium and zero magnesium")
        n = complex_spec.n_nucleotides
        if constraints is not None and len(constraints.partners) != n:
            raise ValueError("Constraints must have one entry per nucleotide")
        chemistry = next(iter(materials))
        # hc_add_bp can enable a noncanonical pair, or silently ignore one
        # shorter than min_loop_size. Such constraints specify an empty subset
        # of our ensemble, not a different pairing model.
        impossible_constraints = False
        if constraints is not None:
            allowed = {"AU", "UA", "CG", "GC", "GU", "UG"} if chemistry == "rna" else {
                "AT", "TA", "CG", "GC", "GT", "TG"}
            strand_ids = [k for k, strand in enumerate(complex_spec.strands)
                          for _ in strand.sequence]
            sequence = complex_spec.sequence
            impossible_constraints = any(
                sequence[i] + sequence[j] not in allowed
                or (strand_ids[i] == strand_ids[j] and j - i <= 3)
                for i, j in enumerate(constraints.partners) if j > i
            )
        RNA = self.RNA
        # Global parameter loading and native work stay serialized within this
        # adapter. Separate application workers may run separate processes.
        with _NATIVE_LOCK:
            if chemistry == "rna":
                RNA.params_load_RNA_Turner2004()
                name, text = "RNA-Turner2004", RNA.parameter_set_rna_turner2004
            else:
                RNA.params_load_DNA_Mathews2004()
                name, text = "DNA-Mathews2004", RNA.parameter_set_dna_mathews2004
            provenance = ModelProvenance(
                backend="ViennaRNA", version=self.runtime_identity.version_label, model=name,
                parameter_digest="sha256:" + hashlib.sha256(text.encode()).hexdigest(),
                ensemble=_ENSEMBLE,
                standard_state="1 M;native association;native cyclic symmetry",
            )
            md = _model_details(RNA, conditions, with_pairs=with_pairs)
            # Stock 2.7.2 omits salt-loop slot 31. A hairpin with 30 free
            # bases (32 nt including its closing pair), or a 29-base
            # interior, then silently disappears even from an otherwise
            # finite unconstrained PF. Checking only PF=+inf misses this.
            if n >= 32 and conditions.sodium != RNA.MODEL_DEFAULT_SALT:
                snapshot = RNA.exp_param(md)
                if RNA.doubleP_getitem(snapshot.expSaltLoop, 31) == 0:
                    raise ArithmeticError(
                        "ViennaRNA 2.7.2 needs the audited salt-loop-31 fix "
                        "before evaluating complexes of 32 or more nucleotides at this sodium concentration")
            fc = RNA.fold_compound("&".join(s.sequence for s in complex_spec.strands), md)
            if constraints is not None and not impossible_constraints:
                for i, j in enumerate(constraints.partners):
                    if i == j:
                        fc.hc_add_up(i + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
                    elif j > i:
                        fc.hc_add_bp(i + 1, j + 1,
                            RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE)
            if impossible_constraints:
                empty = True
            else:
                mfe_structure, mfe = fc.mfe()
                if not math.isfinite(mfe):
                    raise ArithmeticError("Native MFE returned a nonfinite numerical result")
                empty = mfe >= RNA.INF / 100.0
            if empty:
                energy = math.inf
                # Native PF parameters are allocated lazily. The empty subset
                # still carries the same thermal convention as a nonempty one.
                if fc.exp_params is None:
                    fc.exp_params_rescale(0.0)
            else:
                fc.exp_params_rescale(mfe)
                _, energy = fc.pf()
                # Native multistrand MFE can return an empty traceback and a
                # finite zero even when the connected PF is empty. Positive
                # infinity is therefore the native zero-partition result.
                if energy == -math.inf or math.isnan(energy):
                    raise ArithmeticError("Native partition calculation failed numerically")
                if energy == math.inf:
                    # A finite, connected MFE witness proves Z > 0. Without
                    # that witness, native multistrand MFE's empty traceback
                    # cannot distinguish a zero partition from a failure.
                    lengths = tuple(len(strand.sequence) for strand in complex_spec.strands)
                    try:
                        witness = parse_structure(mfe_structure, lengths)
                    except ValueError:
                        witness = None
                    if (witness is not None and is_connected(witness, lengths)
                            and (constraints is None or all(
                                j < 0 or witness.partners[i] == j
                                for i, j in enumerate(constraints.partners)))):
                        raise ArithmeticError("Native partition is zero despite a finite connected MFE structure")
            thermal = fc.exp_params.kT / 1000.0
            partition = PartitionResult(
                complex_spec=complex_spec, conditions=conditions, provenance=provenance,
                log_partition=-energy / thermal, free_energy_kcal_mol=energy,
                thermal_energy_kcal_mol=thermal,
                constraints=constraints,
            )
            if not with_pairs:
                return partition
            partition.require_nonempty()
            upper = fc.bpp()
            if len(upper) != n + 1 or any(len(row) != n + 1 for row in upper):
                raise ArithmeticError("Native pair-probability dimensions are inconsistent")
            matrix = [[0.0] * n for _ in range(n)]
            tolerance = 1e-10
            for i in range(n):
                for j in range(i + 1, n):
                    p = upper[i + 1][j + 1]
                    if not math.isfinite(p) or p < -tolerance or p > 1 + tolerance:
                        raise ArithmeticError("Native pair probability is outside [0, 1]")
                    # Only floating-point endpoint overshoot is adjusted.
                    matrix[i][j] = matrix[j][i] = min(1.0, max(0.0, p))
            for i, row in enumerate(matrix):
                unpaired = 1 - math.fsum(row)
                if unpaired < -tolerance:
                    raise ArithmeticError("Native pair probabilities sum to more than one")
                row[i] = max(0.0, unpaired)
            result = PairProbabilities(partition=partition, matrix=tuple(tuple(row) for row in matrix))
            if constraints is not None:
                for i, j in enumerate(constraints.partners):
                    if j >= 0 and abs(result.matrix[i][j] - 1) > 1e-8:
                        raise ArithmeticError("Native partition did not enforce a hard constraint")
            return result

    def partition(self, complex_spec: ComplexSpec, conditions: Conditions,
                  constraints: PairingConstraints | None = None) -> PartitionResult:
        return self._calculate(complex_spec, conditions, constraints, with_pairs=False)

    def pairs(self, complex_spec: ComplexSpec, conditions: Conditions,
              constraints: PairingConstraints | None = None) -> PairProbabilities:
        return self._calculate(complex_spec, conditions, constraints, with_pairs=True)

    def sampler(self, complex_spec: ComplexSpec, conditions: Conditions, *, seed: int,
                cache_size: int = 4096) -> ConditionalSampler:
        return ConditionalSampler(complex_spec, conditions,
            lambda constraints: self.pairs(complex_spec, conditions, constraints),
            seed=seed, cache_size=cache_size)
