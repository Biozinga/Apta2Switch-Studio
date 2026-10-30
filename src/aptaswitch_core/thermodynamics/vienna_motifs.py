"""Homogeneous structure energies assembled from ViennaRNA's public loop weights.

This is an independent loop traversal, not a native structure evaluator or a
mixed-material model. Energies use continuous PF parameters at the requested
temperature. The uniform cyclic correction belongs to labelled-state weights;
it must not be described as a structure-specific physical MFE correction.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import importlib
import math

from .loops import Loop, decompose_structure
from .runtime_identity import identify_vienna_runtime
from .structures import PairingConstraints, format_structure, is_connected, parse_structure
from .types import ComplexSpec, Conditions, ModelProvenance
from .vienna import UnsupportedChemistryError, _ENSEMBLE, _NATIVE_LOCK, _model_details


_BASE_CODE = {"A": 1, "C": 2, "G": 3, "U": 4, "T": 4}
_PAIR_CODE = {"CG": 1, "GC": 2, "GU": 3, "GT": 3,
              "UG": 4, "TG": 4, "AU": 5, "AT": 5, "UA": 6, "TA": 6}


def validate_loop_face(spec: ComplexSpec, loop: Loop, pair_allowed) -> bool:
    """Check local coordinates before native table access; false means zero weight.

    Full planar topology/connectivity is validated when the matching is
    decomposed. Both homogeneous and mixed oracles use the same local guards.
    """
    if not isinstance(loop, Loop) or loop.kind not in (
            "hairpin", "stack", "bulge", "interior", "multiloop", "exterior"):
        raise ValueError("A recognized decomposed loop is required")
    if any(not 0 <= a <= b <= spec.n_nucleotides for a, b in loop.unpaired_runs):
        raise ValueError("Loop free-base intervals are outside the complex")
    pairs = loop.child_pairs + (() if loop.closing_pair is None else (loop.closing_pair,))
    if any(i >= j or not pair_allowed(spec, i, j) for i, j in pairs):
        return False
    if loop.kind != "exterior" and loop.closing_pair is None:
        raise ValueError("A non-exterior face needs a closing pair")
    if loop.kind in ("stack", "bulge", "interior") and len(loop.child_pairs) != 1:
        raise ValueError("An interior face must have exactly one child pair")
    if loop.kind in ("stack", "bulge", "interior"):
        i, j = loop.closing_pair
        k, l = loop.child_pairs[0]
        if not i < k < l < j or loop.unpaired_count != k-i-1+j-l-1:
            raise ValueError("Interior face coordinates are inconsistent")
        left, right = k-i-1, j-l-1
        kind = "stack" if left == right == 0 else "bulge" if min(left, right) == 0 else "interior"
        if loop.kind != kind:
            raise ValueError("Interior face kind does not match its coordinates")
    if loop.kind == "hairpin" and (loop.child_pairs or
            loop.unpaired_count != loop.closing_pair[1]-loop.closing_pair[0]-1):
        raise ValueError("Hairpin face coordinates are inconsistent")
    if ((loop.kind == "hairpin" and loop.unpaired_count < 3)
            or (loop.kind in ("interior", "bulge") and loop.unpaired_count > 30)):
        return False
    return True


class ViennaMotifModel:
    """One RNA or DNA model suitable for ``ExhaustiveBackend`` callbacks.

    DNA GT/TG pairs are admitted deliberately, matching the homogeneous
    ``ViennaBackend`` convention ``noGU=False``. Other chemistry/pairing models
    need a separate oracle and parameter provenance.
    """

    def __init__(self, material: str = "rna", *, build_manifest=None):
        if material not in ("rna", "dna"):
            raise ValueError("material must be 'rna' or 'dna'")
        try:
            self.RNA = importlib.import_module("RNA")
        except ImportError as error:
            raise RuntimeError("ViennaRNA 2.7.2 Python bindings are required") from error
        if self.RNA.__version__ != "2.7.2":
            raise RuntimeError("This motif oracle is validated for ViennaRNA 2.7.2 only")
        self.runtime_identity = identify_vienna_runtime(self.RNA, build_manifest)
        self.material = material
        parameter_text = (self.RNA.parameter_set_rna_turner2004 if material == "rna"
                          else self.RNA.parameter_set_dna_mathews2004)
        self.provenance = ModelProvenance(
            backend="ViennaRNA-loop-weight-oracle", version=self.runtime_identity.version_label,
            model="RNA-Turner2004" if material == "rna" else "DNA-Mathews2004",
            parameter_digest="sha256:" + hashlib.sha256(parameter_text.encode()).hexdigest(),
            ensemble=_ENSEMBLE + ";max_interior=30",
            standard_state="1 M;native association;uniform cyclic correction in labelled-state weights",
        )

    def _validate_spec(self, spec: ComplexSpec) -> None:
        if not isinstance(spec, ComplexSpec):
            raise ValueError("complex_spec must be a ComplexSpec")
        if any(strand.material != self.material for strand in spec.strands):
            raise UnsupportedChemistryError("The motif oracle requires one matching homogeneous material")

    @staticmethod
    def _validate_conditions(conditions: Conditions) -> None:
        if not isinstance(conditions, Conditions):
            raise ValueError("conditions must be Conditions")
        if conditions.magnesium != 0 or conditions.sodium <= 0:
            raise ValueError("This oracle requires positive sodium and zero magnesium")

    @lru_cache(maxsize=32)
    def _parameters(self, conditions: Conditions):
        self._validate_conditions(conditions)
        with _NATIVE_LOCK:
            if self.material == "rna":
                self.RNA.params_load_RNA_Turner2004()
            else:
                self.RNA.params_load_DNA_Mathews2004()
            md = _model_details(self.RNA, conditions, with_pairs=True)
            return self.RNA.exp_param(md)

    def thermal_energy(self, conditions: Conditions) -> float:
        """Return the native RT, in kcal/mol, used by this parameter snapshot."""
        self._validate_conditions(conditions)
        return float(self._parameters(conditions).kT) / 1000.0

    def pair_allowed(self, spec: ComplexSpec, i: int, j: int) -> bool:
        self._validate_spec(spec)
        if any(type(k) is not int or not 0 <= k < spec.n_nucleotides for k in (i, j)):
            raise ValueError("Pair indices must be valid zero-based nucleotide positions")
        return i != j and spec.sequence[i] + spec.sequence[j] in _PAIR_CODE

    @staticmethod
    def _energy_of_weight(weight: float, thermal: float) -> float:
        weight = float(weight)
        if weight == 0:
            return math.inf
        if not math.isfinite(weight) or weight < 0:
            raise ArithmeticError("Native loop weight is negative or nonfinite")
        return -thermal * math.log(weight)

    def energy(self, spec: ComplexSpec, conditions: Conditions,
               state: PairingConstraints | str) -> float:
        """Evaluate a complete connected matching, or return +inf outside its grammar.

        Invalid syntax/partial matchings raise ValueError. Disconnected states,
        disallowed pairs, short hairpins and interiors with >30 unpaired bases
        have zero ensemble weight. No unsafe native interior-table lookup is
        attempted for those large loops.
        """
        self._validate_spec(spec)
        self._validate_conditions(conditions)
        lengths = tuple(len(strand.sequence) for strand in spec.strands)
        if isinstance(state, str):
            state = parse_structure(state, lengths)
        if not isinstance(state, PairingConstraints) or not state.complete or len(state.partners) != spec.n_nucleotides:
            raise ValueError("A complete matching of the supplied complex is required")
        if not is_connected(state, lengths):
            return math.inf
        if any(not self.pair_allowed(spec, i, j) for i, j in enumerate(state.partners) if j > i):
            return math.inf
        loops = decompose_structure(spec, format_structure(state, lengths))
        if any((loop.kind == "hairpin" and loop.unpaired_count < 3)
               or (loop.kind in ("interior", "bulge") and loop.unpaired_count > 30) for loop in loops):
            return math.inf
        with _NATIVE_LOCK:
            params = self._parameters(conditions)
            thermal = float(params.kT) / 1000.0
            contributions = [thermal * math.log(spec.symmetry_factor)]
            if len(spec.strands) > 1:
                contributions.append((len(spec.strands) - 1) * self._energy_of_weight(params.expDuplexInit, thermal))
            contributions.extend(self.loop_energy(spec, conditions, loop) for loop in loops)
            return math.fsum(contributions)

    def loop_energy(self, spec: ComplexSpec, conditions: Conditions, loop: Loop) -> float:
        """One decomposed face, excluding association and cyclic corrections.

        ``loop`` must be a face obtained from ``decompose_structure`` (or an
        equivalent validated planar grammar). Its boundary coordinates refer
        to ``spec``. Local geometry and array bounds are checked before native
        calls; connectivity is a property of the complete structure.
        """
        self._validate_spec(spec)
        self._validate_conditions(conditions)
        if not validate_loop_face(spec, loop, self.pair_allowed):
            return math.inf
        sequence = spec.sequence
        codes = [_BASE_CODE[base] for base in sequence]
        pair_type = lambda i, j: _PAIR_CODE[sequence[i] + sequence[j]]
        with _NATIVE_LOCK:
            params = self._parameters(conditions)
            thermal = float(params.kT) / 1000.0
            e = lambda weight: self._energy_of_weight(weight, thermal)
            contributions = []
            # The unmodified 2.7.2 constructor leaves salt-loop slot 31
            # zero although the weight APIs read it for these motifs.
            # Report this native failure instead of inventing an empty
            # physical ensemble. A separately audited native fix is needed.
            salt_backbones = (loop.unpaired_count + 1 if loop.kind == "hairpin" else
                              loop.unpaired_count + 2 if loop.kind in ("bulge", "interior") else 0)
            if (salt_backbones == 31 and conditions.sodium != self.RNA.MODEL_DEFAULT_SALT
                    and self.RNA.doubleP_getitem(params.expSaltLoop, 31) == 0):
                raise ArithmeticError("ViennaRNA 2.7.2 leaves the 31-backbone salt-loop weight uninitialized")
            if loop.kind == "hairpin":
                i, j = loop.closing_pair
                contributions.append(e(self.RNA.exp_E_Hairpin(j-i-1, pair_type(i, j),
                    codes[i+1], codes[j-1], sequence[i:j+1].replace("T", "U"), params)))
            elif loop.kind in ("stack", "bulge", "interior"):
                i, j = loop.closing_pair
                k, l = loop.child_pairs[0]
                contributions.append(e(self.RNA.exp_E_IntLoop(k-i-1, j-l-1,
                    pair_type(i, j), pair_type(l, k), codes[i+1], codes[j-1],
                    codes[k-1], codes[l+1], params)))
            else:
                stem_function = self.RNA.exp_E_MLstem if loop.kind == "multiloop" else self.RNA.exp_E_ExtLoop
                if loop.kind == "multiloop":
                    contributions.append(e(params.expMLclosing))
                    if loop.unpaired_count:
                        contributions.append(loop.unpaired_count * e(params.expMLbase))
                if loop.closing_pair is not None:
                    i, j = loop.closing_pair
                    contributions.append(e(stem_function(pair_type(j, i), -1, -1, params)))
                for i, j in loop.child_pairs:
                    contributions.append(e(stem_function(pair_type(i, j), -1, -1, params)))
            return math.fsum(contributions)
