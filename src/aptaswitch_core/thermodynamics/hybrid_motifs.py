"""Experimental RNA/DNA motif oracle with explicit homogeneous-loop approximations.

Version 1 is a nominal 1 M Na, Mg-free reference model, 25--55 Celsius.
Vienna's monomaterial parameter tables are used at their reference salt
(1.021 M in its API), disabling its additional sodium correction. Mixed
interiors use complete RNA/DNA loop-energy means, inspired by Lorenz,
Hofacker and Bernhart (2012), not the empirical hybrid mismatch model of 2026.
Mixed multiloops average all three coefficients, including free-base cost.

Hybrid WC stacks/initiation: Sugimoto et al., Biochemistry 34 (1995),
11211--11216; public R_DNA_NN1 table in Biopython 1.86. Hybrid terminal cost
is zero, consistently with that fit. Association composition interpolation
is the published S93 rule of DOI 10.1021/acssynbio.5c00470, applied to these
different monomaterial tables: this is an explicitly experimental extension.
No hybrid wobble values, dangles or coaxial stacks are fabricated.

Optional version 2 (``sodium_correction='published-pair'``) adds the common
per-pair entropy correction published for DNA (S69--S70), RNA (S73--S75),
and mixed RNA/DNA (S90--S92) in that supplement, restricted here to 0.12--1 M.
This does not validate the experimental fusion of its different loop tables.
https://acs.figshare.com/articles/journal_contribution/31017451

The small H/S table below is from Biopython's MeltingTemp.py:
Copyright 2004-2008 Sebastian Bassi; Copyright 2013-2018 Markus Piotrowski.
Copyright (c) 1999-2024, The Biopython Contributors. All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math

from .loops import Loop, decompose_structure
from .structures import PairingConstraints, format_structure, is_connected, parse_structure
from .types import ComplexSpec, Conditions, ModelProvenance, Strand
from .vienna import _NATIVE_LOCK
from .vienna_motifs import ViennaMotifModel, _PAIR_CODE, validate_loop_face


# Keys: RNA 5'->3' (U encoded T) / DNA 3'->5'. Values: kcal/mol, cal/mol/K.
_SUGIMOTO_HS = {
    "init": (1.9, -3.9),
    "TT/AA": (-11.5, -36.4), "GT/CA": (-7.8, -21.6), "CT/GA": (-7.0, -19.7),
    "AT/TA": (-8.3, -23.9), "TG/AC": (-10.4, -28.4), "GG/CC": (-12.8, -31.9),
    "CG/GC": (-16.3, -47.1), "AG/TC": (-9.1, -23.5), "TC/AG": (-8.6, -22.9),
    "GC/CG": (-8.0, -17.1), "CC/GG": (-9.3, -23.2), "AC/TG": (-5.9, -12.3),
    "TA/AT": (-7.8, -23.2), "GA/CT": (-5.5, -13.5), "CA/GT": (-9.0, -26.1),
    "AA/TT": (-7.8, -21.9),
}
_WC = frozenset(("AT", "TA", "GC", "CG"))


class HybridMotifModel:
    """Complete connected-state energies within the declared restricted grammar.

    Each input strand retains its real chemistry. RNA permits WC/GU pairs;
    DNA and hybrid pairs permit WC only. The local ``loop_energy`` interface
    excludes association and cyclic factors, for independent grammar use.
    """

    def __init__(self, *, sodium_correction: str = "reference", build_manifest=None):
        if sodium_correction not in ("reference", "published-pair"):
            raise ValueError("sodium_correction must be 'reference' or 'published-pair'")
        self._sodium_correction = sodium_correction
        version = 1 if sodium_correction == "reference" else 2
        self._mono = {material: ViennaMotifModel(material, build_manifest=build_manifest)
                      for material in ("rna", "dna")}
        self.runtime_identity = self._mono["rna"].runtime_identity
        self.RNA = self._mono["rna"].RNA
        policy = {
            "hybrid_hs": _SUGIMOTO_HS,
            "monomaterial_parameters": {m: obj.provenance.parameter_digest for m, obj in self._mono.items()},
            "version": version, "native_salt_reference_m": 1.021,
            "sodium_correction": sodium_correction,
            "pair_salt_kcal_mol": "-0.114*ln(Na_M)*T_K/310.15;all-pairs" if version == 2 else "none",
            "hybrid_terminal": 0, "interior": "complete-projection-energy-mean",
            "bulge1": "mean-size-initiation-plus-hybrid-stack",
            "multiloop": "mean-initiation-branch-unpaired;actual-pair-terminal",
            "association": "composition-S93-around-Sugimoto-initiation",
        }
        self.provenance = ModelProvenance(
            backend="experimental-hybrid-motif-oracle", version=f"{version};" + self.runtime_identity.version_label,
            model="experimental-Sugimoto95-Vienna04-loopmean",
            parameter_digest="sha256:" + hashlib.sha256(json.dumps(policy, sort_keys=True).encode()).hexdigest(),
            ensemble="connected-pseudoknot-free;dangles=0;coaxial=0;RNA-WC-GU;DNA-WC;hybrid-WC;min_loop=3;max_interior=30",
            standard_state="1 M;" + ("Na=1M-reference" if version == 1 else "Na=0.12..1M;published-pair-S69-S75-S90")
                + ";Mg=0;T=25..55C;native-salt=1.021M;association-S93;uniform-cyclic-labelled-state",
        )

    @property
    def sodium_correction(self) -> str:
        return self._sodium_correction

    @staticmethod
    def _validate_spec(spec: ComplexSpec) -> None:
        if not isinstance(spec, ComplexSpec):
            raise ValueError("complex_spec must be a ComplexSpec")

    def _validate_conditions(self, conditions: Conditions) -> None:
        if not isinstance(conditions, Conditions):
            raise ValueError("conditions must be Conditions")
        valid_sodium = (conditions.sodium == 1.0 if self._sodium_correction == "reference"
                        else .12 <= conditions.sodium <= 1.0)
        if not valid_sodium or conditions.magnesium != 0 or not 25 <= conditions.celsius <= 55:
            sodium_range = "1 M" if self._sodium_correction == "reference" else "0.12..1 M"
            raise ValueError(f"Experimental hybrid {self._sodium_correction} requires Na={sodium_range}, Mg=0 and 25..55 Celsius")

    @staticmethod
    @lru_cache(maxsize=128)
    def _materials(spec: ComplexSpec) -> tuple[str, ...]:
        return tuple(strand.material for strand in spec.strands for _ in strand.sequence)

    @staticmethod
    @lru_cache(maxsize=128)
    def _projection(spec: ComplexSpec, material: str) -> ComplexSpec:
        # These are parameter lookup projections, never returned as candidates.
        return ComplexSpec(tuple(Strand(s.sequence.replace("T", "U") if material == "rna"
            else s.sequence.replace("U", "T"), material, s.name) for s in spec.strands))

    def _reference(self, conditions: Conditions) -> Conditions:
        self._validate_conditions(conditions)
        return Conditions(celsius=conditions.celsius, sodium=1.021, magnesium=0)

    def thermal_energy(self, conditions: Conditions) -> float:
        return self._mono["rna"].thermal_energy(self._reference(conditions))

    def pair_salt_energy(self, conditions: Conditions) -> float:
        """Per-pair kcal/mol correction, H=0 and S=114/310.15*ln(Na) cal/mol/K.

        This applies equally to RNA, DNA and hybrid pairs. The reference
        parameter snapshots never receive the requested non-reference salt.
        """
        self._validate_conditions(conditions)
        if self._sodium_correction == "reference" or conditions.sodium == 1.0:
            return 0.0
        return -.114 * math.log(conditions.sodium) * conditions.kelvin / 310.15

    def pair_allowed(self, spec: ComplexSpec, i: int, j: int) -> bool:
        self._validate_spec(spec)
        if any(type(k) is not int or not 0 <= k < spec.n_nucleotides for k in (i, j)):
            raise ValueError("Pair indices must be valid zero-based nucleotide positions")
        if i == j:
            return False
        pair = (spec.sequence[i] + spec.sequence[j]).replace("U", "T")
        materials = self._materials(spec)
        return pair in _WC or (materials[i] == materials[j] == "rna" and pair in ("GT", "TG"))

    @staticmethod
    def _hs_energy(key: str, conditions: Conditions) -> float:
        enthalpy, entropy = _SUGIMOTO_HS[key]
        return enthalpy - conditions.kelvin * entropy / 1000

    def association_energy(self, spec: ComplexSpec, conditions: Conditions) -> float:
        """Total (L-1) association penalty with the declared composition rule."""
        self._validate_spec(spec)
        reference = self._reference(conditions)
        if len(spec.strands) == 1:
            return 0.0
        n_rna = sum(len(s.sequence) for s in spec.strands if s.material == "rna")
        fraction = n_rna / spec.n_nucleotides
        majority = "rna" if fraction >= .5 else "dna"
        proportion = fraction if majority == "rna" else 1-fraction
        mono = self._mono[majority]
        with _NATIVE_LOCK:
            params = mono._parameters(reference)
            homogeneous = mono._energy_of_weight(params.expDuplexInit, params.kT/1000)
        hybrid = self._hs_energy("init", conditions)
        return (len(spec.strands)-1) * (hybrid + (2*proportion-1)*(homogeneous-hybrid))

    def _terminal_energy(self, spec: ComplexSpec, conditions: Conditions, pair: tuple[int, int]) -> float:
        i, j = pair
        materials = self._materials(spec)
        if materials[i] != materials[j]:
            return 0.0
        mono = self._mono[materials[i]]
        params = mono._parameters(self._reference(conditions))
        pair_type = _PAIR_CODE[spec.sequence[i] + spec.sequence[j]]
        return mono._energy_of_weight(self.RNA.exp_E_ExtLoop(pair_type, -1, -1, params), params.kT/1000)

    def _hybrid_stack_energy(self, spec: ComplexSpec, conditions: Conditions, loop: Loop) -> float:
        i, j = loop.closing_pair
        k, l = loop.child_pairs[0]
        materials, sequence = self._materials(spec), spec.sequence
        if not (materials[i] == materials[k] and materials[j] == materials[l] and materials[i] != materials[j]):
            raise ValueError("Hybrid stacks require two homogeneous, different backbones")
        upper, lower = ((sequence[i]+sequence[k], sequence[j]+sequence[l]) if materials[i] == "rna"
                        else (sequence[l]+sequence[j], sequence[k]+sequence[i]))
        return self._hs_energy(upper.replace("U", "T") + "/" + lower, conditions)

    def loop_energy(self, spec: ComplexSpec, conditions: Conditions, loop: Loop) -> float:
        """Evaluate one validated planar face; exclude association and symmetry."""
        energy = self._reference_loop_energy(spec, conditions, loop)
        # Every pair closes exactly one face, including a nicked exterior
        # face. Child stems and the unbounded root never duplicate this term.
        return energy + (self.pair_salt_energy(conditions) if loop.closing_pair is not None else 0.0)

    def _reference_loop_energy(self, spec: ComplexSpec, conditions: Conditions, loop: Loop) -> float:
        self._validate_spec(spec)
        reference = self._reference(conditions)
        if not validate_loop_face(spec, loop, self.pair_allowed):
            return math.inf
        pairs = loop.child_pairs + (() if loop.closing_pair is None else (loop.closing_pair,))
        with _NATIVE_LOCK:
            if len(loop.materials) == 1:
                material = next(iter(loop.materials))
                return self._mono[material].loop_energy(self._projection(spec, material), reference, loop)
            if loop.kind == "stack":
                return self._hybrid_stack_energy(spec, conditions, loop)
            if loop.kind == "interior":
                return math.fsum(m.loop_energy(self._projection(spec, material), reference, loop)
                                 for material, m in self._mono.items()) / 2
            params = {material: m._parameters(reference) for material, m in self._mono.items()}
            thermal = self.thermal_energy(conditions)
            e = lambda weight: ViennaMotifModel._energy_of_weight(weight, thermal)
            if loop.kind == "bulge":
                size = loop.unpaired_count
                initiation = math.fsum(e(self.RNA.doubleP_getitem(p.expbulge, size)) for p in params.values()) / 2
                if size == 1:
                    return initiation + self._hybrid_stack_energy(spec, conditions, loop)
                return initiation + math.fsum(self._terminal_energy(spec, conditions, pair) for pair in pairs)
            if loop.kind in ("exterior", "multiloop"):
                terminal = math.fsum(self._terminal_energy(spec, conditions, pair) for pair in pairs)
                if loop.kind == "exterior":
                    return terminal
                initiation = math.fsum(e(p.expMLclosing) for p in params.values()) / 2
                branch = math.fsum(e(self.RNA.doubleP_getitem(p.expMLintern, 1)) for p in params.values()) / 2
                unpaired = math.fsum(e(p.expMLbase) for p in params.values()) / 2 if loop.unpaired_count else 0.0
                return initiation + len(pairs)*branch + loop.unpaired_count*unpaired + terminal
            raise ValueError("Mixed hairpins/chimeric backbones are outside this model")

    def energy(self, spec: ComplexSpec, conditions: Conditions, state: PairingConstraints | str) -> float:
        self._validate_spec(spec)
        self._validate_conditions(conditions)
        lengths = tuple(len(s.sequence) for s in spec.strands)
        if isinstance(state, str):
            state = parse_structure(state, lengths)
        if not isinstance(state, PairingConstraints) or not state.complete or len(state.partners) != spec.n_nucleotides:
            raise ValueError("A complete matching of the supplied complex is required")
        if not is_connected(state, lengths):
            return math.inf
        if any(not self.pair_allowed(spec, i, j) for i, j in enumerate(state.partners) if j > i):
            return math.inf
        loops = decompose_structure(spec, format_structure(state, lengths))
        return math.fsum([self.association_energy(spec, conditions),
                          self.thermal_energy(conditions)*math.log(spec.symmetry_factor),
                          *(self.loop_energy(spec, conditions, loop) for loop in loops)])
