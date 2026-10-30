"""Automatic oxNA box sizing and periodic-image geometry checks.

The initial box uses the *nominal* fully extended backbone, not the extreme
FENE extension limit. The accompanying geometry check must therefore remain
active during sampling. No box-size cap is imposed for long strands.
"""

from __future__ import annotations

import math
from collections.abc import Sequence


OXDNA_LENGTH_NM = .8518
_AVOGADRO = 6.02214076e23
# oxDNA src/model.h and src/Interactions/rna_model.h, as used by NA.
_BACKBONE_SPACING = {"dna": .7564, "rna": .761070781051}
_BACKBONE_OFFSET = {"dna": math.hypot(.3400, .3408), "rna": math.hypot(.4, .2)}
# DNA2Interaction.cpp and RNAInteraction2.cpp; T is converted to T(K)/300.
_DEBYE_LAMBDA = {"dna": .3616455, "rna": .3667258}


def _finite(value: float, name: str, minimum: float, *, inclusive: bool = False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number.")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number.") from exc
    if not math.isfinite(result) or (result < minimum if inclusive else result <= minimum):
        relation = "at least" if inclusive else "greater than"
        raise ValueError(f"{name} must be finite and {relation} {minimum}.")
    return result


def _length(value: int, name: str) -> int:
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive nucleotide count.")
    return value


def estimate_ffs_box(
    switch_length: int,
    trigger_length: int,
    switch_material: str,
    trigger_material: str,
    temperature_c: float,
    salt_m: float,
    initial_separation: float = 10.,
) -> dict:
    """Size a cubic box from the longer nominal strand plus a salt-aware margin.

    All returned distances except ``box_size_nm`` use oxDNA length units.
    ``periodic_margin`` is the extra clearance *beyond* the interaction cutoff.
    Preparation must additionally check that its actual two-strand placement
    fits: no assumed sum of both full contour lengths is imposed here.
    """
    lengths = {"switch": _length(switch_length, "switch_length"),
               "trigger": _length(trigger_length, "trigger_length")}
    materials = {"switch": switch_material, "trigger": trigger_material}
    if any(material not in _BACKBONE_SPACING for material in materials.values()):
        raise ValueError("The oxNA box requires explicit dna or rna strand materials.")
    temperature = _finite(temperature_c, "temperature_c", -273.15)
    salt = _finite(salt_m, "salt_m", 0.)
    separation = _finite(initial_separation, "initial_separation", 0.)
    contour = {name: (length - 1) * _BACKBONE_SPACING[materials[name]]
               for name, length in lengths.items()}
    diameters = {name: contour[name] + 2 * _BACKBONE_OFFSET[material]
                 for name, material in materials.items()}

    # NA dispatches hybrid electrostatics to RNA2. Taking independently the
    # largest backbone offset and Debye length covers DNA, RNA and mixed pairs.
    radius = max(_BACKBONE_OFFSET[material] for material in materials.values())
    debye = max(_DEBYE_LAMBDA[material] for material in materials.values())
    debye *= math.sqrt((temperature + 273.15) / 300.) / math.sqrt(salt)
    # RHIGH = 3 lambda; the smoothed Debye-Huckel tail reaches zero at 4.5 lambda.
    # Keep the non-electrostatic cutoff even if temperature/salt makes DH tiny.
    cutoff = max(2 * radius + 4.5 * debye, 2 * radius + .711879214356, .8 + .783775)
    margin = 2.
    clearance = cutoff + margin
    box = float(math.ceil(max(max(diameters.values()) + clearance,
                              2 * clearance, 2 * separation + margin)))
    return {
        "sizing_method": "nominal_unfolded_contour_with_periodic_clearance_check",
        "box_size": box,
        "box_size_nm": box * OXDNA_LENGTH_NM,
        "concentration_each_m": 1 / (_AVOGADRO * (box * OXDNA_LENGTH_NM * 1e-9) ** 3 * 1000),
        "contour_lengths": contour,
        "nominal_diameters": diameters,
        "interaction_cutoff": cutoff,
        "periodic_margin": margin,
        "minimum_periodic_clearance": clearance,
        "initial_separation": separation,
        "geometry_monitor_required": True,
    }


def check_periodic_self_clearance(
    positions: Sequence[Sequence[float]],
    strand_lengths: Sequence[int],
    box_sides: Sequence[float],
    cutoff: float,
    safety_margin: float = 1.,
) -> dict:
    """Reject a strand approaching its own periodic image.

    First unwrap each covalent chain using consecutive minimum-image steps.
    For every axis, box side minus the strand's unwrapped span bounds its
    clearance to its closest copy. This sufficient condition is conservative
    for tilted strands and does not confuse an ordinary boundary crossing
    with a molecule filling the box. Separate chains are checked separately.
    """
    lengths = [_length(value, "strand length") for value in strand_lengths]
    if not lengths or sum(lengths) != len(positions):
        raise ValueError("Strand lengths must cover exactly all nucleotide positions.")
    if len(box_sides) != 3:
        raise ValueError("The periodic box must have three positive sides.")
    box = [_finite(value, "box side", 0.) for value in box_sides]
    distance = _finite(cutoff, "interaction cutoff", 0.)
    margin = _finite(safety_margin, "safety margin", 0., inclusive=True)
    required = distance + margin
    rows = []
    for position in positions:
        if len(position) != 3:
            raise ValueError("Each nucleotide position must have three coordinates.")
        try:
            row = [float(value) for value in position]
        except (TypeError, ValueError) as exc:
            raise ValueError("Nucleotide coordinates must be finite numbers.") from exc
        if not all(math.isfinite(value) for value in row):
            raise ValueError("Nucleotide coordinates must be finite numbers.")
        rows.append(row)

    strand_spans, strand_clearances = [], []
    start = 0
    for strand_index, length in enumerate(lengths, 1):
        strand = rows[start:start + length]
        unwrapped = [strand[0][:]]
        for previous, current in zip(strand, strand[1:]):
            unwrapped.append([
                unwrapped[-1][axis] + (current[axis] - previous[axis] + box[axis] / 2) % box[axis]
                - box[axis] / 2 for axis in range(3)
            ])
        spans = [max(row[axis] for row in unwrapped) - min(row[axis] for row in unwrapped)
                 for axis in range(3)]
        clearances = [box[axis] - spans[axis] for axis in range(3)]
        for axis, clearance in enumerate(clearances):
            if clearance + 1e-9 < required:
                raise ValueError(
                    f"The oxDNA periodic box is too small for strand {strand_index}: "
                    f"clearance to its periodic image on {'xyz'[axis]} is {clearance:.3f} "
                    f"but at least {required:.3f} oxDNA units are required. "
                    "This configuration must not be retained for FFS."
                )
        strand_spans.append(spans)
        strand_clearances.append(clearances)
        start += length
    return {"minimum_clearance": min(min(values) for values in strand_clearances),
            "required_clearance": required, "strand_spans": strand_spans,
            "strand_clearances": strand_clearances}
