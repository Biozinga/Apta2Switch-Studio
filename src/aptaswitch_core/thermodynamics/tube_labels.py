"""Separate chemical species at equilibrium from their supplied origin labels."""
from __future__ import annotations

from collections import defaultdict
from itertools import product
import math
from numbers import Real

from .types import ComplexSpec, Strand


def merge_species(species, totals):
    """Sum identical chemistry/sequence inputs before enumerating complexes."""
    species, totals = tuple(species), tuple(totals)
    if (not species or len(species) != len(totals) or any(not isinstance(s, Strand) for s in species)
            or len({s.name for s in species}) != len(species)):
        raise ValueError("Each input species needs a distinct origin label and one concentration")
    groups = {}
    for strand, total in zip(species, totals):
        if isinstance(total, bool) or not isinstance(total, Real) or not math.isfinite(total) or total < 0:
            raise ValueError("Input concentrations must be finite and nonnegative")
        representative, values = groups.setdefault(strand.chemical_identity, (strand, []))
        values.append(float(total))
    return tuple(group[0] for group in groups.values()), tuple(math.fsum(group[1]) for group in groups.values())


def labelled_populations(tube, species, totals):
    """Distribute origin labels without altering any equilibrium constant.

    For every physical ordered complex, label assignments have the product of
    their pool fractions. Only cyclic rotations of a labelled order are merged;
    mirrors remain distinct. Repeated identical positions automatically produce
    the binomial/multinomial factors. No division by symmetry is added here.
    """
    species, totals = tuple(species), tuple(totals)
    merged, amounts = merge_species(species, totals)
    physical = {s.chemical_identity: float(t) for s, t in zip(tube.species, tube.total_concentrations_molar)}
    expected = {s.chemical_identity: t for s, t in zip(merged, amounts)}
    if physical.keys() != expected.keys() or any(not math.isclose(physical[k], v, rel_tol=1e-12)
                                                for k, v in expected.items()):
        raise ValueError("Origin-label totals do not match the solved physical tube")
    groups = defaultdict(list)
    for s, amount in zip(species, totals):
        groups[s.chemical_identity].append((s, float(amount)))
    choices = {key: [(s, t / physical[key] if physical[key] else 1 / len(rows)) for s, t in rows]
               for key, rows in groups.items()}
    values = defaultdict(list)
    for record, concentration in zip(tube.partitions, tube.concentrations_molar):
        for assignment in product(*(choices[s.chemical_identity] for s in record.complex_spec.strands)):
            labels = tuple(row[0] for row in assignment)
            cycle = min((labels[i:] + labels[:i] for i in range(len(labels))),
                        key=lambda order: tuple(s.name for s in order))
            values[ComplexSpec(cycle)].append(float(concentration) * math.prod(row[1] for row in assignment))
    return {spec: math.fsum(concentrations) for spec, concentrations in values.items()}
