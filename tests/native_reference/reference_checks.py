"""Small, independent combinatorial references; standard library only.

These helpers do not calculate molecular energies. They enumerate linear,
pseudoknot-free dot-bracket structures in one fixed concatenated strand order.
``cutpoints`` are zero-based boundaries: (2,) splits ``GGCC`` into ``GG&CC``.
No strand permutation, symmetry factor, association penalty or concentration
correction is applied. Callers must select their ensemble and supply energies
with those conventions explicitly. Enumeration grows exponentially.

Allowed pair types are AU, UA, CG, GC, GU and UG, including lonely pairs.
T is mapped to U ONLY for this combinatorial alphabet; this is not a DNA or
mixed RNA/DNA energy model. Comparisons must use noGU=False, noLP=False and
consistent loop constraints. Use dangles=0 when one structure means one state.

ViennaRNA's model documentation defines min_loop_size as the number of
unpaired nucleotides in a hairpin:
https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/grammar/model_details.html
The cut-loop exception is described in section 3.1.1 of:
https://www.tbi.univie.ac.at/newpapers/pdfs/TBI-t-2016-1.pdf
Interstrand pairs therefore have no minimum separation. An intrastrand pair
still requires the configured separation, even when other strand breaks exist.
"""
from __future__ import annotations

from collections.abc import Iterable, Iterator
import math
from numbers import Real


_PAIRS = frozenset(("AU", "UA", "CG", "GC", "GU", "UG"))


def _strand_ids(length: int, cutpoints: tuple[int, ...]) -> tuple[int, ...]:
    if (not isinstance(cutpoints, tuple)
            or any(type(point) is not int or not 0 < point < length for point in cutpoints)
            or tuple(sorted(set(cutpoints))) != cutpoints):
        raise ValueError("cutpoints must be strictly increasing integer boundaries inside the sequence")
    strand, result = 0, []
    for index in range(length):
        if strand < len(cutpoints) and index == cutpoints[strand]:
            strand += 1
        result.append(strand)
    return tuple(result)


def _partners(structure: str) -> tuple[int, ...]:
    if not isinstance(structure, str) or not structure:
        raise ValueError("structure must be a nonempty dot-bracket string without separators")
    partners = list(range(len(structure)))
    stack = []
    for index, character in enumerate(structure):
        if character == "(":
            stack.append(index)
        elif character == ")":
            if not stack:
                raise ValueError("structure contains an unmatched closing parenthesis")
            other = stack.pop()
            partners[index], partners[other] = other, index
        elif character != ".":
            raise ValueError("structure may contain only '.', '(' and ')'; no separators or pseudoknots")
    if stack:
        raise ValueError("structure contains an unmatched opening parenthesis")
    return tuple(partners)


def enum_structures(sequence: str, cutpoints: tuple[int, ...] = (),
                    min_loop_size: int = 3) -> Iterator[str]:
    """Enumerate each admissible dot-bracket string exactly once, including free strands.

    Inputs are nonempty A/C/G/U/T strings (lowercase is accepted), with separate
    cutpoints rather than '&'. Filter with ``connected`` for a connected complex.
    No energy cutoff, maximal internal-loop size or coaxial/dangle state is used.
    """
    if not isinstance(sequence, str) or not sequence or set(sequence.upper()) - set("ACGUT"):
        raise ValueError("sequence must contain only A, C, G, U or T and cannot be empty")
    if type(min_loop_size) is not int or min_loop_size < 0:
        raise ValueError("min_loop_size must be a nonnegative integer")
    normalized = sequence.upper().replace("T", "U")
    strands = _strand_ids(len(sequence), cutpoints)

    def interval(first: int, stop: int) -> Iterator[str]:
        if first == stop:
            yield ""
            return
        for suffix in interval(first + 1, stop):
            yield "." + suffix
        for partner in range(first + 1, stop):
            if normalized[first] + normalized[partner] not in _PAIRS:
                continue
            if strands[first] == strands[partner] and partner - first - 1 < min_loop_size:
                continue
            for inside in interval(first + 1, partner):
                for outside in interval(partner + 1, stop):
                    yield "(" + inside + ")" + outside

    return interval(0, len(sequence))


def connected(structure: str, cutpoints: tuple[int, ...]) -> bool:
    """Whether every strand belongs to one component of the interstrand-pair graph.

    Single strands are connected even if fully unpaired. Dot-bracket syntax and
    boundaries are checked; nucleotide compatibility requires the sequence and
    is deliberately outside this helper's contract.
    """
    partners = _partners(structure)
    strands = _strand_ids(len(structure), cutpoints)
    neighbors = [set() for _ in range(len(cutpoints) + 1)]
    for index, partner in enumerate(partners):
        first, second = strands[index], strands[partner]
        if first != second:
            neighbors[first].add(second)
    reached, pending = {0}, [0]
    while pending:
        for neighbor in neighbors[pending.pop()]:
            if neighbor not in reached:
                reached.add(neighbor)
                pending.append(neighbor)
    return len(reached) == len(neighbors)


def pair_matrix(structure: str) -> list[list[float]]:
    """Binary symmetric pairing matrix with an unpaired probability on the diagonal."""
    partners = _partners(structure)
    matrix = [[0.0] * len(partners) for _ in partners]
    for index, partner in enumerate(partners):
        matrix[index][partner] = 1.0
    return matrix


def _finite(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite real number")
    return number


def boltzmann_summary(structures: Iterable[str], energies: Iterable[float],
                      kT: float) -> dict:
    """Return -kT*log(Z) and pair/unpaired marginals using shifted Boltzmann weights.

    Energies and kT must share units; kT must be positive. Equal-energy DISTINCT
    structures retain their degeneracy. Duplicate dot-bracket strings are rejected
    to prevent accidental double counting: dangling-end microstates, strand
    permutations and symmetry factors must be handled explicitly by the caller.
    The returned free energy covers exactly the supplied structures, with no
    automatic connectedness filter or concentration correction.
    """
    temperature = _finite(kT, "kT")
    if temperature <= 0:
        raise ValueError("kT must be positive")
    states, values = list(structures), list(energies)
    if not states or len(states) != len(values):
        raise ValueError("structures and energies must be nonempty and have equal length")
    partners = [_partners(state) for state in states]
    size = len(partners[0])
    if any(len(state) != size for state in partners):
        raise ValueError("all structures must have the same nucleotide count")
    if len(set(states)) != len(states):
        raise ValueError("duplicate structures would count the same state more than once")
    values = [_finite(value, "energy") for value in values]
    minimum = min(values)
    weights = [math.exp(-(value - minimum) / temperature) for value in values]
    partition = math.fsum(weights)
    contributions: dict[tuple[int, int], list[float]] = {}
    for state, weight in zip(partners, weights):
        for index, partner in enumerate(state):
            contributions.setdefault((index, partner), []).append(weight)
    matrix = [[0.0] * size for _ in range(size)]
    for (index, partner), terms in contributions.items():
        matrix[index][partner] = math.fsum(terms) / partition
    return {"free_energy": minimum - temperature * math.log(partition),
            "pair_probabilities": matrix, "n_structures": len(states)}
