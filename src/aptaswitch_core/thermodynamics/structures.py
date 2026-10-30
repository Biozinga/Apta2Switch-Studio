"""Structural indices shared by thermodynamic backends; no energy model here."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import accumulate


@dataclass(frozen=True)
class PairingConstraints:
    """Partial matching: -1 unresolved, own index unpaired, otherwise partner."""

    partners: tuple[int, ...]

    def __post_init__(self):
        if not isinstance(self.partners, tuple) or not self.partners:
            raise ValueError("A nonempty immutable tuple of pairing constraints is required")
        n = len(self.partners)
        for i, j in enumerate(self.partners):
            if type(j) is not int or not -1 <= j < n:
                raise ValueError("Constraint indices must be -1 or zero-based nucleotide indices")
            if j >= 0 and self.partners[j] != i:
                raise ValueError("Pairing constraints must be reciprocal")
        pairs = [(i, j) for i, j in enumerate(self.partners) if j > i]
        for a, (i, j) in enumerate(pairs):
            if any(i < k < j < l for k, l in pairs[a + 1:]):
                raise ValueError("Pseudoknots are outside the secondary-structure ensemble")

    @classmethod
    def unconstrained(cls, length: int) -> PairingConstraints:
        if type(length) is not int or length < 1:
            raise ValueError("Nucleotide count must be a positive integer")
        return cls((-1,) * length)

    @property
    def complete(self) -> bool:
        return -1 not in self.partners

    def assign(self, index: int, partner: int) -> PairingConstraints:
        if any(type(k) is not int or not 0 <= k < len(self.partners)
               for k in (index, partner)):
            raise ValueError("A pairing decision requires valid nucleotide indices")
        if self.partners[index] != -1 or self.partners[partner] != -1:
            raise ValueError("Cannot overwrite an existing pairing decision")
        result = list(self.partners)
        result[index], result[partner] = partner, index
        return PairingConstraints(tuple(result))


def _validated_lengths(lengths) -> tuple[int, ...]:
    lengths = tuple(lengths)
    if not lengths or any(type(n) is not int or n < 1 for n in lengths):
        raise ValueError("Each strand must have a positive integer length")
    return lengths


def parse_structure(structure: str, strand_lengths) -> PairingConstraints:
    """Read flat dot-bracket or correctly positioned '+'/'&' strand breaks."""
    lengths = _validated_lengths(strand_lengths)
    if not isinstance(structure, str):
        raise ValueError("Structure must be a dot-bracket string")
    separators = set(structure) & {"+", "&"}
    if len(separators) > 1:
        raise ValueError("Use one consistent strand separator")
    if separators:
        chunks = structure.split(next(iter(separators)))
        if tuple(map(len, chunks)) != lengths:
            raise ValueError("Structure strand breaks do not match the supplied strands")
        structure = "".join(chunks)
    if len(structure) != sum(lengths) or set(structure) - set(".()"):
        raise ValueError("Structure length or dot-bracket symbols are invalid")
    state = list(range(len(structure)))
    stack = []
    for i, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(i)
        elif symbol == ")":
            if not stack:
                raise ValueError("Unmatched closing base pair")
            j = stack.pop()
            state[i], state[j] = j, i
    if stack:
        raise ValueError("Unmatched opening base pair")
    return PairingConstraints(tuple(state))


def format_structure(state: PairingConstraints, strand_lengths) -> str:
    lengths = _validated_lengths(strand_lengths)
    if len(state.partners) != sum(lengths) or not state.complete:
        raise ValueError("A complete matching of the supplied strands is required")
    flat = "".join("." if j == i else "(" if j > i else ")"
                   for i, j in enumerate(state.partners))
    boundaries = (0,) + tuple(accumulate(lengths))
    return "+".join(flat[a:b] for a, b in zip(boundaries, boundaries[1:]))


def is_connected(state: PairingConstraints, strand_lengths) -> bool:
    """Connectivity of the interstrand pairing graph, ignoring unpaired ends."""
    lengths = _validated_lengths(strand_lengths)
    if len(state.partners) != sum(lengths) or not state.complete:
        raise ValueError("Connectivity requires a complete matching")
    ids = [strand for strand, n in enumerate(lengths) for _ in range(n)]
    neighbors = [set() for _ in lengths]
    for i, j in enumerate(state.partners):
        if j > i and ids[i] != ids[j]:
            neighbors[ids[i]].add(ids[j])
            neighbors[ids[j]].add(ids[i])
    visited, queue = {0}, [0]
    while queue:
        new = neighbors[queue.pop()] - visited
        visited.update(new)
        queue.extend(new)
    return len(visited) == len(lengths)
