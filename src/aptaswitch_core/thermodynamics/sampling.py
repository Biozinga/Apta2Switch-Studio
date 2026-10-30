"""Seeded Boltzmann sampling from an oracle of conditional pair probabilities.

The chain rule, not independent draws from unconditional marginals, defines
the joint law. Every oracle call must describe the same connected ensemble,
conditioned on the complete history of resolved nucleotide states.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import math
import random
from typing import Callable

from .structures import PairingConstraints, format_structure, is_connected
from .types import ComplexSpec, Conditions, PairProbabilities


class ProbabilityConsistencyError(RuntimeError):
    """A probability oracle violated its scientific contract."""


@dataclass(frozen=True)
class SamplingBranch:
    index: int
    choices: tuple[tuple[int, float], ...]


class ConditionalSampler:
    """One sampler per ordered complex, conditions and immutable energy model.

    Cached values contain only the next O(n) conditional distribution; they
    never retain native fold compounds or full O(n²) matrices. The callback
    is responsible for applying every hard constraint to its partition.
    """

    def __init__(self, complex_spec: ComplexSpec, conditions: Conditions,
                 oracle: Callable[[PairingConstraints], PairProbabilities], *,
                 seed: int, cache_size: int = 4096, tolerance: float = 1e-8):
        if type(seed) is not int:
            raise ValueError("An integer seed is required for reproducible sampling")
        if type(cache_size) is not int or cache_size < 0:
            raise ValueError("Cache size must be a nonnegative integer")
        if not math.isfinite(tolerance) or not 0 < tolerance <= 1e-6:
            raise ValueError("Probability tolerance must be in (0, 1e-6]")
        self.complex_spec, self.conditions, self.oracle = complex_spec, conditions, oracle
        self.lengths = tuple(len(s.sequence) for s in complex_spec.strands)
        self.length = sum(self.lengths)
        self.tolerance, self.cache_size = tolerance, cache_size
        self.random = random.Random(seed)
        self.cache: OrderedDict[PairingConstraints, SamplingBranch] = OrderedDict()
        self.oracle_calls = self.cache_hits = 0
        self.maximum_roundoff_adjustment = 0.0
        self._model_signature = None

    def conditional(self, state: PairingConstraints) -> SamplingBranch:
        if len(state.partners) != self.length or state.complete:
            raise ValueError("A partial matching of this complex is required")
        if state in self.cache:
            self.cache_hits += 1
            self.cache.move_to_end(state)
            return self.cache[state]
        result = self.oracle(state)
        self.oracle_calls += 1
        context = result.partition
        if (context.complex_spec != self.complex_spec or context.conditions != self.conditions):
            raise ProbabilityConsistencyError("The probability oracle changed strands, their order or conditions")
        expected_constraints = None if all(j == -1 for j in state.partners) else state
        if context.constraints != expected_constraints:
            raise ProbabilityConsistencyError("The probability oracle returned a different constrained ensemble")
        signature = (context.provenance, context.thermal_energy_kcal_mol)
        if self._model_signature is not None and signature != self._model_signature:
            raise ProbabilityConsistencyError("The probability oracle changed its energy model")
        self._model_signature = signature
        for i, partner in enumerate(state.partners):
            if partner >= 0 and abs(result.matrix[i][partner] - 1) > self.tolerance:
                raise ProbabilityConsistencyError("A previously selected pairing constraint was not enforced")
        index = state.partners.index(-1)
        row = result.matrix[index]
        forbidden = math.fsum(p for j, p in enumerate(row) if state.partners[j] != -1)
        choices = [(j, p) for j, p in enumerate(row) if state.partners[j] == -1 and p > 0]
        mass = math.fsum(p for _, p in choices)
        if forbidden > self.tolerance or not choices or abs(mass - 1) > self.tolerance:
            raise ProbabilityConsistencyError("Conditional probability mass is inconsistent with the matching")
        self.maximum_roundoff_adjustment = max(self.maximum_roundoff_adjustment, forbidden, abs(mass - 1))
        branch = SamplingBranch(index, tuple((j, p / mass) for j, p in choices))
        if self.cache_size:
            self.cache[state] = branch
            if len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)
        return branch

    def sample_one(self) -> str:
        state = PairingConstraints.unconstrained(self.length)
        while not state.complete:
            branch = self.conditional(state)
            cutoff, cumulative = self.random.random(), 0.0
            chosen = branch.choices[-1][0]
            for partner, probability in branch.choices:
                cumulative += probability
                if cutoff < cumulative:
                    chosen = partner
                    break
            try:
                state = state.assign(branch.index, chosen)
            except ValueError as error:
                raise ProbabilityConsistencyError("The probability oracle proposed an inadmissible pair") from error
        if not is_connected(state, self.lengths):
            raise ProbabilityConsistencyError("The probability oracle sampled a disconnected complex")
        return format_structure(state, self.lengths)

    def sample(self, count: int) -> list[str]:
        if type(count) is not int or count < 0:
            raise ValueError("Sample count must be a nonnegative integer")
        return [self.sample_one() for _ in range(count)]
