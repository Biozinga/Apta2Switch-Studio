"""Engine-independent sequence constraints for the current switch designs.

Only sequence admissibility is represented here: generating a sequence does
not optimize its ensemble, pass the MFE screening, or predict function. Target
pairs impose complementarity. OFF RNA pairs may be GU;
trigger recognition is strictly Watson--Crick, including inverse extensions.

Finite relational libraries handle multi-state couplings without treating GU
as a deterministic complement. Small coupled components are sampled uniformly
from their admissible assignments; generic larger components use bounded
backtracking. Neither proposal mechanism samples a thermodynamic ensemble.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from itertools import product
import math
from numbers import Real
import random
from types import MappingProxyType
from typing import Mapping

from ..architecture import architecture_layout, target_switch_structures, trigger_binding_slices
from ..localization import tr
from ..models import ReporterContext, RunRequest
from ..sequences import OPTIMIZABLE_BASE, STOP_CODONS, normalize_rna, reading_frame_report, reverse_complement_rna_pattern
from .structures import parse_structure
from .types import Conditions, Strand


_DESIGN_DNA = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
              "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
              "D": "AGT", "H": "ACT", "V": "ACG", OPTIMIZABLE_BASE: "ACGT"}


class UnsatisfiableDesignError(ValueError):
    """No assignment satisfies the supplied patterns and relational libraries."""


class NoMutableSitesError(ValueError):
    """No different admissible assignment exists in the requested mutation scope."""


def _name(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("A nonempty name is required")
    return value.strip()


def design_bases(symbol: str, material: str) -> tuple[str, ...]:
    """Resolve one design code (X is free); chemistry is explicit, with no T/U conversion."""
    if material not in {"rna", "dna"} or not isinstance(symbol, str) or len(symbol) != 1:
        raise ValueError("One design code and an explicit RNA/DNA material are required")
    code = symbol.upper()
    if (material == "rna" and code == "T") or (material == "dna" and code == "U"):
        raise ValueError("The nucleotide code belongs to the other strand chemistry")
    key = "T" if code == "U" else code
    if key not in _DESIGN_DNA:
        raise ValueError(f"Unsupported design code: {symbol!r}")
    bases = _DESIGN_DNA[key]
    return tuple(bases.replace("T", "U") if material == "rna" else bases)


@dataclass(frozen=True, slots=True)
class StrandPattern:
    name: str
    material: str
    pattern: str

    def __post_init__(self):
        object.__setattr__(self, "name", _name(self.name))
        if not isinstance(self.pattern, str) or not self.pattern:
            raise ValueError("A nonempty explicit design pattern is required")
        pattern = self.pattern.upper()
        for symbol in pattern:
            design_bases(symbol, self.material)
        object.__setattr__(self, "pattern", pattern)


@dataclass(frozen=True, slots=True, order=True)
class BaseSite:
    strand: str
    index: int

    def __post_init__(self):
        object.__setattr__(self, "strand", _name(self.strand))
        if type(self.index) is not int or self.index < 0:
            raise ValueError("Base indices must be nonnegative integers")


@dataclass(frozen=True, slots=True)
class LibraryConstraint:
    """Allowed concrete words in the exact order of ``sites`` (not reversed)."""
    name: str
    sites: tuple[BaseSite, ...]
    allowed_words: tuple[str, ...]

    def __post_init__(self):
        object.__setattr__(self, "name", _name(self.name))
        sites, words = tuple(self.sites), tuple(self.allowed_words)
        if (not sites or any(not isinstance(site, BaseSite) for site in sites)
                or len(set(sites)) != len(sites)):
            raise ValueError("Library sites must be distinct BaseSite objects")
        if not words or any(not isinstance(word, str) or len(word) != len(sites) for word in words):
            raise ValueError("Every library word must have exactly one base per site")
        object.__setattr__(self, "sites", sites)
        object.__setattr__(self, "allowed_words", tuple(sorted(set(words))))


@dataclass(frozen=True, slots=True)
class DesignTarget:
    name: str
    strands: tuple[str, ...]
    structure: str

    def __post_init__(self):
        object.__setattr__(self, "name", _name(self.name))
        names = tuple(_name(name) for name in self.strands)
        if not names:
            raise ValueError("A target needs at least one strand")
        object.__setattr__(self, "strands", names)


@dataclass(frozen=True, slots=True)
class DesignTube:
    name: str
    on_targets: tuple[tuple[str, float], ...]
    max_complex_size: int

    def __post_init__(self):
        object.__setattr__(self, "name", _name(self.name))
        rows = tuple((_name(name), value) for name, value in self.on_targets)
        if not rows or len({name for name, _ in rows}) != len(rows):
            raise ValueError("Tube on-target names must be nonempty and unique")
        if any(isinstance(value, bool) or not isinstance(value, Real)
               or not math.isfinite(value) or value <= 0 for _, value in rows):
            raise ValueError("Tube target concentrations must be finite positive molar values")
        if type(self.max_complex_size) is not int or self.max_complex_size < 1:
            raise ValueError("Tube maximum complex size must be a positive integer")
        object.__setattr__(self, "on_targets", tuple((name, float(value)) for name, value in rows))


@dataclass(frozen=True, slots=True)
class StrandRegion:
    name: str
    strand: str
    start: int
    stop: int

    def __post_init__(self):
        object.__setattr__(self, "name", _name(self.name))
        object.__setattr__(self, "strand", _name(self.strand))
        if any(type(value) is not int for value in (self.start, self.stop)) or not 0 <= self.start < self.stop:
            raise ValueError("A region needs valid nonempty half-open nucleotide offsets")


@dataclass(frozen=True, slots=True)
class ScreeningPolicy:
    """Sequence policy and the exact region used for the ON structural screen."""
    aug_start: int
    reporter_start: int
    exclude_stop_candidates: bool
    exclude_non_linear_candidates: bool
    on_unpaired_start: int
    on_unpaired_stop: int


@dataclass(frozen=True, slots=True)
class DesignAssignment:
    strands: tuple[Strand, ...]

    def __post_init__(self):
        strands = tuple(self.strands)
        if not strands or any(not isinstance(strand, Strand) for strand in strands):
            raise ValueError("An assignment requires concrete, explicitly typed strands")
        if len({strand.name for strand in strands}) != len(strands):
            raise ValueError("Assignment strand names must be unique")
        object.__setattr__(self, "strands", strands)

    def sequence(self, name: str) -> str:
        return next(strand.sequence for strand in self.strands if strand.name == name)


@dataclass(frozen=True, slots=True)
class _Component:
    sites: tuple[BaseSite, ...]
    libraries: tuple[LibraryConstraint, ...]


def _solve(domains: dict[BaseSite, tuple[str, ...]], libraries: tuple[LibraryConstraint, ...],
           rng: random.Random) -> dict[BaseSite, str] | None:
    """Finite backtracking with relational propagation; no retry-until-success."""
    while True:
        changed = False
        for library in libraries:
            words = [word for word in library.allowed_words
                     if all(base in domains[site] for site, base in zip(library.sites, word))]
            if not words:
                return None
            for column, site in enumerate(library.sites):
                allowed = {word[column] for word in words}
                updated = tuple(base for base in domains[site] if base in allowed)
                if updated != domains[site]:
                    domains[site], changed = updated, True
        if not changed:
            break
    if any(not values for values in domains.values()):
        return None
    unresolved = [site for site, values in domains.items() if len(values) > 1]
    if not unresolved:
        return {site: values[0] for site, values in domains.items()}
    site = min(unresolved, key=lambda key: (len(domains[key]), key))
    choices = list(domains[site])
    rng.shuffle(choices)
    for base in choices:
        result = _solve({**domains, site: (base,)}, libraries, rng)
        if result is not None:
            return result
    return None


@dataclass(frozen=True, slots=True)
class DesignProblem:
    strands: tuple[StrandPattern, ...]
    libraries: tuple[LibraryConstraint, ...] = ()
    targets: tuple[DesignTarget, ...] = ()
    tubes: tuple[DesignTube, ...] = ()
    regions: tuple[StrandRegion, ...] = ()
    layout: Mapping[str, int | str] = field(default_factory=dict)
    conditions: Conditions | None = None
    screening: ScreeningPolicy | None = None
    reporter: ReporterContext | None = None
    _domains: Mapping[BaseSite, tuple[str, ...]] = field(init=False, repr=False, compare=False)
    _components: tuple[_Component, ...] = field(init=False, repr=False, compare=False)
    _component_choices: tuple[tuple[tuple[str, ...], ...] | None, ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        for field_name, cls in (("strands", StrandPattern), ("libraries", LibraryConstraint),
                                ("targets", DesignTarget), ("tubes", DesignTube), ("regions", StrandRegion)):
            values = tuple(getattr(self, field_name))
            if any(not isinstance(value, cls) for value in values):
                raise ValueError(f"Invalid {field_name} entries")
            if len({value.name for value in values}) != len(values):
                raise ValueError(f"Duplicate names in {field_name}")
            object.__setattr__(self, field_name, values)
        if not self.strands:
            raise ValueError("A design requires at least one strand pattern")
        patterns = {strand.name: strand for strand in self.strands}
        domains = {BaseSite(strand.name, i): design_bases(code, strand.material)
                   for strand in self.strands for i, code in enumerate(strand.pattern)}
        neighbors = {site: set() for site in domains}
        for library in self.libraries:
            if any(site not in domains for site in library.sites):
                raise ValueError("A library references an absent nucleotide")
            for word in library.allowed_words:
                for site, base in zip(library.sites, word):
                    alphabet = "ACGU" if patterns[site.strand].material == "rna" else "ACGT"
                    if base not in alphabet:
                        raise ValueError("A library word has a nonconcrete base or incompatible chemistry")
            for site in library.sites:
                neighbors[site].update(library.sites)
        components, unseen = [], set(domains)
        # Pairing joins trigger, recognition arm and OFF partner. Sense-codon
        # constraints can join three such groups, rather than the whole strand.
        while unseen:
            pending, connected = [min(unseen)], set()
            while pending:
                site = pending.pop()
                if site not in connected:
                    connected.add(site)
                    pending.extend(neighbors[site] - connected)
            unseen.difference_update(connected)
            components.append(_Component(tuple(sorted(connected)), tuple(
                lib for lib in self.libraries if lib.sites[0] in connected)))
        for target in self.targets:
            if any(name not in patterns for name in target.strands):
                raise ValueError("A target references an absent strand")
            parse_structure(target.structure, [len(patterns[name].pattern) for name in target.strands])
        target_names = {target.name for target in self.targets}
        if any(name not in target_names for tube in self.tubes for name, _ in tube.on_targets):
            raise ValueError("A tube references an absent target")
        if any(region.strand not in patterns or region.stop > len(patterns[region.strand].pattern)
               for region in self.regions):
            raise ValueError("A region extends outside its strand")
        object.__setattr__(self, "layout", MappingProxyType(dict(self.layout)))
        object.__setattr__(self, "_domains", MappingProxyType(domains))
        object.__setattr__(self, "_components", tuple(components))
        choices = []
        for component in components:
            # Enumerating small components avoids the bias of taking the first
            # CSP solution; larger codon couplings retain finite backtracking.
            if math.prod(len(domains[s]) for s in component.sites) > 4096:
                choices.append(None)
                continue
            valid = []
            for word in product(*(domains[s] for s in component.sites)):
                bases = dict(zip(component.sites, word))
                if all("".join(bases[s] for s in lib.sites) in lib.allowed_words
                       for lib in component.libraries):
                    valid.append(word)
            choices.append(tuple(valid))
        object.__setattr__(self, "_component_choices", tuple(choices))

    def validate(self, assignment: DesignAssignment) -> None:
        if not isinstance(assignment, DesignAssignment):
            raise ValueError("An explicit DesignAssignment is required")
        if tuple(s.name for s in assignment.strands) != tuple(s.name for s in self.strands):
            raise ValueError("Assignment strand order/names differ from the design")
        bases = {}
        for strand, pattern in zip(assignment.strands, self.strands):
            if strand.material != pattern.material or len(strand.sequence) != len(pattern.pattern):
                raise ValueError("Assignment strand length or chemistry differs from its pattern")
            for i, base in enumerate(strand.sequence):
                site = BaseSite(strand.name, i)
                if base not in self._domains[site]:
                    raise ValueError("Assignment changed a fixed base or violates a design pattern")
                bases[site] = base
        if any("".join(bases[site] for site in lib.sites) not in lib.allowed_words for lib in self.libraries):
            raise ValueError("Assignment violates a coupled sequence library")

    def _assignment(self, bases: Mapping[BaseSite, str]) -> DesignAssignment:
        result = DesignAssignment(tuple(Strand(
            "".join(bases[BaseSite(s.name, i)] for i in range(len(s.pattern))), s.material, s.name)
            for s in self.strands))
        self.validate(result)
        return result

    def generate(self, rng: random.Random) -> DesignAssignment:
        if not isinstance(rng, random.Random):
            raise ValueError("Pass an explicit random.Random instance for reproducible proposals")
        bases = {}
        for component, choices in zip(self._components, self._component_choices):
            if choices is not None:
                if not choices:
                    raise UnsatisfiableDesignError(f"Incompatible constraints near {component.sites[0]}")
                bases.update(zip(component.sites, rng.choice(choices)))
                continue
            solution = _solve({site: self._domains[site] for site in component.sites}, component.libraries, rng)
            if solution is None:
                raise UnsatisfiableDesignError(f"Incompatible constraints near {component.sites[0]}")
            bases.update(solution)
        return self._assignment(bases)

    def mutable_components(self) -> tuple[tuple[BaseSite, ...], ...]:
        """Independent coupled groups; fixed groups never inflate diversity."""
        return tuple(c.sites for c, choices in zip(self._components, self._component_choices)
                     if (len(choices) > 1 if choices is not None
                         else any(len(self._domains[s]) > 1 for s in c.sites)))

    def enumerate_assignments(self, limit: int) -> tuple[DesignAssignment, ...] | None:
        """Enumerate only when the entire constrained space fits the budget."""
        if type(limit) is not int or limit < 1:
            raise ValueError("An enumeration limit must be a positive integer")
        if any(c is None for c in self._component_choices):
            return None
        if any(not c for c in self._component_choices):
            raise UnsatisfiableDesignError("Incompatible coupled sequence constraints")
        if math.prod(len(c) for c in self._component_choices) > limit:
            return None
        return tuple(self._assignment({s: b for c, word in zip(self._components, words)
                                       for s, b in zip(c.sites, word)})
                     for words in product(*self._component_choices))

    def design_distance(self, first: DesignAssignment, second: DesignAssignment) -> float:
        """Hamming distance over mutable components, not fixed sequence length."""
        self.validate(first)
        self.validate(second)
        a, b = {s.name: s.sequence for s in first.strands}, {s.name: s.sequence for s in second.strands}
        groups = self.mutable_components()
        return (sum(any(a[s.strand][s.index] != b[s.strand][s.index] for s in group)
                    for group in groups) / len(groups)) if groups else 0.0

    def mutate_guided(self, assignment: DesignAssignment, rng: random.Random,
                      site_weights, *, exploration_fraction: float = .2,
                      multi_fraction: float = .15) -> DesignAssignment:
        """Resample one or several complete groups without breaking OFF/ON pairs."""
        self.validate(assignment)
        weights = dict(site_weights)
        if any(not math.isfinite(v) or v < 0 for v in weights.values()):
            raise ValueError("Mutation weights must be finite and nonnegative")
        if not 0 <= exploration_fraction <= 1 or not 0 <= multi_fraction <= 1:
            raise ValueError("Mutation fractions must lie between zero and one")
        groups = list(self.mutable_components())
        if not groups:
            raise NoMutableSitesError("No mutable coupled component")
        count = min(len(groups), rng.choice((2, 3))) if rng.random() < multi_fraction else 1
        current = assignment
        changed = False
        for _ in range(count):
            scores = [math.fsum(weights.get(s, 0.0) for s in group) for group in groups]
            if rng.random() < exploration_fraction or not any(scores):
                group = rng.choice(groups)
            else:
                group = rng.choices(groups, weights=scores, k=1)[0]
            groups.remove(group)
            index = next(i for i, c in enumerate(self._components) if c.sites == group)
            choices = self._component_choices[index]
            if choices is not None:
                bases = {BaseSite(s.name, i): b for s in current.strands for i, b in enumerate(s.sequence)}
                word = tuple(bases[s] for s in group)
                alternatives = [value for value in choices if value != word]
                if alternatives:
                    bases.update(zip(group, rng.choice(alternatives)))
                    current, changed = self._assignment(bases), True
            else:
                candidates = list(group)
                rng.shuffle(candidates)
                for site in candidates:
                    try:
                        current = self.mutate(current, rng, site=site)
                        changed = True
                        break
                    except NoMutableSitesError:
                        continue
        if not changed:
            raise NoMutableSitesError("No different assignment in the selected groups")
        return current

    def mutate(self, assignment: DesignAssignment, rng: random.Random, *, site: BaseSite | None = None) -> DesignAssignment:
        """Change one feasible site and repair its entire coupled component.

        Other components and all fixed bases stay unchanged. Exhausting the
        finite candidate sites raises NoMutableSitesError instead of silently
        returning the input or retrying indefinitely.
        """
        self.validate(assignment)
        if not isinstance(rng, random.Random):
            raise ValueError("Pass an explicit random.Random instance for reproducible proposals")
        if site is not None and (not isinstance(site, BaseSite) or site not in self._domains):
            raise ValueError("The requested mutation site is absent from the design")
        bases = {BaseSite(s.name, i): base for s in assignment.strands for i, base in enumerate(s.sequence)}
        candidates = [site] if site is not None else [s for s in self._domains if len(self._domains[s]) > 1]
        rng.shuffle(candidates)
        for selected in candidates:
            alternatives = tuple(base for base in self._domains[selected] if base != bases[selected])
            if not alternatives:
                continue
            component = next(c for c in self._components if selected in c.sites)
            domains = {s: self._domains[s] for s in component.sites}
            domains[selected] = alternatives
            solution = _solve(domains, component.libraries, rng)
            if solution is not None:
                return self._assignment({**bases, **solution})
        raise NoMutableSitesError("No different assignment satisfies the constraints in this mutation scope")

    def reading_frame_report(self, assignment: DesignAssignment) -> dict:
        self.validate(assignment)
        if self.screening is None:
            raise ValueError("This design has no switch reading-frame policy")
        return reading_frame_report(assignment.sequence("Switch"), self.screening.aug_start,
                                    reporter_start=self.screening.reporter_start)

    def passes_sequence_screen(self, assignment: DesignAssignment) -> bool:
        """Defensive STOP check; MFE screening still needs an engine."""
        report = self.reading_frame_report(assignment)
        return not (self.screening.exclude_stop_candidates and report["has_premature_stop"])


def _pair_library(name: str, left: BaseSite, right: BaseSite,
                  patterns: Mapping[str, StrandPattern], *, wobble: bool) -> LibraryConstraint:
    lm, rm = patterns[left.strand].material, patterns[right.strand].material
    lu, ru = ("U" if lm == "rna" else "T"), ("U" if rm == "rna" else "T")
    words = ("A" + ru, lu + "A", "CG", "GC")
    if wobble and lm == rm == "rna":
        words += ("GU", "UG")
    return LibraryConstraint(name, (left, right), words)


def _with_stop_free_codons(problem: DesignProblem) -> DesignProblem:
    """Constrain only coding triplets that could become STOPs after coupling.

    Pairing libraries narrow the available bases before codons are joined into
    components. Fixed or effectively fixed sense codons therefore add no work.
    Codons span domain boundaries; the reporter's own codons remain unchanged.
    """
    policy = problem.screening
    if policy is None or not policy.exclude_stop_candidates:
        return problem
    domains = dict(problem._domains)
    for component, choices in zip(problem._components, problem._component_choices):
        if choices is not None:
            if not choices:
                raise UnsatisfiableDesignError(tr("Les bases fixées sont incompatibles avec les contraintes d’appariement."))
            for column, site in enumerate(component.sites):
                domains[site] = tuple(sorted({word[column] for word in choices}))
    libraries = []
    for start in range(policy.aug_start, policy.reporter_start, 3):
        sites = tuple(BaseSite("Switch", start + offset) for offset in range(3))
        if not any(all(base in domains[site] for site, base in zip(sites, stop)) for stop in STOP_CODONS):
            continue
        words = tuple("".join(word) for word in product(*(domains[site] for site in sites))
                      if "".join(word) not in STOP_CODONS)
        if not words:
            raise UnsatisfiableDesignError(tr(
                "Les bases fixées et les contraintes d’appariement imposent un STOP prématuré au codon {codon}.",
                codon=(start-policy.aug_start)//3+1))
        libraries.append(LibraryConstraint(f"sense_codon_{start}", sites, words))
    if not libraries:
        return problem
    constrained = replace(problem, libraries=problem.libraries + tuple(libraries))
    # Prove feasibility once, before any trial or expensive thermodynamics.
    # Application codon components are small; generic larger components use
    # the same complete finite constraint solver as sequence generation.
    for component, choices in zip(constrained._components, constrained._component_choices):
        if choices == () or (choices is None and _solve(
                {site: constrained._domains[site] for site in component.sites},
                component.libraries, random.Random(0)) is None):
            raise UnsatisfiableDesignError(tr(
                "Les bases fixées et les contraintes d’appariement ne permettent aucune séquence sans STOP prématuré avant le rapporteur."))
    return constrained


def build_design_problem(request: RunRequest) -> DesignProblem:
    """Snapshot the current request without loading a thermodynamic runtime."""
    request = request.normalized()
    arch, thermo = request.architecture, request.thermo
    layout = architecture_layout(arch)
    trigger = request.design_trigger_pattern()
    binding = reverse_complement_rna_pattern(trigger)
    slices = trigger_binding_slices(arch)

    def left(name, size):
        return binding[slices[name]] if name in slices else OPTIMIZABLE_BASE * size

    reporter = normalize_rna(request.reporter.sequence_after_start_rna, field_name="reporter")
    segments = [("t7_leader", arch.t7_leader), ("toehold", left("toehold", arch.toehold_length)),
                ("stem_l", left("stem_l", arch.stem_length)),
                ("upper2_l", left("upper2_l", arch.upper_stem2_length)),
                ("bulge_l", left("bulge_l", arch.bulge_length)),
                ("upper_l", left("upper_l", arch.upper_stem_length)),
                ("rbs_prefix", OPTIMIZABLE_BASE * arch.rbs_prefix_length), ("rbs", arch.rbs_loop),
                ("upper_r", OPTIMIZABLE_BASE * arch.upper_stem_length), ("aug_spacer", OPTIMIZABLE_BASE * arch.aug_spacer),
                ("aug", "AUG"), ("upper2_r", OPTIMIZABLE_BASE * arch.upper_stem2_length),
                ("stem_r", OPTIMIZABLE_BASE * arch.stem_length), ("frame_linker", layout["frame_linker"]),
                ("reporter_after_start", reporter)]
    regions, sequence = [], ""
    for name, pattern in segments:
        if pattern:
            pattern = "".join(pattern.split()).upper().replace("T", "U")
            regions.append(StrandRegion(name, "Switch", len(sequence), len(sequence) + len(pattern)))
            sequence += pattern
    patterns = (StrandPattern("Trigger", "dna", trigger),
                StrandPattern("Switch", "rna", request.fixed_switch_sequence or sequence))
    by_name = {s.name: s for s in patterns}
    off, on = target_switch_structures(arch, len(reporter))
    targets = (DesignTarget("switch_off", ("Switch",), off),
               DesignTarget("switch_on", ("Trigger", "Switch"), on),
               DesignTarget("trigger_free", ("Trigger",), "." * len(trigger)))
    libraries = []
    for target in (() if request.fixed_switch_sequence is not None else targets[:2]):
        coordinates = tuple(BaseSite(name, i) for name in target.strands
                            for i in range(len(by_name[name].pattern)))
        pairs = parse_structure(target.structure, [len(by_name[name].pattern) for name in target.strands])
        for i, j in enumerate(pairs.partners):
            if i < j:
                libraries.append(_pair_library(f"{target.name}_{i}_{j}", coordinates[i], coordinates[j],
                                               by_name, wobble=target.name == "switch_off"))
    # Assert the shared layout still describes the same domain construction.
    region_map = {region.name: region for region in regions}
    if (region_map["aug"].start != layout["aug_start"]
            or region_map["rbs"].start != layout["rbs_start"]
            or region_map["reporter_after_start"].start != layout["reporter_start"]):
        raise ValueError("Design sequence domains disagree with the shared architecture offsets")
    t, s = thermo.resolved_trigger_concentration_m, thermo.resolved_switch_concentration_m
    duplex = min(t, s)
    on_targets = [("switch_on", duplex)]
    if t > duplex:
        on_targets.append(("trigger_free", t - duplex))
    if s > duplex:
        on_targets.append(("switch_off", s - duplex))
    problem = DesignProblem(patterns, tuple(libraries), targets,
        (DesignTube("tube_off", (("switch_off", s),), thermo.screening_max_complex_size),
         DesignTube("tube_on", tuple(on_targets), thermo.screening_max_complex_size)),
        tuple(regions), layout, Conditions(thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m),
        ScreeningPolicy(layout["aug_start"], layout["reporter_start"],
                        request.exclude_stop_candidates, request.exclude_non_linear_candidates,
                        layout["rbs_linker_start"], layout["rbs_linker_end"]), request.reporter)
    return problem if request.fixed_switch_sequence is not None else _with_stop_free_codons(problem)
