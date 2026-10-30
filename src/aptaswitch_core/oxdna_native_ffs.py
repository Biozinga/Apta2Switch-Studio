"""Documented CPU ``FFS_MD`` execution for sequence-dependent oxNA.

oxDNA evaluates the order parameters and AND/OR stop conditions at *every*
integration step. Python only schedules shots and retrieves progress in chunks;
a chunk boundary is never interpreted as a physical failure. Configuration
exports retain orientations, velocities, angular momenta and the native clock.

References: https://lorenzo-rovigatti.github.io/oxDNA/ffs.html and
https://lorenzo-rovigatti.github.io/oxDNA/input.html#forward-flux-sampling-ffs-options
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Sequence

try:
    from .oxdna_ffs_worker import Observation, OxpyAdapter, _pairs, distance_pair_groups
except ImportError:  # External oxpy interpreter executes the worker directly.
    from oxdna_ffs_worker import Observation, OxpyAdapter, _pairs, distance_pair_groups


_EXPRESSION = re.compile(r"^(binding|stem|opening|associations|distance|separation)\s*(<=|>=|==|<|>)\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$")
_STOP = re.compile(r"Reached condition(\d+)\b")
_COMPARE = {"<=": lambda a, b: a <= b, ">=": lambda a, b: a >= b,
            "==": lambda a, b: a == b, "<": lambda a, b: a < b,
            ">": lambda a, b: a > b}


def _condition(expression: str) -> tuple[str, str, float]:
    match = _EXPRESSION.fullmatch(expression.strip())
    if match is None:
        raise ValueError(f"Invalid native FFS expression: {expression!r}.")
    name, operator, raw = match.groups()
    value = float(raw)
    if not math.isfinite(value) or value < 0 or (name not in {"distance", "separation"} and not value.is_integer()):
        raise ValueError("FFS thresholds must be nonnegative; bond counts must be integers.")
    return name, operator, value


def conditions_text(conditions: Sequence[Sequence[str]]) -> str:
    """Each condition is an AND block; oxDNA stops at the first true block."""
    blocks = []
    for number, expressions in enumerate(conditions, 1):
        if isinstance(expressions, str) or not expressions:
            raise ValueError("Each FFS stop condition needs a nonempty list of expressions.")
        lines = []
        for expression in expressions:
            name, operator, value = _condition(expression)
            lines.append(f"    {name} {operator} {value:.17g}")
        blocks.append(f"condition{number} = {{\n" + "\n".join(lines) + "\n}\n")
    return "\n".join(blocks)


def first_condition(observation: Observation, conditions: Sequence[Sequence[str]]) -> int | None:
    """Classify a phase point already at a boundary without advancing past it.

    Only the initial state uses this guard. Subsequent crossings are identified
    by the native backend's log, including crossings inside a requested chunk.
    """
    for index, expressions in enumerate(conditions):
        if not expressions:
            raise ValueError("Empty FFS stop condition.")
        matches = True
        for expression in expressions:
            name, operator, threshold = _condition(expression)
            value = getattr(observation, name)
            if value is None:
                raise ValueError(f"Native observation is missing {name}.")
            matches = matches and _COMPARE[operator](value, threshold)
        if matches:
            return index
    return None


def order_parameters_text(plan: dict) -> str:
    """Create native bond counts and explicit nucleotide-centre distances."""
    n, total = plan["switch_length"], plan["nucleotide_count"]
    if type(n) is not int or type(total) is not int or not 0 < n < total:
        raise ValueError("Native FFS requires two explicit nonempty strands.")
    cutoff = float(plan.get("hb_threshold", -.1))
    if cutoff != -.1:
        # HBList and native FFS share this documented cutoff. Other values
        # require a different observation reader, not silently mixed criteria.
        raise ValueError("Native FFS currently requires the oxDNA HB cutoff -0.1.")
    groups = {name: _pairs(plan.get(f"{name}_pairs", []), name)
              for name in ("binding", "stem", "opening")}
    if not groups["binding"]:
        raise ValueError("Native FFS needs target binding pairs.")
    if any(max(pair) >= total for pairs in groups.values() for pair in pairs):
        raise ValueError("Native FFS pair index exceeds the topology.")
    groups["associations"] = tuple((i, j) for i in range(n) for j in range(n, total))
    blocks = []
    for name, pairs in groups.items():
        lines = ["{", "    order_parameter = bond", f"    name = {name}", "    cutoff = -0.1"]
        # oxDNA requires at least one pair key. A self-pair never appears in
        # its nonbonded neighbor list, so this represents the constant zero
        # count exactly when a monitored region contains no candidate pairs.
        lines += [f"    pair{k} = {i}, {j}" for k, (i, j) in enumerate(pairs or ((0, 0),), 1)]
        blocks.append("\n".join([*lines, "}"]))
    for name, pairs in distance_pair_groups(plan).items():
        lines = ["{", "    order_parameter = mindistance", f"    name = {name}", "    use_COM = true"]
        lines += [f"    pair{k} = {i}, {j}" for k, (i, j) in enumerate(pairs, 1)]
        blocks.append("\n".join([*lines, "}"]))
    return "\n\n".join(blocks) + "\n"


@dataclass(frozen=True)
class NativeChunk:
    steps: int
    condition: int | None
    observation: Observation


class NativeFFSAdapter(OxpyAdapter):
    """One real CPU FFS_MD trajectory, with native first-crossing stops.

    Managers cannot coexist inside one oxpy process. Parallel trajectories
    therefore run in separate worker processes. ``initial=False`` preserves
    momenta for shooting and resume; a fresh seed supplies stochastic dynamics.
    """

    def __init__(self, request: dict, configuration: Path, seed: int, initial: bool,
                 *, directory: Path, conditions: Sequence[Sequence[str]]):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.conditions = tuple(tuple(expressions) for expressions in conditions)
        self.op_file = self.directory / "order_parameters.txt"
        self.ffs_file = self.directory / "conditions.txt"
        self.log_file = self.directory / "native_ffs.log"
        self.op_file.write_text(order_parameters_text(request["plan"]), encoding="utf-8")
        self.ffs_file.write_text(conditions_text(self.conditions), encoding="utf-8")
        self.log_file.write_text("", encoding="utf-8")
        self._stopped = None
        self._log_offset = 0
        super().__init__(request, configuration, seed, initial)
        self._log_offset = self.log_file.stat().st_size
        self._initial_checked = False

    def _configure_input(self, values: dict[str, str]) -> None:
        values.update(sim_type="FFS_MD", order_parameters_file=str(self.op_file),
                      ffs_file=str(self.ffs_file), log_file=str(self.log_file))
        # Keep a readable record of exactly the input supplied to oxDNA.
        (self.directory / "native_input.txt").write_text(
            "\n".join(f"{key} = {value}" for key, value in values.items()) + "\n",
            encoding="utf-8")

    def advance(self, steps: int = 1000) -> NativeChunk:
        if type(steps) is not int or steps < 1:
            raise ValueError("The native reporting chunk must be a positive integer.")
        if self._stopped is not None:
            return NativeChunk(0, self._stopped, self.observe())
        if not self._initial_checked:
            self._initial_checked = True
            observation = self.observe()
            self._stopped = first_condition(observation, self.conditions)
            if self._stopped is not None:
                return NativeChunk(0, self._stopped, observation)
        before = self.manager.current_step
        self.manager.run(steps, print_output=False)
        completed = self.manager.current_step - before
        with self.log_file.open(encoding="utf-8") as handle:
            handle.seek(self._log_offset)
            new_log = handle.read()
            self._log_offset = handle.tell()
        stop = _STOP.search(new_log)
        if stop:
            self._stopped = int(stop.group(1)) - 1
            if not 0 <= self._stopped < len(self.conditions):
                raise RuntimeError("oxDNA reported an unknown FFS condition.")
        if completed < steps and self._stopped is None:
            raise RuntimeError("Native FFS_MD stopped without a documented stop condition.")
        return NativeChunk(completed, self._stopped, self.observe())

    def step(self) -> Observation:
        """Run one native integration step; production normally uses chunks."""
        return self.advance(1).observation


def probe_runtime(runtime: dict) -> dict:
    """Check actual native FFS stops and NA forces against the MD backend.

    This rejects older oxDNA FFS builds that fail to reset forces each step.
    An import/version check cannot detect that physical integration defect.
    No experimental calculation is created and no user input is modified.
    """
    import tempfile
    try:
        from .oxdna_ffs_worker import validate_request, _same_checkpoint
    except ImportError:
        from oxdna_ffs_worker import validate_request, _same_checkpoint
    with tempfile.TemporaryDirectory(prefix="aptaswitch_ffs_probe_") as temporary:
        root = Path(temporary)
        topology, configuration = root / "probe.top", root / "probe.dat"
        topology.write_text("2 2 5->3\nG type=DNA\nC type=RNA\n")
        configuration.write_text(
            "t = 1\nb = 40 40 40\nE = 0 0 0\n"
            "5 5 5 1 0 0 0 0 1 0.03 0.02 0.01 0.01 0.02 0.03\n"
            "6.5 5 5 -1 0 0 0 0 -1 -0.03 -0.02 -0.01 0.03 0.02 0.01\n")
        request = validate_request({
            "runtime": runtime,
            "prepared": {"topology_file": str(topology), "configuration_file": str(configuration)},
            "sampling": {"temperature_c": 37., "salt_m": .12},
            "plan": {"switch_length": 1, "nucleotide_count": 2,
                     "binding_pairs": [[0, 1]], "stem_pairs": [], "opening_pairs": [],
                     "interfaces": [{"name": "bound", "min_binding": 1}]},
        })
        plain = OxpyAdapter(request, configuration, 31, False)
        try:
            plain.manager.run(10, print_output=False)
            expected = plain.frame()
        finally:
            plain.close()
        native = NativeFFSAdapter(request, configuration, 31, False,
                                  directory=root / "equivalence", conditions=[])
        try:
            chunk = native.advance(10)
            if chunk.steps != 10 or not _same_checkpoint(expected, native.frame()):
                raise ValueError("Native FFS_MD does not reproduce NA MD forces. Update the oxDNA FFS CPU backend.")
        finally:
            native.close()
        # A deterministic small approach crosses between two observations;
        # the C++ backend must interrupt a large requested chunk exactly there.
        plain = OxpyAdapter(request, configuration, 31, False)
        try:
            distances = [plain.observe().distance]
            for _ in range(10):
                distances.append(plain.step().distance)
        finally:
            plain.close()
        decreasing = min(distances[1:]) < distances[0]
        bound = (distances[0] + (min(distances[1:]) if decreasing else max(distances[1:]))) / 2
        operator = "<=" if decreasing else ">="
        expected_step = next(i for i, value in enumerate(distances)
                             if (value <= bound if decreasing else value >= bound))
        native = NativeFFSAdapter(request, configuration, 31, False,
                                  directory=root / "stop", conditions=[[f"distance {operator} {bound:.17g}"]])
        try:
            chunk = native.advance(100)
            if chunk.condition != 0 or chunk.steps != expected_step or expected_step == 0:
                raise ValueError("The native FFS_MD backend failed its first-crossing stop check.")
            checkpoint = native.frame()
        finally:
            native.close()
        resumed = root / "checkpoint.dat"
        resumed.write_text(checkpoint)
        native = NativeFFSAdapter(request, resumed, 41, False,
                                  directory=root / "resume", conditions=[])
        try:
            if not _same_checkpoint(checkpoint, native.frame()):
                raise ValueError("The native FFS restart altered the stored phase point.")
            chunk = native.advance(2)
            if chunk.steps != 2 or chunk.observation.step != 1 + expected_step + 2:
                raise ValueError("The native FFS restart did not preserve the physical step counter.")
        finally:
            native.close()
    return {"backend": "CPU", "sim_type": "FFS_MD", "interaction_type": "NA",
            "force_equivalence": True, "native_first_crossing": True,
            "phase_point_restart": True, "probe_version": 1}


if __name__ == "__main__":
    import json
    import sys
    if len(sys.argv) != 3 or sys.argv[1] != "--probe":
        raise SystemExit("Usage: oxdna_native_ffs.py --probe runtime-json")
    print(json.dumps(probe_runtime(json.loads(sys.argv[2])), allow_nan=False))
