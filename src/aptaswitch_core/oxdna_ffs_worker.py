"""Native oxNA primitives and entry point for resumable direct FFS.

The controller launches independent native FFS_MD replica processes.
Rate estimates are expressed in oxDNA model time, not experimental seconds.
"""

from __future__ import annotations

import gc
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from dataclasses import dataclass
from typing import Any, Callable

try:
    from .oxdna_box import check_periodic_self_clearance
except ImportError:  # Direct invocation in the external oxpy interpreter.
    from oxdna_box import check_periodic_self_clearance


@dataclass(frozen=True)
class Observation:
    binding: int
    stem: int
    opening: int
    distance: float
    step: int
    associations: int | None = None
    separation: float | None = None


def _pairs(values: Any, name: str) -> tuple[tuple[int, int], ...]:
    pairs = set()
    for pair in values:
        if len(pair) != 2 or any(type(i) is not int or i < 0 for i in pair):
            raise ValueError(f"Invalid nucleotide indices in {name}.")
        i, j = sorted(pair)
        if i == j:
            raise ValueError(f"A nucleotide cannot pair with itself in {name}.")
        pairs.add((i, j))
    return tuple(sorted(pairs))


def distance_pair_groups(plan: dict) -> dict[str, tuple[tuple[int, int], ...]]:
    """One definition for native and observed centre-to-centre PBC distances.

    Complete plans follow the trigger 3′–toehold entry pair; artificial
    proximity plans follow the entire trigger against the exposed toehold.
    Frozen requests retain their original coordinate. Separation monitors both whole
    strands; only ``binding_pairs`` count the intended hydrogen bonds.
    """
    binding = _pairs(plan.get("binding_pairs", []), "binding_pairs")
    scope = plan.get("distance_scope")
    if scope is None and not any(key in plan for key in ("distance_pairs", "separation_pairs")):
        return {"distance": binding}
    if scope not in {"trigger_to_toehold", "trigger_3prime_to_toehold_start"}:
        raise ValueError("Explicit distance pairs require a supported toehold distance scope.")
    n, total = plan.get("switch_length"), plan.get("nucleotide_count")
    if type(n) is not int or type(total) is not int or not 0 < n < total:
        raise ValueError("Distance groups require two nonempty strands.")
    toehold = _pairs(plan.get("toehold_pairs", []), "toehold_pairs")
    if (not toehold or not set(toehold).issubset(binding)
            or any(not 0 <= i < n <= j < total for i, j in toehold)):
        raise ValueError("Toehold pairs must be nonempty intended switch/trigger pairs.")
    exposed = sorted({i for i, _ in toehold})
    if scope == "trigger_3prime_to_toehold_start":
        approach = ((exposed[0], total - 1),)
        if approach[0] not in toehold:
            raise ValueError("The toehold start must target the trigger's 3′ terminal nucleotide.")
    else:
        approach = tuple((i, j) for i in exposed for j in range(n, total))
    expected = {
        "distance": approach,
        "separation": tuple((i, j) for i in range(n) for j in range(n, total)),
    }
    for name, pairs in expected.items():
        if _pairs(plan.get(f"{name}_pairs", []), f"{name}_pairs") != pairs:
            raise ValueError(f"Invalid {name}_pairs for the declared toehold distance scope.")
    return expected


def validate_request(request: dict[str, Any]) -> dict[str, Any]:
    """Validate native geometry and MD settings without importing oxpy.

    Input fields are retained unchanged unless their native numeric type needs
    normalization. The identity therefore remains tied to the exact frozen
    scientific request while the controller owns its production allocation.
    """
    request = json.loads(json.dumps(request))
    plan = request["plan"]
    if plan.setdefault("distance_metric", "nucleotide_centres") != "nucleotide_centres":
        raise ValueError("The native distance criterion uses nucleotide centres.")
    n, total = plan.get("switch_length"), plan.get("nucleotide_count")
    if type(n) is not int or type(total) is not int or not 0 < n < total:
        raise ValueError("The native plan must contain two nonempty strands.")
    for key in ("binding_pairs", "stem_pairs", "opening_pairs"):
        plan[key] = _pairs(plan.get(key, []), key)
        if any(max(pair) >= total for pair in plan[key]):
            raise ValueError(f"Pair index exceeds the topology in {key}.")
    if not plan["binding_pairs"] or any(not i < n <= j for i,j in plan["binding_pairs"]):
        raise ValueError("Binding pairs must connect the switch and trigger.")
    distances = distance_pair_groups(plan)
    if "separation" in distances:
        plan["toehold_pairs"] = _pairs(plan["toehold_pairs"], "toehold_pairs")
        for name, pairs in distances.items():
            plan[f"{name}_pairs"] = pairs
        if "switch_sequence" in request or "trigger_sequence" in request:
            switch = request.get("switch_sequence", "").upper().replace("U", "T")
            trigger = request.get("trigger_sequence", "").upper().replace("U", "T")
            bases = switch + trigger
            complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
            if (len(switch) != n or len(bases) != total
                    or any(complement.get(bases[i]) != bases[j] for i, j in plan["binding_pairs"])):
                raise ValueError("The intended binding pairs do not match the RNA/DNA sequence register.")
    minimum = plan.setdefault("basin_min_stem", 0)
    if type(minimum) is not int or not 0 <= minimum <= len(plan["stem_pairs"]):
        raise ValueError("The preparation stem requirement must be attainable.")
    for key, default in (("basin_distance", 4.), ("hb_threshold", -.1)):
        plan[key] = float(plan.get(key, default))
        if not math.isfinite(plan[key]):
            raise ValueError(f"Nonfinite {key}.")
    if plan["basin_distance"] <= 0 or plan["hb_threshold"] >= 0:
        raise ValueError("Preparation separation must be positive and the hydrogen-bond threshold negative.")
    if not plan.get("interfaces"):
        raise ValueError("A terminal opening criterion is required.")
    names = set()
    for interface in plan["interfaces"]:
        name = interface["name"]
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Interface names must be nonempty and unique.")
        names.add(name)
        binding = interface.setdefault("min_binding", 0)
        if type(binding) is not int or not 0 <= binding <= len(plan["binding_pairs"]):
            raise ValueError("An interface binding count must be attainable.")
        if interface.get("distance_max") is not None:
            distance = float(interface["distance_max"])
            if not math.isfinite(distance) or distance <= 0:
                raise ValueError("An interface distance must be finite and positive.")
            interface["distance_max"] = distance
        for field, pairs in (("max_stem", "stem_pairs"), ("max_opening", "opening_pairs")):
            value = interface.get(field)
            if value is not None and (type(value) is not int or not 0 <= value <= len(plan[pairs])):
                raise ValueError(f"An interface {field} limit must be attainable.")
    settings = request.setdefault("sampling", {})
    if settings.get("diversity") is not None:
        raise ValueError("Adaptive exploration banks are not a kinetic production protocol.")
    if any(settings.get(key) is not None for key in ("budget_seconds", "max_steps_per_trial")):
        raise ValueError("Kinetic trajectories cannot be censored by time or step limits.")
    for key, default in (("frame_interval", 1000), ("max_movie_frames", 600), ("newtonian_steps", 103)):
        value = settings.setdefault(key, default)
        if type(value) is not int or value < 1:
            raise ValueError(f"{key} must be a positive integer.")
    if settings["max_movie_frames"] < 2:
        raise ValueError("At least two movie frames are required.")
    for key, default in (("dt", .003), ("diff_coeff", 2.5)):
        value = float(settings.setdefault(key, default))
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{key} must be finite and positive.")
        settings[key] = value
    for key in ("temperature_c", "salt_m"):
        value = float(settings[key])
        if not math.isfinite(value) or (key == "salt_m" and value <= 0) or (key == "temperature_c" and value <= -273.15):
            raise ValueError(f"Invalid {key}.")
        settings[key] = value
    if type(settings.get("seed", 1)) is not int:
        raise ValueError("The random seed must be an integer.")
    settings.setdefault("seed", 1)
    if request.get("use_average_seq", False) is not False:
        raise ValueError("This protocol requires all three sequence-dependent NA parameter files.")
    return request


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def request_identity(request: dict[str, Any]) -> str:
    """Bind checkpoints to the scientific request and actual input contents."""
    identity = json.loads(json.dumps(validate_request(request)))
    identity.pop("resume", None)
    runtime = identity.pop("runtime", {})
    files = {**identity.pop("prepared", {}), **{
        key: value for key, value in runtime.items() if key.startswith("seq_dep_file_")}}
    identity["input_files"] = {
        key: hashlib.sha256(Path(value).read_bytes()).hexdigest() if Path(value).is_file() else str(Path(value).resolve())
        for key, value in files.items() if key.endswith("_file") or key.startswith("seq_dep_file_")}
    identity["parameter_sha256"] = runtime.get("parameter_sha256")
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


def _tuples(value):
    return tuple(_tuples(item) for item in value) if isinstance(value, list) else value


class FrameBuffer:
    """Thin long segments without losing their first or last actual frame."""

    def __init__(self, interval: int, capacity: int):
        self.interval = interval
        self.capacity = max(2, capacity)
        self.frames: list[tuple[int, str]] = []
        self.next_step = 0

    def add(self, step: int, frame: Callable[[], str], *, force: bool = False) -> None:
        if self.frames and self.frames[-1][0] == step:
            return
        if not force and self.frames and step < self.next_step:
            return
        self.frames.append((step, frame()))
        self.next_step = step + self.interval
        if len(self.frames) > self.capacity:
            self.frames = [self.frames[0], *self.frames[2:-1:2], self.frames[-1]]
            self.interval *= 2


class OxpyAdapter:
    """One native NA MD trajectory; no external forces or velocity refresh in shots."""

    def __init__(self, request: dict[str, Any], configuration: Path, seed: int, initial: bool):
        import oxpy  # Supplied by the external runtime, never imported by the UI process.

        self.context = oxpy.Context()
        self.context.__enter__()
        self.manager = None
        self.info = None
        self.observable = None
        self.hb_observable = None
        try:
            settings = request["sampling"]
            input_file = oxpy.InputFile()
            values = {
                "backend": "CPU", "sim_type": "MD", "interaction_type": "NA",
                "topology": str(Path(request["prepared"]["topology_file"]).resolve()),
                "conf_file": str(configuration.resolve()), "seed": str(seed),
                "T": f"{settings['temperature_c']}C", "salt_concentration": str(settings["salt_m"]),
                "dt": str(settings["dt"]), "thermostat": "brownian",
                "newtonian_steps": str(settings["newtonian_steps"]),
                "diff_coeff": str(settings["diff_coeff"]), "verlet_skin": "0.2",
                "refresh_vel": "true" if initial else "false", "restart_step_counter": "false",
                "steps": "0", "time_scale": "linear", "print_conf_interval": "1000000000",
                "print_energy_every": "1000000000", "external_forces": "false",
                "use_average_seq": "false", "trajectory_file": os.devnull,
                "lastconf_file": os.devnull, "energy_file": os.devnull, "log_file": os.devnull,
                "data_output_1": "{\nname = " + os.devnull + "\nprint_every = 1000000000\nonly_last = true\ncol_1 = {\ntype = Configuration\nid = ffs_checkpoint\n}\ncol_2 = {\ntype = hb_list\nid = ffs_hb\n}\n}",
            }
            for field in ("seq_dep_file_DNA", "seq_dep_file_RNA", "seq_dep_file_NA"):
                parameter_file = Path(request["runtime"][field]).resolve()
                if not parameter_file.is_file():
                    raise ValueError(f"Missing NA parameter file: {field}.")
                values[field] = str(parameter_file)
            self._configure_input(values)
            for key, value in values.items():
                input_file[key] = value
            self.manager = oxpy.OxpyManager(input_file)
            self.info = self.manager.config_info()
            self.observable = self.info.get_observable_by_id("ffs_checkpoint")
            self.hb_observable = self.info.get_observable_by_id("ffs_hb")
            particles = self.info.particles()
            self.switch_length = request["plan"].get("switch_length")
            if not isinstance(self.switch_length, int) or not 0 < self.switch_length < len(particles):
                raise ValueError("Native sampling requires the switch_length strand boundary.")
            if request["plan"].get("nucleotide_count", len(particles)) != len(particles):
                raise ValueError("The planned nucleotide count does not match the topology.")
            self.particles = particles
            self.pairs = {}
            self.pair_indices = {}
            for key in ("binding_pairs", "stem_pairs", "opening_pairs"):
                indices = request["plan"][key]
                if any(max(pair) >= len(particles) for pair in indices):
                    raise ValueError(f"Pair index exceeds the topology in {key}.")
                if key == "binding_pairs" and any(not i < self.switch_length <= j for i, j in indices):
                    raise ValueError("Intended binding pairs must connect the switch and trigger.")
                self.pairs[key] = [(particles[i], particles[j]) for i, j in indices]
                self.pair_indices[key] = set(indices)
            self.distance_pairs = {
                name: [(particles[i], particles[j]) for i, j in indices]
                for name, indices in distance_pair_groups(request["plan"]).items()
            }
            self.hb_threshold = request["plan"]["hb_threshold"]
            self.box_policy = request.get("box_policy")
            self.box_sides = list(self.info.box.box_sides)
            if self.box_policy and any(not math.isclose(side, request['box_size'], abs_tol=1e-8)
                                       for side in self.box_sides):
                raise ValueError("The native box does not match the fixed automatic FFS volume.")
            self.periodic_clearance = None
        except BaseException:
            self.close()
            raise

    def _configure_input(self, values: dict[str, str]) -> None:
        """Allow a native backend to add its documented input settings."""

    def observe(self) -> Observation:
        if self.box_policy:
            # The contour estimate is not a hard extension bound. Check the
            # actual unwrapped chains each step, before accepting a crossing,
            # with a margin outside the native nonbonded cutoff.
            self.periodic_clearance = check_periodic_self_clearance(
                [list(particle.pos) for particle in self.particles],
                [self.switch_length, len(self.particles)-self.switch_length],
                self.box_sides, self.box_policy["interaction_cutoff"])
        self.info.interaction.begin_energy_computation()
        if self.hb_threshold == -0.1:
            # Native HBList evaluates the NA hydrogen-bond term on the neighbor
            # list. Its official HB_CUTOFF is -0.1 oxDNA energy units. Intersect
            # indices in Python rather than evaluating every possible RBS pair.
            bonded = {tuple(sorted(pair)) for pair in self.hb_observable.hb_list()}
            counts = {key: len(indices & bonded) for key, indices in self.pair_indices.items()}
            associations = sum(i < self.switch_length <= j for i, j in bonded)
        else:
            counts = {
                key: sum(self.info.interaction.pair_interaction_term(4, p, q) < self.hb_threshold for p, q in pairs)
                for key, pairs in self.pairs.items()
            }
            associations = sum(
                self.info.interaction.pair_interaction_term(4, p, q) < self.hb_threshold
                for p in self.particles[:self.switch_length] for q in self.particles[self.switch_length:]
            )
        distances = {
            name: math.sqrt(min(self.info.box.sqr_min_image_distance(p, q) for p, q in pairs))
            for name, pairs in self.distance_pairs.items()
        }
        if any(not math.isfinite(distance) for distance in distances.values()):
            raise ValueError("Nonfinite native distance.")
        return Observation(counts["binding_pairs"], counts["stem_pairs"], counts["opening_pairs"],
                           distances["distance"], self.manager.current_step, associations,
                           distances.get("separation"))

    def step(self) -> Observation:
        before = self.manager.current_step
        self.manager.run(1, print_output=False)
        if self.manager.current_step != before + 1:
            raise RuntimeError("The native MD engine stopped before completing a step.")
        return self.observe()

    def frame(self) -> str:
        self.manager.update_observable_data()
        return self.observable.get_output_string(self.manager.current_step).rstrip() + "\n"

    def close(self) -> None:
        # oxpy permits only one live manager. Drop all dependent references first.
        self.pairs = {}
        self.distance_pairs = {}
        self.particles = []
        self.observable = None
        self.hb_observable = None
        self.info = None
        self.manager = None
        gc.collect()
        context, self.context = self.context, None
        if context is not None:
            context.__exit__(None, None, None)


class SamplingStopped(Exception):
    def __init__(self, reason: str):
        self.reason = reason
        super().__init__(reason)


def _same_checkpoint(first: str, second: str) -> bool:
    """Compare physical state across a native restart, allowing text roundoff.

    oxpy recomputes the energy header and orthonormalizes axes on import, so
    byte equality is too strict. No coordinates or momenta are interpolated.
    """
    left, right = first.splitlines(), second.splitlines()
    if len(left) != len(right) or len(left) < 4:
        return False
    try:
        steps = [float(lines[0].split("=", 1)[1]) for lines in (left, right)]
        if not all(map(math.isfinite, steps)) or steps[0] != steps[1]:
            return False
        left_values = [left[1].split("=", 1)[1], *left[3:]]
        right_values = [right[1].split("=", 1)[1], *right[3:]]
        for left_row, right_row in zip(left_values, right_values):
            a, b = list(map(float, left_row.split())), list(map(float, right_row.split()))
            if len(a) != len(b) or not a:
                return False
            if any(not math.isfinite(x) or not math.isfinite(y)
                   or not math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-10)
                   for x, y in zip(a, b)):
                return False
    except (ValueError, IndexError):
        return False
    return True


def main() -> int:
    import fcntl  # Supported runtime platforms: macOS and Linux.

    started = time.monotonic()
    lock = None
    acquired = False
    sampler_started = False
    try:
        if len(sys.argv) != 2:
            raise ValueError("Usage: python oxdna_ffs_worker.py request.json")
        request = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
        lock = Path("worker.lock").open("a+")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("An oxDNA worker is already running in this simulation directory.") from None
        acquired = True
        atomic_json(Path("worker_process.json"), {"pid": os.getpid(), "started_at": time.time()})
        if __package__:
            from .oxdna_kinetics_worker import KineticSampler
        else:
            from oxdna_kinetics_worker import KineticSampler
        sampler = KineticSampler(request, Path.cwd(), started=started)
        sampler_started = True
        result = sampler.run()
    except Exception as error:
        result = {"schema_version": 1, "status": "failed", "stop_reason": "invalid_request",
                  "error": str(error), "kinetic_rate": None, "trajectory_file": None}
        # A rejected resume or duplicate worker must not destroy the actual run.
        if acquired and (sampler_started or not any(Path(name).exists() for name in ("resume_state.json", "banks.json", "result.json", "progress.json"))):
            atomic_json(Path("result.json"), result)
            atomic_json(Path("progress.json"), result)
        elif acquired:
            atomic_json(Path("resume_error.json"), result)
        print(str(error), file=sys.stderr)
    finally:
        if lock is not None:
            lock.close()
    return 1 if result["status"] == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
