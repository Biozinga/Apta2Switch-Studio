"""Runtime data models for Apta2Switch-Studio."""

from __future__ import annotations

from .localization import tr

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable

from .run_storage import default_export_dir
from .sequences import normalize_dna


ProgressCallback = Callable[[int, int, str], None]
DEFAULT_RBS_SEQUENCE = "AGAGGAGA"
STOP_FILTER_DISABLED_WARNING = (
    "L’exclusion des STOP prématurés est désactivée : ces candidats seront analysés "
    "et pourront être bien classés, même si leur traduction s’arrête avant le reporter."
)


@dataclass(frozen=True)
class SequenceInput:
    aptamer_dna: str
    trigger_dna: str
    molecule_name: str = "custom target"
    binding_start: int = 0
    binding_length: int | None = None

    def _clamped_window(self, trigger: str) -> tuple[int, int]:
        start = max(0, min(self.binding_start, max(0, len(trigger) - 1)))
        if self.binding_length is None:
            length = len(trigger) - start
        else:
            length = self.binding_length
        length = max(1, min(length, len(trigger) - start))
        return start, length

    def normalized(self) -> "SequenceInput":
        trigger = normalize_dna(self.trigger_dna, field_name="trigger")
        start, length = self._clamped_window(trigger)
        return SequenceInput(
            aptamer_dna=normalize_dna(self.aptamer_dna, field_name="aptamer"),
            trigger_dna=trigger,
            molecule_name=(self.molecule_name or "custom target").strip(),
            binding_start=start,
            binding_length=length,
        )

    def effective_trigger(self) -> str:
        """The trigger sub-sequence used to build the switch toehold/stem."""
        trigger = normalize_dna(self.trigger_dna, field_name="trigger")
        start, length = self._clamped_window(trigger)
        return trigger[start : start + length]


@dataclass(frozen=True)
class SwitchArchitecture:
    toehold_length: int
    stem_length: int
    upper_stem_length: int = 6
    upper_stem2_length: int = 4
    bulge_length: int = 3
    rbs_loop: str = DEFAULT_RBS_SEQUENCE
    t7_leader: str = "GGG"
    frame_linker: str | None = None
    aug_spacer: int = 0
    architecture_name: str = "custom_runtime_architecture"
    # Three optimizable, unpaired bases precede the fixed, user-supplied RBS.
    # Appended to preserve positional compatibility with saved integrations.
    rbs_prefix_length: int = 3

    @property
    def trigger_region_length(self) -> int:
        return self.toehold_length + self.stem_length

    @property
    def rbs_to_aug_distance(self) -> int:
        """Bases strictly between the end of the RBS and the start of AUG."""
        return self.upper_stem_length + self.aug_spacer

    def resolved_frame_linker(self) -> str:
        """Return a design pattern with the padding needed for the reading frame.

        Automatic bases are free ``N`` domains for NUPACK to optimize against
        the target structures. Only an explicit custom pattern fixes bases.
        """
        if self.frame_linker is not None:
            return self.frame_linker
        post_aug = self.upper_stem2_length + self.stem_length
        return "N" * ((3 - (post_aug % 3)) % 3)

    def validate_for_trigger(self, trigger_length: int) -> None:
        if self.toehold_length <= 0 or self.stem_length <= 0:
            raise ValueError(tr('Les longueurs du toehold et du stem doivent être positives.'))
        if self.trigger_region_length != trigger_length:
            raise ValueError(
                tr('La longueur toehold + stem doit correspondre à celle du trigger ({v0} != {v1}).', v0=self.trigger_region_length, v1=trigger_length)
            )
        if self.upper_stem_length <= 0:
            raise ValueError(tr('La longueur du stem supérieur doit être positive.'))
        if self.upper_stem2_length < 0:
            raise ValueError(tr('La longueur du stem supérieur facultatif ne peut pas être négative.'))
        if self.bulge_length <= 0:
            raise ValueError(tr('La longueur du bulge doit être positive.'))
        if self.rbs_prefix_length < 0 or self.aug_spacer < 0:
            raise ValueError(tr('Les longueurs du préfixe RBS et de l’espace avant AUG ne peuvent pas être négatives.'))


@dataclass(frozen=True)
class ThermoConditions:
    temperature_c: float = 29.0
    sodium_m: float = 0.12
    magnesium_m: float = 0.0
    # Kept as the legacy shared concentration for projects created before
    # trigger/switch concentrations could be configured independently.
    concentration_m: float = 5e-6
    trigger_concentration_m: float | None = None
    switch_concentration_m: float | None = None
    material: str = "rna-dna06"
    switch_material: str = "rna"
    trigger_material: str = "dna"
    screening_max_complex_size: int = 2
    validation_max_complex_size: int = 5

    @property
    def is_mixed_material(self) -> bool:
        return "-" in self.material

    @property
    def resolved_trigger_concentration_m(self) -> float:
        """Trigger concentration, with compatibility for older projects."""

        value = self.trigger_concentration_m
        return float(self.concentration_m if value is None else value)

    @property
    def resolved_switch_concentration_m(self) -> float:
        """Toehold-switch concentration, with compatibility for older projects."""

        value = self.switch_concentration_m
        return float(self.concentration_m if value is None else value)


@dataclass(frozen=True)
class ReporterContext:
    key: str = "sfgfp"
    label: str = "sfGFP"
    sequence_after_start_rna: str = "CGUAAAGGCGAGGAGCUGUUC"


REPORTERS = {
    "sfgfp": ReporterContext(),
    "mng": ReporterContext(
        key="mng",
        label="mNeonGreen",
        sequence_after_start_rna="GUGAGCAAGGGCGAGGAGGAU",
    ),
}


@dataclass(frozen=True)
class ScoringWeights:
    on_yield: float = 3.0
    leak: float = 0.05
    defect: float = 6.0
    ddg: float = 0.25
    rbs: float = 6.0


@dataclass(frozen=True)
class RunRequest:
    sequences: SequenceInput
    architecture: SwitchArchitecture
    thermo: ThermoConditions = field(default_factory=ThermoConditions)
    reporter: ReporterContext = field(default_factory=ReporterContext)
    scoring: ScoringWeights = field(default_factory=ScoringWeights)
    trials: int = 50
    output_dir: Path = field(default_factory=default_export_dir)
    engine: str = "nupack"
    nupack_path: str | None = None
    exclude_stop_candidates: bool = True
    exclude_non_linear_candidates: bool = True

    def normalized(self) -> "RunRequest":
        sequences = self.sequences.normalized()
        self.architecture.validate_for_trigger(len(sequences.effective_trigger()))
        if self.trials <= 0:
            raise ValueError(tr('Le nombre d’essais doit être positif.'))
        if self.thermo.resolved_trigger_concentration_m <= 0:
            raise ValueError(tr('La concentration du trigger doit être positive.'))
        if self.thermo.resolved_switch_concentration_m <= 0:
            raise ValueError(tr('La concentration du toehold switch doit être positive.'))
        return RunRequest(
            sequences=sequences,
            architecture=self.architecture,
            thermo=self.thermo,
            reporter=self.reporter,
            scoring=self.scoring,
            trials=self.trials,
            output_dir=Path(self.output_dir),
            engine=self.engine,
            nupack_path=self.nupack_path,
            exclude_stop_candidates=self.exclude_stop_candidates,
            exclude_non_linear_candidates=self.exclude_non_linear_candidates,
        )

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["output_dir"] = str(self.output_dir)
        return data


@dataclass
class DesignCandidate:
    trial: int
    switch_seq: str
    trigger_seq: str
    score: float | None
    defect: float | None
    on_yield_pct: float | None
    leak_pct: float | None
    delta_g_off: float | None
    delta_g_on: float | None
    ddg_activation: float | None
    delta_g_rbs_linker: float | None
    has_bad_rbs_linker: bool | None
    has_stop_codon: bool
    first_stop: int | None = None
    structure_off: str = ""
    structure_on: str = ""
    reporter_label: str = ""
    excluded: bool = False
    exclusion_reason: str = ""

    @classmethod
    def excluded_for_stop(
        cls, *, trial: int, switch_seq: str, trigger_seq: str,
        first_stop: int, reporter_label: str,
    ) -> "DesignCandidate":
        """Retain the generated sequence without inventing uncomputed metrics."""
        return cls(
            trial=trial, switch_seq=switch_seq, trigger_seq=trigger_seq,
            score=None, defect=None, on_yield_pct=None, leak_pct=None,
            delta_g_off=None, delta_g_on=None, ddg_activation=None,
            delta_g_rbs_linker=None, has_bad_rbs_linker=None,
            has_stop_codon=True, first_stop=first_stop, reporter_label=reporter_label,
            excluded=True, exclusion_reason=f"Codon STOP prématuré — codon {first_stop}",
        )

    @classmethod
    def excluded_for_non_linear(
        cls, *, trial: int, switch_seq: str, trigger_seq: str,
        structure_on: str, delta_g_on: float, has_stop_codon: bool,
        first_stop: int | None, reporter_label: str,
    ) -> "DesignCandidate":
        """Keep the failing ON prediction, leaving skipped analyses unset."""
        return cls(
            trial=trial, switch_seq=switch_seq, trigger_seq=trigger_seq,
            score=None, defect=None, on_yield_pct=None, leak_pct=None,
            delta_g_off=None, delta_g_on=delta_g_on, ddg_activation=None,
            delta_g_rbs_linker=None, has_bad_rbs_linker=None,
            has_stop_codon=has_stop_codon, first_stop=first_stop,
            structure_on=structure_on, reporter_label=reporter_label,
            excluded=True,
            exclusion_reason=(
                "RBS–linker non linéaire en ON — au moins une base appariée "
                "dans une structure MFE renvoyée"
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class RunResult:
    run_id: str
    request: RunRequest
    candidates: list[DesignCandidate]
    output_dir: Path
    started_at: str
    completed_at: str
    status: str = "completed"
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "request": self.request.to_dict(),
            "results": [candidate.to_dict() for candidate in self.candidates],
            "output_dir": str(self.output_dir),
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "status": self.status,
            "warnings": list(self.warnings),
        }
