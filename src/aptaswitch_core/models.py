"""Runtime data models for Apta2Switch-Studio."""

from __future__ import annotations

from .localization import tr

from dataclasses import asdict, dataclass, field
import math
from numbers import Real
from pathlib import Path
from typing import Any, Callable

from .run_storage import default_export_dir
from .sequences import OPTIMIZABLE_BASE, normalize_dna, normalize_rna
from .reporter_sequences import (
    SFGFP_CDS_DNA, SFGFP_ACCESSION, SFGFP_SOURCE_URL,
    MNEONGREEN_CDS_DNA, MNEONGREEN_ACCESSION, MNEONGREEN_SOURCE_URL,
    MCHERRY_CDS_DNA, MCHERRY_ACCESSION, MCHERRY_SOURCE_URL,
    MRFP1_CDS_DNA, MRFP1_ACCESSION, MRFP1_SOURCE_URL,
    LACZ_CDS_DNA, LACZ_ACCESSION, LACZ_SOURCE_URL,
)


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
    # Historical designs recognize only the toehold and lower stem. The
    # corrected 30-nt activator recognizes the entire contiguous left arm.
    binding_layout: str = "legacy_lower_stem"

    @property
    def trigger_region_length(self) -> int:
        length = self.toehold_length + self.stem_length
        if self.binding_layout == "full_left_arm":
            length += self.upper_stem2_length + self.bulge_length + self.upper_stem_length
        return length

    @property
    def rbs_to_aug_distance(self) -> int:
        """Bases strictly between the end of the RBS and the start of AUG."""
        return self.upper_stem_length + self.aug_spacer

    def resolved_frame_linker(self) -> str:
        """Return a design pattern with the padding needed for the reading frame.

        Automatic bases are free ``X`` domains for the sequence optimizer to evaluate against
        the target structures. Only an explicit custom pattern fixes bases.
        """
        if self.frame_linker is not None:
            if not isinstance(self.frame_linker, str) or set(self.frame_linker.upper()) - set("ACGTUX"):
                raise ValueError(tr("Le linker ne peut contenir que A, C, G, T, U ou X (base à optimiser)."))
            return self.frame_linker.upper().replace("T", "U")
        post_aug = self.upper_stem2_length + self.stem_length
        return OPTIMIZABLE_BASE * ((3 - (post_aug % 3)) % 3)

    def validate_sequence_patterns(self) -> None:
        """Fixed leader/RBS and an optional linker with explicit free sites."""
        for sequence in (self.t7_leader, self.rbs_loop):
            if not isinstance(sequence, str) or set(sequence.upper()) - set("ACGTU"):
                raise ValueError(tr("Le leader et le RBS doivent contenir uniquement des bases déterminées A, C, G, T ou U."))
        self.resolved_frame_linker()

    def validate_for_trigger(self, trigger_length: int) -> None:
        self.validate_sequence_patterns()
        if self.binding_layout not in {"legacy_lower_stem", "full_left_arm"}:
            raise ValueError(tr('La disposition de la région reconnue du trigger est invalide.'))
        if self.toehold_length <= 0 or self.stem_length <= 0:
            raise ValueError(tr('Les longueurs du toehold et du stem doivent être positives.'))
        if self.trigger_region_length != trigger_length:
            raise ValueError(
                tr('La région reconnue doit correspondre à la longueur du trigger ({v0} != {v1}).', v0=self.trigger_region_length, v1=trigger_length)
            )
        if self.upper_stem_length <= 0:
            raise ValueError(tr('La longueur du stem supérieur doit être positive.'))
        if self.upper_stem2_length < 0:
            raise ValueError(tr('La longueur du stem supérieur facultatif ne peut pas être négative.'))
        if self.bulge_length <= 0:
            raise ValueError(tr('La longueur du bulge doit être positive.'))
        if self.rbs_prefix_length < 0 or self.aug_spacer < 0:
            raise ValueError(tr('Les longueurs du préfixe RBS et de l’espace avant AUG ne peuvent pas être négatives.'))
        if (self.upper_stem2_length + self.stem_length + len(self.resolved_frame_linker())) % 3:
            raise ValueError(tr('Le linker personnalisé place le rapporteur hors du cadre de lecture : la longueur après AUG avant le rapporteur doit être un multiple de 3.'))


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
    material: str = "experimental-Sugimoto95-Vienna04-loopmean"
    switch_material: str = "rna"
    trigger_material: str = "dna"
    screening_max_complex_size: int = 2
    validation_max_complex_size: int = 3

    def validate_model(self) -> None:
        """Validate a new calculation without rewriting archived model identities."""
        if self.material != "experimental-Sugimoto95-Vienna04-loopmean":
            raise ValueError("Unsupported thermodynamic model; use experimental-Sugimoto95-Vienna04-loopmean for a new calculation.")
        if self.switch_material not in {"dna", "rna"} or self.trigger_material not in {"dna", "rna"}:
            raise ValueError(tr("La chimie de chaque brin doit être ADN ou ARN."))
        values = (self.temperature_c, self.sodium_m, self.magnesium_m)
        if any(isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) for value in values):
            raise ValueError(tr("Les conditions thermodynamiques doivent être finies."))
        if self.magnesium_m != 0:
            raise ValueError(tr("Le modèle ADN/ARN impose Mg²⁺ = 0."))
        low, high = 25, 55
        if not low <= self.temperature_c <= high or not .12 <= self.sodium_m <= 1:
            raise ValueError(f"{self.material} requires {low}–{high} °C, Na⁺ 0.12–1 M and Mg²⁺ = 0.")

    def validate_switch_chemistry(self) -> None:
        """The designed switch is RNA and its trigger is DNA."""
        if (self.switch_material, self.trigger_material) != ("rna", "dna"):
            raise ValueError(tr("La conception nécessite un toehold switch ARN et un trigger ADN."))

    def validate_system_chemistry(self) -> None:
        """The initial aptamer and trigger are both DNA."""
        if (self.switch_material, self.trigger_material) != ("dna", "dna"):
            raise ValueError(tr("Le système initial nécessite un aptamère ADN et un trigger ADN."))

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
class InverseExtensionSettings:
    """Joint switch/extension design, followed by an independent structure check."""

    target_length: int = 30
    extension_side: str = "5prime"
    top_candidates: int = 20
    thermo: ThermoConditions = field(default_factory=lambda: ThermoConditions(switch_material="dna", trigger_material="dna", temperature_c=37))

    def extension_lengths(self, core_length: int) -> tuple[int, int]:
        if type(self.target_length) is not int or self.target_length <= core_length:
            raise ValueError(tr("La longueur cible de l’extension inverse doit dépasser celle du trigger original."))
        if type(self.top_candidates) is not int or self.top_candidates <= 0:
            raise ValueError(tr("Le nombre de switches à comparer doit être un entier positif."))
        added = self.target_length - core_length
        if self.extension_side == "5prime":
            return added, 0
        if self.extension_side == "3prime":
            return 0, added
        if self.extension_side == "both":
            return (added + 1) // 2, added // 2
        raise ValueError(tr("La position de l’extension inverse est invalide."))

    def validate(self, core_length: int) -> None:
        self.extension_lengths(core_length)
        self.thermo.validate_model()
        self.thermo.validate_system_chemistry()


@dataclass(frozen=True)
class ReporterContext:
    key: str = "sfgfp"
    label: str = "sfGFP"
    sequence_after_start_rna: str = "CGUAAAGGCGAGGAGCUGUUC"
    # New reporters store their full CDS after the initial AUG independently of
    # the short region used during sequence design. None preserves historical prefixes.
    full_sequence_after_start_rna: str | None = None
    source_accession: str | None = None
    source_url: str | None = None

    def validate(self) -> None:
        if not isinstance(self.key, str) or not self.key or not isinstance(self.label, str):
            raise ValueError(tr('Le rapporteur doit avoir un nom et un identifiant valides.'))
        prefix = self.sequence_after_start_rna
        if not isinstance(prefix, str) or set(prefix.upper()) - set("ACGTU"):
            raise ValueError(tr('La séquence du rapporteur est invalide.'))
        full = self.full_sequence_after_start_rna
        if full is not None and (not isinstance(full, str) or not full
                                 or set(full.upper()) - set("ACGTU")
                                 or not full.upper().replace("T", "U").startswith(prefix.upper().replace("T", "U"))):
            raise ValueError(tr('Le rapporteur complet doit commencer exactement par la séquence utilisée pour le design.'))
        for source in (self.source_accession, self.source_url):
            if source is not None and not isinstance(source, str):
                raise ValueError(tr('La provenance du rapporteur est invalide.'))


REPORTERS = {
    "sfgfp": ReporterContext(
        sequence_after_start_rna=SFGFP_CDS_DNA[3:33].replace("T", "U"),
        full_sequence_after_start_rna=SFGFP_CDS_DNA[3:].replace("T", "U"),
        source_accession=SFGFP_ACCESSION,
        source_url=SFGFP_SOURCE_URL,
    ),
    "mng": ReporterContext(
        key="mng",
        label="mNeonGreen",
        sequence_after_start_rna=MNEONGREEN_CDS_DNA[3:33].replace("T", "U"),
        full_sequence_after_start_rna=MNEONGREEN_CDS_DNA[3:].replace("T", "U"),
        source_accession=MNEONGREEN_ACCESSION,
        source_url=MNEONGREEN_SOURCE_URL,
    ),
    "mcherry": ReporterContext(
        key="mcherry", label="mCherry",
        sequence_after_start_rna=MCHERRY_CDS_DNA[3:33].replace("T", "U"),
        full_sequence_after_start_rna=MCHERRY_CDS_DNA[3:].replace("T", "U"),
        source_accession=MCHERRY_ACCESSION, source_url=MCHERRY_SOURCE_URL,
    ),
    "mrfp1": ReporterContext(
        key="mrfp1", label="mRFP1",
        sequence_after_start_rna=MRFP1_CDS_DNA[3:33].replace("T", "U"),
        full_sequence_after_start_rna=MRFP1_CDS_DNA[3:].replace("T", "U"),
        source_accession=MRFP1_ACCESSION, source_url=MRFP1_SOURCE_URL,
    ),
    "lacz": ReporterContext(
        key="lacz", label="LacZ (E. coli K-12 MG1655)",
        sequence_after_start_rna=LACZ_CDS_DNA[3:33].replace("T", "U"),
        full_sequence_after_start_rna=LACZ_CDS_DNA[3:].replace("T", "U"),
        source_accession=LACZ_ACCESSION, source_url=LACZ_SOURCE_URL,
    ),
}


@dataclass(frozen=True)
class ScoringWeights:
    on_yield: float = 3.0
    leak: float = 0.0
    defect: float = 6.0
    ddg: float = 0.25
    rbs: float = 6.0
    trigger_linearity: float = 3.0
    rbs_calculator: float = 3.0




@dataclass(frozen=True)
class FunctionalDefectWeights:
    """Relative defects of OFF accessibility/lock and two ON open regions."""

    off_toehold: float = 1.0
    off_stem: float = 1.0
    on_rbs_aug: float = 1.0
    on_rbs_linker: float = 0.75

    def validate(self) -> None:
        values = (self.off_toehold, self.off_stem, self.on_rbs_aug, self.on_rbs_linker)
        for value in values:
            try:
                valid = (not isinstance(value, bool) and isinstance(value, Real)
                         and math.isfinite(value) and value >= 0)
            except (OverflowError, TypeError, ValueError):
                valid = False
            if not valid:
                raise ValueError(tr('Les poids des quatre zones doivent être des nombres finis positifs ou nuls.'))
        if not any(value > 0 for value in values):
            raise ValueError(tr('Au moins un poids de zone doit être strictement positif. Pour désactiver ce critère, mettez le poids global du défaut à zéro.'))

    def fractions(self) -> dict[str, float]:
        """Normalize finite weights without overflowing their sum."""
        self.validate()
        values = {name: float(value) for name, value in asdict(self).items()}
        scale = max(values.values())
        scaled = {name: value / scale for name, value in values.items()}
        total = math.fsum(scaled.values())
        return {name: value / total for name, value in scaled.items()}


@dataclass(frozen=True)
class RunRequest:
    sequences: SequenceInput
    architecture: SwitchArchitecture
    thermo: ThermoConditions = field(default_factory=ThermoConditions)
    reporter: ReporterContext = field(default_factory=lambda: REPORTERS["sfgfp"])
    scoring: ScoringWeights = field(default_factory=ScoringWeights)
    trials: int = 50
    output_dir: Path = field(default_factory=default_export_dir)
    engine: str = "native-hybrid"
    exclude_stop_candidates: bool = True
    exclude_non_linear_candidates: bool = True
    inverse_extension: InverseExtensionSettings | None = None
    defect_metric: str = "functional"
    functional_defect_weights: FunctionalDefectWeights = field(default_factory=FunctionalDefectWeights)
    score_scale: str = "bounded_v1"
    thermodynamics_library_path: str | None = None
    optimization_options: dict[str, Any] = field(default_factory=dict)
    fixed_switch_sequence: str | None = None

    def normalized(self) -> "RunRequest":
        if self.engine != "native-hybrid":
            raise ValueError("This calculation requires the native-hybrid engine.")
        sequences = self.sequences.normalized()
        self.thermo.validate_model()
        self.thermo.validate_switch_chemistry()
        self.reporter.validate()
        if self.inverse_extension is not None:
            self.inverse_extension.validate(len(sequences.trigger_dna))
            if sequences.binding_start != 0 or sequences.effective_trigger() != sequences.trigger_dna:
                raise ValueError(tr("L’extension inverse conserve le trigger original entier ; une fenêtre partielle n’est pas compatible."))
            self.architecture.validate_for_trigger(self.inverse_extension.target_length)
        else:
            self.architecture.validate_for_trigger(len(sequences.effective_trigger()))
        if type(self.trials) is not int or self.trials <= 0:
            raise ValueError(tr('Le nombre d’essais doit être positif.'))
        if not math.isfinite(self.thermo.resolved_trigger_concentration_m) or self.thermo.resolved_trigger_concentration_m <= 0:
            raise ValueError(tr('La concentration du trigger doit être positive.'))
        if not math.isfinite(self.thermo.resolved_switch_concentration_m) or self.thermo.resolved_switch_concentration_m <= 0:
            raise ValueError(tr('La concentration du toehold switch doit être positive.'))
        for size in (self.thermo.screening_max_complex_size, self.thermo.validation_max_complex_size):
            if type(size) is not int or not 2 <= size <= 3:
                raise ValueError("Switch design tubes must include complexes of two to three strands.")
        if self.defect_metric != "functional":
            raise ValueError("The defect scoring metric must be functional.")
        if self.score_scale != "bounded_v1":
            raise ValueError("The selection score scale must be bounded_v1.")
        value = self.scoring.trigger_linearity
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value < 0:
            raise ValueError(tr("Le poids de linéarité du trigger doit être un nombre fini positif ou nul."))
        value = self.scoring.rbs
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value < 0:
            raise ValueError(tr("Le poids RBS–linker doit être un nombre fini positif ou nul."))
        value = self.scoring.rbs_calculator
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value < 0:
            raise ValueError("The RBS Calculator weight must be finite and nonnegative.")
        if not isinstance(self.functional_defect_weights, FunctionalDefectWeights):
            raise ValueError(tr('Les poids des quatre zones sont invalides.'))
        self.functional_defect_weights.validate()
        from .thermodynamics.optimization import TrialOptions
        if not isinstance(self.optimization_options, dict):
            raise ValueError("Optimization options must be a JSON object.")
        try:
            optimization_options = asdict(TrialOptions(**self.optimization_options))
        except TypeError as exc:
            raise ValueError("Unknown optimization options.") from exc
        if self.thermodynamics_library_path is not None and not isinstance(self.thermodynamics_library_path, str):
            raise ValueError("The thermodynamics library path must be text.")
        fixed_switch = None
        if self.fixed_switch_sequence is not None:
            if (self.trials != 1 or self.inverse_extension is not None
                    or sequences.binding_start != 0 or sequences.effective_trigger() != sequences.trigger_dna):
                raise ValueError(tr("Le mode lecture analyse un seul switch fixe avec le trigger entier, sans extension inverse."))
            if self.engine != "native-hybrid":
                raise ValueError("Fixed-sequence reading requires the native-hybrid engine.")
            fixed_switch = normalize_rna(self.fixed_switch_sequence, field_name="Switch en lecture")
            from .architecture import architecture_layout
            layout = architecture_layout(self.architecture)
            reporter = normalize_rna(self.reporter.sequence_after_start_rna, field_name="reporter")
            if len(fixed_switch) != layout["reporter_start"] + len(reporter):
                raise ValueError(tr("La longueur du switch fourni ne correspond pas à l’architecture et au rapporteur."))
            motifs = (
                (0, self.architecture.t7_leader, "leader"),
                (layout["rbs_start"], self.architecture.rbs_loop, "RBS"),
                (layout["aug_start"], "AUG", "AUG"),
                (layout["reporter_start"] - len(layout["frame_linker"]), layout["frame_linker"], "linker"),
                (layout["reporter_start"], reporter, "rapporteur"),
            )
            for start, motif, label in motifs:
                motif = motif.upper().replace("T", "U")
                actual = fixed_switch[start:start + len(motif)]
                if any(expected != OPTIMIZABLE_BASE and expected != base for expected, base in zip(motif, actual)):
                    raise ValueError(tr("Le motif {motif} du switch fourni ne correspond pas à l’architecture enregistrée.", motif=label))
        return RunRequest(
            sequences=sequences,
            architecture=self.architecture,
            thermo=self.thermo,
            reporter=self.reporter,
            scoring=self.scoring,
            trials=self.trials,
            output_dir=Path(self.output_dir),
            engine=self.engine,
            exclude_stop_candidates=self.exclude_stop_candidates if fixed_switch is None else False,
            exclude_non_linear_candidates=self.exclude_non_linear_candidates if fixed_switch is None else False,
            inverse_extension=self.inverse_extension,
            defect_metric=self.defect_metric,
            functional_defect_weights=self.functional_defect_weights,
            score_scale=self.score_scale,
            thermodynamics_library_path=self.thermodynamics_library_path,
            optimization_options=optimization_options,
            fixed_switch_sequence=fixed_switch,
        )

    def design_trigger_pattern(self) -> str:
        core = self.sequences.effective_trigger()
        if self.inverse_extension is None:
            return core
        prefix, suffix = self.inverse_extension.extension_lengths(len(core))
        return OPTIMIZABLE_BASE * prefix + core + OPTIMIZABLE_BASE * suffix

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["output_dir"] = str(self.output_dir)
        if data["fixed_switch_sequence"] is None:
            data.pop("fixed_switch_sequence")
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
    inverse_extension_metrics: dict[str, Any] = field(default_factory=dict)
    defect_structural: float | None = None
    defect_concentration: float | None = None
    defect_breakdown: dict[str, Any] = field(default_factory=dict)
    regional_defects: dict[str, Any] = field(default_factory=dict)
    trigger_linearity: dict[str, Any] = field(default_factory=dict)
    rbs_accessibility: dict[str, Any] = field(default_factory=dict)
    rbs_calculator: dict[str, Any] = field(default_factory=dict)
    score_raw: float | None = None

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
    inverse_extension_analysis: dict[str, Any] = field(default_factory=dict)
    thermodynamics: dict[str, Any] = field(default_factory=dict)
    optimization: dict[str, Any] = field(default_factory=dict)

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
            "inverse_extension_analysis": self.inverse_extension_analysis,
            "thermodynamics": self.thermodynamics,
            "optimization": self.optimization,
        }
