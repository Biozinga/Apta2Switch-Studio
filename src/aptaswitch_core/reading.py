"""Prepare an exact, user-supplied switch for the existing metric pipeline."""

from __future__ import annotations

from dataclasses import replace

from .architecture import CLASSIC_30NT_LINKER, architecture_layout, default_architecture_for_trigger
from .localization import tr
from .models import REPORTERS, ReporterContext, RunRequest, SwitchArchitecture
from .sequences import normalize_dna, normalize_rna, reverse_complement_rna


def _reporter_starts(switch: str, minimum: int, aug_start: int, full: str) -> list[int]:
    """Locate an exact reporter prefix with enough bases to identify its start.

    Twenty-one bases retain compatibility with the application's earliest
    saved reporter fragments. A shorter complete custom reporter must match
    in its entirety. A short coincidence with a longer CDS is not used as an
    annotation.
    """
    anchor_length = min(21, len(full))
    return [start for start in range(minimum, len(switch) - anchor_length + 1)
            if (start - aug_start - 3) % 3 == 0
            and full.startswith(switch[start:])]


def infer_reading_architecture(
    switch_sequence: str, trigger_sequence: str, reporter: ReporterContext, *,
    reporter_is_complete: bool = False,
) -> SwitchArchitecture:
    """Recognize a supported activator family using only supplied sequences.

    Recognition requires the whole trigger's reverse complement, the family's
    RBS and initiation AUG. The leader is the actual preceding sequence. The
    reporter identifies the exact linker boundary, so optimized linker bases
    and lengths are retained. No parameter from an old design form is read.
    """
    switch = normalize_rna(switch_sequence, field_name="Switch en lecture")
    trigger = normalize_dna(trigger_sequence, field_name="trigger")
    if len(trigger) < 2:
        raise ValueError(tr("La région du trigger utilisée pour le switch doit contenir au moins 2 bases."))
    recognized = reverse_complement_rna(trigger)
    family = default_architecture_for_trigger(len(trigger))
    full_reporters = []
    reporter.validate()
    if reporter.full_sequence_after_start_rna is not None:
        full_reporters.append(normalize_rna(reporter.full_sequence_after_start_rna))
    if reporter_is_complete:
        # A complete supplied reporter may be recognized independently of the
        # disabled reporter selector. This identifies bases, not provenance.
        full_reporters.extend(normalize_rna(value.full_sequence_after_start_rna)
                              for value in REPORTERS.values()
                              if value.full_sequence_after_start_rna is not None)
    candidates = {}
    anchor = switch.find(recognized)
    found_motifs = False
    while anchor >= 0:
        architecture = replace(
            family, t7_leader=switch[:anchor],
            architecture_name=f"reading_inferred_{len(trigger)}nt_{family.binding_layout}",
        )
        layout = architecture_layout(architecture)
        aug_start = layout["aug_start"]
        if (switch[layout["rbs_start"]:layout["rbs_start"] + len(architecture.rbs_loop)] == architecture.rbs_loop
                and switch[aug_start:aug_start + 3] == "AUG"):
            found_motifs = True
            linker_start = aug_start + 3 + architecture.upper_stem2_length + architecture.stem_length
            starts = {start for full in full_reporters
                      for start in _reporter_starts(switch, linker_start, aug_start, full)}
            if (reporter_is_complete and family.binding_layout == "full_left_arm"
                    and switch[linker_start:linker_start + len(CLASSIC_30NT_LINKER)] == CLASSIC_30NT_LINKER):
                start = linker_start + len(CLASSIC_30NT_LINKER)
                if start < len(switch):
                    starts.add(start)
            for reporter_start in starts:
                inferred = replace(architecture, frame_linker=switch[linker_start:reporter_start])
                inferred.validate_for_trigger(len(trigger))
                # Offsets fully define the annotation, independent of names.
                key = (anchor, aug_start, reporter_start)
                candidates[key] = inferred
        anchor = switch.find(recognized, anchor + 1)
    if len(candidates) > 1:
        raise ValueError(tr(
            "Plusieurs découpages du switch correspondent aux séquences fournies. "
            "Le leader, le linker ou le début du rapporteur est ambigu ; fournissez le rapporteur correspondant pour lever l’ambiguïté."
        ))
    if candidates:
        return next(iter(candidates.values()))
    if found_motifs:
        if reporter_is_complete:
            raise ValueError(tr(
                "Le trigger, le RBS et AUG sont reconnus, mais la frontière entre le linker et le rapporteur ne peut pas être déterminée. "
                "Sélectionnez ou ajoutez le rapporteur correspondant dans Réglages ; un rapporteur inconnu nécessite un linker classique identifiable."
            ))
        raise ValueError(tr(
            "Le début du rapporteur n’est pas reconnu dans le switch fourni. "
            "Sélectionnez le rapporteur correspondant et fournissez au moins 21 bases de son début, ou sa séquence complète si elle est plus courte."
        ))
    raise ValueError(tr(
        "L’architecture du switch n’est pas reconnue automatiquement. "
        "Vérifiez le trigger entier, sa séquence complémentaire dans le switch, le RBS et AUG ; seules les familles de toehold activateurs prises en charge peuvent être annotées."
    ))


def prepare_reading_request(
    request: RunRequest, switch_sequence: str, *, reporter_is_complete: bool = False,
    infer_architecture: bool = True,
) -> RunRequest:
    """Freeze a supplied switch without redesigning or replacing its bases.

    Architecture is recognized from the switch, whole trigger and reporter.
    Explicit ``infer_architecture=False`` is reserved for API callers that
    already have exact annotations for an unsupported or imperfect activator.
    A known reporter contributes its full CDS only on an exact prefix match.
    """
    request.thermo.validate_model()
    request.thermo.validate_switch_chemistry()
    switch = normalize_rna(switch_sequence, field_name="Switch en lecture")
    if infer_architecture:
        request = replace(request, architecture=infer_reading_architecture(
            switch, request.sequences.trigger_dna, request.reporter,
            reporter_is_complete=reporter_is_complete,
        ))
    layout = architecture_layout(request.architecture)
    reporter_start = layout["reporter_start"]
    if len(switch) <= reporter_start:
        raise ValueError(tr(
            "Le switch fourni doit contenir une séquence de rapporteur après le linker. "
            "Vérifiez la séquence complète et l’architecture sélectionnée."
        ))
    actual_suffix = switch[reporter_start:]
    if reporter_is_complete:
        reporter = ReporterContext(
            key="reading_reporter", label="Rapporteur fourni",
            sequence_after_start_rna=actual_suffix,
            full_sequence_after_start_rna=actual_suffix,
        )
    else:
        request.reporter.validate()
        full = request.reporter.full_sequence_after_start_rna
        if full is None or not normalize_rna(full, field_name="Rapporteur complet").startswith(actual_suffix):
            raise ValueError(tr(
                "La séquence du rapporteur dans le switch fourni ne correspond pas au rapporteur complet sélectionné. "
                "Sélectionnez le rapporteur correspondant, ou indiquez que le suffixe fourni est le rapporteur complet."
            ))
        reporter = replace(request.reporter, sequence_after_start_rna=actual_suffix)

    # Validate against the original motifs first, including every fixed base
    # of an explicit linker. Only X bases and automatic padding can vary.
    prepared = replace(
        request, reporter=reporter, fixed_switch_sequence=switch,
        sequences=replace(request.sequences, binding_start=0, binding_length=None),
        trials=1, inverse_extension=None, optimization_options={},
        exclude_stop_candidates=False, exclude_non_linear_candidates=False,
    ).normalized()
    linker_length = len(layout["frame_linker"])
    actual_linker = switch[reporter_start - linker_length:reporter_start]
    return replace(
        prepared,
        architecture=replace(prepared.architecture, frame_linker=actual_linker),
    )
