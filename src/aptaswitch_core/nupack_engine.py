"""Generic NUPACK engine adapter.

NUPACK is intentionally not bundled. The caller must provide an import path that
belongs to a user-provided installation.
"""

from __future__ import annotations

from .localization import tr

import importlib
import sys
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

from .architecture import architecture_layout
from .models import DesignCandidate, RunRequest, RunResult, STOP_FILTER_DISABLED_WARNING
from .nupack_setup import configured_nupack_import_path
from .run_storage import new_run_directory
from .scoring import compute_selection_score, sort_candidates
from .sequences import (
    as_nupack_dna,
    as_nupack_rna,
    normalize_rna,
    reading_frame_report,
    reverse_complement_rna,
)


@contextmanager
def nupack_import_context(import_path: str):
    path = str(Path(import_path).expanduser().resolve())
    sys.path.insert(0, path)
    try:
        yield importlib.import_module("nupack")
    finally:
        try:
            sys.path.remove(path)
        except ValueError:
            pass


def _domain_seq(sequence: str, chem: str) -> str:
    """Material-prefixed sequence for a NUPACK domain (allows ``N`` for design)."""
    seq = "".join(sequence.split()).upper()
    seq = seq.replace("T", "U") if chem == "rna" else seq.replace("U", "T")
    prefix = "r" if chem == "rna" else "d"
    return prefix + seq


def _to_chem(sequence: str, chem: str) -> str:
    """Strip material prefixes and normalize letters to the given chemistry."""
    seq = sequence.replace("d", "").replace("r", "").replace("m", "").upper()
    return seq.replace("T", "U") if chem == "rna" else seq.replace("U", "T")


def _parse_cross_pairs(structure: str) -> tuple[list[tuple[int, int]], int]:
    """Return base pairs and the index of the strand break for a 2-strand MFE."""
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    break_index = 0
    index = 0
    for char in structure:
        if char == "+":
            break_index = index
            continue
        if char in "([{<":
            stack.append(index)
        elif char in ")]}>":
            if stack:
                pairs.append((stack.pop(), index))
        index += 1
    return pairs, break_index


def analyze_aptamer_trigger(
    import_path: str,
    aptamer_dna: str,
    trigger_dna: str,
    thermo,
    trigger_material: str = "dna",
) -> dict:
    """Compute the aptamer+trigger MFE and the trigger bases bound to the aptamer.

    Returns a dict with ``paired_trigger_indices`` (0-based indices into the
    trigger that pair with the aptamer), plus the structure and free energy.
    Uses a mixed RNA/DNA model so a DNA aptamer can bind a DNA or RNA trigger.
    """
    with nupack_import_context(import_path) as nupack:
        model = nupack.Model(
            material="rna-dna06",
            celsius=thermo.temperature_c,
            sodium=max(thermo.sodium_m, 0.12),
            magnesium=0.0,
        )
        result = nupack.mfe(
            strands=[
                _domain_seq(aptamer_dna, "dna"),
                _domain_seq(trigger_dna, trigger_material),
            ],
            model=model,
        )
        structure = str(result[0].structure)
        energy = float(result[0].energy)

    aptamer_len = len(aptamer_dna)
    pairs, _ = _parse_cross_pairs(structure)
    paired: set[int] = set()
    for left, right in pairs:
        # A trigger base is bound to the aptamer when exactly one end is in the
        # aptamer segment ([0, aptamer_len)) and the other is in the trigger.
        left_apt = left < aptamer_len
        right_apt = right < aptamer_len
        if left_apt != right_apt:
            trigger_index = (right if left_apt else left) - aptamer_len
            if 0 <= trigger_index < len(trigger_dna):
                paired.add(trigger_index)
    return {
        "structure": structure,
        "energy": energy,
        "paired_trigger_indices": sorted(paired),
        "aptamer_len": aptamer_len,
        "trigger_len": len(trigger_dna),
    }


def _layout(request: RunRequest) -> dict:
    return architecture_layout(request.architecture)


def _create_model(nupack, request: RunRequest):
    thermo = request.thermo
    return nupack.Model(
        material=thermo.material,
        celsius=thermo.temperature_c,
        sodium=thermo.sodium_m,
        magnesium=thermo.magnesium_m,
    )


def _build_design_spec(nupack, request: RunRequest, model):
    arch = request.architecture
    layout = _layout(request)
    trigger_binding = reverse_complement_rna(request.sequences.effective_trigger())
    toehold_seq = trigger_binding[: arch.toehold_length]
    stem_l_seq = trigger_binding[arch.toehold_length :]
    reporter_seq = normalize_rna(request.reporter.sequence_after_start_rna, field_name="reporter")
    switch_chem = request.thermo.switch_material
    trigger_chem = request.thermo.trigger_material

    Domain = nupack.Domain
    TargetComplex = nupack.TargetComplex
    TargetStrand = nupack.TargetStrand
    TargetTube = nupack.TargetTube
    SetSpec = nupack.SetSpec
    DesignOptions = nupack.DesignOptions

    t7_leader = (
        Domain(_domain_seq(arch.t7_leader, switch_chem), name="t7_leader", material=switch_chem)
        if arch.t7_leader
        else None
    )
    toehold = Domain(_domain_seq(toehold_seq, switch_chem), name="toehold", material=switch_chem)
    stem_l = Domain(_domain_seq(stem_l_seq, switch_chem), name="stem_l", material=switch_chem)
    aug = Domain(_domain_seq("AUG", switch_chem), name="aug", material=switch_chem)
    rbs = Domain(_domain_seq(arch.rbs_loop, switch_chem), name="rbs", material=switch_chem)
    rbs_prefix = (
        Domain(_domain_seq("N" * arch.rbs_prefix_length, switch_chem), name="rbs_prefix", material=switch_chem)
        if arch.rbs_prefix_length
        else None
    )
    frame_linker = (
        Domain(_domain_seq(layout["frame_linker"], switch_chem), name="frame_linker", material=switch_chem)
        if layout["frame_linker"]
        else None
    )
    reporter = Domain(_domain_seq(reporter_seq, switch_chem), name="reporter_after_start", material=switch_chem)
    trigger_domain = Domain(_domain_seq(request.sequences.effective_trigger(), trigger_chem), name="trigger", material=trigger_chem)
    aptamer_domain = Domain(_domain_seq(request.sequences.aptamer_dna, "dna"), name="aptamer", material="dna")

    bulge_l = Domain(_domain_seq("N" * arch.bulge_length, switch_chem), name="bulge_l", material=switch_chem)
    upper_l = Domain(_domain_seq("N" * arch.upper_stem_length, switch_chem), name="upper_l", material=switch_chem)
    upper_r = Domain(_domain_seq("N" * arch.upper_stem_length, switch_chem), name="upper_r", material=switch_chem)
    upper2_l = (
        Domain(_domain_seq("N" * arch.upper_stem2_length, switch_chem), name="upper2_l", material=switch_chem)
        if arch.upper_stem2_length
        else None
    )
    upper2_r = (
        Domain(_domain_seq("N" * arch.upper_stem2_length, switch_chem), name="upper2_r", material=switch_chem)
        if arch.upper_stem2_length
        else None
    )
    stem_r = Domain(_domain_seq("N" * arch.stem_length, switch_chem), name="stem_r", material=switch_chem)
    aug_spacer = (
        Domain(_domain_seq("N" * layout["len_aug_spacer"], switch_chem), name="aug_spacer", material=switch_chem)
        if layout["len_aug_spacer"]
        else None
    )

    switch_domains = [
        toehold,
        stem_l,
        bulge_l,
        upper_l,
        rbs,
        upper_r,
        aug,
        stem_r,
    ]
    if upper2_l is not None:
        switch_domains.insert(switch_domains.index(bulge_l), upper2_l)
        switch_domains.insert(switch_domains.index(stem_r), upper2_r)
    if rbs_prefix is not None:
        switch_domains.insert(switch_domains.index(rbs), rbs_prefix)
    if aug_spacer is not None:
        switch_domains.insert(switch_domains.index(aug), aug_spacer)
    if t7_leader is not None:
        switch_domains.insert(0, t7_leader)
    if frame_linker is not None:
        switch_domains.append(frame_linker)
    switch_domains.append(reporter)

    trigger = TargetStrand([trigger_domain], name="Trigger")
    aptamer = TargetStrand([aptamer_domain], name="Aptamer")
    switch = TargetStrand(switch_domains, name="Switch")

    structure_off = (
        "." * layout["len_t7"]
        + "." * layout["len_toehold"]
        + "(" * layout["len_stem"]
        + "(" * layout["len_upper2"]
        + "." * layout["len_bulge"]
        + "(" * layout["len_upper"]
        + "." * layout["len_rbs_prefix"]
        + "." * layout["len_rbs"]
        + ")" * layout["len_upper"]
        + "." * layout["len_aug_spacer"]
        + "." * len("AUG")
        + ")" * layout["len_upper2"]
        + ")" * layout["len_stem"]
        + "." * len(layout["frame_linker"])
        + "." * len(reporter_seq)
    )
    structure_switch_on = (
        "." * layout["len_t7"]
        + ")" * layout["len_toehold"]
        + ")" * layout["len_stem"]
        + "." * layout["len_upper2"]
        + "." * layout["len_bulge"]
        + "." * layout["len_upper"]
        + "." * layout["len_rbs_prefix"]
        + "." * layout["len_rbs"]
        + "." * layout["len_upper"]
        + "." * layout["len_aug_spacer"]
        + "." * len("AUG")
        + "." * layout["len_upper2"]
        + "." * layout["len_stem"]
        + "." * len(layout["frame_linker"])
        + "." * len(reporter_seq)
    )
    structure_on = "(" * len(request.sequences.effective_trigger()) + "+" + structure_switch_on

    switch_off = TargetComplex([switch], structure_off, name="switch_off")
    switch_on = TargetComplex([trigger, switch], structure_on, name="switch_on")
    trigger_free = TargetComplex(
        [trigger],
        "." * len(request.sequences.effective_trigger()),
        name="trigger_free",
    )
    trigger_concentration = request.thermo.resolved_trigger_concentration_m
    switch_concentration = request.thermo.resolved_switch_concentration_m
    duplex_concentration = min(trigger_concentration, switch_concentration)
    on_targets = {switch_on: duplex_concentration}
    if trigger_concentration > duplex_concentration:
        on_targets[trigger_free] = trigger_concentration - duplex_concentration
    if switch_concentration > duplex_concentration:
        on_targets[switch_off] = switch_concentration - duplex_concentration
    tube_off = TargetTube(
        on_targets={switch_off: switch_concentration},
        off_targets=SetSpec(max_size=request.thermo.screening_max_complex_size),
        name="tube_off",
    )
    tube_on = TargetTube(
        on_targets=on_targets,
        off_targets=SetSpec(max_size=request.thermo.screening_max_complex_size),
        name="tube_on",
    )
    spec = nupack.tube_design(
        tubes=[tube_off, tube_on],
        model=model,
        options=DesignOptions(f_stop=0.005, wobble_mutations=True),
    )
    return {
        "spec": spec,
        "trigger": trigger,
        "aptamer": aptamer,
        "switch": switch,
        "layout": layout,
    }


def _complex_pct(nupack, *, strands, model, denominator_m, max_size, expected_names):
    Strand = nupack.Strand
    Tube = nupack.Tube
    SetSpec = nupack.SetSpec
    tube_strands = {
        Strand(sequence, name=name): float(concentration_m)
        for name, sequence, concentration_m in strands
    }
    tube = Tube(
        strands=tube_strands,
        complexes=SetSpec(max_size=max_size),
        name="aptaswitch_analysis",
    )
    result = nupack.tube_analysis(tubes=[tube], model=model)
    expected_conc = 0.0
    for complex_obj, conc in result.tubes[tube].complex_concentrations.items():
        names = [strand.name for strand in complex_obj.strands]
        if sorted(names) == sorted(expected_names) and len(names) == len(expected_names):
            expected_conc += float(conc)
    return round((expected_conc / denominator_m * 100.0) if denominator_m > 0 else 0.0, 4)


def _first_non_linear_on_mfe(
    records, *, trigger_length: int, switch_length: int, start: int, end: int,
):
    """Return a failing ON MFE record, or None when the region is all unpaired.

    The interval is the switch's RBS-loop prefix through the last linker base,
    excluding the reporter. Every MFE proxy structure returned by NUPACK must
    pass; this makes no claim about other structures in the ensemble. Pairing
    with any base, including the trigger or the reporter, fails the screen.
    """
    if not (0 <= start < end <= switch_length):
        raise ValueError(tr('Criblage RBS–linker : bornes de région invalides.'))
    if not records:
        raise ValueError(tr('Criblage RBS–linker : aucune structure MFE ON renvoyée.'))
    first_failure = None
    for record in records:
        structure = str(record.structure)
        strands = structure.split("+")
        if len(strands) != 2 or list(map(len, strands)) != [trigger_length, switch_length]:
            raise ValueError(tr('Criblage RBS–linker : les longueurs de la structure MFE ON ne correspondent pas aux brins.'))
        balance = 0
        for symbol in structure:
            if symbol == "(":
                balance += 1
            elif symbol == ")":
                balance -= 1
            elif symbol not in ".+":
                raise ValueError(tr('Criblage RBS–linker : structure MFE ON mal formée.'))
            if balance < 0:
                raise ValueError(tr('Criblage RBS–linker : structure MFE ON mal formée.'))
        if balance:
            raise ValueError(tr('Criblage RBS–linker : structure MFE ON mal formée.'))
        if first_failure is None and any(symbol != "." for symbol in strands[1][start:end]):
            first_failure = record
    return first_failure


def run_nupack_engine(request: RunRequest, progress_callback=None, cancel_check=None) -> RunResult:
    request = request.normalized()
    import_path = request.nupack_path or configured_nupack_import_path()
    if not import_path:
        raise RuntimeError(tr('NUPACK n’est pas configuré. Sélectionnez d’abord votre chemin NUPACK.'))

    started = datetime.now().astimezone()
    candidates: list[DesignCandidate] = []
    if progress_callback:
        progress_callback(0, request.trials, tr('NUPACK · préparation du modèle et des tubes…'))
    with nupack_import_context(import_path) as nupack:
        model = _create_model(nupack, request)
        design = _build_design_spec(nupack, request, model)
        spec = design["spec"]
        trigger = design["trigger"]
        switch = design["switch"]
        layout = design["layout"]

        for trial in range(1, request.trials + 1):
            if cancel_check and cancel_check():
                break
            if progress_callback:
                progress_callback(
                    trial - 1,
                    request.trials,
                    tr('NUPACK · essai {v0}/{v1} : optimisation…', v0=trial, v1=request.trials),
                )

            result = spec.run(trials=1)[0]
            switch_chem = request.thermo.switch_material
            trigger_chem = request.thermo.trigger_material
            trigger_seq = _to_chem(str(result.to_analysis(trigger)), trigger_chem)
            switch_seq = _to_chem(str(result.to_analysis(switch)), switch_chem)
            frame = reading_frame_report(switch_seq, layout["aug_start"], layout["reporter_start_codon"])
            if request.exclude_stop_candidates and frame["has_premature_stop"]:
                candidates.append(DesignCandidate.excluded_for_stop(
                    trial=trial, switch_seq=switch_seq, trigger_seq=trigger_seq,
                    first_stop=frame["first_stop"], reporter_label=request.reporter.label,
                ))
                if progress_callback:
                    progress_callback(trial, request.trials,
                                      tr('NUPACK · essai {v0}/{v1} exclu : STOP prématuré au codon {v2}. MFE et analyses en tube ignorées.', v0=trial, v1=request.trials, v2=frame['first_stop']))
                continue
            trigger_nupack = _domain_seq(trigger_seq, trigger_chem)
            switch_nupack = _domain_seq(switch_seq, switch_chem)
            if progress_callback:
                phase = tr('criblage RBS–linker dans les structures MFE ON') if request.exclude_non_linear_candidates else tr('analyses thermodynamiques')
                progress_callback(trial - 1, request.trials,
                                  tr('NUPACK · essai {v0}/{v1} : {v2}…', v0=trial, v1=request.trials, v2=phase))
            # Screen before OFF, free-trigger, fragment and tube calculations.
            # The retained ON result is reused below, without another MFE call.
            mfe_on = list(nupack.mfe(strands=[trigger_nupack, switch_nupack], model=model))
            if not mfe_on:
                raise ValueError(tr('NUPACK n’a renvoyé aucune structure MFE ON.'))
            if request.exclude_non_linear_candidates:
                failing_on = _first_non_linear_on_mfe(
                    mfe_on, trigger_length=len(trigger_seq), switch_length=len(switch_seq),
                    start=layout["rbs_linker_start"], end=layout["rbs_linker_end"],
                )
                if failing_on is not None:
                    candidates.append(DesignCandidate.excluded_for_non_linear(
                        trial=trial, switch_seq=switch_seq, trigger_seq=trigger_seq,
                        structure_on=str(failing_on.structure), delta_g_on=float(failing_on.energy),
                        has_stop_codon=bool(frame["has_premature_stop"]), first_stop=frame["first_stop"],
                        reporter_label=request.reporter.label,
                    ))
                    if progress_callback:
                        progress_callback(trial, request.trials,
                                          tr('NUPACK · essai {v0}/{v1} exclu : RBS–linker non linéaire en ON. Autres MFE et analyses en tube ignorées.', v0=trial, v1=request.trials))
                    continue
                if progress_callback:
                    progress_callback(trial - 1, request.trials,
                                      tr('NUPACK · essai {v0}/{v1} : criblage RBS–linker validé, analyses thermodynamiques…', v0=trial, v1=request.trials))
            defect = float(result.defects.ensemble_defect)

            mfe_off = nupack.mfe(strands=[switch_nupack], model=model)
            delta_g_off = float(mfe_off[0].energy)
            structure_off = str(mfe_off[0].structure)
            delta_g_on = float(mfe_on[0].energy)
            structure_on = str(mfe_on[0].structure)
            mfe_trigger = nupack.mfe(strands=[trigger_nupack], model=model)
            delta_g_trigger = float(mfe_trigger[0].energy)
            ddg = delta_g_on - delta_g_off - delta_g_trigger

            rbs_linker_seq = switch_seq[layout["rbs_linker_start"]:layout["rbs_linker_end"]]
            mfe_rbs = nupack.mfe(strands=[_domain_seq(rbs_linker_seq, switch_chem)], model=model)
            delta_g_rbs = float(mfe_rbs[0].energy)

            on_yield = _complex_pct(
                nupack,
                strands=[
                    ("Trigger", trigger_nupack, request.thermo.resolved_trigger_concentration_m),
                    ("Switch", switch_nupack, request.thermo.resolved_switch_concentration_m),
                ],
                model=model,
                denominator_m=request.thermo.resolved_switch_concentration_m,
                max_size=request.thermo.validation_max_complex_size,
                expected_names=["Trigger", "Switch"],
            )
            leak = _complex_pct(
                nupack,
                strands=[
                    (
                        "Aptamer",
                        _domain_seq(request.sequences.aptamer_dna, "dna"),
                        request.thermo.resolved_trigger_concentration_m,
                    ),
                    ("Trigger", trigger_nupack, request.thermo.resolved_trigger_concentration_m),
                    ("Switch", switch_nupack, request.thermo.resolved_switch_concentration_m),
                ],
                model=model,
                denominator_m=request.thermo.resolved_switch_concentration_m,
                max_size=request.thermo.validation_max_complex_size,
                expected_names=["Trigger", "Switch"],
            )
            score = compute_selection_score(
                on_yield_pct=on_yield,
                leak_pct=leak,
                defect=defect,
                ddg_activation=ddg,
                delta_g_rbs_linker=delta_g_rbs,
                weights=request.scoring,
            )
            candidates.append(
                DesignCandidate(
                    trial=trial,
                    switch_seq=switch_seq,
                    trigger_seq=trigger_seq,
                    score=score,
                    defect=round(defect, 6),
                    on_yield_pct=on_yield,
                    leak_pct=leak,
                    delta_g_off=round(delta_g_off, 6),
                    delta_g_on=round(delta_g_on, 6),
                    ddg_activation=round(ddg, 6),
                    delta_g_rbs_linker=round(delta_g_rbs, 6),
                    has_bad_rbs_linker=delta_g_rbs < -8.0,
                    has_stop_codon=bool(frame["has_premature_stop"]),
                    first_stop=frame["first_stop"],
                    structure_off=structure_off,
                    structure_on=structure_on,
                    reporter_label=request.reporter.label,
                )
            )
            if progress_callback:
                progress_callback(
                    trial,
                    request.trials,
                    tr('NUPACK · essai {v0}/{v1} terminé.', v0=trial, v1=request.trials),
                )

    completed = datetime.now().astimezone()
    run_id, output_dir = new_run_directory(request.output_dir, "design", started)
    return RunResult(
        run_id=run_id,
        request=request,
        candidates=sort_candidates(candidates),
        output_dir=output_dir,
        started_at=started.isoformat(),
        completed_at=completed.isoformat(),
        status="cancelled" if len(candidates) < request.trials else "completed",
        warnings=[] if request.exclude_stop_candidates else [tr(STOP_FILTER_DISABLED_WARNING)],
    )
