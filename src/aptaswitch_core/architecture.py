"""Switch architecture helpers."""

from __future__ import annotations

from .localization import tr
from .models import SwitchArchitecture


DEFAULT_TRIGGER_LENGTH = 24
DEFAULT_TOEHOLD_LENGTH = 15
DEFAULT_RBS_AUG_DISTANCE = 6
DEFAULT_RBS_PREFIX_LENGTH = 3


def default_architecture_for_trigger(trigger_length: int) -> SwitchArchitecture:
    """Choose a conservative toehold/stem split for the trigger length."""
    if trigger_length == DEFAULT_TRIGGER_LENGTH:
        toehold = DEFAULT_TOEHOLD_LENGTH
    elif trigger_length <= 18:
        toehold = max(8, trigger_length // 2)
    elif trigger_length <= 24:
        toehold = 13
    else:
        toehold = 12
    toehold = min(max(1, toehold), trigger_length - 1)
    stem = trigger_length - toehold
    return SwitchArchitecture(
        toehold_length=toehold,
        stem_length=stem,
        upper_stem_length=DEFAULT_RBS_AUG_DISTANCE,
        # The 30-nt activator ends its coding stem exactly 18 nt after AUG.
        # No extra upper stem or reading-frame padding is needed in this format.
        upper_stem2_length=0 if trigger_length == 30 else 4,
        rbs_prefix_length=DEFAULT_RBS_PREFIX_LENGTH,
        architecture_name=f"runtime_{trigger_length}nt_toehold{toehold}_stem{stem}",
    )


def architecture_layout(architecture: SwitchArchitecture) -> dict[str, int | str]:
    """Shared sequence offsets for design, frame scanning and visualizations.

    The upper stem's right arm lies between the fixed RBS and AUG. It is
    designed by NUPACK, giving the standard activator six optimized bases
    in that interval. The separate prefix forms the unpaired RBS loop.
    ``aug_spacer`` remains supported for explicitly customized legacy runs.
    """
    arch = architecture
    # Results already open during a source reload can still contain the old
    # architecture class: those sequences did not contain an RBS prefix.
    rbs_prefix_length = getattr(arch, "rbs_prefix_length", 0)
    aug_spacer = getattr(arch, "aug_spacer", 0)
    rbs_to_aug_distance = arch.upper_stem_length + aug_spacer
    frame_linker = arch.resolved_frame_linker()
    rbs_prefix_start = (
        len(arch.t7_leader) + arch.trigger_region_length
        + arch.upper_stem2_length + arch.bulge_length + arch.upper_stem_length
    )
    rbs_start = rbs_prefix_start + rbs_prefix_length
    aug_start = rbs_start + len(arch.rbs_loop) + rbs_to_aug_distance
    reporter_start = (
        aug_start + 3 + arch.upper_stem2_length + arch.stem_length + len(frame_linker)
    )
    return {
        "len_t7": len(arch.t7_leader),
        "len_toehold": arch.toehold_length,
        "len_stem": arch.stem_length,
        "len_upper": arch.upper_stem_length,
        "len_upper2": arch.upper_stem2_length,
        "len_bulge": arch.bulge_length,
        "len_rbs": len(arch.rbs_loop),
        "len_rbs_prefix": rbs_prefix_length,
        "len_aug_spacer": aug_spacer,
        "frame_linker": frame_linker,
        "reporter_start_codon": (
            3 + arch.upper_stem2_length + arch.stem_length + len(frame_linker)
        ) // 3 + 1,
        "rbs_prefix_start": rbs_prefix_start,
        "rbs_start": rbs_start,
        "aug_start": aug_start,
        "reporter_start": reporter_start,
        # Green's RBS–linker region starts at the first base of the RBS loop,
        # including its optimized prefix, and ends immediately before reporter.
        # Nucleotide offsets must stay exact even for custom, unaligned linkers.
        "rbs_linker_start": rbs_prefix_start,
        "rbs_linker_end": reporter_start,
        "rbs_linker_length": reporter_start - rbs_prefix_start,
    }


def architecture_summary(architecture: SwitchArchitecture) -> str:
    return tr(
        "Région reconnue du trigger : {trigger_length} nt · toehold : {toehold_length} nt · tige : {stem_length} nt",
        trigger_length=architecture.trigger_region_length,
        toehold_length=architecture.toehold_length,
        stem_length=architecture.stem_length,
    )
