"""Switch architecture helpers."""

from __future__ import annotations

from .localization import tr
from .models import SwitchArchitecture


DEFAULT_TRIGGER_LENGTH = 24
DEFAULT_TOEHOLD_LENGTH = 15
DEFAULT_RBS_AUG_DISTANCE = 6
DEFAULT_RBS_PREFIX_LENGTH = 3
CLASSIC_30NT_LINKER = "AACCUGGCGGCAGCGCAAAAG"


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
    full_left_arm = trigger_length == 30
    # The 18 recognized bases after the 12-nt toehold comprise the 9-nt
    # lower arm, 3 unpaired bulge bases and the 6-nt upper arm.
    stem = 9 if full_left_arm else trigger_length - toehold
    return SwitchArchitecture(
        toehold_length=toehold,
        stem_length=stem,
        upper_stem_length=DEFAULT_RBS_AUG_DISTANCE,
        upper_stem2_length=0 if full_left_arm else 4,
        rbs_prefix_length=DEFAULT_RBS_PREFIX_LENGTH,
        binding_layout="full_left_arm" if full_left_arm else "legacy_lower_stem",
        # Set the new default explicitly: old saved architectures with None
        # must retain their original automatic reading-frame padding.
        frame_linker=CLASSIC_30NT_LINKER if full_left_arm else None,
        architecture_name=f"runtime_{trigger_length}nt_toehold{toehold}_stem{stem}",
    )


def architecture_layout(architecture: SwitchArchitecture) -> dict[str, int | str]:
    """Shared sequence offsets for design, frame scanning and visualizations.

    The upper stem's right arm lies between the fixed RBS and AUG. It is
    designed by the sequence optimizer, giving the standard activator six optimized bases
    in that interval. The separate prefix forms the unpaired RBS loop.
    ``aug_spacer`` remains supported for explicitly customized legacy runs.
    """
    arch = architecture
    rbs_prefix_length = arch.rbs_prefix_length
    aug_spacer = arch.aug_spacer
    rbs_to_aug_distance = arch.upper_stem_length + aug_spacer
    frame_linker = arch.resolved_frame_linker()
    rbs_prefix_start = (
        len(arch.t7_leader) + arch.toehold_length + arch.stem_length
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
        "len_trigger": sum(part.stop - part.start for part in trigger_binding_slices(arch).values()),
        "len_rbs": len(arch.rbs_loop),
        "len_rbs_prefix": rbs_prefix_length,
        "len_aug_spacer": aug_spacer,
        "frame_linker": frame_linker,
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


def trigger_binding_slices(architecture: SwitchArchitecture) -> dict[str, slice]:
    """Return fixed-domain slices into the trigger's reverse complement.

    Every recognized domain is contiguous on the switch's left arm. Keeping
    these slices shared prevents the design and its preview from assigning
    different bases to the bulge or the upper stem.
    """
    arch = architecture
    lengths = [("toehold", arch.toehold_length), ("stem_l", arch.stem_length)]
    if arch.binding_layout == "full_left_arm":
        lengths.extend([
            ("upper2_l", arch.upper_stem2_length),
            ("bulge_l", arch.bulge_length),
            ("upper_l", arch.upper_stem_length),
        ])
    result = {}
    offset = 0
    for name, length in lengths:
        if length:
            result[name] = slice(offset, offset + length)
        offset += length
    return result


def target_switch_structures(
    architecture: SwitchArchitecture, reporter_length: int,
) -> tuple[str, str]:
    """Shared OFF target and ON target (trigger + switch), in dot-bracket form."""
    layout = architecture_layout(architecture)
    off = (
        "." * (layout["len_t7"] + layout["len_toehold"])
        + "(" * (layout["len_stem"] + layout["len_upper2"])
        + "." * layout["len_bulge"]
        + "(" * layout["len_upper"]
        + "." * (layout["len_rbs_prefix"] + layout["len_rbs"])
        + ")" * layout["len_upper"]
        + "." * (layout["len_aug_spacer"] + 3)
        + ")" * (layout["len_upper2"] + layout["len_stem"])
        + "." * (len(layout["frame_linker"]) + reporter_length)
    )
    recognized = layout["len_trigger"]
    on_switch = (
        "." * layout["len_t7"] + ")" * recognized
        + "." * (len(off) - layout["len_t7"] - recognized)
    )
    return off, "(" * recognized + "+" + on_switch


def architecture_summary(architecture: SwitchArchitecture) -> str:
    if architecture.binding_layout == "full_left_arm":
        return tr(
            "Région reconnue du trigger : {trigger_length} nt · toehold : {toehold_length} nt · avant AUG : {upper_length} paires · AUG : 3 nt non appariés · après AUG : {lower_length} paires",
            trigger_length=architecture.trigger_region_length,
            toehold_length=architecture.toehold_length,
            upper_length=architecture.upper_stem_length,
            lower_length=architecture.stem_length + architecture.upper_stem2_length,
        )
    return tr(
        "Région reconnue du trigger : {trigger_length} nt · toehold : {toehold_length} nt · tige : {stem_length} nt",
        trigger_length=architecture.trigger_region_length,
        toehold_length=architecture.toehold_length,
        stem_length=architecture.stem_length,
    )
