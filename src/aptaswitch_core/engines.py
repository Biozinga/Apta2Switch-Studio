"""Engine dispatch for Apta2Switch-Studio."""

from __future__ import annotations

from .localization import tr

from .nupack_engine import run_nupack_engine


def run_design(request, progress_callback=None, cancel_check=None):
    if request.engine != "nupack":
        raise ValueError(tr('Seul le moteur NUPACK est pris en charge. Les runs à blanc ont été supprimés.'))
    return run_nupack_engine(request, progress_callback=progress_callback, cancel_check=cancel_check)
