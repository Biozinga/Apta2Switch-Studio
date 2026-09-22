"""High-level run workflow."""

from __future__ import annotations

from .localization import tr

from .engines import run_design
from .exports import export_run


def run_and_export(request, progress_callback=None, cancel_check=None):
    result = run_design(request, progress_callback=progress_callback, cancel_check=cancel_check)
    if progress_callback:
        progress_callback(len(result.candidates), request.trials, tr('Écriture des exports…'))
    exported = export_run(result)
    if progress_callback:
        progress_callback(len(result.candidates), request.trials, tr('Exports terminés.'))
    return result, exported
