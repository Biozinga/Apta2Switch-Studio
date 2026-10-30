"""Native equilibrium populations in trigger-only and aptamer/trigger tubes.

Complex particle fractions and strand fractions are reported separately: one
trigger dimer is one complex but contains two trigger molecules.
"""
from __future__ import annotations

from typing import Callable

from .models import ThermoConditions
from .thermodynamics.tube_reporting import (
    AVOGADRO_CONSTANT, MASS_BALANCE_RELATIVE_TOLERANCE,
    NEGATIVE_ROUNDOFF_RELATIVE_TOLERANCE, TubeAnalysisConditions, _population,
)


def analyze_parasite_tubes(
    import_path: str | None,
    aptamer_dna: str,
    trigger_dna: str,
    thermo: ThermoConditions,
    conditions: TubeAnalysisConditions = TubeAnalysisConditions(),
    *,
    cancel_check: Callable[[], bool] | None = None,
    progress_callback: Callable[[int, int, str], None] | None = None,
    backend=None,
    thermodynamics_library_path: str | None = None,
) -> dict:
    """Solve complete connected-complex catalogs with native partition functions.

    The positional ``import_path`` is retained for compatibility and ignored.
    Select a native library explicitly with ``thermodynamics_library_path`` or
    configure it for the application. Cancellation is checked between calls.
    """
    from .thermodynamics.application_analysis import analyze_parasite_tubes as analyze_native
    from .thermodynamics.setup import get_backend

    if backend is None:
        backend = get_backend(thermodynamics_library_path)
    return analyze_native(
        backend, aptamer_dna, trigger_dna, thermo, conditions,
        standard_state=backend.standard_state,
        cancel_check=cancel_check, progress_callback=progress_callback,
    )
