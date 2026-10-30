"""Native MFE analysis of the supplied aptamer and trigger before extension."""
from __future__ import annotations

from .models import ThermoConditions


def analyze_system(
    import_path: str | None,
    aptamer_dna: str,
    trigger_dna: str,
    thermo: ThermoConditions,
    *,
    backend=None,
    thermodynamics_library_path: str | None = None,
) -> dict:
    """Compute the initial two-strand minimum with the configured native model.

    ``import_path`` remains a positional compatibility placeholder. Historical
    import paths never select an engine or become native library paths.
    """
    from .thermodynamics.application_analysis import analyze_system as analyze_native
    from .thermodynamics.setup import get_backend

    thermo.validate_model()
    if backend is None:
        backend = get_backend(thermodynamics_library_path)
    return analyze_native(backend, aptamer_dna, trigger_dna, thermo)
