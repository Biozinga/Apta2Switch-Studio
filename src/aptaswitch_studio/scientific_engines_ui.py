"""Scientific engine identity and availability messages for the interface."""
from __future__ import annotations

from dataclasses import asdict
import json


from aptaswitch_core.localization import tr
from aptaswitch_core.thermodynamics.setup import (
    configured_library_path, get_backend,
)


NATIVE_MODEL_HELP = (
    "ARN isolé : algorithmes ViennaRNA 2.7.2. Complexes ADN–ARN : moteur hybride "
    "expérimental Apta2Switch. Les paramètres ARN et la correction saline sont "
    "communs ; ce ne sont pas les réglages par défaut de RNAfold."
)
NATIVE_UNAVAILABLE = "Le moteur thermodynamique intégré est indisponible. Réinstallez le paquet complet de l’application."


_DEFAULT_PATH = object()


def native_backend_identity(library_path=_DEFAULT_PATH) -> str | None:
    """Include binary, parameters and ensemble conventions, not just a path."""
    path = configured_library_path() if library_path is _DEFAULT_PATH else library_path
    if not path:
        return None
    try:
        return json.dumps(asdict(get_backend(path).provenance), sort_keys=True, separators=(",", ":"))
    except (OSError, ValueError, RuntimeError):
        return None


def result_engine_label(result: dict) -> str:
    """Attribute archived calculations to their recorded producer."""
    model = result.get("model") if isinstance(result.get("model"), dict) else {}
    provenance = result.get("thermodynamics") or {}
    if not provenance and isinstance(result.get("model"), dict):
        model = result["model"]
        if model.get("backend"):
            provenance = model
    if provenance:
        return str(provenance.get("backend") or provenance.get("engine") or tr("Moteur non renseigné"))
    request = result.get("request") if isinstance(result.get("request"), dict) else {}
    if request.get("engine") == "native-hybrid":
        return tr("Apta2Switch hybride expérimental")
    return str(request.get("engine") or tr("Moteur non renseigné"))


