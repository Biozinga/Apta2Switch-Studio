"""Pinned offline molecular-viewer JavaScript shared by oxDNA movie rendering."""
from __future__ import annotations

from functools import lru_cache
import hashlib
from pathlib import Path
import re


VENDOR_DIR = Path(__file__).resolve().parent / "assets" / "3dmol"
VENDOR_SHA256 = "f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016"


@lru_cache(maxsize=1)
def bundled_viewer_javascript() -> str:
    """Read one pinned runtime, including in frozen apps; never fetch a CDN."""
    data = (VENDOR_DIR / "3Dmol-min.js").read_bytes()
    if hashlib.sha256(data).hexdigest() != VENDOR_SHA256:
        raise RuntimeError("The bundled 3Dmol.js asset failed its integrity check.")
    # Avoid turning script content into HTML; retain the upstream file unchanged.
    return re.sub(r"</script", r"<\/script", data.decode("utf-8"), flags=re.IGNORECASE)


