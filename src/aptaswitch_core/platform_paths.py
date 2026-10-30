"""Translate Windows paths at the boundary of the Linux backend on WSL."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import sys

from .localization import tr


def is_wsl() -> bool:
    if not sys.platform.startswith("linux"):
        return False
    if os.environ.get("WSL_DISTRO_NAME"):
        return True
    try:
        return "microsoft" in Path("/proc/sys/kernel/osrelease").read_text().lower()
    except OSError:
        return False


def wsl_path(path: str | os.PathLike[str], *, to_windows: bool = False) -> str:
    """Convert a path using WSL's own mount configuration, without a shell."""
    environment = os.environ.copy()
    if getattr(sys, "frozen", False):
        # wslpath is a system binary, independent of PyInstaller's libraries.
        environment["LD_LIBRARY_PATH"] = environment.get("LD_LIBRARY_PATH_ORIG", "")
    try:
        result = subprocess.run(
            ["wslpath", "-w" if to_windows else "-u", os.fspath(path)],
            check=True, capture_output=True, text=True, timeout=5, env=environment,
        )
        converted = result.stdout.rstrip("\r\n")
        if not converted or "\n" in converted or "\r" in converted or "\x00" in converted:
            raise ValueError("Invalid wslpath result")
        if not to_windows and not converted.startswith("/"):
            raise ValueError("wslpath did not return an absolute Linux path")
        return converted
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        raise ValueError(tr("Impossible de convertir ce chemin Windows. Vérifiez que son lecteur est accessible dans WSL.")) from exc


def normalize_user_path(raw_path: str | os.PathLike[str]) -> Path:
    """Accept pasted Windows drive/WSL paths only when running inside WSL.

    Network shares and another distribution's private files are intentionally
    rejected: external engines must be accessible in this app's distribution.
    Normal POSIX paths retain their meaning on macOS and ordinary Linux.
    """
    value = os.fspath(raw_path)
    if not is_wsl():
        return Path(value).expanduser()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1]
    if "\x00" in value or "\n" in value or "\r" in value:
        raise ValueError(tr("Le chemin ne doit pas contenir de retour à la ligne ni de caractère nul."))
    if value.startswith(("\\\\", "//")):
        parts = value.replace("\\", "/").split("/")
        distro = os.environ.get("WSL_DISTRO_NAME", "")
        if (len(parts) < 4 or parts[2].lower() not in {"wsl$", "wsl.localhost"}
                or not distro or parts[3].casefold() != distro.casefold()):
            raise ValueError(tr("Utilisez un lecteur Windows local ou un chemin de la distribution WSL de cette application. Les autres partages réseau ne sont pas pris en charge."))
        # Windows treats the share/distro name case-insensitively, whereas
        # older wslpath versions require the registered spelling.
        canonical = "\\\\" + parts[2].lower() + "\\" + distro + "\\" + "\\".join(parts[4:])
        return Path(wsl_path(canonical))
    if re.match(r"^[A-Za-z]:[\\/]", value):
        return Path(wsl_path(value))
    if re.match(r"^[A-Za-z]:", value):
        raise ValueError(tr("Indiquez un chemin Windows absolu, par exemple C:\\Users\\nom\\Downloads."))
    return Path(value).expanduser()
