"""Open system applications without exposing bundled Linux shared libraries."""

from __future__ import annotations

import os
import subprocess
import sys
import webbrowser


def external_process_environment() -> dict[str, str]:
    """Restore the library search path from before the PyInstaller bootloader.

    Only the child's environment is changed. NUPACK and other bundled helpers
    continue to inherit the app's own library paths.
    """
    environment = os.environ.copy()
    if sys.platform.startswith("linux") and getattr(sys, "frozen", False):
        original = environment.get("LD_LIBRARY_PATH_ORIG")
        if original is None:
            environment.pop("LD_LIBRARY_PATH", None)
        else:
            environment["LD_LIBRARY_PATH"] = original
    return environment


def open_browser(url: str) -> bool:
    """Open a URL with the system browser; never alter the app's environment."""
    if sys.platform.startswith("linux"):
        for command in (["xdg-open", url], ["gio", "open", url]):
            try:
                subprocess.Popen(
                    command, env=external_process_environment(),
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                )
                return True
            except FileNotFoundError:
                continue
            except OSError:
                return False
        return False
    return webbrowser.open(url)
