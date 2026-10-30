"""Open system applications without exposing bundled Linux shared libraries."""

from __future__ import annotations

import os
import subprocess
import sys
import webbrowser

from aptaswitch_core.platform_paths import is_wsl, wsl_path


def external_process_environment() -> dict[str, str]:
    """Restore the library search path from before the PyInstaller bootloader.

    Only the child's environment is changed. Bundled scientific helpers
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


def open_directory(path: str | os.PathLike[str]) -> None:
    """Open a folder with the platform launcher and the system environment."""
    if is_wsl():
        subprocess.Popen(
            ["explorer.exe", wsl_path(os.path.abspath(path), to_windows=True)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            env=external_process_environment(),
        )
    elif sys.platform == "win32":
        os.startfile(os.fspath(path))
    else:
        subprocess.run(
            ["open" if sys.platform == "darwin" else "xdg-open", os.fspath(path)],
            check=True, capture_output=True, env=external_process_environment(),
        )


def open_browser(url: str) -> bool:
    """Open a URL with the system browser; never alter the app's environment."""
    if is_wsl():
        try:
            subprocess.Popen(
                ["explorer.exe", url], env=external_process_environment(),
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            return True
        except OSError:
            return False
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
