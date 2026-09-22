"""Shared locations and names for design and trigger-extension runs."""

from __future__ import annotations

from .localization import tr

import sys
import uuid
from datetime import datetime
from pathlib import Path


RUN_FOLDERS = {"design": "designs", "extension": "extensions"}


def _project_directory(path: Path) -> Path | None:
    """Identify this app's checkout, rather than an unrelated Python project."""

    for parent in path.parents:
        if (parent / "pyproject.toml").is_file() and (parent / "src" / "aptaswitch_studio").is_dir():
            return parent
    return None


def application_directory() -> Path:
    """Locate the project or the folder containing the installed application.

    A macOS bundle is signed application code, so its writable run library lives
    beside the ``.app``, never inside ``Contents`` or PyInstaller's resources.
    """

    if getattr(sys, "frozen", False):
        executable = Path(sys.executable).resolve()
        project = _project_directory(executable)
        if project is not None:
            return project
        for parent in executable.parents:
            if parent.suffix.lower() == ".app":
                return parent.parent
        return executable.parent
    return _project_directory(Path(__file__).resolve()) or Path.cwd()


def default_export_dir(app_dir: Path | None = None) -> Path:
    """Return the shared run library in the application's directory."""

    return (Path(app_dir) if app_dir is not None else application_directory()) / "runs"


def legacy_export_dir(home_dir: Path | None = None) -> Path:
    """Return the previous location, also used when the app folder is read-only."""

    return (Path(home_dir) if home_dir is not None else Path.home()) / "Apta2Switch-Studio Runs"


def run_output_directory(root: Path, kind: str) -> Path:
    """Select a run category without nesting existing category directories."""

    try:
        folder = RUN_FOLDERS[kind]
    except KeyError as exc:
        raise ValueError(tr('Type de run inconnu : {kind}', kind=kind)) from exc
    root = Path(root).expanduser()
    if root.name in RUN_FOLDERS.values():
        root = root.parent
    return root / folder


def new_run_id(kind: str, started: datetime | None = None) -> str:
    """Create a readable ID that also separates runs started in one second."""

    if kind not in RUN_FOLDERS:
        raise ValueError(tr('Type de run inconnu : {kind}', kind=kind))
    started = started or datetime.now().astimezone()
    return f"{started.strftime('%Y%m%d_%H%M%S')}_{kind}_{uuid.uuid4().hex[:6]}"


def new_run_directory(
    root: Path, kind: str, started: datetime | None = None,
) -> tuple[str, Path]:
    """Return the common run ID and final export folder without writing files."""

    run_id = new_run_id(kind, started)
    return run_id, run_output_directory(root, kind) / run_id
