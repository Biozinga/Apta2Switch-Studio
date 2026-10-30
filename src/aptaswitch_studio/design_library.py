"""Portable, named design folders, independent of the scientific run format.

A design owns its inputs and calculation artifacts. Its small JSON manifest is
replaced atomically; imports are fully staged and validated before publication.
Overwriting a design keeps the previous folder in the library's ``.backups``.
"""

from __future__ import annotations

import io
import json
import os
import shutil
import stat
import tempfile
import unicodedata
import uuid
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Iterator

from aptaswitch_core.localization import tr
from aptaswitch_core.settings import load_settings
from aptaswitch_core.run_storage import default_export_dir


FORMAT = "aptaswitch-design"
VERSION = 1
MANIFEST = "design.json"
MAX_MANIFEST_BYTES = 4 * 1024 * 1024
MAX_FILE_BYTES = 128 * 1024 * 1024
MAX_TOTAL_BYTES = 512 * 1024 * 1024
MAX_FILES = 10000
# Scientific data is portable; user-supplied engines and development files are
# deliberately outside a design archive, including when copied into it manually.
_EXCLUDED = frozenset({
    ".backups", ".git", ".venv", "__pycache__", "__macosx", ".ds_store", ".library.lock", ".saved-switch-library.lock",
    ".local_deps", "site-packages", "dist-packages", "oxdna-runtime", "rna-score-runtime",
})


@dataclass(frozen=True)
class DesignProject:
    path: Path
    name: str
    created_at: str
    updated_at: str
    state: dict[str, Any]


class DesignExistsError(ValueError):
    """The caller must obtain an explicit overwrite choice before retrying."""

    def __init__(self, path: Path):
        self.path = path
        super().__init__(tr('Le design « {name} » existe déjà.', name=path.name))


def default_design_root() -> Path:
    settings = load_settings()
    configured = settings.get("design_library_dir") if isinstance(settings, dict) else None
    if isinstance(configured, str) and configured.strip():
        return Path(configured).expanduser().resolve()
    return default_export_dir().parent / "designs"


def design_output_dir(project: DesignProject | Path) -> Path:
    return (project.path if isinstance(project, DesignProject) else Path(project)) / "calculations"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _fold(value: str) -> str:
    return unicodedata.normalize("NFC", value).casefold()


def _name(value: Any) -> str:
    if not isinstance(value, str):
        raise ValueError(tr('Le nom du design doit être un texte.'))
    value = unicodedata.normalize("NFC", value.strip())
    if (not value or value.startswith(".") or value.endswith(".")
            or any(c in '/\\<>:"|?*' or unicodedata.category(c).startswith("C") for c in value)
            or len(value.encode("utf-8")) > 240):
        raise ValueError(tr('Choisissez un nom de design sans séparateur de dossier ni caractère spécial.'))
    if _fold(value.split(".")[0]) in {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)), *(f"lpt{i}" for i in range(1, 10))}:
        raise ValueError(tr('Ce nom de design est réservé au système.'))
    return value


def _state(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(tr('Les données du design doivent être un objet JSON.'))
    try:
        encoded = json.dumps(value, ensure_ascii=False, allow_nan=False)
        if len(encoded.encode("utf-8")) > MAX_MANIFEST_BYTES - 2048:
            raise ValueError(tr('Les données du design sont trop volumineuses.'))
        return json.loads(encoded)
    except (TypeError, OverflowError, RecursionError, UnicodeError) as exc:
        raise ValueError(tr('Les données du design ne peuvent pas être enregistrées en JSON.')) from exc


def _manifest(name: str, state: dict, created_at: str | None = None) -> dict:
    now = _now()
    return {"format": FORMAT, "version": VERSION, "name": _name(name),
            "created_at": created_at or now, "updated_at": now, "state": _state(state)}


def _read_manifest(path: Path) -> dict:
    if path.is_symlink() or not path.is_file():
        raise ValueError(tr('Le dossier ne contient pas de manifeste design.json valide.'))
    try:
        if path.stat().st_size > MAX_MANIFEST_BYTES:
            raise ValueError(tr('Le manifeste du design est trop volumineux.'))
        with path.open("rb") as source:
            raw = source.read(MAX_MANIFEST_BYTES + 1)
        if len(raw) > MAX_MANIFEST_BYTES:
            raise ValueError(tr('Le manifeste du design est trop volumineux.'))
        data = json.loads(raw.decode("utf-8-sig"))
        if not isinstance(data, dict) or data.get("format") != FORMAT or type(data.get("version")) is not int or data["version"] != VERSION:
            raise ValueError(tr('Le format ou la version de ce design est inconnu.'))
        data["name"] = _name(data.get("name"))
        for key in ("created_at", "updated_at"):
            value = data.get(key)
            if not isinstance(value, str) or datetime.fromisoformat(value).tzinfo is None:
                raise ValueError(tr('Les dates du design sont invalides.'))
        data["state"] = _state(data.get("state"))
        return data
    except (OSError, UnicodeError, json.JSONDecodeError, RecursionError, TypeError) as exc:
        raise ValueError(tr('Impossible de lire le manifeste du design : {error}', error=str(exc))) from exc


def _folder(path: Path) -> Path:
    path = Path(path).expanduser()
    if path.name == MANIFEST and path.is_file():
        if path.is_symlink():
            raise ValueError(tr('Le manifeste du design ne peut pas être un lien symbolique.'))
        path = path.parent
    if path.is_symlink() or not path.is_dir():
        raise ValueError(tr('Sélectionnez un dossier de design valide.'))
    return path.resolve()


def load_design(path: Path) -> DesignProject:
    folder = _folder(path)
    data = _read_manifest(folder / MANIFEST)
    return DesignProject(folder, data["name"], data["created_at"], data["updated_at"], data["state"])


def list_designs(root: Path) -> list[DesignProject]:
    root = Path(root).expanduser()
    if not root.is_dir():
        return []
    projects = []
    for folder in root.iterdir():
        if folder.name.startswith(".") or folder.is_symlink() or not folder.is_dir():
            continue
        try:
            projects.append(load_design(folder))
        except (ValueError, OSError):
            continue
    return sorted(projects, key=lambda project: (project.updated_at, _fold(project.name)), reverse=True)


@contextmanager
def _library_lock(root: Path) -> Iterator[None]:
    # macOS and Linux are the supported desktop targets. A library lock prevents
    # concurrent create/import/save requests from bypassing collision checks.
    import fcntl

    flags = os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(root / ".library.lock", flags, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def _write_manifest(folder: Path, data: dict) -> None:
    content = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    descriptor, temporary = tempfile.mkstemp(prefix=".design-", suffix=".json", dir=folder)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as target:
            target.write(content)
            target.flush()
            os.fsync(target.fileno())
        os.replace(temporary, folder / MANIFEST)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _collision(root: Path, name: str) -> Path | None:
    matches = [path for path in root.iterdir() if _fold(path.name) == _fold(name)]
    if len(matches) > 1:
        raise ValueError(tr('Plusieurs dossiers portent ce nom de design. Renommez-les avant de continuer.'))
    return matches[0] if matches else None


def _publish(staging: Path, root: Path, name: str, overwrite: bool) -> DesignProject:
    with _library_lock(root):
        existing = _collision(root, name)
        if existing is not None and not overwrite:
            raise DesignExistsError(existing)
        if existing is not None and (existing.is_symlink() or not existing.is_dir()):
            raise ValueError(tr('Un fichier ou un lien porte déjà ce nom. Choisissez un autre nom de design.'))
        destination = root / name
        backup = None
        if existing is not None:
            backups = root / ".backups"
            if backups.is_symlink():
                raise ValueError(tr('Le dossier de sauvegarde ne peut pas être un lien symbolique.'))
            backups.mkdir(exist_ok=True)
            backup = backups / f"{existing.name[:40]}_{datetime.now(timezone.utc):%Y%m%dT%H%M%S}_{uuid.uuid4().hex[:8]}"
            existing.rename(backup)
        try:
            staging.rename(destination)
        except OSError:
            if backup is not None:
                backup.rename(existing)
            raise
    return load_design(destination)


def create_design(root: Path, name: str, *, overwrite: bool = False, state: dict | None = None) -> DesignProject:
    data = _manifest(name, {} if state is None else state)
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".design-staging-", dir=root))
    try:
        (staging / "sequences").mkdir()
        (staging / "calculations").mkdir()
        _write_manifest(staging, data)
        return _publish(staging, root, data["name"], overwrite)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def save_design(path: Path, state: dict) -> DesignProject:
    # Validate before touching the filesystem; a bad UI snapshot cannot replace
    # the last readable manifest.
    state = _state(state)
    folder = _folder(path)
    with _library_lock(folder.parent):
        previous = _read_manifest(folder / MANIFEST)
        if previous["state"] == state:
            return DesignProject(folder, previous["name"], previous["created_at"], previous["updated_at"], previous["state"])
        _write_manifest(folder, _manifest(previous["name"], state, previous["created_at"]))
    return load_design(folder)


def _excluded(parts: tuple[str, ...]) -> bool:
    return any(_fold(part) in _EXCLUDED
               or _fold(part).endswith(("_user_site", "-runtime", "_runtime"))
               for part in parts)


def _tree_files(folder: Path) -> list[tuple[Path, str]]:
    files = []
    total = 0
    for current, directories, filenames in os.walk(folder, followlinks=False):
        current_path = Path(current)
        directories[:] = sorted(directory for directory in directories
                                if not (current_path / directory).is_symlink()
                                and not _excluded((directory,)))
        for filename in sorted(filenames):
            source = current_path / filename
            relative = source.relative_to(folder).as_posix()
            if _excluded(tuple(PurePosixPath(relative).parts)) or source.is_symlink():
                continue
            details = source.stat()
            if not stat.S_ISREG(details.st_mode):
                continue
            total += details.st_size
            if details.st_size > MAX_FILE_BYTES or total > MAX_TOTAL_BYTES or len(files) >= MAX_FILES:
                raise ValueError(tr('Le design dépasse la taille maximale autorisée pour un import ou un export.'))
            files.append((source, relative))
    return files


def export_design(path: Path) -> bytes:
    project = load_design(path)
    output = io.BytesIO()
    with _library_lock(project.path.parent):
        files = _tree_files(project.path)
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for source, relative in files:
                archive.write(source, relative)
    return output.getvalue()


def _archive_members(archive: zipfile.ZipFile) -> list[tuple[zipfile.ZipInfo, Path]]:
    infos = archive.infolist()
    if len(infos) > MAX_FILES:
        raise ValueError(tr('L’archive contient trop de fichiers.'))
    members = []
    seen: set[str] = set()
    file_paths: set[str] = set()
    total = 0
    for info in infos:
        name = info.orig_filename
        parts = name.rstrip("/").split("/")
        if (not name or any(part in ("", ".", "..") for part in parts)
                or "\\" in name or ":" in name or "\x00" in name
                or any(any(unicodedata.category(c).startswith("C") for c in part) for part in parts)):
            raise ValueError(tr('L’archive contient un chemin de fichier invalide.'))
        mode = info.external_attr >> 16
        kind = stat.S_IFMT(mode)
        if stat.S_ISLNK(mode) or kind not in (0, stat.S_IFREG, stat.S_IFDIR):
            raise ValueError(tr('Les liens symboliques et fichiers spéciaux ne sont pas acceptés dans un design.'))
        if info.flag_bits & 1:
            raise ValueError(tr('Les archives chiffrées ne sont pas acceptées.'))
        key = _fold("/".join(parts))
        if key in seen:
            raise ValueError(tr('L’archive contient plusieurs fichiers portant le même nom.'))
        seen.add(key)
        if not info.is_dir():
            file_paths.add(key)
        if info.file_size > MAX_FILE_BYTES:
            raise ValueError(tr('Un fichier de l’archive est trop volumineux.'))
        total += info.file_size
        if total > MAX_TOTAL_BYTES:
            raise ValueError(tr('L’archive décompressée est trop volumineuse.'))
        members.append((info, PurePosixPath(*parts)))
    for _, relative in members:
        if any(_fold(str(parent)) in file_paths for parent in relative.parents if str(parent) != "."):
            raise ValueError(tr('L’archive contient un conflit entre un fichier et un dossier.'))
    manifests = [relative for info, relative in members if not info.is_dir() and relative.name == MANIFEST]
    prefix = None
    if PurePosixPath(MANIFEST) not in manifests:
        candidates = [relative for relative in manifests if len(relative.parts) == 2]
        if len(candidates) != 1:
            raise ValueError(tr('L’archive ne contient pas un seul dossier de design valide.'))
        prefix = candidates[0].parts[0]
        if any(relative.parts[0] != prefix and not _excluded(relative.parts) for _, relative in members):
            raise ValueError(tr('L’archive doit contenir uniquement le dossier du design.'))
    selected = []
    for info, relative in members:
        parts = relative.parts
        if prefix:
            if parts[0] != prefix:
                continue
            parts = parts[1:]
        if parts and not _excluded(parts):
            selected.append((info, Path(*parts)))
    return selected


def _extract_archive(source: Path | bytes, staging: Path) -> None:
    if (len(source) if isinstance(source, bytes) else source.stat().st_size) > MAX_TOTAL_BYTES:
        raise ValueError(tr('L’archive du design est trop volumineuse.'))
    stream = io.BytesIO(source) if isinstance(source, bytes) else source
    try:
        with zipfile.ZipFile(stream) as archive:
            members = _archive_members(archive)
            for info, relative in members:
                destination = staging / relative
                if info.is_dir():
                    destination.mkdir(parents=True, exist_ok=True)
                    continue
                destination.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(info) as incoming, destination.open("xb") as outgoing:
                    copied = 0
                    while chunk := incoming.read(1024 * 1024):
                        copied += len(chunk)
                        if copied > info.file_size or copied > MAX_FILE_BYTES:
                            raise ValueError(tr('La taille d’un fichier de l’archive est invalide.'))
                        outgoing.write(chunk)
    except (zipfile.BadZipFile, NotImplementedError, RuntimeError, EOFError) as exc:
        raise ValueError(tr('Impossible d’ouvrir cette archive de design : {error}', error=str(exc))) from exc


def import_design(source: Path | bytes, root: Path, name: str | None = None, *, overwrite: bool = False) -> DesignProject:
    """Validate a complete copy before touching any existing named design."""
    if name is not None:
        name = _name(name)
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".design-staging-", dir=root))
    try:
        if not isinstance(source, bytes):
            source = Path(source).expanduser()
        if not isinstance(source, bytes) and (source.is_dir() or source.name == MANIFEST):
            project = load_design(source)
            # Enumerate first: even importing a folder containing the library
            # must not recursively copy the new staging folder into itself.
            if staging.is_relative_to(project.path):
                raise ValueError(tr('Importez ce design dans une bibliothèque située hors de son dossier.'))
            files = _tree_files(project.path)
            for incoming, relative in files:
                destination = staging / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(incoming, destination, follow_symlinks=False)
        else:
            _extract_archive(source, staging)
        data = _read_manifest(staging / MANIFEST)
        # Importing is a transaction: reject malformed known UI fields before
        # a confirmed replacement moves the previous design into its backup.
        from .design_session import validate_design_state
        validate_design_state(data["state"])
        if name is not None:
            data["name"] = name
        data["updated_at"] = _now()
        for folder in ("calculations", "sequences"):
            (staging / folder).mkdir(exist_ok=True)
        _write_manifest(staging, data)
        return _publish(staging, root, data["name"], overwrite)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
