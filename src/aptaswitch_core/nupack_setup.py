"""Validation and settings helpers for user-provided NUPACK paths."""

from __future__ import annotations

from .localization import tr

import importlib.util
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from packaging.tags import sys_tags
from packaging.utils import InvalidWheelFilename, parse_wheel_filename

from .run_storage import application_directory


SETTINGS_DIR = Path.home() / ".aptaswitch-studio"
SETTINGS_PATH = SETTINGS_DIR / "settings.json"


@dataclass(frozen=True)
class NupackValidationResult:
    ok: bool
    path: str
    message: str
    version: str = ""
    import_path: str = ""


PROBE_FLAG = "--aptaswitch-nupack-probe"


def handle_nupack_probe(argv: list[str] | None = None) -> bool:
    """Handle the bundled app's private, isolated NUPACK validation process.

    A frozen executable is not a Python CLI: using ``-c`` or ``-m pip`` would
    launch the GUI again. The launcher dispatches this flag before Streamlit.
    A result file works even for windowed builds without stdout/stderr.
    """
    arguments = sys.argv[1:] if argv is None else argv
    if not arguments or arguments[0] != PROBE_FLAG:
        return False
    if len(arguments) != 3:
        raise SystemExit(2)
    import_path, result_file = map(Path, arguments[1:])
    try:
        sys.path.insert(0, str(import_path))
        import nupack

        # Do not accidentally validate some unrelated installation on sys.path.
        Path(nupack.__file__).resolve().relative_to(import_path.resolve())
        model = nupack.Model(material="rna", celsius=37)
        nupack.config.threads = 1
        result = nupack.mfe(strands=["GGGAAACCC"], model=model)
        if not result:
            raise RuntimeError("NUPACK did not return a structure.")
        payload = {"ok": True, "version": str(getattr(nupack, "__version__", "unknown"))}
    except Exception as exc:
        payload = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
    result_file.write_text(json.dumps(payload), encoding="utf-8")
    return True


def _python_import_check(import_path: Path) -> NupackValidationResult:
    with tempfile.TemporaryDirectory(prefix="aptaswitch-nupack-probe-") as temporary:
        result_file = Path(temporary) / "result.json"
        arguments = [PROBE_FLAG, str(import_path), str(result_file)]
        if getattr(sys, "frozen", False):
            command = [sys.executable, *arguments]
        else:
            source_root = str(Path(__file__).resolve().parents[1])
            code = (
                "import sys; "
                f"sys.path.insert(0, {source_root!r}); "
                "from aptaswitch_core.nupack_setup import handle_nupack_probe; "
                "handle_nupack_probe()"
            )
            command = [sys.executable, "-I", "-c", code, *arguments]
        try:
            completed = subprocess.run(
                command, text=True, capture_output=True, timeout=45, check=False,
            )
            data = json.loads(result_file.read_text(encoding="utf-8")) if result_file.is_file() else {}
        except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
            return NupackValidationResult(False, str(import_path), str(exc), import_path=str(import_path))
        if completed.returncode != 0 or not data.get("ok"):
            message = data.get("error") or completed.stderr or completed.stdout or tr('Impossible d’importer NUPACK.')
            return NupackValidationResult(False, str(import_path), message.strip(), import_path=str(import_path))
    return NupackValidationResult(
        True,
        str(import_path),
        tr('Import NUPACK validé depuis {v0}.', v0=import_path),
        version=str(data.get("version", "unknown")),
        import_path=str(import_path),
    )


def candidate_import_paths(path: Path) -> list[Path]:
    if path.is_file():
        return []
    candidates = [path] if (path / "nupack" / "__init__.py").is_file() else []
    if (path / ".local_deps").is_dir():
        candidates.append(path / ".local_deps")
    if path.name == "nupack" and (path / "__init__.py").is_file():
        candidates.append(path.parent)
    return list(dict.fromkeys(candidates))


def compatible_nupack_wheels(path: Path) -> list[Path]:
    """Find official wheels in a download, its package folder or the drop folder."""
    ranks = {tag: index for index, tag in enumerate(sys_tags())}
    found = set()
    for pattern in ("*.whl", "package/*.whl", "nupack*/*.whl", "nupack*/package/*.whl"):
        found.update(path.glob(pattern))
    candidates = []
    for wheel in found:
        try:
            name, version, _, tags = parse_wheel_filename(wheel.name)
        except InvalidWheelFilename:
            continue
        matching = tags & ranks.keys()
        if name == "nupack" and matching:
            candidates.append((version, -min(ranks[tag] for tag in matching), wheel))
    return [item[2] for item in sorted(candidates, reverse=True)]


def _incompatible_wheel_message() -> str:
    return tr(
        'Aucun fichier NUPACK compatible avec cette application (Python {python}, {system}, {machine}). '
        'Téléchargez le paquet NUPACK pour votre système et déposez son dossier décompressé ici.',
        python=f"{sys.version_info.major}.{sys.version_info.minor}",
        system=platform.system(), machine=platform.machine(),
    )


def validate_nupack_path(raw_path: str) -> NupackValidationResult:
    path = Path(raw_path).expanduser().resolve()
    if not path.exists():
        return NupackValidationResult(False, str(path), tr('Ce chemin n’existe pas.'))

    if sys.platform == "win32":
        return NupackValidationResult(False, str(path), tr('NUPACK ne prend pas en charge Windows en natif. Les calculs de cette bêta sont disponibles sur macOS et Linux.'))

    if path.is_file() and path.suffix == ".whl":
        return validate_nupack_wheel(path)

    if path.is_file():
        return NupackValidationResult(False, str(path), tr('Sélectionnez un dossier, un dossier .local_deps ou un fichier .whl.'))

    last_error = tr('Aucun paquet NUPACK importable trouvé dans ce dossier.')
    for import_path in candidate_import_paths(path):
        result = _python_import_check(import_path)
        if result.ok:
            return result
        last_error = result.message
    for wheel in compatible_nupack_wheels(path):
        result = validate_nupack_wheel(wheel)
        if result.ok:
            return NupackValidationResult(True, str(path), result.message, result.version, result.import_path)
        last_error = result.message
    if last_error == tr('Aucun paquet NUPACK importable trouvé dans ce dossier.'):
        last_error = _incompatible_wheel_message()
    return NupackValidationResult(False, str(path), last_error)


def validate_nupack_wheel(path: Path) -> NupackValidationResult:
    if path.suffix != ".whl":
        return NupackValidationResult(False, str(path), tr('Un fichier .whl est attendu.'))
    try:
        name, _, _, tags = parse_wheel_filename(path.name)
    except InvalidWheelFilename:
        return NupackValidationResult(False, str(path), tr('Le nom du fichier wheel NUPACK est invalide.'))
    if name != "nupack":
        return NupackValidationResult(False, str(path), tr('Le nom du fichier wheel ne semble pas correspondre à NUPACK.'))
    if not tags.intersection(sys_tags()):
        return NupackValidationResult(False, str(path), _incompatible_wheel_message())
    try:
        with path.open("rb") as source:
            digest = hashlib.file_digest(source, "sha256").hexdigest()[:16]
        cache_root = SETTINGS_DIR / "nupack_user_site"
        cache_root.mkdir(parents=True, exist_ok=True)
        wheel_site = cache_root / f"{path.stem}-{digest}"
        if wheel_site.is_dir():
            result = _python_import_check(wheel_site)
            if result.ok:
                return NupackValidationResult(True, str(path), result.message, result.version, str(wheel_site))
        with tempfile.TemporaryDirectory(prefix=".install-", dir=cache_root) as temporary:
            target = Path(temporary) / "site"
            target.mkdir()
            _unpack_wheel(path, target)
            result = _python_import_check(target)
            if not result.ok:
                return NupackValidationResult(False, str(path), result.message)
            if wheel_site.exists():
                shutil.rmtree(wheel_site)
            target.rename(wheel_site)
        return NupackValidationResult(
            True, str(path), tr('Fichier wheel NUPACK validé : {v0}.', v0=path.name),
            version=result.version, import_path=str(wheel_site),
        )
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        return NupackValidationResult(False, str(path), str(exc))


def _unpack_wheel(wheel: Path, target: Path) -> None:
    """Install a binary wheel locally, without pip, networking or host Python."""
    with zipfile.ZipFile(wheel) as archive:
        for info in archive.infolist():
            relative = PurePosixPath(info.filename)
            if relative.is_absolute() or ".." in relative.parts or "\\" in info.filename:
                raise ValueError(tr('Le fichier wheel contient un chemin non valide.'))
            parts = relative.parts
            if parts and parts[0].endswith(".data"):
                if len(parts) < 3 or parts[1] not in {"purelib", "platlib"}:
                    continue
                relative = PurePosixPath(*parts[2:])
            output = target.joinpath(*relative.parts)
            if info.is_dir():
                output.mkdir(parents=True, exist_ok=True)
                continue
            output.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as source, output.open("wb") as destination:
                shutil.copyfileobj(source, destination)
        if not (target / "nupack" / "__init__.py").is_file():
            raise ValueError(tr('Ce fichier wheel ne contient pas le paquet NUPACK.'))


def load_settings() -> dict:
    if not SETTINGS_PATH.exists():
        return {}
    try:
        return json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def save_settings(settings: dict) -> None:
    SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    SETTINGS_PATH.write_text(json.dumps(settings, indent=2), encoding="utf-8")


def save_nupack_settings(validation: NupackValidationResult, license_confirmed: bool) -> None:
    settings = load_settings()
    settings["nupack"] = {
        "path": validation.path,
        "import_path": validation.import_path,
        "version": validation.version,
        "license_confirmed": bool(license_confirmed),
    }
    save_settings(settings)


def configured_nupack_import_path() -> str | None:
    settings = load_settings().get("nupack", {})
    if not settings.get("license_confirmed"):
        return None
    return settings.get("import_path") or settings.get("path")


def nupack_drop_directory(create: bool = False) -> Path:
    """The visible portable drop folder, with a writable user-folder fallback."""
    folder = application_directory() / "nupack"
    fallback = SETTINGS_DIR / "nupack"
    if create:
        try:
            folder.mkdir(parents=True, exist_ok=True)
        except OSError:
            fallback.mkdir(parents=True, exist_ok=True)
            return fallback
    elif not os.access(folder if folder.exists() else folder.parent, os.W_OK):
        return fallback
    return folder


def default_nupack_search_paths() -> list[Path]:
    """Folders that may already contain a user-provided NUPACK install."""
    candidates = []
    for base in (application_directory(), SETTINGS_DIR, Path.cwd()):
        for name in (".local_deps", "nupack"):
            folder = base / name
            if folder.exists():
                candidates.append(folder)
    return list(dict.fromkeys(candidates))


def detect_nupack() -> NupackValidationResult | None:
    """Return an already-available NUPACK install, if any, without prompting.

    Detection order: NUPACK importable in the current interpreter, then the
    project-local ``.local_deps`` / ``nupack`` folders. Returns ``None`` when no
    NUPACK is present, so the manual selection flow stays available.
    """
    try:
        spec = importlib.util.find_spec("nupack")
    except (ImportError, ValueError):
        spec = None
    if spec is not None and spec.submodule_search_locations:
        package_dir = Path(list(spec.submodule_search_locations)[0])
        result = _python_import_check(package_dir.parent)
        if result.ok:
            return result

    for candidate in default_nupack_search_paths():
        result = validate_nupack_path(str(candidate))
        if result.ok:
            return result
    return None
