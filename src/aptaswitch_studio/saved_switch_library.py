"""Portable, named selections of actual switch-design candidates.

A selection is a detached snapshot, not a pointer into a run or today's form.
The original snapshot never changes when context analyses or 3D models are
added. Replacing a named selection requires an explicit overwrite and keeps
the previous *entire* folder, including its analyses, in ``.backups``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import tempfile
import unicodedata
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from aptaswitch_core.localization import tr
from aptaswitch_core.models import DesignCandidate, RunResult
from aptaswitch_core.run_storage import default_export_dir


FORMAT = "aptaswitch-saved-switch"
VERSION = 1
LIBRARY_FOLDER = "tswitch_enregistres"
MANIFEST = "snapshot.json"
MAX_SNAPSHOT_BYTES = 16 * 1024 * 1024


@dataclass(frozen=True)
class SavedSwitch:
    path: Path
    name: str
    saved_at: str
    snapshot: dict[str, Any]

    @property
    def manifest_path(self) -> Path:
        return self.path / MANIFEST

    @property
    def fasta_path(self) -> Path:
        return self.path / "sequences.fasta"

    @property
    def context_path(self) -> Path:
        return self.path / "context.json"

    @property
    def analyses_dir(self) -> Path:
        return self.path / "calculations" / "context"


class SavedSwitchExistsError(ValueError):
    """An existing named selection needs explicit user consent to replace."""

    def __init__(self, path: Path):
        self.path = path
        super().__init__(tr('Le switch enregistré « {name} » existe déjà.', name=path.name))


def _fold(value: str) -> str:
    return unicodedata.normalize("NFC", value).casefold()


def _name(value: Any) -> str:
    if not isinstance(value, str):
        raise ValueError(tr('Le nom du switch doit être un texte.'))
    value = unicodedata.normalize("NFC", value.strip())
    if (not value or value.startswith(".") or value.endswith(".")
            or any(c in '/\\<>:"|?*' or unicodedata.category(c).startswith("C") for c in value)
            or len(value.encode("utf-8")) > 180):
        raise ValueError(tr('Choisissez un nom de switch sans séparateur de dossier ni caractère spécial.'))
    if _fold(value.split(".")[0]) in {
        "con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)),
        *(f"lpt{i}" for i in range(1, 10)),
    }:
        raise ValueError(tr('Ce nom de switch est réservé au système.'))
    return value


def _checked_path(path: Path) -> Path:
    """Reject links before resolution, including a linked library ancestor.

    Canonicalize a trusted temporary-directory location before calling this
    helper on macOS (where ``/var`` and ``/tmp`` themselves are system links).
    No path contained in a snapshot is ever used for filesystem operations.
    """
    path = Path(path).expanduser().absolute()
    if ".." in path.parts:
        raise ValueError(tr('Le chemin du switch enregistré est invalide.'))
    for item in (path, *path.parents):
        if item.is_symlink():
            raise ValueError(tr('Un dossier de switches enregistrés ne peut pas traverser un lien symbolique.'))
    return path


def saved_switch_root(
    design_path: Path | None = None, export_root: Path | None = None,
) -> Path:
    """Keep selections within their design; use the export library otherwise."""
    base = Path(design_path) if design_path is not None else Path(
        default_export_dir() if export_root is None else export_root
    )
    if base.name == LIBRARY_FOLDER:
        return _checked_path(base)
    return _checked_path(base / LIBRARY_FOLDER)


def _json_copy(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(tr('Le switch enregistré doit être un objet JSON.'))
    try:
        encoded = json.dumps(value, ensure_ascii=False, allow_nan=False)
        if len(encoded.encode("utf-8")) > MAX_SNAPSHOT_BYTES:
            raise ValueError(tr('Le switch enregistré est trop volumineux.'))
        return json.loads(encoded)
    except (TypeError, OverflowError, RecursionError, UnicodeError) as exc:
        raise ValueError(tr('Les données du switch ne peuvent pas être enregistrées en JSON.')) from exc


def _sequence(value: Any, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or (not value and not allow_empty) or set(value.upper()) - set("ACGTU"):
        raise ValueError(tr('Le switch enregistré contient une séquence absente ou invalide.'))
    return value


def _rna(sequence: str) -> str:
    return sequence.upper().replace("T", "U")


def _validate_snapshot(value: Any) -> dict[str, Any]:
    # This clone also rejects nonfinite metrics anywhere in extensible metadata.
    data = _json_copy(value)
    if (data.get("format") != FORMAT or type(data.get("version")) is not int
            or data["version"] != VERSION):
        raise ValueError(tr('Le format ou la version du switch enregistré est inconnu.'))
    data["name"] = _name(data.get("name"))
    try:
        identity = data["snapshot_id"]
        if not isinstance(identity, str) or uuid.UUID(hex=identity).hex != identity:
            raise ValueError
        if (not isinstance(data["saved_at"], str)
                or datetime.fromisoformat(data["saved_at"]).tzinfo is None):
            raise ValueError
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(tr('L’identité ou la date du switch enregistré est invalide.')) from exc
    source = data.get("source")
    sequences = data.get("sequences")
    if not isinstance(source, dict) or not isinstance(sequences, dict):
        raise ValueError(tr('Les séquences ou la provenance du switch enregistré sont absentes.'))
    # Reuse the established scientific-run validator without writing any file,
    # reranking candidates, or replacing their saved request with current UI.
    from .run_library import _load_design

    try:
        result = _load_design({
            **source, "request": data.get("request"), "results": [data.get("candidate")],
        }, Path("/saved-switch/snapshot.json"), retained_candidate=True)
    except (KeyError, TypeError, AttributeError, OverflowError, RecursionError) as exc:
        raise ValueError(tr('Les paramètres scientifiques du switch enregistré sont invalides.')) from exc
    candidate = result.candidates[0]
    expected = {
        "switch": candidate.switch_seq,
        "trigger": candidate.trigger_seq,
        "aptamer": result.request.sequences.aptamer_dna,
        "input_trigger": result.request.sequences.trigger_dna,
        "reporter_after_start": result.request.reporter.sequence_after_start_rna,
    }
    for key, sequence in expected.items():
        if sequences.get(key) != sequence:
            raise ValueError(tr('Les séquences du switch enregistré ne correspondent pas au candidat conservé.'))
        _sequence(sequence, allow_empty=key == "reporter_after_start")
    full_reporter = result.request.reporter.full_sequence_after_start_rna
    if sequences.get("reporter_is_complete") is not (full_reporter is not None):
        raise ValueError(tr('Le contexte complet du rapporteur ne correspond pas au design conservé.'))
    if full_reporter is not None and sequences.get("reporter_full_after_start") != full_reporter:
        raise ValueError(tr('Le rapporteur complet enregistré ne correspond pas au design conservé.'))
    original = source.get("original_sequences")
    if original is not None:
        if not isinstance(original, dict):
            raise ValueError(tr('Les séquences originales du design sont invalides.'))
        for key in ("aptamer_dna", "trigger_dna"):
            _sequence(original.get(key))
        if (sequences.get("original_trigger") != original["trigger_dna"]
                or sequences.get("original_trigger_source") != "design_original_sequences"):
            raise ValueError(tr('Le trigger original ne correspond pas au design source.'))
    elif (sequences.get("original_trigger") != expected["input_trigger"]
          or sequences.get("original_trigger_source") != "request_input"):
        raise ValueError(tr('La provenance du trigger original est invalide.'))
    # candidate.switch_seq already contains the designed reporter segment.
    # Never append that segment a second time when exporting the selection.
    reporter = _rna(expected["reporter_after_start"])
    if not _rna(candidate.switch_seq).endswith(reporter):
        raise ValueError(tr('La fin du switch ne correspond pas au rapporteur du design.'))
    if full_reporter is not None:
        expected_transcript = candidate.switch_seq[:len(candidate.switch_seq) - len(reporter)] + _rna(full_reporter)
        if sequences.get("switch_full_transcript") != expected_transcript:
            raise ValueError(tr('Le transcrit complet enregistré ne correspond pas au switch et au rapporteur.'))
    return data


def candidate_snapshot(
    result: RunResult, candidate: DesignCandidate, name: str, *,
    original_sequences: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Freeze one candidate and all the conditions that actually produced it."""
    if candidate not in result.candidates:
        raise ValueError(tr('Ce candidat ne fait pas partie du résultat de design sélectionné.'))
    request = result.request.to_dict()
    source = {
        "run_id": result.run_id, "started_at": result.started_at,
        "completed_at": result.completed_at, "status": result.status,
        "warnings": list(result.warnings),
        "inverse_extension_analysis": result.inverse_extension_analysis,
        "thermodynamics": result.thermodynamics,
        "optimization": dict(result.optimization),
    }
    if result.request.engine == "native-hybrid":
        # A selection retains its own measured trial, not every unrelated
        # optimizer trace from a potentially hundred-thousand-candidate run.
        # The original complete run and its outcome list remain unchanged.
        reports = result.optimization.get("trials", [])
        selected = [report for report in reports if isinstance(report, dict)
                    and report.get("trial") == candidate.trial]
        source["optimization"] = {
            **result.optimization,
            "scope": "retained_candidate",
            "source_run_id": result.run_id,
            "source_trial_count": result.request.trials,
            "retained_trial": candidate.trial,
            "trials": selected,
        }
    # A candidate snapshot is a selection, never a resumable full design run.
    # Keep its optimizer settings and measured outcome without the run cursor
    # or integrity digest, which belong to the original checkpoint only.
    source["optimization"].pop("checkpoint", None)
    if original_sequences is not None:
        if not isinstance(original_sequences, dict):
            raise ValueError(tr('Les séquences originales du design sont invalides.'))
        original_sequences = {
            **original_sequences,
            "aptamer_dna": original_sequences.get("aptamer_dna", original_sequences.get("aptamer")),
            "trigger_dna": original_sequences.get("trigger_dna", original_sequences.get("trigger")),
        }
        source["original_sequences"] = original_sequences
    original_trigger = (original_sequences.get("trigger_dna") if original_sequences is not None
                        else result.request.sequences.trigger_dna)
    full_reporter = result.request.reporter.full_sequence_after_start_rna
    full_sequences = {}
    if full_reporter is not None:
        prefix_length = len(result.request.reporter.sequence_after_start_rna)
        full_sequences = {
            "reporter_full_after_start": full_reporter,
            "switch_full_transcript": candidate.switch_seq[:len(candidate.switch_seq) - prefix_length]
            + _rna(full_reporter),
        }
    return _validate_snapshot({
        "format": FORMAT, "version": VERSION,
        "snapshot_id": uuid.uuid4().hex, "name": _name(name),
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "candidate": candidate.to_dict(), "request": request, "source": source,
        "sequences": {
            "switch": candidate.switch_seq, "trigger": candidate.trigger_seq,
            "aptamer": result.request.sequences.aptamer_dna,
            "input_trigger": result.request.sequences.trigger_dna,
            "original_trigger": original_trigger,
            "original_trigger_source": ("design_original_sequences" if original_sequences is not None
                                        else "request_input"),
            "reporter_after_start": result.request.reporter.sequence_after_start_rna,
            "reporter_is_complete": full_reporter is not None,
            **full_sequences,
        },
    })


def _fasta(snapshot: dict[str, Any]) -> str:
    sequences = snapshot["sequences"]
    entries = [
        ("switch_designed_including_reporter_segment", sequences["switch"]),
        ("trigger_candidate", sequences["trigger"]),
        ("aptamer_design_input", sequences["aptamer"]),
        ("trigger_design_input", sequences["input_trigger"]),
        ("trigger_original", sequences["original_trigger"]),
        ("reporter_segment_after_start_completeness_unspecified", sequences["reporter_after_start"]),
    ]
    if sequences["reporter_is_complete"]:
        entries.extend([
            ("reporter_full_after_start", sequences["reporter_full_after_start"]),
            ("switch_full_transcript", sequences["switch_full_transcript"]),
        ])
    return "".join(
        f">{label}\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"
        for label, sequence in entries if sequence
    )


def _write_file(path: Path, text: str) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as target:
        target.write(text)
        target.flush()
        os.fsync(target.fileno())


@contextmanager
def _library_lock(root: Path) -> Iterator[None]:
    # Both supported desktop platforms provide flock and O_NOFOLLOW.
    import fcntl

    descriptor = os.open(root / ".saved-switch-library.lock", os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK, 0o600)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(tr('Le verrou de la bibliothèque est invalide.'))
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def _collision(root: Path, name: str) -> Path | None:
    matches = [path for path in root.iterdir() if _fold(path.name) == _fold(name)]
    if len(matches) > 1:
        raise ValueError(tr('Plusieurs dossiers portent ce nom de switch. Renommez-les avant de continuer.'))
    return matches[0] if matches else None


def save_switch(
    root: Path, name: str, snapshot: dict[str, Any], *, overwrite: bool = False,
    source_output_dir: Path | None = None,
) -> SavedSwitch:
    """Publish a fully staged snapshot; duplicate names never replace silently."""
    name = _name(name)
    data = _validate_snapshot(snapshot)
    if data["name"] != name:
        raise ValueError(tr('Le nom du dossier doit correspondre au nom du switch enregistré.'))
    content = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if len(content.encode("utf-8")) > MAX_SNAPSHOT_BYTES:
        raise ValueError(tr('Le switch enregistré est trop volumineux.'))
    root = _checked_path(root)
    root.mkdir(parents=True, exist_ok=True)
    _checked_path(root)
    staging = Path(tempfile.mkdtemp(prefix=".switch-staging-", dir=root))
    try:
        _write_file(staging / MANIFEST, content)
        _write_file(staging / "sequences.fasta", _fasta(data))
        (staging / "calculations" / "context").mkdir(parents=True)
        if source_output_dir is not None:
            from .saved_switch_design import copy_saved_design_ensembles

            copy_saved_design_ensembles(data, source_output_dir, staging / "calculations" / "design")
        with _library_lock(root):
            existing = _collision(root, name)
            if existing is not None and not overwrite:
                raise SavedSwitchExistsError(existing)
            if existing is not None and (existing.is_symlink() or not existing.is_dir()):
                raise ValueError(tr('Un fichier ou un lien porte déjà ce nom de switch.'))
            if existing is not None:
                # Candidate workers share this lock during registration so an
                # overwrite cannot move their destination while they run.
                from .saved_switch_ffs_runtime import active_ffs
                from .saved_switch_free_md_runtime import active_free_md
                from .saved_switch_ensemble_runtime import active_context_ensemble

                if active_ffs(existing) is not None or active_free_md(existing) is not None:
                    raise ValueError(tr('Attendez ou annulez la simulation avant de remplacer ce switch enregistré.'))
                if active_context_ensemble(existing) is not None:
                    raise ValueError(tr('Attendez ou annulez l’analyse des familles avant de remplacer ce switch enregistré.'))
            destination = root / name
            backup = None
            if existing is not None:
                backups = _checked_path(root / ".backups")
                backups.mkdir(exist_ok=True)
                backup = backups / f"{existing.name[:40]}_{datetime.now(timezone.utc):%Y%m%dT%H%M%S}_{uuid.uuid4().hex[:8]}"
                existing.rename(backup)
            try:
                staging.rename(destination)
            except OSError:
                if backup is not None:
                    backup.rename(existing)
                raise
            return load_switch(destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    data = {}
    for key, value in pairs:
        if key in data:
            raise ValueError(tr('Le manifeste contient des clés JSON en double.'))
        data[key] = value
    return data


def _read_snapshot(path: Path) -> tuple[Path, Any]:
    """Read the immutable manifest without importing scientific or UI modules."""
    path = _checked_path(path)
    folder = path.parent if path.name == MANIFEST else path
    if not folder.is_dir():
        raise ValueError(tr('Le dossier du switch enregistré est absent.'))
    manifest = _checked_path(folder / MANIFEST)
    descriptor = os.open(manifest, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        details = os.fstat(descriptor)
        if not stat.S_ISREG(details.st_mode) or details.st_size > MAX_SNAPSHOT_BYTES:
            raise ValueError(tr('Le manifeste du switch enregistré est invalide ou trop volumineux.'))
        with os.fdopen(descriptor, "rb", closefd=False) as source:
            raw = source.read(MAX_SNAPSHOT_BYTES + 1)
        if len(raw) > MAX_SNAPSHOT_BYTES:
            raise ValueError(tr('Le manifeste du switch enregistré est trop volumineux.'))
        try:
            data = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique_object)
        except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
            raise ValueError(tr('Le manifeste du switch enregistré contient un JSON invalide.')) from exc
    finally:
        os.close(descriptor)
    return folder, data


def saved_switch_identity(saved: Any) -> str:
    """Bind artifacts to immutable sequence/configuration, independent of location."""
    record = saved if isinstance(saved, dict) else getattr(saved, "snapshot", None)
    if not isinstance(record, dict):
        raise ValueError(tr("Le candidat enregistré ne contient pas de séquence exploitable."))
    selected = {key: record.get(key) for key in ("snapshot_id", "saved_at", "candidate", "request", "sequences")}
    return hashlib.sha256(json.dumps(selected, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def assert_switch_current(saved: SavedSwitch) -> None:
    """Guard a running job against replacement, using its already validated snapshot.

    Progress callbacks must not reimport validators while Streamlit evicts local
    modules on reload. Comparing the whole immutable record also catches edits
    that leave the snapshot ID unchanged.
    """
    _, current = _read_snapshot(saved.path)
    if current != saved.snapshot:
        raise ValueError("The saved candidate was replaced; this calculation belongs to its previous revision")


def load_switch(path: Path) -> SavedSwitch:
    """Read and validate a candidate, even after its design was moved."""
    folder, raw = _read_snapshot(path)
    data = _validate_snapshot(raw)
    return SavedSwitch(folder, data["name"], data["saved_at"], data)


def list_saved_switches(root: Path) -> list[SavedSwitch]:
    """List valid direct children; never traverse backup, staging or link trees."""
    root = _checked_path(root)
    if not root.exists():
        return []
    if not root.is_dir():
        raise ValueError(tr('La bibliothèque de switches enregistrés doit être un dossier.'))
    entries = []
    for folder in root.iterdir():
        if folder.name.startswith(".") or folder.is_symlink() or not folder.is_dir():
            continue
        try:
            entries.append(load_switch(folder))
        except (ValueError, OSError):
            continue
    return sorted(entries, key=lambda entry: (entry.saved_at, _fold(entry.name)), reverse=True)
