#!/usr/bin/env python3
"""Install the bundled Linux application in its dedicated WSL account.

Run through ``wsl -d Ubuntu-22.04 -u root --exec python3``. This standalone
Python 3.10 helper installs no system packages and changes no WSL settings.
Only the final JSON record is written to stdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import posixpath
import pwd
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath


USER = "aptaswitch-studio"
HOME = Path("/home") / USER
APP = "Apta2Switch-Studio"
OWNER_MARKER = Path("/var/lib/aptaswitch-studio/windows-backend.json")
SCHEMA = 1


class InstallError(RuntimeError):
    """An installation precondition failed without deleting user data."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _archive_path(name: str) -> str:
    path = PurePosixPath(name)
    if (
        not name
        or "\\" in name
        or "\x00" in name
        or path.is_absolute()
        or ".." in path.parts
        or not path.parts
        or path.parts[0] != APP
    ):
        raise InstallError(f"Unsafe payload path: {name!r}")
    return str(path)


def validate_archive(archive: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    """Validate every member and link before any payload file is extracted."""
    members: dict[str, tarfile.TarInfo] = {}
    links: dict[str, str] = {}
    for member in archive.getmembers():
        name = _archive_path(member.name)
        if name in members:
            raise InstallError(f"Duplicate payload path: {name}")
        if not (member.isdir() or member.isreg() or member.issym() or member.islnk()):
            raise InstallError(f"Unsupported payload entry: {name}")
        members[name] = member
        if member.issym() or member.islnk():
            target = member.linkname
            if not target or target.startswith("/") or "\\" in target or "\x00" in target:
                raise InstallError(f"Unsafe payload link: {name}")
            if member.issym():
                # Reject even a temporary escape followed by a return into the
                # payload: intermediate external directories can be symlinks.
                parts = list(PurePosixPath(name).parent.parts)
                for part in PurePosixPath(target).parts:
                    if part == "..":
                        if len(parts) <= 1:
                            raise InstallError(f"Unsafe payload link: {name}")
                        parts.pop()
                    elif part != ".":
                        parts.append(part)
                target = posixpath.normpath(posixpath.join(posixpath.dirname(name), target))
            links[name] = _archive_path(target)

    root = members.get(APP)
    if root is not None and not root.isdir():
        raise InstallError("The payload root must be a directory.")
    for name, member in members.items():
        for parent in PurePosixPath(name).parents:
            parent_member = members.get(str(parent))
            if parent_member is not None and not parent_member.isdir():
                raise InstallError(f"A payload entry has a non-directory parent: {name}")
        if name in links:
            target = links[name]
            visited = {name}
            while target in links:
                if target in visited:
                    raise InstallError(f"Cyclic payload link: {name}")
                visited.add(target)
                target = links[target]
            target_member = members.get(target)
            if target_member is None:
                raise InstallError(f"Missing payload link target: {name}")
            if member.islnk() and not target_member.isreg():
                raise InstallError(f"A payload hard link must target a file: {name}")

    executable = members.get(f"{APP}/{APP}")
    if executable is None or not executable.isreg() or not executable.mode & 0o111:
        raise InstallError("The Linux payload does not contain the executable application.")
    return members


def _real_directory(path: Path, *, create: bool = False, mode: int = 0o755) -> None:
    """Reject symlinked installation paths, including their existing parents."""
    for part in [*reversed(path.parents), path]:
        try:
            info = part.lstat()
        except FileNotFoundError:
            if not create:
                raise InstallError(f"Required directory is missing: {part}")
            part.mkdir(mode=mode)
        else:
            if not stat.S_ISDIR(info.st_mode):
                raise InstallError(f"The installation path is not a real directory: {part}")


def _write_json(path: Path, value: dict) -> None:
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=".install-", delete=False) as output:
        temporary = Path(output.name)
        try:
            json.dump(value, output, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class AccountManager:
    """Provision only an account created and marked by this installer."""

    def __init__(self, home: Path = HOME, marker: Path = OWNER_MARKER, owner_uid: int = 0):
        self.home = home
        self.marker = marker
        self.owner_uid = owner_uid

    def _read_marker(self) -> dict | None:
        try:
            info = self.marker.lstat()
        except FileNotFoundError:
            return None
        _real_directory(self.marker.parent)
        parent = self.marker.parent.stat()
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_uid != self.owner_uid
            or info.st_mode & 0o022
            or parent.st_uid != self.owner_uid
            or parent.st_mode & 0o022
        ):
            raise InstallError("The backend ownership marker is not trusted.")
        try:
            marker = json.loads(self.marker.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise InstallError("The backend ownership marker is unreadable.") from error
        if not isinstance(marker, dict) or any(
            marker.get(key) != value
            for key, value in {"schema": SCHEMA, "user": USER, "home": str(self.home)}.items()
        ):
            raise InstallError("The backend ownership marker does not identify this application account.")
        return marker

    def _prepare_marker_directory(self) -> None:
        _real_directory(self.marker.parent, create=True)
        info = self.marker.parent.stat()
        if info.st_uid != self.owner_uid or info.st_mode & 0o022:
            raise InstallError("The backend ownership directory is not trusted.")

    def prepare(self):
        marker = self._read_marker()
        try:
            account = pwd.getpwnam(USER)
        except KeyError:
            account = None
        if account is not None:
            if marker is None:
                raise InstallError(f"The Linux account {USER!r} already exists and is not owned by this application.")
            if (
                account.pw_dir != str(self.home)
                or account.pw_uid < 1000
                or (marker.get("uid") is not None and marker["uid"] != account.pw_uid)
            ):
                raise InstallError("The existing backend account does not match its ownership marker.")
        else:
            if marker is not None and marker.get("uid") is not None:
                raise InstallError("The managed Linux account is missing; its stored data has been preserved.")
            if self.home.exists() or self.home.is_symlink():
                raise InstallError(f"The backend home already exists without its managed account: {self.home}")
            self._prepare_marker_directory()
            marker = {"schema": SCHEMA, "user": USER, "home": str(self.home), "uid": None}
            _write_json(self.marker, marker)
            try:
                subprocess.run(
                    ["useradd", "--create-home", "--user-group", "--shell", "/bin/bash",
                     "--home-dir", str(self.home), "--comment", "Apta2Switch Studio Windows backend", USER],
                    check=True, capture_output=True, text=True,
                )
            except (OSError, subprocess.CalledProcessError) as error:
                raise InstallError("Could not create the dedicated Linux account. The bundled Ubuntu installation must provide useradd.") from error
            account = pwd.getpwnam(USER)
        _real_directory(self.home)
        if self.home.stat().st_uid != account.pw_uid:
            raise InstallError("The backend home is not owned by its managed account.")
        _write_json(self.marker, {**marker, "uid": account.pw_uid})
        return account


def _extract(archive: tarfile.TarFile, members: dict[str, tarfile.TarInfo], destination: Path) -> None:
    for name, member in sorted(members.items(), key=lambda item: len(PurePosixPath(item[0]).parts)):
        path = destination / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if member.isdir():
            path.mkdir(exist_ok=True)
        elif member.isreg():
            source = archive.extractfile(member)
            if source is None:
                raise InstallError(f"Cannot read payload file: {name}")
            with source, path.open("xb") as output:
                shutil.copyfileobj(source, output, 1024 * 1024)
            path.chmod(member.mode & 0o777)
    # Symlinks are created last so no extraction can traverse an archive link.
    for name, member in members.items():
        path = destination / name
        if member.issym():
            path.symlink_to(member.linkname)
        elif member.islnk():
            target = _archive_path(member.linkname)
            while members[target].islnk() or members[target].issym():
                link = members[target]
                target = _archive_path(
                    posixpath.normpath(posixpath.join(posixpath.dirname(target), link.linkname))
                    if link.issym() else link.linkname
                )
            os.link(destination / target, path)


def _chown_tree(path: Path, uid: int, gid: int) -> None:
    for current, directories, files in os.walk(path, followlinks=False):
        os.chown(current, uid, gid, follow_symlinks=False)
        for name in directories + files:
            os.chown(Path(current) / name, uid, gid, follow_symlinks=False)


def install(payload: Path, expected_sha256: str, *, accounts: AccountManager | None = None) -> dict:
    if not re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256):
        raise InstallError("The payload SHA-256 must contain exactly 64 hexadecimal characters.")
    expected_sha256 = expected_sha256.lower()
    if sha256_file(payload) != expected_sha256:
        raise InstallError("The Linux payload checksum does not match. Download the complete Windows package again.")

    with tarfile.open(payload, "r:*") as archive:
        members = validate_archive(archive)
        accounts = accounts or AccountManager()
        account = accounts.prepare()
        home = accounts.home
        data = home / APP
        state = home / ".aptaswitch-studio"
        releases = home / ".local/share/aptaswitch/releases"
        for directory in (data, state, releases):
            _real_directory(directory, create=True)
            os.chown(directory, account.pw_uid, account.pw_gid)
        release = releases / expected_sha256
        executable = release / APP / APP
        marker_path = release / "payload.json"
        if release.exists() or release.is_symlink():
            _real_directory(release)
            _real_directory(executable.parent)
            try:
                if marker_path.is_symlink():
                    raise InstallError("The installed payload marker must not be a symbolic link.")
                marker = json.loads(marker_path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as error:
                raise InstallError(f"An incomplete installation is present at {release}; no saved data was changed.") from error
            if (
                marker != {"schema": SCHEMA, "sha256": expected_sha256}
                or not executable.is_file()
                or executable.is_symlink()
                or not executable.stat().st_mode & 0o111
            ):
                raise InstallError(f"The installed payload is incomplete at {release}; no saved data was changed.")
        else:
            staging = Path(tempfile.mkdtemp(prefix=".install-", dir=releases))
            try:
                _extract(archive, members, staging)
                _write_json(staging / "payload.json", {"schema": SCHEMA, "sha256": expected_sha256})
                _chown_tree(staging, account.pw_uid, account.pw_gid)
                staging.chmod(0o755)
                # A complete release appears atomically; older releases stay usable.
                staging.rename(release)
            finally:
                if staging.exists():
                    shutil.rmtree(staging)
    return {"executable": str(executable), "data_dir": str(data), "state_dir": str(state),
            "user": USER, "payload_sha256": expected_sha256}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    arguments = parser.parse_args()
    try:
        if sys.platform != "linux" or platform.machine().lower() not in {"x86_64", "amd64"}:
            raise InstallError("This Windows package requires an x86-64 Ubuntu WSL distribution.")
        if os.geteuid() != 0:
            raise InstallError("Run the backend installer through WSL with '-u root'.")
        result = install(arguments.payload, arguments.sha256)
    except (InstallError, OSError, tarfile.TarError) as error:
        print(f"Apta2Switch backend installation failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
