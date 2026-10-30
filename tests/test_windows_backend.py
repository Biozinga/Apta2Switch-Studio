"""The WSL helper is standalone; tests never create real system accounts."""

from __future__ import annotations

import importlib.util
import io
import json
import os
import tarfile
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch


if os.name == "nt":
    raise unittest.SkipTest("The Linux WSL backend installer requires POSIX ownership and links.")


_SCRIPT = Path(__file__).resolve().parents[1] / "packaging/windows/install_backend.py"
_SPEC = importlib.util.spec_from_file_location("windows_backend_installer", _SCRIPT)
backend = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(backend)


class WindowsBackendTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        # macOS /var is a symlink. The helper intentionally refuses such paths.
        self.root = Path(self.temporary.name).resolve()
        self.home = self.root / "home" / backend.USER
        self.home.mkdir(parents=True)
        self.account = SimpleNamespace(pw_uid=os.getuid(), pw_gid=os.getgid(), pw_dir=str(self.home))
        self.accounts = Mock(home=self.home)
        self.accounts.prepare.return_value = self.account
        self.chown = self.enterContext(patch.object(backend.os, "chown"))

    def make_payload(self, name="backend.tar.gz", *, content=b"#!/bin/sh\necho backend\n", extras=()):
        payload = self.root / name
        with tarfile.open(payload, "w:gz") as archive:
            root = tarfile.TarInfo(backend.APP)
            root.type = tarfile.DIRTYPE
            root.mode = 0o755
            archive.addfile(root)
            executable = tarfile.TarInfo(f"{backend.APP}/{backend.APP}")
            executable.mode = 0o755
            executable.size = len(content)
            archive.addfile(executable, io.BytesIO(content))
            for member, data in extras:
                archive.addfile(member, io.BytesIO(data) if data is not None else None)
        return payload, backend.sha256_file(payload)

    def file(self, name, data=b"library"):
        member = tarfile.TarInfo(name)
        member.size = len(data)
        member.mode = 0o644
        return member, data

    def link(self, name, target, *, hard=False):
        member = tarfile.TarInfo(name)
        member.type = tarfile.LNKTYPE if hard else tarfile.SYMTYPE
        member.linkname = target
        return member, None

    def test_install_update_and_reuse_preserve_saved_data_and_prior_release(self):
        payload, digest = self.make_payload()
        first = backend.install(payload, digest, accounts=self.accounts)
        self.assertEqual(set(first), {"executable", "data_dir", "state_dir", "user", "payload_sha256"})
        executable = Path(first["executable"])
        self.assertTrue(executable.stat().st_mode & 0o111)
        data = Path(first["data_dir"]) / "my-design.txt"
        data.write_text("experiment")
        settings = Path(first["state_dir"]) / "settings.json"
        settings.write_text('{"external_engine": "/mine"}')
        with patch.object(backend, "_extract", side_effect=AssertionError("Must reuse the existing release")):
            self.assertEqual(first, backend.install(payload, digest, accounts=self.accounts))
        newer_payload, newer_digest = self.make_payload("new.tar.gz", content=b"#!/bin/sh\necho v2\n")
        second = backend.install(newer_payload, newer_digest, accounts=self.accounts)
        self.assertNotEqual(first["executable"], second["executable"])
        self.assertTrue(executable.exists())
        self.assertEqual(data.read_text(), "experiment")
        self.assertEqual(settings.read_text(), '{"external_engine": "/mine"}')
        self.assertEqual(first["data_dir"], second["data_dir"])
        self.assertEqual(first["state_dir"], second["state_dir"])
        self.chown.assert_any_call(Path(first["data_dir"]), self.account.pw_uid, self.account.pw_gid)

    def test_checksum_failure_precedes_account_or_filesystem_changes(self):
        payload, _ = self.make_payload()
        for digest in ("0" * 64, "not-a-hash"):
            with self.subTest(digest=digest), self.assertRaises(backend.InstallError):
                backend.install(payload, digest, accounts=self.accounts)
        self.accounts.prepare.assert_not_called()
        self.assertEqual(list(self.home.iterdir()), [])

    def test_traversal_absolute_and_other_top_level_paths_are_rejected(self):
        for bad_name in ("../escape", "/tmp/escape", f"{backend.APP}/../escape", "another-app/file", f"{backend.APP}\\escape"):
            with self.subTest(path=bad_name):
                payload, digest = self.make_payload(extras=[self.file(bad_name)])
                with self.assertRaises(backend.InstallError):
                    backend.install(payload, digest, accounts=self.accounts)
        self.accounts.prepare.assert_not_called()
        self.assertEqual(list(self.home.iterdir()), [])

    def test_unsafe_links_and_entries_nested_under_links_are_rejected(self):
        cases = [
            [self.link(f"{backend.APP}/escape", "../../outside")],
            [self.link(f"{backend.APP}/escape", f"../{backend.APP}/{backend.APP}")],
            [self.link(f"{backend.APP}/escape", "/etc/passwd")],
            [self.link(f"{backend.APP}/escape", "missing")],
            [self.link(f"{backend.APP}/a", "b"), self.link(f"{backend.APP}/b", "a")],
            [self.link(f"{backend.APP}/escape", "../outside", hard=True)],
            [self.link(f"{backend.APP}/folder", backend.APP), self.file(f"{backend.APP}/folder/file")],
        ]
        for extras in cases:
            with self.subTest(extras=[member.name for member, _ in extras]):
                payload, digest = self.make_payload(extras=extras)
                with self.assertRaises(backend.InstallError):
                    backend.install(payload, digest, accounts=self.accounts)
        self.accounts.prepare.assert_not_called()

    def test_safe_relative_symlinks_and_hardlinks_are_preserved(self):
        extras = [
            self.file(f"{backend.APP}/_internal/libreal.so"),
            self.link(f"{backend.APP}/_internal/libalias.so", "libreal.so"),
            self.link(f"{backend.APP}/_internal/libcopy.so", f"{backend.APP}/_internal/libreal.so", hard=True),
        ]
        payload, digest = self.make_payload(extras=extras)
        result = backend.install(payload, digest, accounts=self.accounts)
        internal = Path(result["executable"]).parent / "_internal"
        self.assertTrue((internal / "libalias.so").is_symlink())
        self.assertEqual((internal / "libalias.so").read_bytes(), b"library")
        self.assertEqual((internal / "libcopy.so").stat().st_ino, (internal / "libreal.so").stat().st_ino)

    def test_devices_and_duplicate_members_are_rejected(self):
        device = tarfile.TarInfo(f"{backend.APP}/device")
        device.type = tarfile.CHRTYPE
        for extras in ([(device, None)], [self.file(f"{backend.APP}/{backend.APP}")]):
            with self.subTest(extras=extras):
                payload, digest = self.make_payload(extras=extras)
                with self.assertRaises(backend.InstallError):
                    backend.install(payload, digest, accounts=self.accounts)
        self.accounts.prepare.assert_not_called()

    def test_failed_extraction_never_publishes_partial_release(self):
        payload, digest = self.make_payload()
        with patch.object(backend, "_extract", side_effect=OSError("disk full")):
            with self.assertRaisesRegex(OSError, "disk full"):
                backend.install(payload, digest, accounts=self.accounts)
        releases = self.home / ".local/share/aptaswitch/releases"
        self.assertEqual(list(releases.iterdir()), [])

    def test_symlinked_data_or_release_directory_is_not_followed(self):
        payload, digest = self.make_payload()
        outside = self.root / "outside"
        outside.mkdir()
        (self.home / backend.APP).symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(backend.InstallError, "real directory"):
            backend.install(payload, digest, accounts=self.accounts)
        self.assertEqual(list(outside.iterdir()), [])

    def test_corrupt_existing_release_is_preserved_and_not_overwritten(self):
        payload, digest = self.make_payload()
        release = self.home / ".local/share/aptaswitch/releases" / digest
        release.mkdir(parents=True)
        diagnostic = release / "diagnostic.txt"
        diagnostic.write_text("partial previous installation")
        with self.assertRaises(backend.InstallError):
            backend.install(payload, digest, accounts=self.accounts)
        self.assertEqual(diagnostic.read_text(), "partial previous installation")


class WindowsBackendAccountTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.home = self.root / "home" / backend.USER
        self.marker = self.root / "managed" / "windows-backend.json"
        self.manager = backend.AccountManager(self.home, self.marker, owner_uid=os.getuid())
        self.account = SimpleNamespace(pw_uid=max(1000, os.getuid()), pw_gid=os.getgid(), pw_dir=str(self.home))

    def write_marker(self, uid=None):
        self.marker.parent.mkdir(parents=True, exist_ok=True)
        self.marker.write_text(json.dumps({"schema": 1, "user": backend.USER, "home": str(self.home), "uid": uid}))
        self.marker.chmod(0o600)

    def test_preexisting_unrelated_account_is_never_adopted(self):
        with patch.object(backend.pwd, "getpwnam", return_value=self.account), patch.object(backend.subprocess, "run") as run:
            with self.assertRaisesRegex(backend.InstallError, "not owned"):
                self.manager.prepare()
        run.assert_not_called()
        self.assertFalse(self.marker.exists())

    def test_mismatched_uid_or_home_is_rejected(self):
        for change in ({"uid": self.account.pw_uid + 1}, {"home": "/someone-else"}):
            self.write_marker(self.account.pw_uid)
            marker = json.loads(self.marker.read_text())
            marker.update(change)
            self.marker.write_text(json.dumps(marker))
            with self.subTest(change=change), patch.object(backend.pwd, "getpwnam", return_value=self.account):
                with self.assertRaises(backend.InstallError):
                    self.manager.prepare()

    def test_untrusted_or_symlinked_marker_is_rejected(self):
        self.write_marker(self.account.pw_uid)
        self.marker.chmod(0o666)
        with self.assertRaisesRegex(backend.InstallError, "not trusted"):
            self.manager.prepare()
        self.marker.unlink()
        elsewhere = self.root / "elsewhere.json"
        elsewhere.write_text("{}")
        self.marker.symlink_to(elsewhere)
        with self.assertRaisesRegex(backend.InstallError, "not trusted"):
            self.manager.prepare()

    def test_new_account_is_marked_before_useradd_and_finalized_afterwards(self):
        def create_account(*args, **kwargs):
            self.assertEqual(json.loads(self.marker.read_text())["uid"], None)
            self.home.mkdir(parents=True)

        real_stat = Path.stat

        def account_owned_home(path, **kwargs):
            value = real_stat(path, **kwargs)
            if path == self.home:
                return SimpleNamespace(st_uid=self.account.pw_uid, st_mode=value.st_mode)
            return value

        with patch.object(backend.pwd, "getpwnam", side_effect=[KeyError(backend.USER), self.account]), \
             patch.object(backend.subprocess, "run", side_effect=create_account) as run, \
             patch.object(Path, "stat", account_owned_home):
            self.assertEqual(self.manager.prepare(), self.account)
        command = run.call_args.args[0]
        self.assertEqual(command[0], "useradd")
        self.assertEqual(command[-1], backend.USER)
        self.assertEqual(json.loads(self.marker.read_text())["uid"], self.account.pw_uid)

    def test_preexisting_home_without_account_is_preserved(self):
        self.home.mkdir(parents=True)
        design = self.home / "my-design"
        design.write_text("keep")
        with patch.object(backend.pwd, "getpwnam", side_effect=KeyError(backend.USER)), patch.object(backend.subprocess, "run") as run:
            with self.assertRaisesRegex(backend.InstallError, "home already exists"):
                self.manager.prepare()
        self.assertEqual(design.read_text(), "keep")
        run.assert_not_called()

    def test_managed_account_is_reused_without_useradd(self):
        self.write_marker(self.account.pw_uid)
        self.home.mkdir(parents=True)
        real_stat = Path.stat

        def account_owned_home(path, **kwargs):
            value = real_stat(path, **kwargs)
            if path == self.home:
                return SimpleNamespace(st_uid=self.account.pw_uid, st_mode=value.st_mode)
            return value

        with patch.object(backend.pwd, "getpwnam", return_value=self.account), \
             patch.object(backend.subprocess, "run") as run, \
             patch.object(Path, "stat", account_owned_home):
            self.assertEqual(self.manager.prepare(), self.account)
        run.assert_not_called()

    def test_missing_managed_account_does_not_reassign_its_existing_data(self):
        self.write_marker(self.account.pw_uid)
        with patch.object(backend.pwd, "getpwnam", side_effect=KeyError(backend.USER)), \
             patch.object(backend.subprocess, "run") as run:
            with self.assertRaisesRegex(backend.InstallError, "managed Linux account is missing"):
                self.manager.prepare()
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
