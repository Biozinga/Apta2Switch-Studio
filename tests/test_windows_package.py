from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from zipfile import ZipFile

from tests import test_distribution_archive as distribution_fixtures

ROOT = Path(__file__).resolve().parents[1]
# Maintainer scripts are intentionally independent from the application package.
sys.path.insert(0, str(ROOT / 'scripts'))
try:
    import windows_package as package
finally:
    sys.path.pop(0)


class WindowsPackageTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        entries, _ = distribution_fixtures.DistributionArchiveTests.entries(None, 'tar')
        backend = distribution_fixtures.DistributionArchiveTests.write(SimpleNamespace(work=self.root), 'tar', entries)
        self.app = self.root / 'Apta2Switch-Studio'
        runtime = self.app / 'runtime'
        runtime.mkdir(parents=True)
        shutil.copyfile(backend, runtime / 'Apta2Switch-Studio-linux-x64.tar.gz')
        (runtime / 'install_backend.py').write_text('print("fixture")')
        for name in ('LICENSE', 'START-HERE.txt', 'Apta2Switch-Studio.exe', '_internal/python.dll'):
            target = self.app / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b'fixture')
        package.write_manifest(self.app)

    def test_archive_includes_native_runtime_and_verifies_embedded_linux_payload(self):
        archive = self.root / 'windows.zip'
        package.write_archive(self.app, archive)
        package.check_archive(archive)
        with ZipFile(archive) as content:
            self.assertIn('Apta2Switch-Studio/_internal/python.dll', content.namelist())
            data = json.loads(content.read('Apta2Switch-Studio/runtime/windows-package.json'))
            self.assertEqual(data['backend']['sha256'], package.sha256(self.app / 'runtime/Apta2Switch-Studio-linux-x64.tar.gz'))

    def test_hash_mismatch_does_not_replace_an_existing_distribution(self):
        archive = self.root / 'windows.zip'
        archive.write_bytes(b'previous distribution')
        (self.app / 'runtime/install_backend.py').write_text('different content')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            package.write_archive(self.app, archive)
        self.assertEqual(archive.read_bytes(), b'previous distribution')

    def test_incomplete_linux_runtime_cannot_be_wrapped_in_valid_windows_package(self):
        entries, runtime = distribution_fixtures.DistributionArchiveTests.entries(None, 'tar')
        entries.pop(runtime + 'aptaswitch_core/oxdna_initialization.py')
        backend = distribution_fixtures.DistributionArchiveTests.write(SimpleNamespace(work=self.root), 'tar', entries)
        shutil.copyfile(backend, self.app / 'runtime/Apta2Switch-Studio-linux-x64.tar.gz')
        with self.assertRaisesRegex(RuntimeError, 'oxDNA workers'):
            package.write_manifest(self.app)

    def test_user_data_or_secrets_are_not_packaged(self):
        for name in ('designs/private.json', '.env', '_internal/.aws/credentials'):
            file = self.app / name
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_bytes(b'private')
            with self.subTest(name=name), self.assertRaises(ValueError):
                package.write_archive(self.app, self.root / 'windows.zip')
            file.unlink()


if __name__ == '__main__':
    unittest.main()
