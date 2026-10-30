"""Assemble and inspect the Windows launcher plus its complete Linux payload."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import tomllib
from zipfile import ZipFile

from check_bundled_app import check_distribution_archive

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_manifest(package: Path) -> None:
    runtime = package / 'runtime'
    backend = runtime / 'Apta2Switch-Studio-linux-x64.tar.gz'
    check_distribution_archive(backend)
    installer = runtime / 'install_backend.py'
    with (ROOT / 'pyproject.toml').open('rb') as handle:
        version = tomllib.load(handle)['project']['version']
    data = {'schema_version': 1, 'version': version,
            'backend': {'file': backend.name, 'sha256': sha256(backend)},
            'installer': {'file': installer.name, 'sha256': sha256(installer)}}
    (runtime / 'windows-package.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')


def check_archive(path: Path) -> None:
    root = 'Apta2Switch-Studio/'
    with ZipFile(path) as archive:
        names = archive.namelist()
        required = {root + p for p in ['Apta2Switch-Studio.exe', 'LICENSE', 'START-HERE.txt',
                    'runtime/windows-package.json', 'runtime/install_backend.py',
                    'runtime/Apta2Switch-Studio-linux-x64.tar.gz']}
        if not required.issubset(names) or not any(n.startswith(root + '_internal/') for n in names):
            raise ValueError('The Windows package is incomplete')
        if len(names) != len(set(names)):
            raise ValueError('Duplicate Windows archive members')
        for name in names:
            parts = name.split('/')
            if not name.startswith(root) or '..' in parts or '\\' in name or ':' in name:
                raise ValueError('Unsafe Windows archive path')
            if set(parts) & {'.git', '.aws', '.ssh', '.local_deps', 'designs', 'runs', 'oxdna', 'oxpy'}:
                raise ValueError('User data or external tools in Windows package')
            if any(p == '.env' or p.startswith('.env.') for p in parts):
                raise ValueError('Local environment in Windows package')
        manifest = json.loads(archive.read(root + 'runtime/windows-package.json'))
        if manifest.get('schema_version') != 1:
            raise ValueError('Invalid Windows package version')
        with tempfile.TemporaryDirectory() as folder:
            for key, filename in [('backend', 'Apta2Switch-Studio-linux-x64.tar.gz'), ('installer', 'install_backend.py')]:
                item = manifest[key]
                if item.get('file') != filename:
                    raise ValueError('Unexpected Windows runtime file')
                target = Path(folder) / filename
                with archive.open(root + 'runtime/' + filename) as source, target.open('wb') as output:
                    shutil.copyfileobj(source, output)
                if sha256(target) != item.get('sha256'):
                    raise ValueError('Windows payload hash mismatch')
                if key == 'backend':
                    check_distribution_archive(target)


def write_archive(package: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Build the ZIP in a sibling temporary path, then publish only after checks.
    with tempfile.TemporaryDirectory(dir=destination.parent) as folder:
        staging = Path(folder) / destination.name
        with ZipFile(staging, 'w') as archive:
            for path in sorted(package.rglob('*')):
                if path.is_symlink():
                    raise ValueError('Unexpected link in native Windows launcher package')
                if path.is_file():
                    archive.write(path, 'Apta2Switch-Studio/' + path.relative_to(package).as_posix())
        check_archive(staging)
        staging.replace(destination)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('manifest', 'archive', 'check'))
    parser.add_argument('path', type=Path)
    parser.add_argument('output', type=Path, nargs='?')
    args = parser.parse_args()
    if args.action == 'manifest':
        write_manifest(args.path)
    elif args.action == 'check':
        check_archive(args.path)
    elif args.output:
        write_archive(args.path, args.output)
    else:
        parser.error('archive requires an output ZIP')
