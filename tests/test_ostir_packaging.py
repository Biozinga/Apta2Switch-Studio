"""Source distribution excludes user data while retaining actual edited code."""
import importlib.util
from pathlib import Path
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ostir_sources", ROOT / "scripts/prepare_ostir_sources.py")
sources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sources)


class OstirPackagingTests(unittest.TestCase):
    def test_source_snapshot_keeps_dirty_build_inputs_but_not_runs_or_caches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in {
                "LICENSE": "original licence", "pyproject.toml": "actual dependency pins",
                "src/app.py": "changed, uncommitted code", "scripts/build.sh": "build instructions",
                "src/__pycache__/app.pyc": "cached", "src/app.egg-info/PKG-INFO": "generated",
                "src/.env": "secret", "designs/private.json": "user sequence",
                ".streamlit/config.toml": "public application configuration",
                ".streamlit/secrets.toml": "private service key",
                "src/.aws/credentials": "private cloud key",
                "src/private/secrets.toml": "private service key",
                "src/validation-results/private.json": "private calculation",
                "runs/ffs/checkpoint.dat": "user checkpoint", ".local_deps/settings.json": "local paths",
            }.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
            target = root / "source.tar.gz"
            sources.write_application_source(root, target)
            with tarfile.open(target) as package:
                names = set(package.getnames())
                self.assertEqual(names, {
                    "Apta2Switch-Studio/LICENSE", "Apta2Switch-Studio/pyproject.toml",
                    "Apta2Switch-Studio/src/app.py", "Apta2Switch-Studio/scripts/build.sh",
                    "Apta2Switch-Studio/.streamlit/config.toml",
                })
                self.assertEqual(package.extractfile("Apta2Switch-Studio/src/app.py").read(),
                                 b"changed, uncommitted code")

    def test_pyinstaller_collects_parameters_native_bindings_and_metadata(self):
        spec = (ROOT / "packaging/aptaswitch-studio.spec").read_text()
        self.assertIn('collect_all("ostir")', spec)
        self.assertIn('collect_all("RNA")', spec)
        self.assertIn('copy_metadata("OSTIR")', spec)
        self.assertIn('copy_metadata("ViennaRNA")', spec)
        self.assertIn('prepare_ostir_sources.py', spec)


if __name__ == "__main__":
    unittest.main()
