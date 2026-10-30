"""The oxDNA movie keeps its pinned offline molecular viewer dependency."""
import hashlib
import json
import unittest
from unittest.mock import patch

from aptaswitch_studio import molecular_viewer_assets as viewer


class MolecularViewerAssetTests(unittest.TestCase):
    def setUp(self):
        viewer.bundled_viewer_javascript.cache_clear()
        self.addCleanup(viewer.bundled_viewer_javascript.cache_clear)

    def test_offline_asset_is_pinned_with_license_and_manifest(self):
        manifest = json.loads((viewer.VENDOR_DIR / "manifest.json").read_text())
        self.assertEqual(manifest["version"], "2.5.5")
        self.assertEqual(manifest["license"], "BSD-3-Clause")
        for filename, info in manifest["files"].items():
            self.assertEqual(hashlib.sha256((viewer.VENDOR_DIR / filename).read_bytes()).hexdigest(), info["sha256"])
        self.assertIn("Redistribution", (viewer.VENDOR_DIR / "LICENSE.txt").read_text())
        self.assertIn("createViewer", viewer.bundled_viewer_javascript())

    def test_modified_asset_is_rejected_before_embedding(self):
        with patch.object(viewer.Path, "read_bytes", return_value=b"altered script"):
            with self.assertRaisesRegex(RuntimeError, "integrity"):
                viewer.bundled_viewer_javascript()
