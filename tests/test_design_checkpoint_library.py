"""Pending designs remain discoverable without stale or duplicate results."""
import json
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio.run_library import discover_runs, load_run
from tests.design_fixtures import design_result_fixture


class DesignCheckpointLibraryTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name).resolve()
        self.result = design_result_fixture(RunRequest(
            SequenceInput('GGTT', 'AACC'), SwitchArchitecture(2, 2),
            trials=3, output_dir=self.root))
        self.result.run_id = 'checkpoint_trial_run'
        self.result.output_dir = self.root / self.result.run_id
        self.result.output_dir.mkdir()
        self.checkpoint = self.result.output_dir / 'design_checkpoint.json'
        self.export = self.result.output_dir / f'{self.result.run_id}_results.json'

    def write(self, path, *, status, count):
        payload = json.loads(json.dumps(self.result.to_dict()))
        payload['status'] = status
        payload['results'] = [row for row in payload['results'] if row['trial'] <= count]
        payload['optimization']['trials'] = payload['optimization']['trials'][:count]
        path.write_text(json.dumps(payload))

    def test_checkpoint_opens_before_any_trial_or_exports_exist(self):
        self.write(self.checkpoint, status='paused', count=0)
        loaded = load_run(self.result.output_dir)
        self.assertEqual(loaded.source, self.checkpoint)
        self.assertEqual(loaded.result.status, 'paused')
        self.assertEqual(loaded.result.candidates, [])
        self.assertEqual([r.path for r in discover_runs([self.root])], [self.checkpoint])

    def test_pending_checkpoint_supersedes_old_pause_exports_without_rewriting_them(self):
        self.write(self.export, status='paused', count=1)
        old_bytes = self.export.read_bytes()
        self.write(self.checkpoint, status='running', count=2)
        loaded = load_run(self.result.output_dir)
        self.assertEqual(len(loaded.result.optimization['trials']), 2)
        rows = discover_runs([self.root, self.result.output_dir, self.export])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].path, self.checkpoint)
        self.assertEqual(rows[0].candidate_count, 2)
        self.assertEqual(self.export.read_bytes(), old_bytes)

    def test_final_archive_wins_after_clean_finish(self):
        for status, count in [('stopped', 2), ('completed', 3)]:
            with self.subTest(status=status):
                self.write(self.export, status=status, count=count)
                self.write(self.checkpoint, status=status, count=count)
                self.assertEqual(load_run(self.result.output_dir).source, self.export)
                self.assertEqual([r.path for r in discover_runs([self.root])], [self.export])

    def test_broken_checkpoint_does_not_hide_intact_final_archive(self):
        self.write(self.export, status='completed', count=3)
        self.checkpoint.write_text('{')
        self.assertEqual(load_run(self.result.output_dir).source, self.export)
        self.assertEqual([r.path for r in discover_runs([self.root])], [self.export])

    def test_checkpoint_symlink_is_not_discovered(self):
        source = self.root / 'external.json'
        self.write(source, status='paused', count=1)
        self.checkpoint.symlink_to(source)
        self.assertEqual(discover_runs([self.result.output_dir]), [])
        with self.assertRaises(ValueError):
            load_run(self.result.output_dir)
