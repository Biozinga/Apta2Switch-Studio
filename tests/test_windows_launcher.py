from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('windows_launcher', ROOT / 'packaging/windows/launcher.py')
launcher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(launcher)


class WindowsLauncherTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        runtime = self.root / 'runtime'
        runtime.mkdir()
        self.manifest = {'schema_version': 1, 'version': '0.2.1b2'}
        for key, name in [('backend', 'Apta2Switch-Studio-linux-x64.tar.gz'), ('installer', 'install_backend.py')]:
            file = runtime / name
            file.write_bytes(b'fixture ' + key.encode())
            self.manifest[key] = {'file': name, 'sha256': launcher.digest(file)}
        self.save_manifest()
        self.backend = launcher.WindowsBackend(self.root)
        self.install = {'executable': '/home/aptaswitch-studio/.local/share/aptaswitch/releases/' + 'a'*64 + '/Apta2Switch-Studio/Apta2Switch-Studio',
                        'data_dir': '/home/aptaswitch-studio/Apta2Switch-Studio',
                        'state_dir': '/home/aptaswitch-studio/.aptaswitch-studio',
                        'user': launcher.USER, 'payload_sha256': self.manifest['backend']['sha256']}

    def save_manifest(self):
        (self.root / 'runtime/windows-package.json').write_text(json.dumps(self.manifest))

    def test_entire_payload_is_checked_before_use(self):
        self.assertEqual(launcher.read_package(self.root), self.manifest)
        (self.root / 'runtime/install_backend.py').write_bytes(b'changed')
        with self.assertRaises(ValueError):
            launcher.read_package(self.root)

    def test_manifest_cannot_redirect_execution_to_another_script(self):
        self.manifest['installer']['file'] = '../unrelated.py'
        self.save_manifest()
        with self.assertRaises(ValueError):
            launcher.read_package(self.root)

    def test_wsl_lists_decode_utf16_and_ignore_localized_column_headers(self):
        text = '  NOM           ÉTAT         VERSION\r\n* Ubuntu-22.04  En cours        2\r\n'
        for encoding in ('utf-16', 'utf-16-le', 'utf-8'):
            self.assertEqual(launcher.distro_version(launcher.decode_output(text.encode(encoding))), 2)
        self.assertIsNone(launcher.distro_version('Ubuntu-22.04-old Stopped 2'))
        self.assertEqual(launcher.distro_version('  Ubuntu-22.04 Stopped 1'), 1)

    def test_missing_or_wsl1_distribution_requires_setup_without_installing(self):
        for output in ('', 'Ubuntu-22.04 Stopped 1'):
            with self.subTest(output=output), patch.object(launcher, 'run_command', return_value=output) as run:
                with self.assertRaises(launcher.SetupRequired):
                    self.backend.check_wsl()
                self.assertEqual(run.call_count, 1)
                self.assertNotIn('--install', run.call_args.args[0])

    def test_setup_elevates_only_machine_activation_not_the_users_distribution(self):
        outputs = ['', RuntimeError('WSL disabled'), '', 'ready', '', 'Ubuntu-22.04 Running 2', 'Ubuntu-22.04 Running 2']
        with patch.object(launcher, 'run_command', side_effect=outputs) as run:
            self.backend.setup()
        commands = [call.args[0] for call in run.call_args_list]
        elevated = [cmd for cmd in commands if cmd[0] == 'powershell.exe']
        self.assertEqual(len(elevated), 1)
        self.assertIn('-Verb RunAs', elevated[0][-1])
        self.assertIn("'--no-distribution'", elevated[0][-1])
        self.assertNotIn(launcher.DISTRO, elevated[0][-1])
        self.assertIn([launcher.wsl_executable(), '--install', '--distribution', launcher.DISTRO, '--no-launch'], commands)

    def test_setup_preserves_ready_environment_and_resumes_after_restart(self):
        with patch.object(launcher, 'run_command', return_value='Ubuntu-22.04 Running 2') as run:
            self.backend.setup()
            self.assertEqual(run.call_count, 1)
        with patch.object(launcher, 'run_command', side_effect=['', RuntimeError('disabled'), '', RuntimeError('reboot')]) as run:
            self.assertIn('Redémarrez', self.backend.setup())
            self.assertFalse(any('--distribution' in call.args[0] for call in run.call_args_list))

    def test_conversion_remains_in_current_user_context(self):
        with patch.object(launcher, 'run_command', side_effect=['Ubuntu-22.04 Stopped 1', 'ready',
                    'Ubuntu-22.04 Stopped 1', '', 'Ubuntu-22.04 Stopped 2']) as run:
            self.backend.setup()
        commands = [call.args[0] for call in run.call_args_list]
        self.assertIn([launcher.wsl_executable(), '--set-version', launcher.DISTRO, '2'], commands)
        self.assertFalse(any(cmd[0] == 'powershell.exe' for cmd in commands))

    def test_install_uses_separate_arguments_and_validates_returned_paths(self):
        with patch.object(self.backend, 'mounted_path', side_effect=['/mnt/c/Éxample folder/install.py', '/mnt/c/a $(x).tar.gz']), \
             patch.object(self.backend, 'linux', return_value=json.dumps(self.install)) as call:
            self.assertEqual(self.backend.install(self.manifest), self.install)
        self.assertEqual(call.call_args.args[0][:3], ['python3', '/mnt/c/Éxample folder/install.py', '--payload'])
        self.assertEqual(call.call_args.kwargs['user'], 'root')
        bad = {**self.install, 'executable': '/home/aptaswitch-studio/../unrelated'}
        with patch.object(self.backend, 'mounted_path', return_value='/mnt/c/path'), \
             patch.object(self.backend, 'linux', return_value=json.dumps(bad)), self.assertRaises(ValueError):
            self.backend.install(self.manifest)

    def test_existing_server_is_reused_without_starting_or_stopping_workers(self):
        with patch.object(self.backend, 'check_wsl'), patch.object(self.backend, 'install', return_value=self.install), \
             patch.object(self.backend, 'running_url', return_value='http://127.0.0.1:8765'), patch.object(self.backend, 'linux') as linux:
            self.assertEqual(self.backend.launch(), 'http://127.0.0.1:8765')
        linux.assert_not_called()

    def test_background_start_preserves_data_and_is_detached_from_launcher(self):
        with patch.object(self.backend, 'check_wsl'), patch.object(self.backend, 'install', return_value=self.install), \
             patch.object(self.backend, 'running_url', side_effect=[None, 'http://127.0.0.1:8765']), patch.object(self.backend, 'linux') as linux:
            self.backend.launch()
        args = linux.call_args.args[0]
        self.assertIn('nohup "$@"', args[2])
        self.assertIn('APTASWITCH_DATA_DIR=' + self.install['data_dir'], args)
        self.assertIn('APTASWITCH_STATE_DIR=' + self.install['state_dir'], args)
        self.assertEqual(args[-1], self.install['executable'])
        self.assertFalse(any(word in args for word in ('--terminate', '--shutdown', '--unregister')))

    def test_server_identity_and_windows_localhost_are_both_checked(self):
        state = {'app_id': 'aptaswitch-studio', 'installation': self.install['executable'], 'pid': 42, 'port': 8765}
        response = Mock(status=200)
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)
        for installation, valid in [(state['installation'], True), (state['installation'].replace('a'*64, 'b'*64), True), ('/other/application', False)]:
            response.read.side_effect = [b'ok', b'<div id="root"></div>']
            with self.subTest(installation=installation), patch.object(self.backend, 'linux', return_value=json.dumps({**state, 'installation': installation})), \
                 patch.object(launcher.urllib.request, 'urlopen', return_value=response) as http:
                self.assertEqual(self.backend.running_url(self.install), 'http://127.0.0.1:8765' if valid else None)
                self.assertEqual(http.call_count, 2 if valid else 0)


if __name__ == '__main__':
    unittest.main()
