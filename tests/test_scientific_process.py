"""Cancellation, deadlines and environment isolation for scientific workers."""
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core import scientific_process as execution
from aptaswitch_core.scientific_process import CalculationCancelled, run_stage


class ScientificStageTests(unittest.TestCase):
    def setUp(self):
        self.work = Path(self.enterContext(tempfile.TemporaryDirectory()))


    def run_process(self, code, **kwargs):
        return run_stage([sys.executable, "-c", code], cwd=self.work, env=os.environ.copy(),
                         stage="fit", log_path=self.work / "stage.log", **kwargs)


    def assert_worker_stopped(self):
        pid = int((self.work / "pid").read_text())
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)


    def test_deadline_stops_native_process(self):
        with self.assertRaises(TimeoutError):
            self.run_process("import os,time;open('pid','w').write(str(os.getpid()));time.sleep(20)",
                             deadline=time.monotonic() + .3)
        self.assert_worker_stopped()


    def test_cancellation_stops_native_process(self):
        cancel = threading.Event()
        timer = threading.Timer(.3, cancel.set)
        timer.start()
        try:
            with self.assertRaises(CalculationCancelled):
                self.run_process("import os,time;open('pid','w').write(str(os.getpid()));time.sleep(20)",
                                 deadline=time.monotonic() + 10, cancel_check=cancel.is_set)
        finally:
            timer.cancel()
        self.assert_worker_stopped()


    def test_silent_native_work_has_progress_heartbeat(self):
        updates = []
        report = self.run_process("import time;time.sleep(.2)", deadline=time.monotonic() + 10,
                                  progress_callback=lambda *args: updates.append(args))
        self.assertEqual(report["return_code"], 0)
        self.assertTrue(updates)



class ScientificProcessTests(unittest.TestCase):
    def test_frozen_linux_restores_loader_path_without_changing_application(self):
        for original in [None, "", "/system/custom-libraries"]:
            environment = {"LD_LIBRARY_PATH": "/app/_internal", "PRESERVED": "yes"}
            if original is not None:
                environment["LD_LIBRARY_PATH_ORIG"] = original
            with self.subTest(original=original), \
                 patch.dict(os.environ, environment, clear=True), \
                 patch.object(execution.sys, "platform", "linux"), \
                 patch.object(execution.sys, "frozen", True, create=True):
                restored = execution.external_process_environment()
                self.assertEqual(dict(os.environ), environment)
            self.assertEqual(restored.get("LD_LIBRARY_PATH"), original)
            self.assertEqual(restored["PRESERVED"], "yes")
        with patch.dict(os.environ, {"LD_LIBRARY_PATH": "/custom/native-libs"}, clear=True), \
             patch.object(execution.sys, "platform", "linux"), \
             patch.object(execution.sys, "frozen", False, create=True):
            self.assertEqual(execution.external_process_environment()["LD_LIBRARY_PATH"], "/custom/native-libs")


    def test_process_groups_are_stopped_on_posix_and_windows(self):
        for system in ["Linux", "Darwin", "Windows"]:
            process = Mock(pid=123)
            process.poll.return_value = None
            process.wait.side_effect = [subprocess.TimeoutExpired("scientific-worker", 2), 0]
            with self.subTest(system=system), \
                 patch.object(execution.platform, "system", return_value=system), \
                 patch.object(execution.os, "killpg", create=True) as killpg, \
                 patch.object(execution.subprocess, "run") as run:
                execution.stop_process(process)
                options = execution.process_group_options()
            if system in {"Linux", "Darwin"}:
                self.assertEqual(killpg.call_args_list[0].args, (123, signal.SIGTERM))
                self.assertEqual(killpg.call_args_list[1].args, (123, signal.SIGKILL))
                self.assertTrue(options["start_new_session"])
            else:
                process.send_signal.assert_called_once()
                self.assertEqual(run.call_args.args[0], ["taskkill", "/PID", "123", "/T", "/F"])
                self.assertIn("creationflags", options)


    def test_stopped_process_is_left_alone(self):
        process = Mock(pid=123)
        process.poll.return_value = 0
        with patch.object(execution.os, "killpg", create=True) as killpg:
            execution.stop_process(process)
        killpg.assert_not_called()
        process.send_signal.assert_not_called()


    def test_windows_stop_falls_back_when_process_signalling_and_taskkill_fail(self):
        process = Mock(pid=123)
        process.poll.return_value = None
        process.send_signal.side_effect = OSError("no console")
        process.wait.side_effect = [subprocess.TimeoutExpired("scientific-worker", 2), 0]
        with patch.object(execution.platform, "system", return_value="Windows"), \
             patch.object(execution.subprocess, "run", side_effect=OSError("taskkill unavailable")):
            execution.stop_process(process)
        process.terminate.assert_called_once_with()
        process.kill.assert_called_once_with()


