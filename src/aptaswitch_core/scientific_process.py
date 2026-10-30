"""Cancellable subprocess execution for external scientific tools."""
from __future__ import annotations

from contextlib import ExitStack
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time


def external_process_environment() -> dict[str, str]:
    """Keep PyInstaller's Linux libraries out of independently supplied tools.

    Restore only the child's loader path. The application and its bundled
    numerical libraries retain their original process environment.
    """
    environment = os.environ.copy()
    if sys.platform.startswith("linux") and getattr(sys, "frozen", False):
        original = environment.get("LD_LIBRARY_PATH_ORIG")
        if original is None:
            environment.pop("LD_LIBRARY_PATH", None)
        else:
            environment["LD_LIBRARY_PATH"] = original
    return environment


def process_group_options() -> dict:
    """Keep worker subprocesses independently cancellable on every host."""
    if platform.system() == "Windows":
        return {"creationflags": getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x200)}
    return {"start_new_session": True}


def stop_process(process: subprocess.Popen) -> None:
    """Stop the worker and its child processes."""
    if process.poll() is not None:
        return
    if platform.system() == "Windows":
        try:
            process.send_signal(getattr(signal, "CTRL_BREAK_EVENT", 1))
        except (OSError, ValueError):
            process.terminate()
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    try:
        process.wait(timeout=2)
        return
    except subprocess.TimeoutExpired:
        pass
    if platform.system() == "Windows":
        try:
            subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           timeout=5, check=False, env=external_process_environment())
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
    else:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()


class CalculationCancelled(RuntimeError):
    """A local scientific stage was cancelled before publishing its output."""


def run_stage(command: list[str], *, cwd: Path, env: dict, stage: str,
              log_path: Path, deadline: float, cancel_check=None,
              progress_callback=None, stdout_path: Path | None = None,
              stdin_path: Path | None = None) -> dict:
    """Run one stage; caller records commands before entry, even on failure.

    Native stages have a wall deadline, cancellation and progress heartbeat.
    """
    started = time.monotonic()
    if cancel_check and cancel_check():
        raise CalculationCancelled("Scientific calculation cancelled")
    if started >= deadline:
        raise TimeoutError(f"Scientific stage {stage} exceeded its time limit")
    process = None
    with ExitStack() as files:
        log = files.enter_context(log_path.open("ab"))
        output = files.enter_context(stdout_path.open("wb")) if stdout_path else log
        source = files.enter_context(stdin_path.open("rb")) if stdin_path else subprocess.DEVNULL
        try:
            process = subprocess.Popen(command, cwd=cwd, env=env, stdin=source,
                                       stdout=output, stderr=log, **process_group_options())
            next_update = started
            while process.poll() is None:
                if cancel_check and cancel_check():
                    raise CalculationCancelled("Scientific calculation cancelled")
                now = time.monotonic()
                if now >= deadline:
                    raise TimeoutError(f"Scientific stage {stage} exceeded its time limit")
                if now >= next_update:
                    if progress_callback:
                        progress_callback(stage, now - started, "")
                    next_update = now + 5.0
                time.sleep(0.1)
            if process.returncode:
                raise RuntimeError(f"Scientific stage {stage} failed (code {process.returncode}); see {log_path.name}")
            return {"return_code": process.returncode, "elapsed_seconds": time.monotonic() - started}
        finally:
            if process is not None and process.poll() is None:
                stop_process(process)
