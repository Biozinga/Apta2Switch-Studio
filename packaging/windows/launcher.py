"""Native Windows entry point for the full Linux scientific application in WSL2.

This small stdlib-only launcher is frozen separately on Windows. It never
imports the scientific application into Windows Python or stops running jobs.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import queue
import re
import subprocess
import sys
import threading
import time
import urllib.request
import webbrowser

DISTRO = "Ubuntu-22.04"
USER = "aptaswitch-studio"
DATA_HOME = f"/home/{USER}"


class SetupRequired(RuntimeError):
    pass


def package_root() -> Path:
    return Path(sys.executable).parent if getattr(sys, "frozen", False) else Path(__file__).parent


def windows_state_dir() -> Path:
    return Path(os.environ.get("LOCALAPPDATA", Path.home())) / "Apta2Switch-Studio" / "windows"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def read_package(root: Path) -> dict:
    runtime = root / "runtime"
    data = json.loads((runtime / "windows-package.json").read_text(encoding="utf-8"))
    if data.get("schema_version") != 1 or not isinstance(data.get("version"), str):
        raise ValueError("Le paquet Windows est incomplet ou incompatible.")
    for key, filename in (("backend", "Apta2Switch-Studio-linux-x64.tar.gz"),
                          ("installer", "install_backend.py")):
        item = data.get(key, {})
        if (item.get("file") != filename
                or not re.fullmatch(r"[0-9a-f]{64}", str(item.get("sha256", "")))):
            raise ValueError("Le manifeste Windows est invalide.")
        path = runtime / filename
        if path.is_symlink() or not path.is_file() or digest(path) != item["sha256"]:
            raise ValueError(f"Le fichier {filename} est incomplet. Extrayez à nouveau toute l’archive.")
    return data


def decode_output(value: bytes) -> str:
    # wsl --list writes UTF-16LE even when redirected; Linux processes use UTF-8.
    if value.startswith(b"\xff\xfe") or b"\0" in value[:100]:
        return value.decode("utf-16", errors="replace").lstrip("\ufeff").strip()
    return value.decode("utf-8", errors="replace").strip()


def wsl_executable() -> str:
    return str(Path(os.environ.get("WINDIR", r"C:\Windows")) / "System32" / "wsl.exe")


def run_command(arguments: list[str], *, timeout: float = 120, check: bool = True) -> str:
    result = subprocess.run(arguments, stdin=subprocess.DEVNULL, capture_output=True,
                            timeout=timeout, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    output = decode_output(result.stdout)
    if check and result.returncode:
        error = decode_output(result.stderr) or output
        raise RuntimeError(error or f"La commande a échoué (code {result.returncode}).")
    return output if result.returncode == 0 else ""


def distro_version(output: str) -> int | None:
    for line in output.splitlines():
        match = re.match(r"^\s*\*?\s*" + re.escape(DISTRO) + r"\s+.+?\s+([12])\s*$", line)
        if match:
            return int(match.group(1))
    return None


class WindowsBackend:
    def __init__(self, root: Path, notify=lambda text: None):
        self.root = root
        self.notify = notify

    def linux(self, arguments: list[str], *, user: str = USER, **kwargs) -> str:
        return run_command([wsl_executable(), "--distribution", DISTRO, "--user", user,
                            "--exec", *arguments], **kwargs)

    def check_wsl(self) -> None:
        try:
            output = run_command([wsl_executable(), "--list", "--verbose"], timeout=20)
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
            raise SetupRequired("Windows doit d’abord activer WSL2. Cliquez sur « Préparer Windows ».") from exc
        version = distro_version(output)
        if version != 2:
            if version == 1:
                raise SetupRequired("Ubuntu 22.04 utilise WSL1. Le bouton « Préparer Windows » propose sa conversion vers WSL2.")
            raise SetupRequired("Ubuntu 22.04 manque. Cliquez sur « Préparer Windows » ; une connexion Internet sera nécessaire.")

    def setup(self) -> str:
        """Called only after an explicit setup confirmation in the Windows UI."""
        output = run_command([wsl_executable(), "--list", "--verbose"], check=False, timeout=20)
        version = distro_version(output)
        if version == 2:
            return "Windows est prêt. Cliquez sur « Ouvrir l’application »."
        # Only fixed application arguments enter PowerShell; no path or user text
        # is interpolated. UAC obtains the user's Windows administrator consent.
        try:
            run_command([wsl_executable(), "--status"], timeout=20)
        except (OSError, RuntimeError, subprocess.TimeoutExpired):
            # Distributions belong to a Windows user. Only enable the machine
            # component under UAC; never create Ubuntu in another admin's account.
            command = ("$ErrorActionPreference='Stop'; "
                       "$p=Start-Process -FilePath \"$env:WINDIR\\System32\\wsl.exe\" "
                       "-ArgumentList @('--install','--no-distribution') -Verb RunAs -Wait -PassThru; "
                       "if ($p.ExitCode -eq 3010) { exit 0 }; exit $p.ExitCode")
            run_command(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", command], timeout=3600)
            try:
                run_command([wsl_executable(), "--status"], timeout=20)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                return "WSL a été activé. Redémarrez Windows, puis cliquez à nouveau sur « Préparer Windows »."
        if version is None:
            run_command([wsl_executable(), "--install", "--distribution", DISTRO, "--no-launch"], timeout=3600)
        # Both installation and conversion run with the original Windows user.
        output = run_command([wsl_executable(), "--list", "--verbose"], timeout=20)
        if distro_version(output) == 1:
            run_command([wsl_executable(), "--set-version", DISTRO, "2"], timeout=3600)
        self.check_wsl()
        return "Windows est prêt. Cliquez sur « Ouvrir l’application »."

    def mounted_path(self, path: Path) -> str:
        return self.linux(["wslpath", "-a", "-u", str(path.resolve())], user="root")

    def install(self, package: dict) -> dict:
        self.notify("Préparation de l’application et des dossiers de travail…")
        installer = self.mounted_path(self.root / "runtime" / package["installer"]["file"])
        archive = self.mounted_path(self.root / "runtime" / package["backend"]["file"])
        try:
            raw = self.linux(["python3", installer, "--payload", archive,
                              "--sha256", package["backend"]["sha256"]], user="root", timeout=900)
        except RuntimeError as exc:
            raise RuntimeError("La préparation Ubuntu a échoué. Python 3 doit être disponible dans Ubuntu 22.04.\n" + str(exc)) from exc
        data = json.loads(raw)
        if data.get("user") != USER or data.get("payload_sha256") != package["backend"]["sha256"]:
            raise ValueError("L’installation Ubuntu ne correspond pas au paquet demandé.")
        for key in ("executable", "data_dir", "state_dir"):
            path = PurePosixPath(data.get(key, ""))
            if not path.is_relative_to(DATA_HOME) or ".." in path.parts:
                raise ValueError("L’installation Ubuntu a renvoyé un emplacement invalide.")
        return data

    def running_url(self, install: dict) -> str | None:
        try:
            state = json.loads(self.linux(["cat", install["state_dir"] + "/web_server.json"],
                                         check=False, timeout=10))
            # Reuse a previous managed version that still has calculations in
            # progress. Installing an update must not replace its live server.
            managed = re.fullmatch(re.escape(DATA_HOME) +
                r"/\.local/share/aptaswitch/releases/[0-9a-f]{64}/Apta2Switch-Studio/Apta2Switch-Studio",
                str(state.get("installation", "")))
            if (state.get("app_id") != "aptaswitch-studio"
                    or not managed
                    or type(state.get("port")) is not int or not 1 <= state["port"] <= 65535
                    or type(state.get("pid")) is not int or state["pid"] <= 0):
                return None
            self.linux(["/bin/kill", "-0", str(state["pid"])], timeout=10)
            url = f"http://127.0.0.1:{state['port']}"
            with urllib.request.urlopen(url + "/_stcore/health", timeout=2) as response:
                if response.status != 200 or response.read(16).strip() != b"ok":
                    return None
            with urllib.request.urlopen(url + "/", timeout=2) as response:
                if response.status != 200 or b'id="root"' not in response.read(65536):
                    return None
            return url
        except (OSError, ValueError, TypeError, RuntimeError, subprocess.TimeoutExpired):
            return None

    def launch(self) -> str:
        self.check_wsl()
        self.notify("Vérification du paquet Windows…")
        package = read_package(self.root)
        install = self.install(package)
        url = self.running_url(install)
        if url:
            return url
        # nohup decouples the server and its scientific workers from this GUI.
        # Closing either the launcher or the browser never kills a calculation.
        self.notify("Démarrage de l’application…")
        self.linux(["/bin/sh", "-c", 'log=$1; shift; nohup "$@" </dev/null >>"$log" 2>&1 &',
                    "aptaswitch-start", install["state_dir"] + "/server.log", "env",
                    "APTASWITCH_NO_BROWSER=1", "APTASWITCH_DATA_DIR=" + install["data_dir"],
                    "APTASWITCH_STATE_DIR=" + install["state_dir"], install["executable"]])
        deadline = time.monotonic() + 120
        while time.monotonic() < deadline:
            url = self.running_url(install)
            if url:
                return url
            time.sleep(1)
        raise RuntimeError("Le serveur n’a pas encore répondu. Les calculs ne sont pas interrompus. "
                           "Consultez le journal dans Ubuntu : " + install["state_dir"] + "/server.log")


def main() -> int:
    if sys.argv[1:2] == ["--check-package"]:
        # Release CI runs this in the actual frozen Windows executable.
        destination = Path(sys.argv[2])
        try:
            data = read_package(package_root())
            import tkinter as tk
            root = tk.Tk()
            root.withdraw()
            root.update()
            root.destroy()
            destination.write_text(json.dumps({"ok": True, "version": data["version"]}), encoding="utf-8")
            return 0
        except Exception as exc:
            destination.write_text(json.dumps({"ok": False, "error": str(exc)}), encoding="utf-8")
            return 1
    if sys.platform != "win32" or platform.machine().lower() not in {"amd64", "x86_64"}:
        raise RuntimeError("Ce lanceur est destiné à Windows 10/11 x64 avec WSL2.")
    import tkinter as tk
    from tkinter import messagebox, ttk
    import msvcrt

    state = windows_state_dir()
    state.mkdir(parents=True, exist_ok=True)
    lock = (state / "launcher.lock").open("a+b")
    if lock.tell() == 0:
        lock.write(b"0")
        lock.flush()
    lock.seek(0)
    try:
        msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
    except OSError:
        messagebox.showinfo("Apta2Switch Studio", "Le lanceur est déjà ouvert.")
        lock.close()
        return 0
    window = tk.Tk()
    window.title("Apta2Switch Studio")
    window.geometry("620x330")
    window.resizable(False, False)
    frame = ttk.Frame(window, padding=24)
    frame.pack(fill="both", expand=True)
    ttk.Label(frame, text="Apta2Switch Studio", font=("Segoe UI", 20, "bold")).pack(anchor="w")
    ttk.Label(frame, text="L’application s’ouvre dans votre navigateur.\nLes calculs restent sur cet ordinateur.",
              padding=(0, 12)).pack(anchor="w")
    status = tk.StringVar(value="Ouvrez l’application, ou préparez Windows au premier lancement.")
    ttk.Label(frame, textvariable=status, wraplength=570).pack(anchor="w", pady=8)
    progress = ttk.Progressbar(frame, mode="indeterminate")
    progress.pack(fill="x", pady=10)
    messages = queue.Queue()
    backend = WindowsBackend(package_root(), lambda text: messages.put(("status", text)))
    busy = False

    def submit(action, kind):
        nonlocal busy
        if busy:
            return
        busy = True
        for button in buttons:
            button.configure(state="disabled")
        progress.start()
        def work():
            try:
                result = action()
                messages.put((kind, result))
            except Exception as exc:
                try:
                    (state / "launcher-error.log").write_text(str(exc), encoding="utf-8")
                except OSError:
                    pass
                messages.put(("error", str(exc)))
        threading.Thread(target=work, daemon=True).start()

    def setup():
        if messagebox.askokcancel("Préparer Windows", "Windows va installer WSL2 et Ubuntu 22.04 via Microsoft. "
                "Internet, une autorisation administrateur et parfois un redémarrage sont nécessaires. "
                "Si Ubuntu 22.04 utilise déjà WSL1, il sera converti vers WSL2 ; sauvegardez ses données avant de continuer."):
            status.set("Préparation Windows en cours… Acceptez la demande Windows si elle apparaît.")
            submit(backend.setup, "setup")

    actions = ttk.Frame(frame)
    actions.pack(fill="x", pady=8)
    buttons = [ttk.Button(actions, text="Ouvrir l’application", command=lambda: submit(backend.launch, "opened")),
               ttk.Button(actions, text="Préparer Windows", command=setup)]
    for button in buttons:
        button.pack(side="left", padx=(0, 12))
    ttk.Label(frame, text="Les moteurs scientifiques sont inclus et prêts à utiliser.",
              wraplength=570).pack(anchor="w", pady=8)

    def poll():
        nonlocal busy
        try:
            while True:
                kind, value = messages.get_nowait()
                if kind == "status":
                    status.set(value)
                    continue
                busy = False
                progress.stop()
                for button in buttons:
                    button.configure(state="normal")
                if kind == "opened":
                    webbrowser.open(value)
                    status.set("Application ouverte : " + value + ". Vous pouvez fermer cette fenêtre ; les calculs continuent.")
                elif kind == "setup":
                    status.set(value)
                else:
                    status.set(value)
                    messagebox.showerror("Apta2Switch Studio", value)
        except queue.Empty:
            pass
        window.after(150, poll)

    def close():
        if busy:
            messagebox.showinfo("Apta2Switch Studio", "Attendez la fin de la préparation avant de fermer cette fenêtre.")
        else:
            window.destroy()
    window.protocol("WM_DELETE_WINDOW", close)
    window.after(150, poll)
    try:
        window.mainloop()
    finally:
        lock.seek(0)
        msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
        lock.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
