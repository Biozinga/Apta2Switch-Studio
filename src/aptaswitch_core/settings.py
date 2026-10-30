"""Persistent application settings shared by independent optional runtimes."""
from __future__ import annotations

import json
from pathlib import Path

SETTINGS_DIR = Path.home() / ".aptaswitch-studio"
SETTINGS_PATH = SETTINGS_DIR / "settings.json"


def load_settings() -> dict:
    if not SETTINGS_PATH.exists():
        return {}
    try:
        return json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def save_settings(settings: dict) -> None:
    SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    SETTINGS_PATH.write_text(json.dumps(settings, indent=2), encoding="utf-8")
