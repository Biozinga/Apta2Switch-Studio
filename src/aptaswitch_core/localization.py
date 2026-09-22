"""Application text translations, independent of scientific data and UI state."""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from functools import lru_cache
import re
from string import Formatter

from .locales.en_core import TRANSLATIONS as CORE_TRANSLATIONS
from .locales.en_components import TRANSLATIONS as COMPONENT_TRANSLATIONS
from .locales.en_main import TRANSLATIONS as MAIN_TRANSLATIONS


LANGUAGES = {"en": "English", "fr": "Français"}
_language = ContextVar("aptaswitch_language", default="en")


def get_language() -> str:
    return _language.get()


def set_language(language: str) -> None:
    _language.set(language if language in LANGUAGES else "en")


@contextmanager
def language_context(language: str):
    token = _language.set(language if language in LANGUAGES else "en")
    try:
        yield
    finally:
        _language.reset(token)


@lru_cache(maxsize=1)
def _catalogs():
    english = {}
    for catalog in (CORE_TRANSLATIONS, COMPONENT_TRANSLATIONS, MAIN_TRANSLATIONS):
        english.update(catalog)
    return english, {value: key for key, value in english.items()}


@lru_cache(maxsize=2)
def _message_patterns(language: str):
    """Recognize formatted messages from saved runs without changing their data."""
    english, french = _catalogs()
    mapping = english if language == "en" else french
    patterns = []
    formatter = Formatter()
    for source, target in mapping.items():
        if source == target:
            continue
        try:
            parts = list(formatter.parse(source))
        except ValueError:
            continue
        if not any(field is not None for _, field, _, _ in parts):
            continue
        if not any(literal.strip() for literal, _, _, _ in parts):
            continue
        fields, pattern = [], ""
        for literal, field, _, _ in parts:
            pattern += re.escape(literal)
            if field is not None:
                if field in fields:
                    pattern += f"(?P=p{fields.index(field)})"
                else:
                    pattern += f"(?P<p{len(fields)}>.*?)"
                    fields.append(field)
        patterns.append((parts[0][0], re.compile(pattern, re.DOTALL), fields, target))
    return patterns


@lru_cache(maxsize=4096)
def _translate_message(text: str, language: str) -> str:
    timestamp = re.match(r"^(\[\d{2}:\d{2}:\d{2}\] )(.*)$", text, re.DOTALL)
    if timestamp:
        return timestamp[1] + _translate_message(timestamp[2], language)
    english, french = _catalogs()
    mapping = english if language == "en" else french
    if text in mapping:
        return mapping[text]
    # Recognize templates whose parameters were formatted by a worker or an
    # older saved run. Captured values are opaque: names, paths and sequences
    # must never be translated or interpreted as templates themselves.
    for prefix, pattern, fields, target in _message_patterns(language):
        if prefix and not text.startswith(prefix):
            continue
        match = pattern.fullmatch(text)
        if match:
            values = {field: match.group(f"p{i}") for i, field in enumerate(fields)}
            return "".join(literal + (values[field] if field is not None else "")
                           for literal, field, _, _ in Formatter().parse(target))
    return text


def tr(text: str, **values) -> str:
    """Translate a display string, formatting only explicitly supplied values."""
    if not isinstance(text, str):
        return text
    translated = _translate_message(text, get_language())
    return translated.format(**values) if values else translated
