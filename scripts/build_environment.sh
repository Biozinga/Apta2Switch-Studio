#!/usr/bin/env bash
# Shared by the native desktop builders; never changes the application runtime.

aptaswitch_check_build_tools() {
  local tool
  for tool in cmake make "${CXX:-c++}" tar "$@"; do
    if ! command -v "$tool" >/dev/null 2>&1; then
      echo "Missing build tool: $tool. See README.md#run-from-source." >&2
      return 1
    fi
  done
}

aptaswitch_check_python_version() {
  "$1" -c 'import platform, sys; sys.exit(0 if platform.python_implementation() == "CPython" and sys.version_info[:2] == (3, 11) else "Use CPython 3.11 for the pinned native build dependencies.")'
}

aptaswitch_prepare_build_python() {
  local check_only="${1:-0}"
  local build_python="${PYTHON_BIN:-.build-venv/bin/python}"
  local bootstrap_python="${APTASWITCH_PYTHON:-python3.11}"

  if command -v "$build_python" >/dev/null 2>&1; then
    aptaswitch_check_python_version "$build_python" || return 1
  elif [[ -n "${PYTHON_BIN:-}" ]]; then
    echo "PYTHON_BIN does not name an executable: $PYTHON_BIN" >&2
    return 1
  else
    if ! command -v "$bootstrap_python" >/dev/null 2>&1; then
      echo "Install CPython 3.11 or set APTASWITCH_PYTHON to its executable." >&2
      return 1
    fi
    aptaswitch_check_python_version "$bootstrap_python" || return 1
    if [[ "$check_only" == "1" ]]; then
      build_python="$bootstrap_python"
    else
      "$bootstrap_python" -m venv .build-venv || return 1
    fi
  fi
  PYTHON_BIN="$build_python"
}
