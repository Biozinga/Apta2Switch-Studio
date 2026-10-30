#!/usr/bin/env bash
# Native build. The portable archive needs neither Python nor FUSE.
set -euo pipefail
cd "$(dirname "$0")/.."

APTASWITCH_CHECK_ENV=0
if [[ "$#" -gt 0 ]]; then
  if [[ "$#" == 1 && "$1" == "--check-env" ]]; then
    APTASWITCH_CHECK_ENV=1
  else
    echo "Usage: bash scripts/build_linux.sh [--check-env]" >&2
    exit 1
  fi
fi

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
  echo "This build must run on Linux x86_64." >&2
  exit 1
fi

source scripts/build_environment.sh
aptaswitch_check_build_tools
aptaswitch_prepare_build_python "$APTASWITCH_CHECK_ENV"
if [[ "$APTASWITCH_CHECK_ENV" == 1 ]]; then
  echo "Linux x86_64 build environment verified (CPython 3.11). No packages installed or binaries built."
  exit 0
fi
"$PYTHON_BIN" -m pip install -e ".[build]"
"$PYTHON_BIN" packaging/make_icons.py
"$PYTHON_BIN" scripts/prepare_thermodynamics_runtime.py
"$PYTHON_BIN" scripts/prepare_bundled_science.py

BUILD_ROOT="$(mktemp -d -t aptaswitch-beta-build.XXXXXXXX)"
trap 'rm -rf "$BUILD_ROOT"' EXIT
"$PYTHON_BIN" -m PyInstaller --noconfirm \
  --workpath "$BUILD_ROOT/build" --distpath "$BUILD_ROOT/dist" \
  packaging/aptaswitch-studio.spec
APP_PATH="$BUILD_ROOT/dist/Apta2Switch-Studio"
if [[ ! -x "$APP_PATH/Apta2Switch-Studio" ]]; then
  echo "The Linux executable was not produced." >&2
  exit 1
fi
mkdir -p "$APP_PATH/runs/designs" "$APP_PATH/runs/extensions"
cp LICENSE "$APP_PATH/LICENSE"
cp packaging/START-HERE.txt "$APP_PATH/START-HERE.txt"
"$PYTHON_BIN" scripts/check_bundled_app.py "$APP_PATH/Apta2Switch-Studio"

# tar preserves the executable bit, unlike downloading a loose AppImage.
mkdir -p dist
tar -C "$BUILD_ROOT/dist" -czf dist/Apta2Switch-Studio-linux-x64.tar.gz Apta2Switch-Studio
"$PYTHON_BIN" scripts/check_bundled_app.py --archive dist/Apta2Switch-Studio-linux-x64.tar.gz
# Refresh application files while preserving saved runs.
mkdir -p dist/Apta2Switch-Studio
rm -rf dist/Apta2Switch-Studio/_internal
cp -R "$APP_PATH"/. dist/Apta2Switch-Studio/

echo "Archive: dist/Apta2Switch-Studio-linux-x64.tar.gz"
