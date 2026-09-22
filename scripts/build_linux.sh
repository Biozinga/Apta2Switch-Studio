#!/usr/bin/env bash
# Maintainer-only build. The portable archive needs neither Python nor FUSE.
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
  echo "This build must run on Linux x86_64." >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  python3 -m venv .venv
fi
"$PYTHON_BIN" -m pip install -e ".[build]"
"$PYTHON_BIN" packaging/make_icons.py

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
mkdir -p "$APP_PATH/nupack" "$APP_PATH/runs/designs" "$APP_PATH/runs/extensions"
cp LICENSE "$APP_PATH/LICENSE"
cp packaging/START-HERE.txt "$APP_PATH/START-HERE.txt"
"$PYTHON_BIN" scripts/check_bundled_app.py "$APP_PATH/Apta2Switch-Studio"

# tar preserves the executable bit, unlike downloading a loose AppImage.
mkdir -p dist
tar -C "$BUILD_ROOT/dist" -czf dist/Apta2Switch-Studio-linux-x64.tar.gz Apta2Switch-Studio
# Refresh application files while preserving any local NUPACK and saved runs.
mkdir -p dist/Apta2Switch-Studio
rm -rf dist/Apta2Switch-Studio/_internal
cp -R "$APP_PATH"/. dist/Apta2Switch-Studio/

echo "Archive: dist/Apta2Switch-Studio-linux-x64.tar.gz"
