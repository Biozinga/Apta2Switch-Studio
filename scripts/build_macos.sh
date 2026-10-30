#!/usr/bin/env bash
# Build a standalone application locally or a signed maintainer release.
set -euo pipefail
cd "$(dirname "$0")/.."

APTASWITCH_CHECK_ENV=0
if [[ "$#" -gt 0 ]]; then
  if [[ "$#" == 1 && "$1" == "--check-env" ]]; then
    APTASWITCH_CHECK_ENV=1
  else
    echo "Usage: bash scripts/build_macos.sh [--check-env]" >&2
    exit 1
  fi
fi

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This build must run on macOS." >&2
  exit 1
fi

source scripts/build_environment.sh
aptaswitch_check_build_tools xcrun codesign xattr ditto
aptaswitch_prepare_build_python "$APTASWITCH_CHECK_ENV"
# Release is the safe default. Local development must opt into an ad-hoc build,
# which remains subject to Gatekeeper when downloaded on another Mac.
"$PYTHON_BIN" scripts/sign_macos_app.py --check-config
if [[ "$APTASWITCH_CHECK_ENV" == 1 ]]; then
  echo "macOS build environment verified (CPython 3.11). No packages installed or binaries built."
  exit 0
fi
"$PYTHON_BIN" -m pip install -e ".[build]"
"$PYTHON_BIN" packaging/make_icons.py
"$PYTHON_BIN" scripts/prepare_thermodynamics_runtime.py
"$PYTHON_BIN" scripts/prepare_bundled_science.py

ARCH="$("$PYTHON_BIN" -c 'import platform; print("arm64" if platform.machine() == "arm64" else "x64")')"
BUILD_ROOT="$(mktemp -d -t aptaswitch-beta-build)"
trap 'rm -rf "$BUILD_ROOT"' EXIT

"$PYTHON_BIN" -m PyInstaller --noconfirm \
  --workpath "$BUILD_ROOT/build" --distpath "$BUILD_ROOT/dist" \
  packaging/aptaswitch-studio.spec
APP_PATH="$BUILD_ROOT/dist/Apta2Switch-Studio.app"
if [[ ! -d "$APP_PATH" ]]; then
  echo "The macOS application was not produced." >&2
  exit 1
fi

# Assemble a fresh distribution folder; never copy user engines or runs.
PACKAGE_PATH="$BUILD_ROOT/package/Apta2Switch-Studio"
mkdir -p "$PACKAGE_PATH/runs/designs" "$PACKAGE_PATH/runs/extensions"
mv "$APP_PATH" "$PACKAGE_PATH/"
APP_PATH="$PACKAGE_PATH/Apta2Switch-Studio.app"
cp LICENSE "$PACKAGE_PATH/LICENSE"
cp packaging/START-HERE.txt "$PACKAGE_PATH/START-HERE.txt"

# Clean build metadata outside Desktop/iCloud before sealing the app. This is a
# maintainer build step, never a user-side quarantine-removal workaround.
xattr -cr "$APP_PATH"
"$PYTHON_BIN" scripts/sign_macos_app.py "$APP_PATH"
"$PYTHON_BIN" scripts/check_bundled_app.py "$APP_PATH/Contents/MacOS/Apta2Switch-Studio"

ARCHIVE_SUFFIX=""
if [[ "${APTASWITCH_MACOS_SIGNING:-release}" == "release" ]]; then
  "$PYTHON_BIN" scripts/sign_macos_app.py --notarize "$APP_PATH"
else
  ARCHIVE_SUFFIX="-local-unsigned"
fi
mkdir -p dist
ARCHIVE="dist/Apta2Switch-Studio-macos-${ARCH}${ARCHIVE_SUFFIX}.zip"
rm -f "$ARCHIVE"
# All app resources are ordinary files; omit Finder/AppleDouble metadata so the
# archive contains only the authenticated payload, including its symlinks.
ditto -c -k --norsrc --keepParent "$PACKAGE_PATH" "$ARCHIVE"

# Verify supporting folders and reject external engines, databases or user runs.
"$PYTHON_BIN" scripts/check_bundled_app.py --archive "$ARCHIVE"

# Keep the same complete layout locally, preserving any existing user data.
LOCAL_PACKAGE="${APTASWITCH_LOCAL_PACKAGE:-dist/Apta2Switch-Studio${ARCHIVE_SUFFIX}}"
mkdir -p "$LOCAL_PACKAGE/runs/designs" "$LOCAL_PACKAGE/runs/extensions"
rm -rf "$LOCAL_PACKAGE/Apta2Switch-Studio.app"
ditto --noextattr "$APP_PATH" "$LOCAL_PACKAGE/Apta2Switch-Studio.app"
cp "$PACKAGE_PATH/LICENSE" "$PACKAGE_PATH/START-HERE.txt" "$LOCAL_PACKAGE/"

echo "Application: $LOCAL_PACKAGE/Apta2Switch-Studio.app"
echo "Archive: $ARCHIVE"
