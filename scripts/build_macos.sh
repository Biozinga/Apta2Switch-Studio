#!/usr/bin/env bash
# Maintainer-only build: users download the resulting standalone application.
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This build must run on macOS." >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  python3 -m venv .venv
fi
"$PYTHON_BIN" -m pip install -e ".[build]"
"$PYTHON_BIN" packaging/make_icons.py

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

# Assemble a fresh distribution folder; never copy a user's NUPACK or runs.
PACKAGE_PATH="$BUILD_ROOT/package/Apta2Switch-Studio"
mkdir -p "$PACKAGE_PATH/nupack" "$PACKAGE_PATH/runs/designs" "$PACKAGE_PATH/runs/extensions"
mv "$APP_PATH" "$PACKAGE_PATH/"
APP_PATH="$PACKAGE_PATH/Apta2Switch-Studio.app"
cp LICENSE "$PACKAGE_PATH/LICENSE"
cp packaging/START-HERE.txt "$PACKAGE_PATH/START-HERE.txt"

# Sign and check outside Desktop/iCloud (Finder can add resource-fork metadata).
xattr -cr "$APP_PATH"
codesign --force --deep -s - "$APP_PATH"
codesign --verify --deep --strict "$APP_PATH"
"$PYTHON_BIN" scripts/check_bundled_app.py "$APP_PATH/Contents/MacOS/Apta2Switch-Studio"

mkdir -p dist
ARCHIVE="dist/Apta2Switch-Studio-macos-${ARCH}.zip"
rm -f "$ARCHIVE"
ditto -c -k --sequesterRsrc --keepParent "$PACKAGE_PATH" "$ARCHIVE"

# Fail the build if the downloadable ZIP loses the visible supporting folders.
"$PYTHON_BIN" - "$ARCHIVE" <<'PY'
import sys
from zipfile import ZipFile

root = "Apta2Switch-Studio/"
with ZipFile(sys.argv[1]) as archive:
    names = set(archive.namelist())
    expected = {
        root + "Apta2Switch-Studio.app/Contents/MacOS/Apta2Switch-Studio",
        root + "nupack/", root + "runs/designs/", root + "runs/extensions/",
        root + "START-HERE.txt", root + "LICENSE",
    }
    missing = expected - names
    if missing:
        raise SystemExit(f"Incomplete application archive: {sorted(missing)}")
    if any(not name.startswith((root, "__MACOSX/")) for name in names):
        raise SystemExit("The archive must contain one Apta2Switch-Studio folder")
    for folder in (root + "nupack/", root + "runs/"):
        if any(name.startswith(folder) and not name.endswith("/") for name in names):
            raise SystemExit(f"User data must not be included in {folder}")
print("Application archive structure verified")
PY

# Keep the same complete layout locally, preserving any existing user data.
LOCAL_PACKAGE="dist/Apta2Switch-Studio"
mkdir -p "$LOCAL_PACKAGE/nupack" "$LOCAL_PACKAGE/runs/designs" "$LOCAL_PACKAGE/runs/extensions"
rm -rf "$LOCAL_PACKAGE/Apta2Switch-Studio.app"
ditto --noextattr "$APP_PATH" "$LOCAL_PACKAGE/Apta2Switch-Studio.app"
cp "$PACKAGE_PATH/LICENSE" "$PACKAGE_PATH/START-HERE.txt" "$LOCAL_PACKAGE/"

echo "Application: $LOCAL_PACKAGE/Apta2Switch-Studio.app"
echo "Archive: $ARCHIVE"
