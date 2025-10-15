#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <path-to-app> <output-dmg>" >&2
  exit 1
fi

APP_PATH="$1"
DMG_PATH="$2"
VOL_NAME="$(basename "${APP_PATH}" .app)"
WORK_DIR="$(mktemp -d)"

cleanup() {
  rm -rf "${WORK_DIR}"
}
trap cleanup EXIT

mkdir -p "${WORK_DIR}/${VOL_NAME}"
cp -R "${APP_PATH}" "${WORK_DIR}/${VOL_NAME}/"
ln -s /Applications "${WORK_DIR}/${VOL_NAME}/Applications"

# Placeholder for optional codesign/notarization steps.
# Example:
# codesign --deep --force --options runtime --sign "Developer ID Application: ..." "${WORK_DIR}/${VOL_NAME}/${VOL_NAME}.app"

hdiutil create \
  -volname "${VOL_NAME}" \
  -srcfolder "${WORK_DIR}/${VOL_NAME}" \
  -format UDZO \
  -fs HFS+ \
  "${DMG_PATH}"

echo "Created DMG at ${DMG_PATH}"
