#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
VERSION=$(awk '/SPACEWIND_VERSION_STRING/ {gsub(/"/, "", $3); print $3; exit}' "$ROOT/include/spacewind/version.h")
DEST=${1:-"$ROOT/dist"}
NAME="spacewind-native-v${VERSION}"
STAGE=$(mktemp -d "${TMPDIR:-/tmp}/spacewind-package.XXXXXX")
trap 'rm -rf "$STAGE"' EXIT HUP INT TERM

mkdir -p "$DEST" "$STAGE/$NAME"
for path in \
  .github .gitignore CHANGELOG.md CMakeLists.txt LICENSE Makefile README.md \
  claims configs docs include scripts src tests; do
  cp -R "$ROOT/$path" "$STAGE/$NAME/"
done
mkdir -p "$STAGE/$NAME/output"
cp "$ROOT/output/.gitkeep" "$STAGE/$NAME/output/.gitkeep"

(
  cd "$STAGE"
  tar -czf "$DEST/$NAME.tar.gz" "$NAME"
  if command -v zip >/dev/null 2>&1; then
    zip -qr "$DEST/$NAME.zip" "$NAME"
  fi
)

if [ -f "$DEST/$NAME.zip" ]; then
  ARCHIVES="$NAME.tar.gz $NAME.zip"
else
  ARCHIVES="$NAME.tar.gz"
fi
if command -v sha256sum >/dev/null 2>&1; then
  (cd "$DEST" && sha256sum $ARCHIVES > "$NAME.sha256")
elif command -v shasum >/dev/null 2>&1; then
  (cd "$DEST" && shasum -a 256 $ARCHIVES > "$NAME.sha256")
fi

printf 'created source archives in %s\n' "$DEST"
