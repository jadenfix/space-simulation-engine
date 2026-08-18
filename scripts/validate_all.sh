#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

STRICT='-std=c11 -O2 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Werror'

printf '%s\n' '[1/7] release build, tests, examples, and receipts'
./scripts/validate.sh

if command -v clang >/dev/null 2>&1; then
  printf '%s\n' '[2/7] Clang warning-as-error build and static analysis'
  make clean
  make CC=clang CFLAGS="$STRICT" all test
  ANALYZE_DIR=${TMPDIR:-/tmp}/spacewind-clang-analyze
  rm -rf "$ANALYZE_DIR"
  mkdir -p "$ANALYZE_DIR"
  for source in src/*.c tests/test_main.c; do
    clang --analyze -std=c11 -Iinclude \
      -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes \
      "$source" -Xanalyzer -analyzer-output=text \
      -o "$ANALYZE_DIR/$(basename "$source" .c).plist"
  done
else
  printf '%s\n' '[2/7] Clang unavailable; skipped'
fi

if command -v gcc >/dev/null 2>&1; then
  printf '%s\n' '[3/7] GCC warning-as-error build and analyzer'
  make clean
  make CC=gcc CFLAGS="$STRICT" all test
  ANALYZE_DIR=${TMPDIR:-/tmp}/spacewind-gcc-analyze
  rm -rf "$ANALYZE_DIR"
  mkdir -p "$ANALYZE_DIR"
  if printf 'int main(void){return 0;}\n' | gcc -x c -fanalyzer -c -o "$ANALYZE_DIR/probe.o" - 2>/dev/null; then
    for source in src/*.c tests/test_main.c; do
      gcc -std=c11 -O0 -g -Iinclude \
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes \
        -fanalyzer -Werror -c "$source" \
        -o "$ANALYZE_DIR/$(basename "$source" .c).o"
    done
  fi
else
  printf '%s\n' '[3/7] GCC unavailable; skipped'
fi

if command -v cmake >/dev/null 2>&1; then
  printf '%s\n' '[4/7] independent CMake/CTest build'
  CMAKE_DIR=${TMPDIR:-/tmp}/spacewind-cmake-validation
  rm -rf "$CMAKE_DIR"
  cmake -S . -B "$CMAKE_DIR" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$CMAKE_DIR"
  ctest --test-dir "$CMAKE_DIR" --output-on-failure
else
  printf '%s\n' '[4/7] CMake unavailable; skipped'
fi

printf '%s\n' '[5/7] independent physical-assurance kernel'
make assurance

if [ "${SPACEWIND_SKIP_SANITIZERS:-0}" = 1 ]; then
  printf '%s\n' '[6/7] sanitizers skipped by SPACEWIND_SKIP_SANITIZERS=1'
else
  printf '%s\n' '[6/7] address and undefined-behavior sanitizers'
  ./scripts/sanitize.sh
fi

printf '%s\n' '[7/7] restore and rerun the release artifacts'
./scripts/validate.sh
printf '%s\n' 'all available validation layers completed'
