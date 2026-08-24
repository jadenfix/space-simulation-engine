#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

make clean
make all
make test
make examples

for cfg in dynamic_soaring_shear cme_encounter magnetic_sail hybrid_glider; do
  ./build/spacewind simulate \
    "configs/${cfg}.cfg" \
    "output/${cfg}.csv" \
    "output/${cfg}.receipt.json"
done

for receipt in output/*.receipt.json; do
  if command -v python3 >/dev/null 2>&1; then
    python3 -m json.tool "$receipt" >/dev/null
  fi
done

printf '%s\n' "release validation completed"
