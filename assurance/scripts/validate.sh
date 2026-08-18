#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${ROOT}/build-assurance-validation"
OUTPUT="${ROOT}/assurance/output"
COMPILER="${CC:-cc}"

rm -rf "${BUILD}"
rm -f "${OUTPUT}"/*.json

cmake -S "${ROOT}/assurance" -B "${BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${COMPILER}" \
  -DSPACEWIND_ASSURANCE_WERROR=ON
cmake --build "${BUILD}" --parallel 2
ctest --test-dir "${BUILD}" --output-on-failure

for receipt in "${OUTPUT}"/*.json; do
  python3 -m json.tool "${receipt}" >/dev/null
done

find "${OUTPUT}" -maxdepth 1 -name '*.json' -print0 | sort -z | \
  xargs -0 sha256sum > "${BUILD}/receipts-first.sha256"

rm -f "${OUTPUT}"/*.json
ctest --test-dir "${BUILD}" --output-on-failure
find "${OUTPUT}" -maxdepth 1 -name '*.json' -print0 | sort -z | \
  xargs -0 sha256sum > "${BUILD}/receipts-second.sha256"
diff -u "${BUILD}/receipts-first.sha256" "${BUILD}/receipts-second.sha256"

python3 - "${OUTPUT}/current_claim_gate.json" <<'PY'
import json
import pathlib
import sys

receipt = json.loads(pathlib.Path(sys.argv[1]).read_text())
if receipt["propulsion_claim_ready"] is not False:
    raise SystemExit("current project must not be promoted to propulsion-ready")
if receipt["flight_propulsion_ready"] is not False:
    raise SystemExit("current project must not be promoted to flight propulsion")
print("claim gate remains fail-closed:", receipt["highest_tier"])
PY

printf 'assurance validation passed with %s\n' "${COMPILER}"
