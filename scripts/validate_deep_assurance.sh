#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_ROOT="${ROOT}/build-deep-assurance"
CC_VALUE="${CC:-cc}"

rm -rf "${BUILD_ROOT}"
mkdir -p "${BUILD_ROOT}"

printf '\n== Original Spacewind engine ==\n'
make -C "${ROOT}" clean
make -C "${ROOT}" CC="${CC_VALUE}"
make -C "${ROOT}" CC="${CC_VALUE}" test
make -C "${ROOT}" CC="${CC_VALUE}" examples

printf '\n== Independent assurance kernel ==\n'
cmake -S "${ROOT}/assurance" -B "${BUILD_ROOT}/assurance" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${CC_VALUE}" \
  -DSPACEWIND_ASSURANCE_WERROR=ON
cmake --build "${BUILD_ROOT}/assurance" --parallel 2
ctest --test-dir "${BUILD_ROOT}/assurance" --output-on-failure

printf '\n== Manufactured PDE and particle checks ==\n'
cmake -S "${ROOT}/assurance/pde" -B "${BUILD_ROOT}/pde" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${CC_VALUE}"
cmake --build "${BUILD_ROOT}/pde" --parallel 2
ctest --test-dir "${BUILD_ROOT}/pde" --output-on-failure

printf '\n== Adversarial robustness checks ==\n'
cmake -S "${ROOT}/assurance/robustness" -B "${BUILD_ROOT}/robustness" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${CC_VALUE}"
cmake --build "${BUILD_ROOT}/robustness" --parallel 2
ctest --test-dir "${BUILD_ROOT}/robustness" --output-on-failure

printf '\n== Receipt validation ==\n'
for receipt in "${ROOT}/assurance/output"/*.json; do
  python3 -m json.tool "${receipt}" >/dev/null
done
python3 -m json.tool "${ROOT}/claims/spacewind-assurance-claims.json" >/dev/null

python3 - "${ROOT}/assurance/output/current_claim_gate.json" <<'PY'
import json
import pathlib
import sys

receipt = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert receipt["propulsion_claim_ready"] is False
assert receipt["flight_propulsion_ready"] is False
print("Current highest claim tier:", receipt["highest_tier"])
print("Propulsion claim remains blocked as required.")
PY

printf '\nDeep assurance validation passed with compiler: %s\n' "${CC_VALUE}"
