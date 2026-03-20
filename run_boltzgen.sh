#!/usr/bin/env bash
set -euo pipefail

NUM_DESIGNS=50
PROTOCOL="nanobody-anything"
RESULTS_DIR="results"

for yaml in designs/graft_*.yaml; do
  name=$(basename "$yaml" .yaml)
  echo "=== Running: $name ==="
  boltzgen run "$yaml" \
    --protocol "$PROTOCOL" \
    --output "$RESULTS_DIR/$name" \
    --num_designs "$NUM_DESIGNS"
done

echo "=== All done ==="
