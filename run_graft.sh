#!/usr/bin/env bash
set -euo pipefail

# === Configuration ===
INPUT_PDB="input/target+motif.pdb"
SCAFFOLD_DIR="example/nanobody_scaffolds"
OUTPUT_DIR="designs"
RESULTS_DIR="results"
CDR=3
FLANK_N="2..5"
FLANK_C="2..5"

# === Run motif grafting for all scaffolds ===
SCAFFOLDS=("$SCAFFOLD_DIR"/*.yaml)

echo "Input:      $INPUT_PDB"
echo "Scaffolds:  ${#SCAFFOLDS[@]} found in $SCAFFOLD_DIR"
echo "Output:     $OUTPUT_DIR/"
echo "CDR target: $CDR"
echo ""

python graft_motif.py \
    --input "$INPUT_PDB" \
    --scaffold "${SCAFFOLDS[@]}" \
    --output "$OUTPUT_DIR/graft.yaml" \
    --cdr "$CDR" \
    --flank-n "$FLANK_N" \
    --flank-c "$FLANK_C"

echo ""
echo "Done. To run BoltzGen on the generated designs:"
echo "  boltzgen run $OUTPUT_DIR/<design>.yaml --protocol nanobody-anything --output $RESULTS_DIR/"
