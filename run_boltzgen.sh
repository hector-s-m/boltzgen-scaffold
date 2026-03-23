#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
#  CONFIGURATION — Edit these parameters before running
# =============================================================================

# --- Input ---
INPUT_PDB="input/target+motif_2.pdb"        # PDB with target + motif chains (already docked)
TARGET_CHAIN=""                              # Target chain ID (leave empty to auto-detect as largest chain)
MOTIF_CHAIN=""                               # Motif chain ID (leave empty to auto-detect as second largest)

# --- Scaffolds ---
SCAFFOLD_DIR="example/nanobody_scaffolds"    # Directory containing nanobody scaffold YAMLs

# --- Output ---
DESIGNS_DIR="designs"                        # Where to write design YAML specs
RESULTS_DIR="results"                        # Where to write BoltzGen results

# --- Grafting ---
CDR=3                                        # Which CDR to graft motif into (1, 2, or 3)
FLANK_N="2..5"                               # N-terminal flank length range
FLANK_C="2..5"                               # C-terminal flank length range
LINKER_LENGTH=""                             # Linker length between disconnected hotspot fragments (default: 1..6)
MOTIF_NOISE=0.3                              # Gaussian noise σ (Å) on motif template coords (0=rigid, 0.3=slight flex, 1.0=loose)

# --- CDR lengths (leave empty for natural defaults: CDR1=5..13, CDR2=3..10, CDR3=10..25) ---
CDR1_LENGTH=""
CDR2_LENGTH=""
CDR3_LENGTH=""                               # Only used if this CDR is NOT the graft target

# --- BoltzGen pipeline ---
NUM_DESIGNS=50                               # Number of designs per scaffold
PROTOCOL="nanobody-anything"                 # BoltzGen protocol

# =============================================================================
#  STEP 1: Generate design YAMLs
# =============================================================================

SCAFFOLDS=("$SCAFFOLD_DIR"/*.yaml)

echo "================================================"
echo "  Motif Grafting Pipeline"
echo "================================================"
echo "Input:      $INPUT_PDB"
echo "Scaffolds:  ${#SCAFFOLDS[@]} found in $SCAFFOLD_DIR"
echo "CDR target: $CDR"
echo "Designs:    $NUM_DESIGNS per scaffold"
echo ""

# Build optional arguments
GRAFT_ARGS=()
[[ -n "$TARGET_CHAIN" ]]  && GRAFT_ARGS+=(--target-chain "$TARGET_CHAIN")
[[ -n "$MOTIF_CHAIN" ]]   && GRAFT_ARGS+=(--motif-chain "$MOTIF_CHAIN")
[[ -n "$CDR1_LENGTH" ]]   && GRAFT_ARGS+=(--cdr1-length "$CDR1_LENGTH")
[[ -n "$CDR2_LENGTH" ]]   && GRAFT_ARGS+=(--cdr2-length "$CDR2_LENGTH")
[[ -n "$CDR3_LENGTH" ]]   && GRAFT_ARGS+=(--cdr3-length "$CDR3_LENGTH")
[[ -n "$LINKER_LENGTH" ]] && GRAFT_ARGS+=(--linker-length "$LINKER_LENGTH")
[[ -n "$MOTIF_NOISE" ]]  && GRAFT_ARGS+=(--motif-noise "$MOTIF_NOISE")

echo "--- Step 1: Generating design YAMLs ---"
python graft_motif.py \
    --input "$INPUT_PDB" \
    --scaffold "${SCAFFOLDS[@]}" \
    --output "$DESIGNS_DIR/graft.yaml" \
    --cdr "$CDR" \
    --flank-n "$FLANK_N" \
    --flank-c "$FLANK_C" \
    "${GRAFT_ARGS[@]}"
echo ""

# =============================================================================
#  STEP 2: Run BoltzGen on each generated design
# =============================================================================

echo "--- Step 2: Running BoltzGen pipeline ---"
for yaml in "$DESIGNS_DIR"/graft_*.yaml; do
    [[ -f "$yaml" ]] || continue
    name=$(basename "$yaml" .yaml)
    echo ""
    echo "=== Running: $name ==="
    boltzgen run "$yaml" \
        --protocol "$PROTOCOL" \
        --output "$RESULTS_DIR/$name" \
        --num_designs "$NUM_DESIGNS"
done

echo ""
echo "================================================"
echo "  All done. Results in $RESULTS_DIR/"
echo "================================================"
