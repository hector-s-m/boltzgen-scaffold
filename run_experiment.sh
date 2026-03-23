#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
#  CONFIGURATION — Edit these parameters before running
# =============================================================================

# --- Input ---
INPUT_PDB="input/target+motif_2.pdb"

# --- Scaffold ---
SCAFFOLD="example/nanobody_scaffolds/7eow.yaml"

# --- Output ---
DESIGNS_DIR="designs/experiment"
RESULTS_DIR="results/experiment"

# --- Design parameters ---
FLANK_N="2..9"
FLANK_C="2..9"
LINKER_LENGTH="1..6"
MOTIF_NOISE=0.2
NUM_DESIGNS=1000                              # Per config (6 configs × 250 = 1500 total)
PROTOCOL="nanobody-anything"

# =============================================================================
#  HELPER FUNCTION
# =============================================================================

COMBINED_REGISTRY="$RESULTS_DIR/design_registry_all.csv"

run_config() {
    local name="$1"; shift
    local yaml="$DESIGNS_DIR/${name}.yaml"

    echo ""
    echo "========================================"
    echo "  Config: $name"
    echo "========================================"

    python graft_motif.py \
        --input "$INPUT_PDB" \
        --scaffold "$SCAFFOLD" \
        --output "$yaml" \
        --flank-n "$FLANK_N" --flank-c "$FLANK_C" \
        --linker-length "$LINKER_LENGTH" \
        --motif-noise "$MOTIF_NOISE" \
        "$@"

    boltzgen run "$yaml" \
        --protocol "$PROTOCOL" \
        --output "$RESULTS_DIR/$name" \
        --num_designs "$NUM_DESIGNS"

    # Append this config's registry to the combined file
    local registry="$RESULTS_DIR/$name/final_ranked_designs/design_registry.csv"
    if [[ -f "$registry" ]]; then
        if [[ ! -f "$COMBINED_REGISTRY" ]]; then
            head -1 "$registry" | sed 's/$/,config/' > "$COMBINED_REGISTRY"
        fi
        tail -n +2 "$registry" | sed "s/$/,$name/" >> "$COMBINED_REGISTRY"
        echo "  -> Registry: $(tail -n +2 "$COMBINED_REGISTRY" | wc -l | tr -d ' ') total designs"
    fi
}

# =============================================================================
#  RUN ALL CONFIGURATIONS
# =============================================================================

mkdir -p "$DESIGNS_DIR" "$RESULTS_DIR"
rm -f "$COMBINED_REGISTRY"

echo "================================================"
echo "  Multi-CDR Grafting Experiment"
echo "================================================"
echo "Input:       $INPUT_PDB"
echo "Scaffold:    $SCAFFOLD"
echo "Designs:     $NUM_DESIGNS per config (6 configs)"
echo "Motif noise: $MOTIF_NOISE Å"

# Config 1: Both residues (PHE + TYR) grafted to CDR3
run_config "both_cdr3" \
    --cdr-motif 3:B,C

# Config 2: Both residues (PHE + TYR) grafted to CDR2
run_config "both_cdr2" \
    --cdr-motif 2:B,C

# Config 3a: PHE → CDR2, TYR → CDR3
run_config "split_B2_C3" \
    --cdr-motif 2:B 3:C

# Config 3b: TYR → CDR2, PHE → CDR3
run_config "split_C2_B3" \
    --cdr-motif 2:C 3:B

# Config 4a: PHE → CDR1, TYR → CDR3
run_config "split_B1_C3" \
    --cdr-motif 1:B 3:C

# Config 4b: TYR → CDR1, PHE → CDR3
run_config "split_C1_B3" \
    --cdr-motif 1:C 3:B

echo ""
echo "================================================"
echo "  All 6 configs complete."
echo "  Results:  $RESULTS_DIR/"
echo "  Registry: $COMBINED_REGISTRY"
echo "================================================"
