#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#  2LGV Hotspot Grafting Experiment
#  Target: chain A  |  Hotspots: B(PHE), C(PHE), D(TRP)
# ============================================================

# --- Tunable parameters ---
INPUT_PDB="input/2LGV_des.pdb"
SCAFFOLD="example/nanobody_scaffolds/7eow.yaml"
NUM_DESIGNS=100             # per configuration
PROTOCOL="nanobody-anything"
FLANK_N="2..5"
FLANK_C="2..5"
LINKER_LENGTH="1..6"
MOTIF_NOISE=0.3             # σ (Å) for side-chain soft flexibility
TARGET_RESIDUES="53,57,60,61,75"  # PHE53, PHE75, SER57, LYS61, LEU60

# --- Configurations ---
# 1. All hotspots → CDR3
# 2. All hotspots → CDR2
# 3. TRP (D) → CDR2, both PHE (B,C) → CDR3
# 4. TRP (D) → CDR1, both PHE (B,C) → CDR3
# 5. B(PHE) → CDR1, C(PHE) → CDR2, D(TRP) → CDR3
# 6. B(PHE) → CDR3, C(PHE) → CDR2, D(TRP) → CDR1

RESULTS_BASE="results/2lgv"

COMMON_ARGS=(
    --input "$INPUT_PDB"
    --scaffold "$SCAFFOLD"
    --flank-n "$FLANK_N" --flank-c "$FLANK_C"
    --linker-length "$LINKER_LENGTH"
    --motif-noise "$MOTIF_NOISE"
    --target-residues "$TARGET_RESIDUES"
)

echo "================================================"
echo "  2LGV Hotspot Grafting Experiment"
echo "================================================"
echo "Input:           $INPUT_PDB"
echo "Scaffold:        $SCAFFOLD"
echo "Designs:         $NUM_DESIGNS per config"
echo "Target residues: $TARGET_RESIDUES"
echo ""

run_config() {
    local config="$1"
    local label="$2"
    shift 2
    echo "--- ${label} ---"
    python graft_motif.py \
        "${COMMON_ARGS[@]}" \
        --output "designs/2lgv_${config}.yaml" \
        "$@"

    boltzgen run "designs/2lgv_${config}.yaml" \
        --protocol "$PROTOCOL" \
        --output "${RESULTS_BASE}/${config}" \
        --num_designs "$NUM_DESIGNS"
    echo ""
}

# Config 1: All hotspots → CDR3
run_config "config1_all_cdr3" "Config 1: B,C,D → CDR3" \
    --cdr-motif 3:B,C,D

# Config 2: All hotspots → CDR2
run_config "config2_all_cdr2" "Config 2: B,C,D → CDR2" \
    --cdr-motif 2:B,C,D

# Config 3: TRP (D) → CDR2, both PHE (B,C) → CDR3
run_config "config3_D_cdr2_BC_cdr3" "Config 3: D → CDR2, B,C → CDR3" \
    --cdr-motif 2:D 3:B,C

# Config 4: TRP (D) → CDR1, both PHE (B,C) → CDR3
run_config "config4_D_cdr1_BC_cdr3" "Config 4: D → CDR1, B,C → CDR3" \
    --cdr-motif 1:D 3:B,C

# Config 5: Each aromatic in its own CDR — B→CDR1, C→CDR2, D→CDR3
run_config "config5_B_cdr1_C_cdr2_D_cdr3" "Config 5: B → CDR1, C → CDR2, D → CDR3" \
    --cdr-motif 1:B 2:C 3:D

# Config 6: Each aromatic in its own CDR — D→CDR1, C→CDR2, B→CDR3
run_config "config6_D_cdr1_C_cdr2_B_cdr3" "Config 6: D → CDR1, C → CDR2, B → CDR3" \
    --cdr-motif 1:D 2:C 3:B

echo "================================================"
echo "  All 6 configurations complete"
echo "  Results in: ${RESULTS_BASE}/"
echo "================================================"
