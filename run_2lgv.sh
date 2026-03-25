#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#  2LGV Hotspot Grafting Experiment
#  Target: chain A  |  Hotspots: B(PHE), C(PHE), D(TRP)
# ============================================================

# --- Tunable parameters ---
INPUT_PDB="input/2LGV_des.pdb"
SCAFFOLD="example/nanobody_scaffolds/7eow.yaml"
NUM_DESIGNS=50              # per configuration
PROTOCOL="nanobody-anything"
FLANK_N="2..5"
FLANK_C="2..5"
LINKER_LENGTH="1..6"
MOTIF_NOISE=0.3             # σ (Å) for side-chain soft flexibility

# --- Configurations ---
# 1. Both PHE (B,C) + TRP (D) all in CDR3
# 2. Both PHE (B,C) + TRP (D) all in CDR2
# 3. TRP (D) in CDR2, both PHE (B,C) in CDR3
# 4. TRP (D) in CDR1, both PHE (B,C) in CDR3

RESULTS_BASE="results/2lgv"

echo "================================================"
echo "  2LGV Hotspot Grafting Experiment"
echo "================================================"
echo "Input:     $INPUT_PDB"
echo "Scaffold:  $SCAFFOLD"
echo "Designs:   $NUM_DESIGNS per config"
echo ""

# ----------------------------------------------------------
# Config 1: All hotspots → CDR3
# ----------------------------------------------------------
CONFIG="config1_all_cdr3"
echo "--- Config 1: B,C,D → CDR3 ---"
python graft_motif.py \
    --input "$INPUT_PDB" \
    --scaffold "$SCAFFOLD" \
    --output "designs/2lgv_${CONFIG}.yaml" \
    --cdr-motif 3:B,C,D \
    --flank-n "$FLANK_N" --flank-c "$FLANK_C" \
    --linker-length "$LINKER_LENGTH" \
    --motif-noise "$MOTIF_NOISE"

boltzgen run "designs/2lgv_${CONFIG}.yaml" \
    --protocol "$PROTOCOL" \
    --output "${RESULTS_BASE}/${CONFIG}" \
    --num_designs "$NUM_DESIGNS"

echo ""

# ----------------------------------------------------------
# Config 2: All hotspots → CDR2
# ----------------------------------------------------------
CONFIG="config2_all_cdr2"
echo "--- Config 2: B,C,D → CDR2 ---"
python graft_motif.py \
    --input "$INPUT_PDB" \
    --scaffold "$SCAFFOLD" \
    --output "designs/2lgv_${CONFIG}.yaml" \
    --cdr-motif 2:B,C,D \
    --flank-n "$FLANK_N" --flank-c "$FLANK_C" \
    --linker-length "$LINKER_LENGTH" \
    --motif-noise "$MOTIF_NOISE"

boltzgen run "designs/2lgv_${CONFIG}.yaml" \
    --protocol "$PROTOCOL" \
    --output "${RESULTS_BASE}/${CONFIG}" \
    --num_designs "$NUM_DESIGNS"

echo ""

# ----------------------------------------------------------
# Config 3: TRP (D) → CDR2, both PHE (B,C) → CDR3
# ----------------------------------------------------------
CONFIG="config3_D_cdr2_BC_cdr3"
echo "--- Config 3: D → CDR2, B,C → CDR3 ---"
python graft_motif.py \
    --input "$INPUT_PDB" \
    --scaffold "$SCAFFOLD" \
    --output "designs/2lgv_${CONFIG}.yaml" \
    --cdr-motif 2:D 3:B,C \
    --flank-n "$FLANK_N" --flank-c "$FLANK_C" \
    --linker-length "$LINKER_LENGTH" \
    --motif-noise "$MOTIF_NOISE"

boltzgen run "designs/2lgv_${CONFIG}.yaml" \
    --protocol "$PROTOCOL" \
    --output "${RESULTS_BASE}/${CONFIG}" \
    --num_designs "$NUM_DESIGNS"

echo ""

# ----------------------------------------------------------
# Config 4: TRP (D) → CDR1, both PHE (B,C) → CDR3
# ----------------------------------------------------------
CONFIG="config4_D_cdr1_BC_cdr3"
echo "--- Config 4: D → CDR1, B,C → CDR3 ---"
python graft_motif.py \
    --input "$INPUT_PDB" \
    --scaffold "$SCAFFOLD" \
    --output "designs/2lgv_${CONFIG}.yaml" \
    --cdr-motif 1:D 3:B,C \
    --flank-n "$FLANK_N" --flank-c "$FLANK_C" \
    --linker-length "$LINKER_LENGTH" \
    --motif-noise "$MOTIF_NOISE"

boltzgen run "designs/2lgv_${CONFIG}.yaml" \
    --protocol "$PROTOCOL" \
    --output "${RESULTS_BASE}/${CONFIG}" \
    --num_designs "$NUM_DESIGNS"

echo ""
echo "================================================"
echo "  All 4 configurations complete"
echo "  Results in: ${RESULTS_BASE}/"
echo "================================================"
