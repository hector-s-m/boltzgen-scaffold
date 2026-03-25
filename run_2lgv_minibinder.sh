#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#  2LGV Minibinder Design (de novo, no scaffold)
#  Target: chain A  |  Hotspots: B(PHE), C(PHE), D(TRP)
#  Fixed side chains with designed backbone and linkers
# ============================================================

# --- Tunable parameters ---
INPUT_PDB="input/2LGV_des.pdb"
NUM_DESIGNS=10              # set to 1000 for production
PROTOCOL="protein-anything"
NTAIL="15..30"              # N-terminal designed segment
CTAIL="15..30"              # C-terminal designed segment
LINKER_LENGTH="3..10"       # linkers between hotspot fragments
MOTIF_NOISE=0.3             # σ (Å) for template coordinate perturbation
TARGET_RESIDUES="53,57,60,61,75"  # PHE53, PHE75, SER57, LYS61, LEU60

RESULTS_BASE="results/2lgv_minibinder"
DESIGNS_DIR="designs"
YAML_OUT="${DESIGNS_DIR}/2lgv_minibinder.yaml"

echo "================================================"
echo "  2LGV Minibinder Design (de novo)"
echo "================================================"
echo "Input:           $INPUT_PDB"
echo "Binder size:     ~${NTAIL%%..* }+${CTAIL%%..* }+2*${LINKER_LENGTH%%..* }+3 to ~${NTAIL##*..}+${CTAIL##*..}+2*${LINKER_LENGTH##*..}+3 residues"
echo "Designs:         $NUM_DESIGNS"
echo "Target residues: $TARGET_RESIDUES"
echo ""

mkdir -p "$DESIGNS_DIR"

# Generate YAML via Python to avoid heredoc/quoting issues
python -c "
import os, yaml

input_rel = os.path.relpath('${INPUT_PDB}', '${DESIGNS_DIR}')

design = {
    'entities': [
        # N-terminal designed segment
        {'protein': {'id': 'B', 'sequence': '${NTAIL}'}},
        # Hotspot 1: PHE from chain B
        {'file': {
            'path': input_rel, 'fuse': 'B',
            'include': [{'chain': {'id': 'B'}}],
            'not_design': [{'chain': {'id': 'B'}}],
            'structure_groups': [{'group': {'id': 'B', 'visibility': 1}}],
        }},
        # Designed linker 1
        {'protein': {'id': 'LK1', 'fuse': 'B', 'sequence': '${LINKER_LENGTH}'}},
        # Hotspot 2: PHE from chain C
        {'file': {
            'path': input_rel, 'fuse': 'B',
            'include': [{'chain': {'id': 'C'}}],
            'not_design': [{'chain': {'id': 'C'}}],
            'structure_groups': [{'group': {'id': 'C', 'visibility': 1}}],
        }},
        # Designed linker 2
        {'protein': {'id': 'LK2', 'fuse': 'B', 'sequence': '${LINKER_LENGTH}'}},
        # Hotspot 3: TRP from chain D
        {'file': {
            'path': input_rel, 'fuse': 'B',
            'include': [{'chain': {'id': 'D'}}],
            'not_design': [{'chain': {'id': 'D'}}],
            'structure_groups': [{'group': {'id': 'D', 'visibility': 1}}],
        }},
        # C-terminal designed segment
        {'protein': {'id': 'CT', 'fuse': 'B', 'sequence': '${CTAIL}'}},
        # Target protein
        {'file': {
            'path': input_rel,
            'include': [{'chain': {'id': 'A'}}],
            'binding_types': [{'chain': {'id': 'A', 'binding': '${TARGET_RESIDUES}'}}],
        }},
    ],
    'motif_noise': ${MOTIF_NOISE},
}

with open('${YAML_OUT}', 'w') as f:
    yaml.dump(design, f, default_flow_style=False, sort_keys=False)
"

echo "Generated: ${YAML_OUT}"
echo ""

# Run BoltzGen
echo "--- Running BoltzGen minibinder design ---"
boltzgen run "${YAML_OUT}" \
    --protocol "$PROTOCOL" \
    --output "${RESULTS_BASE}" \
    --num_designs "$NUM_DESIGNS"

echo ""
echo "================================================"
echo "  Minibinder design complete"
echo "  Results in: ${RESULTS_BASE}/"
echo "================================================"
