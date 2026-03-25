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
BINDER_MIN=50               # min binder length
BINDER_MAX=80               # max binder length
LINKER_LENGTH="3..12"       # linkers between hotspot fragments
MOTIF_NOISE=0.3             # σ (Å) for template coordinate perturbation
TARGET_RESIDUES="53,57,60,61,75"  # PHE53, PHE75, SER57, LYS61, LEU60

RESULTS_BASE="results/2lgv_minibinder"
DESIGNS_DIR="designs"
YAML_OUT="${DESIGNS_DIR}/2lgv_minibinder.yaml"

echo "================================================"
echo "  2LGV Minibinder Design (de novo)"
echo "================================================"
echo "Input:           $INPUT_PDB"
echo "Binder size:     ${BINDER_MIN}..${BINDER_MAX} residues"
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
        {'protein': {'id': 'B', 'sequence': '10..25'}},
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
        {'protein': {'id': 'CT', 'fuse': 'B', 'sequence': '10..25'}},
        # Target protein
        {'file': {
            'path': input_rel,
            'include': [{'chain': {'id': 'A'}}],
            'binding_types': [{'chain': {'id': 'A', 'binding': '${TARGET_RESIDUES}'}}],
        }},
    ],
    'constraints': [
        {'total_len': {'id': 'B', 'min': ${BINDER_MIN}, 'max': ${BINDER_MAX}}},
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
