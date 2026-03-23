#!/usr/bin/env python3
"""Generate BoltzGen design YAML for motif grafting onto nanobody CDRs.

Takes an input structure containing a target protein and a motif (already docked),
and generates a BoltzGen design YAML that grafts the motif into a nanobody CDR
(default: CDR3) with designed flanking regions. All other CDRs are designed
freely to develop additional complementarity against the target.

Usage:
    python graft_motif.py \\
        --input complex.pdb \\
        --scaffold example/nanobody_scaffolds/7eow.yaml \\
        --output designs/my_graft.yaml

    # Then run with BoltzGen:
    boltzgen run designs/my_graft.yaml --protocol nanobody-anything --output results/
"""

import argparse
import os
import re
import sys
from pathlib import Path

import yaml

# Natural CDR length ranges for nanobodies (approximate, IMGT-like)
DEFAULT_CDR_LENGTHS = {
    1: "5..13",
    2: "3..10",
    3: "10..25",
}


def parse_res_ranges(range_str):
    """Parse a comma-separated residue range string into list of (start, end) tuples.

    Examples:
        '26..34,52..59,98..118' -> [(26, 34), (52, 59), (98, 118)]
        '26..28' -> [(26, 28)]
        '5' -> [(5, 5)]
    """
    ranges = []
    for part in str(range_str).split(","):
        part = part.strip()
        if ".." in part:
            tokens = part.split("..")
            start = int(tokens[0]) if tokens[0] else None
            end = int(tokens[1]) if tokens[1] else None
            ranges.append((start, end))
        else:
            idx = int(part)
            ranges.append((idx, idx))
    return ranges


def range_overlaps(r1, r2):
    """Check if two (start, end) ranges overlap."""
    return r1[0] <= r2[1] and r2[0] <= r1[1]


def format_range(start, end):
    """Format a (start, end) tuple as a range string."""
    if start == end:
        return str(start)
    return f"{start}..{end}"


def format_ranges(ranges):
    """Format a list of (start, end) tuples as a comma-separated range string."""
    return ",".join(format_range(s, e) for s, e in ranges)


def parse_length_spec(spec):
    """Parse a CDR length specification like '5..13' or '10' into (min, max)."""
    spec = str(spec).strip()
    if ".." in spec:
        parts = spec.split("..")
        return int(parts[0]), int(parts[1])
    else:
        n = int(spec)
        return n, n


def parse_scaffold_yaml(scaffold_path):
    """Parse a BoltzGen nanobody scaffold YAML to extract CDR info.

    Returns a dict with:
        - pdb_path: absolute path to the scaffold PDB file
        - chain_id: chain identifier in the PDB
        - cdrs: list of 3 (start, end) tuples for CDR1, CDR2, CDR3
        - excludes: list of lists of (start, end) tuples per CDR
        - insertions: list of dicts with id, res_index, num_residues per CDR
        - extra_vis0: any extra residue ranges with visibility 0 (beyond CDRs)
    """
    scaffold_path = Path(scaffold_path)
    with open(scaffold_path) as f:
        data = yaml.safe_load(f)

    chain_id = data["include"][0]["chain"]["id"]
    pdb_path = (scaffold_path.parent / data["path"]).resolve()
    if not pdb_path.exists() and pdb_path.suffix == ".cif":
        # Try .pdb fallback if .cif referenced but not found
        pdb_path = pdb_path.with_suffix(".pdb")

    # Parse CDR design ranges (expect 3: CDR1, CDR2, CDR3)
    design_str = data["design"][0]["chain"]["res_index"]
    cdrs = parse_res_ranges(design_str)
    if len(cdrs) != 3:
        print(f"Warning: expected 3 CDR ranges, got {len(cdrs)}", file=sys.stderr)

    # Parse exclude entries - match each to a CDR
    excludes_per_cdr = [[] for _ in cdrs]
    for exc_entry in data.get("exclude", []):
        exc_ranges = parse_res_ranges(exc_entry["chain"]["res_index"])
        for exc_range in exc_ranges:
            matched = False
            for i, cdr in enumerate(cdrs):
                if range_overlaps(exc_range, cdr):
                    excludes_per_cdr[i].append(exc_range)
                    matched = True
                    break
            if not matched:
                print(
                    f"Warning: exclude range {exc_range} doesn't match any CDR",
                    file=sys.stderr,
                )

    # Parse insertions (expect 3, in CDR order)
    insertions = []
    for ins_entry in data.get("design_insertions", []):
        ins = ins_entry["insertion"]
        insertions.append(
            {
                "id": ins["id"],
                "res_index": ins["res_index"],
                "num_residues": ins["num_residues"],
            }
        )

    # Parse extra visibility 0 ranges (beyond CDRs)
    extra_vis0 = []
    for sg in data.get("structure_groups", []):
        group = sg["group"]
        if group.get("visibility") == 0 and "res_index" in group:
            vis0_ranges = parse_res_ranges(group["res_index"])
            for vr in vis0_ranges:
                # Check if this range is NOT one of the CDRs
                is_cdr = any(vr == cdr for cdr in cdrs)
                if not is_cdr:
                    extra_vis0.append(vr)

    return {
        "pdb_path": pdb_path,
        "chain_id": chain_id,
        "cdrs": cdrs,
        "excludes": excludes_per_cdr,
        "insertions": insertions,
        "extra_vis0": extra_vis0,
    }


def detect_chains(pdb_path):
    """Auto-detect target and motif chains.

    The largest chain is treated as the target. ALL remaining chains are
    treated as motif/hotspot chains (supports multi-chain disconnected motifs
    where each residue may be on its own chain, e.g. PHE on chain B, TYR on
    chain C).

    Returns
    -------
    target : str
        Chain ID of the target.
    motif_chains : list[str]
        Chain IDs of all motif chains (one or more).
    """
    chains = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain_id = line[21]
                if chain_id not in chains:
                    chains[chain_id] = set()
                try:
                    res_num = int(line[22:26].strip())
                    chains[chain_id].add(res_num)
                except ValueError:
                    pass
    chains = [(cid, len(residues)) for cid, residues in chains.items()]

    if len(chains) < 2:
        print(
            "Error: input structure must have at least 2 chains (target + motif)",
            file=sys.stderr,
        )
        sys.exit(1)

    chains.sort(key=lambda x: x[1], reverse=True)
    target = chains[0][0]
    motif_chains = [cid for cid, _ in chains[1:]]
    return target, motif_chains


def detect_fragments(pdb_path, chain_id, bond_cutoff=2.0):
    """Detect contiguous fragments in a PDB chain by peptide bond distance.

    Groups residues into fragments where consecutive residues are connected
    by a peptide bond (C of residue i to N of residue i+1 within cutoff).
    Disconnected residues start new fragments.

    Parameters
    ----------
    pdb_path : str or Path
        Path to the PDB file.
    chain_id : str
        Chain identifier to analyze.
    bond_cutoff : float
        Maximum C→N distance (Å) to consider residues connected. Default 2.0.

    Returns
    -------
    list of list of int
        Each inner list is a contiguous fragment of residue numbers.
        E.g., [[88, 89], [91]] means residues 88-89 are connected, 91 is isolated.
    """
    import math as _math

    # Parse C and N backbone atom positions per residue
    residues = {}  # res_num -> {"N": (x,y,z), "C": (x,y,z)}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain_id:
                continue
            atom_name = line[12:16].strip()
            if atom_name not in ("C", "N"):
                continue
            try:
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            if res_num not in residues:
                residues[res_num] = {}
            residues[res_num][atom_name] = (x, y, z)

    if not residues:
        return []

    sorted_res = sorted(residues.keys())
    fragments = [[sorted_res[0]]]

    for i in range(1, len(sorted_res)):
        prev_res = sorted_res[i - 1]
        curr_res = sorted_res[i]

        # Check if C of prev and N of curr are within bond distance
        connected = False
        if "C" in residues.get(prev_res, {}) and "N" in residues.get(curr_res, {}):
            c_xyz = residues[prev_res]["C"]
            n_xyz = residues[curr_res]["N"]
            dist = _math.sqrt(sum((a - b) ** 2 for a, b in zip(c_xyz, n_xyz)))
            connected = dist <= bond_cutoff

        if connected:
            fragments[-1].append(curr_res)
        else:
            fragments.append([curr_res])

    return fragments


def make_relative(path, relative_to):
    """Make a path relative to a directory."""
    try:
        return os.path.relpath(path, relative_to)
    except ValueError:
        return str(path)


def build_cdr_design(chain_id, cdr_range, scaffold_excludes, scaffold_insertion, cdr_length_override):
    """Build exclude/insertion entries for a single non-graft CDR.

    If cdr_length_override is given (e.g. '5..13'), the entire CDR is excluded
    and replaced with the specified number of designed residues.
    Otherwise, uses the scaffold's own exclude/insertion settings.
    """
    cdr_start, cdr_end = cdr_range

    if cdr_length_override is not None:
        # Override: exclude all scaffold CDR residues, insert desired length
        min_len, max_len = parse_length_spec(cdr_length_override)
        exclude_entries = [
            {"chain": {"id": chain_id, "res_index": format_range(cdr_start, cdr_end)}}
        ]
        insertion_entries = [
            {
                "insertion": {
                    "id": chain_id,
                    "res_index": cdr_start,
                    "num_residues": f"{min_len}..{max_len}" if min_len != max_len else min_len,
                }
            }
        ]
        return exclude_entries, insertion_entries

    # Use scaffold defaults
    exclude_entries = []
    for exc in scaffold_excludes:
        exclude_entries.append(
            {"chain": {"id": chain_id, "res_index": format_ranges([exc])}}
        )
    insertion_entries = []
    if scaffold_insertion is not None:
        insertion_entries.append(
            {
                "insertion": {
                    "id": scaffold_insertion["id"],
                    "res_index": scaffold_insertion["res_index"],
                    "num_residues": scaffold_insertion["num_residues"],
                }
            }
        )
    return exclude_entries, insertion_entries


def generate_graft_yaml(
    input_pdb,
    target_chain,
    motif_chains,
    scaffold_info,
    output_dir,
    cdr_target=3,
    flank_n="2..5",
    flank_c="2..5",
    cdr_lengths=None,
    linker_length="1..6",
    motif_noise=0.0,
):
    """Generate a BoltzGen design YAML for motif grafting.

    The generated YAML constructs the nanobody chain by fusing:
    1. Framework pre-graft CDR (with other CDRs designed normally)
    2. N-terminal designed flank
    3. Grafted motif (fixed structure, anchored to target) — or alternating
       hotspot fragments + designed linkers if the motif is disconnected
    4. C-terminal designed flank
    5. Framework post-graft CDR

    The motif has visibility 1 (same coordinate frame as target), anchoring the
    nanobody's position. The framework has visibility 2 (internal structure known,
    global position determined by the model). All non-motif CDR residues are
    designed freely to develop complementarity against the target.

    Motif residues can span multiple chains (e.g. PHE on chain B, TYR on chain C).
    Each chain is treated as a separate fragment connected by designed linkers.
    Within a single chain, disconnected fragments (no peptide bond) are also
    detected and linked automatically.

    Parameters
    ----------
    motif_chains : list[str]
        One or more chain IDs containing motif/hotspot residues.
    cdr_lengths : dict, optional
        Override CDR lengths for non-graft CDRs. Keys are CDR numbers (1, 2, 3),
        values are length specs like '5..13' or '10'. If None, uses natural defaults.
    linker_length : str
        Length range for designed linkers between hotspot fragments (default: '1..6').
        Only used in hotspot mode (disconnected motif).
    """
    if cdr_lengths is None:
        cdr_lengths = {}

    cdr_idx = cdr_target - 1
    cdr_start, cdr_end = scaffold_info["cdrs"][cdr_idx]
    chain_id = scaffold_info["chain_id"]

    # Check for chain ID collision between scaffold and target/motif
    if chain_id == target_chain:
        print(
            f"Warning: scaffold chain '{chain_id}' collides with target chain '{target_chain}'. "
            f"BoltzGen will auto-rename, but verify the output with 'boltzgen check'.",
            file=sys.stderr,
        )

    # Resolve paths relative to output directory
    input_rel = make_relative(Path(input_pdb).resolve(), output_dir)
    scaffold_pdb_rel = make_relative(scaffold_info["pdb_path"], output_dir)

    # Split CDRs into pre-graft and post-graft groups
    pre_cdrs = [(i, scaffold_info["cdrs"][i]) for i in range(len(scaffold_info["cdrs"])) if i < cdr_idx]
    post_cdrs = [(i, scaffold_info["cdrs"][i]) for i in range(len(scaffold_info["cdrs"])) if i > cdr_idx]

    # --- Entity 1: Target protein ---
    target_entity = {
        "file": {
            "path": input_rel,
            "include": [{"chain": {"id": target_chain}}],
        }
    }

    # --- Entity 2: Nanobody framework pre-graft CDR ---
    pre_design_ranges = [cdr for _, cdr in pre_cdrs]
    pre_exclude_entries = []
    pre_insertions = []
    pre_vis0_ranges = list(pre_design_ranges)

    for orig_idx, cdr_range in pre_cdrs:
        cdr_num = orig_idx + 1
        length_override = cdr_lengths.get(cdr_num, DEFAULT_CDR_LENGTHS[cdr_num])
        scaffold_ins = scaffold_info["insertions"][orig_idx] if orig_idx < len(scaffold_info["insertions"]) else None
        excl, ins = build_cdr_design(
            chain_id, cdr_range, scaffold_info["excludes"][orig_idx],
            scaffold_ins, length_override,
        )
        pre_exclude_entries.extend(excl)
        pre_insertions.extend(ins)

    # Include extra visibility 0 ranges that fall before the graft CDR
    for vr in scaffold_info["extra_vis0"]:
        if vr[1] < cdr_start:
            pre_vis0_ranges.append(vr)

    pre_structure_groups = [{"group": {"id": chain_id, "visibility": 2}}]
    if pre_vis0_ranges:
        pre_structure_groups.append(
            {
                "group": {
                    "id": chain_id,
                    "visibility": 0,
                    "res_index": format_ranges(pre_vis0_ranges),
                }
            }
        )

    pre_framework = {
        "file": {
            "path": scaffold_pdb_rel,
            "include": [
                {"chain": {"id": chain_id, "res_index": f"..{cdr_start - 1}"}}
            ],
            "structure_groups": pre_structure_groups,
            "reset_res_index": [{"chain": {"id": chain_id}}],
        }
    }
    if pre_design_ranges:
        pre_framework["file"]["design"] = [
            {
                "chain": {
                    "id": chain_id,
                    "res_index": format_ranges(pre_design_ranges),
                }
            }
        ]
    if pre_exclude_entries:
        pre_framework["file"]["exclude"] = pre_exclude_entries
    if pre_insertions:
        pre_framework["file"]["design_insertions"] = pre_insertions

    # --- Entity 3: N-terminal CDR flank (designed) ---
    n_flank = {"protein": {"id": "NF", "fuse": chain_id, "sequence": flank_n}}

    # --- Entity 4: Grafted motif / hotspot fragments ---
    # Collect all fragments across all motif chains.  Each chain may contain
    # one or more contiguous fragments; residues on different chains are always
    # separate fragments.
    all_fragments = []  # list of (chain_id, [res_nums])
    for mc in motif_chains:
        chain_frags = detect_fragments(input_pdb, mc)
        for frag in chain_frags:
            all_fragments.append((mc, frag))

    if len(all_fragments) == 1:
        # Single contiguous motif on one chain
        mc, frag = all_fragments[0]
        motif_entities = [
            {
                "file": {
                    "path": input_rel,
                    "fuse": chain_id,
                    "include": [{"chain": {"id": mc}}],
                    "not_design": [{"chain": {"id": mc}}],
                    "structure_groups": [
                        {"group": {"id": mc, "visibility": 1}},
                    ],
                }
            }
        ]
        print(
            f"  Motif: contiguous ({len(frag)} residues on chain {mc})",
            file=sys.stderr,
        )
    else:
        # Hotspot mode: alternating fragment + designed linker entities
        motif_entities = []
        total_hotspot_res = sum(len(f) for _, f in all_fragments)
        print(
            f"  Motif: {len(all_fragments)} disconnected fragments across "
            f"chain(s) {','.join(motif_chains)} "
            f"({total_hotspot_res} hotspot residues, linkers: {linker_length})",
            file=sys.stderr,
        )
        # Count total residues per motif chain to detect whole-chain fragments
        chain_res_counts = {}
        for mc, frag in all_fragments:
            chain_res_counts.setdefault(mc, 0)
            chain_res_counts[mc] += len(frag)

        for i, (mc, fragment) in enumerate(all_fragments):
            # Fragment entity: anchored, not designed
            # If the fragment spans the entire chain, include the whole chain
            # (BoltzGen's res_index uses 1-indexed positions within the chain,
            # not PDB residue numbers, so large PDB numbers would be out of bounds)
            chain_frags_for_mc = [f for c, f in all_fragments if c == mc]
            is_whole_chain = len(chain_frags_for_mc) == 1
            if is_whole_chain:
                include_entry = {"chain": {"id": mc}}
            else:
                res_index = format_range(fragment[0], fragment[-1])
                include_entry = {"chain": {"id": mc, "res_index": res_index}}
            frag_entity = {
                "file": {
                    "path": input_rel,
                    "fuse": chain_id,
                    "include": [include_entry],
                    "not_design": [{"chain": {"id": mc}}],
                    "structure_groups": [
                        {"group": {"id": mc, "visibility": 1}},
                    ],
                }
            }
            motif_entities.append(frag_entity)

            # Designed linker between fragments (not after last)
            if i < len(all_fragments) - 1:
                linker_id = f"LK{i + 1}"
                linker_entity = {
                    "protein": {
                        "id": linker_id,
                        "fuse": chain_id,
                        "sequence": linker_length,
                    }
                }
                motif_entities.append(linker_entity)

    # --- Entity 5: C-terminal CDR flank (designed) ---
    c_flank = {"protein": {"id": "CF", "fuse": chain_id, "sequence": flank_c}}

    # --- Entity 6: Nanobody framework post-graft CDR ---
    post_design_ranges = [cdr for _, cdr in post_cdrs]
    post_exclude_entries = []
    post_insertions = []
    post_vis0_ranges = list(post_design_ranges)

    for orig_idx, cdr_range in post_cdrs:
        cdr_num = orig_idx + 1
        length_override = cdr_lengths.get(cdr_num, DEFAULT_CDR_LENGTHS[cdr_num])
        scaffold_ins = scaffold_info["insertions"][orig_idx] if orig_idx < len(scaffold_info["insertions"]) else None
        excl, ins = build_cdr_design(
            chain_id, cdr_range, scaffold_info["excludes"][orig_idx],
            scaffold_ins, length_override,
        )
        post_exclude_entries.extend(excl)
        post_insertions.extend(ins)

    for vr in scaffold_info["extra_vis0"]:
        if vr[0] > cdr_end:
            post_vis0_ranges.append(vr)

    post_structure_groups = [{"group": {"id": chain_id, "visibility": 2}}]
    if post_vis0_ranges:
        post_structure_groups.append(
            {
                "group": {
                    "id": chain_id,
                    "visibility": 0,
                    "res_index": format_ranges(post_vis0_ranges),
                }
            }
        )

    post_framework = {
        "file": {
            "path": scaffold_pdb_rel,
            "fuse": chain_id,
            "include": [
                {"chain": {"id": chain_id, "res_index": f"{cdr_end + 1}.."}}
            ],
            "structure_groups": post_structure_groups,
        }
    }
    if post_design_ranges:
        post_framework["file"]["design"] = [
            {
                "chain": {
                    "id": chain_id,
                    "res_index": format_ranges(post_design_ranges),
                }
            }
        ]
    if post_exclude_entries:
        post_framework["file"]["exclude"] = post_exclude_entries
    if post_insertions:
        post_framework["file"]["design_insertions"] = post_insertions

    # Assemble the full design YAML
    entities = [target_entity, pre_framework, n_flank]
    entities.extend(motif_entities)
    entities.extend([c_flank, post_framework])
    design = {"entities": entities}

    # Add motif noise if specified (σ in Å for template coordinate perturbation)
    if motif_noise > 0:
        design["motif_noise"] = motif_noise

    return design


def list_scaffolds(scaffolds_dir):
    """List available nanobody scaffold YAML files."""
    scaffolds_dir = Path(scaffolds_dir)
    if not scaffolds_dir.exists():
        return []
    return sorted(scaffolds_dir.glob("*.yaml"))


def main():
    default_scaffolds = Path(__file__).parent / "example" / "nanobody_scaffolds"

    parser = argparse.ArgumentParser(
        description="Generate BoltzGen YAML for motif grafting onto nanobody CDRs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Example:
    python graft_motif.py \\
        --input complex.pdb \\
        --scaffold example/nanobody_scaffolds/7eow.yaml \\
        --output designs/my_graft.yaml

    # Tune CDR lengths:
    python graft_motif.py \\
        --input complex.pdb \\
        --scaffold example/nanobody_scaffolds/7eow.yaml \\
        --output designs/my_graft.yaml \\
        --cdr1-length 6..10 --cdr2-length 4..8

    # Then run:
    boltzgen run designs/my_graft.yaml --protocol nanobody-anything --output results/
""",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input PDB file with target + motif chains (already docked)",
    )
    parser.add_argument(
        "--target-chain",
        default=None,
        help="Target chain ID (default: auto-detect as largest chain)",
    )
    parser.add_argument(
        "--motif-chain",
        default=None,
        help="Motif chain ID(s), comma-separated for multi-chain motifs "
             "(e.g. 'B,C'). Default: auto-detect all non-target chains.",
    )
    parser.add_argument(
        "--scaffold",
        required=True,
        nargs="+",
        help="Path(s) to nanobody scaffold YAML file(s)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output YAML file path",
    )
    parser.add_argument(
        "--cdr",
        type=int,
        default=3,
        choices=[1, 2, 3],
        help="Which CDR to graft the motif into (default: 3)",
    )
    parser.add_argument(
        "--flank-n",
        default="2..5",
        help="N-terminal flank length range (default: 2..5)",
    )
    parser.add_argument(
        "--flank-c",
        default="2..5",
        help="C-terminal flank length range (default: 2..5)",
    )
    parser.add_argument(
        "--cdr1-length",
        default=None,
        help=f"CDR1 length range, e.g. '6..10' (default: {DEFAULT_CDR_LENGTHS[1]})",
    )
    parser.add_argument(
        "--cdr2-length",
        default=None,
        help=f"CDR2 length range, e.g. '4..8' (default: {DEFAULT_CDR_LENGTHS[2]})",
    )
    parser.add_argument(
        "--cdr3-length",
        default=None,
        help=f"CDR3 length range, e.g. '8..20' (default: {DEFAULT_CDR_LENGTHS[3]})",
    )
    parser.add_argument(
        "--linker-length",
        default="1..6",
        help="Linker length range between hotspot fragments (default: 1..6). "
             "Only used when the motif chain has disconnected residues.",
    )
    parser.add_argument(
        "--motif-noise",
        type=float,
        default=0.0,
        help="Gaussian noise σ (Å) added to motif template coordinates for soft rigidity. "
             "0 = rigid (default), 0.3 = slight flex, 1.0 = loose.",
    )
    parser.add_argument(
        "--list-scaffolds",
        action="store_true",
        help="List available nanobody scaffolds and exit",
    )

    args = parser.parse_args()

    if args.list_scaffolds:
        scaffolds = list_scaffolds(default_scaffolds)
        if scaffolds:
            print("Available nanobody scaffolds:")
            for s in scaffolds:
                print(f"  {s}")
        else:
            print(f"No scaffolds found in {default_scaffolds}")
        return

    # Auto-detect chains if not specified
    if args.target_chain and args.motif_chain:
        target_chain = args.target_chain
        motif_chains = [c.strip() for c in args.motif_chain.split(",")]
    else:
        detected_target, detected_motif_chains = detect_chains(args.input)
        target_chain = args.target_chain or detected_target
        if args.motif_chain:
            motif_chains = [c.strip() for c in args.motif_chain.split(",")]
        else:
            motif_chains = detected_motif_chains
        print(
            f"Auto-detected chains: target={target_chain}, motif={','.join(motif_chains)}",
            file=sys.stderr,
        )

    # Build CDR length overrides (only for non-graft CDRs)
    cdr_lengths = {}
    cdr_length_args = {1: args.cdr1_length, 2: args.cdr2_length, 3: args.cdr3_length}
    for cdr_num, length_spec in cdr_length_args.items():
        if cdr_num == args.cdr:
            if length_spec is not None:
                print(
                    f"Warning: --cdr{cdr_num}-length is ignored because CDR{cdr_num} "
                    f"is the graft target (replaced by motif + flanks)",
                    file=sys.stderr,
                )
            continue
        if length_spec is not None:
            cdr_lengths[cdr_num] = length_spec

    # Parse scaffold(s)
    scaffold_infos = []
    for scaffold_path in args.scaffold:
        scaffold_infos.append(parse_scaffold_yaml(scaffold_path))

    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if len(scaffold_infos) == 1:
        # Single scaffold: generate one YAML
        design = generate_graft_yaml(
            input_pdb=args.input,
            target_chain=target_chain,
            motif_chains=motif_chains,
            scaffold_info=scaffold_infos[0],
            output_dir=output_path.parent,
            cdr_target=args.cdr,
            flank_n=args.flank_n,
            flank_c=args.flank_c,
            cdr_lengths=cdr_lengths,
            linker_length=args.linker_length,
            motif_noise=args.motif_noise,
        )
        with open(output_path, "w") as f:
            yaml.dump(design, f, default_flow_style=False, sort_keys=False)
        print(f"Generated: {output_path}", file=sys.stderr)
    else:
        # Multiple scaffolds: generate one YAML per scaffold
        stem = output_path.stem
        suffix = output_path.suffix
        for i, scaffold_info in enumerate(scaffold_infos):
            scaffold_name = Path(args.scaffold[i]).stem
            out_file = output_path.parent / f"{stem}_{scaffold_name}{suffix}"
            design = generate_graft_yaml(
                input_pdb=args.input,
                target_chain=target_chain,
                motif_chains=motif_chains,
                scaffold_info=scaffold_info,
                output_dir=output_path.parent,
                cdr_target=args.cdr,
                flank_n=args.flank_n,
                flank_c=args.flank_c,
                cdr_lengths=cdr_lengths,
            )
            with open(out_file, "w") as f:
                yaml.dump(design, f, default_flow_style=False, sort_keys=False)
            print(f"Generated: {out_file}", file=sys.stderr)

    print("\nRun with:", file=sys.stderr)
    print(
        f"  boltzgen run {output_path} --protocol nanobody-anything --output results/",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
