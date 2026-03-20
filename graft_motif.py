#!/usr/bin/env python3
"""Generate BoltzGen design YAML for motif grafting onto nanobody CDRs.

Takes an input structure containing a target protein and a motif (already docked),
and generates a BoltzGen design YAML that grafts the motif into a nanobody CDR
(default: CDR3) with designed flanking regions. All other CDRs are designed
freely to develop additional complementarity against the target.

Usage:
    python graft_motif.py \\
        --input complex.cif \\
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


def parse_scaffold_yaml(scaffold_path):
    """Parse a BoltzGen nanobody scaffold YAML to extract CDR info.

    Returns a dict with:
        - cif_path: absolute path to the scaffold CIF file
        - chain_id: chain identifier in the CIF
        - cdrs: list of 3 (start, end) tuples for CDR1, CDR2, CDR3
        - excludes: list of lists of (start, end) tuples per CDR
        - insertions: list of dicts with id, res_index, num_residues per CDR
        - extra_vis0: any extra residue ranges with visibility 0 (beyond CDRs)
    """
    scaffold_path = Path(scaffold_path)
    with open(scaffold_path) as f:
        data = yaml.safe_load(f)

    chain_id = data["include"][0]["chain"]["id"]
    cif_path = (scaffold_path.parent / data["path"]).resolve()

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
        "cif_path": cif_path,
        "chain_id": chain_id,
        "cdrs": cdrs,
        "excludes": excludes_per_cdr,
        "insertions": insertions,
        "extra_vis0": extra_vis0,
    }


def detect_chains(cif_path):
    """Auto-detect target and motif chains. Larger chain = target, smaller = motif."""
    try:
        import gemmi

        st = gemmi.read_structure(str(cif_path))
        model = st[0]
        chains = []
        for chain in model:
            polymer = chain.get_polymer()
            n_res = sum(1 for _ in polymer)
            chains.append((chain.name, n_res))
    except ImportError:
        # Fallback: parse ATOM records manually
        chains = {}
        with open(cif_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    # PDB format
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
    motif = chains[1][0]
    return target, motif


def make_relative(path, relative_to):
    """Make a path relative to a directory."""
    try:
        return os.path.relpath(path, relative_to)
    except ValueError:
        return str(path)


def generate_graft_yaml(
    input_cif,
    target_chain,
    motif_chain,
    scaffold_info,
    output_dir,
    cdr_target=3,
    flank_n="2..5",
    flank_c="2..5",
):
    """Generate a BoltzGen design YAML for motif grafting.

    The generated YAML constructs the nanobody chain by fusing:
    1. Framework pre-graft CDR (with other CDRs designed normally)
    2. N-terminal designed flank
    3. Grafted motif (fixed structure, anchored to target)
    4. C-terminal designed flank
    5. Framework post-graft CDR

    The motif has visibility 1 (same coordinate frame as target), anchoring the
    nanobody's position. The framework has visibility 2 (internal structure known,
    global position determined by the model). All non-motif CDR residues are
    designed freely to develop complementarity against the target.
    """
    cdr_idx = cdr_target - 1
    cdr_start, cdr_end = scaffold_info["cdrs"][cdr_idx]
    chain_id = scaffold_info["chain_id"]

    # Resolve paths relative to output directory
    input_rel = make_relative(Path(input_cif).resolve(), output_dir)
    scaffold_cif_rel = make_relative(scaffold_info["cif_path"], output_dir)

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

    for orig_idx, _ in pre_cdrs:
        for exc in scaffold_info["excludes"][orig_idx]:
            pre_exclude_entries.append(
                {"chain": {"id": chain_id, "res_index": format_ranges([exc])}}
            )
        if orig_idx < len(scaffold_info["insertions"]):
            ins = scaffold_info["insertions"][orig_idx]
            pre_insertions.append(
                {
                    "insertion": {
                        "id": ins["id"],
                        "res_index": ins["res_index"],
                        "num_residues": ins["num_residues"],
                    }
                }
            )

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
            "path": scaffold_cif_rel,
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

    # --- Entity 4: Grafted motif (fixed structure) ---
    motif_entity = {
        "file": {
            "path": input_rel,
            "fuse": chain_id,
            "include": [{"chain": {"id": motif_chain}}],
        }
    }

    # --- Entity 5: C-terminal CDR flank (designed) ---
    c_flank = {"protein": {"id": "CF", "fuse": chain_id, "sequence": flank_c}}

    # --- Entity 6: Nanobody framework post-graft CDR ---
    post_design_ranges = [cdr for _, cdr in post_cdrs]
    post_exclude_entries = []
    post_insertions = []
    post_vis0_ranges = list(post_design_ranges)

    for orig_idx, _ in post_cdrs:
        for exc in scaffold_info["excludes"][orig_idx]:
            post_exclude_entries.append(
                {"chain": {"id": chain_id, "res_index": format_ranges([exc])}}
            )
        if orig_idx < len(scaffold_info["insertions"]):
            ins = scaffold_info["insertions"][orig_idx]
            post_insertions.append(
                {
                    "insertion": {
                        "id": ins["id"],
                        "res_index": ins["res_index"],
                        "num_residues": ins["num_residues"],
                    }
                }
            )

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
            "path": scaffold_cif_rel,
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
    design = {
        "entities": [
            target_entity,
            pre_framework,
            n_flank,
            motif_entity,
            c_flank,
            post_framework,
        ]
    }

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
        --input complex.cif \\
        --scaffold example/nanobody_scaffolds/7eow.yaml \\
        --output designs/my_graft.yaml

    # Then run:
    boltzgen run designs/my_graft.yaml --protocol nanobody-anything --output results/
""",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input CIF/PDB with target + motif chains (already docked)",
    )
    parser.add_argument(
        "--target-chain",
        default=None,
        help="Target chain ID (default: auto-detect as largest chain)",
    )
    parser.add_argument(
        "--motif-chain",
        default=None,
        help="Motif chain ID (default: auto-detect as second largest chain)",
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
        motif_chain = args.motif_chain
    else:
        detected_target, detected_motif = detect_chains(args.input)
        target_chain = args.target_chain or detected_target
        motif_chain = args.motif_chain or detected_motif
        print(
            f"Auto-detected chains: target={target_chain}, motif={motif_chain}",
            file=sys.stderr,
        )

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
            input_cif=args.input,
            target_chain=target_chain,
            motif_chain=motif_chain,
            scaffold_info=scaffold_infos[0],
            output_dir=output_path.parent,
            cdr_target=args.cdr,
            flank_n=args.flank_n,
            flank_c=args.flank_c,
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
                input_cif=args.input,
                target_chain=target_chain,
                motif_chain=motif_chain,
                scaffold_info=scaffold_info,
                output_dir=output_path.parent,
                cdr_target=args.cdr,
                flank_n=args.flank_n,
                flank_c=args.flank_c,
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
