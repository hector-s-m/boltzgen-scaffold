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


def _build_motif_entities(input_pdb, input_rel, motif_chains, chain_id, linker_length, cdr_num):
    """Build motif fragment + linker entities for one graft CDR.

    Returns list of entities and a description string for logging.
    """
    all_fragments = []
    for mc in motif_chains:
        chain_frags = detect_fragments(input_pdb, mc)
        for frag in chain_frags:
            all_fragments.append((mc, frag))

    if len(all_fragments) == 1:
        mc, frag = all_fragments[0]
        desc = f"contiguous ({len(frag)} residues on chain {mc})"
        return [
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
        ], desc

    entities = []
    total_res = sum(len(f) for _, f in all_fragments)
    desc = (
        f"{len(all_fragments)} fragments across chain(s) "
        f"{','.join(motif_chains)} ({total_res} residues, linkers: {linker_length})"
    )
    for i, (mc, fragment) in enumerate(all_fragments):
        chain_frags_for_mc = [f for c, f in all_fragments if c == mc]
        is_whole_chain = len(chain_frags_for_mc) == 1
        if is_whole_chain:
            include_entry = {"chain": {"id": mc}}
        else:
            res_index = format_range(fragment[0], fragment[-1])
            include_entry = {"chain": {"id": mc, "res_index": res_index}}
        entities.append(
            {
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
        )
        if i < len(all_fragments) - 1:
            entities.append(
                {
                    "protein": {
                        "id": f"LK{cdr_num}_{i + 1}",
                        "fuse": chain_id,
                        "sequence": linker_length,
                    }
                }
            )
    return entities, desc


def _build_framework_section(
    scaffold_pdb_rel, chain_id, res_range_str, non_graft_cdrs,
    scaffold_info, cdr_lengths, extra_vis0, is_first,
):
    """Build a framework section entity.

    Parameters
    ----------
    res_range_str : str
        Residue range string for this section (e.g., '..25', '35..97', '119..').
    non_graft_cdrs : list of (orig_idx, cdr_range)
        CDRs that fall within this section and should be designed.
    is_first : bool
        If True, includes reset_res_index and no fuse.
    """
    design_ranges = [cdr for _, cdr in non_graft_cdrs]
    exclude_entries = []
    insertions = []

    for orig_idx, cdr_range in non_graft_cdrs:
        cdr_num = orig_idx + 1
        length_override = cdr_lengths.get(cdr_num, DEFAULT_CDR_LENGTHS.get(cdr_num))
        scaffold_ins = (
            scaffold_info["insertions"][orig_idx]
            if orig_idx < len(scaffold_info["insertions"])
            else None
        )
        excl, ins = build_cdr_design(
            chain_id, cdr_range, scaffold_info["excludes"][orig_idx],
            scaffold_ins, length_override,
        )
        exclude_entries.extend(excl)
        insertions.extend(ins)

    vis0_ranges = list(design_ranges) + list(extra_vis0)
    structure_groups = [{"group": {"id": chain_id, "visibility": 2}}]
    if vis0_ranges:
        structure_groups.append(
            {
                "group": {
                    "id": chain_id,
                    "visibility": 0,
                    "res_index": format_ranges(vis0_ranges),
                }
            }
        )

    section = {
        "file": {
            "path": scaffold_pdb_rel,
            "include": [
                {"chain": {"id": chain_id, "res_index": res_range_str}}
            ],
            "structure_groups": structure_groups,
        }
    }
    if is_first:
        section["file"]["reset_res_index"] = [{"chain": {"id": chain_id}}]
    else:
        section["file"]["fuse"] = chain_id

    if design_ranges:
        section["file"]["design"] = [
            {"chain": {"id": chain_id, "res_index": format_ranges(design_ranges)}}
        ]
    if exclude_entries:
        section["file"]["exclude"] = exclude_entries
    if insertions:
        section["file"]["design_insertions"] = insertions

    return section


def generate_graft_yaml(
    input_pdb,
    target_chain,
    scaffold_info,
    output_dir,
    cdr_motif_map,
    flank_n="2..5",
    flank_c="2..5",
    cdr_lengths=None,
    linker_length="1..6",
    motif_noise=0.0,
    target_residues=None,
):
    """Generate a BoltzGen design YAML for motif grafting.

    Supports grafting motif residues into one or multiple CDRs simultaneously.
    The generated YAML constructs the nanobody chain by fusing framework sections
    with graft blocks (flank + motif fragments + flank) at each graft CDR.

    Parameters
    ----------
    cdr_motif_map : dict[int, list[str]]
        Mapping of CDR numbers to motif chain IDs.
        Single CDR: ``{3: ["B", "C"]}`` — both chains grafted into CDR3.
        Multi CDR: ``{2: ["B"], 3: ["C"]}`` — chain B → CDR2, chain C → CDR3.
    cdr_lengths : dict, optional
        Override CDR lengths for non-graft CDRs.
    linker_length : str
        Length range for designed linkers between hotspot fragments (default: '1..6').
    motif_noise : float
        Gaussian noise σ (Å) on motif side-chain atoms for soft flexibility.
    target_residues : str, optional
        Comma-separated residue numbers on the target chain to mark as binding
        site (e.g. '53,57,60,61,75'). Tells BoltzGen where the nanobody should
        make contacts.
    """
    if cdr_lengths is None:
        cdr_lengths = {}

    chain_id = scaffold_info["chain_id"]
    cdrs = scaffold_info["cdrs"]
    input_rel = make_relative(Path(input_pdb).resolve(), output_dir)
    scaffold_pdb_rel = make_relative(scaffold_info["pdb_path"], output_dir)

    if chain_id == target_chain:
        print(
            f"Warning: scaffold chain '{chain_id}' collides with target chain "
            f"'{target_chain}'. BoltzGen will auto-rename, but verify with "
            f"'boltzgen check'.",
            file=sys.stderr,
        )

    # Sort graft CDRs by position
    graft_cdr_nums = sorted(cdr_motif_map.keys())
    graft_cdr_indices = [n - 1 for n in graft_cdr_nums]  # 0-indexed

    # Determine non-graft CDRs and which framework section they belong to.
    # Sections are numbered 0..N where N = len(graft_cdr_nums).
    # Section boundaries are defined by graft CDR start/end positions.
    all_cdr_indices = list(range(len(cdrs)))
    non_graft_indices = [i for i in all_cdr_indices if i not in graft_cdr_indices]

    # Compute section boundaries: list of (section_start, section_end) in scaffold residue numbers
    # Section 0: start of chain .. first graft CDR start - 1
    # Section k (mid): prev graft CDR end + 1 .. next graft CDR start - 1
    # Section N: last graft CDR end + 1 .. end of chain
    section_boundaries = []
    for s in range(len(graft_cdr_nums) + 1):
        if s == 0:
            sec_end = cdrs[graft_cdr_indices[0]][0] - 1
            section_boundaries.append((None, sec_end))  # None = start of chain
        elif s == len(graft_cdr_nums):
            sec_start = cdrs[graft_cdr_indices[-1]][1] + 1
            section_boundaries.append((sec_start, None))  # None = end of chain
        else:
            sec_start = cdrs[graft_cdr_indices[s - 1]][1] + 1
            sec_end = cdrs[graft_cdr_indices[s]][0] - 1
            section_boundaries.append((sec_start, sec_end))

    # Assign non-graft CDRs to sections based on position
    section_cdrs = [[] for _ in range(len(graft_cdr_nums) + 1)]
    for idx in non_graft_indices:
        cdr_start, cdr_end = cdrs[idx]
        for s, (sec_s, sec_e) in enumerate(section_boundaries):
            s_start = sec_s if sec_s is not None else 0
            s_end = sec_e if sec_e is not None else 99999
            if cdr_start >= s_start and cdr_end <= s_end:
                section_cdrs[s].append((idx, cdrs[idx]))
                break

    # Assign extra_vis0 ranges to sections
    section_vis0 = [[] for _ in range(len(graft_cdr_nums) + 1)]
    for vr in scaffold_info["extra_vis0"]:
        for s, (sec_s, sec_e) in enumerate(section_boundaries):
            s_start = sec_s if sec_s is not None else 0
            s_end = sec_e if sec_e is not None else 99999
            if vr[0] >= s_start and vr[1] <= s_end:
                section_vis0[s].append(vr)
                break

    # --- Build entities ---
    # Target protein
    target_entity = {
        "file": {
            "path": input_rel,
            "include": [{"chain": {"id": target_chain}}],
        }
    }
    if target_residues:
        target_entity["file"]["binding_types"] = [
            {"chain": {"id": target_chain, "binding": target_residues}},
        ]
    entities = [target_entity]

    # Interleave: framework section → [flank_N, motif, flank_C] → next section
    for s in range(len(graft_cdr_nums) + 1):
        sec_s, sec_e = section_boundaries[s]

        # Skip empty mid-sections (adjacent graft CDRs with no residues between)
        if sec_s is not None and sec_e is not None and sec_s > sec_e:
            pass  # no framework residues in this section
        else:
            if sec_s is None:
                res_range = f"..{sec_e}"
            elif sec_e is None:
                res_range = f"{sec_s}.."
            else:
                res_range = f"{sec_s}..{sec_e}"

            section = _build_framework_section(
                scaffold_pdb_rel, chain_id, res_range,
                section_cdrs[s], scaffold_info, cdr_lengths,
                section_vis0[s], is_first=(s == 0),
            )
            entities.append(section)

        # After each section (except the last), insert the graft block
        if s < len(graft_cdr_nums):
            cdr_num = graft_cdr_nums[s]
            mc_list = cdr_motif_map[cdr_num]

            # N-terminal flank
            entities.append(
                {"protein": {"id": f"NF{cdr_num}", "fuse": chain_id, "sequence": flank_n}}
            )

            # Motif fragments + linkers
            motif_ents, desc = _build_motif_entities(
                input_pdb, input_rel, mc_list, chain_id, linker_length, cdr_num,
            )
            print(f"  CDR{cdr_num} motif: {desc}", file=sys.stderr)
            entities.extend(motif_ents)

            # C-terminal flank
            entities.append(
                {"protein": {"id": f"CF{cdr_num}", "fuse": chain_id, "sequence": flank_c}}
            )

    design = {"entities": entities}
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
        help="Which CDR to graft the motif into (default: 3). "
             "Ignored when --cdr-motif is used.",
    )
    parser.add_argument(
        "--cdr-motif",
        nargs="+",
        metavar="CDR:CHAINS",
        default=None,
        help="Map CDR numbers to motif chains for multi-CDR grafting. "
             "E.g. --cdr-motif 2:B 3:C (CDR2 gets chain B, CDR3 gets chain C). "
             "Overrides --cdr and --motif-chain.",
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
        "--target-residues",
        default=None,
        help="Comma-separated residue numbers on the target chain to mark as binding site "
             "(e.g. '53,57,60,61,75'). Guides the nanobody to interact with these positions.",
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

    # Build cdr_motif_map: {cdr_num: [chain_ids]}
    if args.cdr_motif:
        # Multi-CDR mode: parse --cdr-motif specs
        cdr_motif_map = {}
        for spec in args.cdr_motif:
            cdr_str, chains_str = spec.split(":", 1)
            cdr_num = int(cdr_str)
            cdr_motif_map[cdr_num] = [c.strip() for c in chains_str.split(",")]
        # Auto-detect target chain
        detected_target, _ = detect_chains(args.input)
        target_chain = args.target_chain or detected_target
        print(
            "Multi-CDR grafting: " + ", ".join(
                f"CDR{k}={','.join(v)}" for k, v in sorted(cdr_motif_map.items())
            ),
            file=sys.stderr,
        )
    else:
        # Single-CDR mode (backward-compatible)
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
        cdr_motif_map = {args.cdr: motif_chains}

    # Build CDR length overrides (only for non-graft CDRs)
    graft_cdrs = set(cdr_motif_map.keys())
    cdr_lengths = {}
    cdr_length_args = {1: args.cdr1_length, 2: args.cdr2_length, 3: args.cdr3_length}
    for cdr_num, length_spec in cdr_length_args.items():
        if cdr_num in graft_cdrs:
            if length_spec is not None:
                print(
                    f"Warning: --cdr{cdr_num}-length is ignored because CDR{cdr_num} "
                    f"is a graft target (replaced by motif + flanks)",
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

    common_kwargs = dict(
        input_pdb=args.input,
        target_chain=target_chain,
        cdr_motif_map=cdr_motif_map,
        flank_n=args.flank_n,
        flank_c=args.flank_c,
        cdr_lengths=cdr_lengths,
        linker_length=args.linker_length,
        motif_noise=args.motif_noise,
        target_residues=args.target_residues,
    )

    if len(scaffold_infos) == 1:
        design = generate_graft_yaml(
            scaffold_info=scaffold_infos[0],
            output_dir=output_path.parent,
            **common_kwargs,
        )
        with open(output_path, "w") as f:
            yaml.dump(design, f, default_flow_style=False, sort_keys=False)
        print(f"Generated: {output_path}", file=sys.stderr)
    else:
        stem = output_path.stem
        suffix = output_path.suffix
        for i, scaffold_info in enumerate(scaffold_infos):
            scaffold_name = Path(args.scaffold[i]).stem
            out_file = output_path.parent / f"{stem}_{scaffold_name}{suffix}"
            design = generate_graft_yaml(
                scaffold_info=scaffold_info,
                output_dir=output_path.parent,
                **common_kwargs,
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
