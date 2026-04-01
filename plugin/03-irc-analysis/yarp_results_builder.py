#!/usr/bin/env python3
"""
Frame-focused IRC/TS processor for text-based geometry files.

For each input directory containing reaction coordinate files, this script:
  1) Verifies required inputs exist: finished_first.xyz, finished_last.xyz, finished_irc.trj
  2) Determines the Transition State (TS) geometry source.
  3) For each frame in each file, it computes the adjacency and bond-electron matrices.
  4) **For finished_irc.trj, it identifies the frame with the highest energy as the TS.**
  5) Writes the results into a structured output folder named 'IRC_Analysis'.
     - Output includes per-frame JSON, labelled CSVs, and a ts_frame.txt for the IRC.
  6) Is idempotent: by default, it will not recompute frames if the output files already exist.
"""
from __future__ import annotations
import argparse
import sys
import json
import tempfile
import csv
from pathlib import Path
from typing import List, Dict, Any, Tuple

import numpy as np

import re

# Your helpers (assumed available in env)
from metal_ligand.yarp_helpers import silent_yarpecule
# --- Configuration ---
REQUIRED_FILES = [
    "finished_first.xyz",
    "finished_last.xyz",
    "finished_irc.trj",
]
OUTPUT_FOLDER_NAME = "IRC_Analysis"

# -------------------------------
# Parsing Utilities
# -------------------------------

def parse_xyz_text_frames(text: str) -> List[Tuple[List[str], np.ndarray, str]]:
    """
    Parses XYZ/TRJ-style text containing one or more frames.
    Returns a list of (elements, coords[n,3], comment_line) per frame.
    """
    lines = text.splitlines()
    i = 0
    frames: List[Tuple[List[str], np.ndarray, str]] = []
    nlines = len(lines)
    while i < nlines:
        # Skip blank lines
        while i < nlines and not lines[i].strip():
            i += 1
        if i >= nlines:
            break

        comment = ""
        # Read atom count from the header
        try:
            num_atoms = int(lines[i].strip())
            i += 1  # Move to the comment line
            # Capture the optional comment line
            if i < nlines:
                comment = lines[i].strip()
                i += 1
        except (ValueError, IndexError):
            # If no valid header, assume the rest of the file is a single frame
            num_atoms = nlines - i
            comment = "" # No comment in fallback mode

        # Read atom coordinate lines
        elems, coords = [], []
        for _ in range(num_atoms):
            if i >= nlines:
                break
            parts = lines[i].split()
            i += 1
            if len(parts) >= 4:
                try:
                    elems.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except (ValueError, IndexError):
                    continue # Skip malformed lines
        
        if elems:
            frames.append((elems, np.asarray(coords, dtype=float), comment))
            
    return frames


def yarpecule_from_frame(elems: List[str], coords: np.ndarray, charge: int | None = None):
    """Creates a yarpecule object from frame data via a temporary XYZ file."""
    with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False, encoding='utf-8') as tf:
        tf.write(f"{len(elems)}\n\n")
        for e, (x, y, z) in zip(elems, coords):
            tf.write(f"{e} {x:.10f} {y:.10f} {z:.10f}\n")
        tf_path = Path(tf.name)
    try:
        if charge is None:
            y = silent_yarpecule(tf_path, canon=False)
        else:
            # requires your wrapper to accept charge= (as shown earlier)
            y = silent_yarpecule(tf_path, charge = int(charge), canon=False)
    finally:
        tf_path.unlink(missing_ok=True)
    return y

# -------------------------------
# File Output Utilities
# -------------------------------

def build_clean_json_payload(source_str: str, labels: List[str], elements: List[str], adj_i: np.ndarray, be_i: np.ndarray) -> Dict[str, Any]:
    """Builds a JSON-serializable dictionary for a molecule's structure."""
    n_atoms = adj_i.shape[0]
    atoms = [{"id": labels[k], "el": str(elements[k]), "e": int(be_i[k, k])} for k in range(n_atoms)]
    
    bonds = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if int(adj_i[i, j]) > 0:
                bonds.append({"i": labels[i], "j": labels[j], "order": int(be_i[i, j])})
                
    return {"source": source_str, "atoms": atoms, "bonds": bonds}

def write_labelled_csv(path: Path, matrix: np.ndarray, labels: List[str]) -> None:
    """Writes a labelled N×N matrix to a CSV file with headers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([""] + labels)
        for i, row in enumerate(matrix.astype(int)):
            writer.writerow([labels[i]] + row.tolist())

# -------------------------------
# Core Processing Logic
# -------------------------------

#def process_directory(source_dir: Path, overwrite: bool = False, verbose: bool = True):
def process_directory(source_dir: Path, overwrite: bool = False, verbose: bool = True, charge: int | None = None):
    """
    Processes a single directory containing reaction coordinate files.
    """
    # 1. --- Validate Inputs ---
    if not source_dir.is_dir():
        print(f"  -> error: Path is not a directory: {source_dir}", file=sys.stderr)
        return False

    missing_files = [f for f in REQUIRED_FILES if not (source_dir / f).exists()]
    if missing_files:
        print(f"  -> error: Missing required files in {source_dir}: {', '.join(missing_files)}", file=sys.stderr)
        return False
        
    # 2. --- Setup Output Directory ---
    output_dir = source_dir / OUTPUT_FOLDER_NAME
    output_dir.mkdir(parents=True, exist_ok=True)

    # 3. --- Determine TS File Sources ---
    ts_sources: List[Tuple[str, Path]] = []  # (label, path_to_file)
    ts_final_geom_path = source_dir / "ts_final_geometry.xyz"
    # NEW: prefer input.xyz in the same folder; fall back to the old rule if missing
    original_ts_path = source_dir / "input.xyz"

    if ts_final_geom_path.exists():
        ts_sources.append(("TS", ts_final_geom_path))
        if original_ts_path.exists():
            ts_sources.append(("initial-TS", original_ts_path))
        elif verbose:
            print(f"  -> [info] Sibling initial-TS file not found: {original_ts_path}")
    elif original_ts_path.exists():
        ts_sources.append(("TS", original_ts_path))
    else:
        print(f"  -> error: Could not find a TS file (neither {ts_final_geom_path.name} nor {original_ts_path.name})", file=sys.stderr)
        return False
        
    # 4. --- Main Processing Loop ---
    all_sources = [
        ("finished_first", source_dir / "finished_first.xyz"),
        ("finished_last", source_dir / "finished_last.xyz"),
        ("finished_irc", source_dir / "finished_irc.trj"),
    ] + ts_sources
    
    total_frames_processed = 0

    for label, file_path in all_sources:
        if verbose:
            print(f"  -> Processing {label} from {file_path.name}")
        
        try:
            text_content = file_path.read_text(encoding='utf-8')
            frames = parse_xyz_text_frames(text_content)
        except Exception as e:
            print(f"  -> error: Failed to read or parse {file_path}: {e}", file=sys.stderr)
            continue

        # --- NEW: IRC Transition State Detection ---
        if label == "finished_irc" and frames:
            max_energy = -float('inf')
            ts_frame_index = -1
            
            for i, (elems, coords, comment) in enumerate(frames):
                try:
                    # Energy is usually the first value in the comment
                    energy = float(comment.split()[0])
                    if energy > max_energy:
                        max_energy = energy
                        ts_frame_index = i
                except (ValueError, IndexError):
                    if verbose:
                        print(f"     [warn] Could not parse energy from comment in frame {i}: '{comment}'")
                    continue
            
            if ts_frame_index != -1:
                ts_info_path = output_dir / "finished_irc" / "ts_frame.txt"
                ts_info_path.parent.mkdir(parents=True, exist_ok=True)
                with ts_info_path.open("w", encoding='utf-8') as f:
                    f.write(f"TS-Frame: {ts_frame_index}, energy: {max_energy}\n")
                if verbose:
                    print(f"     ✨ IRC Transition State identified at frame {ts_frame_index} (Energy: {max_energy})")
        # --- END NEW LOGIC ---

        for i, (elements, coords, comment) in enumerate(frames):
            frame_tag = f"frame_{i:05d}"
            frame_out_dir = output_dir / label
            
            # Define expected output files for idempotency check
            json_path = frame_out_dir / f"{frame_tag}.json"
            adj_csv_path = frame_out_dir / f"{frame_tag}_adjacency.csv"
            be_csv_path = frame_out_dir / f"{frame_tag}_bond_electrons.csv"
            
            # Skip if all output files exist and overwrite is false
            if not overwrite and all(p.exists() for p in [json_path, adj_csv_path, be_csv_path]):
                continue

            # Perform the calculation
            try:
                eff_charge = charge if charge is not None else parse_charge_from_comment(comment)
                y = yarpecule_from_frame(elements, coords, charge = eff_charge)
                elements_full = [e.capitalize() for e in y.elements]
                labels = [f"{sym}{k}" for k, sym in enumerate(elements_full)]
                adj_i = y.adj_mat.astype(np.int32)
                be_i = y.bond_mats[0].astype(np.int32)
            except Exception as e:
                print(f"  -> error: YARP calculation failed for {label}/{frame_tag}: {e}", file=sys.stderr)
                continue

            # Write the output files
            frame_out_dir.mkdir(parents=True, exist_ok=True)
            source_str = f"{label}/{frame_tag}"

            # Write JSON
            payload = build_clean_json_payload(source_str, labels, elements_full, adj_i, be_i)
            with json_path.open("w", encoding='utf-8') as f:
                json.dump(payload, f, separators=(",", ":"))

            # Write CSVs
            write_labelled_csv(adj_csv_path, adj_i, labels)
            write_labelled_csv(be_csv_path, be_i, labels)
            
            total_frames_processed += 1
    # === BEGIN: write reaction/deterministic_{forward,reverse}.clean.json ===
    # We assume finished_first.json and finished_last.json live under output_dir
    reaction_dir = (output_dir / "reaction")
    reaction_dir.mkdir(parents=True, exist_ok=True)

    def _load_graph_from_endpoint_json(p: Path) -> dict:
        """
        Accepts an endpoint JSON you just wrote for finished_first/finished_last and
        returns a dict with keys {"atoms": [...], "bonds": [...]}.
        We support a few common shapes:
          - top-level {"atoms": [...], "bonds": [...]}
          - {"graph": {"atoms": [...], "bonds": [...]}}
        """
        with p.open("r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict) and "atoms" in data and "bonds" in data:
            return {"atoms": data["atoms"], "bonds": data["bonds"]}
        if isinstance(data, dict) and "graph" in data \
           and isinstance(data["graph"], dict) \
           and "atoms" in data["graph"] and "bonds" in data["graph"]:
            return {"atoms": data["graph"]["atoms"], "bonds": data["graph"]["bonds"]}
        raise ValueError(f"Could not find atoms/bonds in {p}")

    def _bkey(i, j):
        # canonical, order-independent pair key (string IDs like "C0", "Rh1", etc.)
        return tuple(sorted((i, j)))

    def _to_pair_map(bonds: list[dict]) -> dict:
        """
        bonds: list of {"i": <id>, "j": <id>, "order": int}
        returns: dict[(i,j)] -> int order; a missing key means 'NaN' (no bond present)
        """
        out = {}
        for b in bonds:
            i, j = b["i"], b["j"]
            o = int(b["order"])
            out[_bkey(i, j)] = o
        return out

    def _collect_reactivity(reactant_bonds: list[dict], product_bonds: list[dict]) -> dict:
        """
        Build the `target` for the reaction clean json:
          - reactive_pairs: only true formations/breaks (one side 'NaN' or 0, the other >0)
          - bond_order_changes: both sides >0 but order changed (e.g., 2 -> 1)
          - reactive_atoms: union of atoms touching the above pairs
        We treat 0 as 'non-bond/dative' for classification per your rule.
        """
        R = _to_pair_map(reactant_bonds)
        P = _to_pair_map(product_bonds)
        all_pairs = set(R.keys()) | set(P.keys())

        def tok(v):  # pretty token for the change string
            return "NaN" if v is None else str(int(v))

        reactive_pairs = []        # [["A","B","NaN -> 1"], ...]
        bond_order_changes = []    # [["A","B","2 -> 1"], ...]
        rattoms = set()

        for pair in sorted(all_pairs):
            i, j = pair
            ro = R.get(pair, None)   # None means 'NaN' (absent)
            po = P.get(pair, None)

            if ro == po:
                continue

            # 'Non-bond' test: treat None or 0 as "no bond" for actual break/form
            r_is_nb = (ro is None) or (ro == 0)
            p_is_nb = (po is None) or (po == 0)

            if r_is_nb != p_is_nb:
                # actual formation or breaking
                reactive_pairs.append([i, j, f"{tok(ro)} -> {tok(po)}"])
                rattoms.update((i, j))
            else:
                # both sides bonds (>0) but order changed
                if (ro is not None and ro > 0) and (po is not None and po > 0):
                    bond_order_changes.append([i, j, f"{tok(ro)} -> {tok(po)}"])
                    rattoms.update((i, j))

        return {
            "reactive_atoms": sorted(rattoms),
            "reactive_pairs": reactive_pairs,
            "bond_order_changes": bond_order_changes,
        }

    # Locate the two endpoint files we just wrote
    ff_path = output_dir / "finished_first" / "frame_00000.json"
    fl_path = output_dir / "finished_last" / "frame_00000.json"
    if ff_path.exists() and fl_path.exists():
        try:
            reactant_graph = _load_graph_from_endpoint_json(ff_path)  # Forward: first → last
            product_graph  = _load_graph_from_endpoint_json(fl_path)

            # Forward file
            forward_payload = {
                "task": {"reaction_class": "Unprocessed"},
                "input_graphs": {
                    "reactant": reactant_graph,
                    "product":  product_graph,
                },
                "target": _collect_reactivity(reactant_graph["bonds"], product_graph["bonds"]),
                "schema_version": 3,
                "direction": "Forward",
            }
            with (reaction_dir / "deterministic_forward.clean.json").open("w", encoding="utf-8") as f:
                json.dump(forward_payload, f, indent=2)

            # Reverse file (swap roles)
            reverse_payload = {
                "task": {"reaction_class": "Unprocessed"},
                "input_graphs": {
                    "reactant": product_graph,
                    "product":  reactant_graph,
                },
                "target": _collect_reactivity(product_graph["bonds"], reactant_graph["bonds"]),
                "schema_version": 3,
                "direction": "Reverse",
            }
            with (reaction_dir / "deterministic_reverse.clean.json").open("w", encoding="utf-8") as f:
                json.dump(reverse_payload, f, indent=2)

            if verbose:
                print(f"  -> wrote reaction/{'deterministic_forward.clean.json'}")
                print(f"  -> wrote reaction/{'deterministic_reverse.clean.json'}")
        except Exception as e:
            print(f"  -> error while writing reaction clean jsons: {e}", file=sys.stderr)
    else:
        if verbose:
            missing = []
            if not ff_path.exists(): missing.append("finished_first.json")
            if not fl_path.exists(): missing.append("finished_last.json")
            print(f"  -> skip reaction clean jsons (missing: {', '.join(missing)})")
    # === END: write reaction/deterministic_{forward,reverse}.clean.json ===
            
    if verbose:
        print(f"  -> ✅ Finished. Processed {total_frames_processed} new frames.")
    return True

# -------------------------------
# Command-Line Interface
# -------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build per-frame analysis (JSON, CSV) from IRC/TS text files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("paths", nargs="+", type=Path,
                        help="One or more directories containing the required .xyz/.trj files.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Force recalculation and overwrite existing output files.")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress informational messages, only show errors.")
    parser.add_argument("--charge", type=int,
                        help="Total molecular charge to use for ALL frames (overrides comment parsing).")
    args = parser.parse_args()

    # Find unique directories to process
    job_dirs = set()
    for p in args.paths:
        resolved_p = p.resolve()
        if resolved_p.is_dir():
            job_dirs.add(resolved_p)
        else:
            print(f"[warn] Skipping non-directory path: {p}", file=sys.stderr)
    
    if not job_dirs:
        print("No valid input directories found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(job_dirs)} unique directories to process...")
    success_count, error_count = 0, 0

    for i, job_dir in enumerate(sorted(list(job_dirs))):
        if not args.quiet:
            print(f"\n--- Processing directory {i+1}/{len(job_dirs)}: {job_dir} ---")
        
        try:
            if process_directory(job_dir, overwrite=args.overwrite, verbose=not args.quiet, charge=args.charge):
                success_count += 1
            else:
                error_count += 1
        except KeyboardInterrupt:
            print("\nInterrupted by user.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"  -> UNHANDLED EXCEPTION in {job_dir}: {e}", file=sys.stderr)
            error_count += 1
            
    if not args.quiet:
        print("\n========================================")
        print(f"🎉 All jobs complete.")
        print(f"   Successful directories: {success_count}")
        print(f"   Failed/skipped directories: {error_count}")
        print("========================================")

if __name__ == "__main__":
    main()

