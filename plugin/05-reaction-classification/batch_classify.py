#!/usr/bin/env python3
"""
batch_classify.py — Batch classification of reactions from tar.gz workflow archives.

Extracts BEL graphs from workflow_status.json files inside tar.gz archives,
deduplicates by canonical graph hash, runs all reaction classifiers, and
appends results to a JSONL registry.

Usage:
    # Process a directory of tar.gz files
    python batch_classify.py /path/to/targz_dir -o results.jsonl

    # Process specific tar.gz files
    python batch_classify.py file1.tar.gz file2.tar.gz -o results.jsonl

    # With custom duplicates log
    python batch_classify.py /path/to/targz_dir -o results.jsonl -d dupes.jsonl

Output (results.jsonl): one JSON line per unique reaction:
    {
        "hash": "abc123...",
        "source_tar": "path/to/file.tar.gz",
        "idx": 273640,
        "direction": "Forward",
        "reactant_bel": "...",
        "product_bel": "...",
        "bond_changes_raw": "N1-C13:1>2;C13-O43:1>n;...",
        "bond_changes_sigma": "C13-O43:1>n;O42-H44:1>n;O43-H44:n>1",
        "positive_classes": ["reductive_elimination"],
        "classifications": { ... }
    }

Duplicates log (dupes.jsonl): one JSON line per skipped duplicate:
    {
        "hash": "abc123...",
        "skipped_tar": "path/to/dup.tar.gz",
        "skipped_idx": 999,
        "original_tar": "path/to/first.tar.gz"
    }
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tarfile
from io import BytesIO
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Checker registry
# ---------------------------------------------------------------------------

# (display_name, module_filename, positive_func_name)
_CHECKER_DEFS: List[Tuple[str, str, str]] = [
    ("OA",        "FailOA_only_checker.py",          "determine_oxidative_addition_inplace"),
    ("RE",        "FailRE_only_checker.py",          "determine_reductive_elimination_inplace"),
    ("MI",        "FailMI_only_checker.py",          "determine_migratory_insertion_inplace"),
    ("BXA",       "FailBXA_only_checker.py",         "determine_beta_atom_elimination_inplace"),
    ("BXA_noH",   "FailBXA-noBHE_only_checker.py",  "determine_beta_atom_elimination_inplace"),
    ("BXE",       "FailBXE_only_checker.py",         "determine_beta_heteroatom_elimination_inplace"),
    ("CH",        "FailCH_only_checker.py",          "determine_ch_activation_inplace"),
    ("MET",       "FailMET_only_checker.py",         "determine_metathesis_inplace"),
    ("TM",        "FailTM_only_checker.py",          "determine_transmetalation_inplace"),
]


def _load_func(module_path: Path, func_name: str) -> Callable:
    """Import a single function from a .py file (handles hyphens in filenames)."""
    spec = importlib.util.spec_from_file_location("_checker", str(module_path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, func_name)


def load_checkers(pkg_dir: Path) -> List[Tuple[str, Callable]]:
    """
    Load all positive-detector functions from the checker scripts.
    Returns [(name, func), ...].  Skips any checkers whose file is missing.
    """
    loaded = []
    for name, filename, funcname in _CHECKER_DEFS:
        fpath = pkg_dir / filename
        if not fpath.exists():
            print(f"  [warn] checker not found, skipping: {fpath}", file=sys.stderr)
            continue
        func = _load_func(fpath, funcname)
        loaded.append((name, func))
    return loaded


# ---------------------------------------------------------------------------
# Tar.gz extraction
# ---------------------------------------------------------------------------

_WF_PATH_IN_TAR = "output/workflow_status.json"


def iter_workflow_jsons(paths: List[Path]):
    """
    Yield (tar_path, workflow_dict) for every valid workflow_status.json
    found inside the given tar.gz files or directories of tar.gz files.
    """
    tar_files: List[Path] = []
    for p in paths:
        if p.is_dir():
            tar_files.extend(sorted(p.rglob("*.tar.gz")))
        elif p.suffix == ".gz" or p.name.endswith(".tar.gz"):
            tar_files.append(p)
        else:
            print(f"  [skip] not a tar.gz: {p}", file=sys.stderr)

    for tf_path in tar_files:
        try:
            with tarfile.open(tf_path, "r:gz") as tf:
                # Try the expected path; also try without leading directory
                candidate_names = [_WF_PATH_IN_TAR]
                for member in tf.getnames():
                    if member.endswith("workflow_status.json"):
                        candidate_names.append(member)

                for name in candidate_names:
                    try:
                        f = tf.extractfile(name)
                    except (KeyError, AttributeError):
                        continue
                    if f is None:
                        continue
                    try:
                        data = json.loads(f.read())
                        yield (tf_path, data)
                        break  # one workflow per tar
                    except json.JSONDecodeError:
                        continue
        except (tarfile.TarError, OSError) as exc:
            print(f"  [error] {tf_path}: {exc}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Registry (seen hashes)
# ---------------------------------------------------------------------------

def load_registry(jsonl_path: Path) -> Dict[str, str]:
    """
    Load existing results JSONL → {hash: source_tar}.
    Returns empty dict if file doesn't exist.
    """
    seen: Dict[str, str] = {}
    if not jsonl_path.exists():
        return seen
    with jsonl_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                h = rec.get("hash")
                if h:
                    seen[h] = rec.get("source_tar", "")
            except json.JSONDecodeError:
                continue
    return seen


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_reaction(
    rec: dict,
    checkers: List[Tuple[str, Callable]],
) -> Tuple[List[str], dict]:
    """
    Run all checkers on a deep copy of rec.
    Returns (positive_class_names, {name: deterministic_block}).
    """
    positive_names: List[str] = []
    classifications: dict = {}
    for name, func in checkers:
        rec_copy = json.loads(json.dumps(rec))
        result = func(rec_copy)
        if result is not None:
            det = result.get("deterministic", {})
            positive_names.append(name)
            classifications[name] = det
    return positive_names, classifications


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------

def process_all(
    input_paths: List[Path],
    output_path: Path,
    dupes_path: Path,
    pkg_dir: Path,
) -> None:
    from bel_parser import parse_bel, reaction_hash, filter_sigma_changes, build_rec_from_bel

    # Load checkers
    checkers = load_checkers(pkg_dir)
    print(f"Loaded {len(checkers)} checkers: {[n for n, _ in checkers]}")

    # Load existing registry
    seen = load_registry(output_path)
    print(f"Registry has {len(seen)} existing entries from {output_path}")

    stats = {"processed": 0, "skipped_no_irc": 0, "skipped_no_sigma": 0, "duplicates": 0, "errors": 0, "classified": 0}
    class_counts: Dict[str, int] = {}

    out_fh = output_path.open("a")
    dup_fh = dupes_path.open("a")

    try:
        for tar_path, ws in iter_workflow_jsons(input_paths):
            # ---- Filter: must have irc_results ----
            irc = ws.get("irc_results")
            status = ws.get("status", "")
            if not irc:
                stats["skipped_no_irc"] += 1
                continue

            r_bel = irc.get("reactant_geom_graph", "")
            p_bel = irc.get("product_geom_graph", "")
            if not r_bel or not p_bel:
                stats["skipped_no_irc"] += 1
                continue

            # ---- Parse & hash ----
            try:
                r_atoms, r_bonds = parse_bel(r_bel)
                p_atoms, p_bonds = parse_bel(p_bel)
                h = reaction_hash(r_atoms, r_bonds, p_atoms, p_bonds)
            except Exception as exc:
                print(f"  [error] parse failed for {tar_path}: {exc}", file=sys.stderr)
                stats["errors"] += 1
                continue

            # ---- Dedup ----
            if h in seen:
                dup_rec = {
                    "hash": h,
                    "skipped_tar": str(tar_path.resolve()),
                    "skipped_idx": (ws.get("input_data") or {}).get("idx"),
                    "original_tar": seen[h],
                }
                dup_fh.write(json.dumps(dup_rec) + "\n")
                stats["duplicates"] += 1
                continue

            # ---- σ-filter bond_changes (skip if no σ-crossing) ----
            bc_raw = irc.get("bond_changes", "")
            bc_sigma = filter_sigma_changes(bc_raw)
            if not bc_sigma:
                stats["skipped_no_sigma"] += 1
                continue

            # ---- Build rec & classify ----
            try:
                source = str(tar_path.resolve())
                rec = build_rec_from_bel(r_bel, p_bel, source=source)
                positive_names, classifications = classify_reaction(rec, checkers)
            except Exception as exc:
                print(f"  [error] classify failed for {tar_path}: {exc}", file=sys.stderr)
                stats["errors"] += 1
                continue

            # ---- Barriers (truncate to 2 decimal places) ----
            def _trunc2(v):
                if v is None:
                    return None
                try:
                    return round(float(v), 2)
                except (ValueError, TypeError):
                    return None

            irc_barrier = _trunc2(irc.get("irc_barrier"))
            ssm_barrier = _trunc2(ws.get("ssm_barrier"))

            # ---- Write result ----
            input_data = ws.get("input_data") or {}
            out_rec = {
                "hash": h,
                "source_tar": str(tar_path.resolve()),
                "idx": input_data.get("idx"),
                "direction": input_data.get("direction"),
                "parent_path": input_data.get("parent_path"),
                "mode": input_data.get("mode"),
                "source": input_data.get("source"),
                "prediction": input_data.get("prediction"),
                "ssm_barrier": ssm_barrier,
                "irc_barrier": irc_barrier,
                "reactant_bel": r_bel,
                "product_bel": p_bel,
                "bond_changes_raw": bc_raw,
                "bond_changes_sigma": bc_sigma,
                "positive_classes": positive_names,
                "classifications": classifications,
            }
            out_fh.write(json.dumps(out_rec) + "\n")
            out_fh.flush()

            seen[h] = str(tar_path.resolve())
            stats["classified"] += 1
            stats["processed"] += 1

            # Track per-class counts
            if positive_names:
                for cname in positive_names:
                    class_counts[cname] = class_counts.get(cname, 0) + 1
            else:
                class_counts["no_class"] = class_counts.get("no_class", 0) + 1

            if stats["processed"] % 100 == 0:
                print(f"  ... processed {stats['processed']} reactions", file=sys.stderr)

    finally:
        out_fh.close()
        dup_fh.close()

    print(f"\nDone.  classified={stats['classified']}  duplicates={stats['duplicates']}  "
          f"skipped_no_irc={stats['skipped_no_irc']}  skipped_no_sigma={stats['skipped_no_sigma']}  "
          f"errors={stats['errors']}")

    # ---- Classification summary ----
    if class_counts:
        print("\nClassification counts:")
        for cname in sorted(class_counts, key=lambda k: (k == "no_class", k)):
            print(f"  {cname} = {class_counts[cname]}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Batch-classify reactions from tar.gz workflow archives.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "inputs", nargs="*", type=Path,
        help="Tar.gz files or directories containing them.",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=Path("results.jsonl"),
        help="Output JSONL registry (append-only).  Default: results.jsonl",
    )
    ap.add_argument(
        "-d", "--duplicates", type=Path, default=Path("duplicates.jsonl"),
        help="Duplicates log (append-only).  Default: duplicates.jsonl",
    )
    ap.add_argument(
        "--pkg-dir", type=Path, default=None,
        help="Directory containing the *_only_checker.py scripts.  "
             "Default: same directory as this script.",
    )
    ap.add_argument(
        "--chunk-file", type=Path, default=None,
        help="Text file listing tar.gz paths (one per line).  "
             "Used by job arrays instead of positional inputs.",
    )
    args = ap.parse_args()

    # Build input list: --chunk-file takes priority, then positional args
    input_paths: List[Path] = []
    if args.chunk_file:
        with args.chunk_file.open() as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    input_paths.append(Path(line))
    if args.inputs:
        input_paths.extend(args.inputs)
    if not input_paths:
        ap.error("No inputs: provide positional paths or --chunk-file.")

    pkg = args.pkg_dir or Path(__file__).resolve().parent
    if not (pkg / "rxnclass_helper.py").exists():
        ap.error(f"Cannot find rxnclass_helper.py in {pkg}")

    # Ensure the package dir is importable
    if str(pkg) not in sys.path:
        sys.path.insert(0, str(pkg))

    process_all(input_paths, args.output, args.duplicates, pkg)


if __name__ == "__main__":
    main()

