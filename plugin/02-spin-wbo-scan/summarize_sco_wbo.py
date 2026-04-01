#!/usr/bin/env python3
"""Summarize spin-crossover and WBO/BEM events from xTB-scan results.

Usage
-----
Run from a directory that contains `finished.txt`:

    python summarize_sco_wbo.py

By default this will:
  * read TS-level folders from finished.txt
  * find their parent folders
  * scan each parent for *.zip files that already contain xTB-scan/
  * for each such zip that has xTB-scan/reactive_summary.csv, emit one
    JSON record into SCO_WBO_summary.jsonl

You can override the finished file and output path with:

    python summarize_sco_wbo.py --finished path/to/finished.txt \
                                --output SCO_WBO_summary.jsonl
"""

from __future__ import annotations

import argparse
import io
import json
import math
import os
import sys
import zipfile
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Summarize TS, WBO changes, and spin crossovers into a JSONL file."
    )
    p.add_argument(
        "--finished",
        default="finished.txt",
        help="Path to finished.txt listing TS-level folders (default: finished.txt in CWD)",
    )
    p.add_argument(
        "--output",
        default="SCO_WBO_summary.jsonl",
        help="Output JSONL path (default: SCO_WBO_summary.jsonl in CWD)",
    )
    return p.parse_args()


def load_parent_dirs(finished_path: str) -> List[str]:
    """Read finished.txt and return unique parent directories (absolute paths)."""
    parents = set()
    with open(finished_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Take the first whitespace-separated token as the TS directory path
            ts_dir = line.split()[0]
            parent = os.path.dirname(ts_dir)
            if parent:
                parents.add(os.path.abspath(parent))
    return sorted(parents)


def find_member(namelist: List[str], suffix: str) -> Optional[str]:
    """Return the first member whose name ends with the given suffix, or None."""
    for name in namelist:
        if name.endswith(suffix):
            return name
    return None

def find_xtb_scan_member(namelist: List[str], filename: str) -> Optional[str]:
    """Return the member inside an xTB-scan/ folder with the given basename."""
    suffix = f"xTB-scan/{filename}"
    for name in namelist:
        if name.endswith(suffix):
            return name
    return None

def parse_ts_from_trj(file_bytes: bytes) -> Tuple[Optional[int], Optional[float], int]:
    """Parse finished_irc.trj content and return (ts_frame_idx, ts_energy, n_frames).

    The file is assumed to contain repeated XYZ blocks:
        line 1: integer number of atoms
        line 2: energy (float, possibly followed by comments)
        next N lines: coordinates

    The TS frame is defined as the frame with the highest energy.
    """
    text = file_bytes.decode("utf-8", errors="replace")
    lines = text.splitlines()

    energies: List[float] = []
    i = 0
    frame_idx = 0

    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            # If we can't parse natoms here, stop parsing further
            break
        if i + 1 >= len(lines):
            break
        energy_line = lines[i + 1].strip()
        tokens = energy_line.split()
        if not tokens:
            energy = float("nan")
        else:
            try:
                energy = float(tokens[0])
            except ValueError:
                energy = float("nan")
        energies.append(energy)
        # Skip coordinates: N atom lines
        i += 2 + natoms
        frame_idx += 1

    n_frames = len(energies)
    if n_frames == 0:
        return None, None, 0

    # Find index of max energy, ignoring NaNs
    valid_indices = [idx for idx, e in enumerate(energies) if not math.isnan(e)]
    if not valid_indices:
        return None, None, n_frames
    ts_idx = max(valid_indices, key=lambda idx: energies[idx])
    ts_energy = energies[ts_idx]
    return ts_idx, ts_energy, n_frames


def compute_bond_slopes(df: pd.DataFrame, top_n: int = 3) -> Dict[str, dict]:
    """Compute first-derivative-like slopes for each bond column.

    We:
      * sort by the frame column
      * for each bond column, compute a slope at each frame:
          - frame 0: forward difference: v[1] - v[0]
          - last frame: backward difference: v[-1] - v[-2]
          - interior frames: central difference: (v[i+1] - v[i-1]) / 2
      * find the top |slope| values (up to top_n) and record their frames

    Returns a dict mapping column name -> {
        "first": float,
        "last": float,
        "delta": float,
        "top_slope_events": [
            {"frame": int, "slope": float, "value": float}, ...
        ]
    }
    """
    if "frame" not in df.columns:
        raise ValueError("reactive_summary.csv is missing a 'frame' column")

    df_sorted = df.sort_values("frame").reset_index(drop=True)
    frames = df_sorted["frame"].to_numpy()
    result: Dict[str, dict] = {}

    for col in df_sorted.columns:
        if col == "frame":
            continue
        vals = df_sorted[col].to_numpy(dtype=float)
        n = len(vals)
        if n == 0:
            continue

        slopes = np.zeros(n, dtype=float)
        if n == 1:
            slopes[0] = 0.0
        elif n == 2:
            slopes[0] = vals[1] - vals[0]
            slopes[1] = slopes[0]
        else:
            slopes[0] = vals[1] - vals[0]
            slopes[-1] = vals[-1] - vals[-2]
            slopes[1:-1] = (vals[2:] - vals[:-2]) / 2.0

        order = np.argsort(-np.abs(slopes))  # descending by |slope|
        k = min(top_n, n)
        top_idx = order[:k]

        events = []
        for idx in top_idx:
            events.append(
                {
                    "frame": int(frames[idx]),
                    "slope": float(slopes[idx]),
                    "value": float(vals[idx]),
                }
            )

        result[col] = {
            "first": float(vals[0]),
            "last": float(vals[-1]),
            "delta": float(vals[-1] - vals[0]),
            "top_slope_events": events,
        }

    return result


def load_spin_crossover_info(
    z: zipfile.ZipFile, namelist: List[str]
) -> Optional[dict]:
    """Load spin crossover info (if present) from the zip.

    Looks for:
      * spin_crossover_frames.txt
      * spin_scan_summary.csv

    Returns a dict of the form:
        {
          "frames": [int, ...],
          "events": [
             {
               "frame": int,
               "from_best_UHF": Optional[int],
               "to_best_UHF": Optional[int],
               "from_multiplicity": Optional[int],
               "to_multiplicity": Optional[int],
             }, ...
          ],
        }
    or None if no spin_crossover_frames.txt is present.
    """
    sco_member = find_xtb_scan_member(namelist, "spin_crossover_frames.txt")
    if sco_member is None:
        return None

    # Read crossover frames
    with z.open(sco_member, "r") as f:
        text = f.read().decode("utf-8", errors="replace")
    frames: List[int] = []
    for tok in text.replace(",", " ").split():
        try:
            frames.append(int(tok))
        except ValueError:
            continue
    frames = sorted(set(frames))
    if not frames:
        return None

    # Load spin_scan_summary.csv to get best_UHF vs frame
    spin_member = find_xtb_scan_member(namelist, "spin_scan_summary.csv")
    if spin_member is None:
        # We know frames but not spins; still record frames
        return {"frames": frames, "events": []}

    with z.open(spin_member, "r") as f:
        df_spin = pd.read_csv(io.TextIOWrapper(f, encoding="utf-8"))

    if "frame" not in df_spin.columns or "best_UHF" not in df_spin.columns:
        return {"frames": frames, "events": []}

    frame_to_best = {
        int(row["frame"]): int(row["best_UHF"])
        for _, row in df_spin.iterrows()
        if not pd.isna(row["frame"]) and not pd.isna(row["best_UHF"])
    }

    events = []
    for fr in frames:
        to_spin = frame_to_best.get(fr)
        from_spin = frame_to_best.get(fr - 1)
        event = {
            "frame": int(fr),
            "from_best_UHF": int(from_spin) if from_spin is not None else None,
            "to_best_UHF": int(to_spin) if to_spin is not None else None,
            "from_multiplicity": int(from_spin + 1) if from_spin is not None else None,
            "to_multiplicity": int(to_spin + 1) if to_spin is not None else None,
        }
        events.append(event)

    return {"frames": frames, "events": events}


def process_zip(zip_path: str) -> Optional[dict]:
    """Process a single .zip archive.

    Returns a JSON-serializable dict summarizing TS, WBO slopes, and spin
    crossovers, or None if the archive should be skipped (e.g., no
    xTB-scan/reactive_summary.csv).
    """
    try:
        z = zipfile.ZipFile(zip_path, "r")
    except Exception as e:
        print(f"[warn] Failed to open zip {zip_path}: {e}", file=sys.stderr)
        return None

    with z:
        namelist = z.namelist()

        # Skip zips that do not contain xTB-scan/ at all
        if not any("xTB-scan/" in name for name in namelist):
            return None

        # TS info from finished_irc.trj
        trj_member = find_member(namelist, "finished_irc.trj")
        ts_frame = None
        ts_energy = None
        n_frames = None
        if trj_member is not None:
            with z.open(trj_member, "r") as f:
                ts_frame, ts_energy, n_frames = parse_ts_from_trj(f.read())

        # Reactive summary (required to emit a record)
        reactive_member = find_xtb_scan_member(namelist, "reactive_summary.csv")
        if reactive_member is None:
            return None

        with z.open(reactive_member, "r") as f:
            df_reactive = pd.read_csv(io.TextIOWrapper(f, encoding="utf-8"))

        bond_summary = compute_bond_slopes(df_reactive, top_n=3)

        rec: dict = {
            "zip_path": os.path.abspath(zip_path),
            "ts_frame": int(ts_frame) if ts_frame is not None else None,
            "ts_energy": float(ts_energy) if ts_energy is not None else None,
            "n_frames": int(n_frames) if n_frames is not None else None,
            "reactive_bonds": bond_summary,
        }

        # Spin crossover (optional)
        sco_info = load_spin_crossover_info(z, namelist)
        if sco_info is not None:
            rec["spin_crossover"] = sco_info

        return rec


def main() -> None:
    args = parse_args()
    finished_path = os.path.abspath(args.finished)

    if not os.path.exists(finished_path):
        print(f"[error] finished file not found: {finished_path}", file=sys.stderr)
        sys.exit(1)

    parent_dirs = load_parent_dirs(finished_path)
    if not parent_dirs:
        print("[warn] No parent directories found from finished.txt", file=sys.stderr)

    out_path = os.path.abspath(args.output)
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    n_zips_scanned = 0
    n_records = 0

    with open(out_path, "w") as out_f:
        for parent in parent_dirs:
            if not os.path.isdir(parent):
                print(f"[warn] Parent directory does not exist: {parent}", file=sys.stderr)
                continue

            for name in sorted(os.listdir(parent)):
                if not name.lower().endswith(".zip"):
                    continue
                zip_path = os.path.join(parent, name)
                n_zips_scanned += 1

                rec = process_zip(zip_path)
                if rec is None:
                    continue

                out_f.write(json.dumps(rec) + "\n")
                n_records += 1
                print(f"[info] Recorded summary for {zip_path}", file=sys.stderr)

    print(
        f"Done. Scanned {n_zips_scanned} zip files, wrote {n_records} records to {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

