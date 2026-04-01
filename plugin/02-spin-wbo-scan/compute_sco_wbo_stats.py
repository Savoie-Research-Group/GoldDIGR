#!/usr/bin/env python3
"""
Compute statistics from SCO_WBO_summary.jsonl.

For each reaction (one JSON line), we look at:
  - ts_frame (transition state frame)
  - reactive_bonds[...]["top_slope_events"] (top-3 bond-change frames)
  - spin_crossover["frames"] (spin-crossover frames, if present)

We then report:
  1. Among reactions with a TS frame + reactive bonds:
       - % where ts_frame lies in top-1 bond frames
       - % where ts_frame lies in top-2 (ranks 1 or 2) bond frames
       - % where ts_frame lies in top-3 (ranks 1, 2, or 3) bond frames

  2. Among reactions with spin crossover + reactive bonds:
       - % where any SCO frame lies in top-1 bond frames
       - % where any SCO frame lies in top-2 bond frames
       - % where any SCO frame lies in top-3 bond frames

Usage:
    python compute_sco_wbo_stats.py SCO_WBO_summary.jsonl
"""

import argparse
import json
from typing import Set, Tuple


def collect_bond_event_frames(reactive_bonds: dict) -> Tuple[Set[int], Set[int], Set[int]]:
    """
    From reactive_bonds, build three sets of frames:

        top1_frames: union of rank-1 frames across all bonds
        top2_frames: union of rank-1 and rank-2 frames
        top3_frames: union of rank-1, rank-2, and rank-3 frames

    We assume that for each bond:
        bond_info["top_slope_events"] is a list sorted by descending |slope|.
    """
    top1_frames: Set[int] = set()
    top2_frames: Set[int] = set()
    top3_frames: Set[int] = set()

    for bond_name, bond_info in reactive_bonds.items():
        events = bond_info.get("top_slope_events", [])
        for rank, ev in enumerate(events):
            frame = ev.get("frame")
            if frame is None:
                continue

            # rank 0 = top1; 1 = top2; 2 = top3
            if rank == 0:
                top1_frames.add(frame)
                top2_frames.add(frame)
                top3_frames.add(frame)
            elif rank == 1:
                top2_frames.add(frame)
                top3_frames.add(frame)
            elif rank == 2:
                top3_frames.add(frame)
            # ignore any events beyond rank-3 if present

    return top1_frames, top2_frames, top3_frames


def main():
    parser = argparse.ArgumentParser(
        description="Compute TS and spin-crossover alignment stats from SCO_WBO_summary.jsonl."
    )
    parser.add_argument(
        "jsonl",
        help="Path to SCO_WBO_summary.jsonl produced by summarize_sco_wbo.py",
    )
    args = parser.parse_args()

    ts_total = ts_top1 = ts_top2 = ts_top3 = 0
    sco_total = sco_top1 = sco_top2 = sco_top3 = 0

    with open(args.jsonl, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            rec = json.loads(line)

            reactive_bonds = rec.get("reactive_bonds", {})
            if not reactive_bonds:
                # No bond events => nothing to compare against
                continue

            top1_frames, top2_frames, top3_frames = collect_bond_event_frames(reactive_bonds)

            # --- TS statistics ---
            ts_frame = rec.get("ts_frame")
            if ts_frame is not None:
                ts_total += 1
                if ts_frame in top1_frames:
                    ts_top1 += 1
                if ts_frame in top2_frames:
                    ts_top2 += 1
                if ts_frame in top3_frames:
                    ts_top3 += 1

            # --- Spin-crossover statistics ---
            sco = rec.get("spin_crossover")
            sco_frames = []
            if sco and isinstance(sco, dict):
                frames = sco.get("frames", [])
                if isinstance(frames, list):
                    # make sure these are ints
                    sco_frames = [int(fr) for fr in frames if fr is not None]

            if sco_frames:
                sco_total += 1
                if any(fr in top1_frames for fr in sco_frames):
                    sco_top1 += 1
                if any(fr in top2_frames for fr in sco_frames):
                    sco_top2 += 1
                if any(fr in top3_frames for fr in sco_frames):
                    sco_top3 += 1

    # --- Report ---
    def pct(num: int, den: int) -> float:
        return 100.0 * num / den if den > 0 else 0.0

    print("=== TS vs bond-event frames ===")
    print(f"Reactions with TS + reactive bonds: {ts_total}")
    if ts_total > 0:
        print(f"TS in top-1 bond frames: {ts_top1} ({pct(ts_top1, ts_total):.2f} %)")
        print(f"TS in top-2 bond frames: {ts_top2} ({pct(ts_top2, ts_total):.2f} %)")
        print(f"TS in top-3 bond frames: {ts_top3} ({pct(ts_top3, ts_total):.2f} %)")
    else:
        print("No reactions with both ts_frame and reactive_bonds.")

    print("\n=== Spin-crossover vs bond-event frames ===")
    print(f"Reactions with spin-crossover + reactive bonds: {sco_total}")
    if sco_total > 0:
        print(f"SCO in top-1 bond frames: {sco_top1} ({pct(sco_top1, sco_total):.2f} %)")
        print(f"SCO in top-2 bond frames: {sco_top2} ({pct(sco_top2, sco_total):.2f} %)")
        print(f"SCO in top-3 bond frames: {sco_top3} ({pct(sco_top3, sco_total):.2f} %)")
    else:
        print("No reactions with both spin_crossover and reactive_bonds.")


if __name__ == "__main__":
    main()

