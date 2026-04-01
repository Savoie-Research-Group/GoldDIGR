#!/usr/bin/env python3
from __future__ import annotations
import argparse
from typing import Optional, Tuple, List, Set
from pathlib import Path

from rxnclass_helper import (
    elt, pairs_from_reactive, metals_in_record,
    replace_task_with_deterministic, process_path,
    metal_before_after_electrons, electron_direction,
    rebuild_bond_changes_with_tags, split_required_optional_pairs_inplace,
)
from metal_ligand.labels import is_metal

import json

def _nonmetal_nonmetal(p: Tuple[str, str]) -> bool:
    return (not is_metal(p[0])) and (not is_metal(p[1]))


# -------------------- POSITIVE OA DETECTION (unchanged) --------------------
def _detect_oxidative_addition(rec: dict) -> Optional[str]:
    formed, broken = pairs_from_reactive(rec)
    # Core feature: a nonmetal–nonmetal σ bond (X–Y) is BROKEN.
    broken_nm: List[Tuple[str, str]] = [p for p in broken if _nonmetal_nonmetal(p)]
    if not broken_nm:
        return None

    formed_set = {tuple(p) for p in formed}
    candidate_M: Set[str] = set(metals_in_record(rec))
    if not candidate_M:
        # Very permissive fallback: any metal that appears in formed pairs
        for a, b in formed:
            if is_metal(a):
                candidate_M.add(a)
            if is_metal(b):
                candidate_M.add(b)

    # Relaxed criterion: it's enough that **one** of (M–X) or (M–Y) formed.
    for x, y in broken_nm:
        for m in candidate_M:
            mx = (m, x) if m <= x else (x, m)
            my = (m, y) if m <= y else (y, m)
            saw_mx = mx in formed_set
            saw_my = my in formed_set
            if saw_mx or saw_my:
                b, a = metal_before_after_electrons(rec, m)
                edir = electron_direction(b, a)
                # Enforce OA direction when data are present:
                if edir is not None and edir != "decrease":
                    continue
                e_msg = ""
                if edir is not None:
                    e_msg = f" Metal electrons: {b} → {a} (decrease)."

                # REQUIRED = {X–Y broken} plus ONE of (M–X) or (M–Y) if observed
                required_pairs: Set[frozenset[str]] = {frozenset((x, y))}
                if saw_mx:
                    required_pairs.add(frozenset((m, x)))
                if saw_my:
                    required_pairs.add(frozenset((m, y)))
                rebuild_bond_changes_with_tags(rec, required_pairs)
                split_required_optional_pairs_inplace(rec, required_pairs)  # optional but recommended

                which = "M–X" if saw_mx else "M–Y"
                rec["deterministic"] = {
                    "reaction_class": "oxidative_addition",
                    "explanation": (
                        f"{elt(x)}–{elt(y)} σ bond cleaves; at least one metal–ligand bond "
                        f"({which}) is observed formed → oxidative addition (relaxed).{e_msg}"
                    ),
                }
                replace_task_with_deterministic(rec)
                return "oxidative_addition"
    return None


def determine_oxidative_addition_inplace(rec: dict) -> Optional[dict]:
    return rec if _detect_oxidative_addition(rec) else None


# -------------------- NEGATIVE OA GENERATION --------------------
def _oa_failure_reasons(rec: dict) -> Tuple[List[str], Set[frozenset[str]]]:
    """
    Collect ALL violated OA rules (independently).
    Returns (reasons, required_pairs_if_any).
    """
    reasons: List[str] = []

    formed, broken = pairs_from_reactive(rec)
    formed_set = {tuple(p) for p in formed}

    # --- Rule A: must break a nonmetal–nonmetal σ bond (X–Y) ---
    broken_nm: List[Tuple[str, str]] = [p for p in broken if _nonmetal_nonmetal(p)]
    required_pairs: Set[frozenset[str]] = {frozenset(p) for p in broken_nm}
    if not broken_nm:
        reasons.append("No nonmetal–nonmetal σ bond (X–Y) is broken")

    # --- Candidate metals (present anywhere in the record or infer from formed) ---
    candidate_M: Set[str] = set(metals_in_record(rec))
    if not candidate_M:
        for a, b in formed:
            if is_metal(a): candidate_M.add(a)
            if is_metal(b): candidate_M.add(b)

    if not candidate_M:
        reasons.append("No metal center is present in the record")
    else:
        # --- Rule B: metal electrons should DECREASE (oxidation); check independently ---
        elec_known = False
        elec_any_decrease = False
        elec_any_violate = False
        bad_e_msg = None

        for m in candidate_M:
            before_e, after_e = metal_before_after_electrons(rec, m)
            edir = electron_direction(before_e, after_e)  # None | 'increase' | 'decrease' | 'no_change'
            if edir is None:
                continue
            elec_known = True
            if edir == "decrease":
                elec_any_decrease = True
            else:
                elec_any_violate = True
                if bad_e_msg is None:
                    dir_txt = "increase" if edir == "increase" else "no change"
                    bad_e_msg = f"Metal {m}: electrons {before_e} → {after_e} ({dir_txt})."

        if elec_known and (not elec_any_decrease) and elec_any_violate:
            reasons.append("Metal center is not oxidized (electron count does not decrease)" +
                           (f" {bad_e_msg}" if bad_e_msg else ""))

    # --- Rule C: for any broken X–Y, at least one of M–X or M–Y must form ---
    if broken_nm and candidate_M:
        any_mxmy = False
        for x, y in broken_nm:
            for m in candidate_M:
                mx = (m, x) if m <= x else (x, m)
                my = (m, y) if m <= y else (y, m)
                if mx in formed_set or my in formed_set:
                    any_mxmy = True
                    break
            if any_mxmy:
                break
        if not any_mxmy:
            reasons.append("Neither M–X nor M–Y σ bond is formed for any broken X–Y bond")

    return reasons, required_pairs

def determine_oxidative_addition_failed_inplace(rec: dict) -> Optional[dict]:
    """
    Return a *negative* OA record if OA was not detected, including rule-by-rule reasons.
    Never emits a failed record if the positive detector would classify as OA.
    """
    # --- Guard: if positive OA would pass, do NOT emit a failure ---
    # _detect_oxidative_addition mutates in-place, so use a deep copy.
    rec_copy = json.loads(json.dumps(rec))
    if _detect_oxidative_addition(rec_copy) is not None:
        return None

    # --- Evaluate all violated rules independently ---
    reasons, required_pairs = _oa_failure_reasons(rec)

    # If no reasons, do not emit a failed record.
    if not reasons:
        return None
    # Tag any broken X–Y σ bonds as 'required' for easier inspection.
    if required_pairs:
        rebuild_bond_changes_with_tags(rec, required_pairs)
        split_required_optional_pairs_inplace(rec, required_pairs)

    rec["deterministic"] = {
        "reaction_class": "not_oxidative_addition",
        "explanation": "; ".join(reasons),
        #"failed_rules": reasons,
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Oxidative addition checker (relaxed) with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()

    # 1) Try positive OA → writes OA-*.json if detected
    process_path(Path(args.path), determine_oxidative_addition_inplace, prefix="OA")

    # 2) Emit negative OA if not detected → writes OA_Failed-*.json
    process_path(Path(args.path), determine_oxidative_addition_failed_inplace, prefix="OA_Failed")

if __name__ == "__main__":
    main()

