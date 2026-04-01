#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from typing import Optional, Dict, Set, Tuple, List, FrozenSet
from pathlib import Path

from rxnclass_helper import (
    elt, pairs_from_reactive, metals_in_record,
    replace_task_with_deterministic, process_path,
    rebuild_bond_changes_with_tags, split_required_optional_pairs_inplace,
)
from metal_ligand.labels import is_metal

# -------------------- POSITIVE: Transmetalation --------------------
def _detect_transmetalation(rec: dict) -> Optional[str]:
    formed, broken = pairs_from_reactive(rec)
    formed_set = {tuple(p) for p in formed}

    # Map ligand L (non-metal) → metals that newly bind it
    lig_to_new_metals: Dict[str, Set[str]] = {}
    for a, b in formed:
        if is_metal(a) and not is_metal(b):
            lig_to_new_metals.setdefault(b, set()).add(a)
        elif is_metal(b) and not is_metal(a):
            lig_to_new_metals.setdefault(a, set()).add(b)

    # For any broken M1–L, if L forms to a different metal M2 → TM
    for a, b in broken:
        if is_metal(a) and not is_metal(b):
            m1, L = a, b
        elif is_metal(b) and not is_metal(a):
            m1, L = b, a
        else:
            continue
        for m2 in lig_to_new_metals.get(L, ()):
            if m2 == m1:
                continue
            required_pairs: Set[FrozenSet[str]] = {
                frozenset((m1, L)),  # broken
                frozenset((m2, L)),  # formed
            }
            rebuild_bond_changes_with_tags(rec, required_pairs)
            split_required_optional_pairs_inplace(rec, required_pairs)

            rec["deterministic"] = {
                "reaction_class": "transmetalation",
                "explanation": f"Ligand {elt(L)} transfers: {elt(m1)}–{elt(L)} breaks and {elt(m2)}–{elt(L)} forms.",
            }
            replace_task_with_deterministic(rec)
            return "transmetalation"
    return None

def determine_transmetalation_inplace(rec: dict) -> Optional[dict]:
    return rec if _detect_transmetalation(rec) else None

# -------------------- NEGATIVE: reasons --------------------
def _tm_failure_reasons(rec: dict) -> Tuple[List[str], Set[FrozenSet[str]]]:
    reasons: List[str] = []
    formed, broken = pairs_from_reactive(rec)
    if not formed and not broken:
        return (["No σ bond changes recorded (formed/broken are empty)"], set())

    cand_metals = set(metals_in_record(rec))
    if not cand_metals:
        # infer from pairs, if any metal appears
        for a, b in formed + broken:
            if is_metal(a): cand_metals.add(a)
            if is_metal(b): cand_metals.add(b)
    if not cand_metals:
        reasons.append("No metal centers present in the record")

    # any ligand L that leaves a metal?
    broken_ML: List[Tuple[str, str]] = []
    for a, b in broken:
        if is_metal(a) and not is_metal(b):
            broken_ML.append((a, b))
        elif is_metal(b) and not is_metal(a):
            broken_ML.append((b, a))
    if not broken_ML:
        reasons.append("No M–L bond is broken (ligand does not leave any metal)")

    # any L that binds a (different) metal?
    formed_ML: List[Tuple[str, str]] = []
    for a, b in formed:
        if is_metal(a) and not is_metal(b):
            formed_ML.append((a, b))
        elif is_metal(b) and not is_metal(a):
            formed_ML.append((b, a))
    if not formed_ML:
        reasons.append("No M–L bond is formed (ligand does not bind any metal)")

    # check “transfer” consistency: same L must appear in a broken M1–L and a formed M2–L with M2≠M1
    transfer_ok = False
    broken_by_L: Dict[str, Set[str]] = {}
    formed_by_L: Dict[str, Set[str]] = {}
    for m, L in broken_ML: broken_by_L.setdefault(L, set()).add(m)
    for m, L in formed_ML: formed_by_L.setdefault(L, set()).add(m)
    for L, m1s in broken_by_L.items():
        m2s = formed_by_L.get(L, set())
        if any(m2 for m2 in m2s if m2 not in m1s):
            transfer_ok = True
            break
    if (broken_ML or formed_ML) and not transfer_ok:
        reasons.append("Broken M–L and formed M–L do not involve the same ligand transferring to a different metal")

    # Required tagging (if we have candidate M–L pairs)
    required_pairs: Set[FrozenSet[str]] = {frozenset(p) for p in (broken_ML + formed_ML)}
    return reasons, required_pairs

def determine_transmetalation_failed_inplace(rec: dict) -> Optional[dict]:
    rec_copy = json.loads(json.dumps(rec))
    if _detect_transmetalation(rec_copy):
        return None
    reasons, req = _tm_failure_reasons(rec)
    if not reasons:
        return None
    if req:
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
    rec["deterministic"] = {
        "reaction_class": "not_transmetalation",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Transmetalation checker with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()
    process_path(Path(args.path), determine_transmetalation_inplace, prefix="TM")
    process_path(Path(args.path), determine_transmetalation_failed_inplace, prefix="TM_Failed")

if __name__ == "__main__":
    main()

