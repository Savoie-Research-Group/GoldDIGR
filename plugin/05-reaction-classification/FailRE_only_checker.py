#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from typing import Optional, Tuple, List, Set, FrozenSet
from pathlib import Path

from rxnclass_helper import (
    elt, pairs_from_reactive, metals_in_record,
    replace_task_with_deterministic, process_path,
    metal_before_after_electrons, electron_direction,
    rebuild_bond_changes_with_tags, split_required_optional_pairs_inplace,
)
from metal_ligand.labels import is_metal

def _nonmetal_nonmetal(p: Tuple[str,str]) -> bool:
    return (not is_metal(p[0])) and (not is_metal(p[1]))

# -------------------- POSITIVE --------------------
def _detect_reductive_elimination(rec: dict) -> Optional[str]:
    formed, broken = pairs_from_reactive(rec)
    formed_nm = [p for p in formed if _nonmetal_nonmetal(p)]
    if not formed_nm:
        return None

    broken_set = {tuple(p) for p in broken}
    cand_M: Set[str] = set(metals_in_record(rec))
    if not cand_M:
        for a, b in broken:
            if is_metal(a): cand_M.add(a)
            if is_metal(b): cand_M.add(b)

    for x, y in formed_nm:
        for m in cand_M:
            mx = (m, x) if m <= x else (x, m)
            my = (m, y) if m <= y else (y, m)
            saw_mx = mx in broken_set
            saw_my = my in broken_set
            if not (saw_mx or saw_my):
                continue
            bE, aE = metal_before_after_electrons(rec, m)
            ed = electron_direction(bE, aE)
            if ed is not None and ed != "increase":
                continue
            e_msg = f" Metal electrons: {bE} → {aE} (increase)." if ed is not None else ""
            required = {frozenset((x,y))}
            if saw_mx: required.add(frozenset((m,x)))
            if saw_my: required.add(frozenset((m,y)))
            rebuild_bond_changes_with_tags(rec, required)
            split_required_optional_pairs_inplace(rec, required)
            which = "M–X" if saw_mx else "M–Y"
            rec["deterministic"] = {
                "reaction_class": "reductive_elimination",
                "explanation": (f"New {elt(x)}–{elt(y)} σ bond forms; at least one metal–ligand "
                                f"bond ({which}) is broken; metal electron count increases if available.{e_msg}")
            }
            replace_task_with_deterministic(rec)
            return "reductive_elimination"
    return None

def determine_reductive_elimination_inplace(rec: dict) -> Optional[dict]:
    return rec if _detect_reductive_elimination(rec) else None

# -------------------- NEGATIVE --------------------
def _re_failure_reasons(rec: dict) -> Tuple[List[str], Set[FrozenSet[str]]]:
    reasons: List[str] = []
    formed, broken = pairs_from_reactive(rec)
    if not formed and not broken:
        return (["No σ bond changes recorded (formed/broken are empty)"], set())

    formed_nm = [p for p in formed if _nonmetal_nonmetal(p)]
    if not formed_nm:
        reasons.append("No nonmetal–nonmetal σ bond (X–Y) is formed")

    cand_M: Set[str] = set(metals_in_record(rec))
    if not cand_M:
        for a,b in broken:
            if is_metal(a): cand_M.add(a)
            if is_metal(b): cand_M.add(b)
    if not cand_M:
        reasons.append("No metal center is present in the record")

    # need at least one of M–X or M–Y broken for a formed X–Y
    if formed_nm and cand_M:
        any_ok = False
        broken_set = {tuple(p) for p in broken}
        for x, y in formed_nm:
            for m in cand_M:
                mx = (m, x) if m <= x else (x, m)
                my = (m, y) if m <= y else (y, m)
                if (mx in broken_set) or (my in broken_set):
                    any_ok = True
                    break
            if any_ok: break
        if not any_ok:
            reasons.append("Neither M–X nor M–Y σ bond is broken for any formed X–Y bond")

    # electron count should increase if present
    elec_known=False; elec_any_increase=False; elec_any_violate=False; bad=None
    for m in cand_M:
        bE,aE = metal_before_after_electrons(rec,m)
        ed = electron_direction(bE,aE)
        if ed is None: continue
        elec_known = True
        if ed=="increase":
            elec_any_increase = True
        else:
            elec_any_violate = True
            if bad is None:
                dir_txt = "decrease" if ed=="decrease" else "no change"
                bad = f"Metal {m}: electrons {bE} → {aE} ({dir_txt})."
    if elec_known and (not elec_any_increase) and elec_any_violate:
        reasons.append("Metal center is not reduced (electron count does not increase) " + (bad or ""))

    # tag pairs
    req: Set[FrozenSet[str]] = set()
    for x,y in formed_nm: req.add(frozenset((x,y)))
    for a,b in broken:
        if is_metal(a) ^ is_metal(b):
            req.add(frozenset((a,b)))
    return reasons, req

def determine_reductive_elimination_failed_inplace(rec: dict) -> Optional[dict]:
    rec_copy = json.loads(json.dumps(rec))
    if _detect_reductive_elimination(rec_copy):
        return None
    reasons, req = _re_failure_reasons(rec)
    if not reasons:
        return None
    if req:
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
    rec["deterministic"] = {
        "reaction_class": "not_reductive_elimination",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Reductive elimination checker with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()
    process_path(Path(args.path), determine_reductive_elimination_inplace, prefix="RE")
    process_path(Path(args.path), determine_reductive_elimination_failed_inplace, prefix="RE_Failed")

if __name__ == "__main__":
    main()

