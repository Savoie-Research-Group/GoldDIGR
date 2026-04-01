#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from typing import Optional, Tuple, List, Set, Dict, FrozenSet
from pathlib import Path

from rxnclass_helper import (
    elt, pairs_from_reactive, replace_task_with_deterministic, process_path,
    _idx_label_maps, max_sigma_order_by_atom, rebuild_bond_changes_with_tags,
    split_required_optional_pairs_inplace,
    metal_before_after_electrons, electron_direction,
)
from metal_ligand.labels import is_metal

def _is_H(x:str)->bool: return x.startswith("H")

# ---- utilities copied (trimmed) from your MI draft ----
def _parse_change(ch:str)->Tuple[Optional[int], Optional[int]]:
    m = re.match(r"\s*(NaN|-?\d+)\s*->\s*(NaN|-?\d+)\s*$", str(ch))
    if not m: return None, None
    c = lambda t: None if t=="NaN" else int(t)
    return c(m.group(1)), c(m.group(2))

def _collect_sigma_weaken(rec: dict) -> Tuple[Set[FrozenSet[str]], Set[FrozenSet[str]]]:
    broken_sigma, sigma_to_dat = set(), set()
    tp = (rec.get("target") or {})
    for item in (tp.get("reactive_pairs") or []):
        if isinstance(item,(list,tuple)) and len(item)>=3: i,j,ch=item[0],item[1],item[2]
        elif isinstance(item,dict): i,j,ch=item.get("i"),item.get("j"),item.get("change")
        else: continue
        frm,to = _parse_change(ch)
        if frm is not None and frm>=1 and to is None:
            broken_sigma.add(frozenset((str(i),str(j))))
        if frm == 1 and to == 0:
            sigma_to_dat.add(frozenset((str(i),str(j))))
    return broken_sigma, sigma_to_dat

def _has_ch_transfer(formed: List[Tuple[str,str]], broken: List[Tuple[str,str]]) -> bool:
    # SAME H appears in broken X–H (X non-metal) and formed C–H
    hx = {}
    for a,b in broken:
        if _is_H(a) and not is_metal(b): hx[a]=b
        if _is_H(b) and not is_metal(a): hx[b]=a
    for a,b in formed:
        if _is_H(a) and not is_metal(b) and a in hx and hx[a]!=b: return True
        if _is_H(b) and not is_metal(a) and b in hx and hx[b]!=a: return True
    return False

def _choose_R(a,b,ra,pa,rb,pb,broken_sigma,sigma_to_dat)->Optional[str]:
    a_unsat2sat, b_unsat2sat = (ra>1 and pa==1), (rb>1 and pb==1)
    if a_unsat2sat ^ b_unsat2sat:
        return b if a_unsat2sat else a
    for R in (a,b):
        for p in (broken_sigma|sigma_to_dat):
            if R in p:
                oth = next(iter(p-{R}))
                if is_metal(oth):
                    return R
    return None

def _find_M_for_R(R:str, broken_sigma:Set[FrozenSet[str]], sigma_to_dat:Set[FrozenSet[str]])->Tuple[Optional[str],Optional[str]]:
    for p in broken_sigma:
        if R in p:
            oth = next(iter(p-{R}))
            if is_metal(oth): return oth, "broken"
    for p in sigma_to_dat:
        if R in p:
            oth = next(iter(p-{R}))
            if is_metal(oth): return oth, "sigma->dative"
    return None, None

# -------------------- POSITIVE --------------------
def _detect_mi(rec: dict) -> Optional[str]:
    formed, broken = pairs_from_reactive(rec)
    if not formed:
        return None
    if _has_ch_transfer(formed, broken):
        return None  # CH/PDM-like — not MI

    react_max = max_sigma_order_by_atom(rec, "reactant")
    prod_max  = max_sigma_order_by_atom(rec, "product")

    for a, b in formed:
        if is_metal(a) or is_metal(b):
            continue  # witness must be X–Y
        ra, rb = react_max.get(a,1), react_max.get(b,1)
        pa, pb = prod_max.get(a,1),  prod_max.get(b,1)
        if not ((ra>1 and pa==1) or (rb>1 and pb==1)):
            continue  # need unsat→sat on this formed pair

        # hydrometallation allowance: if H in formed pair, require broken M–H
        if _is_H(a) or _is_H(b):
            if not any(((_is_H(x) and is_metal(y)) or (_is_H(y) and is_metal(x))) for x,y in broken):
                continue

        broken_sigma, sigma_to_dat = _collect_sigma_weaken(rec)
        R = _choose_R(a,b,ra,pa,rb,pb,broken_sigma,sigma_to_dat)
        if R is None: 
            continue
        M, weaken_kind = _find_M_for_R(R, broken_sigma, sigma_to_dat)
        if M is None:
            continue

        # redox-neutrality at that M if e data exist
        bE, aE = metal_before_after_electrons(rec, M)
        if bE is not None and aE is not None and bE != aE:
            continue

        req = {frozenset((a,b)), frozenset((M,R))}
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
        wtxt = "M–R broken" if weaken_kind=="broken" else "M–R σ→dative"
        rec["deterministic"] = {
            "reaction_class": "migratory_insertion",
            "explanation": (
                f"{elt(a)}–{elt(b)} σ forms with unsaturated→saturated change on one partner; "
                f"{wtxt} on {M}–{R}; metal redox-neutral if known."
            )
        }
        replace_task_with_deterministic(rec)
        return "migratory_insertion"
    return None

def determine_migratory_insertion_inplace(rec: dict) -> Optional[dict]:
    return rec if _detect_mi(rec) else None

# -------------------- NEGATIVE --------------------
def _mi_failure_reasons(rec: dict) -> Tuple[List[str], Set[FrozenSet[str]]]:
    reasons: List[str] = []
    formed, broken = pairs_from_reactive(rec)
    if not formed and not broken:
        return (["No σ bond changes recorded (formed/broken are empty)"], set())

    if _has_ch_transfer(formed, broken):
        reasons.append("Pattern matches proton/CH transfer (CMD/PDM-like), not migratory insertion")

    react_max = max_sigma_order_by_atom(rec, "reactant")
    prod_max  = max_sigma_order_by_atom(rec, "product")

    # witness candidates among formed X–Y (no metal)
    candidates = []
    for a,b in formed:
        if is_metal(a) or is_metal(b): 
            continue
        ra, rb = react_max.get(a,1), react_max.get(b,1)
        pa, pb = prod_max.get(a,1),  prod_max.get(b,1)
        if (ra>1 and pa==1) or (rb>1 and pb==1):
            candidates.append((a,b))
    if not candidates:
        reasons.append("No formed X–Y bond shows unsaturated→saturated change on either atom")

    # require a weakened M–R (broken or 1→0) consistent with migration
    broken_sigma, sigma_to_dat = _collect_sigma_weaken(rec)
    if not (broken_sigma or sigma_to_dat):
        reasons.append("No weakened metal–group pair detected (no M–R broken or σ→dative change)")

    # hydrometallation guard
    if any(_is_H(a) or _is_H(b) for a,b in formed):
        if not any(((_is_H(x) and is_metal(y)) or (_is_H(y) and is_metal(x))) for x,y in broken):
            reasons.append("Formed pair contains H but no broken M–H (hydrometallation inconsistency)")

    # electron neutrality (if known) for any metal adjacent to weaken pairs
    viol: List[str] = []
    metals: Set[str] = set()
    for p in (broken_sigma|sigma_to_dat):
        for t in p:
            if is_metal(t): metals.add(t)
    for M in sorted(metals):
        bE, aE = metal_before_after_electrons(rec, M)
        ed = electron_direction(bE, aE)
        if ed in ("increase","decrease"):
            viol.append(f"{M}: {bE} → {aE} ({ed})")
    if viol:
        reasons.append("Metal electron count changes (MI is typically redox-neutral): " + "; ".join(viol))

    # required tagging
    req: Set[FrozenSet[str]] = set()
    for a,b in formed:
        if not is_metal(a) and not is_metal(b):
            req.add(frozenset((a,b)))
    req |= broken_sigma | sigma_to_dat
    return reasons, req

def determine_migratory_insertion_failed_inplace(rec: dict) -> Optional[dict]:
    rec_copy = json.loads(json.dumps(rec))
    if _detect_mi(rec_copy):
        return None
    reasons, req = _mi_failure_reasons(rec)
    if not reasons:
        return None
    if req:
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
    rec["deterministic"] = {
        "reaction_class": "not_migratory_insertion",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Migratory insertion checker with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()
    process_path(Path(args.path), determine_migratory_insertion_inplace, prefix="MI")
    process_path(Path(args.path), determine_migratory_insertion_failed_inplace, prefix="MI_Failed")

if __name__ == "__main__":
    main()

