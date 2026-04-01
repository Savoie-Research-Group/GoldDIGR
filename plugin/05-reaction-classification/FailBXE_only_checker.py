#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from typing import Optional, Set, Tuple, List, Dict, FrozenSet, Iterable
from pathlib import Path
from collections import defaultdict

from rxnclass_helper import (
    elt, pairs_from_reactive, replace_task_with_deterministic, process_path,
    rebuild_bond_changes_with_tags, split_required_optional_pairs_inplace,
    metal_before_after_electrons, electron_direction,
)
from metal_ligand.labels import is_metal

# Allowed heteroatoms X (no H here)
HETERO = {"F","Cl","Br","I","O","N","S","P","Se","C"}  # you included C in your draft
ALLOW_DATIVE_ALPHA = True  # allow M–Cα dative for late metals
LATE = {"Au","Pt","Pd","Ag","Cu"}

def _why_not_bxe(rec: dict) -> List[str]:
    reasons, _ = _bxe_failure_reasons(rec)
    return reasons or ["no explicit reason collected (unexpected)"]

def _bond_order_map(rec: dict, which: str) -> dict[FrozenSet[str], int]:
    bonds = (((rec.get("input_graphs") or {}).get(which) or {}).get("bonds") or [])
    m: Dict[FrozenSet[str], int] = {}
    for b in bonds:
        i, j = b.get("i"), b.get("j")
        try:
            o = int(round(float(b.get("order", 0))))
        except Exception:
            o = 0
        if i is None or j is None: 
            continue
        m[frozenset((str(i), str(j)))] = o
    return m

#def _bond_order_increase(rec: dict, a: str, b: str) -> bool:
def _bond_order_increase(rec: dict, a: str, b: str, min_delta: int = 1) -> bool:
    r, p = _bond_order_map(rec, "reactant"), _bond_order_map(rec, "product")
    key = frozenset((a, b))
    before, after = r.get(key, 0), p.get(key, 0)
    #return (before >= 1) and (after >= 2) and (after - before) >= 1
    # require any increase by ≥ min_delta; do NOT insist on “becomes 2”
    return (before >= 1) and (after - before) >= min_delta

def _sigma_graph(rec: dict, which: str="reactant"):
    atoms = (((rec.get("input_graphs") or {}).get(which) or {}).get("atoms") or [])
    bonds = (((rec.get("input_graphs") or {}).get(which) or {}).get("bonds") or [])
    G = defaultdict(set)
    for b in bonds:
        try:
            i, j = b.get("i"), b.get("j")
            order = int(round(float(b.get("order", 0))))
        except Exception:
            continue
        if order >= 1:
            G[str(i)].add(str(j)); G[str(j)].add(str(i))
    return G

def _reactant_m_ca_cb_x(rec: dict, hetero: Set[str]) -> Iterable[Tuple[str,str,str,str]]:
    G = _sigma_graph(rec, "reactant")
    # allow M–Cα as dative only for late metals if requested
    def neighbors(ep):
        return G.get(ep, set())
    # metals present (by label presence in graph keys)
    candidates = set(G.keys())
    for M in list(candidates):
        if not is_metal(M): 
            continue
        for Ca in neighbors(M):
            if elt(Ca) != "C":
                continue
            for Cb in neighbors(Ca):
                if Cb == M or elt(Cb) != "C":
                    continue
                for X in neighbors(Cb):
                    if X == Ca or elt(X) not in hetero:
                        continue
                    yield (M, Ca, Cb, X)

def _has_same_H_CtoM_transfer(formed: List[Tuple[str,str]], broken: List[Tuple[str,str]]) -> bool:
    # same as your draft guard — if C–H broken and M–H formed for SAME H → that’s BHE (β-H)
    h_left_carbon = set()
    for a, b in broken:
        ea, eb = elt(a), elt(b)
        if ea == "C" and eb == "H": h_left_carbon.add(b)
        if eb == "C" and ea == "H": h_left_carbon.add(a)
    for a, b in formed:
        ea, eb = elt(a), elt(b)
        if ea == "H" and is_metal(b) and a in h_left_carbon: return True
        if eb == "H" and is_metal(a) and b in h_left_carbon: return True
    return False

# -------------------- POSITIVE --------------------
def _detect_bxe(rec: dict) -> Optional[str]:
    formed, broken = pairs_from_reactive(rec)
    if _has_same_H_CtoM_transfer(formed, broken):
        return None  # let BHE checker own β-H

    formed_set = {frozenset(p) for p in formed}
    broken_set = {frozenset(p) for p in broken}

    # Need a reactant motif M–Cα–Cβ–X and SAME X moves: Cβ–X broken & M–X formed
    for M, Ca, Cb, X in _reactant_m_ca_cb_x(rec, HETERO):
        if frozenset((Cb, X)) not in broken_set:
            continue
        '''
        if (frozenset((M, X)) not in formed_set):
            # try dative formed in product only when absent in reactant (order 0)
            pmap, rmap = _bond_order_map(rec, "product"), _bond_order_map(rec, "reactant")
            key = frozenset((M, X))
            if not (pmap.get(key, -1) == 0 and key not in rmap):
                continue
        '''
        # Accept either σ formation OR product-only dative (order 0 in product, absent in reactant)
        if (frozenset((M, X)) not in formed_set):
            pmap, rmap = _bond_order_map(rec, "product"), _bond_order_map(rec, "reactant")
            key = frozenset((M, X))
            if not (pmap.get(key, -1) == 0 and key not in rmap):
                continue

        # α–β becomes more unsaturated
        #if not _bond_order_increase_to_double(rec, Ca, Cb):
        if not _bond_order_increase(rec, Ca, Cb, min_delta=1):
            continue

        # metal redox neutrality (prefer)
        b, a = metal_before_after_electrons(rec, M)
        if b is not None and a is not None and b != a:
            continue

        required = {frozenset((Cb, X)), frozenset((M, X))}
        rebuild_bond_changes_with_tags(rec, required)
        split_required_optional_pairs_inplace(rec, required)
        e_msg = (f" Metal electrons: {b} → {a} (no change)." if (b is not None and a is not None) else "")
        rec["deterministic"] = {
            "reaction_class": "beta_heteroatom_elimination",
            "explanation": (f"Reactant motif {M}–Cα({Ca})–Cβ({Cb})–{X}; Cβ–{X} breaks and {M}–{X} forms; "
                            f"Cα–Cβ gains unsaturation.{e_msg}")
        }
        replace_task_with_deterministic(rec)
        return "beta_heteroatom_elimination"
    return None

def determine_beta_heteroatom_elimination_inplace(rec: dict) -> Optional[dict]:
    return rec if _detect_bxe(rec) else None

# -------------------- NEGATIVE --------------------
def _bxe_failure_reasons(rec: dict) -> Tuple[List[str], Set[FrozenSet[str]]]:
    reasons: List[str] = []
    formed, broken = pairs_from_reactive(rec)
    if not formed and not broken:
        return (["No σ bond changes recorded (formed/broken are empty)"], set())

    # β-H elimination exclusion
    if _has_same_H_CtoM_transfer(formed, broken):
        reasons.append("Pattern matches β-hydride elimination (same H transfers C→M), not β-heteroatom elimination")

    # need some Cβ–X broken with X∈HETERO
    has_CX_broken = False
    for a, b in broken:
        if "C" in (elt(a), elt(b)):
            X = b if elt(a) == "C" else a
            if elt(X) in HETERO:
                has_CX_broken = True
                break
    if not has_CX_broken:
        reasons.append("No C–X σ bond (X in allowed hetero set) is broken")

    # need M–X formed (or dative formed only in product) for same X
    formed_set = {frozenset(p) for p in formed}
    broken_set = {frozenset(p) for p in broken}
    # coarse pairing test
    pair_ok = False
    # build a bag of candidate X from broken C–X
    cand_X: Set[str] = set()
    for a, b in broken:
        if "C" in (elt(a), elt(b)):
            X = b if elt(a) == "C" else a
            if elt(X) in HETERO:
                cand_X.add(X)
    rmap, pmap = _bond_order_map(rec, "reactant"), _bond_order_map(rec, "product")
    for X in cand_X:
        # any metal M such that M–X formed or dative formed in product only?
        if any((frozenset((M, X)) in formed_set) for M in list(cand_X) + list({a for a, _ in formed} | {b for _, b in formed})):
            pair_ok = True
            break
        for M in {a for a, _ in formed} | {b for _, b in formed}:
            if not is_metal(M): 
                continue
            k = frozenset((M, X))
            if pmap.get(k, -1) == 0 and k not in rmap:
                pair_ok = True
                break
    if cand_X and not pair_ok:
        reasons.append("No matching M–X formation for the broken Cβ–X (same X)")

    # α–β unsaturation increase
    motif_any = False
    unsat_any = False
    for M, Ca, Cb, X in _reactant_m_ca_cb_x(rec, HETERO):
        motif_any = True
        #if _bond_order_increase_to_double(rec, Ca, Cb):
        if _bond_order_increase(rec, Ca, Cb, min_delta=1):
            unsat_any = True
            break
    if motif_any and not unsat_any:
        reasons.append("Cα–Cβ bond does not gain unsaturation")

    # redox neutrality (if e data present)
    # if any metal with known electrons changes, mention it
    e_viol: List[str] = []
    igm = set()
    for M, *_ in _reactant_m_ca_cb_x(rec, HETERO):
        igm.add(M)
    for M in sorted(igm):
        b, a = metal_before_after_electrons(rec, M)
        ed = electron_direction(b, a)
        if ed in ("increase", "decrease"):
            e_viol.append(f"{M}: {b} → {a} ({ed})")
    if e_viol:
        reasons.append("Metal electron count changes (β-heteroatom elimination is typically electron-neutral): " + "; ".join(e_viol))

    # Tag all potentially relevant pairs
    req: Set[FrozenSet[str]] = set()
    for a, b in broken:
        if "C" in (elt(a), elt(b)):
            X = b if elt(a) == "C" else a
            if elt(X) in HETERO:
                req.add(frozenset((a, b)))
    for a, b in formed:
        if is_metal(a) and elt(b) in HETERO: req.add(frozenset((a, b)))
        if is_metal(b) and elt(a) in HETERO: req.add(frozenset((a, b)))
    return reasons, req

def determine_beta_heteroatom_elimination_failed_inplace(rec: dict) -> Optional[dict]:
    rec_copy = json.loads(json.dumps(rec))
    if _detect_bxe(rec_copy):
        return None
    reasons, req = _bxe_failure_reasons(rec)
    if not reasons:
        reasons = ["Record does not satisfy BXE, but no specific rule fired. "
                   "Likely causes: (i) no Cβ–X broken matching the M–X change; "
                   "(ii) α–β unsaturation test not met; (iii) M–X only broken (reverse anchor) "
                   "without accepting that pattern."]
    if req:
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
    rec["deterministic"] = {
        "reaction_class": "not_beta_heteroatom_elimination",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="β-Heteroatom elimination checker with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()
    process_path(Path(args.path), determine_beta_heteroatom_elimination_inplace, prefix="BXE")
    process_path(Path(args.path), determine_beta_heteroatom_elimination_failed_inplace, prefix="BXE_Failed")

if __name__ == "__main__":
    main()

