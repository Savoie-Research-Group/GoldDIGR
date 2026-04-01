#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from typing import Dict, Set, Tuple, List, FrozenSet, Iterable
from collections import defaultdict

from rxnclass_helper import (
    elt, pairs_from_reactive, replace_task_with_deterministic, process_path,
    rebuild_bond_changes_with_tags, split_required_optional_pairs_inplace,
    metal_before_after_electrons, electron_direction,
)
from metal_ligand.labels import is_metal

# -------------------- Config --------------------
# We handle β-hetero (X in HETERO_NO_H), β-alkyl (X=C), and (optionally) β-hydride (X=H) here.
HETERO_NO_H = {"F","Cl","Br","I","O","N","S","P","Se","Si","B"}
INCLUDE_HYDRIDE = True  # set True to unify β-H inside BXA; set False to keep β-H in BHE only
HETERO = (HETERO_NO_H | {"C"}) | ({"H"} if INCLUDE_HYDRIDE else set())

ALLOW_PRODUCT_ONLY_DATIVE_MX = True   # accept product-only M–X dative (order 0) as "formed"
REQUIRE_REDOX_NEUTRAL = True          # if electron counts are present, require neutrality at the involved metal

# -------------------- Small helpers --------------------
def _bond_order_map(rec: dict, which: str) -> Dict[FrozenSet[str], int]:
    ig = (rec.get("input_graphs") or {}).get(which) or {}
    bonds = ig.get("bonds") or []
    m: Dict[FrozenSet[str], int] = {}
    for b in bonds:
        i, j = str(b.get("i")), str(b.get("j"))
        try:
            o = int(round(float(b.get("order", 0))))
        except Exception:
            o = 0
        if i is None or j is None:
            continue
        m[frozenset((i, j))] = o
    return m

def _sigma_graph(rec: dict, which: str="reactant"):
    ig = (rec.get("input_graphs") or {}).get(which) or {}
    bonds = ig.get("bonds") or []
    G = defaultdict(set)
    for b in bonds:
        i, j = str(b.get("i")), str(b.get("j"))
        try:
            o = int(round(float(b.get("order", 0))))
        except Exception:
            continue
        if o >= 1:
            G[i].add(j); G[j].add(i)
    return G

def _bond_order_increase(rec: dict, a: str, b: str, min_delta: int = 1) -> bool:
    r, p = _bond_order_map(rec, "reactant"), _bond_order_map(rec, "product")
    key = frozenset((a, b))
    before, after = r.get(key, 0), p.get(key, 0)
    return (before >= 1) and (after - before) >= min_delta

def _alpha_unsaturation_increase_via_neighbor(rec: dict, Ca: str, exclude: str) -> bool:
    r, p = _bond_order_map(rec, "reactant"), _bond_order_map(rec, "product")
    # Look for any Ca–Cγ (γ ≠ exclude) with Δorder >= 1
    for key in list(r.keys()):
        if Ca in key:
            Cg = next(iter(key - {Ca}))
            if Cg == exclude:
                continue
            before = r.get(frozenset((Ca, Cg)), 0)
            after  = p.get(frozenset((Ca, Cg)), 0)
            if before >= 1 and (after - before) >= 1:
                return True
    return False

def _reactant_m_ca_cb_x(rec: dict, allow_X: Set[str]) -> Iterable[Tuple[str,str,str,str]]:
    """Enumerate M–Cα–Cβ–X motifs in the reactant σ-graph, with X element in allow_X (no H)."""
    G = _sigma_graph(rec, "reactant")
    for M in list(G.keys()):
        if not is_metal(M):
            continue
        for Ca in G[M]:
            if elt(Ca) != "C":
                continue
            for Cb in G[Ca]:
                if Cb == M or elt(Cb) != "C":
                    continue
                for X in G[Cb]:
                    if X == Ca or elt(X) not in allow_X:
                        continue
                    yield (M, Ca, Cb, X)

def _is_H(x: str) -> bool: return elt(x) == "H"
def _same_H_transfer_present(formed: List[Tuple[str,str]], broken: List[Tuple[str,str]]) -> bool:
    """
    β-H requires: break Cβ–H and form M–H with the SAME H atom.
    """
    h_from_c: Set[str] = set()
    for a, b in broken:
        ea, eb = elt(a), elt(b)
        if ea == "C" and eb == "H":
            h_from_c.add(b)
        elif eb == "C" and ea == "H":
            h_from_c.add(a)
    for a, b in formed:
        ea, eb = elt(a), elt(b)
        if ea == "H" and is_metal(b) and a in h_from_c:
            return True
        if eb == "H" and is_metal(a) and b in h_from_c:
            return True
    return False

# -------------------- POSITIVE DETECTOR (general β-atom) --------------------
def _detect_bxa(rec: dict) -> str | None:
    formed, broken = pairs_from_reactive(rec)
    formed_set = {frozenset(p) for p in formed}
    broken_set = {frozenset(p) for p in broken}

    # ---- Branch A: β-heteroatom (X ≠ H) ----
    # Require: (M–X formed) + (Cβ–X broken) with same X; and Δorder(Cα–Cβ) ≥ 1
    # Also allow product-only dative M–X (order 0 in product, absent in reactant)
    rmap, pmap = _bond_order_map(rec, "reactant"), _bond_order_map(rec, "product")

    # Try every M–Cα–Cβ–X motif in reactant with X in the unified set (hetero, carbon, and optional H)
    for M, Ca, Cb, X in _reactant_m_ca_cb_x(rec, HETERO):
        # ---- β-hetero/hydride arm (X in hetero or X==H if enabled) ----
        if elt(X) in HETERO_NO_H or (INCLUDE_HYDRIDE and elt(X) == "H"):
            mx = frozenset((M, X))
            cbx = frozenset((Cb, X))

            mx_ok = (mx in formed_set)
            if not mx_ok and ALLOW_PRODUCT_ONLY_DATIVE_MX:
                # product-only dative: order 0 in product and absent in reactant
                if pmap.get(mx, -1) == 0 and mx not in rmap:
                    mx_ok = True

            # If X is H (β-H), require the SAME H to appear in M–H (formed) and Cβ–H (broken)
            if elt(X) == "H" and INCLUDE_HYDRIDE and not _same_H_transfer_present(formed, broken):
                mx_ok = False

            if mx_ok and (cbx in broken_set) and _bond_order_increase(rec, Ca, Cb, 1):
                subtype = "hydride" if (elt(X) == "H") else "hetero"
                if REQUIRE_REDOX_NEUTRAL:
                    bE, aE = metal_before_after_electrons(rec, M)
                    if bE is not None and aE is not None and bE != aE:
                        pass  # reject
                    else:
                        required = {mx, cbx}
                        rebuild_bond_changes_with_tags(rec, required)
                        split_required_optional_pairs_inplace(rec, required)
                        rec["deterministic"] = {
                            "reaction_class": "beta_atom_elimination",
                            "subtype": subtype,
                            "explanation": (
                                ("β-hydride: Cβ–H breaks and {M}–H forms; " if subtype=="hydride"
                                 else "β-heteroatom: Cβ–X breaks and {M}–X forms; ")
                                + "Cα–Cβ gains unsaturation; metal redox-neutral if known."
                            )

                        }
                        replace_task_with_deterministic(rec)
                        return "beta_atom_elimination"
                else:
                    required = {mx, cbx}
                    rebuild_bond_changes_with_tags(rec, required)
                    split_required_optional_pairs_inplace(rec, required)
                    rec["deterministic"] = {
                        "reaction_class": "beta_atom_elimination",
                        "subtype": ("hydride" if (elt(X)=="H" and INCLUDE_HYDRIDE) else "hetero"),
                        "explanation": (
                            "β-hydride: Cβ–H breaks and M–H forms; Cα–Cβ gains unsaturation."
                            if (elt(X)=="H" and INCLUDE_HYDRIDE)
                            else "β-heteroatom: Cβ–X breaks and M–X forms; Cα–Cβ gains unsaturation."
                        )
                    }
                    replace_task_with_deterministic(rec)
                    return "beta_atom_elimination"

    # ---- Branch B: β-alkyl (X = C) ----
    # Require: (M–Cβ formed) + (Cα–Cβ broken) + α’s other bond gains unsaturation (Δ≥1)
    for a, b in formed:
        M, Cb = (a, b) if is_metal(a) and elt(b) == "C" else ((b, a) if is_metal(b) and elt(a) == "C" else (None, None))
        if M is None:
            continue

        # find any broken pair Cα–Cβ
        for x, y in broken:
            if Cb not in (x, y):
                continue
            Ca = x if y == Cb else y
            if elt(Ca) != "C":
                continue

            # α gains unsaturation with a neighbor other than β
            if not _alpha_unsaturation_increase_via_neighbor(rec, Ca, exclude=Cb):
                continue

            if REQUIRE_REDOX_NEUTRAL:
                bE, aE = metal_before_after_electrons(rec, M)
                if bE is not None and aE is not None and bE != aE:
                    continue

            required = {frozenset((M, Cb)), frozenset((Ca, Cb))}
            rebuild_bond_changes_with_tags(rec, required)
            split_required_optional_pairs_inplace(rec, required)
            rec["deterministic"] = {
                "reaction_class": "beta_atom_elimination",
                "subtype": "alkyl",
                "explanation": (f"β-alkyl: Cα–Cβ breaks and {M}–Cβ forms; "
                                f"Cα gains unsaturation with another neighbor; metal redox-neutral if known.")
            }
            replace_task_with_deterministic(rec)
            return "beta_atom_elimination"

    return None

def determine_beta_atom_elimination_inplace(rec: dict) -> dict | None:
    return rec if _detect_bxa(rec) else None

# -------------------- NEGATIVE with reasons --------------------
def _bxa_failure_reasons(rec: dict) -> Tuple[List[str], Set[FrozenSet[str]]]:
    reasons: List[str] = []
    formed, broken = pairs_from_reactive(rec)
    if not formed and not broken:
        return (["No σ bond changes recorded (formed/broken are empty)"], set())

    formed_set = {frozenset(p) for p in formed}
    broken_set = {frozenset(p) for p in broken}
    rmap, pmap = _bond_order_map(rec, "reactant"), _bond_order_map(rec, "product")

    # H guard: if the only candidate is β-H, that belongs to BHE
    if any(_is_H(x) or _is_H(y) for x, y in formed + broken):
        # Don't blanket-reject; just note if the only viable X is H
        pass

    # Check presence of hetero/hydride pattern pieces
    hetero_any = False
    hetero_anchor_ok = False
    unsat_ok_hetero = False

    hydride_anchor_seen = False
    hydride_sameH_ok = False
    for M, Ca, Cb, X in _reactant_m_ca_cb_x(rec, HETERO):
        if elt(X) in HETERO_NO_H or (INCLUDE_HYDRIDE and elt(X) == "H"):
            hetero_any = True
            mx = frozenset((M, X)); cbx = frozenset((Cb, X))
            mx_ok = (mx in formed_set)
            if not mx_ok and ALLOW_PRODUCT_ONLY_DATIVE_MX:
                if pmap.get(mx, -1) == 0 and mx not in rmap:
                    mx_ok = True
            if mx_ok and (cbx in broken_set):
                if elt(X) == "H" and INCLUDE_HYDRIDE:
                    hydride_anchor_seen = True
                    if _same_H_transfer_present(formed, broken):
                        hydride_sameH_ok = True
                        hetero_anchor_ok = True   # treat hydride as a valid anchor for aggregation
                        if _bond_order_increase(rec, Ca, Cb, 1):
                            unsat_ok_hetero = True
                else:
                    hetero_anchor_ok = True
                    if _bond_order_increase(rec, Ca, Cb, 1):
                        unsat_ok_hetero = True

    # Check presence of alkyl pieces
    alkyl_anchor_ok = False
    unsat_ok_alkyl = False
    for a, b in formed:
        M, Cb = (a, b) if is_metal(a) and elt(b) == "C" else ((b, a) if is_metal(b) and elt(a) == "C" else (None, None))
        if M is None:
            continue
        for x, y in broken:
            if Cb not in (x, y): 
                continue
            Ca = x if y == Cb else y
            if elt(Ca) != "C":
                continue
            alkyl_anchor_ok = True
            if _alpha_unsaturation_increase_via_neighbor(rec, Ca, exclude=Cb):
                unsat_ok_alkyl = True

    # Reasons aggregation
    if not (hetero_anchor_ok or alkyl_anchor_ok):
        if INCLUDE_HYDRIDE:
            reasons.append("No matching anchors: need (M–X formed + Cβ–X broken) or (M–Cβ formed + Cα–Cβ broken); X may be H only if the same H forms M–H")
        else:
            reasons.append("No matching anchors: need (M–X formed + Cβ–X broken) or (M–Cβ formed + Cα–Cβ broken)")

    if hetero_anchor_ok and not unsat_ok_hetero:
        reasons.append("β-hetero/hydride branch: Cα–Cβ does not gain unsaturation")
    if INCLUDE_HYDRIDE and hydride_anchor_seen and not hydride_sameH_ok:
        reasons.append("β-hydride branch: Cβ–H and M–H must involve the same H atom")

    if alkyl_anchor_ok and not unsat_ok_alkyl:
        reasons.append("β-alkyl branch: Cα does not gain unsaturation with a neighbor other than Cβ")

    # Redox neutrality (if known)
    if REQUIRE_REDOX_NEUTRAL:
        # collect metals from anchors
        metals: Set[str] = set()
        for a, b in formed:
            if is_metal(a): metals.add(a)
            if is_metal(b): metals.add(b)
        redox_viol = []
        for M in sorted(metals):
            bE, aE = metal_before_after_electrons(rec, M)
            ed = electron_direction(bE, aE)
            if ed in ("increase", "decrease"):
                redox_viol.append(f"{M}: {bE} → {aE} ({ed})")
        if redox_viol:
            reasons.append("Metal electron count changes (β-atom eliminations are typically redox-neutral): " + "; ".join(redox_viol))

    # Tag the obvious required pairs for context
    req: Set[FrozenSet[str]] = set()
    for a, b in formed: 
        if is_metal(a) or is_metal(b):
            req.add(frozenset((a, b)))
    for a, b in broken:
        req.add(frozenset((a, b)))

    return reasons, req

def determine_beta_atom_elimination_failed_inplace(rec: dict) -> dict | None:
    rec_copy = json.loads(json.dumps(rec))
    if _detect_bxa(rec_copy):
        return None
    reasons, req = _bxa_failure_reasons(rec)
    if not reasons:
        reasons = ["Record does not satisfy β-atom elimination (hetero or alkyl), but no specific rule fired."]
    if req:
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
    rec["deterministic"] = {
        "reaction_class": "not_beta_atom_elimination",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="β-atom elimination (hetero + alkyl + optional hydride) with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()
    process_path(Path(args.path), determine_beta_atom_elimination_inplace, prefix="BXA")
    process_path(Path(args.path), determine_beta_atom_elimination_failed_inplace, prefix="BXA_Failed")

if __name__ == "__main__":
    main()

