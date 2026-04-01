#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from typing import Optional, Tuple, List, Set, FrozenSet, Dict
from pathlib import Path

from rxnclass_helper import (
    elt, pairs_from_reactive, metals_in_record,
    replace_task_with_deterministic, process_path,
    metal_before_after_electrons, electron_direction,
    rebuild_bond_changes_with_tags, split_required_optional_pairs_inplace,
)
from metal_ligand.labels import is_metal

_E = {"B","C","N","O","F","Si","P","S","Cl","Se","Br","I","Ge","H"}

def _fro(a,b): return frozenset((a,b))
def _norm(a,b): return (a,b) if a<=b else (b,a)

# ---------- helpers shared with positive/negative ----------
def _bond_order_pairs(rec: dict, which: str, order_val: int, elems: Optional[Set[str]]=None) -> Set[FrozenSet[str]]:
    out: Set[FrozenSet[str]] = set()
    ig = (rec.get("input_graphs") or {}).get(which) or {}
    idx2lab = {i: str(a.get("id")) for i, a in enumerate(ig.get("atoms") or [])}
    for b in ig.get("bonds") or []:
        i, j, o = b.get("i"), b.get("j"), b.get("order")
        try:
            o = int(round(float(o)))
        except Exception:
            continue
        if o != order_val:
            continue
        u = idx2lab.get(i) if isinstance(i, int) else str(i)
        v = idx2lab.get(j) if isinstance(j, int) else str(j)
        if not u or not v:
            continue
        if elems and (elt(u) not in elems or elt(v) not in elems):
            continue
        out.add(_fro(u, v))
    return out

def _all_metals_redox_neutral(rec: dict) -> bool:
    for m in metals_in_record(rec):
        b, a = metal_before_after_electrons(rec, m)
        if b is not None and a is not None and b != a:
            return False
    return True

# -------------------- POSITIVE --------------------
def _detect_sigma_metathesis(rec: dict) -> Optional[str]:
    formed, broken = pairs_from_reactive(rec)
    fs, bs = { _norm(a,b) for a,b in formed }, { _norm(a,b) for a,b in broken }
    if not fs or not bs:
        return None

    broken_MX, broken_EH = [], []
    for a, b in broken:
        ea, eb = elt(a), elt(b)
        if is_metal(a) and not is_metal(b): broken_MX.append((a,b))
        elif is_metal(b) and not is_metal(a): broken_MX.append((b,a))
        if (ea in _E and eb=="H"): broken_EH.append((a,b))
        elif (eb in _E and ea=="H"): broken_EH.append((b,a))

    for (M, X) in broken_MX:
        for (E, H) in broken_EH:
            if _norm(M, E) not in fs: continue
            if _norm(X, H) not in fs: continue
            b, a = metal_before_after_electrons(rec, M)
            if b is not None and a is not None and b != a:
                continue
            req = {_fro(M,E), _fro(X,H), _fro(M,X), _fro(E,H)}
            rebuild_bond_changes_with_tags(rec, req)
            split_required_optional_pairs_inplace(rec, req)
            e_msg = (f" Metal electrons: {b} → {a} (no change)." if b is not None and a is not None else "")
            rec["deterministic"] = {
                "reaction_class": "metathesis",
                "subtype": "sigma",
                "explanation": f"Break {elt(M)}–{elt(X)} and {elt(E)}–H; form {elt(M)}–{elt(E)} and {elt(X)}–H; redox-neutral at metal.{e_msg}",
            }
            replace_task_with_deterministic(rec)
            return "sigma"
    return None

def _detect_olefin_metathesis(rec: dict) -> Optional[str]:
    dbl_r = _bond_order_pairs(rec, "reactant", 2, {"C"})
    dbl_p = _bond_order_pairs(rec, "product",  2, {"C"})
    removed, added = dbl_r - dbl_p, dbl_p - dbl_r
    if len(removed)==2 and len(added)==2:
        # partner swap if degree=1 for all four atoms in removed and in added; also same 4 atoms
        rem_atoms = set().union(*removed); add_atoms = set().union(*added)
        if rem_atoms == add_atoms and len(rem_atoms) == 4:
            if not _all_metals_redox_neutral(rec):
                return None
            req = set(removed) | set(added)
            rebuild_bond_changes_with_tags(rec, req)
            split_required_optional_pairs_inplace(rec, req)
            rec["deterministic"] = {
                "reaction_class": "metathesis",
                "subtype": "olefin",
                "explanation": "Olefin metathesis: partner swap among two C=C bonds (AB+CD→AD+CB) with redox-neutral metals.",
            }
            replace_task_with_deterministic(rec)
            return "olefin"
    return None

def determine_metathesis_inplace(rec: dict) -> Optional[dict]:
    if _detect_sigma_metathesis(rec): return rec
    if _detect_olefin_metathesis(rec): return rec
    return None

# -------------------- NEGATIVE --------------------
def _met_failure_reasons(rec: dict) -> Tuple[List[str], Set[FrozenSet[str]]]:
    reasons: List[str] = []
    formed, broken = pairs_from_reactive(rec)
    if not formed and not broken:
        return (["No σ bond changes recorded (formed/broken are empty)"], set())

    # σ/σ-CAM endpoint check
    fs, bs = { _norm(a,b) for a,b in formed }, { _norm(a,b) for a,b in broken }
    broken_MX = [(a,b) if is_metal(a) and not is_metal(b) else (b,a)
                 for a,b in broken if (is_metal(a) ^ is_metal(b))]
    broken_EH = []
    for a,b in broken:
        ea, eb = elt(a), elt(b)
        if (ea in _E and eb=="H"): broken_EH.append((a,b))
        if (eb in _E and ea=="H"): broken_EH.append((b,a))
    if not broken_MX:
        reasons.append("σ/σ-CAM: no M–X bond is broken")
    if not broken_EH:
        reasons.append("σ/σ-CAM: no E–H bond is broken")
    if broken_MX and broken_EH:
        # require M–E and X–H formed
        mx_eh_ok = False
        for (M,X) in broken_MX:
            for (E,H) in broken_EH:
                if _norm(M,E) in fs and _norm(X,H) in fs:
                    mx_eh_ok = True
                    break
            if mx_eh_ok: break
        if not mx_eh_ok:
            reasons.append("σ/σ-CAM: missing formation of M–E and X–H corresponding to broken M–X and E–H")

    # olefin partner-swap check
    dbl_r = _bond_order_pairs(rec, "reactant", 2, {"C"})
    dbl_p = _bond_order_pairs(rec, "product",  2, {"C"})
    removed, added = dbl_r - dbl_p, dbl_p - dbl_r
    if (removed or added) and not (len(removed)==2 and len(added)==2 and set().union(*removed)==set().union(*added)):
        reasons.append("Olefin: C=C changes do not match a clean 2-bond partner swap")

    # redox neutrality
    redox_viol = []
    for m in metals_in_record(rec):
        b, a = metal_before_after_electrons(rec, m)
        ed = electron_direction(b, a)
        if ed in ("increase", "decrease"):
            redox_viol.append(f"{m}: {b} → {a} ({ed})")
    if redox_viol:
        reasons.append("Metal electron count changes (metathesis is redox-neutral): " + "; ".join(redox_viol))

    # required tagging (collect obvious pairs)
    req: Set[FrozenSet[str]] = set()
    for (M,X) in broken_MX: req.add(_fro(M,X))
    for (E,H) in broken_EH: req.add(_fro(E,H))
    # collect plausible formed counterparts
    for a,b in formed:
        req.add(_fro(a,b))
    for p in removed | added: req.add(p)
    return reasons, req

def determine_metathesis_failed_inplace(rec: dict) -> Optional[dict]:
    rec_copy = json.loads(json.dumps(rec))
    if determine_metathesis_inplace(rec_copy):
        return None
    reasons, req = _met_failure_reasons(rec)
    if not reasons:
        return None
    if req:
        rebuild_bond_changes_with_tags(rec, req)
        split_required_optional_pairs_inplace(rec, req)
    rec["deterministic"] = {
        "reaction_class": "not_metathesis",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Metathesis checker (σ/σ-CAM + olefin) with negative outputs.")
    ap.add_argument("path", type=str)
    args = ap.parse_args()
    process_path(Path(args.path), determine_metathesis_inplace, prefix="MET")
    process_path(Path(args.path), determine_metathesis_failed_inplace, prefix="MET_Failed")

if __name__ == "__main__":
    main()

