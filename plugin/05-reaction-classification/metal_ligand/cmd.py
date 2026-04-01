from __future__ import annotations
from typing import List
from .labels import elem_of, is_metal, normalize_pair, build_index

HETERO_BASES = {"O","N"}

def metals_hint_from_raw(raw: dict | None) -> List[str]:
    metals: List[str] = []
    if not raw: return metals
    mec = raw.get("metal_electron_change", {})
    if isinstance(mec, dict):
        metals.extend(mec.keys())
    lig = raw.get("ligand_changes", {})
    if isinstance(lig, dict):
        for k in lig.keys():
            if elem_of(k) not in HETERO_BASES and k not in metals:
                metals.append(k)
    return metals

def detect_cmd_forward(broken, formed, include_sbm=False, raw=None, optional_metal_bonds=True):
    broken_set = {normalize_pair(*p) for p in broken}
    formed_set = {normalize_pair(*p) for p in formed}
    broken_idx = build_index(list(broken_set))
    formed_idx = build_index(list(formed_set))
    # forbid formed M–H unless SBM allowed
    if any((is_metal(a) and elem_of(b)=="H") or (is_metal(b) and elem_of(a)=="H") for a,b in formed_set):
        if not include_sbm: return None
    # need a broken C–H with the same H now bound to X∈{O,N}
    for a,b in broken_set:
        if {elem_of(a), elem_of(b)} != {"C","H"}: continue
        c_label = a if elem_of(a)=="C" else b
        h_label = a if elem_of(a)=="H" else b
        h_new = formed_idx.get(h_label, set())
        if len(h_new)!=1: continue
        x_label = next(iter(h_new))
        if elem_of(x_label) not in HETERO_BASES: continue
        # identify metal
        cand = [n for n in broken_idx.get(x_label,set()) if is_metal(n)]
        for n in (n for n in formed_idx.get(c_label,set()) if is_metal(n)):
            if n not in cand: cand.append(n)
        if not cand and raw:
            hint = metals_hint_from_raw(raw)
            if len(hint)==1 and is_metal(hint[0]): cand=[hint[0]]
        if not cand: continue
        for m in cand:
            has_mx_broken = (min(m,x_label), max(m,x_label)) in broken_set
            has_mc_formed = (min(m,c_label), max(m,c_label)) in formed_set
            if not include_sbm and (min(m,h_label), max(m,h_label)) in formed_set: continue
            if (has_mx_broken and has_mc_formed) or optional_metal_bonds:
                return {
                    "metal": m, "carbon": c_label, "hydrogen": h_label, "base_atom": x_label,
                    "explanation": "CMD C–H activation: C–H broken, X–H formed; intramolecular base."
                }
    return None

