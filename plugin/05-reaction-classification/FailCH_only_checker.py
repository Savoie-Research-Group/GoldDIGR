#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
from typing import Dict, List, Tuple, Set, Optional
from pathlib import Path

from rxnclass_helper import (
    elt, is_H, is_C, sigma_val, split_change,
    replace_task_with_deterministic, process_path,
    is_metal,
    pairs_from_reactive, metals_in_record,
    metal_before_after_electrons, electron_direction,
    split_required_optional_pairs_inplace,  # reuse helper for splitting pairs
)

# =============================================================================
# -------------------- CH-activation (CMD) logic (shared) ---------------------
# =============================================================================

def _collect_by_H(rec: dict):
    """
    For each H atom, bucket reactive_pairs into:
      HC_broken, HX_broken, HC_formed, HX_formed  (σ semantics)
    Each entry: (a, b, before, after)
    HX_* lists only count X that is a NON-METAL (to exclude β-H elimination).
    """
    rps = (rec.get("target", {}) or {}).get("reactive_pairs", [])
    HC_broken_by_H: Dict[str, List[Tuple[str, str, str, str]]] = {}
    HX_broken_by_H: Dict[str, List[Tuple[str, str, str, str]]] = {}
    HC_formed_by_H: Dict[str, List[Tuple[str, str, str, str]]] = {}
    HX_formed_by_H: Dict[str, List[Tuple[str, str, str, str]]] = {}

    for trip in rps:
        if not isinstance(trip, (list, tuple)) or len(trip) < 3:
            continue
        a, b, change = trip[0], trip[1], str(trip[2])
        before, after = split_change(change)
        bsv, asv = sigma_val(before), sigma_val(after)

        # only consider σ-broken / σ-formed (ignore σ→σ order changes)
        if not (is_H(a) or is_H(b)):
            continue
        H = a if is_H(a) else b
        X = b if H == a else a

        # broken
        if bsv >= 1 and asv < 1:
            if is_C(X):
                HC_broken_by_H.setdefault(H, []).append((a, b, before, after))
            elif not is_metal(X):  # only non-metal acceptors considered for CMD
                HX_broken_by_H.setdefault(H, []).append((a, b, before, after))
            # else: X is a metal → β-H elimination pattern; ignore for CMD

        # formed
        elif bsv < 1 and asv >= 1:
            if is_C(X):
                HC_formed_by_H.setdefault(H, []).append((a, b, before, after))
            elif not is_metal(X):  # only non-metal acceptors considered for CMD
                HX_formed_by_H.setdefault(H, []).append((a, b, before, after))
            # else: X is a metal → H–M formation; exclude from CMD

    return HC_broken_by_H, HX_broken_by_H, HC_formed_by_H, HX_formed_by_H


def _decide_orientation(rec: dict) -> tuple[Optional[str], Optional[str]]:
    """
    Decide CMD orientation on this single record:
      - 'C_H_activation' if (H–C broken) & (H–X formed) exist for same H (X non-metal)
      - 'C_H_activation_reversed' if (H–X broken) & (H–C formed) exist for same H
    Picks the H with the strongest minimal count; tie → forward.
    """
    HC_broken_by_H, HX_broken_by_H, HC_formed_by_H, HX_formed_by_H = _collect_by_H(rec)

    f_scores = {H: min(len(HC_broken_by_H[H]), len(HX_formed_by_H[H]))
                for H in set(HC_broken_by_H) & set(HX_formed_by_H)}
    r_scores = {H: min(len(HC_formed_by_H[H]), len(HX_broken_by_H[H]))
                for H in set(HC_formed_by_H) & set(HX_broken_by_H)}

    fmax = max(f_scores.values()) if f_scores else 0
    rmax = max(r_scores.values()) if r_scores else 0
    if fmax == 0 and rmax == 0:
        return None, None
    if fmax >= rmax:
        H = max(f_scores, key=f_scores.get)
        return "C_H_activation", H
    else:
        H = max(r_scores, key=r_scores.get)
        return "C_H_activation_reversed", H


def _rebuild_bond_changes_with_tags(rec: dict, required_pairs: Set[frozenset[str]]) -> None:
    """
    Build/replace top-level rec['bond_changes'] with σ-formed/σ-broken entries only:
      [A, B, "before -> after", "required|optional"]
    """
    formed_out, broken_out = [], []
    for trip in (rec.get("target", {}) or {}).get("reactive_pairs", []):
        if not isinstance(trip, (list, tuple)) or len(trip) < 3:
            continue
        a, b, change = trip[0], trip[1], str(trip[2])
        before, after = split_change(change)
        bsv, asv = sigma_val(before), sigma_val(after)
        tag = "required" if frozenset((a, b)) in required_pairs else "optional"
        if bsv < 1 and asv >= 1:
            formed_out.append([a, b, f"{before} -> {after}", tag])
        elif bsv >= 1 and asv < 1:
            broken_out.append([a, b, f"{before} -> {after}", tag])
    rec["bond_changes"] = {"formed": formed_out, "broken": broken_out}


def _reactant_has_CH_sigma(rec: dict) -> bool:
    """
    Check if reactant graph contains any C–H σ bond (order >= 1).
    Supports bonds with endpoints as indices or labels.
    """
    side = "reactant"
    atoms = (((rec.get("input_graphs") or {}).get(side) or {}).get("atoms")) or []
    bonds = (((rec.get("input_graphs") or {}).get(side) or {}).get("bonds")) or []
    idx2lab = {i: str(a.get("id", "")) for i, a in enumerate(atoms)}
    for b in bonds:
        i_ep = b.get("i")
        j_ep = b.get("j")
        try:
            order = int(b.get("order", 0))
        except Exception:
            continue
        if order < 1:
            continue
        li = idx2lab.get(i_ep) if isinstance(i_ep, int) else str(i_ep)
        lj = idx2lab.get(j_ep) if isinstance(j_ep, int) else str(j_ep)
        if (li and lj) and ((is_H(li) and is_C(lj)) or (is_H(lj) and is_C(li))):
            return True
    return False

# =============================================================================
# --------------------- POSITIVE: forward CMD only ----------------------------
# =============================================================================

def determine_ch_activation_inplace(rec: dict) -> Optional[dict]:
    """
    Single-record CH activation (CMD) checker.
    STRICT: only forward CMD ('C_H_activation') is accepted as positive.
      - If not forward CMD → return None (caller writes no positive output).
      - If forward CMD:
          * set rec['deterministic'] with reaction_class + explanation
          * remove/replace any top-level 'task'
          * rebuild 'bond_changes' with required/optional tags
          * keep only required pairs in target.reactive_pairs and move others → target.optional_pairs
    """
    label, chosen_H = _decide_orientation(rec)
    if label != "C_H_activation":
        return None  # strict: do not accept reverse as positive

    # Minimal required set: one (H–C) broken + one (H–X non-metal) formed for the same H
    HC_broken_by_H, _, _, HX_formed_by_H = _collect_by_H(rec)
    req_pairs: Set[frozenset[str]] = {
        frozenset(HC_broken_by_H[chosen_H][0][0:2]),
        frozenset(HX_formed_by_H[chosen_H][0][0:2]),
    }

    _rebuild_bond_changes_with_tags(rec, req_pairs)

    # Split required vs optional in target
    tgt = rec.setdefault("target", {})
    all_rps = list(tgt.get("reactive_pairs", []))
    req_list, opt_list = [], []
    for trip in all_rps:
        if not isinstance(trip, (list, tuple)) or len(trip) < 2:
            continue
        a, b = trip[0], trip[1]
        if frozenset((a, b)) in req_pairs:
            req_list.append(trip)
        else:
            opt_list.append(trip)
    tgt["reactive_pairs"] = req_list
    if opt_list:
        tgt["optional_pairs"] = opt_list
    else:
        tgt.pop("optional_pairs", None)

    expl = f"Detected CMD (forward): H {chosen_H} has H–C σ-bond broken and H–X σ-bond formed."
    rec["deterministic"] = {"reaction_class": "C_H_activation", "explanation": expl}
    replace_task_with_deterministic(rec)
    return rec

# =============================================================================
# --------------------- NEGATIVE: rule-by-rule failure ------------------------
# =============================================================================

def _ch_failure_reasons(rec: dict) -> Tuple[List[str], Set[frozenset[str]]]:
    """
    Evaluate all CH rules independently, collect ALL violated reasons (no early returns).
    Returns (reasons, required_pairs_for_tagging).
    """
    reasons: List[str] = []

    # Orientation (to detect explicit reverse CMD)
    label, Hpicked = _decide_orientation(rec)

    # Bucket by H to test subset rules
    HC_broken_by_H, HX_broken_by_H, HC_formed_by_H, HX_formed_by_H = _collect_by_H(rec)

    # --- Rule 0: reactant must contain at least one C–H σ bond
    if not _reactant_has_CH_sigma(rec):
        reasons.append("No C–H σ bond exists on the reactant side")

    # --- Rule 1: need an H–C σ BOND BROKEN (some H)
    if not HC_broken_by_H:
        reasons.append("No H–C σ bond is broken in reactive_pairs")

    # --- Rule 2: need an H–X σ BOND FORMED with non-metal X and X ≠ C
    if not HX_formed_by_H:
        reasons.append("No H–X σ bond (X non-metal, X≠C) is formed in reactive_pairs")

    # --- Rule 3: both events must occur on the SAME hydrogen
    if HC_broken_by_H and HX_formed_by_H:
        if not (set(HC_broken_by_H.keys()) & set(HX_formed_by_H.keys())):
            reasons.append("H–C broken and H–X formed occur on different hydrogens (must involve the same H)")

    # --- Rule 4: exclude β-H elimination signature (H–M σ formed)
    formed_pairs, _ = pairs_from_reactive(rec)  # formed, broken (σ-only) from target.reactive_pairs
    if any((is_H(a) and is_metal(b)) or (is_H(b) and is_metal(a)) for a, b in formed_pairs):
        reasons.append("H–M σ bond is formed (acceptor is a metal), indicative of β-hydride elimination, not CMD")

    # --- Rule 5: metal electron count should not change for CMD
    cand_M = metals_in_record(rec)
    elec_msgs: List[str] = []
    for m in sorted(cand_M):
        before_e, after_e = metal_before_after_electrons(rec, m)
        edir = electron_direction(before_e, after_e)  # None | increase | decrease | no_change
        if edir in ("increase", "decrease"):
            tag = "reduction" if edir == "increase" else "oxidation"
            elec_msgs.append(f"{m}: {before_e} → {after_e} ({edir}; {tag})")
    if elec_msgs:
        reasons.append("Metal center electron count changes (CMD should not change metal electrons): " + "; ".join(elec_msgs))

    # --- Rule 6: reverse CMD is NOT C–H activation (strict request)
    if label == "C_H_activation_reversed":
        reasons.append(
            f"Reverse CMD detected on H {Hpicked}: H–X σ bond is broken and H–C σ bond is formed; "
            "this is not C–H activation (treat as negative)"
        )

    # For tagging/inspection, mark all relevant σ pairs (if any) as 'required'
    required_pairs: Set[frozenset[str]] = set()
    for d in (HC_broken_by_H, HX_formed_by_H, HX_broken_by_H, HC_formed_by_H):
        for entries in d.values():
            for a, b, _, _ in entries:
                required_pairs.add(frozenset((a, b)))

    return reasons, required_pairs


def determine_ch_activation_failed_inplace(rec: dict) -> Optional[dict]:
    """
    Emit a *negative* CH record if forward CMD was NOT detected, including a list of
    all violated rules in 'explanation'. Never emits a failed record when the
    positive detector would classify as forward CMD.
    """
    # Guard: if positive forward CMD would pass, do NOT emit failure
    rec_copy = json.loads(json.dumps(rec))
    if determine_ch_activation_inplace(rec_copy) is not None:
        return None

    # Collect all violated rules
    reasons, required_pairs = _ch_failure_reasons(rec)

    if not reasons:
        # If no specific reason was collected (should be rare), avoid writing noise.
        return None

    if required_pairs:
        _rebuild_bond_changes_with_tags(rec, required_pairs)
        split_required_optional_pairs_inplace(rec, required_pairs)

    rec["deterministic"] = {
        "reaction_class": "not_c_h_activation",
        "explanation": "; ".join(reasons),
    }
    replace_task_with_deterministic(rec)
    return rec

# =============================================================================
# --------------------------------- CLI ---------------------------------------
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="C–H activation (CMD) checker with negative outputs. "
                    "Writes CH-*.json for forward CMD and CH_Failed-*.json otherwise."
    )
    ap.add_argument("path", type=str, help="Path to deterministic_*.clean.json OR a directory containing them.")
    args = ap.parse_args()

    # 1) Positive: forward CMD only → CH-*.json
    process_path(Path(args.path), determine_ch_activation_inplace, prefix="CH")

    # 2) Negative: anything not passing forward CMD (including reverse CMD) → CH_Failed-*.json
    process_path(Path(args.path), determine_ch_activation_failed_inplace, prefix="CH_Failed")


if __name__ == "__main__":
    main()

