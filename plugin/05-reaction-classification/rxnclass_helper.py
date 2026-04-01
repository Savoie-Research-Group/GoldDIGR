#!/usr/bin/env python3
from __future__ import annotations
import json
import re
from pathlib import Path
from typing import Callable, Iterable, Optional, Dict, Tuple

from metal_ligand.labels import is_metal
from metal_ligand.labels import normalize_pair

# -------------------- Lightweight shared helpers --------------------

_RX_EL = re.compile(r"[A-Za-z]+")

def elt(label: str) -> str:
    m = _RX_EL.match(label or "")
    return m.group(0) if m else ""

def is_H(x: str) -> bool: return elt(x) == "H"
def is_C(x: str) -> bool: return elt(x) == "C"

def sigma_val(s: object) -> int:
    """
    Convert a 'before/after' token into σ semantics:
      - NaN or 0  → 0 (no σ-bond; allow '0' to mean dative)
      - >=1       → 1 (σ-bond present)
    """
    if s is None:
        return 0
    t = str(s).strip()
    if t == "NaN":
        return 0
    try:
        v = int(t)
    except Exception:
        try:
            v = int(round(float(t)))
        except Exception:
            v = 0
    return 0 if v <= 0 else 1

def split_change(change: str) -> Tuple[str, str]:
    if not isinstance(change, str) or "->" not in change:
        return ("NaN", "NaN")
    before, after = [t.strip() for t in change.split("->", 1)]
    return before, after

def replace_task_with_deterministic(rec: dict) -> None:
    """
    Ensure there is NO top-level 'task' in outputs.
    If 'deterministic' missing but 'task' exists, migrate content; then drop 'task'.
    """
    if "task" in rec:
        if "deterministic" not in rec and isinstance(rec["task"], dict):
            rec["deterministic"] = rec["task"]
        rec.pop("task", None)

# -------------------- Output naming helpers --------------------

REACTION_CLASS_TO_PREFIX: Dict[str, str] = {
    "C_H_activation": "CH",
    "oxidative_addition": "OA",
    "reductive_elimination": "RE",
    "migratory_insertion": "MI",
    "beta_hydride_elimination": "BHE",
    "transmetalation": "TM",
}

def prefix_for_class(reaction_class: str) -> str:
    """
    Map 'reaction_class' → short prefix. Adds '-rev' suffix for *_reversed.
    Unknown classes → 'RXN'.
    """
    if not isinstance(reaction_class, str) or not reaction_class:
        return "RXN"
    base = reaction_class.replace("_reversed", "")
    pref = REACTION_CLASS_TO_PREFIX.get(base, "RXN")
    return f"{pref}-rev" if reaction_class.endswith("_reversed") else pref

def write_with_prefix(path: Path, rec: dict, prefix: str) -> Path:
    """
    Write JSON next to 'path' as '<PREFIX>-<original-filename>'.
    Returns the output Path.
    """
    out = path.with_name(f"{prefix}-{path.name}")
    out.write_text(json.dumps(rec, indent=2))
    return out

# -------------------- File/dir dispatch --------------------

DEFAULT_INPUT_FILES = (
    "deterministic_forward.clean.json",
    "deterministic_reverse.clean.json",
)

DEFAULT_BEL_FILES = {
    "reactant": "BEL-reactant.txt",
    "product":  "BEL-product.txt",
}

def _try_bel_directory(
    p: Path,
    determine_inplace_func: Callable[[dict], Optional[dict]],
    prefix: str,
    bel_files: dict = DEFAULT_BEL_FILES,
) -> bool:
    """
    If 'p' contains a BEL reactant/product pair, parse them into a rec,
    run the classifier, and write the result.  Returns True if an output
    was written.
    """
    r_bel = p / bel_files["reactant"]
    p_bel = p / bel_files["product"]
    if not (r_bel.exists() and p_bel.exists()):
        return False

    from bel_parser import build_rec_from_bel_files
    rec = build_rec_from_bel_files(r_bel, p_bel)
    rec = determine_inplace_func(rec)
    if rec is None:
        return False

    # Use a synthetic reference path so write_with_prefix produces a
    # reasonable filename:  <dir>/<PREFIX>-bel_result.json
    ref = p / "bel_result.json"
    out = write_with_prefix(ref, rec, prefix)
    print(f"✓ Wrote {out}")
    return True


def process_path(
    path: Path | str,
    determine_inplace_func: Callable[[dict], Optional[dict]],
    prefix: str,
    candidates: Iterable[str] = DEFAULT_INPUT_FILES,
) -> None:
    """
    Generic dispatcher:
      - If 'path' is a JSON file: load → determine_inplace_func → if not None, write_with_prefix(prefix)
      - If 'path' is a dir: try JSON 'candidates' first, then fall back to
        BEL-reactant.txt / BEL-product.txt if present.
      - If determine_inplace_func returns None → NO OUTPUT is written.
    """
    p = Path(path).resolve()

    # Single file
    if p.is_file():
        if p.suffix.lower() != ".json":
            raise ValueError(f"Unsupported file: {p}")
        rec = json.loads(p.read_text())
        rec = determine_inplace_func(rec)
        if rec is not None:
            out = write_with_prefix(p, rec, prefix)
            print(f"✓ Wrote {out}")
        return

    # Directory
    if p.is_dir():
        wrote_any = False

        # --- Phase 1: try JSON candidates ---
        for name in candidates:
            f = p / name
            if not f.exists():
                continue
            rec = json.loads(f.read_text())
            rec = determine_inplace_func(rec)
            if rec is not None:
                out = write_with_prefix(f, rec, prefix)
                print(f"✓ Wrote {out}")
                wrote_any = True

        # --- Phase 2: try BEL pair (only if no JSON was found) ---
        if not wrote_any:
            wrote_any = _try_bel_directory(p, determine_inplace_func, prefix)

        if not wrote_any:
            print("No outputs written (no matching inputs or no classified records).")
        return

    raise FileNotFoundError(f"Path not found: {p}")

from typing import List, Tuple, Set

def pairs_from_reactive(rec: dict) -> Tuple[List[Tuple[str,str]], List[Tuple[str,str]]]:
    """
    Extract σ-formed / σ-broken pairs from target.reactive_pairs.
    Return (formed_pairs, broken_pairs) as normalized label tuples.
    """
    formed, broken = [], []
    for trip in ((rec.get("target") or {}).get("reactive_pairs") or []):
        if not isinstance(trip, (list, tuple)) or len(trip) < 2:
            continue
        a, b = str(trip[0]), str(trip[1])
        before, after = split_change(str(trip[2]) if len(trip) > 2 else "")
        bsv, asv = sigma_val(before), sigma_val(after)
        if bsv < 1 and asv >= 1:
            formed.append(normalize_pair(a, b))
        elif bsv >= 1 and asv < 1:
            broken.append(normalize_pair(a, b))
    return formed, broken

def metals_in_record(rec: dict) -> Set[str]:
    """
    Collect atom labels that look like metals from input_graphs.reactant/product.atoms.
    """
    out: Set[str] = set()
    for side in ("reactant", "product"):
        atoms = (((rec.get("input_graphs") or {}).get(side) or {}).get("atoms")) or []
        for a in atoms:
            lab = str(a.get("id", ""))
            if lab and is_metal(lab):
                out.add(lab)
    return out

# --- Metal electron utilities ---

def _atom_e_map(rec: dict, side: str) -> dict[str, float]:
    """Map label -> electron count for one side ('reactant' or 'product')."""
    atoms = (((rec.get("input_graphs") or {}).get(side) or {}).get("atoms")) or []
    out = {}
    for a in atoms:
        lab = str(a.get("id", ""))
        if not lab:
            continue
        # Prefer 'e' if present; fall back to 'electrons' if your older JSON used that key
        if "e" in a:
            out[lab] = float(a["e"])
        elif "electrons" in a:
            out[lab] = float(a["electrons"])
    return out

def metal_before_after_electrons(rec: dict, m_label: str) -> tuple[float | None, float | None]:
    """
    Return (before, after) electron counts for metal m_label.
    Priority 1: top-level 'metal_electron_change'[m_label] = {before, after}
    Fallback:   input_graphs.reactant/product.atoms[].e matched by 'id' == m_label
    """
    mec = rec.get("metal_electron_change") or {}
    if m_label in mec:
        b = mec[m_label].get("before")
        a = mec[m_label].get("after")
        try:
            return (float(b), float(a))
        except Exception:
            pass  # fall through to atom maps

    react = _atom_e_map(rec, "reactant")
    prod  = _atom_e_map(rec, "product")
    return (react.get(m_label), prod.get(m_label))

def electron_direction(before: float | None, after: float | None) -> str | None:
    """
    'increase' if after>before, 'decrease' if after<before, 'no_change' if equal,
    None if either is missing.
    """
    if before is None or after is None:
        return None
    if after > before:
        return "increase"
    if after < before:
        return "decrease"
    return "no_change"

def rebuild_bond_changes_with_tags(rec: dict, required_pairs: Set[frozenset[str]]) -> None:
    """
    Rebuild top-level rec['bond_changes'] from target.reactive_pairs,
    tagging entries as 'required' vs 'optional'. Only σ-formed/σ-broken included.
    """
    formed_out, broken_out = [], []
    for trip in ((rec.get("target") or {}).get("reactive_pairs") or []):
        if not isinstance(trip, (list, tuple)) or len(trip) < 3:
            continue
        a, b, change = str(trip[0]), str(trip[1]), str(trip[2])
        before, after = split_change(change)
        bsv, asv = sigma_val(before), sigma_val(after)
        tag = "required" if frozenset((a, b)) in required_pairs else "optional"
        if bsv < 1 and asv >= 1:
            formed_out.append([a, b, f"{before} -> {after}", tag])
        elif bsv >= 1 and asv < 1:
            broken_out.append([a, b, f"{before} -> {after}", tag])
    rec["bond_changes"] = {"formed": formed_out, "broken": broken_out}

def split_required_optional_pairs_inplace(rec: dict, required_pairs: Set[frozenset[str]]) -> None:
    """
    Keep REQUIRED pairs in target.reactive_pairs, move the rest to target.optional_pairs.
    """
    tgt = rec.setdefault("target", {})
    rps = list(tgt.get("reactive_pairs", []))
    req_list, opt_list = [], []
    for trip in rps:
        if not isinstance(trip, (list, tuple)) or len(trip) < 2:
            continue
        a, b = str(trip[0]), str(trip[1])
        if frozenset((a, b)) in required_pairs:
            req_list.append(trip)
        else:
            opt_list.append(trip)
    tgt["reactive_pairs"] = req_list
    if opt_list:
        tgt["optional_pairs"] = opt_list
    else:
        tgt.pop("optional_pairs", None)

def _idx_label_maps(rec: dict, side: str) -> Tuple[Dict[int,str], Dict[str,int]]:
    atoms = (((rec.get("input_graphs") or {}).get(side) or {}).get("atoms")) or []
    idx2lab: Dict[int,str] = {}
    lab2idx: Dict[str,int] = {}
    for i, a in enumerate(atoms):
        lab = str(a.get("id",""))
        idx2lab[i] = lab
        if lab:
            lab2idx[lab] = i
    return idx2lab, lab2idx

def max_sigma_order_by_atom(rec: dict, side: str) -> Optional[Dict[str,int]]:
    bonds = (((rec.get("input_graphs") or {}).get(side) or {}).get("bonds")) or []
    idx2lab, _ = _idx_label_maps(rec, side)
    if not bonds or not idx2lab:
        return None
    out: Dict[str,int] = {}
    for b in bonds:
        i, j = int(b.get("i", -1)), int(b.get("j", -1))
        order = int(b.get("order", 0))
        if order < 1:  # σ-only
            continue
        ai, aj = idx2lab.get(i), idx2lab.get(j)
        if ai: out[ai] = max(out.get(ai, 0), order)
        if aj: out[aj] = max(out.get(aj, 0), order)
    return out

# --- Max σ order by atom from bond list (labels or indices supported) ---
from typing import Dict, Tuple

def _idx2lab(rec: dict, side: str) -> Dict[int, str]:
    atoms = (((rec.get("input_graphs") or {}).get(side) or {}).get("atoms")) or []
    out: Dict[int, str] = {}
    for i, a in enumerate(atoms):
        lab = str(a.get("id", ""))
        out[i] = lab
    return out

def _label_for(endpoint, idx2lab: Dict[int, str]) -> str | None:
    # Bonds can store endpoints as integers (indices) or strings (labels).
    if isinstance(endpoint, int):
        return idx2lab.get(endpoint)
    try:
        # if it’s a string, use it as-is (e.g., "Au0")
        return str(endpoint)
    except Exception:
        return None

def max_sigma_order_by_atom(rec: dict, side: str) -> Dict[str, int]:
    """
    Build {label -> max σ-bond order} from input_graphs[side].bonds.
    σ means order >= 1. Dative (0) is ignored. Works for i/j as labels or indices.
    """
    bonds = (((rec.get("input_graphs") or {}).get(side) or {}).get("bonds")) or []
    idx2lab = _idx2lab(rec, side)
    out: Dict[str, int] = {}
    for b in bonds:
        i_ep = b.get("i")
        j_ep = b.get("j")
        try:
            order = int(b.get("order", 0))
        except Exception:
            continue
        if order < 1:
            continue  # ignore non-σ (e.g., dative)
        ai = _label_for(i_ep, idx2lab)
        aj = _label_for(j_ep, idx2lab)
        if ai:
            out[ai] = max(out.get(ai, 0), order)
        if aj:
            out[aj] = max(out.get(aj, 0), order)
    return out

def _get_e(atoms, atom_id):
    for a in atoms:
        if str(a.get("id")) == str(atom_id):
            return a.get("e")
    return None

