#!/usr/bin/env python3
"""
bel_parser.py — Parse Bond-Electron List (BEL) files into the JSON-style
record dicts consumed by the *_only_checker.py classification scripts.

BEL format (from BEL-rules.txt):
  Schema: ID-neighbors bond:1=default,^=dative,==dbl,#=trpl /e=valence(default 0).
  Geometry: Metal{trans-pairs}/geom or Metal[ligands]/geom. A~B=trans.

A BEL string has two sections separated by a newline (real or literal '\\n'):
  1) Geometry header (optional):  Metal{trans}/geom-Metal[ligs]/geom  or  @Metal[ligs]/geom-...
  2) Bond list:  space-separated tokens of the form  CENTER-NBR1,NBR2^,NBR3=,...

Neighbor modifiers:
  (none) → σ single bond  (order 1)
  ^      → dative          (order 0)
  =      → double bond     (order 2)
  #      → triple bond     (order 3)

Non-bonding electrons can be annotated in two ways:
  - On a neighbor token:      Cs26-Ni25^/7   → Ni25 has 7 e
  - On the center before dash: Co0/6-C3^,...  → Co0 has 6 e
  - As a trailing /e on the last neighbor:    Co0-...,O42^/4/6  → O42 has 4 e, Co0 has 6 e
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from metal_ligand.labels import is_metal

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_EL_RE = re.compile(r"^([A-Za-z]+)")


def _element_of(atom_id: str) -> str:
    """Extract element symbol from an atom label like 'Ni25' → 'Ni'."""
    m = _EL_RE.match(atom_id)
    if not m:
        return atom_id
    raw = m.group(1)
    return raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()


def _sort_key(atom_id: str):
    """Sort atoms by numeric index within element groups."""
    m = re.search(r"(\d+)$", atom_id)
    return (int(m.group(1)) if m else 0, atom_id)


def _try_int(s: str) -> Optional[int]:
    try:
        return int(s)
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# BEL correction — fix misplaced metal electron annotations
# ---------------------------------------------------------------------------

def correct_bel(bel_text: str) -> str:
    """
    Fix a known BEL generator bug: when a metal center has non-bonding
    electrons, the generator appends ``/e`` at the end of the neighbor chain
    instead of placing it on the center atom before the dash.

    Bad:      ``Fe0-C55^,N1^/2,N2^/2,N3^/2,N4^/2,H56/6``
    Correct:  ``Fe0/6-C55^,N1^/2,N2^/2,N3^/2,N4^/2,H56``

    Bad:      ``Fe0-N1^/2,N2^/2,N3^/2,N4^/2/8``
    Correct:  ``Fe0/8-N1^/2,N2^/2,N3^/2,N4^/2``

    Only metal center tokens are corrected.  Non-metal tokens (where ``/e``
    on a neighbor is legitimately that neighbor's electrons) are left alone.

    Disambiguation rules for the last neighbor's trailing ``/e``:

    * **Double-slash** (``NBR/e1/e2``): always unambiguous — neighbor keeps
      ``/e1``, last ``/e2`` is relocated to the metal center.
    * **Single-slash** (``NBR/e``): depends on the neighbor's element:

      - H or C never carry non-bonding electrons → ``/e`` belongs to the
        metal center (relocate).
      - N, O, S, P, halogens, etc. can carry lone pairs → ``/e`` belongs
        to the neighbor (leave as-is, metal has 0 electrons).
    """
    # Elements that never carry non-bonding electrons in BEL encoding
    _NO_LONE_PAIR = {"H", "C"}

    text = bel_text.strip()

    # ---- Split header from bond section ----
    if "\n" in text:
        header, bond_section = text.split("\n", 1)
        sep = "\n"
    elif "\\n" in text:
        header, bond_section = text.split("\\n", 1)
        sep = "\\n"
    else:
        header = None
        bond_section = text
        sep = None

    corrected_tokens = []
    for token in bond_section.strip().split():
        if "-" not in token:
            # Standalone atom — leave as-is
            corrected_tokens.append(token)
            continue

        dash = token.index("-")
        center_raw = token[:dash]
        rest = token[dash + 1:]

        # ---- Only fix metal centers without /e before dash ----
        if "/" in center_raw:
            # Center already has /e — no fix needed
            corrected_tokens.append(token)
            continue

        if not is_metal(center_raw):
            # Non-metal center — leave as-is
            corrected_tokens.append(token)
            continue

        # ---- Metal center with no /e — check last neighbor ----
        neighbors = rest.split(",")
        last_nb = neighbors[-1]
        slash_parts = last_nb.split("/")

        if len(slash_parts) < 2:
            # No /e on last neighbor — nothing to relocate
            corrected_tokens.append(token)
            continue

        if len(slash_parts) >= 3:
            # Double+ slash: e.g. "O42^/4/6" or "N4^/2/8"
            # Always unambiguous: last /e belongs to the metal center
            center_e = slash_parts[-1]
            neighbors[-1] = "/".join(slash_parts[:-1])
            new_token = f"{center_raw}/{center_e}-{','.join(neighbors)}"
            corrected_tokens.append(new_token)
            continue

        # Single slash: e.g. "H56/6" or "O42^/4"
        # Determine the neighbor's element to decide ownership
        nb_base = slash_parts[0]
        # Strip bond-type suffix (^, =, #) to get atom id
        nb_id = nb_base.rstrip("^=#")
        nb_el = _element_of(nb_id)

        if nb_el in _NO_LONE_PAIR:
            # H or C cannot have lone pairs → /e belongs to the metal
            center_e = slash_parts[1]
            neighbors[-1] = nb_base
            new_token = f"{center_raw}/{center_e}-{','.join(neighbors)}"
            corrected_tokens.append(new_token)
        else:
            # Heteroatom (N, O, S, ...) → /e legitimately belongs to
            # the neighbor; metal has 0 electrons.  Leave as-is.
            corrected_tokens.append(token)

    corrected_bond_section = " ".join(corrected_tokens)
    if header is not None:
        return f"{header}{sep}{corrected_bond_section}"
    else:
        return corrected_bond_section


# ---------------------------------------------------------------------------
# Core BEL parser
# ---------------------------------------------------------------------------

def _split_header(text: str) -> str:
    """
    Strip geometry header, return only the bond-list section.
    Handles both real newlines (from JSON) and literal two-char '\\n' (from .txt files).
    """
    # Real newline: take everything after first line
    if "\n" in text:
        lines = text.split("\n")
        return " ".join(lines[1:]).strip()
    # Literal \n (two chars: backslash + n)
    if "\\n" in text:
        return text.split("\\n", 1)[1].strip()
    # No header at all
    return text.strip()


def parse_bel(bel_text: str) -> Tuple[List[dict], List[dict]]:
    """
    Parse a single BEL string into *atoms* and *bonds* lists that are
    structurally identical to the JSON ``{"atoms": [...], "bonds": [...]}``
    format used by the checker scripts.

    Parameters
    ----------
    bel_text : str
        One BEL line (may contain a geometry header before a newline).

    Returns
    -------
    atoms : list[dict]
        ``[{"id": "Ni25", "el": "Ni", "e": 7}, ...]``
    bonds : list[dict]
        ``[{"i": "C0", "j": "C1", "order": 1}, ...]``
    """
    bond_section = _split_header(bel_text.strip())

    atoms: Dict[str, dict] = {}          # id → {id, el, e}
    bonds: Dict[frozenset, int] = {}     # frozenset(i,j) → order

    def _ensure(aid: str):
        if aid not in atoms:
            atoms[aid] = {"id": aid, "el": _element_of(aid), "e": 0}

    for token in bond_section.split():
        # ---- Standalone atom (no bonds, e.g. detached "H44") ----
        if "-" not in token:
            parts = token.split("/")
            atom_id = parts[0].strip()
            if atom_id and _EL_RE.match(atom_id):
                _ensure(atom_id)
                if len(parts) > 1:
                    v = _try_int(parts[1])
                    if v is not None:
                        atoms[atom_id]["e"] = v
            continue

        # ---- CENTER-NBR1,NBR2,... token ----
        dash = token.index("-")
        center_raw = token[:dash]
        rest = token[dash + 1:]

        # Center may carry /e before dash:  Co0/6-C3^,...
        center_e: Optional[int] = None
        if "/" in center_raw:
            cparts = center_raw.split("/")
            center = cparts[0]
            center_e = _try_int(cparts[1])
        else:
            center = center_raw

        _ensure(center)
        if center_e is not None and center_e > 0:
            atoms[center]["e"] = center_e

        # ---- Parse comma-separated neighbor tokens ----
        for nb_raw in rest.split(","):
            nb_raw = nb_raw.strip()
            if not nb_raw:
                continue

            # Split on '/' to separate atom+modifier from electron annotations.
            # Patterns:  "C7^/2"    → [C7^, 2]           neighbor e=2
            #            "O42^/4/6" → [O42^, 4, 6]       neighbor e=4, center e=6
            #            "C3^"      → [C3^]              no electrons
            slash_parts = nb_raw.split("/")
            main_part = slash_parts[0]

            neighbor_e: Optional[int] = None
            trailing_center_e: Optional[int] = None

            if len(slash_parts) == 2:
                neighbor_e = _try_int(slash_parts[1])
            elif len(slash_parts) >= 3:
                neighbor_e = _try_int(slash_parts[1])
                trailing_center_e = _try_int(slash_parts[-1])

            # ---- Bond-type suffix ----
            if main_part.endswith("^"):
                atom_id, order = main_part[:-1], 0   # dative
            elif main_part.endswith("="):
                atom_id, order = main_part[:-1], 2   # double
            elif main_part.endswith("#"):
                atom_id, order = main_part[:-1], 3   # triple
            else:
                atom_id, order = main_part, 1         # σ single

            _ensure(atom_id)
            if neighbor_e is not None and neighbor_e > 0:
                atoms[atom_id]["e"] = neighbor_e
            if trailing_center_e is not None and trailing_center_e > 0:
                atoms[center]["e"] = trailing_center_e

            pair = frozenset((center, atom_id))
            if pair not in bonds:
                bonds[pair] = order

    # ---- Convert to sorted lists ----
    atoms_list = [atoms[k] for k in sorted(atoms, key=_sort_key)]
    bonds_list = []
    for pair, order in bonds.items():
        a, b = sorted(pair)
        bonds_list.append({"i": a, "j": b, "order": order})
    bonds_list.sort(key=lambda b: (_sort_key(b["i"]), _sort_key(b["j"])))

    return atoms_list, bonds_list


# ---------------------------------------------------------------------------
# Reactive-pair computation (diff two bond graphs)
# ---------------------------------------------------------------------------

def compute_reactive_pairs(
    r_bonds: List[dict],
    p_bonds: List[dict],
) -> List[list]:
    """
    Diff reactant and product bond lists to produce ``target.reactive_pairs``.

    Returns a list of ``[atom_a, atom_b, "before -> after"]`` for every bond
    whose order changed.  Absent bonds are encoded as ``"NaN"``.
    """
    def _to_map(blist):
        m: Dict[frozenset, int] = {}
        for b in blist:
            m[frozenset((b["i"], b["j"]))] = b["order"]
        return m

    r_map = _to_map(r_bonds)
    p_map = _to_map(p_bonds)

    reactive: List[list] = []
    for pair in sorted(r_map.keys() | p_map.keys(),
                       key=lambda p: tuple(sorted(p))):
        r_o = r_map.get(pair)
        p_o = p_map.get(pair)
        if r_o == p_o:
            continue
        a, b = sorted(pair)
        before = str(r_o) if r_o is not None else "NaN"
        after  = str(p_o) if p_o is not None else "NaN"
        reactive.append([a, b, f"{before} -> {after}"])

    return reactive


# ---------------------------------------------------------------------------
# Canonical hashing
# ---------------------------------------------------------------------------

def _canonical_graph_str(atoms: List[dict], bonds: List[dict]) -> str:
    """
    Produce a deterministic string for a parsed graph (atoms + bonds),
    independent of BEL neighbor ordering or formatting.
    """
    atom_tuples = tuple((a["id"], a["el"], a["e"]) for a in atoms)
    bond_tuples = tuple((b["i"], b["j"], b["order"]) for b in bonds)
    return repr((atom_tuples, bond_tuples))


def reaction_hash(
    r_atoms: List[dict], r_bonds: List[dict],
    p_atoms: List[dict], p_bonds: List[dict],
) -> str:
    """
    SHA-256 hex digest of the canonical (reactant, product) graph pair.
    Deterministic regardless of BEL formatting or neighbor ordering.
    """
    combined = (
        "R:" + _canonical_graph_str(r_atoms, r_bonds) +
        "|P:" + _canonical_graph_str(p_atoms, p_bonds)
    )
    return hashlib.sha256(combined.encode()).hexdigest()


# ---------------------------------------------------------------------------
# σ-crossing filter for bond_changes strings
# ---------------------------------------------------------------------------

def normalize_bond_changes(bond_changes: str) -> str:
    """
    Canonicalize a bond-changes string by:
      1) Sorting the two atom IDs within each entry (``H29-Fe1:1>n`` → ``Fe1-H29:1>n``)
      2) Sorting entries alphabetically (``Fe1-H29:1>n;B0-C2:n>1`` → ``B0-C2:n>1;Fe1-H29:1>n``)

    Empty/whitespace input returns ``""``.
    """
    if not bond_changes or not bond_changes.strip():
        return ""
    normalized = []
    for entry in bond_changes.split(";"):
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue
        pair_str, change = entry.rsplit(":", 1)
        if "-" in pair_str:
            a, b = pair_str.split("-", 1)
            pair_str = f"{a}-{b}" if a <= b else f"{b}-{a}"
        normalized.append(f"{pair_str}:{change}")
    return ";".join(sorted(normalized))


def reaction_hash_from_bel(reactant_bel: str, bond_changes_sigma: str) -> str:
    """
    SHA-256 hex digest of (canonical reactant graph, normalized σ bond changes).

    This identifies a reaction by *what molecule reacts* and *which σ bonds
    change*, independent of BEL formatting, neighbor ordering, or entry order
    in the bond-changes string.
    """
    r_atoms, r_bonds = parse_bel(reactant_bel)
    bc_norm = normalize_bond_changes(bond_changes_sigma)
    combined = "R:" + _canonical_graph_str(r_atoms, r_bonds) + "|BC:" + bc_norm
    return hashlib.sha256(combined.encode()).hexdigest()

def _sigma_class(tok: str) -> int:
    """0 for {n, NaN, 0}, else 1."""
    tok = tok.strip().lower()
    if tok in ("n", "nan", "0"):
        return 0
    try:
        return 0 if int(tok) <= 0 else 1
    except ValueError:
        return 0


def filter_sigma_changes(bond_changes_raw: str) -> str:
    """
    Keep only σ-crossing entries from a bond_changes string.

    σ-crossing: one side in {n, 0} and the other in {1, 2, 3}.
    Exclude: n↔0  and  ≥1↔≥1.
    """
    if not bond_changes_raw:
        return ""
    kept = []
    for entry in bond_changes_raw.split(";"):
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue
        _, change = entry.rsplit(":", 1)
        if ">" not in change:
            continue
        before_tok, after_tok = change.split(">", 1)
        bv, av = _sigma_class(before_tok), _sigma_class(after_tok)
        # Keep only when exactly one side is 0 (σ-crossing)
        if (bv == 0) ^ (av == 0):
            kept.append(entry)
    return ";".join(kept)


# ---------------------------------------------------------------------------
# Build the full ``rec`` dict consumed by checker scripts
# ---------------------------------------------------------------------------

def build_rec_from_bel(
    reactant_bel: str,
    product_bel: str,
    source: str = "bel_input",
) -> dict:
    """
    Construct a record dict from reactant and product BEL strings.

    The returned dict has exactly the structure the *_only_checker.py
    scripts expect::

        {
          "source": ...,
          "input_graphs": {
              "reactant": {"atoms": [...], "bonds": [...]},
              "product":  {"atoms": [...], "bonds": [...]},
          },
          "target": {
              "reactive_pairs": [...],
          },
          "metal_electron_change": { ... },
        }
    """
    r_atoms, r_bonds = parse_bel(reactant_bel)
    p_atoms, p_bonds = parse_bel(product_bel)

    reactive_pairs = compute_reactive_pairs(r_bonds, p_bonds)

    # ---- metal_electron_change (mirrors JSON convention) ----
    r_e = {a["id"]: a["e"] for a in r_atoms}
    p_e = {a["id"]: a["e"] for a in p_atoms}
    mec: Dict[str, dict] = {}
    for aid in sorted(set(r_e) | set(p_e)):
        if not is_metal(aid):
            continue
        be, ae = r_e.get(aid, 0), p_e.get(aid, 0)
        if be != 0 or ae != 0:
            mec[aid] = {"before": be, "after": ae}

    rec: dict = {
        "source": source,
        "input_graphs": {
            "reactant": {"atoms": r_atoms, "bonds": r_bonds},
            "product":  {"atoms": p_atoms, "bonds": p_bonds},
        },
        "target": {
            "reactive_pairs": reactive_pairs,
        },
    }
    if mec:
        rec["metal_electron_change"] = mec

    return rec


def build_rec_from_bel_files(
    reactant_path: Path,
    product_path: Path,
    source: Optional[str] = None,
) -> dict:
    """Convenience: build a rec dict from BEL file paths."""
    r_text = Path(reactant_path).read_text().strip()
    p_text = Path(product_path).read_text().strip()
    src = source or str(Path(reactant_path).resolve().parent.name)
    return build_rec_from_bel(r_text, p_text, source=src)

