import json
from typing import Dict, List
from metal_ligand.labels import is_metal, elem_of, build_index  # uses your helpers

def check_18e(json_path: str) -> Dict[str, List[Dict]]:
    """
    Return per-side electron counts around each metal center.

    Output shape:
      {
        "reactant": [{"metal": "Ir11", "neighbors": [...], "electron_sum": 18, "is_18e": True}, ...],
        "product":  [{"metal": "Ir11", "neighbors": [...], "electron_sum": 18, "is_18e": True}, ...]
      }
    """
    with open(json_path, "r") as fh:
        rec = json.load(fh)

    out: Dict[str, List[Dict]] = {"reactant": [], "product": []}
    ig = rec.get("input_graphs") or {}

    for side in ("reactant", "product"):
        g = ig.get(side) or {}
        atoms = g.get("atoms") or []
        bonds = g.get("bonds") or []

        # Map atom id -> (element, electrons)
        e_by_id = {a["id"]: int(a.get("e") or 0) for a in atoms}
        el_by_id = {a["id"]: (a.get("el") or elem_of(a["id"])) for a in atoms}

        # Build adjacency (ignore bond order)
        pairs = [(b.get("i"), b.get("j")) for b in bonds if b.get("i") and b.get("j")]
        idx = build_index(pairs)  # id -> set(neighbors)

        # Find metal centers and compute 18e counts
        for atom_id, el in el_by_id.items():
            if not is_metal(el):
                continue
            neighbors = sorted(idx.get(atom_id, set()))
            metal_e = e_by_id.get(atom_id, 0)
            total_e = e_by_id.get(atom_id, 0) + sum(e_by_id.get(n, 0) for n in neighbors)
            out[side].append({
                "metal": atom_id,
                "metal_e": metal_e,
                "neighbors": neighbors,
                "electron_sum": total_e,
                "is_18e": (total_e == 18),
            })

    return out

