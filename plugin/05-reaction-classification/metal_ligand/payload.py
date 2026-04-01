from __future__ import annotations
from .matrices import graph_from_section

def make_schema_comments(direction_note: str | None = None) -> dict:
    return {
        "_about": "Deterministic-style payload (classification may be added by downstream tools).",
        "_direction": direction_note or "",
        "deterministic": {
            "reaction_class": "Deterministic label for this direction.",
            "explanation": "Short rationale."
        },
        "yarp_graph": {
            "reactant_product_note": "Graphs aligned with this direction.",
            "elements": "Labels with running indices (Pd0, C1, ...).",
            "electrons": "Per-atom electrons (BE diagonal).",
            "bond_list": {"legend": ["atom_A", "atom_B", "bond_order", "bond_type"],
                          "bond_type_rule": "bond_type = 'dative' if bond_order == 0 else 'sigma"}
        }
    }

def build_output_record(base_rec: dict,
                        rxn_class: str,
                        explanation: str,
                        adj_sections: dict | None,
                        be_sections: dict | None,
                        add_schema_comments: bool = True,
                        direction_note: str = "") -> dict:
    out = dict(base_rec)
    out["deterministic"] = {"reaction_class": rxn_class, "explanation": explanation}
    graphs = {}
    if adj_sections or be_sections:
        r_adj = (adj_sections or {}).get("reactant")
        p_adj = (adj_sections or {}).get("product")
        r_be  = (be_sections  or {}).get("reactant")
        p_be  = (be_sections  or {}).get("product")
        r_graph = graph_from_section(r_adj, r_be)
        p_graph = graph_from_section(p_adj, p_be)
        if r_graph: graphs["reactant"] = r_graph
        if p_graph: graphs["product"]  = p_graph
    if graphs:
        out["yarp_graph"] = graphs
    if add_schema_comments:
        out["_schema_comments"] = make_schema_comments(direction_note)
    return out
'''
def _graph_clean_from_section(section_adj: dict | None, section_be: dict | None) -> dict:
    if not section_adj or not section_be:
        return {}
    labels = section_adj.get("labels") or section_be.get("labels") or []
    adj    = section_adj.get("matrix", [])
    be     = section_be.get("matrix",  [])
    if not labels or not adj or not be:
        return {}
    n = min(len(labels), len(adj), len(be))
    atoms = []
    bonds = []
    for i in range(n):
        el = str(labels[i])
        j = len(el)
        while j > 0 and el[j-1].isdigit():
            j -= 1
        el_sym = el[:j] if j > 0 else el
        e_diag = 0
        if i < len(be) and i < len(be[i]):
            try:
                e_diag = int(round(float(be[i][i])))
            except Exception:
                e_diag = 0
        atoms.append({"id": i, "el": el_sym, "e": e_diag})
    for i in range(n):
        for j in range(i+1, n):
            try:
                aij = int(round(float(adj[i][j])))
            except Exception:
                aij = 0
            if aij != 0:
                try:
                    bo = int(round(float(be[i][j])))
                except Exception:
                    bo = 0
                bonds.append({"i": i, "j": j, "order": bo})
    return {"atoms": atoms, "bonds": bonds}

def build_clean_payload(bond_changes_idx: dict,
                        adj_sections: dict | None,
                        be_sections: dict | None,
                        direction_note: str = "") -> dict:
    out = {
        "schema_version": 3,
        "direction": direction_note,
        "input_graphs": {},
        "bond_changes": {
            "formed": bond_changes_idx.get("formed", []),
            "broken": bond_changes_idx.get("broken", []),
        },
        "deterministic": {
            "reaction_class": "Unprocessed",
            "explanation": "Unprocessed"
        }
    }
    if adj_sections or be_sections:
        r_adj = (adj_sections or {}).get("reactant")
        p_adj = (adj_sections or {}).get("product")
        r_be  = (be_sections  or {}).get("reactant")
        p_be  = (be_sections  or {}).get("product")
        r_graph = _graph_clean_from_section(r_adj, r_be)
        p_graph = _graph_clean_from_section(p_adj, p_be)
        if r_graph: out["input_graphs"]["reactant"] = r_graph
        if p_graph: out["input_graphs"]["product"]  = p_graph
    return out
'''
def _graph_clean_from_section(section_adj: dict | None, section_be: dict | None) -> dict:
    """Return {"atoms":[{id,el,e}], "bonds":[{i,j,order}]} with *label-string* ids (e.g., 'C0')."""
    if not section_adj or not section_be:
        return {}
    labels = (section_adj.get("labels") or section_be.get("labels") or [])
    adj    = section_adj.get("matrix", [])
    be     = section_be.get("matrix",  [])
    if not labels or not adj or not be:
        return {}

    n = min(len(labels), len(adj), len(be))
    atoms = []
    bonds = []

    # Atoms: id uses the label string directly (e.g. "C0"), el extracted, e from diagonal
    for i in range(n):
        lab = str(labels[i])  # e.g. "C0", "Au0"
        # element = leading letters only
        j = len(lab)
        while j > 0 and lab[j-1].isdigit():
            j -= 1
        el_sym = lab[:j] if j > 0 else lab
        e_diag = 0
        if i < len(be) and i < len(be[i]):
            try:
                e_diag = int(round(float(be[i][i])))
            except Exception:
                e_diag = 0
        atoms.append({"id": lab, "el": el_sym, "e": e_diag})

    # Bonds: i/j are label-strings, order from BE off-diagonal
    for i in range(n):
        for j in range(i+1, n):
            try:
                aij = int(round(float(adj[i][j])))
            except Exception:
                aij = 0
            if aij != 0:
                try:
                    bo = int(round(float(be[i][j])))
                except Exception:
                    bo = 0
                bonds.append({"i": labels[i], "j": labels[j], "order": bo})

    return {"atoms": atoms, "bonds": bonds}


def build_clean_payload(
    source: str,
    reaction_class: str,
    adj_sections: dict | None,
    be_sections: dict | None,
    direction_note: str = ""
) -> dict:
    """
    Build the *new* deterministic_*.clean.json:
      - top-level: source, task:{reaction_class}, input_graphs{reactant/product}, target{...}, schema_version
      - ids are label-strings like 'C0' everywhere
      - target.reactive_pairs carry 'before -> after' with 'NaN' for no bond
      - target.bond_order_changes is the subset where adjacency stayed bonded (no break/form)
    """
    out = {
        "source": source,
        "task": {"reaction_class": reaction_class},
        "input_graphs": {},
        "target": {"reactive_atoms": [], "reactive_pairs": [], "bond_order_changes": []},
        "schema_version": 3,
        "direction": direction_note,
    }

    if not (adj_sections or be_sections):
        return out

    # Pull reactant/product labelled sections
    r_adj = (adj_sections or {}).get("reactant");  r_be = (be_sections or {}).get("reactant")
    p_adj = (adj_sections or {}).get("product");   p_be = (be_sections or {}).get("product")
    if not (r_adj and r_be and p_adj and p_be):
        return out

    # Build reactant/product graphs with label-string ids
    r_graph = _graph_clean_from_section(r_adj, r_be)
    p_graph = _graph_clean_from_section(p_adj, p_be)
    if r_graph: out["input_graphs"]["reactant"] = r_graph
    if p_graph: out["input_graphs"]["product"]  = p_graph

    # Compute reactive pairs and bond-order changes
    r_labels = r_adj.get("labels", []); p_labels = p_adj.get("labels", [])
    # assume labels align between R and P (as produced by analyzer)
    labels = r_labels if r_labels else p_labels
    R_adj, R_be = r_adj.get("matrix", []), r_be.get("matrix", [])
    P_adj, P_be = p_adj.get("matrix", []), p_be.get("matrix", [])

    '''
    reactive_atoms = set()
    reactive_pairs = []
    bond_order_changes = []

    n = min(len(labels), len(R_adj), len(R_be), len(P_adj), len(P_be))
    for i in range(n):
        for j in range(i+1, n):
            before = _state(R_adj[i][j], R_be[i][j])
            after  = _state(P_adj[i][j], P_be[i][j])
            if before != after:
                ai, aj = labels[i], labels[j]
                reactive_pairs.append([ai, aj, f"{before} -> {after}"])
                reactive_atoms.update([ai, aj])
                # bond-order-only change: both bonded in R and P
                if before != "NaN" and after != "NaN":
                    bond_order_changes.append([ai, aj, f"{before} -> {after}"])

    out["target"]["reactive_atoms"] = sorted(reactive_atoms)
    out["target"]["reactive_pairs"] = reactive_pairs
    out["target"]["bond_order_changes"] = bond_order_changes
    '''
    reactive_atoms = set()
    reactive_pairs = []
    bond_order_changes = []

    n = min(len(labels), len(R_adj), len(R_be), len(P_adj), len(P_be))

    def _parse_state(adj_val, be_val):
        # returns (state_str, sigma_order_int) where state_str is "NaN" or "0"/"1"/"2"/...
        try:
            a = int(round(float(adj_val)))
        except Exception:
            a = 0
        if a == 0:
            return "NaN", 0  # no bond
        try:
            bo = int(round(float(be_val)))
        except Exception:
            bo = 0
        return str(bo), bo  # bonded; bo==0 means dative, bo>=1 means sigma

    for i in range(n):
        for j in range(i+1, n):
            before_str, before_bo = _parse_state(R_adj[i][j], R_be[i][j])
            after_str,  after_bo  = _parse_state(P_adj[i][j], P_be[i][j])

            if before_str != after_str:
                ai, aj = labels[i], labels[j]
                change_txt = f"{before_str} -> {after_str}"
                # only document sigma bond change
                if ((before_bo < 1) and (after_bo >= 1)) or ((before_bo >= 1) and (after_bo < 1)):
                    reactive_pairs.append([ai, aj, change_txt])
                    reactive_atoms.update([ai, aj])

                # NEW rule: bond_order_changes only when σ exists on both sides (>=1) and changes
                if (before_bo >= 1) and (after_bo >= 1) and (before_bo != after_bo):
                    bond_order_changes.append([ai, aj, change_txt])

    out["target"]["reactive_atoms"] = sorted(reactive_atoms)
    out["target"]["reactive_pairs"] = reactive_pairs
    out["target"]["bond_order_changes"] = bond_order_changes
    return out

