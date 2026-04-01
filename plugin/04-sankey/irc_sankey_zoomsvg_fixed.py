#!/usr/bin/env python3
"""
irc_sankey.py (fixed with SVG zoom + transport-aggregate switch)
===============================================================

Builds two Sankey diagrams from an IRC:
  1) TRANSPORT (atom → atom electron ownership flow R_i)
  2) TWO-RAIL (atom ⇄ bond-storage (i-j) ⇄ atom)

Adds:
- `--export-svg` normal+zoomed SVGs (uses python-kaleido if available)
- Always-on HTML modebar button “Download SVG (zoomed)”
- `--transport-aggregate` toggle to show a single "Other transport" donut
  instead of detailed per-node "Other from/to …" links
- Taller export canvas and tight margins

"""
from __future__ import annotations
import argparse
import numpy as np
from pathlib import Path
import json, re, csv
from typing import List, Dict, Tuple, Optional

# --------------------------- helpers ---------------------------

def _read_text_number(path: Path) -> int:
    txt = Path(path).read_text().strip()
    m = re.search(r"(\d+)", txt)
    if not m:
        raise ValueError(f"Could not parse integer from {path}")
    return int(m.group(1))


def _load_labelled_matrix(path: Path):
    """Load a square matrix that may include row/column labels and metadata.
    Returns (M, labels) where labels is a list of strings (or None)."""
    import pandas as pd
    df = pd.read_csv(path, index_col=0, low_memory=False)

    # Drop obvious metadata/string columns
    bad_cols = [c for c in df.columns if isinstance(c, str) and re.search(r'(source|path|file|meta|json)', c, re.I)]
    if bad_cols:
        df = df.drop(columns=bad_cols, errors="ignore")

    num = df.apply(pd.to_numeric, errors="coerce")
    common = [c for c in num.columns if c in num.index]
    if len(common) >= 2:
        sub = num.loc[common, common]
    else:
        row_mask = num.notna().any(axis=1)
        col_mask = num.notna().any(axis=0)
        sub = num.loc[row_mask, col_mask]
        if sub.shape[0] != sub.shape[1]:
            n = min(sub.shape[0], sub.shape[1])
            sub = sub.iloc[:n, :n]

    if sub.shape[0] != sub.shape[1]:
        raise ValueError(f"{path}: could not extract a square numeric matrix")
    sub = sub.fillna(0.0)

    labels = [str(x) for x in sub.index.tolist()]
    M = sub.values.astype(float)
    return M, labels


def _auto_detect_V_from_BEM(B: np.ndarray) -> np.ndarray:
    """Extract bond order V from BEM matrix B.
    Assumes diag(B)=U (valence electrons), offdiag(B)= V or 2*V.
    Heuristic: if median positive offdiagonal ≈ 2 or 4, divide by 2.
    """
    n = B.shape[0]
    off = B.copy()
    np.fill_diagonal(off, 0.0)
    vals = off[np.triu_indices(n, k=1)]
    pos = vals[vals > 1e-12]
    if pos.size == 0:
        return off
    med = np.median(pos)
    if 1.5 < med < 2.5 or 3.0 < med < 5.0:
        return off / 2.0
    return off


def _load_frame(frame_dir: Path):
    """Return dict with keys: U (N,), V (N,N), AM (N,N) or None, labels (N,)"""
    candidates_bem = (
        list(frame_dir.glob("*BEM*.*")) + list(frame_dir.glob("*bem*.*")) +
        list(frame_dir.glob("*bond*electron*.*")) + list(frame_dir.glob("*bond*electrons*.*")) +
        list(frame_dir.glob("*bond*_*electron*.*")) + list(frame_dir.glob("*bond*.*")) + list(frame_dir.glob("*electron*.*")) +
        list(frame_dir.glob("BEM.*"))
    )
    candidates_am  = list(frame_dir.glob("*AM*.*")) + list(frame_dir.glob("*adj*.*")) + list(frame_dir.glob("AM.*"))
    bem_path = candidates_bem[0] if candidates_bem else None
    am_path  = candidates_am[0]  if candidates_am else None
    if bem_path is None:
        raise FileNotFoundError(f"No BEM-like file found in {frame_dir}")
    B, labels = _load_labelled_matrix(bem_path)
    U = np.diag(B).copy()
    V = _auto_detect_V_from_BEM(B)
    AM = None
    if am_path:
        AM, _ = _load_labelled_matrix(am_path)
    return dict(U=U, V=V, AM=AM, labels=labels)


def _node_labels_from_any(labels, N):
    if labels and len(labels)==N:
        return [str(x) for x in labels]
    return [f"A{i}" for i in range(N)]


def _graph_edges_from_AM_or_V(AM: Optional[np.ndarray], V: np.ndarray, thresh=1e-8) -> List[Tuple[int,int]]:
    N = V.shape[0]
    edges = set()
    if AM is not None:
        if isinstance(AM, tuple):
            AM = AM[0]
        AM = np.asarray(AM, dtype=float)
        if AM.ndim != 2 or AM.shape[0] != AM.shape[1]:
            AM = None
    if AM is not None:
        idx = np.argwhere(AM > thresh)
        for i,j in idx:
            if i<j:
                edges.add((int(i),int(j)))
    else:
        idx = np.argwhere(np.triu(V,1) > thresh)
        for i,j in idx:
            edges.add((int(i), int(j)))
    return sorted(edges)


def _incidence_matrix(N:int, edges: List[Tuple[int,int]]):
    oriented = [(i,j) if i<j else (j,i) for (i,j) in edges]
    E = len(oriented)
    B = np.zeros((N, E), dtype=float)
    for e,(i,j) in enumerate(oriented):
        B[i,e] = -1.0
        B[j,e] = +1.0
    return B, oriented


def _min_norm_flow(B: np.ndarray, dR: np.ndarray, ridge: float=1e-12) -> np.ndarray:
    BBt = B @ B.T
    N = BBt.shape[0]
    A = BBt + ridge*np.eye(N)
    x = np.linalg.solve(A, dR)
    f = B.T @ x
    return f


def _aggregate_window_flows(U_seq, V_seq, AM_seq, ts_idx: int, W: int, v_thresh: float = 1e-8):
    Kp1 = len(U_seq)
    N = U_seq[0].shape[0]
    k0 = max(0, ts_idx - W)
    k1 = min(Kp1-2, ts_idx + W - 1)

    F = {}  # transport (i->j)
    S_to_bond = {}
    S_from_bond = {}

    for k in range(k0, k1+1):
        U0, U1 = U_seq[k], U_seq[k+1]
        V0, V1 = V_seq[k], V_seq[k+1]
        AM = AM_seq[k] if AM_seq[k] is not None else AM_seq[k+1]
        R0 = U0 + V0.sum(axis=1)
        R1 = U1 + V1.sum(axis=1)
        dR = R1 - R0
        dE = 2.0*(V1 - V0)
        edges = _graph_edges_from_AM_or_V(AM, np.maximum(V0,V1), thresh=v_thresh)
        if not edges:
            edges = _graph_edges_from_AM_or_V(None, np.maximum(V0,V1), thresh=v_thresh)
        B, oriented = _incidence_matrix(N, edges)
        f = _min_norm_flow(B, dR)
        for e,(i,j) in enumerate(oriented):
            val = f[e]
            tol_abs=1e-5
            if abs(val) < tol_abs: continue
            if val > 0:
                F[(i,j)] = F.get((i,j), 0.0) + val
            elif val < 0:
                F[(j,i)] = F.get((j,i), 0.0) + (-val)
        for (i,j) in edges:
            dEij = dE[i,j]
            if abs(dEij) < 1e-12:
                continue
            if dEij > 0:
                S_to_bond[(i,(i,j))]  = S_to_bond.get((i,(i,j)),0.0) + 0.5*dEij
                S_to_bond[(j,(i,j))]  = S_to_bond.get((j,(i,j)),0.0) + 0.5*dEij
            else:
                S_from_bond[((i,j),i)] = S_from_bond.get(((i,j),i),0.0) + 0.5*(-dEij)
                S_from_bond[((i,j),j)] = S_from_bond.get(((i,j),j),0.0) + 0.5*(-dEij)
    return F, S_to_bond, S_from_bond


def _find_bond_events_from_V(V_seq: List[np.ndarray], v_min: float, ts_idx: int):
    Kp1 = len(V_seq)
    N = V_seq[0].shape[0]
    events = []
    V0 = V_seq[0]; VK = V_seq[-1]
    changed = np.argwhere(np.triu(np.abs(VK - V0), 1) >= v_min)
    for (i, j) in changed:
        dV = np.diff([V[i, j] for V in V_seq])
        if dV.size == 0:
            continue
        k_star = int(np.argmax(np.abs(dV)))
        k_event = k_star + 1
        delta_total = float(VK[i, j] - V0[i, j])
        kind = "form" if delta_total > 0 else "break" if delta_total < 0 else "unchanged"
        if kind == "unchanged":
            continue
        events.append(dict(i=int(i), j=int(j), kind=kind,
                           k_event=k_event, delta_bo=delta_total,
                           dist_to_ts=int(k_event - ts_idx)))
    return events


def _classify_concerted_sequential(V_seq, ts_idx: int, window: int, v_min: float):
    events = _find_bond_events_from_V(V_seq, v_min=v_min, ts_idx=ts_idx)
    if not events:
        return dict(label="no-rearrangements", synchrony_index=0.0,
                    skew_break_minus_form=None, n_events=0, n_break=0, n_form=0,
                    within_ts_frac=1.0, events=[])
    k_list = [ev["k_event"] for ev in events]
    within = [abs(k - ts_idx) <= window for k in k_list]
    label = "concerted" if all(within) else "sequential"
    span = max(k_list) - min(k_list)
    S = span / (2.0 * max(window, 1))
    krel_break = [ev["k_event"] - ts_idx for ev in events if ev["kind"] == "break"]
    krel_form  = [ev["k_event"] - ts_idx for ev in events if ev["kind"] == "form"]
    skew = None
    if krel_break and krel_form:
        import numpy as _np
        skew = float(_np.median(krel_break) - _np.median(krel_form))
    return dict(label=label, synchrony_index=float(S),
                skew_break_minus_form=skew,
                n_events=len(events), n_break=sum(1 for ev in events if ev["kind"]=="break"),
                n_form=sum(1 for ev in events if ev["kind"]=="form"),
                within_ts_frac=float(sum(within)/len(within)), events=events)


def _write_mechanism_summary(outdir: Path, summary: Dict, ts_idx: int, window: int, v_min: float):
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "mechanism_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label","synchrony_index","skew_break_minus_form",
                    "n_events","n_break","n_form","within_ts_frac","ts_idx","window","event_thresh_bo"])
        w.writerow([
            summary["label"], f"{summary['synchrony_index']:.6f}",
            "" if summary["skew_break_minus_form"] is None else f"{summary['skew_break_minus_form']:.6f}",
            summary["n_events"], summary["n_break"], summary["n_form"],
            f"{summary['within_ts_frac']:.6f}", ts_idx, window, v_min
        ])
    with open(outdir / "bond_events.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i","j","kind","k_event","dist_to_ts","delta_bo"])
        for ev in summary["events"]:
            w.writerow([ev["i"], ev["j"], ev["kind"], ev["k_event"], ev["dist_to_ts"], f"{ev['delta_bo']:.6f}"])

# ------------------- Sankey node/link builders -------------------

def _build_nodes_and_links_transport(args, F, labels, topK: int, min_frac: float):
    tot = sum(F.values()) if F else 0.0
    if tot <= 0:
        return [], []
    items = sorted(F.items(), key=lambda kv: kv[1], reverse=True)
    kept: List[Tuple[Tuple[int,int], float]] = []
    other_in: Dict[int,float] = {}
    other_out: Dict[int,float] = {}
    agg_trans_T = 0.0

    for (i,j), v in items:
        if len(kept) < topK and v >= min_frac*tot:
            kept.append(((i,j), v))
        else:
            if args.transport_aggregate:
                agg_trans_T += v
            else:
                other_out[i] = other_out.get(i, 0.0) + v
                other_in[j]  = other_in.get(j,  0.0) + v

    # Build nodes FIRST
    node_ids: Dict[Tuple[str,object], int] = {}
    nodes: List[Dict] = []

    def add_node(kind: str, idx: object, label: str):
        key = (kind, idx)
        if key not in node_ids:
            node_ids[key] = len(nodes)
            nodes.append(dict(id=len(nodes), key=key, label=label, kind=kind))
        return node_ids[key]

    # Atoms
    for i, lab in enumerate(labels):
        add_node("atom", i, lab)

    # Aggregation endpoints
    if args.transport_aggregate:
        add_node("other", "trans", "Other transport")
    else:
        for i in sorted(other_out):
            add_node("other_out", i, f"Other from {labels[i]}")
        for j in sorted(other_in):
            add_node("other_in", j,  f"Other to {labels[j]}")

    # Now LINKS
    links: List[Dict] = []
    for (i,j), v in kept:
        links.append(dict(source=node_ids[("atom", i)], target=node_ids[("atom", j)],
                          value=float(v), kind="transport", edge=f"{labels[i]}->{labels[j]}"))

    if args.transport_aggregate:
        if agg_trans_T > 0:
            links.append(dict(source=node_ids[("other","trans")], target=node_ids[("other","trans")],
                              value=float(agg_trans_T), kind="transport-agg", edge="Other transport"))
    else:
        for i, v in other_out.items():
            links.append(dict(source=node_ids[("atom", i)], target=node_ids[("other_out", i)],
                              value=float(v), kind="transport-agg-out", edge=f"{labels[i]}→Other"))
        for j, v in other_in.items():
            links.append(dict(source=node_ids[("other_in", j)], target=node_ids[("atom", j)],
                              value=float(v), kind="transport-agg-in", edge=f"Other→{labels[j]}"))

    return nodes, links


def _build_nodes_and_links_two_rail(args, F, S_to_bond, S_from_bond, labels, topK: int, min_frac: float):
    tot = sum(F.values()) + sum(S_to_bond.values()) + sum(S_from_bond.values())
    if tot <= 0:
        return [], []

    all_links = []
    for (i,j), v in F.items():
        all_links.append(("trans", (i,j), v))
    for (i,ij), v in S_to_bond.items():
        all_links.append(("toBond", (i,ij), v))
    for (ij,i), v in S_from_bond.items():
        all_links.append(("fromBond", (ij,i), v))
    all_links.sort(key=lambda x: x[2], reverse=True)

    other_out: Dict[int,float] = {}
    other_in:  Dict[int,float] = {}

    kept = []
    rest = []
    agg_trans = 0.0
    for t, key, v in all_links:
        if len(kept) < topK and v >= min_frac*tot:
            kept.append((t, key, v))
        else:
            if t == "trans":
                if args.transport_aggregate:
                    agg_trans += v
                else:
                    i, j = key
                    other_out[i] = other_out.get(i, 0.0) + v
                    other_in[j]  = other_in.get(j,  0.0) + v
            else:
                rest.append((t, key, v))

    # Build nodes
    node_ids: Dict[Tuple[object,...], int] = {}
    nodes: List[Dict] = []
    '''
    def add_node(key, label, kind):
        if key not in node_ids:
            node_ids[key] = len(nodes)
            nodes.append(dict(id=len(nodes), key=str(key), label=label, kind=kind))
        return node_ids[key]
    '''
    def add_node(key, label, kind):
        if key not in node_ids:
            node_ids[key] = len(nodes)
            # set default x unless overridden by kind
            x = None
            if kind == "atom":
                x = 0.12         # left rail (atoms)
            elif kind == "bond":
                x = 0.82         # right rail (bonds)
            elif kind == "other_trans":
                x = 0.20         # << keep the donut CLOSE to atoms
            elif kind == "other_trans_out":
                x = 0.08
            elif kind == "other_trans_in":
                x = 0.92
            elif kind == "other_toBond":
                x = 0.76
            elif kind == "other_fromBond":
                x = 0.88
            nodes.append(dict(id=len(nodes), key=str(key), label=label, kind=kind, x=x))
        return node_ids[key]

    # Atoms
    for i, lab in enumerate(labels):
        add_node(("atom", i), lab, "atom")

    # Transport aggregation endpoints
    if args.transport_aggregate:
        add_node(("other","trans"), "Other transport", "other_trans")
        nid = add_node("other", "trans", "Other transport")
        nodes[nid]["x"] = 0.20   # << keep it near the atoms
    else:
        for i in sorted(other_out):
            add_node(("other_out", i), f"Other from {labels[i]}", "other_trans_out")
        for j in sorted(other_in):
            add_node(("other_in",  j), f"Other to {labels[j]}",   "other_trans_in")

    # Bond nodes for kept storage links
    kept_bonds = set()
    for t,key,v in kept:
        if t == "toBond":
            _, (a,b) = key
            kept_bonds.add(tuple(sorted((a,b))))
        elif t == "fromBond":
            (a,b), _ = key
            kept_bonds.add(tuple(sorted((a,b))))
    for (a,b) in sorted(kept_bonds):
        add_node(("bond", a, b), f"({labels[a]}–{labels[b]})", "bond")

    # Global storage aggregators
    add_node(("other","toBond"),   "Other -> bond storage",  "other_toBond")
    add_node(("other","fromBond"), "Other <- bond storage",  "other_fromBond")

    # LINKS
    links: List[Dict] = []
    for t,key,v in kept:
        if t == "trans":
            i,j = key
            links.append(dict(source=node_ids[("atom", i)], target=node_ids[("atom", j)],
                              value=float(v), kind="transport", edge=f"{labels[i]}->{labels[j]}"))
        elif t == "toBond":
            i,(a,b) = key; a,b = tuple(sorted((a,b)))
            links.append(dict(source=node_ids[("atom", i)], target=node_ids[("bond", a, b)],
                              value=float(v), kind="storage_in",
                              edge=f"{labels[i]}->({labels[a]}-{labels[b]})"))
        else:  # fromBond
            (a,b), i = key; a,b = tuple(sorted((a,b)))
            links.append(dict(source=node_ids[("bond", a, b)], target=node_ids[("atom", i)],
                              value=float(v), kind="storage_out",
                              edge=f"({labels[a]}-{labels[b]})->{labels[i]}"))

    # Aggregate remainder: TRANSPORT
    if args.transport_aggregate:
        if agg_trans > 0:
            links.append(dict(source=node_ids[("other","trans")], target=node_ids[("other","trans")],
                              value=float(agg_trans), kind="transport-agg", edge="Other transport"))
    else:
        for i, v in other_out.items():
            links.append(dict(source=node_ids[("atom", i)], target=node_ids[("other_out", i)],
                              value=float(v), kind="transport-agg-out", edge=f"{labels[i]}→Other"))
        for j, v in other_in.items():
            links.append(dict(source=node_ids[("other_in", j)], target=node_ids[("atom", j)],
                              value=float(v), kind="transport-agg-in", edge=f"Other→{labels[j]}"))

    # Aggregate remainder: STORAGE
    agg = {"toBond":0.0, "fromBond":0.0}
    for t,key,v in rest:
        agg[t] += v
    if agg["toBond"] > 0:
        links.append(dict(source=node_ids[("other","toBond")], target=node_ids[("other","toBond")],
                          value=float(agg["toBond"]), kind="storage_in-agg", edge="Other -> bond storage"))
    if agg["fromBond"] > 0:
        links.append(dict(source=node_ids[("other","fromBond")], target=node_ids[("other","fromBond")],
                          value=float(agg["fromBond"]), kind="storage_out-agg", edge="Other <- bond storage"))

    return nodes, links


# ------------------ IO for CSV nodes/links ------------------

def _write_csv_nodes_links(outdir: Path, nodes, links, prefix):
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / f"{prefix}_nodes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id","label","kind","key"])
        for n in nodes:
            w.writerow([n["id"], n["label"], n["kind"], n["key"]])
    with open(outdir / f"{prefix}_links.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source","target","value","kind","edge"])
        for l in links:
            w.writerow([l["source"], l["target"], f"{l['value']:.6f}", l["kind"], l["edge"]])

# --------------------------- Plotly ---------------------------

def _try_plotly_sankey(out_html: Path, nodes, links, title,
                       font_size=12, node_thickness=18, node_pad=12,
                       export_svg=False, svg_basename=None,
                       zoom_font_mult=1.8, zoom_node_thickness=None, zoom_node_pad=None):
    try:
        import plotly.graph_objects as go
        from plotly.io import write_image
    except Exception as e:
        try:
            import plotly.graph_objects as go  # type: ignore
        except Exception as e2:
            print(f"[info] plotly not available ({e2}); skipping HTML {title}")
            return False

    source = [l["source"] for l in links]
    target = [l["target"] for l in links]
    value  = [l["value"]  for l in links]
    label  = [n["label"]  for n in nodes]

    # NEW: collect x-coordinates if any
    x = [n.get("x") for n in nodes]
    use_fixed = any(v is not None for v in x)

    node_dict = dict(label=label, pad=node_pad, thickness=node_thickness)
    if use_fixed:
        node_dict["x"] = x

    fig = go.Figure(data=[go.Sankey(
        arrangement=("fixed" if use_fixed else "snap"),   # << honor our x if present
        node=node_dict,
        link=dict(source=source, target=target, value=value)
    )])

    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(label=label, pad=node_pad, thickness=node_thickness),
        link=dict(source=source, target=target, value=value)
    )])
    fig.update_layout(title=title, font_size=font_size)#, width=1400, height=1200)

    # HTML with camera defaulting to SVG
    fig.write_html(str(out_html), include_plotlyjs="cdn",
                   config={"toImageButtonOptions": {"format": "svg"}})

    # Always-on "SVG (zoomed)" modebar button
    try:
        html = out_html.read_text()
        base = (svg_basename or out_html.with_suffix('').name)
        button_js = f'''
<script>
(function(){{
  var gd = document.querySelector('div.js-plotly-plot');
  if(!gd || !window.Plotly) return;
  function addButton(){{
    var mb = gd._fullLayout && gd._fullLayout._modeBar;
    if(!mb || !mb.addButton) return;
    mb.addButton({{name:'svg-zoom', title:'Download SVG (zoomed labels)',
      icon: Plotly.Icons.camera,
      click: function(){{
        var traceIndex = 0;
        var origFont = (gd.layout && gd.layout.font && gd.layout.font.size) || {font_size};
        var origThk  = (gd.data[traceIndex] && gd.data[traceIndex].node && gd.data[traceIndex].node.thickness) || {node_thickness};
        var origPad  = (gd.data[traceIndex] && gd.data[traceIndex].node && gd.data[traceIndex].node.pad) || {node_pad};
        var zFont = Math.round(origFont * {zoom_font_mult});
        var zThk  = Math.round(origThk * 1.2);
        var zPad  = origPad;
        Promise.all([
          Plotly.relayout(gd, {{'font.size': zFont}}),
          Plotly.restyle(gd, {{'node.thickness':[zThk], 'node.pad':[zPad]}}, [traceIndex])
        ]).then(function(){{
          return Plotly.downloadImage(gd, {{format:'svg', filename:'{base}_zoom'}});
        }}).then(function(){{
          return Promise.all([
            Plotly.relayout(gd, {{'font.size': origFont}}),
            Plotly.restyle(gd, {{'node.thickness':[origThk], 'node.pad':[origPad]}}, [traceIndex])
          ]);
        }});
      }}
    }});
  }}
  if(document.readyState==='loading'){{
    document.addEventListener('DOMContentLoaded', function(){{ setTimeout(addButton, 100); }});
  }} else {{ setTimeout(addButton, 100); }}
}})();
</script>
'''
        out_html.write_text(html + "\n" + button_js)
    except Exception as _e:
        print("[warn] Could not inject persistent SVG (zoomed) button:", _e)

    # Tight export copy (no title, small margins, tall canvas)
    fig_tight = go.Figure(fig)
    fig_tight.update_layout(margin=dict(l=4, r=4, t=4, b=4), title=None)
                            #width=1400, height=1200

    if export_svg:
        base = svg_basename or out_html.with_suffix('').name
        try:
            write_image(fig_tight, str(out_html.with_name(f"{base}.svg")), format="svg")
            z_font  = int(round(font_size * zoom_font_mult))
            z_thick = zoom_node_thickness if zoom_node_thickness is not None else int(round(node_thickness*1.2))
            z_pad   = zoom_node_pad if zoom_node_pad is not None else node_pad
            fig_zoom = go.Figure(fig_tight)
            fig_zoom.update_layout(font_size=z_font)
            try:
                fig_zoom.data[0]['node']['thickness'] = z_thick
                fig_zoom.data[0]['node']['pad'] = z_pad
            except Exception:
                pass
            write_image(fig_zoom, str(out_html.with_name(f"{base}_zoom.svg")), format="svg")
            print(f"[ok] SVGs: {out_html.with_name(f'{base}.svg')} and {out_html.with_name(f'{base}_zoom.svg')}")
        except Exception as e:
            print(f"[warn] SVG export failed (install conda-forge 'python-kaleido'): {e}")
            # Fallback: inject a zoomed-SVG button if not already present (best-effort)
            try:
                html2 = out_html.read_text()
                if 'name:\'svg-zoom\'' not in html2:
                    out_html.write_text(html2 + "\n" + button_js)
                    print("[info] Added a 'Download SVG (zoomed)' button to the HTML as a fallback.")
            except Exception as e2:
                print(f"[warn] Could not inject zoom button: {e2}")

    return True


# --------------------------- Frame loading ---------------------------

def load_irc_frames(root: Path):
    root = Path(root / "finished_irc")
    ts_path = root / "ts_frame.txt"
    if not ts_path.exists():
        raise FileNotFoundError(f"Missing ts_frame.txt at {ts_path}")
    ts_idx = _read_text_number(ts_path)

    subdirs = [p for p in root.iterdir() if p.is_dir()]
    if subdirs:
        def keydir(p):
            m = re.search(r"(\d+)", p.name)
            return (int(m.group(1)) if m else 10**9, p.name)
        subdirs = sorted(subdirs, key=keydir)
        U_seq, V_seq, AM_seq = [], [], []
        labels = None
        for d in subdirs:
            data = _load_frame(d)
            if labels is None:
                labels = _node_labels_from_any(data.get("labels"), len(data["U"]))
            U_seq.append(data["U"]); V_seq.append(data["V"]); AM_seq.append(data["AM"])
        return U_seq, V_seq, AM_seq, labels, ts_idx

    files = list(root.glob("*"))
    buckets = {}
    for p in files:
        if p.name == "ts_frame.txt":
            continue
        m = re.search(r"(\d+)", p.stem)
        if not m:
            continue
        k = int(m.group(1))
        buckets.setdefault(k, []).append(p)
    if not buckets:
        raise FileNotFoundError(f"No per-frame files or subfolders found in {root}")

    U_seq, V_seq, AM_seq = [], [], []
    labels = None
    for k in sorted(buckets.keys()):
        U = V = AM = None
        bems = [p for p in buckets[k] if re.search(r'bem|bond.*electron', p.name, re.I)]
        bem_path = bems[0] if bems else None
        if bem_path is None:
            squares = []
            for p in buckets[k]:
                try:
                    M, lab = _load_labelled_matrix(p)
                    if M.shape[0] == M.shape[1]:
                        squares.append((p, M, lab))
                except Exception:
                    pass
            if squares:
                bem_path, B, _labels = squares[0]
            else:
                raise FileNotFoundError(f"No BEM-like file for frame {k}")
        else:
            B, _labels = _load_labelled_matrix(bem_path)
        U = np.diag(B).copy()
        V = _auto_detect_V_from_BEM(B)
        ams = [p for p in buckets[k] if re.search(r'adj|am', p.name, re.I)]
        AM = None
        if ams:
            AM = _load_labelled_matrix(ams[0])
        if labels is None:
            labels = _node_labels_from_any(_labels, len(U))
        U_seq.append(U); V_seq.append(V); AM_seq.append(AM)
    return U_seq, V_seq, AM_seq, labels, ts_idx


# ------------------------------ CLI ------------------------------

def main():
    ap = argparse.ArgumentParser()
    # SVG/zoom
    ap.add_argument('--export-svg', action='store_true', help='Export SVG next to HTML (requires python-kaleido).')
    ap.add_argument('--font-size', type=int, default=12, help='Base font size for labels (default 12).')
    ap.add_argument('--node-thickness', type=int, default=18, help='Sankey node thickness (default 18).')
    ap.add_argument('--node-pad', type=int, default=12, help='Sankey node pad (default 12).')
    ap.add_argument('--zoom-font-mult', type=float, default=1.8, help='Multiplier for zoomed-in SVG labels (default 1.8).')
    ap.add_argument('--zoom-node-thickness', type=int, default=0, help='Optional thickness for zoomed-in SVG (0 => auto 1.2x).')
    ap.add_argument('--zoom-node-pad', type=int, default=0, help='Optional pad for zoomed-in SVG (0 => keep).')

    # Core
    ap.add_argument('--root', required=True, help='Path to IRC folder containing per-frame data and ts_frame.txt')
    ap.add_argument('--window', type=int, default=20, help='Half-window in frames around TS (default: 20)')
    ap.add_argument('--topK', type=int, default=30, help='Max number of explicit links to draw (rest aggregated)')
    ap.add_argument('--min-frac', type=float, default=0.03, help='Drop links < min_frac of total (unless within topK)')
    ap.add_argument('--out', default='out_sankey', help='Output directory')

    # Mechanism timing
    ap.add_argument('--event-thresh', type=float, default=0.3, help='Minimum |ΔBO| to flag rearrangement (default 0.3 BO).')
    ap.add_argument('--hist-bins', default='-50:50:1', help="Bins for TS histograms as 'start:end:step' in frames")
    ap.add_argument('--cum-N', default='5,10,20', help='Comma-separated absolute frame windows for cumulative fractions')
    ap.add_argument('--write-event-times', action='store_true', help='Write raw per-event times CSV')

    # Toggle: donut vs. detailed per-node Others for TRANSPORT remainder
    ap.add_argument('--transport-aggregate', action='store_true',
                    help='Show transport remainder as a single "Other transport" donut (default: detailed per-node)')

    args = ap.parse_args()

    root = Path(args.root)
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    print(f"[info] Loading IRC frames from: {root}")
    U_seq, V_seq, AM_seq, labels, ts_idx = load_irc_frames(root)
    print(f"[info] Frames: {len(U_seq)}, Atoms: {len(labels)}, TS frame index: {ts_idx}")
    # Mechanism summary
    mech = _classify_concerted_sequential(V_seq, ts_idx=ts_idx, window=args.window, v_min=args.event_thresh)
    def _parse_bins(spec: str):
        import numpy as _np
        start, end, step = [int(x) for x in spec.split(":")]
        edges = _np.arange(start, end + step, step, dtype=int)
        centers = (edges[:-1] + edges[1:]) / 2.0
        return edges, centers
    def _compute_ts_histograms(events, edges):
        import numpy as _np
        dists = _np.array([e["dist_to_ts"] for e in events], dtype=int) if events else _np.empty((0,), dtype=int)
        any_counts, _ = _np.histogram(dists, bins=edges)
        bd = _np.array([e["dist_to_ts"] for e in events if e["kind"] == "break"], dtype=int)
        fd = _np.array([e["dist_to_ts"] for e in events if e["kind"] == "form"], dtype=int)
        br_counts, _ = _np.histogram(bd, bins=edges)
        fo_counts, _ = _np.histogram(fd, bins=edges)
        return {"break": br_counts, "form": fo_counts, "any": any_counts}
    def _compute_cumulative(events, Ns):
        import numpy as _np
        if not events:
            return [[int(N), 0.0, 0.0, 0.0] for N in Ns]
        d_all = _np.array([e["dist_to_ts"] for e in events], dtype=int)
        d_br  = _np.array([e["dist_to_ts"] for e in events if e["kind"] == "break"], dtype=int)
        d_fo  = _np.array([e["dist_to_ts"] for e in events if e["kind"] == "form"], dtype=int)
        rows = []
        for N in Ns:
            N = abs(int(N))
            any_frac = (abs(d_all) <= N).sum() / max(len(d_all), 1)
            br_frac  = (abs(d_br)  <= N).sum() / max(len(d_br),  1) if len(d_br) else 0.0
            fo_frac  = (abs(d_fo)  <= N).sum() / max(len(d_fo),  1) if len(d_fo) else 0.0
            rows.append([int(N), float(br_frac), float(fo_frac), float(any_frac)])
        return rows

    # Write mechanism summary files
    def _write_ts_hist_outputs(outdir, centers, histo, cum_rows):
        outdir.mkdir(parents=True, exist_ok=True)
        with open(outdir / "ts_histograms.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["bin_center", "count_break", "count_form", "count_any"])
            for i in range(len(centers)):
                w.writerow([f"{centers[i]:.1f}", int(histo["break"][i]), int(histo["form"][i]), int(histo["any"][i])])
        with open(outdir / "ts_cumulative.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["N_frames", "frac_break", "frac_form", "frac_any"])
            for N, fb, ff, fa in cum_rows:
                w.writerow([N, f"{fb:.6f}", f"{ff:.6f}", f"{fa:.6f}"])

    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "mechanism_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label","synchrony_index","skew_break_minus_form","n_events","n_break","n_form","within_ts_frac","ts_idx","window","event_thresh_bo"])
        w.writerow([
            mech["label"], f"{mech['synchrony_index']:.6f}",
            "" if mech.get("skew_break_minus_form") is None else f"{mech['skew_break_minus_form']:.6f}",
            mech["n_events"], mech["n_break"], mech["n_form"], f"{mech['within_ts_frac']:.6f}",
            ts_idx, args.window, args.event_thresh
        ])

    # Histograms/cumulative
    try:
        edges, centers = _parse_bins(args.hist_bins)
        histo = _compute_ts_histograms(mech["events"], edges)
        Ns = [int(x) for x in args.cum_N.split(',') if x.strip()]
        cum_rows = _compute_cumulative(mech["events"], Ns)
        _write_ts_hist_outputs(outdir, centers, histo, cum_rows)
        print(f"[mech] wrote ts_histograms.csv and ts_cumulative.csv  (bins={args.hist_bins}, cum={args.cum_N})")
    except Exception as e:
        print(f"[warn] histogram/cumulative generation failed: {e}")

    if args.write_event_times:
        with open(outdir / "event_times.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["i","j","kind","dist_to_ts","k_event","delta_bo"])
            for ev in mech["events"]:
                w.writerow([ev["i"], ev["j"], ev["kind"], ev["dist_to_ts"], ev["k_event"], f"{ev['delta_bo']:.6f}"])
        print("[mech] wrote event_times.csv")

    # Aggregate flows over TS window
    F, S_to_bond, S_from_bond = _aggregate_window_flows(U_seq, V_seq, AM_seq, ts_idx, args.window)

    # TRANSPORT Sankey
    nodes_T, links_T = _build_nodes_and_links_transport(args, F, labels, args.topK, args.min_frac)
    _write_csv_nodes_links(outdir, nodes_T, links_T, prefix="transport")
    html_T = outdir / "sankey_transport.html"
    ok_T = _try_plotly_sankey(
        html_T, nodes_T, links_T, title="",
        font_size=args.font_size, node_thickness=args.node_thickness, node_pad=args.node_pad,
        export_svg=args.export_svg, svg_basename='sankey_transport',
        zoom_font_mult=args.zoom_font_mult,
        zoom_node_thickness=(None if args.zoom_node_thickness==0 else args.zoom_node_thickness),
        zoom_node_pad=(None if args.zoom_node_pad==0 else args.zoom_node_pad)
    )
    if ok_T:
        print(f"[ok] Wrote {html_T} and CSVs under {outdir}")

    # TWO-RAIL Sankey
    nodes_B, links_B = _build_nodes_and_links_two_rail(args, F, S_to_bond, S_from_bond, labels, args.topK, args.min_frac)
    _write_csv_nodes_links(outdir, nodes_B, links_B, prefix="two_rail")
    html_B = outdir / "sankey_two_rail.html"
    ok_B = _try_plotly_sankey(
        html_B, nodes_B, links_B, title="",
        font_size=args.font_size, node_thickness=args.node_thickness, node_pad=args.node_pad,
        export_svg=args.export_svg, svg_basename='sankey_two_rail',
        zoom_font_mult=args.zoom_font_mult,
        zoom_node_thickness=(None if args.zoom_node_thickness==0 else args.zoom_node_thickness),
        zoom_node_pad=(None if args.zoom_node_pad==0 else args.zoom_node_pad)
    )
    if ok_B:
        print(f"[ok] Wrote {html_B} and CSVs under {outdir}")

    print("[done]")

if __name__ == "__main__":
    main()

