from __future__ import annotations
import csv, re
from typing import List, Dict

def _parse_matrix_lines(lines: list[str]) -> tuple[list[str], list[list[int]]]:
    lines = [ln for ln in lines if ln.strip()]
    if not lines: return [], []
    header_tokens = next(csv.reader([lines[0]]))
    if header_tokens and header_tokens[0] == "":
        labels = [t.strip() for t in header_tokens[1:] if t.strip() != ""]
    else:
        labels = [t.strip() for t in header_tokens if t.strip() != ""]
    n = len(labels)
    if n == 0: return [], []
    matrix: list[list[int]] = []
    rdr = csv.reader(lines[1:])
    for row in rdr:
        if not row: continue
        offset = len(row) - n
        if offset >= 1: tokens = row[offset:offset+n]
        elif offset == 0: tokens = row[:n]
        else: tokens = (row + ["0"]*(n-len(row)))[:n]
        vals: list[int] = []
        for tok in tokens:
            t = tok.strip()
            if t == "": vals.append(0); continue
            try: vals.append(int(round(float(t))))
            except ValueError: vals.append(0)
        matrix.append(vals)
    return labels, matrix

def load_labelled_sections(csv_path) -> dict:
    from pathlib import Path
    csv_path = Path(csv_path)
    if not csv_path.exists(): return {}
    text = csv_path.read_text()
    sections = {}; current = None; buf: list[str] = []
    def _flush():
        nonlocal sections, current, buf
        if current and buf:
            labels, mat = _parse_matrix_lines(buf)
            sections[current] = {"labels": labels, "matrix": mat}
        buf = []
    for line in text.splitlines():
        tag = line.strip().lower()
        if tag in ("reactant","product"):
            _flush(); current = tag; continue
        buf.append(line)
    _flush()
    return sections

_ELEMENT_DIGITS = re.compile(r"^[A-Za-z]+[0-9]+$")

def indexed_labels(labels: list[str]) -> list[str]:
    if labels and all(_ELEMENT_DIGITS.match(lab) for lab in labels):
        return labels
    return [f"{sym}{i}" for i, sym in enumerate(labels)]

def coerce_int(x, default=0) -> int:
    try: return int(round(float(x)))
    except Exception: return default

def graph_from_section(section_adj: dict | None, section_be: dict | None) -> dict:
    if not section_adj or not section_be: return {}
    labels = section_adj.get("labels") or section_be.get("labels") or []
    adj    = section_adj.get("matrix", [])
    be     = section_be.get("matrix",  [])
    if not labels or not adj or not be: return {}
    n = min(len(labels), len(adj), len(be))
    idx_labels = indexed_labels(labels[:n])
    electrons = [coerce_int(be[i][i], 0) if i < len(be[i]) else 0 for i in range(n)]
    bond_list = []
    for i in range(n):
        for j in range(i+1, n):
            aij = coerce_int(adj[i][j], 0)
            if aij != 0:
                bo = coerce_int(be[i][j], 0) if j < len(be[i]) else 0
                btype = "dative" if bo == 0 else "sigma"
                bond_list.append([idx_labels[i], idx_labels[j], bo, btype])
    return {"elements": idx_labels, "electrons": electrons, "bond_list": bond_list}

def swap_reactant_product(sections: dict) -> dict:
    out = {}
    if "product" in sections: out["reactant"] = sections["product"]
    if "reactant" in sections: out["product"]  = sections["reactant"]
    return out

