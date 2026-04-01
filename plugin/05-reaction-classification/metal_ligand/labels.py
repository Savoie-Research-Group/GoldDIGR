import re
from typing import List, Tuple, Dict, Set
from .constants import METALS

def elem_of(label: str) -> str:
    m = re.match(r"([A-Za-z]+)", str(label))
    if not m: return str(label)
    raw = m.group(1)
    return raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()

def is_metal(label: str) -> bool:
    return elem_of(label) in METALS

def normalize_pair(a: str, b: str) -> Tuple[str,str]:
    return (a,b) if a <= b else (b,a)

def to_pairs(items) -> List[Tuple[str,str]]:
    out = []
    for rec in items or []:
        if not rec: continue
        a = rec[0]; b = rec[1] if len(rec) > 1 else None
        if b is None: continue
        out.append((a,b))
    return out

def build_index(pairs: List[Tuple[str,str]]) -> Dict[str, Set[str]]:
    idx: Dict[str, Set[str]] = {}
    for a,b in pairs:
        idx.setdefault(a,set()).add(b)
        idx.setdefault(b,set()).add(a)
    return idx

