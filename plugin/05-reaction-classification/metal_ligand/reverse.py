from typing import Optional

def invert_change_string(s: Optional[str]) -> Optional[str]:
    if not s: return s
    txt = str(s).strip()
    if "->" not in txt: return txt
    left, right = [t.strip() for t in txt.split("->", 1)]
    return f"{right} -> {left}"

def reverse_record(rec: dict) -> dict:
    out = {k: v for k, v in rec.items()
           if k not in ("bond_changes", "reactive_atom_saturation", "metal_electron_change")}
    bc = rec.get("bond_changes", {})
    def _copy_and_invert(entries: list[list[str]]) -> list[list[str]]:
        flipped: list[list[str]] = []
        for e in entries or []:
            if not e: continue
            a = e[0]; b = e[1] if len(e) > 1 else None
            if b is None: continue
            if len(e) >= 3 and e[2] is not None:
                flipped.append([a, b, invert_change_string(e[2])])
            else:
                flipped.append([a, b])
        return flipped
    out["bond_changes"] = {
        "formed": _copy_and_invert(bc.get("broken", [])),
        "broken": _copy_and_invert(bc.get("formed", [])),
    }
    ras = {}
    for atom, info in (rec.get("reactive_atom_saturation", {}) or {}).items():
        info = dict(info or {})
        if "before" in info and "after" in info:
            info["before"], info["after"] = info["after"], info["before"]
        if "max_bond_order_before" in info and "max_bond_order_after" in info:
            info["max_bond_order_before"], info["max_bond_order_after"] = (
                info["max_bond_order_after"],
                info["max_bond_order_before"]
            )
        ras[atom] = info
    out["reactive_atom_saturation"] = ras
    mec = {}
    for m, info in (rec.get("metal_electron_change", {}) or {}).items():
        info = dict(info or {})
        if "before" in info and "after" in info:
            info["before"], info["after"] = info["after"], info["before"]
        mec[m] = info
    out["metal_electron_change"] = mec
    return out

