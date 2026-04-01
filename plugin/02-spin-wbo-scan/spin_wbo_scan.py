#!/usr/bin/env python3
import argparse, os, re, json, subprocess, sys
from pathlib import Path
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

XYZ_HEADER_INT = re.compile(r"^\s*(\d+)\s*$")

# ---- BEM CSV helpers ----

# Nominal valence electrons for common elements (extend as needed)
# Nominal valence electrons (H–Rn), keys capitalized to match XYZ symbols
NOMINAL_VALENCE = {
    "H": 1, "He": 2,
    "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
    "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Ga": 3, "Ge": 4, "As": 5, "Se": 6, "Br": 7, "Kr": 8,
    "Rb": 1, "Sr": 2, "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "In": 3, "Sn": 4, "Sb": 5, "Te": 6, "I": 7, "Xe": 8,
    "Cs": 1, "Ba": 2, "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
    "Tl": 3, "Pb": 4, "Bi": 5, "Po": 6, "At": 7, "Rn": 8, 

    # La group (Lanthanides)
    "Ce": 3,"Pr": 3,"Nd": 3,"Pm": 3,"Sm": 3,"Eu": 3,"Gd": 3,
    "Tb": 3,"Dy": 3,"Ho": 3,"Er": 3,"Tm": 3,"Yb": 3,"Lu": 3,

    # Ac group (Actinides)
    "Ac": 3,"Th": 3,"Pa": 3,"U": 3, "Np": 3,"Pu": 3,"Am": 3,"Cm": 3,
    "Bk": 3,"Cf": 3,"Es": 3,"Fm": 3,"Md": 3,"No": 3,"Lr": 3
}

import csv
from glob import glob

def _load_bem_e_csv(path):
    """
    Read frame_#####_bond_electrons.csv → (labels, matrix)
    labels: list like ['C0','H1',...]
    matrix: n×n float matrix (off-diags 0.0 if missing; diagonal may be blank→0.0)
    """
    with open(path, "r", newline="") as f:
        rows = list(csv.reader(f))
    if not rows or not rows[0]:
        raise ValueError(f"Empty/invalid BEM CSV: {path}")
    labels = rows[0][1:]
    n = len(labels)
    mat = [[0.0]*n for _ in range(n)]
    for i in range(n):
        # row i+1: [row_label, v1, v2, ...]
        row = rows[i+1]
        for j in range(n):
            tok = row[j+1] if j+1 < len(row) else ""
            try:
                mat[i][j] = float(tok)
            except Exception:
                mat[i][j] = 0.0
    return labels, mat

def _series_from_frames(workdir):
    """
    Load all frame_#####_bond_electrons.csv in workdir, sorted by frame index.
    Returns: labels, [mat0, mat1, ...], [frame_indices]
    """
    files = sorted(glob(str(Path(workdir)/"frame_*_bond_electrons.csv")))
    if not files:
        return [], [], []
    # sort by numeric frame
    def _frame_idx(p):
        name = Path(p).name
        # frame_00088_bond_electrons.csv → 00088
        return int(name.split("_")[1])
    files.sort(key=_frame_idx)
    labels0, mat0 = _load_bem_e_csv(files[0])
    mats = [mat0]
    frames = [_frame_idx(files[0])]
    for fp in files[1:]:
        lab, m = _load_bem_e_csv(fp)
        # ensure label consistency; if not, intersect in order
        if lab != labels0:
            # conservative: map by name
            idx_map = {lab: i for i, lab in enumerate(lab)}
            reorder = [idx_map.get(L, None) for L in labels0]
            mm = [[0.0]*len(labels0) for _ in range(len(labels0))]
            for i, src_i in enumerate(reorder):
                if src_i is None: continue
                for j, src_j in enumerate(reorder):
                    if src_j is None: continue
                    mm[i][j] = m[src_i][src_j]
            m = mm
        mats.append(m)
        frames.append(_frame_idx(fp))
    return labels0, mats, frames

def summarize_reactive_elements(workdir, thresh=0.5, out_csv="reactive_summary.csv", include_meta=False):
    """
    Compare first and last electron-BEM; pick entries with |Δ| >= thresh.
    For each flagged entry (diag or off-diag), write a row with its full time series.
    Row label convention:
      - diagonal: 'Rh1'
      - off-diagonal: 'C0-Rh1' with alphabetical order as produced in CSV (i<j)
    Columns: key, type, first, last, delta, frame_00000, frame_00001, ...
    """
    labels, mats, frames = _series_from_frames(workdir)
    if not mats:
        print("[reactive] No bond_electrons CSVs found; skipping.")
        return None
    n = len(labels)
    first = mats[0]
    last  = mats[-1]

    # numeric frame column names ("0","1",...)
    frame_cols = [str(k) for k in frames]

    # Build list of candidate keys
    candidates = []
    # diagonals → valence
    for i in range(n):
        dv = last[i][i] - first[i][i]
        if abs(dv) >= thresh:
            candidates.append(("valence", (i, i), labels[i]))
    # off-diagonals → bond (i<j)
    for i in range(n):
        for j in range(i+1, n):
            dv = last[i][j] - first[i][j]
            if abs(dv) >= thresh:
                candidates.append(("bond", (i, j), f"{labels[i]}-{labels[j]}"))

    # --- TRANSPOSED OUTPUT ---
    # Build a dict of {key -> full time series}
    key_to_series = {}
    for typ, (i, j), key in candidates:
        series = [m[i][j] for m in mats]
        key_to_series[key] = series

    # Sort keys (columns) for stable output
    keys_sorted = sorted(key_to_series.keys())

    # Header: frame index numbers only ("0","1","2",...)
    # Rows: one per frame; columns are the keys
    outp = Path(workdir) / out_csv
    with open(outp, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frame"] + keys_sorted)
        for t, fr in enumerate(frames):
            row = [fr]
            for k in keys_sorted:
                row.append(key_to_series[k][t])
            writer.writerow(row)

    print(f"[reactive] Wrote {outp}")
    return str(outp)

def write_bem_csvs(symbols, charges, wbo_pairs, out_dir, frame_idx):
    """
    Write two CSVs per frame:
      frame_#####_bond_electrons.csv  (off-diag = 2*WBO; diag ~ valence - Mulliken charge if available)
      frame_#####_bond_orders.csv     (off-diag = WBO;  diag blank)
    Header row/col are atom labels like C0, H1, ...
    """
    from csv import writer

    n = len(symbols)
    # Sum of WBOs incident to each atom (for conservation-consistent diagonal)
    sum_wbo = [0.0] * n
    for (a, b), w in wbo_pairs.items():
        w = float(w)
        sum_wbo[a] += w
        sum_wbo[b] += w

    labels = [f"{symbols[i]}{i}" for i in range(n)]
    # Build dense matrices
    M_e = [[""] * (n + 1) for _ in range(n + 1)]
    M_o = [[""] * (n + 1) for _ in range(n + 1)]

    # headers
    M_e[0][0] = ""
    M_o[0][0] = ""
    for j in range(n):
        M_e[0][j+1] = labels[j]
        M_o[0][j+1] = labels[j]

    # Initialize all off-diagonals to 0 (so missing WBO → solid 0)
    for i in range(n):
        for j in range(n):
            if i != j:
                M_e[i+1][j+1] = 0.0
                # If you ALSO want zeros in bond_orders.csv off-diagonals, keep the next line:
                M_o[i+1][j+1] = 0.0

    # fill rows
    for i in range(n):
        M_e[i+1][0] = labels[i]
        M_o[i+1][0] = labels[i]
        # diagonals (valence-only): d_i = V0 - sum_j WBO_ij
        sym = symbols[i]
        V0  = NOMINAL_VALENCE.get(sym)
        if V0 is not None:
            di = float(V0) - float(sum_wbo[i])
            # clean tiny negatives due to numerical noise
            if abs(di) < 1e-8:
                di = 0.0
            elif di < 0.0 and di > -1e-3:
                di = 0.0
            M_e[i+1][i+1] = di
        else:
            M_e[i+1][i+1] = ""
        # order-BEM keeps blank diagonal
        M_o[i+1][i+1] = ""
    # off-diagonals from WBO
    for (a,b), w in wbo_pairs.items():
        # electrons matrix gets 2*WBO, order matrix gets WBO
        val_e = 2.0 * float(w)
        val_o = float(w)
        ia, ib = a+1, b+1
        M_e[ia][ib] = val_e
        M_e[ib][ia] = val_e
        M_o[ia][ib] = val_o
        M_o[ib][ia] = val_o

    # write files (use 5-digit frame index to match your YARP naming)
    p_e = Path(out_dir) / f"frame_{frame_idx:05d}_bond_electrons.csv"
    p_o = Path(out_dir) / f"frame_{frame_idx:05d}_bond_orders.csv"
    with open(p_e, "w", newline="") as f:
        w = writer(f)
        for row in M_e: w.writerow(row)
    with open(p_o, "w", newline="") as f:
        w = writer(f)
        for row in M_o: w.writerow(row)

# Element → Atomic Number (H–Rn), symbols capitalized
EL_TO_Z = {
    "H": 1,  "He": 2,
    "Li": 3, "Be": 4,  "B": 5,  "C": 6,  "N": 7,  "O": 8,  "F": 9,  "Ne": 10,
    "Na": 11,"Mg": 12, "Al": 13,"Si": 14,"P": 15, "S": 16,"Cl": 17,"Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21,"Ti": 22,"V": 23, "Cr": 24,"Mn": 25,"Fe": 26,"Co": 27,"Ni": 28,
    "Cu": 29,"Zn": 30,"Ga": 31,"Ge": 32,"As": 33,"Se": 34,"Br": 35,"Kr": 36,
    "Rb": 37,"Sr": 38,"Y": 39, "Zr": 40,"Nb": 41,"Mo": 42,"Tc": 43,"Ru": 44,"Rh": 45,"Pd": 46,
    "Ag": 47,"Cd": 48,"In": 49,"Sn": 50,"Sb": 51,"Te": 52,"I": 53, "Xe": 54,
    "Cs": 55,"Ba": 56,          "Hf": 72,"Ta": 73,"W": 74, "Re": 75,"Os": 76,"Ir": 77,"Pt": 78,
    "Au": 79,"Hg": 80,"Tl": 81,"Pb": 82,"Bi": 83,"Po": 84,"At": 85,"Rn": 86, 
    # La group (Lanthanides)
    "La": 57,"Ce": 58,"Pr": 59,"Nd": 60,"Pm": 61,"Sm": 62,"Eu": 63,"Gd": 64,
    "Tb": 65,"Dy": 66,"Ho": 67,"Er": 68,"Tm": 69,"Yb": 70,"Lu": 71,

    # Ac group (Actinides)
    "Ac": 89,"Th": 90,"Pa": 91,"U": 92, "Np": 93,"Pu": 94,"Am": 95,"Cm": 96,
    "Bk": 97,"Cf": 98,"Es": 99,"Fm": 100,"Md": 101,"No": 102,"Lr": 103
}

import itertools

WBO_PAIRLINE = re.compile(r"^\s*(\d+)\s+(\d+)\s+([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s*$")

def parse_wbo_file(wbo_path):
    """
    Parse xTB .wbo file.
    Supports:
      (A) Pair list format: 'i  j  value' (1-based indices)
      (B) Lower/upper triangular matrix (n lines of floats); we detect by row lengths.
    Returns: dict {(i,j): wbo} with 0-based indices and i<j.
    """
    wbo_pairs = {}
    if not Path(wbo_path).exists():
        return wbo_pairs

    with open(wbo_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.strip().startswith("#")]

    # Try (A) pair list
    hits = 0
    for ln in lines:
        m = WBO_PAIRLINE.match(ln)
        if m:
            i = int(m.group(1)) - 1
            j = int(m.group(2)) - 1
            v = float(m.group(3))
            if i != j:
                a, b = (i, j) if i < j else (j, i)
                prev = wbo_pairs.get((a, b))
                wbo_pairs[(a, b)] = v if (prev is None or abs(v) > abs(prev)) else prev
            hits += 1
    if hits > 0:
        return wbo_pairs

    # Try (B) triangular / square matrix
    # Detect row token counts
    token_rows = [ln.split() for ln in lines]
    lens = [len(r) for r in token_rows]
    n = max(lens) if lens else 0
    if n > 1 and all(all(re.match(r"[+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?$", t) for t in row) for row in token_rows):
        # Accept either n rows or >= n rows; use first n rows of width <= n
        mat = [[0.0]*n for _ in range(n)]
        for i, row in enumerate(token_rows[:n]):
            for j, tok in enumerate(row[:n]):
                try:
                    mat[i][j] = float(tok)
                except:
                    mat[i][j] = 0.0
        # Symmetrize and extract upper triangle
        for i in range(n):
            for j in range(i+1, n):
                v = 0.5*(mat[i][j] + mat[j][i])
                if abs(v) > 0.0:
                    wbo_pairs[(i, j)] = v
        return wbo_pairs

    # Fallback: nothing recognized
    return wbo_pairs

def read_multi_xyz_like_trj(path):
    recs = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        m = XYZ_HEADER_INT.match(lines[i])
        if not m: i += 1; continue
        nat = int(m.group(1))
        if i+1+nat >= n: break
        comment = lines[i+1].rstrip("\n")
        block = lines[i+2:i+2+nat]
        symbols, coords, ok = [], [], True
        for ln in block:
            ps = ln.strip().split()
            if len(ps) < 4: ok=False; break
            sym = ps[0]
            try: x,y,z = map(float, ps[1:4])
            except: ok=False; break
            symbols.append(sym); coords.append((x,y,z))
        if ok:
            recs.append((symbols, coords, comment))
            i = i + 2 + nat
        else:
            i += 1
    return recs

def write_xyz(symbols, coords, outpath, comment=""):
    with open(outpath, "w") as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        for s,(x,y,z) in zip(symbols, coords):
            f.write(f"{s:2s} {x: .8f} {y: .8f} {z: .8f}\n")

def _norm_sym(sym: str) -> str:
    # Handle 'fe'/'FE' → 'Fe', etc.
    return sym[:1].upper() + sym[1:].lower() if sym else sym

def total_electrons(symbols, charge):
    """Sum atomic numbers from EL_TO_Z and subtract net charge."""
    zsum = 0
    for s in symbols:
        zs = EL_TO_Z.get(_norm_sym(s))
        if zs is None:
            raise ValueError(f"Unknown element symbol '{s}'. Please extend EL_TO_Z.")
        zsum += zs
    return zsum - charge

def base_uhf_from_parity(symbols, charge):
    return 0 if (total_electrons(symbols, charge) % 2 == 0) else 1

# ---------------- Pass 1: spin-polarized energy (tblite) ----------------
def run_xtb_energy_tblite(xyz_path, charge, uhf, workdir, acc=0.2, maxiter=500):
    ns = f"{Path(xyz_path).stem}_UHF{uhf}_SPIN"
    cmd = [
        "xtb", xyz_path, "--scc",
        "--spinpol", "--tblite",
        "--chrg", str(charge), "--uhf", str(uhf),
        "--acc", str(acc), "--iterations", str(maxiter),
        "--namespace", ns
    ]
    out_file = Path(workdir) / f"{ns}.out"
    with open(out_file, "w") as fout:
        proc = subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, cwd=workdir)
    if proc.returncode != 0:
        raise RuntimeError(f"xTB(spin) failed for {xyz_path} UHF={uhf} (see {out_file})")

    energy_Eh, gap_eV = None, None
    with open(out_file, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    for ln in lines:
        if "HOMO-LUMO gap" in ln:
            m = re.search(r"HOMO-LUMO gap\s+([-\d\.Ee+]+)\s*eV", ln)
            if m: gap_eV = float(m.group(1)); break
    for ln in reversed(lines):
        if "TOTAL ENERGY" in ln.upper():
            m = re.search(r"TOTAL ENERGY\s+([-\d\.Ee+]+)\s*Eh", ln, re.I)
            if m: energy_Eh = float(m.group(1)); break
    return {"energy_Eh": energy_Eh, "gap_eV": gap_eV, "stdout_file": str(out_file)}

# ---------------- Pass 2: properties (standard xTB) ----------------
def run_xtb_props_standard(xyz_path, charge, uhf, workdir, acc=0.2, maxiter=500):
    ns = f"{Path(xyz_path).stem}_UHF{uhf}_PROP"
    cmd = [
        "xtb", xyz_path, "--scc",
        "--chrg", str(charge), "--uhf", str(uhf),
        "--acc", str(acc), "--iterations", str(maxiter),
        "--pop", "--wbo",
        "--namespace", ns
    ]
    out_file = Path(workdir) / f"{ns}.out"
    with open(out_file, "w") as fout:
        proc = subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, cwd=workdir)
    if proc.returncode != 0:
        raise RuntimeError(f"xTB(props) failed: {xyz_path} UHF={uhf} (see {out_file})")

    with open(out_file, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    charges = []
    in_charge = False
    for ln in lines:
        if re.match(r"\s*#\s+Z\s+.*\sq\s+", ln):
            in_charge = True; continue
        if in_charge:
            if not ln.strip(): in_charge = False; continue
            parts = ln.split()
            if len(parts) >= 5 and parts[0].isdigit():
                try: charges.append(float(parts[3]))
                except: pass
            else:
                in_charge = False

    # ---- NEW: parse WBO from the dedicated *.wbo file ----
    wbo_file = Path(out_file).with_suffix(".wbo")
    wbo_pairs = parse_wbo_file(wbo_file)
    return {"charges": charges, "wbo_pairs": wbo_pairs, "stdout_file": str(out_file)}

# ---------------- Driver ----------------
def main():
    ap = argparse.ArgumentParser(description="Two-pass spin scan (tblite energies) + properties (standard xTB).")
    ap.add_argument("--trj", required=True)
    ap.add_argument("--charge", type=int, required=True)
    ap.add_argument("--mult", type=int, default=-1, help="-1=auto from parity; else odd multiplicity (1,3,5,...)")
    ap.add_argument("--uhf-max", type=int, default=6)
    ap.add_argument("--acc", type=float, default=0.2)
    ap.add_argument("--iterations", type=int, default=500)
    ap.add_argument("--workdir", default="xTB-scan")
    ap.add_argument("--nprocs", type=int, default=1)

    ap.add_argument("--reactive-thresh", type=float, default=0.5,
                help="|Δ| threshold between first and last electron-BEM (>0.5 flagged)")
    ap.add_argument("--reactive-out", default="reactive_summary.csv",
                help="Output CSV summarizing time series of flagged entries")
    ap.add_argument("--reactive-include-meta", action="store_true",
                help="Include type/first/last/delta columns (default: omit)")

    args = ap.parse_args()

    workdir = Path(args.workdir).absolute(); workdir.mkdir(parents=True, exist_ok=True)
    xyzdir = workdir / "xyz_frames"; xyzdir.mkdir(exist_ok=True)
    bemdir = workdir / "bem_snapshots"; bemdir.mkdir(exist_ok=True)

    frames = read_multi_xyz_like_trj(args.trj)
    if not frames:
        print("No frames detected."); sys.exit(1)

    xyz_paths = []
    for idx, (symbols, coords, comment) in enumerate(frames):
        outp = xyzdir / f"step_{idx:04d}.xyz"
        write_xyz(symbols, coords, outp, comment or f"frame {idx}")
        xyz_paths.append((idx, outp, symbols))

    ladder = []
    for idx, xyzp, symbols in xyz_paths:
        if args.mult == -1:
            base = base_uhf_from_parity(symbols, args.charge)
        else:
            if args.mult < 1 or args.mult % 2 == 0:
                raise ValueError("--mult must be odd (1,3,5,...) or -1 for auto.")
            base = max(0, args.mult - 1)
        uhf_vals = sorted(set([u for u in range(base, args.uhf_max+1, 2) if u >= 0]))
        ladder.append((idx, xyzp, symbols, uhf_vals))

    # ---- Pass 1: spin-polarized energies (tblite)
    energy_results = {}  # idx -> {uhf: dict}
    jobs = []
    with ProcessPoolExecutor(max_workers=args.nprocs) as ex:
        futs = []
        keys = []
        for idx, xyzp, symbols, uhf_vals in ladder:
            for u in uhf_vals:
                futs.append(ex.submit(run_xtb_energy_tblite, str(xyzp), args.charge, u, str(workdir), args.acc, args.iterations))
                keys.append((idx,u))
        for (idx,u), fut in tqdm(zip(keys, as_completed(futs)), total=len(keys), desc="xTB spin energy"):
            try:
                res = fut.result()
            except Exception as e:
                print(f"[error] frame {idx} UHF={u} (spin energy): {e}")
                res = None
            energy_results.setdefault(idx, {})[u] = res

    # choose best UHF per frame
    best_by_frame = {}
    summary_rows = []
    for idx in sorted(energy_results.keys()):
        per = energy_results[idx]
        best_u, best_e = None, None
        row = OrderedDict(); row["frame"]=idx
        for u in sorted(per.keys()):
            e = per[u]["energy_Eh"] if per[u] else None
            row[f"E_UHF{u}_Eh"] = e
            if e is not None and (best_e is None or e < best_e):
                best_e, best_u = e, u
        row["best_UHF"]=best_u; row["best_E_Eh"]=best_e
        summary_rows.append(row)
        best_by_frame[idx] = best_u

    # ---- Pass 2: properties for best UHF (standard xTB)
    props_by_frame = {}
    with ProcessPoolExecutor(max_workers=args.nprocs) as ex:
        futs = []
        keys = []
        for idx, xyzp, symbols, _ in ladder:
            bu = best_by_frame.get(idx)
            if bu is None: continue
            futs.append(ex.submit(run_xtb_props_standard, str(xyzp), args.charge, bu, str(workdir), args.acc, args.iterations))
            keys.append((idx, bu, xyzp, symbols))
        for (idx, bu, xyzp, symbols), fut in tqdm(zip(keys, as_completed(futs)), total=len(keys), desc="xTB props"):
            try:
                res = fut.result()
            except Exception as e:
                print(f"[error] frame {idx} props UHF={bu}: {e}")
                res = None
            props_by_frame[idx] = {"best_uhf": bu, "props": res, "symbols": symbols}

    # collect global pair list
    observed_pairs = set()
    for idx, rec in props_by_frame.items():
        if rec["props"] and rec["props"]["wbo_pairs"]:
            observed_pairs.update(rec["props"]["wbo_pairs"].keys())
    pair_list = sorted(observed_pairs)
    pair_labels = [f"{i}-{j}" for (i,j) in pair_list]

    # write BEM-like JSONs and wbo_timeseries
    import csv
    bemdir.mkdir(exist_ok=True)
    ts_rows = []
    for row in summary_rows:
        idx = row["frame"]
        rec = props_by_frame.get(idx)
        if not rec or not rec["props"]:
            continue
        symbols = rec["symbols"]
        charges = rec["props"]["charges"]
        wbo = rec["props"]["wbo_pairs"]
        atoms = [{"id": k, "el": symbols[k],
                  "charge": float(charges[k]) if k < len(charges) else None}
                 for k in range(len(symbols))]
        bonds = [{"i": i, "j": j, "order": float(v)} for (i,j),v in sorted(wbo.items())]
        bem = {"schema":"BEM-lite:v1", "frame": idx,
               "best_uhf": int(rec["best_uhf"]), "energy_Eh": row["best_E_Eh"],
               "note":"Diagonal valences left null; off-diagonals=Wiberg/Mayer from standard xTB.",
               "atoms": atoms, "bonds": bonds}
        out_json = (Path(args.workdir) / "bem_snapshots" / f"step_{idx:04d}.json")
        with open(out_json, "w") as f:
            json.dump(bem, f, indent=2)

        # Also write BEM CSVs in the top-level workdir (to mirror your YARP CSV layout)
        write_bem_csvs(
            symbols=symbols,
            charges=charges,
            wbo_pairs=wbo,
            out_dir=Path(args.workdir),   # CSVs go next to spin_scan_summary.csv
            frame_idx=idx
        )

        ts = {"frame": idx, "best_UHF": int(rec["best_uhf"])}
        for (i,j), lab in zip(pair_list, pair_labels):
            ts[lab] = float(wbo.get((i,j), 0.0))
        ts_rows.append(ts)

    # write summaries
    sum_cols = sorted({k for r in summary_rows for k in r.keys()})
    with open(workdir / "spin_scan_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sum_cols)
        w.writeheader(); [w.writerow(r) for r in summary_rows]

    ts_cols = ["frame","best_UHF"] + pair_labels
    with open(workdir / "wbo_timeseries.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ts_cols)
        w.writeheader()
        for r in sorted(ts_rows, key=lambda x:x["frame"]):
            w.writerow(r)

    # naive crossover detection from best_UHF track
    xovers, prev = [], None
    for r in sorted(summary_rows, key=lambda x:x["frame"]):
        cur = r["best_UHF"]
        if prev is not None and cur is not None and prev != cur:
            xovers.append(r["frame"])
        prev = cur
    with open(workdir / "spin_crossover_frames.txt", "w") as f:
        for fr in xovers: f.write(f"{fr}\n")

    print(f"[done] Frames: {len(frames)}")
    print(f"[done] Summary: {workdir/'spin_scan_summary.csv'}")
    print(f"[done] WBO timeseries: {workdir/'wbo_timeseries.csv'}")
    print(f"[done] BEM-like JSONs: {bemdir}")
    print(f"[done] Spin crossovers (by best UHF change): {xovers}")

    # --- New: reactive elements summary from float-BEMs ---
    summarize_reactive_elements(
        workdir=args.workdir,
        thresh=args.reactive_thresh,
        out_csv=args.reactive_out,
        include_meta=args.reactive_include_meta
    )

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] != "--help":
        main()

