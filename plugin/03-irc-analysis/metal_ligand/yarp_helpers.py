from __future__ import annotations
import contextlib, os, tempfile, time
from pathlib import Path
import numpy as np
import yarp as yp

def silent_yarpecule(xyz_or_tuple, charge=None, **kwargs):
    """
    Accepts:
      - str path: "file.xyz"  (charge auto-parsed)
      - tuple: ("file.xyz", q)  (explicit total charge)
      - str path + charge kwarg: charge=-1
    Forwards any extra kwargs (e.g., canon=False) to yarpecule.
    """
    with open(os.devnull, "w") as dn, \
         contextlib.redirect_stdout(dn), \
         contextlib.redirect_stderr(dn):

        # explicit charge via kwarg
        if charge is not None:
            return yp.yarpecule((str(xyz_or_tuple), int(charge)), **kwargs)

        # tuple branch: (xyz_path, q)
        if isinstance(xyz_or_tuple, (tuple, list)) and len(xyz_or_tuple) == 2:
            return yp.yarpecule((str(xyz_or_tuple[0]), int(xyz_or_tuple[1])), **kwargs)

        # default: string path, auto-parse charge
        return yp.yarpecule(str(xyz_or_tuple), **kwargs)


def bond_state(adj_val, order_val):
    if adj_val == 0: return "none"
    return "dative bond" if order_val == 0 else "sigma bond"

def labelled_matrix(mat: np.ndarray, elements):
    labels = [f"{str(e).capitalize()}{i}" for i, e in enumerate(elements)]
    header = [""] + labels
    labelled = [header]
    for elem, row in zip(labels, mat):
        labelled.append([elem.capitalize()] + [int(v) for v in row])
    return labelled

