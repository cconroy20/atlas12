#!/usr/bin/env python3
"""Build the lookup index that lets SYNTHE interpolate the extracted NLTE grid
directly, with no per-model sidecar.

WHY AN INDEX RATHER THAN A DENSE TABLE
--------------------------------------
The store from nlte_extract_grid.py is a flat array of records; what it lacks
is a way to get from stellar parameters to a record number.  Resampling it
onto a dense rectangular table would be the obvious fix, but it bakes in the
axis choices of whatever grid is current -- and the requirement here is to
build models at ANY Teff, log g, metallicity, Na abundance and microturbulence.
So the sparse store is kept as it is and this writes a small integer index over
the full axis product: cell -> record number, or 0 where the grid has nothing.
SYNTHE then binary-searches each axis and reads only the corner records it
needs (~79 kB), never the 1 GB.

THE AXES, AND ONE THAT ISN'T
----------------------------
    Teff, log g, [Fe/H], v_turb, dNa          five axes
    dNa = A(Na) - [Fe/H] is a fixed 31-value ladder, identical for every model,
    so Na abundance is a proper interpolation axis (about +-1.5 dex on solar).

[alpha/Fe] is NOT an axis.  These are MARCS '_st_' standard-composition models,
in which alpha is a deterministic function of [Fe/H] (0 at [Fe/H] >= 0, ramping
to +0.4 by [Fe/H] <= -1): 15 ([Fe/H], [alpha/Fe]) pairs exist out of 75
possible.  An alpha-enhanced model therefore receives b computed at whatever
alpha the standard relation assigns to its metallicity.  Na's departures depend
on alpha only indirectly, through T(tau) and n_e, so this is expected to be
minor -- but it is a constraint the grid imposes, not a choice, and it is
recorded in the output metadata so it cannot be forgotten.

GEOMETRY
--------
MARCS is plane-parallel only at log g >= 3 and spherical below; ATLAS12 is
plane-parallel throughout.  Policy here: use plane-parallel where it exists,
otherwise spherical at M = 1 Msun.  Every cell carries a geometry flag, because
the switch falls exactly on the log g 2.5/3.0 boundary and interpolating across
it mixes the two -- a seam that has to be measurable rather than assumed small.

OUTPUT
------
  <out>.idx   little-endian binary, readable directly by Fortran:
                 int32   nT, nG, nF, nV, nD
                 float32 T(nT), G(nG), F(nF), V(nV), D(nD)
                 int32   rec(nT,nG,nF,nV,nD)   1-based record, 0 = absent
                 int8    geom(...)             0 absent, 1 plane-par, 2 spherical
              Fortran order on the index arrays, so rec(i,j,k,l,m) reads
              naturally with the first axis fastest.
  <out>.idx.json  axis values, fill statistics and the policy record.

USAGE
  python3 tools/nlte_build_index.py ~/kurucz/nlte_grids/na/Na_MARCS_Jul-14-2023
"""
import argparse
import json
import re
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="store prefix (expects .npz and .f32)")
    ap.add_argument("--mass", type=float, default=1.0,
                    help="spherical model mass to use, Msun (default 1.0)")
    a = ap.parse_args()

    z = np.load(a.prefix + ".npz", allow_pickle=True)
    tf, lg, fe, vt, ab = (z[k] for k in ("teff", "logg", "feh", "vturb", "abund"))
    mod = z["model"]
    geom = np.array([m[0] for m in mod])
    mass = np.array([float(re.search(r"_m([0-9.]+)_", m).group(1)) for m in mod])
    # dNa is the ladder offset; rounded because it is a difference of two
    # values printed to 2 dp in the index file.
    dna = np.round(ab - fe, 3)

    sel = (geom == "p") | ((geom == "s") & (mass == a.mass))
    if not sel.any():
        sys.exit(f"ERROR: no records with geometry p or (s, M={a.mass})")

    T, G, F, V, D = (np.unique(x[sel]) for x in (tf, lg, fe, vt, dna))
    shape = (len(T), len(G), len(F), len(V), len(D))
    rec = np.zeros(shape, np.int32)
    gm = np.zeros(shape, np.int8)

    ix = lambda u, x: int(np.searchsorted(u, x))
    nover = 0
    for i in np.where(sel)[0]:
        c = (ix(T, tf[i]), ix(G, lg[i]), ix(F, fe[i]), ix(V, vt[i]), ix(D, dna[i]))
        if rec[c] == 0:
            rec[c] = i + 1
            gm[c] = 1 if geom[i] == "p" else 2
        else:
            nover += 1
            # plane-parallel wins where both geometries exist
            if geom[i] == "p" and gm[c] == 2:
                rec[c] = i + 1
                gm[c] = 1

    filled = int((rec > 0).sum())
    npp = int((gm == 1).sum())
    nsp = int((gm == 2).sum())

    with open(a.prefix + ".idx", "wb") as f:
        np.array(shape, dtype="<i4").tofile(f)
        for u in (T, G, F, V, D):
            u.astype("<f4").tofile(f)
        # Fortran order: first axis fastest, so SYNTHE indexes rec(i,j,k,l,m)
        f.write(np.asfortranarray(rec).astype("<i4").tobytes(order="F"))
        f.write(np.asfortranarray(gm).astype("<i1").tobytes(order="F"))

    json.dump({
        "shape": list(shape),
        "axes": {"teff": T.tolist(), "logg": G.tolist(), "feh": F.tolist(),
                 "vturb": V.tolist(), "dNa": D.tolist()},
        "cells": int(rec.size), "filled": filled,
        "fill_fraction": filled / rec.size,
        "plane_parallel": npp, "spherical": nsp,
        "spherical_mass_Msun": a.mass,
        "geometry_policy": "plane-parallel where present, else spherical at "
                           f"M={a.mass}; flag in geom array",
        "alpha_note": "[alpha/Fe] is NOT an axis: MARCS _st_ models tie it to "
                      "[Fe/H] by the standard relation (15 pairs of 75). An "
                      "alpha-enhanced model gets b at the standard alpha for "
                      "its [Fe/H].",
        "dNa_definition": "A(Na) - [Fe/H]; a fixed ladder common to all models",
        "overlapping_records_resolved": nover,
    }, open(a.prefix + ".idx.json", "w"), indent=2)

    print(f"axes      : {len(T)} Teff x {len(G)} logg x {len(F)} feh x "
          f"{len(V)} vturb x {len(D)} dNa = {rec.size} cells")
    print(f"filled    : {filled} ({100*filled/rec.size:.1f} %)  "
          f"{npp} plane-parallel, {nsp} spherical M={a.mass}")
    print(f"overlaps  : {nover} resolved in favour of plane-parallel")
    print(f"wrote     : {a.prefix}.idx  "
          f"({(20 + 4*sum(shape) + 4*rec.size + rec.size)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
