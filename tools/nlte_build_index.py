#!/usr/bin/env python3
"""Build the lookup index that lets SYNTHE interpolate the extracted NLTE grid
directly, with no per-model sidecar.

WHY AN INDEX RATHER THAN A DENSE TABLE
--------------------------------------
The store from nlte_extract_grid.py is a flat array of records; what it lacks
is a way to get from stellar parameters to a record number.  Resampling it
onto a dense rectangular table would be the obvious fix, but it bakes in the
axis choices of whatever grid is current -- and the requirement here is to
build models at ANY Teff, log g, metallicity, element abundance and vturb.
So the sparse store is kept as it is and this writes a small integer index over
the full axis product: cell -> record number, or 0 where the grid has nothing.
SYNTHE then binary-searches each axis and reads only the corner records it
needs (~79 kB), never the 1 GB.

THE AXES, AND ONE THAT ISN'T
----------------------------
    Teff, log g, [Fe/H], v_turb, dX           five axes
    dX = A(X) - [Fe/H] is a fixed ladder (31 values for Na, 21 for Mg and Ca),
    identical for every model, so the element's abundance is a proper
    interpolation axis -- roughly +-1.5 dex about the scaled-solar value.

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
    # Geometry letter and mass come from the MARCS model name.  Not every
    # record has one: the Mg grid carries 21 records named simply "sun",
    # which are extra models rather than nodes of the rectangular grid --
    # and whose leading "s" would otherwise be read as "spherical".  Anything
    # that does not match the MARCS pattern is dropped, with a count.
    pat = re.compile(r"^([ps])\d+_g[+-][\d.]+_m([0-9.]+)_")
    mm = [pat.match(m) for m in mod]
    named = np.array([x is not None for x in mm])
    geom = np.array([x.group(1) if x else "?" for x in mm])
    mass = np.array([float(x.group(2)) if x else -1.0 for x in mm])
    # dX = A(X) - [Fe/H] is the ladder offset; rounded because it is a
    # difference of two values printed to 2 dp in the index file.
    dna = np.round(ab - fe, 3)

    sel = named & ((geom == "p") | ((geom == "s") & (mass == a.mass)))
    if (~named).any():
        print(f"unnamed   : {int((~named).sum())} records are not MARCS grid "
              f"nodes (e.g. {mod[~named][0]!r}) and are excluded")
    if not sel.any():
        sys.exit(f"ERROR: no records with geometry p or (s, M={a.mass})")

    # DOES A(X) VARY AT ALL AT FIXED ATMOSPHERE?  For a trace species the
    # grid carries a ladder of 21-31 abundances per atmosphere and dX is a
    # genuine fifth axis.  For two elements it is not:
    #
    #   H   A(H) = 12.000 always -- it is the definition of the scale, not a
    #       free parameter, so dX = 12 - [Fe/H] takes 15 values that are
    #       PERFECTLY DEGENERATE with the [Fe/H] axis.  Left alone that builds
    #       a 15x15 ([Fe/H], dX) plane with only the anti-diagonal filled, and
    #       a model landing between two [Fe/H] nodes then finds 2 of its 4
    #       corners missing, drops half the interpolation weight and
    #       renormalises -- a silent collapse to something far cruder than an
    #       interpolation in metallicity.
    #   Fe  A(Fe) = [Fe/H] + 7.500 exactly, since in a scaled-solar MARCS
    #       model the iron abundance IS the metallicity.  dX is single-valued
    #       already, so this costs nothing but is recorded for the same reason.
    #
    # Detected from the data rather than hard-coded per element: if no
    # atmosphere carries more than one A(X), the axis is collapsed to length 1
    # and the runtime lookup lands on it unconditionally.
    atm = np.stack([tf[sel], lg[sel], fe[sel], z["alpha"][sel], vt[sel]], 1)
    nkey = len(np.unique(atm, axis=0))
    nkeyab = len(np.unique(np.column_stack([atm, dna[sel]]), axis=0))
    varies = nkeyab > nkey
    if not varies:
        dfix = float(np.median(dna[sel]))
        print(f"abundance : A(X) does not vary at fixed atmosphere; the dX "
              f"axis is COLLAPSED to length 1 rather than left degenerate "
              f"with [Fe/H].  The stored value {dfix:+.3f} is a placeholder "
              f"-- a length-1 axis always brackets to itself, so it is never "
              f"compared against anything")
        dna = np.full_like(dna, dfix)

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
                 "vturb": V.tolist(), "dX": D.tolist()},
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
        "dX_definition": "A(X) - [Fe/H]; a fixed ladder common to all models",
        "dX_is_an_axis": bool(varies),
        "dX_note": ("A(X) varies at fixed atmosphere, so dX is a genuine fifth "
                    "axis" if varies else
                    "A(X) does not vary at fixed atmosphere (H: fixed at 12; "
                    "Fe: equal to [Fe/H] by construction), so the dX axis was "
                    "collapsed to length 1; leaving it would make it "
                    "degenerate with [Fe/H] and lose half the interpolation "
                    "weight between metallicity nodes"),
        "overlapping_records_resolved": nover,
    }, open(a.prefix + ".idx.json", "w"), indent=2)

    print(f"axes      : {len(T)} Teff x {len(G)} logg x {len(F)} feh x "
          f"{len(V)} vturb x {len(D)} dX = {rec.size} cells")
    print(f"filled    : {filled} ({100*filled/rec.size:.1f} %)  "
          f"{npp} plane-parallel, {nsp} spherical M={a.mass}")
    print(f"overlaps  : {nover} resolved in favour of plane-parallel")
    print(f"wrote     : {a.prefix}.idx  "
          f"({(20 + 4*sum(shape) + 4*rec.size + rec.size)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
