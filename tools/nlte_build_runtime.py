#!/usr/bin/env python3
"""Build the single self-contained runtime file SYNTHE reads for NLTE_MODE = 3.

WHY A DERIVED FILE RATHER THAN THE EXTRACTED STORE
--------------------------------------------------
tools/nlte_extract_grid.py produces the master store: every node the published
grid contains, with as many low-lying levels as were worth keeping while the
15.9 GB download was open.  That is the archive, and it belongs outside the
repository.  Most of it is unreachable at run time:

only records the index points at are reachable, and the index fixes the
spherical-model mass at 1 Msun, so the 0.5/2/5/15 Msun models are ballast --
27 per cent of records are never referenced.

Dropping them, with the index embedded, gives ONE file to install, no way to
pair a store with a stale index, and no dead weight.  ALL TEN LEVELS ARE KEPT
BY DEFAULT: they cost 500 MB over the three Na D needs, but they are what
makes 8183/8195 (3p-3d) and 1.14 um (3p-4s) reachable later without touching
the grid again.  Narrowing is a local rebuild; widening past ten would need
another 15.9 GB download, so the asymmetry argues for keeping them.

FILE FORMAT (little-endian, read directly by mod_nlte)
-----------------------------------------------------
    int32    nT, nG, nF, nV, nD          axis lengths
    int32    nlev, ndep, nrec            levels kept, depths, records
    float32  T(nT), G(nG), F(nF), V(nV), D(nD)
    int32    level_id(nlev)              original atom level numbers
    int32    rec(nT,nG,nF,nV,nD)         1-based record, 0 = absent, Fortran order
    int8     geom(...)                   0 absent, 1 plane-parallel, 2 spherical
    float32  data(nrec, ndep + nlev*ndep)   tau then b(level, depth), C order

USAGE
  python3 tools/nlte_build_runtime.py MASTER_PREFIX OUT.nlte --levels 1,2,3
"""
import argparse
import os
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master", help="master store prefix (.f32 + .idx)")
    ap.add_argument("out", help="output runtime file")
    ap.add_argument("--levels", default="1,2,3,4,5,6,7,8,9,10",
                    help="atom level numbers to keep, 1-based.  Default is "
                         "all ten from the master: Na D needs only 1,2,3, but "
                         "keeping the rest preserves 8183/8195 and 1.14 um "
                         "without re-downloading the grid.")
    ap.add_argument("--ndep", type=int, default=56)
    ap.add_argument("--nlev-master", type=int, default=10)
    a = ap.parse_args()

    want = [int(v) for v in a.levels.split(",")]
    if any(v < 1 or v > a.nlev_master for v in want):
        sys.exit(f"ERROR: levels must be within 1..{a.nlev_master}")

    # --- read the master index -------------------------------------------
    with open(a.master + ".idx", "rb") as f:
        nT, nG, nF, nV, nD = np.fromfile(f, "<i4", 5)
        T = np.fromfile(f, "<f4", nT); G = np.fromfile(f, "<f4", nG)
        F = np.fromfile(f, "<f4", nF); V = np.fromfile(f, "<f4", nV)
        D = np.fromfile(f, "<f4", nD)
        n = nT*nG*nF*nV*nD
        rec = np.fromfile(f, "<i4", n).reshape((nT,nG,nF,nV,nD), order="F")
        geo = np.fromfile(f, "<i1", n).reshape((nT,nG,nF,nV,nD), order="F")

    used = np.unique(rec[rec > 0])
    nrec = len(used)
    print(f"master index : {n} cells, {nrec} distinct records referenced")

    # renumber: master record -> 1-based slot in the new data block
    remap = np.zeros(used.max() + 1, np.int32)
    remap[used] = np.arange(1, nrec + 1, dtype=np.int32)
    newrec = np.where(rec > 0, remap[np.clip(rec, 0, None)], 0).astype(np.int32)

    # --- copy the wanted levels of the wanted records ---------------------
    W = a.ndep + a.nlev_master*a.ndep
    master = np.memmap(a.master + ".f32", dtype="<f4", mode="r")
    if master.size % W:
        sys.exit(f"ERROR: {a.master}.f32 is not a whole number of records")
    master = master.reshape(-1, W)
    print(f"master store : {master.shape[0]} records x {W} floats")

    nlev = len(want)
    out = np.empty((nrec, a.ndep + nlev*a.ndep), dtype="<f4")
    for slot, r in enumerate(used):
        row = master[r-1]
        out[slot, :a.ndep] = row[:a.ndep]
        for li, lev in enumerate(want):
            s = a.ndep + (lev-1)*a.ndep
            out[slot, a.ndep + li*a.ndep : a.ndep + (li+1)*a.ndep] = row[s:s+a.ndep]

    if not np.isfinite(out).all():
        sys.exit("ERROR: non-finite values in the selected levels")
    if (out[:, a.ndep:] <= 0).any():
        sys.exit("ERROR: non-positive b in the selected levels")

    # --- write ------------------------------------------------------------
    with open(a.out, "wb") as f:
        np.array([nT,nG,nF,nV,nD], "<i4").tofile(f)
        np.array([nlev, a.ndep, nrec], "<i4").tofile(f)
        for u in (T,G,F,V,D):
            u.astype("<f4").tofile(f)
        np.array(want, "<i4").tofile(f)
        f.write(np.asfortranarray(newrec).astype("<i4").tobytes(order="F"))
        f.write(np.asfortranarray(geo).astype("<i1").tobytes(order="F"))
        out.tofile(f)

    sz = os.path.getsize(a.out)
    print(f"levels kept  : {want}")
    print(f"wrote        : {a.out}  {sz/1e6:.0f} MB "
          f"(master was {os.path.getsize(a.master + '.f32')/1e6:.0f} MB "
          f"+ {os.path.getsize(a.master + '.idx')/1e6:.0f} MB)")


if __name__ == "__main__":
    main()
