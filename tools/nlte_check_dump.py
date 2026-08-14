#!/usr/bin/env python3
"""Verify a <model>.nlte dump against the profile it was supposed to come from.

This is the test the constant null test cannot do.  With b identical at every
depth and every transition, a wrong (transition, depth) index, a reversed depth
axis, or coefficients attributed to the wrong line all produce exactly the right
answer.  Here b varies along both axes, so the dump can be checked element by
element against an independent evaluation of the same profile.

Two things are checked, and they fail independently:

  1. INSTALLED b.  Recompute the structured-test ramp from the dumped column
     mass and compare against the dumped b_lo/b_up.  This catches fill_departures
     writing to the wrong slot, and a depth axis that runs the wrong way.

  2. THE FACTORS.  Recompute fkappa, fdev and S_l/B_nu from the dumped b and x
     using the equations in mod_nlte, and compare against what nlte_factors
     actually returned.  Because the dump records the returned values rather
     than re-deriving them, this catches an indexing error inside nlte_factors
     itself -- b read from one depth while the opacity is applied at another.

Usage:  python3 tools/nlte_check_dump.py MODEL.nlte
"""
import sys
import numpy as np

# Must match BLO_DEEP/BLO_TOP/BUP_DEEP/BUP_TOP in mod_nlte.f90.
BLO_DEEP = np.array([1.00, 1.00])
BLO_TOP  = np.array([1.20, 1.10])
BUP_DEEP = np.array([1.00, 1.00])
BUP_TOP  = np.array([0.60, 0.80])

# The dump is written at full double precision, so the only floor here is
# the round trip through decimal text: ~1e-14 relative.
TOL = 1e-10

d = np.loadtxt(sys.argv[1])
k, j = d[:, 0].astype(int), d[:, 1].astype(int)
rhox, T = d[:, 2], d[:, 3]
blo, bup, x = d[:, 4], d[:, 5], d[:, 6]
fkap, fdev, sratio = d[:, 7], d[:, 8], d[:, 9]
seen, inv = d[:, 10].astype(int), d[:, 11].astype(int)

ntr, nj = k.max(), j.max()
print(f"{len(d)} rows: {ntr} transitions x {nj} depths")
print(f"column mass {rhox.min():.3e} -> {rhox.max():.3e} g/cm2, "
      f"T {T.min():.1f} -> {T.max():.1f} K")
print(f"seen at all pairs: {bool((seen == 1).all())}   inversions: {(inv == 1).sum()}")

# --- 1. installed b against an independent evaluation of the ramp ------------
lr = np.log10(rhox)
lr_top, lr_bot = lr[j == 1][0], lr[j == nj][0]
f = (lr_bot - lr) / (lr_bot - lr_top)
blo_exp = BLO_DEEP[k - 1] + (BLO_TOP[k - 1] - BLO_DEEP[k - 1]) * f
bup_exp = BUP_DEEP[k - 1] + (BUP_TOP[k - 1] - BUP_DEEP[k - 1]) * f
e1 = np.abs(blo - blo_exp).max()
e2 = np.abs(bup - bup_exp).max()
print(f"\n1. installed b vs analytic ramp")
print(f"   max |db_lo| = {e1:.3e}   max |db_up| = {e2:.3e}"
      f"   [{'PASS' if max(e1, e2) < TOL else 'FAIL'}]")
print(f"   surface (j=1): b_lo={blo[j==1]}, b_up={bup[j==1]}  -> want top values")
print(f"   base  (j={nj}): b_lo={blo[j==nj]}, b_up={bup[j==nj]}  -> want deep values")

# --- 2. returned factors against the mod_nlte equations ---------------------
ex = np.exp(-x)
fkap_exp = blo * (1.0 - (bup / blo) * ex) / (1.0 - ex)
fdev_exp = (bup - fkap_exp) / fkap_exp
srat_exp = (1.0 - ex) / (blo / bup - ex)
e3 = np.abs(fkap - fkap_exp).max()
e4 = np.abs(fdev - fdev_exp).max()
e5 = np.abs(sratio - srat_exp).max()
print(f"\n2. returned factors vs the equations")
print(f"   max |dfkappa| = {e3:.3e}   max |dfdev| = {e4:.3e}"
      f"   max |dS/B| = {e5:.3e}   [{'PASS' if max(e3,e4,e5) < TOL else 'FAIL'}]")

# The identity the whole design rests on: emissivity = b_up x LTE, i.e.
# fkappa * (1 + fdev) must equal b_up exactly.
e6 = np.abs(fkap * (1.0 + fdev) - bup).max()
print(f"\n3. identity  fkappa*(1+fdev) = b_up   max residual = {e6:.3e}"
      f"   [{'PASS' if e6 < TOL else 'FAIL'}]")

# --- 4. the two transitions must be distinguishable -------------------------
if ntr >= 2:
    a = sratio[k == 1]
    b = sratio[k == 2]
    print(f"\n4. transitions distinguishable: max |S/B(D2) - S/B(D1)| = "
          f"{np.abs(a - b).max():.4f}   [{'PASS' if np.abs(a-b).max() > 0.01 else 'FAIL'}]")
    print(f"   D2 S/B: {a.min():.4f} (surface) -> {a.max():.4f} (base)")
    print(f"   D1 S/B: {b.min():.4f} (surface) -> {b.max():.4f} (base)")
