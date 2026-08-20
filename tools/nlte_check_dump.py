#!/usr/bin/env python3
"""Check a <model>.nlte diagnostic dump for internal consistency.

SYNTHE writes the dump whenever NLTE is on: one row per (transition, depth)
carrying the model's depth scale, the departure coefficients that were
installed, and the factors mod_nlte actually returned.  Because it records
what was RETURNED rather than re-deriving it, differencing it against the
equations catches an indexing error inside nlte_factors, not merely one in
the grid read.

Run this first when an NLTE run looks wrong: it separates "the plumbing is
broken" from "the grid data is not what I expected", which otherwise look
identical in a spectrum.

Usage:  python3 tools/nlte_check_dump.py MODEL.nlte
"""
import sys

import numpy as np

# Column layout of the dump; keep in step with nlte_dump in mod_nlte.f90.
K, J, RHOX, LGTAU5, T, BLO, BUP, X, FKAP, FDEV, SRAT, SEEN, INV, CLM = range(14)

d = np.loadtxt(sys.argv[1])
if d.shape[1] != 14:
    sys.exit(f"ERROR: expected 14 columns, got {d.shape[1]} -- the dump layout "
             "has changed and this script needs updating with it.")

ntr, nj = int(d[:, K].max()), int(d[:, J].max())
print(f"{len(d)} rows: {ntr} transitions x {nj} depths")
print(f"column mass {d[:, RHOX].min():.3e} .. {d[:, RHOX].max():.3e} g/cm2")
print(f"log tau_5000 {d[:, LGTAU5].min():.2f} .. {d[:, LGTAU5].max():.2f}")
print(f"T {d[:, T].min():.1f} .. {d[:, T].max():.1f} K")
print(f"seen at every pair: {bool((d[:, SEEN] == 1).all())}   "
      f"inversions: {int((d[:, INV] == 1).sum())}   "
      f"depth-clamped rows: {int((d[:, CLM] == 1).sum())}")

TOL = 1e-10          # the dump is full double precision; only decimal
                     # round-trip stands between these and exact
ok = True

# 1. the factors mod_nlte returned, against the equations they come from.
#    Rows forced back to LTE by the inversion guard are exempt by design.
live = d[:, INV] == 0
ex = np.exp(-d[live, X])
fkap_e = d[live, BLO] * (1.0 - (d[live, BUP] / d[live, BLO]) * ex) / (1.0 - ex)
fdev_e = (d[live, BUP] - fkap_e) / fkap_e
srat_e = (1.0 - ex) / (d[live, BLO] / d[live, BUP] - ex)
for name, got, want in (("fkappa", d[live, FKAP], fkap_e),
                        ("fdev", d[live, FDEV], fdev_e),
                        ("S_l/B_nu", d[live, SRAT], srat_e)):
    e = np.abs(got - want).max()
    ok &= e < TOL
    print(f"\n{name:<9} vs its equation : max |diff| = {e:.3e}  "
          f"[{'PASS' if e < TOL else 'FAIL'}]")

# 2. the identity the whole design rests on: emissivity = b_up x LTE.
e = np.abs(d[live, FKAP] * (1.0 + d[live, FDEV]) - d[live, BUP]).max()
ok &= e < TOL
print(f"identity  fkappa*(1+fdev) = b_up : max residual = {e:.3e}  "
      f"[{'PASS' if e < TOL else 'FAIL'}]")

# 3. rows the guard rejected must be exactly LTE, not approximately.
bad = d[~live]
if len(bad):
    clean = np.all(bad[:, FKAP] == 1.0) and np.all(bad[:, FDEV] == 0.0)
    ok &= clean
    print(f"{len(bad)} inverted rows forced to exact LTE : "
          f"[{'PASS' if clean else 'FAIL'}]")

# 4. the transitions must be distinguishable, or coefficients have been
#    attributed to the wrong line.
if ntr >= 2:
    a, b = d[d[:, K] == 1][:, SRAT], d[d[:, K] == 2][:, SRAT]
    if len(a) == len(b):
        print(f"\ntransitions distinguishable: max |S/B(1) - S/B(2)| = "
              f"{np.abs(a - b).max():.4f}")
        for t in (1, 2):
            s = d[d[:, K] == t][:, SRAT]
            print(f"   transition {t} S_l/B_nu: {s.min():.4f} .. {s.max():.4f}")

print("\nVERDICT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
