#!/usr/bin/env python3
"""Audit the continuum-opacity budget from a SYNTHE <model>.cont dump.

The .cont file (written when synthe is run with more_output=yes) carries the
per-source continuum opacity in cm^2/g at every continuum frequency and every
depth, exactly as KAPP assembled it.  This script integrates the continuum
optical depth down the column, finds the layer where tau_cont = 1 at each
wavelength, and reports which absorber carries the opacity there.

That is the question that matters for a model-vs-model continuum discrepancy:
not "what is the opacity somewhere" but "what sets the flux at the depth the
photons actually escape from".

Usage:
    python3 tools/cont_audit.py <model>.cont <model>.atm [--windows]
"""
import os
import sys

import numpy as np

SOURCES = ["aHyd", "aH2plus", "aHmin", "aH2min", "aH2coll",
           "aHe1", "aHe2", "aHemin", "aMetal"]
SCATTER = ["sigH", "sigHe", "sigEl", "sigH2", "sigX"]
COLS = (["nu", "wave_nm", "j", "T", "rho"] + SOURCES
        + ["ACONT"] + SCATTER + ["SIGMAC"])

WINDOWS = [
    ("optical 5500-5700", 550., 570.),
    ("TiO gamma 7053-7270", 705.3, 727.0),
    ("I  9000-9500", 900., 950.),
    ("J  11500-13500", 1150., 1350.),
    ("H  15500-17000", 1550., 1700.),
    ("K  21000-23000", 2100., 2300.),
]


def read_atm_rhox(path):
    rows = []
    started = False
    for line in open(path):
        if line.startswith("READ DECK"):
            started = True
            continue
        if started:
            p = line.split()
            if len(p) < 11:
                break
            rows.append(float(p[0]))
    return np.array(rows)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    cont, atm = sys.argv[1], sys.argv[2]
    rhox = read_atm_rhox(atm)

    # Fortran E-format drops the 'E' when the exponent needs three digits
    # ("9.28485-275").  The writer now uses E13.5E3, but stay tolerant of
    # older dumps rather than dying on one underflowed He II entry.
    def fix(tok):
        t = tok.decode() if isinstance(tok, bytes) else tok
        if "E" not in t and "e" not in t:
            i = max(t.rfind("+"), t.rfind("-"))
            if i > 0:
                t = t[:i] + "E" + t[i:]
        return float(t)

    d = np.loadtxt(cont, converters={i: fix for i in range(len(COLS))},
                   encoding=None)
    idx = {c: i for i, c in enumerate(COLS)}
    wl = d[:, idx["wave_nm"]]
    jj = d[:, idx["j"]].astype(int)
    nlay = jj.max()
    if nlay != len(rhox):
        print(f"warning: .cont has {nlay} layers, .atm has {len(rhox)}")
    nlay = min(nlay, len(rhox))

    nus = np.unique(d[:, idx["nu"]].astype(int))
    print(f"{os.path.basename(cont)}: {len(nus)} continuum points "
          f"{wl.min():.2f}-{wl.max():.1f} nm, {nlay} layers\n")

    # ---- per-wavelength: integrate tau_cont, locate tau=1, decompose there
    recs = {}
    nucol = d[:, idx["nu"]].astype(int)
    for nu in nus:
        m = (nucol == nu) & (jj <= nlay)
        if m.sum() < nlay:
            continue
        o = np.argsort(jj[m])
        blk = d[m][o]
        kap = blk[:, idx["ACONT"]] + blk[:, idx["SIGMAC"]]
        # tau = integral kappa dm, trapezoid on the column-mass grid
        tau = np.concatenate([[0.0], np.cumsum(
            0.5 * (kap[1:] + kap[:-1]) * np.diff(rhox[:nlay]))])
        if tau[-1] < 1.0:                      # never reaches tau=1
            j1 = nlay - 1
        else:
            j1 = int(np.searchsorted(tau, 1.0))
            j1 = min(j1, nlay - 1)
        tot = kap[j1]
        frac = {s: blk[j1, idx[s]] / tot for s in SOURCES + SCATTER}
        recs[blk[0, idx["wave_nm"]]] = (j1, blk[j1, idx["T"]], tot, frac)

    # ---- window medians
    print(f"{'window':22s} {'T(tau=1)':>9} {'kappa':>10}   "
          + " ".join(f"{s.replace('a','').replace('sig','s.'):>8s}"
                     for s in SOURCES + SCATTER))
    print("-" * 150)
    for lab, lo, hi in WINDOWS:
        ws = [w for w in recs if lo <= w <= hi]
        if not ws:
            continue
        T = np.median([recs[w][1] for w in ws if w in recs])
        kp = np.median([recs[w][2] for w in ws if w in recs])
        row = f"{lab:22s} {T:9.0f} {kp:10.3e}   "
        for s in SOURCES + SCATTER:
            f = np.median([recs[w][3][s] for w in ws if w in recs])
            row += f"{100 * f:8.2f}"
        print(row)
    print("\n(columns are percent of total continuum opacity at the "
          "tau_cont = 1 layer)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
