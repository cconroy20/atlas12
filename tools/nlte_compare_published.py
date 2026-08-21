#!/usr/bin/env python3
"""Compare our NLTE line results with the published EW grids of
Lind et al. (2022, A&A 665, A33; VizieR J/A+A/665/A33), which tabulate
log10 EW in LTE and in NLTE for 35 Na and 134 Mg lines over 2158 MARCS
models -- the same quantity we compute, from the same lineage of model
atoms, applied by a different code (PySME) on different atmospheres.

WHY ABSOLUTE dEW IS THE RIGHT THING TO COMPARE, and fractional dEW is not:
our equivalent widths are integrated over a fixed window of the full
synthetic spectrum, so they include every blend in that window and are much
larger than the isolated line's EW.  But only the TAGGED line changes between
our LTE and NLTE runs, so dEW = EW_NLTE - EW_LTE is a property of that line
alone, in mA, and is directly comparable with theirs.  Dividing by our
blended EW_LTE would not be.

The abundance correction is the other robust comparison.  Ours comes from
-dEW / (dEW/dlogA) with the slope measured by a second synthesis at +0.1 dex;
theirs is computed here the standard way from their own curve of growth --
find the LTE abundance that reproduces the NLTE equivalent width:

    A_NLTE - A_LTE  =  A0 - x   where  EW_LTE(x) = EW_NLTE(A0)

  python3 tools/nlte_compare_published.py --data <dir with mg.dat/na.dat>
"""
import argparse
import os

import numpy as np

D = "workdir/nlte_lines"

# our line -> (published line id, element, window, vacuum centre, EW half-width)
LINES = [
    ("Mg I b4 5167", "5167", "Mg", "Mgb", 5168.761, 1.0),
    ("Mg I b2 5173", "5173", "Mg", "Mgb", 5174.125, 1.0),
    ("Mg I b1 5184", "5184", "Mg", "Mgb", 5185.048, 1.0),
    ("Na I D2 5890", "5890", "Na", "NaD", 5891.583, 1.5),
    ("Na I D1 5896", "5896", "Na", "NaD", 5897.558, 1.5),
]
# our models -> the (Teff, logg, [Fe/H], vturb) to look up in their grid
MODELS = {
    "sun": (5777.0, 4.44, 0.00, 2.0, "Sun 5777/4.44/0.0"),
    "mp2": (5750.0, 4.50, -2.00, 2.0, "metal-poor 5750/4.5/-2.0"),
}
XFE = 0.0          # our models are scaled-solar in Mg and Na


def read_grid(path, lineid):
    rows = []
    with open(path) as f:
        for ln in f:
            t = ln.split()
            if t[1] != lineid:
                continue
            rows.append((float(t[2]), float(t[3]), float(t[4]), float(t[5]),
                         float(t[6]), float(t[7]), float(t[8])))
    return np.array(rows)


def cog_correction(rows, teff, logg, feh, vturb, xfe=XFE):
    """Their NLTE abundance correction at these parameters, from their own
    curve of growth, with bilinear interpolation in Teff and log g onto the
    nearest bracketing nodes (their [Fe/H] and vturb match ours exactly)."""
    m = (np.abs(rows[:, 2] - feh) < 1e-6) & (np.abs(rows[:, 3] - vturb) < 1e-6)
    r = rows[m]
    if not len(r):
        return None
    Ts = np.unique(r[:, 0]); Gs = np.unique(r[:, 1])
    def brack(ax, x):
        if x <= ax[0]:  return ax[0], ax[0], 0.0
        if x >= ax[-1]: return ax[-1], ax[-1], 0.0
        i = np.searchsorted(ax, x) - 1
        return ax[i], ax[i+1], (x - ax[i])/(ax[i+1] - ax[i])
    t0, t1, wt = brack(Ts, teff)
    g0, g1, wg = brack(Gs, logg)
    acc_d, acc_dew, acc_ewl, wsum = 0.0, 0.0, 0.0, 0.0
    for tv, wtv in ((t0, 1-wt), (t1, wt)):
        for gv, wgv in ((g0, 1-wg), (g1, wg)):
            w = wtv*wgv
            if w <= 0:
                continue
            c = r[(np.abs(r[:, 0]-tv) < 1e-6) & (np.abs(r[:, 1]-gv) < 1e-6)]
            if not len(c):
                return None
            o = np.argsort(c[:, 4])
            xf, lwl, lwn = c[o, 4], c[o, 5], c[o, 6]
            j = np.argmin(np.abs(xf - xfe))
            ewl, ewn = 10**lwl[j], 10**lwn[j]
            # LTE abundance reproducing the NLTE equivalent width
            x = np.interp(lwn[j], lwl, xf)
            acc_d += w*(xf[j] - x); acc_dew += w*(ewn - ewl); acc_ewl += w*ewl
            wsum += w
    return acc_d/wsum, acc_dew/wsum, acc_ewl/wsum


def spec(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]/np.maximum(d[:, 2], 1e-300)


def ew(w, r, wc, half):
    m = np.abs(w - wc) <= half
    return np.trapz(1.0 - r[m], w[m])*1e3


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="dir holding mg.dat, na.dat")
    a = ap.parse_args()

    cache = {}
    print(f"{'model':10s}{'line':15s}{'dEW ours':>10s}{'dEW Lind':>10s}"
          f"{'  ':2s}{'dlogA ours':>11s}{'dlogA Lind':>11s}   {'EW_L Lind':>10s}")
    for mkey, (teff, logg, feh, vturb, mlab) in MODELS.items():
        for lab, lid, el, win, wc, half in LINES:
            f0 = f"{D}/spec_{mkey}_{win}_lte.txt"
            f1 = f"{D}/spec_{mkey}_{win}_nlte.txt"
            fp = f"{D}/cog/spec_{mkey}_{win}_{el}p1.txt"
            if not all(os.path.exists(x) for x in (f0, f1, fp)):
                continue
            w, r0 = spec(f0); _, r1 = spec(f1); _, rp = spec(fp)
            e0, e1, ep = (ew(w, x, wc, half) for x in (r0, r1, rp))
            slope = (ep - e0)/0.1
            ours_dew = e1 - e0
            ours_d = -ours_dew/slope if slope > 5 else np.nan
            key = (el, lid)
            if key not in cache:
                cache[key] = read_grid(os.path.join(a.data,
                                       "mg.dat" if el == "Mg" else "na.dat"), lid)
            got = cog_correction(cache[key], teff, logg, feh, vturb)
            if got is None:
                print(f"{mkey:10s}{lab:15s}  (not in the published grid)")
                continue
            th_d, th_dew, th_ewl = got
            print(f"{mkey:10s}{lab:15s}{ours_dew:10.2f}{th_dew:10.2f}"
                  f"{'':2s}{ours_d:+11.3f}{th_d:+11.3f}   {th_ewl:10.1f}")


if __name__ == "__main__":
    main()
