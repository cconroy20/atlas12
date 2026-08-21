#!/usr/bin/env python3
"""Two independent checks on the NLTE line results, for the Sun.

(1) NLTE ABUNDANCE CORRECTIONS, so the effect can be compared with the
    published literature, which is quoted in dex and not in per cent.  The
    local curve-of-growth slope is measured with the code itself, from an
    otherwise identical synthesis with A(X) raised by 0.1 dex at FIXED
    structure -- which is exactly what an abundance analysis does:

        dlogA = -(EW_NLTE - EW_LTE) / [(EW_+0.1 - EW_LTE)/0.1]

    i.e. the abundance change that would undo what NLTE did.  Sign convention
    is the usual one, A_NLTE - A_LTE.

(2) COMPARISON WITH THE OBSERVED SOLAR FLUX SPECTRUM (IAG atlas, Reiners
    et al. 2016, 405-1065 nm, normalised).  Equivalent width is the right
    observable here because it is invariant under convolution: the observed
    spectrum carries rotation and macroturbulence that the synthesis does
    not, so line DEPTHS cannot be compared without assuming a broadening
    kernel, while EW can.  The question asked is not "does the model match"
    -- a 1D LTE solar spectrum has larger errors than this from granulation,
    log gf and van der Waals -- but the differential one: does NLTE move the
    model TOWARD the observation or away from it?

  python3 tools/nlte_validate.py
"""
import argparse
import os
import numpy as np

D = "workdir/nlte_lines"
# line, window, vacuum centre, EW half-width [A], element
LINES = [
    ("Fe I 3931.4",    "HK",  3931.409, 0.6, "Fe"),
    ("Fe I 3924.0",    "HK",  3924.022, 0.6, "Fe"),
    ("Fe I 3921.4",    "HK",  3921.368, 0.6, "Fe"),
    ("Mg I b4 5167.7", "Mgb", 5168.761, 1.0, "Mg"),
    ("Mg I b2 5172.7", "Mgb", 5174.125, 1.0, "Mg"),
    ("Mg I b1 5183.6", "Mgb", 5185.048, 1.0, "Mg"),
    ("Fe I 5167.7",    "Mgb", 5167.721, 0.4, "Fe"),
    ("Fe I 5170.3",    "Mgb", 5170.337, 0.4, "Fe"),
    ("Ca II 8498",     "CaT", 8500.355, 4.0, "Ca"),
    ("Ca II 8542",     "CaT", 8544.451, 4.0, "Ca"),
    ("Ca II 8662",     "CaT", 8664.520, 4.0, "Ca"),
    ("Fe I 8517.4",    "CaT", 8517.449, 0.4, "Fe"),
    ("Fe I 8584.6",    "CaT", 8584.616, 0.4, "Fe"),
]


def spec(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]/np.maximum(d[:, 2], 1e-300)


def atlas(win):
    d = np.loadtxt(f"{D}/iag_{win}.txt")
    return d[:, 0], d[:, 1]


def ew(w, r, wc, half):
    m = np.abs(w - wc) <= half
    return np.trapz(1.0 - r[m], w[m])*1e3      # mA


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--specdir", default="workdir/nlte_lines")
    ap.add_argument("--model", default="sun",
                    help="which model's spectra to analyse: sun (5777/4.44, "
                         "[Fe/H]=0) or mp2 (5750/4.5, [Fe/H]=-2)")
    a = ap.parse_args()
    globals()["D"] = a.specdir
    print(__doc__.split("\n")[0])
    print(f"\n(1) NLTE ABUNDANCE CORRECTIONS ({a.model}), from our own curve "
          f"of growth\n")
    print(f"{'line':16s}{'EW_LTE':>9s}{'dEW_NLTE':>10s}{'dEW/dlogA':>11s}"
          f"{'dlogA_NLTE':>12s}")
    corr = {}
    for lab, win, wc, half, el in LINES:
        paths = [f"{D}/spec_{a.model}_{win}_lte.txt",
                 f"{D}/spec_{a.model}_{win}_nlte.txt",
                 f"{D}/cog/spec_{a.model}_{win}_{el}p1.txt"]
        if not all(os.path.exists(x) for x in paths):
            print(f"{lab:16s}   (no spectra for this model/window)")
            continue
        w, r0 = spec(f"{D}/spec_{a.model}_{win}_lte.txt")
        _, r1 = spec(f"{D}/spec_{a.model}_{win}_nlte.txt")
        _, rp = spec(f"{D}/cog/spec_{a.model}_{win}_{el}p1.txt")
        e0, e1, ep = (ew(w, x, wc, half) for x in (r0, r1, rp))
        slope = (ep - e0)/0.1
        # A window whose EW barely responds to the element cannot yield an
        # abundance correction: solar Fe I 3931.4 sits in the Ca II K wing,
        # which dominates its EW and does not care about iron, so its slope
        # comes out NEGATIVE and the division explodes to -4 dex.  Refuse
        # rather than print a number that looks like a measurement.
        if slope < 5.0:
            corr[lab] = (e0, e1, np.nan)
            print(f"{lab:16s}{e0:9.1f}{e1-e0:10.2f}{slope:11.1f}"
                  f"{'  n/a':>12s}   (window EW does not respond to "
                  f"{el}; blended)")
            continue
        dlogA = -(e1 - e0)/slope
        corr[lab] = (e0, e1, dlogA)
        print(f"{lab:16s}{e0:9.1f}{e1-e0:10.2f}{slope:11.1f}{dlogA:+12.3f}")

    if a.model != "sun":
        print("\n(2) skipped: the IAG atlas is the Sun, and there is no "
              "observed counterpart for this model\n")
        return
    print("\n(2) AGAINST THE OBSERVED SOLAR FLUX (IAG atlas)\n")
    print(f"{'line':16s}{'EW_obs':>9s}{'EW_LTE':>9s}{'EW_NLTE':>9s}"
          f"{'LTE-obs':>9s}{'NLTE-obs':>10s}   verdict")
    for lab, win, wc, half, el in LINES:
        if not os.path.exists(f"{D}/iag_{win}.txt"):
            continue
        w, r0 = spec(f"{D}/spec_{a.model}_{win}_lte.txt")
        _, r1 = spec(f"{D}/spec_{a.model}_{win}_nlte.txt")
        wo, ro = atlas(win)
        eo = ew(wo, ro, wc, half)
        e0, e1 = ew(w, r0, wc, half), ew(w, r1, wc, half)
        d0, d1 = e0 - eo, e1 - eo
        if abs(e1 - e0) < 0.02*abs(d0):
            v = "NLTE change negligible vs the LTE residual"
        elif abs(d1) < abs(d0):
            v = "TOWARD the observation"
        else:
            v = "AWAY from the observation"
        print(f"{lab:16s}{eo:9.1f}{e0:9.1f}{e1:9.1f}{d0:9.1f}{d1:10.1f}   {v}")


if __name__ == "__main__":
    main()
