#!/usr/bin/env python3
"""LTE vs NLTE zoom panels for the Ca II and Mg I diagnostics.

Reads the .spec pairs written by workdir/nlte_lines/run_ab.sh -- one run with
NLTE_MODE = 0 and one with NLTE_MODE = 1, otherwise identical -- and draws one
column per multiplet and one row per model, plus a per-line summary table.

Normalised flux (flux / continuum flux) is plotted rather than flux: the
continuum is untouched by line-level departure coefficients, so the ratio
isolates the line effect and lets three very different models share an axis.

The core depth is measured within +-0.3 A of the NOMINAL line centre, not at
the deepest pixel of the panel: in an M dwarf the deepest pixel of the Mg b
window belongs to a blend, and a panel-wide minimum silently reports that
instead.  Equivalent width is likewise integrated over a fixed +-2 A about
each centre, so the two models' numbers refer to the same interval.

  python3 tools/nlte_lines_plot.py OUT.png [--specdir workdir/nlte_lines]
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Vacuum wavelengths [A] as they appear in gfall; the runs are in vacuum.
# Vacuum wavelengths [A] as they appear in gfall; the runs are in vacuum.
# One column per LINE, not per multiplet: at multiplet width the two curves
# lie on top of each other and the whole effect hides inside the cores.
LINE_SETS = {
    # Ca II and Mg I diagnostics
    "camg": [
        ("Ca II K",     "HK",  3934.777, 2.5),
        ("Ca II H",     "HK",  3969.591, 2.5),
        (r"Mg I b$_4$", "Mgb", 5168.761, 1.2),
        (r"Mg I b$_2$", "Mgb", 5174.125, 1.2),
        (r"Mg I b$_1$", "Mgb", 5185.048, 1.2),
        ("Ca II 8498",  "CaT", 8500.355, 2.5),
        ("Ca II 8542",  "CaT", 8544.451, 2.5),
        ("Ca II 8662",  "CaT", 8664.520, 2.5),
    ],
    # Fe I.  The first three are multiplet 4 (a5D - z5D), the strongest
    # low-excitation Fe I in the H&K window; 5167.7 is from the ground state
    # and sits beside Mg b4; 5170.7 is high-excitation, included for contrast
    # because overionization should touch it far less.
    "fe": [
        ("Fe I 3931.4", "HK",  3931.409, 0.5),
        ("Fe I 3924.0", "HK",  3924.022, 0.5),
        ("Fe I 3921.4", "HK",  3921.368, 0.5),
        ("Fe I 5167.7", "Mgb", 5167.721, 0.5),
        ("Fe I 5170.3", "Mgb", 5170.337, 0.5),
        ("Fe I 5170.7", "Mgb", 5170.738, 0.5),
        ("Fe I 8517.4", "CaT", 8517.449, 0.6),
        ("Fe I 8584.6", "CaT", 8584.616, 0.6),
    ],
}
MODEL_LABELS = {
    "sun": "Sun\n5777 K / 4.44\n[Fe/H] = 0",
    "k40": "K dwarf\n4000 K / 4.5\n[Fe/H] = 0",
    "m30": "M dwarf\n3000 K / 5.0\n[Fe/H] = 0",
    "mp2": "metal-poor dwarf\n5750 K / 4.5\n[Fe/H] = $-$2",
}
DEFAULT_MODELS = ["sun", "k40", "m30"]

# Okabe-Ito: black for LTE, vermillion for NLTE -- distinguished by lightness
# as well as hue, so it survives red-green colour blindness and greyscale.
C_LTE, C_NLTE = "#000000", "#D55E00"

CORE_HALF = 0.30      # A, window for the core depth
EW_HALF   = 2.00      # A, window for the equivalent width


def read_spec(path):
    d = np.loadtxt(path)
    wl, f, fc = d[:, 0], d[:, 1], d[:, 2]
    return wl, np.where(fc > 0, f / np.maximum(fc, 1e-300), np.nan)


def metrics(w, r0, r1, wc):
    m = np.abs(w - wc) <= CORE_HALF
    c0, c1 = np.nanmin(r0[m]), np.nanmin(r1[m])
    m = np.abs(w - wc) <= EW_HALF
    dw = np.gradient(w[m])
    ew0 = np.nansum((1.0 - r0[m])*dw)*1e3      # mA
    ew1 = np.nansum((1.0 - r1[m])*dw)*1e3
    return c0, c1, ew0, ew1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--specdir", default="workdir/nlte_lines")
    ap.add_argument("--set", default="camg", choices=sorted(LINE_SETS),
                    help="which line set to draw (default camg)")
    ap.add_argument("--models", default=",".join(DEFAULT_MODELS),
                    help="comma-separated model keys, one row each")
    ap.add_argument("--tag", default="nlte",
                    help="tag of the NLTE-side runs to compare against")
    a = ap.parse_args()

    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix",
        "font.serif": ["STIXGeneral"], "font.size": 9,
        "axes.grid": False, "axes.linewidth": 0.8,
        "xtick.direction": "in", "ytick.direction": "in",
    })

    LINES = LINE_SETS[a.set]
    MODELS = [(k, MODEL_LABELS[k]) for k in a.models.split(',')]
    nr, nc = len(MODELS), len(LINES)
    fig, axes = plt.subplots(nr, nc, figsize=(1.85*nc, 2.4*nr))

    print(f"{'model':6s}{'line':13s}{'core LTE':>10s}{'core NLTE':>11s}"
          f"{'d core %':>10s}{'EW LTE mA':>11s}{'d EW %':>9s}")
    for ri, (mkey, mlab) in enumerate(MODELS):
        for ci, (llab, win, wc, half) in enumerate(LINES):
            ax = axes[ri, ci]
            p0 = os.path.join(a.specdir, f"spec_{mkey}_{win}_lte.txt")
            p1 = os.path.join(a.specdir, f"spec_{mkey}_{win}_{a.tag}.txt")
            if not (os.path.exists(p0) and os.path.exists(p1)):
                ax.text(0.5, 0.5, "missing", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            w, r0 = read_spec(p0)
            _, r1 = read_spec(p1)
            m = np.abs(w - wc) <= half
            # Offset from line centre, not absolute lambda: at +-0.5 A
            # matplotlib renders absolute wavelengths as an unreadable
            # "+3.931e3" offset label.
            ax.plot(w[m] - wc, r0[m], color=C_LTE,  lw=0.8, label="LTE")
            ax.plot(w[m] - wc, r1[m], color=C_NLTE, lw=0.8, label="NLTE")

            c0, c1, ew0, ew1 = metrics(w, r0, r1, wc)
            dcore = 100.0*(c1 - c0)/c0 if c0 > 0 else np.nan
            dew = 100.0*(ew1 - ew0)/ew0 if ew0 != 0 else np.nan
            plain = llab.replace("$", "").replace("_", "")
            print(f"{mkey:6s}{plain:13s}{c0:10.4f}{c1:11.4f}{dcore:10.2f}"
                  f"{ew0:11.1f}{dew:9.2f}")
            ax.annotate(f"{dcore:+.0f}%", xy=(0.04, 0.06),
                        xycoords="axes fraction", fontsize=7.5, color=C_NLTE)

            top = np.nanmax([np.nanmax(r0[m]), np.nanmax(r1[m])])
            ax.set_ylim(0, max(1.05*top, 1e-3))
            ax.set_xlim(-half, half)
            ax.spines[["top", "right"]].set_visible(False)
            if ri == 0:
                ax.set_title(f"{llab}\n{wc:.2f} " + r"$\rm\AA$ vac", fontsize=9)
            if ci == 0:
                ax.set_ylabel(mlab + "\n\nnormalised flux", fontsize=8)
            if ri == nr - 1:
                ax.set_xlabel(r"$\lambda-\lambda_0$ [$\rm\AA$]", fontsize=7.5)
            ax.tick_params(labelsize=6.5)
            if ri == 0 and ci == 0:
                ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(a.out, dpi=300)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
