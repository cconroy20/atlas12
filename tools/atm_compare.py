#!/usr/bin/env python3
"""Overlay the T(depth) structure of several Kurucz-format .atm files.

Usage:  python3 tools/atm_compare.py out.png label=file.atm [label=file.atm ...]

Plots T against log column mass and against log tau_5000, and the difference
of every curve from the first one AT FIXED COLUMN MASS.  Used to see whether
an opacity change fed back into the structure or only into the emergent
spectrum.

COMPARE ON COLUMN MASS, NOT ON TAU.  Column mass is the independent variable
of the hydrostatic solution and means the same thing in every code; the tau
column does not.  ATLAS12 writes tau_5000, PHOENIX restart files carry tau at
their own standard wavelength (12000 A for NewEra) on a prescribed grid, and
at M-dwarf temperatures those two scales diverge by ~1 dex in column mass by
log tau = -4 (5000 A sits inside heavy molecular absorption, 12000 A does
not).  Reading T at equal tau across the two codes therefore compares layers a
decade apart and manufactures differences of several hundred K.  The tau panel
is kept only for within-code comparisons.
"""
import os
import sys

import numpy as np

COLORS = ["#0072B2", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "0.3"]


def read_atm(path):
    """Return (rhox, T, tau5000) from the DECK6 block."""
    rows = []
    with open(path) as fh:
        started = False
        for line in fh:
            if line.startswith("READ DECK"):
                started = True
                continue
            if started:
                p = line.split()
                if len(p) < 11:
                    break
                rows.append([float(p[0]), float(p[1]), float(p[10])])
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    out = sys.argv[1]
    models = []
    for a in sys.argv[2:]:
        label, path = a.split("=", 1)
        models.append((label, *read_atm(path)))
        print(f"{label:12s} {len(models[-1][1]):3d} layers  "
              f"T {models[-1][2].min():.0f}-{models[-1][2].max():.0f} K  "
              f"{path}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 14,
        "axes.grid": False, "figure.dpi": 300,
    })
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ref = models[0]
    for k, (lab, rhox, T, tau) in enumerate(models):
        c = COLORS[k % len(COLORS)]
        axes[0].plot(np.log10(rhox), T, color=c, lw=1.6, label=lab)
        axes[1].plot(np.log10(tau), T, color=c, lw=1.6)
        if k:
            # Differences on column mass only -- see the module docstring on
            # why equal-tau differences are meaningless across codes.
            # Restrict to the column-mass range both models actually cover;
            # np.interp clamps outside it and would draw a spurious spike.
            lm, lr = np.log10(rhox), np.log10(ref[1])
            ov = (lm >= lr.min()) & (lm <= lr.max())
            axes[2].plot(lm[ov], T[ov] - np.interp(lm[ov], lr, ref[2]),
                         color=c, lw=1.6, label=f"{lab} - {ref[0]}")
    axes[2].axhline(0.0, color="0.2", lw=0.8)
    axes[0].set_xlabel(r"$\log\,m$  [g cm$^{-2}$]")
    axes[1].set_xlabel(r"$\log\,\tau_{5000}$ (within-code only)")
    axes[2].set_xlabel(r"$\log\,m$  [g cm$^{-2}$]")
    axes[0].set_ylabel("T  [K]")
    axes[2].set_ylabel(r"$\Delta T$  [K]")
    axes[0].legend(frameon=False, fontsize=12)
    if len(models) > 1:
        axes[2].legend(frameon=False, fontsize=12)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
