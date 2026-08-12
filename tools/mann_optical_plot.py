#!/usr/bin/env python3
"""The house optical comparison figure for a Mann star.

One wide panel over 4500-7800 A: data, the model, and PHOENIX NewEra, with no
shaded windows.  This is the standing format for optical model-vs-data work --
the molecular bands are broad and overlapping there, so a single wide panel
reads better than sub-panels, and shading band definitions over the top hides
the very structure being judged.

A second panel gives the ratio to the data, which is where the blue excess
below ~5500 A is obvious: at 2700 K both ATLAS12+SYNTHE and PHOENIX run far
above the observed spectrum there, a shared residual that is NOT a
model-vs-model difference.

Usage:
    python3 tools/mann_optical_plot.py GJ644C out.png [tag]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml                                        # noqa: E402
from mann_compare_plot import _model_on_obs, _phoenix_on_obs   # noqa: E402

RSYN = 300000.0
WLO, WHI = 4500.0, 7800.0


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    s = ml.resolve(sys.argv[1])
    out = sys.argv[2]
    tag = sys.argv[3] if len(sys.argv) > 3 else None
    spec = os.path.join(ml.rundir_for(s, tag), s.name + ".spec")

    wobs, fobs, _ = ml.read_spectrum(s)
    data = fobs / s.dilute
    w, f, _ = ml.read_spec_file(spec)
    mod = _model_on_obs(wobs, w, f, RSYN, s.obs_split)
    phx = _phoenix_on_obs(wobs, s)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 13,
        "axes.grid": False, "figure.dpi": 300,
    })
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(15, 8), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]})

    m = (wobs >= WLO) & (wobs <= WHI)
    sc = 10.0 ** np.floor(np.log10(np.nanmax(data[m])))
    ax1.plot(wobs[m], data[m] / sc, color="0.15", lw=0.8, zorder=0,
             label="Mann data")
    ax1.plot(wobs[m], phx[m] / sc, color="#E69F00", lw=0.8, zorder=1,
             label="PHOENIX NewEra")
    ax1.plot(wobs[m], mod[m] / sc, color="#0072B2", lw=0.8, zorder=2,
             label="ATLAS12+SYNTHE")
    ax1.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                   r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
    ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc="upper left", fontsize=12)
    ax1.set_title(f"{s.name}   $T_{{\\rm eff}}$={s.teff:.0f} K, "
                  f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}", fontsize=13)

    ax2.axhline(1.0, color="0.15", lw=0.8)
    ax2.plot(wobs[m], phx[m] / data[m], color="#E69F00", lw=0.7, alpha=0.9)
    ax2.plot(wobs[m], mod[m] / data[m], color="#0072B2", lw=0.7, alpha=0.9)
    ax2.set_ylim(0.3, 2.2)
    ax2.set_ylabel("model / data")
    ax2.set_xlabel(r"wavelength  [$\mathrm{\AA}$]")
    ax2.set_xlim(WLO, WHI)
    for ax in (ax1, ax2):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
