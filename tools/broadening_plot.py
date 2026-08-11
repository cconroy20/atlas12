#!/usr/bin/env python3
"""Figure for a Mann star: the data, PHOENIX NewEra and one model.

Kept deliberately spare -- one model curve at a time.  Extra labels are
accepted for A/B work but the default use is a single latest model.

Shows the near-IR, where over-broadened H2O lines were deleting the weak-line
haze between the water bands, against the Mann spectrum and PHOENIX NewEra.

Usage:
    python3 tools/broadening_plot.py GJ644C out.png label=tag [label=tag ...]

Each label=tag names a run directory workdir/mann/<STAR>_<tag> ('.' = the
untagged run).  Colours are Okabe-Ito with no red/green pairing.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml                                       # noqa: E402
from mann_compare_plot import _model_on_obs, _phoenix_on_obs  # noqa: E402

RSYN = 300000.0
COLORS = ["#0072B2", "#CC79A7", "#56B4E9", "#F0E442"]
BANDS = [("J", 11500, 13500), ("H", 15500, 17000), ("K", 21000, 23000)]


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 1
    s = ml.resolve(sys.argv[1])
    out = sys.argv[2]
    models = []
    for a in sys.argv[3:]:
        label, tag = a.split("=", 1)
        tag = None if tag in (".", "none", "") else tag
        w, f, _ = ml.read_spec_file(
            os.path.join(ml.rundir_for(s, tag), s.name + ".spec"))
        models.append((label, w, f))

    wobs, fobs, _ = ml.read_spectrum(s)
    data = fobs / s.dilute
    curves = [(lab, _model_on_obs(wobs, w, f, RSYN, s.obs_split))
              for lab, w, f in models]
    phx = _phoenix_on_obs(wobs, s)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 13,
        "axes.grid": False, "figure.dpi": 300,
    })
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9))
    wum = wobs / 1.0e4

    # --- top: near-IR flux, where the effect lives
    m = (wum >= 0.95) & (wum <= 2.45)
    sc = 10.0 ** np.floor(np.log10(np.nanmax(data[m])))
    ax1.plot(wum[m], data[m] / sc, color="0.15", lw=0.7, zorder=0,
             label="Mann data")
    ax1.plot(wum[m], phx[m] / sc, color="#E69F00", lw=0.7, zorder=1,
             label="PHOENIX NewEra")
    for k, (lab, fc) in enumerate(curves):
        ax1.plot(wum[m], fc[m] / sc, color=COLORS[k % len(COLORS)], lw=0.7,
                 zorder=5 - k, label=lab)
    ax1.set_xlim(0.95, 2.45)
    ax1.set_ylim(bottom=0)
    # Band markers last, in axes coordinates, so they cannot be thrown off
    # by a later set_ylim (and sit low, clear of the legend).
    for lab, lo, hi in BANDS:
        ax1.annotate(lab, xy=((lo + hi) / 2e4, 0.04),
                     xycoords=("data", "axes fraction"),
                     ha="center", fontsize=11, color="0.45")
    ax1.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                   r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
    ax1.set_xlabel(r"wavelength  [$\mu$m]")
    ax1.legend(frameon=False, loc="lower left", fontsize=11)
    ax1.set_title(f"{s.name}   $T_{{\\rm eff}}$={s.teff:.0f} K, "
                  f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}", fontsize=13)

    # --- bottom: ratio to data over the full range, to show the optical
    #     is untouched while the near-IR moves
    ax2.axhline(1.0, color="0.15", lw=0.8)
    ax2.plot(wum, phx / data, color="#E69F00", lw=0.5, alpha=0.8)
    for k, (lab, fc) in enumerate(curves):
        ax2.plot(wum, fc / data, color=COLORS[k % len(COLORS)], lw=0.5,
                 alpha=0.9)
    ax2.set_xscale("log")
    ax2.set_xlim(0.4, 2.45)
    ax2.set_xticks([0.4, 0.6, 0.8, 1.0, 1.5, 2.0])
    ax2.set_xticklabels(["0.4", "0.6", "0.8", "1.0", "1.5", "2.0"])
    ax2.set_ylim(0.4, 1.7)
    ax2.set_ylabel("model / data")
    ax2.set_xlabel(r"wavelength  [$\mu$m]")
    for ax in (ax1, ax2):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
