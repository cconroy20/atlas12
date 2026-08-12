#!/usr/bin/env python3
"""Score Ca I 4227 profile variants against the Mann spectrum.

Prints model/data in windows across the blue depression for every tagged
run given on the command line, plus a chi2-like score on the 4000-4500 and
4500-5000 A bands.  Optionally writes a comparison figure.

Usage:  python3 ca4227_compare.py [--plot OUT.png] tag [tag ...]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml
import mann_compare_plot as mcp

WINDOWS = [(3900, 4000), (4000, 4100), (4100, 4200), (4200, 4300),
           (4300, 4400), (4400, 4500), (4500, 4700), (4700, 5000),
           (5000, 5150)]
SCORE = [(4000, 4200), (4200, 4400), (4400, 4600), (4600, 4800),
         (4800, 5000)]


def load(star, tag):
    """Model smoothed to the Mann LSF, kept on the NATIVE R=300k sampling.

    Interpolating onto the observed wavelengths would downsample the model
    to SNIFS pixels (~2.5 A) on top of the LSF convolution, which is a
    second, unwanted degradation.
    """
    rundir = ml.rundir_for(star, tag)
    w, f, _ = ml.read_spec_file(os.path.join(rundir, star.name + ".spec"))
    return w, ml.smooth_mann(w, f, mcp.RSYN, split=star.obs_split)


def win_mean(w, f, lo, hi):
    m = (w > lo) & (w < hi)
    return np.mean(f[m]) if m.sum() else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("tags", nargs="+",
                    help="run-dir tag, or tag=Plot label")
    ap.add_argument("--star", default="GJ644C")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--no-phoenix", action="store_true",
                    help="omit the PHOENIX NewEra curve from the figure")
    a = ap.parse_args()

    s = ml.resolve(a.star)
    wobs, fobs, eobs = ml.read_spectrum(s)
    fsurf = fobs / s.dilute

    labels = {}
    tags = []
    for spec in a.tags:
        tag, _, lab = spec.partition("=")
        tags.append(tag)
        labels[tag] = lab or tag
    a.tags = tags

    models = {t: load(s, t) for t in a.tags}
    wphx, fphx = mcp.ml.phoenix_newera(s.teff, s.logg, s.feh)
    n = int(np.log(wphx[-1] / wphx[0]) * mcp.PHX_RP) + 1
    wphx_log = wphx[0] * np.exp(np.arange(n) / mcp.PHX_RP)
    phx = (wphx_log, ml.smooth_mann(wphx_log, np.interp(wphx_log, wphx, fphx),
                                    mcp.PHX_RP, split=s.obs_split))

    with np.errstate(divide="ignore", invalid="ignore"):
        snr = fobs / eobs
    ok = np.isfinite(snr) & (snr > 3)

    # Window means are taken on each curve's OWN grid; only the ratio is
    # formed between them.  No resampling of the model anywhere.
    hdr = f"{'window':>12}" + "".join(f"{t:>12}" for t in a.tags) + f"{'PHOENIX':>10}"
    print(hdr)
    print("-" * len(hdr))
    for lo, hi in WINDOWS:
        m = ok & (wobs > lo) & (wobs < hi)
        if m.sum() < 3:
            continue
        d = np.mean(fsurf[m])
        row = f"{lo:5d}-{hi:5d}" + "".join(
            f"{win_mean(*models[t], lo, hi)/d:12.2f}" for t in a.tags)
        print(row + f"{win_mean(*phx, lo, hi)/d:10.2f}")

    print()
    print(f"{'score':>12}" + "".join(f"{t:>12}" for t in a.tags) + f"{'PHOENIX':>10}")
    scores = {}
    for t in list(a.tags) + ["PHOENIX"]:
        wm, fm = phx if t == "PHOENIX" else models[t]
        v = []
        for lo, hi in SCORE:
            m = ok & (wobs > lo) & (wobs < hi)
            if m.sum() >= 3:
                v.append(np.log(win_mean(wm, fm, lo, hi)
                                / np.mean(fsurf[m])) ** 2)
        scores[t] = np.sqrt(np.mean(v))
    print(f"{'rms ln ratio':>12}" + "".join(f"{scores[t]:12.3f}" for t in a.tags)
          + f"{scores['PHOENIX']:10.3f}")

    if not a.plot:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                         "font.size": 14, "axes.grid": False, "figure.dpi": 600})
    colors = ["#0072B2", "#CC79A7", "#56B4E9", "#E69F00", "#F0E442"]
    curves = [("Mann data", wobs, fsurf, "0.15", 0)]
    if not a.no_phoenix:
        curves.append(("PHOENIX NewEra", phx[0], phx[1], "#E69F00", 1))
    curves += [(labels[t], models[t][0], models[t][1],
                colors[k % len(colors)], 2 + k) for k, t in enumerate(a.tags)]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 10), sharex=True)
    for ax in (ax1, ax2):
        for lab, wc, fc, col, z in curves:
            m = (wc > 3850) & (wc < 5200)
            ax.plot(wc[m] / 1e4, fc[m], color=col, lw=1.2,
                    label=lab if ax is ax1 else None, zorder=z)
        ax.set_xlim(0.385, 0.52)
        ax.set_ylabel(r"$f_\lambda^{\rm surf}$"
                      r"  [erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
        ax.minorticks_on()
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.tick_params(which="major", length=7)
        ax.tick_params(which="minor", length=3.5)
    ax1.set_ylim(0, 9000)
    ax2.set_yscale("log")
    ax2.set_ylim(30, 2e4)
    ax2.set_xlabel(r"wavelength  [$\mu$m]")
    ax1.legend(frameon=False, fontsize=12, loc="upper left")
    ax1.set_title(f"{s.name}   $T_{{\\rm eff}}$={s.teff:.0f} K, "
                  f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}", fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.04)
    fig.savefig(a.plot, bbox_inches="tight")
    print("wrote", a.plot)


if __name__ == "__main__":
    main()
