#!/usr/bin/env python3
"""Three-panel Mann point-comparison figure (house format, 2026-08-09):
top 0.35-2.5 um (log-lambda); middle 4000-7500 A; bottom 8000-10000 A
(autoscaled).  Data plotted first; models overlaid; no shaded windows;
no title parentheticals.  Surface-flux scale.  Smoothing = measured Mann
LSF (SNIFS constant 11 A FWHM below the library splice, SpeX R=2000
above; mann_lib.smooth_mann with the star's library-correct splice).

Curves: Mann data (black), ATLAS12+SYNTHE RE (blue #0072B2), PHOENIX
NewEra interpolated at the star's params (orange #E69F00) -- fetched
automatically from the LowRes grid via mann_lib.phoenix_newera.
Optional extra model curves in #CC79A7 / #56B4E9 (Okabe-Ito).

The model spectrum must be a FULL-RANGE (wlbeg=350 wlend=2500) synthesis
at R >= 300,000 (see CHANGELOG 2026-08-09 on resolution convergence).

Usage:  python3 mann_compare_plot.py STAR [TAG]
        (rundir workdir/mann/<STAR>[_TAG], model <STAR>.spec inside)
As API: plot_star("GJ887") or plot_star(star, extra=[(label, specpath)]).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml

RSYN = 300000.0
PHX_RP = 500000.0            # log-grid resampling R for the PHOENIX curve
EXTRA_COLORS = ["#CC79A7", "#56B4E9", "#009E73"]


def _model_on_obs(wobs, w, f, model_R, split):
    fs = ml.smooth_mann(w, f, model_R, split=split)
    return np.interp(wobs, w, fs, left=np.nan, right=np.nan)


def _phoenix_on_obs(wobs, s):
    """PHOENIX NewEra LowRes at the star's params, smoothed to the Mann LSF.

    The LowRes files are on a LINEAR 0.1 A grid; resample onto a log grid
    at R=PHX_RP (far finer than the LSF) before the constant-R smoothing.
    """
    wl, flam = ml.phoenix_newera(s.teff, s.logg, s.feh)
    n = int(np.log(wl[-1] / wl[0]) * PHX_RP) + 1
    wlog = wl[0] * np.exp(np.arange(n) / PHX_RP)
    plog = np.interp(wlog, wl, flam)
    return _model_on_obs(wobs, wlog, plog, PHX_RP, s.obs_split)


def plot_star(star, tag=None, spec=None, model_R=RSYN, extra=(),
              phoenix=True, out=None):
    """Write the 3-panel figure for one star; returns the png path.

    star  : name or resolved mann_lib.Star
    spec  : model .spec path (default workdir/mann/<name>[_tag]/<name>.spec)
    extra : iterable of (label, spec_path) additional model curves
    """
    s = ml.resolve(star) if isinstance(star, str) else star
    rundir = ml.rundir_for(s, tag)
    spec = spec or os.path.join(rundir, s.name + ".spec")

    wobs, fobs, _ = ml.read_spectrum(s)
    fsurf = fobs / s.dilute                       # data -> surface flux

    curves = [("Mann data", fsurf, "0.2", 0, 0.9)]
    w, f, _ = ml.read_spec_file(spec)             # surface flux
    curves.append(("ATLAS12+SYNTHE",
                   _model_on_obs(wobs, w, f, model_R, s.obs_split),
                   "#0072B2", 5, 0.9))
    if phoenix:
        curves.append(("PHOENIX NewEra", _phoenix_on_obs(wobs, s),
                       "#E69F00", 1, 0.8))
    for k, (lab, path) in enumerate(extra):
        we, fe, _ = ml.read_spec_file(path)
        curves.append((lab, _model_on_obs(wobs, we, fe, model_R, s.obs_split),
                       EXTRA_COLORS[k % len(EXTRA_COLORS)], 4, 0.9))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 14,
        "axes.grid": False, "figure.dpi": 600,
    })

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(13, 11))
    sc = 10.0 ** np.floor(np.log10(np.nanmax(fsurf)))
    wum = wobs / 1.0e4                              # all x-axes in um
    for lo, hi, ax, lw in [(0.35, 2.5, ax1, 0.6), (0.40, 0.75, ax2, 0.8),
                           (0.80, 1.00, ax3, 0.8)]:
        m = (wum > lo) & (wum < hi)
        for lab, fc, c, z, al in curves:
            ax.plot(wum[m], fc[m] / sc, color=c, lw=lw, alpha=al, zorder=z,
                    label=lab if ax is ax1 else None)
        ax.set_xlim(lo, hi)
        if ax is ax3:
            vals = [fc[m] for _, fc, *_ in curves]
            vlo = np.nanmin([np.nanmin(v) for v in vals])
            vhi = np.nanmax([np.nanmax(v) for v in vals])
            pad = 0.06 * (vhi - vlo)
            ax.set_ylim((vlo - pad) / sc, (vhi + pad) / sc)
        else:
            ax.set_ylim(bottom=0)
        ax.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                      r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
        ax.set_xlabel(r"wavelength  [$\mu$m]")
    ax1.set_xscale("log")
    ax1.set_xticks([0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5])
    ax1.set_xticklabels(["0.4", "0.6", "0.8", "1.0", "1.5", "2.0", "2.5"])
    ax1.legend(frameon=False, loc="upper right", fontsize=12)
    ax1.set_title(f"{s.name}   $T_{{\\rm eff}}$={s.teff:.0f} K, "
                  f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}", fontsize=14)
    fig.tight_layout()
    out = out or os.path.join(rundir, s.name + "_full_compare.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)
    return out


if __name__ == "__main__":
    star = sys.argv[1] if len(sys.argv) > 1 else "GJ887"
    tag = sys.argv[2] if len(sys.argv) > 2 else None
    plot_star(star, tag=tag)
