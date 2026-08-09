#!/usr/bin/env python3
"""Three-panel Mann-comparison figure (house format, 2026-08-09):
top 0.35-2.5 um (log-lambda); middle 4500-7500 A; bottom 8000-10000 A
(autoscaled).  Data plotted first; models overlaid; no shaded windows;
no title parentheticals.  Smoothing = measured Mann LSF (SNIFS constant
11 A FWHM below 9500 A, SpeX R=2000 above; validate_mann.smooth_mann).

The model spectrum must be a FULL-RANGE (wlbeg=350 wlend=2500) synthesis
at R >= 300,000 (see CHANGELOG 2026-08-09 on resolution convergence).

Usage:  python3 mann_compare_plot.py [STAR]        (default GJ887)
Env:    PHOENIX_NPZ  overlay spectrum {wl[A], flam} (skipped if absent)
        PHX_LABEL    its legend label
"""
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAR = sys.argv[1] if len(sys.argv) > 1 else "GJ887"
MANN_DIR = os.path.expanduser("~/sps/SPECTRA/Mann")
RSUN, PC, C_ANG = 6.957e10, 3.0856776e18, 2.99792458e18
RSYN = 300000.0
SNIFS_DLAM = 11.0            # A FWHM (measured; see reference plot)
SPLIT = 9500.0               # 2016-library SNIFS->SpeX splice [A]

d = fits.open(f"{MANN_DIR}/M_params.fits")[1].data[0]
names = np.array([n.strip() for n in d["NAME"]])
i = int(np.where(names == STAR)[0][0])
teff, feh = d["TEFF"][i], d["FEH"][i]
logg = np.log10(6.67430e-8 * d["MASS"][i] * 1.98892e33
                / (d["RADIUS"][i] * RSUN) ** 2)
dilute = (d["RADIUS"][i] * RSUN / (d["DISTANCE"][i] * PC)) ** 2
wobs_um, fobs, _ = np.loadtxt(f"{MANN_DIR}/{STAR}.ascii", unpack=True)
wobs = wobs_um * 1e4
ok = fobs > 0
WOBS, FOBS = wobs[ok], fobs[ok] / dilute

def smooth_to_R(f, mR, oR):
    return gaussian_filter1d(f, mR / (oR * 2.3548), mode="nearest")

def smooth_mann(w, f, mR):
    step = w[0] / mR
    wl = np.arange(w[0], min(w[-1], SPLIT + 200.0), step)
    fl = np.interp(wl, w, f)
    fsn = gaussian_filter1d(fl, SNIFS_DLAM / 2.3548 / step, mode="nearest")
    out = smooth_to_R(f, mR, 2000.)
    m = w < SPLIT
    out[m] = np.interp(w[m], wl, fsn)
    return out

def model_on_obs(spec_path, mR=RSYN):
    w, hnu, _ = np.loadtxt(spec_path, unpack=True)
    f = 4 * np.pi * hnu * C_ANG / w ** 2
    fs = smooth_mann(w, f, mR)
    return np.interp(WOBS, w, fs, left=np.nan, right=np.nan)

curves = [("Mann+15", FOBS, "0.2", 0, 0.9)]
curves.append(("ATLAS12 RE",
               model_on_obs(f"{REPO}/workdir/mann/{STAR}/{STAR}.spec"),
               "#0072B2", 5, 0.9))
phx_npz = os.environ.get("PHOENIX_NPZ")
if phx_npz and os.path.exists(os.path.expanduser(phx_npz)):
    ph = np.load(os.path.expanduser(phx_npz))
    RP = 500000.0
    n = int(np.log(ph["wl"][-1] / ph["wl"][0]) * RP) + 1
    wlog = ph["wl"][0] * np.exp(np.arange(n) / RP)
    plog = np.interp(wlog, ph["wl"], ph["flam"])
    psm = smooth_mann(wlog, plog, RP)
    curves.append((os.environ.get("PHX_LABEL", "PHOENIX NewEra"),
                   np.interp(WOBS, wlog, psm, left=np.nan, right=np.nan),
                   "#E69F00", 1, 0.8))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "stix", "font.size": 14,
    "axes.grid": False, "figure.dpi": 300,
})

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(13, 11))
sc = 10.0 ** np.floor(np.log10(np.nanmax(FOBS)))
for lo, hi, ax, lw in [(3500, 25000, ax1, 0.6), (4500, 7500, ax2, 0.8),
                       (8000, 10000, ax3, 0.8)]:
    m = (WOBS > lo) & (WOBS < hi)
    for lab, f, c, z, al in curves:
        ax.plot(WOBS[m], f[m] / sc, color=c, lw=lw, alpha=al, zorder=z,
                label=lab if ax is ax1 else None)
    ax.set_xlim(lo, hi)
    if ax is ax3:
        vals = [f[m] for _, f, *_ in curves]
        vlo = np.nanmin([np.nanmin(v) for v in vals])
        vhi = np.nanmax([np.nanmax(v) for v in vals])
        pad = 0.06 * (vhi - vlo)
        ax.set_ylim((vlo - pad) / sc, (vhi + pad) / sc)
    else:
        ax.set_ylim(bottom=0)
    ax.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                  r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
ax1.set_xscale("log")
ax1.set_xticks([4000, 6000, 8000, 10000, 15000, 20000, 25000])
ax1.set_xticklabels(["0.4", "0.6", "0.8", "1.0", "1.5", "2.0", "2.5"])
ax1.set_xlabel(r"wavelength  [$\mu$m]")
ax1.legend(frameon=False, loc="upper right", fontsize=12)
ax1.set_title(f"{STAR}   $T_{{\\rm eff}}$={teff:.0f} K, "
              f"log g={logg:.2f}, [Fe/H]={feh:+.2f}", fontsize=14)
ax2.set_xlabel(r"wavelength  [$\mathrm{\AA}$]")
ax3.set_xlabel(r"wavelength  [$\mathrm{\AA}$]")
fig.tight_layout()
out = f"{REPO}/workdir/mann/{STAR}/{STAR}_full_compare.png"
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
