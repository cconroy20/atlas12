#!/usr/bin/env python3
"""Detailed optical plots for the GJ887 T-min fit.

Four wavelength panels across 540-770 nm at the observed resolution:
data, RE baseline, best-fit T-min model, and two bracketing models.
Second figure: model/data ratios in the same panels, with the scored
TiO windows and CaH index bands shaded.
"""
import json
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

REPO = "/Users/cconroy/kurucz/atlas12"
import os as _os
STAR = _os.environ.get("TMIN_STAR", "GJ887")
FIT = f"{REPO}/workdir/mann/{STAR}/tmin_fit"
MANN_DIR = os.path.expanduser("~/sps/SPECTRA/Mann")
RSUN, PC, C_ANG = 6.957e10, 3.0856776e18, 2.99792458e18

d = fits.open(f"{MANN_DIR}/M_params.fits")[1].data[0]
names = np.array([n.strip() for n in d["NAME"]])
i = int(np.where(names == STAR)[0][0])
dilute = (d["RADIUS"][i] * RSUN / (d["DISTANCE"][i] * PC)) ** 2
wobs_um, fobs, _ = np.loadtxt(f"{MANN_DIR}/{STAR}.ascii", unpack=True)
wobs = wobs_um * 1e4
ok = (fobs > 0) & (wobs > 4450) & (wobs < 7850)
WOBS, FOBS = wobs[ok], fobs[ok] / dilute

def smooth_to_R(f, mR, oR):
    return gaussian_filter1d(f, mR / (oR * 2.3548), mode="nearest")

def model_on_obs(spec_path):
    w, hnu, _ = np.loadtxt(spec_path, unpack=True)
    f = 4 * np.pi * hnu * C_ANG / w ** 2
    fs = smooth_to_R(f, 50000., 1000.)
    return np.interp(WOBS, w, fs, left=np.nan, right=np.nan)

scores = json.load(open(f"{FIT}/scores.json"))
scores.update(json.load(open(f"{FIT}/scores3.json")))
if os.path.exists(f"{FIT}/scoresm.json"):
    scores.update(json.load(open(f"{FIT}/scoresm.json")))
best_tag = min(scores, key=lambda k: scores[k]["chi2"])
bs = scores[best_tag]
print(f"best grid point: {best_tag}  (tm={bs['tm']}, s={bs['s']}) "
      f"chi2={bs['chi2']:.1f}")
order = sorted(scores.values(), key=lambda v: v["chi2"])[:6]
print("top 6:", [(v['tm'], v['s'], round(v['chi2'], 1)) for v in order])

def tag_of(v):
    coord = v.get("coord", "tau")
    pre = "cm" if coord == "cmass" else "tm"
    fmt = f"{v['tm']:+.2f}" if coord == "cmass" else f"{v['tm']:+.1f}"
    t = f"{pre}{fmt}_s{v['s']:.0f}".replace("+", "p").replace("-", "m")
    if v.get("dt"):
        t += f"_dt{v['dt']:.0f}"
    return t

def label_of(v):
    coord = v.get("coord", "tau")
    dep = ("log m" if coord == "cmass" else "tau_min")
    lab = f"{dep}={v['tm']}, s={v['s']:.0f}"
    if v.get("dt"):
        lab += f", dT={v['dt']:.0f}"
    return lab

# wide (450-780 nm) syntheses where available; fall back to the fit window
def spec_path(tag):
    for p in (f"{FIT}/mgh_ab_old/{tag}.spec", f"{FIT}/wide_{tag}/{tag}.spec",
              f"{FIT}/{tag}_opt/{tag}.spec"):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(tag)

models = [("RE SYNTHE", f"{REPO}/workdir/mann/{STAR}/{STAR}.spec", "#0072B2"),
          ("best-fit chromo (" + label_of(bs) + ")", spec_path(best_tag),
           "#CC79A7")]

curves = [(lab, model_on_obs(p), col) for lab, p, col in models]

# PHOENIX NewEra 3700/5.0 (surface f_lambda, high-res) -> obs resolution
PHX = _os.environ.get("PHOENIX_NPZ", _os.path.expanduser(
    "~/kurucz/grids/NewEra_HSR/phoenix_hires_3700g5.0.npz"))
ph = np.load(PHX)
RP = 500000.0
n = int(np.log(ph["wl"][-1] / ph["wl"][0]) * RP) + 1
wlog = ph["wl"][0] * np.exp(np.arange(n) / RP)
plog = np.interp(wlog, ph["wl"], ph["flam"])
p1000 = smooth_to_R(plog, RP, 1000.0)
curves.append(("PHOENIX NewEra 3700/5.0",
               np.interp(WOBS, wlog, p1000), "#E69F00"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "stix", "font.size": 10,
    "axes.spines.top": True, "axes.spines.right": True,
    "axes.grid": False, "figure.dpi": 300,
})

LO, HI = 4500, 7800
TIO_W = [(5500, 5700), (5900, 6100), (7053, 7270), (7450, 7650)]
CAH_W = [(6814, 6846), (6960, 6990)]
CONT_W = [(7042, 7046)]

for figname, ratio in [("tmin_fit_spectra.png", False),
                       ("tmin_fit_ratios.png", True)]:
    fig, ax = plt.subplots(figsize=(18, 5.5))
    m = (WOBS > LO) & (WOBS < HI)
    if ratio:
        ax.axhline(1, color="0.5", lw=0.6)
        for lab, f, col in curves:
            ax.plot(WOBS[m], f[m] / FOBS[m], color=col, lw=1.1, label=lab)
        ax.set_ylim(0.55, 1.45)
        ax.set_ylabel("model / data")
    else:
        ax.plot(WOBS[m], FOBS[m] / 1e5, color="0.1", lw=1.6, label="Mann+15")
        for lab, f, col in curves:
            ax.plot(WOBS[m], f[m] / 1e5, color=col, lw=1.1, alpha=0.95,
                    label=lab)
        ax.set_ylabel(r"$F_\lambda^{\rm surf}$ [$10^5$ erg s$^{-1}$ cm$^{-2}$"
                      r" $\mathrm{\AA}^{-1}$]")
    ax.set_xlim(LO, HI)
    ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
    ax.set_xlabel(r"wavelength [$\mathrm{\AA}$]")
    ax.set_title("GJ887 T-min fit, optical")
    fig.tight_layout()
    fig.savefig(f"{FIT}/{figname}", bbox_inches="tight")
    print("wrote", f"{FIT}/{figname}")
