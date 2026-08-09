#!/usr/bin/env python3
"""Plots for the tmin_rf verified inversion: wide-panel spectra comparison
(data, RE baseline, verified signed-tilt model, PHOENIX NewEra) and the
T(log tau) structures.  Star-generic via TMIN_STAR / PHOENIX_NPZ / PHX_LABEL.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tmin_rf import (REPO, STAR, RF, BASE, WOBS, FOBS, smooth_to_R)

TILT = os.environ.get("TMIN_RF_TAG", "vsigned30")
PHX = os.environ.get("PHOENIX_NPZ", os.path.expanduser(
    "~/kurucz/grids/NewEra_HSR/phoenix_hires_3700g5.0.npz"))
PHX_LABEL = os.environ.get("PHX_LABEL", "PHOENIX NewEra 3700/5.0")


def model_on_obs_wide(tag):
    p = f"{RF}/{tag}_wide/{tag}.spec"
    if not os.path.exists(p):
        p = f"{RF}/{tag}_opt/{tag}.spec"
    w, hnu, _ = np.loadtxt(p, unpack=True)
    f = 4 * np.pi * hnu * 2.99792458e18 / w ** 2
    fs = np.where(w < 9500., smooth_to_R(f, 50000., 1000.),
                  smooth_to_R(f, 50000., 2000.))
    return np.interp(WOBS, w, fs, left=np.nan, right=np.nan)


def deck_T(atm):
    L = open(atm).readlines()
    i = next(k for k, l in enumerate(L) if l.startswith('READ DECK')) + 1
    nd = int(L[i - 1].split()[2])
    rows = [(float(l.split()[1]), float(l.split()[-1])) for l in L[i:i + nd]]
    T = np.array([r[0] for r in rows])
    lt = np.log10([r[1] for r in rows])
    return lt, T


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "stix", "font.size": 10,
    "axes.grid": False, "figure.dpi": 300,
})

RF_LABEL = os.environ.get("TMIN_RF_LABEL", f"RF inversion ({TILT})")
curves = [("RE SYNTHE", model_on_obs_wide("base"), "#0072B2"),
          (RF_LABEL, model_on_obs_wide(TILT), "#CC79A7")]
# optional extra model curve, e.g. a Teff+50 K RE model: "tag:label"
EXTRA = os.environ.get("TMIN_RF_EXTRA")
if EXTRA:
    etag, elabel = EXTRA.split(":", 1)
    curves.append((elabel, model_on_obs_wide(etag), "#56B4E9"))

ph = np.load(PHX)
RP = 500000.0
n = int(np.log(ph["wl"][-1] / ph["wl"][0]) * RP) + 1
wlog = ph["wl"][0] * np.exp(np.arange(n) / RP)
plog = np.interp(wlog, ph["wl"], ph["flam"])
psm = np.where(wlog < 9500., smooth_to_R(plog, RP, 1000.),
               smooth_to_R(plog, RP, 2000.))
curves.append((PHX_LABEL, np.interp(WOBS, wlog, psm), "#E69F00"))

LO, HI = 4500, 12000
m = (WOBS > LO) & (WOBS < HI)
fig, ax = plt.subplots(figsize=(18, 5.5))
ax.plot(WOBS[m], FOBS[m] / 1e5, color="0.1", lw=1.6, label="Mann+15")
for lab, f, col in curves:
    ax.plot(WOBS[m], f[m] / 1e5, color=col, lw=1.1, alpha=0.95, label=lab)
ax.set_xlim(LO, HI)
ax.set_ylabel(r"$F_\lambda^{\rm surf}$ [$10^5$ erg s$^{-1}$ cm$^{-2}$"
              r" $\mathrm{\AA}^{-1}$]")
ax.set_xlabel(r"wavelength [$\mathrm{\AA}$]")
ax.set_title(f"{STAR}: {RF_LABEL}")
ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
fig.tight_layout()
fig.savefig(f"{RF}/rf_bestfit_spectra.png", bbox_inches="tight")
print("wrote", f"{RF}/rf_bestfit_spectra.png")

fig, ax = plt.subplots(figsize=(7, 5))
struct = [("RE", f"{RF}/base.atm", "#0072B2", 1.6),
          (RF_LABEL, f"{RF}/{TILT}.atm", "#CC79A7", 1.6)]
if EXTRA:
    struct.append((elabel, f"{RF}/{etag}.atm", "#56B4E9", 1.2))
for lab, atm, col, lw in struct:
    if os.path.exists(atm):
        lt, T = deck_T(atm)
        ax.plot(lt, T, color=col, lw=lw, label=lab)
ax.set_xlim(2, -4.5)
ax.set_xlabel(r"$\log\,\tau_{5000}$")
ax.set_ylabel("T [K]")
ax.set_title(f"{STAR}: temperature structures")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(f"{RF}/rf_structure.png", bbox_inches="tight")
print("wrote", f"{RF}/rf_structure.png")
