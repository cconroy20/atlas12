#!/usr/bin/env python3
"""High-resolution model-vs-UVES comparison (line-profile showcase).

Compares the full-range R=300k SYNTHE spectrum of a ladder star against
ESO UVES Phase-3 spectra (~/sps/SPECTRA/UVES/<STAR>_UVES_{uv,blue,red}.fits,
R ~ 41k/46k, AIR wavelengths, stellar RV NOT removed, flux-calibrated at
Earth in 1e-16 erg/s/cm^2/A).

Alignment: UVES air -> vacuum (mann_lib.air_to_vac), then one global
velocity fitted on the Ca II 8542 region shifts UVES to the model rest
frame (this absorbs the stellar RV + barycentric details).  Fluxes are
ABSOLUTE surface fluxes: UVES / dilution with the catalog (R/d)^2 --
UVES slit losses make the absolute scale good to only ~10-20%, so judge
shapes/depths, not the zero point.  Model smoothed to the per-arm
SPEC_RES from the UVES headers.  UVES is NOT telluric-corrected; O2/H2O
bands are marked in the panel titles where relevant.

Usage: python3 uves_compare.py [STAR]        (default GJ411)
"""
import glob
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml

STAR = sys.argv[1] if len(sys.argv) > 1 else "GJ411"
UVES_DIR = os.path.expanduser("~/sps/SPECTRA/UVES")
RSYN = 300000.0
C_KMS = 2.99792458e5

# showcase windows [A, vacuum]: (label, lo, hi).  NOTE UVES GJ411 red arm:
# flux calibration valid only from 6708 A (air), chip gap 8530-8667 A (air)
# -- Ca II 8542/8662 both fall IN the gap; 8498 is the usable triplet line.
PANELS = [
    ("CH G band", 4295, 4335),
    (r"H$\beta$", 4845, 4885),
    ("CaH A-X + TiO", 6712, 6762),
    (r"TiO $\gamma$(0,0) heads", 7050, 7140),
    ("K I 7699", 7690, 7720),
    ("Na I 8183/8195", 8175, 8215),
    (r"TiO $\epsilon$ head", 8430, 8470),
    ("Ca II 8498", 8485, 8520),
]

s = ml.resolve(STAR)
spec = os.path.join(ml.rundir_for(s), s.name + ".spec")
wm, fm, _ = ml.read_spec_file(spec)                    # vacuum, surface

# ---- UVES arms: air->vac, flux -> surface
arms = []
for path in sorted(glob.glob(f"{UVES_DIR}/{STAR}_UVES_*.fits")):
    hdul = fits.open(path)
    d = hdul[1].data[0]
    res = float(hdul[0].header.get("SPEC_RES", 45000.0))
    wl = ml.air_to_vac(d["WAVE"])
    fl = d["FLUX"].astype(float) * 1e-16 / s.dilute
    ok = np.isfinite(fl) & (fl > 0)
    arms.append((wl[ok], fl[ok], res))
if not arms:
    sys.exit(f"no UVES spectra for {STAR} in {UVES_DIR}")


def model_at(res):
    return gaussian_filter1d(fm, RSYN / (res * 2.3548), mode="nearest")


def uves_window(lo, hi, pad=5.0):
    for wl, fl, res in arms:
        m = (wl > lo - pad) & (wl < hi + pad)
        if m.sum() > 50:
            return wl[m], fl[m], res
    return None


# ---- global velocity: fit on the Na I 8183/8195 region (gap-free)
win = uves_window(8150, 8230, pad=0)
vgrid = np.arange(-150.0, 151.0, 0.5)
best = (0.0, np.inf)
wl_u, fl_u, res_u = win
fsm = model_at(res_u)
for v in vgrid:
    fmi = np.interp(wl_u, wm * (1.0 + v / C_KMS), fsm)
    sc = np.sum(fmi * fl_u) / np.sum(fmi * fmi)
    chi = float(np.sum((fl_u - sc * fmi) ** 2))
    if chi < best[1]:
        best = (v, chi)
VSHIFT = best[0]
print(f"{STAR}: UVES-vs-model velocity = {VSHIFT:+.1f} km/s "
      f"(stellar RV + conventions; UVES shifted to model frame)")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "stix", "font.size": 12,
    "axes.grid": False, "figure.dpi": 600,
})

fig, axes = plt.subplots(4, 2, figsize=(14, 14))
for (lab, lo, hi), ax in zip(PANELS, axes.ravel()):
    win = uves_window(lo, hi)
    if win is None:
        ax.set_visible(False)
        continue
    wl_u, fl_u, res_u = win
    wl_shift = wl_u / (1.0 + VSHIFT / C_KMS)           # -> model frame
    fmi = np.interp(wl_shift, wm, model_at(res_u))
    m = (wl_shift > lo) & (wl_shift < hi)
    sc = 10.0 ** np.floor(np.log10(np.nanmax(fl_u[m])))
    ax.plot(wl_shift[m] / 1e4, fl_u[m] / sc, color="0.2", lw=0.9,
            label="UVES")
    ax.plot(wl_shift[m] / 1e4, fmi[m] / sc, color="#0072B2", lw=0.9,
            alpha=0.9, label="ATLAS12+SYNTHE")
    r = np.median(fmi[m] / fl_u[m])
    ax.set_title(f"{lab}   (model/UVES = {r:.2f})", fontsize=12)
    ax.set_xlim(lo / 1e4, hi / 1e4)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"wavelength  [$\mu$m]")
    ax.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                  r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
    print(f"  {lab:<24} {lo}-{hi} A  model/UVES median = {r:.3f}")
axes[0, 0].legend(frameon=False, fontsize=11, loc="lower right")
fig.suptitle(f"{s.name}   $T_{{\\rm eff}}$={s.teff:.0f} K, "
             f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}"
             f"   --   UVES R$\\approx$41-46k", y=0.995, fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.985))
out = os.path.join(ml.rundir_for(s), s.name + "_uves_compare.png")
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
