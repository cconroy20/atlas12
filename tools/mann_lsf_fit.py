#!/usr/bin/env python3
"""R_eff(lambda) of the Mann spectra, measured by chunk-fitting reference
models smoothed to trial R.  Star-generic: run once per star; overlays
ATLAS R=300k reference, PHOENIX reference (GJ887 only), strong lines, and
the assumed R=1000/2000 convention.  Writes <star>_data_resolution.png and
a JSON of the measurements for the cross-star overlay."""
import json
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml

REPO = ml.REPO
_S = ml.resolve(sys.argv[1] if len(sys.argv) > 1 else "GJ887")
STAR = _S.name
# reference model: full-range R=300k synthesis from the standard run dir
# (override with a path as argv[2])
SPEC = (sys.argv[2] if len(sys.argv) > 2
        else os.path.join(ml.rundir_for(_S), STAR + ".spec"))
OUTDIR = ml.rundir_for(_S)
MANN_DIR = ml.MANN_DIR
C_ANG = ml.C_ANG
RSYN = 300000.0

# dilution irrelevant in principle (every fit has a free scale) but kept
# so overlays sit on the data
_wobs, _fobs, _ = ml.read_spectrum(_S)
WOBS, FOBS = _wobs, _fobs / _S.dilute

w, hnu, _ = np.loadtxt(SPEC, unpack=True)
fmod = 4 * np.pi * hnu * C_ANG / w ** 2

R_TRIALS = np.unique(np.round(np.geomspace(150, 6000, 44)))
SM = {R: gaussian_filter1d(fmod, RSYN / (R * 2.3548), mode="nearest")
      for R in R_TRIALS}

# velocity nuisance per chunk: per-star rest-frame residual (+-20 km/s,
# Mann's cross-correlation accuracy) + SNIFS wavelength-solution wobble.
# Without it, registration errors bias the fitted width BROAD (the
# 2026-08-09 air-vs-vac episode: the original 11 A SNIFS FWHM).
C_KMS = 2.99792458e5
V_NUIS = np.arange(-100.0, 101.0, 4.0)


def chunk_fit(edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (WOBS > lo) & (WOBS < hi)
        if m.sum() < 40:
            continue
        wd, fd = WOBS[m], FOBS[m]
        chi = np.empty((len(R_TRIALS), len(V_NUIS)))
        for i, R in enumerate(R_TRIALS):
            for k, v in enumerate(V_NUIS):
                fm = np.interp(wd, w * (1.0 + v / C_KMS), SM[R])
                s = np.sum(fd * fm) / np.sum(fm * fm)
                chi[i, k] = float(np.sum((fd - s * fm) ** 2))
        chiR = chi.min(axis=1)                    # profile out velocity
        j = int(np.argmin(chiR))
        kbest = int(np.argmin(chi[j]))
        csc = chiR / chiR[j] * len(wd)
        okr = np.flatnonzero(csc < csc[j] + 2 * np.sqrt(2 * len(wd)))
        out.append([float(np.sqrt(lo * hi)), float(R_TRIALS[j]),
                    float(R_TRIALS[okr[0]]), float(R_TRIALS[okr[-1]]),
                    float(V_NUIS[kbest])])
    return out

rows = chunk_fit(np.geomspace(4200, 24000, 30))

LINES = {"Na D": (5820, 5975), "K I 7699": (7660, 7740),
         "Ca II 8542": (8480, 8610), "Na I 8190": (8150, 8240)}
lrows = []
for lab, (lo, hi) in LINES.items():
    m = (WOBS > lo) & (WOBS < hi)
    wd, fd = WOBS[m], FOBS[m]
    best = (None, np.inf)
    for R in R_TRIALS:
        for v in V_NUIS:
            fm = np.interp(wd, w * (1.0 + v / C_KMS), SM[R])
            s = np.sum(fd * fm) / np.sum(fm * fm)
            rms = float(np.sqrt(np.mean((fd - s * fm) ** 2)))
            if rms < best[1]:
                best = (float(R), rms)
    lrows.append([float(np.sqrt(lo * hi)), best[0], lab])

json.dump({"star": STAR, "chunks": rows, "lines": lrows},
          open(f"{OUTDIR}/resolution_fit.json", "w"),
          indent=1)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "stix", "font.size": 12,
    "axes.grid": False, "figure.dpi": 300,
})
rows = np.array([r for r in rows])
fig, ax = plt.subplots(figsize=(10.5, 4.8))
ax.errorbar(rows[:, 0], rows[:, 1],
            yerr=[rows[:, 1] - rows[:, 2], rows[:, 3] - rows[:, 1]],
            fmt="o", ms=4.5, color="#0072B2", lw=1, capsize=2,
            label=f"chunk fit vs ATLAS R=300k ({STAR})")
for x, R, lab in lrows:
    ax.plot(x, R, marker="*", ms=13, color="#CC79A7", ls="none",
            zorder=5)
    ax.annotate(lab, (x, R), textcoords="offset points", xytext=(4, -14),
                fontsize=9, color="#CC79A7")
ax.plot([], [], marker="*", ms=11, color="#CC79A7", ls="none",
        label="strong-line fits")
ax.plot([4200, 9500], [1000, 1000], color="0.4", lw=1.3, ls="--",
        label="assumed R=1000 / R=2000")
ax.plot([9500, 24000], [2000, 2000], color="0.4", lw=1.3, ls="--")
ax.axvline(9500, color="0.8", lw=0.8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([5000, 7000, 9000, 12000, 16000, 21000])
ax.set_xticklabels(["5000", "7000", "9000", "12000", "16000", "21000"])
ax.set_yticks([200, 300, 500, 1000, 2000, 4000])
ax.set_yticklabels(["200", "300", "500", "1000", "2000", "4000"])
ax.set_xlabel(r"wavelength  [$\mathrm{\AA}$]")
ax.set_ylabel(r"$R_{\rm eff}$")
ax.legend(frameon=False, fontsize=10, loc="upper left")
ax.set_title(f"Mann spectrum effective resolution: {STAR}", fontsize=12)
out = f"{OUTDIR}/{STAR}_data_resolution.png"
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
