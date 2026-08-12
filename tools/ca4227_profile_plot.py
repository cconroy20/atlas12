#!/usr/bin/env python3
"""What the modified Ca I 4227 profile does to the far wing.

The enhancement over the impact-approximation Lorentzian,
    H_mod / H_Lorentz = ps * sqrt(pi) * x^(2-px),   x = dlam/dlam_Doppler,
is independent of the damping constant, so it can be drawn without
reference to a particular atmosphere.  Only the conversion from x to
Angstrom needs the Doppler width, taken at 2700 K with xi = 1 km/s.

Usage:  python3 ca4227_profile_plot.py [OUT.png]
"""
import sys

import numpy as np

LAM0 = 4227.918          # A, vacuum
TEFF = 2700.0            # K, line-forming layers of GJ644C
XI = 1.0e5               # cm/s microturbulence
MCA = 40.078 * 1.66054e-24
KB = 1.380649e-16
CLIGHT = 2.99792458e10

VARIANTS = [
    ("Jones+23 as published:  $p_s$=1e-3, $p_x$=0.5", 1.0e-3, 0.5, 500.0,
     "#56B4E9"),
    ("fitted to GJ644C:  $p_s$=1.2e-5, $p_x$=0.5", 1.2e-5, 0.5, 750.0,
     "#CC79A7"),
]

out = sys.argv[1] if len(sys.argv) > 1 else "ca4227_profile.png"

vdop = np.sqrt(2.0 * KB * TEFF / MCA + XI ** 2)
dlam_D = LAM0 * vdop / CLIGHT            # Angstrom
print(f"Doppler width at {TEFF:.0f} K, xi={XI/1e5:.1f} km/s: "
      f"{dlam_D:.4f} A  ({vdop/1e5:.2f} km/s)")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                     "font.size": 14, "axes.grid": False, "figure.dpi": 600})

dl = np.logspace(-1, 3.1, 600)           # Angstrom
x = dl / dlam_D

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.axhline(1.0, color="0.15", lw=1.6, label="Lorentzian far wing (baseline)")
for lab, ps, px, dlmax, col in VARIANTS:
    enh = ps * np.sqrt(np.pi) * x ** (2.0 - px)
    enh = np.maximum(enh, 1.0)           # SYNTHE takes MAX(Voigt, modified)
    enh[dl > dlmax] = np.nan
    ax.plot(dl, enh, color=col, lw=1.8, label=lab)
    k = np.nanargmax(np.where(dl <= dlmax, dl, np.nan))
    ax.plot([dl[k]], [enh[k]], marker="o", ms=5, color=col)
    ax.annotate(f"cut at {dlmax:.0f} " + r"$\mathrm{\AA}$",
                xy=(dl[k], enh[k]), xytext=(6, 4), textcoords="offset points",
                fontsize=10, color=col)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.5, 1300)
ax.set_ylim(0.5, 3e4)
ax.set_xlabel(r"$|\lambda - 4227.9|$  [$\mathrm{\AA}$]")
ax.set_ylabel("opacity relative to Lorentzian wing")
ax.set_title("Ca I 4227: what the modified profile adds", fontsize=14)
ax.legend(frameon=False, fontsize=11, loc="upper left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

secax = ax.secondary_xaxis("top", functions=(lambda d: d / dlam_D,
                                             lambda v: v * dlam_D))
secax.set_xlabel(r"$x = \delta\lambda / \delta\lambda_{\rm Doppler}$",
                 fontsize=12)

fig.tight_layout()
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
