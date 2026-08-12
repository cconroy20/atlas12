#!/usr/bin/env python3
"""GJ644C (vB 8, M7V, 2700 K): Mann+15 data vs ATLAS12+SYNTHE vs PHOENIX
NewEra, focused on the classic 0.40-0.45 um blue mismatch.

Panel 1: 0.35-0.75 um, log f_lambda, with the model continuum for scale.
Panel 2: 0.38-0.50 um zoom, linear, with the data 1-sigma band.
Panel 3: model/data ratio, log scale, 0.38-0.75 um.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/Users/cconroy/kurucz/atlas12/tools")
import mann_lib as ml
import mann_compare_plot as mcp

STAR = sys.argv[1] if len(sys.argv) > 1 else "GJ644C"
TAG = sys.argv[2] if len(sys.argv) > 2 else None

s = ml.resolve(STAR)
rundir = ml.rundir_for(s, TAG)
SPEC = os.path.join(rundir, s.name + ".spec")
OUT = os.path.join(os.path.dirname(rundir.rstrip("/")),
                   s.name + "_blue_compare.png")
wobs, fobs, eobs = ml.read_spectrum(s)
fsurf, esurf = fobs / s.dilute, eobs / s.dilute

w, f, c = ml.read_spec_file(SPEC)
fmod = mcp._model_on_obs(wobs, w, f, mcp.RSYN, s.obs_split)
cmod = mcp._model_on_obs(wobs, w, c, mcp.RSYN, s.obs_split)
fphx = mcp._phoenix_on_obs(wobs, s)

# ------------------------------------------------------------- band table
print(f"{s}")
print(f"{'window [A]':>13} {'S/N':>6} {'ours/data':>10} {'PHX/data':>9} "
      f"{'ours/PHX':>9} {'blanketing':>11}")
for lo, hi in [(3800, 4000), (4000, 4500), (4500, 5000), (5000, 5500),
               (5500, 6000), (6000, 7000), (7000, 8000)]:
    m = (wobs > lo) & (wobs < hi)
    d, ph = np.mean(fsurf[m]), np.mean(fphx[m])
    print(f"{lo:6d}-{hi:6d} {np.mean(fobs[m]/eobs[m]):6.1f} "
          f"{np.mean(fmod[m])/d:10.2f} {ph/d:9.2f} "
          f"{np.mean(fmod[m])/ph:9.2f} {np.mean(cmod[m]/fmod[m]):11.1f}")

# ------------------------------------------------------------------ figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                     "font.size": 14, "axes.grid": False, "figure.dpi": 600})

C_DATA, C_OURS, C_PHX, C_CONT = "0.15", "#0072B2", "#E69F00", "#56B4E9"
wum = wobs / 1e4

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(13, 12))

m = (wum > 0.35) & (wum < 0.75)
ax1.plot(wum[m], cmod[m], color=C_CONT, lw=0.9, ls="--", zorder=2,
         label="ATLAS12 continuum")
ax1.plot(wum[m], fsurf[m], color=C_DATA, lw=0.8, zorder=0, label="Mann data")
ax1.plot(wum[m], fphx[m], color=C_PHX, lw=0.9, zorder=1,
         label="PHOENIX NewEra")
ax1.plot(wum[m], fmod[m], color=C_OURS, lw=0.9, zorder=3,
         label="ATLAS12+SYNTHE")
ax1.set_yscale("log")
ax1.set_xlim(0.35, 0.75)
ax1.set_ylim(1e2, 3e5)
ax1.legend(frameon=False, loc="lower right", fontsize=12, ncol=2)
ax1.set_title(f"{s.name} (vB 8, M7V)   $T_{{\\rm eff}}$={s.teff:.0f} K, "
              f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}", fontsize=14)

m = (wum > 0.385) & (wum < 0.47)
ax2.fill_between(wum[m], fsurf[m] - esurf[m], fsurf[m] + esurf[m],
                 color="0.7", lw=0, zorder=0)
ax2.plot(wum[m], fsurf[m], color=C_DATA, lw=1.2, zorder=1, label="Mann data")
ax2.plot(wum[m], fphx[m], color=C_PHX, lw=1.2, zorder=2,
         label="PHOENIX NewEra")
ax2.plot(wum[m], fmod[m], color=C_OURS, lw=1.2, zorder=3,
         label="ATLAS12+SYNTHE")
ax2.set_xlim(0.385, 0.47)
ax2.set_ylim(0, 6500)
# chromospheric emission in this active M7Ve -- data spikes, not continuum
for lab, wl in [("Ca II K", 3934.8), ("Ca II H", 3969.6), ("H$\\delta$", 4102.9),
                ("H$\\gamma$", 4341.7)]:
    k = (wobs > wl - 15) & (wobs < wl + 15)
    ax2.annotate(lab, xy=(wl / 1e4, fsurf[k].max() + 150), ha="center",
                 va="bottom", fontsize=10, color="0.15", rotation=90)

for ax in (ax1, ax2):
    ax.set_ylabel(r"$f_\lambda^{\rm surf}$"
                  r"  [erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")


def binned(x, y, n=200):
    e = np.linspace(x.min(), x.max(), n + 1)
    i = np.clip(np.digitize(x, e) - 1, 0, n - 1)
    xb = np.array([np.mean(x[i == k]) if (i == k).any() else np.nan
                   for k in range(n)])
    yb = np.array([np.mean(y[i == k]) if (i == k).any() else np.nan
                   for k in range(n)])
    return xb, yb


mr = (wum > 0.38) & (wum < 0.75)
for lab, fc, col in [("PHOENIX NewEra / data", fphx, C_PHX),
                     ("ATLAS12+SYNTHE / data", fmod, C_OURS)]:
    xb, yb = binned(wum[mr], fc[mr] / fsurf[mr])
    ax3.plot(xb, yb, color=col, lw=1.5, label=lab)
ax3.axhline(1.0, color="0.4", lw=0.8, ls=":")
ax3.set_yscale("log")
ax3.set_xlim(0.38, 0.75)
ax3.set_ylim(0.4, 20)
ax3.set_yticks([0.5, 1, 2, 5, 10, 20])
ax3.set_yticklabels(["0.5", "1", "2", "5", "10", "20"])
ax3.legend(frameon=False, loc="upper right", fontsize=12)
ax3.set_ylabel("model / data")

for ax in (ax1, ax2, ax3):
    ax.set_xlabel(r"wavelength  [$\mu$m]")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, bbox_inches="tight")
print("wrote", OUT)
