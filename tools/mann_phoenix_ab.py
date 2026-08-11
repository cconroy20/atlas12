#!/usr/bin/env python3
"""Band-by-band comparison of one or more ATLAS12+SYNTHE models against the
Mann spectrum AND the PHOENIX NewEra model at the star's parameters.

Unlike tools/cia_ab.py (same-structure, PHOENIX high-res, opacity-only A/B)
this compares FULLY SELF-CONSISTENT models: each spectrum comes from its own
converged structure, and both chains are put on the observed grid through the
measured Mann LSF.  That is the comparison that matters once the structure is
allowed to respond to an opacity change.

Usage:
    python3 tools/mann_phoenix_ab.py GJ644C                    # untagged run
    python3 tools/mann_phoenix_ab.py GJ644C new=. old=preopacity

Each extra argument is label=tag, where tag selects workdir/mann/<STAR>_<tag>
('.' or 'none' means the untagged run directory).  With no such arguments the
untagged run is used alone under the label 'ATLAS12'.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml                                      # noqa: E402
from mann_compare_plot import _model_on_obs, _phoenix_on_obs   # noqa: E402

RSYN = 300000.0
COLORS = ["#0072B2", "#CC79A7", "#56B4E9", "#009E73"]

# (label, lo, hi) in Angstroms.  Optical windows are the molecular-band set
# used throughout the campaign; the near-IR set is the one the CIA / H2- ff
# work is scored on.
WINDOWS = [
    ("TiO blue  5500-5700", 5500., 5700.),
    ("TiO gamma' 5900-6100", 5900., 6100.),
    ("CaH       6800-7000", 6800., 7000.),
    ("TiO gamma 7053-7270", 7053., 7270.),
    ("TiO       7450-7650", 7450., 7650.),
    ("TiO eps   8400-8800", 8400., 8800.),
    ("I/cont    9000-9500", 9000., 9500.),
    ("J         11500-13500", 11500., 13500.),
    ("H2O 1.4um 13800-14600", 13800., 14600.),
    ("H         15500-17000", 15500., 17000.),
    ("H2O 1.9um 18200-19200", 18200., 19200.),
    ("K         21000-23000", 21000., 23000.),
]


def median_ratio(mask, a, b):
    good = mask & np.isfinite(a) & np.isfinite(b) & (b > 0)
    return float(np.median(a[good] / b[good])) if good.sum() else np.nan


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    s = ml.resolve(sys.argv[1])
    specs = []
    for a in sys.argv[2:]:
        label, tag = a.split("=", 1)
        tag = None if tag in (".", "none", "") else tag
        specs.append((label, os.path.join(ml.rundir_for(s, tag),
                                          s.name + ".spec")))
    if not specs:
        specs = [("ATLAS12", os.path.join(ml.rundir_for(s),
                                          s.name + ".spec"))]

    wobs, fobs, _ = ml.read_spectrum(s)
    data = fobs / s.dilute                       # -> surface flux
    print(f"{s.name}  Teff={s.teff:.0f}  logg={s.logg:.2f}  "
          f"[Fe/H]={s.feh:+.2f}   ({len(wobs)} obs points "
          f"{wobs[0]:.0f}-{wobs[-1]:.0f} A)")

    models = []
    for label, path in specs:
        w, f, _ = ml.read_spec_file(path)
        models.append((label, _model_on_obs(wobs, w, f, RSYN, s.obs_split)))
        print(f"  {label:10s} {path}")
    phx = _phoenix_on_obs(wobs, s)

    # ------------------------------------------------------------------ table
    head = (f"{'window':22s}" + "".join(f"{lab[:9]+'/data':>16s}"
                                        for lab, _ in models)
            + f"{'PHOENIX/data':>16s}"
            + "".join(f"{lab[:9]+'/PHX':>16s}" for lab, _ in models))
    print()
    print(head)
    print("-" * len(head))
    for lab, lo, hi in WINDOWS:
        m = (wobs >= lo) & (wobs <= hi)
        if not m.any():
            continue
        row = f"{lab:22s}"
        for _, fm in models:
            row += f"{median_ratio(m, fm, data):16.3f}"
        row += f"{median_ratio(m, phx, data):16.3f}"
        for _, fm in models:
            row += f"{median_ratio(m, fm, phx):16.3f}"
        print(row)

    # broad integrals
    print()
    for lo, hi, lab in [(4000, 9500, "optical 0.40-0.95"),
                        (9500, 24000, "NIR     0.95-2.40")]:
        m = (wobs >= lo) & (wobs <= hi)
        if not m.any():
            continue
        parts = [f"{lab}: "]
        for mlab, fm in models:
            parts.append(f"{mlab} {median_ratio(m, fm, data):.3f}")
        parts.append(f"PHOENIX {median_ratio(m, phx, data):.3f}")
        print("  median model/data,  " + "   ".join(parts))
    for mlab, fm in models:
        m = np.isfinite(fm)
        print(f"  integral {mlab}/data = "
              f"{np.trapz(fm[m], wobs[m]) / np.trapz(data[m], wobs[m]):.3f}")
    m = np.isfinite(phx)
    print("  integral PHOENIX/data = "
          f"{np.trapz(phx[m], wobs[m]) / np.trapz(data[m], wobs[m]):.3f}")

    # ----------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 14,
        "axes.grid": False, "figure.dpi": 300,
    })
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    wum = wobs / 1.0e4
    sc = 10.0 ** np.floor(np.log10(np.nanmax(data)))
    ax1.plot(wum, data / sc, color="0.2", lw=0.6, zorder=0, label="Mann data")
    ax1.plot(wum, phx / sc, color="#E69F00", lw=0.6, zorder=1,
             label="PHOENIX NewEra")
    for k, (lab, fm) in enumerate(models):
        ax1.plot(wum, fm / sc, color=COLORS[k % len(COLORS)], lw=0.6,
                 zorder=5 - k, label=lab)
    ax1.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                   r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
    ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc="upper right", fontsize=12)
    ax1.set_title(f"{s.name}   $T_{{\\rm eff}}$={s.teff:.0f} K, "
                  f"log g={s.logg:.2f}, [Fe/H]={s.feh:+.2f}", fontsize=14)

    ax2.axhline(1.0, color="0.2", lw=0.8)
    ax2.plot(wum, phx / data, color="#E69F00", lw=0.5)
    for k, (lab, fm) in enumerate(models):
        ax2.plot(wum, fm / data, color=COLORS[k % len(COLORS)], lw=0.5)
    ax2.set_ylim(0.3, 1.7)
    ax2.set_ylabel("model / data")
    ax2.set_xlabel(r"wavelength  [$\mu$m]")
    ax2.set_xscale("log")
    ax2.set_xlim(0.4, 2.45)
    ax2.set_xticks([0.4, 0.6, 0.8, 1.0, 1.5, 2.0])
    ax2.set_xticklabels(["0.4", "0.6", "0.8", "1.0", "1.5", "2.0"])
    fig.tight_layout()
    out = os.path.join(ml.rundir_for(s), s.name + "_phoenix_ab.png")
    fig.savefig(out, bbox_inches="tight")
    print("\nwrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
