#!/usr/bin/env python3
"""Na I D: stacked panels comparing LTE against NLTE (or against each other).

Each panel is one star.  The first curve in a panel is the reference and is
drawn in blue; the rest step away through the Okabe-Ito sequence.

PANEL SYNTAX
------------
    --panel "TITLE|opt|opt|LABEL:FILE|LABEL:FILE|..."

Fields are separated by '|'.  The first is the panel title.  A field of the
form key=value (or a bare 'log') is an option; anything else is a curve given
as LABEL:PATH.  Options:

    ylim=LO,HI    y limits.  Set them per panel: a cool dwarf's Na D window
                  never comes near the continuum, so the 0-1 range that suits
                  the Sun renders it as a flat black floor.
    log           logarithmic y axis
    note=TEXT     small grey annotation inside the panel

Example:
    python3 tools/nad_nlte_plot.py out.png \\
      --panel "Sun|ylim=0,1.09|LTE:sun_lte.spec|NLTE:sun_nlte.spec" \\
      --panel "M dwarf|ylim=0,0.075|LTE:m_lte.spec|NLTE:m_nlte.spec" \\
      --title "..." --subtitle "..."

Files are SYNTHE .spec (vacuum wavelength in A, flux, continuum flux); the
plotted quantity is flux/continuum.
"""
import argparse
import textwrap

import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Vacuum wavelengths as they appear in gfall (nm x 10).  The doublet is
# universally quoted in air (5889.95 / 5895.92); SYNTHE works in vacuum
# throughout and .spec is vacuum, so the vacuum values are used here and the
# axis is labelled accordingly.
D2_VAC = 5891.583
D1_VAC = 5897.558

# Okabe-Ito.  Blue is always the reference curve; anything plotted against it
# steps away through orange, reddish purple and bluish green, all of which
# stay distinguishable under the common colour-vision deficiencies.
PALETTE = ["#0072B2", "#E69F00", "#CC79A7", "#009E73", "#56B4E9"]


def load(path, smooth_to=None, native_r=300000.0):
    """Wavelength (A, vacuum) and F/Fcont, optionally at reduced resolution.

    SYNTHE writes onto a LOGARITHMIC wavelength grid with ratio 1 + 1/R, so a
    fixed resolving power is a fixed number of PIXELS everywhere -- the kernel
    width does not have to track lambda.  Degrading R_native to R_out needs a
    Gaussian of FWHM = R_native/R_out pixels.

    Flux and continuum are convolved separately and then divided.  The
    continuum is smooth so this barely differs from smoothing the ratio, but
    it is what an instrument actually does to each.

    Note this is the correct way round: synthesising directly on a coarse grid
    is not the same thing and biases cool-star spectra (it over-absorbs the
    M-dwarf optical by ~10 per cent), which is why the syntheses are always run
    at R >= 3e5 and only smoothed afterwards.
    """
    d = np.loadtxt(path)
    w, f, c = d[:, 0], d[:, 1], d[:, 2]
    if smooth_to:
        sigma_pix = (native_r / smooth_to) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        f = gaussian_filter1d(f, sigma_pix, mode="nearest")
        c = gaussian_filter1d(c, sigma_pix, mode="nearest")
    return w, f / c


def parse_panel(spec):
    """'TITLE|opt|LABEL:FILE|...' -> (title, curves, opts)."""
    fields = [f.strip() for f in spec.split("|") if f.strip()]
    title, rest = fields[0], fields[1:]
    curves, opts = [], {"ylim": None, "log": False, "note": None}
    for f in rest:
        if f == "log":
            opts["log"] = True
        elif f.startswith("ylim="):
            opts["ylim"] = tuple(float(v) for v in f[5:].split(","))
        elif f.startswith("note="):
            opts["note"] = f[5:]
        else:
            label, path = f.split(":", 1)
            curves.append((label, path))
    if not curves:
        raise SystemExit(f"panel '{title}' has no curves")
    return title, curves, opts


def draw(ax, title, curves, opts, wlo, whi, label_lines, smooth_to, native_r):
    for i, (label, path) in enumerate(curves):
        w, r = load(path, smooth_to, native_r)
        m = (w >= wlo) & (w <= whi)
        ax.plot(w[m], r[m], color=PALETTE[i % len(PALETTE)],
                lw=1.9 if i == 0 else 1.5, label=label)

    # The guide lines run through every panel, but only the top one is
    # labelled -- lower down the text would land on the panel title.
    for wl, name in ((D2_VAC, "D$_2$"), (D1_VAC, "D$_1$")):
        ax.axvline(wl, color="0.75", lw=0.7, zorder=0)
        if label_lines:
            ax.annotate(name, xy=(wl, 1.005), xycoords=("data", "axes fraction"),
                        ha="center", va="bottom", fontsize=11, color="0.35",
                        annotation_clip=False)

    ax.set_xlim(wlo, whi)
    if opts["log"]:
        ax.set_yscale("log")
    if opts["ylim"]:
        ax.set_ylim(*opts["ylim"])
    ax.set_ylabel("$F_\\lambda / F_\\mathrm{cont}$")
    ax.set_title(title, loc="left", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if opts["note"]:
        ax.text(0.985, 0.06, opts["note"], transform=ax.transAxes, ha="right",
                va="bottom", fontsize=10.5, color="0.35")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--panel", action="append", required=True,
                    metavar="TITLE|opt|LABEL:FILE|...")
    ap.add_argument("--wlo", type=float, default=5886.0)
    ap.add_argument("--whi", type=float, default=5903.0)
    ap.add_argument("--title", default="")
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--panel-height", type=float, default=3.7)
    ap.add_argument("--smooth-to", type=float, default=None,
                    help="degrade to this resolving power (FWHM)")
    ap.add_argument("--native-r", type=float, default=300000.0,
                    help="resolving power the .spec files were synthesised at")
    ap.add_argument("--no-line-labels", action="store_true",
                    help="draw the D1/D2 guide lines without their labels")
    args = ap.parse_args()

    panels = [parse_panel(p) for p in args.panel]

    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 13,
        "axes.grid": False, "figure.dpi": 300,
    })

    n = len(panels)
    fig, axes = plt.subplots(n, 1, figsize=(9.5, args.panel_height * n),
                             sharex=True)
    if n == 1:
        axes = [axes]

    for i, (ax, (title, curves, opts)) in enumerate(zip(axes, panels)):
        draw(ax, title, curves, opts, args.wlo, args.whi,
             label_lines=(i == 0 and not args.no_line_labels),
             smooth_to=args.smooth_to, native_r=args.native_r)

    axes[0].legend(frameon=False, loc="lower left", fontsize=11)
    axes[-1].set_xlabel("vacuum wavelength [$\\mathrm{\\AA}$]")

    if args.title:
        fig.suptitle(args.title, fontsize=13.5, y=0.985)
    # Wrap rather than trusting the caller to pre-break it: a long caption
    # passed on one line silently runs off both edges of the figure.
    if args.subtitle:
        sub = "\n".join(textwrap.wrap(args.subtitle, 105))
        nl = sub.count("\n") + 1
        fig.text(0.5, 0.004, sub, ha="center", va="bottom", fontsize=10,
                 color="0.35", linespacing=1.5)
        bot = 0.012 + 0.020 * nl * (2.0 / n)
    else:
        bot = 0.0
    fig.tight_layout(rect=[0, bot, 1, 0.97 if args.title else 1.0])
    fig.savefig(args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
