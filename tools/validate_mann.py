#!/usr/bin/env python3
"""
External validation of ATLAS12+SYNTHE models against the Mann et al. (2015)
M-dwarf sample (flux-calibrated 0.3-2.4 um spectra + interferometric-quality
Teff/R/M/[Fe/H]).

For a chosen star this script:
  1. pulls Teff, [Fe/H], M, R, d from M_params.fits (logg from M and R),
  2. runs atlas12c.exe from the nearest cool-dwarf starting model,
  3. runs synthe.exe over the observed wavelength range,
  4. converts the model to f_lambda at Earth  [ 4*pi*H_nu * c/lambda^2 * (R/d)^2 ],
     convolves to the SNIFS/SpeX resolution, and plots model vs. data.

The comparison is ABSOLUTE (no renormalization): the (R/d)^2 dilution uses the
measured radius and distance, so the plot tests both the spectral shape and
the overall flux scale.  An FBOL cross-check is printed.

Usage:
  python3 validate_mann.py --list                    # show the sample
  python3 validate_mann.py --star PM_I11054+4331     # run one star
  python3 validate_mann.py --star ... --force        # rerun atlas12+synthe

Runs land in workdir/mann/<star>/ ; the plot is <star>_compare.png there.
"""

import argparse
import os
import subprocess
import sys

import numpy as np
from astropy.io import fits

REPO      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANN_DIR  = os.path.expanduser("~/sps/SPECTRA/Mann")
PARAMS    = os.path.join(MANN_DIR, "M_params.fits")
RUN_ROOT  = os.path.join(REPO, "workdir", "mann")

# Starting models for atlas12 (SCALE_MODEL regrids to the target Teff/logg).
START_MODELS = [
    (2500.0, os.path.join(REPO, "workdir", "czc17_2500.atm")),
    (2800.0, os.path.join(REPO, "workdir", "czc19_2800c.atm")),
    (3500.0, os.path.join(REPO, "workdir", "czc19_3500.atm")),
    (5777.0, os.path.join(REPO, "workdir", "czc17_solar.atm")),
]

# Physical constants (cgs)
GRAV   = 6.67430e-8
MSUN   = 1.98892e33
RSUN   = 6.957e10
PC     = 3.0856776e18
C_ANG  = 2.99792458e18       # speed of light in Angstrom/s

# Observed resolution: SNIFS constant-dlam (see SNIFS_DLAM_FWHM below),
# SpeX SXD R=2000 in the NIR; 2016-library splice at 9500 A
OBS_R_RED   = 2000.0
OBS_R_SPLIT = 9500.0         # Angstrom

# Reference grid computed with the original Kurucz f77 programs (same line
# lists): solar [Fe/H], logg 4.50 only, native R=300000.
C3K_DIR = os.path.expanduser(
    "~/kurucz/grids/c3k_v2.3/at12_feh+0.00_afe+0.0/spec")
C3K_R   = 300000.0


def load_mann():
    """Return the Mann table as a dict of 1-D arrays (183 stars)."""
    d = fits.open(PARAMS)[1].data[0]
    name = np.array([n.strip() for n in d["NAME"]])
    mass, radius = d["MASS"], d["RADIUS"]
    logg = np.log10(GRAV * mass * MSUN / (radius * RSUN) ** 2)
    return dict(name=name, spt=np.array([s.strip() for s in d["SPT"]]),
                teff=d["TEFF"], feh=d["FEH"], mh=d["MH"], mass=mass,
                radius=radius, dist=d["DISTANCE"], fbol=d["FBOL"],
                lum=d["LUMINOSITY"], logg=logg)


def list_sample(t):
    order = np.argsort(t["teff"])[::-1]
    print(f"{'name':<18}{'SpT':<8}{'Teff':>6}{'logg':>7}{'[Fe/H]':>8}"
          f"{'R/Rsun':>8}{'d(pc)':>8}  spectrum")
    for i in order:
        has = os.path.isfile(os.path.join(MANN_DIR, t["name"][i] + ".ascii"))
        print(f"{t['name'][i]:<18}{t['spt'][i]:<8}{t['teff'][i]:>6.0f}"
              f"{t['logg'][i]:>7.2f}{t['feh'][i]:>8.2f}{t['radius'][i]:>8.3f}"
              f"{t['dist'][i]:>8.2f}  {'yes' if has else 'NO'}")


def run(cmd, cwd, logfile):
    env = dict(os.environ, ATLAS12=REPO)
    with open(logfile, "w") as log:
        proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=log,
                              stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        sys.exit(f"error: {' '.join(cmd)} failed (rc={proc.returncode}); "
                 f"see {logfile}")


def iter_convergence(iter_file):
    """Return (rms, max_abs) of the flux ERROR column, final iteration."""
    blocks, cur = [], []
    with open(iter_file) as fh:
        for line in fh:
            if "log10TAU" in line:
                if cur:
                    blocks.append(cur)
                cur = []
                continue
            parts = line.split()
            if not parts or not parts[0].lstrip("-").isdigit():
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            if len(vals) >= 13:
                cur.append(vals[7])
    if cur:
        blocks.append(cur)
    if not blocks:
        return np.nan, np.nan
    err = np.array(blocks[-1])
    return float(np.sqrt(np.mean(err ** 2))), float(np.max(np.abs(err)))


def load_c3k(teff):
    """Load the C3K v2.3 grid spectrum closest in Teff.

    Returns (teff_grid, wave_A, hnu) or None if the grid is absent.
    """
    import re
    if not os.path.isdir(C3K_DIR):
        return None
    pts = []
    for f in sorted(os.listdir(C3K_DIR)):
        m = re.search(r"t(\d{5})g", f)
        if m:
            pts.append((float(m.group(1)), os.path.join(C3K_DIR, f)))
    if not pts:
        return None
    tgrid, path = min(pts, key=lambda p: abs(p[0] - teff))
    w, hnu, _ = np.loadtxt(path, unpack=True)
    return tgrid, w, hnu


def smooth_to_R(flux, model_R, obs_R):
    """Gaussian-smooth a constant-R log-lambda spectrum to resolving power obs_R."""
    from scipy.ndimage import gaussian_filter1d
    sigma_pix = model_R / (obs_R * 2.3548)
    return gaussian_filter1d(flux, sigma_pix, mode="nearest")


# SNIFS effective LSF: constant-FWHM Gaussian, measured 2026-08-09 by
# chunk-fitting R=300k syntheses to three stars (GJ887, PM_I10113+4927,
# GJ105A): dlam ~ 11 A FWHM, i.e. R rising ~400 (blue) -> ~850 (9500 A).
# The old flat R=1000 was ~2x too sharp in the blue.
SNIFS_DLAM_FWHM = 11.0     # Angstrom


def smooth_mann(wave, flux, model_R):
    """Smooth a constant-R model to the Mann-library LSF: SNIFS
    constant-dlam Gaussian below OBS_R_SPLIT, SpeX R=2000 above."""
    from scipy.ndimage import gaussian_filter1d
    step = wave[0] / model_R                    # finest linear spacing
    wl = np.arange(wave[0], min(wave[-1], OBS_R_SPLIT + 200.0), step)
    fl = np.interp(wl, wave, flux)
    fsnifs = gaussian_filter1d(fl, SNIFS_DLAM_FWHM / 2.3548 / step,
                               mode="nearest")
    out = smooth_to_R(flux, model_R, OBS_R_RED)
    m = wave < OBS_R_SPLIT
    out[m] = np.interp(wave[m], wl, fsnifs)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--star", help="Mann star name, e.g. PM_I11054+4331")
    ap.add_argument("--list", action="store_true", help="list the sample and exit")
    ap.add_argument("--vturb", type=float, default=1.0, help="microturbulence km/s")
    ap.add_argument("--numit", type=int, default=30, help="atlas12 iterations")
    ap.add_argument("--resolu", type=float, default=300000.0,
                    help="synthe resolving power before smoothing "
                         "(>=300000: R=50k over-absorbs the molecular haze "
                         "by ~10 percent median in the optical, see "
                         "GJ887_resolution_test.png 2026-08-09)")
    ap.add_argument("--wlbeg", type=float, default=300.0, help="nm")
    ap.add_argument("--wlend", type=float, default=2500.0, help="nm")
    ap.add_argument("--use-mh", action="store_true",
                    help="scale metals by [M/H] instead of [Fe/H]")
    ap.add_argument("--solar", default="berg25",
                    help="solar reference pattern for atlas12 (ag89, agss09, "
                         "berg25; default berg25). Pass an empty string to "
                         "keep the table carried by the input model.")
    ap.add_argument("--theta-mas", type=float, default=None,
                    help="interferometric limb-darkened angular diameter "
                         "[mas].  Overrides Teff (from table FBOL + theta via "
                         "Stefan-Boltzmann), radius (theta/2 * d), logg "
                         "(table mass + that radius), and hence the flux "
                         "dilution -- a fully model-free comparison.  "
                         "Default tag becomes 'int'.")
    ap.add_argument("--tag", default=None,
                    help="suffix for the run directory (workdir/mann/<star>_<tag>)")
    ap.add_argument("--force", action="store_true", help="rerun even if outputs exist")
    args = ap.parse_args()

    t = load_mann()
    if args.list or not args.star:
        list_sample(t)
        return

    idx = np.where(t["name"] == args.star)[0]
    if len(idx) == 0:
        sys.exit(f"error: star {args.star!r} not in M_params.fits (see --list)")
    i = idx[0]

    teff, logg, feh = t["teff"][i], t["logg"][i], t["feh"][i]
    mh = t["mh"][i]
    zmet = mh if args.use_mh else feh
    radius, dist = t["radius"][i], t["dist"][i]
    if not np.isfinite([teff, logg, zmet, radius, dist]).all():
        sys.exit(f"error: {args.star} has missing parameters "
                 f"(Teff={teff}, logg={logg}, z={zmet}, R={radius}, d={dist})")

    obs_file = os.path.join(MANN_DIR, args.star + ".ascii")
    if not os.path.isfile(obs_file):
        sys.exit(f"error: no observed spectrum {obs_file}")

    print(f"{args.star} ({t['spt'][i]}):  Teff={teff:.0f} K  logg={logg:.2f}  "
          f"[Fe/H]={feh:+.2f}  [M/H]={mh:+.2f}  R={radius:.3f} Rsun  d={dist:.2f} pc")

    if args.theta_mas:
        SIGMA_SB = 5.670374419e-5                       # cgs
        theta_rad = args.theta_mas * 1.0e-3 / 206264.806
        fbol = t["fbol"][i] * 1.0e-8                    # Mann convention
        teff = (4.0 * fbol / (SIGMA_SB * theta_rad ** 2)) ** 0.25
        radius = (theta_rad / 2.0) * dist * PC / RSUN
        logg = np.log10(GRAV * t["mass"][i] * MSUN / (radius * RSUN) ** 2)
        print(f"INTERFEROMETRIC MODE: theta_LD={args.theta_mas:.3f} mas  ->  "
              f"Teff={teff:.0f} K (table {t['teff'][i]:.0f})  "
              f"R={radius:.3f} Rsun (table {t['radius'][i]:.3f})  "
              f"logg={logg:.2f}")
        if args.tag is None:
            args.tag = "int"

    # --------------- run atlas12 ---------------
    rundir = os.path.join(RUN_ROOT,
                          args.star + (f"_{args.tag}" if args.tag else ""))
    os.makedirs(rundir, exist_ok=True)
    base = args.star
    atm = os.path.join(rundir, base + ".atm")

    start_teff, start_atm = min(START_MODELS,
                                key=lambda m: abs(np.log(m[0] / teff)))
    if args.force or not os.path.isfile(atm):
        print(f"running atlas12 from {os.path.basename(start_atm)} "
              f"(start Teff {start_teff:.0f} K) ...")
        cmd = [os.path.join(REPO, "bin", "atlas12c.exe"), start_atm, base,
               f"teff={teff:.0f}", f"logg={logg:.2f}",
               f"zscale={10.0 ** zmet:.4f}", f"vturb={args.vturb}",
               f"numit={args.numit}"]
        if args.solar:
            cmd.append(f"solar={args.solar}")
        run(cmd, rundir, os.path.join(rundir, base + ".stdout"))
    else:
        print("atlas12: using existing .atm (use --force to rerun)")
    rms, emax = iter_convergence(os.path.join(rundir, base + ".iter"))
    print(f"atlas12 convergence (final iter): flux-error RMS={rms:.2f}%  "
          f"max={emax:.2f}%")

    # --------------- run synthe ---------------
    spec_file = os.path.join(rundir, base + ".spec")
    if args.force or not os.path.isfile(spec_file):
        print(f"running synthe {args.wlbeg:.0f}-{args.wlend:.0f} nm at "
              f"R={args.resolu:.0f} ...")
        cmd = [os.path.join(REPO, "bin", "synthe.exe"), atm,
               f"wlbeg={args.wlbeg}", f"wlend={args.wlend}",
               f"resolu={args.resolu:.0f}"]
        run(cmd, rundir, os.path.join(rundir, base + ".synthe.stdout"))
    else:
        print("synthe: using existing .spec (use --force to rerun)")

    # --------------- model -> f_lambda at Earth ---------------
    wmod, hnu, hcont = np.loadtxt(spec_file, unpack=True)     # Angstrom, H_nu
    dilute = (radius * RSUN / (dist * PC)) ** 2
    flam = 4.0 * np.pi * hnu * C_ANG / wmod ** 2 * dilute     # erg/s/cm^2/A

    fsm = smooth_mann(wmod, flam, args.resolu)

    # --------------- observed spectrum ---------------
    wobs_um, fobs, eobs = np.loadtxt(obs_file, unpack=True)
    wobs = wobs_um * 1.0e4                                     # Angstrom
    ok = (wobs >= wmod[0]) & (wobs <= wmod[-1]) & (fobs > 0)
    fmod_i = np.interp(wobs, wmod, fsm)

    # --------------- reference C3K grid spectrum (f77 Kurucz) ---------------
    c3k = load_c3k(teff)
    if c3k is not None:
        c3k_teff, wc3k, hnu_c3k = c3k
        flam_c3k = 4.0 * np.pi * hnu_c3k * C_ANG / wc3k ** 2 * dilute
        csm = smooth_mann(wc3k, flam_c3k, C3K_R)
        okc = ok & (wobs >= wc3k[0]) & (wobs <= wc3k[-1])
        fc3k_i = np.interp(wobs, wc3k, csm)
        c3k_label = f"C3K v2.3 f77 (t{c3k_teff:.0f} g4.5 [Fe/H]=0)"
        print(f"reference grid point: {c3k_label}")

    # --------------- metrics ---------------
    fbol_obs = np.trapz(fobs, wobs)
    fbol_mod = np.trapz(np.interp(wobs, wmod, flam), wobs)
    print(f"integral over obs range:  data={fbol_obs:.3e}  model={fbol_mod:.3e} "
          f"erg/s/cm^2  (model/data = {fbol_mod / fbol_obs:.3f})")
    print(f"table FBOL = {t['fbol'][i]:.3e}")
    for lo, hi, lab in [(4000, 9500, "optical 0.40-0.95 um"),
                        (9500, 24000, "NIR     0.95-2.40 um")]:
        m = ok & (wobs > lo) & (wobs < hi)
        if m.any():
            r = np.median(fmod_i[m] / fobs[m])
            line = f"median model/data, {lab}: {r:.3f}"
            if c3k is not None:
                mc = okc & (wobs > lo) & (wobs < hi)
                line += f"   (c3k: {np.median(fc3k_i[mc] / fobs[mc]):.3f})"
            print(line)

    # --------------- plot ---------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "figure.dpi": 300,
    })

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 6.5), sharex=True,
        gridspec_kw={"height_ratios": [2.6, 1], "hspace": 0.06})

    sc = 10.0 ** np.floor(np.log10(np.nanmax(fobs[ok])))
    ax1.plot(wobs[ok] / 1e4, fobs[ok] / sc, color="0.25", lw=0.8, label="Mann+15")
    ax1.plot(wobs[ok] / 1e4, fmod_i[ok] / sc, color="crimson", lw=0.8, alpha=0.85,
             label="ATLAS12+SYNTHE")
    if c3k is not None:
        ax1.plot(wobs[okc] / 1e4, fc3k_i[okc] / sc, color="steelblue", lw=0.8,
                 alpha=0.75, label=c3k_label)
    ax1.set_ylabel(rf"$f_\lambda$  [$10^{{{int(np.log10(sc))}}}$"
                   r" erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
    ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc="upper right")
    ax1.set_title(f"{args.star}  ({t['spt'][i]})   "
                  f"$T_{{\\rm eff}}$={teff:.0f} K, log g={logg:.2f}, "
                  f"[Fe/H]={feh:+.2f}", fontsize=11)

    ax2.axhline(1.0, color="0.6", lw=0.7)
    if c3k is not None:
        ax2.plot(wobs[okc] / 1e4, fc3k_i[okc] / fobs[okc], color="steelblue",
                 lw=0.7, alpha=0.75)
    ax2.plot(wobs[ok] / 1e4, fmod_i[ok] / fobs[ok], color="crimson", lw=0.7)
    ax2.set_ylim(0.5, 1.5)
    ax2.set_xlabel(r"wavelength  [$\mu$m]")
    ax2.set_ylabel("model / data")

    out = os.path.join(rundir, base + "_compare.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"plot: {out}")


if __name__ == "__main__":
    main()
