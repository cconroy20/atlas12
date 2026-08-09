#!/usr/bin/env python3
"""Model-free validation against the Mann+13 interferometric CALIBRATOR
spectra (public Calibrators.zip; SNIFS+SpeX absolute f_lambda, 0.33-3 um)
-- the K-dwarf extension of the interferometric campaign.

All comparison inputs are measured: Teff from --fbol + --theta-mas via
Stefan-Boltzmann, dilution = (theta/2)^2, logg supplied (Gaia benchmark /
literature).  Runs land in workdir/mann/<name>_int/ like the M-dwarf
campaign.

Usage:
  python3 validate_kdwarf.py --name GJ820A --theta-mas 1.775 \
      --fbol 3.90e-7 --logg 4.67 --feh -0.33 [--spec path.fits]
"""
import argparse
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_mann import (REPO, RUN_ROOT, START_MODELS, run,
                           iter_convergence, smooth_to_R,
                           OBS_R_BLUE, OBS_R_RED, OBS_R_SPLIT, C_ANG)

SIGMA_SB = 5.670374419e-5
CAL_DIR = os.environ.get(
    "CALIB_DIR",
    os.path.expanduser("~/sps/SPECTRA/Mann/Calibrators2013"))

WINDOWS = [("5500-5700", 5500, 5700), ("5900-6100", 5900, 6100),
           ("CaH 6800-7000", 6800, 7000), ("gam 7053-7270", 7053, 7270),
           ("7450-7650", 7450, 7650), ("8500-9000", 8500, 9000),
           ("NIR 15000-18000", 15000, 18000), ("9800-11000", 9800, 11000)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--theta-mas", type=float, required=True)
    ap.add_argument("--fbol", type=float, required=True,
                    help="bolometric flux at Earth [erg/s/cm^2]")
    ap.add_argument("--logg", type=float, required=True)
    ap.add_argument("--feh", type=float, default=0.0)
    ap.add_argument("--vturb", type=float, default=1.0)
    ap.add_argument("--numit", type=int, default=30)
    ap.add_argument("--resolu", type=float, default=50000.0)
    ap.add_argument("--wlbeg", type=float, default=300.0)
    ap.add_argument("--wlend", type=float, default=2500.0)
    ap.add_argument("--solar", default="berg25")
    ap.add_argument("--spec", default=None,
                    help="calibrator FITS (default CAL_DIR/<name>.fits)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    theta_rad = args.theta_mas * 1.0e-3 / 206264.806
    teff = (4.0 * args.fbol / (SIGMA_SB * theta_rad ** 2)) ** 0.25
    dilute = (theta_rad / 2.0) ** 2
    print(f"{args.name}: theta_LD={args.theta_mas:.3f} mas  "
          f"Fbol={args.fbol:.3e}  ->  Teff={teff:.0f} K   "
          f"logg={args.logg:.2f}  [Fe/H]={args.feh:+.2f}")

    spec_fits = args.spec or os.path.join(CAL_DIR, args.name + ".fits")
    d = fits.open(spec_fits)[0].data
    wobs, fobs, eobs = d[0] * 1.0e4, d[1], d[2]
    ok = np.isfinite(fobs) & (fobs > 0)
    wobs, fobs = wobs[ok], fobs[ok]

    rundir = os.path.join(RUN_ROOT, args.name + "_int")
    os.makedirs(rundir, exist_ok=True)
    base = args.name
    atm = os.path.join(rundir, base + ".atm")
    start_teff, start_atm = min(START_MODELS,
                                key=lambda m: abs(np.log(m[0] / teff)))
    if args.force or not os.path.isfile(atm):
        print(f"running atlas12 from {os.path.basename(start_atm)} ...")
        cmd = [os.path.join(REPO, "bin", "atlas12c.exe"), start_atm, base,
               f"teff={teff:.0f}", f"logg={args.logg:.2f}",
               f"zscale={10.0 ** args.feh:.4f}", f"vturb={args.vturb}",
               f"numit={args.numit}"]
        if args.solar:
            cmd.append(f"solar={args.solar}")
        run(cmd, rundir, os.path.join(rundir, base + ".stdout"))
    rms, emax = iter_convergence(os.path.join(rundir, base + ".iter"))
    print(f"atlas12 convergence: RMS={rms:.2f}%  max={emax:.2f}%")

    spec_file = os.path.join(rundir, base + ".spec")
    if args.force or not os.path.isfile(spec_file):
        print(f"running synthe {args.wlbeg:.0f}-{args.wlend:.0f} nm ...")
        run([os.path.join(REPO, "bin", "synthe.exe"), atm,
             f"wlbeg={args.wlbeg}", f"wlend={args.wlend}",
             f"resolu={args.resolu:.0f}"],
            rundir, os.path.join(rundir, base + ".synthe.stdout"))

    wmod, hnu, hcont = np.loadtxt(spec_file, unpack=True)
    flam = 4.0 * np.pi * hnu * C_ANG / wmod ** 2 * dilute
    clam = 4.0 * np.pi * hcont * C_ANG / wmod ** 2 * dilute
    fb = smooth_to_R(flam, args.resolu, OBS_R_BLUE)
    fr = smooth_to_R(flam, args.resolu, OBS_R_RED)
    fsm = np.where(wmod < OBS_R_SPLIT, fb, fr)
    cb = smooth_to_R(clam, args.resolu, OBS_R_BLUE)
    cr = smooth_to_R(clam, args.resolu, OBS_R_RED)
    csm = np.where(wmod < OBS_R_SPLIT, cb, cr)

    sel = (wobs >= wmod[0]) & (wobs <= wmod[-1])
    wo, do = wobs[sel], fobs[sel]
    fm = np.interp(wo, wmod, fsm)
    cm = np.interp(wo, wmod, csm)

    print(f"integral over obs range: model/data = "
          f"{np.trapz(fm, wo) / np.trapz(do, wo):.3f}")
    print("window            model/data   data/cont")
    for lab, a, b in WINDOWS:
        s = (wo > a) & (wo < b)
        if s.sum() < 5:
            continue
        print(f"  {lab:16s} {np.median(fm[s] / do[s]):8.3f}   "
              f"{np.median(do[s] / cm[s]):8.3f}")


if __name__ == "__main__":
    main()
