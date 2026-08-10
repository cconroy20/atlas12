#!/usr/bin/env python3
"""External validation of ATLAS12+SYNTHE models against the Mann M/K-dwarf
spectra (flux-calibrated 0.3-2.4 um SNIFS+SpeX).

Handles BOTH spectral libraries transparently (dispatch in mann_lib):
  - Mann+15 M dwarfs (private 2016 library): catalog params from
    M_params.fits (Teff/R/M/[Fe/H]; logg from M and R).
  - Mann+13 K-dwarf calibrators (public Calibrators2013): model-free
    curated params (theta_LD + Fbol -> Teff; dilution = (theta/2)^2).

Stars may be named by GJ designation (GJ887, "Gl 15A") or Mann PM_I name;
run directories use the GJ name where one exists: workdir/mann/<name>[_tag].

For a chosen star this script:
  1. resolves parameters (mann_lib.resolve),
  2. runs atlas12c.exe from the nearest starting model,
  3. runs synthe.exe (default full range 350-2500 nm at R=300000),
  4. converts the model to f_lambda at Earth, convolves with the measured
     Mann LSF, prints absolute metrics, and writes a quick 2-panel plot.

The comparison is ABSOLUTE (no renormalization).  For the standing
3-panel point-comparison figure (incl. PHOENIX NewEra) use
mann_compare_plot.py; for the ladder batch use mann_pointcomp.py.

Usage:
  python3 validate_mann.py --list                # show the sample
  python3 validate_mann.py --star GJ887          # run one star
  python3 validate_mann.py --star GJ887 --force  # rerun atlas12+synthe
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml
from mann_lib import (REPO, MANN_DIR, RUN_ROOT, C3K_DIR, C3K_R,
                      smooth_mann, smooth_to_R, run, iter_convergence)


def list_sample():
    t = ml.load_mann()
    _, i2gj = ml.gj_maps()
    order = np.argsort(t["teff"])[::-1]
    print(f"{'name':<10}{'PM name':<18}{'SpT':<8}{'Teff':>6}{'logg':>7}"
          f"{'[Fe/H]':>8}{'R/Rsun':>8}{'d(pc)':>8}  spectrum")
    for k in ml.KDWARF_PARAMS:
        s = ml.resolve(k)
        print(f"{s.name:<10}{'-':<18}{s.spt:<8}{s.teff:>6.0f}{s.logg:>7.2f}"
              f"{s.feh:>8.2f}{'-':>8}{'-':>8}  calib2013")
    for i in order:
        has = os.path.isfile(os.path.join(MANN_DIR, t["name"][i] + ".ascii"))
        gj = i2gj.get(i, "-")
        print(f"{gj:<10}{t['name'][i]:<18}{t['spt'][i]:<8}{t['teff'][i]:>6.0f}"
              f"{t['logg'][i]:>7.2f}{t['feh'][i]:>8.2f}{t['radius'][i]:>8.3f}"
              f"{t['dist'][i]:>8.2f}  {'yes' if has else 'NO'}")


def load_c3k(teff):
    """C3K v2.3 grid spectrum closest in Teff (solar, logg 4.5 only)."""
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


def run_star(s, args):
    """Run atlas12 + synthe for a resolved Star; return (rundir, metrics)."""
    zmet = s.mh if args.use_mh else s.feh
    if not np.isfinite([s.teff, s.logg, zmet]).all():
        sys.exit(f"error: {s.name} has missing parameters")

    mode = ("model-free (theta+Fbol)" if s.library == "calib2013"
            or s.theta_mas else "catalog table")
    print(f"{s.name} ({s.spt}, {s.library}, {mode}):  Teff={s.teff:.0f} K  "
          f"logg={s.logg:.2f}  [Fe/H]={s.feh:+.2f}  [M/H]={s.mh:+.2f}")

    rundir = ml.rundir_for(s, args.tag)
    os.makedirs(rundir, exist_ok=True)
    base = s.name
    atm = os.path.join(rundir, base + ".atm")

    start_teff, start_atm = ml.start_model(s.teff)
    if args.force or not os.path.isfile(atm):
        print(f"running atlas12 from {os.path.basename(start_atm)} "
              f"(start Teff {start_teff:.0f} K) ...")
        cmd = [os.path.join(REPO, "bin", "atlas12c.exe"), start_atm, base,
               f"teff={s.teff:.0f}", f"logg={s.logg:.2f}",
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

    spec_file = os.path.join(rundir, base + ".spec")
    if args.force or not os.path.isfile(spec_file):
        print(f"running synthe {args.wlbeg:.0f}-{args.wlend:.0f} nm at "
              f"R={args.resolu:.0f} ...")
        run([os.path.join(REPO, "bin", "synthe.exe"), atm,
             f"wlbeg={args.wlbeg}", f"wlend={args.wlend}",
             f"resolu={args.resolu:.0f}"],
            rundir, os.path.join(rundir, base + ".synthe.stdout"))
    else:
        print("synthe: using existing .spec (use --force to rerun)")

    # ---------------- model -> f_lambda at Earth ----------------
    wmod, flam, _ = ml.read_spec_file(spec_file, dilute=s.dilute)
    fsm = smooth_mann(wmod, flam, args.resolu, split=s.obs_split)

    # ---------------- observed spectrum ----------------
    wobs, fobs, _ = ml.read_spectrum(s)
    ok = (wobs >= wmod[0]) & (wobs <= wmod[-1])
    fmod_i = np.interp(wobs, wmod, fsm)

    # ---------------- reference C3K grid spectrum (f77 Kurucz) ------------
    c3k = load_c3k(s.teff)
    if c3k is not None:
        c3k_teff, wc3k, hnu_c3k = c3k
        flam_c3k = ml.hnu_to_flam(wc3k, hnu_c3k) * s.dilute
        csm = smooth_mann(wc3k, flam_c3k, C3K_R, split=s.obs_split)
        okc = ok & (wobs >= wc3k[0]) & (wobs <= wc3k[-1])
        fc3k_i = np.interp(wobs, wc3k, csm)
        c3k_label = f"C3K v2.3 f77 (t{c3k_teff:.0f} g4.5 [Fe/H]=0)"
        print(f"reference grid point: {c3k_label}")

    # ---------------- metrics ----------------
    met = ml.band_metrics(wobs[ok], fobs[ok], fmod_i[ok])
    met["rms"], met["emax"] = rms, emax
    met["teff"], met["logg"], met["feh"] = s.teff, s.logg, s.feh
    fbol_obs = np.trapz(fobs[ok], wobs[ok])
    fbol_mod = np.trapz(np.interp(wobs[ok], wmod, flam), wobs[ok])
    print(f"integral over obs range:  data={fbol_obs:.3e}  "
          f"model={fbol_mod:.3e} erg/s/cm^2 "
          f"(model/data = {met['integral']:.3f})")
    print(f"star FBOL = {s.fbol:.3e}")
    for key, lab in [("opt_med", "optical 0.40-0.95 um"),
                     ("nir_med", "NIR     0.95-2.40 um")]:
        if key in met:
            line = f"median model/data, {lab}: {met[key]:.3f}"
            if c3k is not None:
                lo, hi = (4000, 9500) if key == "opt_med" else (9500, 24000)
                mc = okc & (wobs > lo) & (wobs < hi)
                if mc.any():
                    line += f"   (c3k: {np.median(fc3k_i[mc] / fobs[mc]):.3f})"
            print(line)
    with open(os.path.join(rundir, base + "_metrics.json"), "w") as fh:
        json.dump(met, fh, indent=1)

    # ---------------- figure: house 3-panel format ----------------
    if not args.tag:               # tagged runs are experiments; skip figure
        from mann_compare_plot import plot_star
        plot_star(s, model_R=args.resolu)
    return rundir, met


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--star", help="GJ or Mann PM_I name, e.g. GJ887")
    ap.add_argument("--list", action="store_true",
                    help="list the sample and exit")
    ap.add_argument("--vturb", type=float, default=1.0,
                    help="microturbulence km/s")
    ap.add_argument("--numit", type=int, default=30,
                    help="atlas12 iterations")
    ap.add_argument("--resolu", type=float, default=300000.0,
                    help="synthe resolving power before smoothing "
                         "(>=300000: R=50k over-absorbs the molecular haze "
                         "by ~10 percent median in the optical)")
    ap.add_argument("--wlbeg", type=float, default=350.0, help="nm")
    ap.add_argument("--wlend", type=float, default=2500.0, help="nm")
    ap.add_argument("--use-mh", action="store_true",
                    help="scale metals by [M/H] instead of [Fe/H]")
    ap.add_argument("--solar", default="berg25",
                    help="solar reference pattern for atlas12 (ag89, agss09, "
                         "berg25; default berg25). Pass an empty string to "
                         "keep the table carried by the input model.")
    ap.add_argument("--theta-mas", type=float, default=None,
                    help="interferometric limb-darkened angular diameter "
                         "[mas] for a mann15 star: overrides Teff/R/logg "
                         "and the flux dilution (model-free comparison). "
                         "Default tag becomes 'int'.")
    ap.add_argument("--tag", default=None,
                    help="run-directory suffix (workdir/mann/<star>_<tag>)")
    ap.add_argument("--force", action="store_true",
                    help="rerun even if outputs exist")
    args = ap.parse_args()

    if args.list or not args.star:
        list_sample()
        return

    try:
        s = ml.resolve(args.star)
    except KeyError as e:
        sys.exit(f"error: {e} (see --list)")

    if args.theta_mas:
        if s.library != "mann15":
            sys.exit("error: --theta-mas applies to mann15 stars "
                     "(calib2013 params are already interferometric)")
        told = s.teff
        ml.apply_theta(s, args.theta_mas)
        print(f"INTERFEROMETRIC MODE: theta_LD={args.theta_mas:.3f} mas -> "
              f"Teff={s.teff:.0f} K (table {told:.0f})  "
              f"R={s.radius:.3f} Rsun  logg={s.logg:.2f}")
        if args.tag is None:
            args.tag = "int"

    run_star(s, args)


if __name__ == "__main__":
    main()
