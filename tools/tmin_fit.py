#!/usr/bin/env python3
"""GJ887 T-min (tau_min, slope) fit driver.

Perturbs the converged GJ887 structure with tools/tmin_perturb.perturb over a
(logtau_min, slope) grid, synthesizes optical (540-770) and NIR (1500-1800)
windows, scores each against the Mann spectrum, and writes scores.json.

Usage:
  python3 tmin_fit.py --smoke          # single point (-1.7, 125), report score
  python3 tmin_fit.py --grid           # full 5x5 grid (skips existing .spec)
"""
import argparse
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, "/Users/cconroy/kurucz/atlas12/tools")
from tmin_perturb import perturb

import os as _os
REPO = "/Users/cconroy/kurucz/atlas12"
STAR = _os.environ.get("TMIN_STAR", "GJ887")
BASE = f"{REPO}/workdir/mann/{STAR}/{STAR}.atm"
FIT = f"{REPO}/workdir/mann/{STAR}/tmin_fit"
MANN_DIR = os.path.expanduser("~/sps/SPECTRA/Mann")
RSUN, PC, C_ANG = 6.957e10, 3.0856776e18, 2.99792458e18

TM_GRID = [float(x) for x in _os.environ.get("TMIN_TM", "-1.3,-1.5,-1.7,-1.9,-2.1").split(",")]
S_GRID = [float(x) for x in _os.environ.get("TMIN_S", "50,100,150,200,300").split(",")]
# VAL-type saturating grid
TM3 = [-1.3, -1.5, -1.7, -1.9]
S3 = [150.0, 300.0]
DT3 = [150.0, 250.0, 350.0]
# column-mass-coordinate grid (anchors in log10 m [g/cm2], slopes K/dex(m))
M_GRID = [0.3, 0.15, 0.0, -0.15, -0.3]
SM_GRID = [100.0, 200.0, 300.0, 450.0, 600.0]

# ---------------- observed spectrum (surface) ----------------
d = fits.open(f"{MANN_DIR}/M_params.fits")[1].data[0]
names = np.array([n.strip() for n in d["NAME"]])
i = int(np.where(names == STAR)[0][0])
dilute = (d["RADIUS"][i] * RSUN / (d["DISTANCE"][i] * PC)) ** 2
wobs_um, fobs, _ = np.loadtxt(f"{MANN_DIR}/{STAR}.ascii", unpack=True)
wobs = wobs_um * 1e4
ok = fobs > 0
WOBS, FOBS = wobs[ok], fobs[ok] / dilute

def smooth_to_R(f, mR, oR):
    return gaussian_filter1d(f, mR / (oR * 2.3548), mode="nearest")

def model_on_obs(spec_path):
    w, hnu, _ = np.loadtxt(spec_path, unpack=True)
    f = 4 * np.pi * hnu * C_ANG / w ** 2
    b = smooth_to_R(f, 50000., 1000.)
    r = smooth_to_R(f, 50000., 2000.)
    fs = np.where(w < 9500., b, r)
    m = (WOBS > w[0] + 30) & (WOBS < w[-1] - 30)
    return m, np.interp(WOBS[m], w, fs)

def run_synthe(atm, cwd, wlbeg, wlend, log):
    env = dict(os.environ, ATLAS12=REPO)
    with open(log, "w") as fh:
        subprocess.run([f"{REPO}/bin/synthe.exe", atm,
                        f"wlbeg={wlbeg}", f"wlend={wlend}", "resolu=50000"],
                       cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT,
                       check=True)

def one_point(tm, s, dt=math.inf, coord="tau"):
    pre = "cm" if coord == "cmass" else "tm"
    tag = f"{pre}{tm:+.2f}_s{s:.0f}".replace("+", "p").replace("-", "m")
    if math.isfinite(dt):
        tag += f"_dt{dt:.0f}"
    atm = f"{FIT}/{tag}.atm"
    if not os.path.exists(atm):
        perturb(BASE, atm, tm, s, dt, coord)
    for sub, (b, e) in [("opt", (540.0, 770.0)), ("nir", (1500.0, 1800.0))]:
        rd = f"{FIT}/{tag}_{sub}"
        spec = f"{rd}/{tag}.spec"
        if os.path.exists(spec):
            continue
        os.makedirs(rd, exist_ok=True)
        run_synthe(os.path.join("..", f"{tag}.atm"), rd, b, e,
                   f"{rd}/synthe.stdout")
    return tag

def index_obs(w, f, num, den=(7042., 7046.)):
    mn = (w > num[0]) & (w < num[1])
    md = (w > den[0]) & (w < den[1])
    return float(np.mean(f[mn]) / np.mean(f[md]))

TIO_BANDS = [(5500, 5700), (5900, 6100), (7053, 7270), (7450, 7650)]

def score(tag):
    mo, fo_m = model_on_obs(f"{FIT}/{tag}_opt/{tag}.spec")
    mn, fn_m = model_on_obs(f"{FIT}/{tag}_nir/{tag}.spec")
    wo, do = WOBS[mo], FOBS[mo]
    wn, dn = WOBS[mn], FOBS[mn]
    tio = [float(np.median(fo_m[(wo > a) & (wo < b)] / do[(wo > a) & (wo < b)]))
           for a, b in TIO_BANDS]
    out = dict(
        tio_bands=tio,
        tio_med=float(np.median(tio)),
        cah2=index_obs(wo, fo_m, (6814., 6846.)),
        cah3=index_obs(wo, fo_m, (6960., 6990.)),
        cah2_data=index_obs(wo, do, (6814., 6846.)),
        cah3_data=index_obs(wo, do, (6960., 6990.)),
        tio5=index_obs(wo, fo_m, (7126., 7135.)),
        tio5_data=index_obs(wo, do, (7126., 7135.)),
        nir=float(np.median(fn_m / dn)),
        opt_med=float(np.median(fo_m / do)),
    )
    # chi2: TiO windows to 1, CaH2/CaH3 indices to data, NIR to 1
    out["chi2"] = (sum((r - 1.0) ** 2 for r in tio) / len(tio) / 0.03 ** 2
                   + (out["cah2"] - out["cah2_data"]) ** 2 / 0.03 ** 2
                   + (out["cah3"] - out["cah3_data"]) ** 2 / 0.03 ** 2
                   + (out["nir"] - 1.0) ** 2 / 0.02 ** 2)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--grid", action="store_true")
    ap.add_argument("--grid3", action="store_true",
                    help="VAL-type saturating grid (tm, s, dTmax)")
    ap.add_argument("--gridm", action="store_true",
                    help="column-mass-coordinate grid (log m anchor, K/dex m)")
    args = ap.parse_args()
    os.makedirs(FIT, exist_ok=True)

    if args.smoke:
        tag = one_point(-1.7, 125.0)
        s = score(tag)
        print(json.dumps({tag: s}, indent=1))
        return

    if args.grid3:
        pts = [(tm, s, dt) for tm in TM3 for s in S3 for dt in DT3]
        outname = "scores3.json"
    elif args.gridm:
        pts = [(m, s, math.inf, "cmass") for m in M_GRID for s in SM_GRID]
        outname = "scoresm.json"
    else:
        pts = [(tm, s) for tm in TM_GRID for s in S_GRID]
        outname = "scores.json"
    with ThreadPoolExecutor(max_workers=5) as ex:
        tags = list(ex.map(lambda p: one_point(*p), pts))
    scores = {}
    for p, tag in zip(pts, tags):
        sc = score(tag)
        sc["tm"], sc["s"] = p[0], p[1]
        sc["dt"] = p[2] if len(p) > 2 and math.isfinite(p[2]) else None
        sc["coord"] = p[3] if len(p) > 3 else "tau"
        scores[tag] = sc
        print(f"{tag}: tio_med={sc['tio_med']:.3f} cah2={sc['cah2']:.3f} "
              f"cah3={sc['cah3']:.3f} nir={sc['nir']:.3f} chi2={sc['chi2']:.1f}")
    json.dump(scores, open(f"{FIT}/{outname}", "w"), indent=1)
    print("wrote", f"{FIT}/{outname}")

if __name__ == "__main__":
    main()
