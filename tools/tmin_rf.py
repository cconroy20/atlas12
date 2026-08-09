#!/usr/bin/env python3
"""tmin_rf.py -- finite-difference response functions of band/line indices
to binned dT perturbations of a converged ATLAS12 model (SIR-style, but FD).

Perturbs T by +DT (default 50 K) in top-hat bins of log10 tau5000, then
synthesizes three windows -- optical (540-770 nm), Ca II triplet (838-872 nm),
NIR (1500-1800 nm) -- scores each against the Mann spectrum with the same
indices as tmin_fit.py, and writes rf.json: per-bin indices + the Jacobian
d(index)/dT [per K].  A linearity check reruns one strong-response bin at
DT/2.  Analysis (SVD, mode count, linear solvability of the data residual)
is in --analyze, which only reads rf.json.

Usage:
  TMIN_STAR=GJ887 python3 tmin_rf.py --smoke    # baseline + one bin, report
  TMIN_STAR=GJ887 python3 tmin_rf.py --run      # full battery (skips done)
  TMIN_STAR=GJ887 python3 tmin_rf.py --analyze  # Jacobian/SVD/solvability
"""
import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

REPO = "/Users/cconroy/kurucz/atlas12"
STAR = os.environ.get("TMIN_STAR", "GJ887")
BASE = f"{REPO}/workdir/mann/{STAR}/{STAR}.atm"
RF = f"{REPO}/workdir/mann/{STAR}/tmin_rf"
MANN_DIR = os.path.expanduser("~/sps/SPECTRA/Mann")
RSUN, PC, C_ANG = 6.957e10, 3.0856776e18, 2.99792458e18

DT = float(os.environ.get("TMIN_RF_DT", "50"))
# bin edges in log10 tau5000, deep -> shallow; bin i = (EDGES[i+1], EDGES[i]]
# TMIN_RF_EDGES overrides (comma list); results then go to TMIN_RF_JSON
# (default rf_fine.json) with 2-decimal tags so runs never collide with the
# default battery.
_CUSTOM = "TMIN_RF_EDGES" in os.environ
EDGES = ([float(x) for x in os.environ["TMIN_RF_EDGES"].split(",")]
         if _CUSTOM else
         [1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -4.0])
OUTJSON = os.environ.get("TMIN_RF_JSON",
                         "rf_fine.json" if _CUSTOM else "rf.json")
TAGFMT = "+.2f" if _CUSTOM else "+.1f"
LIN_BIN = 5  # (-2.0, -1.5], expected strong TiO response; rerun at DT/2
WINDOWS = {"opt": (540.0, 770.0), "ca": (838.0, 872.0), "nir": (1500.0, 1800.0)}

TIO_BANDS = [(5500, 5700), (5900, 6100), (7053, 7270), (7450, 7650)]
CAH = {"cah2": (6814., 6846.), "cah3": (6960., 6990.), "tio5": (7126., 7135.)}
CAH_DEN = (7042., 7046.)
CA_CORES = {"ca8542": (8538., 8551.), "ca8662": (8658., 8671.)}
CA_SIDE = (8600., 8640.)

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
    b = smooth_to_R(f, 300000., 1000.)
    r = smooth_to_R(f, 300000., 2000.)
    fs = np.where(w < 9500., b, r)
    m = (WOBS > w[0] + 30) & (WOBS < w[-1] - 30)
    return m, np.interp(WOBS[m], w, fs)


def index_win(w, f, num, den):
    mn = (w > num[0]) & (w < num[1])
    md = (w > den[0]) & (w < den[1])
    return float(np.mean(f[mn]) / np.mean(f[md]))


# ---------------- perturbation ----------------
def perturb_bin(in_atm, out_atm, lo, hi, dt):
    """T += dt for layers with lo < log10(tau5000) <= hi; all else untouched."""
    lines = open(in_atm).readlines()
    istart = next(k for k, l in enumerate(lines)
                  if l.startswith('READ DECK')) + 1
    ndep = int(lines[istart - 1].split()[2])
    nhit = 0
    for k in range(ndep):
        p = lines[istart + k].split()
        rhox, T = float(p[0]), float(p[1])
        rest = [float(x) for x in p[2:]]
        lt = np.log10(rest[-1])
        if lo < lt <= hi:
            T += dt
            nhit += 1
        lines[istart + k] = (f'{rhox:25.5E}{T:10.2f}' +
                             ''.join(f'{v:10.3E}' for v in rest) + '\n')
    open(out_atm, 'w').write(''.join(lines))
    return nhit


# ---------------- synthesis + scoring ----------------
def run_synthe(atm, cwd, wlbeg, wlend, log):
    env = dict(os.environ, ATLAS12=REPO)
    with open(log, "w") as fh:
        subprocess.run([f"{REPO}/bin/synthe.exe", atm,
                        f"wlbeg={wlbeg}", f"wlend={wlend}", "resolu=300000"],
                       cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT,
                       check=True)


def one_model(tag, windows=WINDOWS):
    atm = f"{RF}/{tag}.atm"
    for sub, (b, e) in windows.items():
        rd = f"{RF}/{tag}_{sub}"
        spec = f"{rd}/{tag}.spec"
        if os.path.exists(spec):
            continue
        os.makedirs(rd, exist_ok=True)
        run_synthe(os.path.join("..", f"{tag}.atm"), rd, b, e,
                   f"{rd}/synthe.stdout")
    return tag


def score(tag):
    mo, fo_m = model_on_obs(f"{RF}/{tag}_opt/{tag}.spec")
    mc, fc_m = model_on_obs(f"{RF}/{tag}_ca/{tag}.spec")
    mn, fn_m = model_on_obs(f"{RF}/{tag}_nir/{tag}.spec")
    wo, do = WOBS[mo], FOBS[mo]
    wc, dc = WOBS[mc], FOBS[mc]
    wn, dn = WOBS[mn], FOBS[mn]
    out = {}
    for j, (a, b) in enumerate(TIO_BANDS):
        s = (wo > a) & (wo < b)
        out[f"tio{j}"] = float(np.median(fo_m[s] / do[s]))
    for k, num in CAH.items():
        out[k] = index_win(wo, fo_m, num, CAH_DEN)
        out[k + "_data"] = index_win(wo, do, num, CAH_DEN)
    for k, num in CA_CORES.items():
        out[k] = index_win(wc, fc_m, num, CA_SIDE)
        out[k + "_data"] = index_win(wc, dc, num, CA_SIDE)
    out["nir"] = float(np.median(fn_m / dn))
    out["opt_med"] = float(np.median(fo_m / do))
    return out


# names of the Jacobian observables (ratios/indices, in fixed order)
OBS = [f"tio{j}" for j in range(4)] + ["cah2", "cah3", "nir",
                                       "ca8542", "ca8662"]
# per-observable sigmas, matching tmin_fit chi2 conventions (Ca cores looser:
# LTE core systematics)
SIG = dict(**{f"tio{j}": 0.03 for j in range(4)}, cah2=0.03, cah3=0.03,
           nir=0.02, ca8542=0.05, ca8662=0.05)


def bin_tag(ib, dt):
    return (f"b{ib}_{EDGES[ib + 1]:{TAGFMT}}_{EDGES[ib]:{TAGFMT}}_dt{dt:.0f}"
            .replace("+", "p").replace("-", "m").replace(".", ""))


def make_atms(bins, dt):
    tags = []
    for ib in bins:
        tag = bin_tag(ib, dt)
        atm = f"{RF}/{tag}.atm"
        if not os.path.exists(atm):
            n = perturb_bin(BASE, atm, EDGES[ib + 1], EDGES[ib], dt)
            print(f"{tag}: {n} layers perturbed by +{dt:.0f} K")
        tags.append((ib, tag))
    return tags


def run_battery(bins, dt, workers=4):
    os.makedirs(RF, exist_ok=True)
    if not os.path.exists(f"{RF}/base.atm"):
        import shutil
        shutil.copy(BASE, f"{RF}/base.atm")
    jobs = ["base"] + [t for _, t in make_atms(bins, dt)]
    # linearity check at half amplitude (default battery only)
    if dt == DT and LIN_BIN in bins and not _CUSTOM:
        jobs += [t for _, t in make_atms([LIN_BIN], DT / 2)]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one_model, jobs))
    return jobs


def collect(bins, dt):
    out = {"star": STAR, "dt": dt, "edges": EDGES,
           "base": score("base"), "bins": {}}
    for ib in bins:
        tag = bin_tag(ib, dt)
        out["bins"][str(ib)] = dict(score(tag), tag=tag,
                                    lo=EDGES[ib + 1], hi=EDGES[ib])
    lt = bin_tag(LIN_BIN, DT / 2)
    if not _CUSTOM and os.path.exists(f"{RF}/{lt}_opt/{lt}.spec"):
        out["linearity"] = dict(score(lt), tag=lt, dt=DT / 2)
    json.dump(out, open(f"{RF}/{OUTJSON}", "w"), indent=1)
    return out


# ---------------- analysis ----------------
def _jacobian():
    rf = json.load(open(f"{RF}/{OUTJSON}"))
    base, dt = rf["base"], rf["dt"]
    E = rf["edges"]
    ibs = sorted(int(k) for k in rf["bins"])
    centers = [(E[i] + E[i + 1]) / 2 for i in ibs]
    # Jacobian: d(obs)/dT  [per K], rows = OBS, cols = bins
    J = np.array([[(rf["bins"][str(i)][o] - base[o]) / dt for i in ibs]
                  for o in OBS])
    # residual target: what the data want minus what the RE model gives
    targ = {**{f"tio{j}": 1.0 for j in range(4)}, "nir": 1.0}
    for k in ("cah2", "cah3", "ca8542", "ca8662"):
        targ[k] = base[k + "_data"]
    r = np.array([targ[o] - base[o] for o in OBS])
    w = np.array([1.0 / SIG[o] for o in OBS])
    return rf, base, centers, J, r, w, targ


def _solve(Jw, rw, sig_s, free=None):
    """curvature-regularized signed lsq; free = bin indices allowed to move
    (others pinned to dT=0; deep-layers-untouched constraint)."""
    n = Jw.shape[1]
    if free is None:
        free = list(range(n))
    Jf = Jw[:, free]
    L = np.diff(np.eye(len(free)), 2, axis=0)
    if len(free) > 2:
        A = np.vstack([Jf, L / sig_s])
        b = np.concatenate([rw, np.zeros(L.shape[0])])
    else:
        A, b = Jf, rw
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    dT = np.zeros(n)
    dT[free] = x
    return dT


def _solve_mono(Jw, rw):
    from scipy.optimize import nnls
    n = Jw.shape[1]
    M = np.tril(np.ones((n, n))).T
    u, _ = nnls(Jw @ M, rw)
    return M @ u


def perturb_profile(in_atm, out_atm, centers, dT):
    """apply per-layer delta-T interpolated from (bin center, dT) pairs;
    flat extrapolation beyond the end bins."""
    lines = open(in_atm).readlines()
    istart = next(k for k, l in enumerate(lines)
                  if l.startswith('READ DECK')) + 1
    ndep = int(lines[istart - 1].split()[2])
    c = np.asarray(centers)
    o = np.argsort(c)
    for k in range(ndep):
        p = lines[istart + k].split()
        rhox, T = float(p[0]), float(p[1])
        rest = [float(x) for x in p[2:]]
        lt = np.log10(rest[-1])
        T += float(np.interp(lt, c[o], np.asarray(dT)[o]))
        lines[istart + k] = (f'{rhox:25.5E}{T:10.2f}' +
                             ''.join(f'{v:10.3E}' for v in rest) + '\n')
    open(out_atm, 'w').write(''.join(lines))


def verify(sig_s=30.0, workers=3, scales=(1.0,), anchor=None):
    """forward-verify the linear-inversion solutions: synthesize the
    regularized signed profile (optionally damped by the given scale
    factors -- Gauss-Newton line search for the strongly nonlinear cool
    stars; optionally with dT=0 pinned deeper than `anchor`) and, at
    scale 1 unanchored, the monotone ceiling; score everything."""
    rf, base, centers, J, r, w, targ = _jacobian()
    Jw, rw = J * w[:, None], r * w
    free = (None if anchor is None
            else [i for i, c in enumerate(centers) if c < anchor])
    dT0 = _solve(Jw, rw, sig_s, free)
    stem = (f"vsigned{sig_s:.0f}" if anchor is None
            else f"vthin{anchor:+.1f}".replace("+", "p")
            .replace("-", "m").replace(".", ""))
    if anchor is not None and sig_s != 30.0:
        stem += f"s{sig_s:.0f}"  # sig_s=30 keeps legacy tag names
    if _CUSTOM:
        stem = "f" + stem
    sols = {}
    for f in scales:
        tag = stem if f == 1.0 else f"{stem}x{f:.2f}".replace(".", "p")
        sols[tag] = f * dT0
    if 1.0 in scales and anchor is None:
        sols["vmono"] = _solve_mono(Jw, rw)
    for tag, dT in sols.items():
        atm = f"{RF}/{tag}.atm"
        if not os.path.exists(atm):
            perturb_profile(BASE, atm, centers, dT)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one_model, sols))
    for tag, dT in sols.items():
        sc = score(tag)
        chi2_lin = float(np.sum((Jw @ dT - rw) ** 2))
        chi2_act = float(np.sum([((sc[o] - targ[o]) / SIG[o]) ** 2
                                 for o in OBS]))
        print(f"\n== {tag}  dT(deep->shallow): " +
              " ".join(f"{v:+.0f}" for v in dT))
        print(f"chi2 linear-predicted {chi2_lin:.1f}  actual {chi2_act:.1f}"
              f"  (RE baseline {float(np.sum(rw ** 2)):.1f})")
        print("obs      target   linpred  actual")
        for k, o in enumerate(OBS):
            lin = base[o] + (J @ dT)[k]
            print(f"  {o:7s} {targ[o]:7.3f}  {lin:7.3f}  {sc[o]:7.3f}")
        json.dump(dict(sc, dT=list(dT), chi2=chi2_act),
                  open(f"{RF}/{tag}_score.json", "w"), indent=1)


def analyze():
    rf, base, centers, J, r, w, targ = _jacobian()
    Jw, rw = J * w[:, None], r * w
    dt = rf["dt"]

    print(f"== {rf['star']}  (dT probe = {dt:.0f} K) ==")
    print("bin centers (log tau5000):",
          " ".join(f"{c:+.2f}" for c in centers))
    print("\nJacobian d(index)/dT [1e-3 per K]  (rows=obs, cols=deep->shallow)")
    for o, row in zip(OBS, J):
        print(f"  {o:7s} " + " ".join(f"{1e3 * v:+7.3f}" for v in row))
    print("\nresidual (target - RE):")
    print("  " + "  ".join(f"{o}={v:+.3f}" for o, v in zip(OBS, r)))

    U, s, Vt = np.linalg.svd(Jw, full_matrices=False)
    # a mode is "constrained" if a 100 K excursion along it exceeds 1 sigma
    print("\nsingular values (sigmas of signal per 100 K of mode amplitude):",
          " ".join(f"{100 * v:.1f}" for v in s))
    print("modes with 100K excursion > 1 sigma: "
          f"{int(np.sum(100 * s > 1))}")
    for k in range(min(4, len(s))):
        print(f"  mode {k + 1} (s100={100 * s[k]:.1f}) dT shape: " +
              " ".join(f"{v:+.2f}" for v in Vt[k]))

    chi2_0 = float(np.sum(rw ** 2))
    print(f"\nchi2(RE baseline) = {chi2_0:.1f}")

    # (a) signed, second-difference-regularized; sigma_s = curvature scale
    # in K per bin^2 (Pelletier-style; sigma_s -> inf is unregularized)
    n = Jw.shape[1]
    L = np.diff(np.eye(n), 2, axis=0)
    print("\nsigned dT, curvature-regularized:")
    print("sig_s[K/bin2]   chi2   dT(deep->shallow)")
    for sig_s in (10., 30., 100., 300., None):
        A = np.vstack([Jw, L / sig_s]) if sig_s else Jw
        b = np.concatenate([rw, np.zeros(L.shape[0])]) if sig_s else rw
        dT, *_ = np.linalg.lstsq(A, b, rcond=None)
        chi2 = float(np.sum((Jw @ dT - rw) ** 2))
        lab = f"{sig_s:8.0f}" if sig_s else "     inf"
        print(f"{lab}      {chi2:6.1f}  " +
              " ".join(f"{v:+5.0f}" for v in dT))

    # (a2) thin-layers-only: deep bins pinned to 0, signed outward
    # (user constraint 2026-08-08: leave the deep layers untouched)
    print("\nthin-layers-only (dT=0 deeper than anchor), signed, sig_s=30:")
    print("anchor   chi2   dT(deep->shallow)")
    for anchor in (0.0, -0.5, -1.0, -1.5):
        free = [i for i, c in enumerate(centers) if c < anchor]
        dTa = _solve(Jw, rw, 30., free)
        chi2 = float(np.sum((Jw @ dTa - rw) ** 2))
        print(f"{anchor:+5.1f}  {chi2:6.1f}  " +
              " ".join(f"{v:+5.0f}" for v in dTa))

    # (b) monotone-outward-warming ceiling: dT = cumsum(u), u >= 0 outward
    # (the family all previous forward fits lived in)
    from scipy.optimize import nnls
    M = np.tril(np.ones((n, n))).T  # dT_j = sum_{i<=j} u_i, outward cumsum
    u, _ = nnls(Jw @ M, rw)
    dT_mono = M @ u
    chi2_mono = float(np.sum((Jw @ dT_mono - rw) ** 2))
    print(f"\nmonotone-outward-warming ceiling: chi2 = {chi2_mono:.1f}")
    print("  dT: " + " ".join(f"{v:+5.0f}" for v in dT_mono))

    # per-observable residuals after the two ceilings
    dT_free, *_ = np.linalg.lstsq(Jw, rw, rcond=None)
    print("\nper-obs residual/sigma:  RE     signed  monotone")
    for k, o in enumerate(OBS):
        print(f"  {o:7s} {rw[k]:+7.2f} {rw[k] - (Jw @ dT_free)[k]:+7.2f} "
              f"{rw[k] - (Jw @ dT_mono)[k]:+7.2f}")
    if "linearity" in rf:
        lin = rf["linearity"]
        print(f"\nlinearity check, bin {LIN_BIN} at dT={lin['dt']:.0f}:")
        for o in OBS:
            full = rf["bins"][str(LIN_BIN)][o] - base[o]
            half = lin[o] - base[o]
            frac = 2 * half / full if abs(full) > 1e-6 else float("nan")
            print(f"  {o:7s} 2*d(half)/d(full) = {frac:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="baseline + one bin, report indices")
    ap.add_argument("--run", action="store_true", help="full battery")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--verify", action="store_true",
                    help="synthesize + score the inversion solutions")
    ap.add_argument("--sig-s", type=float, default=30.0,
                    help="curvature scale [K/bin^2] for the verified profile")
    ap.add_argument("--scales", type=str, default="1.0",
                    help="comma list of damping factors for the signed "
                         "profile (Gauss-Newton line search)")
    ap.add_argument("--anchor", type=float, default=None,
                    help="pin dT=0 for bins deeper than this log tau5000 "
                         "(thin-layers-only constraint)")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    if args.smoke:
        run_battery([LIN_BIN], DT, args.workers)
        base, pert = score("base"), score(bin_tag(LIN_BIN, DT))
        print(json.dumps({"base": base, "perturbed": pert}, indent=1))
        print("\ndeltas (perturbed - base):")
        for o in OBS:
            print(f"  {o:7s} {pert[o] - base[o]:+.4f}")
        return
    if args.run:
        run_battery(range(len(EDGES) - 1), DT, args.workers)
        out = collect(range(len(EDGES) - 1), DT)
        print("wrote", f"{RF}/rf.json")
        return
    if args.analyze:
        analyze()
    if args.verify:
        verify(args.sig_s, args.workers,
               tuple(float(x) for x in args.scales.split(",")), args.anchor)


if __name__ == "__main__":
    main()
