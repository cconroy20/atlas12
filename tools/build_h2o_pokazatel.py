#!/usr/bin/env python3
"""Build a Kurucz-format H2O line file from ExoMol POKAZATEL super-lines.

Input: ExoMol pre-computed super-line histograms (R=1e6 binning, 0-41200
cm^-1) at a grid of temperatures, ~/kurucz/upgrade/raw_data/h2o_super/
<T>K.super.bz2 (two columns: nu [cm^-1], S [cm/molecule] at T).

Method: the bins are IDENTICAL across temperature files, and for a bin
holding lines of a single lower energy the intensity obeys exactly
    ln[ S(T) * Q(T) / (1 - exp(-c2 nu / T)) ] = ln(C1 * gf) - c2 * E_low / T
so a linear fit of the left side against 1/T over the FIT_TEMPS recovers
an effective (gf, E_low) per bin -- a pseudo-line whose Boltzmann scaling
matches the true bin contents at the fit temperatures (exact for
single-E bins; least-squares-effective for mixtures, with the residual
kept as a diagnostic).  T=2500 K is deliberately excluded from the fit
and used as a blind reconstruction test (--validate).

Conventions: Q(T) = ExoMol POKAZATEL .pf (full, nuclear-spin included);
the fitted gf is therefore HITRAN-full-convention and is divided by
g_ns = 4 for the Kurucz bare-gf convention (matching the spin-free EOS
row in molecules.dat, refit 2026-08-09, and exomol_to_kurucz.py).
C1 = pi e^2 / (m_e c^2) = 8.85282e-13 cm per (cm/molecule).

Output (--write): h2opokazatel.bin in the exact h2ofastfix.bin record
format (read by both mod_mklinelist.read_h2o and the ATLAS12 unit-51
path): direct-access RECL=8, int32 wavelength index iwl with
wl_vac[nm] = exp(iwl * ln(1 + 1/2e6)), and a packed pair of int16s --
low 16 bits E_low [integer cm^-1], high 16 bits round(1000*log10 gf) +
16384 -- both kept positive (isotopologue 1, 1H2-16O; the POKAZATEL
list is that isotopologue only, and the reader multiplies by 0.9976).

Usage:
  python3 build_h2o_pokazatel.py --cache      # parse .super.bz2 -> npz
  python3 build_h2o_pokazatel.py --fit        # per-bin (gf, E) inversion
  python3 build_h2o_pokazatel.py --validate   # blind 2500 K + stats
  python3 build_h2o_pokazatel.py --write      # emit data/h2opokazatel.bin
"""
import argparse
import bz2
import os
import sys

import numpy as np

RAW = os.path.expanduser("~/kurucz/upgrade/raw_data/h2o_super")
CACHE = os.path.join(RAW, "cache")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PF_FILE = os.path.expanduser("~/kurucz/upgrade/raw_data/exomol_pf/H2O.pf")
OUT_BIN = os.path.join(REPO, "data", "h2opokazatel.bin")

FIT_TEMPS = [1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800,
             3000, 3200, 3400, 3600, 4000, 4400, 4500, 5000]
HOLDOUT_T = 2500
C2 = 1.4387769            # cm K
C1 = 8.85282e-13          # cm/molecule per gf (pi e^2 / m_e c^2)
GNS = 4.0                 # 1H2-16O nuclear-spin degeneracy
XISO1 = 0.9976            # reader-applied isotopic abundance (iso=1)

# keep bins contributing to opacity: S(T_keep)/max > CUT for any T_keep
CUT = 1e-12
NU_MIN, NU_MAX = 400.0, 41200.0        # 243 nm .. 25 um


def _npz(T):
    return os.path.join(CACHE, f"{T}K.npz")


def cmd_cache():
    import subprocess
    import warnings
    os.makedirs(CACHE, exist_ok=True)
    nu_ref = None
    for T in FIT_TEMPS + [HOLDOUT_T]:
        if os.path.isfile(_npz(T)):
            continue
        src = os.path.join(RAW, f"{T}K.super.bz2")
        print(f"parsing {src} ...", flush=True)
        raw = subprocess.run(["bzcat", src], capture_output=True,
                             check=True).stdout
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            arr = np.fromstring(raw, dtype=np.float64, sep=" ")
        arr = arr.reshape(-1, 2)
        nu = arr[:, 0].copy()
        S = arr[:, 1].copy()
        if nu_ref is None:
            nu_ref = nu
        else:
            assert len(nu) == len(nu_ref) and \
                np.abs(nu[::100000] - nu_ref[::100000]).max() < 1e-9, \
                f"{T}K grid differs"
        np.savez(_npz(T), nu=nu if T == FIT_TEMPS[0] else np.array([]),
                 S=S.astype(np.float64))
        print(f"  {len(S)} bins")
    print("cache complete")


def _load_grid():
    nu = np.load(_npz(FIT_TEMPS[0]))["nu"]
    S = np.empty((len(FIT_TEMPS), len(nu)))
    for i, T in enumerate(FIT_TEMPS):
        S[i] = np.load(_npz(T))["S"]
    return nu, S


def _qpf():
    d = np.loadtxt(PF_FILE)
    return lambda T: np.interp(T, d[:, 0], d[:, 1])


# E_low ladder for the NNLS decomposition [cm^-1]: dense at low E where
# the Boltzmann factor discriminates best at 1000-5000 K; capped below
# the int16 storage limit.  Each bin's S(T) is decomposed as a
# non-negative sum over these bands -- the multi-T generalization that
# fixes the single-(gf,E) convexity bias (blind 2500 K median ratio 1.45).
E_LADDER = np.array([0., 700., 1600., 2800., 4300., 6200., 8600., 11600.,
                     15300., 19800., 25200., 31000.])


def cmd_fit_nnls(n_iter=3000, chunk=1_000_000):
    """Batched non-negative decomposition: per bin solve
    min ||W (A x - y)||^2, x >= 0, via normal equations M x = b with
    multiplicative updates x <- x * b / (M x)  (all quantities positive,
    so x stays nonnegative; NMF-style convergence).  Fully vectorized --
    no per-bin Python."""
    nu, S = _load_grid()
    Q = _qpf()
    Tv = np.array(FIT_TEMPS, float)

    # data-side opacity-relevance mask (independent of any fit)
    keep = np.zeros(len(nu), bool)
    for i, T in enumerate(Tv):
        keep |= S[i] > S[i].max() * CUT
    keep &= (nu >= NU_MIN) & (nu <= NU_MAX)
    nu_k = nu[keep]
    nk = keep.sum()
    print(f"NNLS decomposition: {nk} bins x {len(Tv)} temps "
          f"x {len(E_LADDER)} E bands", flush=True)

    # y_j = S_j Q_j / stim_j = C1 * sum_m gf_m exp(-c2 E_m / T_j)
    y = np.empty((nk, len(Tv)))
    for j, T in enumerate(Tv):
        stim = -np.expm1(-C2 * nu_k / T)
        y[:, j] = S[j][keep] * Q(T) / np.maximum(stim, 1e-300)
    A0 = np.exp(-C2 * E_LADDER[None, :] / Tv[:, None])   # (10, 8)

    gf = np.zeros((nk, len(E_LADDER)), np.float64)
    resid = np.zeros(nk, np.float32)
    eps = 1e-300
    for lo in range(0, nk, chunk):
        hi = min(lo + chunk, nk)
        yc = y[lo:hi]                                    # (n, 10)
        ymax = np.maximum(yc.max(axis=1, keepdims=True), eps)
        w2 = 1.0 / np.maximum(yc, ymax * 1e-6) ** 2      # relative weights^2
        # normal equations per bin
        M = np.einsum('nj,jk,jl->nkl', w2, A0, A0)       # (n, 8, 8)
        b = np.einsum('nj,jk->nk', w2 * yc, A0)          # (n, 8)
        x = b / np.maximum(np.einsum('nkk->nk', M), eps)  # positive start
        for it in range(n_iter):
            Mx = np.einsum('nkl,nl->nk', M, x)
            x *= b / np.maximum(Mx, eps)
        gf[lo:hi] = x / C1
        pred = x @ A0.T                                  # (n, 10)
        resid[lo:hi] = np.sqrt(
            (w2 * (pred - yc) ** 2).sum(axis=1) / len(Tv))
        print(f"  {hi}/{nk}", flush=True)
    # refinement: exact NNLS on the bins the multiplicative updates left
    # with the worst weighted residuals (slow-converging strong bins)
    bad = np.flatnonzero(resid > 0.02)
    print(f"refining {len(bad)} bins with exact NNLS ...", flush=True)
    from scipy.optimize import nnls as _nnls
    nref = 0
    for i in bad:
        yi = y[i]
        w = 1.0 / np.maximum(yi, yi.max() * 1e-6)
        try:
            x, rn = _nnls(A0 * w[:, None], yi * w)
        except Exception:
            continue
        rn /= np.sqrt(len(Tv))
        if rn < resid[i]:
            gf[i] = x / C1
            resid[i] = rn
            nref += 1
    print(f"  improved {nref}/{len(bad)}")
    np.savez_compressed(os.path.join(CACHE, "fit_nnls.npz"),
                        nu=nu_k, gf_full=gf, resid=resid,
                        e_ladder=E_LADDER, idx=np.flatnonzero(keep))
    nz = (gf > gf.max(axis=1, keepdims=True) * 1e-6).sum(axis=1)
    print(f"done; mean active bands/bin {nz.mean():.2f}")


def cmd_fit():
    nu, S = _load_grid()
    Q = _qpf()
    Tv = np.array(FIT_TEMPS, float)
    print(f"{len(nu)} bins x {len(Tv)} temperatures")

    # y = ln[S Q / (1-exp(-c2 nu/T))] ; guard S=0
    y = np.empty_like(S)
    for i, T in enumerate(Tv):
        stim = -np.expm1(-C2 * nu / T)
        stim[stim <= 0] = 1e-300
        y[i] = np.log(np.maximum(S[i], 1e-300) * Q(T) / stim)

    # closed-form 2-parameter LSQ per bin: y = b - a * x, x = 1/T
    x = 1.0 / Tv
    n = len(Tv)
    sx, sxx = x.sum(), (x * x).sum()
    sy = y.sum(axis=0)
    sxy = (x[:, None] * y).sum(axis=0)
    denom = n * sxx - sx * sx
    a = -(n * sxy - sx * sy) / denom          # slope -> c2 * E_low
    b = (sy * sxx - sx * sxy) / denom         # intercept -> ln(C1 gf)
    E = a / C2                                # cm^-1
    gf_full = np.exp(b) / C1
    resid = y - (b[None, :] - a[None, :] * x[:, None])
    rms = np.sqrt((resid ** 2).mean(axis=0))

    np.savez(os.path.join(CACHE, "fit.npz"),
             nu=nu, E=E, gf_full=gf_full, rms=rms)
    print(f"fit written; E range {E.min():.0f}..{E.max():.0f} cm-1; "
          f"median ln-rms {np.median(rms):.3f}")


def _strength(nu, gf_full, E, T, Q):
    stim = -np.expm1(-C2 * nu / T)
    return C1 * gf_full * np.exp(-C2 * E / T) * stim / Q(T)


def _keep_mask(nu, gf_full, E, Q):
    """Opacity-relevant bins: S/S_max > CUT at any of 1500/2500/4000 K,
    within the wavenumber window."""
    keep = np.zeros(len(nu), bool)
    for T in (1500.0, 2500.0, 4000.0):
        S = _strength(nu, gf_full, E, T, Q)
        keep |= S > S.max() * CUT
    keep &= (nu >= NU_MIN) & (nu <= NU_MAX) & (gf_full > 0) & (E >= 0)
    return keep


def _recon_nnls(nu, gf, e_ladder, T, Q):
    stim = -np.expm1(-C2 * nu / T)
    boltz = np.exp(-C2 * e_ladder[None, :] / T)
    return C1 * (gf * boltz).sum(axis=1) * stim / Q(T)


def cmd_validate():
    p = os.path.join(CACHE, "fit_nnls.npz")
    T = float(HOLDOUT_T)
    Q = _qpf()
    if os.path.isfile(p):
        d = np.load(p)
        nu, gf, el, idx = d["nu"], d["gf_full"], d["e_ladder"], d["idx"]
        Sh = np.load(_npz(HOLDOUT_T))["S"][idx]
        Sr = _recon_nnls(nu, gf, el, T, Q)
        label = f"NNLS ({len(el)} E bands)"
    else:
        d = np.load(os.path.join(CACHE, "fit.npz"))
        nu, E, gf_full = d["nu"], d["E"], d["gf_full"]
        keep = _keep_mask(nu, gf_full, E, Q)
        nu, Sh = nu[keep], np.load(_npz(HOLDOUT_T))["S"][keep]
        Sr = _strength(nu, gf_full[keep], E[keep], T, Q)
        label = "single-line"
    print(f"{label}: {len(nu)} bins; blind test at {HOLDOUT_T} K:")
    print(f"  total S: recon/actual = {Sr.sum() / Sh.sum():.4f}")
    r = Sr / np.maximum(Sh, 1e-300)
    w = Sh / Sh.sum()
    print(f"  S-weighted mean ratio  = {(r * w).sum():.4f}")
    q = np.percentile(r, [2, 16, 50, 84, 98])
    print(f"  ratio percentiles 2/16/50/84/98: "
          + " ".join(f"{v:.3f}" for v in q))
    # S-weighted percentiles (what the opacity actually feels)
    o = np.argsort(r)
    cw = np.cumsum(w[o])
    qw = [r[o][np.searchsorted(cw, f)] for f in (0.02, 0.16, 0.5, 0.84, 0.98)]
    print(f"  S-WEIGHTED percentiles:          "
          + " ".join(f"{v:.3f}" for v in qw))
    for lo, hi, lab in [(3900, 7150, "1.4-2.6 um"), (7150, 10500, "0.95-1.4"),
                        (10500, 15000, "0.67-0.95"), (15000, 41200, "<0.67um")]:
        m = (nu > lo) & (nu < hi)
        if m.any():
            print(f"  band {lab:<10}: recon/actual = "
                  f"{Sr[m].sum() / Sh[m].sum():.4f}   "
                  f"(share of total S: {Sh[m].sum()/Sh.sum():.2e})")


def cmd_write():
    d = np.load(os.path.join(CACHE, "fit_nnls.npz"))
    nu_b, gf_mat, el = d["nu"], d["gf_full"], d["e_ladder"]
    # expand components: one pseudo-line per (bin, active E band);
    # multiplicative updates only decay dead bands, so apply a relative
    # floor (1e-4 of the bin's strongest band = 0.01% strength)
    bi, mi = np.nonzero(gf_mat > gf_mat.max(axis=1, keepdims=True) * 1e-4)
    nu = nu_b[bi]
    E = el[mi]
    gf_full = gf_mat[bi, mi]
    print(f"{len(nu_b)} bins -> {len(nu)} pseudo-lines")

    gf_bare = gf_full / GNS / XISO1     # reader multiplies xiso back on
    wl_nm = 1.0e7 / nu                  # vacuum nm
    ratiolog = np.log(1.0 + 1.0 / 2.0e6)
    iwl = np.round(np.log(wl_nm) / ratiolog).astype(np.int64)

    igf = np.round(1000.0 * np.log10(gf_bare)).astype(np.int64) + 16384
    ielo = np.round(E).astype(np.int64)
    ok = (igf >= 1) & (igf <= 32767)
    ielo = np.clip(ielo, 1, 32767)      # int16-positive; E>32767 harmless
    iwl, igf, ielo = iwl[ok], igf[ok], ielo[ok]
    print(f"writing {ok.sum()} lines ({(~ok).sum()} below gf floor)")

    order = np.argsort(iwl, kind="stable")
    iwl, igf, ielo = iwl[order], igf[order], ielo[order]

    rec = np.empty((len(iwl), 2), dtype="<i4")
    rec[:, 0] = iwl
    rec[:, 1] = (igf.astype(np.int64) << 16) | ielo
    rec.tofile(OUT_BIN)
    print(f"wrote {OUT_BIN}: {len(iwl)} records, "
          f"{os.path.getsize(OUT_BIN)/1e6:.1f} MB, "
          f"{1e7/np.exp(iwl[-1]*np.log(1+0.5e-6)*1):.0f}..")
    print(f"wl coverage {np.exp(iwl[0]*ratiolog):.1f} - "
          f"{np.exp(iwl[-1]*ratiolog):.1f} nm (vacuum)")


# ===================================================================
# RAW-TRANS ROUTE (production): bin the full 5.7e9-transition POKAZATEL
# list directly onto (iwl on the R=2e6 storage grid) x (E_low bands),
# accumulating exact gf sums and gf-weighted mean E_low per cell.  This
# removes the NNLS approximation floor entirely: temperature scaling is
# exact up to the within-band E spread, which enters only at second
# order once the first moment (mean E) is stored.
RATIOLOG = np.log(1.0 + 1.0 / 2.0e6)
E_EDGES = np.array([0., 300., 700., 1200., 1800., 2600., 3600., 4800.,
                    6200., 8000., 10200., 13000., 16500., 21000., 26500.,
                    33000., 41000., 1.0e9])
TRANS_URL = ("https://exomol.com/db/H2O/1H2-16O/POKAZATEL/"
             "1H2-16O__POKAZATEL__%05d-%05d.trans.bz2")
NCHUNK, CHUNK_CM = 413, 100
STATES_BZ2 = os.path.join(RAW, "H2O.states.bz2")
CKPT = os.path.join(CACHE, "rawbin")
GF_CONST = 1.4992e-16          # gf = GF_CONST * g_u * A * lambda_A^2


def _iwl_bounds():
    lo = int(np.round(np.log(1.0e7 / NU_MAX) / RATIOLOG))
    hi = int(np.round(np.log(1.0e7 / NU_MIN) / RATIOLOG))
    return lo, hi - lo + 1


def _load_states():
    import subprocess
    import warnings
    raw = subprocess.run(["bzcat", STATES_BZ2], capture_output=True,
                         check=True).stdout
    E = {}
    ids, Ev, gv = [], [], []
    for line in raw.splitlines():
        p = line.split(None, 3)
        ids.append(int(p[0])); Ev.append(float(p[1])); gv.append(float(p[2]))
    n = max(ids)
    Ea = np.zeros(n + 1); ga = np.zeros(n + 1)
    Ea[np.array(ids)] = Ev
    ga[np.array(ids)] = gv
    return Ea, ga


def cmd_raw():
    import subprocess
    import warnings
    os.makedirs(CKPT, exist_ok=True)
    done_file = os.path.join(CKPT, "done.txt")
    done = set()
    if os.path.isfile(done_file):
        done = set(open(done_file).read().split())
    iwl0, niwl = _iwl_bounds()
    nb = len(E_EDGES) - 1
    sgf_p, sge_p = os.path.join(CKPT, "sgf.npy"), os.path.join(CKPT, "sge.npy")
    if os.path.isfile(sgf_p):
        sgf = np.load(sgf_p); sge = np.load(sge_p)
    else:
        sgf = np.zeros(nb * niwl); sge = np.zeros(nb * niwl)
    print(f"states ...", flush=True)
    Ea, ga = _load_states()
    print(f"{len(Ea)-1} states; grid {niwl} iwl x {nb} E bands", flush=True)

    for ic in range(NCHUNK):
        a, b = ic * CHUNK_CM, (ic + 1) * CHUNK_CM
        tag = f"{a:05d}-{b:05d}"
        if tag in done:
            continue
        fpath = os.path.join(RAW, f"trans_{tag}.bz2")
        if not os.path.isfile(fpath):
            url = TRANS_URL % (a, b)
            r = subprocess.run(["curl", "-s", "-f", "-o", fpath, url])
            if r.returncode != 0:
                print(f"  {tag}: download FAILED, skipping", flush=True)
                continue
        raw = subprocess.run(["bzcat", fpath], capture_output=True,
                             check=True).stdout
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            first = raw.split(b"\n", 1)[0]
            ncol = len(first.split())
            arr = np.fromstring(raw, dtype=np.float64, sep=" ")
        del raw
        arr = arr.reshape(-1, ncol)
        u = arr[:, 0].astype(np.int64)
        l = arr[:, 1].astype(np.int64)
        A = arr[:, 2]
        del arr
        El = Ea[l]
        nu = Ea[u] - El
        ok = (nu >= NU_MIN) & (nu <= NU_MAX) & (A > 0)
        if ok.any():
            nu, El, A, u = nu[ok], El[ok], A[ok], u[ok]
            gf = GF_CONST * ga[u] * A * (1.0e8 / nu) ** 2
            iwl = np.round(np.log(1.0e7 / nu) / RATIOLOG).astype(np.int64) \
                - iwl0
            band = np.searchsorted(E_EDGES, El, "right") - 1
            idx = band * niwl + iwl
            sgf += np.bincount(idx, weights=gf, minlength=nb * niwl)
            sge += np.bincount(idx, weights=gf * El, minlength=nb * niwl)
        os.remove(fpath)
        done.add(tag)
        with open(done_file, "a") as fh:
            fh.write(tag + "\n")
        if ic % 10 == 0 or ic == NCHUNK - 1:
            np.save(sgf_p, sgf); np.save(sge_p, sge)
        print(f"  {tag}: {ok.sum()} lines binned "
              f"(running total gf {sgf.sum():.4e})", flush=True)
    np.save(sgf_p, sgf); np.save(sge_p, sge)
    print("raw binning complete")


def _load_raw_lines():
    iwl0, niwl = _iwl_bounds()
    nb = len(E_EDGES) - 1
    sgf = np.load(os.path.join(CKPT, "sgf.npy"))
    sge = np.load(os.path.join(CKPT, "sge.npy"))
    nz = np.flatnonzero(sgf > 0)
    gf_full = sgf[nz]
    E = sge[nz] / gf_full                     # gf-weighted mean E_low
    iwl = (nz % niwl) + iwl0
    nu = 1.0e7 / np.exp(iwl * RATIOLOG)
    return iwl, nu, gf_full, E


def cmd_validate_raw():
    """Reconstruct S(T) from the raw-binned pseudo-lines and compare to
    the ExoMol super-line histograms on coarse 10 cm^-1 bins."""
    Q = _qpf()
    iwl, nu, gf_full, E = _load_raw_lines()
    print(f"{len(nu)} raw pseudo-lines")
    nu_e = np.load(_npz(FIT_TEMPS[0]))["nu"]
    edges = np.arange(NU_MIN, NU_MAX + 10.0, 10.0)
    for T in (1400, 2500, 4000):
        if not os.path.isfile(_npz(T)):
            continue
        Se = np.load(_npz(T))["S"]
        Sr = _strength(nu, gf_full, E, float(T), Q)
        hr, _ = np.histogram(nu, bins=edges, weights=Sr)
        he, _ = np.histogram(nu_e, bins=edges, weights=Se)
        m = he > he.max() * 1e-10
        r = hr[m] / he[m]
        q = np.percentile(r, [2, 16, 50, 84, 98])
        print(f"T={T}: total recon/actual = {hr.sum()/he.sum():.4f}   "
              f"10cm-1-bin pct 2/16/50/84/98: "
              + " ".join(f"{v:.3f}" for v in q))
        for lo, hi, lab in [(3900, 7150, "1.4-2.6um"),
                            (7150, 10500, "0.95-1.4"),
                            (10500, 15000, "0.67-0.95")]:
            mm = (nu > lo) & (nu < hi)
            me = (nu_e > lo) & (nu_e < hi)
            print(f"    band {lab:<9}: "
                  f"{Sr[mm].sum()/Se[me].sum():.4f}")


def cmd_write_raw():
    Q = _qpf()
    iwl, nu, gf_full, E = _load_raw_lines()
    keep = np.zeros(len(nu), bool)
    for T in (1500.0, 2500.0, 4000.0):
        S = _strength(nu, gf_full, E, T, Q)
        keep |= S > S.max() * CUT
    iwl, gf_full, E = iwl[keep], gf_full[keep], E[keep]
    print(f"{keep.sum()}/{len(keep)} pseudo-lines pass the strength cut")

    gf_bare = gf_full / GNS / XISO1
    igf = np.round(1000.0 * np.log10(gf_bare)).astype(np.int64) + 16384
    ielo = np.clip(np.round(E).astype(np.int64), 1, 32767)
    ok = (igf >= 1) & (igf <= 32767)
    iwl, igf, ielo = iwl[ok], igf[ok], ielo[ok]
    order = np.argsort(iwl, kind="stable")
    iwl, igf, ielo = iwl[order], igf[order], ielo[order]
    rec = np.empty((len(iwl), 2), dtype="<i4")
    rec[:, 0] = iwl
    rec[:, 1] = (igf.astype(np.int64) << 16) | ielo
    rec.tofile(OUT_BIN)
    print(f"wrote {OUT_BIN}: {len(iwl)} records, "
          f"{os.path.getsize(OUT_BIN)/1e6:.1f} MB; wl "
          f"{np.exp((iwl[0]) * RATIOLOG):.1f}-"
          f"{np.exp((iwl[-1]) * RATIOLOG):.1f} nm vacuum")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", action="store_true")
    ap.add_argument("--fit", action="store_true",
                    help="single-(gf,E) fit (superseded diagnostic)")
    ap.add_argument("--fit-nnls", action="store_true",
                    help="E-ladder NNLS decomposition (production)")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--write", action="store_true",
                    help="write .bin from the NNLS fit (superseded)")
    ap.add_argument("--raw", action="store_true",
                    help="download+bin the full raw trans set (production)")
    ap.add_argument("--validate-raw", action="store_true")
    ap.add_argument("--write-raw", action="store_true",
                    help="write .bin from the raw binning (production)")
    args = ap.parse_args()
    if args.raw:
        cmd_raw()
    if args.validate_raw:
        cmd_validate_raw()
    if args.write_raw:
        cmd_write_raw()
    if args.cache:
        cmd_cache()
    if args.fit:
        cmd_fit()
    if args.fit_nnls:
        cmd_fit_nnls()
    if args.validate:
        cmd_validate()
    if args.write:
        cmd_write()
    if not any(vars(args).values()):
        ap.print_help()


if __name__ == "__main__":
    main()
