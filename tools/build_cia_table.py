#!/usr/bin/env python3
"""Build data/h2collop.dat: collision-induced absorption (CIA) tables.

Replaces the Borysow, Jorgensen & Zheng (1997) H2-H2 + H2-He pair of
7 x 81 tables that the Kurucz lineage shipped, and adds the H2-H and
H-He pairs, which were absent entirely.

SOURCES (all fetched from HITRAN's CIA section, https://hitran.org/cia/,
"Main" folder = the HITRAN-recommended set; see CIA_Readme_2024.pdf
Table 1 for the reference indices quoted below):

  H2-H2  [6]  Abel, Frommhold, Li & Hunt 2011, J.Phys.Chem.A 115, 6805
              T = 200-3000 K, nu = 20-10000 cm^-1
  H2-He  [7]  Abel, Frommhold, Li & Hunt 2012, J.Phys.Chem.A 116, 3068
              T = 200-9900 K, nu = 20-20000 cm^-1
  H2-H  [17]  Gustafsson & Frommhold 2003, A&A 400, 1161
              T = 1000-2500 K (4 points), nu = 100-10000 cm^-1
  H-He  [18]  Gustafsson & Frommhold 2001, ApJ 546, 1168
              T = 1500-10000 K (10 points), nu = 50-11000 cm^-1

Abel et al. supersede the Borysow-lineage tabulations: they are a single
uniform ab initio calculation over the whole (T, nu) domain, where the
older tables "used different approximations and potential surfaces for
different temperature and frequency ranges" (Saumon, Marley, Abel,
Frommhold & Freedman 2012, ApJ 750, 74, sec. 2.1).  Relative to
Borysow, Jorgensen & Fu (2001) the H2-H2 fundamental band is weaker by
24% at 1000 K and 44% at 2000 K -- reproduced by this script's readers
to 3 significant figures, which is how the unit conversion below was
validated.

The H2-H and H-He files on HITRAN are the same Gustafsson & Frommhold
calculations that Turbospectrum ships in COM-v19.1/DATA, but tabulated
at 1 cm^-1 instead of ~50 points, so HITRAN is used for all four pairs
and there is a single provenance.

H2-H2 CONTINUATION.  Abel's H2-H2 box stops at 3000 K and 10^4 cm^-1,
but ATLAS12 runs to 25000 K and this table is used down to 0.5 um.
Outside the box we fall back to the only calculation that covers it,
Borysow, Jorgensen & Fu 2001 (JQSRT 68, 235; Turbospectrum's
CIA.H2-H2.Yi.dat, T = 1000-7000 K, nu = 20-20000 cm^-1), rescaled by
the Abel/BJF01 ratio evaluated along the boundary being crossed.  That
keeps BJF01's *shape* -- the only information available out there --
while carrying Abel's *normalisation* across the seam, so the table is
continuous.  Both continuation regions are physically unimportant:
above 3000 K H2 is dissociating and n(H2)^2 collapses, and above
10^4 cm^-1 CIA is a negligible opacity next to H- and the metals.

UNITS AND CONVENTIONS -- both verified, not assumed:

  * Stored quantity is log10 of the binary absorption coefficient in
    cm^5, i.e. kappa[cm^2/g] = 10**stored * n_1 * n_2 / rho with the
    number densities in cm^-3.  HITRAN tabulates cm^5 molecule^-2
    directly, so no conversion is needed; the Turbospectrum files used
    for the continuation are in cm^-1 amagat^-2 and are divided by
    n_L^2 = (2.6867811e19)^2, which is exactly the 1.38529e-39
    multiplier quoted in their own headers.

  * Stimulated emission IS ALREADY INCLUDED in these coefficients and
    must NOT be applied again by the consuming code.  Turbospectrum's
    detabs.f is explicit -- "Aleksandra Borysow me signale que stim
    deja inclu dans CIA" -- and divides every CIA term by its global
    (1 - exp(-hnu/kT)) factor to cancel it.  Independently, the low
    frequency limit of every one of these tables goes as nu^2 (measured
    slope 2.00 below ~50 cm^-1), which is the signature of an
    absorption coefficient carrying the (1 - exp(-hc nu / kT)) factor;
    without it the limit would be nu^1.

Usage:
    python3 tools/build_cia_table.py --raw <dir> [--out data/h2collop.dat]
    python3 tools/build_cia_table.py --raw <dir> --validate

<dir> holds the five source files, none of which are kept in the repo
(~180 MB).  Fetch them with:

    H=https://hitran.org/data/CIA/main
    for f in H2-H2_2011 H2-He_2011 H2-H_2011 He-H_2011; do
        curl -sSLO $H/$f.cia
    done
    curl -sSLO https://raw.githubusercontent.com/bertrandplez/\\
Turbospectrum2019/master/COM-v19.1/DATA/CIA.H2-H2.Yi.dat

The HITRAN four are the "Main" (recommended) sets; the Turbospectrum
file supplies only the H2-H2 continuation above Abel's range.  Run
--validate after any re-fetch: it re-derives Saumon et al. (2012)'s
published Abel/BJF01 comparison from the raw files, so a changed or
truncated download shows up as a failed check rather than a silently
different table.
"""

import argparse
import os
import sys

import numpy as np

# Loschmidt constant [cm^-3]; TS headers quote 1/n_L^2 = 1.38529e-39
N_LOSCHMIDT = 2.6867811e19

# Common wavenumber grid.  25 cm^-1 keeps the interpolation error under
# 0.4% peak / 0.03% rms against the native 1 cm^-1 tables over the
# synthesised range (lambda < 50 um); the shipped 250 cm^-1 grid was
# 82% peak / 11% rms, dominated by the sharp H2 fundamental at 4160.
NU_MIN, NU_MAX, D_NU = 25.0, 20000.0, 25.0

LOG_FLOOR = -99.0          # log10 cm^5 written where a pair has no opacity

# Temperature at which the H2-H2 seam-correction extrapolation stops and
# is held constant: two seam-steps above Abel's 3000 K ceiling.
HOT_FREEZE = 5000.0


def temperature_grid(pair):
    """Per-pair temperature grid.

    100 K spacing over the range where CIA matters holds the temperature
    interpolation error to ~0.2% rms / 1.2% peak; it is coarsened where
    n(H2) has collapsed.  H2-H and H-He keep their native points -- the
    underlying calculations provide only 4 and 10 temperatures, and
    resampling would invent structure.
    """
    if pair == 'H2-H2':
        return np.concatenate([np.arange(200., 3001., 100.),
                               np.array([3250., 3500., 3750., 4000., 4500.,
                                         5000., 5500., 6000., 6500., 7000.])])
    if pair == 'H2-He':
        return np.concatenate([np.arange(200., 4001., 100.),
                               np.arange(4250., 6001., 250.),
                               np.arange(6500., 9501., 500.),
                               np.array([9900.])])
    if pair == 'H2-H':
        return np.array([1000., 1500., 2000., 2500.])
    if pair == 'H-He':
        return np.array([1500., 2250., 3000., 4000., 5000., 6000.,
                         7000., 8000., 9000., 10000.])
    raise ValueError(pair)


def read_hitran_cia(path):
    """HITRAN .cia block format -> (T[nt], nu[nnu], k[nt, nnu]) in cm^5.

    Every block of these four files shares one wavenumber grid; asserted
    rather than assumed.
    """
    temps, nus, ks = [], None, []
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            npts = int(header[40:47])
            temps.append(float(header[47:54]))
            block = np.array([fh.readline().split() for _ in range(npts)],
                             dtype=float)
            if nus is None:
                nus = block[:, 0]
            elif not np.allclose(nus, block[:, 0]):
                raise ValueError(f'{path}: block wavenumber grids differ')
            ks.append(block[:, 1])
    order = np.argsort(temps)
    return np.array(temps)[order], nus, np.array(ks)[order]


def read_turbospectrum_cia(path):
    """Turbospectrum COM-v19.1/DATA/CIA.*.dat -> (T, nu, k) in cm^5."""
    with open(path) as fh:
        lines = fh.readlines()
    nwave = int(lines[1].split()[0])
    ntemp = int(lines[2].split()[0])
    temps = np.array([float(x) for x in lines[3].split()[:ntemp]])
    rows = []
    for line in lines[5:]:
        parts = line.split()
        if len(parts) < ntemp + 1:
            continue
        rows.append([float(p) for p in parts[:ntemp + 1]])
        if len(rows) == nwave:
            break
    rows = np.array(rows)
    if rows.shape != (nwave, ntemp + 1):
        raise ValueError(f'{path}: expected {nwave}x{ntemp + 1}, got {rows.shape}')
    return temps, rows[:, 0], rows[:, 1:].T / N_LOSCHMIDT**2


def resample(T_src, nu_src, logk_src, T_out, nu_out):
    """Bilinear in (nu, T) on log10 k, clamped at every edge.

    Clamping is what the consuming code does too, so building the table
    this way means the file edge and the runtime edge agree.
    """
    it = np.clip(np.searchsorted(T_src, T_out) - 1, 0, len(T_src) - 2)
    wt = np.clip((T_out - T_src[it]) / (T_src[it + 1] - T_src[it]), 0.0, 1.0)
    inu = np.clip(np.searchsorted(nu_src, nu_out) - 1, 0, len(nu_src) - 2)
    wn = np.clip((nu_out - nu_src[inu]) / (nu_src[inu + 1] - nu_src[inu]), 0.0, 1.0)
    lo = logk_src[np.ix_(it, inu)] * (1 - wn) + logk_src[np.ix_(it, inu + 1)] * wn
    hi = logk_src[np.ix_(it + 1, inu)] * (1 - wn) + logk_src[np.ix_(it + 1, inu + 1)] * wn
    return lo * (1 - wt[:, None]) + hi * wt[:, None]


def extend_low_nu(nu_out, logk, nu_data_min):
    """Below the tabulated range use the exact nu^2 low-frequency limit.

    alpha(nu) = (4 pi^2 nu / 3 hbar c)[1 - exp(-hc nu / kT)] n1 n2 g(nu)
    -> nu^2 g(0) as nu -> 0, and g is smooth there.  Only reached below
    ~50 cm^-1 (lambda > 200 um), where the Planck function is dead.
    """
    below = nu_out < nu_data_min
    if not below.any():
        return logk
    anchor = np.argmin(np.abs(nu_out - nu_data_min))
    scale = 2.0 * np.log10(nu_out[below] / nu_out[anchor])
    logk[:, below] = logk[:, [anchor]] + scale[None, :]
    return logk


def extend_high_nu(nu_out, logk, nu_data_max, fit_span=1000.0):
    """Above the tabulated range continue the exponential band decay.

    Fits log10 k linearly in nu over the last `fit_span` cm^-1 of real
    data and continues that slope, forced non-positive so the
    extrapolation can only decay.  Used for H2-H and H-He, whose
    calculations stop at 10^4 and 1.1e4 cm^-1; a hard cut to zero there
    would put a step in the continuum inside the J band.
    """
    above = nu_out > nu_data_max
    if not above.any():
        return logk
    fit = (nu_out >= nu_data_max - fit_span) & (nu_out <= nu_data_max)
    edge = np.argmin(np.abs(nu_out - nu_data_max))
    x = nu_out[fit]
    for i in range(logk.shape[0]):
        slope = np.polyfit(x, logk[i, fit], 1)[0]
        logk[i, above] = logk[i, edge] + min(slope, 0.0) * (nu_out[above] - nu_out[edge])
    return logk


def fit_seam_slope(T_a, nu_a, la, T_b, nu_b, lb, nu_out, floor=200.0):
    """d log10(Abel/BJF01) / dT, in dex per 1000 K, averaged over nu.

    Fitted only at the temperatures that are exact grid points of BOTH
    tables (1000, 2000, 3000 K), so no interpolation of either source
    enters the slope, and only where BJF01 carries real signal.
    """
    shared = [t for t in T_b if t <= T_a.max() and t >= T_a.min()]
    ratio = []
    for t in shared:
        ra = resample(T_a, nu_a, la, np.array([t]), nu_out)[0]
        rb = resample(T_b, nu_b, lb, np.array([t]), nu_out)[0]
        ratio.append(ra - rb)
    ratio = np.array(ratio)
    good = (nu_out >= floor) & (nu_out <= min(nu_a.max(), nu_b.max()))
    slopes = np.polyfit(np.array(shared), ratio[:, good], 1)[0] * 1000.0
    return float(np.median(slopes))


def build_h2h2(raw, nu_out, hot_scale=1.0):
    """Abel 2011 inside its box, ratio-continued BJF01 outside.

    `hot_scale` multiplies the T > 3000 K continuation only, and exists
    to sweep the one genuinely uncertain choice in this table.  The
    Abel/BJF01 ratio is not constant -- it falls from ~0.85 at 1000 K to
    ~0.62 at 2000 K to ~0.44 at 3000 K -- and holding it at its 3000 K
    value above the box therefore likely overestimates, since the trend
    would carry it lower.  Any model whose answer depends on this knob
    is resting on an extrapolation, not on Abel.
    """
    T_out = temperature_grid('H2-H2')
    T_a, nu_a, k_a = read_hitran_cia(os.path.join(raw, 'H2-H2_2011.cia'))
    T_b, nu_b, k_b = read_turbospectrum_cia(os.path.join(raw, 'CIA.H2-H2.Yi.dat'))
    la, lb = np.log10(k_a), np.log10(k_b)

    abel = resample(T_a, nu_a, la, T_out, nu_out)
    bjf = resample(T_b, nu_b, lb, T_out, nu_out)

    # Seam offsets, in log space, along each boundary of Abel's box.
    hot = T_out > T_a.max()
    blue = nu_out > nu_a.max()
    # ratio along the T = 3000 K face, as a function of nu
    face_T = resample(T_a, nu_a, la, np.array([T_a.max()]), nu_out)[0] \
        - resample(T_b, nu_b, lb, np.array([T_a.max()]), nu_out)[0]
    # ratio along the nu = 10^4 face, as a function of T
    face_nu = resample(T_a, nu_a, la, T_out, np.array([nu_a.max()]))[:, 0] \
        - resample(T_b, nu_b, lb, T_out, np.array([nu_a.max()]))[:, 0]
    # the corner offset governs the doubly-outside quadrant
    corner = face_T[np.argmin(np.abs(nu_out - nu_a.max()))]

    # The seam correction is not constant in T: log10(Abel/BJF01) falls
    # almost linearly, and with a nearly wavenumber-independent slope
    # (-0.142 dex per 1000 K, spread 0.029 over 300-9500 cm^-1, three
    # point fits good to 0.01 dex).  Freezing it above the seam would
    # therefore bias the continuation high.  Continue the measured trend
    # instead, using one global slope -- per-nu slopes would inject the
    # noise of the poorly determined fits near 8000 cm^-1 -- and stop
    # extrapolating at HOT_FREEZE, two seam-steps out, beyond which
    # n(H2)^2 is negligible and the claim would outrun the evidence.
    slope = fit_seam_slope(T_a, nu_a, la, T_b, nu_b, lb, nu_out)
    dT = np.clip(T_out - T_a.max(), 0.0, HOT_FREEZE - T_a.max()) / 1000.0
    trend = slope * dT + np.log10(hot_scale)

    out = abel.copy()
    out[np.ix_(hot, ~blue)] = bjf[np.ix_(hot, ~blue)] + face_T[None, ~blue] \
        + trend[hot, None]
    out[np.ix_(~hot, blue)] = bjf[np.ix_(~hot, blue)] + face_nu[~hot, None]
    out[np.ix_(hot, blue)] = bjf[np.ix_(hot, blue)] + corner + trend[hot, None]

    # BJF01 itself starts at 1000 K, so the blue quadrant below 1000 K is
    # clamped by `resample`.  Decay it in log T instead of holding it flat:
    # the overtone bands that live above 10^4 cm^-1 weaken as the gas
    # cools, and a clamp would leave a cold-model overestimate.
    cold = T_out < T_b.min()
    if cold.any() and blue.any():
        ref = np.argmin(np.abs(T_out - T_b.min()))
        nxt = np.argmin(np.abs(T_out - 2000.0))
        slope = (out[nxt, blue] - out[ref, blue]) / \
                np.log10(T_out[nxt] / T_out[ref])
        slope = np.maximum(slope, 0.0)          # only ever decay downwards
        dlog = np.log10(T_out[cold] / T_out[ref])
        out[np.ix_(cold, blue)] = out[ref, blue][None, :] + slope[None, :] * dlog[:, None]

    out = extend_low_nu(nu_out, out, min(nu_a.min(), nu_b.min()))
    return T_out, out


def build_simple(raw, nu_out, pair, filename):
    """Pairs taken straight from one HITRAN file."""
    T_out = temperature_grid(pair)
    T_s, nu_s, k_s = read_hitran_cia(os.path.join(raw, filename))
    logk = resample(T_s, nu_s, np.log10(k_s), T_out, nu_out)
    logk = extend_low_nu(nu_out, logk, nu_s.min())
    logk = extend_high_nu(nu_out, logk, nu_s.max())
    return T_out, logk


PAIRS = [
    ('H2-H2', 'H2 number density (squared)',
     'Abel, Frommhold, Li & Hunt 2011, J.Phys.Chem.A 115, 6805'
     ' [+ Borysow, Jorgensen & Fu 2001 continuation]'),
    ('H2-He', 'He I number density',
     'Abel, Frommhold, Li & Hunt 2012, J.Phys.Chem.A 116, 3068'),
    ('H2-H', 'H I number density',
     'Gustafsson & Frommhold 2003, A&A 400, 1161'),
    ('H-He', 'H I and He I number densities',
     'Gustafsson & Frommhold 2001, ApJ 546, 1168'),
]


def build_all(raw, hot_scale=1.0):
    nu_out = np.arange(NU_MIN, NU_MAX + 0.5 * D_NU, D_NU)
    tables = {}
    T, k = build_h2h2(raw, nu_out, hot_scale=hot_scale)
    tables['H2-H2'] = (T, k)
    tables['H2-He'] = build_simple(raw, nu_out, 'H2-He', 'H2-He_2011.cia')
    tables['H2-H'] = build_simple(raw, nu_out, 'H2-H', 'H2-H_2011.cia')
    tables['H-He'] = build_simple(raw, nu_out, 'H-He', 'He-H_2011.cia')
    for name, (T, k) in tables.items():
        if not np.all(np.isfinite(k)):
            raise ValueError(f'{name}: non-finite entries')
        tables[name] = (T, np.maximum(k, LOG_FLOOR))
    return nu_out, tables


def write_table(path, nu_out, tables):
    with open(path, 'w') as fh:
        fh.write('# Collision-induced absorption (CIA) tables.\n')
        fh.write('# Generated by tools/build_cia_table.py -- see that file for\n')
        fh.write('# full provenance, unit and convention notes.  Do not hand-edit.\n')
        fh.write('#\n')
        fh.write('# Stored value is log10 of the binary absorption coefficient in\n')
        fh.write('# cm^5:  kappa [cm^2/g] = 10**value * n_1 * n_2 / rho,  with n in\n')
        fh.write('# cm^-3.  STIMULATED EMISSION IS ALREADY INCLUDED -- the consuming\n')
        fh.write('# code must not apply a (1 - exp(-h nu / kT)) factor to these.\n')
        fh.write('#\n')
        fh.write('# Sources:\n')
        for name, _, ref in PAIRS:
            fh.write(f'#   {name:6s} {ref}\n')
        fh.write('#\n')
        fh.write('# Layout: NPAIR and the largest per-pair temperature count, then\n')
        fh.write('# the shared wavenumber grid (NNU, NU_MIN, D_NU), then per pair:\n')
        fh.write('# name, NTEMP, the temperature list, and NTEMP records of NNU\n')
        fh.write('# values (one record per temperature, ascending nu).  Read\n')
        fh.write('# list-directed.  Below NU_MIN the nu^2 low frequency limit\n')
        fh.write('# applies; above the last point the coefficient is zero.\n')
        fh.write('# Temperatures clamp at both ends.\n')
        ntmax = max(len(T) for T, _ in tables.values())
        fh.write(f'{len(tables)} {ntmax}\n')
        fh.write(f'{len(nu_out)} {NU_MIN:.1f} {D_NU:.1f}\n')
        for name, _, _ in PAIRS:
            T, k = tables[name]
            fh.write(f'{name} {len(T)}\n')
            for i in range(0, len(T), 8):
                fh.write(' '.join(f'{t:9.1f}' for t in T[i:i + 8]) + '\n')
            for row in k:
                for i in range(0, len(row), 8):
                    fh.write(' '.join(f'{v:10.4f}' for v in row[i:i + 8]) + '\n')


def validate(raw, nu_out, tables):
    """Checks that would have caught a units or convention slip."""
    ok = True

    def check(label, got, want, tol):
        nonlocal ok
        good = abs(got - want) <= tol
        ok &= good
        print(f'  [{"ok" if good else "FAIL"}] {label}: {got:.4f} (want {want:.4f} +/- {tol})')

    print('Saumon et al. 2012 sec 2.1: H2-H2 fundamental band weaker than')
    print('BJF01 by 24% at 1000 K and 44% at 2000 K')
    T_b, nu_b, k_b = read_turbospectrum_cia(os.path.join(raw, 'CIA.H2-H2.Yi.dat'))
    T_a, nu_a, k_a = read_hitran_cia(os.path.join(raw, 'H2-H2_2011.cia'))
    for T0, want in ((1000., 0.76), (2000., 0.56)):
        ia, ib = np.argmin(np.abs(T_a - T0)), np.argmin(np.abs(T_b - T0))
        ja, jb = np.argmin(np.abs(nu_a - 4160.)), np.argmin(np.abs(nu_b - 4160.))
        check(f'  Abel/BJF01 at 4160 cm^-1, {T0:.0f} K', k_a[ia, ja] / k_b[ib, jb],
              want, 0.02)

    print('low frequency limit is nu^2 (stimulated emission included)')
    for name, (T, k) in tables.items():
        i = len(T) // 2
        j1 = np.argmin(np.abs(nu_out - 50.)), np.argmin(np.abs(nu_out - 100.))
        slope = (k[i, j1[1]] - k[i, j1[0]]) / np.log10(nu_out[j1[1]] / nu_out[j1[0]])
        check(f'  {name} slope at T={T[i]:.0f} K', slope, 2.0, 0.35)

    print('round trip through the written file')
    for name, (T, k) in tables.items():
        print(f'  {name:6s} {len(T):3d} T x {k.shape[1]:4d} nu, '
              f'log10 k in [{k.min():.1f}, {k.max():.1f}]')

    print('continuity across the H2-H2 continuation seams')
    T, k = tables['H2-H2']
    j = np.argmin(np.abs(nu_out - 10000.))
    for i in (np.argmin(np.abs(T - 2000.)), np.argmin(np.abs(T - 3000.))):
        jump = k[i, j + 1] - k[i, j]
        check(f'  d log10 k across nu=10^4 at {T[i]:.0f} K', jump, 0.0, 0.05)
    i = np.argmin(np.abs(T - 3000.))
    for jj in (np.argmin(np.abs(nu_out - 4160.)), np.argmin(np.abs(nu_out - 6000.))):
        jump = k[i + 1, jj] - k[i, jj]
        check(f'  d log10 k across T=3000 at {nu_out[jj]:.0f} cm^-1', jump, 0.0, 0.12)

    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--raw', required=True, help='directory of downloaded source tables')
    ap.add_argument('--out', default=None, help='output path (default: no write)')
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--hot-scale', type=float, default=1.0,
                    help='sensitivity knob on the H2-H2 T>3000 K continuation')
    args = ap.parse_args()

    nu_out, tables = build_all(args.raw, hot_scale=args.hot_scale)
    if args.hot_scale != 1.0:
        print(f'*** H2-H2 continuation above 3000 K scaled by {args.hot_scale} ***')
    print(f'wavenumber grid: {len(nu_out)} points, '
          f'{nu_out[0]:.0f}-{nu_out[-1]:.0f} cm^-1 at {D_NU:.0f} cm^-1')

    ok = True
    if args.validate:
        ok = validate(args.raw, nu_out, tables)

    if args.out:
        write_table(args.out, nu_out, tables)
        print(f'wrote {args.out} '
              f'({os.path.getsize(args.out) / 1e6:.2f} MB)')

    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
