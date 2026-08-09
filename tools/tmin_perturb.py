#!/usr/bin/env python3
"""
tmin_perturb.py -- attach a two-parameter lower-chromosphere/T-min structure
to a converged ATLAS12 model (.atm).

Parameterization (user-specified, 2026-08-06; --tmax added 2026-08-07):
  tau >= tau_min :  T = T_model(tau)                              (untouched)
  tau <  tau_min :  T = T_model(tau_min)
                        + min(s * (log10(tau_min) - log10(tau)), dT_max)

dT_max = inf reproduces the original two-parameter linear ramp; finite
dT_max gives a VAL-like morphology (rise above the T-min, then a plateau).

i.e. the model is exact down to the temperature-minimum location tau_min
(in tau_5000, the TAU5000 column of DECK6), and outward of it T itself is
a straight line in log tau anchored at T(tau_min) with slope s [K/dex] --
strictly monotonic rise, VAL-like morphology.  s=0 gives an isothermal
outer atmosphere at T(tau_min); the perturbation never cools below the
anchor temperature.

Only the T column is modified; the pressure/density structure is left in
place (SYNTHE recomputes XNE from T,P).  The perturbed model is a
DIAGNOSTIC, not a flux-conserving solution.

Usage:
  python3 tmin_perturb.py in.atm out.atm --logtau-min -2.0 --slope 200
"""

import argparse
import numpy as np


def _read_deck(in_atm, coord):
    lines = open(in_atm).readlines()
    istart = next(i for i, l in enumerate(lines)
                  if l.startswith('READ DECK')) + 1
    ndep = int(lines[istart - 1].split()[2])
    lts, Ts = [], []
    for k in range(ndep):
        p = lines[istart + k].split()
        depth = float(p[0]) if coord == 'cmass' else float(p[-1])
        lts.append(np.log10(depth))
        Ts.append(float(p[1]))
    return lines, istart, ndep, np.array(lts), np.array(Ts)


def _write_deck(lines, istart, ndep, Tnew, out_atm):
    for k in range(ndep):
        p = lines[istart + k].split()
        rhox = float(p[0])
        rest = [float(x) for x in p[2:]]
        lines[istart + k] = (f'{rhox:25.5E}{Tnew[k]:10.2f}' +
                             ''.join(f'{v:10.3E}' for v in rest) + '\n')
    open(out_atm, 'w').write(''.join(lines))


def perturb_nodes(in_atm, out_atm, knots, smooth_dex=0.0, coord='tau'):
    """Anchored monotone node profile (Piette/Molliere morphology, 2026-08-08).

    knots: [(log_depth, val), ...] ordered inward->outward (log depth
    decreasing).  The FIRST knot is the anchor: val = SIGNED offset dT_a
    applied at the anchor depth (allows a true T-min dipping below RE).
    Each subsequent val is a cumulative increment dT_i >= 0, i.e.
    T(knot_k) = T_model(anchor) + dT_a + sum(dT_1..dT_k).  Between knots T
    is PCHIP-interpolated in log depth (shape-preserving -- no spline
    overshoot); outward of the last knot the atmosphere is isothermal at
    the last knot's T (declared, not fit); inward of the anchor the model
    is untouched.  smooth_dex > 0 gaussian-smooths the resulting deltaT in
    log depth (sigma in dex) to soften the anchor kink/step (Piette+20 use
    0.3 dex).

    The old 2-parameter ramp is the special case knots=[(a,0),(b,s*(a-b))]
    with the last knot at the top of the grid and smooth_dex=0.
    """
    from scipy.interpolate import PchipInterpolator
    lines, istart, ndep, lts, Ts = _read_deck(in_atm, coord)
    kx = np.array([k[0] for k in knots], float)
    kv = np.array([k[1] for k in knots], float)
    if not np.all(np.diff(kx) < 0):
        raise ValueError('knots must be ordered inward->outward '
                         '(log depth strictly decreasing)')
    if np.any(kv[1:] < 0):
        raise ValueError('increments dT_1.. must be >= 0 (only the anchor '
                         'offset dT_a may be negative)')
    T_anchor = float(np.interp(kx[0], lts[np.argsort(lts)],
                               Ts[np.argsort(lts)]))
    kT = T_anchor + np.cumsum(kv)
    # PCHIP in outward coordinate x = -log depth (increasing outward)
    pch = PchipInterpolator(-kx, kT)
    Tnew = Ts.copy()
    inside = (lts <= kx[0]) & (lts >= kx[-1])
    Tnew[inside] = pch(-lts[inside])
    Tnew[lts < kx[-1]] = kT[-1]
    if smooth_dex > 0:
        from scipy.ndimage import gaussian_filter1d
        dT = Tnew - Ts
        step = float(np.median(np.abs(np.diff(lts))))
        dT = gaussian_filter1d(dT, smooth_dex / step, mode='nearest')
        Tnew = Ts + dT
    _write_deck(lines, istart, ndep, Tnew, out_atm)


def perturb(in_atm, out_atm, logtau_min, slope, dtmax=np.inf, coord='tau'):
    """coord='tau': anchor/ramp in log10 TAU5000; coord='cmass': in log10 RHOX
    (column mass, the VAL-native depth variable) — logtau_min is then the
    anchor in log10 m [g cm^-2] and slope is K per dex of column mass."""
    lines = open(in_atm).readlines()
    istart = next(i for i, l in enumerate(lines) if l.startswith('READ DECK')) + 1
    ndep = int(lines[istart - 1].split()[2])
    # first pass: T at the anchor depth (interpolate T_model in log depth)
    lts, Ts = [], []
    for k in range(ndep):
        p = lines[istart + k].split()
        depth = float(p[0]) if coord == 'cmass' else float(p[-1])
        lts.append(np.log10(depth)); Ts.append(float(p[1]))
    T_anchor = float(np.interp(logtau_min, lts, Ts))
    for k in range(ndep):
        p = lines[istart + k].split()
        rhox, T = float(p[0]), float(p[1])
        rest = [float(x) for x in p[2:]]
        lt = np.log10(rhox) if coord == 'cmass' else np.log10(rest[-1])
        if lt < logtau_min:
            T = T_anchor + min(slope * (logtau_min - lt), dtmax)
        lines[istart + k] = (f'{rhox:25.5E}{T:10.2f}' +
                             ''.join(f'{v:10.3E}' for v in rest) + '\n')
    open(out_atm, 'w').write(''.join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('in_atm')
    ap.add_argument('out_atm')
    ap.add_argument('--logtau-min', type=float,
                    help='log10 tau_5000 of the temperature minimum')
    ap.add_argument('--slope', type=float,
                    help='outward temperature rise, K per dex of log tau')
    ap.add_argument('--knots', type=str, default=None,
                    help='node-profile mode: "lt:dTa,lt:dT1,lt:dT2,..." -- '
                         'first pair = anchor depth + signed offset, later '
                         'pairs = outward depths + increments >= 0')
    ap.add_argument('--smooth-dex', type=float, default=0.0,
                    help='gaussian smoothing of deltaT in dex (node mode)')
    ap.add_argument('--tmax', type=float, default=np.inf,
                    help='saturation ceiling dT_max in K (default: none)')
    ap.add_argument('--coord', choices=['tau', 'cmass'], default='tau',
                    help='depth variable for anchor+ramp (default tau5000)')
    args = ap.parse_args()
    if args.knots:
        knots = [tuple(float(v) for v in pair.split(':'))
                 for pair in args.knots.split(',')]
        perturb_nodes(args.in_atm, args.out_atm, knots, args.smooth_dex,
                      args.coord)
        print(f'{args.out_atm}: node profile {knots}, '
              f'smooth={args.smooth_dex} dex')
    else:
        if args.logtau_min is None or args.slope is None:
            ap.error('either --knots or both --logtau-min and --slope')
        perturb(args.in_atm, args.out_atm, args.logtau_min, args.slope,
                args.tmax, args.coord)
        print(f'{args.out_atm}: T-min at logtau={args.logtau_min}, '
              f'slope={args.slope} K/dex')


if __name__ == '__main__':
    main()
