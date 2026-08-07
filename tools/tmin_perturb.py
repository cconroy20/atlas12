#!/usr/bin/env python3
"""
tmin_perturb.py -- attach a two-parameter lower-chromosphere/T-min structure
to a converged ATLAS12 model (.atm).

Parameterization (user-specified, 2026-08-06):
  tau >= tau_min :  T = T_model(tau)                              (untouched)
  tau <  tau_min :  T = T_model(tau_min) + s * (log10(tau_min) - log10(tau))

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


def perturb(in_atm, out_atm, logtau_min, slope):
    lines = open(in_atm).readlines()
    istart = next(i for i, l in enumerate(lines) if l.startswith('READ DECK')) + 1
    ndep = int(lines[istart - 1].split()[2])
    # first pass: T at the anchor depth (interpolate T_model in log tau)
    lts, Ts = [], []
    for k in range(ndep):
        p = lines[istart + k].split()
        lts.append(np.log10(float(p[-1]))); Ts.append(float(p[1]))
    T_anchor = float(np.interp(logtau_min, lts, Ts))
    for k in range(ndep):
        p = lines[istart + k].split()
        rhox, T = float(p[0]), float(p[1])
        rest = [float(x) for x in p[2:]]
        lt = np.log10(rest[-1])            # TAU5000 is the last DECK6 column
        if lt < logtau_min:
            T = T_anchor + slope * (logtau_min - lt)
        lines[istart + k] = (f'{rhox:25.5E}{T:10.2f}' +
                             ''.join(f'{v:10.3E}' for v in rest) + '\n')
    open(out_atm, 'w').write(''.join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('in_atm')
    ap.add_argument('out_atm')
    ap.add_argument('--logtau-min', type=float, required=True,
                    help='log10 tau_5000 of the temperature minimum')
    ap.add_argument('--slope', type=float, required=True,
                    help='outward temperature rise, K per dex of log tau')
    args = ap.parse_args()
    perturb(args.in_atm, args.out_atm, args.logtau_min, args.slope)
    print(f'{args.out_atm}: T-min at logtau={args.logtau_min}, '
          f'slope={args.slope} K/dex')


if __name__ == '__main__':
    main()
