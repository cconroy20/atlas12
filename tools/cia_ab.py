#!/usr/bin/env python3
"""Same-structure A/B for the CIA upgrade at 2700 K.

Compares one or more ATLAS12/SYNTHE spectra computed on the PHOENIX
NewEra 2700/5.0/0.0 structure (so only the opacity differs) against each
other and against the PHOENIX spectrum computed on that same structure.

Usage:
    python3 tools/cia_ab.py label=path [label=path ...] [--win a,b ...]
"""
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mann_lib import read_spec_file        # noqa: E402

PHX = os.path.expanduser('~/kurucz/grids/NewEra_HSR/phoenix_hires_2700g5.0.npz')

# diagnostic windows, Angstroms
WINDOWS = [
    ('optical 5000-7800', 5000., 7800.),
    ('I  8000-9000', 8000., 9000.),
    ('J  11500-13500', 11500., 13500.),
    ('H  15500-17000', 15500., 17000.),
    ('K  21000-23000', 21000., 23000.),
]


def band_ratio(wl, fl, wref, fref, lo, hi):
    """Median of model/PHOENIX over a window, on the PHOENIX grid."""
    m = (wref >= lo) & (wref <= hi)
    if m.sum() == 0:
        return np.nan
    interp = np.interp(wref[m], wl, fl)
    good = np.isfinite(interp) & (fref[m] > 0)
    if good.sum() == 0:
        return np.nan
    return float(np.median(interp[good] / fref[m][good]))


def main():
    args = [a for a in sys.argv[1:] if '=' in a]
    if not args:
        print(__doc__)
        return 1

    phx = np.load(PHX)
    wref, fref = phx['wl'], phx['flam']

    models = {}
    for a in args:
        label, path = a.split('=', 1)
        wl, fl, fc = read_spec_file(path)
        models[label] = (np.asarray(wl), np.asarray(fl), np.asarray(fc))
        print(f'{label:22s} {len(wl):8d} points, '
              f'{wl.min():.1f}-{wl.max():.1f} A   {path}')
    print()

    hdr = f'{"window":22s}' + ''.join(f'{k:>14s}' for k in models)
    print(hdr)
    print('-' * len(hdr))
    for name, lo, hi in WINDOWS:
        cells = []
        for label, (wl, fl, _) in models.items():
            if wl.max() < lo or wl.min() > hi:
                cells.append(f'{"-":>14s}')
            else:
                cells.append(f'{band_ratio(wl, fl, wref, fref, lo, hi):14.4f}')
        print(f'{name:22s}' + ''.join(cells))
    print('\n(values are median ours/PHOENIX on the same structure)')

    # pairwise ratio between models, where they overlap.  The continuum
    # column is the one to read for a continuum-opacity change: it is
    # free of the line blanketing that dilutes the flux ratio.
    labels = list(models)
    if len(labels) > 1:
        print()
        base = labels[0]
        for label in labels[1:]:
            wl0, fl0, fc0 = models[base]
            wl1, fl1, fc1 = models[label]
            lo, hi = max(wl0.min(), wl1.min()), min(wl0.max(), wl1.max())
            print(f'{label} / {base}:')
            print(f'   {"window":22s} {"flux":>18s} {"continuum":>10s}')
            for name, a, b in WINDOWS:
                a, b = max(a, lo), min(b, hi)
                if a >= b:
                    continue
                m = (wl0 >= a) & (wl0 <= b)
                r = np.interp(wl0[m], wl1, fl1) / fl0[m]
                rc = np.interp(wl0[m], wl1, fc1) / fc0[m]
                print(f'   {name:22s} {np.median(r):8.4f} '
                      f'[{np.percentile(r, 5):.4f},{np.percentile(r, 95):.4f}]'
                      f' {np.median(rc):9.4f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
