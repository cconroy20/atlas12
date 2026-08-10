#!/usr/bin/env python3
"""Plot T(tau) for the ExoMol Toto vs Schwenke TiO line lists at 2800 K.

Reads either converged .atm models or, while a run is still iterating,
the last complete block of its .iter file (same T and log tau columns),
so the comparison can be previewed before convergence.

Usage:
    python3 tools/tio_ttau_plot.py new=<dir> old=<dir> [--out fig.png]

Each <dir> holds t2800.atm and/or t2800.iter.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Okabe-Ito, colour-blind safe.  Blue is reserved for our new/RE model.
BLUE = '#0072B2'
ORANGE = '#E69F00'
GREY = '#666666'

plt.rcParams.update({
    'font.family': 'serif',
    'mathtext.fontset': 'stix',
    'font.serif': ['STIXGeneral'],
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.grid': False,
})


def read_atm(path):
    """Converged model -> (log10 tau5000, T)."""
    lines = open(path).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith('READ DECK')][0]
    rows = []
    for l in lines[i + 1:]:
        p = l.split()
        if len(p) < 11:
            break
        try:
            rows.append([float(x) for x in p[:11]])
        except ValueError:
            break
    d = np.array(rows)
    return np.log10(d[:, 10]), d[:, 1]


def read_iter(path):
    """Last complete iteration block of a .iter file -> (log10 tau, T).

    Blocks restart at 'J log10TAU'; a run still in progress may have a
    truncated final block, so the last COMPLETE one is used and the
    iteration index is returned for honest labelling.
    """
    lines = open(path).read().splitlines()
    starts = [k for k, l in enumerate(lines) if l.lstrip().startswith('J log10TAU')]
    blocks = []
    for n, s in enumerate(starts):
        end = starts[n + 1] if n + 1 < len(starts) else len(lines)
        rows = []
        for l in lines[s + 2:end]:
            p = l.split()
            if len(p) < 3:
                continue
            try:
                rows.append([float(p[1]), float(p[2])])
            except ValueError:
                continue
        if rows:
            blocks.append(np.array(rows))
    if not blocks:
        raise ValueError(f'{path}: no complete iteration block')
    nfull = max(len(b) for b in blocks)
    full = [b for b in blocks if len(b) == nfull]
    return full[-1][:, 0], full[-1][:, 1], len(full)


def load(d):
    atm = os.path.join(d, 't2800.atm')
    if os.path.exists(atm) and os.path.getsize(atm) > 0:
        tau, T = read_atm(atm)
        return tau, T, 'converged'
    tau, T, n = read_iter(os.path.join(d, 't2800.iter'))
    return tau, T, f'iteration {n}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+', help='new=<dir> old=<dir>')
    ap.add_argument('--out', default='tio_ttau.png')
    args = ap.parse_args()

    data = {}
    for a in args.dirs:
        label, d = a.split('=', 1)
        data[label] = load(d)
        print(f'{label:5s} {data[label][2]:>14s}  {len(data[label][0])} layers  '
              f'T {data[label][1].min():.0f}-{data[label][1].max():.0f} K')

    tn, Tn, sn = data['new']
    to, To, so = data['old']
    n = min(len(tn), len(to))
    dT = Tn[:n] - To[:n]

    fig, (ax, axd) = plt.subplots(
        2, 1, figsize=(7.0, 6.4), sharex=True,
        gridspec_kw={'height_ratios': [2.4, 1], 'hspace': 0.08})

    ax.plot(to, To, color=ORANGE, lw=1.8, label=f'Schwenke 1998 TiO ({so})')
    ax.plot(tn, Tn, color=BLUE, lw=1.8, label=f'ExoMol Toto 2024 TiO ({sn})')
    ax.set_ylabel(r'$T$  [K]')
    ax.legend(frameon=False, loc='upper left', fontsize=9)
    ax.set_title('2800 K dwarf: ATLAS12 structure vs TiO line list', fontsize=11)

    axd.axhline(0.0, color=GREY, lw=0.8, ls='--')
    axd.plot(tn[:n], dT, color=BLUE, lw=1.8)
    axd.set_xlabel(r'$\log_{10}\,\tau_{5000}$')
    axd.set_ylabel(r'$\Delta T$  [K]')

    for a in (ax, axd):
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.tick_params(direction='in', which='both')

    lo = max(dT.max(), -dT.min(), 1.0) * 1.25
    axd.set_ylim(-lo, lo)
    axd.text(0.98, 0.06,
             f'Toto $-$ Schwenke:  max {dT[np.argmax(np.abs(dT))]:+.0f} K, '
             f'rms {np.sqrt((dT**2).mean()):.0f} K',
             transform=axd.transAxes, ha='right', fontsize=8.5, color=GREY)

    fig.savefig(args.out, bbox_inches='tight')
    print(f'wrote {args.out}')
    j = np.argmax(np.abs(dT))
    print(f'max |dT| {dT[j]:+.1f} K at log tau = {tn[j]:.2f}; '
          f'rms {np.sqrt((dT**2).mean()):.2f} K')
    return 0


if __name__ == '__main__':
    sys.exit(main())
