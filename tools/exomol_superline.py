#!/usr/bin/env python3
"""
exomol_superline.py -- stream huge ExoMol .trans files into Kurucz-format
SUPER-LINES for the polymol reader.

For lists too large to ingest line-by-line (e.g. CaOH OYT6, ~1e10
transitions), bin gf on a fixed (wavenumber x lower-state-energy) grid and
emit one effective line per populated cell:

  gf_cell  = sum of nuclear-spin-free gf of all transitions in the cell
  nu_cell  = gf-weighted mean wavenumber within the cell
  Elo_cell = gf-weighted mean lower energy within the cell

With dnu << thermal Doppler width and dElo small against kT, synthesis
with super-lines is indistinguishable from the full list (this is the
same construction as Kurucz's 'h2ofast').  Defaults: dnu = 0.01 cm^-1
(R ~ 1.5e6 at 700 nm), dElo = 250 cm^-1 (Boltzmann error < 7% at 3000 K).

Conventions match exomol_to_kurucz.py: gf = 1.4991938e-16*(g_u/g_ns)*A*
lambda_A^2, nuclear-spin-free (--gns required).  Output is the polymol
I6-code ascii format, sorted by wavelength.  Labels are generic ('X00'
lower, 'B00' upper) -- super-lines carry no level identity; damping is
the classical default (loggr=0, non-X upper => molecular Stark/vdW
defaults in the reader).

Example (CaOH OYT6 optical chunks):
  python3 exomol_superline.py --states 40Ca-16O-1H__OYT6.states.bz2 \\
      --trans 40Ca-16O-1H__OYT6__1[46]000-1[68]000.trans.bz2 --gns 2 \\
      --out caoh_oyt6_super.dat
"""

import argparse
import re
import subprocess
import sys

import numpy as np

GF_CONST = 1.4991938e-16
sys.path.insert(0, __path__[0] if '__path__' in dir() else '.')


def zstream(path):
    """Line-block stream via external bzcat (faster than python bz2)."""
    if path.endswith('.bz2'):
        proc = subprocess.Popen(['bzcat', path], stdout=subprocess.PIPE,
                                bufsize=1024 * 1024)
        return proc.stdout
    return open(path, 'rb')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--states', required=True)
    ap.add_argument('--trans', required=True, nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--gns', type=float, required=True)
    ap.add_argument('--icode', type=int, default=None)
    ap.add_argument('--iso', type=int, default=None)
    ap.add_argument('--dnu', type=float, default=0.01, help='cm^-1 bin')
    ap.add_argument('--delo', type=float, default=250.0, help='cm^-1 bin')
    ap.add_argument('--numin', type=float, default=None)
    ap.add_argument('--numax', type=float, default=None)
    ap.add_argument('--chunk', type=int, default=5_000_000)
    args = ap.parse_args()

    from exomol_to_kurucz import parse_isotopologue, ELEMENTS  # noqa: F401
    m = re.search(r'((?:[\d]+[A-Z][a-z]?-?)+)__', args.states.split('/')[-1])
    parsed = parse_isotopologue(m.group(1)) if m else None
    icode = args.icode if args.icode is not None else (parsed and parsed[1])
    iso = args.iso if args.iso is not None else (parsed and parsed[2])
    if icode is None or iso is None:
        sys.exit('error: pass --icode/--iso (could not infer from filename)')
    print(f'species: icode={icode} iso={iso} gns={args.gns}')

    # --- states -> arrays -------------------------------------------------
    import bz2
    ids, Es, Gs = [], [], []
    op = bz2.open(args.states, 'rt') if args.states.endswith('.bz2') else open(args.states)
    for ln in op:
        p = ln.split()
        ids.append(int(p[0])); Es.append(float(p[1])); Gs.append(float(p[2]))
    nmax = max(ids)
    E = np.zeros(nmax + 2); G = np.zeros(nmax + 2)
    E[np.array(ids)] = Es; G[np.array(ids)] = Gs
    print(f'states: {len(ids)}')

    # --- accumulation grids ----------------------------------------------
    numin = args.numin; numax = args.numax
    if numin is None or numax is None:
        # infer from trans filenames like __14000-16000
        lo, hi = 1e30, -1e30
        for tf in args.trans:
            mm = re.search(r'__(\d+)-(\d+)\.trans', tf)
            if mm:
                lo = min(lo, float(mm.group(1))); hi = max(hi, float(mm.group(2)))
        if hi < 0:
            sys.exit('error: pass --numin/--numax (not inferable from filenames)')
        numin, numax = lo, hi
    nnu = int(np.ceil((numax - numin) / args.dnu))
    nel = int(np.ceil(35000.0 / args.delo))
    print(f'grid: {nnu} nu-bins x {nel} Elo-bins over {numin}-{numax} cm-1')
    SGF  = np.zeros(nnu * nel)           # sum gf
    SNU  = np.zeros(nnu * nel)           # sum gf*nu
    SEL  = np.zeros(nnu * nel)           # sum gf*Elo

    # ------------------------------------------------------------------
    # Fixed-width byte parser for the ExoMol 3-column trans format:
    #   cols  1-12: upper id (right-justified int)
    #   cols 13-25: lower id
    #   cols 26-36: A in the rigid form ' d.dddde+dd'
    # Wavenumber is NOT in the file; it is E[u]-E[l] from the states.
    # ------------------------------------------------------------------
    LINE = 37                                  # 36 chars + newline
    P10 = 10.0 ** np.arange(11, -1, -1)
    P10L = 10.0 ** np.arange(12, -1, -1)

    def parse_block(buf):
        a = np.frombuffer(buf, dtype=np.uint8).reshape(-1, LINE)
        d = np.where((a >= 48) & (a <= 57), a - 48, 0).astype(np.float64)
        u = (d[:, 0:12] @ P10).astype(np.int64)
        l = (d[:, 12:25] @ P10L).astype(np.int64)
        # ' d.dddde+dd' at cols 25..35 (0-based): digit 26, '.', 4 digits 28-31,
        # 'e' 32, sign 33, 2-digit exponent 34-35
        mant = d[:, 26] + d[:, 28] * 1e-1 + d[:, 29] * 1e-2 + d[:, 30] * 1e-3 + d[:, 31] * 1e-4
        ex = d[:, 34] * 10 + d[:, 35]
        sgn = np.where(a[:, 33] == 45, -1.0, 1.0)        # '-' = 45
        A = mant * 10.0 ** (sgn * ex)
        return u, l, A

    ntot = 0
    blocklines = args.chunk
    for tf in args.trans:
        print(f'streaming {tf} ...', flush=True)
        fh = zstream(tf)
        pending = b''
        while True:
            buf = fh.read(LINE * blocklines - len(pending))
            if not buf:
                break
            buf = pending + buf
            nfull = len(buf) // LINE
            pending = buf[nfull * LINE:]
            if nfull == 0:
                continue
            u, l, A = parse_block(buf[:nfull * LINE])
            nu = np.abs(E[u] - E[l])
            ok = (nu > max(numin, 1.0)) & (nu < numax)
            ntot += nfull
            if ok.any():
                u, l, A, nu = u[ok], l[ok], A[ok], nu[ok]
                lam_A = 1.0e8 / nu
                gf = GF_CONST * (G[u] / args.gns) * A * lam_A**2
                elo = np.minimum(E[u], E[l])
                inu = np.minimum(((nu - numin) / args.dnu).astype(np.int64), nnu - 1)
                iel = np.minimum((elo / args.delo).astype(np.int64), nel - 1)
                idx = inu * nel + iel
                np.add.at(SGF, idx, gf)
                np.add.at(SNU, idx, gf * nu)
                np.add.at(SEL, idx, gf * elo)
            if ntot % 200_000_000 < blocklines:
                print(f'  ... {ntot/1e9:.2f}e9 transitions', flush=True)
        fh.close()
    print(f'total transitions streamed: {ntot}')

    # --- emit super-lines -------------------------------------------------
    nz = np.nonzero(SGF > 0)[0]
    gf  = SGF[nz]
    nuc = SNU[nz] / gf
    elc = SEL[nz] / gf
    wl_nm = 1.0e7 / nuc
    gflog = np.log10(gf)
    keep = gflog > -99.0
    nz, gf, nuc, elc, wl_nm, gflog = nz[keep], gf[keep], nuc[keep], elc[keep], wl_nm[keep], gflog[keep]
    order = np.argsort(wl_nm)
    print(f'super-lines: {len(order)}')

    with open(args.out, 'w') as out:
        for k in order:
            e_lo = elc[k]
            e_up = elc[k] + nuc[k]
            out.write(f'{wl_nm[k]:10.4f}{gflog[k]:7.3f}{0.0:5.1f}{e_lo:10.3f}'
                      f'{0.0:5.1f}{e_up:11.3f}{icode:6d}{"X00":<8s}{"B00":<8s}'
                      f'{iso:2d}{0:4d}\n')
    print(f'wrote {args.out} ({wl_nm[order[0]]:.1f}-{wl_nm[order[-1]]:.1f} nm)')


if __name__ == '__main__':
    main()
