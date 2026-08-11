#!/usr/bin/env python3
"""
exomol_to_kurucz.py -- convert an ExoMol line list (.states + .trans) to the
Kurucz molecular ascii format read by read_molec_ascii (mod_mklinelist.f90).

Output record (FORMAT F10.4,F7.3,F5.1,F10.3,F5.1,F11.3,I4,A8,A8,I2,I4):
  wl        vacuum wavelength [nm] (reader recomputes from energies; the
            column is used for windowing, so the FILE MUST BE SORTED
            ascending in wl -- the reader stops at the first line beyond
            wlend)
  gflog     log10(g_lower * f_lu), NUCLEAR-SPIN-FREE and with NO isotope
            abundance weighting (molec_dispatch applies the isotopic
            x1+x2+fudge correction at ingestion)
  xj,  e    J and energy [cm^-1] of the LOWER level
  xjp, ep   J and energy [cm^-1] of the UPPER level; negative energy marks
            a computed (non-measured) level, Kurucz convention -- the
            reader uses MIN(ABS(e),ABS(ep)), so only cosmetics downstream
  icode     Kurucz species code: Z_small*100 + Z_big (CaH=120, TiO=822)
  labels    lower/upper level labels '<State><vv><e/f><F-comp>', e.g.
            'X00e1'; the reader tests clabelp(1:1)=='X' to pick default
            Stark/vdW damping, so the leading electronic-state letter of
            the UPPER label matters
  iso       isotope tag keyed by molec_dispatch (heavy-atom mass number
            for hydrides/oxides, e.g. 40 for 40CaH, 48 for 48Ti16O)
  loggr     100*log10(Gamma_rad [1/s]) of the line, 0 => classical.
            Computed from the states-file lifetime column (Gamma = 1/tau
            of the upper level) when available.

Conventions (must not drift):
  gf = 1.4991938e-16 * (g_upper/g_ns) * A_ul [1/s] * (lambda [A])^2
  ExoMol g_tot includes the nuclear-spin factor g_ns; BC16/Kurucz
  partition functions exclude it, so gf here must exclude it too --
  pass the correct --gns (product of (2I+1) over nuclei; 2 for CaH).

ExoMol .states columns 1-4 are always (id, E [cm^-1], g_tot, J).  The
optional columns (e/f parity, State, v, Sigma, lifetime, source tag) vary
by dataset; they are auto-detected from the first data row and can be
overridden with --*-col flags (1-based indices).  Levels whose source tag
is NOT in --measured-tags get negative (computed) energies when a source
column is present.

Example (Owens et al. 2022 CaH XAB):
  python3 exomol_to_kurucz.py --states 40Ca-1H__XAB.states.bz2 \\
      --trans 40Ca-1H__XAB.trans.bz2 --gns 2 --out cah_owens22.dat
icode/iso default from the isotopologue filename ('40Ca-1H' -> 120 / 40).
"""

import argparse
import bz2
import re
import sys

import numpy as np

GF_CONST = 1.4991938e-16     # m_e*c/(8*pi^2*e^2) in s/A^2  (gf = C*g_u*A*lam_A^2)
WL_FMT_MAX = 99999.9999      # nm; largest wavelength the F10.4 field holds

ELEMENTS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Y': 39, 'Zr': 40, 'La': 57,
}


def zopen(path):
    return bz2.open(path, 'rt') if path.endswith('.bz2') else open(path)


def parse_isotopologue(name):
    """'40Ca-1H' or '40Ca-16O-1H' -> (atom list, icode, iso).

    icode is the Kurucz species code: atomic numbers sorted ascending and
    concatenated as 2-digit fields after the first (CaH: 1,20 -> 120;
    H2O: 1,1,8 -> 10108; CaOH: 1,8,20 -> 10820).  iso = mass number of
    the highest-Z atom.  Codes > 9999 need the 'polymol' manifest tag.
    """
    pairs = re.findall(r'(\d+)([A-Z][a-z]?)', name)
    if len(pairs) < 2:
        return None
    atoms = []
    for m, el in pairs:
        z = ELEMENTS.get(el)
        if z is None:
            return None
        atoms.append((int(m), el, z))
    zs = sorted(a[2] for a in atoms)
    icode = int(str(zs[0]) + ''.join(f'{z:02d}' for z in zs[1:]))
    iso = max(atoms, key=lambda a: a[2])[0]
    return atoms, icode, iso


def autodetect_columns(tokens):
    """Guess optional state-file columns from one data row (1-based)."""
    guess = {}
    for i, tok in enumerate(tokens[4:], start=5):
        if tok in ('e', 'f') and 'ef' not in guess:
            guess['ef'] = i
        elif re.match(r'^[A-Za-z]+\d?\(', tok) or re.match(r"^[A-Za-z]['`]?p?$", tok):
            guess.setdefault('state', i)
        elif re.match(r'^[A-Z][a-z]$', tok) and i > 4:
            guess.setdefault('source', i)
    # v: first pure-integer column after the state label
    if 'state' in guess:
        for i in range(guess['state'] + 1, len(tokens) + 1):
            if re.match(r'^-?\d+$', tokens[i - 1]):
                guess['v'] = i
                break
        # Sigma: first half-integer-valued float after v (Sigma precedes Omega)
        if 'v' in guess:
            for i in range(guess['v'] + 1, len(tokens) + 1):
                if re.match(r'^-?\d+\.5$', tokens[i - 1]) or tokens[i - 1] in ('0.5', '-0.5'):
                    guess['sigma'] = i
                    break
    # lifetime: a column that parses as float with values spanning many
    # decades or 'Inf' -- detect by literal 'Inf' or exponent format
    for i, tok in enumerate(tokens[4:], start=5):
        if tok in ('Inf', 'inf', 'INF') or re.match(r'^\d\.\d+E[+-]\d+$', tok):
            guess['tau'] = i
            break
    return guess


def build_label(state, v, ef, fcomp):
    lab = (state or 'U')[:1].upper()
    if v is not None:
        lab += f'{min(v, 99):02d}'
    if ef:
        lab += ef
    if fcomp:
        lab += fcomp
    return lab[:8]


def main():
    ap = argparse.ArgumentParser(
        description='Convert ExoMol .states/.trans to Kurucz molecular ascii')
    ap.add_argument('--states', required=True)
    ap.add_argument('--trans', required=True, nargs='+',
                    help='one or more .trans files (large lists come split)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--gns', type=float, required=True,
                    help='nuclear-spin degeneracy included in ExoMol g_tot '
                         '(product of 2I+1; 2 for CaH, 1 for 48Ti16O)')
    ap.add_argument('--icode', type=int, default=None,
                    help='Kurucz species code (default: from filename)')
    ap.add_argument('--iso', type=int, default=None,
                    help='isotope tag for molec_dispatch (default: from filename)')
    ap.add_argument('--wlmin', type=float, default=0.0, help='nm')
    ap.add_argument('--wlmax', type=float, default=1.0e6, help='nm')
    ap.add_argument('--gflog-min', type=float, default=-99.0,
                    help='drop lines weaker than this log gf (format floor)')
    for key in ('ef', 'state', 'v', 'sigma', 'tau', 'source'):
        ap.add_argument(f'--{key}-col', type=int, default=None,
                        help=f'1-based states-file column for {key} (overrides autodetect)')
    ap.add_argument('--measured-tags', default='Ma',
                    help='comma-separated source tags treated as measured '
                         '(others get negative energies)')
    args = ap.parse_args()

    # --- species identification -----------------------------------------
    m = re.search(r'([\d]+[A-Z][a-z]?-[\d]+[A-Z][a-z]?)', args.states)
    parsed = parse_isotopologue(m.group(1)) if m else None
    icode = args.icode if args.icode is not None else (parsed and parsed[1])
    iso = args.iso if args.iso is not None else (parsed and parsed[2])
    if icode is None or iso is None:
        sys.exit('error: could not infer icode/iso from filename; '
                 'pass --icode and --iso explicitly')
    print(f'species: icode={icode}  iso={iso}  gns={args.gns}')

    # --- read states -----------------------------------------------------
    with zopen(args.states) as fh:
        rows = [ln.split() for ln in fh if ln.strip()]
    cols = autodetect_columns(rows[0])
    for key in ('ef', 'state', 'v', 'sigma', 'tau', 'source'):
        ovr = getattr(args, f'{key}_col')
        if ovr is not None:
            cols[key] = ovr
    print(f'states columns (1-based): id=1 E=2 g=3 J=4  optional={cols}')

    nst = max(int(r[0]) for r in rows) + 1
    E = np.full(nst, np.nan)
    g = np.zeros(nst)
    J = np.zeros(nst)
    lab = np.empty(nst, dtype=object)
    measured = np.ones(nst, dtype=bool)
    gamrad = np.zeros(nst)
    tags = set(args.measured_tags.split(','))
    for r in rows:
        i = int(r[0])
        E[i] = float(r[1]); g[i] = float(r[2]); J[i] = float(r[3])
        state = r[cols['state'] - 1] if 'state' in cols else None
        v = int(r[cols['v'] - 1]) if 'v' in cols else None
        ef = r[cols['ef'] - 1] if 'ef' in cols else ''
        fcomp = ''
        if 'sigma' in cols:
            try:
                fcomp = '1' if float(r[cols['sigma'] - 1]) >= 0 else '2'
            except ValueError:
                pass
        lab[i] = build_label(state, v, ef, fcomp)
        if 'source' in cols:
            measured[i] = r[cols['source'] - 1] in tags
        if 'tau' in cols:
            try:
                tau = float(r[cols['tau'] - 1])
                if np.isfinite(tau) and tau > 0:
                    gamrad[i] = 1.0 / tau
            except ValueError:
                pass
    print(f'states: {len(rows)}  (measured: {measured[[int(r[0]) for r in rows]].sum()})')

    # --- read transitions ------------------------------------------------
    up, lo, A = [], [], []
    for tf in args.trans:
        with zopen(tf) as fh:
            for ln in fh:
                p = ln.split()
                if len(p) < 3:
                    continue
                up.append(int(p[0])); lo.append(int(p[1])); A.append(float(p[2]))
    up = np.array(up); lo = np.array(lo); A = np.array(A)
    print(f'transitions: {len(up)}')

    # --- physics ---------------------------------------------------------
    nu = np.abs(E[up] - E[lo])                      # cm^-1
    ok = nu > 1.0e-6
    up, lo, A, nu = up[ok], lo[ok], A[ok], nu[ok]
    wl_nm = 1.0e7 / nu
    # The Kurucz wavelength field is F10.4, so 100000 nm (100 um) does not
    # fit: such a line shifts every later field and the Fortran reader
    # rejects the whole record.  Drop them here rather than write records
    # that only look valid.  They are pure-rotational far-IR transitions,
    # far outside any photospheric flux this code models.
    ok = (wl_nm >= args.wlmin) & (wl_nm <= min(args.wlmax, WL_FMT_MAX))
    n_farir = int(((wl_nm > WL_FMT_MAX) & (wl_nm <= args.wlmax)).sum())
    up, lo, A, nu, wl_nm = up[ok], lo[ok], A[ok], nu[ok], wl_nm[ok]
    if n_farir:
        print(f'dropped {n_farir} lines beyond {WL_FMT_MAX:.0f} nm '
              f'(F10.4 wavelength field limit)')

    gf = GF_CONST * (g[up] / args.gns) * A * (wl_nm * 10.0) ** 2
    gflog = np.log10(np.maximum(gf, 1.0e-99))
    ok = gflog >= args.gflog_min
    n_weak = (~ok).sum()
    up, lo, nu, wl_nm, gflog = up[ok], lo[ok], nu[ok], wl_nm[ok], gflog[ok]

    order = np.argsort(wl_nm)
    up, lo, nu, wl_nm, gflog = up[order], lo[order], nu[order], wl_nm[order], gflog[order]

    loggr = np.where(gamrad[up] > 0,
                     np.rint(100.0 * np.log10(np.maximum(gamrad[up], 1.0))),
                     0).astype(int)

    # --- write -----------------------------------------------------------
    # Diatomics use the classic I4 code field (read_molec_ascii); codes
    # over 9999 get the I6 'polymol' format (read_molec_poly).
    poly = icode > 9999
    cfmt = f'{icode:6d}' if poly else f'{icode:4d}'
    with open(args.out, 'w') as out:
        for k in range(len(up)):
            u, l = up[k], lo[k]
            e_lo = E[l] if measured[l] else -E[l]
            e_up = E[u] if measured[u] else -E[u]
            out.write(f'{wl_nm[k]:10.4f}{gflog[k]:7.3f}{J[l]:5.1f}{e_lo:10.3f}'
                      f'{J[u]:5.1f}{e_up:11.3f}{cfmt}{lab[l]:<8s}{lab[u]:<8s}'
                      f'{iso:2d}{loggr[k]:4d}\n')

    print(f'wrote {len(up)} lines to {args.out} '
          f'({wl_nm[0]:.1f}-{wl_nm[-1]:.1f} nm; dropped {n_weak} below '
          f'gflog={args.gflog_min})')
    if poly:
        print(f'NOTE: code {icode} > 9999 -> POLYATOMIC format written; '
              f'list it in lines.list with the "polymol" tag (and make sure '
              f'poly_nelion in mod_mklinelist.f90 has a case for {icode})')


if __name__ == '__main__':
    main()
