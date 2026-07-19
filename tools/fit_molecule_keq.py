"""
fit_molecule_keq.py -- fit molecules.dat physical-form K_eq rows from BC16.

Reconstructs the diatomic path of the April 2026 refit pipeline
(~/kurucz/upgrade/build_molecules_physical_dat.py, commit bb386c9) that
produced the "BC16 fit" rows in data/molecules.dat.  Kept in the repo so
the fit provenance survives with the data it generates; the original
pipeline's helper modules (bc16_loader.py, matcher.py, ...) were cleaned
up after that commit and no longer exist.

Conventions (must not drift -- see --validate):

  BC16 Table 7 stores log10(K_p) for DISSOCIATION in Pa
  (Barklem & Collet 2016, A&A 588, A96; Nov 2022 bug-fixed version).

  molecules.dat stores the ASSOCIATION constant in number density,
  n_mol = K * n_A * n_B, evaluated by NMOLEC as

    ln K_eq(T) = D0/(k_B T) - 1.5 * n_trans * ln T
               + a0 + a1/T + a2 ln T + a3 T + a4/T^2

  with D0 in eV, k_B in eV/K, and n_trans = N_atoms - ION - 1.

  Conversion:  log10 K = log10(k_B[Pa cm^3/K] * T) - log10 K_p(Table 7)

  Fit: subtract the fixed part D0/(k_B T) - 1.5 n_trans ln T, then
  unweighted lstsq of the residual onto {1, 1/T, ln T, T, 1/T^2} on the
  BC16 native grid restricted to T >= 100 K (19 points, 100..10^4 K).
  The D0 subtracted is the D0 WRITTEN TO THE ROW: E1 and a0..a4 are a
  matched set, so a row's D0 cannot be changed without refitting.

Scope: NEUTRAL diatomics only.  Molecular ions use a different D0 sign
convention and the polyatomic path needs the (lost) ExoMol machinery.

Usage:
  python3 tools/fit_molecule_keq.py --validate
      Refit every neutral-diatomic "BC16 fit" row in data/molecules.dat
      and require exact agreement with the filed coefficients at the
      printed E12.4 precision.  Run this before trusting new fits.

  python3 tools/fit_molecule_keq.py --fit CrO CrS TiN TiS TiH
      Print formatted molecules.dat rows for the named species, taking
      D0 from BC16 Table 1 (adopted column).

  python3 tools/fit_molecule_keq.py --fit ... --write
      Also insert the rows into data/molecules.dat, keeping the file's
      species-code ordering.  Refuses to overwrite an existing code.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

KB_PA_CM3_PER_K = 1.380649e-17     # k_B [Pa cm^3 / K]
K_BOLTZMANN_EV = 8.617333262e-5    # k_B [eV / K]
T_FIT_MIN = 100.0                  # BC16 diatomic fit mask

REPO_ROOT = Path(__file__).resolve().parent.parent
MOLECULES_DAT = REPO_ROOT / 'data' / 'molecules.dat'
BC16_DIR = Path.home() / 'kurucz' / 'upgrade' / 'raw_data'

ELEMENT_Z = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Ba': 56, 'La': 57,
}


def parse_formula(name):
    """'TiO' -> [22, 8]; 'C2' -> [6, 6].  Neutral diatomics only."""
    if name.endswith('+') or name.endswith('-'):
        raise ValueError(f'{name}: ions not supported (D0 sign convention)')
    atoms = []
    for sym, count in re.findall(r'([A-Z][a-z]?)(\d?)', name):
        if not sym:
            continue
        if sym not in ELEMENT_Z:
            raise ValueError(f'{name}: unknown element {sym}')
        atoms += [ELEMENT_Z[sym]] * (int(count) if count else 1)
    if len(atoms) != 2:
        raise ValueError(f'{name}: not a diatomic ({len(atoms)} atoms)')
    return atoms


def species_code(atoms):
    """Kurucz code: atomic-number digit pairs in ascending order."""
    a, b = sorted(atoms)
    return a * 100 + b


def load_table7():
    """BC16 Table 7: name -> (T_grid, log10 K_p [Pa], dissociation)."""
    lines = (BC16_DIR / 'bc16_table7_vNov2022.dat').read_text().split('\n')
    tgrid = np.array([float(x) for x in lines[2].split()[2:]])
    table = {}
    for line in lines[3:]:
        parts = line.split()
        if len(parts) == len(tgrid) + 1:
            table[parts[0]] = np.array([float(x) for x in parts[1:]])
    return tgrid, table


def load_table1_d0():
    """BC16 Table 1: name -> adopted D0 [eV]."""
    d0 = {}
    for line in (BC16_DIR / 'bc16_table1.dat').read_text().split('\n')[4:]:
        parts = line.split()
        if len(parts) >= 10 and parts[9] != '.':
            d0[parts[0]] = float(parts[9])
    return d0


def bc16_name_for(atoms, table):
    """Match an atom multiset against BC16 species names (handles SH/HS)."""
    want = sorted(atoms)
    for name in table:
        if name.endswith('+') or name.endswith('-'):
            continue
        try:
            if sorted(parse_formula(name)) == want:
                return name
        except ValueError:
            continue
    return None


def fit_diatomic(d0_row, log10_kp, tgrid, n_trans=1):
    """The _do_fit kernel from build_molecules_physical_dat.py."""
    mask = tgrid >= T_FIT_MIN
    t = tgrid[mask]
    ln_k = (np.log10(KB_PA_CM3_PER_K * t) - log10_kp[mask]) * np.log(10.0)
    fixed = d0_row / (K_BOLTZMANN_EV * t) - 1.5 * n_trans * np.log(t)
    basis = np.column_stack(
        [np.ones_like(t), 1.0 / t, np.log(t), t, 1.0 / t**2])
    coeffs, *_ = np.linalg.lstsq(basis, ln_k - fixed, rcond=None)
    return coeffs


def format_row(label, code, d0, coeffs, comment):
    body = f'{label:<10s}{code:12.2f} {d0:7.3f}' \
           + ''.join(f'{c:12.4E}' for c in coeffs)
    return f'{body:<90s}   ! {comment}'


def parse_molecules_rows():
    """Yield (line_index, label, code, d0, coeffs, comment) for data rows."""
    lines = MOLECULES_DAT.read_text().split('\n')
    for i, line in enumerate(lines):
        if not line.strip() or line.startswith('!'):
            continue
        label = line[:10].strip()
        code = float(line[10:22])
        d0 = float(line[23:30]) if line[23:30].strip() else 0.0
        coeffs = None
        if len(line) >= 90 and line[30:90].strip():
            coeffs = [float(line[30 + 12 * k:42 + 12 * k]) for k in range(5)]
        bang = line.find('!')
        comment = line[bang + 1:].strip() if bang > 0 else ''
        yield i, label, code, d0, coeffs, comment
    return


# Saha channel (the atom that carries the charge) for the molecular-ion
# rows, per BC16 Table 1 dissociation products: the lower-IP partner.
ION_SAHA_CHANNEL = {'H2+': 'H', 'CH+': 'C', 'OH+': 'O', 'SiH+': 'Si',
                    'CN+': 'C', 'N2+': 'N', 'NO+': 'O', 'O2+': 'O'}


def fit_diatomic_ion(e1, log10_kp, tgrid, saha_atom):
    """
    Molecular-cation fit in the Kurucz e- convention:
    K = n(M+) n_e / prod n(atoms) = [BC16 association K] * Saha(atom).
    n_trans = 0 (M+ + e- vs two atoms).  E1 = D0(BC16) - IE(atom).
    """
    from atomic_saha import AtomicData
    atomic = _ion_atomic.setdefault('a', AtomicData())
    mask = tgrid >= T_FIT_MIN
    t = tgrid[mask]
    ln_k = (np.log10(KB_PA_CM3_PER_K * t) - log10_kp[mask]) * np.log(10.0)
    ln_k = ln_k + atomic.log10_K_Saha(saha_atom, t) * np.log(10.0)
    fixed = e1 / (K_BOLTZMANN_EV * t)          # n_trans = 0
    basis = np.column_stack(
        [np.ones_like(t), 1.0 / t, np.log(t), t, 1.0 / t**2])
    coeffs, *_ = np.linalg.lstsq(basis, ln_k - fixed, rcond=None)
    return coeffs


_ion_atomic = {}


def cmd_validate():
    tgrid, table7 = load_table7()
    n_ok = n_bad = n_skip = 0
    for _, label, code, d0, coeffs, comment in parse_molecules_rows():
        if ('BC16 fit' not in comment and 'BC16+Saha' not in comment) \
                or coeffs is None:
            continue
        ion = int(round(code * 100)) % 100
        atoms = []
        n = int(round(code * 100)) // 100
        while n > 0:
            atoms.append(n % 100)
            n //= 100
        if len(atoms) != 2 or 100 in atoms:
            n_skip += 1
            continue
        if ion == 0:
            name = bc16_name_for(atoms, table7)
            if name is None:
                print(f'  {label:<6s} NO BC16 MATCH')
                n_bad += 1
                continue
            refit = fit_diatomic(d0, table7[name], tgrid)
        elif label in ION_SAHA_CHANNEL and label in table7:
            refit = fit_diatomic_ion(d0, table7[label], tgrid,
                                     ION_SAHA_CHANNEL[label])
        else:
            n_skip += 1
            continue
        same = all(f'{a:12.4E}' == f'{b:12.4E}'
                   for a, b in zip(refit, coeffs))
        if same:
            n_ok += 1
        else:
            n_bad += 1
            print(f'  {label:<6s} MISMATCH:')
            for a, b in zip(refit, coeffs):
                print(f'     refit {a:12.4E}   filed {b:12.4E}')
    print(f'validate: {n_ok} exact, {n_bad} mismatched, '
          f'{n_skip} skipped')
    return 1 if n_bad else 0


def cmd_fit(names, write):
    tgrid, table7 = load_table7()
    table1 = load_table1_d0()
    existing = {code for _, _, code, _, _, _ in parse_molecules_rows()}

    new_rows = []
    for name in names:
        atoms = parse_formula(name)
        code = species_code(atoms)
        if float(code) in existing:
            sys.exit(f'{name}: code {code}.00 already in molecules.dat')
        bc_name = bc16_name_for(atoms, table7)
        if bc_name is None:
            sys.exit(f'{name}: not in BC16 Table 7')
        if bc_name not in table1:
            sys.exit(f'{name}: no adopted D0 in BC16 Table 1')
        # Round D0 to the F7.3 file precision BEFORE fitting: E1 and a0..a4
        # are a matched set (a1 is degenerate with D0/k_B), so the fit must
        # use exactly the D0 that lands in the row.  Table 1 carries up to
        # six decimals; fitting with those and writing three shifted a1 by
        # up to ~5 in its 4th significant figure.
        d0 = round(table1[bc_name], 3)
        coeffs = fit_diatomic(d0, table7[bc_name], tgrid)
        row = format_row(name, float(code), d0, coeffs,
                         'BC16 fit; D0 = BC16 Table 1; new species')
        new_rows.append((float(code), row))
        print(row)

    if not write:
        return 0

    lines = MOLECULES_DAT.read_text().split('\n')
    # Data-row line numbers and codes, for sorted insertion
    positions = [(i, code) for i, _, code, _, _, _ in parse_molecules_rows()]
    for code, row in sorted(new_rows, reverse=True):
        after = max(i for i, c in positions if c < code)
        lines.insert(after + 1, row)
    MOLECULES_DAT.write_text('\n'.join(lines))
    print(f'wrote {len(new_rows)} rows to {MOLECULES_DAT}')
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--fit', nargs='+', metavar='SPECIES')
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    if args.validate:
        sys.exit(cmd_validate())
    if args.fit:
        sys.exit(cmd_fit(args.fit, args.write))
    ap.print_help()


if __name__ == '__main__':
    main()
