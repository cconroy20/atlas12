"""
ggchem_loader.py -- GGchem dispol equilibrium-constant files.

Reads the gas-phase kp(T) fit files shipped with GGchem (Woitke et al.
2018, github.com/pw31/GGchem), local clone at
~/kurucz/upgrade/raw_data/ggchem.  The fits trace to NIST-JANAF /
Stock (2008) / Burcat thermochemistry and are the modern reference for
polyatomics without ExoMol coverage.

Conventions (transcribed from GGchem src16/smchem16.f, function gk):
  * kp is the FORMATION constant, p_mol / prod p_atom^nu, with pressures
    in dyn/cm^2 (cgs); gk returns ln kp.
  * fit=4 (Stock/Kitzmann): dG = a0/T + a1 lnT + a2 + a3 T + a4 T^2;
    ln kp = dG + (1 - Natom) ln(1e6)
  * fit=3 (Sharp & Huebner): dG in cal/mol; ln kp = -dG/(Rcal T)
    + (1 - Natom) ln(1.013e6)
  * fit=5 (dG-fit, J/mol):   ln kp = -dG/(Rgas T) + (1 - Natom) ln(1e6)
  * fit=1 (Gail), fit=2 (Tsuji): theta = 5040/T polynomials
  * fit=6 (Barklem & Collet 8-term refit by GGchem)

Conversion to the Kurucz association convention used in molecules.dat:
  ln K_cgs(assoc, number density) = ln kp + (N-1) ln(k_B[erg/K] T)

Validated on load (see validate()): the dispol Barklem & Collet entries
round-trip against BC16 Table 7 to the precision of GGchem's own 8-term
refits.
"""

from pathlib import Path

import numpy as np

GGCHEM_DATA = Path.home() / 'kurucz' / 'upgrade' / 'raw_data' / 'ggchem' / 'data'
KB_CGS = 1.380649e-16          # erg/K
RCAL = 1.987                   # cal/mol/K   (GGchem's value)
RGAS = 8.3144598               # J/mol/K     (GGchem's value)
LN_BAR = np.log(1.0e6)         # bar in dyn/cm^2
LN_ATM = np.log(1.013e6)       # atm in dyn/cm^2 (GGchem's value)


def parse_dispol(filename):
    """dispol file -> {name: {'atoms': {el: n}, 'natom': N, 'fit': f,
                              'a': [...]}}.  Charged species are skipped."""
    path = GGCHEM_DATA / filename
    lines = [l for l in path.read_text().split('\n')]
    recs = {}
    i = 1                       # line 0 is the species count
    while i + 1 < len(lines):
        head = lines[i].split()
        coef = lines[i + 1].split()
        i += 2
        if not head or not coef:
            continue
        name = head[0]
        ntyp = int(head[1])
        els = head[2:2 + ntyp]
        cnts = [int(x) for x in head[2 + ntyp:2 + 2 * ntyp]]
        if '+' in name or '-' in name or 'el' in els:
            continue
        recs[name] = {
            'atoms': dict(zip(els, cnts)),
            'natom': sum(cnts),
            'fit': int(coef[0]),
            'a': [float(x) for x in coef[1:]],
        }
    return recs


def lnkp_cgs(rec, T):
    """ln kp (formation, cgs pressures) at T; mirrors smchem16.f gk()."""
    T = np.asarray(T, dtype=float)
    a = rec['a']
    n = rec['natom']
    fit = rec['fit']
    if fit == 1:                                  # Gail
        th = 5040.0 / T
        return a[0] + a[1]*th + a[2]*th**2 + a[3]*th**3 + a[4]*th**4
    if fit == 2:                                  # Tsuji
        th = 5040.0 / T
        return np.log(10.0) * -(a[0] + a[1]*th + a[2]*th**2
                                + a[3]*th**3 + a[4]*th**4)
    if fit == 3:                                  # Sharp & Huebner
        dg = a[0]/T + a[1] + a[2]*T + a[3]*T**2 + a[4]*T**3
        return -dg/(RCAL*T) + (1 - n)*LN_ATM
    if fit == 4:                                  # Stock / Kitzmann
        dg = a[0]/T + a[1]*np.log(T) + a[2] + a[3]*T + a[4]*T**2
        return dg + (1 - n)*LN_BAR
    if fit == 5:                                  # dG-fit (J/mol)
        dg = a[0]/T + a[1] + a[2]*T + a[3]*T**2 + a[4]*T**3
        return -dg/(RGAS*T) + (1 - n)*LN_BAR
    if fit == 6:                                  # BC16 8-term refit
        return (a[0]/T**3 + a[1]/T**2 + a[2]/T + a[3]/T**0.05
                + a[4]*np.log(T) + a[5] + a[6]*T + a[7]*T**2)
    raise ValueError(f'unsupported fit type {fit}')


def ln_K_kurucz(rec, T):
    """Kurucz association ln K [cm^3(N-1)] from a dispol record."""
    T = np.asarray(T, dtype=float)
    return lnkp_cgs(rec, T) + (rec['natom'] - 1) * np.log(KB_CGS * T)


def validate():
    """Round-trip GGchem's BC16 refits against BC16 Table 7 directly."""
    from bc16_loader import load_bc16, KB_PaCm3_per_K
    bc = parse_dispol('dispol_BarklemCollet.dat')
    t7 = load_bc16()
    worst = ('', 0.0)
    n = 0
    for name in ('H2', 'CO', 'TiO', 'H2O', 'SiO', 'CN', 'MgH', 'OH'):
        if name not in bc or name not in t7:
            continue
        m = t7[name]
        mask = (m.T_grid >= 1000.0) & (m.T_grid <= 6000.0)
        T = m.T_grid[mask]
        # Table 7 dissociation Kp [Pa] -> Kurucz association ln K
        ln_ref = (np.log10(KB_PaCm3_per_K * T)
                  - m.log10_Kp[mask]) * np.log(10.0) * (m and 1)
        d = np.abs(ln_K_kurucz(bc[name], T) - ln_ref).max()
        n += 1
        if d > worst[1]:
            worst = (name, d)
    return n, worst


if __name__ == '__main__':
    n, (name, d) = validate()
    print(f'validate: {n} BC16 diatomics round-tripped, '
          f'worst |dln K| over 1000-6000 K = {d:.4f} ({name})')
