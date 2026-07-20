"""
janaf_loader.py -- raw NIST-JANAF thermochemical tables (Chase 1998).

Reads the per-species JANAF .txt tables shipped with GGchem
(~/kurucz/upgrade/raw_data/ggchem/data/JANAF/) and assembles
formation-from-atoms equilibrium constants directly from the log Kf
columns -- fit-free, the primary source itself.

  JanafTable(name)     one species table: .T, .logKf, .dHf (kJ/mol),
                       .dHf0 (kJ/mol at 0 K)
  ln_K_kurucz_janaf(mol, atoms, T)
                       Kurucz association ln K [cm^3(N-1)] on grid T,
                       from log Kf(mol) - sum log Kf(atoms); JANAF Kf
                       is referenced to p_std = 1 bar
  atomization_D0(mol, atoms)
                       D0 [eV] = sum dHf0(atoms) - dHf0(mol), the
                       0 K atomization energy consistent with the same
                       tables
"""

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

JANAF_DIR = (Path.home() / 'kurucz' / 'upgrade' / 'raw_data'
             / 'ggchem' / 'data' / 'JANAF')
KB_CGS = 1.380649e-16          # erg/K
LN_BAR = np.log(1.0e6)         # bar in dyn/cm^2
KJMOL_TO_EV = 1.0 / 96.4853321233


# GGchem's JANAF folder names files by formula, which collides on a
# case-insensitive filesystem (CO carbon monoxide vs Co cobalt -- the
# only such pair).  The cobalt table is renamed in the local clone.
FILENAME_ALIAS = {'Co': 'Co_atom'}


class JanafTable:

    def __init__(self, name):
        path = JANAF_DIR / f'{FILENAME_ALIAS.get(name, name)}.txt'
        text = path.read_text()
        # Guard against filename case collisions (CO vs Co): the header's
        # formula field (e.g. "C1O1(g)") must contain the same element
        # multiset as the requested name, compared case-sensitively.
        # Element order differs between filename and header (JANAF writes
        # MgOH as H1Mg1O1), so compare sorted element tokens.
        import re as _re

        def _elements(s):
            return sorted(_re.findall(r'[A-Z][a-z]?', s.split('(')[0]))

        header_formula = text.split('\n')[0].split('\t')[-1]
        if _elements(header_formula) != _elements(name):
            raise ValueError(
                f"JANAF file {path.name} header says '{header_formula}', "
                f"not '{name}' -- filename case collision?")
        T, logkf, dhf = [], [], []
        self.dHf0 = None
        for line in text.split('\n')[2:]:
            parts = line.split('\t')
            if len(parts) < 8 or not parts[0].strip():
                continue
            t = float(parts[0])
            if t == 0.0:
                self.dHf0 = float(parts[5])
                continue
            if parts[7].strip() in ('INFINITE', ''):
                continue
            T.append(t)
            dhf.append(float(parts[5]))
            logkf.append(float(parts[7]))
        self.T = np.array(T)
        self.logKf = np.array(logkf)
        self.dHf = np.array(dhf)
        # spline in 1/T (log Kf is near-linear in 1/T); reverse so the
        # abscissa is strictly increasing
        self._cs = CubicSpline((1.0 / self.T)[::-1], self.logKf[::-1])

    def logKf_at(self, T):
        """log10 Kf interpolated in 1/T (log Kf is near-linear in 1/T)."""
        T = np.asarray(T, dtype=float)
        out = self._cs(1.0 / T)
        out = np.where((T >= self.T.min()) & (T <= self.T.max()),
                       out, np.nan)
        return out


def ln_K_kurucz_janaf(mol, atoms, T):
    """
    Kurucz association ln K for `mol` (JANAF table name) built from
    atoms (list of JANAF table names, with repetition).
    """
    T = np.asarray(T, dtype=float)
    m = JanafTable(mol) if isinstance(mol, str) else mol
    logk = m.logKf_at(T).copy()
    n = len(atoms)
    for a in atoms:
        at = JanafTable(a) if isinstance(a, str) else a
        logk = logk - at.logKf_at(T)
    # formation-from-atoms, p_std = 1 bar -> cgs -> association (number
    # density): ln K = ln kp_cgs + (N-1) ln(k_B T)
    ln_kp_cgs = logk * np.log(10.0) + (1 - n) * LN_BAR
    return ln_kp_cgs + (n - 1) * np.log(KB_CGS * T)


def atomization_D0(mol, atoms):
    """0 K atomization energy [eV] from the same tables' dHf(0 K)."""
    m = JanafTable(mol) if isinstance(mol, str) else mol
    d = -m.dHf0
    for a in atoms:
        at = JanafTable(a) if isinstance(a, str) else a
        d += at.dHf0
    return d * KJMOL_TO_EV
