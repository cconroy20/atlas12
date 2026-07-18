"""
atomic_saha.py -- atomic partition functions, masses, and Saha factors.

Reconstruction of the April 2026 module deleted after commit bb386c9.
Data: BC16 Table 8 (Nov 2022) atomic partition functions and Table 4
ionization energies, both in ~/kurucz/upgrade/raw_data/ -- the same
files the April pipeline downloaded.

Interface (fixed by comp_pf.py / polyatomic_assembly.py call sites):
  AtomicData()
    .Q(symbol, T, ion=0)      partition function of X I/II/III at T
    .mass(symbol)             atomic mass [amu]
    .log10_K_Saha(symbol, T)  log10[ n(X+) n_e / n(X) ] in cm^-3

Interpolation of the Table 8 grid is a reconstruction unknown; the
scheme below (cubic spline of log10 Q vs log10 T) is selected by
INTERP_MODE and validated against the bb386c9 coefficient reproduction.
"""

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

K_BOLTZMANN_EV = 8.617333262e-5

# Electron translational prefactor (2 pi m_e k / h^2)^{3/2} [cm^-3 K^-1.5],
# CODATA 2018
SAHA_TR_E = 2.4147187e15

BC16_DIR = Path.home() / 'kurucz' / 'upgrade' / 'raw_data'
TABLE8 = BC16_DIR / 'bc16_table8_vNov2022.dat'
TABLE4 = BC16_DIR / 'bc16_table4.dat'

# Interpolation scheme for the 42-point Table 8 grid.  Candidates tried
# during reconstruction: 'loglog' (spline of log10 Q vs log10 T),
# 'loglin' (spline of Q vs log10 T), 'linear' (np.interp of Q vs log10 T).
INTERP_MODE = 'loglog'

# Standard atomic weights [amu] for species appearing in assemblies
ATOMIC_MASS = {
    'H': 1.008, 'He': 4.003, 'Li': 6.94, 'Be': 9.012, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'K': 39.098, 'Ca': 40.078, 'Sc': 44.956,
    'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845,
    'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38, 'Sr': 87.62,
    'Y': 88.906, 'Zr': 91.224, 'Ba': 137.33, 'La': 138.91,
}

_ROMAN = {0: 'I', 1: 'II', 2: 'III'}


class AtomicData:

    def __init__(self):
        lines = TABLE8.read_text().split('\n')
        # Row 3 (index 2) is the T grid; data rows follow, labelled X_I etc.
        self._T = np.array([float(x) for x in lines[2].split()[2:]])
        self._logT = np.log10(self._T)
        self._Q = {}
        self._interp = {}
        for line in lines[3:]:
            parts = line.split()
            if len(parts) == len(self._T) + 1 and '_' in parts[0]:
                self._Q[parts[0]] = np.array([float(x) for x in parts[1:]])

        # Table 4: first three ionization energies [eV]
        self._ie = {}
        for line in TABLE4.read_text().split('\n'):
            if line.startswith('*') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 5:
                self._ie[parts[1]] = [float(parts[2]), float(parts[3]),
                                      float(parts[4])]

    def _interpolator(self, key):
        if key not in self._interp:
            q = self._Q[key]
            if INTERP_MODE == 'loglog':
                cs = CubicSpline(self._logT, np.log10(q))
                self._interp[key] = lambda lt: 10.0 ** cs(lt)
            elif INTERP_MODE == 'pchip':
                from scipy.interpolate import PchipInterpolator
                cs = PchipInterpolator(self._logT, np.log10(q))
                self._interp[key] = lambda lt: 10.0 ** cs(lt)
            elif INTERP_MODE == 'spline_T':
                cs = CubicSpline(self._T, q)
                self._interp[key] = lambda lt: cs(10.0 ** lt)
            elif INTERP_MODE == 'loglin':
                cs = CubicSpline(self._logT, q)
                self._interp[key] = cs
            else:
                self._interp[key] = lambda lt, q=q: np.interp(
                    lt, self._logT, q)
        return self._interp[key]

    def Q(self, symbol, T, ion=0):
        """Partition function of `symbol` at ionization stage `ion`."""
        key = f'{symbol}_{_ROMAN[ion]}'
        T = np.asarray(T, dtype=float)
        return self._interpolator(key)(np.log10(T))

    def mass(self, symbol):
        return ATOMIC_MASS[symbol]

    def ionization_energy(self, symbol, stage=1):
        """IE of the `stage`-th ionization [eV] (1 = neutral -> +)."""
        return self._ie[symbol][stage - 1]

    def log10_K_Saha(self, symbol, T):
        """
        log10 of the Saha constant n(X+) n_e / n(X) in cm^-3:
        2 * (2 pi m_e k T / h^2)^{3/2} * U_II/U_I * exp(-IE/kT)
        """
        T = np.asarray(T, dtype=float)
        u1 = self.Q(symbol, T, 0)
        u2 = self.Q(symbol, T, 1)
        ie = self.ionization_energy(symbol, 1)
        return (np.log10(2.0 * SAHA_TR_E) + 1.5 * np.log10(T)
                + np.log10(u2 / u1)
                - ie / (K_BOLTZMANN_EV * T) / np.log(10.0))
