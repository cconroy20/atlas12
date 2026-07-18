"""
exomol_loader.py -- local ExoMol partition function files.

Reconstruction of the April 2026 module deleted after commit bb386c9.
Reads the .pf files in <root>/raw_data/exomol_pf/ (two columns: T [K],
Q(T), typically on a 1 K grid).

  ExoMolData(root)
    .has(name)                name = .pf basename ('H2O', 'CaOH', 'H3p')
    .molecules[name]          {'T': grid, 'Q': values, 'Tmax': float}
    .Q(name, T)               interpolated Q; NaN outside native range

  assemble_log10_Kp_exomol(atom_a, atom_b, D0, Qfunc, atomic, T)
    diatomic wrapper around the polyatomic assembly (used for ExoMol
    overlays on diatomic pages when a diatomic .pf is present).
"""

from pathlib import Path

import numpy as np

from polyatomic_assembly import assemble_log10_Kp_polyatomic


class ExoMolData:

    def __init__(self, root='.'):
        self.molecules = {}
        pf_dir = Path(root) / 'raw_data' / 'exomol_pf'
        for pf in sorted(pf_dir.glob('*.pf')):
            data = np.loadtxt(pf)
            self.molecules[pf.stem] = {
                'T': data[:, 0], 'Q': data[:, 1],
                'Tmax': float(data[:, 0].max()),
            }

    def has(self, name):
        return name in self.molecules

    def Q(self, name, T):
        m = self.molecules[name]
        T = np.asarray(T, dtype=float)
        q = np.interp(T, m['T'], m['Q'], left=np.nan, right=np.nan)
        return q


def assemble_log10_Kp_exomol(atom_a, atom_b, D0_eV, Qfunc, atomic, T):
    """Diatomic assembly: strip any ion tag (BC16 Table 1 products can be
    e.g. 'H+'); the neutral-channel assembly is all comp_pf uses."""
    atoms = [atom_a.rstrip('+'), atom_b.rstrip('+')]
    return assemble_log10_Kp_polyatomic(atoms, D0_eV, Qfunc, atomic, T)
