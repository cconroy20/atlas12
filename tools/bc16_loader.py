"""
bc16_loader.py -- Barklem & Collet (2016) molecular tables.

Reconstruction of the April 2026 module deleted after commit bb386c9.
Interface fixed by the surviving comp_pf.py and
build_molecules_physical_dat.py call sites; conventions validated by
reproducing the bb386c9 molecules.dat (see tools/rebuild notes).

Provides:
  KB_PaCm3_per_K   Boltzmann constant in Pa cm^3 / K
  load_bc16()      dict name -> BC16Mol with .name, .T_grid, .log10_Kp,
                   .D0_adopted, .atom_a, .atom_b
"""

from pathlib import Path

import numpy as np

KB_PaCm3_per_K = 1.380649e-17     # k_B [Pa cm^3 / K]

BC16_DIR = Path.home() / 'kurucz' / 'upgrade' / 'raw_data'
TABLE7 = BC16_DIR / 'bc16_table7_vNov2022.dat'
TABLE1 = BC16_DIR / 'bc16_table1.dat'


class BC16Mol:
    """One BC16 molecule: Table 7 log10 K_p(T) [Pa] plus Table 1 D0."""

    def __init__(self, name, T_grid, log10_Kp, D0_adopted, atom_a, atom_b):
        self.name = name
        self.T_grid = T_grid
        self.log10_Kp = log10_Kp
        self.D0_adopted = D0_adopted
        self.atom_a = atom_a          # dissociation products from Table 1,
        self.atom_b = atom_b          # e.g. ('H', 'H+') for H2+


def load_bc16():
    """name -> BC16Mol for every species in Table 7."""
    lines = TABLE7.read_text().split('\n')
    T_grid = np.array([float(x) for x in lines[2].split()[2:]])

    # Table 1: adopted D0 (10th column) and the two dissociation products
    d0 = {}
    products = {}
    for line in TABLE1.read_text().split('\n')[4:]:
        parts = line.split()
        if len(parts) >= 10:
            name = parts[0]
            products[name] = (parts[1], parts[2])
            d0[name] = float(parts[9]) if parts[9] != '.' else float('nan')

    bc16 = {}
    for line in lines[3:]:
        parts = line.split()
        if len(parts) != len(T_grid) + 1:
            continue
        name = parts[0]
        a, b = products.get(name, (None, None))
        bc16[name] = BC16Mol(
            name, T_grid, np.array([float(x) for x in parts[1:]]),
            d0.get(name, float('nan')), a, b)
    return bc16
