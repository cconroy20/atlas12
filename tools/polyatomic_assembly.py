"""
polyatomic_assembly.py -- assemble dissociation K_p from partition functions.

Reconstruction of the April 2026 module deleted after commit bb386c9;
conventions validated by reproducing bb386c9's 16 "ExoMol fit" rows.

Physics: for molecule M with atoms {A_i}, the dissociation constant in
number density is

  K_n = prod_i n(A_i) / n(M)
      = exp(-D0/kT) * [prod_i q_tr,i U_i] / [q_tr,M U_M]      [cm^-3]^(N-1)

with q_tr,X = C_TR * (m_X[amu] * T)^{3/2} cm^-3, and the BC16 pressure
convention K_p = K_n * (k_B T)^{N-1} with k_B in Pa cm^3/K.

Conventions (reconstruction unknowns, fixed by validation):
  * C_TR from CODATA 2018: (2 pi m_u k / h^2)^{3/2} = 1.8793e20
  * ExoMol partition functions carry HITRAN-convention nuclear-spin
    degeneracy; BC16 atomic partition functions do not.  NS_DIVIDE=True
    divides the molecular Q by prod_i (2 I_i + 1) of the constituent
    nuclei (main isotopes) to put both on the astrophysical convention.
"""

import numpy as np

from bc16_loader import KB_PaCm3_per_K

K_BOLTZMANN_EV = 8.617333262e-5
C_TR = 1.8793e20            # (2 pi m_u k / h^2)^{3/2}  [cm^-3 amu^-1.5 K^-1.5]
NS_DIVIDE = True

# Nuclear spin degeneracy (2I+1) of the main isotope
GNS = {'H': 2, 'He': 1, 'Li': 4, 'B': 4, 'C': 1, 'N': 3, 'O': 1,
       'F': 2, 'Na': 4, 'Mg': 1, 'Al': 6, 'Si': 1, 'P': 2, 'S': 1,
       'Cl': 4, 'K': 4, 'Ca': 1, 'Ti': 1, 'Fe': 1}


def _gns_total(atoms_sym):
    g = 1
    for s in atoms_sym:
        g *= GNS[s]
    return g


def assemble_log10_Kp_polyatomic(atoms_sym, D0_eV, Qmol, atomic, T):
    """
    log10 K_p [Pa^(N-1)] for M -> sum(atoms), on temperatures T.

    atoms_sym : list of element symbols (with repetition), e.g. ['H','H','O']
    D0_eV     : atomization energy [eV]
    Qmol      : callable T -> molecular partition function (ExoMol, HITRAN
                nuclear-spin convention; NaN outside its native range)
    atomic    : AtomicData (BC16 Table 8 partition functions, masses)
    """
    T = np.asarray(T, dtype=float)
    N = len(atoms_sym)

    q_mol = np.asarray(Qmol(T), dtype=float)
    if NS_DIVIDE:
        q_mol = q_mol / _gns_total(atoms_sym)

    m_mol = sum(atomic.mass(s) for s in atoms_sym)

    ln_kn = -D0_eV / (K_BOLTZMANN_EV * T)
    ln_kn = ln_kn + (N - 1) * np.log(C_TR) + 1.5 * (N - 1) * np.log(T)
    ln_kn = ln_kn - 1.5 * np.log(m_mol) - np.log(q_mol)
    for s in atoms_sym:
        ln_kn = ln_kn + 1.5 * np.log(atomic.mass(s)) \
                      + np.log(atomic.Q(s, T, 0))

    log10_kp = (ln_kn + (N - 1) * np.log(KB_PaCm3_per_K * T)) / np.log(10.0)
    return log10_kp


def kurucz_log10_K_polyatomic(entry, T):
    """Kurucz polynomial ln K -> log10 K (cgs association convention)."""
    from kurucz_molec import kurucz_lnKeq
    return kurucz_lnKeq(entry, T) / np.log(10.0)


def kurucz_to_bc16_pressure_convention(log10_K_cgs, T, N):
    """
    Kurucz cgs association K [cm^{3(N-1)}] -> BC16 dissociation K_p
    [Pa^(N-1)]:  log10 K_p = (N-1) log10(k_B T) - log10 K_cgs.
    """
    T = np.asarray(T, dtype=float)
    return (N - 1) * np.log10(KB_PaCm3_per_K * T) - log10_K_cgs
