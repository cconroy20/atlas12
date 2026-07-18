"""
build_molecules_physical_dat.py
================================

Generates molecules_physical.dat: same content as molecules_updated.dat
except that every "standard molecule" row now has its 5 polynomial
coefficients (E2..E6) replaced by the 5 physical-form coefficients
(a0..a4) describing the SAME equilibrium constant via:

    ln K_eq(T) = D0/(k_B T) - 1.5 * n_trans * ln T
               + a0 + a1/T + a2 ln T + a3 T + a4/T^2

The basis {1, 1/T, ln T, T, 1/T^2} replaces Kurucz's polynomial
basis {1, T, T^2, T^3, T^4}, with the same number of free parameters
(5) and the same total information content per row. The new form has
physically motivated terms (rotational/vibrational growth, electronic
excitation, vibrational equipartition, etc.) and is well-behaved at
high T (no T^4 blowup).

D0 (column 1 of the 6-coefficient block) is unchanged from the input.
n_trans is encoded implicitly in the molecule's atom count and is
computed by the F90 reader (NMOLEC) from the species code.

Fit method varies by species:
  - Diatomics matched to BC16:        fit to BC16 native grid (T >= 1000 K)
  - Polyatomics matched to ExoMol:    fit to ExoMol-assembled K_p
                                       on ExoMol native grid (T >= 1000 K)
  - All other true molecules:         fit to the Kurucz polynomial
                                       sampled at 1000-6000 K (self-fit)

Atom-only rows (E1=0, ionization bookkeeping) and auxiliary-equation
rows (atom-pool marker) pass through unchanged.

Usage: python build_molecules_physical_dat.py
"""

import sys
import numpy as np

from kurucz_molec import parse_molecules_dat, kurucz_lnKeq, ELEMENT_SYMBOLS
from bc16_loader import load_bc16, KB_PaCm3_per_K
from atomic_saha import AtomicData
from exomol_loader import ExoMolData
from polyatomic_d0 import POLYATOMIC_D0
from polyatomic_assembly import assemble_log10_Kp_polyatomic
from matcher import match_kurucz_to_bc16

K_BOLTZMANN_EV = 8.617333262e-5     # eV/K, matches kurucz_molec


# ---------- The three fitters --------------------------------------------

def fit_to_bc16(entry, bc16_mol):
    """Fit physical form to BC16 native grid (T >= 100 K)."""
    Tmask = bc16_mol.T_grid >= 100.0
    T_fit = bc16_mol.T_grid[Tmask]

    log10_K_kurucz = (np.log10(KB_PaCm3_per_K * T_fit)
                      - bc16_mol.log10_Kp[Tmask])
    lnK = log10_K_kurucz * np.log(10.0)

    return _do_fit(entry, T_fit, lnK)


def fit_to_exomol(entry, atoms, d0_info, atomic, exomol):
    """Fit physical form to ExoMol-assembled K_p on ExoMol's grid."""
    ex_name = d0_info['exomol_name']
    em_meta = exomol.molecules[ex_name]
    T_native = em_meta['T']
    Tmask = T_native >= 1000.0
    T_fit = T_native[Tmask]

    N = len(atoms)
    log10_Kp_em = assemble_log10_Kp_polyatomic(
        atoms, d0_info['D0'],
        lambda Tx: exomol.Q(ex_name, Tx), atomic, T_fit)

    # Convert from BC16 K_p convention to Kurucz cgs ln K_eq
    log10_K_kurucz = (N - 1) * np.log10(KB_PaCm3_per_K * T_fit) - log10_Kp_em
    lnK = log10_K_kurucz * np.log(10.0)

    finite = np.isfinite(lnK)
    return _do_fit(entry, T_fit[finite], lnK[finite])


def fit_to_kurucz_self(entry, T_lo=1000.0, T_hi=6000.0, n_pts=50):
    """Fit physical form to Kurucz's own polynomial K_eq."""
    T_fit = np.linspace(T_lo, T_hi, n_pts)
    lnK = kurucz_lnKeq(entry, T_fit)
    return _do_fit(entry, T_fit, lnK)


def _do_fit(entry, T_fit, lnK):
    """Common LSQ kernel: fit physical-form a0..a4 against lnK on T_fit."""
    D0 = entry['E'][0]
    n_trans = entry['n_trans']
    fixed = D0 / (K_BOLTZMANN_EV * T_fit) - 1.5 * n_trans * np.log(T_fit)
    R = lnK - fixed
    A = np.column_stack([
        np.ones_like(T_fit),
        1.0 / T_fit,
        np.log(T_fit),
        T_fit,
        1.0 / T_fit**2,
    ])
    coeffs, *_ = np.linalg.lstsq(A, R, rcond=None)
    a0, a1, a2, a3, a4 = coeffs
    return D0, a0, a1, a2, a3, a4


# ---------- Header text --------------------------------------------------

NEW_HEADER = """! ===========================================================================
! molecules_physical.dat -- Kurucz molecules.dat with physical-form fits
! ---------------------------------------------------------------------------
! Format (90 columns):
!   cols 1-10  (A10):    species label (informational; ignored on read)
!   cols 11-22 (F12.2):  Kurucz species code (encodes atoms + ion stage)
!   col  23    (1X):     space
!   cols 24-30 (F7.3):   D0 [eV]  (atomization energy; signed)
!   cols 31-90 (5E12.4): physical-form coefficients a0, a1, a2, a3, a4
!   cols 91+:            trailing comment with provenance (ignored on read)
!
! Equilibrium-constant evaluation:
!   ln K_eq(T) = D0/(k_B T) - 1.5 * n_trans * ln T
!              + a0 + a1/T + a2 ln T + a3 T + a4/T^2
!
! where n_trans is determined by the species code (number of dissociation
! products minus one) and k_B is in eV/K.  The basis {1, 1/T, ln T, T, 1/T^2}
! replaces Kurucz's earlier polynomial basis {1, T, T^2, T^3, T^4} with the
! same number of free parameters (5) but with physically motivated terms:
! constant, electronic-excitation 1/T term, rotational/vibrational ln T
! term, vibrational-equipartition T term, and electronic-correction 1/T^2.
!
! Fit methods (per row, identified in trailing comment):
!   "BC16 fit"         - diatomic refit to BC16 Table-1 native grid (T>=100 K).
!   "ExoMol fit"       - polyatomic refit to ExoMol-assembled K_p on ExoMol's
!                        native partition-function grid (T>=1000 K).
!   "Kurucz poly fit"  - true molecule with no external reference (BC16 has
!                        no entry, no ExoMol .pf available).  Fitted to
!                        Kurucz's original polynomial K_eq sampled in
!                        1000-6000 K, so D0 and global behavior are
!                        preserved while replacing the polynomial basis
!                        with the well-behaved physical basis.
!
! D0 source (after "; ") is the same as in the previous molecules_updated.dat:
!   BC16     = Barklem & Collet 2016 (A&A 588, A96), Table 1.
!   ATcT     = Active Thermochemical Tables (Ruscic et al., Argonne).
!              Squib codes like '2006Rus/Pin' identify primary sources
!              within the CCCBDB chain.
!   ATcT1.220 = ATcT version 1.220 (December 2025), with silicon compounds
!               added (Bross & Ruscic, to be published).  Used for SiO2,
!               SiH4, C3, and CH3.
!   CCCBDB   = NIST Computational Chemistry Comparison and Benchmark
!              Database (cccbdb.nist.gov/hf0kx.asp).
!   CODATA   = CODATA atomic enthalpies of formation at 0 K.
!   Gurvich  = Gurvich et al. Thermodynamic Properties of Individual
!              Substances (1989-1996).
!   JANAF    = NIST-JANAF Thermochemical Tables (Chase 1998).
!   Koput2023 = Koput, J. Mol. Spectrosc. 395, 111805 (2023).  State-of-the-
!               art CCSD(T) ab initio with septuple-zeta basis (used for MgOH).
!   ST84     = Sauval & Tatum (1984) - source of many original Kurucz values.
!
! Atom rows (E1 = 0.0) pass through unchanged from molecules_updated.dat
! -- they have no polynomial to refit (E2..E6 are all zero placeholders).
! Every row with non-zero K_eq coefficients has been refit to the
! physical form using one of the methods above.
! ===========================================================================
"""


# ---------- Main build ---------------------------------------------------

def main():
    bc16   = load_bc16()
    atomic = AtomicData()
    exomol = ExoMolData('.')

    entries = parse_molecules_dat('molecules_updated.dat')
    matches = match_kurucz_to_bc16(entries, bc16)
    matched_bc16_name_by_code = {
        e['code']: name for e, name, _ in matches if name is not None}

    # Read raw lines so we can preserve per-row formatting / comments.
    # entry['line'] is the 1-indexed line number in the source file.
    with open('molecules_updated.dat') as f:
        raw_lines = f.readlines()

    n_bc16 = n_em = n_self = n_passthru = 0

    out_lines = []
    out_lines.append(NEW_HEADER.rstrip('\n'))

    for e in entries:
        line_idx = e['line'] - 1   # convert to 0-indexed
        line_orig = raw_lines[line_idx]
        # Strip the original trailing comment (after the first '!') if any
        bang = line_orig.find('!')
        body = line_orig[:bang].rstrip() if bang >= 0 else line_orig.rstrip()

        # Atom-only rows (E1=0): pass through unchanged.  These hold no
        # equilibrium constant -- their K_eq is computed via Saha rather
        # than from polynomial coefficients, and slots E2..E6 are zero
        # placeholders.
        if e['E'][0] == 0.0:
            out_lines.append(line_orig.rstrip('\n'))
            n_passthru += 1
            continue

        # Any row with a non-zero E1 has K_eq coefficients that the F90
        # evaluator will USE (the runtime check is `IF (EQUIL(1) .NE. 0)`),
        # so it must be refit to the physical form.  This includes:
        #  - normal multi-atom molecules
        #  - anion species (e.g. H-, OH-, C2-, ...) which have an atom-pool
        #    marker but are real species with polynomial K_eq
        #  - atom + electron equilibria (e.g. H-, C-, O-, F-, ...) which
        #    parse as "not is_true_molecule" but still carry polynomial
        #    coefficients used by the runtime evaluator

        # Standard molecule: pick a fitter
        atoms_sym = [ELEMENT_SYMBOLS[a] for a in e['real_atoms']]
        N = e['n_real_atoms']

        method = None
        coeffs = None

        # 1. Diatomic with BC16 match
        if N == 2 and e['code'] in matched_bc16_name_by_code:
            bc_name = matched_bc16_name_by_code[e['code']]
            try:
                coeffs = fit_to_bc16(e, bc16[bc_name])
                method = 'BC16 fit'
                n_bc16 += 1
            except Exception as exc:
                print(f'  WARN: BC16 fit failed for {e["formula"]}: {exc}',
                      file=sys.stderr)
                coeffs = None

        # 2. Polyatomic with ExoMol match
        if coeffs is None and N >= 3 and e['ion'] == 0:
            for ex_name, d0_info in POLYATOMIC_D0.items():
                if (sorted(d0_info['atoms']) == sorted(atoms_sym)
                        and exomol.has(ex_name)):
                    try:
                        coeffs = fit_to_exomol(e, atoms_sym, d0_info,
                                                atomic, exomol)
                        method = 'ExoMol fit'
                        n_em += 1
                    except Exception as exc:
                        print(f'  WARN: ExoMol fit failed for {e["formula"]}: '
                              f'{exc}', file=sys.stderr)
                    break

        # 3. Fallback: self-fit to Kurucz polynomial
        if coeffs is None:
            coeffs = fit_to_kurucz_self(e)
            method = 'Kurucz poly fit'
            n_self += 1

        # Encode the new row.  We need to preserve cols 1-22 (label+code)
        # and col 23 (space), then write D0 (cols 24-30, F7.3) and
        # a0..a4 (cols 31-90, 5E12.4).
        # The original line already has cols 1-23 in the right format; we
        # rebuild cols 24-90.
        D0, a0, a1, a2, a3, a4 = coeffs
        prefix = body[:23]  # cols 1-23

        # D0 in F7.3
        d0_str = f'{D0:7.3f}'
        if len(d0_str) > 7:
            d0_str = f'{D0:7.2f}'

        # 5 x E12.4 for the coefficients
        coef_str = ''.join(f'{c:12.4E}' for c in (a0, a1, a2, a3, a4))

        new_body = prefix + d0_str + coef_str
        # Pad/truncate to width 90
        new_body = f'{new_body:<90s}'

        # Build new comment (preserve any original D0 provenance)
        # The original line had a comment like:
        #    ! BC16 Table 1; was 4.392
        # or  ! Kurucz original
        # or  ! ATcT1.220; was 13.957
        original_comment = ''
        if bang >= 0:
            original_comment = line_orig[bang+1:].rstrip().lstrip()
        # Replace the leading "BC16 Table 1" / "Kurucz original" / etc.
        # with the new method tag, preserving the "; was X.XXX" or
        # "; D0 from <source>" pieces if present.
        if 'was' in original_comment:
            # extract "was X.XXX" piece
            was_idx = original_comment.find('was')
            d0_provenance = original_comment[:was_idx].rstrip(' ;')
            was_piece = original_comment[was_idx:]
            new_comment = f'{method}; D0 = {d0_provenance}; {was_piece}'
        elif 'Kurucz original' in original_comment:
            new_comment = f'{method}; D0 = Kurucz original'
        else:
            # Fall back on the original wording
            new_comment = (f'{method}; D0 = {original_comment}'
                           if original_comment else method)

        out_lines.append(f'{new_body}   ! {new_comment}')

    out_path = 'molecules_physical.dat'
    with open(out_path, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')

    print(f'Wrote {out_path}')
    print(f'  Atom-only rows passed through:           {n_passthru}')
    print(f'  Diatomic refit to BC16:                  {n_bc16}')
    print(f'  Polyatomic refit to ExoMol:              {n_em}')
    print(f'  Fallback self-fit to Kurucz polynomial:  {n_self}')
    print(f'  Total rows:                              '
          f'{n_passthru + n_bc16 + n_em + n_self}')


if __name__ == '__main__':
    main()
