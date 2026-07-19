"""
Unified Kurucz-vs-modern comparison driver.

Produces a single PDF report covering:
  - 49 diatomics matched to BC16 (41 neutrals + 8 ions Saha-corrected)
  - polyatomics matched to ExoMol .pf files (currently 13)

Each per-molecule page shows up to 4 curves on the same axis:
  1. Kurucz original (from molecules.dat)
  2. Kurucz updated  (from molecules_updated.dat - same polynomial,
                      new E1)
  3. BC16            (diatomics: tabulated K_p; polyatomics: not shown)
  4. ExoMol assembled (from ExoMol Q + BC16 atomic Q + modern D_0,
                       in BC16 pressure convention).

Convention: ALL curves on a single page are converted to a common
axis. For diatomics, that axis is BC16 cgs (n_AB / (n_A * n_B), units
cm^3). For polyatomics, that axis is BC16 pressure (K_p in Pa^(N-1)).

Sort order:
  cover page,
  then diatomic neutrals (by descending max|Δ Kurucz original - BC16|),
  then diatomic ions (same),
  then polyatomics (by descending max|Δ Kurucz original - reference|).

Outputs:
  kurucz_unified_report.pdf       all per-molecule pages
  kurucz_unified_summary.csv      one row per molecule
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import CubicSpline

from kurucz_molec import parse_molecules_dat, kurucz_lnKeq, ELEMENT_SYMBOLS
from bc16_loader import load_bc16, KB_PaCm3_per_K
from matcher import match_kurucz_to_bc16
from atomic_saha import AtomicData
from exomol_loader import (
    ExoMolData, assemble_log10_Kp_exomol)
from polyatomic_assembly import (
    assemble_log10_Kp_polyatomic, kurucz_log10_K_polyatomic,
    kurucz_to_bc16_pressure_convention)
from polyatomic_d0 import POLYATOMIC_D0


# T grids
T_DIATOMIC  = np.logspace(np.log10(1000.0), np.log10(20000.0), 250)
T_POLY      = np.logspace(np.log10(1000.0), np.log10(20000.0), 250)
T_STAT_LO_D, T_STAT_HI_D = 1500.0, 8000.0
T_STAT_LO_P, T_STAT_HI_P = 1500.0, 6000.0
K_BOLTZMANN_EV = 8.617333262e-5     # eV/K, matches kurucz_molec


# ============================================================
# Physical-form fit (Option C-lite)
#
# The functional form is:
#   ln K_eq(T) = D0/(k_B T) - 1.5 * n_trans * ln T
#              + a0 + a1/T + a2 ln T + a3 T + a4/T^2
#
# This matches Kurucz's polynomial form's count (5 free params: a0..a4)
# and sign convention exactly, but uses {1, 1/T, ln T, T, 1/T^2} as the
# basis instead of {1, T, T^2, T^3, T^4}. Each new basis function has a
# physical interpretation (translational + electronic excitation +
# rotational/vibrational growth + vibrational equipartition + electronic
# excited-state correction).
# ============================================================

def fit_physical_diatomic(entry, bc16_mol):
    """
    Fit the physical form to BC16 data (sampled at BC16's actual grid).

    Returns:
        coeffs : tuple (D0, a0, a1, a2, a3, a4)  -- same length as Kurucz E
        T_fit  : array of temperatures used in fit (BC16's actual grid)
        lnK_bc : reference ln K at those temperatures (Kurucz units)
    """
    # Use BC16's actual temperature grid (no spline interpolation). Restrict
    # to T in (0, T_max] to avoid the unphysical T<<1 K entries BC16 carries.
    Tmask = bc16_mol.T_grid >= 100.0
    T_fit = bc16_mol.T_grid[Tmask]

    # Convert BC16 log10 K_p (Pa) to Kurucz's natural log K_eq (cgs)
    log10_K_kurucz = (np.log10(KB_PaCm3_per_K * T_fit)
                      - bc16_mol.log10_Kp[Tmask])
    lnK_bc = log10_K_kurucz * np.log(10.0)

    # Fixed pieces: D0 from Kurucz entry (which has been pinned to BC16)
    D0 = entry['E'][0]
    n_trans = entry['n_trans']
    fixed = D0 / (K_BOLTZMANN_EV * T_fit) - 1.5 * n_trans * np.log(T_fit)

    # Residual to be matched by physical-basis terms
    R = lnK_bc - fixed

    # Linear LSQ design matrix
    A = np.column_stack([
        np.ones_like(T_fit),    # a0   (constant)
        1.0 / T_fit,            # a1   (1/T)
        np.log(T_fit),          # a2   (ln T)
        T_fit,                  # a3   (T)
        1.0 / T_fit**2,         # a4   (1/T^2)
    ])

    coeffs, *_ = np.linalg.lstsq(A, R, rcond=None)
    a0, a1, a2, a3, a4 = coeffs

    return (D0, a0, a1, a2, a3, a4), T_fit, lnK_bc


def fit_physical_polyatomic(entry, atoms, d0_info, atomic, exomol):
    """
    Fit the physical form to ExoMol-assembled K_p for a polyatomic.

    Same form as for diatomics, with n_trans coming from the entry. We
    fit in Kurucz's natural ln K_eq (cgs) convention so the resulting
    coefficients can later be evaluated in the same way as the existing
    Kurucz polynomial.

    Returns:
        coeffs : tuple (D0, a0, a1, a2, a3, a4)
        T_fit  : ExoMol's native T grid restricted to T >= 1000 K
        lnK_em : reference ln K at those temperatures (Kurucz units)
    """
    ex_name = d0_info['exomol_name']
    em_meta = exomol.molecules[ex_name]

    # Native ExoMol grid, restricted to T >= 1000 K (stellar-relevant)
    T_native = em_meta['T']
    Tmask = T_native >= 1000.0
    T_fit = T_native[Tmask]

    # Assemble log10 K_p in BC16 pressure convention at these T points,
    # using the verified D0 from d0_info.
    N = len(atoms)
    log10_Kp_em = assemble_log10_Kp_polyatomic(
        atoms, d0_info['D0'],
        lambda Tx: exomol.Q(ex_name, Tx), atomic, T_fit)

    # Convert from BC16 pressure convention (Pa^(N-1)) to Kurucz's natural
    # log K_eq (cgs, cm^{3*(N-1)}). The conversion factor matches the
    # one used to send Kurucz cgs -> BC16 Pa: K_p [Pa^(N-1)] =
    # (k_B T)^(N-1) / K_cgs [cm^{3(N-1)}], so
    # log10 K_cgs = (N-1)*log10(kB*T) - log10 K_p.
    log10_K_kurucz = (N - 1) * np.log10(KB_PaCm3_per_K * T_fit) - log10_Kp_em
    lnK_em = log10_K_kurucz * np.log(10.0)

    # Drop any NaN points (shouldn't occur within ExoMol range, but guard)
    finite = np.isfinite(lnK_em)
    T_fit = T_fit[finite]
    lnK_em = lnK_em[finite]

    # Fixed pieces
    D0 = entry['E'][0]
    n_trans = entry['n_trans']
    fixed = D0 / (K_BOLTZMANN_EV * T_fit) - 1.5 * n_trans * np.log(T_fit)

    # Residual to be matched by physical-basis terms
    R = lnK_em - fixed

    # Linear LSQ design matrix
    A = np.column_stack([
        np.ones_like(T_fit),
        1.0 / T_fit,
        np.log(T_fit),
        T_fit,
        1.0 / T_fit**2,
    ])
    coeffs, *_ = np.linalg.lstsq(A, R, rcond=None)
    a0, a1, a2, a3, a4 = coeffs

    return (D0, a0, a1, a2, a3, a4), T_fit, lnK_em


def evaluate_physical(coeffs, n_trans, T):
    """Evaluate the physical form ln K_eq at temperatures T."""
    D0, a0, a1, a2, a3, a4 = coeffs
    return (D0 / (K_BOLTZMANN_EV * T) - 1.5 * n_trans * np.log(T)
            + a0 + a1 / T + a2 * np.log(T) + a3 * T + a4 / T**2)


# Canonical display names for species whose internal formula string uses
# the code-ordered (ascending-Z) convention rather than chemical convention.
_DISPLAY_OVERRIDE = {
    'O2Ti': 'TiO2', 'O2V': 'VO2', 'O2Si': 'SiO2', 'O2S': 'SO2',
    'HOAl': 'AlOH', 'OAl2': 'Al2O', 'HO2Al': 'AlO2H', 'HOMg': 'MgOH',
    'HOCa': 'CaOH', 'HONa': 'NaOH', 'O2H': 'HO2', 'H2C': 'CH2',
    'H2N': 'NH2', 'C2Si': 'SiC2', 'CSi2': 'Si2C', 'COS': 'OCS',
    # the generator's H-first branch mishandles these (ExoMol convention)
    'H3C': 'CH3', 'H3N': 'NH3', 'H2C2': 'C2H2', 'H4C': 'CH4',
    'H4Si': 'SiH4', 'H3P': 'PH3', 'H2Si': 'SiH2', 'HLiO': 'LiOH',
    'HOK': 'KOH',
}


def display_formula(name):
    """'O2Ti' -> 'TiO$_2$', 'H3+' -> 'H$_3^+$': chemical-convention name
    with mathtext subscripts for counts and superscripts for charge."""
    import re as _re
    name = _DISPLAY_OVERRIDE.get(name, name)
    m = _re.match(r'^(.*?)([+-]+)$', name)
    base, charge = (m.group(1), m.group(2)) if m else (name, '')
    out = _re.sub(r'(\d+)', r'$_\1$', base)
    if charge:
        if out.endswith('$'):                    # merge H$_3$ + ^+ -> H$_3^+$
            out = out[:-1] + f'^{charge}$'
        else:
            out += f'$^{charge}$'
    return out


# Diatomic ion channels (same as elsewhere)
ION_DISSOCIATION_CHANNEL = {
    'H2+':  ('H', 'H'),
    'CH+':  ('H', 'C'),
    'OH+':  ('O', 'H'),
    'SiH+': ('H', 'Si'),
    'CN+':  ('N', 'C'),
    'N2+':  ('N', 'N'),
    'NO+':  ('N', 'O'),
    'O2+':  ('O', 'O'),
}


# ============================================================
# Per-molecule comparison computation
# ============================================================

def compute_diatomic(entry_orig, entry_upd, bc16_mol, atomic, exomol):
    """
    Build all comparison curves for a diatomic.

    Returns a dict containing T, log10_K curves for Kurucz original,
    Kurucz updated, BC16, and ExoMol-assembled (when applicable),
    plus residuals and stats.

    All log10 K's are in BC16 cgs convention (units cm^3) so they
    can be plotted on a single axis.
    """
    T = T_DIATOMIC
    # entry_orig is None for species added after the April 2026 refit:
    # they have no Kurucz polynomial ancestry, so the orig/updated curves
    # (and their residuals) are absent and only the BC16 reference plus
    # the physical-form fit are shown.
    is_new = entry_orig is None
    base = entry_upd if is_new else entry_orig
    ion = base['ion']
    formula = base['formula']

    # --- Kurucz original and updated (cgs from polynomial)
    if is_new:
        log10_K_kur = None
        log10_K_upd = None
    else:
        log10_K_kur = kurucz_lnKeq(entry_orig, T) / np.log(10.0)
        log10_K_upd = kurucz_lnKeq(entry_upd,  T) / np.log(10.0)

    # --- BC16 reference (convert from K_p [Pa] to cgs).
    # BC16 data extend only to T_max = 10000 K; mask the curve outside
    # that range so the plotted line reflects actual data, not spline
    # extrapolation. Below 1000 K the polyatomic K_p values get extreme
    # but BC16 does not always tabulate them; mask similarly.
    cs = CubicSpline(np.log10(bc16_mol.T_grid), bc16_mol.log10_Kp)
    logKp_bc = cs(np.log10(T))
    log10_K_bc = np.log10(KB_PaCm3_per_K * T) - logKp_bc
    bc16_T_max = float(bc16_mol.T_grid.max())
    bc16_T_min = max(1000.0, float(bc16_mol.T_grid[bc16_mol.T_grid > 0].min()))
    bc16_mask = (T >= bc16_T_min) & (T <= bc16_T_max)
    log10_K_bc = np.where(bc16_mask, log10_K_bc, np.nan)

    # --- BC16 raw data points (no interpolation) for the scatter overlay.
    # Restrict to T >= 1000 K so the symbols line up with the plotted axis;
    # BC16 carries unphysical T<<1 K entries we don't want to display.
    bc16_T_pts = bc16_mol.T_grid[bc16_mol.T_grid >= 1000.0]
    bc16_log10K_pts = (np.log10(KB_PaCm3_per_K * bc16_T_pts)
                       - bc16_mol.log10_Kp[bc16_mol.T_grid >= 1000.0])

    # --- Physical-form fit (refit a0..a4 against BC16 with D0 fixed)
    phys_coeffs, _, _ = fit_physical_diatomic(entry_upd, bc16_mol)
    log10_K_phys = evaluate_physical(phys_coeffs, entry_upd['n_trans'], T) / np.log(10.0)

    # --- Saha-correction for ions: convert Kurucz cgs to BC16 cgs
    # NOTE: this applies ONLY to the polynomial-form (Kurucz original/updated)
    # curves. The physical-form fit is fitted to BC16 directly and is
    # therefore already on BC16's axis, so no Saha correction is needed
    # for it. The raw BC16 points are also already on BC16's axis.
    saha_log10 = None
    channel = None
    if ion == 1 and formula in ION_DISSOCIATION_CHANNEL and not is_new:
        neutral_atom, ionized_atom = ION_DISSOCIATION_CHANNEL[formula]
        saha_log10 = atomic.log10_K_Saha(ionized_atom, T)
        log10_K_kur = log10_K_kur - saha_log10
        log10_K_upd = log10_K_upd - saha_log10
        # log10_K_phys: NO Saha subtraction (already on BC16 axis)
        channel = (neutral_atom, ionized_atom)

    # --- ExoMol-assembled (neutrals only)
    log10_K_em = None
    if ion == 0 and exomol.has(bc16_mol.name):
        try:
            log10_Kp_em = assemble_log10_Kp_exomol(
                bc16_mol.atom_a, bc16_mol.atom_b, bc16_mol.D0_adopted,
                lambda Tx: exomol.Q(bc16_mol.name, Tx), atomic, T)
            log10_K_em = np.log10(KB_PaCm3_per_K * T) - log10_Kp_em
        except (KeyError, ValueError):
            pass

    # --- Residuals (vs BC16, the diatomic reference)
    delta_phys = log10_K_phys - log10_K_bc
    mask = (T >= T_STAT_LO_D) & (T <= T_STAT_HI_D)
    rms_phys = float(np.sqrt(np.nanmean(delta_phys[mask]**2)))
    max_phys = float(np.nanmax(np.abs(delta_phys[mask])))
    if is_new:
        delta_kur = None
        delta_upd = None
        rms_kur = max_kur = rms_upd = max_upd = float('nan')
    else:
        delta_kur = log10_K_kur - log10_K_bc
        delta_upd = log10_K_upd - log10_K_bc
        rms_kur  = float(np.sqrt(np.nanmean(delta_kur[mask]**2)))
        max_kur  = float(np.nanmax(np.abs(delta_kur[mask])))
        rms_upd  = float(np.sqrt(np.nanmean(delta_upd[mask]**2)))
        max_upd  = float(np.nanmax(np.abs(delta_upd[mask])))

    return {
        'kind':        'diatomic',
        'is_new':      is_new,
        'T':           T,
        'log10_K_kurucz_orig':    log10_K_kur,
        'log10_K_kurucz_updated': log10_K_upd,
        'log10_K_bc16':           log10_K_bc,
        'log10_K_exomol':         log10_K_em,
        'log10_K_physical':       log10_K_phys,
        'phys_coeffs':            phys_coeffs,
        'bc16_T_pts':             bc16_T_pts,
        'bc16_log10K_pts':        bc16_log10K_pts,
        'delta_kurucz':           delta_kur,
        'delta_updated':          delta_upd,
        'delta_physical':         delta_phys,
        'rms_kurucz':  rms_kur,
        'max_kurucz':  max_kur,
        'rms_updated': rms_upd,
        'max_updated': max_upd,
        'rms_physical': rms_phys,
        'max_physical': max_phys,
        'channel':     channel,
        'ion':         ion,
        'D0_kurucz_orig': float('nan') if is_new else entry_orig['E'][0],
        'D0_kurucz_upd':  entry_upd['E'][0],
        'D0_BC16':        bc16_mol.D0_adopted,
        'bc16_name':      bc16_mol.name,
        'formula':        formula,
        'Tmax_em':        exomol.molecules[bc16_mol.name]['Tmax']
                          if exomol.has(bc16_mol.name) else None,
    }


def compute_polyatomic_self_fit(entry_orig, entry_upd, atoms_sym,
                                T_lo=1000.0, T_hi=6000.0, janaf_rec=None,
                                cur_entry=None):
    """
    Polyatomics where we have no external reference (BC16 doesn't tabulate,
    no ExoMol .pf available, or no POLYATOMIC_D0 entry). The "data" here
    is Kurucz's own polynomial K_eq, sampled in the 1000-6000 K range.

    The result has the same dict structure as compute_polyatomic so that
    _plot_one can render it directly. There is no separate updated curve,
    no log10_K_bc16, no log10_K_exomol -- just Kurucz original (==updated)
    plus the physical-form fit and the data points used as the fit target.
    """
    T = T_POLY
    N = entry_orig['n_real_atoms']

    # Sample Kurucz polynomial on the fit grid (Kurucz's natural cgs units)
    n_pts = 50
    T_fit = np.linspace(T_lo, T_hi, n_pts)
    lnK_kur_fit = kurucz_lnKeq(entry_orig, T_fit)

    # Fit physical form. D0 fixed at entry's E1 (Kurucz's own atomization).
    D0 = entry_orig['E'][0]
    n_trans = entry_orig['n_trans']
    fixed_fit = D0 / (K_BOLTZMANN_EV * T_fit) - 1.5 * n_trans * np.log(T_fit)
    R = lnK_kur_fit - fixed_fit
    A = np.column_stack([
        np.ones_like(T_fit),
        1.0 / T_fit,
        np.log(T_fit),
        T_fit,
        1.0 / T_fit**2,
    ])
    coeffs_a, *_ = np.linalg.lstsq(A, R, rcond=None)
    a0, a1, a2, a3, a4 = coeffs_a
    phys_coeffs = (D0, a0, a1, a2, a3, a4)

    # Evaluate Kurucz and physical fit on the dense plot grid.
    # Convert log10_K (cgs) -> log10 K_p (Pa^(N-1)) for plotting (matches
    # the rest of the polyatomic panels' axis convention).
    log10_K_kur_cgs = kurucz_lnKeq(entry_orig, T) / np.log(10.0)
    log10_K_phys_cgs = (evaluate_physical(phys_coeffs, n_trans, T)
                        / np.log(10.0))
    log10_Kp_kur  = kurucz_to_bc16_pressure_convention(log10_K_kur_cgs, T, N)
    log10_Kp_phys = kurucz_to_bc16_pressure_convention(log10_K_phys_cgs, T, N)

    # The "data" points in the same plot units
    log10_K_fit_cgs = lnK_kur_fit / np.log(10.0)
    log10_Kp_fit    = kurucz_to_bc16_pressure_convention(
        log10_K_fit_cgs, T_fit, N)

    # Residuals are physical-vs-Kurucz on the dense plot grid.
    delta_phys = log10_Kp_phys - log10_Kp_kur

    # Stats: residuals over the FIT range only (where the curves should
    # agree by construction). Outside the fit range the physical form
    # diverges from the polynomial, but that divergence is the *point*
    # (it represents the well-behaved extrapolation).
    mask = (T >= T_lo) & (T <= T_hi)
    if mask.any():
        rms_phys = float(np.sqrt(np.mean(delta_phys[mask]**2)))
        max_phys = float(np.max(np.abs(delta_phys[mask])))
    else:
        rms_phys = max_phys = float('nan')

    # --- Optional NIST-JANAF reference (GGchem dispol fit) ---
    # For self-fit species with no ExoMol coverage, the Stock/Kitzmann
    # kp(T) compilation (traceable to NIST-JANAF/Burcat) provides the
    # modern external reference.  Overlaid for assessment; the physical
    # fit itself remains anchored to the Kurucz polynomial.
    log10_Kp_janaf = None
    janaf_stats = {}
    if janaf_rec is not None:
        from ggchem_loader import ln_K_kurucz as janaf_lnK
        log10_K_janaf_cgs = janaf_lnK(janaf_rec, T) / np.log(10.0)
        log10_Kp_janaf = kurucz_to_bc16_pressure_convention(
            log10_K_janaf_cgs, T, N)
        jmask = (T >= 1500.0) & (T <= T_hi)
        d_kur_j = log10_Kp_kur - log10_Kp_janaf
        d_phys_j = log10_Kp_phys - log10_Kp_janaf
        janaf_stats = {
            'rms_kurucz_janaf': float(np.sqrt(np.mean(d_kur_j[jmask]**2))),
            'max_kurucz_janaf': float(np.max(np.abs(d_kur_j[jmask]))),
            'rms_physical_janaf': float(np.sqrt(np.mean(d_phys_j[jmask]**2))),
            'max_physical_janaf': float(np.max(np.abs(d_phys_j[jmask]))),
            'delta_kurucz_janaf': d_kur_j,
        }

    # --- Overlay of the CURRENT data/molecules.dat row when it has been
    # refit since April (e.g. the MgOH/HO2 D0-swap refits): the page's
    # "physical form fit" reproduces the April state, so a materially
    # different filed row is shown as its own curve, with its residual
    # vs the JANAF reference when one is present.
    log10_Kp_cur = None
    cur_stats = {}
    if cur_entry is not None and cur_entry['E'][0] != 0.0:
        lnK_cur = (cur_entry['E'][0] / (K_BOLTZMANN_EV * T)
                   - 1.5 * cur_entry['n_trans'] * np.log(T)
                   + cur_entry['E'][1] + cur_entry['E'][2] / T
                   + cur_entry['E'][3] * np.log(T) + cur_entry['E'][4] * T
                   + cur_entry['E'][5] / T**2)
        kp_cur = kurucz_to_bc16_pressure_convention(
            lnK_cur / np.log(10.0), T, N)
        wmask = (T >= T_lo) & (T <= T_hi)
        if np.max(np.abs((kp_cur - log10_Kp_phys)[wmask])) > 0.02:
            log10_Kp_cur = kp_cur
            if log10_Kp_janaf is not None:
                d_cur_j = kp_cur - log10_Kp_janaf
                jm = (T >= 1500.0) & (T <= T_hi)
                cur_stats = {
                    'rms_current_janaf':
                        float(np.sqrt(np.mean(d_cur_j[jm]**2))),
                    'max_current_janaf':
                        float(np.max(np.abs(d_cur_j[jm]))),
                    'delta_current_janaf': d_cur_j,
                    'D0_current': cur_entry['E'][0],
                }

    return {
        'kind':        'polyatomic',
        'self_fit':    True,    # marker so _plot_one knows to label differently
        'T':           T,
        'N':           N,
        'log10_K_kurucz_orig':    log10_Kp_kur,
        'log10_K_kurucz_updated': log10_Kp_kur,   # no E1 update for these
        'log10_K_exomol':         None,
        'log10_K_physical':       log10_Kp_phys,
        'log10_K_bc16':           None,
        'log10_K_janaf':          log10_Kp_janaf,
        'log10_K_current':        log10_Kp_cur,
        'phys_coeffs':            phys_coeffs,
        **cur_stats,
        # Fit-target scatter points
        'em_T_pts':               T_fit,
        'em_log10K_pts':          log10_Kp_fit,
        'delta_kurucz':           np.zeros_like(T),    # by definition
        'delta_updated':          np.zeros_like(T),
        'delta_physical':         delta_phys,
        'rms_kurucz':  0.0,
        'max_kurucz':  0.0,
        'rms_updated': 0.0,
        'max_updated': 0.0,
        'rms_physical': rms_phys,
        'max_physical': max_phys,
        **janaf_stats,
        'ion':         entry_orig['ion'],
        'channel':     None,
        'D0_kurucz_orig': entry_orig['E'][0],
        'D0_kurucz_upd':  entry_upd['E'][0],
        'D0_BC16':        float('nan'),    # no external reference
        'D0_uncert':      float('nan'),
        'D0_source':      '',
        'D0_verified':    False,
        'formula':        entry_orig['formula'],
        'bc16_name':      entry_orig['formula'],
        'atoms':          atoms_sym,
        'Tmax_em':        None,
        'fit_T_lo':       T_lo,
        'fit_T_hi':       T_hi,
    }


def compute_polyatomic_new(entry_cur, janaf_rec, em_ref=None):
    """
    Page for a polyatomic added after the April refit: no Kurucz
    ancestry, so the curves are the CURRENT filed row and the reference
    it was fit to -- the ExoMol assembly when em_ref = (d0_info, exomol,
    atomic) is given (with any JANAF/SK curve kept as an overlay), else
    the NIST-JANAF (Stock/Kitzmann) fit.
    """
    T = T_POLY
    N = entry_cur['n_real_atoms']
    from ggchem_loader import ln_K_kurucz as janaf_lnK
    log10_Kp_janaf = None
    if janaf_rec is not None:
        log10_Kp_janaf = kurucz_to_bc16_pressure_convention(
            janaf_lnK(janaf_rec, T) / np.log(10.0), T, N)
    E = entry_cur['E']
    lnK_cur = (E[0] / (K_BOLTZMANN_EV * T)
               - 1.5 * entry_cur['n_trans'] * np.log(T)
               + E[1] + E[2] / T + E[3] * np.log(T) + E[4] * T + E[5] / T**2)
    log10_Kp_cur = kurucz_to_bc16_pressure_convention(
        lnK_cur / np.log(10.0), T, N)

    if em_ref is not None:
        # ExoMol-assembled reference on its native range
        d0i, exomol, atomic = em_ref
        name = d0i['exomol_name']
        ref_curve = assemble_log10_Kp_polyatomic(
            d0i['atoms'], d0i['D0'],
            lambda x: exomol.Q(name, x), atomic, T)
        Tn = exomol.molecules[name]['T']
        T_pts = Tn[(Tn >= 1000.0)][::max(1, len(Tn) // 25)]
        pts = assemble_log10_Kp_polyatomic(
            d0i['atoms'], d0i['D0'],
            lambda x: exomol.Q(name, x), atomic, T_pts)
        ref_label, pts_label = 'ExoMol', 'ExoMol data'
        tmax_em = exomol.molecules[name]['Tmax']
    else:
        ref_curve = log10_Kp_janaf
        T_pts = np.arange(100.0, 6001.0, 200.0)
        pts = kurucz_to_bc16_pressure_convention(
            janaf_lnK(janaf_rec, T_pts) / np.log(10.0), T_pts, N)
        ref_label, pts_label = 'JANAF/SK', 'JANAF/SK fit points'
        tmax_em = None

    delta = log10_Kp_cur - ref_curve
    mask = ((T >= T_STAT_LO_P) & (T <= T_STAT_HI_P)
            & np.isfinite(delta))
    return {
        'kind': 'polyatomic', 'is_new': True, 'self_fit': False,
        'ref_label': ref_label, 'ref_points_label': pts_label,
        'T': T, 'N': N,
        'log10_K_kurucz_orig': None, 'log10_K_kurucz_updated': None,
        'log10_K_exomol': ref_curve if em_ref is not None else None,
        'log10_K_bc16': None,
        'log10_K_janaf': log10_Kp_janaf, 'log10_K_physical': log10_Kp_cur,
        'phys_coeffs': (E[0],) + tuple(E[1:]),
        'em_T_pts': T_pts, 'em_log10K_pts': pts,
        'delta_kurucz': None, 'delta_updated': None,
        'delta_physical': delta,
        'rms_kurucz': float('nan'), 'max_kurucz': float('nan'),
        'rms_updated': float('nan'), 'max_updated': float('nan'),
        'rms_physical': float(np.sqrt(np.mean(delta[mask]**2))),
        'max_physical': float(np.max(np.abs(delta[mask]))),
        'ion': 0, 'channel': None,
        'D0_kurucz_orig': float('nan'), 'D0_kurucz_upd': E[0],
        'D0_BC16': float('nan'), 'D0_uncert': float('nan'),
        'D0_source': 'SK/JANAF', 'D0_verified': True,
        'formula': entry_cur['formula'], 'bc16_name': entry_cur['formula'],
        'atoms': [ELEMENT_SYMBOLS[a] for a in entry_cur['real_atoms']],
        'Tmax_em': tmax_em,
    }


def compute_polyatomic(entry_orig, entry_upd, atoms, d0_info,
                       atomic, exomol):
    """
    Build comparison curves for a polyatomic.

    Reference axis: BC16 pressure convention, K_p in Pa^(N-1).
    Curves: Kurucz orig, Kurucz updated, ExoMol-assembled, and
    physical-form fit (refit a0..a4 to ExoMol-assembled K_p).
    No BC16 reference for polyatomics.
    """
    T = T_POLY
    N = len(atoms)
    ex_name = d0_info['exomol_name']

    # Kurucz original / updated (in cgs, then convert to pressure)
    log10_K_kur_cgs = kurucz_log10_K_polyatomic(entry_orig, T)
    log10_K_upd_cgs = kurucz_log10_K_polyatomic(entry_upd, T)
    log10_Kp_kur = kurucz_to_bc16_pressure_convention(log10_K_kur_cgs, T, N)
    log10_Kp_upd = kurucz_to_bc16_pressure_convention(log10_K_upd_cgs, T, N)

    # ExoMol-assembled (using verified D_0 from d0_info)
    log10_Kp_em = assemble_log10_Kp_polyatomic(
        atoms, d0_info['D0'],
        lambda Tx: exomol.Q(ex_name, Tx), atomic, T)

    # Physical-form fit (refit a0..a4 against ExoMol-assembled K_p,
    # with D0 fixed at entry_upd's E1)
    phys_coeffs, _, _ = fit_physical_polyatomic(
        entry_upd, atoms, d0_info, atomic, exomol)
    # Evaluate the fit on the dense plot grid (in Kurucz's cgs convention),
    # then convert to BC16 pressure for plotting and residuals.
    log10_K_phys_cgs = (evaluate_physical(phys_coeffs, entry_upd['n_trans'], T)
                        / np.log(10.0))
    log10_Kp_phys = kurucz_to_bc16_pressure_convention(log10_K_phys_cgs, T, N)

    # ExoMol native data points for the scatter overlay (T >= 1000 K)
    T_native = exomol.molecules[ex_name]['T']
    pts_mask = T_native >= 1000.0
    em_T_pts = T_native[pts_mask]
    em_log10K_pts = assemble_log10_Kp_polyatomic(
        atoms, d0_info['D0'],
        lambda Tx: exomol.Q(ex_name, Tx), atomic, em_T_pts)
    pts_finite = np.isfinite(em_log10K_pts)
    em_T_pts = em_T_pts[pts_finite]
    em_log10K_pts = em_log10K_pts[pts_finite]

    # Residuals (vs ExoMol, the polyatomic reference)
    delta_kur = log10_Kp_kur - log10_Kp_em
    delta_upd = log10_Kp_upd - log10_Kp_em
    delta_phys = log10_Kp_phys - log10_Kp_em

    em_mask = ~np.isnan(log10_Kp_em)
    mask = (T >= T_STAT_LO_P) & (T <= T_STAT_HI_P) & em_mask
    if mask.any():
        rms_kur  = float(np.sqrt(np.mean(delta_kur[mask]**2)))
        max_kur  = float(np.max(np.abs(delta_kur[mask])))
        rms_upd  = float(np.sqrt(np.mean(delta_upd[mask]**2)))
        max_upd  = float(np.max(np.abs(delta_upd[mask])))
        rms_phys = float(np.sqrt(np.mean(delta_phys[mask]**2)))
        max_phys = float(np.max(np.abs(delta_phys[mask])))
    else:
        rms_kur = max_kur = rms_upd = max_upd = float('nan')
        rms_phys = max_phys = float('nan')

    return {
        'kind':        'polyatomic',
        'T':           T,
        'N':           N,
        'log10_K_kurucz_orig':    log10_Kp_kur,
        'log10_K_kurucz_updated': log10_Kp_upd,
        'log10_K_exomol':         log10_Kp_em,
        'log10_K_physical':       log10_Kp_phys,
        'log10_K_bc16':           None,
        'phys_coeffs':            phys_coeffs,
        'em_T_pts':               em_T_pts,
        'em_log10K_pts':          em_log10K_pts,
        'delta_kurucz':           delta_kur,
        'delta_updated':          delta_upd,
        'delta_physical':         delta_phys,
        'rms_kurucz':  rms_kur,
        'max_kurucz':  max_kur,
        'rms_updated': rms_upd,
        'max_updated': max_upd,
        'rms_physical': rms_phys,
        'max_physical': max_phys,
        'ion':         0,
        'channel':     None,
        'D0_kurucz_orig': entry_orig['E'][0],
        'D0_kurucz_upd':  entry_upd['E'][0],
        'D0_BC16':        d0_info['D0'],   # CCCBDB-verified value
        'D0_uncert':      d0_info.get('uncert', float('nan')),
        'D0_source':      d0_info.get('source', ''),
        'D0_verified':    d0_info.get('verified', False),
        'formula':        entry_orig['formula'],
        'bc16_name':      ex_name,    # using ExoMol name as label
        'atoms':          atoms,
        'Tmax_em':        exomol.molecules[ex_name]['Tmax'],
    }


# ============================================================
# Plotting
# ============================================================

def _plot_one(pdf, r):
    """Plot one molecule onto a fresh page in pdf."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.5, 6.5), sharex=True,
        gridspec_kw={'height_ratios': [2.2, 1.0], 'hspace': 0.08})

    T = r['T']
    is_poly = (r['kind'] == 'polyatomic')

    # --- Top panel: log10 K curves
    # (orig/updated absent for species added after the April 2026 refit)
    if r['log10_K_kurucz_orig'] is not None:
        ax_top.plot(T, r['log10_K_kurucz_orig'], '-', color='C0', lw=1.5,
                    label=f'Kurucz original ($D_0$={r["D0_kurucz_orig"]:.3f})')
    if (r['log10_K_kurucz_updated'] is not None
            and abs(r['D0_kurucz_upd'] - r['D0_kurucz_orig']) > 1e-4):
        ax_top.plot(T, r['log10_K_kurucz_updated'], '--', color='C4', lw=1.5,
                    label=f'Kurucz updated ($D_0$={r["D0_kurucz_upd"]:.3f})')

    # New physical-form fit (both diatomics and polyatomics now)
    if r.get('log10_K_physical') is not None:
        ax_top.plot(T, r['log10_K_physical'], '-', color='C2', lw=1.8,
                    label='Physical form fit')

    # Reference data: scatter symbols at actual data points (no line)
    if not is_poly:
        # Diatomics: BC16 raw points
        if r.get('bc16_T_pts') is not None:
            ax_top.plot(r['bc16_T_pts'], r['bc16_log10K_pts'],
                        'o', color='k', ms=4, label='BC16 data',
                        markerfacecolor='white', markeredgewidth=1.0,
                        zorder=5)
    else:
        # Polyatomics: ExoMol-assembled points (no line connecting them).
        # The BC16 line is intentionally absent because BC16 doesn't
        # tabulate K_p for these polyatomics. For self-fit species (no
        # external reference), the points are Kurucz polynomial samples
        # over the fit window, marked with a different label.
        if r.get('em_T_pts') is not None:
            scatter_label = r.get('ref_points_label') or (
                'Kurucz poly samples (fit target)'
                if r.get('self_fit') else 'ExoMol data')
            ax_top.plot(r['em_T_pts'], r['em_log10K_pts'],
                        'o', color='k', ms=4, label=scatter_label,
                        markerfacecolor='white', markeredgewidth=1.0,
                        zorder=5)
        if r.get('log10_K_janaf') is not None:
            ax_top.plot(T, r['log10_K_janaf'], '-', color='C3', lw=1.4,
                        label='NIST-JANAF (GGchem/Stock-Kitzmann)',
                        zorder=4)
        if r.get('log10_K_current') is not None:
            lbl = 'current molecules.dat row'
            if r.get('D0_current') is not None:
                lbl += f' ($D_0$={r["D0_current"]:.3f})'
            ax_top.plot(T, r['log10_K_current'], '--', color='C5', lw=1.7,
                        label=lbl, zorder=6)

    ax_top.set_xscale('log')
    if is_poly:
        ax_top.set_ylabel(rf'$\log_{{10}}\,K_p$  [Pa$^{r["N"]-1}$]')
    else:
        ax_top.set_ylabel(r'$\log_{10}\,K_{\mathrm{cgs}}$  [cm$^3$]')
    ax_top.grid(True, alpha=0.3, which='both')
    ax_top.legend(loc='best', fontsize=8)

    # Anchor y-limits to the reference curve's range so the polynomial
    # endpoint blowup beyond ~10000 K (where the T^4 term dominates and
    # can drive log K to +/- 100) doesn't squash the meaningful T<10000 K
    # region. Pad by 5 dex on each side.
    #
    # For polyatomics, use the physical-form fit as the anchor (when
    # available) because it is well-behaved across the full plot range
    # 1000-20000 K, whereas log10_K_exomol is masked outside ExoMol's
    # T_max which can be << 10000 K. The physical-form fit is also a
    # reasonable extrapolation guide.
    if is_poly and r.get('log10_K_physical') is not None:
        ref_curve = r['log10_K_physical']
    elif r['log10_K_bc16'] is not None:
        ref_curve = r['log10_K_bc16']
    else:
        ref_curve = r.get('log10_K_exomol')
    if ref_curve is not None:
        finite = ref_curve[np.isfinite(ref_curve)]
        if finite.size > 0:
            ax_top.set_ylim(float(np.min(finite)) - 5.0,
                            float(np.max(finite)) + 5.0)

    if r.get('Tmax_em') is not None and r['Tmax_em'] < T[-1]:
        ax_top.axvline(r['Tmax_em'], color='C2', ls=':', alpha=0.4, lw=0.8)

    # Title -- always use BC16 name as the primary species label for diatomics
    # (e.g. "AlO" not Kurucz's "OAl"); polyatomics use Kurucz formula since
    # there's no BC16 entry.
    # One-line title: the species name in chemical convention, plus a
    # status tag where applicable.  D0 values live in the legend labels.
    if is_poly:
        disp = display_formula(r['formula'])
    else:
        disp = display_formula(r['bc16_name'])
    if not is_poly and r['channel'] is not None:
        ch = r['channel']
        title = (f'{disp} $\\rightarrow$ {display_formula(ch[0])} + '
                 f'{display_formula(ch[1] + "+")}')
    else:
        title = disp
    ax_top.set_title(title, fontsize=12)

    # --- Bottom panel: residuals vs reference
    ax_bot.axhline(0, color='k', lw=0.6)
    if r.get('self_fit'):
        # Kurucz orig/updated residuals are identically zero by
        # construction here; only plot the physical-form residual which
        # represents the difference between the new fit and the
        # polynomial it was fit to.  When a JANAF reference is present,
        # additionally show how far the legacy polynomial sits from it.
        if r.get('delta_kurucz_janaf') is not None:
            ax_bot.plot(T, r['delta_kurucz_janaf'], '-', color='C3',
                        lw=1.4, label='Kurucz poly - JANAF')
        if r.get('delta_current_janaf') is not None:
            ax_bot.plot(T, r['delta_current_janaf'], '--', color='C5',
                        lw=1.7, label='current row - JANAF')
        if r.get('delta_physical') is not None:
            ax_bot.plot(T, r['delta_physical'], '-', color='C2', lw=1.6,
                        label='Physical form - Kurucz poly')
    else:
        if r.get('delta_kurucz') is not None:
            ax_bot.plot(T, r['delta_kurucz'], '-', color='C0', lw=1.3,
                        label='Kurucz original')
        if (r.get('delta_updated') is not None
                and abs(r['D0_kurucz_upd'] - r['D0_kurucz_orig']) > 1e-4):
            ax_bot.plot(T, r['delta_updated'], '--', color='C4', lw=1.3,
                        label='Kurucz updated')
        if r.get('delta_physical') is not None:
            ax_bot.plot(T, r['delta_physical'], '-', color='C2', lw=1.6,
                        label='Physical form')
    ax_bot.set_xscale('log')
    ax_bot.set_xlabel('Temperature [K]')
    if r.get('ref_label'):
        ref_label = r['ref_label']
    elif r.get('self_fit'):
        ref_label = 'Kurucz poly'
    elif is_poly:
        ref_label = 'ExoMol'
    else:
        ref_label = 'BC16'
    ax_bot.set_ylabel(rf'$\Delta\log_{{10}}\,K$ vs {ref_label}')
    ax_bot.grid(True, alpha=0.3, which='both')
    ax_bot.legend(loc='best', fontsize=8)

    # Auto-scale: be generous if endpoint blowup, but anchored to stat range
    if is_poly:
        mask = (T >= T_STAT_LO_P) & (T <= T_STAT_HI_P)
    else:
        mask = (T >= T_STAT_LO_D) & (T <= T_STAT_HI_D)
    if r.get('delta_kurucz_janaf') is not None:
        d_ref = r['delta_kurucz_janaf']
    elif r.get('delta_kurucz') is not None:
        d_ref = r['delta_kurucz']
    else:
        d_ref = r['delta_physical']
    d_all = np.concatenate([
        d_ref[mask],
        r['delta_updated'][mask]
        if (r.get('delta_updated') is not None
            and abs(r['D0_kurucz_upd']-r['D0_kurucz_orig']) > 1e-4)
        else d_ref[mask],
    ])
    d_all = d_all[~np.isnan(d_all)]
    if len(d_all) > 0:
        stat_max = max(0.1, 1.5 * np.max(np.abs(d_all)))
    else:
        stat_max = 0.5
    full_max = max(1e-9, np.nanmax(np.abs(d_ref)))
    if full_max <= 2.0 * stat_max:
        ylim = max(stat_max, full_max * 1.1)
    else:
        ylim = 2.0 * stat_max
    ax_bot.set_ylim(-ylim, ylim)

    # Stats annotation
    if r.get('self_fit'):
        stats = (f'In fit window {r["fit_T_lo"]:.0f}-{r["fit_T_hi"]:.0f} K:  '
                 f'physical-vs-Kurucz: RMS={r["rms_physical"]:.4f}, '
                 f'max={r["max_physical"]:.4f}')
        if r.get('rms_kurucz_janaf') is not None:
            stats += (f'\n1500-{r["fit_T_hi"]:.0f} K vs JANAF:  '
                      f'Kurucz: RMS={r["rms_kurucz_janaf"]:.3f}, '
                      f'max={r["max_kurucz_janaf"]:.3f}  |  '
                      f'physical: RMS={r["rms_physical_janaf"]:.3f}, '
                      f'max={r["max_physical_janaf"]:.3f}')
        if r.get('rms_current_janaf') is not None:
            stats += (f'  |  current row: RMS={r["rms_current_janaf"]:.3f}, '
                      f'max={r["max_current_janaf"]:.3f}')
    elif r.get('is_new'):
        ref = r.get('ref_label', 'BC16')
        stats = (f'In stat range:  filed row vs {ref}: '
                 f'RMS={r["rms_physical"]:.3f}, max={r["max_physical"]:.3f}')
    else:
        stats = (f'In stat range:  '
                 f'orig: RMS={r["rms_kurucz"]:.3f}, max={r["max_kurucz"]:.3f}')
        if abs(r['D0_kurucz_upd']-r['D0_kurucz_orig']) > 1e-4:
            stats += (f'  |  updated: RMS={r["rms_updated"]:.3f}, '
                      f'max={r["max_updated"]:.3f}')
        if r.get('rms_physical') is not None:
            stats += (f'  |  physical: RMS={r["rms_physical"]:.3f}, '
                      f'max={r["max_physical"]:.3f}')
    ax_bot.text(0.02, 0.97, stats, transform=ax_bot.transAxes,
                fontsize=7.5, va='top',
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'))

    # Enforce common x-axis range across all pages, with explicit
    # round-number tick labels (the default log locator labels only
    # 10^3 and 10^4 over this range)
    ax_top.set_xlim(1000.0, 20000.0)
    ax_bot.set_xlim(1000.0, 20000.0)
    xticks = [1000, 2000, 3000, 5000, 10000, 20000]
    ax_bot.set_xticks(xticks)
    ax_bot.set_xticklabels([f'{t:,d}' for t in xticks])
    ax_bot.xaxis.set_minor_formatter(plt.NullFormatter())

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _plot_cover(pdf, results_diat_n, results_diat_i, results_poly,
                results_new=()):
    """Cover page summarizing the contents."""
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.95, 'Kurucz molecules.dat: original / updated / physical-form fit',
             ha='center', fontsize=14, fontweight='bold')
    counts = (f'{len(results_diat_n)} diatomic neutrals, '
              f'{len(results_diat_i)} diatomic ions, '
              f'{len(results_poly)} polyatomics')
    if results_new:
        counts += (f', {len(results_new)} diatomics added July 2026 '
                   f'(BC16 fit only; no Kurucz ancestry)')
    fig.text(0.5, 0.92, counts, ha='center', fontsize=11)

    lines = []
    lines.append('')
    lines.append('SOURCES (Kurucz updated):')
    lines.append('  Diatomics:     BC16 Table 1 ($D_0$ adopted)')
    lines.append('  Diatomic ions: $E_1 = D_0(\\mathrm{AB}^+)_\\mathrm{BC16}$ '
                 '$- \\mathrm{IE}(B)$')
    lines.append('  Polyatomics:   ATcT 1.220 / Koput 2023 / CCCBDB')
    lines.append('')
    lines.append('PHYSICAL FORM (diatomics + polyatomics):')
    lines.append('  $\\ln K = D_0/(k_B T) - 1.5 n_{trans}\\ln T$')
    lines.append('         $+ a_0 + a_1/T + a_2 \\ln T + a_3 T + a_4/T^2$')
    lines.append('  Replaces Kurucz polynomial $\\{T, T^2, T^3, T^4\\}$')
    lines.append('  basis with physically-motivated basis')
    lines.append('  $\\{1, 1/T, \\ln T, T, 1/T^2\\}$. Same parameter count.')
    lines.append('  Diatomic fits: refit to BC16 native grid (T >= 1000 K).')
    lines.append('  Polyatomic fits: refit to ExoMol-assembled K_p on')
    lines.append('    ExoMol native grid (T >= 1000 K).')
    lines.append('')
    lines.append('REFERENCES (in residual plots):')
    lines.append('  Diatomics:    BC16 (tabulated $K_p$ in Pa)')
    lines.append('  Polyatomics:  ExoMol $Q$ + verified $D_0$ '
                 '(assembled $K_p$ in Pa$^{N-1}$)')
    lines.append('')
    lines.append('STAT RANGES:')
    lines.append('  Diatomics:    1500 - 8000 K')
    lines.append('  Polyatomics:  1500 - 6000 K (and within ExoMol $Q$ range)')

    text = '\n'.join(lines)
    fig.text(0.05, 0.85, text, fontsize=9, va='top', family='monospace')

    # --- Top changes table (all categories combined, by |delta D0|) ---
    all_results = results_diat_n + results_diat_i + results_poly
    sortable = [
        (r, abs(r['D0_kurucz_upd'] - r['D0_kurucz_orig']))
        for r in all_results
        if abs(r['D0_kurucz_upd'] - r['D0_kurucz_orig']) > 1e-3
    ]
    sortable.sort(key=lambda x: -x[1])
    top_changes = sortable[:25]

    rows = []
    rows.append(f'{"formula":<8s}  {"orig":>8s}  {"upd":>8s}  '
                f'{"delta":>7s}  {"RMS_o":>6s} {"RMS_u":>6s}  category')
    rows.append('-' * 80)
    for r, _ in top_changes:
        d_old = r['D0_kurucz_orig']
        d_new = r['D0_kurucz_upd']
        cat = ('diat-ion' if (r['kind']=='diatomic' and r['ion']) else
               'diat'     if r['kind']=='diatomic' else
               'poly')
        rows.append(f'{r["formula"]:<8s}  {d_old:>8.3f}  {d_new:>8.3f}  '
                    f'{d_new-d_old:>+7.3f}  '
                    f'{r["rms_kurucz"]:>6.3f} {r["rms_updated"]:>6.3f}  {cat}')
    fig.text(0.05, 0.55, '\n'.join(rows),
             fontsize=8.5, va='top', family='monospace')

    fig.text(0.5, 0.04,
             'Generated by compare_unified.py',
             ha='center', fontsize=8, style='italic')
    pdf.savefig(fig)
    plt.close(fig)


# ============================================================
# Top-level driver
# ============================================================

def run(orig_dat='molecules.dat',
        updated_dat='molecules_updated.dat',
        current_dat=None,
        out_pdf='kurucz_unified_report.pdf',
        out_csv='kurucz_unified_summary.csv'):
    bc16   = load_bc16()
    atomic = AtomicData()
    exomol = ExoMolData('.')

    k_orig = parse_molecules_dat(orig_dat)
    k_upd  = parse_molecules_dat(updated_dat)
    upd_by_code = {e['code']: e for e in k_upd}

    matches = match_kurucz_to_bc16(k_orig, bc16)

    # ----- DIATOMICS (matched against BC16) -----
    diat_neutrals = []
    diat_ions     = []
    for entry_orig, bc_name, _ in matches:
        if bc_name is None:
            continue
        if entry_orig['n_real_atoms'] != 2:
            continue
        entry_upd = upd_by_code[entry_orig['code']]
        bc_mol = bc16[bc_name]
        try:
            r = compute_diatomic(entry_orig, entry_upd, bc_mol, atomic, exomol)
        except Exception as e:
            print(f'  diatomic {entry_orig["formula"]}: skipped ({e})')
            continue
        if r['ion'] == 0:
            diat_neutrals.append(r)
        else:
            diat_ions.append(r)

    diat_neutrals.sort(key=lambda r: -r['max_kurucz'])
    diat_ions.sort(key=lambda r: -r['max_kurucz'])

    # ----- NEW DIATOMICS (post-April species from the current file) -----
    # Species in current_dat with a BC16 match but no entry in orig_dat
    # have no Kurucz-polynomial ancestry: their pages show the BC16
    # reference and the physical-form fit only.
    diat_new = []
    cur_by_code = {}
    if current_dat is not None:
        k_cur = parse_molecules_dat(current_dat)
        cur_by_code = {e['code']: e for e in k_cur}
        lineage_codes = {e['code'] for e in k_orig}
        for e_cur, bc_name, _ in match_kurucz_to_bc16(k_cur, bc16):
            if bc_name is None or e_cur['code'] in lineage_codes:
                continue
            try:
                r = compute_diatomic(None, e_cur, bc16[bc_name],
                                     atomic, exomol)
            except Exception as e:
                print(f'  new diatomic {e_cur["formula"]}: skipped ({e})')
                continue
            diat_new.append(r)
        diat_new.sort(key=lambda r: r['bc16_name'])

    # ----- POLYATOMICS (matched against POLYATOMIC_D0 + .pf availability) -----
    poly_results = []
    # Build a Kurucz-by-atom-composition lookup
    def kurucz_by_atoms(target_atoms, n_real):
        target = sorted(target_atoms)
        for e in k_orig:
            if (e['is_standard_molecule'] and e['n_real_atoms'] == n_real
                    and sorted(ELEMENT_SYMBOLS[a] for a in e['real_atoms']) == target):
                return e
        return None

    for ex_name, d0_info in POLYATOMIC_D0.items():
        atoms = d0_info['atoms']
        if not exomol.has(ex_name):
            continue
        entry_orig = kurucz_by_atoms(atoms, len(atoms))
        if entry_orig is None:
            continue
        if entry_orig['ion'] != 0:
            continue   # skip H3+ etc.
        entry_upd = upd_by_code[entry_orig['code']]
        try:
            r = compute_polyatomic(entry_orig, entry_upd, atoms, d0_info,
                                   atomic, exomol)
        except Exception as e:
            print(f'  polyatomic {ex_name}: skipped ({e})')
            continue
        poly_results.append(r)

    poly_results.sort(key=lambda r: -r['max_kurucz'])

    # ----- POLYATOMICS without an external reference (self-fit) -----
    # These are the polyatomics in molecules.dat that either (a) lack a
    # POLYATOMIC_D0 entry, (b) are ions, or (c) have a D0 entry but no
    # ExoMol .pf available. For these we fit the new physical form to
    # Kurucz's own polynomial K_eq sampled in 1000-6000 K, so that every
    # species in molecules.dat ends up with a physical-form set of
    # coefficients.
    poly_selffit_results = []
    # Build a set of (sorted-atom-string, ion) keys already covered above
    covered = set()
    for r in poly_results:
        atoms_key = tuple(sorted(r['atoms']))
        covered.add((atoms_key, r['ion']))

    # NIST-JANAF reference fits (GGchem's dispol_StockKitzmann.dat) for
    # self-fit species, matched by atom multiset (neutrals only).
    try:
        from ggchem_loader import parse_dispol, ln_K_kurucz as _janaf_lnK
        janaf_by_atoms = {}
        for name, rec in parse_dispol('dispol_StockKitzmann.dat').items():
            key = tuple(sorted(
                s for s, n in rec['atoms'].items() for _ in range(n)))
            # Isomer collisions (AlOH vs HAlO, FOO vs OFO, ...): keep the
            # more-bound (ground) isomer -- the larger K at 2000 K -- which
            # dominates the equilibrium abundance.
            old = janaf_by_atoms.get(key)
            if (old is None or float(_janaf_lnK(rec, np.array([2000.0]))[0])
                    > float(_janaf_lnK(old, np.array([2000.0]))[0])):
                janaf_by_atoms[key] = rec
    except Exception as exc:
        print(f'  JANAF reference unavailable ({exc})')
        janaf_by_atoms = {}

    for entry_orig in k_orig:
        if not entry_orig['is_standard_molecule']:
            continue
        if entry_orig['n_real_atoms'] < 3:
            continue   # diatomics handled elsewhere
        atoms_sym = [ELEMENT_SYMBOLS[a] for a in entry_orig['real_atoms']]
        key = (tuple(sorted(atoms_sym)), entry_orig['ion'])
        if key in covered:
            continue
        entry_upd = upd_by_code[entry_orig['code']]
        janaf_rec = (janaf_by_atoms.get(tuple(sorted(atoms_sym)))
                     if entry_orig['ion'] == 0 else None)
        try:
            r = compute_polyatomic_self_fit(entry_orig, entry_upd, atoms_sym,
                                            janaf_rec=janaf_rec,
                                            cur_entry=cur_by_code.get(
                                                entry_orig['code']))
        except Exception as exc:
            print(f'  poly self-fit {entry_orig["formula"]}: skipped ({exc})')
            continue
        poly_selffit_results.append(r)

    poly_selffit_results.sort(key=lambda r: r['formula'])

    # ----- NEW POLYATOMICS (post-April, fit to JANAF/SK) -----
    poly_new = []
    if current_dat is not None and janaf_by_atoms:
        lineage_codes = {e['code'] for e in k_orig}
        for e_cur in k_cur:
            if (not e_cur['is_standard_molecule']
                    or e_cur['n_real_atoms'] < 3
                    or e_cur['ion'] != 0
                    or e_cur['code'] in lineage_codes):
                continue
            atoms_sym = [ELEMENT_SYMBOLS[a] for a in e_cur['real_atoms']]
            rec = janaf_by_atoms.get(tuple(sorted(atoms_sym)))
            # ExoMol-fit rows reference the ExoMol assembly instead
            em_ref = None
            for d0i in POLYATOMIC_D0.values():
                if (sorted(d0i['atoms']) == sorted(atoms_sym)
                        and exomol.has(d0i['exomol_name'])):
                    em_ref = (d0i, exomol, atomic)
                    break
            if rec is None and em_ref is None:
                print(f'  new polyatomic {e_cur["formula"]}: no reference')
                continue
            try:
                poly_new.append(
                    compute_polyatomic_new(e_cur, rec, em_ref=em_ref))
            except Exception as exc:
                print(f'  new polyatomic {e_cur["formula"]}: skipped ({exc})')
        poly_new.sort(key=lambda r: r['formula'])

    # ----- Write PDF -----
    with PdfPages(out_pdf) as pdf:
        _plot_cover(pdf, diat_neutrals, diat_ions, poly_results, diat_new)
        for r in diat_neutrals:
            _plot_one(pdf, r)
        for r in diat_ions:
            _plot_one(pdf, r)
        for r in diat_new:
            _plot_one(pdf, r)
        for r in poly_results:
            _plot_one(pdf, r)
        for r in poly_new:
            _plot_one(pdf, r)
        for r in poly_selffit_results:
            _plot_one(pdf, r)

    # ----- Write CSV -----
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['molecule', 'kind', 'ion', 'N',
                    'D0_orig', 'D0_updated', 'D0_modern_ref', 'deltaD0_upd',
                    'rms_orig', 'max_orig', 'rms_updated', 'max_updated',
                    'rms_physical', 'max_physical',
                    'reference', 'D0_source'])
        all_results = (
            [(r, 'diatomic-neutral')   for r in diat_neutrals] +
            [(r, 'diatomic-ion')       for r in diat_ions]     +
            [(r, 'diatomic-new')       for r in diat_new]      +
            [(r, 'polyatomic')         for r in poly_results]  +
            [(r, 'polyatomic-new')     for r in poly_new]      +
            [(r, 'polyatomic-selffit') for r in poly_selffit_results])
        for r, label in all_results:
            if label == 'polyatomic-selffit':
                ref = 'Kurucz polynomial (self-fit)'
            elif r['kind']=='diatomic':
                ref = 'BC16'
            else:
                ref = 'ExoMol-assembled'
            d0_src = (r.get('D0_source', 'BC16 Table 1')
                      if r['kind']=='polyatomic' else 'BC16 Table 1')
            n = r.get('N', 2)
            rms_phys_str = (f'{r["rms_physical"]:.4f}'
                            if r.get('rms_physical') is not None else 'NA')
            max_phys_str = (f'{r["max_physical"]:.4f}'
                            if r.get('max_physical') is not None else 'NA')
            d0_modern = (f'{r["D0_BC16"]:.4f}'
                         if not np.isnan(r["D0_BC16"]) else 'NA')
            w.writerow([
                r['formula'], label, r['ion'], n,
                f'{r["D0_kurucz_orig"]:.4f}', f'{r["D0_kurucz_upd"]:.4f}',
                d0_modern,
                f'{r["D0_kurucz_upd"]-r["D0_kurucz_orig"]:+.4f}',
                f'{r["rms_kurucz"]:.4f}', f'{r["max_kurucz"]:.4f}',
                f'{r["rms_updated"]:.4f}', f'{r["max_updated"]:.4f}',
                rms_phys_str, max_phys_str,
                ref, d0_src,
            ])

    print(f'Wrote {out_pdf}')
    print(f'Wrote {out_csv}')
    print()
    print(f'Summary by category:')
    print(f'{"category":<22s}  count  '
          f'{"avg RMS orig":>13s}  {"avg RMS upd":>13s}  {"avg RMS phys":>13s}')
    for label, results in [('diatomic neutrals',     diat_neutrals),
                           ('diatomic ions',         diat_ions),
                           ('diatomics new (2026)',  diat_new),
                           ('polyatomics (ref)',     poly_results),
                           ('polyatomics new (2026)', poly_new),
                           ('polyatomics (self-fit)', poly_selffit_results)]:
        if not results:
            continue
        rms_orig = np.mean([r['rms_kurucz']  for r in results])
        rms_upd  = np.mean([r['rms_updated'] for r in results])
        rms_phys_vals = [r['rms_physical'] for r in results
                         if r.get('rms_physical') is not None]
        rms_phys = np.mean(rms_phys_vals) if rms_phys_vals else float('nan')
        rms_phys_str = (f'{rms_phys:>13.4f}' if not np.isnan(rms_phys)
                        else f'{"NA":>13s}')
        print(f'{label:<22s}  {len(results):>5d}  '
              f'{rms_orig:>13.3f}  {rms_upd:>13.3f}  {rms_phys_str}')

    return diat_neutrals, diat_ions, poly_results, poly_selffit_results


if __name__ == '__main__':
    run()
