"""
plot_keq_fits.py -- regenerate the K_eq fit-validation atlas (comp_pf.pdf).

One page per neutral-diatomic "BC16 fit" row in data/molecules.dat:
  top    log10 K_cgs(T): BC16 Table 7 reference points, the filed
         physical-form fit, and (where the species predates the April 2026
         refit) the old Kurucz polynomial for comparison
  bottom residual (filed fit - BC16) at the reference points, with the
         stellar band 1000-6000 K shaded and its max |dln K| annotated

Successor to ~/kurucz/upgrade/comp_pf.py (April 2026), whose data-loader
modules no longer exist; the original PDF with its ExoMol polyatomic pages
is preserved as comp_pf_apr2026.pdf.  This version reads the fits directly
from molecules.dat (ground truth, not a refit), so it validates what the
code actually uses.  Molecular-ion rows (8) are not plotted -- their D0
sign convention has no reconstructed reference path.

Usage:
  python3 tools/plot_keq_fits.py [-o out.pdf] [--old-ref GITREV]
      defaults: -o ~/kurucz/upgrade/comp_pf.pdf, --old-ref bb386c9^
"""

import argparse
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from fit_molecule_keq import (KB_PA_CM3_PER_K, K_BOLTZMANN_EV, REPO_ROOT,
                              MOLECULES_DAT, load_table7, load_table1_d0,
                              bc16_name_for, parse_molecules_rows)

plt.rcParams.update({
    'font.family': 'serif', 'mathtext.fontset': 'stix',
    'font.serif': ['STIXGeneral'], 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.dpi': 300,
})

T_BAND = (1000.0, 6000.0)          # stellar band highlighted in residuals


def decode_atoms(code):
    """Species code -> (atom list, ion stage); mirrors READMOL."""
    icode = round(code * 100)
    ion = icode % 100
    n = icode // 100
    atoms = []
    while n > 0:
        atoms.append(n % 100)
        n //= 100
    return atoms, ion


def lnk_physical(d0, coeffs, t, n_trans=1):
    """NMOLEC's physical-form ln K_eq."""
    a0, a1, a2, a3, a4 = coeffs
    return (d0 / (K_BOLTZMANN_EV * t) - 1.5 * n_trans * np.log(t)
            + a0 + a1 / t + a2 * np.log(t) + a3 * t + a4 / t**2)


def lnk_kurucz_poly(d0, coeffs, t, n_trans=1):
    """Pre-April-2026 Kurucz polynomial ln K_eq (NMOLEC's old expression)."""
    e2, e3, e4, e5, e6 = coeffs
    return (d0 / (K_BOLTZMANN_EV * t) - 1.5 * n_trans * np.log(t)
            - e2 + e3 * t - e4 * t**2 + e5 * t**3 - e6 * t**4)


def load_old_rows(gitrev):
    """code -> (D0, poly coeffs) from the pre-refit molecules.dat."""
    try:
        txt = subprocess.run(
            ['git', 'show', f'{gitrev}:data/molecules.dat'],
            capture_output=True, text=True, check=True,
            cwd=REPO_ROOT).stdout
    except subprocess.CalledProcessError:
        print(f'warning: git show {gitrev} failed; no old-poly overlays')
        return {}
    old = {}
    for line in txt.split('\n'):
        if not line.strip() or line.startswith('!'):
            continue
        code = float(line[10:22])
        if len(line) >= 90 and line[30:90].strip() and line[23:30].strip():
            d0 = float(line[23:30])
            coeffs = [float(line[30 + 12 * k:42 + 12 * k]) for k in range(5)]
            old[round(code * 100)] = (d0, coeffs)
    return old


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', type=Path,
                    default=Path.home() / 'kurucz' / 'upgrade' / 'comp_pf.pdf')
    ap.add_argument('--old-ref', default='bb386c9^',
                    help='git rev of the pre-refit molecules.dat')
    args = ap.parse_args()

    tgrid, table7 = load_table7()
    old_rows = load_old_rows(args.old_ref)

    # Collect plottable rows: neutral diatomic BC16 fits
    species = []
    for _, label, code, d0, coeffs, comment in parse_molecules_rows():
        if 'BC16 fit' not in comment or coeffs is None:
            continue
        atoms, ion = decode_atoms(code)
        if ion != 0 or len(atoms) != 2 or 100 in atoms:
            continue
        name = bc16_name_for(atoms, table7)
        if name is None:
            print(f'warning: {label} has no BC16 match; skipped')
            continue
        species.append((code, label, d0, coeffs, comment, name))
    species.sort()

    ln10 = np.log(10.0)
    mask = tgrid >= 100.0
    t_pts = tgrid[mask]
    t_dense = np.logspace(2, 4, 300)
    n_new = sum('new species' in s[4] for s in species)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output) as pdf:
        # --- Cover page
        fig = plt.figure(figsize=(8.5, 6.5))
        fig.text(0.5, 0.72, 'molecules.dat equilibrium-constant fits vs. '
                 'Barklem & Collet (2016)', ha='center', fontsize=15)
        fig.text(0.5, 0.60,
                 f'{len(species)} neutral diatomics ("BC16 fit" rows), '
                 f'of which {n_new} added July 2026\n'
                 'curves are the coefficients as filed in data/molecules.dat',
                 ha='center', fontsize=11)
        fig.text(0.5, 0.42,
                 'reference: BC16 Table 7 log10 Kp converted to the Kurucz\n'
                 'association convention, '
                 r'$\log_{10}K = \log_{10}(k_B T) - \log_{10}K_p$'
                 '\n\nold Kurucz polynomial overlaid where the species '
                 f'predates the April 2026 refit ({args.old_ref})\n'
                 'molecular ions (8 rows) and ExoMol polyatomics not shown; '
                 'see comp_pf_apr2026.pdf for the April polyatomic pages',
                 ha='center', fontsize=9)
        fig.text(0.5, 0.18, 'generated by tools/plot_keq_fits.py',
                 ha='center', fontsize=8, style='italic')
        pdf.savefig(fig)
        plt.close(fig)

        # --- One page per species
        for code, label, d0, coeffs, comment, bc_name in species:
            is_new = 'new species' in comment
            ref_pts = (np.log10(KB_PA_CM3_PER_K * t_pts)
                       - table7[bc_name][mask])
            fit_dense = lnk_physical(d0, coeffs, t_dense) / ln10
            fit_pts = lnk_physical(d0, coeffs, t_pts) / ln10

            fig, (ax, axr) = plt.subplots(
                2, 1, figsize=(8.5, 6.5), sharex=True,
                gridspec_kw={'height_ratios': [2.2, 1.0], 'hspace': 0.08})

            ax.plot(t_dense, fit_dense, '-', color='C2', lw=1.8,
                    label=f'physical-form fit as filed ($D_0$={d0:.3f} eV)')
            old = old_rows.get(round(code * 100))
            if old is not None:
                old_dense = lnk_kurucz_poly(old[0], old[1], t_dense) / ln10
                ax.plot(t_dense, old_dense, '--', color='C0', lw=1.3,
                        label=f'Kurucz polynomial, pre-refit '
                              f'($D_0$={old[0]:.3f})')
            ax.plot(t_pts, ref_pts, 'o', color='k', ms=4,
                    markerfacecolor='white', markeredgewidth=1.0, zorder=5,
                    label='BC16 Table 7')
            ax.set_xscale('log')
            ax.set_ylabel(r'$\log_{10}\,K_{\mathrm{cgs}}$  [cm$^3$]')
            lo, hi = ref_pts.min(), ref_pts.max()
            ax.set_ylim(lo - 5.0, hi + 5.0)
            ax.legend(loc='best', fontsize=8, frameon=False)
            tag = '   [new species, July 2026]' if is_new else ''
            ax.set_title(f'{label}  (code {code:.2f}, BC16 "{bc_name}")'
                         f'{tag}', fontsize=11)

            res = (fit_pts - ref_pts) * ln10          # residual in ln K
            band = (t_pts >= T_BAND[0]) & (t_pts <= T_BAND[1])
            axr.axhline(0.0, color='k', lw=0.6)
            axr.axvspan(*T_BAND, color='C2', alpha=0.08)
            axr.plot(t_pts, res, 'o-', color='C3', ms=3, lw=0.9)
            axr.set_xscale('log')
            axr.set_xlabel('T  [K]')
            axr.set_ylabel(r'fit $-$ BC16  [$\Delta \ln K$]')
            axr.annotate(
                f'max $|\\Delta\\ln K|$ in {T_BAND[0]:.0f}-{T_BAND[1]:.0f} K:'
                f' {np.abs(res[band]).max():.3f}',
                xy=(0.98, 0.92), xycoords='axes fraction',
                ha='right', va='top', fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)

    print(f'wrote {args.output}: cover + {len(species)} species pages '
          f'({n_new} new)')


if __name__ == '__main__':
    main()
