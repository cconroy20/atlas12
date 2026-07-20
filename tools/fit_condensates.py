#!/usr/bin/env python3
"""Build data/condensates.dat for the ATLAS12 Cond-limit (Stage 1).

Filed convention (one canonical form, evaluated by READCOND/NMOLEC):

  ln S = sum_i nu_i * ln(n_i * kB * T)  +  c0/T + c1*lnT + c2 + c3*T + c4*T^2

with n_i the number density [cm^-3] of reference species i (a free
neutral atom, or a gas molecule carried in molecules.dat), and
kB = 1.380662e-16 erg/K (GGchem's value, kept for source consistency).
The c-polynomial is ln K(T) = -dG(ref species -> solid)/RT with all
standard-pressure constants absorbed into c2.

Sources: GGchem DustChem.dat (Woitke+ 2018) via tools/dustchem_loader.
Per-species selection policy:
  - active GGchem fit when it is an atom-referenced dG form (1/2/5);
  - vapor-pressure (fit-3) species keep that reference: monatomic ones
    (Fe, Ni) file the atom; polyatomic ones file the gas molecule IF it
    exists in our network (TiO2, SiO2, SiO, VO), else fall back to the
    best atom-referenced alternative (ZrO2 -> SUPCRTBL Stock fit).
Liquids are excluded: no [l] phase forms above 1000 K at any pressure
in the Stage-0 reference sweeps.

Usage:
  python3 fit_condensates.py            # refit + write .dat + atlas
  python3 fit_condensates.py --validate # reread .dat, verify vs sources
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dustchem_loader import load_dustchem, lnk_atoms, ln_psat  # noqa: E402

REPO = os.path.expanduser("~/kurucz/atlas12")
OUT_DAT = os.path.join(REPO, "data", "condensates.dat")
OUT_PDF = os.path.expanduser("~/kurucz/upgrade/comp_cond.pdf")
SUMMARY_CSV = os.path.expanduser(
    "~/kurucz/upgrade/cond_ref/cond_ref_summary.csv")

TLO, THI = 700.0, 2600.0
TGRID = np.arange(TLO, THI + 1.0, 10.0)
LN10 = math.log(10.0)

# name-in-DustChem, mode, source selector, note
#   mode 'A'  = atoms via dG form (fit 1/2/5)
#   mode 'PA' = vapor pressure, monatomic (files the atom)
#   mode 'M:<label>' = vapor pressure, files gas molecule <label>
#   selector: 'active' or ('alt', code, substring-of-source)
SPECIES = [
    ("Al2O3[s]",      "A",       "active", "corundum"),
    ("CaTiO3[s]",     "A",       "active", "perovskite; SH90 only, T>=1000K"),
    ("Ti2O3[s]",      "A",       "active", "Kitzmann+2024 Stock fit"),
    ("Ti3O5[s]",      "A",       "active", ""),
    ("Ti4O7[s]",      "A",       "active", ""),
    ("TiO2[s]",       "M:TiO2",  "active", "rutile"),
    ("ZrO2[s]",       "A",       ("alt", 5, "SUPCRTBL"),
     "baddeleyite; no ZrO2 gas in network"),
    ("Ca2Al2SiO7[s]", "A",       "active", "gehlenite"),
    ("CaAl2Si2O8[s]", "A",       "active", "anorthite (SH90 known bad)"),
    ("MgAl2O4[s]",    "A",       "active", "spinel"),
    ("Ca2MgSi2O7[s]", "A",       "active", "akermanite"),
    ("CaMgSi2O6[s]",  "A",       "active", "diopside"),
    ("CaSiO3[s]",     "A",       "active", "wollastonite"),
    ("Mg2SiO4[s]",    "A",       "active", "forsterite"),
    ("MgSiO3[s]",     "A",       "active", "enstatite"),
    ("SiO2[s]",       "M:SiO2",  "active", "quartz"),
    ("SiO[s]",        "M:SiO",   "active",
     "amorphous SiO; GGchem use_SiO default, revisitable"),
    ("Fe[s]",         "PA",      "active", "iron"),
    ("Ni[s]",         "PA",      "active", "nickel"),
    ("VO[s]",         "M:VO",    "active", ""),
    ("V2O3[s]",       "A",       "active", ""),
]


def pick_fit(cond, selector):
    if selector == "active":
        fit = cond.active_fit
        if fit is None:
            raise ValueError(f"{cond.name}: no active fit")
        return fit
    _, code, key = selector
    for f in cond.fits:
        if f.code == code and key in f.source:
            return f
    raise ValueError(f"{cond.name}: no alt fit{code} matching '{key}'")


def target_lnk(cond, fit, mode, T):
    """ln K'(T) in the filed convention (per-point)."""
    if mode.startswith("M") or mode == "PA":
        if fit.code != 3:
            raise ValueError(f"{cond.name}: mode {mode} needs fit-3")
        return -ln_psat(fit, T)
    lnk, pst = lnk_atoms(cond, fit, T)
    nutot = cond.natoms()
    return lnk - nutot * math.log(pst)


def refs_of(cond, mode):
    """Reference species list [(nu, label)] in the filed convention."""
    if mode.startswith("M:"):
        return [(1, mode.split(":")[1])]
    if mode == "PA":
        return [(1, cond.stoich[0][1])]
    return list(cond.stoich)


def refit(cond, fit, mode):
    y = np.array([target_lnk(cond, fit, mode, t) for t in TGRID])
    A = np.column_stack([1.0/TGRID, np.log(TGRID), np.ones_like(TGRID),
                         TGRID, TGRID**2])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = (A @ c - y)/LN10          # dex
    return c, y, resid


def build():
    db = load_dustchem()
    rows = []
    print(f"{'species':<16}{'mode':<9}{'src':>5}  {'prec':>6}  "
          f"{'max|resid| dex':>14}")
    for name, mode, selector, note in SPECIES:
        cond = db[name]
        fit = pick_fit(cond, selector)
        c, y, resid = refit(cond, fit, mode)
        wr = np.max(np.abs(resid))
        prec = f"{fit.precision:.3f}" if fit.precision else "  --"
        print(f"{name:<16}{mode:<9}fit{fit.code:>2}  {prec:>6}  {wr:14.5f}")
        rows.append(dict(name=name.replace("[s]", ""), cond=cond, fit=fit,
                         mode=mode, note=note, coeff=c, target=y,
                         resid=resid))
    return rows


def write_dat(rows):
    hdr = f"""\
# ATLAS12 equilibrium-condensate data (Cond-limit; no dust opacity)
#
# ln S = sum_i nu_i * ln(n_i * kB * T) + c0/T + c1*lnT + c2 + c3*T + c4*T^2
#   n_i  : number density [cm^-3] of reference species i
#          (element symbol = free neutral atom; else a molecules.dat label)
#   kB   : 1.380662e-16 erg/K
#   valid: {TLO:.0f}-{THI:.0f} K refit range (source validity varies; see atlas)
#
# Source: GGchem DustChem.dat (Woitke+ 2018) fits, re-referenced and refit
# by tools/fit_condensates.py.  Atlas: ~/kurucz/upgrade/comp_cond.pdf
#
"""
    with open(OUT_DAT, "w") as f:
        f.write(hdr)
        f.write(f"{len(rows)}\n")
        for r in rows:
            refs = refs_of(r["cond"], r["mode"])
            reftxt = f"{len(refs)}  " + "  ".join(
                f"{nu:1d} {sp:<5}" for nu, sp in refs)
            coefs = "".join(f"{v: .7E} " for v in r["coeff"])
            f.write(f"{r['name']:<13} {reftxt:<32} {coefs}  "
                    f"! fit{r['fit'].code} {r['fit'].source[:48]}\n")
    print(f"\nwrote {OUT_DAT}")


def read_dat(path=OUT_DAT):
    rows = []
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#") and ln.strip()]
    n = int(lines[0])
    for ln in lines[1:n+1]:
        body = ln.split("!")[0].split()
        name = body[0]
        nref = int(body[1])
        refs = [(int(body[2+2*k]), body[3+2*k]) for k in range(nref)]
        coeff = np.array([float(x) for x in body[2+2*nref:7+2*nref]])
        rows.append((name, refs, coeff))
    assert len(rows) == n
    return rows


def validate():
    db = load_dustchem()
    filed = {name: (refs, coeff) for name, refs, coeff in read_dat()}
    A = np.column_stack([1.0/TGRID, np.log(TGRID), np.ones_like(TGRID),
                         TGRID, TGRID**2])
    nbad = 0
    for name, mode, selector, note in SPECIES:
        short = name.replace("[s]", "")
        if short not in filed:
            print(f"{short:<14} MISSING from {OUT_DAT}")
            nbad += 1
            continue
        refs, coeff = filed[short]
        cond = db[name]
        fit = pick_fit(cond, selector)
        want = refs_of(cond, mode)
        y = np.array([target_lnk(cond, fit, mode, t) for t in TGRID])
        d = np.max(np.abs(A @ coeff - y))/LN10
        ok = (d < 0.02) and (want == refs)
        print(f"{short:<14} refs={'OK ' if want == refs else 'BAD'} "
              f"max|filed-source|={d:8.5f} dex  {'PASS' if ok else 'FAIL'}")
        nbad += 0 if ok else 1
    print(f"\n{len(SPECIES)-nbad}/{len(SPECIES)} PASS")
    return nbad == 0


def appearance_table():
    """condensate -> appearance-T row from the Stage-0 summary CSV."""
    import csv
    out = {}
    if not os.path.exists(SUMMARY_CSV):
        return out, []
    with open(SUMMARY_CSV) as f:
        rd = list(csv.reader(f))
    header = rd[0]
    for row in rd[1:]:
        if not row or not row[0]:
            break
        out[row[0]] = row[1:]
    return out, header[1:]


def atlas(rows):
    """comp_cond.pdf in the comp_pf.py house style: 8.5x6.5 two-panel
    pages (curve + residual), log-T axis with round-number comma ticks,
    black open circles for the fit target, C-cycle lines, grid alpha 0.3,
    stats box, and a monospace cover page."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from comp_pf import display_formula

    XTICKS = [700, 1000, 1500, 2000, 2500]
    ALT_STYLES = [("C3", "-", 1.4), ("C5", "--", 1.7), ("C4", "--", 1.5)]

    app, appcols = appearance_table()

    def t_app(nm):
        if nm in app and appcols:
            try:
                idx = [i for i, c in enumerate(appcols)
                       if c.endswith("logP-3")][0]
                v = app[nm][idx]
                return f"{float(v):5.0f}" if v else "   --"
            except (IndexError, ValueError):
                pass
        return "   --"

    with PdfPages(OUT_PDF) as pdf:
        # ---------------- cover ----------------
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.95, "ATLAS12 condensates.dat: "
                 "GGchem/DustChem sources and filed fits",
                 ha="center", fontsize=14, fontweight="bold")
        n_a = sum(1 for r in rows if r["mode"] == "A")
        n_m = sum(1 for r in rows if r["mode"].startswith("M:"))
        n_pa = sum(1 for r in rows if r["mode"] == "PA")
        fig.text(0.5, 0.92,
                 f"{len(rows)} solids: {n_a} atom-referenced, "
                 f"{n_m} molecule-referenced, {n_pa} vapor-pressure metals",
                 ha="center", fontsize=11)
        conv = [
            "Filed convention (data/condensates.dat):",
            "",
            "  ln S = sum_i nu_i ln(n_i kB T) + c0/T + c1 lnT + c2"
            " + c3 T + c4 T^2",
            "",
            f"  kB = 1.380662e-16 erg/K (GGchem value); cgs; refit range"
            f" {TLO:.0f}-{THI:.0f} K.",
            "  Reference species n_i: element symbol = free neutral atom;",
            "  otherwise a molecules.dat gas-species label (mode M).",
            "  Standard-pressure constants are absorbed into c2.",
            "",
            "  mode A  = atom-referenced dG fit (GGchem fit 1/2/5)",
            "  mode PA = vapor pressure, monatomic (files the atom)",
            "  mode M  = vapor pressure, files the named gas molecule",
            "",
            "  Sources: GGchem DustChem.dat (Woitke+ 2018); liquids"
            " excluded (none",
            "  form above 1000 K in the Stage-0 reference sweeps).",
        ]
        fig.text(0.05, 0.87, "\n".join(conv), fontsize=9, va="top",
                 family="monospace")
        hdr = (f"{'species':<13}{'mode':<9}{'fit':<4}{'prec':>6}"
               f"{'RMS':>8}{'max':>8}{'T_app':>7}   source")
        tbl = [hdr, "-" * len(hdr)]
        for r in rows:
            prec = (f"{r['fit'].precision:6.3f}"
                    if r["fit"].precision else "    --")
            rms = float(np.sqrt(np.mean(r["resid"]**2)))
            tbl.append(
                f"{r['name']:<13}{r['mode']:<9}{r['fit'].code:<4}{prec:>6}"
                f"{rms:8.4f}{np.max(np.abs(r['resid'])):8.4f}"
                f"{t_app(r['name']):>7}   {r['fit'].source[:34]}")
        fig.text(0.05, 0.60, "\n".join(tbl), fontsize=8.5, va="top",
                 family="monospace")
        fig.text(0.5, 0.04,
                 "RMS/max = filed-vs-source residual [dex] over the refit "
                 "range.  T_app [K] from Stage-0 GGchem eqcond sweep at "
                 "log P[bar] = -3 (solar).\nGenerated by "
                 "tools/fit_condensates.py; rerun --validate after any "
                 "condensates.dat edit.",
                 ha="center", fontsize=8, style="italic")
        pdf.savefig(fig)
        plt.close(fig)

        # ---------------- per-species pages ----------------
        A = np.column_stack([1.0/TGRID, np.log(TGRID),
                             np.ones_like(TGRID), TGRID, TGRID**2])
        for r in rows:
            cond, fit, mode = r["cond"], r["fit"], r["mode"]
            fig, (ax_top, ax_bot) = plt.subplots(
                2, 1, figsize=(8.5, 6.5), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.08})

            y10 = r["target"]/LN10
            fit10 = (A @ r["coeff"])/LN10
            src = (fit.source.split(", +/-")[0].split(",  +/-")[0]
           .rstrip(": ") or "DustChem")[:40]
            ax_top.plot(TGRID[::5], y10[::5], "o", color="k", ms=4,
                        markerfacecolor="white", markeredgewidth=1.0,
                        label=f"source: fit{fit.code} ({src})", zorder=5)
            ax_top.plot(TGRID, fit10, "-", color="C2", lw=1.8,
                        label="filed fit (condensates.dat)", zorder=6)

            # same-reference alternative fits (atom-referenced pages only)
            alts = []
            if mode == "A":
                k = 0
                for alt in cond.fits:
                    if alt is fit or alt.code not in (1, 2, 5):
                        continue
                    try:
                        ya = np.array([target_lnk(cond, alt, mode, t)
                                       for t in TGRID])/LN10
                    except Exception:
                        continue
                    color, ls, lw = ALT_STYLES[k % len(ALT_STYLES)]
                    ax_top.plot(TGRID, ya, ls, color=color, lw=lw,
                                label=f"alt fit{alt.code} "
                                      f"({alt.source.split(chr(44))[0].rstrip(chr(58)+chr(32))[:30]})")
                    alts.append((ya, color, ls, lw,
                                 f"alt fit{alt.code} - source"))
                    k += 1

            ax_top.set_xscale("log")
            refs = refs_of(cond, mode)
            reftxt = " + ".join(
                (f"{nu} " if nu > 1 else "") + display_formula(sp)
                for nu, sp in refs)
            ax_top.set_ylabel(r"$\log_{10}\,K_{\mathrm{cond}}$  (cgs)")
            ax_top.grid(True, alpha=0.3, which="both")
            ax_top.legend(loc="best", fontsize=8)
            finite = y10[np.isfinite(y10)]
            ax_top.set_ylim(float(np.min(finite)) - 5.0,
                            float(np.max(finite)) + 5.0)
            if fit.tmax is not None and fit.tmax < THI:
                ax_top.axvline(fit.tmax, color="C2", ls=":", alpha=0.4,
                               lw=0.8)
            title = (display_formula(r["name"]) + "[s]  "
                     + rf"$\leftarrow$ {reftxt}")
            if cond.trivial:
                title += f"   ({cond.trivial.lower()})"
            ax_top.set_title(title, fontsize=12)

            ax_bot.axhline(0, color="k", lw=0.6)
            ax_bot.plot(TGRID, r["resid"], "-", color="C2", lw=1.6,
                        label="filed - source")
            for ya, color, ls, lw, lbl in alts:
                ax_bot.plot(TGRID, ya - y10, ls, color=color, lw=lw*0.8,
                            label=lbl)
            ax_bot.set_xscale("log")
            ax_bot.set_xlabel("Temperature [K]")
            ax_bot.set_ylabel(r"$\Delta\log_{10}\,K$ vs source")
            ax_bot.grid(True, alpha=0.3, which="both")
            ax_bot.legend(loc="upper right", fontsize=8)
            dmax = max(0.02, 1.5*float(np.max(np.abs(r["resid"]))))
            if alts:
                amax = max(float(np.nanmax(np.abs(ya - y10)))
                           for ya, *_ in alts)
                dmax = min(max(dmax, 1.1*amax), 20.0*dmax) if amax > dmax \
                    else dmax
            ax_bot.set_ylim(-dmax, dmax)

            rms = float(np.sqrt(np.mean(r["resid"]**2)))
            stats = (f"In fit window {TLO:.0f}-{THI:.0f} K:  "
                     f"filed-vs-source: RMS={rms:.4f}, "
                     f"max={np.max(np.abs(r['resid'])):.4f}")
            if fit.precision is not None:
                stats += f"   |   source precision +/-{fit.precision}"
            if r["note"]:
                stats += f"\n{r['note']}"
            ax_bot.text(0.02, 0.97, stats, transform=ax_bot.transAxes,
                        fontsize=7.5, va="top",
                        bbox=dict(facecolor="white", alpha=0.85,
                                  edgecolor="none"))

            for ax in (ax_top, ax_bot):
                ax.set_xlim(TLO, THI)
            ax_bot.set_xticks(XTICKS)
            ax_bot.set_xticklabels([f"{t:,d}" for t in XTICKS])
            ax_bot.xaxis.set_minor_formatter(plt.NullFormatter())

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    if "--validate" in sys.argv:
        sys.exit(0 if validate() else 1)
    rows = build()
    write_dat(rows)
    atlas(rows)
