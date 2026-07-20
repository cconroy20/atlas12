#!/usr/bin/env python3
"""Stage-3 condensation validation: ATLAS12 Cond-limit vs GGchem eqcond,
point-matched per layer.

For every condensing layer of a flag-ON ATLAS12 run (plus boundary
buffers), runs GGchem equilibrium condensation at the layer's exact
(T, P) with the model's own abundances, and compares:
  - the active phase assemblage,
  - the depleted gas-phase element abundances (dex).

Inputs: the run's stdout (COND:/CONDEPS: lines; last iteration block),
its .iter file (T, P per layer), and the input .atm (abundance table).

Usage:
  python3 validate_condensation.py --stdout run.stdout --iter run.iter \
      --atm model.atm [--out report.txt]

Known systematics (see condensation memory / comp_cond.pdf):
  - ZrO2 is filed from SUPCRTBL while GGchem's default is the
    vapor-pressure fit: ~0.1 dex in S, ~10-15 K in appearance T.
  - W[s] is in GGchem but deliberately not in ATLAS12 (trace, no
    opacity); it is ignored in the phase comparison.
"""
import argparse
import math
import os
import re
import shutil
import subprocess

GGCHEM_DIR = os.path.expanduser("~/kurucz/upgrade/raw_data/ggchem")
ELEMENTS_LINE = ("H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li F P V Cr "
                 "Mn Ni Zr W el")
COMPARE_ELS = ["Ti", "Al", "Ca", "Mg", "Si", "Fe", "V", "Zr", "Ni"]
IGNORE_PHASES = {"W"}


def parse_atlas_stdout(path):
    """Last COND block: {J: {'phases': {name: dens}, 'eps': {el: frac}}}."""
    cond, eps = {}, {}
    for ln in open(path):
        m = re.match(r"\s*COND: J=\s*(\d+)\s+T=\s*([\d.]+)\s*(.*)", ln)
        if m:
            j = int(m.group(1))
            if j in cond and j == min(cond):   # new block starts
                cond, eps = {}, {}
            phases = dict((p, float(v)) for p, v in
                          re.findall(r"(\S+)=\s*([\dEe.+-]+)", m.group(3)))
            cond[j] = {"T": float(m.group(2)), "phases": phases}
        m = re.match(r"\s*CONDEPS: J=\s*(\d+)\s*(.*)", ln)
        if m:
            j = int(m.group(1))
            eps[j] = dict((e, float(v)) for e, v in
                          re.findall(r"(\S+)=\s*([\dEe.+-]+)", m.group(2)))
    for j in cond:
        cond[j]["eps"] = eps.get(j, {})
    return cond


def parse_iter(path):
    """Final block: {J: (T, P_dyn)}."""
    lines = open(path).readlines()
    idx = [i for i, ln in enumerate(lines) if "log10TAU" in ln]
    out = {}
    for ln in lines[idx[-1] + 2:]:
        p = ln.split()
        if len(p) < 15 or not p[0].isdigit():
            continue
        out[int(p[0])] = (float(p[2]), float(p[13]))
    return out


def parse_atm_abund(path):
    """{el: log10 number fraction}; H/He converted from linear."""
    out = {}
    for ln in open(path):
        if "ABUNDANCE" in ln or ln.startswith("TEFF") or "TITLE" in ln:
            continue
        for z, el, val in re.findall(
                r"(\d+)([A-Za-z ]{2})\s*(-?\d+\.\d+)", ln):
            el = el.strip()
            v = float(val)
            out[el] = math.log10(v) if v > 0 else v
        if "READ DECK" in ln or "TAU5000" in ln:
            break
    return out


def write_ggchem_inputs(abund, workdir):
    """Custom abundance file (12-scale) from the model's abundances."""
    fh = 10 ** abund["H"]
    apath = os.path.join(workdir, "abund_atlas.dat")
    with open(apath, "w") as f:
        for el in ELEMENTS_LINE.split():
            if el == "el":
                continue
            la = abund.get(el)
            if la is None:
                continue
            f.write(f"{el:<3} {12.0 + la - math.log10(fh):8.3f}\n")
    return apath


def ggchem_point(T, P_dyn, apath, tag):
    """Run one eqcond point; return (phases set, {el: log10 eps_gas})."""
    pbar = P_dyn / 1.0e6
    inp = os.path.join(GGCHEM_DIR, "input", f"atlas_val_{tag}.in")
    with open(inp, "w") as f:
        f.write(f"""# selected elements
{ELEMENTS_LINE}

dispol_BarklemCollet.dat                ! dispol_file
dispol_StockKitzmann_withoutTsuji.dat   ! dispol_file2
dispol_WoitkeRefit.dat                  ! dispol_file3

0                     ! abund_pick
{apath}

.true.                ! model_eqcond

1                     ! model_dim  (0,1,2)
.true.                ! model_pconst
{T + 0.02:.2f}                ! Tmax [K]
{T:.2f}                ! Tmin [K]      (if model_dim>0)
{pbar:.6E}          ! pmax [bar]    (if pconst=.true.)
{pbar:.6E}          ! pmin [bar]
4.E+19                ! nHmax [cm-3]  (if pconst=.false.)
4.E+19                ! nHmin [cm-3]
2                     ! Npoints
5                     ! NewBackIt
1000.0                ! Tfast
""")
    db = os.path.join(GGCHEM_DIR, "database.dat")
    if os.path.exists(db):
        os.remove(db)
    r = subprocess.run(["./ggchem", f"input/atlas_val_{tag}.in"],
                       cwd=GGCHEM_DIR, capture_output=True, text=True,
                       timeout=600)
    if r.returncode != 0:
        return None, None
    ls = open(os.path.join(GGCHEM_DIR, "Static_Conc.dat")).readlines()
    nout, nmole, ndust, _ = map(int, ls[1].split())
    hdr = ls[2].split()
    row = [float(x) for x in ls[-1].split()]
    n0 = 3 + 1 + nout + nmole + ndust
    phases = set(hdr[i][1:] for i in range(n0, n0 + ndust)
                 if row[i] > -299) - IGNORE_PHASES
    epsc = {}
    for k, c in enumerate(hdr):
        if c.startswith("eps"):
            epsc[c[3:]] = row[k]          # log10 n_el,gas / n<H>
    return phases, epsc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stdout", required=True)
    ap.add_argument("--iter", dest="iterfile", required=True)
    ap.add_argument("--atm", required=True)
    ap.add_argument("--out", default="s3_report.txt")
    args = ap.parse_args()

    cond = parse_atlas_stdout(args.stdout)
    tp = parse_iter(args.iterfile)
    abund = parse_atm_abund(args.atm)
    fh = 10 ** abund["H"]

    jmin, jmax = min(cond), max(cond)
    layers = [j for j in sorted(tp) if jmin - 2 <= j <= jmax + 2]

    lines = [f"{'J':>3} {'T':>7} {'logP':>6}  "
             f"{'phases: atlas | ggchem':<58} {'set':>5}"]
    devs = []
    for j in layers:
        T, P = tp[j]
        ours = set(cond.get(j, {}).get("phases", {}))
        gg, epsg = ggchem_point(T, P, write_ggchem_inputs(
            abund, os.path.dirname(os.path.abspath(args.out)) or "."),
            f"J{j}")
        if gg is None:
            lines.append(f"{j:>3} {T:7.1f} {'':>6}  GGCHEM FAILED")
            continue
        setok = "OK" if ours == gg else "DIFF"
        lines.append(f"{j:>3} {T:7.1f} {math.log10(P/1e6):6.2f}  "
                     f"{'+'.join(sorted(ours)) or '-':<28}| "
                     f"{'+'.join(sorted(gg)) or '-':<28} {setok:>5}")
        gasfrac = cond.get(j, {}).get("eps", {})
        row = []
        for el in COMPARE_ELS:
            if el not in gasfrac and el not in (epsg or {}):
                continue
            la = abund.get(el)
            if la is None:
                continue
            ours_log = (la - math.log10(fh)
                        + math.log10(max(gasfrac.get(el, 1.0), 1e-30)))
            d = ours_log - epsg[el]
            if abs(d) > 0.005 or el in gasfrac:
                row.append(f"{el} {d:+6.2f}")
                devs.append((abs(d), j, el))
        if row:
            lines.append(f"     d(log eps_gas): {'  '.join(row)}")

    devs.sort(reverse=True)
    lines.append("")
    lines.append(f"worst eps deviations: " + "  ".join(
        f"{el}@J{j}:{d:.2f}" for d, j, el in devs[:6]))
    report = "\n".join(lines)
    open(args.out, "w").write(report + "\n")
    print(report)


if __name__ == "__main__":
    main()
