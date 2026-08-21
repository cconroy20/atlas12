#!/usr/bin/env python3
"""Find gfall lines whose BOTH level energies map unambiguously onto a model
atom, so they can be added to mod_mklinelist's NLTE transition table.

Hand-writing that table is fine for an alkali doublet and hopeless for iron:
atom.fe607a has 548 Fe I levels, many within a few cm^-1 of one another, and a
transition assigned to the wrong level gets a plausible-looking but wrong b.
This does the matching against the atom the grid was actually solved with, and
REFUSES anything ambiguous rather than picking the nearest.

Ranking is by a crude line strength, log gf - Elo*5040/(T*8065.54) at a stated
temperature -- enough to put the lines that will actually show up first.

  python3 tools/nlte_pick_transitions.py --atom atom.fe607a --code 26.00 \\
      --window 5164 5190 [--tol 1.0] [--teff 5000] [--top 20]
"""
import argparse
import re
import sys


def read_atom(path):
    lev, started = [], False
    for ln in open(path, errors="replace"):
        s = ln.strip()
        if not s or s.startswith("*"):
            continue
        t = s.split()
        if not started:
            if len(t) >= 3 and all(x.lstrip("+-").isdigit() for x in t[:3]):
                started = True
            continue
        m = re.match(r"\s*([-\d.eEdD+]+)\s+([-\d.eEdD+]+)\s+'([^']*)'\s+(\d+)", ln)
        if not m:
            break
        lev.append((len(lev) + 1, float(m.group(1).replace("D", "E")),
                    float(m.group(2).replace("D", "E")), m.group(3).strip(),
                    int(m.group(4))))
    return lev


def match(lev, E, ion, tol, zero):
    """Atom levels of this ion within tol of E, on gfall's per-stage zero."""
    hits = [l for l in lev if l[4] == ion and abs((l[1] - zero) - E) <= tol]
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atom", required=True)
    ap.add_argument("--gfall", default="data/gfallvac08oct17.dat")
    ap.add_argument("--code", required=True, help="Kurucz species code, e.g. 26.00")
    ap.add_argument("--window", nargs=2, type=float, required=True,
                    help="vacuum wavelength range [A]")
    ap.add_argument("--tol", type=float, default=1.0, help="cm^-1 on each level")
    ap.add_argument("--teff", type=float, default=5000.0)
    ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()

    lev = read_atom(a.atom)
    ion = int(round((float(a.code) % 1) * 100)) + 1     # 26.00 -> 1, 26.01 -> 2
    # gfall zeroes each ionization stage separately; the atom runs one scale
    # across stages, so the stage's own ground level is the offset.
    zero = min(l[1] for l in lev if l[4] == ion)
    print(f"atom      : {len(lev)} levels, {sum(1 for l in lev if l[4]==ion)} "
          f"in stage {ion}; stage zero at {zero:.3f} cm^-1")

    out = []
    theta = 5040.0 / a.teff / 8065.54
    with open(a.gfall, errors="replace") as f:
        for ln in f:
            if ln[18:24].strip() != a.code:
                continue
            try:
                wl = float(ln[0:11])*10.0                # nm -> A
                gf = float(ln[11:18])
                e1 = abs(float(ln[24:36])); e2 = abs(float(ln[51:63]))
            except ValueError:
                continue
            if not (a.window[0] <= wl <= a.window[1]):
                continue
            lo, hi = min(e1, e2), max(e1, e2)
            ml = match(lev, lo, ion, a.tol, zero)
            mu = match(lev, hi, ion, a.tol, zero)
            out.append((gf - lo*theta, wl, gf, lo, hi, ml, mu))

    out.sort(reverse=True)
    print(f"lines     : {len(out)} of species {a.code} in "
          f"{a.window[0]:.1f}-{a.window[1]:.1f} A\n")
    print(f"{'wl_vac':>10}{'loggf':>8}{'Elo':>11}{'Eup':>11}  lower -> upper")
    shown = 0
    for _, wl, gf, lo, hi, ml, mu in out:
        if shown >= a.top:
            break
        if len(ml) == 1 and len(mu) == 1:
            tag = f"{ml[0][0]:4d} -> {mu[0][0]:4d}   {ml[0][3]} -> {mu[0][3]}"
        elif not ml or not mu:
            tag = f"NO MATCH  (lower {len(ml)}, upper {len(mu)})"
        else:
            tag = f"AMBIGUOUS (lower {len(ml)}, upper {len(mu)})"
        print(f"{wl:10.3f}{gf:8.3f}{lo:11.3f}{hi:11.3f}  {tag}")
        shown += 1


if __name__ == "__main__":
    main()
