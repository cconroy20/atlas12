#!/usr/bin/env python3
"""Build data/mol_broad.dat -- per-species molecular van der Waals broadening.

Replaces the two hardcoded constants (1e-8 if the upper-state label starts
with 'X', else 1e-7) that every molecular line in the codes used to get.  The
1e-7 branch -- which catches all the electronic bands that make the optical
spectrum, and all of H2O -- was 10-40x too large, giving Voigt damping
parameters of 27-44 at M-dwarf photospheric conditions where ~0.5 is physical.

Source: ExoMol `.broad` files, https://www.exomol.com/db/<SP>/<ISO>/<ISO>__H2.broad
(and __He).  Format is `code gamma n [J...]` with gamma the Lorentz HWHM in
cm^-1/bar at 296 K and n the temperature exponent, gamma(T) = gamma * (296/T)^n
per bar of perturber.  Three parameterisations appear: `a0` (keyed on J_low),
`a1` (J_up, J_low) and `m0` (keyed on |m|).  Most species publish a SINGLE a0
record -- one constant, n fixed at 0.500 -- so a per-species constant is the
state of the art rather than a compromise.  H2O, CO and AlH are J/m-resolved;
we take the median, because our H2O comes from super-line histograms that
discarded J, and because CO and AlH vary by only ~12% across J.

Usage:
    python3 tools/build_mol_broad.py --fetch     # download into a cache dir
    python3 tools/build_mol_broad.py --write     # emit data/mol_broad.dat
"""
import argparse
import os
import sys
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(REPO, "workdir", "exomol_broad")
OUT = os.path.join(REPO, "data", "mol_broad.dat")
BASE = "https://www.exomol.com/db/{sp}/{iso}/{iso}__{p}.broad"

# Kurucz molecule code -> (ExoMol species, isotopologue).  The code is what
# molec_dispatch / the readers carry per line.
SPECIES = [
    (106,   "CH",   "12C-1H"),
    (107,   "NH",   "14N-1H"),
    (108,   "OH",   "16O-1H"),
    (111,   "NaH",  "23Na-1H"),
    (112,   "MgH",  "24Mg-1H"),
    (113,   "AlH",  "27Al-1H"),
    (114,   "SiH",  "28Si-1H"),
    (120,   "CaH",  "40Ca-1H"),
    (124,   "CrH",  "52Cr-1H"),
    (126,   "FeH",  "56Fe-1H"),
    (101,   "H2",   "1H2"),
    (606,   "C2",   "12C2"),
    (607,   "CN",   "12C-14N"),
    (608,   "CO",   "12C-16O"),
    (812,   "MgO",  "24Mg-16O"),
    (814,   "SiO",  "28Si-16O"),
    (822,   "TiO",  "48Ti-16O"),
    (823,   "VO",   "51V-16O"),
    (10108, "H2O",  "1H2-16O"),
    (10820, "CaOH", "40Ca-16O-1H"),
]

# Species with no ExoMol .broad file, mapped to their nearest measured
# chemical analogue.  Each substitution is recorded in the output file so it
# is never mistaken for a measurement.
FALLBACK = {
    "CH":  ("NH", "first-row hydride"),
    "OH":  ("NH", "first-row hydride"),
    "FeH": ("CrH", "transition-metal hydride"),
    "CN":  ("CO", "first-row C-bearing diatomic"),
    "C2":  ("CO", "first-row C-bearing diatomic"),
    "H2":  ("CO", "no analogue; H2 quadrupole lines are negligible here"),
}


def fetch():
    os.makedirs(CACHE, exist_ok=True)
    for _, sp, iso in SPECIES:
        for p in ("H2", "He"):
            dst = os.path.join(CACHE, f"{sp}_{p}.broad")
            if os.path.exists(dst):
                continue
            try:
                urllib.request.urlretrieve(BASE.format(sp=sp, iso=iso, p=p), dst)
                print(f"  fetched {sp}__{p}")
            except Exception as e:
                if os.path.exists(dst):
                    os.remove(dst)
                print(f"  no {sp}__{p} ({e.__class__.__name__})")


def parse(path):
    """Median (gamma, n) over whatever parameterisation the file uses."""
    if not os.path.exists(path):
        return None
    g, n = [], []
    for line in open(path):
        r = line.split()
        if len(r) >= 3 and r[0] in ("a0", "a1", "m0"):
            g.append(float(r[1]))
            n.append(float(r[2]))
    if not g:
        return None
    g.sort()
    n.sort()
    mid = len(g) // 2
    return g[mid], n[mid], len(g)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    if args.fetch:
        fetch()
    if not args.write:
        return 0

    byname = {}
    for _, sp, _ in SPECIES:
        r2 = parse(os.path.join(CACHE, f"{sp}_H2.broad"))
        re_ = parse(os.path.join(CACHE, f"{sp}_He.broad"))
        if r2:
            byname[sp] = (r2, re_)

    with open(OUT, "w") as fh:
        fh.write(
            "# Molecular van der Waals broadening, from ExoMol .broad files.\n"
            "# Built by tools/build_mol_broad.py; see its header for provenance.\n"
            "#\n"
            "# gamma is the Lorentz HWHM in cm^-1/bar at 296 K, n the exponent\n"
            "# in gamma(T) = gamma*(296/T)^n per bar of perturber.  nJ is the\n"
            "# number of records in the source file (1 = a single constant is\n"
            "# all ExoMol publishes; >1 = J- or m-resolved, median taken).\n"
            "#\n"
            "# 'src' names the species the numbers came from: its own name for a\n"
            "# measurement, or the analogue used when ExoMol has no file for it.\n"
            "#\n"
            "# code species   gamma_H2    n_H2  gamma_He    n_He   nJ  src\n")
        for code, sp, _ in SPECIES:
            src, note = sp, ""
            if sp not in byname:
                src, note = FALLBACK[sp]
                if src not in byname:
                    print(f"  WARNING: no data and no usable fallback for {sp}")
                    continue
            (g2, n2, nj), he = byname[src]
            ge, ne = (he[0], he[1]) if he else (g2 * 0.55, n2)
            tag = src if not note else f"{src}  # fallback: {note}"
            fh.write(f"{code:6d} {sp:<8s} {g2:9.4f} {n2:7.3f} "
                     f"{ge:9.4f} {ne:7.3f} {nj:4d}  {tag}\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
