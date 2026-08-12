#!/usr/bin/env python3
"""Build data/mol_broad.dat -- per-species molecular van der Waals broadening.

Replaces the two hardcoded constants (1e-8 if the upper-state label starts
with 'X', else 1e-7) that every molecular line in the codes used to get.  The
1e-7 branch -- which catches all the electronic bands that make the optical
spectrum, and all of H2O -- was 10-40x too large, giving Voigt damping
parameters of 27-44 at M-dwarf photospheric conditions where ~0.5 is physical.

Source: ExoMol `.broad` files, https://www.exomol.com/db/<SP>/<ISO>/<ISO>__H2.broad
(and __He).  Format is `code gamma n [J...]` with gamma the Lorentz HWHM in
cm^-1/ATM at 296 K and n the temperature exponent, gamma(T) = gamma * (296/T)^n
per atm of perturber (2024 release paper Sect. 4: P_ref = 1 atm, T_ref = 296 K).

Provenance matters and is recorded per species: ExoMol has EMPIRICAL data for
only seven molecules (H2O, NH3, SO2, CH4, PH3, HCN, H2CO).  Everything else --
including TiO, CaH, VO and CaOH -- carries THEORETICAL values, and for the
'exotic' pairs those are simple semi-classical, rotationally-independent
estimates from molecular masses and kinetic diameters, with n = 0.5 assumed.
That is why most of our species show a single record with n exactly 0.500 --
it is an assumed default, not a fitted exponent.

Three parameterisations appear: `a0` (keyed on J_low), `a1` (J_up, J_low) and
`m0` (keyed on |m|).  H2O, CO and AlH are J/m-resolved; we take the median,
because our H2O comes from super-line histograms that discarded J, and because
CO and AlH vary by only ~12% across J.

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

# Gharib-Nezhad et al. (2021, ApJS 254, 34; the EXOPLINES database) Table 1:
# gamma_L(J) = gamma_0 - (dgamma/dJ) * J_lower, for J_lower the LOWER-state
# total angular momentum, at 296 K.  Where a species appears here we adopt
# their gamma_0 and slope in preference to ExoMol's rotationally-independent
# estimate, because a single value is simply wrong for a molecule whose
# populated J is high.  Their J=0 values agree with ExoMol's to ~13% anyway
# (TiO 0.10 vs 0.0882; H2O 0.09/0.02 vs 0.0916/0.0219), so this changes the
# J-dependence, not the scale.
#
# They floor the decline to avoid negative widths at high J (their Eq. 5,
# gamma_L(J=40) ~ 0.1 gamma_L(J=0)); we apply the same 0.1 floor.
EXOPLINES = {           # species: (g0_H2, dgdJ_H2, g0_He, dgdJ_He)
    "TiO":  (0.100, 0.00200, 0.060, 0.00120),
    "VO":   (0.100, 0.00200, 0.060, 0.00120),
    "FeH":  (0.075, 0.00100, 0.045, 0.00060),
    "CrH":  (0.075, 0.00100, 0.045, 0.00060),
}
JFLOOR = 0.1            # gamma_L(J) >= JFLOOR * gamma_L(0)

# Rotational constants B [cm^-1], used ONLY to population-weight gamma(J) for
# the species whose line lists carry no per-line J (TiO, H2O, CaOH -- all
# binary/super-line formats).  H2O is an asymmetric top; 14.5 is an effective
# value adequate for a Boltzmann weight.  The rest are spectroscopic B_e.
BROT = {
    "CH": 14.19, "NH": 16.34, "OH": 18.55, "NaH": 4.90, "MgH": 5.82,
    "AlH": 6.39, "SiH": 7.50, "CaH": 4.23, "CrH": 6.22, "FeH": 6.55,
    "H2": 60.85, "C2": 1.82, "CN": 1.90, "CO": 1.93, "MgO": 0.574,
    "SiO": 0.727, "TiO": 0.5355, "VO": 0.5463, "H2O": 14.5, "CaOH": 0.334,
}

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
    """Return (gamma[], n[], J[]) for whatever parameterisation the file uses."""
    if not os.path.exists(path):
        return None
    g, n, J = [], [], []
    for line in open(path):
        r = line.split()
        if len(r) >= 3 and r[0] in ("a0", "a1", "m0"):
            g.append(float(r[1]))
            n.append(float(r[2]))
            # a0/m0: J (or |m|) is the last column; a1: J_lower is the last
            J.append(float(r[-1]))
    if not g:
        return None
    return g, n, J


def collapse(rec, B, tref=3000.0):
    """Reduce a .broad file to (g0, n, dgdJ, nrec).

    For a single-record file there is nothing to do.  For a J-resolved file we
    fit the linear g0 - dgdJ*J form over the J range that actually matters at
    stellar temperatures -- out to twice the Boltzmann peak -- rather than
    taking a median over the J grid.  Medianing the grid is what we did first
    and it is wrong: it returns the value at the middle of the tabulated J
    range (J ~ 25 for water) when the populated J is ~8, and for H2O that
    understated gamma by a factor of two.
    """
    g, n, J = rec
    if len(g) == 1:
        return g[0], n[0], 0.0, 1
    import math
    g = [float(x) for x in g]
    J = [float(x) for x in J]
    jpk = math.sqrt(tref / (2.0 * B * 1.4388)) - 0.5 if B > 0 else 10.0
    jmax = max(2.0 * jpk, 4.0)
    sel = [(jj, gg) for jj, gg in zip(J, g) if jj <= jmax]
    if len(sel) < 2:
        sel = list(zip(J, g))[:2]
    sx = sum(j for j, _ in sel); sy = sum(v for _, v in sel)
    sxx = sum(j * j for j, _ in sel); sxy = sum(j * v for j, v in sel)
    nn = len(sel)
    den = nn * sxx - sx * sx
    slope = (nn * sxy - sx * sy) / den if den else 0.0
    icept = (sy - slope * sx) / nn
    nmed = sorted(n)[len(n) // 2]
    return icept, nmed, max(-slope, 0.0), len(g)


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
        B = BROT.get(sp, 1.0)
        r2 = parse(os.path.join(CACHE, f"{sp}_H2.broad"))
        re_ = parse(os.path.join(CACHE, f"{sp}_He.broad"))
        if r2:
            byname[sp] = (collapse(r2, B), collapse(re_, B) if re_ else None)

    with open(OUT, "w") as fh:
        fh.write(
            "# Molecular van der Waals broadening, with J-dependence.\n"
            "# Built by tools/build_mol_broad.py; see its header for provenance.\n"
            "#\n"
            "# gamma_L(J) = max(g0 - dgdJ*J_lower, %.2f*g0), in cm^-1/ATM HWHM at\n"
            "# 296 K, with gamma(T) = gamma*(296/T)^n per atm of perturber\n"
            "# (ExoMol 2024 release Sect. 4: P_ref = 1 atm, T_ref = 296 K).\n"
            "# J_lower is the LOWER-state total angular momentum.\n"
            "#\n"
            "# The J-dependence and the g0 that goes with it are Gharib-Nezhad\n"
            "# et al. (2021, ApJS 254, 34) Table 1 where that paper covers the\n"
            "# species; dgdJ = 0 means no J-dependence is known and the ExoMol\n"
            "# rotationally-independent estimate is used unchanged.\n"
            "#\n"
            "# B is the rotational constant, used ONLY to population-weight\n"
            "# gamma(J) for line lists that carry no per-line J (TiO, H2O, CaOH).\n"
            "#\n"
            "# PROVENANCE: ExoMol has EMPIRICAL H2/He broadening for only seven\n"
            "# molecules (H2O, NH3, SO2, CH4, PH3, HCN, H2CO) -- of ours, only\n"
            "# H2O.  Everything else is theoretical: ExoMol's are semi-classical\n"
            "# estimates from molecular masses and kinetic diameters with n = 0.5\n"
            "# assumed, and Gharib-Nezhad et al. note there is not a single\n"
            "# laboratory measurement for any of the metal oxides/hydrides here.\n"
            "# Their J-slopes are likewise assumed, by analogy with HCl for\n"
            "# large-dipole species.  Treat all of it as good to a factor of a\n"
            "# few, not to percent.\n"
            "#\n"
            "# 'src' names where the numbers came from.\n"
            "#\n"
            "# code species   g0_H2    n_H2   g0_He    n_He  dgdJ_H2  dgdJ_He"
            "       B   nJ  src\n" % JFLOOR)
        for code, sp, _ in SPECIES:
            src, note = sp, ""
            if sp not in byname:
                src, note = FALLBACK[sp]
                if src not in byname:
                    print(f"  WARNING: no data and no usable fallback for {sp}")
                    continue
            (g2, n2, d2, nj), he = byname[src]
            ge, ne, de = (he[0], he[1], he[2]) if he else (g2 * 0.55, n2, d2 * 0.55)
            if sp in EXOPLINES:                 # prefer the J-resolved source
                g2, d2, ge, de = EXOPLINES[sp]
                note = "g0+dgdJ from GN+21 Tab.1"
                src = sp        # GN+21 covers it directly; no analogue needed
            b = BROT.get(sp, 0.0)
            tag = src if not note else f"{src}  # {note}"
            fh.write(f"{code:6d} {sp:<8s} {g2:8.4f} {n2:7.3f} {ge:8.4f} {ne:7.3f}"
                     f" {d2:8.5f} {de:8.5f} {b:7.3f} {nj:4d}  {tag}\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
