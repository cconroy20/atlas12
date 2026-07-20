#!/usr/bin/env python3
"""Loader for GGchem's condensate thermochemistry (data/DustChem.dat).

Parses every entry (active fit + commented alternatives) and evaluates
ln S contributions in GGchem's own conventions (supersat.f, Woitke+ 2018):

  S = exp( sum_j nu_j * ln(n_at,j * kT / p_std) + lnK(T) )

where lnK(T) = -dG_f/RT for formation of the solid from FREE NEUTRAL GAS
ATOMS at standard pressure p_std.  Fit forms:

  1  Sharp & Huebner 1990:  dG[cal/mol] = a0/T + a1 + a2 T + a3 T^2 + a4 T^3,
     p_std = 1 atm
  2  NIST-JANAF dG:         dG[J/mol]   = a0/T + a1 + a2 T + a3 T^2 + a4 T^3,
     p_std = 1 bar
  3  ln p_vap fit:          ln psat[dyn/cm2] = a0/T + a1 + a2 T + a3 T^2
     + a4 T^3;  S = n1*kT/psat with n1 = atom (monatomic) or gas molecule
  5  Stock:                 -dG/RT = a0/T + a1 lnT + a2 + a3 T + a4 T^2,
     p_std = 1 bar (linear dG extrapolation above Tmax, as in supersat.f)

Fits 4/6/7/99 are parsed but not evaluated (none are needed for the warm
condensate set).  Units: cgs; kB in erg/K.
"""
import math
import os
import re

KB = 1.380662e-16          # erg/K   (GGchem datamod.f value)
RGAS = 8.3144598           # J/mol/K (GGchem value)
CAL = 4.184                # J/cal
BAR = 1.0e6                # dyn/cm2
ATM = 1.013e6              # dyn/cm2 (GGchem value)

DEFAULT_PATH = os.path.expanduser(
    "~/kurucz/upgrade/raw_data/ggchem/data/DustChem.dat")


class CondFit:
    def __init__(self, code, coeff, source, precision, tmax, active):
        self.code = code            # GGchem fit form (int)
        self.coeff = coeff          # list of floats (5 for codes != 7)
        self.source = source        # provenance text scraped from comments
        self.precision = precision  # +/- dex, if stated
        self.tmax = tmax            # stated Tmax of fit validity (K)
        self.active = active        # True if this is GGchem's active fit

    def __repr__(self):
        tag = "active" if self.active else "alt"
        return f"<fit{self.code} {tag} src='{self.source}'>"


class Condensate:
    def __init__(self, name, trivial, rho, stoich, tmelt=None):
        self.name = name            # e.g. 'Al2O3[s]'
        self.trivial = trivial      # e.g. 'CORUNDUM(alpha)'
        self.rho = rho              # g/cm3
        self.stoich = stoich        # list of (count, element) tuples
        self.tmelt = tmelt          # K, liquids only
        self.fits = []              # CondFit list, active first

    @property
    def active_fit(self):
        return self.fits[0] if self.fits and self.fits[0].active else None

    def natoms(self):
        return sum(n for n, _ in self.stoich)

    def __repr__(self):
        return f"<{self.name} {self.stoich} {len(self.fits)} fits>"


def _parse_fit_line(line, active):
    """Parse a '<code> <c0> <c1> ...' fit line (comment marker stripped)."""
    parts = line.split("!")[0].split()
    try:
        code = int(parts[0])
        coeff = [float(x) for x in parts[1:]]
    except (ValueError, IndexError):
        return None
    if not coeff:
        return None
    return CondFit(code, coeff, "", None, None, active)


def load_dustchem(path=DEFAULT_PATH):
    """Return {name: Condensate} for every entry in DustChem.dat."""
    with open(path) as f:
        lines = f.read().splitlines()
    imax = int(lines[2].split()[0])
    out = {}
    i = 4
    nread = 0
    while nread < imax and i < len(lines):
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break
        head = lines[i].split()
        name = head[0]
        trivial = head[1] if len(head) > 1 else ""
        tmelt = None
        if name.endswith("[l]") and len(head) > 2:
            try:
                tmelt = float(head[2])
            except ValueError:
                pass
        i += 1
        rho = float(lines[i].split()[0]); i += 1
        nel = int(lines[i].split()[0]); i += 1
        stoich = []
        for _ in range(nel):
            nu = int(lines[i][0:2])
            el = lines[i][3:5].strip()
            stoich.append((nu, el))
            i += 1
        cond = Condensate(name, trivial, rho, stoich, tmelt)

        # comment block: alternative fits + provenance; first non-# line
        # is the active fit; blank line ends the entry
        pending_note = ""
        got_active = False
        while i < len(lines) and lines[i].strip():
            ln = lines[i]
            if ln.lstrip().startswith("#"):
                body = ln.lstrip()[1:].strip()
                fit = _parse_fit_line(body, active=False)
                if fit is not None:
                    fit.source = pending_note
                    _scrape_meta(fit, pending_note)
                    cond.fits.append(fit)
                else:
                    pending_note = body
            elif not got_active:
                fit = _parse_fit_line(ln, active=True)
                if fit is not None:
                    fit.source = pending_note
                    _scrape_meta(fit, pending_note)
                    cond.fits.insert(0, fit)
                    got_active = True
            i += 1
        out[name] = cond
        nread += 1
    return out


def _scrape_meta(fit, note):
    m = re.search(r"\+/-\s*([0-9.]+)", note)
    if m:
        fit.precision = float(m.group(1))
    m = re.search(r"Tmax\s*=\s*([0-9.]+)", note)
    if m:
        fit.tmax = float(m.group(1))


# ---------------------------------------------------------------------------
# evaluation, GGchem conventions

def lnk_atoms(cond, fit, T):
    """ln K(T) = -dG_f(atoms->solid)/RT with GGchem's p_std convention
    folded OUT (returned separately): returns (lnK, p_std).

    S = exp( sum nu ln(n_at kT / p_std) + lnK ).  Only fit codes 1,2,5.
    """
    a = fit.coeff
    if fit.code == 1:
        dG = a[0]/T + a[1] + a[2]*T + a[3]*T**2 + a[4]*T**3   # cal/mol
        return -dG/(RGAS/CAL*T), ATM
    if fit.code == 2:
        dG = a[0]/T + a[1] + a[2]*T + a[3]*T**2 + a[4]*T**3   # J/mol
        return -dG/(RGAS*T), BAR
    if fit.code == 5:
        tmax = fit.tmax or 6000.0
        if T > tmax:
            # linear dG extrapolation beyond Tmax (supersat.f:119-133)
            t1 = tmax
            f0 = a[0]/t1 + a[1]*math.log(t1) + a[2] + a[3]*t1 + a[4]*t1**2
            df = -a[0]/t1**2 + a[1]/t1 + a[3] + 2.0*a[4]*t1
            dg0 = f0*t1
            dgdt = f0 + t1*df
            return (dg0 + dgdt*(T - t1))/T, BAR
        return a[0]/T + a[1]*math.log(T) + a[2] + a[3]*T + a[4]*T**2, BAR
    raise ValueError(f"fit code {fit.code} is not an atom-referenced dG form")


def ln_psat(fit, T):
    """ln p_sat [dyn/cm2] for fit code 3."""
    if fit.code != 3:
        raise ValueError("ln_psat needs fit code 3")
    a = fit.coeff
    return a[0]/T + a[1] + a[2]*T + a[3]*T**2 + a[4]*T**3


def log10_S(cond, fit, T, nat, nmol=None):
    """log10 supersaturation, exactly as GGchem's SUPERSAT computes it.

    nat  : {element: number density [cm^-3]} of free neutral atoms
    nmol : {molecule: number density} -- required for fit-3 polyatomics
           (key = condensate name without the [s]/[l] suffix, e.g. 'ZrO2')
    """
    if fit.code == 3:
        if len(cond.stoich) == 1 and cond.stoich[0][0] == 1:
            n1 = nat[cond.stoich[0][1]]
        else:
            key = cond.name.split("[")[0]
            if nmol is None or key not in nmol:
                raise KeyError(f"fit-3 species {cond.name} needs gas "
                               f"molecule '{key}'")
            n1 = nmol[key]
        return (math.log(n1*KB*T) - ln_psat(fit, T))/math.log(10.0)
    lnk, pst = lnk_atoms(cond, fit, T)
    s = lnk
    for nu, el in cond.stoich:
        s += nu*math.log(nat[el]*KB*T/pst)
    return s/math.log(10.0)


if __name__ == "__main__":
    db = load_dustchem()
    print(f"{len(db)} condensates parsed")
    for nm in ("Al2O3[s]", "CaTiO3[s]", "ZrO2[s]", "Fe[s]", "Mg2SiO4[s]",
               "VO[s]", "V2O3[s]", "SiO[s]", "Ti3O5[s]"):
        c = db.get(nm)
        if c is None:
            print(f"{nm:<14} MISSING")
            continue
        af = c.active_fit
        print(f"{nm:<14} rho={c.rho:5.2f} {str(c.stoich):<40} "
              f"active=fit{af.code} alts={len(c.fits)-1} "
              f"prec={af.precision} Tmax={af.tmax}")
