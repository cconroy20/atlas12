# SYNTHE/ATLAS F90 Upgrade To-Do List

---

## 1. H2 Collision-Induced Absorption

**Status: open.**

**Current implementation (H2COLLOP):** Tabulated log(alpha) on a
7-temperature x 81-wavenumber grid from Borysow, Jorgensen & Zheng (1997),
interpolated bilinearly.

**Modern alternative:** Abel et al. (2011, J. Phys. Chem. A 115, 6805;
2012, J. Chem. Phys. 136, 044319) provide updated H2-H2 and H2-He CIA
coefficients for 200-9900 K.  The HITRAN CIA database (Richard et al. 2012)
compiles these on a fine grid.  Differences from Borysow et al. can reach
10-30% in specific spectral windows.

**Impact:** Critical for cool stars (Teff < 5000 K), brown dwarfs, and
giant planets, where H2 CIA dominates the IR pseudo-continuum.  For
solar-type stars the impact is small.

---

## 2. vdW -> ABO Transition for Line Broadening

**Status: open.**

**Current implementation:** Classical Unsold formula for van der Waals
broadening of neutral-atom lines.

**Modern alternative:** Anstee, Barklem & O'Mara (ABO) theory provides
tabulated sigma and alpha coefficients that correctly describe the
broadening of neutral-atom lines by H collisions, fitting to ab initio
atomic-structure calculations.  Where ABO data exist they should be used
in preference to Unsold, which is known to systematically under-predict
the broadening of strong metal lines by factors of 2-3.

**Impact:** Most pronounced for strong lines of neutral metals (e.g.,
Na I D, Ca II H and K, Fe I lines) used in abundance analyses.  Moving
from Unsold to ABO can shift derived abundances by 0.1-0.3 dex for
saturated lines.

---

## 3. Molecular Equilibrium Network Overhaul

**Status: COMPLETED July 18-19, 2026** (commits fecddd0..0056192).

Network grew 190 rows / 23 equations -> 297 rows / 33 equations.  All 96
diatomics validate against BC16 in the correct conventions
(tools/fit_molecule_keq.py --validate); polyatomics are ExoMol-fit,
JANAF/SK-referenced, or deliberately reconciled (MgOH via Koput D0, HO2
via ATcT).  Fixes along the way: READMOL integer decoder (La exposed it),
electron budget closed (Co/Ni/Cu/Zn -- solar photospheres cool 1-3 K vs
all prior models), Y/Zr/Co ionization potentials, molecular-ion Saha
convention (8 rows were 4-13 dex low since April), H3+ (Kurucz original
was 4.5 dex low, missing one T^3/2).  Headline physics: TiO2 takes 23%
of Ti at 2500 K (TiO bands -12%), HBO holds ~99% of B, s-process oxides
lock Zr/Y/La.  Full pipeline + fit atlas (comp_pf.pdf, 138 pp) live in
tools/; regenerate after any molecules.dat change.

---

## 4. Equilibrium Condensation (Cond-limit)

**Status: COMPLETED July 20-21, 2026** (commits d60069d..5e3a493);
**ON by default** (`USE_CONDENSATION = .TRUE.`).

Cond-limit solver inside NMOLEC (depletion only, no grain opacity):
21 condensates in data/condensates.dat refit from the GGchem
compilation (tools/fit_condensates.py --validate), outer quasi-Newton
active-set iteration wrapped around the untouched legacy gas solve.
GGchem point-matched validation: phase assemblages exact, major-element
gas abundances to <=0.04 dex down to 1200 K surfaces.  A/B: 2500-2800 K
dwarf surfaces cool 25-42 K over the condensing layers (same class as
the MSG ~25 K benchmark); solar/3500 K bit-identical to flag-OFF.
Full ledger in CHANGELOG.md.

---

## 5. Condensate Species Extension

**Status: open.**

The two species-list gaps identified by the GGchem validation at 1200 K
surfaces: MgCr2O4 (no Cr condensates are filed at all) and ZrSiO4
(zircon, which supersedes baddeleyite ZrO2 at low T).  They account for
the Si (-0.18 dex) and trace-Zr deviations vs GGchem.  Add both to
data/condensates.dat via tools/fit_condensates.py and rerun --validate.

---

## 6. Low-Temperature Gas-Solve Hardening

**Status: open; blocks lowering TFLOOR_ATM below 1200 K.**

The molecular-equilibrium Newton diverges when a layer reaches
T ~ 1050 K: the equilibrium n_e (~1e2 cm^-3) against 1e12+ atom
densities and K_eq ~ e^100+ spans ~90 orders of magnitude and the
charge row is lost in double precision (the NaN then spreads through
the depth warm-start).  GGchem handles this regime in quadruple
precision below ~1000 K.  Candidate fixes: log-space formulation,
extended precision for the charge row, or charge-row preconditioning.

---

## 7. Dust Opacity (Dusty Mode)

**Status: open; only relevant if the grid floor drops below ~2400 K.**

The Cond-limit solver depletes the gas phase but adds no grain opacity
(the AMES-Cond limit).  Below ~2400 K the condensate column becomes
optically significant and a Dusty-mode treatment (grain opacity from
the condensed fractions) would be needed.  Settling/microphysics
(BT-Settl/StaticWeather style) is brown-dwarf territory -- out of scope.

---
