# SYNTHE/ATLAS F90 Upgrade To-Do List

Open items only; completed work is recorded in
[CHANGELOG.md](CHANGELOG.md).

---

## 1. H2 Collision-Induced Absorption

**Status: DONE** (see CHANGELOG, "Numerical and physical").  H2COLLOP now
carries Abel et al. (2011/2012) H2-H2 and H2-He from HITRAN's recommended
set, plus the H2-H and H-He pairs that were never in the Kurucz lineage,
on a 25 cm^-1 / ~100 K grid, with the stimulated-emission double count
removed.  Rebuild the table with tools/build_cia_table.py.

Two residual limitations, neither blocking:

- Abel's H2-H2 calculation stops at 3000 K and 10^4 cm^-1.  Beyond it the
  table continues Borysow, Jorgensen & Fu (2001) rescaled by the measured
  Abel/BJF01 trend, frozen at 5000 K.  No public calculation covers that
  region; if one appears, drop the continuation.  Measured sensitivity at
  2700 K: 0.11% on the K continuum between the two defensible choices,
  0.71% for a deliberate factor of two.
- H2-H is tabulated only over 1000-2500 K and clamps outside, which
  underestimates it above 2500 K where n(H2)*n(H) actually peaks.  The
  pair is worth 0.17% of the K continuum at 2700 K, so this is small, but
  it would grow toward 3000-4000 K.

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

## 3. Low-Temperature Gas-Solve Hardening

**Status: open; blocks lowering TFLOOR_ATM below 1200 K.**

The molecular-equilibrium Newton diverges when a layer reaches
T ~ 1050 K: the equilibrium n_e (~1e2 cm^-3) against 1e12+ atom
densities and K_eq ~ e^100+ spans ~90 orders of magnitude and the
charge row is lost in double precision (the NaN then spreads through
the depth warm-start).  GGchem handles this regime in quadruple
precision below ~1000 K.  Candidate fixes: log-space formulation,
extended precision for the charge row, or charge-row preconditioning.

---

## 4. Dust Opacity (Dusty Mode)

**Status: open; only relevant if the grid floor drops below ~2400 K.**

The Cond-limit solver depletes the gas phase but adds no grain opacity
(the AMES-Cond limit).  Below ~2400 K the condensate column becomes
optically significant and a Dusty-mode treatment (grain opacity from
the condensed fractions) would be needed.  Settling/microphysics
(BT-Settl/StaticWeather style) is brown-dwarf territory -- out of scope.

---

## 5. Mann-Validation Open Items (post D0-audit + resolution fixes, 2026-08-09)

**Status: open.**  At converged synthesis resolution (R=300k) with the
corrected molecular network, the RE GJ887 model matches the Mann spectrum
at the percent level; what survives:

- **CaH bands ~10% too strong** (CaH2 index 0.587 vs 0.659): no plausible
  T-structure fixes it without breaking TiO; candidates are the Owens+22
  CaH gf scale and the D0 tail (1.699 +/- ~0.02 eV = 8% lever).  Test:
  A/B the old CaH line list and D0 = 1.679.
- **4700-5400 A: both ATLAS12 and PHOENIX ~10-15% BRIGHTER than data**
  (shared -> common missing opacity or SNIFS blue calibration).
- **8000-8800 A: models few % low; 9100-9600 few % high** (shared with
  PHOENIX; SNIFS red edge / telluric splice caveats).
- **Ca II IRT cores too deep** (LTE; known).
- **~4100 K limit cycle worsened**: PM_I10113+4927 reconverged at 18% RMS
  (was ~2.9%) on the corrected network -- rerun the limit-cycle diagnosis.
- **Solar absolute (+2.1%) recheck at R=300k**; then the sample-wide
  rerun (all medians shift with the resolution fix).
