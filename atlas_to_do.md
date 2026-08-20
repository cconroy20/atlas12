# SYNTHE/ATLAS F90 Upgrade To-Do List

Open items only; completed work is recorded in
[CHANGELOG.md](CHANGELOG.md).

---

## 1. vdW -> ABO Transition for Line Broadening

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

## 2. Dust Opacity (Dusty Mode)

**Status: open; only relevant if the grid floor drops below ~2400 K.**

The Cond-limit solver depletes the gas phase but adds no grain opacity
(the AMES-Cond limit).  Below ~2400 K the condensate column becomes
optically significant and a Dusty-mode treatment (grain opacity from
the condensed fractions) would be needed.  Settling/microphysics
(BT-Settl/StaticWeather style) is brown-dwarf territory -- out of scope.

---

## 3. Mann-Validation Open Items (post D0-audit + resolution fixes, 2026-08-09)

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
- **H/K windows ~4-7% bright vs PHOENIX at 2700 K on THEIR structure.**
  The two leading data hypotheses are now both excluded by measurement:
  modern CIA (Abel/HITRAN, + H2-H/H-He, + the stimulated-emission double
  count) moved K by +0.4% on the continuum, and unifying TiO with SYNTHE's
  ExoMol Toto moved the 2800 K structure by <10 K.  Look elsewhere.
- **Ca II IRT cores too deep** (LTE; known).
- **~4100 K limit cycle worsened**: PM_I10113+4927 reconverged at 18% RMS
  (was ~2.9%) on the corrected network -- rerun the limit-cycle diagnosis.
- **Solar absolute (+2.1%) recheck at R=300k**; then the sample-wide
  rerun (all medians shift with the resolution fix).

