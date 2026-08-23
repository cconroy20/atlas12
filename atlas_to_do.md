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

## 3. Ion stages absent from molecules.dat

**Status: the silent-zero trap is FIXED (2026-08-23); one judgement call
is recorded below so it is not reopened.**

`MOLEC`'s atomic lookup has three outcomes, and the middle one was the
problem: an exact code match uses the NMOLEC population, an element absent
from `molecules.dat` falls through to `PFSAHA` (correct, every stage), and
an element that IS in the table but lacks the requested stage got a silent
`NUMBER(:, ION) = 0`.  `SELECTLINES` then discarded every line of that
species on its zero-population test -- correct behaviour for a species that
genuinely is not there, and therefore indistinguishable from it.  Partial
coverage was worse than none, and the failure was silent at the lookup,
silent at the discard, and reported as `0 lines from ...`, which reads as a
physical result.  That is how 10.3M lines of Ca-Ni VI-IX stayed invisible
from Kurucz's original through to 2026-08.

`MOLZERO_RECORD` now notes every species/stage that gets zeroed and
`MOLZERO_REPORT`, called once at the end of `COMPUTE_ALL_POPS`, names them
-- but only those whose population would not have been negligible, judged
by a Saha estimate at the model's hottest layer against a 0.1% threshold,
and in one line rather than a per-species block.  Cool models zero four
stages per iron-group element on every run and stay silent, as they should;
a 10000 K model prints a single line naming the count and the worst offender
(12 stages, worst B IV at 223x the stage below), with the full list under
IDEBUG.

**Known residual, deliberately left alone.**  In the Teff = 8000-10000 K
band the deepest layers (41,000-52,000 K) sit at a ~7% flux error that does
not iterate away -- the warning above now fires there, which is the honest
outcome.  It is pre-existing (the code before the 2026-08 line list work
gave 7.71% at 10,000 K where it now gives 6.52%) and sub-photospheric, with
the photosphere itself converging to 0.06-0.19%.  Molecules-off cures it
(7.05% -> 0.134% at 9000 K) but costs a real photospheric error there
(rms 38.5 K at 9000 K, 692 K at 8000 K), so `TEFF_MOLEC_LIMIT` stays at
10000 K.  If this is ever worth fixing properly the target is the
temperature correction in the diffusion regime, not the line lists.

**Not to do: extending `molecules.dat` to stages VI-X.**  Settled by
measurement, recorded here so it is not reopened.  It would only be worth
it if some regime needed molecules and stage-VI ions at once, and none
does.  Below the gate the hottest layer any model reaches is 52,156 K (a
Teff = 10000 K, log g 4.5 model), where Fe VI is 5.5% of iron -- but that
is the deepest layer, at log tau ~ 3, and at its tau = 1 layer the Fe VI
fraction is 2e-31.  The direct test is stronger: at 25,000 K, where stage
VI+ genuinely dominates the deep layers (131,000 K), switching those lines
on moved the emergent flux by +0.002% and the photosphere by 0.16 K rms.
Adding stages would also widen the window in which the zeroing branch can
fire, rather than closing it.
