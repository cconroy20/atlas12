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

## 1b. Reconverge the cool-star structures

**Status: GJ644C done 2026-08-10; GJ905 and GJ699 open — and ALL of them
now need redoing again for a different reason: the molecular van der Waals
widths were corrected the same day (see CHANGELOG), which changes molecular
blanketing enough that every cool-star structure on disk is built with the
wrong line widths.  The GJ644C numbers below are syntheses on the old
structure and are not self-consistent.**  Every cool-star
ATLAS12 structure on disk predates the opacity fixes and none has CIA
feedback: the dropped-absorber fix (106d8b7), H2- free-free (dff6224),
and the CIA modernization + TiO unification.  Band-level conclusions
drawn on the old structures are not safe.  Cost is the limiting factor,
not correctness -- and it is now much lower than feared: with the
block-stream binary readers the 2700 K full 350-2500 nm synthesis takes
15 min, not the 2.5 h the pre-`be077cb` run needed, so a whole ladder
star is well under an hour.

The GJ644C (2700 K) rerun sets expectations for the rest.  The near-IR
excess vs the Mann spectrum collapsed -- H 1.328 -> 1.212, K 1.209 ->
1.078 -- while the optical did not move at all (TiO/CaH windows shift by
<0.02, and the model/data ratio curves lie on top of each other below
1 um).  The structural response is modest and confined to the optically
thin layers: at fixed column mass the new model is 30-55 K warmer for
log m < -1 and within 5-20 K through the photosphere.  So the opacity
changes are near-additive on the emergent spectrum, and reconverging the
remaining stars is bookkeeping rather than a physics unknown.

Two caveats found doing it: the reconverged model landed at 1.96% RMS /
11.7% max (old 1.61% / 6.91%) -- the usual cool-dwarf floor -- and it
carries a one-layer +87 K step in T(m) at J=25 (log m = -1.14, T ~ 1750 K)
that the old structure does not have.  It sits well outside the tau ~ 1
layers where the H/K continuum forms and the optical is unchanged, so it
does not affect the result above, but it should be understood before this
structure is used for fine T(m) work.

---

## 1c. Retired line lists staged for deletion

**Status: awaiting a decision.**  4.2 GB moved out of `data/` on
2026-08-10 to `~/kurucz/superseded_linelists/`, unreferenced by either
code: `schwenke.bin` (Schwenke 1997 TiO, replaced by ExoMol Toto),
`mol/tiototo.bin` (previous Toto), `h2ofastfix.bin` (P&S H2O, replaced by
POKAZATEL).  All were gitignored, so deleting them is unrecoverable from
the repo -- hence staged rather than removed.  Delete when satisfied.

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

---

## 6. NLTE departure coefficients (Na I D shipped, default off)

**Status: machinery DONE and committed (854228b), `NLTE_MODE = 0`.**  Four
modes; see CHANGELOG and `data/nlte/README.txt`.  Open items below, roughly
in the order they would have to be settled before a C3K grid rollout.

- **The grid's atmospheres do not reach the top of ours.**  MARCS stops
  1.1-2.0 dex of column mass short, so the outermost 15-27 layers hold an
  endpoint `b`.  |S_l/B - 1| is largest exactly at that edge and still
  growing, and the D-line cores form above it, so **every effect measured
  so far is a lower limit**.  This cannot be fixed by better interpolation
  -- it needs `b` computed on our own structures.
- **Giants use spherical `b` on plane-parallel structures.**  MARCS has no
  plane-parallel model at log g <= 2 (it is 100% spherical there) while
  ATLAS12 is plane-parallel throughout.  Sphericity changes the geometry of
  the radiation field, which is what sets `b` for a resonance line.
  Unmeasured.  Worth checking on one cool giant before committing the giant
  half of any grid; it may argue for shipping dwarfs first.
- **[alpha/Fe] cannot be varied.**  MARCS `_st_` models tie it to [Fe/H] by
  the standard relation, so an alpha-enhanced model necessarily gets `b` at
  the standard alpha for its metallicity.  Expected minor (Na's departures
  depend on alpha only through T(tau) and n_e) but unquantified.
- **The log g 2.5/3.0 geometry seam.**  Interpolating across it mixes
  spherical and plane-parallel corners; SYNTHE warns.  C3K's 0.5-dex log g
  grid never lands strictly between, so this only bites for arbitrary log g.
- **QA sweep before any rollout.**  Run mode 3 over a coarse subset spanning
  the HR diagram and map the clamped fraction, the surviving interpolation
  weight and S_l/B on the HR diagram.  That table is the real go/no-go, and
  it is cheap now that the grid is local.
- **Verify `b -> 1` at the 8000 K grid ceiling.**  Above it there is no NLTE
  data and Na I is ionised away, so LTE is right -- but if `b` has not
  converged to 1 by the boundary, switching there puts a seam into the grid.
- **Other elements.**  The extraction, index and runtime tooling is
  element-agnostic; only the transition tagging is per-element.  The same
  grid family covers ~20 elements, and Fe I would matter far more broadly
  for C3K than Na does.  Decide whether this is "NLTE Na D" or "the NLTE
  pipeline" before building more of it.
- **Confounders for any comparison with data**, both larger than the NLTE
  effect in the cool regime: Na D wings in an M dwarf are broadened by H2
  and He rather than H I, so ABO does not apply and the Unsold mix in
  `txnxn` is doing real work (see item 2); and the D cores are
  chromospherically filled, which no radiative-equilibrium model reproduces.

## 7. Second out-of-bounds read blocking `-fcheck=all`

`number` in `atlas12_modules.f90` is indexed at 2 in a dimension of extent 1
(around line 8713).  Found while trying to run a bounds-checked build during
the NLTE work; the companion violation in `interp_logu_cached` was an
`.AND.` short-circuit assumption and is fixed, but this one is not obviously
the same class and was left alone.  Until it is resolved no whole-program
`-fcheck=all` build will run, which costs a useful debugging tool.
