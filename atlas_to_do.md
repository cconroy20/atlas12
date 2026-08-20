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

## 3. Retire ATLAS12's private atomic line files

**Status: open.**

ATLAS12's `SELECTLINES` still reads Kurucz-era binaries that SYNTHE does not
use, while SYNTHE resolves every list through the `lines.list` manifest and
`mod_mklinelist`:

| source | ATLAS12 (`SELECTLINES`) | SYNTHE |
|--------|-------------------------|--------|
| predicted atomic | `gfpred29dec2014.bin` (unit 11) | same file, via `lines.list` |
| observed atomic  | `lowobsat12.bin` (unit 111) + `hilines.bin` (unit 21) | `gfallvac08oct17.dat` |
| complex profiles | `nltelinobsat12.bin` (unit 19) | derived from `gfall` by `read_gfall`'s TYPE dispatch |
| diatomics, TiO, H2O, polyatomics | already via `lines.list` | via `lines.list` |

The molecular half was unified in 2026-08 (ATLAS12 calls
`read_diatomics_for_atlas` in `mod_mklinelist`), so the pattern and the
machinery both exist; the atomic half was left behind.

**Why it matters, and it is not hypothetical.**  `lowobsat12.bin` +
`hilines.bin` are a pre-split binary packaging of nominally the same Kurucz
atomic data that `gfallvac08oct17.dat` carries, so the two codes can silently
disagree about what lines exist and with what parameters -- structure built on
one list, spectrum synthesised from another.  That exact failure has already
happened once here with TiO (ATLAS12 on Schwenke 1997 while SYNTHE had moved to
ExoMol Toto, the case hard-wired past the manifest), and it went unnoticed
until someone went looking.  Nothing prevents a repeat on the atomic side.

**What to do.**  Extend `mod_mklinelist` with an ATLAS12-facing atomic reader
(the counterpart of `read_diatomics_for_atlas`), have `SELECTLINES` take its
atomic lines from `lines.list`, and delete the three private binaries.  Any
difference in the resulting structure is itself the interesting result: it is a
direct measure of how far the two lists had drifted.

**Watch for:** `SELECTLINES` applies its own `XNFDOPMAX`/`TABCONT` rejection
before storing a line, which is not the same filter as SYNTHE's `LINE_CUTOFF`
retention -- the goal is a shared line *source*, not identical selection.  The
predicted list is also the one where `gfpred26apr2018.bin` was measured as not
worth adopting (see CHANGELOG), so leave that choice alone while moving it.
