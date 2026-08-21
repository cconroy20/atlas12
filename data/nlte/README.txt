NLTE departure-coefficient grids, read by SYNTHE when NLTE_MODE = 1.

One self-contained file per element.  Flat directory, no subdirectories: the
element and the source grid's release date are in the filename, so several
elements coexist here.

  Na_MARCS_Jul-14-2023.nlte   789 MB   Na I,   levels 1-10
  Mg_MARCS_Nov-13-2024.nlte   570 MB   Mg I,   levels 1-10
  Ca_MARCS_Jun-02-2021.nlte   825 MB   Ca I+II, levels 1-10 and 68-72
  Fe_MARCS_May-07-2021.nlte    55 MB   Fe I,   21 selected levels

Fe is the cheapest element in the collection despite having the MOST levels
(607), because it carries no A(X) abundance ladder -- A(Fe) is the metallicity
by construction -- so it has ~1 record per atmosphere instead of ~21.  Its
runtime file holds only the 21 levels the eleven tagged Fe I transitions need;
widening that is a seconds-long rebuild from the master store.  Two limits:
levels 549-606 are Fe II SUPERLEVELS (whole terms merged, g = 30, 28, ...) and
cannot represent an individual Fe II line, and levels above 566 are
contaminated in the source and must never be selected.

Ca's model atom runs both ionization stages in ONE level list (Ca II 4s sits
at 49305.952 cm^-1, the Ca I ionization energy), which is why its level
selection is in two blocks: 1-10 are Ca I, 68-72 are the Ca II levels that
H&K and the infrared triplet need.  gfall zeroes each stage separately, so
the transition table in mod_mklinelist carries E_low = 0 for both Ca I 4226
and Ca II H&K while their level indices differ.

Which files are needed depends on the synthesis window: mod_nlte opens only
the elements that own a tagged transition inside it, so a Na D run does not
require Mg or Ca to be installed.  Each has its own override -- $NLTE_GRID_NA,
$NLTE_GRID_MG, $NLTE_GRID_CA, $NLTE_GRID_FE.

Each file carries its own axes, its parameter -> record index, and the records,
so there is nothing to pair up and no way to combine an index with the wrong
data.  SYNTHE reads only the interpolation corners (~79 kB per model), never
the whole file, so size here costs disk and not run time.

PROVENANCE
  Amarsi et al. (2020, A&A 642, A62), 1D MARCS, repackaged for Turbospectrum
  by Gerber et al. (2023, A&A 669, A43), from the MPDL keeper server; links in
  TSFitPy's utilities/nlte_grids_links.cfg.  That collection holds 17 elements
  (H, O, Mg, Al, Si, Ca, Ti, Cr, Mn, Fe, Co, Ni, Na, Sr, Y, Ba, Eu), ALL on the
  same MARCS base: Teff 2500-8000 K, log g -0.5..5.5, [Fe/H] -5..+1, vturb
  0/1/2/5, ~11200 atmospheres, ~60 per cent of them spherical.  Sources here:

    Na  NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin  15.9 GB zip, 436255 rec,
        atom.na_qmh, 290 levels
    Mg  NLTEgrid4TS_Mg_MARCS_Nov-13-2024.bin   4.6 GB zip, 319260 rec,
        atom.mg86d,  86 levels
    Ca  NLTEgrid4TS_Ca_MARCS_Jun-02-2021.bin   6.8 GB zip, 314482 rec,
        atom.ca105b, 105 levels
    Fe  NLTEgrid4TS_Fe_MARCS_May-07-2021.bin   1.7 GB zip,  15229 rec,
        atom.fe607a, 607 levels

  The master stores keeping EVERY level of these three (plus H, O and Fe,
  which are not wired into the code) are in ~/kurucz/nlte_grids/; see the
  README there.  Rebuilding a runtime file with a different level selection
  is a local operation against those and needs no download.

HOW IT WAS BUILT, AND HOW TO REBUILD
  1. tools/nlte_extract_grid.py   one streaming pass over the 15.9 GB zip,
                                  keeping tau_5000 + the 10 lowest levels of
                                  every record.  32.6 min, 15.90 GB
                                  transferred, a few hundred kB of memory.
  2. tools/nlte_build_index.py    the parameter -> record lookup, with the
                                  geometry policy (plane-parallel where MARCS
                                  has it, otherwise spherical at 1 Msun).
  3. tools/nlte_build_runtime.py  merges the two, drops the 27 per cent of
                                  records the index can never reach, and
                                  writes the single file here.

  Steps 1-2 produce the master archive, kept OUTSIDE the repository at
  ~/kurucz/nlte_grids/na/ (1.0 GB) along with the provenance JSON and its own
  README.  Step 3 is cheap and local, so narrowing the level set, changing the
  spherical mass, or regenerating after a policy change costs seconds -- only
  widening past the 10 extracted levels would need the download again.

WHAT IS IN IT
  b = b(level, tau_5000 ; Teff, log g, [Fe/H], v_turb, A(Na))

  Axes: Teff 32 (2500-8000 K), log g 13 (-0.5..5.5), [Fe/H] 15 (-5..+1),
  v_turb 4 (0,1,2,5 km/s), dNa = A(Na)-[Fe/H] 31 (about +-1.5 dex on solar).
  773760 cells, 41 % filled -- the empty ones are the unphysical corners of
  the HR diagram, so the grid is HR-shaped rather than randomly holey.
  tau_5000 is stored PER RECORD: every MARCS model has its own depth grid.

  Levels kept (of the atom's 290), and what they reach:
    1 3s 2S    2 3p 2P1/2   3 3p 2P3/2      Na D 5891.6 / 5897.6 = 1-3 / 1-2
    4 4s 2S    5 3d 2D      6 3d 2D         8183 / 8195          = 2-5 / 3-5
    7 4p 2P1/2 8 4p 2P3/2   9 5s 2S         1.14 um              = 2-4 / 3-4
   10 4d 2D
  Only 1, 2 and 3 are read for Na D; the rest cost 500 MB and are what makes
  the other features reachable without touching the published grid again.

TWO CONSTRAINTS THE GRID IMPOSES
  [alpha/Fe] is NOT an axis.  These are MARCS '_st_' models, in which alpha is
  a deterministic function of [Fe/H] (0 at [Fe/H] >= 0, ramping to +0.4 by
  <= -1): 15 ([Fe/H],[alpha/Fe]) pairs exist of 75 possible.  An alpha-enhanced
  model therefore gets b computed at the standard alpha for its metallicity.

  MARCS is plane-parallel only at log g >= 3 and 100 % spherical below, while
  ATLAS12 is plane-parallel throughout.  Giants use spherical models at
  1 Msun; every cell carries a geometry flag and SYNTHE warns when an
  interpolation mixes the two, which can only happen for log g strictly
  between 2.5 and 3.0.
