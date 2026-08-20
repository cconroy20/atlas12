# ATLAS12 + SYNTHE

Fortran 90 translation of Robert L. Kurucz's stellar atmosphere and spectral
synthesis codes, originally developed at the Harvard-Smithsonian Center for
Astrophysics.

This translation was written by Charlie Conroy and Claude.AI.

Companion documents: [CHANGELOG.md](CHANGELOG.md) — the full modernization
ledger, with design rationale and validation numbers for every change —
and [atlas_to_do.md](atlas_to_do.md) — open items.

## Programs

### ATLAS12 — Opacity-Sampling Stellar Atmosphere Code

Iteratively solves for the structure of a plane-parallel, LTE stellar
atmosphere in radiative and hydrostatic equilibrium.  Opacity sampling
over ~30,000 wavelength points replaces the older opacity-distribution
function approach of ATLAS9.

At each iteration the code solves hydrostatic equilibrium, computes
Saha-Boltzmann ionization and molecular equilibrium for all species
(including equilibrium condensation at cool photospheric temperatures),
builds continuous and line opacity tables, loops over the full
wavelength grid solving the Feautrier radiative transfer equation at
each point, and applies a temperature correction to enforce radiative
plus convective flux conservation.

### SYNTHE — Spectral Synthesis

Computes the emergent spectrum from a converged ATLAS12 model atmosphere.
Integrates three formerly separate programs (XNFPELSYN, SYNTHE, SPECTRV)
into a single executable with all data flow in memory:

- **XNFPELSYN**: computes continuum opacities and ion populations
- **SYNTHE**: builds the line opacity vector at each wavelength point
- **SPECTRV**: solves radiative transfer for the emergent surface flux

The output is an ASCII spectrum file with wavelength (Angstroms), flux,
and continuum flux at each point.

## Quick Start

```
git clone https://github.com/cconroy20/atlas12.git
cd atlas12

# Download the eight large data files from the Google Drive folder
# linked in the Input Data section, and place them in data/.
# Unpack the molecular archive in place:
cd data && tar xzf mol.tar.gz && cd ..

# Build both executables
cd src && make && cd ..

# Point the code at the data directory
export ATLAS12=$(pwd)

# Move into the work directory, which ships with a solar starting atmosphere
cd workdir

# Run ATLAS12 on the solar model (first arg = input atmosphere, second = output basename)
../bin/atlas12c.exe ap00t5777g4.44at12.dat sun

# Synthesize a visible spectrum from the converged model
../bin/synthe.exe sun.atm wlbeg=400 wlend=700
```

See [Running ATLAS12](#running-atlas12) and [Running SYNTHE](#running-synthe)
for the full argument syntax, and [Input Data](#input-data) for the data
files required.

## Source Files

| File | Description |
|------|-------------|
| `atlas12_modules.f90` | Shared modules, all ATLAS subroutines (EOS, opacity, transfer, READIN, JOSH, ...) |
| `atlas12c.f90`        | ATLAS12 main program (iteration driver) |
| `synthe_module.f90`   | SYNTHE shared data and procedures (hydrogen/He profiles, line opacity, `run_xnfpelsyn`) |
| `mod_mklinelist.f90`  | In-memory line-list preprocessor (replaces standalone mklinelist + `synbeg`/`rgfall`/`rpredict`/`rmolecasc`/`rh2ofast` pipeline) |
| `mod_nlte.f90`        | NLTE departure coefficients for named transitions (SYNTHE only; `NLTE_MODE`, default off) |
| `synthe.f90`          | SYNTHE main program (spectral synthesis driver) |

## Building

Requires `gfortran`.  From the source directory:

```
make              # build both atlas12c.exe and synthe.exe
make atlas        # build atlas12c.exe only
make synthe       # build synthe.exe only
```

Executables are installed to `../bin/`.  The source directory is left
clean (no `.o` or `.mod` files).

## Running ATLAS12

```
export ATLAS12=/path/to/atlas12/
atlas12c.exe <input_atm> [basename] [key=value ...]
```

The first positional argument is the input atmosphere file and is
required.  This file contains the standard Kurucz keyword cards (TEFF,
GRAVITY, ABUNDANCE, TURBULENCE, CONVECTION, READ DECK, ...) terminated
by a `BEGIN` card.  The optional second positional `basename` argument
sets the prefix for all output files (default `mystar`); keyword
arguments may appear in any order, before or after the positionals.
Pass `--help` (or `-h`, `help`) to print usage and exit.

Output files:

| File | Contents |
|------|----------|
| `<basename>.atm`   | Converged model atmosphere. The `READ DECK6` block has one line per depth with columns: `RHOX`, `T`, `P`, `XNE`, `ABROSS` (Rosseland mean opacity κ), `ACCRAD`, `VTURB`, `FLXCNV`, `VCONV`, `RHO` (mass density, g/cm³), `TAU5000` (continuum optical depth at 5000 Å — absorption + continuum scattering, no lines — the reference depth scale tabulated by MARCS/PHOENIX) |
| `<basename>.flux`  | Emergent flux vs. wavelength |
| `<basename>.iter`  | Per-iteration summary, including the temperature-correction diagnostics.  When the deep-CZ polish runs, its post-polish per-layer state is appended as a final block |

ATLAS12 writes exactly these three files.  Earlier versions also emitted
separate `.taunu` and `.tcorr` files; `.taunu` is no longer written, and the
`.tcorr` temperature-correction diagnostics are now columns in `.iter`.

Only the first seven DECK6 columns (`RHOX`…`VTURB`) are read back when a
model is used as input; `FLXCNV`, `VCONV`, `RHO`, and `TAU5000` are
write-only diagnostics for downstream use.

Command-line options (keyword=value):

| Option       | Default     | Description |
|--------------|-------------|-------------|
| `numit=N`    | 30          | Number of iterations |
| `vturb=X`    | from model  | Microturbulence (km/s) |
| `mlt=X`      | from model  | Mixing length parameter |
| `teff=X`     | from model  | Rescale model to this Teff (K) |
| `logg=X`     | from model  | Rescale model to this log g (cgs) |
| `solar=name` | from model  | Solar reference abundance pattern (`ag89`, `agss09`, `berg25`); replaces the abundance table carried by the model file |
| `zscale=X`   | 1.0         | Metal abundance scale factor (multiplicative on Z≥3) |
| `heabnd=X`   | from model  | He number fraction Y; H is recomputed as X = 1 − Y − Z |
| `abund=file` | none        | Individual element overrides (see below) |

Numerics and physics switches are deliberately not CLI options.  They
live as developer flags at their declarations in `mod_atlas_data` (set
the default value and recompile): `USE_CONDENSATION` (equilibrium
condensation, default on), `USE_TOPBASE_MBF` (TOPbase metal continuum
vs. legacy analytic fits, default on), `USE_CZ_CONSTRUCTOR` (deep-CZ
temperature constructor, default on), `USE_CZC_POLISH` (terminal
deep-CZ flux-closure polish, default on), `USE_FLXCNV_SMOOTH` (interior
1-2-1 convective-flux smoothing, default on), `TFLOOR_ATM` (temperature-
correction floor, 1200 K), `IROSSTAB` (Rosseland-table interpolation:
1=bilinear, 2=Shepard, 3=moving least squares), and `IQUAD` (`INTEG`
quadrature: 0=legacy blended-parabola, 1=Steffen monotone cubic).
`USE_KP_HYDROGEN` (Kurucz–Peterson hydrogen Stark profiles instead of
Stehlé–Hutcheon) lives in `synthe_module.f90`, and `CA4227_MODE`
(Ca I 4227 resonance profile) and `NLTE_MODE` (departure coefficients
for selected lines, see "Departure coefficients" below) live in
`mod_parameters`.

Abundance override file format: one element per line with two
whitespace-separated columns, `Z  log10(number_fraction)`.  Lines
starting with `#` or `!` are treated as comments.  Example:

```
# carbon and iron abundances
 6  -3.52
26  -4.54
```

The abundance overrides are applied in order after the model file is
read: `solar=` replaces the reference pattern, `zscale=` shifts all
metals, `abund=` overrides individual elements, and `heabnd=` sets Y.
`solar=` leaves the relative offsets (`XRELATIVE`) untouched, so a
model's [M/H] scaling is preserved across a change of solar zero-point.
When any of these is used, ATLAS12 renormalizes so that X + Y + Z = 1
(aborting if the specified Y + Z would drive X negative), then
recomputes all abundance-dependent quantities before the iteration
loop.  If both `teff=` and `logg=` are given, the model is regridded
via `SCALE_MODEL` before iteration begins.

## Running SYNTHE

```
export ATLAS12=/path/to/atlas12/
synthe.exe <model_file> wlbeg=<nm> wlend=<nm> [resolu=<R>] [turbv=<kms>] [more_output=<yes|no>]
```

The merged executable performs line-list construction, continuum opacity
computation, line opacity accumulation, and radiative transfer in a
single run — there is no longer a separate SYNBEG / line-reader /
SYNTHE pipeline, and no intermediate `fort.*` files are written or
read.  Line lists are built in memory by `run_mklinelist`, which reads
`lines.list` from `$ATLAS12/data/` and dispatches internally to the
appropriate readers (gfall, predict, mol, h2o).

### NLTE departure coefficients (default off)

`NLTE_MODE` in `mod_parameters` is **0 = off** (pure LTE; nothing allocated,
no NLTE code runs) or **1 = on**.  When on, `b_l` and `b_u` for named
transitions rescale the line opacity and replace `SLINE`:

```
kappa = kappa_LTE * b_l * [1 - (b_u/b_l) e^-x] / [1 - e^-x]
S_l   = (2h nu^3/c^2) / [ (b_l/b_u) e^x - 1 ]        x = h nu / kT
```

These multiply to `kappa*S = b_u * kappa_LTE * B_nu * (1-e^-x)`, so only one
extra accumulator is needed, carrying the opacity-weighted *deviation* of the
emissivity from LTE; it stays identically zero for every LTE line.  Eligible
transitions are declared in `mod_mklinelist` and matched on species plus both
level energies, so every hyperfine component is tagged at once — currently the
Na I D doublet.

SYNTHE reads one self-contained file from `data/nlte/` (`$NLTE_GRID`
overrides) and interpolates `b` over Teff, log g, [Fe/H], v_turb and Na
abundance, reading only the interpolation corners.  Each corner carries its own
τ₅₀₀₀ grid, so corners are placed on this model's τ₅₀₀₀ before being combined;
missing corners have their weight dropped and the remainder renormalised, with
the surviving fraction reported.  Two constraints the grid imposes: `[α/Fe]` is
not an axis (MARCS `_st_` models tie it to `[Fe/H]`), and MARCS is
plane-parallel only at log g ≥ 3, so giants use spherical models at 1 M⊙.

The file is derived — `tools/nlte_extract_grid.py`, `nlte_build_index.py`,
`nlte_build_runtime.py` — from a published grid whose master copy lives outside
the repository.  See `data/nlte/README.txt` for provenance and the rebuild, and
CHANGELOG for the physics and its caveats.

The electron density (with its consistent `XNATOM` and `RHO`) is
recomputed self-consistently from the model structure rather than taken
from the `.atm` file (compile-time toggle `RECOMPUTE_XNE` in
`run_xnfpelsyn`); for atmospheres converged with the current ATLAS12 the
stored and recomputed values agree to solver tolerance.

Arguments:

| Argument              | Required | Description |
|-----------------------|:--------:|-------------|
| `<model_file>`        | yes      | ATLAS12 `.atm` model atmosphere (positional, 1st) |
| `wlbeg=<nm>`          | yes      | Start wavelength in nanometers |
| `wlend=<nm>`          | yes      | End wavelength in nanometers (> wlbeg) |
| `resolu=<R>`          | no       | Resolving power λ/Δλ (default 300 000) |
| `turbv=<kms>`         | no       | Microturbulence in km/s. If > 0, **replaces** the model atmosphere's microturbulence at all depths; if omitted or ≤ 0 (default 0.0), the per-layer value from the input model is used |
| `more_output=<yes\|no>` | no     | If yes, also write `.linform` and `.mol` diagnostic files (default no).  Accepted truthy values: `yes`, `true`, `1`, `on`, `y` (case-insensitive); falsy: `no`, `false`, `0`, `off`, `n`. |

The output basename is derived from the model filename by stripping the
leading directory and trailing extension.  For example,
`synthe.exe models/sun.atm wlbeg=400 wlend=700` produces `sun.spec` in
the current directory.  Adding `more_output=yes` additionally produces
`sun.linform`, `sun.mol`, and `sun.lines`.

Output files:

| File              | Written                 | Contents |
|-------------------|-------------------------|----------|
| `<base>.spec`     | always                  | ASCII spectrum: wavelength (Å, F11.4), flux (E15.6), continuum flux (E15.6) |
| `<base>.linform`  | only if `more_output=yes` | Per-wavelength diagnostic: wavelength, emergent H, surface H, monochromatic optical depth at each atmospheric layer |
| `<base>.mol`      | only if `more_output=yes` | Molecular number-density diagnostics vs. depth for all species tracked by the equation of state |
| `<base>.nlte`       | only if `NLTE_MODE` /= 0   | Per-(transition, depth) departure coefficients and the opacity/source-function factors derived from them; see "Departure coefficients" above |
| `<base>.lines`      | only if `more_output=yes` | Line list for all lines used in synthesis. Columns are: LTE/NLTE, vacuum wavelengths (Å), species code, nelion (internal species identifier), ELO, cgf (strength indicator), gamma_rad, gamma_stark, gamma_vdW |

Wavelengths are handled internally in nanometers on a logarithmic grid
with spacing `ratio = 1 + 1/resolu`; vacuum wavelengths are used
throughout.  The `.spec` file reports wavelengths in Angstroms for
compatibility with legacy post-processing.

## Input Data

Set the `$ATLAS12` environment variable to the installation root.  All
data files are read from `$ATLAS12/data/` (falling back to `./data/`
if unset).

Most files in the data directory are tracked in the GitHub repository.
Eight large files are **not** in the repository and must be downloaded
separately from

> https://drive.google.com/drive/u/0/folders/1vzl0j_aUIpOQpz480vwhUCsmWKNR2WB9

and placed in `$ATLAS12/data/` before running ATLAS12 or SYNTHE.  These
are marked with † in the tables below.  After downloading, unpack
`mol.tar.gz` in place — `lines.list` references the individual molecular
sub-lists that the archive expands to.

The full contents of the data directory, organized by purpose:

**Equation of state and partition functions**

| File | Used by | Contents |
|------|---------|----------|
| `ionpots.dat`       | `IONPOTS`  | Ionization potentials, all species |
| `isotopes.dat`      | `ISOTOPES` | Isotope mass fractions |
| `molecules.dat`     | `READMOL`  | Molecular equilibrium constants, 297 rows (see [Tools](#tools)) |
| `condensates.dat`   | `READCOND` | Condensate saturation ln K(T) fits, 21 solids, for equilibrium condensation (see [Tools](#tools)) |
| `partfn_bc2016.dat` | B&C partition-function module (lazy-loaded) | Barklem & Collet (2016) atomic partition functions, Z = 1–92, ion stages I–III — the production U(T) source |
| `pfsaha.dat`        | `PFSAHA`   | Legacy Kurucz atomic partition-function data (retained as an internal safety net) |
| `pfiron.dat`        | `PFIRON`   | Iron-group partition functions (pressure-lowering correction) |
| `partfnh2.dat`      | `PARTFNH2` | H₂ partition function vs. temperature |

**Continuous opacities**

| File | Used by | Contents |
|------|---------|----------|
| `crossch.dat`        | `CHOP`       | CH bound-free + bound-bound cross-section table |
| `crossoh.dat`        | `OHOP`       | OH bound-free + bound-bound cross-section table |
| `h2collop.dat`       | `H2COLLOP`   | Collision-induced absorption, H₂–H₂/H₂–He/H₂–H/H–He (Abel et al. 2011, 2012; Gustafsson & Frommhold 2001, 2003) |
| `hotop.dat`          | `HOTOP`, `MBF_HIGH_ION` | Hot-star opacities (high-ionization species) |
| `mbf/` (33 files)    | `MBF_TOPBASE` | TOPbase resonance-averaged photoionization grids, 30 species (Allende Prieto et al. 2003) |
| `op_fe1.dat`, `op_fe2.dat` | `FELO_OPACITY` | Fe I / Fe II R-matrix bound-free cross sections (Bautista 1997; Nahar & Pradhan 1994) |
| `gauntff_vanhoof.dat` | `read_gauntff_table` | Free-free Gaunt factors (van Hoof et al. 2014) |
| `karzas_ekarzas.dat` | `read_karzas_tables` | Karzas–Latter tabulated Gaunt-factor energy grid |
| `karzas_freqn.dat`   | `read_karzas_tables` | Karzas–Latter frequency grid |
| `karzas_xl.dat`      | `read_karzas_tables` | Karzas–Latter ℓ-resolved cross sections |
| `karzas_xn.dat`      | `read_karzas_tables` | Karzas–Latter total-n cross sections |

**Radiative transfer and line profiles**

| File | Used by | Contents |
|------|---------|----------|
| `blockj.dat`          | `BLOCKJ`            | Feautrier J-operator coefficient matrices |
| `blockh.dat`          | `BLOCKH`            | Feautrier H-operator coefficient matrices |
| `stark_profile.dat`   | `XLINOP`            | Legacy Kurucz–Peterson Stark broadening profiles |
| `stehle_lyman.bin`    | `INIT_STARK_TABLES` | Stehlé–Hutcheon (1999) Stark tables, Lyman series |
| `stehle_balmer.bin`   | `INIT_STARK_TABLES` | Stehlé–Hutcheon (1999) Stark tables, Balmer series |
| `stehle_paschen.bin`  | `INIT_STARK_TABLES` | Stehlé–Hutcheon (1999) Stark tables, Paschen series |
| `stehle_brackett.bin` | `INIT_STARK_TABLES` | Stehlé–Hutcheon (1999) Stark tables, Brackett series |
| `he1tables.dat`       | `read_he1_stark_tables` | He I Stark broadening tables (port incomplete — currently unused) |

**Line lists**

| File | Used by | Contents |
|------|---------|----------|
| `lines.list`             | `run_mklinelist` | Plain-text manifest pointing at the line-list files below |
| `gfallvac08oct17.dat` †  | `read_gfall`     | Kurucz atomic line list (vacuum wavelengths, Oct 2017) |
| `gfpred29dec2014.bin` †  | `read_predict`, SELECTLINES | Kurucz predicted atomic lines (Dec 2014) |
| `hilines.bin` †          | SELECTLINES      | Hydrogen/helium line table |
| `lowobsat12.bin` †       | SELECTLINES      | Low-excitation observed atomic lines |
| `nltelinobsat12.bin` †   | ATLAS12 / XLINOP | NLTE line data |
| `mol.tar.gz` †           | —                | Archive of molecular sub-lists referenced from `lines.list`; unpack in place |
| `h2opokazatel.bin` ‡     | ATLAS12 / SYNTHE | H₂O pseudo-line list (51.3M records) built from ExoMol POKAZATEL; replaces `h2ofastfix.bin` (P&S 1997).  Rebuild with `tools/build_h2o_pokazatel.py --raw --write-raw` |
| `mol/tiototo2024.bin` ‡  | ATLAS12 / SYNTHE | TiO line list (131.6M records, ⁴⁶Ti–⁵⁰Ti) from ExoMol Toto; replaces `schwenke.bin` (Schwenke 1997).  Both codes resolve it through `lines.list`, so they cannot diverge |
| `nlte/*.nlte` ‡          | SYNTHE           | NLTE departure-coefficient grids, one self-contained file per element, for `NLTE_MODE = 3`.  Ships with Na I (789 MB).  Derived — see `data/nlte/README.txt` for provenance and the three-stage rebuild |
| `mol/alo_atp.dat` ‡      | ATLAS12 / SYNTHE | AlO line list (4.93M records) from ExoMol ATP.  The B–X bands at 4842 and 4648 Å reach 60% of the local extinction at the τ(4500 Å) = 1 layer of a 2900 K dwarf; on by default.  Rebuild with `tools/exomol_to_kurucz.py --gns 6 --icode 813 --iso 16` |

† Not tracked in the repository; download from the Google Drive folder above.
‡ Not tracked; regenerated from ExoMol by the named tool.

## Tools

`tools/` holds the Python pipeline that generates and validates the
fitted data files.  None of it is needed to build or run the code; it
keeps the provenance of every fitted row in the repository.  (The raw
reference tables it consumes — BC16, ExoMol, the GGchem clone — are
not in the repository.)

| Tool | Purpose |
|------|---------|
| `fit_molecule_keq.py` | Fit molecular equilibrium constants for `molecules.dat`; `--validate` regenerates all 96 diatomic rows (including the 8 molecular ions) from source and must stay exact |
| `fit_condensates.py`  | Build `data/condensates.dat` from the GGchem condensate compilation; `--validate` round-trips all 21 solids (<0.02 dex) |
| `validate_condensation.py` | Point-matched GGchem eqcond comparison at each condensing layer of an ATLAS12 run (consumes the `COND:`/`CONDEPS:` run-log diagnostics) |
| `comp_pf.py` (+ `bc16_loader.py`, `exomol_loader.py`, `matcher.py`, `atomic_saha.py`, `polyatomic_d0.py`, `polyatomic_assembly.py`, `kurucz_molec.py`) | Regenerate the per-species fit atlas `comp_pf.pdf` — filed fits against BC16 / ExoMol / JANAF references |
| `dustchem_loader.py`, `ggchem_loader.py`, `janaf_loader.py` | Parsers/evaluators for the GGchem DustChem, dispol, and raw NIST-JANAF data in their native conventions |
| `build_molecules_physical_dat.py` | Faithful reproduction of the April 2026 refit (regenerates the pre-correction molecular-ion rows by design — historical record, not production) |
| `mann_lib.py` | Shared library for the Mann external-validation family: star resolution by GJ name, both Mann spectral libraries (registration-corrected to vacuum/rest), measured LSF, band windows, curated K-dwarf params, the point-comparison star ladder, and interpolated PHOENIX NewEra LowRes access |
| `validate_mann.py` | Run the full ATLAS12+SYNTHE chain for one Mann star and compare absolutely (no renormalization); writes metrics JSON + the standing three-panel figure |
| `mann_pointcomp.py` | Batch driver for the ~10-star point-comparison ladder; writes `workdir/mann/summary.md` |
| `mann_compare_plot.py` | The house three-panel comparison figure (data / ATLAS12+SYNTHE / PHOENIX NewEra), importable as `plot_star()` |
| `mann_lsf_fit.py` | Per-star effective-resolution measurement of the Mann spectra (chunk fits with velocity nuisance) |
| `uves_compare.py` | High-resolution line-profile comparison of a ladder star's R=300k synthesis against ESO UVES spectra |
| `nlte_extract_grid.py` | One streaming pass over a published NLTE grid (15.9 GB zipped for Na I), keeping τ₅₀₀₀ and the low-lying levels of every record.  Pays the download once instead of per model; writes the master store, which lives outside the repository |
| `nlte_build_index.py` | Parameter→record lookup over the master store (Teff, log g, [Fe/H], v_turb, ΔNa), including the geometry policy: plane-parallel where MARCS has it, otherwise spherical at 1 M⊙ |
| `nlte_build_runtime.py` | Merges store and index into the single self-contained `data/nlte/*.nlte` file SYNTHE reads, dropping the records the index cannot reach.  Seconds to run, so policy changes are a local rebuild |
| `nlte_check_dump.py` | Checks a `<model>.nlte` diagnostic dump for internal consistency (returned factors against their equations, the `κ·S = b_u·κ_LTE·B_ν` identity) — separates a plumbing bug from bad grid data |
| `nad_nlte_plot.py` | Stacked LTE-vs-NLTE panels for the Na D region, with per-panel axis control and optional instrumental smoothing (`--smooth-to`) |
| `tmin_perturb.py`, `tmin_fit.py`, `tmin_rf.py` (+ plot drivers) | T(τ) perturbation / T-min fitting / response-function machinery on converged models |
| `build_h2o_pokazatel.py` | Build `data/h2opokazatel.bin` from ExoMol POKAZATEL: exact raw-transition binning (`--raw`/`--validate-raw`/`--write-raw`) plus the super-line NNLS cross-check route |
| `build_cia_table.py` | Build `data/h2collop.dat` from the HITRAN CIA sets (+ BJF01 for the H₂–H₂ continuation above Abel's range); `--validate` re-derives the published Abel/Borysow comparison, the ν² low-frequency slope, and seam continuity from the raw files.  Source URLs in the docstring |
| `cia_ab.py` | Same-structure A/B of CIA variants at 2700 K against the PHOENIX NewEra spectrum: per-window flux and continuum ratios |
| `tio_ttau_plot.py` | T(τ) comparison of two ATLAS12 runs (used for the TiO line-list swap); reads a partial `.iter` file to preview a run still converging |
| `build_mol_broad.py` | Build `data/mol_broad.dat` from ExoMol `.broad` files — per-species molecular van der Waals widths, replacing the two hardcoded constants that carried no molecular data.  Handles the `a0`/`a1`/`m0` parameterisations and records, per species, whether the numbers are its own or a named chemical analogue.  Source URLs in the docstring |
| `cont_audit.py` | Decompose the continuum-opacity budget from a `<model>.cont` dump (`DUMP_CONTINUUM` in `mod_atlas_data`): integrates τ_cont and reports which absorber carries the opacity at the τ = 1 layer, per window |
| `mann_phoenix_ab.py` | Band-by-band comparison of one or more self-consistent models against both the Mann spectrum and PHOENIX NewEra (the counterpart to `cia_ab.py`, which is same-structure) |
| `atm_compare.py` | Overlay the T(depth) structure of several `.atm` files.  Differences are taken on COLUMN MASS: the τ columns of different codes are not the same quantity, and comparing at equal τ manufactures hundreds of K |
| `broadening_plot.py` | Spare data / PHOENIX / one-model figure for a Mann star, near-IR panel plus a full-range ratio panel |

Discipline: after **any** edit to `molecules.dat` or `condensates.dat`,
rerun the corresponding `--validate` and regenerate the fit atlas.

`data/mol_broad.dat` (per-species molecular van der Waals widths) is
small enough to live in the repository; the ExoMol `.broad` files it is
built from are not, and `build_mol_broad.py --fetch` re-downloads them.

## Translation from Fortran 77

The original ATLAS12 was ~23,000 lines of Fortran 77 fixed-format code.
The original SYNTHE suite comprised three separate programs (XNFPELSYN,
SYNTHE, SPECTRV) communicating through intermediate files.

The port modernizes both the engineering and the physics — free-format
Fortran 90 with modules in place of COMMON blocks and EQUIVALENCE
statements, full GOTO elimination, and the SYNTHE pipeline merged into
one executable with all data flow in memory; Barklem & Collet (2016)
partition functions, TOPbase/Iron Project metal continua, McLaughlin
(2017) H⁻, Stehlé & Hutcheon (1999) hydrogen Stark profiles, a 297-row
molecular equilibrium network, equilibrium condensation, and rebuilt
deep convection-zone convergence machinery.  The complete change
ledger — every change, with design rationale and validation numbers —
is in [CHANGELOG.md](CHANGELOG.md).

## References

- Kurucz, R. L. 1970, SAO Special Report 309
- Kurucz, R. L. 1993, ATLAS9 Stellar Atmosphere Programs, CD-ROM No. 13
- Kurucz, R. L. 2005, Memorie della Società Astronomica Italiana Supplementi, 8, 14
- Sbordone, L., Bonifacio, P., Castelli, F., & Kurucz, R. L. 2004, MSAIS, 5, 93
- Castelli, F., & Kurucz, R. L. 2004, astro-ph/0405087 (new grids of ATLAS9 model atmospheres)
- Castelli, F. 2005, MSAIS, 8, 25 (ATLAS12: how to use it)
- Allende Prieto, C., Hubeny, I., & Lambert, D. L. 2003, ApJS, 147, 363 (TOPbase metal photoionization grids)
- Anders, E., & Grevesse, N. 1989, Geochim. Cosmochim. Acta, 53, 197 (default solar abundance scale, `solar=ag89`)
- Asplund, M., Grevesse, N., Sauval, A. J., & Scott, P. 2009, ARA&A, 47, 481 (solar abundance scale, `solar=agss09`)
- Barklem, P. S., & Collet, R. 2016, A&A, 588, A96 (atomic and molecular partition functions)
- Bautista, M. A. 1997, A&AS, 122, 167 (Fe I R-matrix bound-free)
- Bell, K. L., & Berrington, K. A. 1987, J. Phys. B, 20, 801 (H⁻ free-free)
- John, T. L. 1975, MNRAS, 172, 305 (H₂⁻ free-free; table as filed for MARCS)
- Bergemann, M., Lodders, K., & Palme, H. 2025, Zenodo record 14988840 (solar abundance scale, `solar=berg25`)
- Abel, M., Frommhold, L., Li, X., & Hunt, K. L. C. 2011, J. Phys. Chem. A, 115, 6805 (H₂–H₂ collision-induced absorption)
- Abel, M., Frommhold, L., Li, X., & Hunt, K. L. C. 2012, J. Phys. Chem. A, 116, 3068 (H₂–He collision-induced absorption)
- Borysow, A., Jørgensen, U. G., & Fu, Y. 2001, JQSRT, 68, 235 (H₂–H₂ CIA above the Abel et al. temperature range)
- Gustafsson, M., & Frommhold, L. 2001, ApJ, 546, 1168 (H–He collision-induced absorption)
- Gustafsson, M., & Frommhold, L. 2003, A&A, 400, 1161 (H₂–H collision-induced absorption)
- Humlíček, J. 1982, JQSRT, 27, 437 (Voigt function asymptotic)
- Karzas, W. J., & Latter, R. 1961, ApJS, 6, 167 (hydrogenic bound-free Gaunt factors)
- Koester, D. 1980, A&AS, 39, 401 (convective flux at the grid bottom)
- McLaughlin, B. M., Stancil, P. C., Sadeghpour, H. R., & Forrey, R. C. 2017, J. Phys. B, 50, 114001 (H⁻ bound-free)
- Nahar, S. N., & Pradhan, A. K. 1994, J. Phys. B, 27, 429 (Fe II R-matrix bound-free)
- Sharp, C. M., & Huebner, W. F. 1990, ApJS, 72, 417 (condensate thermochemistry fits)
- Steffen, M. 1990, A&A, 239, 443 (monotone cubic quadrature, `IQUAD=1`)
- Stehlé, C., & Hutcheon, R. 1999, A&AS, 140, 93 (hydrogen Stark profiles)
- van Hoof, P. A. M., Williams, R. J. R., Volk, K., et al. 2014, MNRAS, 444, 420 (free-free Gaunt factors)
- Weideman, J. A. C. 1994, SIAM J. Numer. Anal., 31, 1497 (Voigt function rational approximation)
- Woitke, P., Helling, Ch., Hunter, G. H., et al. 2018, A&A, 614, A1 (GGchem; condensate thermochemistry compilation)
