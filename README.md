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
Stehlé–Hutcheon) lives in `synthe_module.f90`.

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
| `h2collop.dat`       | `H2COLLOP`   | H₂ collision-induced absorption (Borysow et al. 1997) |
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
| `h2ofastfix.bin` †       | SELECTLINES      | H₂O line list (Partridge & Schwenke) |
| `schwenke.bin` †         | SELECTLINES      | Schwenke diatomic/metal-hydride line list |
| `hilines.bin` †          | SELECTLINES      | Hydrogen/helium line table |
| `lowobsat12.bin` †       | SELECTLINES      | Low-excitation observed atomic lines |
| `nltelinobsat12.bin` †   | ATLAS12 / XLINOP | NLTE line data |
| `mol.tar.gz` †           | —                | Archive of molecular sub-lists referenced from `lines.list`; unpack in place |

† Not tracked in the repository; download from the Google Drive folder above.

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

Discipline: after **any** edit to `molecules.dat` or `condensates.dat`,
rerun the corresponding `--validate` and regenerate the fit atlas.

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
- Bergemann, M., Lodders, K., & Palme, H. 2025, Zenodo record 14988840 (solar abundance scale, `solar=berg25`)
- Borysow, A., Jørgensen, U. G., & Zheng, C. 1997, A&A, 324, 185 (H₂ collision-induced absorption)
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
