# ATLAS12 + SYNTHE

Fortran 90 translation of Robert L. Kurucz's stellar atmosphere and spectral
synthesis codes, originally developed at the Harvard-Smithsonian Center for
Astrophysics.

This translation was written by Charlie Conroy and Claude.AI.

## Programs

### ATLAS12 — Opacity-Sampling Stellar Atmosphere Code

Iteratively solves for the structure of a plane-parallel, LTE stellar
atmosphere in radiative and hydrostatic equilibrium.  Opacity sampling
over ~30,000 wavelength points replaces the older opacity-distribution
function approach of ATLAS9.

At each iteration the code solves hydrostatic equilibrium, computes
Saha-Boltzmann ionization and molecular equilibrium for all species,
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

# Move into the work directory, which ships with a solar input_model.dat
cd workdir

# Run ATLAS12 on the solar model
../bin/atlas12c.exe sun

# Synthesize a visible spectrum from the converged model
../bin/synthe.exe sun.atm wlbeg=400 wlend=700
```

See [Running ATLAS12](#running-atlas12) and [Running SYNTHE](#running-synthe)
for the full argument syntax, and [Input Data](#input-data) for the data
files required.

## Source Files

| File | Lines | Description |
|------|------:|-------------|
| `atlas12_modules.f90` | 16,449 | Shared modules, all ATLAS subroutines (EOS, opacity, transfer, READIN, JOSH, ...) |
| `atlas12c.f90`        |    623 | ATLAS12 main program (iteration driver) |
| `synthe_module.f90`   |  3,970 | SYNTHE shared data and procedures (hydrogen/He profiles, line opacity, `run_xnfpelsyn`) |
| `mod_mklinelist.f90`  |  1,645 | In-memory line-list preprocessor (replaces standalone mklinelist + `synbeg`/`rgfall`/`rpredict`/`rmolecasc`/`rh2ofast` pipeline) |
| `synthe.f90`          |  1,197 | SYNTHE main program (spectral synthesis driver) |

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
atlas12c.exe [basename] [key=value ...]
```

The input model is read from a file named `input_model.dat` in the
current working directory.  This file contains the standard Kurucz
keyword cards (TEFF, GRAVITY, ABUNDANCE, TURBULENCE, CONVECTION,
READ DECK, ...) terminated by a `BEGIN` card.  The positional
`basename` argument sets the prefix for all output files (default
`mystar`); keyword arguments may appear in any order, before or after
it.

Output files:

| File | Contents |
|------|----------|
| `<basename>.atm`   | Converged model atmosphere (T, P, κ vs. depth) |
| `<basename>.flux`  | Emergent flux vs. wavelength |
| `<basename>.taunu` | Monochromatic optical depth profiles |
| `<basename>.iter`  | Per-iteration summary |
| `<basename>.tcorr` | Temperature correction diagnostics |

Command-line options (keyword=value):

| Option       | Default     | Description |
|--------------|-------------|-------------|
| `numit=N`    | 30          | Number of iterations |
| `vturb=X`    | from model  | Microturbulence (km/s) |
| `mlt=X`      | from model  | Mixing length parameter |
| `teff=X`     | from model  | Rescale model to this Teff (K) |
| `logg=X`     | from model  | Rescale model to this log g (cgs) |
| `zscale=X`   | 1.0         | Metal abundance scale factor (multiplicative on Z≥3) |
| `heabnd=X`   | from model  | He number fraction Y; H is recomputed as X = 1 − Y − Z |
| `abund=file` | none        | Individual element overrides (see below) |

Abundance override file format: one element per line with two
whitespace-separated columns, `Z  log10(number_fraction)`.  Lines
starting with `#` or `!` are treated as comments.  Example:

```
# carbon and iron abundances
 6  -3.52
26  -4.54
```

When `zscale`, `heabnd`, or `abund=` is used, ATLAS12 renormalizes so
that X + Y + Z = 1 (aborting if the specified Y + Z would drive X
negative), then recomputes all abundance-dependent quantities before
the iteration loop.  If both `teff=` and `logg=` are given, the model
is regridded via `SCALE_MODEL` before iteration begins.

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

Arguments:

| Argument              | Required | Description |
|-----------------------|:--------:|-------------|
| `<model_file>`        | yes      | ATLAS12 `.atm` model atmosphere (positional, 1st) |
| `wlbeg=<nm>`          | yes      | Start wavelength in nanometers |
| `wlend=<nm>`          | yes      | End wavelength in nanometers (> wlbeg) |
| `resolu=<R>`          | no       | Resolving power λ/Δλ (default 300 000) |
| `turbv=<kms>`         | no       | Extra microturbulence in km/s added in quadrature to the model value (default 0.0) |
| `more_output=<yes\|no>` | no     | If yes, also write `.linform` and `.mol` diagnostic files (default no).  Accepted truthy values: `yes`, `true`, `1`, `on`, `y` (case-insensitive); falsy: `no`, `false`, `0`, `off`, `n`. |

The output basename is derived from the model filename by stripping the
leading directory and trailing extension.  For example,
`synthe.exe models/sun.atm wlbeg=400 wlend=700` produces `sun.spec` in
the current directory.  Adding `more_output=yes` additionally produces
`sun.linform` and `sun.mol`.

Output files:

| File              | Written                 | Contents |
|-------------------|-------------------------|----------|
| `<base>.spec`     | always                  | ASCII spectrum: wavelength (Å, F11.4), flux (E15.6), continuum flux (E15.6) |
| `<base>.linform`  | only if `more_output=yes` | Per-wavelength diagnostic: wavelength, emergent H, surface H, monochromatic optical depth at each atmospheric layer |
| `<base>.mol`      | only if `more_output=yes` | Molecular number-density diagnostics vs. depth for all species tracked by the equation of state |

Wavelengths are handled internally in nanometers on a logarithmic grid
with spacing `ratio = 1 + 1/resolu`; vacuum wavelengths are used
throughout.  The `.spec` file reports wavelengths in Angstroms for
compatibility with legacy post-processing.

## Example Usage

*TBD — end-to-end worked example (e.g. a solar-model `input_model.dat`,
the ATLAS12 invocation that converges it, and the SYNTHE call that
synthesizes the visible spectrum from the resulting `.atm` file).*

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
| `ionpots.dat`  | `IONPOTS`  | Ionization potentials, all species |
| `isotopes.dat` | `ISOTOPES` | Isotope mass fractions |
| `molecules.dat`| `READMOL`  | Molecular equilibrium constants |
| `pfsaha.dat`   | `PFSAHA`   | Atomic partition functions |
| `pfiron.dat`   | `PFIRON`   | Iron-group partition functions |
| `partfnh2.dat` | `PARTFNH2` | H₂ partition function vs. temperature |

**Continuous opacities**

| File | Used by | Contents |
|------|---------|----------|
| `continua.dat`       | SYNTHE continuum setup | Continuum edge frequency list |
| `crossch.dat`        | `CHOP`       | CH bound-free + bound-bound cross-section table |
| `crossoh.dat`        | `OHOP`       | OH bound-free + bound-bound cross-section table |
| `h2collop.dat`       | `H2COLLOP`   | H₂ collision-induced absorption |
| `hotop.dat`          | `HOTOP`      | Hot-star opacities (high-ionization species) |
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
| `gfpred29dec2014.bin` †  | `read_predict`   | Kurucz predicted atomic lines (Dec 2014) |
| `h2ofastfix.bin` †       | SELECTLINES      | H₂O line list (Partridge & Schwenke) |
| `schwenke.bin` †         | SELECTLINES      | Schwenke diatomic/metal-hydride line list |
| `diatomicspacksrt.bin` † | SELECTLINES      | Packed diatomic molecule line list (sorted) |
| `hilines.bin` †          | SELECTLINES      | Hydrogen/helium line table |
| `lowlines_obs.bin`       | SELECTLINES      | Low-excitation observed atomic lines (symlink → `lowobsat12.bin` †) |
| `lowlines_pl.bin`        | SELECTLINES      | Low-excitation predicted atomic lines (symlink) |
| `nltelines_obs.bin`      | ATLAS12 / XLINOP | NLTE line data (symlink → `nltelinobsat12.bin`) |
| `mol.tar.gz` †           | —                | Archive of molecular sub-lists referenced from `lines.list`; unpack in place |

† Not tracked in the repository; download from the Google Drive folder above.

Several files in this directory are symbolic links that provide stable
"logical" names pointing at versioned data (for example
`lowlines_obs.bin → lowobsat12.bin`,
`nltelines_obs.bin → nltelinobsat12.bin`).  The code always opens the
logical name; versioned targets may be swapped without source changes.

## Translation from Fortran 77

The original ATLAS12 was ~23,000 lines of Fortran 77 fixed-format code.
The original SYNTHE suite comprised three separate programs (XNFPELSYN,
SYNTHE, SPECTRV) communicating through intermediate files.

Key changes in the modernization:

**Structural**

- Free-format Fortran 90 source; explicit typing throughout; `-r8` flag removed
- 58 COMMON blocks consolidated into two modules (`mod_parameters`, `mod_atlas_data`)
- 829 EQUIVALENCE statements eliminated
- Full GOTO elimination pass: 38 → 0 active GOTOs across all four files, including a refactor of the `JOSH` radiative transfer solver
- Command-line argument parsing replaces card-based input for ATLAS12 driver options
- The `mklinelist` line-list preprocessor (formerly the standalone `synbeg | rgfall | rpredict | rmolecasc | rh2ofast` shell pipeline) absorbed into the SYNTHE executable as `mod_mklinelist`, dispatching from a single plain-text `lines.list` manifest
- SYNTHE three-program pipeline (XNFPELSYN, SYNTHE, SPECTRV) merged into a single executable with all inter-stage data flow in memory; intermediate binary scratch files (`fort.9`, `fort.10`, `fort.13`, `fort.14`, `fort.15`, `fort.93`, `spectrv.input`) eliminated
- ASCII spectrum output written directly by SYNTHE (`.spec`), eliminating the standalone `syntoascanga` post-processing step

**Numerical and physical**

- Atmosphere state arrays (T, P, ρ, χ, populations) promoted from REAL(4) to REAL(8)
- Physical constants consolidated into `mod_constants` with CODATA 2018 values, replacing ~130 scattered literals
- Hydrogen line profiles now use the Stehlé & Hutcheon (1999) MMM tabulated Stark broadening profiles via `hydrogen_line_profile`, replacing the Kurucz–Peterson (1982) analytic approximation; the K–P path is retained as a fallback for atmospheric layers exceeding the Inglis–Teller density limit, controlled by the `USE_KP_HYDROGEN` flag
- H₂ partition function read from an external data file rather than hardcoded
- Convergence improvements for cool dwarf models (~2800 K): adiabatic sweep in `DTCONV`, iteration damping in convective layers, gap-filling in `CONVEC`

**Data and code hygiene**

- Hardcoded data arrays (partition functions, Feautrier matrices, Karzas–Latter tables, ionization potentials, isotope fractions) externalized to readable data files
- Dead-code audit with static call-graph analysis: ~233 lines removed from `atlas12_modules.f90`
- He I line-profile island (`HE1_GENERIC_PROFILE` and 11 helpers, ~904 lines) retained but marked `PORT INCOMPLETE` and not currently called; `stark_quasistatic_profile` handles He I until that port is completed
- All remaining Fortran 2008 `BLOCK...END BLOCK` constructs eliminated
- All routines carry descriptive header comments

## References

- Kurucz, R. L. 1970, SAO Special Report 309
- Kurucz, R. L. 1993, ATLAS9 Stellar Atmosphere Programs, CD-ROM No. 13
- Kurucz, R. L. 2005, Memorie della Società Astronomica Italiana Supplementi, 8, 14
- Sbordone, L., Bonifacio, P., Castelli, F., & Kurucz, R. L. 2004, MSAIS, 5, 93
- Castelli, F., & Kurucz, R. L. 2004, astro-ph/0405087 (new grids of ATLAS9 model atmospheres)
- Castelli, F. 2005, MSAIS, 8, 25 (ATLAS12: how to use it)
