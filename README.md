# ATLAS12 + SYNTHE

Fortran 90 translation of Robert L. Kurucz's stellar atmosphere and spectral
synthesis codes, originally developed at the Harvard-Smithsonian Center for
Astrophysics.

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

## Source Files

| File | Lines | Description |
|------|------:|-------------|
| `atlas12_modules.f90` | 14,789 | Shared modules, all subroutines (EOS, opacity, transfer, etc.) |
| `atlas12c.f90` | 623 | ATLAS12 main program (iteration driver) |
| `synthe_module.f90` | 3,498 | SYNTHE shared data and procedures (line profiles, abundances) |
| `synthe.f90` | 1,105 | SYNTHE main program (spectral synthesis driver) |

## Building

Requires `gfortran`.  From the source directory:

```
make              # build both atlas12c.exe and synthe_spectrv.exe
make atlas        # build atlas12c.exe only
make synthe       # build synthe_spectrv.exe only
```

Executables are installed to `../bin/`.  The source directory is left
clean (no `.o` or `.mod` files).

## Running ATLAS12

```
export ATLAS12=/path/to/atlas12/
atlas12c.exe mystar [options]
```

The model atmosphere is read from `input_model.dat` in the working
directory.  Output files use the base name given on the command line
(default `mystar`):

| File | Contents |
|------|----------|
| `mystar.atm` | Converged model atmosphere (T, P, κ vs. depth) |
| `mystar.flux` | Emergent flux vs. wavelength |
| `mystar.taunu` | Monochromatic optical depth profiles |
| `mystar.iter` | Per-iteration summary |
| `mystar.tcorr` | Temperature correction diagnostics |

Command-line options (keyword=value):

| Option | Default | Description |
|--------|---------|-------------|
| `numit=N` | 30 | Number of iterations |
| `vturb=X` | from model | Microturbulence (km/s) |
| `mlt=X` | 2.0 | Mixing length parameter |
| `teff=X` | from model | Rescale to this Teff |
| `logg=X` | from model | Rescale to this log g |
| `zscale=X` | 1.0 | Metal abundance scale factor |
| `heabnd=X` | from model | He number fraction |
| `abund=file` | none | Individual element overrides (Z  log_abund) |

## Running SYNTHE

The synthesis workflow has three steps:

1. **SYNBEG** — sets wavelength range and resolution, writes `fort.93`
2. **Line readers** (RGFALLLINELIST, etc.) — select lines, write `fort.12`/`fort.14`
3. **SYNTHE_SPECTRV** — reads the model and line data, computes the spectrum

```
synthe_spectrv.exe <model.atm> <output_basename>
```

Input files:

| File | Description |
|------|-------------|
| `model.atm` | ATLAS12 model atmosphere (unit 5) |
| `fort.93` | Run parameters from SYNBEG |
| `fort.12` | Preprocessed line data |
| `spectrv.input` | SPECTRV parameters (in `$ATLAS12/data/`) |
| `continua.dat` | Continuum edge frequencies (in `$ATLAS12/data/`) |

Output files:

| File | Contents |
|------|----------|
| `<basename>.spec` | ASCII spectrum: wavelength, flux, continuum flux |
| `<basename>.linform` | Line identifications |

## Data Directory

Set the `$ATLAS12` environment variable to the installation root.
Data files are read from `$ATLAS12/data/`.  If unset, defaults to `./data/`.

Key data files:

| File | Used by | Contents |
|------|---------|----------|
| `molecules.dat` | READMOL | Molecular equilibrium constants |
| `pfiron.dat` | PFIRON | Iron-group partition functions |
| `blockj.dat`, `blockh.dat` | BLOCKJ, BLOCKH | Feautrier coefficient matrices |
| `ionpots.dat` | IONPOTS | Ionization potentials |
| `isotopes.dat` | ISOTOPES | Isotope mass fractions |
| `karzas_*.dat` | XKARZAS | Hydrogen bound-free cross sections |
| `nltelines_obs.bin` | XLINOP | NLTE line data |
| `spectrv.input` | SPECTRV | Synthesis run parameters |
| `continua.dat` | SPECTRV | Continuum edge frequency list |

## Translation from Fortran 77

The original ATLAS12 was ~23,000 lines of Fortran 77 fixed-format code.
The original SYNTHE suite comprised three separate programs (XNFPELSYN,
SYNTHE, SPECTRV) communicating through intermediate files.

Key changes in the modernization:

- 58 COMMON blocks consolidated into two modules (`mod_parameters`, `mod_atlas_data`)
- 829 EQUIVALENCE statements eliminated
- Data arrays (partition functions, coefficients) externalized to readable text files
- SYNTHE three-program pipeline merged into a single executable with in-memory data flow
- All intermediate file I/O between SYNTHE stages eliminated
- Explicit typing throughout; `-r8` flag removed
- Free-format source; descriptive comments on all routines
- Command-line argument parsing replaces card-based input for ATLAS12

## References

- Kurucz, R. L. 1970, SAO Special Report 309
- Kurucz, R. L. 1993, ATLAS9 Stellar Atmosphere Programs, CD-ROM No. 13
- Kurucz, R. L. 2005, Memorie della Società Astronomica Italiana Supplementi, 8, 14
- Sbordone, L., et al. 2004, MSAIS, 5, 93
