# ATLAS12

Opacity-sampling stellar atmosphere code, originally written by Robert L. Kurucz at the Harvard-Smithsonian Center for Astrophysics. This version has been modernized from Fortran 77 to Fortran 90 by Charlie Conroy.

ATLAS12 iteratively solves for the structure of a plane-parallel, LTE stellar atmosphere in radiative and hydrostatic equilibrium, using opacity sampling over ~30,000 wavelength points.

## Building

```bash
cd $ATLAS12
make
```

The executable is built as `bin/atlas12c.exe`. The `$ATLAS12` environment variable should point to the installation root directory, with atomic/molecular data files in `$ATLAS12/data/`.

## Quick Start

```bash
# Run with default parameters (reads model from input_model.dat)
atlas12c.exe mystar

# Specify Teff, logg, and number of iterations
atlas12c.exe mystar teff=5500 logg=4.50 numit=30

# Scale metallicity to [Fe/H] = -1.0
atlas12c.exe mystar teff=5500 logg=4.50 zscale=0.1

# Full example with all options
atlas12c.exe mystar teff=5500 logg=4.50 numit=30 vturb=2.0 mlt=1.5 \
    zscale=0.1 heabnd=0.078 abund=my_elements.dat
```

## Command-Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `basename` | Output file base name (positional) | `mystar` |
| `numit=N` | Number of iterations | `30` |
| `teff=X` | Effective temperature (K); rescales model | from model |
| `logg=X` | Log surface gravity (cgs); rescales model | from model |
| `vturb=X` | Microturbulent velocity (km/s) | from model |
| `mlt=X` | Mixing length parameter (l/H_p) | `2.0` |
| `zscale=X` | Metal abundance scale factor (linear) | no scaling |
| `heabnd=X` | Helium number fraction | from model |
| `abund=file` | File with individual element abundance overrides | none |

When `teff=` and/or `logg=` are specified, the input model is regridded onto a standard 80-point log(tau_Ross) grid (from log tau = -6.875 to +2.875 in steps of 0.125) and the temperature/pressure structure is rescaled to the new parameters.

When `zscale=` is specified, all metal abundances (Z >= 3) are scaled by the given factor. The hydrogen abundance is automatically recomputed as X = 1 - Y - Z to maintain consistency. If `heabnd=` is also given, it sets the helium number fraction before the hydrogen recomputation.

## Input Files

### `input_model.dat`

The starting model atmosphere, read automatically at startup. This file uses the Kurucz card format and must contain at minimum:

- `TEFF` and `GRAVITY` cards
- `ABUNDANCE TABLE` with all 99 elements
- `READ DECK6` with the model structure (RHOX, T, P, XNE, ABROSS, ACCRAD, VTURB, FLXCNV, VCONV)
- `PRADK` surface radiation pressure
- `BEGIN` to signal end of model data

A previously converged `.atm` output file can be used directly as `input_model.dat`.

### Abundance override file (`abund=`)

Optional plain text file with individual element abundance overrides. Format is one element per line with two columns:

```
# Z   log10(number_fraction)
6    -3.61
8    -3.31
26   -4.54
```

Lines that cannot be parsed (e.g., comments starting with `#`) are skipped. Elements specified here override the metallicity scaling from `zscale=`.

### Data files (`$ATLAS12/data/`)

Atomic, molecular, and opacity data files located in `$ATLAS12/data/`. Files marked **(required)** must be present; files marked **(optional)** add additional opacity sources if available.

**Radiative transfer and structure:**

- `blockj.dat` — Feautrier J-moment coefficient matrices **(required)**
- `blockh.dat` — Feautrier H-moment coefficient matrices **(required)**
- `ionpots.dat` — ionization potentials for all species **(required)**
- `isotopes.dat` — isotope mass fractions **(required)**
- `molecules.dat` — molecular equilibrium constants and species codes **(required)**
- `pfsaha.dat` — partition functions for Saha ionization equilibrium **(required)**
- `pfiron.dat` — iron partition function data **(required)**
- `partfnh2.dat` — H₂ partition function (Barklem & Collet 2016) **(required)**

**Continuous opacity:**

- `karzas_ekarzas.dat` — Karzas-Latter bound-free gaunt factors **(required)**
- `karzas_xn.dat` — Karzas-Latter cross-section data **(required)**
- `karzas_freqn.dat` — Karzas-Latter frequency grid **(required)**
- `karzas_xl.dat` — Karzas-Latter angular momentum data **(required)**
- `crossch.dat` — CH molecule photo-dissociation cross-sections **(required)**
- `crossoh.dat` — OH molecule photo-dissociation cross-sections **(required)**
- `h2collop.dat` — H₂ collision-induced opacity **(required)**
- `hotop.dat` — high-ionization opacity data **(required)**
- `stark_profile.dat` — Stark broadening profile tables **(required)**

**Line opacity (optional, cascading):**

These files are opened in sequence; missing files are skipped. Each adds an additional class of line opacity.

- `lowlines_pl.dat` — predicted atomic lines (low excitation) **(optional)**
- `lowlines_obs.bin` — observed atomic lines (low excitation) **(optional)**
- `nltelines_obs.bin` — pre-computed NLTE line data **(optional)**
- `hilines.bin` — high-excitation atomic lines **(optional)**
- `diatomicspacksrt.bin` — diatomic molecular lines **(optional)**
- `schwenke.bin` — TiO/H₂O lines (Schwenke) **(optional)**
- `h2ofastfix.bin` — H₂O lines (fast computation) **(optional)**
- `h3plus.dat` — H₃⁺ lines **(optional)**

## Output Files

All output files use the base name specified on the command line:

| File | Description |
|------|-------------|
| `basename.atm` | Converged model atmosphere (can be re-used as input) |
| `basename.flux` | Surface flux spectrum |
| `basename.taunu` | Monochromatic optical depths |
| `basename.iter` | Full model structure at each printed iteration |
| `basename.tcorr` | Temperature correction diagnostics |

### Model atmosphere format (`.atm`)

The `.atm` file contains the converged model in a self-describing format:

```
TEFF  5500.  GRAVITY  4.5000 LTE
TITLE ATLAS12 l/H= 2.00
 ABUNDANCE TABLE
    1H   0.920440       2He  0.078340
    3Li-10.990 0.000    4Be-10.620 0.000    5B  -9.250 0.000  ...
READ DECK6 80        RHOX         T         P       XNE    ABROSS    ACCRAD     VTURB    FLXCNV     VCONV
              5.26119E-04   3449.70 1.663E+01 2.256E+09 2.533E-04 9.034E-02 2.000E+05 0.000E+00 0.000E+00
              ...
PRADK 1.2144E+00
BEGIN                    ITERATION  30 COMPLETED
```

Columns in the model structure:

| Column | Variable | Description | Units |
|--------|----------|-------------|-------|
| 1 | RHOX | Column mass density | g/cm² |
| 2 | T | Temperature | K |
| 3 | P | Gas pressure | dyn/cm² |
| 4 | XNE | Electron number density | cm⁻³ |
| 5 | ABROSS | Rosseland mean opacity | cm²/g |
| 6 | ACCRAD | Radiative acceleration | cm/s² |
| 7 | VTURB | Microturbulent velocity | cm/s |
| 8 | FLXCNV | Convective flux | erg/cm²/s |
| 9 | VCONV | Convective velocity | cm/s |

## Physics

ATLAS12 computes model atmospheres under the following assumptions:

- **Plane-parallel geometry** — suitable for dwarf and giant stars; not appropriate for supergiants or extended atmospheres.
- **Hydrostatic equilibrium** — gas pressure + radiation pressure + turbulent pressure balance gravity.
- **Radiative equilibrium** — enforced via temperature corrections (Lambda iteration on the flux integral).
- **Local thermodynamic equilibrium (LTE)** — source function equals the Planck function; populations follow Saha-Boltzmann.
- **Opacity sampling** — continuous and line opacities are computed at ~30,000 wavelength points spanning the UV through IR.
- **Mixing-length convection** — convective energy transport using the Böhm-Vitense formulation with optional overshooting.
- **Molecular equilibrium** — simultaneous solution for ionization and molecular equilibrium including H₂, CO, TiO, H₂O, and other species.

## History

ATLAS12 was written by R.L. Kurucz as a successor to ATLAS9, replacing opacity distribution functions with direct opacity sampling for improved accuracy in cool star atmospheres where molecular line blanketing is important. This version has been ported from Fortran 77 to Fortran 90 with modernized input handling, named constants, and explicit variable typing.

## References

- Kurucz, R.L. 1970, SAO Special Report 309 (original ATLAS code)
- Kurucz, R.L. 1993, ATLAS9 Stellar Atmosphere Programs (CD-ROM 13)
- Kurucz, R.L. 2005, Memorie della Società Astronomica Italiana Supplementi, 8, 14 (ATLAS12 description)
- Castelli, F. & Kurucz, R.L. 2003, IAU Symposium 210, A20 (New grids of ATLAS9 models)
