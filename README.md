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
| `<basename>.iter`  | Per-iteration summary, including the temperature-correction diagnostics (the former `.tcorr` columns are merged into this file) |

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

Numerics A/B switches are deliberately not CLI options.  They live as
developer flags at their declarations in `mod_atlas_data` (set the
default value and recompile): `USE_CZ_CONSTRUCTOR` (deep-CZ temperature
constructor), `USE_CZC_POLISH` (terminal deep-CZ flux-closure polish),
`USE_FLXCNV_SMOOTH` (interior 1-2-1 convective-flux smoothing),
`IROSSTAB` (Rosseland-table interpolation: 1=bilinear, 2=Shepard,
3=moving least squares), and `IQUAD` (`INTEG` quadrature: 0=legacy
blended-parabola, 1=Steffen monotone cubic).

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

## Example Usage

*TBD — end-to-end worked example (e.g. a solar input atmosphere file,
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
| `gfpred29dec2014.bin` †  | `read_predict`, SELECTLINES | Kurucz predicted atomic lines (Dec 2014) |
| `h2ofastfix.bin` †       | SELECTLINES      | H₂O line list (Partridge & Schwenke) |
| `schwenke.bin` †         | SELECTLINES      | Schwenke diatomic/metal-hydride line list |
| `diatomicspacksrt.bin` † | SELECTLINES      | Packed diatomic molecule line list (sorted) |
| `hilines.bin` †          | SELECTLINES      | Hydrogen/helium line table |
| `lowobsat12.bin` †       | SELECTLINES      | Low-excitation observed atomic lines |
| `nltelinobsat12.bin` †   | ATLAS12 / XLINOP | NLTE line data |
| `mol.tar.gz` †           | —                | Archive of molecular sub-lists referenced from `lines.list`; unpack in place |

† Not tracked in the repository; download from the Google Drive folder above.

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
- Atomic partition functions replaced with the Barklem & Collet (2016) tabulated U(T) calculated by direct NIST level summation, covering Z = 1..92 in ionization stages I–III; the legacy Kurucz hand-selected few-term sums used by `PFGROUND` (the low-T floor in `PFSAHA`) and the generic `NNN` interpolation table are retained as an internal safety net but are no longer the production data source.  Iron-group partition functions (Z = 20..28) still flow through `PFIRON` for its pressure-lowering correction, but the unperturbed POTLOW = 0 baseline is now anchored to B&C
- Metal continuum bound-free opacity overhauled: a new `CONT_METAL_OPACITY_TOPBASE` dispatcher replaces the legacy collection of analytic Seaton/Peach fits (Li1OP, C1OP, MG1OP, AL1OP, SI1OP, FE1OP, …).  `MBF_TOPBASE` supplies 30 species (neutrals and first ions of Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, S, Ar, Ca) from R-matrix Opacity Project cross sections processed by Allende Prieto, Hubeny & Lambert (2003, ApJS 147, 363) into resonance-averaged 3%-log-spaced grids; `MBF_HIGH_ION` supplies higher ionization stages and Coulomb free-free from a filtered subset of `hotop.dat`; `FELO_OPACITY` supplies Fe I and Fe II from Iron Project / Opacity Project R-matrix calculations (Bautista 1997 for Fe I, Nahar & Pradhan 1994 for Fe II, sourced from the TLUSTY model atoms).  Toggleable via `USE_TOPBASE_MBF`; the legacy path is preserved for reference.  Known limitation: TOPbase identifies levels by LS term without resolving J, so per-species opacity carries a U_TB/U_BC factor of a few percent for most species (≲13% for the worst case, Ar II)
- Solar reference abundance pattern is now selectable at the command line (`solar=name`).  The compiled-in default remains Anders & Grevesse (1989), defined once in `mod_atlas_data` and shared by ATLAS12 and SYNTHE; `solar=` replaces the abundance table carried by the input model file while preserving its relative [M/H] offsets, and `zscale=`/`abund=`/`heabnd=` overrides apply on top.  Additional named scales are stored as log-ε tables and converted internally (`LOGEPS_TO_ABUND`); registered scales: `ag89` (Anders & Grevesse 1989), `agss09` (Asplund et al. 2009), `berg25` (Bergemann, Lodders & Palme 2025, with A(He) = 10.922 from the published compilation).  (An earlier changelog entry here claimed the default had been switched to Asplund et al. 2009 — that change was never made; the code default has always been AG89.)
- Hydrogen line profiles now use the Stehlé & Hutcheon (1999) MMM tabulated Stark broadening profiles via `hydrogen_line_profile`, replacing the Kurucz–Peterson (1982) analytic approximation; the K–P path is retained as a fallback for atmospheric layers exceeding the Inglis–Teller density limit, controlled by the `USE_KP_HYDROGEN` flag
- Voigt function evaluated via Weideman (1994) N=32 rational approximation of w(z) plus Humlíček (1982) R₁,₂ asymptotic, replacing the Kurucz three-regime table approximation that had ~4% error at the a = 0.2 regime boundary.  Accuracy is now ~10⁻⁵ relative across the full (a, v) plane with no regime boundaries
- H₂ partition function read from an external data file rather than hardcoded
- The `JOSH` scattering solution is now a direct dense solve of the 51-point integral-operator system (I − diag(α)·Λ_J)·S = (1−α)·S̄ via `SOLVIT`, replacing the Gauss–Seidel Λ-iteration.  The iteration stalled against its 50-sweep cap once the scattering albedo α exceeded ~0.99 across the grid — 0.3% error in S at α = 0.99, 10% at α = 0.999, accepted after a warning — and its 10⁻⁵ early-exit tolerance quantized the response of every scattering solve, leaving discontinuous noise in H_ν at the tolerance level.  The direct solve is exact at comparable cost (one 51×51 LU per frequency) and reproduces the iteration to its own tolerance wherever the iteration converged (`workdir/pfverify/josh_verify.py`, Test D)
- Optional Steffen (1990) monotone-cubic quadrature for `INTEG` (developer flag `IQUAD = 1` in `mod_atlas_data`), the routine that builds τ_Ross, the per-frequency τ_ν, P_rad, and the height scale.  The legacy blended-parabola quadrature (`PARCOE`) integrates a 5-decade exponential on real model grids with 1.3% (solar) to 11% (2800 K) maximum error — only ~2× better than a trapezoid — while the Steffen cubic is 3–6× more accurate and its monotone interpolant keeps τ scales monotone with no overshoot (`workdir/pfverify/josh_verify.py`, Test C).  The quadrature error is concentrated in the deepest layers, so switching moves converged structures mainly at the grid bottom: ≤1.5 K in cool/warm stars (convection pins the deep adiabat) but 0.1–0.5% of T in hot-star radiative envelopes (47–105 K at 8000–40000 K), with equal convergence quality and negligible emergent-flux change.  The default remains the legacy quadrature
- Partition-function tables (H₂ and Barklem & Collet) evaluated with natural cubic splines instead of (log-)linear interpolation.  The interpolation kinks at table nodes, sampled by the EOS finite-difference stencil, imprinted spurious structure on ∇_ad that F_conv ∝ (∇−∇_ad)^{3/2} amplified into large local flux errors in cool-dwarf deep convection zones.  Node values, boundary behavior, and hot-star fallbacks are unchanged (verified by the `workdir/pfverify/` battery); converged structures shift by ≲0.5 K
- Deep convection-zone convergence overhauled ("CZ constructor", `CZ_CONSTRUCT`): where convection is efficient the classical flux-error correction is structurally dead, so the required gradient is obtained by inverting the Böhm–Vitense relation in closed form for the superadiabatic excess carrying F_conv = F_tot − F_rad, applied as incremental interval-space corrections blended into the Kurucz correction by convective flux fraction.  Engagement is error-gated; release is on step size (ΔT < 1.5 K), standard practice for efficient convection zones.  Cool dwarfs (2500–3500 K) now converge from cold starts in ~10–30 iterations to few-percent flux balance — grid points that previously never converged; solar and hot-star models never engage it.  Developer flags: `USE_CZ_CONSTRUCTOR`, `USE_FLXCNV_SMOOTH` in `mod_atlas_data`
- Terminal deep-CZ polish (`CZC_POLISH`, developer flag `USE_CZC_POLISH`, default on): on the final iteration, the CZ-constructor block's temperatures are re-placed by a Newton solve of the node-gradient system ∇_j(T) = ∇_ad + ∇_required (tridiagonal finite-difference Jacobian, banded LU; the interval-sweep construction is retained as a fallback), closing the MLT flux relation against the final radiative flux to ≲0.1% where the iterative floor was 5–30%.  Two structural facts force this design.  First, deep-CZ flux balance is stiff: F_conv ∝ (∇−∇_ad)^{3/2} at ∇−∇_ad ≈ 10⁻³ demands milli-Kelvin temperature placement — three orders below the constructor's release threshold, and unreachable during the iterations because each TAUSTD remap re-perturbs T at the Kelvin level (which is also why deep-CZ flux errors of tens of percent are endemic to published cool-dwarf models).  Second, ∇_ad's rapidly varying curvature through the H₂-dissociation region (the EOS chain itself is smooth to ~10⁻⁷, verified by the fine-T probe `GRDADB_SCAN_MAYBE`, env `ATLAS_GRDADB_SCAN`) defeats interval-space constructions: the interval-slope vs node-average operator mismatch has an alternating component invisible to averaged relaxation, which DEL^{3/2} stiffness amplifies to ~16% flux error.  Two further constraints shape where and how the placement is applied.  The polish closes the *raw* MLT flux relation, which equals the model's actual (1-2-1 smoothed, gap-filled) convective flux only where convection is efficient, so its authority is a smoothstep in convective flux fraction over `CZC_POL_W_LO`…`CZC_POL_W_HI` (0.85…0.95): full in the deep interior — which contains the entire flux-balance sawtooth — tapering to zero at the CZ-top transition, which is left to the normal iteration.  (Applying it across the transition instead degrades those layers from <1% to several percent, and the final consistency refresh must likewise use the full production `CONVEC`, not the MLT-only path, or the `.atm` convective-flux columns are overwritten with raw MLT that switches convection off near τ≈1.)  And because re-placing the block in isolation shifts its temperatures — by several K where the constructor left the deep CZ far from closure — against a still-unrelaxed atmosphere above, the polish runs on the last `CZC_POL_NHEAL` (8) iterations rather than only the final one: between calls the normal flux correction, which is alive in the transition, relaxes those layers toward the polished deep boundary, healing the boundary seam (3500 K: 7% → a converged 1.8%, with the final polish then applying ≈0 ΔT).  Verified on 2500/2800/3500 K: the deep convection zone closes from ±11–34% to <0.04%, with the CZ-top transition left at a normal 1–2%; solar and hot-star models are untouched (the polish only runs where the constructor engaged).  A closure summary goes to the run log, the post-polish per-layer state is appended to `.iter` as a final block in the standard table layout (`ERROR` is the flux residual at frozen F_rad, `T1` the polish ΔT, correction columns zero by construction), and the final `.atm` carries the polished structure
- Convective boundary fixes: marginally subadiabatic layers at the grid bottom take the flux-conservation value F_c = F_tot − F_rad instead of zero (Koester 1980), and the deepest-layer gradient uses the log-space one-sided slope, removing a ~13% systematic overestimate of ∇ — the long-standing "deepest layer never converges" artifact

- The `.atm` deck gains a `TAU5000` column: the continuum optical depth at 5000 Å (continuum absorption + continuum scattering at 500 nm from `KAPP`, integrated on the same quadrature as τ_Ross), the reference depth scale tabulated by MARCS and PHOENIX; computed once at output time from the final structure.  The deck-count field is widened `I3`→`I4`, fixing the round-trip failure of models with ≥ 100 depth points (the count printed glued to `DECK6` and could not be parsed back by `READIN`)

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
- Allende Prieto, C., Hubeny, I., & Lambert, D. L. 2003, ApJS, 147, 363 (TOPbase metal photoionization grids)
- Anders, E., & Grevesse, N. 1989, Geochim. Cosmochim. Acta, 53, 197 (default solar abundance scale, `solar=ag89`)
- Asplund, M., Grevesse, N., Sauval, A. J., & Scott, P. 2009, ARA&A, 47, 481 (solar abundance scale, `solar=agss09`)
- Bergemann, M., Lodders, K., & Palme, H. 2025, Zenodo record 14988840 (solar abundance scale, `solar=berg25`)
- Barklem, P. S., & Collet, R. 2016, A&A, 588, A96 (atomic and molecular partition functions)
- Bautista, M. A. 1997, A&AS, 122, 167 (Fe I R-matrix bound-free)
- Humlíček, J. 1982, JQSRT, 27, 437 (Voigt function asymptotic)
- Nahar, S. N., & Pradhan, A. K. 1994, J. Phys. B, 27, 429 (Fe II R-matrix bound-free)
- Stehlé, C., & Hutcheon, R. 1999, A&AS, 140, 93 (hydrogen Stark profiles)
- Weideman, J. A. C. 1994, SIAM J. Numer. Anal., 31, 1497 (Voigt function rational approximation)
