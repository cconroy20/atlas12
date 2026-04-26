!=========================================================================
!  ATLAS12  —  Opacity-Sampling Stellar Atmosphere Code
!
!  Robert L. Kurucz, Harvard-Smithsonian Center for Astrophysics
!  Translated from Fortran 77 (fixed-form) to Fortran 90 (free-form),
!
!  Data files are located via the environment variable
!  ATLAS12.  If unset, defaults to ./data/
!       export ATLAS12=/path/to/atlas12/
!
!  LIMITATIONS:
!    - The iron PF and CH, OH cross sections have a low-T limit of 2000K
!    - The OS grid only extends to ~14.5um, which might not be red enough
!      for very cool stars
!
!  CHANGE LOG:
!    2/26: initial translation from F77 to F90 (C. Conroy w/ Claude.ai)
!
!=========================================================================
!
!  STRUCTURE OUTLINE
!  -----------------
!
!  MODULES
!    mod_parameters              Compile-time constants (kw, mion, maxmol, ...)
!    mod_constants               Physical constants (CODATA 2018) and
!                                dissolved-level support data (Holtsmark Q
!                                table, Hummer-Mihalas beta coefficient)
!    mod_partition_functions     Barklem & Collet (2016) atomic partition
!                                functions; U_BC accessor + lazy file load
!    mod_atlas_data              All shared state (replaces 58 COMMON blocks)
!                                plus the main physics/numerics routines
!
!  MAIN PROGRAM  (in atlas12c.f90)
!    ATLAS12                     Iteration driver: read input, iterate atmosphere
!
!  --- Procedures in mod_partition_functions -----------------------------
!    set_bc_data_dir             Set directory for partfn_bc2016.dat
!    init_bc_partition_functions Lazy one-time load of the B&C table
!    parse_t_grid                Parser: T grid header line
!    parse_species_record        Parser: one species line (label + 42 U(T))
!    classify_label              Parser: "Fe_I"/"H-"/... -> (Z, ION) or neg-ion
!    element_Z                   Atomic symbol -> Z
!    roman_to_int                Roman numeral -> integer
!    U_BC                        Positive-ion partition function U(Z,ION,T)
!    U_BC_NEG                    Negative-ion partition function
!    BC_HAS                      Coverage query (Z, ION)
!    BC_HAS_NEG                  Coverage query (Z, negative ion)
!    interp_log_linear_cached    Log-linear T interpolation with bracket hint
!
!  --- Procedures in mod_atlas_data --------------------------------------
!
!  MODEL OUTPUT
!    PUTOUT                      Write model atmosphere, fluxes, intensities
!
!  TEMPERATURE CORRECTION & STRUCTURE
!    TCORR                       Temperature correction via flux conservation
!    STATEQ                      Statistical equilibrium (H NLTE rates)
!    RADIAP                      Radiative acceleration iteration
!    ROSS                        Rosseland mean opacity computation
!    CONVEC                      Convective flux by mixing-length theory
!    COMPUTE_PTURB               Turbulent pressure computation
!    VTURB_VARYDEPTH             Standard microturbulent velocity profile
!    COMPUTE_HEIGHT              Geometric depth from hydrostatic integration
!    TTAUP                       Starting model: T-tau-P relation
!
!  EQUATION OF STATE
!    NELECT                      Electron density by charge conservation
!    COMPUTE_ONE_POP             Single-species population (Saha-Boltzmann)
!    COMPUTE_ALL_POPS            Populations for all species
!    PFSAHA                      Partition functions and Saha ionization
!    pfsaha_highlevels           Hydrogenic high-n Rydberg tail correction
!    PFGROUND_HYBRID             Atomic partition function dispatcher (B&C/Kurucz)
!    PFGROUND_KURUCZ             Kurucz hand-coded partition functions (fallback)
!    PFIRON                      Iron-group partition function interpolation
!    MOLEC                       Single-molecule equilibrium
!    NMOLEC                      Full molecular equilibrium solution
!    READMOL                     Read molecular equilibrium input data
!    EQUILH2                     H2 equilibrium constant
!    PARTFNH2                    H2 molecular partition function
!
!  RADIATIVE TRANSFER
!    JOSH                        Source function solver (Feautrier method)
!    ENERGY_DENSITY              Radiation energy density integration
!    BLOCKJ                      J-coefficient matrix (loaded from file)
!    BLOCKH                      H-coefficient matrix (loaded from file)
!
!  CONTINUOUS OPACITY  — Hydrogen
!    HOP                         H I bound-free + free-free
!    XKARZAS                     Karzas-Latter bound-free cross sections
!    COULX                       Hydrogenic photoionization cross section
!    COULBF1S                    Ground-state Coulomb bound-free
!    COULFF                      Coulomb free-free gaunt factor
!    HMINOP                      H- bound-free + free-free
!    H2PLOP                      H2+ opacity
!    H2RAOP                      H2 Rayleigh scattering
!    H2COLLOP                    H2 collision-induced absorption
!    HRAYOP                      H I Rayleigh scattering
!    occupation_prob             Hummer-Mihalas (1988) level occupation w_n
!    holtsmark_Q                 Holtsmark microfield cumulative Q(beta)
!
!  CONTINUOUS OPACITY  — Helium
!    HE1OP                       He I bound-free (with satellites)
!    CROSSHE                     He I photoionization cross section
!    HE111S, HE12S1S, ...        He I level-specific cross sections
!    HE2OP                       He II bound-free + free-free
!    HEMIOP                      He- free-free
!    HERAOP                      He I Rayleigh scattering
!
!  CONTINUOUS OPACITY  — Metals
!    C1OP, C2OP                  Carbon  I–II
!    O1OP                        Oxygen   I
!    MG1OP, MG2OP                Magnesium I–II
!    AL1OP                       Aluminum I
!    SI1OP, SI2OP                Silicon  I–II
!    CA2OP                       Calcium  II
!    FE1OP                       Iron I
!    CHOP                        CH molecular opacity
!    OHOP                        OH molecular opacity
!    SEATON                      Seaton formula for photoionization x-section
!
!  CONTINUOUS OPACITY  — Assembly
!    KAPP                        Total continuous opacity from all sources
!    KAPCONT                     Continuous opacity tabulation
!    CONT_METAL_OPACITY_LEGACY   Consolidated metal/molecular/CIA continuum
!                                (replaces former COOLOP + WARMOP + HOTOP)
!    ELECOP                      Electron scattering
!    XCONOP                      Tabulated continuum opacity lookup
!    ROSSTAB                     Rosseland opacity table interpolation
!
!  LINE OPACITY
!    SELECTLINES                 Line selection for current frequency
!    LINOP1                      Line opacity computation (single frequency)
!    XLINOP                      Extended line opacity computation
!    IONPOTS                     Load ionization potentials
!    ISOTOPES                    Load isotope data
!
!  LINE PROFILES
!    HLINOP                      Hydrogen line opacity (Balmer, Lyman, ...)
!    STARK                       Kurucz analytic Stark broadening (VCS)
!    STARK_MMM                   Stehle (1999) MMM tabulated Stark broadening
!    INIT_STARK_TABLES           Lazy load of the Stehle binary tables
!    hydrogen_f_value            Hydrogen oscillator strength f(n,m)
!    TABVOIGT                    Voigt profile tabulation
!    HFNM                        Hydrogen f(n,m) (Kurucz convention)
!    VCSE1F                      VCS electric field distribution
!    STARK_PROFILE               Stark profile wing function (formerly SOFBET)
!    EXPI                        Exponential integral E_n(x)
!
!  NUMERICAL UTILITIES
!    DERIV                       Differentiation by parabolic interpolation
!    INTEG                       Integration using parabolic coefficients
!    PARCOE                      Parabolic interpolation coefficients
!    MAP1                        Linear interpolation/remapping
!    MAP4                        Cubic interpolation/remapping
!    SOLVIT                      LU decomposition linear equation solver
!    LINTER                      Linear interpolation utility
!
!  I/O & PARSING
!    READIN                      Read and parse all input control cards
!    SCALE_MODEL                 Rescale existing model to new Teff/logg
!    FREEFF                      Free-format floating point reader
!    FREEFR                      Free-format real number reader
!    NEXTWORD                    Free-format word reader (returns character string)
!    DUMP_ARRAY                  Debug array print utility
!
!=========================================================================

!=========================================================================
! mod_parameters: Compile-time constants for ATLAS12
!=========================================================================

MODULE mod_parameters

  IMPLICIT NONE
  INTEGER, PARAMETER :: kw = 80        ! Number of atmospheric depth points
  INTEGER, PARAMETER :: mion = 1006    ! Number of ion species
  INTEGER, PARAMETER :: maxmol = 200   ! Maximum number of molecules
  INTEGER, PARAMETER :: max1 = maxmol + 1
  INTEGER, PARAMETER :: maxeq = 35     ! Maximum number of equilibrium equations
  INTEGER, PARAMETER :: maxloc = 3 * maxmol

  ! Debug flag: set to 1 to print subroutine entry tracing to unit 6
  INTEGER, PARAMETER :: IDEBUG = 0

END MODULE mod_parameters

!=========================================================================
! mod_constants: Physical constants and derived quantities
!
!   Fundamental constants: CODATA 2022 / 2019 SI redefinition.
!   k_B, h, c, e_charge (in C), and 1 eV are exact by the 2019 SI
!   redefinition.  The remaining quantities (m_u, m_e, sigma_T, R_inf,
!   e in esu) carry small experimental uncertainties that are
!   negligible for stellar atmosphere work but kept to full tabulated
!   precision so future recompilations match the literature.
!
!   Derived constants: pre-computed combinations that appear throughout
!   ATLAS12/SYNTHE (Planck function prefactor, Saha prefactor, free-free
!   coefficient, etc.).  Defining them here once eliminates scattered
!   magic numbers and ensures self-consistency.
!
!   Spectroscopic constants: measured quantities (Rydberg, H⁻ electron
!   affinity, H₂ dissociation energy) that are not fundamental but appear
!   in many routines.  Kept in a separate section with references.
!
!   CHANGE LOG:
!     Initial creation consolidating CODATA-1963 literals scattered
!     across ~100 locations in the codebase.  Values updated to
!     CODATA 2018 (Tiesinga et al. 2021, Rev. Mod. Phys. 93, 025010)
!     and then to CODATA 2022 (Mohr et al. 2024, Rev. Mod. Phys. 97,
!     025002).  Differences between 2018 and 2022 are all <10 ppm.
!
!   Previously used (CODATA 1963):
!     k_B = 1.38054E-16,  h = 6.6256E-27,  m_u = 1.660E-24
!=========================================================================

MODULE mod_constants

  IMPLICIT NONE

  ! =====================================================================
  !  FUNDAMENTAL CONSTANTS  (CGS-Gaussian units)
  !  Source: CODATA 2022 (Mohr et al. 2024, RMP 97, 025002).
  !  Values marked (exact) are exact under the 2019 SI redefinition.
  !  Precision is written to the last digit REAL(8) can store.
  ! =====================================================================

  REAL(8), PARAMETER :: PI        = 3.141592653589793D0
  REAL(8), PARAMETER :: FOURPI    = 4.0D0 * PI
  REAL(8), PARAMETER :: SQRTPI    = 1.7724538509055159D0
  REAL(8), PARAMETER :: INVSQRTPI = 0.5641895835477563D0

  REAL(8), PARAMETER :: KBOL    = 1.380649D-16          ! Boltzmann constant [erg/K] (exact)
  REAL(8), PARAMETER :: HPLANCK = 6.62607015D-27        ! Planck constant [erg s]    (exact)
  REAL(8), PARAMETER :: CLIGHT  = 2.99792458D10         ! speed of light [cm/s]      (exact)

  ! Speed of light in other wavelength-per-time units
  REAL(8), PARAMETER :: CLIGHT_NMS  = 2.99792458D17     ! [nm/s]   freq = CLIGHT_NMS  / lambda_nm
  REAL(8), PARAMETER :: CLIGHT_ANGS = 2.99792458D18     ! [Å/s]    freq = CLIGHT_ANGS / lambda_Å
  REAL(8), PARAMETER :: CLIGHT_KMS  = 2.99792458D5      ! [km/s]

  REAL(8), PARAMETER :: AMU     = 1.66053906892D-24     ! atomic mass unit [g]
  REAL(8), PARAMETER :: EMASS   = 9.1093837139D-28      ! electron mass [g]
  REAL(8), PARAMETER :: ECHARGE = 4.80320471D-10        ! elementary charge [esu]

  REAL(8), PARAMETER :: SIGMA_SB      = 5.670374419D-5   ! Stefan-Boltzmann [erg/cm²/s/K⁴]
  REAL(8), PARAMETER :: SIGMA_THOMSON = 6.6524587051D-25 ! Thomson cross section [cm²]
  REAL(8), PARAMETER :: EV_ERG        = 1.602176634D-12  ! 1 eV in erg (exact)

  ! =====================================================================
  !  DERIVED CONSTANTS
  !  Pre-computed combinations of the fundamentals, used throughout the code.
  ! =====================================================================

  REAL(8), PARAMETER :: HCK     = HPLANCK * CLIGHT / KBOL     ! hc/k  [cm K]  exp(-E_cm * HCK / T)
  REAL(8), PARAMETER :: HOVERK  = HPLANCK / KBOL              ! h/k   [s K]   HKT = HOVERK / T
  REAL(8), PARAMETER :: KBOL_EV = KBOL / EV_ERG               ! k/eV  [1/K]   TKEV = T * KBOL_EV
  REAL(8), PARAMETER :: LN10    = 2.30258509299405D0          ! ln(10); log(x) = LN10 * log10(x)

  ! log₁₀(e) / k_B(eV) -- Boltzmann θ parameter:  θ = THETA_COEFF / T
  REAL(8), PARAMETER :: THETA_COEFF = 0.4342944819032518D0 / KBOL_EV

  REAL(8), PARAMETER :: FOURPI_OVER_C      = FOURPI / CLIGHT        ! 4π/c  [cm⁻¹ s]
  REAL(8), PARAMETER :: FOUR_SIGMA_OVER_PI = 4.0D0 * SIGMA_SB / PI  ! 4σ/π  (convection)

  ! Planck function prefactor: 2h/c² × (10¹⁵)³  [cgs], with ν in units of 10¹⁵ Hz.
  !   B_ν = BNU_PREFAC * (ν/10¹⁵)³ * exp(-hν/kT) / (1 - exp(-hν/kT))
  REAL(8), PARAMETER :: BNU_PREFAC = 2.0D0 * HPLANCK / CLIGHT**2 * 1.0D45

  ! Saha equation prefactor: 2 × (2π m_e k / h²)^(3/2)  [cm⁻³ K⁻³/²].
  ! The leading 2 is the electron spin weight g_e, so SAHA_PREFAC is the
  ! FULL prefactor; no extra factor of 2 is needed at call sites.
  !   n_e × n⁺/n₀ = SAHA_PREFAC × T^(3/2) × (U⁺/U₀) × exp(-χ/kT)
  REAL(8), PARAMETER :: SAHA_PREFAC = 2.0D0 &
    * (2.0D0 * PI * EMASS * KBOL / HPLANCK**2)**1.5D0

  ! Free-free opacity coefficient [cgs].  Mihalas (1978) eq. 4-64.
  !   κ_ff = COEFF_FF / ν³ × g_ff × n_e × n_ion / (ρ √T)
  !   COEFF_FF = 32π² e⁶ / (3√3 m_e c h) × √(2π / (3 m_e k))
  ! Kept at the ATLAS canonical value; e appears to the 6th power, so
  ! the CODATA 2022 shift in e (~4 ppm) would move this by ~25 ppm,
  ! well below the precision at which the constant is quoted.
  REAL(8), PARAMETER :: COEFF_FF = 3.69196D8

  ! =====================================================================
  !  SPECTROSCOPIC CONSTANTS  (measured, not fundamental)
  ! =====================================================================

  ! Rydberg constants [cm⁻¹].  CODATA 2022.
  !   R_∞ = 10 973 731.568 157 m⁻¹
  !   R_H  = R_∞ × 1/(1 + m_e/m_p)  with m_p/m_e = 1836.152 673 426
  !   R_He = R_∞ × M_He/(m_e + M_He) with ⁴He mass = 4.002 603 254 u
  REAL(8), PARAMETER :: RYDBERG_INF = 109737.31568157D0
  REAL(8), PARAMETER :: RYDBERG_H   = 109677.583403D0
  REAL(8), PARAMETER :: RYDBERG_HE  = 109722.277609D0

  ! H I ionization limit [cm⁻¹]  (observed, includes quantum defect)
  REAL(8), PARAMETER :: ELIM_HI = 109678.764D0

  ! Hydrogen series ionization limits [cm⁻¹]: H_SERIES_LIMITS(n) is the
  ! energy from level n to the continuum.  For n=1..8 these are the
  ! observed values ELIM_HI - E_obs(n) (matching the hand-tabulated
  ! Kurucz values).  For n=9..15 they agree with R_H/n² to <1 part in 10⁵
  ! since Lamb shift and fine structure become negligible.
  REAL(8), PARAMETER :: H_SERIES_LIMITS(15) = [ &
    109678.764D0, 27419.659D0, 12186.462D0, 6854.871D0, 4387.113D0, &
      3046.604D0,  2238.320D0,  1713.713D0, 1354.044D0, 1096.776D0, &
       906.426D0,   761.650D0,   648.980D0,  559.579D0,  487.456D0 ]

  ! Hydrogen ionization frequency [Hz]
  REAL(8), PARAMETER :: FREQ_RYDH = RYDBERG_H * CLIGHT

  ! Rydberg wavelength for hydrogen [Å] (= 1e8 / R_H_cm⁻¹).  Used to
  ! compute hydrogen transition wavelengths: λ = RYD_ANG × (n·m)² / (m²-n²).
  REAL(8), PARAMETER :: RYD_ANG   = 1.0D8 / RYDBERG_H

  ! H⁻ electron affinity [eV]
  !   Lykke, Murray & Lineberger (1991, PRA 43, 6104):
  !   E_A(H⁻) = 6082.99 ± 0.15 cm⁻¹ = 0.75419(2) eV
  REAL(8), PARAMETER :: HMINUS_EA = 0.75419D0


  ! =====================================================================
  !  HOLTSMARK MICROFIELD DISTRIBUTION TABLE
  !
  !  Q(β) = Prob(F < β F₀) for the Holtsmark microfield in a plasma,
  !  used in the Hummer & Mihalas (1988) occupation-probability formalism:
  !  the fraction of hydrogen level n that survives dissolution is
  !  w_n = Q(β_crit) with β_crit = BETA_COEFF_HM88 / (n⁵ × N_e^(2/3)).
  !
  !  Table: 150 points on a uniform log₁₀(β) grid from 0.01 to 50.
  !
  !  References:
  !    Holtsmark 1919, Ann. Phys. 363, 577
  !    Hummer & Mihalas 1988, ApJ 331, 794
  !    Nayfonov et al. 1999, ApJ 526, 451
  ! =====================================================================
  INTEGER, PARAMETER :: NQ_HOLTSMARK  = 150
  REAL(8), PARAMETER :: LOG_BETA_MIN  = -2.0000000000000000D+00
  REAL(8), PARAMETER :: LOG_BETA_MAX  =  1.6989700043360187D+00
  REAL(8), PARAMETER :: LOG_BETA_STEP =  2.4825302042523617D-02

  REAL(8), PARAMETER :: Q_HOLTSMARK(150) = [ &
    1.414671303101461D-07, 1.679306576215498D-07, 1.993444994252022D-07, 2.366346434794731D-07, 2.809002805410456D-07, &
    3.334461984807552D-07, 3.958212340900693D-07, 4.698639150426803D-07, 5.577566360795190D-07, 6.620899645927724D-07, &
    7.859389687582573D-07, 9.329538149374408D-07, 1.107467300593689D-06, 1.314622486716282D-06, 1.560524184270007D-06, &
    1.852418749731333D-06, 2.198907475766682D-06, 2.610199848759794D-06, 3.098414113872032D-06, 3.677933974566057D-06, &
    4.365831897218437D-06, 5.182371440131865D-06, 6.151603336163127D-06, 7.302072795786128D-06, 8.667658741262937D-06, &
    1.028856952546623D-05, 1.221252424026525D-05, 1.449615410836965D-05, 1.720666483126171D-05, 2.042380831343906D-05, &
    2.424222111026057D-05, 2.877419750061379D-05, 3.415297755655442D-05, 4.053664530978627D-05, 4.811274949652059D-05, &
    5.710377986116004D-05, 6.777365615449882D-05, 8.043541539930427D-05, 9.546031643883215D-05, 1.132886200658068D-04, &
    1.344423491071419D-04, 1.595403868046069D-04, 1.893163349213092D-04, 2.246396266107989D-04, 2.665404747620448D-04, &
    3.162393359894821D-04, 3.751816855286531D-04, 4.450790309992992D-04, 5.279572453545271D-04, 6.262134733849464D-04, &
    7.426830638017916D-04, 8.807182017925553D-04, 1.044280166083684D-03, 1.238047410111971D-03, 1.467541967671182D-03, &
    1.739277006102283D-03, 2.060928688577077D-03, 2.441535851096500D-03, 2.891731333794821D-03, 3.424009106941961D-03, &
    4.053031566870220D-03, 4.795981500493777D-03, 5.672963167595142D-03, 6.707456645870922D-03, 7.926828918285955D-03, &
    9.362904019250405D-03, 1.105259272456321D-02, 1.303857956034104D-02, 1.537006106786488D-02, 1.810352400480589D-02, &
    2.130354516765621D-02, 2.504358545069891D-02, 2.940673929734126D-02, 3.448638660516894D-02, 4.038667732307838D-02, &
    4.722275959985061D-02, 5.512064100258095D-02, 6.421655023415057D-02, 7.465564600958080D-02, 8.658990347532045D-02, &
    1.001750012765402D-01, 1.155660400501043D-01, 1.329119530555128D-01, 1.523485300474136D-01, 1.739900743830450D-01, &
    1.979198568471881D-01, 2.241797192777529D-01, 2.527594102861463D-01, 2.835864859104609D-01, 3.165178467055500D-01, &
    3.513341604156768D-01, 3.877384740815520D-01, 4.253601858765986D-01, 4.637651729352414D-01, 5.024722401913855D-01, &
    5.409752086892451D-01, 5.787690154213567D-01, 6.153773349921458D-01, 6.503786770405923D-01, 6.834278586173648D-01, &
    7.142702953716408D-01, 7.427476403869326D-01, 7.687946986070335D-01, 7.924289205182540D-01, 8.137347922132799D-01, &
    8.328458692792989D-01, 8.499270198328847D-01, 8.651588016097912D-01, 8.787250560031417D-01, 8.908040093261350D-01, &
    9.015625910920573D-01, 9.111533622086631D-01, 9.197133565854327D-01, 9.273641961484557D-01, 9.342129682962188D-01, &
    9.403534976718694D-01, 9.458677646981941D-01, 9.508273225009358D-01, 9.552946289851881D-01, 9.593242514728686D-01, &
    9.629639386838167D-01, 9.662555465715050D-01, 9.692358425479507D-01, 9.719372021216177D-01, 9.743882111568676D-01, &
    9.766141278943968D-01, 9.786373698849190D-01, 9.804778420731699D-01, 9.821532811516551D-01, 9.836795130160503D-01, &
    9.850706740272015D-01, 9.863394492399232D-01, 9.874972091274948D-01, 9.885541671623754D-01, 9.895195463975206D-01, &
    9.904015671946884D-01, 9.912077851485863D-01, 9.919449707318418D-01, 9.926192297837247D-01, 9.932361028219120D-01, &
    9.938006248722731D-01, 9.943173688128598D-01, 9.947904964619882D-01, 9.952237635494758D-01, 9.956206036156809D-01, &
    9.959840685797559D-01, 9.963172215001466D-01, 9.966224531295311D-01, 9.969021568157113D-01, 9.971585495409699D-01 ]

  ! β_crit(n, N_e) = BETA_COEFF_HM88 / (n⁵ × N_e^(2/3))
  ! Derived from Inglis-Teller critical field and Holtsmark normal field.
  REAL(8), PARAMETER :: BETA_COEFF_HM88 = 8.798905203085208D+14

END MODULE mod_constants


!=========================================================================
! MODULE mod_partition_functions
!
! Atomic partition functions from Barklem & Collet (2016), A&A 588, A96.
!
! Provides U(T) = sum_i g_i exp(-E_i / kT) evaluated by direct NIST-level
! summation, for Z = 1..92 and ionization stages I, II, III (284 species
! total), plus seven negative ions (H-, C-, O-, F-, Si-, S-, Cl-).
!
! These are UNPERTURBED (isolated-atom) partition functions.  They do not
! include pressure-lowering effects; those remain the province of the
! PFIRON table (for Z=20..28), which is treated via a multiplicative
! hybrid in PFIRON() itself (see that routine for details).
!
! Intended usage pattern:
!   1. At program startup, call set_bc_data_dir('/path/to/data/')
!   2. Query with U_BC(Z, ION, T) -- returns negative sentinel if the
!      requested (Z, ION) is outside B&C's coverage.
!
! On first call to any accessor, partfn_bc2016.dat is read lazily from
! the configured data directory.  Subsequent calls hit the cached table.
!
! Interpolation is log-linear in T (equivalent to piecewise power-law in
! U).  A bisection with a persistent bracket-hint makes lookups near-O(1)
! after warmup.  Below T_GRID(1) = 1e-5 K (unreachable in practice), the
! value at the lowest grid point is returned.  ABOVE T_GRID(42) = 10,000 K,
! lookups return the out-of-coverage sentinel -- B&C is NEVER extrapolated.
! Callers are expected to have their own high-T fallback (e.g. Kurucz's
! hand-coded partition function expressions, invoked by the PFGROUND_HYBRID
! dispatcher).  This policy keeps the model smooth and predictable in
! hot-star atmospheres, where log-linear extrapolation from the 5600/10000 K
! endpoints would give badly wrong partition functions.
!
! Sentinel convention: out-of-coverage accesses return -1.0d0.  Callers
! branch on U < 0 to fall back to whatever previous scheme they used.
! "Out of coverage" means any of: Z outside [1, 92], ION outside [1, 3],
! (Z, ION) not tabulated by B&C, OR T > 10,000 K.
!
! Data file format (partfn_bc2016.dat): upstream B&C Table 8, renamed
! from the published "table8_vNov2022.dat" to make its provenance clear
! in this codebase.
!=========================================================================

MODULE mod_partition_functions

  USE mod_parameters
  IMPLICIT NONE
  SAVE
  PRIVATE

  ! --- Public API ---
  PUBLIC :: U_BC, U_BC_NEG, BC_HAS, BC_HAS_NEG
  PUBLIC :: set_bc_data_dir, init_bc_partition_functions

  ! --- Table dimensions (upstream file is fixed-format) ---
  INTEGER, PARAMETER :: NT_BC     = 42    ! number of temperature grid points
  INTEGER, PARAMETER :: NZ_MAX    = 92    ! highest Z covered (U)
  INTEGER, PARAMETER :: NION_MAX  = 3     ! stages I, II, III
  INTEGER, PARAMETER :: NNEG_MAX  = 7     ! H-, C-, O-, F-, Si-, S-, Cl-

  ! --- Cached data (loaded once) ---
  REAL(8) :: T_GRID(NT_BC)                        ! temperature grid [K]
  REAL(8) :: LOGT_GRID(NT_BC)                     ! ln(T) for interpolation
  REAL(8) :: U_TAB(NT_BC, NZ_MAX, NION_MAX)       ! U(T, Z, ION)
  REAL(8) :: LOGU_TAB(NT_BC, NZ_MAX, NION_MAX)    ! ln(U) precomputed
  LOGICAL :: HAS_TAB(NZ_MAX, NION_MAX) = .FALSE. ! coverage flag

  REAL(8) :: U_NEG_TAB(NT_BC, NNEG_MAX)           ! U(T) for negative ions
  REAL(8) :: LOGU_NEG_TAB(NT_BC, NNEG_MAX)        ! ln(U) for negative ions
  INTEGER :: NEG_Z(NNEG_MAX) = 0                 ! parent atom Z per slot
  LOGICAL :: HAS_NEG(NZ_MAX) = .FALSE.           ! coverage flag by Z
  INTEGER :: NEG_SLOT(NZ_MAX) = 0                ! Z -> slot index

  LOGICAL :: INITIALIZED = .FALSE.

  ! --- Data directory (set by the caller before first use) ---
  CHARACTER(len=512) :: BC_DATADIR = './'

CONTAINS

!-----------------------------------------------------------------------
! set_bc_data_dir(dir): stash the path where partfn_bc2016.dat lives.
! Should be called once at startup, before the first partition-function
! lookup.  If not called, defaults to './'.
!-----------------------------------------------------------------------
SUBROUTINE set_bc_data_dir(dir)
  CHARACTER(len=*), INTENT(IN) :: dir
  BC_DATADIR = dir
END SUBROUTINE set_bc_data_dir

!-----------------------------------------------------------------------
! Lazy initialization: read the B&C table from DATADIR/partfn_bc2016.dat
! on first need.  Parses comments, the title line, the species count,
! the 42-point T grid, and 284 species records.
!-----------------------------------------------------------------------
SUBROUTINE init_bc_partition_functions()
  INTEGER            :: LUN, IOS, IT, IZ, ION, NEG_IDX
  INTEGER            :: N_SPECIES, N_READ
  CHARACTER(len=4096) :: LINE
  CHARACTER(len=16)   :: LABEL
  REAL(8)             :: VALS(NT_BC)

  IF (INITIALIZED) RETURN

  LUN = 89
  OPEN(LUN, FILE=trim(BC_DATADIR)//'partfn_bc2016.dat', &
       STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS .NE. 0) THEN
    WRITE(6,'(A,A)') ' INIT_BC_PARTITION_FUNCTIONS ERROR: cannot open ', &
         trim(BC_DATADIR)//'partfn_bc2016.dat'
    STOP 1
  END IF

  ! --- Skip header comments and blank lines to find the title line ---
  DO
    READ(LUN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A)') ' INIT_BC_PARTITION_FUNCTIONS ERROR: no title line found'
      STOP 1
    END IF
    LINE = adjustl(LINE)
    IF (len_trim(LINE) .EQ. 0) CYCLE
    IF (LINE(1:1) .EQ. '#')   CYCLE
    ! Title line is cosmetic; we just discard LINE here.
    EXIT
  END DO

  ! --- Read species count (next non-comment, non-blank line) ---
  N_SPECIES = 0
  DO
    READ(LUN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A)') ' INIT_BC_PARTITION_FUNCTIONS ERROR: missing species count'
      STOP 1
    END IF
    LINE = adjustl(LINE)
    IF (len_trim(LINE) .EQ. 0) CYCLE
    IF (LINE(1:1) .EQ. '#')   CYCLE
    READ(LINE, *, IOSTAT=IOS) N_SPECIES
    IF (IOS .NE. 0 .OR. N_SPECIES .LE. 0) THEN
      WRITE(6,'(A,A)') ' INIT_BC_PARTITION_FUNCTIONS ERROR: bad species count line: ', trim(LINE)
      STOP 1
    END IF
    EXIT
  END DO

  ! --- Read T grid line: first two tokens are "T" and "[K]" ---
  DO
    READ(LUN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A)') ' INIT_BC_PARTITION_FUNCTIONS ERROR: missing T grid'
      STOP 1
    END IF
    LINE = adjustl(LINE)
    IF (len_trim(LINE) .EQ. 0) CYCLE
    IF (LINE(1:1) .EQ. '#')   CYCLE
    CALL parse_t_grid(LINE, T_GRID)
    EXIT
  END DO
  DO IT = 1, NT_BC
    LOGT_GRID(IT) = log(T_GRID(IT))
  END DO

  ! --- Now read species records.  Each non-blank, non-comment line is a
  !     record: <label> <v1> <v2> ... <v42>.  Label is e.g. "Fe_I", "Cr_II",
  !     "H-", "D_I".  Deuterium is ignored (not distinguished from H).  ---
  NEG_IDX = 0
  N_READ  = 0
  DO
    READ(LUN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) EXIT
    LINE = adjustl(LINE)
    IF (len_trim(LINE) .EQ. 0) CYCLE
    IF (LINE(1:1) .EQ. '#')   CYCLE
    CALL parse_species_record(LINE, LABEL, VALS, IOS)
    IF (IOS .NE. 0) CYCLE
    N_READ = N_READ + 1
    ! Classify the label and stash the data.
    CALL classify_label(LABEL, IZ, ION)
    IF (IZ .LT. 0) CYCLE                   ! deuterium or unrecognized -- skip
    IF (ION .GT. 0) THEN
      ! Positive ion or neutral
      IF (IZ .LE. NZ_MAX .AND. ION .LE. NION_MAX) THEN
        U_TAB(:, IZ, ION) = VALS(:)
        ! Precompute ln(U) for speed.  B&C values are strictly positive
        ! for every tabulated species at every T, but guard anyway.
        DO IT = 1, NT_BC
          IF (VALS(IT) .GT. 0.0d0) THEN
            LOGU_TAB(IT, IZ, ION) = log(VALS(IT))
          ELSE
            LOGU_TAB(IT, IZ, ION) = -huge(0.0d0)
          END IF
        END DO
        HAS_TAB(IZ, ION) = .TRUE.
      END IF
    ELSE
      ! Negative ion
      NEG_IDX = NEG_IDX + 1
      IF (NEG_IDX .LE. NNEG_MAX) THEN
        U_NEG_TAB(:, NEG_IDX) = VALS(:)
        DO IT = 1, NT_BC
          IF (VALS(IT) .GT. 0.0d0) THEN
            LOGU_NEG_TAB(IT, NEG_IDX) = log(VALS(IT))
          ELSE
            LOGU_NEG_TAB(IT, NEG_IDX) = -huge(0.0d0)
          END IF
        END DO
        NEG_Z(NEG_IDX) = IZ
        IF (IZ .GE. 1 .AND. IZ .LE. NZ_MAX) THEN
          HAS_NEG(IZ) = .TRUE.
          NEG_SLOT(IZ) = NEG_IDX
        END IF
      END IF
    END IF
  END DO
  CLOSE(LUN)

  IF (N_READ .NE. N_SPECIES) THEN
    WRITE(6,'(A,I5,A,I5)') ' INIT_BC_PARTITION_FUNCTIONS WARNING: read ', &
         N_READ, ' species but file header claimed ', N_SPECIES
  END IF

  INITIALIZED = .TRUE.
END SUBROUTINE init_bc_partition_functions

!-----------------------------------------------------------------------
! Parse the temperature grid line.  Expected form:
!   "  T [K]   1.00000e-05   1.00000e-04   ... 1.00000e+04"
! That is, two leading tokens "T" and "[K]" followed by NT_BC values.
!-----------------------------------------------------------------------
SUBROUTINE parse_t_grid(LINE, TGRID)
  CHARACTER(len=*), INTENT(IN)  :: LINE
  REAL(8),           INTENT(OUT) :: TGRID(NT_BC)

  CHARACTER(len=len(LINE)) :: BUF
  INTEGER :: I, IOS, POS

  BUF = LINE
  ! Strip leading "T" and "[K]" tokens, whatever the exact spacing.
  ! Locate the "]" of "[K]" and start reading numbers after it.
  POS = index(BUF, ']')
  IF (POS .GT. 0) BUF = BUF(POS+1:)
  READ(BUF, *, IOSTAT=IOS) (TGRID(I), I=1,NT_BC)
  IF (IOS .NE. 0) THEN
    WRITE(6,'(A,A)') ' PARSE_T_GRID ERROR on line: ', trim(LINE)
    STOP 1
  END IF
END SUBROUTINE parse_t_grid

!-----------------------------------------------------------------------
! Parse a species record: "<label> <v1> <v2> ... <v42>".
! Returns IOS = 0 on success.
!-----------------------------------------------------------------------
SUBROUTINE parse_species_record(LINE, LABEL, VALS, IOS)
  CHARACTER(len=*),  INTENT(IN)  :: LINE
  CHARACTER(len=*),  INTENT(OUT) :: LABEL
  REAL(8),            INTENT(OUT) :: VALS(NT_BC)
  INTEGER,           INTENT(OUT) :: IOS

  CHARACTER(len=len(LINE)) :: BUF
  INTEGER :: I, POS

  BUF = adjustl(LINE)
  ! Label is the first whitespace-delimited token
  POS = scan(BUF, ' '//char(9))
  IF (POS .LE. 0) THEN
    IOS = 1
    RETURN
  END IF
  LABEL = BUF(1:POS-1)
  BUF   = adjustl(BUF(POS+1:))
  READ(BUF, *, IOSTAT=IOS) (VALS(I), I=1,NT_BC)
END SUBROUTINE parse_species_record

!-----------------------------------------------------------------------
! Classify a species label into (IZ, ION) or negative-ion form.
! Recognized label formats:
!   "H_I", "He_II", "Fe_III", ...  -> IZ > 0, ION in {1,2,3}
!   "D_I", "D_II"                   -> IZ = -1 (skip)
!   "H-", "C-", "O-", ...           -> IZ > 0, ION = 0 (negative ion)
!   anything else                   -> IZ = -2 (skip, unrecognized)
!-----------------------------------------------------------------------
SUBROUTINE classify_label(LABEL, IZ, ION)
  CHARACTER(len=*), INTENT(IN)  :: LABEL
  INTEGER,          INTENT(OUT) :: IZ
  INTEGER,          INTENT(OUT) :: ION

  INTEGER :: POS_U, POS_M
  CHARACTER(len=8) :: SYM, ROMAN

  POS_U = index(LABEL, '_')
  POS_M = index(LABEL, '-')

  IF (POS_U .GT. 0) THEN
    ! Positive ion or neutral: "<sym>_<roman>"
    SYM   = LABEL(1:POS_U-1)
    ROMAN = LABEL(POS_U+1:)
    IF (trim(SYM) .EQ. 'D') THEN
      IZ = -1   ! deuterium: skip
      ION = 0
      RETURN
    END IF
    IZ  = element_Z(trim(SYM))
    ION = roman_to_int(trim(ROMAN))
    IF (IZ .LE. 0 .OR. ION .LE. 0) THEN
      IZ = -2
      ION = 0
    END IF
  ELSE IF (POS_M .GT. 0 .AND. POS_M .EQ. len_trim(LABEL)) THEN
    ! Negative ion: "<sym>-"
    SYM = LABEL(1:POS_M-1)
    IZ  = element_Z(trim(SYM))
    ION = 0     ! negative-ion marker
    IF (IZ .LE. 0) THEN
      IZ = -2
    END IF
  ELSE
    IZ  = -2
    ION = 0
  END IF
END SUBROUTINE classify_label

!-----------------------------------------------------------------------
! Map a chemical symbol to atomic number Z.  Case-sensitive (B&C uses
! standard mixed-case symbols: H, He, Li, ... U).
!-----------------------------------------------------------------------
FUNCTION element_Z(SYM) RESULT(IZ)
  CHARACTER(len=*), INTENT(IN) :: SYM
  INTEGER :: IZ

  CHARACTER(len=3), PARAMETER :: SYMBOLS(92) = [ &
    'H  ', 'He ', 'Li ', 'Be ', 'B  ', 'C  ', 'N  ', 'O  ', 'F  ', 'Ne ', &
    'Na ', 'Mg ', 'Al ', 'Si ', 'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ', &
    'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ', &
    'Ga ', 'Ge ', 'As ', 'Se ', 'Br ', 'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ', &
    'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ', &
    'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ', 'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ', &
    'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ', &
    'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ', &
    'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ', 'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', &
    'Pa ', 'U  ' ]
  INTEGER :: I

  IZ = 0
  DO I = 1, 92
    IF (trim(SYMBOLS(I)) .EQ. SYM) THEN
      IZ = I
      RETURN
    END IF
  END DO
END FUNCTION element_Z

!-----------------------------------------------------------------------
! Convert a Roman numeral string "I"/"II"/"III" to integer 1/2/3.
! Returns 0 on unrecognized input.
!-----------------------------------------------------------------------
FUNCTION roman_to_int(ROMAN) RESULT(ION)
  CHARACTER(len=*), INTENT(IN) :: ROMAN
  INTEGER :: ION
  SELECT CASE (trim(ROMAN))
  CASE ('I');   ION = 1
  CASE ('II');  ION = 2
  CASE ('III'); ION = 3
  CASE DEFAULT; ION = 0
  END SELECT
END FUNCTION roman_to_int

!-----------------------------------------------------------------------
! U_BC(IZ, ION, T): partition function for atom Z in ionization
! stage ION (1=neutral, 2=singly ionized, 3=doubly ionized) at
! temperature T.  Returns -1.0d0 if (IZ, ION) is not in B&C coverage.
!-----------------------------------------------------------------------
FUNCTION U_BC(IZ, ION, T) RESULT(U)
  INTEGER, INTENT(IN) :: IZ, ION
  REAL(8),  INTENT(IN) :: T
  REAL(8) :: U

  IF (.NOT. INITIALIZED) CALL init_bc_partition_functions()

  U = -1.0d0
  IF (IZ .LT. 1 .OR. IZ .GT. NZ_MAX) RETURN
  IF (ION .LT. 1 .OR. ION .GT. NION_MAX) RETURN
  IF (.NOT. HAS_TAB(IZ, ION)) RETURN
  U = interp_log_linear_cached(T, U_TAB(:, IZ, ION), LOGU_TAB(:, IZ, ION))
END FUNCTION U_BC

!-----------------------------------------------------------------------
! U_BC_NEG(IZ, T): negative-ion partition function for the negative
! ion of atom Z (e.g., H- for IZ=1).  Returns -1.0d0 if not tabulated.
!-----------------------------------------------------------------------
FUNCTION U_BC_NEG(IZ, T) RESULT(U)
  INTEGER, INTENT(IN) :: IZ
  REAL(8),  INTENT(IN) :: T
  REAL(8) :: U
  INTEGER :: SLOT

  IF (.NOT. INITIALIZED) CALL init_bc_partition_functions()

  U = -1.0d0
  IF (IZ .LT. 1 .OR. IZ .GT. NZ_MAX) RETURN
  IF (.NOT. HAS_NEG(IZ)) RETURN
  SLOT = NEG_SLOT(IZ)
  IF (SLOT .LT. 1 .OR. SLOT .GT. NNEG_MAX) RETURN
  U = interp_log_linear_cached(T, U_NEG_TAB(:, SLOT), LOGU_NEG_TAB(:, SLOT))
END FUNCTION U_BC_NEG

!-----------------------------------------------------------------------
! Coverage predicates.
!-----------------------------------------------------------------------
FUNCTION BC_HAS(IZ, ION) RESULT(YES)
  INTEGER, INTENT(IN) :: IZ, ION
  LOGICAL :: YES
  IF (.NOT. INITIALIZED) CALL init_bc_partition_functions()
  YES = .FALSE.
  IF (IZ .LT. 1 .OR. IZ .GT. NZ_MAX) RETURN
  IF (ION .LT. 1 .OR. ION .GT. NION_MAX) RETURN
  YES = HAS_TAB(IZ, ION)
END FUNCTION BC_HAS

FUNCTION BC_HAS_NEG(IZ) RESULT(YES)
  INTEGER, INTENT(IN) :: IZ
  LOGICAL :: YES
  IF (.NOT. INITIALIZED) CALL init_bc_partition_functions()
  YES = .FALSE.
  IF (IZ .LT. 1 .OR. IZ .GT. NZ_MAX) RETURN
  YES = HAS_NEG(IZ)
END FUNCTION BC_HAS_NEG

!-----------------------------------------------------------------------
! Log-linear interpolation in T, using precomputed ln(U) values.
!
! Performance notes:
!   (1) ln(U) is cached at load time -- no log() calls at runtime.
!   (2) Only one exp() per call (to convert back from ln space).
!   (3) The bracket is located via a persistent index hint LAST_LO:
!       subsequent calls with similar T re-use the bracket without
!       searching.  PFSAHA calls this thousands of times with slowly
!       varying T (different depth points in the atmosphere), so the
!       hint hits frequently.  On miss we fall back to bisection over
!       the 42-point grid (~6 iterations).
!
! Arguments:
!   T        -- temperature [K]
!   U_ARR    -- raw partition function values on the 42-point grid.
!               Used only for the low-T clamp and the fallback linear-
!               in-U path when a grid endpoint is zero (never occurs
!               in practice for B&C data).
!   LOGU_ARR -- precomputed ln(U_ARR), cached at init time.
!-----------------------------------------------------------------------
FUNCTION interp_log_linear_cached(T, U_ARR, LOGU_ARR) RESULT(U)
  REAL(8), INTENT(IN) :: T
  REAL(8), INTENT(IN) :: U_ARR(NT_BC)
  REAL(8), INTENT(IN) :: LOGU_ARR(NT_BC)
  REAL(8) :: U

  INTEGER :: LO, HI, MID
  REAL(8)  :: LT, LT_LO, LT_HI, F, LU_LO, LU_HI

  ! Persistent hint: remember the bracket from the last call.
  INTEGER :: LAST_LO = 10
  SAVE    :: LAST_LO

  ! --- Clamp below the grid (harmless; T_GRID(1) = 1e-5 K is unreachable
  !     for any stellar atmosphere, but cool-atmosphere callers occasionally
  !     probe a few hundred K outer boundary layers). ---
  IF (T .LE. T_GRID(1)) THEN
    U = U_ARR(1)
    RETURN
  END IF

  ! --- Hard ceiling at the top of the grid: NO extrapolation.
  !     Return the out-of-coverage sentinel so the caller can fall back
  !     to whatever scheme it uses for high T (typically a hand-coded
  !     Kurucz expression).  This is important for ATLAS convergence:
  !     blind log-linear extrapolation above 10,000 K gives wildly wrong
  !     partition functions for hot-star work. ---
  IF (T .GE. T_GRID(NT_BC)) THEN
    U = -1.0d0
    RETURN
  END IF

  ! --- Try the hint first: does T lie in the last bracket? ---
  LO = LAST_LO
  IF (LO .LT. 1 .OR. LO .GE. NT_BC) LO = 10
  IF (T .GE. T_GRID(LO) .AND. T .LT. T_GRID(LO+1)) THEN
    HI = LO + 1
  ELSE IF (LO+1 .LT. NT_BC .AND. T .GE. T_GRID(LO+1) .AND. T .LT. T_GRID(LO+2)) THEN
    ! Common case: T moved up by one grid cell (next depth point)
    LO = LO + 1
    HI = LO + 1
    LAST_LO = LO
  ELSE IF (LO .GT. 1 .AND. T .GE. T_GRID(LO-1) .AND. T .LT. T_GRID(LO)) THEN
    ! T moved down by one grid cell
    LO = LO - 1
    HI = LO + 1
    LAST_LO = LO
  ELSE
    ! Miss: bisect
    LO = 1
    HI = NT_BC
    DO WHILE (HI - LO .GT. 1)
      MID = (LO + HI) / 2
      IF (T_GRID(MID) .GT. T) THEN
        HI = MID
      ELSE
        LO = MID
      END IF
    END DO
    LAST_LO = LO
  END IF

  ! --- Log-linear interpolation using cached LOGU_ARR ---
  LT    = log(T)
  LT_LO = LOGT_GRID(LO)
  LT_HI = LOGT_GRID(HI)
  F     = (LT - LT_LO) / (LT_HI - LT_LO)

  IF (U_ARR(LO) .GT. 0.0d0 .AND. U_ARR(HI) .GT. 0.0d0) THEN
    LU_LO = LOGU_ARR(LO)
    LU_HI = LOGU_ARR(HI)
    U = exp(LU_LO + F * (LU_HI - LU_LO))
  ELSE
    ! Fallback linear-in-U (never occurs for B&C data in practice)
    U = U_ARR(LO) + F * (U_ARR(HI) - U_ARR(LO))
  END IF
END FUNCTION interp_log_linear_cached

END MODULE mod_partition_functions

!=========================================================================
! mod_atlas_data: Shared arrays replacing all COMMON blocks
!   Original ATLAS12 used 58 COMMON blocks for shared state.
!   All are consolidated here as module variables with SAVE attribute.
!=========================================================================

MODULE mod_atlas_data

  USE mod_parameters
  USE mod_constants
  USE mod_partition_functions
  IMPLICIT NONE
  SAVE

  ! --- Rosseland optical depth scale ---
  ! COMMON /ABROSS/
  REAL(8)  :: ABROSS(kw), TAUROS(kw)

  ! --- Total absorption and scattering ---
  ! COMMON /ABTOT/
  REAL(8)  :: ABTOT(kw), ALPHA(kw)

  ! --- Cool-star continuous opacities (C I, Mg I, Al I, Si I, Fe I) ---
  ! COMMON /ACOOL/
  REAL(8)  :: AC1(kw), AMG1(kw), AAL1(kw), ASI1(kw), AFE1(kw)
  ! Per-species source functions: computed by C1OP/AL1OP/FE1OP/HE2OP
  ! but NOT consumed by KAPP, which follows the atlas12.for convention
  ! of weighting all metal/He continua by BNU.  Kept for potential
  ! future per-species source function treatment (e.g. NLTE departure
  ! coefficient hooks).
  REAL(8)  :: SAL1(kw), SFE1(kw), SHE2(kw), SC1(kw)

  ! --- Lukewarm opacities (C II, N I, O I, Mg II, Si II, Ca II) ---
  ! COMMON /ALUKE/
  REAL(8)  :: AC2(kw), AN1(kw), AO1(kw), AMG2(kw), ASI2(kw), ACA2(kw)

  ! --- Convection parameters ---
  ! COMMON /CONV/
  REAL(8)  :: DLTDLP(kw) = 0.0d0, HEATCP(kw) = 0.0d0, DLRDLT(kw) = 0.0d0, VELSND(kw) = 0.0d0
  REAL(8)  :: GRDADB(kw) = 0.0d0, HSCALE(kw) = 0.0d0, FLXCNV(kw) = 0.0d0, VCONV(kw) = 0.0d0
  REAL(8)  :: MIXLTH = 2.0d0, OVERWT = 0.0d0
  REAL(8)  :: FLXCNV0(kw), FLXCNV1(kw)
  INTEGER :: IFCONV = 1, NCONV = 30

  ! --- NLTE departure coefficients ---
  ! COMMON /DEPART/
  REAL(8)  :: BHYD(kw, 6) = 1.0d0, BMIN(kw) = 1.0d0
  INTEGER :: NLTEON = 0

  ! --- Electron density ---
  ! COMMON /EDENS/
  REAL(8)  :: EDENS(kw)
  INTEGER :: IFEDNS

  ! --- Element abundances, atomic masses, labels ---
  ! COMMON /ELEM/
  REAL(8)       :: YABUND(99)
  ! Default solar abundances: Anders & Grevesse (1989, Geochim. Cosmochim.
  ! Acta, 53, 197).  H and He are fractional number densities (H+He=1);
  ! Z >= 3 are log10(N_Z / N_total).  Elements with -20.00 have no
  ! astrophysically relevant abundance.
  REAL(8)       :: ABUND(99) = (/ &
  !   1 H        2 He
       0.911D0,  0.089D0, &
  !   3 Li       4 Be       5 B        6 C        7 N        8 O        9 F       10 Ne
     -10.88D0, -10.89D0,  -9.44D0,  -3.48D0,  -3.99D0,  -3.11D0,  -7.48D0,  -3.95D0, &
  !  11 Na      12 Mg      13 Al      14 Si      15 P       16 S       17 Cl      18 Ar
      -5.71D0,  -4.46D0,  -5.57D0,  -4.49D0,  -6.59D0,  -4.83D0,  -6.54D0,  -5.48D0, &
  !  19 K       20 Ca      21 Sc      22 Ti      23 V       24 Cr      25 Mn      26 Fe
      -6.92D0,  -5.68D0,  -8.94D0,  -7.05D0,  -8.04D0,  -6.37D0,  -6.65D0,  -4.37D0, &
  !  27 Co      28 Ni      29 Cu      30 Zn      31 Ga      32 Ge      33 As      34 Se
      -7.12D0,  -5.79D0,  -7.83D0,  -7.44D0,  -9.16D0,  -8.63D0,  -9.67D0,  -8.69D0, &
  !  35 Br      36 Kr      37 Rb      38 Sr      39 Y       40 Zr      41 Nb      42 Mo
      -9.41D0,  -8.81D0,  -9.44D0,  -9.14D0,  -9.80D0,  -9.44D0, -10.62D0, -10.12D0, &
  !  43 Tc      44 Ru      45 Rh      46 Pd      47 Ag      48 Cd      49 In      50 Sn
     -20.00D0, -10.20D0, -10.92D0, -10.35D0, -11.10D0, -10.18D0, -10.38D0, -10.04D0, &
  !  51 Sb      52 Te      53 I       54 Xe      55 Cs      56 Ba      57 La      58 Ce
     -11.04D0,  -9.80D0, -10.53D0,  -9.81D0, -10.92D0,  -9.91D0, -10.82D0, -10.49D0, &
  !  59 Pr      60 Nd      61 Pm      62 Sm      63 Eu      64 Gd      65 Tb      66 Dy
     -11.33D0, -10.54D0, -20.00D0, -11.04D0, -11.53D0, -10.92D0, -12.14D0, -10.94D0, &
  !  67 Ho      68 Er      69 Tm      70 Yb      71 Lu      72 Hf      73 Ta      74 W
     -11.78D0, -11.11D0, -12.04D0, -10.96D0, -11.28D0, -11.16D0, -11.91D0, -10.93D0, &
  !  75 Re      76 Os      77 Ir      78 Pt      79 Au      80 Hg      81 Tl      82 Pb
     -11.77D0, -10.59D0, -10.69D0, -10.24D0, -11.03D0, -10.95D0, -11.14D0, -10.19D0, &
  !  83 Bi      84 Po      85 At      86 Rn      87 Fr      88 Ra      89 Ac      90 Th
     -11.33D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0, -11.92D0, &
  !  91 Pa      92 U       93 Np      94 Pu      95 Am      96 Cm      97 Bk      98 Cf      99 Es
     -20.00D0, -12.51D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0, -20.00D0/)
  REAL(8)       :: ATMASS(99) = (/ 1.008D0, 4.003D0, &
   6.939D0,9.013D0,10.81D0,12.01D0,14.01D0,16.00D0,19.00D0,20.18D0,22.99D0,24.31D0, &
  26.98D0,28.09D0,30.98D0,32.07D0,35.45D0,39.95D0,39.10D0,40.08D0,44.96D0,47.90D0, &
  50.94D0,52.00D0,54.94D0,55.85D0,58.94D0,58.71D0,63.55D0,65.37D0,69.72D0,72.60D0, &
  74.92D0,78.96D0,79.91D0,83.80D0,85.48D0,87.63D0,88.91D0,91.22D0,92.91D0,95.95D0, &
  99.00D0,101.1D0,102.9D0,106.4D0,107.9D0,112.4D0,114.8D0,118.7D0,121.8D0,127.6D0, &
  126.9D0,131.3D0,132.9D0,137.4D0,138.9D0,140.1D0,140.9D0,144.3D0,147.0D0,150.4D0, &
  152.0D0,157.3D0,158.9D0,162.5D0,164.9D0,167.3D0,168.9D0,173.0D0,175.0D0,178.5D0, &
  181.0D0,183.9D0,186.3D0,190.2D0,192.2D0,195.1D0,197.0D0,200.6D0,204.4D0,207.2D0, &
  209.0D0,210.0D0,211.0D0,222.0D0,223.0D0,226.1D0,227.1D0,232.0D0,231.0D0,238.0D0, &
  237.0D0,244.0D0,243.0D0,247.0D0,247.0D0,251.0D0,254.0D0/)
  CHARACTER(2) :: ELEM(99) = (/ 'H ','He', &
  'Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg', &
  'Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti', &
  'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge', &
  'As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo', &
  'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te', &
  'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
  'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
  'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb', &
  'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ', &
  'NP','Pu','Am','Cm','Bk','Cf','Es'/)

  ! --- Radiative flux and flux derivatives ---
  ! COMMON /FLUX/
  REAL(8)  :: FLUX = 0.0d0, FLXERR(kw) = 0.0d0, FLXDRV(kw) = 0.0d0, FLXRAD(kw)

  ! --- Free-format input parser state ---
  ! COMMON /FREE/
  INTEGER :: NUMCOL, LETCOL, LAST, MORE, IFFAIL, MAXPOW

  ! --- Current frequency point ---
  ! COMMON /FREQ/
  REAL(8)  :: FREQ, FREQLG, EHVKT(kw), STIM(kw), BNU(kw), WAVE, WAVENO

  ! --- Frequency set for opacity distribution functions ---
  ! COMMON /FRESET/
  REAL(8)  :: WAVESET(30000), RCOSET(30000)
  INTEGER :: NULO, NUHI, NUMNU = 0, NUSTEP

  ! --- Hydrogen bound-free opacity tables ---
  ! COMMON /H1TAB/
  REAL(8)  :: H0TAB(2001), H1TAB(2001), H2TAB(2001)

  ! --- Geometric height scale ---
  ! COMMON /HEIGHT/
  REAL(8)  :: HEIGHT(kw)

  ! --- Control flags ---
  ! COMMON /IF/
  INTEGER :: IFCORR = 1, IFPRES = 1, IFSURF = 0, IFSCAT = 1, IFMOL = 1, IFREADLINES = 1

  ! ROSSTAB interpolation mode:
  !   1 = original bilinear (4-quadrant nearest neighbor)
  !   2 = Shepard (K-nearest, inverse-distance weighted, p=3)
  INTEGER :: IROSSTAB = 1

  ! --- Molecular equilibrium ---
  ! COMMON /IFEQUA/
  REAL(8)  :: XNMOLCODE(maxmol), EQUIL(6, maxmol)
  INTEGER :: IFEQUA(101), KCOMPS(maxloc), LOCJ(max1), IDEQUA(maxeq)
  INTEGER :: NEQUA, NEQUA1, NEQNEQ, NUMMOL, NLOC

  ! --- Population control flag ---
  ! COMMON /IFPOP/
  INTEGER :: IFPOP

  ! --- Line data integer header words ---
  ! COMMON /IIIIIII/ (renamed LINEREC) — packed binary line-data record
  INTEGER    :: IWL
  INTEGER(2)  :: IELION, IELO, IGFLOG, IGR, IGS, IGW

  ! --- Isotope fractional abundances ---
  ! COMMON /ISOTOPE/
  REAL(8)  :: ISOTOPE(10, 2, mion)

  ! --- Iteration control ---
  ! COMMON /ITER/
  INTEGER :: ITER, ifprnt(60) = 2, ifpnch(60) = 0, NUMITS = 0

  ! --- Title and header data ---
  ! COMMON /JUNK/
  CHARACTER(1) :: TITLE(74) = ' '
  REAL(8)  :: XSCALE = 1.0d0
  CHARACTER(4) :: WLTE = 'LTE '
  CHARACTER(256) :: DATADIR    ! Path to data files (from $ATLAS12)
  INTEGER :: INPUTDATA

  ! --- Flag set if SYNTHE is running ---
  INTEGER :: IFSYNTHE = 0

  ! --- Flag set to use TOPbase metal continuum opacities ---
  LOGICAL :: USE_TOPBASE_MBF = .TRUE.
  
  ! --- Flag set if SYNTHE should emit .mol/.linform output files ---
  ! Gated by the more_output CLI argument.  NMOLEC checks IFMOLOUT to
  ! decide whether to write its molecular-density table to unit 35.
  INTEGER :: IFMOLOUT = 0

  ! --- Radiative transfer J-coefficient matrix ---
  ! COMMON /MATXJ/
  REAL(8)  :: COEFJ(51, 51)

  ! --- Radiative transfer H-coefficient matrix ---
  ! COMMON /MATXH/
  REAL(8)  :: COEFH(51, 51)

  ! --- Angular quadrature ---
  ! COMMON /MUS/
  REAL(8)  :: ANGLE(20), SURFI(20) = 0.0d0
  INTEGER :: NMU = 1

  ! --- Individual continuous opacity sources ---
  ! COMMON /OPS/
  REAL(8)  :: AHYD(kw), AH2P(kw), AHMIN(kw), SIGH(kw)
  REAL(8)  :: AHE1(kw), AHE2(kw), AHEMIN(kw), SIGHE(kw)
  REAL(8)  :: ACOOL(kw), ALUKE(kw), AHOT(kw)
  REAL(8)  :: AH2COLL(kw)           ! H2-H2 collision-induced absorption
  REAL(8)  :: ACONT_METAL(kw)       ! consolidated metal+molecular continuum
  REAL(8)  :: SIGEL(kw), SIGH2(kw), AHLINE(kw), ALINES(kw), SIGLIN(kw)
  REAL(8)  :: AXLINE(kw), SIGXL(kw), SIGX(kw)
  REAL(8)  :: SHYD(kw), SHMIN(kw), SHLINE(kw), SXLINE(kw), SXCONT(kw)

  ! --- Total opacity ---
  ! COMMON /OPTOT/
  REAL(8)  :: ACONT(kw), SCONT(kw), ALINE(kw), SLINE(kw)
  REAL(8)  :: SIGMAC(kw), SIGMAL(kw)

  ! --- TOPbase / Allende-Prieto metal bound-free data ---
  !
  ! Each TOPbase level stores pre-smoothed photoionization cross section
  ! on a uniform 3% log-spaced energy grid (Allende Prieto et al. 2003,
  ! ApJS 147, 363).  Threshold energies are shifted at init time so the
  ! ground-state ionization potential matches NIST exactly; the same
  ! shift is applied to every level's photon-energy grid.
  !
  TYPE :: TOPBASE_LEVEL
    INTEGER :: ISLP                    ! term: 100*(2S+1) + 10*L + parity
    INTEGER :: ILV                     ! level index within term
    INTEGER :: G_STAT                  ! stat. weight (2S+1)*(2L+1)
    INTEGER :: NP                      ! number of (E, sigma) points
    REAL(8) :: E_BIND_RY               ! binding energy (negative, Ry)
    REAL(8) :: E_EXC_RY                ! excitation energy from ground (Ry)
    REAL(8) :: LOG_E_MIN               ! log10 of first photon-energy point
    REAL(8) :: LOG_E_STRIDE            ! log10 stride between points
    REAL(8), ALLOCATABLE :: E_PHOTON_RY(:)   ! photon energy grid (Ry)
    REAL(8), ALLOCATABLE :: SIGMA_MB(:)      ! cross section (Mb = 1e-18 cm^2)
  END TYPE TOPBASE_LEVEL

  TYPE :: TOPBASE_SPECIES
    CHARACTER(len=16) :: LABEL        ! e.g. "CaI", "FeII"
    INTEGER :: XNFP_INDEX             ! XNFP/XNF slot for this species
    INTEGER :: NLEVELS                ! # levels (after init-time truncation)
    REAL(8) :: E_IP_NIST_RY           ! NIST ionization potential (Ry)
    REAL(8) :: E_SHIFT_RY             ! applied E shift = IP_NIST - |E_ground_raw|
    REAL(8) :: LOWEST_THR_RY          ! min |E_BIND| across retained levels (Ry)
    TYPE(TOPBASE_LEVEL), ALLOCATABLE :: LEVELS(:)
  END TYPE TOPBASE_SPECIES

  INTEGER, PARAMETER :: N_TB_SPECIES_MAX = 30
  TYPE(TOPBASE_SPECIES) :: TB_SPECIES(N_TB_SPECIES_MAX)
  INTEGER :: N_TB_SPECIES = 0
  LOGICAL :: TB_INITIALIZED = .FALSE.

  ! Per-depth Boltzmann cache for MBF_TOPBASE.  Indexed by depth J and
  ! species S; each slot holds sum over levels of g_i * exp(-E_exc_i/kT_J)
  ! times the per-level sigma(nu).  Rebuilt when either T(J) or FREQ
  ! changes.  We use two cache levels: Boltzmann factors alone (cached
  ! until T(J) changes) plus sigma at the current FREQ (recomputed every
  ! frequency call).
  REAL(8), ALLOCATABLE :: TB_BOLTZ(:,:)    ! (level_flat_idx, J)
  REAL(8), ALLOCATABLE :: TB_T_LAST(:)     ! last T(J) at which TB_BOLTZ was built
  INTEGER, ALLOCATABLE :: TB_LEVEL_SPECIES(:)  ! flat -> species index
  INTEGER :: N_TB_LEVELS_TOTAL = 0

  ! Filtered hotop.dat data (excludes species covered by TOPbase)
  INTEGER :: N_HIGH_ION_EDGES = 0
  REAL(8), ALLOCATABLE :: HI_EDGE_DATA(:,:)   ! (NPAR=7, N_HIGH_ION_EDGES)
  LOGICAL :: HIGH_ION_INITIALIZED = .FALSE.

  ! --- Low-ionization iron (Fe I / Fe II) bound-free data ---
  ! Opacity Project / Iron Project R-matrix cross sections:
  !   Fe I  : Bautista (1997), A&AS 122, 167
  !   Fe II : Nahar & Pradhan (1994), J. Phys. B 27, 429
  ! Source data from TLUSTY model atoms (Hubeny & Lanz 2017), reformatted
  ! to op_fe1.dat and op_fe2.dat in the TLUSTY OP fit-point convention:
  !     x    = log10(nu / nu_threshold)
  !     sig  = 10^(logsig) * 1e-18 cm^2     (logsig = log10(sigma in Mb))
  ! Used by FELO_OPACITY when USE_TOPBASE_MBF = .TRUE.
  !
  ! Filtering note: hotop.dat contains NO iron edges at any ionization
  ! stage (verified by inventory -- only C/N/O/Ne, IDs 22-60), so no
  ! double-counting occurs with felo regardless of filter configuration.
  INTEGER, PARAMETER :: FELO_N_ION    = 2      ! 1 = Fe I, 2 = Fe II
  INTEGER, PARAMETER :: FELO_NLEV_MAX = 60     ! Fe I has 45, Fe II has 33
  INTEGER, PARAMETER :: FELO_NPTS_MAX = 200    ! Fe II max NF = 118
  INTEGER, PARAMETER :: FELO_XNFP_IDX(FELO_N_ION) = (/ 351, 352 /)
  INTEGER :: FELO_NLEV  (FELO_N_ION)                                 = 0
  REAL(8) :: FELO_NUTH  (FELO_NLEV_MAX, FELO_N_ION)                  = 0.0D0
  REAL(8) :: FELO_G     (FELO_NLEV_MAX, FELO_N_ION)                  = 0.0D0
  REAL(8) :: FELO_E     (FELO_NLEV_MAX, FELO_N_ION)                  = 0.0D0
  INTEGER :: FELO_IFANCY(FELO_NLEV_MAX, FELO_N_ION)                  = 0
  INTEGER :: FELO_NPTS  (FELO_NLEV_MAX, FELO_N_ION)                  = 0
  REAL(8) :: FELO_X     (FELO_NPTS_MAX, FELO_NLEV_MAX, FELO_N_ION)   = 0.0D0
  REAL(8) :: FELO_LOGSIG(FELO_NPTS_MAX, FELO_NLEV_MAX, FELO_N_ION)   = 0.0D0
  REAL(8) :: FELO_SEATON(4, FELO_NLEV_MAX, FELO_N_ION)               = 0.0D0
  LOGICAL :: FELO_INITIALIZED = .FALSE.

  ! --- Ionization potentials ---
  ! COMMON /POTION/
  REAL(8)  :: POTION(999), POTIONSUM(999)

  ! --- Total pressure ---
  ! COMMON /PTOTAL/
  REAL(8)  :: PTOTAL(kw)

  ! --- Output control ---
  ! COMMON /PUT/
  REAL(8)  :: PUT
  INTEGER :: IPUT

  ! --- Zero-point pressure and radiation field ---
  ! COMMON /PZERO/
  REAL(8)  :: PZERO, PCON, PRADK0, PTURB0
  REAL(8)  :: KNU(kw), PRADK(kw), RADEN(kw)

  ! --- Radiation pressure ---
  ! COMMON /RAD/
  REAL(8)  :: ACCRAD(kw) = 0.0d0, PRAD(kw) = 0.0d0

  ! --- Column mass depth scale ---
  ! COMMON /RHOX/
  REAL(8)  :: RHOX(kw)
  INTEGER :: NRHOX = 0

  ! --- Mean intensity diagnostic ---
  ! COMMON /RR/
  REAL(8)  :: RJMINSNU(kw), RDIAGJNU(kw)

  ! --- Thermodynamic state (pressure, electron density, etc.) ---
  ! COMMON /STATE/
  REAL(8)  :: P(kw), XNE(kw), XNATOM(kw), RHO(kw), CHARGESQ(kw)

  ! --- Depth integration step parameters ---
  ! COMMON /STEPLG/
  REAL(8)  :: STEPLG = 0.125d0, TAU1LG = -6.875d0
  INTEGER :: KRHOX = 0

  ! --- Continuous opacity table ---
  ! COMMON /TABCONT/
  REAL(8)  :: TABCONT(kw, 344)
  REAL(8)  :: WAVETAB(344)
  INTEGER :: IWAVETAB(344)

  ! --- Log lookup table (replaced array with function) ---
  ! COMMON /TABLOG/ removed — now computed inline via TABLOG()

  ! --- Optical depth, source function, flux moments ---
  ! COMMON /TAUSHJ/
  REAL(8)  :: TAUNU(kw), SNU(kw), HNU(kw), JNU(kw), JMINS(kw)

  ! --- Standard optical depth scale ---
  ! COMMON /TAUSTD/
  REAL(8)  :: TAUSTD(kw)

  ! --- Effective temperature and gravity ---
  ! COMMON /TEFF/
  REAL(8)  :: TEFF = 0.0d0, GRAV = 0.0d0, GLOG

  ! --- Temperature and derived quantities ---
  ! COMMON /TEMP/
  REAL(8)  :: T(kw), TKEV(kw), TK(kw), HKT(kw), HCKT(kw), TLOG(kw)
  INTEGER :: ITEMP

  ! --- Temperature smoothing ---
  ! COMMON /TSMOOTH/
  INTEGER :: J1SMOOTH = 0, J2SMOOTH = 0
  REAL(8)  :: WTJM1 = 0.3d0, WTJ = 0.4d0, WTJP1 = 0.3d0, TSMOOTH(kw)

  ! --- Turbulent pressure ---
  ! COMMON /TURBPR/
  REAL(8)  :: VTURB(kw) = 0.0d0, PTURB(kw) = 0.0d0, TRBFDG = 0.0d0, TRBCON = 0.0d0, TRBPOW = 0.0d0, TRBSND = 0.0d0
  INTEGER :: IFTURB = 0

  ! --- Wavelength grid control ---
  ! COMMON /WAVEY/
  REAL(8)  :: WBEGIN, DELTAW
  INTEGER :: IFWAVE = 0

  ! --- Current line data ---
  ! COMMON /WWWWWWW/ (renamed LINEPARAM) — packed line-parameter record
  REAL(8)  :: WLVAC
  INTEGER :: NELION
  REAL(8)  :: CGF, ELO, GAMMAR, GAMMAS, GAMMAW

  ! --- Depth-dependent abundances ---
  ! COMMON /XABUND/
  REAL(8)  :: XABUND(kw, 99), WTMOLE(kw), XRELATIVE(99) = 0.D0

  ! --- Isotope fractions and masses ---
  ! COMMON /XISO/
  REAL(8)  :: XISO(10, mion), AMASSISO(10, mion)

  ! --- Line opacity distribution ---
  ! COMMON /XLINES/
  REAL(8)  :: XLINES(kw, 30000)

  ! --- Ion number densities ---
  ! COMMON /XNF/
  REAL(8)  :: XNF(kw, mion), XNFP(kw, mion), XNH2(kw)

  ! --- Doppler-broadened line opacity ---
  ! COMMON /XNFDOP/
  REAL(8)  :: XNFDOP(kw, mion), DOPPLE(kw, mion)

  ! --- Molecular number densities ---
  ! COMMON /XNMOL/
  REAL(8)  :: XNMOL(kw, maxmol), XNFPMOL(kw, maxmol)

  ! Flag: set to 1 after MOLEC has read molecular data from INPUTDATA
  ! (or when NMOLEC has populated the arrays directly, skipping the read).
  INTEGER :: MOLEC_IREAD = 0

  ! --- Saved number densities for equilibrium ---
  ! COMMON /XNSAVE/
  REAL(8)  :: XNSAVE(kw, maxeq)

  ! --- In-memory line data storage (replaces fort.12 I/O) ---
  INTEGER(4), ALLOCATABLE :: LINEDATA(:,:)   ! (4, NLINES_STORED)
  INTEGER :: NLINES_STORED = 0

  ! --- Stehlé MMM hydrogen Stark broadening tables ---
  ! Preprocessed from Stehlé & Hutcheon (1999) and Stehlé & Fouquet (2010).
  ! Loaded by INIT_STARK_TABLES; used by STARK_MMM.
  INTEGER, PARAMETER :: NSTARK_SERIES = 4         ! Ly, Ba, Pa, Br
  INTEGER, PARAMETER :: NSTARK_DALPHA = 60        ! Δα grid points
  INTEGER, PARAMETER :: NSTARK_TEMPS  = 10        ! temperature grid points
  INTEGER, PARAMETER :: NSTARK_DENS_MAX = 20      ! max density grid points

  TYPE :: stark_series_t
    INTEGER :: n_lower                             ! lower quantum number
    INTEGER :: n_upper_min, n_upper_max            ! upper quantum number range
    INTEGER :: n_transitions                       ! = n_upper_max - n_upper_min + 1
    INTEGER :: n_dens                              ! actual number of density points
    REAL(8)  :: density_grid(NSTARK_DENS_MAX)       ! Ne [cm^-3]
    REAL(8)  :: temp_grid(NSTARK_TEMPS)             ! T [K]
    REAL(8)  :: log_dalpha_grid(NSTARK_DALPHA)      ! log10(Δα)
    INTEGER, ALLOCATABLE :: max_dens_idx(:)        ! (n_transitions) highest valid density
    REAL(8),  ALLOCATABLE :: k_alpha(:)             ! (n_transitions) asymptotic wing constant
    REAL(8),  ALLOCATABLE :: profiles(:,:,:,:)      ! (NSTARK_DALPHA, NSTARK_TEMPS, n_dens, n_transitions)
    LOGICAL :: loaded = .FALSE.
  END TYPE stark_series_t

  TYPE(stark_series_t), TARGET :: STEHLE_DATA(NSTARK_SERIES)
  LOGICAL :: STEHLE_TABLES_LOADED = .FALSE.

CONTAINS

  ! Unpack LINEREC(4) record into module line-data variables
  ! Replaces the removed EQUIVALENCE (LINEREC(1),IWL)
  SUBROUTINE UNPACK_LINEDATA(III)
    INTEGER(4), INTENT(IN) :: III(4)
    INTEGER(2) :: pair(2)
    IWL = III(1)
    pair = TRANSFER(III(2), pair)
    IELION = pair(1)
    IELO   = pair(2)
    pair = TRANSFER(III(3), pair)
    IGFLOG = pair(1)
    IGR    = pair(2)
    pair = TRANSFER(III(4), pair)
    IGS    = pair(1)
    IGW    = pair(2)
  END SUBROUTINE UNPACK_LINEDATA

  ! Pack module line-data variables back into LINEREC(4) record
  ! For use before WRITE(12) LINEREC when variables have been modified
  SUBROUTINE PACK_LINEDATA(III)
    INTEGER(4), INTENT(OUT) :: III(4)
    INTEGER(2) :: pair(2)
    III(1) = IWL
    pair(1) = IELION; pair(2) = IELO
    III(2) = TRANSFER(pair, III(2))
    pair(1) = IGFLOG; pair(2) = IGR
    III(3) = TRANSFER(pair, III(3))
    pair(1) = IGS; pair(2) = IGW
    III(4) = TRANSFER(pair, III(4))
  END SUBROUTINE PACK_LINEDATA

  ! Pack variables into LINEPARAM(4) matching F77 COMMON /LINEPARAM/ layout:
  ! LINEPARAM(1) = WLVAC (REAL*8)
  ! LINEPARAM(2) = [NELION (INTEGER*4), CGF (REAL*4)]
  ! LINEPARAM(3) = [ELO (REAL*4), GAMMAR (REAL*4)]
  ! LINEPARAM(4) = [GAMMAS (REAL*4), GAMMAW (REAL*4)]
  SUBROUTINE PACK_LINEPARAM(W, WLVAC_in, NELION_in, CGF_in, ELO_in, GAMMAR_in, GAMMAS_in, GAMMAW_in)
    REAL(8), INTENT(OUT) :: W(4)
    REAL(8), INTENT(IN) :: WLVAC_in
    INTEGER, INTENT(IN) :: NELION_in
    REAL(4), INTENT(IN) :: CGF_in, ELO_in, GAMMAR_in, GAMMAS_in, GAMMAW_in
    REAL(4) :: pair(2)
    W(1) = WLVAC_in
    pair(1) = TRANSFER(NELION_in, 1.0)
    pair(2) = CGF_in
    W(2) = TRANSFER(pair, 1.0D0)
    pair(1) = ELO_in
    pair(2) = GAMMAR_in
    W(3) = TRANSFER(pair, 1.0D0)
    pair(1) = GAMMAS_in
    pair(2) = GAMMAW_in
    W(4) = TRANSFER(pair, 1.0D0)
  END SUBROUTINE PACK_LINEPARAM

  ! Unpack LINEPARAM(4) into variables matching F77 COMMON /LINEPARAM/ layout
  SUBROUTINE UNPACK_LINEPARAM(W, WLVAC_out, NELION_out, CGF_out, ELO_out, GAMMAR_out, GAMMAS_out, GAMMAW_out)
    REAL(8), INTENT(IN) :: W(4)
    REAL(8), INTENT(OUT) :: WLVAC_out
    INTEGER, INTENT(OUT) :: NELION_out
    REAL(4), INTENT(OUT) :: CGF_out, ELO_out, GAMMAR_out, GAMMAS_out, GAMMAW_out
    REAL(4) :: pair(2)
    WLVAC_out = W(1)
    pair = TRANSFER(W(2), pair)
    NELION_out = TRANSFER(pair(1), 1)
    CGF_out = pair(2)
    pair = TRANSFER(W(3), pair)
    ELO_out = pair(1)
    GAMMAR_out = pair(2)
    pair = TRANSFER(W(4), pair)
    GAMMAS_out = pair(1)
    GAMMAW_out = pair(2)
  END SUBROUTINE UNPACK_LINEPARAM

  ELEMENTAL REAL(8) FUNCTION TABLOG(I)
    INTEGER(2), INTENT(IN) :: I
    TABLOG = 10.D0**((I-16384)*.001D0)
  END FUNCTION TABLOG


!=========================================================================
! SUBROUTINE PUTOUT(MODE)
!
!   Write model atmosphere results to screen (unit 6), punch file
!   (unit 7), surface flux file (unit 8), and diagnostic files.
!
!   MODE = 1 : Write column headings; initialize frequency counter
!   MODE = 2,3,4 : Print flux or intensity at current wavelength step
!   MODE = 5 : Print end-of-iteration summary; write model to punch file
!
!   Output units:
!     6  — Standard output (screen): iteration tables, diagnostics
!     7  — Punch file: model structure for restart or post-processing
!     8  — Surface flux or intensity vs. wavelength
!     50 — Optical depth tau(nu) vs. depth at each wavelength
!     66 — Compact model atmosphere summary
!=========================================================================

SUBROUTINE PUTOUT(MODE)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: MODE

  ! --- Local variables ---
  REAL(8)  :: HSURF, TAUEND
  REAL(8)  :: HLAM, HLAMLG, HLAMMG, HNULG, HNUMG
  REAL(8)  :: FLXCNVRATIO(kw)
  REAL(8)  :: SURFIN(20)       ! Note: SURFIN accumulation was disabled in original
                               ! code (opacity-sampling replaced frequency groups).
                               ! SURFI from JOSH holds the actual surface intensity.
  INTEGER :: J, I, IZ, MU, JTAU1

  ! Persistent state across calls (NU and IFHEAD survive between MODE=1 and MODE=2-4)
  INTEGER, SAVE :: NU = 0
  INTEGER, SAVE :: IFHEAD = 0

  ! =================================================================

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING PUTOUT'
  SELECT CASE (MODE)

  ! =================================================================
  ! MODE 1 — Write column headings to output and punch files
  ! =================================================================
  CASE (1)

    IF (IFPRNT(ITER) .EQ. 0) RETURN
    IFHEAD = 0
    NU = NULO - NUSTEP

    ! Write header to surface flux file (unit 8)
    IF (IFPNCH(ITER) .LT. 2) RETURN
    WRITE(8, '("TEFF ",F7.0,"  GRAVITY",F8.4,1X,A4 / "TITLE ",74A1)') &
      TEFF, GLOG, WLTE, TITLE

    ! For intensity mode (IFSURF=2), also write angle grid
    IF (IFSURF .NE. 2) RETURN
    WRITE(8, '(I3," ANGLES",10F7.4 / 10X,10F7.4)') NMU, (ANGLE(MU), MU=1,NMU)

  ! =================================================================
  ! MODE 2,3,4 — Print flux/intensity at current wavelength point
  ! =================================================================
  CASE (2, 3, 4)

    NU = NU + NUSTEP
    HSURF = HNU(1)
    TAUEND = log10(TAUNU(NRHOX))
    IF (HSURF .LE. 0.0D0) HSURF = 1.0D-50

    ! --- Detailed frequency-by-frequency output (print level >= 2) ---
    IF (IFPRNT(ITER) .GE. 2.AND.IDEBUG .EQ. 1) THEN

      ! Find first depth where optical depth exceeds unity
      JTAU1 = NRHOX
      DO J = 1, NRHOX
        IF (TAUNU(J) .GT. 1.0D0) THEN
          JTAU1 = J
          EXIT
        END IF
      END DO
      TAUEND = log10(TAUNU(NRHOX))

      ! --- Flux output format (IFSURF = 0 or 1) ---
      IF (IFSURF .EQ. 0 .OR. IFSURF .EQ. 1) THEN
        IF (IFHEAD .EQ. 0) THEN
          WRITE(6, 101)
101       FORMAT('1', ///// 10X, 'WAVE', 7X, 'HLAMBDA', 7X, 'LOG H', 7X, 'MAG', &
            10X, 'FREQUENCY', 8X, 'HNU', 10X, 'LOG H', 7X, 'MAG', 2X, 'TAUONE', &
            ' TAUNU')
        END IF
        IFHEAD = 1

        ! Convert H_nu to H_lambda and compute magnitudes
        HLAM   = HSURF * FREQ / WAVE
        HNULG  = log10(HSURF)
        HLAMLG = log10(HLAM)
        HLAMMG = -2.5D0 * HLAMLG
        HNUMG  = -2.5D0 * HNULG

        WRITE(6, 401) NU, WAVE, HLAM, HLAMLG, HLAMMG, &
                       FREQ, HSURF, HNULG, HNUMG, JTAU1, TAUEND, NU
401     FORMAT(I6, F10.3, 1PE13.4, 0PF12.5, F10.3, 1PE20.6, E13.4, &
          0PF12.5, F10.3, I6, F6.2, I6)

        ! Write tau(nu) profile to diagnostic file (unit 50)
        WRITE(50, "(I6,F12.3,(100F10.4))") NU, WAVE, &
          (log10(TAUNU(J)), J=1,NRHOX)
      END IF

      ! --- Intensity output format (IFSURF = 2) ---
      IF (IFSURF .EQ. 2) THEN
        IF (IFHEAD .EQ. 0) THEN
          WRITE(6, 102)
102       FORMAT('1', ///// 10X, 'WAVE', 5X, 'FREQUENCY', 3X, 'TAUONE TAUNU', &
            5('   MU  INTENSITY '))
        END IF
        IFHEAD = 1

        WRITE(6, 406) NU, WAVE, FREQ, JTAU1, TAUEND, &
                       (ANGLE(MU), SURFIN(MU), MU=1,NMU)
406     FORMAT(I6, F9.3, 1PE15.6, I6, 0PF6.2, 5(0PF7.4, 1PE10.3) / &
          (42X, 5(0PF7.4, 1PE10.3)))
      END IF

    END IF  ! IFPRNT >= 2

    ! --- Write to punch file (unit 8) ---
    IF (IFPNCH(ITER) .LT. 2) RETURN
    IF (IFSURF .GT. 2) RETURN

    IF (IFSURF .EQ. 2) THEN
      ! Intensity mode: write surface intensity at all angles
      WRITE(8, 416) NU, WAVE, FREQ, (SURFIN(MU), MU=1,NMU)
416   FORMAT('INTENSITY', I5, F9.2, 1PE15.6 / (1P8E10.3))
      IF (NU .EQ. NUHI) WRITE(8, 416)
    ELSE
      ! Flux mode: write wavelength (in nm) and Eddington flux H_nu
      WRITE(8, '(F13.2,E13.4)') WAVE * 10.0D0, HSURF
    END IF

  ! =================================================================
  ! MODE 5 — End-of-iteration summary and model punch output
  ! =================================================================
  CASE (5)

    ! --- Print iteration summary tables to screen ---
    IF (IFPRNT(ITER) .NE. 0) THEN

       IF (IDEBUG .EQ. 1) THEN
          ! Convection and structure parameters vs. depth
          WRITE(6, 501) (J, RHOX(J), PTOTAL(J), PTURB(J), GRDADB(J), DLTDLP(J), &
               VELSND(J), DLRDLT(J), HEATCP(J), HSCALE(J), VCONV(J), FLXCNV(J), &
               J=1,NRHOX)
501       FORMAT(///, '       RHOX       PTOTAL     PTURB      GRDADB', &
               '     DLTDLP     VELSND     DLRDLT     HEATCP     HSCALE     VCONV', &
               '     FLXCNV            ', / (I3, 1P11E11.3))
          
          ! Number densities: H I, H II, He I, He II, He III + turbulent velocity
          WRITE(6, 503) (J, XNATOM(J), RADEN(J), PRADK(J), XNFP(J,1), XNFP(J,2), &
               XNFP(J,3), XNFP(J,4), XNFP(J,5), VTURB(J), &
               FLXCNV0(J), FLXCNV1(J), J=1,NRHOX)
503       FORMAT(///, '      XNATOM      RADEN      PRADK     XNFPH1', &
               '    XNFPH2     XNFPHE1    XNFPHE2    XNFPHE3     VTURB', &
               / (I3, 1P11E11.3))
       
       ENDIF

      ! Compute convective flux fraction at each depth
      DO J = 1, NRHOX
        IF (IFCORR .EQ. 0) FLXRAD(J) = FLUX - FLXCNV(J)
        FLXCNVRATIO(J) = FLXCNV(J) / (FLXCNV(J) + FLXRAD(J))
      END DO

      ! Full model atmosphere table to unit 66
      WRITE(66, 542) (J, RHOX(J), T(J), P(J), XNE(J), RHO(J), ABROSS(J), &
        HEIGHT(J), TAUROS(J), FLXCNVRATIO(J), ACCRAD(J), FLXERR(J), FLXDRV(J), &
        J=1,NRHOX)
542   FORMAT('0', 35X, 'ELECTRON', 11X, &
        'ROSSELAND    HEIGHT   ROSSELAND   FRACTION  RADIATIVE', &
        '         PERCENT FLUX', / &
        '       RHOX      TEMP    PRESSURE    NUMBER', &
        '    DENSITY      MEAN       (KM)      DEPTH', &
        '    CONV FLUX  ACCELERATION', &
        '   ERROR   DERIV', / (I3, 1PE10.3, 0PF9.1, 1P8E11.3, 2P2E11.3))

    END IF  ! IFPRNT

    ! --- Write model to punch file (unit 7) ---
    IF (IFPNCH(ITER) .EQ. 0) RETURN

    ! Model parameters and abundances
    WRITE(7, '("TEFF ",F7.0,"  GRAVITY",F8.4,1X,A4)') TEFF, GLOG, WLTE
    WRITE(7, '("TITLE ",74A1)') TITLE

    ! Abundance table with relative offsets
    WRITE(7, 553) ELEM(1), ABUND(1), ELEM(2), ABUND(2), &
      (IZ, ELEM(IZ), ABUND(IZ), XRELATIVE(IZ), IZ=3,99)
553 FORMAT(' ABUNDANCE TABLE' / '    1', A2, F10.6, '       2', A2, F10.6 &
      / (5(I5, A2, F7.3, F6.3)))

    ! Model structure: column mass, T, P, N_e, kappa_Ross, g_rad, v_turb,
    !                  convective flux, convective velocity
    WRITE(7, 554) NRHOX, (RHOX(J), T(J), P(J), XNE(J), ABROSS(J), ACCRAD(J), &
      VTURB(J), FLXCNV(J), VCONV(J), J=1,NRHOX)
554 FORMAT('READ DECK6', I3, &
      '     RHOX         T         P       XNE', &
      '     ABROSS    ACCRAD     VTURB    FLXCNV     VCONV' &
      / (13X, 1PE12.5, 0PF10.2, 1P7E10.3))

    ! Surface radiation pressure constant
    WRITE(7, '("PRADK",1PE11.4)') PRADK0

    ! NLTE departure coefficients (if NLTE mode is on)
    IF (NLTEON .NE. 0) THEN
      WRITE(7, 556) NRHOX, (RHOX(J), (BHYD(J,I), I=1,6), BMIN(J), J=1,NRHOX)
556   FORMAT('READ DEPARTURE COEFFICIENTS', I3, ' RHOX  BHYD 1-6  BMIN' &
        / (1PE11.4, 0P7F9.4))
    END IF

    ! Signal completion
    WRITE(7, '("BEGIN",20X,"ITERATION ",I3," COMPLETED")') ITER
    CLOSE(UNIT=7)

  END SELECT

END SUBROUTINE PUTOUT


!=========================================================================
! SUBROUTINE TCORR(MODE, RCOWT)
!
! Temperature correction by the modified Avrett-Krook method.
!
! Drives the model atmosphere toward radiative (+convective) equilibrium
! by computing temperature corrections from three components:
!   DTFLUX  — flux-constancy correction (integral constraint)
!   DTLAMB  — Lambda-iteration correction (local J = S constraint)
!   DTSURF  — surface boundary correction
!
! After applying T corrections, the atmosphere is remapped onto the
! standard Rosseland optical depth grid (TAUSTD) to maintain the
! prescribed depth scale.
!
! Called in three modes per iteration:
!   MODE = 1: Initialize (zero frequency-integral accumulators)
!   MODE = 2: Accumulate frequency integrals (called per wavelength)
!   MODE = 3: Compute and apply temperature + structure corrections
!
! Key references:
!   Avrett & Krook (1963), ApJ 137, 874
!   Kurucz (1970), SAO Special Report 309
!=========================================================================

SUBROUTINE TCORR(MODE, RCOWT)

  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN)  :: MODE
  REAL(8),  INTENT(IN)  :: RCOWT   ! quadrature weight for current frequency

  ! --- Named constants ---
  ! SIGMA_SB and FOURPI now from mod_constants
  REAL(8), PARAMETER :: EULER_M1  = 0.922784335098467D0  ! 1 - gamma_Euler + ln(2)
  REAL(8), PARAMETER :: RDIAGJ_FLOOR = 1.0D-30    ! prevent division by zero
  REAL(8), PARAMETER :: ABROSS_FLOOR = 1.0D-30    ! prevent division by zero
  REAL(8), PARAMETER :: DEL_FLOOR    = 1.0D-10    ! prevent division by zero in superadiabatic excess

  ! --- Persistent locals (accumulate across MODE 1→2→3 calls) ---
  REAL(8), SAVE :: RJMINS(kw)   ! integrated opacity-weighted (J - S)
  REAL(8), SAVE :: RDABH(kw)    ! integrated flux divergence term
  REAL(8), SAVE :: RDIAGJ(kw)   ! integrated diagonal Lambda-operator response
  REAL(8), SAVE :: OLDT1(kw)    ! previous iteration's total correction (for damping)

  ! --- MODE 2 locals ---
  REAL(8)  :: DABTOT(kw)        ! d(kappa_tot)/d(RHOX)
  REAL(8)  :: TERM1, TERM2, D, EX, DIAGJ, DBDT

  ! --- MODE 3 locals: derivatives and gradients ---
  REAL(8)  :: DTDRHX(kw)        ! dT/d(RHOX)
  REAL(8)  :: DDLT(kw)          ! d(nabla)/d(RHOX)
  REAL(8)  :: DABROS(kw)        ! d(kappa_Ross)/d(RHOX)

  ! --- MODE 3 locals: convection ---
  REAL(8)  :: CNVFLX(kw)        ! local copy of convective flux (smoothed)
  REAL(8)  :: SMOOTH(kw)        ! smoothing buffer for convective flux
  REAL(8)  :: DDEL(kw)          ! convective response factor (1 + D/(D+DEL))/DEL
  REAL(8)  :: HRATIO(kw)        ! convective-to-total flux ratio
  REAL(8)  :: DEL, VCO, FLUXCO, TAUB, CNVFL

  ! --- MODE 3 locals: flux correction ---
  REAL(8)  :: CODRHX(kw)        ! integrand for G(tau) integrating factor
  REAL(8)  :: G(kw)             ! integrating factor exp(integral of CODRHX)
  REAL(8)  :: GFLUX(kw)         ! G * (flux error) / (flux response)
  REAL(8)  :: DTAU(kw)          ! integrated tau correction
  REAL(8)  :: DTFLUX(kw)        ! flux-constancy temperature correction

  ! --- MODE 3 locals: Lambda correction ---
  REAL(8)  :: DTLAMB(kw)        ! Lambda-iteration temperature correction
  REAL(8)  :: DTSURF(kw)        ! surface boundary correction
  REAL(8)  :: DTCONV(kw)        ! convective flux correction (Crivellari-Simonneau)
  REAL(8)  :: T1(kw)            ! total correction
  REAL(8)  :: TEFF25            ! clamp: Teff/25
  REAL(8)  :: DTSUR             ! surface correction magnitude
  REAL(8)  :: DUM(kw), TINTEG(kw), TAV
  REAL(8)  :: TONE_ARR(1), TTWO_ARR(1), XNEW_TMP(1)

  ! --- MODE 3 locals: convective adiabatic sweep ---
  REAL(8)  :: T_SWEEP(kw)       ! corrected T for adiabatic integration
  REAL(8)  :: DT_RAD            ! radiative correction at a layer
  REAL(8)  :: DT_ADIAB          ! adiabatic sweep correction at a layer
  REAL(8)  :: DLNP              ! Δln P between adjacent layers
  REAL(8)  :: T_ADIAB           ! adiabatic extrapolation temperature
  REAL(8)  :: FCONV_RATIO       ! convective flux fraction at a depth
  INTEGER :: JANCHOR           ! shallowest convective layer (anchor)

  ! --- MODE 3 locals: RHOX correction / remapping ---
  REAL(8)  :: TPLUS(kw), TNEW1(kw), TNEW2(kw), PRDNEW(kw)
  REAL(8)  :: AB1(kw), PTOT1(kw), P1(kw)
  REAL(8)  :: AB2(kw), PTOT2(kw), P2(kw)
  REAL(8)  :: PPP(kw), RRR(kw), DRHOX(kw)
  REAL(8)  :: REMAP(kw, 10)     ! buffer for MAP1 remapping

  ! --- Scalar temps ---
  INTEGER :: J, I, K, IDUM, IFUDGE, JCONV
  REAL(8)  :: ABROSS_safe

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING TCORR'

  !=====================================================================
  ! MODE 1: Zero frequency-integral accumulators
  !=====================================================================
  IF (MODE .EQ. 1) THEN

    RJMINS = 0.0D0
    RDABH  = 0.0D0
    RDIAGJ = 0.0D0
    FLXRAD = 0.0D0
    RETURN

  !=====================================================================
  ! MODE 2: Accumulate frequency integrals at current wavelength
  !=====================================================================
  ELSE IF (MODE .EQ. 2) THEN

    ! --- Opacity gradient term: d(ln kappa)/d(RHOX) * H_nu ---
    CALL DERIV(RHOX, ABTOT, DABTOT, NRHOX)
    RDABH  = RDABH + DABTOT / ABTOT * HNU * RCOWT
    RJMINSNU = ABTOT * JMINS * RCOWT
    RJMINS = RJMINS + RJMINSNU
    FLXRAD = FLXRAD + HNU * RCOWT

    ! --- Diagonal Lambda operator: sum kappa * (Lambda_diag - 1) * (1-eps) * dB/dT ---
    !     Lambda_diag is computed from E_3 exponential integrals of the
    !     monochromatic optical depth increments between adjacent layers.
    TERM2 = 0.0D0
    DO J = 1, NRHOX
      TERM1 = TERM2

      ! Optical depth increment to next layer
      IF (J .NE. NRHOX) THEN
        D = TAUNU(J+1) - TAUNU(J)
      END IF
      D = max(1.0D-10, D)

      IF (D .LE. 0.01D0) THEN
        ! Small optical depth: Taylor series for E_3 integral
        ! TERM2 ≈ (1 - gamma + ln2 - lnD)*D/4 + D^2/12 - D^3/96 + D^4/720
        TERM2 = (EULER_M1 - log(D)) * D / 4.0D0 &
              + D**2 / 12.0D0 - D**3 / 96.0D0 + D**4 / 720.0D0
      ELSE
        ! Standard E_3 path
        EX = 0.0D0
        IF (D .LT. 10.0D0) EX = EXPI(3, D)
        ! Cool-star stability patch: suppress E_3 for narrow Dtau range
        IF (TEFF .LE. 4250.0D0 .AND. D .GT. 0.005D0 .AND. D .LT. 0.02D0) EX = 0.0D0
        TERM2 = 0.5D0 * (D + EX - 0.5D0) / D
      END IF

      DIAGJ = TERM1 + TERM2

      ! dB_nu/dT = B_nu * h*nu / (kT^2) / (1 - e^{-h*nu/kT})
      IF (NUMNU .EQ. 1) THEN
        DBDT = FLUX * 16.0D0 / T(J)
      ELSE
        ! Guard: STIM = 1 - exp(-hv/kT) can vanish in Rayleigh-Jeans limit
        DBDT = BNU(J) * FREQ * HKT(J) / T(J) / max(STIM(J), 1.0D-30)
      END IF

      ! Accumulate: kappa * (Lambda_diag - 1) / (1 - eps*Lambda_diag) * (1-eps) * dB/dT
      RDIAGJNU(J) = ABTOT(J) * (DIAGJ - 1.0D0) &
                   / (1.0D0 - ALPHA(J) * DIAGJ) &
                   * (1.0D0 - ALPHA(J)) * DBDT * RCOWT
      RDIAGJ(J) = RDIAGJ(J) + RDIAGJNU(J)
    END DO
    RETURN

  END IF

  !=====================================================================
  ! MODE 3: Compute and apply temperature corrections
  !=====================================================================

  ! --- Compute needed derivatives ---
  CALL DERIV(RHOX, T, DTDRHX, NRHOX)
  CALL DERIV(RHOX, DLTDLP, DDLT, NRHOX)
  CALL DERIV(RHOX, ABROSS, DABROS, NRHOX)

  !---------------------------------------------------------------------
  ! (A) Prepare smoothed convective flux
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    CNVFLX(J) = 0.0D0
    IF (IFCONV .EQ. 1 .AND. J .GE. 3) CNVFLX(J) = FLXCNV(J)
  END DO

  ! 1-2-1 smoothing filter on convective flux (interior points)
  DO J = 2, NRHOX - 1
    SMOOTH(J) = 0.25D0 * CNVFLX(J-1) + 0.50D0 * CNVFLX(J) + 0.25D0 * CNVFLX(J+1)
  END DO
  ! Asymmetric boundary kernel at deepest layer: 75-25 split
  SMOOTH(NRHOX) = 0.75D0 * CNVFLX(NRHOX) + 0.25D0 * CNVFLX(NRHOX-1)
  DO J = 2, NRHOX
    CNVFLX(J) = SMOOTH(J)
    FLXCNV(J) = CNVFLX(J)
  END DO

  !---------------------------------------------------------------------
  ! (B) Build integrating factor for flux correction
  !     CODRHX = d(ln G)/d(RHOX) incorporates opacity gradient and
  !     convective flux response to temperature perturbations
  !---------------------------------------------------------------------

  ! Initialize DDEL for all layers (safe default for non-convective case)
  DDEL(:) = 1.0D0

  DO J = 1, NRHOX
    ABROSS_safe = max(ABROSS(J), ABROSS_FLOOR)
    RDABH(J) = RDABH(J) - FLXRAD(J) * DABROS(J) / ABROSS_safe

    ! Convective response: DEL = superadiabatic excess, D = radiative leak parameter
    DEL  = 1.0D0
    D    = 0.0D0

    IF (CNVFLX(J) .GT. 0.0D0 .AND. FLXCNV0(J) .GT. 0.0D0) THEN
      DEL = DLTDLP(J) - GRDADB(J)
      DEL = max(DEL, DEL_FLOOR)

      VCO = 0.5D0 * MIXLTH * sqrt(max(-0.5D0 * PTOTAL(J) / RHO(J) * DLRDLT(J), 0.0D0))
      FLUXCO = 0.5D0 * RHO(J) * HEATCP(J) * T(J) * MIXLTH / FOURPI

      IF (MIXLTH .GT. 0.0D0 .AND. VCO .GT. 0.0D0) THEN
        D = 8.0D0 * SIGMA_SB * T(J)**4 &
          / (ABROSS_safe * HSCALE(J) * RHO(J)) / (FLUXCO * FOURPI) / VCO
      END IF

      TAUB = ABROSS_safe * RHO(J) * MIXLTH * HSCALE(J)
      D = D * TAUB**2 / (2.0D0 + TAUB**2)
      D = D**2 / 2.0D0
      DDEL(J) = (1.0D0 + D / (D + DEL)) / DEL
    END IF

    ! Only include convective coupling when it's significant
    CNVFL = 0.0D0
    IF (max(FLXRAD(J), 1.0D-30) .GT. 0.0D0) THEN
      IF (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) .GT. 1.0D-3 .AND. &
          FLXCNV0(J) / max(FLXRAD(J), 1.0D-30) .GT. 1.0D-3) THEN
        CNVFL = CNVFLX(J)
      END IF
    END IF

    CODRHX(J) = (RDABH(J) &
               + CNVFL * (DTDRHX(J) / T(J) * (1.0D0 - 9.0D0 * D / (D + DEL)) &
                        + 1.5D0 * DDLT(J) / DEL * (1.0D0 + D / (D + DEL)))) &
              / (FLXRAD(J) + CNVFLX(J) * 1.5D0 * DLTDLP(J) / DEL &
                           * (1.0D0 + D / (D + DEL)))
  END DO

  ! Force zero at boundary (no correction at outermost layers)
  CODRHX(1) = 0.0D0
  CODRHX(2) = 0.0D0

  !---------------------------------------------------------------------
  ! (C) DTFLUX: flux-constancy correction
  !     Integrate CODRHX to get G(RHOX), then integrate G*(flux error)
  !     over tau_Ross and convert to temperature correction.
  !
  !     Crivellari-Simonneau modification: suppress the GFLUX integrand
  !     in convection-dominated layers. A persistent flux error at depth
  !     (from MLT or opacity sampling limitations) accumulates in the
  !     DTAU integral and corrupts DTFLUX at all depths. Weight GFLUX
  !     by (1 - f_conv) so deep convective layers contribute negligibly.
  !---------------------------------------------------------------------
  CALL INTEG(RHOX, CODRHX, G, NRHOX, 0.0D0)
  DO J = 1, NRHOX
    G(J) = exp(G(J))
    GFLUX(J) = G(J) * (FLXRAD(J) + CNVFLX(J) - FLUX) &
             / (FLXRAD(J) + CNVFLX(J) * 1.5D0 * DLTDLP(J) * DDEL(J))

    ! Suppress GFLUX in convection-dominated layers (with hysteresis)
    IF (IFCONV .EQ. 1 .AND. &
        (CNVFLX(J) .GT. 0.0D0 .OR. FLXCNV0(J) .GT. 0.0D0)) THEN
      FCONV_RATIO = max(CNVFLX(J), FLXCNV0(J)) &
                  / (max(CNVFLX(J), FLXCNV0(J)) + max(FLXRAD(J), 1.0D-30))
      GFLUX(J) = GFLUX(J) * (1.0D0 - FCONV_RATIO)
    END IF
  END DO

  CALL INTEG(TAUROS, GFLUX, DTAU, NRHOX, 0.0D0)
  DO J = 1, NRHOX
    DTAU(J) = DTAU(J) / G(J)
    ! Clamp tau correction to ±tau/3 to prevent overshooting
    DTAU(J) = max(-TAUROS(J) / 3.0D0, min(TAUROS(J) / 3.0D0, DTAU(J)))
    ABROSS_safe = max(ABROSS(J), ABROSS_FLOOR)
    DTFLUX(J) = -DTAU(J) * DTDRHX(J) / ABROSS_safe
  END DO

  TEFF25 = TEFF / 25.0D0

  !---------------------------------------------------------------------
  ! (C2) Suppress DTFLUX in convective layers.
  !      The Avrett-Krook integral handles radiative flux errors only;
  !      in convective layers, DTFLUX is attenuated by (1 - f_conv).
  !      The convective correction (adiabatic sweep) is applied later,
  !      after DTLAMB and DTSURF have been computed.
  !---------------------------------------------------------------------
  DTCONV(:) = 0.0D0

  IF (IFCONV .EQ. 1) THEN
    DO J = 1, NRHOX
      IF (CNVFLX(J) .GT. 0.0D0 .OR. FLXCNV0(J) .GT. 0.0D0) THEN
        FCONV_RATIO = max(CNVFLX(J), FLXCNV0(J)) &
                    / (max(CNVFLX(J), FLXCNV0(J)) + max(FLXRAD(J), 1.0D-30))
        DTFLUX(J) = DTFLUX(J) * (1.0D0 - FCONV_RATIO)
      END IF
    END DO
  END IF

  !---------------------------------------------------------------------
  ! (D) DTLAMB: Lambda-iteration correction (optically thin layers)
  !     Uses thermal imbalance (J - S) and diagonal Lambda operator
  !     to drive local radiative balance.
  !---------------------------------------------------------------------

  ! Flux error as percentage
  FLXERR = (FLXRAD + CNVFLX - FLUX) / FLUX * 100.0D0
  CALL DERIV(TAUROS, FLXERR, FLXDRV, NRHOX)

  ! Find first layer where convection becomes significant
  ! (used to safely bound the DTLAMB backward-damping reach)
  JCONV = NRHOX + 1
  DO J = 1, NRHOX
    IF (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) .GE. 1.0D-5 .OR. TAUROS(J) .GE. 1.0D0) THEN
      JCONV = J
      EXIT
    END IF
  END DO

  DO J = 1, NRHOX
    ! In radiative layers, replace flux derivative with thermal imbalance
    IF (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) .LT. 1.0D-5) THEN
      ABROSS_safe = max(ABROSS(J), ABROSS_FLOOR)
      FLXDRV(J) = RJMINS(J) / ABROSS_safe / FLUX * 100.0D0
    END IF

    ! Lambda correction: DT = -(dF/dtau) / (dRDIAGJ) * kappa_Ross
    DTLAMB(J) = -FLXDRV(J) * FLUX / 100.0D0 &
              / max(abs(RDIAGJ(J)), RDIAGJ_FLOOR) * sign(1.0D0, RDIAGJ(J)) &
              * max(ABROSS(J), ABROSS_FLOOR)

    ! Zero DTLAMB in convective / optically thick layers, and
    ! damp preceding layers to ensure smooth transition
    IF (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) .GE. 1.0D-5 .OR. TAUROS(J) .GE. 1.0D0) THEN
      DTLAMB(J) = 0.0D0
      ! Damp 5 preceding layers (with bounds check)
      DO K = 1, 5
        IF (J - K .GE. 1) DTLAMB(J - K) = DTLAMB(J - K) / 2.0D0
      END DO
    END IF

    ! Clamp to ±Teff/25
    DTLAMB(J) = max(-TEFF25, min(TEFF25, DTLAMB(J)))
  END DO

  !---------------------------------------------------------------------
  ! (E) DTSURF: surface boundary correction
  !     Uniform shift from flux error at tau=0, adjusted to not fight
  !     the integral of (DTFLUX + DTLAMB) between tau=0.1 and tau=2.
  !---------------------------------------------------------------------
  DTSUR = (FLUX - FLXRAD(1)) / FLUX * 0.25D0 * T(1)
  DTSUR = max(-TEFF25, min(TEFF25, DTSUR))

  DUM = DTFLUX + DTLAMB
  CALL INTEG(TAUROS, DUM, TINTEG, NRHOX, 0.0D0)
  XNEW_TMP(1) = 0.1D0
  IDUM = MAP1(TAUROS, TINTEG, NRHOX, XNEW_TMP, TONE_ARR, 1)
  XNEW_TMP(1) = 2.0D0
  IDUM = MAP1(TAUROS, TINTEG, NRHOX, XNEW_TMP, TTWO_ARR, 1)

  TAV = (TTWO_ARR(1) - TONE_ARR(1)) / 2.0D0
  IF (DTSUR * TAV .LE. 0.0D0) TAV = 0.0D0
  IF (abs(TAV) .GT. abs(DTSUR)) TAV = DTSUR
  DTSUR = DTSUR - TAV

  DTSURF = DTSUR
  HRATIO = CNVFLX / (CNVFLX + max(FLXRAD, 1.0D-30))

  !---------------------------------------------------------------------
  ! (E2) DTCONV: convective temperature correction via adiabatic sweep
  !
  !      Now that DTFLUX, DTLAMB, and DTSURF are all computed, we can
  !      build the adiabatic sweep. Strategy:
  !        1. Find the anchor (shallowest convective layer)
  !        2. Anchor T = current T + radiative corrections
  !        3. Sweep downward: each layer's T follows the adiabat from
  !           the corrected layer above
  !        4. Blend using f_conv: radiative corrections at low f_conv,
  !           adiabatic sweep at high f_conv
  !---------------------------------------------------------------------

  IF (IFCONV .EQ. 1) THEN

    ! --- Find anchor: shallowest layer with convective flux ---
    JANCHOR = 0
    DO J = 1, NRHOX
      IF (CNVFLX(J) .GT. 0.0D0 .OR. FLXCNV0(J) .GT. 0.0D0) THEN
        JANCHOR = J
        EXIT
      END IF
    END DO

    ! --- Adiabatic sweep from anchor downward ---
    IF (JANCHOR .GT. 0) THEN
      ! Anchor gets normal radiative corrections only (DTCONV=0 there).
      ! T_SWEEP tracks the corrected temperature for integration.
      T_SWEEP(JANCHOR) = T(JANCHOR) + DTFLUX(JANCHOR) &
                       + DTLAMB(JANCHOR) + DTSURF(JANCHOR)

      DO J = JANCHOR + 1, NRHOX
        ! Convective fraction at this layer (with hysteresis)
        IF (CNVFLX(J) .GT. 0.0D0 .OR. FLXCNV0(J) .GT. 0.0D0) THEN
          FCONV_RATIO = max(CNVFLX(J), FLXCNV0(J)) &
                      / (max(CNVFLX(J), FLXCNV0(J)) + max(FLXRAD(J), 1.0D-30))
        ELSE
          FCONV_RATIO = 0.0D0
        END IF

        IF (FCONV_RATIO .GT. 0.0D0 .AND. PTOTAL(J) .GT. 0.0D0 &
            .AND. PTOTAL(J-1) .GT. 0.0D0) THEN
          ! Adiabatic extrapolation from corrected layer above
          DLNP = log(PTOTAL(J) / PTOTAL(J-1))
          T_ADIAB = T_SWEEP(J-1) * exp(GRDADB(J) * DLNP)

          ! Normal (radiative) correction for this layer
          DT_RAD = DTFLUX(J) + DTLAMB(J) + DTSURF(J)

          ! Adiabatic correction: what sweep says T should be minus current T
          DT_ADIAB = T_ADIAB - T(J)

          ! Blend: f_conv=0 → pure radiative, f_conv=1 → pure adiabat
          DTCONV(J) = FCONV_RATIO * DT_ADIAB &
                    + (1.0D0 - FCONV_RATIO) * DT_RAD - DT_RAD
          ! Note: total correction T1 = DT_RAD + DTCONV, so
          !   T1 = DT_RAD + f*DT_ADIAB + (1-f)*DT_RAD - DT_RAD
          !      = f*DT_ADIAB + (1-f)*DT_RAD   ← the desired blend

          ! Clamp to ±Teff/25
          DTCONV(J) = max(-TEFF25, min(TEFF25, DTCONV(J)))

          ! Update sweep temperature to reflect what we actually apply
          ! so next layer integrates from the right place
          T_SWEEP(J) = T(J) + DT_RAD + DTCONV(J)
        ELSE
          ! Non-convective layer below anchor: no DTCONV, pass through
          T_SWEEP(J) = T(J) + DTFLUX(J) + DTLAMB(J) + DTSURF(J)
        END IF
      END DO
    END IF
  END IF

  !---------------------------------------------------------------------
  ! (F) Total correction and iteration damping
  !---------------------------------------------------------------------
  T1 = DTFLUX + DTLAMB + DTSURF + DTCONV

  ! Iteration damping: compare current correction T1 against previous
  ! iteration's correction OLDT1 to detect convergence behavior.
  !   Same sign, shrinking  → accelerate by 1.25x (monotone convergence)
  !   Same sign, growing    → cap at previous magnitude (prevent runaway)
  !   Sign flip             → damp by 0.5x (oscillation control)
  ! Skip damping on first iteration only (no history).
  DO J = 1, NRHOX
    IF (ITER .EQ. 1) THEN
      ! No damping — accept raw correction
    ELSE IF (OLDT1(J) * T1(J) .GT. 0.0D0) THEN
      ! Same sign as previous iteration
      IF (abs(T1(J)) .LT. abs(OLDT1(J))) THEN
        ! Shrinking: accelerate by 25%
        T1(J) = T1(J) * 1.25D0
      ELSE
        ! Growing: cap at previous magnitude to prevent runaway
        T1(J) = sign(abs(OLDT1(J)), T1(J))
      END IF
    ELSE IF (OLDT1(J) * T1(J) .LT. 0.0D0) THEN
      ! Sign flip: damp by 50%
      T1(J) = T1(J) * 0.5D0
    END IF
    OLDT1(J) = T1(J)
  END DO

  ! Diagnostic output (after damping, so T1 reflects what is actually applied)
  IF (IFPRNT(ITER) .NE. 0) THEN
    WRITE(67, 100) (J, log10(max(TAUROS(J),1.0D-30)), T(J), DTLAMB(J), &
                    DTSURF(J), DTFLUX(J), DTCONV(J), T1(J), HRATIO(J), &
                    FLXERR(J), FLXDRV(J), DLTDLP(J), GRDADB(J), J=1,NRHOX)
100 FORMAT('0', 2X, 'lgTAUROS', 6X, 'T', 6X, &
      'DTLAMB   DTSURF   DTFLUX   DTCONV', 5X, &
      'T1   CONV/TOTAL      ERROR     DERIV   NABLA  NABLA_AD', / &
      (I3, F8.3, F10.1, 5F9.1, 1X, 1PE11.3, 1X, &
       0P, 2F10.3, 2F8.4))
    flush(67)
  END IF

  !---------------------------------------------------------------------
  ! (G) Compute RHOX correction to maintain constant TAUROS grid
  !     Uses finite-difference: run TTAUP with T and T+DT, compare
  !     total pressures to infer needed RHOX adjustment.
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    TPLUS(J)  = T(J) + T1(J)
    TAUSTD(J) = 10.0D0**(TAU1LG + (J - 1) * STEPLG)
  END DO

  IDUM = MAP1(TAUROS, T,    NRHOX, TAUSTD, TNEW1,  NRHOX)
  IDUM = MAP1(TAUROS, PRAD, NRHOX, TAUSTD, PRDNEW, NRHOX)
  CALL TTAUP(TNEW1, TAUSTD, AB1, PTOT1, P1, PRDNEW, PTURB, VTURB, GRAV, NRHOX)

  IDUM = MAP1(TAUROS, TPLUS, NRHOX, TAUSTD, TNEW2, NRHOX)
  CALL TTAUP(TNEW2, TAUSTD, AB2, PTOT2, P2, PRDNEW, PTURB, VTURB, GRAV, NRHOX)

  PPP = (PTOT2 - PTOT1) / PTOT1
  IDUM = MAP1(TAUSTD, PPP, NRHOX, TAUROS, RRR, NRHOX)
  DRHOX = RRR * RHOX

  !---------------------------------------------------------------------
  ! (H) Apply temperature correction
  !---------------------------------------------------------------------
   !apply damping to the temperature correction
  ! T1(J) = sign(min(abs(T1(J)), 0.02D0 * T(J)), T1(J))
   T = T + T1

  ! Optional smoothing
  IF (J1SMOOTH .GT. 0) THEN
    DO J = J1SMOOTH, J2SMOOTH
      TSMOOTH(J) = WTJM1 * T(J-1) + WTJ * T(J) + WTJP1 * T(J+1)
    END DO
    DO J = J1SMOOTH, J2SMOOTH
      T(J) = TSMOOTH(J)
    END DO
  END IF

  ! Force monotonicity: T must increase inward
  DO I = 2, NRHOX
    J = NRHOX + 1 - I
    T(J) = min(T(J), T(J+1) - 1.0D0)
  END DO

  ! Floor temperature options
  ! various opacity tables only extend to 2000K
   T = max(T, 2000.0D0)
  ! T(J) = max(T(J),0.7*TEFF)

  !---------------------------------------------------------------------
  ! (I) Recompute thermodynamic quantities from corrected T
  !---------------------------------------------------------------------
  IFUDGE = 0

  TK   = KBOL * T
  HKT  = HPLANCK / TK
  HCKT = HKT * CLIGHT
  TKEV = KBOL_EV * T
  TLOG = log(T)

  IF (IFUDGE .EQ. 1) RETURN

  !---------------------------------------------------------------------
  ! (J) Apply RHOX correction and remap atmosphere onto standard
  !     Rosseland optical depth grid TAUSTD
  !---------------------------------------------------------------------
  RHOX = RHOX + DRHOX

  ! Remap all atmospheric quantities from old TAUROS onto TAUSTD
  IDUM = MAP1(TAUROS, RHOX,   NRHOX, TAUSTD, REMAP(1,1),  NRHOX)
  IDUM = MAP1(TAUROS, T,      NRHOX, TAUSTD, REMAP(1,2),  NRHOX)
  IDUM = MAP1(TAUROS, P,      NRHOX, TAUSTD, REMAP(1,3),  NRHOX)
  IDUM = MAP1(TAUROS, XNE,    NRHOX, TAUSTD, REMAP(1,4),  NRHOX)
  IDUM = MAP1(TAUROS, ABROSS, NRHOX, TAUSTD, REMAP(1,5),  NRHOX)
  IDUM = MAP1(TAUROS, PRAD,   NRHOX, TAUSTD, REMAP(1,6),  NRHOX)
  IDUM = MAP1(TAUROS, VTURB,  NRHOX, TAUSTD, REMAP(1,7),  NRHOX)
  IDUM = MAP1(TAUROS, BMIN,   NRHOX, TAUSTD, REMAP(1,8),  NRHOX)
  IDUM = MAP1(TAUROS, PTURB,  NRHOX, TAUSTD, REMAP(1,9),  NRHOX)
  IDUM = MAP1(TAUROS, ACCRAD, NRHOX, TAUSTD, REMAP(1,10), NRHOX)

  ! Guard: if TAUROS grid starts deeper than TAUSTD, extrapolated
  ! ABROSS can go negative. Fix by clamping to surface value.
  DO J = 1, NRHOX
    IF (TAUROS(1) .GT. TAUSTD(J) .AND. REMAP(J,5) .LT. 0.0D0) THEN
      REMAP(J,5) = ABROSS(1)
    END IF
  END DO

  ! Copy remapped values back to module arrays
  RHOX   = REMAP(:, 1)
  T      = REMAP(:, 2)
  TK     = KBOL * T
  HKT    = HPLANCK / TK
  HCKT   = HKT * CLIGHT
  TKEV   = KBOL_EV * T
  TLOG   = log(T)
  P      = REMAP(:, 3)
  XNE    = REMAP(:, 4)
  ABROSS = REMAP(:, 5)
  PRAD   = REMAP(:, 6)
  PRADK  = PRAD + PRADK0
  VTURB  = REMAP(:, 7)
  BMIN   = REMAP(:, 8)
  PTURB  = REMAP(:, 9)

  ! Remap hydrogen NLTE departure coefficients
  DO I = 1, 6
    IDUM = MAP1(TAUROS, BHYD(1,I), NRHOX, TAUSTD, REMAP(1,1), NRHOX)
    BHYD(:, I) = REMAP(:, 1)
  END DO

  ! Replace old TAUROS with the standard grid
  TAUROS = TAUSTD

  RETURN

END SUBROUTINE TCORR

!=========================================================================
! SUBROUTINE STATEQ(MODE, RCOWT)
!
! NLTE statistical equilibrium for hydrogen and H-.
!
! Solves for departure coefficients b_n (= n_n / n_n^*) for the first
! 6 bound levels of hydrogen plus the H- ion, using accumulated
! radiative rates from the frequency integration and collisional rates
! computed from analytic fits.
!
! Called in three modes per iteration:
!   MODE = 1: Initialize (zero radiative rate accumulators, save T)
!   MODE = 2: Accumulate radiative rates (called per wavelength)
!   MODE = 3: Compute collision rates, solve rate equations, apply
!
! Hydrogen model atom: 8 levels (6 bound + level 7 = H+ + e, level 8 = H-)
! Rate matrix solved: 6x6 system for levels 1-6, with levels 7-8 as
!   sinks feeding the right-hand side.
!
! H- departure coefficient: computed separately from detailed balance
!   of radiative and collisional attachment/detachment rates.
!
! References:
!   Collision rates: Burke, Ormonde & Whitaker (1968), Proc.Phys.Soc. 92, 319
!   Cross-section formula: Q_ij = 4*f_ij*(E_H/E_0)^2 * [ln(E/E_0)/(E/E_0) + 0.148/(E/E_0)^6]
!   Implementation: D.M. Peterson, May 1968
!=========================================================================

SUBROUTINE STATEQ(MODE, RCOWT)
  
  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN) :: MODE
  REAL(8),  INTENT(IN) :: RCOWT    ! quadrature weight for current frequency

  ! --- Number of bound levels in the model atom ---
  INTEGER, PARAMETER :: NLEV = 6   ! bound levels solved explicitly
  INTEGER, PARAMETER :: NATOM = 8  ! total levels (6 bound + H+ + H-)

  ! --- Persistent locals (accumulate across MODE 1 → 2 → 3) ---
  REAL(8), SAVE :: QRADIK(kw, NLEV)  ! radiative photoionization rates (level i → continuum)
  REAL(8), SAVE :: QRADKI(kw, NLEV)  ! radiative recombination rates (continuum → level i)
  REAL(8), SAVE :: DQRAD(kw, NLEV)   ! dQ_recomb/dT (for first-order T correction)
  REAL(8), SAVE :: QRDHMK(kw)        ! H- radiative detachment rate (H- → H + e)
  REAL(8), SAVE :: QRDKHM(kw)        ! H- radiative attachment rate (H + e → H-)
  REAL(8), SAVE :: DQRD(kw)          ! dQ_Hminus_attach/dT
  REAL(8), SAVE :: TOLD(kw)          ! temperature at MODE=1 (for dT correction)

  ! --- MODE 2 locals ---
  REAL(8)  :: HCONT(NLEV)         ! hydrogen bound-free cross-sections at current freq
  REAL(8)  :: HMINBF              ! H- bound-free cross-section at current freq
  REAL(8)  :: RFRWT               ! 4*pi / (h*nu) * RCOWT — rate prefactor
  REAL(8)  :: HVC                 ! 2*h*nu*(nu/c)^2 — stimulated emission correction
  REAL(8)  :: RJ, RJE, RJEDT     ! rate building blocks at each depth

  ! --- MODE 3 locals ---
  REAL(8)  :: QCOLL(NATOM, NATOM) ! electron collision rate matrix
  REAL(8)  :: A(NLEV, NLEV)       ! rate equation matrix (destroyed by SOLVIT)
  REAL(8)  :: RIGHT(NLEV)         ! RHS vector → solution vector
  INTEGER :: IPIVOT(NLEV)        ! pivot scratch for SOLVIT
  REAL(8)  :: DT                  ! temperature change since MODE 1
  REAL(8)  :: THETA               ! 5040/T (excitation temperature parameter)
  REAL(8)  :: TH                  ! ionization potential / kT = 13.595 / T_eV
  REAL(8)  :: Y, Z_lev            ! real-valued level indices
  REAL(8)  :: GIK, X0, Q          ! collision rate intermediates
  REAL(8)  :: QELECT              ! H- collisional detachment by electrons
  REAL(8)  :: QASSOC              ! H- associative detachment (H + H → H2 + e)
  REAL(8)  :: QCHARG              ! H- charge exchange (H- + H+ → 2H)
  REAL(8)  :: DENOM               ! denominator for BMIN (with guard)

  ! --- Loop indices ---
  INTEGER :: I, J, K, L, N

  ! --- Oscillator strengths for bound-bound transitions ---
  ! F(I,K) = oscillator strength for transition I → K
  ! From Peterson (1968); upper triangle only (F is symmetric by convention here)
  REAL(8), PARAMETER :: F(NATOM, NATOM) = reshape( [ &
    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,   &
    0.4162D0, 0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,   &
    0.07910D0,0.6408D0, 0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,   &
    0.02899D0,0.1193D0, 0.8420D0, 0.0D0,    0.0D0,    0.0D0,    0.0D0,    0.0D0,   &
    0.01394D0,0.04467D0,0.1506D0, 1.038D0,  0.0D0,    0.0D0,    0.0D0,    0.0D0,   &
    0.007800D0,0.02209D0,0.05585D0,0.1794D0,1.231D0,  0.0D0,    0.0D0,    0.0D0,   &
    0.004814D0,0.01271D0,0.02768D0,0.06551D0,0.2070D0,1.425D0,  0.0D0,    0.0D0,   &
    0.003184D0,0.008037D0,0.01604D0,0.03229D0,0.07455D0,0.2340D0,1.615D0,0.0D0     &
    ], shape=[NATOM, NATOM] )

  ! --- Physical constants now from mod_constants ---

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING STATEQ'

  !=====================================================================
  ! MODE 1: Zero radiative rate accumulators and save current temperature
  !=====================================================================
  IF (MODE .EQ. 1) THEN

    QRADIK = 0.0D0
    QRADKI = 0.0D0
    DQRAD  = 0.0D0
    TOLD   = T
    QRDHMK = 0.0D0
    QRDKHM = 0.0D0
    DQRD   = 0.0D0
    RETURN

  !=====================================================================
  ! MODE 2: Accumulate radiative rates at current frequency
  !=====================================================================
  ELSE IF (MODE .EQ. 2) THEN

    ! Rate prefactor: (4*pi / h*nu) * quadrature weight
    RFRWT = FOURPI / HPLANCK * RCOWT / FREQ

    ! Stimulated emission term: 2*h*nu * (nu/c)^2
    HVC = 2.0D0 * HPLANCK * FREQ * (FREQ / CLIGHT)**2

    ! Hydrogen bound-free cross-sections for levels 2-6
    DO N = 2, NLEV
      HCONT(N) = COULX(N, FREQ, 1.0D0)
    END DO

    ! H- bound-free cross-section (polynomial fit)
    HMINBF = 0.0D0
    IF (FREQ .GT. 1.8259D14 .AND. FREQ .LT. 2.111D14) THEN
      HMINBF = 3.695D-16 + (-1.251D-1 + 1.052D13 / FREQ) / FREQ
    END IF
    IF (FREQ .GE. 2.111D14) THEN
      HMINBF = 6.801D-20 + (5.358D-3 + (1.481D13 &
             + (-5.519D27 + 4.808D41 / FREQ) / FREQ) / FREQ) / FREQ
    END IF

    ! Accumulate rates at each depth
    DO J = 1, NRHOX
      ! RJ  = (4*pi*J_nu / h*nu) * weight — photoionization driver
      ! RJE = (4*pi*e^{-hv/kT} * (J_nu + 2hv^3/c^2) / h*nu) * weight — recombination
      ! RJEDT = d(RJE)/dT — temperature derivative for MODE 3 correction
      RJ    = RFRWT * JNU(J)
      RJE   = RFRWT * EHVKT(J) * (JNU(J) + HVC)
      RJEDT = RJE * HKT(J) * FREQ / T(J)

      ! Hydrogen levels 2-6 (level 1 has no photoionization at these frequencies)
      DO I = 2, NLEV
        QRADIK(J, I) = QRADIK(J, I) + HCONT(I) * RJ
        QRADKI(J, I) = QRADKI(J, I) + HCONT(I) * RJE
        DQRAD(J, I)  = DQRAD(J, I)  + HCONT(I) * RJEDT
      END DO

      ! H- rates
      QRDHMK(J) = QRDHMK(J) + HMINBF * RJ
      QRDKHM(J) = QRDKHM(J) + HMINBF * RJE
      DQRD(J)   = DQRD(J)   + HMINBF * RJEDT
    END DO
    RETURN

  END IF

  !=====================================================================
  ! MODE 3: Compute collision rates, solve rate equations, update BHYD/BMIN
  !=====================================================================

  !-------------------------------------------------------------------
  ! (A) H- departure coefficient: detailed balance of formation/destruction
  !-------------------------------------------------------------------
  IF (IFPRNT(ITER) .GT. 0) THEN
    WRITE(6, 201)
201 FORMAT('1', /////, 36X, 'HMINUS STATISTICAL EQUILIBRIUM' / &
      10X, 'RHOX', 7X, 'QELECT', 6X, 'QASSOC', 6X, 'QCHARG', &
      6X, 'QRDKHM', 6X, 'QRDHMK', 7X, 'BMIN')
  END IF

  DO J = 1, NRHOX
    DT    = T(J) - TOLD(J)
    THETA = THETA_COEFF / T(J)

    ! Collisional detachment by electrons: rate ∝ theta^{3/2} * n_e
    QELECT = 10.0D0**(-8.7D0) * THETA**1.5D0 * XNE(J)

    ! Associative detachment: H + H → H2 + e  (rate ∝ b_1 * n(H))
    QASSOC = 10.0D0**(-8.7D0) * 2.0D0 * BHYD(J,1) * XNFP(J,1)

    ! Charge exchange: H- + H+ → 2H  (rate ∝ theta^{1/3} * n(H+))
    QCHARG = 10.0D0**(-7.4D0) * THETA**0.333333D0 * XNFP(J,2)

    ! First-order correction of recombination rate for T change
    QRDKHM(J) = QRDKHM(J) + DQRD(J) * DT

    ! b(H-) = total formation / total destruction
    DENOM = QRDHMK(J) + QELECT + QASSOC + QCHARG
    BMIN(J) = (QRDKHM(J) + QELECT + QASSOC + QCHARG) / max(DENOM, 1.0D-30)

    IF (IFPRNT(ITER) .GT. 0) THEN
      WRITE(6, 211) J, RHOX(J), QELECT, QASSOC, QCHARG, &
                     QRDKHM(J), QRDHMK(J), BMIN(J)
211   FORMAT(I5, 1P6E12.3, 0PF10.4)
    END IF
  END DO

  !-------------------------------------------------------------------
  ! (B) Hydrogen 6-level atom: build and solve rate equations
  !-------------------------------------------------------------------
  IF (IFPRNT(ITER) .GT. 0) THEN
    WRITE(6, 31)
31  FORMAT('1', /////, 30X, &
      'STATISTICAL EQUILIBRIUM RATES    RATE=SIGN(ALOG10(MAX(ABS(RATE*1.E20),1.)),RATE)', /, &
      '0 RAD   1-K   K-1   2-K   K-2   3-K   K-3   4-K   K-4   5-K', &
      '   K-5   6-K   K-6   COLL  1-K   2-K   3-K   4-K   5-K   6-K   5-8', &
      '   6-8  ', /, &
      '  COLL  1-2   1-3   1-4   1-5   1-6   1-7   2-3   2-4   2-5', &
      '   2-6   2-7   3-4   3-5   3-6   3-7   4-5   4-6   4-7   5-6   5-7', &
      '   6-7  ')
  END IF

  DO J = 1, NRHOX
    DT = T(J) - TOLD(J)
    TH = 13.595D0 / TKEV(J)    ! chi_H / kT

    !--- Electron collision rates ---
    ! Diagonal: bound-free collisional ionization
    ! Off-diagonal: bound-bound collisional excitation/de-excitation
    DO I = 1, NATOM
      Y = dble(I)

      ! Bound-free collisional ionization rate (Seaton-like)
      QCOLL(I, I) = 2.2D-8 * Y**3 / sqrt(TH) * exp(-TH / Y**2) * XNE(J)

      IF (I .LT. NATOM) THEN
        DO K = I + 1, NATOM
          Z_lev = dble(K)
          GIK = 1.0D0 / Y**2 - 1.0D0 / Z_lev**2
          X0  = TH * GIK

          ! Burke, Ormonde & Whitaker cross-section with E_1 and E_5 integrals
          Q = 2.186D-10 * F(I, K) / GIK**2 * X0 * sqrt(T(J)) &
            * (EXPI(1, X0) + 0.148D0 * X0 * EXPI(5, X0))

          QCOLL(I, K) = Q * XNE(J)
          ! Detailed balance: reverse rate includes Boltzmann factor
          QCOLL(K, I) = QCOLL(I, K) * (Y / Z_lev)**2 * exp(X0)
        END DO
      END IF
    END DO

    !--- Assemble 6x6 rate matrix ---
    ! A(I,I) = total destruction rate for level I (radiative + collisional)
    ! A(I,K) = -collisional rate I → K (off-diagonal, K ≠ I)
    ! RIGHT(I) = total formation rate for level I from continuum/H-
    DO I = 1, NLEV
      ! Start diagonal with radiative ionization rate
      A(I, I) = QRADIK(J, I)

      ! First-order T correction to recombination rate
      QRADKI(J, I) = QRADKI(J, I) + DQRAD(J, I) * DT

      ! RHS: recombination + collisional ionization + coupling to levels 7,8
      RIGHT(I) = QRADKI(J, I) + QCOLL(I, I) + QCOLL(I, 7) + QCOLL(I, 8)

      ! Add all collisional rates to diagonal (total destruction)
      DO K = 1, NATOM
        A(I, I) = A(I, I) + QCOLL(I, K)
      END DO

      ! Off-diagonal: collisional coupling between bound levels
      IF (I .LT. NLEV) THEN
        DO K = I + 1, NLEV
          A(I, K) = -QCOLL(I, K)
          A(K, I) = -QCOLL(K, I)
        END DO
      END IF
    END DO

    !--- Solve for departure coefficients ---
    CALL SOLVIT(A, NLEV, RIGHT, IPIVOT)
    DO L = 1, NLEV
      BHYD(J, L) = RIGHT(L)
    END DO

    !--- Diagnostic output (log-scaled rates) ---
    IF (IFPRNT(ITER) .GT. 1) THEN
      ! Convert rates to sign-preserving log scale for display
      DO I = 1, NLEV
        QRADKI(J, I) = sign(log10(max(abs(QRADKI(J, I) * 1.0D20), 1.0D0)), &
                             QRADKI(J, I))
        QRADIK(J, I) = sign(log10(max(abs(QRADIK(J, I) * 1.0D20), 1.0D0)), &
                             QRADIK(J, I))
      END DO
      DO I = 1, NATOM
        DO K = 1, NATOM
          QCOLL(I, K) = sign(log10(max(abs(QCOLL(I, K) * 1.0D20), 1.0D0)), &
                              QCOLL(I, K))
        END DO
      END DO

      WRITE(6, 100) J, (QRADIK(J, I), QRADKI(J, I), I=1,NLEV), &
                        (QCOLL(I, I), I=1,NLEV), QCOLL(5, 8), QCOLL(6, 8)
100   FORMAT('0', I5, 12F6.2, 6X, 8F6.2)
      WRITE(6, 110) (QCOLL(1, K), K=2,7), (QCOLL(2, K), K=3,7), &
                    (QCOLL(3, K), K=4,7), (QCOLL(4, K), K=5,7), &
                    (QCOLL(5, K), K=6,7), QCOLL(6, 7)
110   FORMAT(6X, 21F6.2)
    END IF

  END DO  ! depth loop

  !-------------------------------------------------------------------
  ! (C) Print final departure coefficients
  !-------------------------------------------------------------------
  WRITE(6, 170) (J, RHOX(J), (BHYD(J, I), I=1,NLEV), J=1,NRHOX)
170 FORMAT('1', /////, 30X, 'STATISTICAL EQUILIBRIUM FOR HYDROGEN' / &
      15X, 'RHOX', 10X, 'B1', 8X, 'B2', 8X, 'B3', 8X, 'B4', 8X, 'B5', 8X, 'B6', / &
      (8X, I2, 1PE11.4, 1X, 0P, 6F10.4))

  RETURN

END SUBROUTINE STATEQ

!=======================================================================
! SUBROUTINE RADIAP(MODE, RCOWT)
!
! Computes the radiation pressure P_rad(depth) by frequency-integrating
! the radiative acceleration (opacity * flux) and then integrating
! inward through the atmosphere via hydrostatic balance.
!
! The radiation pressure enters the equation of hydrostatic equilibrium:
!   P_gas = g * RHOX - P_rad - P_turb - P_continuous
!
! Also accumulates the mean radiation energy density (RADEN) for
! diagnostic output, and the surface K-integral (PRADK0) which
! provides the radiation pressure boundary condition.
!
! Called in three modes per iteration:
!   MODE = 1: Initialize (zero accumulators)
!   MODE = 2: Accumulate frequency integrals (called per wavelength)
!   MODE = 3: Finalize — convert units, integrate to get P_rad(depth)
!
! Physical quantities computed:
!   RADEN(J)  — mean radiation energy density (4*pi/c) * integral(J_nu * dnu)
!   ACCRAD(J) — radiative acceleration = (4*pi/c) * integral(kappa * H_nu * dnu)
!   PRAD(J)   — radiation pressure = integral(ACCRAD * dRHOX) from surface
!   PRADK0    — surface radiation pressure = (4*pi/c) * integral(K_nu * dnu)
!   PRADK(J)  — total radiation pressure = PRAD(J) + PRADK0
!=======================================================================

SUBROUTINE RADIAP(MODE, RCOWT)
  
  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN) :: MODE
  REAL(8),  INTENT(IN) :: RCOWT    ! quadrature weight for current frequency

  ! --- Named constants: FOURPI_OVER_C now from mod_constants ---

  ! --- Persistent local: frequency-integrated flux (for error scaling) ---
  REAL(8), SAVE :: H_TOTAL(kw)    ! integrated Eddington flux at each depth

  ! --- Local variables ---
  REAL(8)  :: ERRORMAX            ! max flux ratio (for stability scaling)
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING RADIAP'

  !=====================================================================
  ! MODE 1: Zero frequency-integral accumulators
  !=====================================================================
  IF (MODE .EQ. 1) THEN

    H_TOTAL = 0.0D0
    RADEN   = 0.0D0
    ACCRAD  = 0.0D0
    PRADK0 = 0.0D0
    RETURN

  !=====================================================================
  ! MODE 2: Accumulate frequency integrals at current wavelength
  !=====================================================================
  ELSE IF (MODE .EQ. 2) THEN

    ! Mean radiation energy density: J_nu * weight
    RADEN = RADEN + JNU * RCOWT

    ! Total Eddington flux: H_nu * weight
    H_TOTAL = H_TOTAL + HNU * RCOWT

    ! Radiative acceleration: kappa_nu * H_nu * weight
    ACCRAD = ACCRAD + ABTOT * HNU * RCOWT

    ! Surface K-integral (second moment, = radiation pressure at tau=0)
    PRADK0 = PRADK0 + KNU(1) * RCOWT
    RETURN

  END IF

  !=====================================================================
  ! MODE 3: Convert to physical units and integrate for P_rad
  !=====================================================================

  ! Convert from frequency sums to physical units: multiply by 4*pi/c
  ! This converts:
  !   RADEN:  sum(J_nu * dnu) → (4*pi/c) * J = radiation energy density
  !   ACCRAD: sum(kappa * H_nu * dnu) → (4*pi/c) * kappa * H = radiative acceleration
  ERRORMAX = 0.0D0
  DO J = 1, NRHOX
    RADEN(J)  = RADEN(J) * FOURPI_OVER_C
    ACCRAD(J) = ACCRAD(J) * FOURPI_OVER_C

    ! Stability safeguard: if the integrated flux exceeds the target
    ! (FLUX = sigma*Teff^4 / 4*pi), the radiative acceleration is
    ! artificially scaled down to prevent runaway radiation pressure
    ! in early iterations with large flux errors.
    IF (H_TOTAL(J) / FLUX .GT. 1.0D0) THEN
      ACCRAD(J) = ACCRAD(J) * FLUX / H_TOTAL(J)
    END IF
    ERRORMAX = max(ERRORMAX, H_TOTAL(J) / FLUX)
  END DO

  ! Same scaling for surface K-integral
  PRADK0 = PRADK0 * FOURPI_OVER_C
  IF (ERRORMAX .GT. 1.0D0) PRADK0 = PRADK0 / ERRORMAX

  ! Integrate radiative acceleration inward to get radiation pressure:
  !   P_rad(RHOX) = integral from 0 to RHOX of ACCRAD * dRHOX
  ! with boundary value ACCRAD(1)*RHOX(1) (linear extrapolation to surface)
  CALL INTEG(RHOX, ACCRAD, PRAD, NRHOX, ACCRAD(1) * RHOX(1))

  ! Total radiation pressure = depth-integrated + surface boundary term
  PRADK = PRAD + PRADK0

  RETURN

END SUBROUTINE RADIAP


!=========================================================================
! SUBROUTINE ROSS(MODE, RCOWT)
!
!   Compute the Rosseland mean opacity and optical depth scale.
!
!   The Rosseland mean is defined as the harmonic mean of opacity
!   weighted by dB/dT (the temperature derivative of the Planck function):
!
!     1/kappa_Ross = integral[ (1/kappa_nu) * (dB_nu/dT) dnu ]
!                    / integral[ (dB_nu/dT) dnu ]
!
!   The denominator equals (4*sigma/pi)*T^3 by the Stefan-Boltzmann law.
!
!   MODE = 1 : Initialize — zero the ABROSS accumulator
!   MODE = 2 : Accumulate — add (dB/dT)/kappa_nu * weight at current freq
!   MODE = 3 : Finalize   — convert sum to kappa_Ross and integrate
!              the Rosseland optical depth scale tau_Ross(RHOX)
!
!   RCOWT = frequency integration weight for current wavelength point
!=========================================================================

SUBROUTINE ROSS(MODE, RCOWT)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: MODE
  REAL(8),  INTENT(IN) :: RCOWT

  ! --- Local variables ---
  REAL(8)  :: DBDT
  INTEGER :: J

  ! FOUR_SIGMA_OVER_PI now from mod_constants

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING ROSS'
  SELECT CASE (MODE)

  ! --- MODE 1: Zero the Rosseland mean accumulator ---
  CASE (1)
    ABROSS = 0.0D0

  ! --- MODE 2: Accumulate (dB_nu/dT) / kappa_nu at current frequency ---
  CASE (2)
    DO J = 1, NRHOX
      ! dB_nu/dT = B_nu * (h*nu/kT) / T / (1 - e^{-h*nu/kT})
      !          = BNU * FREQ * HKT / T / STIM
      DBDT = BNU(J) * FREQ * HKT(J) / T(J) / STIM(J)
      ! Single-frequency fallback: use total dB/dT = (4*sigma/pi)*T^3
      IF (NUMNU .EQ. 1) DBDT = FOUR_SIGMA_OVER_PI * T(J)**3
      ABROSS(J) = ABROSS(J) + DBDT / ABTOT(J) * RCOWT
    END DO

  ! --- MODE 3: Finalize — compute kappa_Ross and tau_Ross scale ---
  CASE (3)
    ! Convert accumulated sum to Rosseland mean opacity:
    !   kappa_Ross = (4*sigma/pi)*T^3 / sum[ (dB/dT)/kappa_nu * dnu ]
    ABROSS = FOUR_SIGMA_OVER_PI * T**3 / ABROSS
    ! Integrate kappa_Ross over column mass to get Rosseland optical depth
    CALL INTEG(RHOX, ABROSS, TAUROS, NRHOX, ABROSS(1) * RHOX(1))

  END SELECT

END SUBROUTINE ROSS


!=========================================================================
! SUBROUTINE DERIV(X, F, DFDX, N)
!
!   Compute the numerical derivative dF/dX at N tabulated points using
!   a smooth tangent-averaging scheme at interior points.
!
!   At endpoints (J=1, J=N): simple one-sided finite differences.
!
!   At interior points (J=2..N-1): the left and right finite-difference
!   slopes are converted to half-angles via
!       t = slope / (S * sqrt(1 + slope^2) + 1)
!   and combined using the tangent addition formula
!       dF/dX = (t_right + t_left) / (1 - t_right * t_left) * SCALE
!   This produces a smoothly varying derivative that avoids the
!   overshooting of simple centered differences near sharp features.
!   S = sign(deltaX) preserves orientation for decreasing X grids.
!
!   The slopes are normalized by SCALE = max(|F|)/|X| at each point
!   to prevent overflow when F and X have very different magnitudes.
!
!   Assumes that any zero in X occurs only at an endpoint.
!=========================================================================

SUBROUTINE DERIV(X, F, DFDX, N)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: N
  REAL(8),  INTENT(IN)  :: X(*), F(*)
  REAL(8),  INTENT(OUT) :: DFDX(*)

  ! --- Local variables ---
  REAL(8)  :: S, SCALE, D_RIGHT, D_LEFT, T_RIGHT, T_LEFT
  INTEGER :: J

  ! --- Endpoint derivatives: one-sided finite differences ---
  DFDX(1) = (F(2) - F(1)) / (X(2) - X(1))
  DFDX(N) = (F(N) - F(N-1)) / (X(N) - X(N-1))
  IF (N .EQ. 2) RETURN

  ! --- Sign of the X spacing (handles decreasing grids) ---
  S = abs(X(2) - X(1)) / (X(2) - X(1))

  ! --- Interior points: tangent half-angle averaging ---
  DO J = 2, N - 1
    ! Normalization scale to prevent overflow
    SCALE = max(abs(F(J-1)), abs(F(J)), abs(F(J+1))) / abs(X(J))
    IF (SCALE .EQ. 0.0D0) SCALE = 1.0D0

    ! Scaled right and left finite-difference slopes
    D_RIGHT = (F(J+1) - F(J)) / (X(J+1) - X(J)) / SCALE
    D_LEFT  = (F(J) - F(J-1)) / (X(J) - X(J-1)) / SCALE

    ! Convert slopes to half-angle tangents: t = d / (S*sqrt(1+d^2) + 1)
    T_RIGHT = D_RIGHT / (S * sqrt(1.0D0 + D_RIGHT**2) + 1.0D0)
    T_LEFT  = D_LEFT  / (S * sqrt(1.0D0 + D_LEFT**2)  + 1.0D0)

    ! Tangent addition formula: tan(a+b) = (tan(a)+tan(b))/(1-tan(a)*tan(b))
    DFDX(J) = (T_RIGHT + T_LEFT) / (1.0D0 - T_RIGHT * T_LEFT) * SCALE
  END DO

END SUBROUTINE DERIV

!=========================================================================
! SUBROUTINE INTEG(X, F, FINT, N, START)
!
!   Integrate F(X) over the tabulated grid X(1..N) using piecewise
!   parabolic interpolation, returning the running integral FINT(1..N).
!
!   The function F is first fit with piecewise parabolas via PARCOE,
!   which returns coefficients A, B, C such that in each interval
!   [X(I), X(I+1)]:
!       F(x) ≈ A(I) + B(I)*x + C(I)*x^2
!
!   The integral over each interval is then computed analytically:
!       integral_{X(I)}^{X(I+1)} (A + Bx + Cx^2) dx
!         = [ A*x + B*x^2/2 + C*x^3/3 ]_{X(I)}^{X(I+1)}
!         = (A + B/2*(X_{I+1}+X_I) + C/3*(X_{I+1}^2+X_{I+1}*X_I+X_I^2))
!           * (X_{I+1} - X_I)
!
!   Arguments:
!     X(N)     — abscissa grid (e.g., column mass RHOX)
!     F(N)     — integrand values at grid points (e.g., kappa_Ross)
!     FINT(N)  — output: running integral; FINT(1) = START
!     N        — number of grid points
!     START    — initial value for FINT(1) (boundary condition)
!=========================================================================

SUBROUTINE INTEG(X, F, FINT, N, START)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: N
  REAL(8),  INTENT(IN)  :: X(*), F(*), START
  REAL(8),  INTENT(OUT) :: FINT(*)

  ! --- Local variables ---
  REAL(8)  :: A(kw), B(kw), C(kw)
  REAL(8)  :: XLO, XHI
  INTEGER :: I

  ! Fit piecewise parabolas to F(X)
  CALL PARCOE(F, X, A, B, C, N)

  ! Integrate by summing the analytic parabolic integrals
  FINT(1) = START
  DO I = 1, N - 1
    XLO = X(I)
    XHI = X(I+1)
    FINT(I+1) = FINT(I) + (A(I) + B(I) * 0.5D0 * (XHI + XLO) &
                + C(I) / 3.0D0 * ((XHI + XLO) * XHI + XLO * XLO)) &
                * (XHI - XLO)
  END DO

END SUBROUTINE INTEG

!=========================================================================
! SUBROUTINE PARCOE(F, X, A, B, C, N)
!
!   Compute piecewise parabolic interpolation coefficients for a
!   tabulated function F(X) on N points.
!
!   Returns coefficients A(I), B(I), C(I) such that in each interval
!   [X(I), X(I+1)]:
!       F(x) ≈ A(I) + B(I)*x + C(I)*x^2
!
!   These coefficients are designed for use with INTEG, which
!   analytically integrates the parabolic segments.
!
!   Algorithm:
!     1. Endpoints (I=1, I=N): linear interpolation (C=0)
!     2. Interior points (N>=4): fit parabola through three consecutive
!        points using divided differences for the curvature C(J)
!     3. Override intervals 2 and 3 with linear fits so that the
!        smoothing pass (step 4) starts from linear near the boundary
!     4. Smoothing pass: blend each interval's coefficients with its
!        right neighbor using curvature-weighted averaging. The weight
!        WT = |C_{J+1}|/(|C_{J+1}|+|C_J|) favors the less curved
!        (smoother) side, suppressing parabolic overshooting.
!     5. Penultimate interval (N-1): inherits endpoint linear fit
!
!   Special cases:
!     N = 2 : Single linear segment (returned after step 1)
!     N = 3 : Two linear segments. The post-loop smoothing code would
!             access out-of-bounds X(4)/F(4), so we return early.
!             (Bug fix: John Lester, 18 May 2007)
!=========================================================================

SUBROUTINE PARCOE(F, X, A, B, C, N)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: N
  REAL(8),  INTENT(IN)  :: F(*), X(*)
  REAL(8),  INTENT(OUT) :: A(*), B(*), C(*)

  ! --- Local variables ---
  REAL(8)  :: D, WT
  INTEGER :: J

  ! --- Step 1: Endpoint linear fits (C = 0) ---
  C(1) = 0.0D0
  B(1) = (F(2) - F(1)) / (X(2) - X(1))
  A(1) = F(1) - X(1) * B(1)

  C(N) = 0.0D0
  B(N) = (F(N) - F(N-1)) / (X(N) - X(N-1))
  A(N) = F(N) - X(N) * B(N)

  IF (N .EQ. 2) RETURN

  ! --- N = 3: use linear interpolation on both intervals ---
  !     Copy last-interval linear fit to the middle point.
  !     (Returning early prevents the post-loop code from accessing
  !     out-of-bounds F(4) and X(4).)
  IF (N .EQ. 3) THEN
    A(2) = A(3)
    B(2) = B(3)
    C(2) = C(3)
    RETURN
  END IF

  ! --- Step 2: Parabolic fits at interior points (N >= 4) ---
  !     For each interior point J, fit a parabola through
  !     (X_{J-1}, F_{J-1}), (X_J, F_J), (X_{J+1}, F_{J+1}).
  !     C(J) is the second divided difference (curvature term).
  DO J = 2, N - 1
    D = (F(J) - F(J-1)) / (X(J) - X(J-1))
    C(J) = F(J+1) / ((X(J+1) - X(J)) * (X(J+1) - X(J-1))) &
         - F(J)   / ((X(J) - X(J-1))   * (X(J+1) - X(J)))   &
         + F(J-1) / ((X(J) - X(J-1))   * (X(J+1) - X(J-1)))
    B(J) = D - (X(J) + X(J-1)) * C(J)
    A(J) = F(J-1) - X(J-1) * D + X(J) * X(J-1) * C(J)
  END DO

  ! --- Step 3: Override first two interior points with linear fits ---
  !     These are replaced before the smoothing pass so that the
  !     blending near the boundaries starts from linear.
  C(2) = 0.0D0
  B(2) = (F(3) - F(2)) / (X(3) - X(2))
  A(2) = F(2) - X(2) * B(2)

  C(3) = 0.0D0
  B(3) = (F(4) - F(3)) / (X(4) - X(3))
  A(3) = F(3) - X(3) * B(3)

  ! --- Step 4: Smoothing pass — curvature-weighted blending ---
  !     Where curvature is non-zero, blend with the right neighbor's
  !     coefficients. The weight WT = |C_{J+1}| / (|C_{J+1}| + |C_J|)
  !     favors the smoother (lower curvature) side.
  DO J = 2, N - 1
    IF (C(J) .NE. 0.0D0) THEN
      WT = abs(C(J+1)) / (abs(C(J+1)) + abs(C(J)))
      A(J) = A(J+1) + WT * (A(J) - A(J+1))
      B(J) = B(J+1) + WT * (B(J) - B(J+1))
      C(J) = C(J+1) + WT * (C(J) - C(J+1))
    END IF
  END DO

  ! --- Step 5: Penultimate point inherits endpoint linear fit ---
  A(N-1) = A(N)
  B(N-1) = B(N)
  C(N-1) = C(N)

END SUBROUTINE PARCOE

!=========================================================================
! FUNCTION MAP1(XOLD, FOLD, NOLD, XNEW, FNEW, NNEW)
!
!   Remap a tabulated function from one grid to another using piecewise
!   parabolic interpolation with curvature-weighted blending.
!
!   Given FOLD(1..NOLD) defined on XOLD(1..NOLD), evaluate the function
!   at the new grid points XNEW(1..NNEW), storing results in FNEW.
!   Both grids must be monotonically increasing.
!
!   At interior points (L = 4..NOLD-1), two parabolas are computed:
!     "backward" through (XOLD(L-2), XOLD(L-1), XOLD(L))
!     "forward"  through (XOLD(L-1), XOLD(L),   XOLD(L+1))
!   and blended using WT = |C_fwd| / (|C_fwd| + |C_bwd|), which
!   favors the less curved fit (same technique as PARCOE).
!
!   At boundaries (L = 2,3 or L = NOLD), linear or single-parabola
!   interpolation is used instead.
!
!   Returns MAP1 = index of the last old-grid interval used.
!=========================================================================

FUNCTION MAP1(XOLD, FOLD, NOLD, XNEW, FNEW, NNEW)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NOLD, NNEW
  REAL(8),  INTENT(IN)  :: XOLD(*), FOLD(*), XNEW(*)
  REAL(8),  INTENT(OUT) :: FNEW(*)
  INTEGER              :: MAP1

  ! --- Local variables ---
  REAL(8)  :: A, B, C                   ! Current interpolation coefficients
  REAL(8)  :: ABAC, BBAC, CBAC          ! Backward parabola coefficients
  REAL(8)  :: AFOR, BFOR, CFOR          ! Forward parabola coefficients
  REAL(8)  :: D, WT
  INTEGER :: K, L, LL, L1, L2
  LOGICAL :: past_right_edge           ! True if XNEW(K) is beyond XOLD(NOLD)

  L  = 2     ! Pointer into XOLD: XNEW(K) is in [XOLD(L-1), XOLD(L)]
  LL = 0     ! Previous L value (0 = no coefficients computed yet)

  DO K = 1, NNEW

    ! -----------------------------------------------------------------
    ! Step 1: Advance L to bracket XNEW(K) in [XOLD(L-1), XOLD(L)]
    ! -----------------------------------------------------------------
    past_right_edge = .FALSE.
    DO WHILE (XNEW(K) .GE. XOLD(L))
      L = L + 1
      IF (L .GT. NOLD) THEN
        ! Beyond right edge of old grid — use boundary linear fit
        past_right_edge = .TRUE.
        EXIT
      END IF
    END DO

    ! -----------------------------------------------------------------
    ! Step 2: If still in the same interval, reuse existing coefficients
    ! -----------------------------------------------------------------
    IF (.NOT. past_right_edge .AND. L .EQ. LL) THEN
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      CYCLE
    END IF

    ! -----------------------------------------------------------------
    ! Step 3: Boundary cases — use linear interpolation
    !         (L=2,3: too few points to the left for a parabola;
    !          past right edge: linear extrapolation from last two points)
    ! -----------------------------------------------------------------
    IF (L .LE. 3 .OR. past_right_edge) THEN
      L = min(NOLD, L)
      IF (.NOT. past_right_edge .AND. L .EQ. LL) THEN
        FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
        CYCLE
      END IF
      C = 0.0D0
      B = (FOLD(L) - FOLD(L-1)) / (XOLD(L) - XOLD(L-1))
      A = FOLD(L) - XOLD(L) * B
      LL = L
      FNEW(K) = A + B * XNEW(K)
      CYCLE
    END IF

    ! -----------------------------------------------------------------
    ! Step 4: Interior — compute backward parabola through (L-2, L-1, L)
    ! -----------------------------------------------------------------
    L1 = L - 1

    ! Check if we can reuse the previous forward parabola as backward
    IF (L .EQ. LL + 1 .AND. L .NE. 3 .AND. L .NE. 4) THEN
      CBAC = CFOR
      BBAC = BFOR
      ABAC = AFOR
    ELSE
      ! Compute fresh backward parabola
      L2 = L - 2
      D = (FOLD(L1) - FOLD(L2)) / (XOLD(L1) - XOLD(L2))
      CBAC = FOLD(L)  / ((XOLD(L) - XOLD(L1)) * (XOLD(L) - XOLD(L2))) &
           + (FOLD(L2) / (XOLD(L) - XOLD(L2)) &
           -  FOLD(L1) / (XOLD(L) - XOLD(L1))) &
           / (XOLD(L1) - XOLD(L2))
      BBAC = D - (XOLD(L1) + XOLD(L2)) * CBAC
      ABAC = FOLD(L2) - XOLD(L2) * D + XOLD(L1) * XOLD(L2) * CBAC
    END IF

    ! -----------------------------------------------------------------
    ! Step 5: At right edge (L=NOLD) — use backward parabola only
    ! -----------------------------------------------------------------
    IF (L .GE. NOLD) THEN
      A = ABAC;  B = BBAC;  C = CBAC
      LL = L
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      CYCLE
    END IF

    ! -----------------------------------------------------------------
    ! Step 6: Compute forward parabola through (L-1, L, L+1) and blend
    ! -----------------------------------------------------------------
    D = (FOLD(L) - FOLD(L1)) / (XOLD(L) - XOLD(L1))
    CFOR = FOLD(L+1) / ((XOLD(L+1) - XOLD(L)) * (XOLD(L+1) - XOLD(L1))) &
         + (FOLD(L1) / (XOLD(L+1) - XOLD(L1)) &
         -  FOLD(L)  / (XOLD(L+1) - XOLD(L))) &
         / (XOLD(L) - XOLD(L1))
    BFOR = D - (XOLD(L) + XOLD(L1)) * CFOR
    AFOR = FOLD(L1) - XOLD(L1) * D + XOLD(L) * XOLD(L1) * CFOR

    ! Curvature-weighted blend: favor the less curved parabola
    WT = 0.0D0
    IF (abs(CFOR) .NE. 0.0D0) WT = abs(CFOR) / (abs(CFOR) + abs(CBAC))
    A = AFOR + WT * (ABAC - AFOR)
    B = BFOR + WT * (BBAC - BFOR)
    C = CFOR + WT * (CBAC - CFOR)
    LL = L

    ! -----------------------------------------------------------------
    ! Evaluate interpolant at XNEW(K)
    ! -----------------------------------------------------------------
    FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)

  END DO

  MAP1 = LL - 1

END FUNCTION MAP1

!=========================================================================
! SUBROUTINE SOLVIT(A, N, B, IPIVOT)
!
!   Solve the N×N linear system  A * x = B  by Gaussian elimination
!   with partial (row) pivoting.  On return, B contains the solution x.
!   A is destroyed (overwritten with the LU factors).
!
!   The algorithm has three phases:
!
!   Phase 1 — LU factorization with partial pivoting (in-place):
!     For each column I = 1..N-1:
!       (a) Find the pivot row M (max |A(K,I)| for K = I..N)
!       (b) Swap rows I and M for columns I+1..N only (column I is
!           handled implicitly — an optimization that avoids N swaps
!           per column while preserving correctness)
!       (c) Store 1/pivot on the diagonal: A(I,I) = 1/A(M,I)
!       (d) Compute multipliers: A(K,I) = A(K,I) / pivot, K = I+1..N
!       (e) Eliminate: A(K,J) -= A(K,I)*A(I,J),  K,J = I+1..N
!     Final: A(N,N) = 1/A(N,N)
!
!   Phase 2 — Forward substitution (apply L^{-1} and row swaps to B):
!     For I = 1..N-1: swap B(I) ↔ B(IPIVOT(I)), then
!       B(K) -= A(K,I)*B(I) for K = I+1..N
!
!   Phase 3 — Back substitution (apply U^{-1} to B):
!     For J = N..2: B(J) *= A(J,J) [= B(J)/U(J,J)], then
!       B(K) -= A(K,J)*B(J) for K = 1..J-1
!     Final: B(1) *= A(1,1)
!
!   Arguments:
!     A(N,N)    — coefficient matrix (destroyed on output)
!     N         — system size
!     B(N)      — right-hand side on input, solution on output
!     IPIVOT(N) — integer scratch array for pivot row indices
!
!   Notes:
!     - Near-singular protection: if a pivot magnitude falls below
!       PIVOT_TOL (1e-30), it is replaced with ±PIVOT_TOL and a
!       warning is printed to unit 6. This prevents division by zero
!       while producing a least-damaging approximate solution.
!     - Callers in ATLAS12 (STATEQ, NMOLEC) pass REAL*8 workspace for
!       IPIVOT. This works through sequence association since SOLVIT
!       treats it as INTEGER and the caller never reads it back.
!=========================================================================

SUBROUTINE SOLVIT(A, N, B, IPIVOT)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: N
  REAL(8),  INTENT(INOUT) :: A(N,N), B(N)
  INTEGER, INTENT(OUT)   :: IPIVOT(N)

  ! --- Local variables ---
  REAL(8)  :: PIVOT, T, C
  INTEGER :: I, J, K, M

  ! Threshold for treating a pivot as effectively zero.
  ! Machine epsilon for REAL*8 is ~2.2e-16; we use a generous threshold
  ! relative to the column's maximum element to catch near-singular cases.
  REAL(8), PARAMETER :: PIVOT_TOL = 1.0D-30

  ! =================================================================
  ! Phase 1: LU factorization with partial pivoting
  ! =================================================================

  DO I = 1, N - 1

    ! --- (a) Find pivot: row with largest |A(K,I)| for K >= I ---
    M = I
    DO K = I + 1, N
      IF (abs(A(K,I)) .GT. abs(A(M,I))) M = K
    END DO
    IPIVOT(I) = M

    ! --- (b) Swap rows I and M for columns I+1..N ---
    !     Column I is not swapped; its entries are handled implicitly
    !     through the multiplier computation below.
    IF (M .NE. I) THEN
      DO K = I + 1, N
        T = A(I,K)
        A(I,K) = A(M,K)
        A(M,K) = T
      END DO
    END IF

    ! --- (c) Store reciprocal pivot on diagonal ---
    IF (abs(A(M,I)) .LT. PIVOT_TOL) THEN
      WRITE(6,'(A,I4,A,1PE12.4)') &
        ' SOLVIT: near-singular matrix at column ', I, &
        ', pivot = ', A(M,I)
      A(M,I) = sign(PIVOT_TOL, A(M,I))
    END IF
    PIVOT = 1.0D0 / A(M,I)
    A(M,I) = A(I,I)
    A(I,I) = PIVOT

    ! --- (d) Compute multipliers L(K,I) = A(K,I) / pivot ---
    DO K = I + 1, N
      A(K,I) = A(K,I) * PIVOT
    END DO

    ! --- (e) Eliminate: update submatrix A(I+1:N, I+1:N) ---
    DO J = I + 1, N
      C = A(I,J)
      IF (C .EQ. 0.0D0) CYCLE
      DO K = I + 1, N
        A(K,J) = A(K,J) - A(K,I) * C
      END DO
    END DO

  END DO

  ! Reciprocal of the last pivot
  IF (abs(A(N,N)) .LT. PIVOT_TOL) THEN
    WRITE(6,'(A,I4,A,1PE12.4)') &
      ' SOLVIT: near-singular matrix at column ', N, &
      ', pivot = ', A(N,N)
    A(N,N) = sign(PIVOT_TOL, A(N,N))
  END IF
  A(N,N) = 1.0D0 / A(N,N)

  ! =================================================================
  ! Phase 2: Forward substitution — apply row swaps and L^{-1} to B
  ! =================================================================

  DO I = 1, N - 1
    M = IPIVOT(I)
    IF (M .NE. I) THEN
      T = B(M)
      B(M) = B(I)
      B(I) = T
    END IF
    C = B(I)
    DO K = I + 1, N
      B(K) = B(K) - A(K,I) * C
    END DO
  END DO

  ! =================================================================
  ! Phase 3: Back substitution — apply U^{-1} to B
  ! =================================================================
  !   A(J,J) stores 1/U(J,J); A(K,J) for K<J stores U(K,J).

  DO J = N, 2, -1
    B(J) = B(J) * A(J,J)
    C = B(J)
    DO K = 1, J - 1
      B(K) = B(K) - A(K,J) * C
    END DO
  END DO
  B(1) = B(1) * A(1,1)

END SUBROUTINE SOLVIT

!=========================================================================
! SUBROUTINE DUMP_ARRAY(LABEL, ARR, N)
!
!   Print a labeled array to standard output (unit 6) for diagnostics.
!   Used for error dumps (e.g., negative pressure in HSE) and iteration
!   summaries.
!
!   Prints LABEL followed by ARR(1..N) in exponential format, 10 values
!   per line.  For scalar diagnostics, pass the scalar directly with N=1
!   (Fortran passes by reference, so W reads it as ARR(1)).
!=========================================================================

SUBROUTINE DUMP_ARRAY(LABEL, ARR, N)

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: LABEL
  REAL(8),       INTENT(IN) :: ARR(*)
  INTEGER,      INTENT(IN) :: N

  INTEGER :: I

  WRITE(6, '(A6,1P10E12.4 / (7X,10E12.4))') LABEL, (ARR(I), I=1,N)

END SUBROUTINE DUMP_ARRAY

!=========================================================================
! FUNCTION ROSSTAB(TEMP, PRESSURE, V)
!
!   Rosseland opacity lookup table: stores and interpolates kappa_Ross
!   as a function of temperature and gas pressure.
!
!   Each time the Rosseland mean opacity is computed (after ROSS mode 3),
!   the current model's kappa_Ross profile is appended to a growing table
!   of (T, P, kappa_Ross) points. Subsequent calls can then interpolate
!   kappa_Ross at any (T, P) — needed when temperature corrections shift
!   the structure away from the original grid points.
!
!   Usage:
!     ROSSTAB(0, 0, 0)       — STORE: append current ABROSS(1..NRHOX)
!     ROSSTAB(T, P, Vturb)   — LOOKUP: interpolate log(kappa_Ross)
!
!   Interpolation mode controlled by module variable IROSSTAB:
!     1 = Bilinear: nearest neighbor in each of 4 quadrants,
!         bilinear interpolation if all populated, inverse-distance
!         weighted fallback otherwise.
!     2 = Shepard: K=12 nearest neighbors, weighted average with
!         w = 1/(d^2 + eps)^(p/2), p=3. Smooth (C1) interpolation.
!=========================================================================

FUNCTION ROSSTAB(TEMP, PRESSURE, V)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: TEMP, PRESSURE, V
  REAL(8)             :: ROSSTAB

  ! --- Table storage (persistent across calls) ---
  INTEGER, PARAMETER :: MAXTAB = kw * 60

  REAL(8),  SAVE :: ROSS(MAXTAB)     ! log10(kappa_Ross)
  REAL(8),  SAVE :: TABT(MAXTAB)     ! normalized log10(T)
  REAL(8),  SAVE :: TABP(MAXTAB)     ! normalized log10(P)
  REAL(8),  SAVE :: ZEROT, ZEROP     ! normalization offsets
  REAL(8),  SAVE :: SLOPET, SLOPEP   ! normalization ranges
  INTEGER, SAVE :: NROSS = 0        ! number of entries in table

  ! --- Shepard parameters ---
  INTEGER, PARAMETER :: KFIT = 12          ! nearest neighbors
  REAL(8),  PARAMETER :: SHEP_PHALF = 1.5D0 ! p/2 where p=3
  REAL(8),  PARAMETER :: SHEP_EPS = 1.0D-12 ! softening

  ! --- Local variables ---
  REAL(8)  :: TEMPLOG, PRESSLOG      ! normalized query coordinates
  REAL(8)  :: DT, DP, RADIUS2        ! displacement and distance
  REAL(8)  :: R, RWT, W              ! interpolated result, weights

  ! Bilinear variables
  REAL(8)  :: RMIN_PP, RMIN_PM, RMIN_MP, RMIN_MM
  INTEGER :: IDX_PP,  IDX_PM,  IDX_MP,  IDX_MM
  REAL(8)  :: DIST_PP, DIST_PM, DIST_MP, DIST_MM
  REAL(8)  :: TPP, PPP, RPP, TPM, PPM, RPM
  REAL(8)  :: TMP, PMP, RMP, TMM, PMM, RMM
  REAL(8)  :: R_HIPRESS, R_LOPRESS
  REAL(8)  :: P_HIPRESS, P_LOPRESS

  ! Shepard variables
  REAL(8)  :: DIST2(KFIT)
  INTEGER :: KIDX(KFIT)
  REAL(8)  :: DMAX2
  INTEGER :: KMAX, KFOUND

  INTEGER :: I, J, K

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING ROSSTAB'

  ! =================================================================
  ! MODE 1: STORE — append current kappa_Ross profile to table
  ! =================================================================

  IF (TEMP .LE. 0.0D0) THEN

    IF (NROSS .EQ. 0) THEN
      ZEROT  = log10(T(1))
      ZEROP  = log10(P(1))
      SLOPET = log10(T(NRHOX)) - ZEROT
      SLOPEP = log10(P(NRHOX)) - ZEROP
    END IF

    DO J = 1, NRHOX
      NROSS = NROSS + 1
      IF (NROSS .GT. MAXTAB) THEN
        WRITE(6,*) 'ROSSTAB ERROR: table overflow, NROSS =', NROSS
        NROSS = MAXTAB
        EXIT
      END IF
      TABT(NROSS) = (log10(T(J)) - ZEROT) / SLOPET
      TABP(NROSS) = (log10(P(J)) - ZEROP) / SLOPEP
      ROSS(NROSS) = log10(ABROSS(J))
      IF (IDEBUG .EQ. 1) &
        WRITE(6, '(" ROSSTAB",I5,F10.1,F10.5,1PE12.3,0PF10.5,F10.5)') &
        NROSS, T(J), TABT(NROSS), P(J), TABP(NROSS), ROSS(NROSS)
    END DO

    ROSSTAB = 0.0D0
    RETURN
  END IF

  ! =================================================================
  ! MODE 2: LOOKUP
  ! =================================================================

  TEMPLOG  = (log10(TEMP) - ZEROT) / SLOPET
  PRESSLOG = (log10(PRESSURE) - ZEROP) / SLOPEP

  ! -----------------------------------------------------------------
  ! IROSSTAB = 1: Original bilinear (4-quadrant nearest neighbor)
  ! -----------------------------------------------------------------
  IF (IROSSTAB .EQ. 1) THEN

    RMIN_PP = 1.0D30;  IDX_PP = 0
    RMIN_PM = 1.0D30;  IDX_PM = 0
    RMIN_MP = 1.0D30;  IDX_MP = 0
    RMIN_MM = 1.0D30;  IDX_MM = 0

    DO I = 1, NROSS
      DT = TABT(I) - TEMPLOG
      DP = TABP(I) - PRESSLOG
      RADIUS2 = DT**2 + DP**2

      IF (DT .GE. 0.0D0) THEN
        IF (DP .GE. 0.0D0) THEN
          IF (RADIUS2 .LT. RMIN_PP) THEN
            RMIN_PP = RADIUS2;  IDX_PP = I
          END IF
        ELSE
          IF (RADIUS2 .LT. RMIN_PM) THEN
            RMIN_PM = RADIUS2;  IDX_PM = I
          END IF
        END IF
      ELSE
        IF (DP .GE. 0.0D0) THEN
          IF (RADIUS2 .LT. RMIN_MP) THEN
            RMIN_MP = RADIUS2;  IDX_MP = I
          END IF
        ELSE
          IF (RADIUS2 .LT. RMIN_MM) THEN
            RMIN_MM = RADIUS2;  IDX_MM = I
          END IF
        END IF
      END IF
    END DO

    IF (IDX_PP .NE. 0 .AND. IDX_PM .NE. 0 .AND. &
        IDX_MP .NE. 0 .AND. IDX_MM .NE. 0) THEN

      TPP = TABT(IDX_PP);  PPP = TABP(IDX_PP);  RPP = ROSS(IDX_PP)
      TPM = TABT(IDX_PM);  PPM = TABP(IDX_PM);  RPM = ROSS(IDX_PM)
      TMP = TABT(IDX_MP);  PMP = TABP(IDX_MP);  RMP = ROSS(IDX_MP)
      TMM = TABT(IDX_MM);  PMM = TABP(IDX_MM);  RMM = ROSS(IDX_MM)

      R_HIPRESS = ((TEMPLOG - TMP) * RPP &
        + (TPP - TEMPLOG) * RMP) / (TPP - TMP)
      P_HIPRESS = ((TEMPLOG - TMP) * PPP &
        + (TPP - TEMPLOG) * PMP) / (TPP - TMP)

      R_LOPRESS = ((TEMPLOG - TMM) * RPM &
        + (TPM - TEMPLOG) * RMM) / (TPM - TMM)
      P_LOPRESS = ((TEMPLOG - TMM) * PPM &
        + (TPM - TEMPLOG) * PMM) / (TPM - TMM)

      R = ((PRESSLOG - P_LOPRESS) * R_HIPRESS &
        + (P_HIPRESS - PRESSLOG) * R_LOPRESS) &
        / (P_HIPRESS - P_LOPRESS)

    ELSE

      DIST_PP = sqrt(RMIN_PP) + 1.0D-5
      DIST_PM = sqrt(RMIN_PM) + 1.0D-5
      DIST_MP = sqrt(RMIN_MP) + 1.0D-5
      DIST_MM = sqrt(RMIN_MM) + 1.0D-5

      RWT = 1.0D0/DIST_PP + 1.0D0/DIST_PM &
          + 1.0D0/DIST_MP + 1.0D0/DIST_MM

      R = (ROSS(max(1, IDX_PP)) / DIST_PP &
         + ROSS(max(1, IDX_PM)) / DIST_PM &
         + ROSS(max(1, IDX_MP)) / DIST_MP &
         + ROSS(max(1, IDX_MM)) / DIST_MM) / RWT

    END IF

  ! -----------------------------------------------------------------
  ! IROSSTAB = 2: Shepard (K-nearest, p=3)
  ! -----------------------------------------------------------------
  ELSE

    KFOUND = 0
    DMAX2  = 0.0D0
    KMAX   = 1

    DO I = 1, NROSS
      DT = TABT(I) - TEMPLOG
      DP = TABP(I) - PRESSLOG
      RADIUS2 = DT * DT + DP * DP

      IF (KFOUND .LT. KFIT) THEN
        KFOUND = KFOUND + 1
        DIST2(KFOUND) = RADIUS2
        KIDX(KFOUND)  = I
        IF (RADIUS2 .GT. DMAX2) THEN
          DMAX2 = RADIUS2
          KMAX  = KFOUND
        END IF
      ELSE IF (RADIUS2 .LT. DMAX2) THEN
        DIST2(KMAX) = RADIUS2
        KIDX(KMAX)  = I
        DMAX2 = DIST2(1)
        KMAX  = 1
        DO K = 2, KFIT
          IF (DIST2(K) .GT. DMAX2) THEN
            DMAX2 = DIST2(K)
            KMAX  = K
          END IF
        END DO
      END IF
    END DO

    R   = 0.0D0
    RWT = 0.0D0
    DO K = 1, KFOUND
      W   = 1.0D0 / (DIST2(K) + SHEP_EPS) ** SHEP_PHALF
      R   = R   + W * ROSS(KIDX(K))
      RWT = RWT + W
    END DO
    R = R / RWT

  END IF

  ROSSTAB = 10.0D0 ** R

END FUNCTION ROSSTAB

!=========================================================================
! SUBROUTINE TTAUP(T_in, TAU, ABSTD, PTOTAL, P_out, PRAD_in,
!                  PTURB_in, VTURB_in, GRAV_in, NUMTAU)
!
! Compute total pressure from T(tau) by integrating hydrostatic
! equilibrium on a log-spaced optical depth grid.
!
! Given a temperature profile T(tau) and a Rosseland opacity lookup
! table (ROSSTAB), this routine iteratively solves for the total
! pressure at each depth point by integrating:
!
!   d(ln P_total) / d(ln tau) = g / (kappa_Ross * P_total) * tau
!
! which is the hydrostatic equilibrium equation in log-tau form:
!   dP/dRHOX = g  and  dTAU = kappa * dRHOX
!   → dP = g * dTAU / kappa  → d(ln P) = g*tau/(kappa*P) * d(ln tau)
!
! The integration uses Adams-Bashforth multistep predictors:
!   J = 1:      direct from surface: P_total = g * tau / kappa
!   J = 2-4:    linear extrapolation from previous point
!   J > 4:      4th-order Adams-Bashforth predictor
! followed by corrector iterations (predictor-corrector) until
! convergence to 5e-5 in ln(P) or MAX_CORRECTOR iterations.
!
! At each depth, the gas pressure is extracted from total pressure:
!   P_gas = P_total + (P_rad(1) - P_rad(J)) + (P_turb(1) - P_turb(J))
! which accounts for the fact that radiation and turbulent pressures
! are defined relative to their surface values.
!
! The Rosseland opacity kappa_Ross(T, P_gas) is obtained from the
! interpolation table ROSSTAB, which accumulates (T, P, kappa) data
! from previous iterations.
!
! Stability features:
!   - Surface bootstrap: re-solves J=1 if initial opacity guess was
!     off by more than a factor of 2
!   - Under-relaxation: detects corrector oscillations (error not
!     decreasing) and reduces step size to ensure convergence
!   - DPLOG clamping: d(ln P)/d(ln tau) clamped to [0, DPLOG_MAX]
!     to prevent unphysical pressure jumps
!   - Convergence failure warning: diagnostic printed if corrector
!     hits iteration limit
!
! Requirements:
!   - TAU array must have uniform log spacing: TAU(J+1)/TAU(J) = const
!   - ROSSTAB table must be populated before calling
!
! Arguments:
!   T_in(NUMTAU)     — temperature at each depth (input)
!   TAU(NUMTAU)      — Rosseland optical depth grid (input, log-spaced)
!   ABSTD(NUMTAU)    — Rosseland opacity at each depth (output)
!   PTOTAL(NUMTAU)   — total pressure at each depth (output)
!   P_out(NUMTAU)    — gas pressure at each depth (output)
!   PRAD_in(NUMTAU)  — radiation pressure at each depth (input)
!   PTURB_in(NUMTAU) — turbulent pressure at each depth (input)
!   VTURB_in(NUMTAU) — microturbulent velocity at each depth (input)
!   GRAV_in          — surface gravity (input)
!   NUMTAU           — number of depth points (input)
!=========================================================================

SUBROUTINE TTAUP(T_in, TAU, ABSTD, PTOTAL, P_out, PRAD_in, PTURB_in, &
                 VTURB_in, GRAV_in, NUMTAU)

  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN)    :: NUMTAU
  REAL(8),  INTENT(IN)    :: T_in(NUMTAU)
  REAL(8),  INTENT(IN)    :: TAU(NUMTAU)
  REAL(8),  INTENT(OUT)   :: ABSTD(NUMTAU)
  REAL(8),  INTENT(OUT)   :: PTOTAL(NUMTAU)
  REAL(8),  INTENT(OUT)   :: P_out(NUMTAU)
  REAL(8),  INTENT(IN)    :: PRAD_in(NUMTAU)
  REAL(8),  INTENT(IN)    :: PTURB_in(NUMTAU)
  REAL(8),  INTENT(IN)    :: VTURB_in(NUMTAU)
  REAL(8),  INTENT(IN)    :: GRAV_in

  ! --- Tuning parameters ---
  INTEGER, PARAMETER :: MAX_CORRECTOR  = 1000
  REAL(8),  PARAMETER :: CONV_TOL       = 5.0D-5   ! convergence in ln(P)
  REAL(8),  PARAMETER :: DPLOG_MAX      = 10.0D0   ! max d(ln P) per step
  REAL(8),  PARAMETER :: OPACITY_REBOOT = 2.0D0    ! opacity ratio triggering surface redo
  INTEGER, PARAMETER :: RELAX_TRIGGER  = 4        ! stall count before under-relaxation
  REAL(8),  PARAMETER :: RELAX_FACTOR   = 0.3D0    ! relaxation weight when oscillating

  ! --- Local variables ---
  REAL(8)  :: DLGTAU              ! log spacing: ln(TAU(2)/TAU(1))
  REAL(8)  :: PLOG                ! current ln(P_total) estimate
  REAL(8)  :: PNEW               ! corrected ln(P_total)
  REAL(8)  :: DPLOG              ! d(ln P)/d(ln tau) * DLGTAU at current point
  REAL(8)  :: ERROR              ! |PNEW - PLOG| convergence measure
  REAL(8)  :: PREV_ERROR         ! previous corrector error
  REAL(8)  :: ABSTD1_INIT        ! initial surface opacity guess (for bootstrap)
  REAL(8)  :: RELAX              ! relaxation weight for current step
  INTEGER :: N_STALL            ! count of stalled corrector steps

  ! Adams-Bashforth history (previous 4 points)
  REAL(8)  :: PLOG1, PLOG2, PLOG3, PLOG4
  REAL(8)  :: DPLOG1, DPLOG2, DPLOG3

  INTEGER :: J, N
  LOGICAL :: CONVERGED

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING TTAUP'

  ! Log spacing of the tau grid (must be uniform)
  DLGTAU = log(TAU(2) / TAU(1))

  ! Initialize Adams-Bashforth history
  PLOG1  = 0.0D0
  PLOG2  = 0.0D0
  PLOG3  = 0.0D0
  DPLOG1 = 0.0D0
  DPLOG2 = 0.0D0

  ! Initial opacity guess: 0.1 cm^2/g, or smaller if radiation pressure
  ! is known (ensures P_total > P_rad at surface)
  ABSTD(1) = 0.1D0
  IF (PRAD_in(1) .GT. 0.0D0) THEN
    ABSTD(1) = min(0.1D0, GRAV_in * TAU(1) / PRAD_in(1) / 2.0D0)
  END IF
  ABSTD1_INIT = ABSTD(1)

  !---------------------------------------------------------------------
  ! March inward through the atmosphere
  !---------------------------------------------------------------------
  DO J = 1, NUMTAU

    ! --- Predictor: initial guess for ln(P_total) ---
    IF (J .EQ. 1) THEN
      ! Surface: P_total = g * tau / kappa (constant-opacity approximation)
      PLOG = log(GRAV_in / ABSTD(1) * TAU(1))
    ELSE IF (J .LE. 4) THEN
      ! Linear extrapolation from previous point
      PLOG = PLOG1 + DPLOG1
    ELSE
      ! 4th-order Adams-Bashforth predictor
      PLOG = (3.0D0 * PLOG4 + 8.0D0 * DPLOG1 &
            - 4.0D0 * DPLOG2 + 8.0D0 * DPLOG3) / 3.0D0
    END IF

    ! --- Corrector iterations ---
    ERROR      = 1.0D0
    PREV_ERROR = 1.0D0
    N          = 0
    CONVERGED  = .FALSE.
    N_STALL    = 0

    DO
      ! Convert ln(P_total) → total pressure → gas pressure → opacity
      PTOTAL(J) = exp(PLOG)

      ! Gas pressure: subtract radiation and turbulent pressure changes
      ! relative to their surface values
      P_out(J) = PTOTAL(J) + (PRAD_in(1) - PRAD_in(J)) &
                            + (PTURB_in(1) - PTURB_in(J))

      IF (P_out(J) .LE. 0.0D0) THEN
        ! Gas pressure went negative — corrector guess for P_total is too low.
        ! Floor P_gas to a small fraction of P_total so the opacity lookup
        ! returns a large kappa, which drives DPLOG upward and self-corrects.
        P_out(J) = PTOTAL(J) * 1.0D-4
        IF (P_out(J) .LE. 0.0D0) P_out(J) = 1.0D-10
      END IF

      ! Look up Rosseland opacity from the interpolation table
      ABSTD(J) = ROSSTAB(T_in(J), P_out(J), VTURB_in(J))

      ! Hydrostatic equilibrium derivative
      DPLOG = GRAV_in / ABSTD(J) * TAU(J) / PTOTAL(J) * DLGTAU

      ! Clamp to prevent unphysical pressure jumps
      DPLOG = max(0.0D0, min(DPLOG_MAX, DPLOG))

      N = N + 1
      IF (N .EQ. 1) CYCLE    ! first pass: go straight to corrector

      ! --- Corrector: update ln(P) ---
      IF (J .EQ. 1) THEN
        PNEW = log(GRAV_in / ABSTD(1) * TAU(1))
      ELSE IF (J .LE. 4) THEN
        ! Trapezoidal corrector
        PNEW = (PLOG + 2.0D0 * PLOG1 + DPLOG + DPLOG1) / 3.0D0
      ELSE
        ! 4th-order Adams-Moulton corrector
        PNEW = (126.0D0 * PLOG1 - 14.0D0 * PLOG3 + 9.0D0 * PLOG4 &
              + 42.0D0 * DPLOG + 108.0D0 * DPLOG1 &
              - 54.0D0 * DPLOG2 + 24.0D0 * DPLOG3) / 121.0D0
      END IF

      PREV_ERROR = ERROR
      ERROR = abs(PNEW - PLOG)

      ! Detect stalling: error not decreasing meaningfully
      IF (N .GT. 3 .AND. ERROR .GT. 0.9D0 * PREV_ERROR) THEN
        N_STALL = N_STALL + 1
      ELSE
        N_STALL = max(N_STALL - 1, 0)
      END IF

      ! Choose relaxation weight:
      !   Normal: average old and new (0.5), same as original code
      !   Oscillating: bias toward old value to damp oscillation
      IF (N_STALL .GE. RELAX_TRIGGER) THEN
        RELAX = RELAX_FACTOR
      ELSE
        RELAX = 0.5D0
      END IF
      PLOG = (1.0D0 - RELAX) * PLOG + RELAX * PNEW

      IF (ERROR .LE. CONV_TOL) THEN
        CONVERGED = .TRUE.
        EXIT
      END IF
      IF (N .GT. MAX_CORRECTOR) EXIT
    END DO

    ! Warn on convergence failure
    IF (.NOT. CONVERGED) THEN
      WRITE(6, '(A,I4,A,ES10.3,A,I5)') &
        ' TTAUP WARNING: corrector did not converge at J=', J, &
        ', error=', ERROR, ', iter=', N
    END IF

    ! --- Surface bootstrap ---
    ! If the converged opacity at J=1 differs substantially from the
    ! initial blind guess, redo J=1 with the better opacity value.
    ! This prevents a poor surface guess from corrupting the entire
    ! Adams-Bashforth history.
    IF (J .EQ. 1 .AND. (ABSTD(1) / ABSTD1_INIT .GT. OPACITY_REBOOT .OR. &
                       ABSTD1_INIT / ABSTD(1) .GT. OPACITY_REBOOT)) THEN
      ABSTD1_INIT = ABSTD(1)   ! prevent infinite re-bootstrapping
      PLOG = log(GRAV_in / ABSTD(1) * TAU(1))
      PTOTAL(1) = exp(PLOG)
      P_out(1)  = PTOTAL(1)    ! at surface: PRAD(1)-PRAD(1) = 0, same for PTURB
      IF (P_out(1) .GT. 0.0D0) THEN
        ABSTD(1) = ROSSTAB(T_in(1), P_out(1), VTURB_in(1))
        PLOG     = log(GRAV_in / ABSTD(1) * TAU(1))
        PTOTAL(1) = exp(PLOG)
        P_out(1)  = PTOTAL(1)
      END IF
      DPLOG = GRAV_in / ABSTD(1) * TAU(1) / PTOTAL(1) * DLGTAU
      DPLOG = max(0.0D0, min(DPLOG_MAX, DPLOG))
    END IF

    ! --- Shift Adams-Bashforth history ---
    PLOG4  = PLOG3
    PLOG3  = PLOG2
    PLOG2  = PLOG1
    PLOG1  = PLOG
    DPLOG3 = DPLOG2
    DPLOG2 = DPLOG1
    DPLOG1 = DPLOG

  END DO

  RETURN

END SUBROUTINE TTAUP



!=======================================================================
! READIN: Read and parse all input data
!=======================================================================

SUBROUTINE READIN(MODE)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: MODE
!     MODE=1  COMPUTE A MODEL (ATLAS12: reads input_model.dat + stdin)
!     MODE=20 READ A MODEL FOR SYNTHESIS (SYNTHE: caller opens file on unit 5)

  ! Local scalars
  REAL(8)  :: XSCALELOG
  INTEGER :: I, IZ
  INTEGER :: J, MU
  INTEGER :: IOS_CARD
  LOGICAL :: FIRST_KW
  CHARACTER(20) :: KEYWORD

  CHARACTER(1) :: CARD(132)

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING READIN'

  IF (MODE .EQ. 1) THEN
    !-------------------------------------------------------------------
    ! ATLAS12 path: read model from input_model.dat, then stdin overrides
    !-------------------------------------------------------------------
    OPEN(UNIT=3, FILE='input_model.dat', STATUS='OLD', ACTION='READ')
    INPUTDATA = 3
    IFPRES = 1
    IFCORR = 1
    CALL READMOL
  ELSE
    !-------------------------------------------------------------------
    ! SYNTHE path: caller has already opened the model file on unit 5.
    ! Force surface flux mode with a single iteration and full output.
    !   SURFACE FLUX
    !   ITERATIONS 1 PRINT 2 PUNCH 2
    !-------------------------------------------------------------------
    INPUTDATA = 5
    IFPRES = 0
    IFSURF = 1
    NUMITS = 1
    IFPRNT(1) = 2
    IFPNCH(1) = 2
    CALL READMOL
  END IF

  LAST=133
  MAXPOW=38

  !---------------------------------------------------------------------
  ! Main card-reading loop
  !---------------------------------------------------------------------
  card_loop: DO
    MORE = 0
    LETCOL = 1
    READ(INPUTDATA, '(132A1)', IOSTAT=IOS_CARD) CARD
    IF (IOS_CARD .NE. 0) EXIT card_loop

    !-------------------------------------------------------------------
    ! Keyword parsing loop (multiple keywords per card)
    !-------------------------------------------------------------------
    KEYWORD = NEXTWORD(CARD)
    NUMCOL = LETCOL
    FIRST_KW = .TRUE.

    keyword_loop: DO
      IF (.NOT. FIRST_KW) THEN
        ! Get next keyword from same card
        LETCOL = MAX(LETCOL, NUMCOL)
        MORE = 1
        KEYWORD = NEXTWORD(CARD)
        IF (IFFAIL .EQ. 1) EXIT keyword_loop
        MORE = 0
        NUMCOL = LETCOL
      END IF
      FIRST_KW = .FALSE.

      !-----------------------------------------------------------------
      ! TEFF
      !-----------------------------------------------------------------
      IF (trim(KEYWORD) .EQ. 'TEFF') THEN
        TEFF = FREEFF(CARD)
        FLUX = SIGMA_SB / FOURPI * TEFF**4
        CYCLE keyword_loop

      !-----------------------------------------------------------------
      ! GRAVITY
      !-----------------------------------------------------------------
      ELSE IF (trim(KEYWORD) .EQ. 'GRAVITY') THEN
        GRAV = FREEFF(CARD)
        IF (GRAV .LT. 10.0D0) GRAV = 10.0D0**(GRAV)
        GLOG = LOG10(GRAV)
        CYCLE keyword_loop

      !-----------------------------------------------------------------
      ! ABUNDANCE (sub-keywords: SCALE, CHANGE, ABSOLUTE, RELATIVE, TABLE)
      !-----------------------------------------------------------------
      ELSE IF (trim(KEYWORD) .EQ. 'ABUNDANCE') THEN
        KEYWORD = NEXTWORD(CARD)
        ! SCALE
        IF (trim(KEYWORD) .EQ. 'SCALE') THEN
          NUMCOL = LETCOL
          XSCALE = FREEFF(CARD)
          IF (XSCALE .GT. 0.0D0) XSCALELOG = LOG10(XSCALE)
          DO IZ = 3, 99
            XRELATIVE(IZ) = XSCALELOG
          END DO
          XSCALE = 1.0D0
          CYCLE keyword_loop
        ! CHANGE
        ELSE IF (trim(KEYWORD) .EQ. 'CHANGE') THEN
          MORE = 1
          DO
            IZ = FREEFF(CARD)
            IF (IFFAIL .EQ. 1) EXIT keyword_loop
            ABUND(IZ) = FREEFF(CARD)
            IF (IZ .GT. 2 .AND. ABUND(IZ) .GT. 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
          END DO
        ! ABSOLUTE
        ELSE IF (trim(KEYWORD) .EQ. 'ABSOLUTE') THEN
          MORE = 1
          DO
            IZ = FREEFF(CARD)
            IF (IFFAIL .EQ. 1) EXIT keyword_loop
            ABUND(IZ) = FREEFF(CARD)
            IF (IZ .GT. 2 .AND. ABUND(IZ) .GT. 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
            XRELATIVE(IZ) = 0.0D0
          END DO
        ! RELATIVE
        ELSE IF (trim(KEYWORD) .EQ. 'RELATIVE') THEN
          MORE = 1
          DO
            IZ = FREEFF(CARD)
            IF (IFFAIL .EQ. 1) EXIT keyword_loop
            XRELATIVE(IZ) = FREEFF(CARD)
          END DO
        ! TABLE
        ELSE IF (trim(KEYWORD) .EQ. 'TABLE') THEN
          READ(INPUTDATA, '(7X,F10.6,10X,F10.6/(5(7X,F7.3,F6.3)))') &
            ABUND(1), ABUND(2), (ABUND(IZ), XRELATIVE(IZ), IZ=3,99)
          DO IZ = 3, 99
            IF (ABUND(IZ) .GT. 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
          END DO
          XSCALE = 1.0D0
          EXIT keyword_loop
        ELSE
          WRITE(6, '(A,A)') ' ABUNDANCE: unknown sub-keyword: ', trim(KEYWORD)
          STOP 1
        END IF

      !-----------------------------------------------------------------
      ! READ (sub-keyword: DECK/DECK6)
      !-----------------------------------------------------------------
      ELSE IF (trim(KEYWORD) .EQ. 'READ') THEN
        KEYWORD = NEXTWORD(CARD)
        NUMCOL = LETCOL

        ! DECK / DECK6
        IF (trim(KEYWORD) .EQ. 'DECK' .OR. trim(KEYWORD) .EQ. 'DECK6') THEN
          NRHOX = FREEFF(CARD)
          DO J = 1, NRHOX
            NUMCOL = 1
            READ(INPUTDATA, '(132A1)') CARD
            RHOX(J) = FREEFF(CARD)
            T(J) = FREEFF(CARD)
            MORE = 1
            P(J) = FREEFF(CARD)
            XNE(J) = FREEFF(CARD)
            ABROSS(J) = FREEFF(CARD)
            ! Column 6 is ACCRAD (labelled PRAD temporarily)
            PRAD(J) = FREEFF(CARD)
            VTURB(J) = FREEFF(CARD)
            MORE = 0
          END DO
          IF (RHOX(1) .LT. 0.0D0) THEN
            RHOX = 10.0D0**(RHOX)
          END IF
          PRADK0 = 0.0D0
          PTURB0 = PTURB(1)
          PCON = 0.0D0
          PZERO = PCON + PRADK0 + PTURB0
          CALL INTEG(RHOX, ABROSS, TAUROS, NRHOX, ABROSS(1)*RHOX(1))
          IF (trim(KEYWORD) .EQ. 'DECK6') THEN
            ! DECK6: read additional PRADK0 card
            READ(INPUTDATA, '(132A1)') CARD
            NUMCOL = 1
            PRADK0 = FREEFF(CARD)
            ACCRAD = PRAD
            CALL INTEG(RHOX, ACCRAD, PRAD, NRHOX, ACCRAD(1)*RHOX(1))
            PRADK = PRAD + PRADK0
          END IF
          EXIT keyword_loop

        ELSE
          WRITE(6, '(A,A)') ' READ: unknown sub-keyword: ', trim(KEYWORD)
          STOP 1
        END IF

      !-----------------------------------------------------------------
      ! BEGIN — end of model file; exit card loop for finalization
      !-----------------------------------------------------------------
      ELSE IF (trim(KEYWORD) .EQ. 'BEGIN') THEN
        IF (INPUTDATA .EQ. 3) CLOSE(UNIT=3)
        EXIT card_loop

      !-----------------------------------------------------------------
      ! TURBULENCE (sub-keywords: ON, OFF)
      !-----------------------------------------------------------------
      ELSE IF (trim(KEYWORD) .EQ. 'TURBULENCE') THEN
        KEYWORD = NEXTWORD(CARD)
        IF (trim(KEYWORD) .EQ. 'ON') THEN
          IFTURB = 1
          NUMCOL = LETCOL
          TRBFDG = FREEFF(CARD)
          TRBPOW = FREEFF(CARD)
          TRBSND = FREEFF(CARD)
          TRBCON = FREEFF(CARD)
        ELSE IF (trim(KEYWORD) .EQ. 'OFF') THEN
          IFTURB = 0
          TRBFDG = 0.0D0
          TRBPOW = 0.0D0
          TRBSND = 0.0D0
          TRBCON = 0.0D0
        ELSE
          WRITE(6, '(A,A)') ' TURBULENCE: unknown sub-keyword: ', trim(KEYWORD)
          STOP 1
        END IF
        CYCLE keyword_loop

      !-----------------------------------------------------------------
      ! SURFACE (sub-keywords: INTENSITY, FLUX, OFF)
      !-----------------------------------------------------------------
      ELSE IF (trim(KEYWORD) .EQ. 'SURFACE') THEN
        KEYWORD = NEXTWORD(CARD)
        IF (trim(KEYWORD) .EQ. 'INTENSITY') THEN
          NMU = FREEFF(CARD)
          DO MU = 1, NMU
            ANGLE(MU) = FREEFR(CARD)
          END DO
          IFSURF = 2
        ELSE IF (trim(KEYWORD) .EQ. 'FLUX') THEN
          IFSURF = 1
        ELSE IF (trim(KEYWORD) .EQ. 'OFF') THEN
          IFSURF = 0
        ELSE
          WRITE(6, '(A,A)') ' SURFACE: unknown sub-keyword: ', trim(KEYWORD)
          STOP 1
        END IF
        CYCLE keyword_loop

      !-----------------------------------------------------------------
      ! Unrecognized keyword — skip to next card
      !-----------------------------------------------------------------
      ELSE
        EXIT keyword_loop
      END IF

    END DO keyword_loop
  END DO card_loop

  ! ===================================================================
  ! FINALIZATION — compute derived quantities
  ! ===================================================================
  IF (ABUND(1) .LT. 0.0D0) ABUND(1) = 10.0D0**ABUND(1)
  IF (ABUND(2) .LT. 0.0D0) ABUND(2) = 10.0D0**ABUND(2)
  DO IZ = 3, 99
    IF (ABUND(IZ) .GT. 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
  END DO
  ! Write abbreviated list of abundances
  WRITE(6, '(/" TEFF",F7.0,"   LOGG",F8.4/' // &
       '"   1",A2,F10.6,"     2",A2,F10.6/(5(I4,A2,F7.2,F5.2)))') &
       TEFF, GLOG, ELEM(1), ABUND(1), ELEM(2), ABUND(2), &
     (IZ, ELEM(IZ), ABUND(IZ), XRELATIVE(IZ), IZ=3,32)
  DO J = 1, NRHOX
    DO IZ = 3, 99
      XABUND(J,IZ) = 10.0D0**(ABUND(IZ) + XRELATIVE(IZ))
    END DO
    XABUND(J,1) = ABUND(1)
    XABUND(J,2) = ABUND(2)
    WTMOLE(J) = 0.0D0
    DO IZ = 1, 99
      WTMOLE(J) = WTMOLE(J) + XABUND(J,IZ) * ATMASS(IZ)
    END DO
  END DO
  YABUND(1) = ABUND(1)
  YABUND(2) = ABUND(2)
  DO IZ = 3, 99
    YABUND(IZ) = ABUND(IZ) + XRELATIVE(IZ)
  END DO
  DO J = 1, NRHOX
    TK(J) = KBOL * T(J)
    HKT(J) = HPLANCK / TK(J)
    HCKT(J) = HKT(J) * CLIGHT
    TKEV(J) = KBOL_EV * T(J)
    TLOG(J) = LOG(T(J))
    XNATOM(J) = P(J) / TK(J) - XNE(J)
    RHO(J) = XNATOM(J) * WTMOLE(J) * AMU
    IF (IFTURB .GT. 0) PTURB(J) = 0.5D0 * RHO(J) * VTURB(J)**2
    CHARGESQ(J) = XNE(J)
  END DO

  IF (IFSYNTHE.EQ.0) THEN
     !only write if running atlas12
     WRITE(6, *)
     WRITE(6, '(" NUMITS",I3,"  IFPRNT",*(I2))') NUMITS, IFPRNT(1:NUMITS)
     WRITE(6, '(10X,"  IFPNCH",*(I2))') IFPNCH(1:NUMITS)
  ENDIF

  ! Auto-generate title from mixing length
  BLOCK
    CHARACTER(74) :: TITLEBUF
    WRITE(TITLEBUF, '(A,F5.2)') 'ATLAS12 l/H=', MIXLTH
    DO I = 1, 74
      TITLE(I) = TITLEBUF(I:I)
    END DO
  END BLOCK

END SUBROUTINE READIN


!=======================================================================
! SUBROUTINE SCALE_MODEL(TEFF_NEW, LOGG_NEW)
!
! Regrid the current model onto a standard uniform log-tau grid and
! optionally rescale to a new Teff and/or log g.
!
! The tau grid is defined by three parameters (hardcoded below):
!   NDEPTHS  — number of depth points (80)
!   TAU1LG   — log10 of smallest Rosseland optical depth (-6.875)
!   STEPLG   — log10 step between successive depth points (0.125)
! This gives tau from 10^{-6.875} to 10^{+2.875} in 80 steps.
!
! The routine:
!   1. Builds the new tau grid
!   2. Integrates ABROSS to get TAUROS from the current model
!   3. Interpolates all structure variables onto the new grid
!   4. If Teff/logg differ from current values, rescales T, Prad,
!      and recomputes hydrostatic equilibrium
!=======================================================================

SUBROUTINE SCALE_MODEL(TEFF_NEW, LOGG_NEW)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: TEFF_NEW, LOGG_NEW

  ! --- Tau grid parameters ---
  ! NDEPTHS: number of depth points in the output model
  ! TAU1LG and STEPLG are module-level variables (mod_atlas_data),
  ! also used by TCORR to rebuild the grid during iterations.
  ! The deepest point is at log10(tau) = TAU1LG + (NDEPTHS-1)*STEPLG
  INTEGER, PARAMETER :: NDEPTHS = 80

  ! Local work arrays for interpolation
  REAL(8)  :: DUM1(kw), DUM2(kw), DUM3(kw), DUM4(kw)
  REAL(8)  :: DUM5(kw), DUM6(kw), DUM7(kw), DUM8(kw)
  REAL(8)  :: GNEW
  INTEGER :: I, J, IDUM

  ! Convert logg to linear gravity
  GNEW = 10.0D0**LOGG_NEW

  ! --- Step 1: Build new uniform log-tau grid ---
  DO J = 1, NDEPTHS
    TAUSTD(J) = 10.0D0**(TAU1LG + (J-1)*STEPLG)
  END DO

  ! --- Step 2: Integrate ABROSS to get TAUROS from current model ---
  CALL INTEG(RHOX, ABROSS, TAUROS, NRHOX, ABROSS(1)*RHOX(1))

  ! --- Step 3: Interpolate all structure variables onto new grid ---
  IDUM = MAP1(TAUROS, RHOX,   NRHOX, TAUSTD, DUM1, NDEPTHS)
  IDUM = MAP1(TAUROS, T,      NRHOX, TAUSTD, DUM2, NDEPTHS)
  IDUM = MAP1(TAUROS, P,      NRHOX, TAUSTD, DUM3, NDEPTHS)
  IDUM = MAP1(TAUROS, XNE,    NRHOX, TAUSTD, DUM4, NDEPTHS)
  IDUM = MAP1(TAUROS, ABROSS, NRHOX, TAUSTD, DUM5, NDEPTHS)
  IDUM = MAP1(TAUROS, PRAD,   NRHOX, TAUSTD, DUM6, NDEPTHS)
  IDUM = MAP1(TAUROS, VTURB,  NRHOX, TAUSTD, DUM7, NDEPTHS)
  IDUM = MAP1(TAUROS, BMIN,   NRHOX, TAUSTD, DUM8, NDEPTHS)
  DO J = 1, NDEPTHS
    RHOX(J)   = DUM1(J)
    T(J)      = DUM2(J)
    P(J)      = DUM3(J)
    XNE(J)    = DUM4(J)
    ABROSS(J) = DUM5(J)
    PRAD(J)   = DUM6(J)
    PRADK(J)  = PRAD(J) + PRADK0
    VTURB(J)  = DUM7(J)
    BMIN(J)   = DUM8(J)
  END DO
  DO I = 1, 6
    IDUM = MAP1(TAUROS, BHYD(1,I), NRHOX, TAUSTD, DUM1, NDEPTHS)
    DO J = 1, NDEPTHS
      BHYD(J,I) = DUM1(J)
    END DO
  END DO
  NRHOX = NDEPTHS

  ! --- Step 4: Rescale to new Teff/logg if different ---
  IF (TEFF_NEW .EQ. 0.0D0) RETURN
  IF (TEFF_NEW .EQ. TEFF .AND. GNEW .EQ. GRAV) RETURN
  IF (TEFF_NEW .LT. TEFF+1.0D0 .AND. TEFF_NEW .GT. TEFF-1.0D0 .AND. &
      GNEW .LT. GRAV*1.001D0 .AND. GNEW .GT. GRAV*0.999D0) RETURN

  TAUROS = TAUSTD
  T      = T * TEFF_NEW / TEFF
  PTURB  = 0.0D0
  PRADK  = PRADK * (TEFF_NEW/TEFF)**4
  PRAD   = PRAD  * (TEFF_NEW/TEFF)**4
  PRADK0 = PRADK0 * (TEFF_NEW/TEFF)**4
  PZERO  = PCON + PRADK0 + PTURB0
  TEFF   = TEFF_NEW
  FLUX   = SIGMA_SB / FOURPI * TEFF**4
  GRAV   = GNEW
  GLOG   = LOG10(GRAV)
  PTOTAL = P + PRAD + PTURB
  RHOX   = PTOTAL / GRAV
  PTOTAL = PTOTAL + PZERO

END SUBROUTINE SCALE_MODEL

!=======================================================================
! FREEFR: Free-format real number reader
!=======================================================================

FUNCTION FREEFR(CARD)

  IMPLICIT NONE
  CHARACTER(1), INTENT(INOUT) :: CARD(*)
  REAL(8) :: FREEFR
  INTEGER :: I, L

  MORE = 1
  FREEFR = FREEFF(CARD)
  IF (IFFAIL .EQ. 0) RETURN
  L = LAST - 1
  READ(INPUTDATA, 1) (CARD(I), I=1, L)
    1 FORMAT(132A1)
  NUMCOL = 1
  FREEFR = FREEFF(CARD)
  RETURN

END FUNCTION FREEFR


!=========================================================================
! FUNCTION FREEFF(CARD)
!
! Free-format floating-point number reader. Scans CARD starting at module
! variable NUMCOL, parses a number (with optional sign, decimal point,
! and E-exponent), and returns it as real*8.
!
! Recognized formats: 123  -45.6  1.23E4  .5E-3  7,  etc.
! Delimiters: blank or comma.
!
! Module variables updated:
!   NUMCOL  — advanced past the number (or past LAST on failure)
!   IFFAIL  — set to 1 if end-of-card reached without finding a number
!=========================================================================

FUNCTION FREEFF(CARD)

  IMPLICIT NONE

  CHARACTER(1), INTENT(IN) :: CARD(*)
  REAL(8) :: FREEFF

  ! Digit lookup table: 0-9
  CHARACTER(1), PARAMETER :: DIGIT(10) = (/ &
    '0','1','2','3','4','5','6','7','8','9' /)

  ! Parser states
  INTEGER, PARAMETER :: ST_INTEGER  = 1  ! before decimal point
  INTEGER, PARAMETER :: ST_FRACTION = 2  ! after decimal point
  INTEGER, PARAMETER :: ST_EXPSIGN  = 3  ! after E (expecting sign or digit)
  INTEGER, PARAMETER :: ST_EXPDIGIT = 4  ! exponent digits

  CHARACTER(1) :: C
  REAL(8)  :: ANSWER, ASIGN
  INTEGER :: I, N, NCOL, NPT, NPOWER, ISIGN, IF0, STATE

  ! Initialize
  IFFAIL = 0
  IF (NUMCOL .GT. LAST) THEN
    IFFAIL = 1
    IF (MORE .GT. 0) THEN
      FREEFF = 0.0D0
      RETURN
    END IF
    WRITE(6, '(A/(1X,131A1))') '1FREEFF HAS READ OFF THE END', &
      (CARD(I), I = 1, LAST)
    STOP 1
  END IF

  ANSWER = 0.0D0
  ASIGN  = 1.0D0
  ISIGN  = 1
  NPT    = 0
  IF0    = 0
  N      = 0
  STATE  = ST_INTEGER

  !---------------------------------------------------------------------
  ! Main scan loop
  !---------------------------------------------------------------------
  scan_loop: DO NCOL = NUMCOL, LAST
    C = CARD(NCOL)

    SELECT CASE (STATE)

    !-----------------------------------------------------------------
    ! State 1: Integer part (before decimal point)
    !-----------------------------------------------------------------
    CASE (ST_INTEGER)
      IF (C .EQ. ' ') THEN
        IF (N .EQ. 0) THEN
          CALL freeff_reset()
          CYCLE scan_loop
        END IF
        ! Blank after digits — return integer value
        FREEFF = ANSWER * ASIGN
        NUMCOL = NCOL + 1
        RETURN
      END IF
      ! Check for digit
      DO I = 1, 10
        IF (C .EQ. DIGIT(I)) THEN
          N = N + 1
          ANSWER = 10.0D0 * ANSWER + dble(I - 1)
          CYCLE scan_loop
        END IF
      END DO
      ! Check for decimal point
      IF (C .EQ. '.') THEN
        STATE = ST_FRACTION
        CYCLE scan_loop
      END IF
      ! Check for comma (delimiter)
      IF (C .EQ. ',') THEN
        IF (N .EQ. 0) THEN
          ! Reset — comma before any digits
          CALL freeff_reset()
          CYCLE scan_loop
        END IF
        FREEFF = ANSWER * ASIGN
        NUMCOL = NCOL + 1
        RETURN
      END IF
      ! Check for minus sign
      IF (C .EQ. '-') THEN
        IF (N .EQ. 0) THEN
          ASIGN = -1.0D0
          CYCLE scan_loop
        END IF
        ! Minus after digits — reset
        CALL freeff_reset()
        CYCLE scan_loop
      END IF
      ! Unrecognized character — reset
      CALL freeff_reset()

    !-----------------------------------------------------------------
    ! State 2: Fractional part (after decimal point)
    !-----------------------------------------------------------------
    CASE (ST_FRACTION)
      ! Check for digit
      DO I = 1, 10
        IF (C .EQ. DIGIT(I)) THEN
          N = N + 1
          NPT = NPT + 1
          ANSWER = 10.0D0 * ANSWER + dble(I - 1)
          CYCLE scan_loop
        END IF
      END DO
      IF (C .EQ. 'E') THEN
        STATE = ST_EXPSIGN
        CYCLE scan_loop
      END IF
      IF (C .EQ. '-') THEN
        ! Minus after decimal digits — treat as exponent sign
        ISIGN = -1
        NPOWER = 0
        STATE = ST_EXPDIGIT
        CYCLE scan_loop
      END IF
      IF (C .EQ. '+') THEN
        NPOWER = 0
        STATE = ST_EXPDIGIT
        CYCLE scan_loop
      END IF
      IF (C .EQ. ' ' .OR. C .EQ. ',') THEN
        IF (N .EQ. 0) THEN
          CALL freeff_reset()
          CYCLE scan_loop
        END IF
        FREEFF = ANSWER * ASIGN / 10.0D0**NPT
        NUMCOL = NCOL + 1
        RETURN
      END IF
      ! Unrecognized — reset
      CALL freeff_reset()

    !-----------------------------------------------------------------
    ! State 3: After E (expecting sign or first exponent digit)
    !-----------------------------------------------------------------
    CASE (ST_EXPSIGN)
      ! Check for digit
      DO I = 1, 10
        IF (C .EQ. DIGIT(I)) THEN
          NPOWER = I - 1
          IF0 = 1
          STATE = ST_EXPDIGIT
          CYCLE scan_loop
        END IF
      END DO
      IF (C .EQ. ' ' .OR. C .EQ. '+') THEN
        NPOWER = 0
        STATE = ST_EXPDIGIT
        CYCLE scan_loop
      END IF
      IF (C .EQ. '-') THEN
        ISIGN = -1
        NPOWER = 0
        STATE = ST_EXPDIGIT
        CYCLE scan_loop
      END IF
      ! Unrecognized — reset
      CALL freeff_reset()

    !-----------------------------------------------------------------
    ! State 4: Exponent digits
    !-----------------------------------------------------------------
    CASE (ST_EXPDIGIT)
      ! Check for digit
      DO I = 1, 10
        IF (C .EQ. DIGIT(I)) THEN
          NPOWER = 10 * NPOWER + I - 1
          IF0 = 1
          IF (NPOWER .GE. MAXPOW) THEN
            CALL freeff_reset()
            CYCLE scan_loop
          END IF
          CYCLE scan_loop
        END IF
      END DO
      IF (C .EQ. ',' .OR. C .EQ. ' ') THEN
        IF (IF0 .EQ. 0) THEN
          CALL freeff_reset()
          CYCLE scan_loop
        END IF
        FREEFF = ANSWER * ASIGN * 10.0D0**(ISIGN * NPOWER - NPT)
        NUMCOL = NCOL + 1
        RETURN
      END IF
      ! Unrecognized — reset
      CALL freeff_reset()

    END SELECT
  END DO scan_loop

  !---------------------------------------------------------------------
  ! Fell off end of card
  !---------------------------------------------------------------------
  NUMCOL = LAST + 1
  IFFAIL = 1
  IF (MORE .GT. 0) THEN
    FREEFF = 0.0D0
    RETURN
  END IF
  WRITE(6, '(A/(1X,131A1))') '1FREEFF HAS READ OFF THE END', &
    (CARD(I), I = 1, LAST)
  STOP 1

CONTAINS

  SUBROUTINE freeff_reset()
    ! Reset parser to initial searching state
    ASIGN  = 1.0D0
    ANSWER = 0.0D0
    ISIGN  = 1
    NPT    = 0
    IF0    = 0
    N      = 0
    STATE  = ST_INTEGER
  END SUBROUTINE freeff_reset

END FUNCTION FREEFF

!=========================================================================
! FUNCTION NEXTWORD(CARD)
!
! Free-format word reader. Scans CARD starting at module variable LETCOL,
! finds the next alphanumeric word, and returns it as an uppercase
! character string (up to 20 characters).
!
! Special case: an "E" preceded by a digit or "." and followed by a digit
! or blank is treated as scientific notation, not a word start.
!
! Module variables updated:
!   LETCOL  — advanced past the word (or past LAST on failure)
!   IFFAIL  — set to 1 if end-of-card reached without finding a word
!=========================================================================

FUNCTION NEXTWORD(CARD)

  IMPLICIT NONE

  CHARACTER(1), INTENT(IN) :: CARD(*)
  CHARACTER(20) :: NEXTWORD

  CHARACTER(1) :: C
  INTEGER :: I, NCOL, N
  LOGICAL :: IN_WORD, SKIP_E, IS_ALPHA

  ! Initialize
  NEXTWORD = ' '
  IFFAIL = 0
  N = 0
  IN_WORD = .FALSE.

  IF (LETCOL .GT. LAST) THEN
    IFFAIL = 1
    IF (MORE .GT. 0) RETURN
    WRITE(6, '(A/(1X,131A1))') 'NEXTWORD HAS READ OFF THE END', &
      (CARD(I), I = 1, LAST)
    STOP 1
  END IF

  !---------------------------------------------------------------------
  ! Main scan loop
  !---------------------------------------------------------------------
  scan_loop: DO NCOL = LETCOL, LAST
    C = CARD(NCOL)

    IF (.NOT. IN_WORD) THEN
      !--- State: searching for word start ---
      IF (C .EQ. ' ') CYCLE scan_loop

      ! Check if C is a letter (A-Z)
      IS_ALPHA = (C .GE. 'A' .AND. C .LE. 'Z')
      IF (.NOT. IS_ALPHA) THEN
        ! Not a letter — reset and continue searching
        N = 0
        CYCLE scan_loop
      END IF

      ! Special case: is this an "E" in scientific notation?
      IF (C .EQ. 'E' .AND. NCOL .GT. 1) THEN
        SKIP_E = .FALSE.
        ! Check preceding character: digit or "."
        IF ((CARD(NCOL-1) .GE. '0' .AND. CARD(NCOL-1) .LE. '9') .OR. &
            CARD(NCOL-1) .EQ. '.') THEN
          ! Check following character: digit or blank
          IF (CARD(NCOL+1) .EQ. ' ' .OR. &
              (CARD(NCOL+1) .GE. '0' .AND. CARD(NCOL+1) .LE. '9')) THEN
            SKIP_E = .TRUE.
          END IF
        END IF
        IF (SKIP_E) THEN
          N = 0
          CYCLE scan_loop
        END IF
      END IF

      ! Start a new word
      N = 1
      NEXTWORD(1:1) = C
      IN_WORD = .TRUE.

    ELSE
      !--- State: inside a word ---
      IF (C .EQ. ' ' .OR. C .EQ. '=' .OR. C .EQ. ',') THEN
        ! Word delimiter — return the word
        LETCOL = NCOL + 1
        RETURN
      END IF

      ! Check if C is alphanumeric (A-Z or 0-9)
      IS_ALPHA = (C .GE. 'A' .AND. C .LE. 'Z') .OR. &
                 (C .GE. '0' .AND. C .LE. '9')
      IF (.NOT. IS_ALPHA) THEN
        ! Not alphanumeric — abandon word, reset to searching
        NEXTWORD = ' '
        N = 0
        IN_WORD = .FALSE.
        CYCLE scan_loop
      END IF

      ! Accumulate character into word (max 20 chars)
      N = N + 1
      IF (N .LE. 20) THEN
        NEXTWORD(N:N) = C
      END IF
    END IF
  END DO scan_loop

  !---------------------------------------------------------------------
  ! Fell off end of card without completing a word
  !---------------------------------------------------------------------
  LETCOL = LAST + 1
  IFFAIL = 1
  IF (MORE .GT. 0) RETURN
  WRITE(6, '(A/(1X,131A1))') 'NEXTWORD HAS READ OFF THE END', &
    (CARD(I), I = 1, LAST)
  STOP 1

END FUNCTION NEXTWORD

!=========================================================================
! SUBROUTINE READMOL
!
! Read molecular equilibrium data and build the component lookup tables.
!
! Reads molecular species definitions from molecules.dat, which lists
! each molecule with its species code and 6 equilibrium constants for
! the dissociation equilibrium calculation in NMOLEC.
!
! Species code encoding:
!   The code is a real number whose digits encode the atomic components.
!   Pairs of digits give atomic numbers, read left to right:
!     60808.    → atoms 6, 8, 8           → CO2
!     100.      → atom 1, charge code 100 → H-
!     26.       → atom 26                 → Fe I
!     101.01    → atoms 1, 1, ion code 01 → H2+
!     14.03     → atom 14, ion code 03    → Si IV (Si 3+)
!   The decimal part (.NN) gives the ionization state (NN electrons
!   removed), with each ion represented as component 101 (= electron).
!
! Outputs (via module variables):
!   NUMMOL              — number of molecules read
!   XNMOLCODE(JMOL)    — species code for molecule JMOL
!   EQUIL(1:6, JMOL)    — equilibrium constants for molecule JMOL
!   KCOMPS(KLOC)        — component list: equation indices for each atom
!                          in each molecule, stored contiguously
!   LOCJ(JMOL)          — index into KCOMPS where molecule JMOL starts
!                          LOCJ(JMOL+1)-LOCJ(JMOL) = number of components
!   IFEQUA(I)           — equation number assigned to element I
!                          (0 if element I not involved in any molecule)
!                          IFEQUA(100) = neutral atom, IFEQUA(101) = electron
!   IDEQUA(K)           — reverse map: element ID for equation K
!   NEQUA               — total number of equilibrium equations
!   NEQUA1              — NEQUA + 1 (convenience: index for inverse XNE)
!   NEQNEQ              — NEQUA^2 (convenience: size of rate matrix)
!
! Equilibrium equation structure:
!   Equation 1:      total particle number (XNATOM)
!   Equations 2..N:  one per element that appears in any molecule
!   Equation NEQUA:  charge conservation (if any ions present, i.e.
!                     if components 100 or 101 appear)
!   Variable NEQUA1: inverse electron density (1/XNE)
!
! File format (unit 2):
!   Each data line: A10, F12.2, 1X, F7.3, 5E12.4  (cols 1-90)
!     Column 1-10:  human-readable label  (e.g. 'OH', 'H2O', 'CaOH')
!                   This is for human readability only; the species code
!                   below is the source of truth for chemistry.
!     Column 11-22: species code (e.g. 60808.00)
!     Column 23:    space gap
!     Column 24-30: E1 (eV)
!     Column 31-90: E2 .. E6
!     Column 91+:   free-form trailing comment (ignored on read)
!   Lines starting with `!` (after optional leading whitespace) are
!     skipped as Fortran-style comments.
!   Blank lines are skipped.
!   End-of-file ends the molecule list. A legacy `0.0` terminator is
!     also still recognized for backward compatibility.
!   See the header at the top of molecules.dat for full file documentation.
!=========================================================================

SUBROUTINE READMOL

  IMPLICIT NONE

  ! --- Code extraction table ---
  ! XCODE(I) is the place value for the I-th pair of digits in the
  ! species code. Dividing the code by XCODE and truncating extracts
  ! the atomic number encoded at that position.
  REAL(8), PARAMETER :: XCODE(8) = &
    (/ 1.0D14, 1.0D12, 1.0D10, 1.0D8, 1.0D6, 1.0D4, 1.0D2, 1.0D0 /)

  ! --- Local variables ---
  REAL(8)  :: C                   ! species code read from file
  REAL(8)  :: E1, E2, E3, E4, E5, E6  ! equilibrium constants
  REAL(8)  :: X                   ! remaining code after extracting components
  INTEGER :: JMOL               ! molecule counter
  INTEGER :: KLOC               ! position in KCOMPS array
  INTEGER :: ID                  ! atomic number extracted from code
  INTEGER :: ION                 ! ionization state from decimal part
  INTEGER :: II                  ! starting position in XCODE for decoding
  INTEGER :: I, IEQUA, IOS, J
  CHARACTER(LEN=256) :: BUFFER  ! one input line, before parsing
  CHARACTER(LEN=10)  :: SPECIES_LABEL  ! human-readable label (cols 1-10)

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING READMOL'

  OPEN(UNIT=2, FILE=trim(DATADIR)//'molecules.dat', &
       STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS .NE. 0) THEN
    STOP 'molecules.dat file not found'
  END IF

  ! --- Initialize: no elements require equations yet ---
  IFEQUA(:) = 0

  ! --- Read molecules and build component table ---
  KLOC    = 1
  LOCJ(1) = 1

  JMOL = 0
  DO

    ! Read one line into a buffer; bail on EOF
    READ(2, '(A)', IOSTAT=IOS) BUFFER
    IF (IOS .NE. 0) EXIT

    ! Skip blank lines: find first non-space character
    J = 1
    DO WHILE (J .LE. LEN_TRIM(BUFFER) .AND. BUFFER(J:J) .EQ. ' ')
      J = J + 1
    END DO
    IF (J .GT. LEN_TRIM(BUFFER)) CYCLE

    ! Skip comment lines (Fortran-style `!`)
    IF (BUFFER(J:J) .EQ. '!') CYCLE

    ! Parse the data fields from the buffer using the fixed-width format.
    ! Layout: A10 label (skipped), F12.2 species code, 1-char gap, F7.3 E1,
    ! 5*E12.4 E2..E6
    READ(BUFFER, '(A10, F12.2, 1X, F7.3, 5E12.4)', IOSTAT=IOS) &
         SPECIES_LABEL, C, E1, E2, E3, E4, E5, E6
    IF (IOS .NE. 0) THEN
      WRITE(6, '(A)') ' READMOL ERROR: malformed data line:'
      WRITE(6, '(A,A)') '   ', trim(BUFFER)
      CLOSE(UNIT=2)
      STOP 1
    END IF

    ! Legacy `0.0` terminator (now optional): treat as end-of-list
    IF (C .EQ. 0.0D0) EXIT

    JMOL = JMOL + 1
    IF (JMOL .GT. maxmol) THEN
      WRITE(6, '(A,I5)') ' READMOL ERROR: molecule count exceeds maxmol =', maxmol
      CLOSE(UNIT=2)
      STOP 1
    END IF

    IF (IDEBUG .EQ. 1) WRITE(6, '(I5, F18.2, F7.3, 1P5E11.4)') &
         JMOL, C, E1, E2, E3, E4, E5, E6

    ! --- Decode species code into atomic components ---
    ! Find the leading non-zero pair of digits
    II = 0
    DO I = 1, 8
      IF (C .GE. XCODE(I)) THEN
        II = I
        EXIT
      END IF
    END DO
    IF (II .EQ. 0) THEN
      WRITE(6, '(A,F18.2)') ' READMOL ERROR: invalid species code ', C
      CLOSE(UNIT=2)
      STOP 1
    END IF

    ! Extract atomic numbers from successive digit pairs
    X = C
    DO I = II, 8
      ID = int(X / XCODE(I) + 0.5D0)
      X  = X - dble(ID) * XCODE(I)
      IF (ID .EQ. 0) ID = 100       ! code 00 → neutral atom (element 100)
      IF (ID .LT. 1 .OR. ID .GT. 101) THEN
        WRITE(6, '(A,I5,A,F18.2)') ' READMOL ERROR: decoded atomic ID =', ID, &
          ' out of range for species code ', C
        CLOSE(UNIT=2)
        STOP 1
      END IF
      IF (KLOC .GT. MAXLOC) THEN
        WRITE(6, '(A,I5)') ' READMOL ERROR: KCOMPS overflow, MAXLOC =', MAXLOC
        CLOSE(UNIT=2)
        STOP 1
      END IF
      IFEQUA(ID) = 1              ! flag: this element needs an equation
      KCOMPS(KLOC) = ID
      KLOC = KLOC + 1
    END DO

    ! Extract ionization state from decimal part: .NN → NN electrons removed
    ! Each removed electron adds component 101 (free electron)
    ION = int(X * 100.0D0 + 0.5D0)
    IF (ION .GE. 1) THEN
      IFEQUA(100) = 1    ! neutral atom appears in ion balance
      IFEQUA(101) = 1    ! electron appears in ion balance
      DO I = 1, ION
        IF (KLOC .GT. MAXLOC) THEN
          WRITE(6, '(A,I5)') ' READMOL ERROR: KCOMPS overflow, MAXLOC =', MAXLOC
          CLOSE(UNIT=2)
          STOP 1
        END IF
        KCOMPS(KLOC) = 101
        KLOC = KLOC + 1
      END DO
    END IF

    ! Store molecule data
    LOCJ(JMOL + 1) = KLOC
    XNMOLCODE(JMOL) = C
    EQUIL(1:6, JMOL) = [E1, E2, E3, E4, E5, E6]

  END DO

  NUMMOL = JMOL
  NLOC   = KLOC - 1

  ! --- Assign equation numbers to each element ---
  ! Equation 1 is for total particle number (XNATOM).
  ! Subsequent equations are assigned in order of element number.
  ! If any ions are present, the last equation is charge conservation
  ! (variable NEQUA = XNE, variable NEQUA1 = 1/XNE).
  IEQUA = 1
  DO I = 1, 100
    IF (IFEQUA(I) .EQ. 0) CYCLE
    IEQUA = IEQUA + 1
    IF (IEQUA .GT. maxeq) THEN
      WRITE(6, '(A,I5)') ' READMOL ERROR: too many equations, maxeq =', maxeq
      CLOSE(UNIT=2)
      STOP 1
    END IF
    IFEQUA(I) = IEQUA        ! element I → equation number IEQUA
    IDEQUA(IEQUA) = I        ! equation IEQUA → element I
  END DO

  NEQUA  = IEQUA
  NEQUA1 = NEQUA + 1
  IFEQUA(101) = NEQUA1       ! electron → variable index NEQUA+1
  NEQNEQ = NEQUA**2

  ! --- Remap KCOMPS from element IDs to equation indices ---
  ! After this, KCOMPS(K) gives the equation number for the K-th
  ! component, ready for direct use in building the rate matrix.
  DO I = 1, NLOC
    ID = KCOMPS(I)
    KCOMPS(I) = IFEQUA(ID)
  END DO

  IF (IDEBUG .EQ. 1) THEN
    WRITE(6,*)
    WRITE(6, '(A,I4,A,I4)') ' MOLECULES  USED', NUMMOL, '  MAX', MAXMOL
    WRITE(6, '(A,I4,A,I4)') ' COMPONENTS USED', NLOC,   '  MAX', MAXLOC
    WRITE(6, '(A,I4,A,I4)') ' EQUATIONS  USED', NEQUA,  '  MAX', MAXEQ
  END IF

  CLOSE(UNIT=2)
  RETURN

END SUBROUTINE READMOL

!=========================================================================
! SUBROUTINE COMPUTE_ONE_POP(CODE, MODE, NUMBER)
!
! Population dispatcher: compute number densities for any atomic or
! molecular species at all depth points.
!
! This is the central interface for obtaining populations throughout
! ATLAS12. It dispatches to the appropriate solver depending on
! whether molecular equilibrium is active (IFMOL) and what species
! is requested (CODE).
!
! Species code encoding (same convention as READMOL):
!   CODE = 0.0       → no species requested; just ensure electron
!                       density is up to date (NELECT or NMOLEC)
!   CODE < 100       → atomic species:
!                       IZ   = INT(CODE)          = atomic number
!                       NION = decimal part * 100  = number of ion stages
!                       e.g. 1.01 = H I (1 stage)
!                            26.05 = Fe I-V (5 stages)
!                       Calls PFSAHA for Saha ionization at each depth,
!                       then multiplies by total number density and
!                       abundance to get absolute number densities.
!   CODE >= 100      → molecular species (requires IFMOL=1):
!                       Dispatched to MOLEC, which looks up the species
!                       in the molecular equilibrium tables built by
!                       NMOLEC.
!                       e.g. 101.00 = H2, 60808.00 = CO2
!
! MODE controls what PFSAHA/MOLEC returns:
!   MODE = 1         → ionization fractions (n_ion / n_total)
!   MODE = 11 or 12  → ionization fractions / partition functions
!                       (= number density per unit partition function,
!                        needed for opacity calculations)
!   For MODE < 10, only the ground ionization stage is converted to
!   number density. For MODE >= 10, all NION stages are converted.
!
! Electron density caching:
!   NELECT (atomic) or NMOLEC (molecular) is called only once per
!   temperature iteration, tracked by comparing ITEMP to ITEMP_PREV.
!   Subsequent COMPUTE_ONE_POP calls within the same iteration reuse the cached
!   electron densities.
!
! Arguments:
!   CODE             — species code (see encoding above)
!   MODE             — output mode (1, 11, or 12)
!   NUMBER(kw,*)     — output: number densities at each depth
!                       For atoms: NUMBER(J,ION) for ION=1..NION
!                       For molecules: NUMBER(J,1)
!=========================================================================

SUBROUTINE COMPUTE_ONE_POP(CODE, MODE, NUMBER)

  IMPLICIT NONE

  ! --- Arguments ---
  REAL(8),  INTENT(IN)    :: CODE
  INTEGER, INTENT(IN)    :: MODE
  REAL(8),  INTENT(INOUT) :: NUMBER(kw, *)

  ! --- Persistent: tracks whether electron density is current ---
  INTEGER, SAVE :: ITEMP_PREV = 0

  ! --- Local variables ---
  INTEGER :: IZ       ! atomic number
  INTEGER :: NION     ! number of ionization stages requested
  INTEGER :: NNNN     ! number of stages to convert to number density
  INTEGER :: J, ION

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COMPUTE_ONE_POP'

  !=====================================================================
  ! Branch: molecular equilibrium ON or OFF
  !=====================================================================
  IF (IFMOL .EQ. 1) THEN

    ! --- Molecular equilibrium path ---
    ! NMOLEC solves the full molecular + ionization equilibrium system,
    ! including electron density, in one shot.
    IF (IFPRES .EQ. 1 .AND. ITEMP .NE. ITEMP_PREV) CALL NMOLEC(MODE)
    ITEMP_PREV = ITEMP

    IF (CODE .EQ. 0.0D0) RETURN

    ! Look up requested species from molecular equilibrium tables
    CALL MOLEC(CODE, MODE, NUMBER)

  ELSE

    ! --- Atomic-only path (no molecules) ---
    ! NELECT solves for electron density from Saha ionization alone.
    IF (IFPRES .EQ. 1 .AND. ITEMP .NE. ITEMP_PREV) CALL NELECT
    ITEMP_PREV = ITEMP

    IF (CODE .EQ. 0.0D0) RETURN

    IF (CODE .GE. 100.0D0) THEN
      ! Molecular species requested but molecules are off
      WRITE(6, '(A)') ' COMPUTE_ONE_POP ERROR: molecular species requested but IFMOL=0'
      STOP 1
    END IF

    ! --- Decode species code ---
    IZ   = int(CODE)                              ! atomic number
    NION = int((CODE - dble(IZ)) * 100.0D0 + 1.5D0)  ! number of ion stages

    ! --- Compute ionization fractions at each depth via Saha equation ---
    DO J = 1, NRHOX
      CALL PFSAHA(J, IZ, NION, MODE, NUMBER)

      ! Convert fractions to absolute number densities:
      !   n_ion = fraction * n_total_atoms * abundance_of_element
      ! For MODE < 10: only the first stage (ground ionization)
      ! For MODE >= 10: all NION stages
      NNNN = 1
      IF (MODE .GE. 10) NNNN = NION

      DO ION = 1, NNNN
        NUMBER(J, ION) = NUMBER(J, ION) * XNATOM(J) * XABUND(J, IZ)
      END DO
    END DO

  END IF

  RETURN

END SUBROUTINE COMPUTE_ONE_POP

!=========================================================================
! SUBROUTINE ENERGY_DENSITY
!
! Compute the internal energy density at each depth point.
!
! The internal energy per unit mass has two contributions:
!   1. Kinetic (translational): (3/2) * n_total * kT
!   2. Ionization + excitation: for each ion stage,
!        n_ion * kT * [ chi_cumul * (hc/kT) + d(ln U)/d(ln T) ]
!      where chi_cumul is the cumulative ionization potential up to
!      that stage, and U is the partition function.
!
! The d(ln U)/d(ln T) term captures the energy stored in excited
! states — as T increases, higher levels become populated, absorbing
! energy even without ionization. It is evaluated by finite differences:
!   d(ln U)/d(ln T) ≈ [U(T+) - U(T-)] / [U(T+) + U(T-)] * (2/dlnT)
! with dlnT = ln(1.001) - ln(0.999) ≈ 0.002, giving a scale factor
! of 2/0.002 = 1000.
!
! The energy density is needed by CONVEC for the adiabatic gradient
! and specific heat in mixing length theory. IFEDNS is set to 1 by
! CONVEC when it needs this quantity.
!
! Note: Currently atoms only (Z=1-99), not molecules.
!=========================================================================

SUBROUTINE ENERGY_DENSITY

  IMPLICIT NONE

  ! --- Temperature perturbation for finite differences ---
  REAL(8), PARAMETER :: DELTA_P     = 1.001D0     ! T+ multiplier
  REAL(8), PARAMETER :: DELTA_M     = 0.999D0     ! T- multiplier
  REAL(8), PARAMETER :: RATIO_PM    = 0.999D0 / 1.001D0  ! T- / T+ (to go from T+ to T-)
  REAL(8), PARAMETER :: DLNUDT_SCALE = 1000.0D0   ! 2 / (ln(1.001) - ln(0.999))

  ! --- Number of ionization stages per element (for PFSAHA calls) ---
  INTEGER, PARAMETER :: NION_Z(30) = (/ &
    2, 3, 4, 5, 5,                        &  ! Z = 1-5
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,     &  ! Z = 6-16
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  &  ! Z = 17-28
    3, 3 /)                                   ! Z = 29-30
  ! Z = 31-99: all use NION = 3

  ! --- Local arrays ---
  REAL(8)  :: PFPLUS(mion)    ! partition functions at T+
  REAL(8)  :: PFMINUS(mion)   ! partition functions at T-

  ! --- Local variables ---
  REAL(8)  :: XNTOT           ! total particle density (atoms + electrons)
  REAL(8)  :: UPF, DLNU       ! partition function sum and log-derivative
  INTEGER :: J, IION, NOFF

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING ENERGY_DENSITY'
  IF (IFEDNS .EQ. 0) RETURN

  DO J = 1, NRHOX

    ! --- Compute partition functions at T+ = T * 1.001 ---
    CALL perturb_T(J, DELTA_P)
    CALL compute_partfcns(J, PFPLUS)

    ! --- Compute partition functions at T- = T * 0.999 ---
    ! T is currently at T*1.001, so multiply by 0.999/1.001
    CALL perturb_T(J, RATIO_PM)
    CALL compute_partfcns(J, PFMINUS)

    ! --- Restore original temperature ---
    ! T is currently at T*0.999, so divide by 0.999
    CALL perturb_T(J, 1.0D0 / DELTA_M)

    ! --- Assemble energy density ---
    ! Kinetic: (3/2) n_total kT
    XNTOT = XNE(J) + XNATOM(J)
    EDENS(J) = 1.5D0 * XNTOT * TK(J)

    ! Ionization + excitation for each ion stage (IION = 1..840)
    DO IION = 1, 840
      UPF  = PFPLUS(IION) + PFMINUS(IION) + 1.0D-30
      DLNU = (PFPLUS(IION) - PFMINUS(IION)) / UPF * DLNUDT_SCALE

      EDENS(J) = EDENS(J) + XNF(J, IION) * TK(J) &
               * (POTIONSUM(IION) * HCKT(J) + DLNU)
    END DO

    ! Convert to energy per unit mass
    EDENS(J) = EDENS(J) / RHO(J)

  END DO

  RETURN

CONTAINS

  !---------------------------------------------------------------------
  ! Multiply temperature and related arrays at depth J by factor FAC.
  !---------------------------------------------------------------------
  SUBROUTINE perturb_T(J_in, FAC)
    INTEGER, INTENT(IN) :: J_in
    REAL(8),  INTENT(IN) :: FAC
    REAL(8) :: FAC_INV
    FAC_INV = 1.0D0 / FAC
    T(J_in)    = T(J_in) * FAC
    TK(J_in)   = TK(J_in) * FAC
    TKEV(J_in) = TKEV(J_in) * FAC
    HCKT(J_in) = HCKT(J_in) * FAC_INV
    HKT(J_in)  = HKT(J_in) * FAC_INV
    TLOG(J_in) = log(T(J_in))
  END SUBROUTINE perturb_T

  !---------------------------------------------------------------------
  ! Compute partition functions for all elements (Z=1-99) at depth J
  ! via PFSAHA MODE 5. Results stored in PF(:) at the standard NELION
  ! offsets: IZ*(IZ+1)/2 for Z=1-30, 496+(IZ-31)*5 for Z=31-99.
  !---------------------------------------------------------------------
  SUBROUTINE compute_partfcns(J_in, PF)
    INTEGER, INTENT(IN)  :: J_in
    REAL(8),  INTENT(OUT) :: PF(mion)
    INTEGER :: IZ_loc

    DO IZ_loc = 1, 30
      NOFF = IZ_loc * (IZ_loc + 1) / 2
      CALL PFSAHA(J_in, IZ_loc, NION_Z(IZ_loc), 5, PF(NOFF))
    END DO
    DO IZ_loc = 31, 99
      NOFF = 496 + (IZ_loc - 31) * 5
      CALL PFSAHA(J_in, IZ_loc, 3, 5, PF(NOFF))
    END DO
  END SUBROUTINE compute_partfcns

END SUBROUTINE ENERGY_DENSITY

!=========================================================================
! SUBROUTINE NELECT
!
! Compute electron density by iterative charge conservation.
!
! At each depth point, solves for the electron density XNE that
! satisfies simultaneous charge conservation and the Saha ionization
! equations for all elements Z=1-99. This is the atomic-only path
! (IFMOL=0); when molecules are active, NMOLEC handles this instead.
!
! Algorithm:
!   Given total pressure P and temperature T, the total particle
!   density is n_total = P / kT. The electron density XNE enters
!   the Saha equation (via PFSAHA), which returns ionization fractions
!   for each element. The charge-weighted sum of all ion populations
!   gives a new electron density estimate. Iteration continues until
!   XNE converges to 0.01%.
!
!   At each iteration:
!   1. n_atom = n_total - n_e  (atoms = total minus electrons)
!   2. For each element, call PFSAHA(MODE=12) to get ionization
!      fractions divided by partition functions
!   3. Multiply by n_atom * abundance to get number densities
!   4. Sum charge-weighted populations: n_e(new) = sum(n_ion * q_ion)
!   5. Average: n_e = (n_e(old) + n_e(new)) / 2
!   6. Check convergence
!
! Also computes:
!   XNATOM(J)   — total atom number density
!   XNF(J,*)    — number densities for all ion stages of all elements
!   RHO(J)      — mass density = n_atom * mean_molecular_weight * m_H
!   CHARGESQ(J) — sum of q^2 * n_ion (for Debye shielding)
!
! Convergence: typically 5-15 iterations; max 200 with warning.
!=========================================================================

SUBROUTINE NELECT

  IMPLICIT NONE

  ! --- Convergence parameters ---
  INTEGER, PARAMETER :: MAX_ITER_XNE = 200
  REAL(8),  PARAMETER :: XNE_TOL      = 1.0D-4   ! relative convergence

  ! --- Number of ionization stages per element for Saha computation ---
  ! These differ slightly from COMPUTE_ALL_POPS because NELECT needs
  ! enough stages to capture charge balance accurately.
  INTEGER, PARAMETER :: NION_NELECT(30) = (/ &
    2, 3, 4, 4, 4,                        &  ! Z = 1-5
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,     &  ! Z = 6-16
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  &  ! Z = 17-28
    3, 3 /)                                   ! Z = 29-30
  ! Z = 31-99: all use NION = 3

  ! --- Local variables ---
  REAL(8)  :: XNTOT           ! total particle density = P / kT
  REAL(8)  :: XNENEW          ! new electron density estimate
  REAL(8)  :: CHARGESQUARE    ! sum of n_ion * q^2
  REAL(8)  :: ERROR           ! relative change in XNE
  REAL(8)  :: DEBYE           ! Debye shielding length
  INTEGER :: J, IZ, ION, IION, NOFF, ITER_XNE

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING NELECT'

  !---------------------------------------------------------------------
  ! Loop over all depth points
  !---------------------------------------------------------------------
  DO J = 1, NRHOX

    ! Total particle density from ideal gas law
    XNTOT = P(J) / TK(J)
    XNATOM(J) = XNTOT - XNE(J)

    !-------------------------------------------------------------------
    ! Iterate to converge electron density
    !-------------------------------------------------------------------
    DO ITER_XNE = 1, MAX_ITER_XNE

      ! Zero all ion populations
      XNF(J, 1:MION) = 0.0D0

      XNENEW = 0.0D0
      CHARGESQUARE = 0.0D0

      ! --- Compute Saha ionization for elements Z=1-30 ---
      DO IZ = 1, 30
        NOFF = IZ * (IZ + 1) / 2
        CALL PFSAHA(J, IZ, NION_NELECT(IZ), 12, XNF(1, NOFF))
      END DO

      ! --- Accumulate electron density and charge sum for Z=1-30 ---
      ! Each element IZ has IZ+1 ion stages in the NELION layout
      IION = 0
      DO IZ = 1, 30
        DO ION = 1, IZ + 1
          IION = IION + 1
          ! Convert fraction to number density
          XNF(J, IION) = XNF(J, IION) * XNATOM(J) * XABUND(J, IZ)
          ! Charge = ION - 1 (neutral = 0, singly ionized = 1, etc.)
          CHARGESQUARE = CHARGESQUARE + XNF(J, IION) * dble((ION - 1)**2)
          XNENEW = XNENEW + XNF(J, IION) * dble(ION - 1)
        END DO
      END DO

      ! --- Elements Z=31-99: 3 stages computed, 5 slots allocated ---
      DO IZ = 31, 99
        NOFF = 496 + (IZ - 31) * 5
        CALL PFSAHA(J, IZ, 3, 12, XNF(1, NOFF))
        DO ION = 1, 5
          IION = IION + 1
          XNF(J, IION) = XNF(J, IION) * XNATOM(J) * XABUND(J, IZ)
          CHARGESQUARE = CHARGESQUARE + XNF(J, IION) * dble((ION - 1)**2)
          XNENEW = XNENEW + XNF(J, IION) * dble(ION - 1)
        END DO
      END DO

      ! Floor to prevent oscillation to zero
      XNENEW = max(XNENEW, XNE(J) * 0.5D0)

      ! Debye length (for diagnostic; also sets ionization potential
      ! lowering inside PFSAHA on the next iteration)
      DEBYE = sqrt(TK(J) / (FOURPI * ECHARGE**2 &
                             * (CHARGESQUARE + XNENEW)))

      ! Damped update: average old and new
      XNENEW = (XNENEW + XNE(J)) * 0.5D0

      ! Check convergence
      ERROR = abs((XNE(J) - XNENEW) / XNENEW)
      XNE(J) = XNENEW
      XNATOM(J) = XNTOT - XNE(J)
      CHARGESQ(J) = CHARGESQUARE + XNE(J)

      IF (ERROR .LT. XNE_TOL) EXIT

    END DO  ! ITER_XNE

    IF (ITER_XNE .GT. MAX_ITER_XNE) THEN
      WRITE(6, '(A,I4,A,ES10.3)') &
        ' NELECT WARNING: XNE did not converge at J=', J, ', error=', ERROR
    END IF

    ! Mass density from atom count and mean molecular weight
    RHO(J) = XNATOM(J) * WTMOLE(J) * AMU

  END DO  ! depth loop

  ! --- Debug output ---
  IF (IDEBUG .EQ. 1) THEN
    WRITE(6, '(3X,4X,A,A,A,A,A,A,A,A)') &
      'RHOX', '     T    ', '     P     ', '     XNE    ', &
      '    XNATOM  ', '    WTMOLE  ', '     RHO    ', '    CHARGESQ'
    DO J = 1, NRHOX
      WRITE(6, '(I3,1PE15.7,0PF10.1,1P6E12.3)') &
        J, RHOX(J), T(J), P(J), XNE(J), XNATOM(J), WTMOLE(J), &
        RHO(J), CHARGESQ(J)
    END DO
  END IF

  RETURN

END SUBROUTINE NELECT

!=======================================================================
! SUBROUTINE PFSAHA
!
! Partition functions and Saha equation solver.
!
! Computes partition functions for each ionization stage by:
!   - Detailed level summation for H, He, B, C, N, O, Na, Mg, Al, Si, K, Ca
!   - NNN interpolation table for generic light elements
!   - PFIRON tabulated partition functions for iron group (Z=20-28)
!
! MODE controls output in ANSWER:
!   MODE 1  or 11: number density / partition function (per-state populations)
!   MODE 2  or 12: ionization fractions F(ION)
!   MODE 3  or 13: partition functions PART(ION) only
!   MODE 4  or 14: mean electron number
!   MODE 5:        partition functions + cumulative IPs (for COMPUTE_ALL_POPS)
!   MODE < 10: returns for all ions; MODE >= 10: returns for NION-th ion only

!=======================================================================

SUBROUTINE PFSAHA(J, IZ, NION, MODE, ANSWER)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)    :: J       ! Depth point index
  INTEGER, INTENT(IN)    :: IZ      ! Atomic number
  INTEGER, INTENT(IN)    :: NION    ! Number of ionization stages requested
  INTEGER, INTENT(IN)    :: MODE    ! Output mode selector
  REAL(8),  INTENT(INOUT) :: ANSWER(kw, *)

  ! External function

  ! NNN interpolation table (loaded from file on first call)
  INTEGER, SAVE :: NNN(6, 365)
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! Local arrays
  REAL(8)  :: F(31), IP(31), PART(31), POTLO(31)
  INTEGER :: LOCZ(29)

  ! Scale factors for NNN unpacking
  REAL(8) :: SCALE(4)

  ! Detailed energy level tables for special elements (cm^-1)
  REAL(8) :: EHYD(6), GHYD(6)
  REAL(8) :: EHE1(29), GHE1(29), EHE2(6), GHE2(6)
  REAL(8) :: EB1(7), GB1(7)
  REAL(8) :: EC1(14), GC1(14), EC2(6), GC2(6)
  REAL(8) :: EO1(13), GO1(13)
  REAL(8) :: ENA1(8), GNA1(8)
  REAL(8) :: EMG1(11), GMG1(11), EMG2(6), GMG2(6)
  REAL(8) :: EAL1(9), GAL1(9)
  REAL(8) :: ESI1(11), GSI1(11), ESI2(6), GSI2(6)
  REAL(8) :: EK1(8), GK1(8)
  REAL(8) :: ECA1(8), GCA1(8), ECA2(5), GCA2(5)

  ! Local scalars
  REAL(8)  :: Z_ion, TV, DEBYE, POTLOW, T2000, DT, P1, P2, PMIN
  REAL(8)  :: G, D1, D2, CF, B_dep
  INTEGER :: MODE1, N, NIONS, NION2, ION, I, L, IT, NNN100, INDEX
  INTEGER :: K1, K2, K3, KSCALE, KP1
  INTEGER :: jj  ! loop variable for data read

  ! Data initialization
  DATA LOCZ /1,3,6,10,14,18,22,27,33,39,45,51,57,63,69,75,81,86,91, &
             96,101,106,111,116,121,126,131,136,141/
  DATA SCALE /0.001d0, 0.01d0, 0.1d0, 1.0d0/

  ! Hydrogen (6 levels)
  DATA EHYD /0.d0, 82259.105d0, 97492.302d0, 102823.893d0, 105291.651d0, 106632.160d0/
  DATA GHYD /2.d0, 8.d0, 18.d0, 32.d0, 50.d0, 72.d0/

  ! Helium I (29 levels)
  DATA EHE1 / 0.d0, 159856.069d0, 166277.546d0, 169087.007d0, 171135.000d0, &
    183236.892d0, 184864.936d0, 185564.694d0, 186101.654d0, 186105.065d0, &
    186209.471d0, 190298.210d0, 190940.331d0, 191217.14d0, 191444.588d0, &
    191446.559d0, 191451.80d0, 191452.08d0, 191492.817d0, 193347.089d0, &
    193663.627d0, 193800.78d0, 193917.245d0, 193918.391d0, 193921.31d0, &
    193921.37d0, 193922.5d0, 193922.5d0, 193942.57d0/
  DATA GHE1 /1.d0,3.d0,1.d0,9.d0,3.d0,3.d0,1.d0,9.d0,15.d0,5.d0, &
    3.d0,3.d0,1.d0,9.d0,15.d0,5.d0,21.d0,7.d0,3.d0,3.d0, &
    1.d0,9.d0,15.d0,5.d0,21.d0,7.d0,27.d0,9.d0,3.d0/

  ! Helium II (6 levels, hydrogenic)
  DATA EHE2 /0.d0, 329182.321d0, 390142.359d0, 411477.925d0, 421353.135d0, 426717.413d0/
  DATA GHE2 /2.d0, 8.d0, 18.d0, 32.d0, 50.d0, 72.d0/

  ! Boron I (7 levels)
  DATA EB1  /10.17d0, 28810.d0, 40039.65d0, 47856.99d0, 48613.01d0, 54767.74d0, 55010.08d0/
  DATA GB1  /6.d0, 12.d0, 2.d0, 10.d0, 6.d0, 10.d0, 2.d0/

  ! Carbon I (14 levels)
  DATA EC1  / 29.60d0, 10192.66d0, 21648.02d0, 33735.20d0, 60373.00d0, 61981.82d0, &
    64088.85d0, 68856.33d0, 69722.00d0, 70743.95d0, 71374.90d0, 72610.72d0, &
    73975.91d0, 75254.93d0/
  DATA GC1  /9.d0,5.d0,1.d0,5.d0,9.d0,3.d0,15.d0,3.d0,15.d0,3.d0, &
    9.d0,5.d0,1.d0,9.d0/

  ! Carbon II (6 levels)
  DATA EC2  /42.48d0, 43035.8d0, 74931.11d0, 96493.74d0, 110652.10d0, 116537.65d0/
  DATA GC2  /6.d0, 12.d0, 10.d0, 2.d0, 6.d0, 2.d0/

  ! Oxygen I (13 levels)
  DATA EO1  / 77.975d0, 15867.862d0, 33792.583d0, 73768.200d0, 76794.978d0, &
    86629.089d0, 88630.977d0, 95476.728d0, 96225.049d0, 97420.748d0, &
    97488.476d0, 99094.065d0, 99681.051d0/
  DATA GO1  /9.d0,5.d0,1.d0,5.d0,3.d0,15.d0,9.d0,5.d0,3.d0,25.d0, &
    15.d0,15.d0,9.d0/

  ! Sodium I (8 levels)
  DATA ENA1 /0.d0, 16956.172d0, 16973.368d0, 25739.991d0, 29172.889d0, &
    29172.839d0, 30266.99d0, 30272.58d0/
  DATA GNA1 /2.d0, 2.d0, 4.d0, 2.d0, 6.d0, 4.d0, 2.d0, 4.d0/

  ! Magnesium I (11 levels)
  DATA EMG1 /0.d0, 21890.854d0, 35051.264d0, 41197.403d0, 43503.333d0, &
    46403.065d0, 47847.797d0, 47957.034d0, 49346.729d0, 51872.526d0, 52556.206d0/
  DATA GMG1 /1.d0,9.d0,3.d0,3.d0,1.d0,5.d0,9.d0,15.d0,3.d0,3.d0,1.d0/

  ! Magnesium II (6 levels)
  DATA EMG2 /0.d0, 35730.36d0, 69804.95d0, 71490.54d0, 80639.85d0, 92790.51d0/
  DATA GMG2 /2.d0, 6.d0, 2.d0, 10.d0, 6.d0, 2.d0/

  ! Aluminum I (9 levels)
  DATA EAL1 /74.707d0, 25347.756d0, 29097.11d0, 32436.241d0, 32960.363d0, &
    37689.413d0, 38932.139d0, 40275.903d0, 41319.377d0/
  DATA GAL1 /6.d0, 2.d0, 12.d0, 10.d0, 6.d0, 2.d0, 10.d0, 6.d0, 14.d0/

  ! Silicon I (11 levels)
  DATA ESI1 /149.681d0, 6298.850d0, 15394.370d0, 33326.053d0, 39859.920d0, &
    40991.884d0, 45303.310d0, 47284.061d0, 47351.554d0, 48161.459d0, 49128.131d0/
  DATA GSI1 /9.d0,5.d0,1.d0,5.d0,9.d0,3.d0,15.d0,3.d0,5.d0,15.d0,9.d0/

  ! Silicon II (6 levels)
  DATA ESI2 /191.55d0, 43002.27d0, 55319.11d0, 65500.73d0, 76665.61d0, 79348.67d0/
  DATA GSI2 /6.d0, 12.d0, 10.d0, 2.d0, 2.d0, 10.d0/

  ! Potassium I (8 levels)
  DATA EK1  /0.d0, 12985.170d0, 13042.876d0, 21026.551d0, 21534.680d0, &
    21536.988d0, 24701.382d0, 24720.139d0/
  DATA GK1  /2.d0, 2.d0, 4.d0, 2.d0, 6.d0, 4.d0, 2.d0, 4.d0/

  ! Calcium I (8 levels)
  DATA ECA1 /0.d0, 15263.089d0, 20356.265d0, 21849.634d0, 23652.304d0, &
    31539.495d0, 33317.264d0, 35831.203d0/
  DATA GCA1 /1.d0, 9.d0, 15.d0, 5.d0, 3.d0, 3.d0, 1.d0, 21.d0/

  ! Calcium II (5 levels)
  DATA ECA2 /0.d0, 13686.60d0, 25340.10d0, 52166.93d0, 56850.78d0/
  DATA GCA2 /2.d0, 10.d0, 6.d0, 2.d0, 10.d0/

  !=====================================================================
  ! Load NNN table on first call
  !=====================================================================
  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING PFSAHA'
  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=TRIM(DATADIR)//'pfsaha.dat', STATUS='OLD', ACTION='READ')
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')  ! skip 3 header lines
    READ(89, *) ((NNN(I,jj), I=1,6), jj=1,365)
    CLOSE(89)
    ! Ensure ionization potentials are loaded (PFSAHA needs POTION for IP)
    CALL IONPOTS
    INITIALIZED = .TRUE.
  ENDIF

  !=====================================================================
  ! Setup: determine NNN row range and number of ions for this element
  !=====================================================================
  MODE1 = MODE
  IF (MODE1 .GT. 10) MODE1 = MODE1 - 10

  ! Debye-Hückel lowering of ionization potential (volts per unit Zeff)
  DEBYE  = SQRT(TK(J) / FOURPI / ECHARGE**2 / CHARGESQ(J))
  POTLOW = MIN(1.D0, 1.44D-7 / DEBYE)
  TV     = TKEV(J)

  ! Locate element in NNN table
  IF (IZ .LE. 28) THEN
    N     = LOCZ(IZ)
    NIONS = LOCZ(IZ+1) - N
  ELSE
    N     = 3*IZ + 54
    NIONS = 3
  ENDIF

  ! Override for C and N (extended tables in NNN67 block)
  IF (IZ .EQ. 6) THEN; N = 354; NIONS = 6; ENDIF
  IF (IZ .EQ. 7) THEN; N = 360; NIONS = 6; ENDIF

  ! Iron group elements have 10 ionization stages in PFIRON
  IF (IZ .GE. 20 .AND. IZ .LT. 29) NIONS = 10

  NION2 = MIN(NION + 2, NIONS)
  N = N - 1   ! will be incremented at start of loop

  !=====================================================================
  ! Phase 1: Compute partition functions PART(ION) for each ion stage
  !=====================================================================
  DO ION = 1, NION2
    Z_ion = dble(ION)
    POTLO(ION) = POTLOW * Z_ion
    N = N + 1

    ! Ionization potential from POTION table
    NNN100 = NNN(6,N) / 100
    IF (IZ .LE. 30) THEN
      INDEX = (IZ*(IZ+1))/2 + ION - 1
    ELSE
      INDEX = IZ*5 + 341 + ION - 1
    ENDIF
    IP(ION) = POTION(INDEX) / 8065.479d0
    IF (IP(ION) .EQ. 0.d0) IP(ION) = POTION(INDEX-1) / 8065.479d0

    ! Iron group: use PFIRON for partition functions
    IF (IZ .GE. 20 .AND. IZ .LT. 29) THEN
      CALL PFIRON(IZ, ION, TLOG(J)/LN10, &
                  POTLO(ION)*8065.479d0, PART(ION))
      CYCLE  ! next ION
    ENDIF

    ! Statistical weight of highest included level
    G = NNN(6,N) - NNN100*100

    !-----------------------------------------------------------------
    ! Detailed level summation for special elements
    ! B_dep = NLTE departure coefficient (always 1.0 in LTE mode)
    !-----------------------------------------------------------------
    B_dep = 1.0d0

    SELECT CASE (N)

    CASE (1)  ! ---- Hydrogen I ----
      ! Partition function with Hummer & Mihalas (1988) occupation
      ! probability weighting.  Each level n is weighted by w_n, the
      ! probability that the level survives Stark dissolution.
      ! For FGK stars (T < 8000 K), the ground state dominates and
      ! this gives PF ≈ 2.000 regardless of density.  For hotter stars
      ! the weighting properly truncates the sum that would otherwise
      ! diverge from unweighted high-n levels.
      !
      ! Replaces the former approach of levels 1-6 at full weight plus
      ! a fixed Euler-Maclaurin integral from n=6.5 to ionization.
      PART(1) = 2.d0 * B_dep       ! n=1: w_1 = 1 always
      DO I = 2, 6
        PART(1) = PART(1) + occupation_prob(I, XNE(J)) &
                * GHYD(I) * B_dep * EXP(-EHYD(I)*HCKT(J))
      END DO
      ! Levels n=7 and above: hydrogenic energies, weighted by w_n.
      ! Sum truncates naturally where w_n -> 0 or Boltzmann factor -> 0.
      DO I = 7, 80
        D1 = occupation_prob(I, XNE(J))
        IF (D1 .LT. 1.0d-6) EXIT   ! remaining levels fully dissolved
        PART(1) = PART(1) + D1 &
                * 2.0d0 * dble(I)**2 * EXP(-(ELIM_HI - RYDBERG_H/dble(I)**2)*HCKT(J))
      END DO
      CYCLE

    CASE (3)  ! ---- Helium I ----
      PART(1) = B_dep
      DO I = 2, 29
        PART(1) = PART(1) + GHE1(I) * B_dep * EXP(-EHE1(I)*HCKT(J))
      END DO
      D1 = RYDBERG_H / 5.5d0 / 5.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (4)  ! ---- Helium II ----
      PART(2) = 2.d0 * B_dep
      DO I = 2, 6
        PART(2) = PART(2) + GHE2(I) * B_dep * EXP(-EHE2(I)*HCKT(J))
      END DO
      D1 = 4.d0 * RYDBERG_HE / 6.5d0 / 6.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (14)  ! ---- Boron I ----
      PART(1) = B_dep * (2.d0 + 4.d0*EXP(-15.25d0*HCKT(J)))
      DO I = 2, 7
        PART(1) = PART(1) + GB1(I) * B_dep * EXP(-EB1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 6.d0*EXP(-57786.80d0*HCKT(J)) &
        + 10.d0*EXP(-59989.d0*HCKT(J)) + 14.d0*EXP(-60031.03d0*HCKT(J)) &
        + 2.d0*EXP(-63561.d0*HCKT(J))
      G = 2.d0
      D1 = 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (45)  ! ---- Sodium I ----
      PART(1) = B_dep * 2.d0
      DO I = 2, 8
        PART(1) = PART(1) + GNA1(I) * B_dep * EXP(-ENA1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 10.d0*EXP(-34548.745d0*HCKT(J)) &
        + 14.d0*EXP(-34586.96d0*HCKT(J))
      G = 2.d0
      D1 = 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (51)  ! ---- Magnesium I ----
      PART(1) = B_dep
      DO I = 2, 11
        PART(1) = PART(1) + GMG1(I) * B_dep * EXP(-EMG1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 5.d0*EXP(-53134.d0*HCKT(J)) &
        + 15.d0*EXP(-54192.d0*HCKT(J)) + 28.d0*EXP(-54676.d0*HCKT(J)) &
        + 9.d0*EXP(-57853.d0*HCKT(J))
      G = 4.d0
      D1 = 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (52)  ! ---- Magnesium II ----
      PART(2) = B_dep * 2.d0
      DO I = 2, 6
        PART(2) = PART(2) + GMG2(I) * B_dep * EXP(-EMG2(I)*HCKT(J))
      END DO
      PART(2) = PART(2) + 10.d0*EXP(-93310.80d0*HCKT(J)) &
        + 14.d0*EXP(-93799.70d0*HCKT(J)) + 6.d0*EXP(-97464.32d0*HCKT(J)) &
        + 10.d0*EXP(-103419.82d0*HCKT(J)) + 14.d0*EXP(-103689.89d0*HCKT(J)) &
        + 18.d0*EXP(-103705.66d0*HCKT(J))
      G = 2.d0
      D1 = 4.d0 * 109734.83d0 / 5.5d0 / 5.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (57)  ! ---- Aluminum I ----
      PART(1) = B_dep * (2.d0 + 4.d0*EXP(-112.061d0*HCKT(J)))
      DO I = 2, 9
        PART(1) = PART(1) + GAL1(I) * B_dep * EXP(-EAL1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 10.d0*EXP(-42235.d0*HCKT(J)) &
        + 14.d0*EXP(-43831.d0*HCKT(J))
      G = 2.d0
      D1 = 109735.08d0 / 5.5d0 / 5.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (63)  ! ---- Silicon I ----
      PART(1) = B_dep * (1.d0 + 3.d0*EXP(-77.115d0*HCKT(J)) &
        + 5.d0*EXP(-223.157d0*HCKT(J)))
      DO I = 2, 11
        PART(1) = PART(1) + GSI1(I) * B_dep * EXP(-ESI1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 76.d0*EXP(-53000.d0*HCKT(J)) &
        + 71.d0*EXP(-57000.d0*HCKT(J)) + 191.d0*EXP(-60000.d0*HCKT(J)) &
        + 240.d0*EXP(-62000.d0*HCKT(J)) + 251.d0*EXP(-63000.d0*HCKT(J)) &
        + 300.d0*EXP(-65000.d0*HCKT(J))
      CYCLE  ! Si I has no high-level correction (goes direct to 18)

    CASE (64)  ! ---- Silicon II ----
      PART(2) = B_dep * (2.d0 + 4.d0*EXP(-287.32d0*HCKT(J)))
      DO I = 2, 6
        PART(2) = PART(2) + GSI2(I) * B_dep * EXP(-ESI2(I)*HCKT(J))
      END DO
      PART(2) = PART(2) + 6.d0*EXP(-81231.59d0*HCKT(J)) &
        + 6.d0*EXP(-83937.08d0*HCKT(J)) + 10.d0*EXP(-101024.09d0*HCKT(J)) &
        + 14.d0*EXP(-103556.35d0*HCKT(J)) + 10.d0*EXP(-108800.d0*HCKT(J)) &
        + 42.d0*EXP(-115000.d0*HCKT(J)) + 6.d0*EXP(-121000.d0*HCKT(J)) &
        + 38.d0*EXP(-125000.d0*HCKT(J)) + 34.d0*EXP(-132000.d0*HCKT(J))
      G = 2.d0
      D1 = 4.d0 * 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (91)  ! ---- Potassium I ----
      PART(1) = B_dep * 2.d0
      DO I = 2, 8
        PART(1) = PART(1) + GK1(I) * B_dep * EXP(-EK1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 10.d0*EXP(-27397.077d0*HCKT(J)) &
        + 14.d0*EXP(-28127.85d0*HCKT(J))
      G = 2.d0
      D1 = 109734.83d0 / 5.5d0 / 5.5d0 * HCKT(J)
      CALL pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      CYCLE

    CASE (354)  ! ---- Carbon I (extended, NNN67) ----
      PART(1) = B_dep * (1.d0 + 3.d0*EXP(-16.42d0*HCKT(J)) &
        + 5.d0*EXP(-43.42d0*HCKT(J)))
      DO I = 2, 14
        PART(1) = PART(1) + GC1(I) * B_dep * EXP(-EC1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 108.d0*EXP(-80000.d0*HCKT(J)) &
        + 189.d0*EXP(-84000.d0*HCKT(J)) + 247.d0*EXP(-87000.d0*HCKT(J)) &
        + 231.d0*EXP(-88000.d0*HCKT(J)) + 190.d0*EXP(-89000.d0*HCKT(J)) &
        + 300.d0*EXP(-90000.d0*HCKT(J))
      CYCLE  ! C I has no high-level correction

    CASE (355)  ! ---- Carbon II (extended) ----
      PART(2) = B_dep * (2.d0 + 4.d0*EXP(-63.42d0*HCKT(J)))
      DO I = 2, 6
        PART(2) = PART(2) + GC2(I) * B_dep * EXP(-EC2(I)*HCKT(J))
      END DO
      PART(2) = PART(2) + 6.d0*EXP(-131731.80d0*HCKT(J)) &
        + 4.d0*EXP(-142027.1d0*HCKT(J)) + 10.d0*EXP(-145550.13d0*HCKT(J)) &
        + 10.d0*EXP(-150463.62d0*HCKT(J)) + 2.d0*EXP(-157234.07d0*HCKT(J)) &
        + 6.d0*EXP(-162500.d0*HCKT(J)) + 42.d0*EXP(-168000.d0*HCKT(J)) &
        + 56.d0*EXP(-178000.d0*HCKT(J)) + 102.d0*EXP(-183000.d0*HCKT(J)) &
        + 400.d0*EXP(-188000.d0*HCKT(J))
      CYCLE  ! C II has no high-level correction

    CASE (367)  ! ---- Oxygen I (extended) ----
      PART(1) = B_dep * (5.d0 + 3.d0*EXP(-158.265d0*HCKT(J)) &
        + EXP(-226.977d0*HCKT(J)))
      DO I = 2, 13
        PART(1) = PART(1) + GO1(I) * B_dep * EXP(-EO1(I)*HCKT(J))
      END DO
      PART(1) = PART(1) + 15.d0*EXP(-101140.d0*HCKT(J)) &
        + 131.d0*EXP(-103000.d0*HCKT(J)) + 128.d0*EXP(-105000.d0*HCKT(J)) &
        + 600.d0*EXP(-107000.d0*HCKT(J))
      CYCLE  ! O I has no high-level correction

    CASE DEFAULT
      !-----------------------------------------------------------------
      ! Generic NNN interpolation for all other elements/ions
      !-----------------------------------------------------------------
      T2000 = IP(ION) * 2000.d0 / 11.d0
      IT = MAX(1, MIN(9, INT(T(J)/T2000 - 0.5d0)))
      DT = T(J)/T2000 - dble(IT) - 0.5d0
      PMIN = 1.d0
      I = (IT + 1) / 2

      K1 = NNN(I,N) / 100000
      K2 = NNN(I,N) - K1*100000
      K3 = K2 / 10
      KSCALE = K2 - K3*10

      IF (MOD(IT,2) .EQ. 0) THEN
        ! Even IT: interpolate between K3 and next column K1
        P1 = K3 * SCALE(KSCALE)
        K1 = NNN(I+1,N) / 100000
        KSCALE = MOD(NNN(I+1,N), 10)
        P2 = K1 * SCALE(KSCALE)
      ELSE
        ! Odd IT: interpolate between K1 and K3
        P1 = K1 * SCALE(KSCALE)
        P2 = K3 * SCALE(KSCALE)
        IF (DT .LT. 0.d0 .AND. KSCALE .LE. 1) THEN
          KP1 = int(P1)
          IF (KP1 .EQ. INT(P2 + 0.5d0)) PMIN = dble(KP1)
        ENDIF
      ENDIF

      PART(ION) = MAX(PMIN, P1 + (P2 - P1)*DT)

      ! Patch (14 Jun 2004, J. Laird): ensure partition function is >=
      ! the ground-state partition function at low T.  PFGROUND_HYBRID
      ! returns the B&C-based value where tabulated (T <= 10000 K), and
      ! otherwise falls back to Kurucz's pristine hand-coded expressions.
      IF (T(J) .LT. T2000*2.0d0) THEN
        PART(ION) = MAX(PFGROUND_HYBRID((IZ-1)*6 + ION, T(J)), PART(ION))
        CYCLE  ! skip high-level correction at low T
      ENDIF

      ! High-level (Rydberg series) correction
      IF (G .EQ. 0.d0 .OR. POTLO(ION) .LT. 0.1d0 .OR. T(J) .LT. T2000*4.d0) CYCLE
      IF (T(J) .GT. T2000*11.d0) TV = (T2000*11.d0) * KBOL_EV
      D1 = 0.1d0 / TV
      D2 = POTLO(ION) / TV
      PART(ION) = PART(ION) + G*EXP(-IP(ION)/TV) * &
        ( SQRT(13.595d0*Z_ion*Z_ion/TV/D2)**3 * &
          (1.d0/3.d0 + (1.d0 - (0.5d0 + (1.d0/18.d0 + D2/120.d0)*D2)*D2)*D2) &
        - SQRT(13.595d0*Z_ion*Z_ion/TV/D1)**3 * &
          (1.d0/3.d0 + (1.d0 - (0.5d0 + (1.d0/18.d0 + D1/120.d0)*D1)*D1)*D1) )
      TV = TKEV(J)
      CYCLE

    END SELECT
  END DO  ! ION

  !=====================================================================
  ! Phase 2: Saha equation and output
  !=====================================================================
  IF (MODE1 .EQ. 5) THEN
    ! MODE 5: partition functions + cumulative ionization potentials
    ANSWER(32,1) = 0.d0
    DO ION = 1, NION
      ANSWER(ION,1)    = PART(ION)
      ANSWER(ION+32,1) = IP(ION) + ANSWER(ION+31,1)
    END DO
    RETURN
  END IF

  IF (MODE1 .NE. 3) THEN
    ! Compute Saha ionization fractions
    N = N - NION2
    CF = SAHA_PREFAC * T(J) * SQRT(T(J)) / XNE(J)

    DO ION = 2, NION2
      N = N + 1
      F(ION) = CF * PART(ION) / PART(ION-1) * EXP(-(IP(ION-1) - POTLO(ION-1))/TV)
    END DO

    ! Normalize fractions
    F(1) = 1.d0
    L = NION2 + 1
    DO ION = 2, NION2
      L = L - 1
      F(1) = 1.d0 + F(L)*F(1)
    END DO
    F(1) = 1.d0 / F(1)

    DO ION = 2, NION2
      F(ION) = F(ION-1) * F(ION)
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Fill ANSWER based on MODE
  !---------------------------------------------------------------------
  IF (MODE .GE. 10) THEN
    ! Multi-ion output modes (MODE = 11, 12, 13, 14)
    SELECT CASE (MODE1)
    CASE (1)
      DO ION = 1, NION
        ANSWER(J,ION) = F(ION) / PART(ION)
      END DO
    CASE (2)
      DO ION = 1, NION
        ANSWER(J,ION) = F(ION)
      END DO
    CASE (3)
      DO ION = 1, NION
        ANSWER(J,ION) = PART(ION)
      END DO
    CASE (4)
      ANSWER(J,1) = 0.d0
      DO ION = 2, NION2
        ANSWER(J,1) = ANSWER(J,1) + F(ION) * dble(ION-1)
      END DO
    END SELECT
  ELSE
    ! Single-ion output modes (MODE = 1, 2, 3, 4)
    SELECT CASE (MODE1)
    CASE (1);  ANSWER(J,1) = F(NION) / PART(NION)
    CASE (2);  ANSWER(J,1) = F(NION)
    CASE (3);  ANSWER(J,1) = PART(NION)
    CASE (4)
      ANSWER(J,1) = 0.d0
      DO ION = 2, NION2
        ANSWER(J,1) = ANSWER(J,1) + F(ION) * dble(ION-1)
      END DO
    END SELECT
  ENDIF
  RETURN

END SUBROUTINE PFSAHA


!-----------------------------------------------------------------------
! Helper: add Rydberg high-level contribution to partition function
!-----------------------------------------------------------------------
SUBROUTINE pfsaha_highlevels(PART_ION, G, IP_ION, Z_ion, TV_in, &
                             POTLO_ION, D1)
!-----------------------------------------------------------------------
! High-level (Rydberg series) correction to partition function.
! Called ONLY by special elements, which jump directly to label 14 in
! the original code, bypassing the T2000 guard, TV cap, and D1=0.1/TV.
! Those guards apply only to the generic NNN interpolation path.
!-----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(8), INTENT(INOUT) :: PART_ION
  REAL(8), INTENT(IN)    :: G, IP_ION, Z_ion, TV_in, POTLO_ION, D1
  REAL(8) :: D2

  D2 = POTLO_ION / TV_in

  PART_ION = PART_ION + G*EXP(-IP_ION/TV_in) * &
    ( SQRT(13.595d0*Z_ion*Z_ion/TV_in/D2)**3 * &
      (1.d0/3.d0 + (1.d0 - (0.5d0 + (1.d0/18.d0 + D2/120.d0)*D2)*D2)*D2) &
    - SQRT(13.595d0*Z_ion*Z_ion/TV_in/D1)**3 * &
      (1.d0/3.d0 + (1.d0 - (0.5d0 + (1.d0/18.d0 + D1/120.d0)*D1)*D1)*D1) )

END SUBROUTINE pfsaha_highlevels


!=======================================================================
! SUBROUTINE MOLEC
!
! Molecular and atomic number density lookup.
!
! Three paths:
!   CODOUT >= 100 : Molecular species — look up by molecule code in
!                   XNMOLCODE, return number density or fractional pop.
!   CODOUT < 100  : Atomic element — first try molecule table, then
!                   fall through to PFSAHA for ionization equilibrium.
!
! MODE controls output:
!   MODE 1 or 11: fractional populations (XNFPMOL)
!   MODE 2 or 12: number densities (XNMOL)
!   MODE 11/12: multi-ion output
!=======================================================================

SUBROUTINE MOLEC(CODOUT, MODE, NUMBER)

  IMPLICIT NONE

  ! Arguments
  REAL(8),  INTENT(IN)    :: CODOUT
  INTEGER, INTENT(IN)    :: MODE
  REAL(8),  INTENT(INOUT) :: NUMBER(kw, 1)

  ! Local variables
  INTEGER :: JMOL, J, NN, ION, ID, I, II
  REAL(8)  :: C
  LOGICAL :: found

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING MOLEC'

  !=====================================================================
  ! Read molecular input data (first call only)
  !=====================================================================
  IF (IFPOP .NE. 2 .AND. MOLEC_IREAD .NE. 1 .AND. IFPRES .NE. 1) THEN
    READ(INPUTDATA, '(I5)') NUMMOL

    DO JMOL = 1, NUMMOL
      READ(INPUTDATA, '(F20.2)') XNMOLCODE(JMOL)
      READ(INPUTDATA, '(1P8E10.3)') (XNMOL(J,JMOL), J=1,NRHOX)
      IF (IDEBUG .EQ. 1) WRITE(6, '(F20.2/(1P8E10.3))') XNMOLCODE(JMOL), (XNMOL(J,JMOL), J=1,NRHOX)
    END DO

    READ(INPUTDATA, '(1P8E10.3)')  ! skip header line
    READ(INPUTDATA, '(1P8E10.3)') (XNATOM(J), RHO(J), J=1,NRHOX)
    IF (IDEBUG .EQ. 1) WRITE(6, '(1P8E10.3)') (XNATOM(J), RHO(J), J=1,NRHOX)

    READ(INPUTDATA, '(1P8E10.3)')  ! skip header line
    READ(INPUTDATA, '(1P8E10.3)') (XNE(J), J=1,NRHOX)
    IF (IDEBUG .EQ. 1) WRITE(6, '(1P8E10.3)') (XNE(J), J=1,NRHOX)

    MOLEC_IREAD = 1
  ENDIF

  !=====================================================================
  ! Path 1: Molecular species lookup (CODOUT >= 100)
  !=====================================================================
  IF (CODOUT .GE. 100.d0) THEN
    DO JMOL = 1, NUMMOL
      IF (XNMOLCODE(JMOL) .EQ. CODOUT) THEN
        DO J = 1, NRHOX
          IF (MODE .EQ. 1 .OR. MODE .EQ. 11) NUMBER(J,1) = XNFPMOL(J,JMOL)
          IF (MODE .EQ. 2 .OR. MODE .EQ. 12) NUMBER(J,1) = XNMOL(J,JMOL)
        END DO
        RETURN
      ENDIF
    END DO
    ! Species not in molecular equilibrium table: return zero populations
    NUMBER(:, 1) = 0.d0
    RETURN
  ENDIF

  !=====================================================================
  ! Path 2: Atomic element lookup (CODOUT < 100)
  !=====================================================================
  C = CODOUT
  NN = 1
  IF (MODE .EQ. 11 .OR. MODE .EQ. 12) NN = INT((C - AINT(C))*100.d0 + 1.5d0)

  DO I = 1, NN
    ION = NN - I + 1

    ! Search molecule table for approximate code match
    found = .FALSE.
    DO JMOL = 1, NUMMOL
      IF (XNMOLCODE(JMOL) + 0.001d0 .GT. C .AND. &
          XNMOLCODE(JMOL) - 0.001d0 .LT. C) THEN
        DO J = 1, NRHOX
          IF (MODE .EQ. 1 .OR. MODE .EQ. 11) NUMBER(J,ION) = XNFPMOL(J,JMOL)
          IF (MODE .EQ. 2 .OR. MODE .EQ. 12) NUMBER(J,ION) = XNMOL(J,JMOL)
        END DO
        found = .TRUE.
        EXIT  ! exit JMOL loop
      ENDIF
    END DO

    IF (.NOT. found) THEN
      ! Try matching by integer element code
      ID = INT(CODOUT)
      found = .FALSE.
      DO JMOL = 1, NUMMOL
        IF (INT(XNMOLCODE(JMOL)) .EQ. ID) THEN
          ! Element exists but this specific ion not tabulated: zero it
          NUMBER(:, ION) = 0.d0
          found = .TRUE.
          EXIT  ! exit JMOL loop
        ENDIF
      END DO

      IF (.NOT. found) THEN
        ! Element not in molecule table at all: use Saha equation
        ION = INT((CODOUT - dble(ID))*100.d0 + 1.5d0)
        NN = ION
        IF (MODE .EQ. 1) NN = 1
        DO J = 1, NRHOX
          CALL PFSAHA(J, ID, ION, MODE, NUMBER)
          DO II = 1, NN
            NUMBER(J,II) = NUMBER(J,II) * XNATOM(J) * XABUND(J,ID)
          END DO
        END DO
        RETURN  ! exit entire subroutine (matches original GO TO 400)
      ENDIF
    ENDIF

    C = C - 0.01d0
  END DO  ! I

  RETURN

END SUBROUTINE MOLEC

!=========================================================================
! SUBROUTINE NMOLEC(MODE)
!
! Molecular chemical equilibrium solver.
!
! Solves the coupled system of chemical equilibrium equations for all
! molecules, atoms, and ions simultaneously at each depth point. This
! is the molecular-equilibrium path (IFMOL=1); when molecules are off,
! NELECT handles the electron density instead.
!
! The system of equations enforces:
!   1. Total particle conservation: sum of all particles = P / kT
!   2. Element conservation: for each element, the total abundance
!      (free atoms + atoms bound in molecules) equals n_atom * A(Z)
!   3. Charge conservation: n_e = sum of ionic charges (if ions present)
!
! Molecules contribute additional terms: each molecule with equilibrium
! constant K_mol consumes its constituent atoms according to:
!   n_mol = K_mol * product(n_component_k)
!
! The nonlinear system is solved by Newton-Raphson iteration:
!   - Build the Jacobian matrix DEQ and residual vector EQ
!   - Solve DEQ * delta = EQ via SOLVIT (Gaussian elimination)
!   - Apply corrections with damping (factor-of-100 limiter plus
!     sign-flip detection for oscillation control)
!   - Iterate until all residuals < 0.01% of the current values
!
! After convergence, the routine also computes:
!   - n/U (number density / partition function) for opacity calculations
!     (when MODE = 1 or 11)
!   - Molecular contribution to the energy density (when IFEDNS = 1,
!     called from CONVEC)
!
! Arguments:
!   MODE = 1 or 11:  compute n/U for opacity (line + continuum)
!   MODE = 2 or 12:  skip n/U computation (populations only)
!=========================================================================

SUBROUTINE NMOLEC(MODE)

  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN) :: MODE

  ! --- Local arrays ---
  REAL(8)  :: EQUILJ(maxmol)          ! equilibrium constants at current depth
  REAL(8)  :: XNZ(kw, maxeq)          ! converged number densities per equation
  REAL(8)  :: EQ(maxeq)               ! residual vector
  REAL(8)  :: XN(maxeq)               ! current number density estimates
  REAL(8)  :: XAB(maxeq)              ! element abundances
  REAL(8)  :: DEQ(maxeq * maxeq)      ! Jacobian matrix (flattened)
  REAL(8)  :: EQOLD(maxeq)            ! previous residuals (for oscillation detection)
  INTEGER :: IPIVOT(maxeq)           ! pivot array for SOLVIT

  REAL(8)  :: FRAC(kw, 6)             ! scratch for PFSAHA output
  REAL(8)  :: PFP(61), PFM(61)        ! partition functions at T+/T- (atoms)
  REAL(8)  :: PFPLUS(kw), PFMIN(kw)   ! partition functions at T+/T- (molecules)

  ! For electron contribution diagnostic
  REAL(8)  :: E_CONTRIB(kw, maxmol)
  REAL(8)  :: XE_CONTRIB(kw, maxmol)
  INTEGER :: IDZ(maxmol), NION_MOL(maxmol)

  ! --- Local scalars ---
  REAL(8)  :: XNTOT          ! total particle density = P / kT
  REAL(8)  :: RATIO          ! pressure ratio for depth extrapolation
  REAL(8)  :: TERM, D        ! molecule contribution to Jacobian
  REAL(8)  :: X, XNEQ, XN100, SCALE
  REAL(8)  :: AMASS          ! molecular mass accumulator
  INTEGER :: J, K, M, KK, K1, MK
  INTEGER :: JMOL, NCOMP, ION, ID
  INTEGER :: LOCK, LOCJ1, LOCJ2, LOCM, NEQUAK
  INTEGER :: IFERR
  INTEGER :: JMOL1, JMOL10, NN, NZ = 0, IZ, IZ1, IZ10

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING NMOLEC'
  
  NEQUA1 = NEQUA + 1
  NEQNEQ = NEQUA**2

  !=====================================================================
  ! Initialize abundances (constant with depth in this version)
  !=====================================================================
  XAB(:) = 0.0D0
  J = 1
  DO K = 2, NEQUA
    ID = IDEQUA(K)
    IF (ID .LT. 100) XAB(K) = max(XABUND(J, ID), 1.0D-20)
  END DO
  IF (ID .EQ. 100) XAB(NEQUA) = 0.0D0

  ! Initial guess for number densities at the surface
  XNTOT = P(1) / TK(1)
  XN(1) = XNTOT / 2.0D0
  IF (T(1) .LT. 4000.0D0) XN(1) = XNTOT
  X = XN(1) / 10.0D0
  DO K = 2, NEQUA
    XN(K) = X * XAB(K)
  END DO
  IF (ID .EQ. 100) XN(NEQUA) = X
  XNE(1) = X
  
  !=====================================================================
  ! Main depth loop: solve molecular equilibrium at each depth
  !=====================================================================
  DO J = 1, NRHOX
     
    XNTOT = P(J) / TK(J)

    ! Extrapolate from previous depth as initial guess
    IF (J .GT. 1) THEN
      RATIO = P(J) / P(J - 1)
      XNE(J) = XNE(J - 1) * RATIO
      DO K = 1, NEQUA
        XN(K) = XN(K) * RATIO
      END DO
    END IF

    ! If computing energy density, restore saved values from previous
    ! normal (non-EDENS) run as the starting point
    IF (IFEDNS .NE. 0) THEN
      DO K = 1, NEQUA
        XN(K) = XNSAVE(J, K)
      END DO
    END IF

    !-------------------------------------------------------------------
    ! Compute equilibrium constants for all molecules at this depth
    !-------------------------------------------------------------------
    DO JMOL = 1, NUMMOL
      NCOMP = LOCJ(JMOL + 1) - LOCJ(JMOL)

      IF (EQUIL(1, JMOL) .NE. 0.0D0) THEN
        ! True molecule: has equilibrium constant polynomial
        ION = int((XNMOLCODE(JMOL) - aint(XNMOLCODE(JMOL))) * 100.0D0 + 0.5D0)
        EQUILJ(JMOL) = 0.0D0

        IF (XNMOLCODE(JMOL) .EQ. 101.0D0) THEN
          ! H2: special equilibrium function (Lester 2005 update)
          IF (T(J) .LE. 20000.0D0) EQUILJ(JMOL) = EQUILH2(T(J))
        ELSE
          ! General molecule: polynomial equilibrium constant
          IF (T(J) .LE. 10000.0D0) THEN
            EQUILJ(JMOL) = exp(EQUIL(1, JMOL) / TKEV(J) - EQUIL(2, JMOL) &
              + (EQUIL(3, JMOL) + (-EQUIL(4, JMOL) + (EQUIL(5, JMOL) &
              - EQUIL(6, JMOL) * T(J)) * T(J)) * T(J)) * T(J) &
              - 1.5D0 * (NCOMP - ION - ION - 1) * TLOG(J))
          END IF
        END IF

      ELSE IF (NCOMP .EQ. 1) THEN
        ! Single atom: trivial equilibrium
        EQUILJ(JMOL) = 1.0D0

      ELSE
        ! Multi-stage atom: get ionization equilibrium from Saha
        ID = int(XNMOLCODE(JMOL))
        ION = NCOMP - 1
        CALL PFSAHA(J, ID, NCOMP, 12, FRAC)
        EQUILJ(JMOL) = FRAC(J, NCOMP) / FRAC(J, 1) * XNE(J)**ION
      END IF
    END DO

    !-------------------------------------------------------------------
    ! Newton-Raphson iteration for chemical equilibrium
    !-------------------------------------------------------------------
    EQOLD(:) = 0.0D0

    newton_loop: DO

      ! --- Build Jacobian (DEQ) and residual (EQ) ---
      DEQ(1:NEQNEQ) = 0.0D0

      ! Equation 1: total particle conservation
      EQ(1) = -XNTOT
      K1 = 1
      KK = 1
      DO K = 2, NEQUA
        EQ(1) = EQ(1) + XN(K)
        K1 = K1 + NEQUA
        DEQ(K1) = 1.0D0              ! d(eq1)/d(XN_k) = 1
        EQ(K) = XN(K) - XAB(K) * XN(1)  ! element conservation
        KK = KK + NEQUA1
        DEQ(KK) = 1.0D0              ! d(eq_k)/d(XN_k) = 1
        DEQ(K) = -XAB(K)             ! d(eq_k)/d(XN_1) = -A_k
      END DO

      ! Charge conservation (if ions present)
      IF (IDEQUA(NEQUA) .GE. 100) THEN
        EQ(NEQUA) = -XN(NEQUA)
        DEQ(NEQNEQ) = -1.0D0
      END IF

      ! --- Add molecular contributions ---
      DO JMOL = 1, NUMMOL
        NCOMP = LOCJ(JMOL + 1) - LOCJ(JMOL)
        IF (NCOMP .EQ. 1) CYCLE         ! single atoms don't contribute

        ! Compute molecule number density: TERM = K * product(n_k)
        TERM = EQUILJ(JMOL)
        LOCJ1 = LOCJ(JMOL)
        LOCJ2 = LOCJ(JMOL + 1) - 1
        DO LOCK = LOCJ1, LOCJ2
          K = KCOMPS(LOCK)
          IF (K .EQ. NEQUA1) THEN
            TERM = TERM / XN(NEQUA)   ! electron: divide by n_e
          ELSE
            TERM = TERM * XN(K)        ! atom: multiply by n_k
          END IF
        END DO

        ! Add to total particle count
        EQ(1) = EQ(1) + TERM

        ! Add to element conservation and Jacobian
        DO LOCK = LOCJ1, LOCJ2
          K = KCOMPS(LOCK)
          IF (K .GE. NEQUA1) THEN
            K = NEQUA
            D = -TERM / XN(K)          ! electron: d(TERM)/d(n_e) = -TERM/n_e
          ELSE
            D = TERM / XN(K)           ! atom: d(TERM)/d(n_k) = TERM/n_k
          END IF
          EQ(K) = EQ(K) + TERM
          NEQUAK = NEQUA * K - NEQUA
          K1 = NEQUAK + 1
          DEQ(K1) = DEQ(K1) + D
          DO LOCM = LOCJ1, LOCJ2
            M = KCOMPS(LOCM)
            IF (M .EQ. NEQUA1) M = NEQUA
            MK = M + NEQUAK
            DEQ(MK) = DEQ(MK) + D
          END DO
        END DO

        ! Correction to charge equation for negative ions
        ! (species whose last component is element 100 = neutral atom)
        K = KCOMPS(LOCJ2)
        IF (IDEQUA(K) .EQ. 100) THEN
          DO LOCK = LOCJ1, LOCJ2
            K = KCOMPS(LOCK)
            D = TERM / XN(K)
            IF (K .EQ. NEQUA) EQ(K) = EQ(K) - TERM - TERM
            NEQUAK = NEQUA * K - NEQUA
            DO LOCM = LOCJ1, LOCJ2
              M = KCOMPS(LOCM)
              IF (M .EQ. NEQUA) THEN
                MK = M + NEQUAK
                DEQ(MK) = DEQ(MK) - D - D
              END IF
            END DO
          END DO
        END IF

      END DO  ! JMOL

      ! --- Solve linearized system ---
      CALL SOLVIT(DEQ, NEQUA, EQ, IPIVOT)

      ! --- Apply corrections with damping ---
      IFERR = 0
      SCALE = 100.0D0
      DO K = 1, NEQUA
        RATIO = abs(EQ(K) / XN(K))
        IF (RATIO .GT. 1.0D-4) IFERR = 1

        ! Reduce correction if oscillating (sign flip)
        IF (EQOLD(K) * EQ(K) .LT. 0.0D0) EQ(K) = EQ(K) * 0.69D0

        XNEQ = XN(K) - EQ(K)
        XN100 = XN(K) / 100.0D0

        IF (abs(XNEQ) .GE. XN100) THEN
          ! Normal correction: accept new value (force positive)
          XN(K) = abs(XNEQ)
        ELSE
          ! Correction too large relative to current value: scale down
          XN(K) = XN(K) / SCALE
          IF (EQOLD(K) * EQ(K) .LT. 0.0D0) SCALE = sqrt(SCALE)
        END IF

        EQOLD(K) = EQ(K)
      END DO

      IF (IFERR .EQ. 0) EXIT newton_loop

    END DO newton_loop

    !-------------------------------------------------------------------
    ! Store converged solution
    !-------------------------------------------------------------------
    DO K = 1, NEQUA
      XNZ(J, K) = XN(K)
    END DO
    XNATOM(J) = XN(1)
    RHO(J) = XNATOM(J) * WTMOLE(J) * AMU
    IF (IDEQUA(NEQUA) .EQ. 100) XNE(J) = XN(NEQUA)

    ! Compute molecular number densities from converged XN
    DO JMOL = 1, NUMMOL
      XNMOL(J, JMOL) = EQUILJ(JMOL)
      LOCJ1 = LOCJ(JMOL)
      LOCJ2 = LOCJ(JMOL + 1) - 1
      DO LOCK = LOCJ1, LOCJ2
        K = KCOMPS(LOCK)
        IF (K .EQ. NEQUA1) THEN
          XNMOL(J, JMOL) = XNMOL(J, JMOL) / XN(NEQUA)
        ELSE
          XNMOL(J, JMOL) = XNMOL(J, JMOL) * XN(K)
        END IF
      END DO
    END DO

  END DO  ! depth loop J

  !=====================================================================
  ! If computing energy density (IFEDNS=1), jump to that section
  !=====================================================================
  IF (IFEDNS .EQ. 1) THEN
    CALL nmolec_energy_density()
    RETURN
  END IF

  
  !=====================================================================
  ! Save solution for next call / EDENS restart
  !=====================================================================

  DO K = 1, NEQUA
    XNSAVE(:, K) = XNZ(:, K)
  END DO

  ! --- Diagnostic output ---
  IF ((ITER .EQ. 1 .OR. ITER .EQ. NUMITS) .AND. IDEBUG .EQ. 1) THEN
    WRITE(6, '(10X,"RHOX",9X,"T",11X,"P",10X,"XNE",8X,"XNATOM",8X,"RHO")') 
    DO J = 1, NRHOX
      WRITE(6, '(I5,1P6E12.3)') J, RHOX(J), T(J), P(J), XNE(J), XNATOM(J), RHO(J)
    END DO
  END IF

  IF (ITER .EQ. NUMITS .AND. IFSYNTHE .EQ. 1 .AND. IFMOLOUT .EQ. 1) THEN
    DO JMOL1 = 1, NUMMOL, 10
      JMOL10 = min(JMOL1 + 9, NUMMOL)
      WRITE(35, '(49X,"MOLECULAR NUMBER DENSITIES"/5X,10F12.2)') &
        (XNMOLCODE(JMOL), JMOL = JMOL1, JMOL10)
      DO J = 1, NRHOX
        WRITE(35, '(I5,1P10E12.3)') J, (XNMOL(J, JMOL), JMOL = JMOL1, JMOL10)
      END DO
    END DO
  END IF

  !=====================================================================
  ! MODE 1 or 11: compute n / partition_function for opacity
  !=====================================================================
  IF (MODE .NE. 2 .AND. MODE .NE. 12) THEN

  ! Convert element number densities to n/U
  DO K = 2, NEQUA
    ID = IDEQUA(K)
    IF (ID .EQ. 100) THEN
      ! Electrons: divide by electron partition function (= 2)
      ! and translational PF: (2*pi*m_e*kT/h^2)^(3/2) = 2.4148e15 * T^(3/2)
      XNZ(:, K) = XNZ(:, K) / SAHA_PREFAC / T / sqrt(T)
    ELSE
      ! Atoms: divide by partition function and translational PF
      DO J = 1, NRHOX
        CALL PFSAHA(J, ID, 1, 3, FRAC)
        XNZ(J, K) = XNZ(J, K) / FRAC(J, 1) / 1.8786D20 &
                   / sqrt((ATMASS(ID) * T(J))**3)
      END DO
    END IF
  END DO

  ! Compute n/U for molecules and track element IDs for electron diagnostics
  IZ = 1
  IDZ(IZ) = 1
  DO JMOL = 1, NUMMOL
    NCOMP = LOCJ(JMOL + 1) - LOCJ(JMOL)

    IF (EQUIL(1, JMOL) .NE. 0.0D0) THEN
      ! True molecule: n_mol/U = exp(D0/kT) * product(n_k/U_k) * translational_PF
      AMASS = 0.0D0
      LOCJ1 = LOCJ(JMOL)
      LOCJ2 = LOCJ(JMOL + 1) - 1

      XNFPMOL(:, JMOL) = exp(EQUIL(1, JMOL) / TKEV)

      DO LOCK = LOCJ1, LOCJ2
        K = KCOMPS(LOCK)
        IF (K .EQ. NEQUA1) THEN
          XNFPMOL(:, JMOL) = XNFPMOL(:, JMOL) / XNZ(:, NEQUA)
        ELSE
          ID = IDEQUA(K)
          IF (ID .LT. 100) AMASS = AMASS + ATMASS(ID)
          XNFPMOL(:, JMOL) = XNFPMOL(:, JMOL) * XNZ(:, K)
        END IF
      END DO

      XNFPMOL(:, JMOL) = XNFPMOL(:, JMOL) * 1.8786D20 * sqrt((AMASS * T)**3)

    ELSE
      ! Atom/ion: use PFSAHA to get n/U
      ID = int(XNMOLCODE(JMOL))
      IF (ID .NE. IDZ(IZ)) THEN
        IF (JMOL .GE. 3) NION_MOL(IZ) = LOCJ(JMOL - 1) - LOCJ(JMOL - 2)
        IZ = IZ + 1
        IDZ(IZ) = ID
      END IF

      DO J = 1, NRHOX
        CALL PFSAHA(J, ID, NCOMP, 3, FRAC)
        XNFPMOL(J, JMOL) = XNMOL(J, JMOL) / FRAC(J, 1)
      END DO
    END IF
  END DO

  ! Compute electron contribution per element
  NZ = IZ
  DO IZ = 1, NZ
    DO J = 1, NRHOX
      CALL PFSAHA(J, IDZ(IZ), NION_MOL(IZ), 4, E_CONTRIB)
      XE_CONTRIB(J, IZ) = E_CONTRIB(J, 1) * XNATOM(J) &
                         * XABUND(J, IDZ(IZ)) / XNE(J)
    END DO
  END DO

  END IF  ! MODE /= 2 and MODE /= 12

  ! --- Print n/U and electron contribution diagnostics ---
  IF (ITER .EQ. NUMITS.AND. IDEBUG .EQ. 1) THEN
    NN = ((NUMMOL / 10) + 1) * 10
    DO JMOL1 = 1, NN, 10
      JMOL10 = JMOL1 + 9
      WRITE(6, '("1",40X,"NUMBER DENSITIES / PARTITION FUNCTIONS"/5X,10F12.2/(I5,1P10E12.3))') &
        (XNMOLCODE(JMOL), JMOL = JMOL1, JMOL10), &
        (J, (XNFPMOL(J, JMOL), JMOL = JMOL1, JMOL10), J = 1, NRHOX)
    END DO

    DO IZ1 = 1, NZ, 10
      IZ10 = IZ1 + 9
      WRITE(6, '("1",40X,"ELECTRON CONTRIBUTION Ne(elem)/Ne_tot"/5X,10I12/(I5,1P10E12.3))') &
        (IDZ(IZ), IZ = IZ1, IZ10), &
        (J, (XE_CONTRIB(J, IZ), IZ = IZ1, IZ10), J = 1, NRHOX)
    END DO
  END IF

  RETURN

CONTAINS

  !-------------------------------------------------------------------
  ! Compute molecular contribution to energy density (IFEDNS=1 path).
  ! Called from CONVEC via COMPUTE_ONE_POP when energy density needed.
  !-------------------------------------------------------------------
  SUBROUTINE nmolec_energy_density()
    INTEGER :: J_e, JMOL_e, NCOMP_e, ION_e, ID_e
    INTEGER :: LOCK_e, LOCJ1_e, LOCJ2_e, K_e
    REAL(8)  :: TPLUS_e, TMINUS_e

    ! Initialize with kinetic energy
    DO J_e = 1, NRHOX
      XNTOT = P(J_e) / TK(J_e)
      EDENS(J_e) = 1.5D0 * XNTOT * TK(J_e)
    END DO

    DO JMOL_e = 1, NUMMOL
      NCOMP_e = LOCJ(JMOL_e + 1) - LOCJ(JMOL_e)

      IF (EQUIL(1, JMOL_e) .NE. 0.0D0) THEN
        !---------------------------------------------------------------
        ! Molecules: partition function derivative by finite differences
        !---------------------------------------------------------------
        DO J_e = 1, NRHOX
          IF (XNMOL(J_e, JMOL_e) .LE. 0.0D0) CYCLE
          TPLUS_e  = T(J_e) * 1.001D0
          TMINUS_e = T(J_e) * 0.999D0
          PFMIN(J_e)  = 0.0D0
          PFPLUS(J_e) = 0.0D0

          IF (XNMOLCODE(JMOL_e) .EQ. 101.0D0) THEN
            ! H2: use dedicated partition function
            PFMIN(J_e)  = PARTFNH2(TMINUS_e)
            PFPLUS(J_e) = PARTFNH2(TPLUS_e)
          ELSE
            ! General molecule: PF from equilibrium constant polynomial
            IF (T(J_e) .LE. 10000.0D0) THEN
              PFPLUS(J_e) = exp(-EQUIL(2, JMOL_e) &
                + (EQUIL(3, JMOL_e) + (-EQUIL(4, JMOL_e) &
                + (EQUIL(5, JMOL_e) - EQUIL(6, JMOL_e) * TPLUS_e) &
                * TPLUS_e) * TPLUS_e) * TPLUS_e) + 1.0D-30
              PFMIN(J_e) = exp(-EQUIL(2, JMOL_e) &
                + (EQUIL(3, JMOL_e) + (-EQUIL(4, JMOL_e) &
                + (EQUIL(5, JMOL_e) - EQUIL(6, JMOL_e) * TMINUS_e) &
                * TMINUS_e) * TMINUS_e) * TMINUS_e) + 1.0D-30
            END IF
          END IF
        END DO

        ! For non-H2 molecules: multiply PF by component atom PFs at T+/T-
        IF (XNMOLCODE(JMOL_e) .NE. 101.0D0) THEN
          LOCJ1_e = LOCJ(JMOL_e)
          LOCJ2_e = LOCJ(JMOL_e + 1) - 1
          DO LOCK_e = LOCJ1_e, LOCJ2_e
            K_e = KCOMPS(LOCK_e)
            IF (K_e .EQ. NEQUA) CYCLE    ! electron: skip
            IF (K_e .GT. NEQUA) EXIT      ! past last equation
            ID_e = IDEQUA(K_e)
            DO J_e = 1, NRHOX
              IF (XNMOL(J_e, JMOL_e) .LE. 0.0D0) CYCLE
              ! PF at T+
              T(J_e) = T(J_e) * 1.001D0
              TK(J_e) = TK(J_e) * 1.001D0
              TKEV(J_e) = TKEV(J_e) * 1.001D0
              CALL PFSAHA(J_e, ID_e, 1, 3, FRAC)
              PFPLUS(J_e) = PFPLUS(J_e) * FRAC(J_e, 1)
              ! PF at T-
              T(J_e) = T(J_e) / 1.001D0 * 0.999D0
              TK(J_e) = TK(J_e) / 1.001D0 * 0.999D0
              TKEV(J_e) = TKEV(J_e) / 1.001D0 * 0.999D0
              CALL PFSAHA(J_e, ID_e, 1, 3, FRAC)
              PFMIN(J_e) = PFMIN(J_e) * FRAC(J_e, 1)
              ! Restore T
              T(J_e) = T(J_e) / 0.999D0
              TK(J_e) = TK(J_e) / 0.999D0
              TKEV(J_e) = TKEV(J_e) / 0.999D0
            END DO
          END DO
        END IF

        ! Accumulate energy density contribution
        IF (XNMOLCODE(JMOL_e) .EQ. 101.0D0) THEN
          ! H2: dissociation energy = 36118.11 cm^-1
          DO J_e = 1, NRHOX
            IF (XNMOL(J_e, JMOL_e) .LE. 0.0D0) CYCLE
            EDENS(J_e) = EDENS(J_e) + XNMOL(J_e, JMOL_e) * TK(J_e) &
              * (-36118.11D0 * HCKT(J_e) &
              + (PFPLUS(J_e) - PFMIN(J_e)) / (PFPLUS(J_e) + PFMIN(J_e) + 1.0D-30) &
              * 1000.0D0)
          END DO
        ELSE
          ! General molecule
          DO J_e = 1, NRHOX
            IF (XNMOL(J_e, JMOL_e) .LE. 0.0D0) CYCLE
            EDENS(J_e) = EDENS(J_e) + XNMOL(J_e, JMOL_e) * TK(J_e) &
              * (-EQUIL(1, JMOL_e) / TKEV(J_e) &
              + (PFPLUS(J_e) - PFMIN(J_e)) / (PFPLUS(J_e) + PFMIN(J_e) + 1.0D-30) &
              * 1000.0D0)
          END DO
        END IF

      ELSE
        !---------------------------------------------------------------
        ! Atoms: partition function derivative by finite differences
        !---------------------------------------------------------------
        ID_e = int(XNMOLCODE(JMOL_e))
        DO J_e = 1, NRHOX
          IF (XNMOL(J_e, JMOL_e) .LE. 0.0D0) CYCLE
          ! PF at T+
          T(J_e) = T(J_e) * 1.001D0
          TK(J_e) = TK(J_e) * 1.001D0
          TKEV(J_e) = TKEV(J_e) * 1.001D0
          CALL PFSAHA(J_e, ID_e, NCOMP_e, 5, PFP)
          ! PF at T-
          T(J_e) = T(J_e) / 1.001D0 * 0.999D0
          TK(J_e) = TK(J_e) / 1.001D0 * 0.999D0
          TKEV(J_e) = TKEV(J_e) / 1.001D0 * 0.999D0
          CALL PFSAHA(J_e, ID_e, NCOMP_e, 5, PFM)
          ! Restore T
          T(J_e) = T(J_e) / 0.999D0
          TK(J_e) = TK(J_e) / 0.999D0
          TKEV(J_e) = TKEV(J_e) / 0.999D0

          ION_e = NCOMP_e
          PFP(ION_e) = max(PFP(ION_e), PFM(ION_e))
          EDENS(J_e) = EDENS(J_e) + XNMOL(J_e, JMOL_e) * TK(J_e) &
            * (PFP(31 + ION_e) / TKEV(J_e) &
            + (PFP(ION_e) - PFM(ION_e)) / (PFP(ION_e) + PFM(ION_e)) &
            * 1000.0D0)
        END DO

      END IF
    END DO  ! JMOL_e

    ! Convert to energy per unit mass
    DO J_e = 1, NRHOX
      EDENS(J_e) = EDENS(J_e) / RHO(J_e)
    END DO

  END SUBROUTINE nmolec_energy_density

END SUBROUTINE NMOLEC

!=========================================================================
! SUBROUTINE COMPUTE_ALL_POPS
!
! Compute number densities and n/U for all species needed by the
! opacity routines.
!
! Populates two master arrays used throughout the opacity calculation:
!   XNF(J, NELION)  — number density * (ion fraction / partition function)
!                      (MODE=12): used for line opacity via LINSOP
!   XNFP(J, NELION) — number density / partition function
!                      (MODE=11): used for continuum opacity via KAPP/HOP
!
! Species layout (NELION index):
!   Atoms Z=1-30:  NELION = IZ*(IZ+1)/2 .. IZ*(IZ+1)/2 + NION - 1
!                  (triangular number layout, IZ+1 slots per element)
!   Atoms Z=31-99: NELION = 496 + (IZ-31)*5 .. +4
!                  (5 slots per element)
!   Molecules:     NELION = 841..1006 (when IFMOL=1)
!                  See species table below for mapping.
!
! The number of ionization stages (NION) differs between the XNF and
! XNFP passes because continuum opacity requires more ion stages than
! line opacity for the transition metals.
!
! For cool stars (T < 9000 K), H2 and CO equilibrium populations are
! computed from analytic fits when molecules are off (IFMOL=0), since
! these are the dominant molecular absorbers.
!
! When IFMOL=1, molecular populations are obtained from the NMOLEC
! equilibrium solution via COMPUTE_ONE_POP.
!
! Bug fixes: Be and B stages reduced from 4,5 to 3,3 (Stift 2009);
! Ar MODE=11 changed from .05 to .04.
!=========================================================================

SUBROUTINE COMPUTE_ALL_POPS

  IMPLICIT NONE

  ! --- Species code tables for Z=1-30 ---
  ! CODE = IZ + NION/100:  IZ = atomic number, NION = number of ion stages
  ! MODE=12 (XNF): for line opacity — fewer stages needed
  ! MODE=11 (XNFP): for continuum opacity — more stages for transition metals

  ! Number of ion stages per element for MODE=12 (XNF, line opacity)
  INTEGER, PARAMETER :: NION_XNF(30) = (/ &
    1, 2, 3, 3, 3,                        &  ! Z = 1-5  (H, He, Li, Be, B)
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,     &  ! Z = 6-16
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  &  ! Z = 17-28
    2, 2 /)                                   ! Z = 29-30

  ! Number of ion stages per element for MODE=11 (XNFP, continuum opacity)
  INTEGER, PARAMETER :: NION_XNFP(30) = (/ &
    1, 2, 3, 3, 3,                         &  ! Z = 1-5
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,      &  ! Z = 6-16
    5, 4, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9,   &  ! Z = 17-28
    2, 2 /)                                    ! Z = 29-30

  ! NELION offsets: IZ*(IZ+1)/2 for Z=1-30
  INTEGER, PARAMETER :: NOFF(30) = (/ &
      1,   3,   6,  10,  15,  21,  28,  36,  45,  55, &
     66,  78,  91, 105, 120, 136, 153, 171, 190, 210, &
    231, 253, 276, 300, 325, 351, 378, 406, 435, 465 /)

  ! --- Molecular species for IFMOL=1 path ---
  ! CODE → NELION mapping for molecules computed by NMOLEC
  ! H2=841, CH=846, NH=847, OH=848, MgH=851, SiH=853,
  ! CaH=858, CrH=862, FeH=864, C2=868, CN=869, CO=870,
  ! SiO=889, TiO=895, VO=896, H2O=940
  INTEGER, PARAMETER :: NMOL_SPECIES = 16
  REAL(8),  PARAMETER :: MOL_CODE(16) = (/ &
    101.0D0, 106.0D0, 107.0D0, 108.0D0, 112.0D0, 114.0D0, &
    120.0D0, 124.0D0, 126.0D0, 606.0D0, 607.0D0, 608.0D0, &
    814.0D0, 822.0D0, 823.0D0, 10108.0D0 /)
  INTEGER, PARAMETER :: MOL_NELION(16) = (/ &
    841, 846, 847, 848, 851, 853, &
    858, 862, 864, 868, 869, 870, &
    889, 895, 896, 940 /)

  ! --- Local variables ---
  REAL(8)  :: CODE
  INTEGER :: IZ, J, IMOL

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COMPUTE_ALL_POPS'

  !=====================================================================
  ! Pass 1: XNF (MODE=12) — number densities for line opacity
  !=====================================================================
  DO IZ = 1, 30
    CODE = dble(IZ) + dble(NION_XNF(IZ)) / 100.0D0
    CALL COMPUTE_ONE_POP(CODE, 12, XNF(1, NOFF(IZ)))
  END DO

  !=====================================================================
  ! Pass 2: XNFP (MODE=11) — n/partition function for continuum opacity
  !=====================================================================
  DO IZ = 1, 30
    CODE = dble(IZ) + dble(NION_XNFP(IZ)) / 100.0D0
    CALL COMPUTE_ONE_POP(CODE, 11, XNFP(1, NOFF(IZ)))
  END DO

  !=====================================================================
  ! Elements Z=31-99: 2 stages for both passes
  !=====================================================================
  DO IZ = 31, 99
    CODE = dble(IZ) + 0.02D0
    CALL COMPUTE_ONE_POP(CODE, 11, XNFP(1, 496 + (IZ - 31) * 5))
    CALL COMPUTE_ONE_POP(CODE, 12, XNF(1, 496 + (IZ - 31) * 5))
  END DO

  !=====================================================================
  ! Approximate H2 and CO for cool stars when molecules are off
  ! (analytic equilibrium fits, valid for T < 9000 K)
  !=====================================================================
  DO J = 1, NRHOX
    IF (T(J) .GT. 9000.0D0) CYCLE

    ! H2 equilibrium: n(H2)/U(H2) from H partition functions
    XNFP(J, 841) = XNFP(J, 1)**2 * exp(4.478D0 / TKEV(J) - 4.64584D1 &
      + (1.63660D-3 + (-4.93992D-7 + (1.11822D-10 + (-1.49567D-14 &
      + (1.06206D-18 - 3.08720D-23 * T(J)) * T(J)) * T(J)) * T(J)) &
      * T(J)) * T(J) - 1.5D0 * TLOG(J))

    ! CO equilibrium: n(CO)/U(CO) from C and O partition functions
    XNFP(J, 870) = XNFP(J, 21) * XNFP(J, 36) &
      * exp(11.091D0 / TKEV(J) - 49.0414D0 &
      + 14.0306D-4 * T(J) - 26.6341D-8 * T(J)**2 &
      + 35.382D-12 * T(J)**3 - 26.5424D-16 * T(J)**4 &
      + 8.32385D-20 * T(J)**5 - 1.5D0 * TLOG(J))
  END DO

  !=====================================================================
  ! Molecular populations from NMOLEC (when molecular equilibrium is on)
  !=====================================================================
  IF (IFMOL .EQ. 0) RETURN

  DO IMOL = 1, NMOL_SPECIES
    CALL COMPUTE_ONE_POP(MOL_CODE(IMOL), 1, XNFP(1, MOL_NELION(IMOL)))
  END DO

  RETURN

END SUBROUTINE COMPUTE_ALL_POPS

!=========================================================================
! SUBROUTINE CONVEC
!
! Convective flux and thermodynamic quantities via mixing length theory.
!
! Computes the convective energy flux at each depth point using the
! Böhm-Vitense (1958) formulation of mixing length theory (MLT).
! The thermodynamic derivatives needed for the adiabatic gradient,
! specific heats, and sound speed are obtained by numerical finite
! differences: perturbing T and P by ±0.1% and recomputing the
! equation of state (via COMPUTE_ONE_POP → ENERGY_DENSITY).
!
! Four perturbation calls give the key thermodynamic derivatives:
!   (dE/dT)_P,  (dρ/dT)_P     from T+ and T- at constant P
!   (dE/dP)_T,  (dρ/dP)_T     from P+ and P- at constant T
!
! From these, the routine computes:
!   HEATCP(J)  — specific heat at constant pressure
!   GRDADB(J)  — adiabatic temperature gradient ∇_ad
!   DLTDLP(J)  — actual temperature gradient d(ln T)/d(ln P)
!   VELSND(J)  — adiabatic sound speed
!   HSCALE(J)  — pressure scale height H_P = P / (ρg)
!   DLRDLT(J)  — d(ln ρ)/d(ln T) at constant P
!
! When convection is present (∇_actual > ∇_ad and MIXLTH > 0):
!   - Iterates on the convective element opacity (up to 30 iterations)
!     using the harmonic mean of opacity at T±ΔT
!   - Includes optically thin correction (Mihalas bubble formulation)
!   - Solves the cubic MLT equation via series expansion for numerical
!     stability when the discriminant is small
!   - Smooths the convective flux with a 3-point filter
!   - Optionally applies overshooting (controlled by OVERWT)
!
! Key outputs stored in module arrays:
!   FLXCNV(J)  — convective flux at each depth
!   VCONV(J)   — convective velocity
!   EDENS(J)   — internal energy density per unit mass
!
! References:
!   Böhm-Vitense (1958), Z. Astrophys. 46, 108
!   Mihalas (1978), Stellar Atmospheres (optically thin correction)
!   Castelli (correction to overshooting weight)
!=========================================================================

SUBROUTINE CONVEC

  IMPLICIT NONE

  ! --- Local arrays: finite-difference energy density and density ---
  REAL(8) :: EDENS1(kw), EDENS2(kw)   ! E at T+, T-
  REAL(8) :: EDENS3(kw), EDENS4(kw)   ! E at P+, P-
  REAL(8) :: RHO1(kw),   RHO2(kw)     ! ρ at T+, T-
  REAL(8) :: RHO3(kw),   RHO4(kw)     ! ρ at P+, P-
  REAL(8) :: SAVXNE(kw),  SAVXNA(kw),  SAVRHO(kw)  ! saved state
  REAL(8) :: DILUT(kw)                 ! dilution factor 1 - exp(-τ_Ross)

  REAL(8) :: DTDRHX(kw)    ! dT/d(RHOX) from DERIV
  REAL(8) :: ABCONV(kw)    ! convective opacity (harmonic mean at T±ΔT)
  REAL(8) :: DELTAT(kw)    ! temperature excess of convective element
  REAL(8) :: ROSST(kw)     ! Rosseland opacity at actual T, P
  REAL(8) :: CNVINT(kw)    ! integrated convective flux (for overshooting)
  REAL(8) :: DELHGT(kw)    ! overshooting height increment

  ! --- Local scalars: thermodynamic derivatives ---
  REAL(8)  :: DEDT          ! (dE/dT)_P
  REAL(8)  :: DRDT          ! (dρ/dT)_P
  REAL(8)  :: DEDPG         ! (dE/dP)_T
  REAL(8)  :: DRDPG         ! (dρ/dP)_T
  REAL(8)  :: DPDPG         ! dP_total/dP_gas (= 1 when ignoring P_turb)
  REAL(8)  :: DPDT          ! dP_rad/dT ≈ 4*P_rad/T * dilution
  REAL(8)  :: HEATCV        ! specific heat at constant volume

  ! --- Local scalars: mixing length theory ---
  REAL(8)  :: DEL           ! superadiabatic excess ∇ - ∇_ad
  REAL(8)  :: VCO           ! convective velocity scale
  REAL(8)  :: FLUXCO        ! convective flux coefficient
  REAL(8)  :: D             ! radiative damping parameter
  REAL(8)  :: TAUB          ! optical thickness of convective bubble
  REAL(8)  :: DDD           ! discriminant for cubic equation
  REAL(8)  :: DELTA         ! convective efficiency parameter
  REAL(8)  :: TERM, UP, DOWN  ! series expansion variables
  REAL(8)  :: DPLUS, DMINUS ! opacity ratios at T±ΔT
  REAL(8)  :: OLDDELT       ! previous iteration's DELTAT
  INTEGER :: ITS30         ! max opacity iterations (30 or 1)
  INTEGER :: IDELTAT       ! opacity iteration counter

  ! --- Local scalars: overshooting ---
  REAL(8)  :: WTCNV         ! overshooting weight
  REAL(8)  :: CNV1(1), CNV2(1)    ! interpolated integrated fluxes
  REAL(8)  :: XNEW_CNV(1)
  INTEGER :: M             ! MAP1 return value

  ! --- Other locals ---
  INTEGER :: J
  INTEGER :: JTOP, JBOT, JA, JB  ! gap-filling indices
  REAL(8)  :: WGHT                 ! interpolation weight

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING CONVEC'

  CALL DERIV(RHOX, T, DTDRHX, NRHOX)

  !=====================================================================
  ! Finite-difference thermodynamic derivatives
  ! Compute E and ρ at T±0.1% and P±0.1% to get dE/dT, dρ/dT, etc.
  !=====================================================================
  IFEDNS = 1

  ! Save current state
  DILUT = 1.0D0 - exp(-TAUROS)
  ! DILUT(J) = PRAD(J) / PRADK(J)
  SAVXNE = XNE
  SAVXNA = XNATOM
  SAVRHO = RHO

  ! --- Perturbation 1: T + 0.1% ---
  TLOG = TLOG + 0.0009995003D0
  T    = T * 1.001D0
  TK   = TK * 1.001D0
  HKT  = HKT / 1.001D0
  HCKT = HCKT / 1.001D0
  TKEV = TKEV * 1.001D0
  ITEMP = ITEMP + 1
  CALL COMPUTE_ONE_POP(0.0D0, 1, XNE)

  ! 3*PRADK is approximately RADEN (radiation energy density).
  ! PRADK is used because it can be reconstructed from model decks
  ! whereas RADEN cannot. Rigorously the radiation field should be
  ! recalculated.
  EDENS1 = EDENS + 3.0D0 * PRADK / RHO &
            * (1.0D0 + DILUT * (1.001D0**4 - 1.0D0))
  RHO1 = RHO

  ! --- Perturbation 2: T - 0.1% ---
  TLOG = TLOG - 0.0009995003D0 - 0.0010005003D0
  T    = T / 1.001D0 * 0.999D0
  TK   = TK / 1.001D0 * 0.999D0
  HKT  = HKT * 1.001D0 / 0.999D0
  HCKT = HCKT * 1.001D0 / 0.999D0
  TKEV = TKEV / 1.001D0 * 0.999D0
  ITEMP = ITEMP + 1
  CALL COMPUTE_ONE_POP(0.0D0, 1, XNE)

  EDENS2 = EDENS + 3.0D0 * PRADK / RHO &
            * (1.0D0 + DILUT * (0.999D0**4 - 1.0D0))
  RHO2 = RHO

  ! --- Perturbation 3: P + 0.1% (restore T first) ---
  TLOG = TLOG + 0.0010005003D0
  T    = T / 0.999D0
  TK   = TK / 0.999D0
  HKT  = HKT * 0.999D0
  HCKT = HCKT * 0.999D0
  TKEV = TKEV / 0.999D0
  P    = P * 1.001D0
  ITEMP = ITEMP + 1
  CALL COMPUTE_ONE_POP(0.0D0, 1, XNE)

  EDENS3 = EDENS + 3.0D0 * PRADK / RHO
  RHO3 = RHO
  P = P / 1.001D0 * 0.999D0

  ! --- Perturbation 4: P - 0.1% ---
  ITEMP = ITEMP + 1
  CALL COMPUTE_ONE_POP(0.0D0, 1, XNE)

  EDENS4 = EDENS + 3.0D0 * PRADK / RHO
  RHO4 = RHO
  ! Restore saved state
  XNE    = SAVXNE
  XNATOM = SAVXNA
  RHO    = SAVRHO
  P      = P / 0.999D0

  !=====================================================================
  ! Compute thermodynamic quantities and convective flux at each depth
  !=====================================================================
  DO J = 1, NRHOX

    ! Finite-difference derivatives (factor 500 = 1/(2*0.001))
    DEDT  = (EDENS1(J) - EDENS2(J)) / T(J) * 500.0D0
    DRDT  = (RHO1(J) - RHO2(J)) / T(J) * 500.0D0
    DEDPG = (EDENS3(J) - EDENS4(J)) / P(J) * 500.0D0
    DRDPG = (RHO3(J) - RHO4(J)) / P(J) * 500.0D0

    ! Thermodynamic quantities, ignoring P_turb and assuming P_rad ∝ T⁴
    DPDPG = 1.0D0
    DPDT  = 4.0D0 * PRADK(J) / T(J) * DILUT(J)
    ! DPDT = 4.0D0 * PRAD(J) / T(J)

    ! Actual temperature gradient
    DLTDLP(J) = PTOTAL(J) / T(J) / GRAV * DTDRHX(J)

    ! Specific heats
    HEATCV = DEDT - DEDPG * DRDT / DRDPG
    HEATCP(J) = DEDT - DEDPG * DPDT / DPDPG &
              - PTOTAL(J) / RHO(J)**2 &
              * (DRDT - DRDPG * DPDT / DPDPG)

    ! Sound speed
    IF (HEATCV .LE. 0.0D0) THEN
      VELSND(J) = 0.0D0
    ELSE
      VELSND(J) = sqrt(max(HEATCP(J) / HEATCV * DPDPG / DRDPG, 0.0D0))
    END IF

    ! Density response and adiabatic gradient
    DLRDLT(J) = T(J) / RHO(J) * (DRDT - DRDPG * DPDT / DPDPG)
    GRDADB(J) = -PTOTAL(J) / RHO(J) / T(J) * DLRDLT(J) / HEATCP(J)
    HSCALE(J) = PTOTAL(J) / RHO(J) / GRAV

    ! Initialize convective quantities to zero
    VCONV(J)   = 0.0D0
    FLXCNV(J)  = 0.0D0
    ABCONV(J)  = ABROSS(J)
    DELTAT(J)  = 0.0D0
    ROSST(J)   = 0.0D0

    ! --- Check if convection can operate ---
    IF (MIXLTH .EQ. 0.0D0) CYCLE
    IF (J .LT. 4) CYCLE

    DEL = DLTDLP(J) - GRDADB(J)
    IF (DEL .LT. 0.0D0) CYCLE    ! subadiabatic: no convection

    VCO = 0.5D0 * MIXLTH &
        * sqrt(max(-0.5D0 * PTOTAL(J) / RHO(J) * DLRDLT(J), 0.0D0))
    IF (VCO .EQ. 0.0D0) CYCLE

    FLUXCO = 0.5D0 * RHO(J) * HEATCP(J) * T(J) * MIXLTH / FOURPI
    ROSST(J) = ROSSTAB(T(J), P(J), VTURB(J))
    OLDDELT = 0.0D0

    ! --- Iterate on the convective opacity ---
    ITS30 = 30
    IF (IFCONV .EQ. 0) ITS30 = 1

    DO IDELTAT = 1, ITS30
      DPLUS  = ROSSTAB(T(J) + DELTAT(J), P(J), VTURB(J)) / ROSST(J)
      DMINUS = ROSSTAB(T(J) - DELTAT(J), P(J), VTURB(J)) / ROSST(J)
      ABCONV(J) = 2.0D0 / (1.0D0 / DPLUS + 1.0D0 / DMINUS) * ABROSS(J)

      ! Radiative damping parameter
      D = 8.0D0 * SIGMA_SB * T(J)**4 &
        / (ABCONV(J) * HSCALE(J) * RHO(J)) &
        / (FLUXCO * FOURPI) / VCO

      ! Correction for optically thin bubbles (Mihalas)
      TAUB = ABCONV(J) * RHO(J) * MIXLTH * HSCALE(J)
      D = D * TAUB**2 / (2.0D0 + TAUB**2)

      D = D**2 / 2.0D0
      DDD = (DEL / (D + DEL))**2

      ! Solve for DELTA = convective efficiency
      ! Uses series expansion when DDD < 0.5 for numerical stability
      ! (binomial series for (1 - sqrt(1-x))/x)
      IF (DDD .GE. 0.5D0) THEN
        DELTA = (1.0D0 - sqrt(1.0D0 - DDD)) / DDD
      ELSE
        DELTA = 0.5D0
        TERM  = 0.5D0
        UP    = -1.0D0
        DOWN  = 2.0D0
        DO
          UP   = UP + 2.0D0
          DOWN = DOWN + 2.0D0
          TERM = UP / DOWN * DDD * TERM
          DELTA = DELTA + TERM
          IF (TERM .LE. 1.0D-6) EXIT
        END DO
      END IF
      DELTA = DELTA * DEL**2 / (D + DEL)

      ! Convective velocity and flux
      VCONV(J)  = VCO * sqrt(DELTA)
      FLXCNV(J) = FLUXCO * VCONV(J) * DELTA
      FLXCNV(J) = max(FLXCNV(J), 0.0D0)

      ! Temperature excess of convective element
      DELTAT(J) = T(J) * MIXLTH * DELTA
      DELTAT(J) = min(DELTAT(J), T(J) * 0.15D0)
      DELTAT(J) = DELTAT(J) * 0.7D0 + OLDDELT * 0.3D0

      IF (DELTAT(J) .LT. OLDDELT + 0.5D0 .AND. &
          DELTAT(J) .GT. OLDDELT - 0.5D0) EXIT
      OLDDELT = DELTAT(J)
    END DO  ! IDELTAT
    IF (IDELTAT .GT. ITS30) THEN
      WRITE(6,'(A,I3,A,F8.1,A,F8.1)') &
        ' CONVEC WARNING: opacity iteration did not converge at layer ', &
        J, '  DELTAT=', DELTAT(J), '  OLDDELT=', OLDDELT
    END IF

  END DO  ! J depth loop

  !=====================================================================
  ! Post-processing: smoothing, artifact removal, overshooting
  !=====================================================================

  !--- Commented-out code: eliminate artifactual convection above
  !    the main convection zone ---
  !      DO 730 J=3,NRHOX
  !      K=NRHOX+1-J
  !      IF(FLXCNV(K).GT.0.)GO TO 731
  !  730 CONTINUE
  !      RETURN
  !  731 DO 732 J=1,K
  !      L=K+1-J
  !      IF(FLXCNV(L).EQ.0.)GO TO 733
  !  732 CONTINUE
  !  733 DO 734 J=1,L
  !      VCONV(J)=0.
  !  734 FLXCNV(J)=0.
  !--- End commented-out code ---

  !--- Commented-out patch: remove numerical artifacts in upper half ---
  !      DO 7735 J=1,NRHOX/2
  ! 7735 FLXCNV(J)=0.
  !--- End commented-out code ---

  ! 1-2-1 smoothing of convective flux
  FLXCNV0 = FLXCNV
  DO J = 2, NRHOX - 1
    FLXCNV(J) = 0.25D0 * FLXCNV0(J - 1) + 0.50D0 * FLXCNV0(J) &
              + 0.25D0 * FLXCNV0(J + 1)
  END DO
  ! Asymmetric boundary kernel at deepest layer: 75-25 split
  FLXCNV(NRHOX) = 0.75D0 * FLXCNV0(NRHOX) + 0.25D0 * FLXCNV0(NRHOX - 1)

  ! Fill radiative gaps inside the convection zone.
  ! If layers above and below are both convective, the layer in between
  ! must be too — a radiative pocket inside a convection zone is unphysical.
  ! Find the shallowest and deepest convective layers, then fill any
  ! zero-flux layers between them by interpolating from the boundaries.
  JTOP = 0
  JBOT = 0
  DO J = 1, NRHOX
    IF (FLXCNV(J) .GT. 0.0D0) THEN
      IF (JTOP .EQ. 0) JTOP = J
      JBOT = J
    END IF
  END DO
  IF (JTOP .GT. 0 .AND. JBOT .GT. JTOP + 1) THEN
    DO J = JTOP + 1, JBOT - 1
      IF (FLXCNV(J) .EQ. 0.0D0) THEN
        ! Linear interpolation in log between nearest convective neighbors
        ! Find nearest convective layer above
        JA = J - 1
        DO WHILE (JA .GT. JTOP .AND. FLXCNV(JA) .EQ. 0.0D0)
          JA = JA - 1
        END DO
        ! Find nearest convective layer below
        JB = J + 1
        DO WHILE (JB .LT. JBOT .AND. FLXCNV(JB) .EQ. 0.0D0)
          JB = JB + 1
        END DO
        IF (FLXCNV(JA) .GT. 0.0D0 .AND. FLXCNV(JB) .GT. 0.0D0) THEN
          WGHT = dble(J - JA) / dble(JB - JA)
          FLXCNV(J) = FLXCNV(JA) * (1.0D0 - WGHT) + FLXCNV(JB) * WGHT
        END IF
      END IF
    END DO
  END IF

  FLXCNV0 = FLXCNV

  !=====================================================================
  ! Overshooting: extend convective flux above the formal boundary
  ! Assumes overshooting by 0.5 H_P if convection is strong,
  ! none if weak. Setting OVERWT=0 turns off overshooting entirely.
  !=====================================================================
  IF (OVERWT .GT. 0.0D0) THEN
    ! Find maximum convective-to-total flux ratio
    !     WTCNV = MIN(FLXCNV(NRHOX)/FLUX, 1.D0) * OVERWT
    ! Correction from Fiorella Castelli:
    WTCNV = 0.0D0
    DO J = 1, NRHOX
      WTCNV = max(WTCNV, FLXCNV(J) / FLUX)
    END DO
    WTCNV = min(WTCNV, 1.0D0) * OVERWT

    !      DELHGT(J) = MIN(HSCALE(J)*MIXLTH*0.5D-5, HEIGHT(NRHOX)-HEIGHT(J),
    DELHGT = min(HSCALE * 0.5D-5 * WTCNV, &
                    HEIGHT(NRHOX) - HEIGHT, &
                    HEIGHT - HEIGHT(1))
    !      WRITE(6,775) J, HEIGHT(J), DELHGT(J), CNVINT(J)
    FLXCNV0 = FLXCNV
    FLXCNV1 = 0.0D0

    CALL INTEG(HEIGHT, FLXCNV, CNVINT, NRHOX, 0.0D0)

    DO J = NRHOX / 2, NRHOX - 1
      IF (DELHGT(J) .EQ. 0.0D0) CYCLE
      XNEW_CNV(1) = HEIGHT(J) - DELHGT(J)
      M = MAP1(HEIGHT, CNVINT, NRHOX, XNEW_CNV, CNV1, 1)
      XNEW_CNV(1) = HEIGHT(J) + DELHGT(J)
      M = MAP1(HEIGHT, CNVINT, NRHOX, XNEW_CNV, CNV2, 1)
      FLXCNV1(J) = FLXCNV1(J) + (CNV2(1) - CNV1(1)) / DELHGT(J) / 2.0D0
    END DO

    FLXCNV = max(FLXCNV0, FLXCNV1)
  END IF

  ! Patch to remove numerical artifacts in outermost layers
  DO J = 1, NCONV
    ! DO 7779 J=1,NRHOX/3
    FLXCNV(J) = 0.0D0
  END DO

  RETURN

END SUBROUTINE CONVEC

!=========================================================================
! SUBROUTINE COMPUTE_HEIGHT
!
! Compute the geometric height scale of the atmosphere.
!
! Integrates HEIGHT(J) = integral(dz) = integral((1/rho) dRHOX)
! from the surface inward, since dRHOX = rho * dz. The factor 1e-5
! is a scaling constant to keep heights in convenient units.
!
! The HEIGHT array is used by CONVEC for the overshooting calculation.
!=========================================================================

SUBROUTINE COMPUTE_HEIGHT

  IMPLICIT NONE

  REAL(8) :: RHOINV(kw)
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COMPUTE_HEIGHT'

  RHOINV = 1.0D-5 / RHO

  CALL INTEG(RHOX, RHOINV, HEIGHT, NRHOX, 0.0D0)

  RETURN

END SUBROUTINE COMPUTE_HEIGHT

!=========================================================================
! SUBROUTINE VTURB_VARYDEPTH(VNEW)
!
! Set up a depth-dependent microturbulent velocity profile.
!
! Constructs VTURB(J) as a function of Rosseland optical depth by
! scaling the solar microturbulence profile from the Avrett Model C
! atmosphere. The profile shape (rising from ~0.5 km/s in the deep
! atmosphere to ~1.83 km/s at the surface) is fixed; only the overall
! amplitude is adjusted.
!
! Two modes controlled by VNEW:
!   VNEW = -99e5:  Look up VMAX from a precomputed table of maximum
!                  convective velocities indexed by (log g, Teff).
!                  Table covers log g = -1.0..5.0, Teff = 3000..9000 K.
!   VNEW = other:  Use |VNEW| directly as VMAX.
!
! The depth profile is then:
!   VTURB(J) = V_solar(tau_J) * VMAX / 1.83e5
! where V_solar is the Avrett Model C turbulence profile interpolated
! onto the model's optical depth grid.
!
! Called from READIN when the user specifies a negative VTURB value.
!=========================================================================

SUBROUTINE VTURB_VARYDEPTH(VNEW)

  IMPLICIT NONE

  ! --- Arguments ---
  REAL(8), INTENT(IN) :: VNEW

  ! --- Avrett Solar Model C microturbulence profile ---
  ! 30-point tabulation of v_turb (cm/s) vs log10(tau_Ross)
  INTEGER, PARAMETER :: NSOLAR = 30
  REAL(8), PARAMETER :: VTURB_SOLAR_MAX = 1.83D5  ! peak value (cm/s)

  REAL(8), PARAMETER :: VSTANDARD(30) = (/ &
    0.50D5, 0.50D5, 0.50D5, 0.51D5, 0.52D5, 0.55D5, 0.63D5, &
    0.80D5, 0.90D5, 1.00D5, 1.10D5, 1.20D5, 1.30D5, 1.40D5, &
    1.46D5, 1.52D5, 1.56D5, 1.60D5, 1.64D5, 1.68D5, 1.71D5, &
    1.74D5, 1.76D5, 1.78D5, 1.80D5, 1.81D5, 1.82D5, 1.83D5, &
    1.83D5, 1.83D5 /)

  REAL(8), PARAMETER :: TAUSTANDARD(30) = (/ &
    -20.0D0, -3.0D0, -2.67313D0, -2.49296D0, -2.31296D0, &
    -1.95636D0, -1.60768D0, -1.26699D0, -1.10007D0, -0.93587D0, &
    -0.77416D0, -0.61500D0, -0.45564D0, -0.29176D0, -0.18673D0, &
    -0.07193D0,  0.01186D0,  0.10342D0,  0.20400D0,  0.31605D0, &
     0.44498D0,  0.58875D0,  0.74365D0,  0.90604D0,  1.07181D0, &
     1.23841D0,  1.39979D0,  1.55300D0,  2.00000D0, 10.00000D0 /)

  ! --- VMAX lookup table: maximum convective velocity (km/s) ---
  ! Rows: log g = -1.0, -0.5, 0.0, 0.5, ..., 5.0  (13 values, step 0.5)
  ! Cols: Teff  = 3000, 3250, 3500, ..., 9000 K    (25 values, step 250)
  ! Entries = 0.0 indicate parameter combinations outside the valid range.
  ! Values are flux-weighted max convective velocities from model grids,
  ! smoothed for plausible behavior (range 0-8 km/s).
  INTEGER, PARAMETER :: NG_TAB = 13, NT_TAB = 25
  REAL(8), PARAMETER :: VMAXSTD(NG_TAB, NT_TAB) = reshape( (/ &
    ! Teff=3000  (log g = -1.0 .. 5.0)
    3.3, 3.0, 2.7, 2.4, 2.1, 1.8, 1.3, 0.9, 0.6, 0.3, 0.2, 0.1, 0.1, &
    ! Teff=3250
    4.1, 3.7, 3.3, 2.9, 2.5, 2.1, 1.6, 1.2, 0.9, 0.6, 0.3, 0.2, 0.1, &
    ! Teff=3500
    5.2, 4.6, 4.0, 3.4, 2.9, 2.4, 1.9, 1.5, 1.2, 0.9, 0.6, 0.4, 0.2, &
    ! Teff=3750
    6.3, 5.5, 4.7, 3.9, 3.3, 2.7, 2.2, 1.8, 1.5, 1.2, 0.9, 0.6, 0.4, &
    ! Teff=4000
    7.3, 6.4, 5.5, 4.6, 3.7, 3.1, 2.6, 2.1, 1.8, 1.4, 1.1, 0.8, 0.6, &
    ! Teff=4250
    8.0, 7.7, 6.4, 5.1, 4.2, 3.5, 2.9, 2.4, 2.0, 1.6, 1.3, 1.0, 0.7, &
    ! Teff=4500
    8.0, 8.0, 7.1, 5.7, 4.7, 3.9, 3.2, 2.7, 2.3, 1.9, 1.5, 1.2, 0.9, &
    ! Teff=4750
    8.0, 8.0, 7.9, 6.3, 5.2, 4.3, 3.6, 3.0, 2.5, 2.1, 1.7, 1.4, 1.1, &
    ! Teff=5000
    8.0, 8.0, 8.0, 6.9, 5.6, 4.7, 4.0, 3.4, 2.8, 2.3, 1.9, 1.5, 1.2, &
    ! Teff=5250
    0.0, 8.0, 8.0, 7.5, 6.1, 5.1, 4.4, 3.7, 3.1, 2.6, 2.1, 1.7, 1.3, &
    ! Teff=5500
    0.0, 8.0, 8.0, 8.0, 6.6, 5.5, 4.7, 4.0, 3.4, 2.8, 2.3, 1.9, 1.5, &
    ! Teff=5750
    0.0, 0.0, 8.0, 8.0, 7.1, 5.9, 5.0, 4.3, 3.6, 3.0, 2.5, 2.0, 1.7, &
    ! Teff=6000
    0.0, 0.0, 8.0, 8.0, 7.6, 6.2, 5.4, 4.6, 3.9, 3.3, 2.7, 2.2, 1.9, &
    ! Teff=6250
    0.0, 0.0, 0.0, 4.6, 8.0, 6.6, 5.7, 4.9, 4.2, 3.5, 2.9, 2.4, 2.0, &
    ! Teff=6500
    0.0, 0.0, 0.0, 0.2, 4.3, 7.0, 6.1, 5.3, 4.4, 3.7, 3.1, 2.6, 2.2, &
    ! Teff=6750
    0.0, 0.0, 0.0, 0.0, 0.3, 4.2, 6.4, 5.6, 4.7, 4.0, 3.3, 2.8, 2.4, &
    ! Teff=7000
    0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 3.9, 5.9, 5.0, 4.2, 3.5, 3.0, 2.6, &
    ! Teff=7250
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 3.7, 5.2, 4.4, 3.7, 3.2, 2.8, &
    ! Teff=7500
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 3.5, 4.7, 3.9, 3.4, 3.0, &
    ! Teff=7750
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 3.4, 4.1, 3.6, 3.1, &
    ! Teff=8000
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 3.6, 3.8, 3.3, &
    ! Teff=8250
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 4.0, 3.5, &
    ! Teff=8500
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6, 3.6, &
    ! Teff=8750
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, &
    ! Teff=9000
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  &
    /), (/ NG_TAB, NT_TAB /) )

  ! --- Local variables ---
  REAL(8)  :: VMAX           ! amplitude scale for the profile (cm/s)
  REAL(8)  :: TAULOG(kw)     ! log10(tau_Ross) at each depth
  REAL(8)  :: DELG, DELT     ! bilinear interpolation weights
  INTEGER :: IG, IT         ! table indices
  INTEGER :: J, MM

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING VTURB_VARYDEPTH'

  !---------------------------------------------------------------------
  ! Determine VMAX: either from table or from user input
  !---------------------------------------------------------------------
  IF (VNEW .EQ. -99.0D5) THEN
    ! Look up from (log g, Teff) table with bilinear interpolation
    IG = int((GLOG + 1.0D0) / 0.5D0) + 1
    IG = min(max(IG, 1), NG_TAB - 1)
    IT = int((TEFF - 3000.0D0) / 250.0D0) + 1
    IT = min(max(IT, 1), NT_TAB - 1)

    DELG = (GLOG - (dble(IG - 1) * 0.5D0 - 1.0D0)) / 0.5D0
    DELT = (TEFF - (dble(IT - 1) * 250.0D0 + 3000.0D0)) / 250.0D0

    VMAX = VMAXSTD(IG,   IT  ) * (1.0D0 - DELG) * (1.0D0 - DELT) &
         + VMAXSTD(IG+1, IT  ) * DELG            * (1.0D0 - DELT) &
         + VMAXSTD(IG,   IT+1) * (1.0D0 - DELG) * DELT            &
         + VMAXSTD(IG+1, IT+1) * DELG            * DELT
    VMAX = VMAX * 1.0D5   ! convert km/s → cm/s
  ELSE
    VMAX = abs(VNEW)
  END IF

  !---------------------------------------------------------------------
  ! Build depth-dependent profile by interpolating the solar template
  ! onto the model's optical depth grid, then scaling by VMAX
  !---------------------------------------------------------------------
  TAULOG = log10(TAUSTD)

  MM = MAP1(TAUSTANDARD, VSTANDARD, NSOLAR, TAULOG, VTURB, NRHOX)

  VTURB = VTURB * VMAX / VTURB_SOLAR_MAX

  RETURN

END SUBROUTINE VTURB_VARYDEPTH

!=========================================================================
! SUBROUTINE COMPUTE_PTURB
!
! Compute microturbulent velocity and turbulent pressure at each depth.
!
! The microturbulent velocity is constructed from three adjustable
! components controlled by input parameters:
!
!   VTURB(J) = TRBFDG * rho^TRBPOW  +  TRBSND * V_sound  +  TRBCON
!                (density term)         (sound speed term)   (constant)
!
! all in km/s, then converted to cm/s. This allows a flexible
! parameterization: density-dependent, sound-speed-scaled, or constant
! microturbulence (or any combination).
!
! The turbulent pressure is then:
!   PTURB(J) = 0.5 * rho * VTURB^2
!
! which enters the hydrostatic equilibrium as an addition to P_gas.
!
! Called from the main iteration loop when IFTURB = 1.
!=========================================================================

SUBROUTINE COMPUTE_PTURB

  IMPLICIT NONE

  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COMPUTE_PTURB'

  VTURB = (TRBFDG * RHO**TRBPOW &
            + TRBSND * VELSND / 1.0D5 &
            + TRBCON) * 1.0D5
  PTURB = RHO * VTURB**2 * 0.5D0

  RETURN

END SUBROUTINE COMPUTE_PTURB

!=========================================================================
! SUBROUTINE KAPP
!
! Assemble total monochromatic opacity from all continuum and line sources.
!
! At the current frequency (set by BNU, FREQ, etc. before calling KAPP),
! this routine:
!   1. Zeros all per-source opacity arrays
!   2. Calls each opacity source
!   3. Assembles the total true absorption (ACONT, ALINE), scattering
!      (SIGMAC, SIGMAL), and source functions (SCONT, SLINE) from
!      the source-weighted contributions.
!
! The total opacity at each depth is:
!   kappa_total = ACONT + ALINE + SIGMAC + SIGMAL
! with source function formed by weighting each source's S_nu by its
! contribution to the absorption coefficient.
!=========================================================================

SUBROUTINE KAPP

  IMPLICIT NONE

  INTEGER  :: J
  REAL(8)  :: A ! sum of non-H, non-Hminus, non-extra continuum absorption

  ! Zero all per-source opacity arrays
  AHYD   = 0.0D0
  AH2P   = 0.0D0
  AHMIN  = 0.0D0
  SIGH   = 0.0D0
  AHE1   = 0.0D0
  AHE2   = 0.0D0
  AHEMIN = 0.0D0
  SIGHE  = 0.0D0
  ACOOL  = 0.0D0
  ALUKE  = 0.0D0
  AHOT   = 0.0D0
  AH2COLL = 0.0D0
  ACONT_METAL = 0.0D0
  SIGEL  = 0.0D0
  SIGH2  = 0.0D0
  AHLINE = 0.0D0
  ALINES = 0.0D0
  SIGLIN = 0.0D0
  AXLINE = 0.0D0
  SIGXL  = 0.0D0
  SIGX   = 0.0D0
  SHYD   = 0.0D0
  SHMIN  = 0.0D0
  SHLINE = 0.0D0
  SXLINE = 0.0D0
  SXCONT = 0.0D0
  SAL1   = 0.0D0
  SFE1   = 0.0D0
  SHE2   = 0.0D0
  SC1    = 0.0D0

  ! Call each opacity source
  CALL HOP
  CALL H2PLOP
  CALL HMINOP
  CALL HRAYOP
  CALL HE1OP
  CALL HE2OP
  CALL HEMIOP
  CALL HERAOP
  IF (USE_TOPBASE_MBF) THEN
    CALL CONT_METAL_OPACITY_TOPBASE
  ELSE
    CALL CONT_METAL_OPACITY_LEGACY
  ENDIF
  CALL ELECOP
  CALL H2RAOP
  CALL HLINOP

  ! Assemble totals: true absorption, source function, and scattering.
  !
  ! The continuum source function follows the F77 atlas12.for KAPP
  ! treatment: only hydrogen (SHYD), H-minus (SHMIN), and the auxiliary
  ! continuum (SXCONT) carry distinct source functions.  All other
  ! continuum opacity sources (HE I, HE II, the metals C/Mg/Al/Si/Fe,
  ! H2+, He-minus, LUKE/WARM, and HOT) are bundled into the "A" term
  ! and weighted by BNU(J).  This is the LTE limit and matches what
  ! atlas12.for has always done.
  !
  ! The atlas7lib.for KAPP carries an elaborate per-species source-
  ! function formula (SHE1, SHE2, SC1, SMG1, SSI1, SAL1, SFE1), but no
  ! F77 program in the SYNTHE/ATLAS12 pipeline ever consumed it: the
  ! F77 SPECTRV unconditionally sets SCONT(J) = BNU(J) at every
  ! wavelength point and never calls KAPP.  The F90 SYNTHE driver
  ! preserves this behaviour via the contabs_sv tabulation path.
  ! Therefore the only consumer of KAPP's SCONT in this codebase is
  ! the ATLAS12 driver, and the atlas12.for treatment is the correct
  ! reference.
  !
  ! Note on the absorption sum: the metal, molecular, and CIA continuum
  ! contributions are consolidated into ACONT_METAL(J) by the routine
  ! CONT_METAL_OPACITY_LEGACY.  This bundle replaces the sum that was
  ! previously assembled inline here (AC1 + AMG1 + AAL1 + ASI1 + AFE1
  ! + ALUKE + AHOT + AH2COLL + CH/OH molecular), and exists so that a
  ! TOPbase-based metal continuum variant can be swapped in as a
  ! drop-in replacement (CONT_METAL_OPACITY_TOPBASE, phase 3 of the
  ! F90 modernization).  The individual per-species arrays (AC1, AMG1,
  ! etc.) remain populated and available for diagnostics.
  DO J = 1, NRHOX
    ! Sources weighted by BNU in the source function
    A = AH2P(J) + AHE1(J) + AHE2(J) + AHEMIN(J) + ACONT_METAL(J)

    ! Total continuum absorption
    ACONT(J) = A + AHYD(J) + AHMIN(J)

    ! Continuum source function (atlas12.for form)
    SCONT(J) = BNU(J)
    IF (ACONT(J) .GT. 0.0D0) THEN
      SCONT(J) = (A * BNU(J) + AHYD(J) * SHYD(J) &
               + AHMIN(J) * SHMIN(J)) / ACONT(J)
    END IF

    ! Line absorption and source function
    ALINE(J) = AHLINE(J) + ALINES(J) + AXLINE(J)
    SLINE(J) = BNU(J)
    IF (ALINE(J) .GT. 0.0D0) THEN
      SLINE(J) = (AHLINE(J) * SHLINE(J) + ALINES(J) * BNU(J) &
               + AXLINE(J) * SXLINE(J)) / ALINE(J)
    END IF

    ! Scattering: continuum and line
    SIGMAC(J) = SIGH(J) + SIGHE(J) + SIGEL(J) + SIGH2(J) + SIGX(J)
    SIGMAL(J) = SIGLIN(J) + SIGXL(J)
  END DO

  RETURN

END SUBROUTINE KAPP

!=========================================================================
! SUBROUTINE HOP
!
! Hydrogen bound-free and free-free opacity.
!
! Computes hydrogen continuum absorption and source function at the
! current frequency for all depth points.  Two contributions:
!
!   1. Bound-free from explicit levels n=1-15, with REAL*8 wavenumber
!      threshold checks for each level.  Each level is weighted by its
!      occupation probability w_n from the Hummer & Mihalas (1988)
!      formalism.  Levels 1-6 include non-LTE departure coefficients
!      BHYD(J,n).  Levels 7-15 are in LTE.
!
!   2. Free-free (bremsstrahlung).
!
! The former "dissolved levels n >= 16" integral has been removed.
! That integral assumed a fixed cutoff at n=16 regardless of density.
! With occupation probabilities, levels 7-15 are smoothly weighted
! and levels above 15 contribute negligibly to the bound-free opacity
! because their w_n is small and their cross-sections are tiny.
!
! The population normalization XNFP(J,1)/RHO(J) is applied once at
! the end (after summing all bound-free terms), matching the F77
! structure exactly.
!
! CRITICAL: The wavenumber threshold checks must be at REAL*8 precision
! (via WAVENO) to produce sharp continuum edges.
!
! References:
!   Hummer, D.G. & Mihalas, D. 1988, ApJ 331, 794
!=========================================================================

SUBROUTINE HOP

  IMPLICIT NONE

  ! Maximum hydrogen level included in HOP.  The explicit-level loop runs
  ! n=1..15 (with full Karzas-Latter cross sections); a high-n extension
  ! loop covers n=16..NMAX_HLEVELS using XKARZAS's hydrogenic rescaling
  ! from the n=15 table; and the dissolved-fraction pseudo-continuum
  ! covers n=7..NMAX_HLEVELS.  All three loops share the same upper bound
  ! so the bookkeeping stays consistent.
  !
  ! NMAX_HLEVELS = 80 matches NMAX_OCC in HLINOP.  The motivation is cool,
  ! low-density atmospheres (M dwarfs, n_e ~ 1e10): at those densities,
  ! levels up to n ~ 40-60 still have non-negligible occupation probability
  ! w_n, and their bound-free contribution near the Paschen/Brackett jumps
  ! is real opacity that should not be dropped.  For hot/dense atmospheres
  ! (B/O stars, WDs), levels with n >> 20 are fully dissolved and the
  ! w_n < 1.0D-6 guard in the high-n loop cycles immediately, so the
  ! larger cap costs essentially nothing there.
  INTEGER, PARAMETER :: NMAX_HLEVELS = 80

  REAL(8)  :: X            ! bound-free cross-section from XKARZAS
  REAL(8)  :: A            ! single-level opacity contribution
  REAL(8)  :: H, S         ! running sums for opacity and source function
  REAL(8)  :: w_n          ! occupation probability for current level
  REAL(8)  :: E_n          ! excitation energy for the current level [cm^-1]
  REAL(8)  :: thresh_n     ! ionization threshold of the current level [cm^-1]
  REAL(8)  :: FREQ3_INV    ! 1 / FREQ^3 (precomputed for pseudo-continuum)
  INTEGER :: J, I_PC, n_hi

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HOP'

  DO J = 1, NRHOX

    H = 0.0D0
    S = 0.0D0

    ! ------------------------------------------------------------------
    ! Explicit levels n=15 down to n=1.
    ! Cascading: exit as soon as WAVENO is below a threshold, skipping
    ! that level and all lower levels (higher thresholds).
    ! Each level is weighted by its occupation probability w_n.
    ! Levels 7-15: LTE (use STIM), weighted by w_n
    ! Levels 1-6:  non-LTE (use BHYD departure coefficients), w_n ~ 1
    ! ------------------------------------------------------------------
    levels: DO

      ! n=15  threshold =    487.456 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(15)) EXIT levels
      w_n = occupation_prob(15, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 15, 15)
      A = w_n * X * 450.0D0 * EXP(-109191.313D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=14  threshold =    559.579 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(14)) EXIT levels
      w_n = occupation_prob(14, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 14, 14)
      A = w_n * X * 392.0D0 * EXP(-109119.188D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=13  threshold =    648.980 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(13)) EXIT levels
      w_n = occupation_prob(13, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 13, 13)
      A = w_n * X * 338.0D0 * EXP(-109029.789D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=12  threshold =    761.649 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(12)) EXIT levels
      w_n = occupation_prob(12, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 12, 12)
      A = w_n * X * 288.0D0 * EXP(-108917.117D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=11  threshold =    906.426 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(11)) EXIT levels
      w_n = occupation_prob(11, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 11, 11)
      A = w_n * X * 242.0D0 * EXP(-108772.336D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=10  threshold =   1096.776 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(10)) EXIT levels
      w_n = occupation_prob(10, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 10, 10)
      A = w_n * X * 200.0D0 * EXP(-108581.992D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=9   threshold =   1354.044 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(9)) EXIT levels
      w_n = occupation_prob(9, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 9, 9)
      A = w_n * X * 162.0D0 * EXP(-108324.719D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=8   threshold =   1713.713 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(8)) EXIT levels
      w_n = occupation_prob(8, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 8, 8)
      A = w_n * X * 128.0D0 * EXP(-107965.051D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=7   threshold =   2238.320 cm^{-1}
      IF (WAVENO .LT. H_SERIES_LIMITS(7)) EXIT levels
      w_n = occupation_prob(7, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 7, 7)
      A = w_n * X * 98.0D0 * EXP(-107440.444D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=6   threshold =   3046.604 cm^{-1}   (non-LTE)
      IF (WAVENO .LT. H_SERIES_LIMITS(6)) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 6, 6)
      A = X * 72.0D0 * EXP(-106632.160D0 * HCKT(J)) * (BHYD(J,6) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,6) - EHVKT(J))

      ! n=5   threshold =   4387.113 cm^{-1}   (non-LTE)
      IF (WAVENO .LT. H_SERIES_LIMITS(5)) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 5, 5)
      A = X * 50.0D0 * EXP(-105291.651D0 * HCKT(J)) * (BHYD(J,5) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,5) - EHVKT(J))

      ! n=4   threshold =   6854.871 cm^{-1}   (non-LTE)
      IF (WAVENO .LT. H_SERIES_LIMITS(4)) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 4, 4)
      A = X * 32.0D0 * EXP(-102823.893D0 * HCKT(J)) * (BHYD(J,4) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,4) - EHVKT(J))

      ! n=3   threshold =  12186.462 cm^{-1}   (non-LTE)
      IF (WAVENO .LT. H_SERIES_LIMITS(3)) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 3, 3)
      A = X * 18.0D0 * EXP(-97492.302D0 * HCKT(J)) * (BHYD(J,3) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,3) - EHVKT(J))

      ! n=2   threshold =  27419.659 cm^{-1}   (non-LTE)  [Balmer limit]
      IF (WAVENO .LT. H_SERIES_LIMITS(2)) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 2, 2)
      A = X * 8.0D0 * EXP(-82259.105D0 * HCKT(J)) * (BHYD(J,2) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,2) - EHVKT(J))

      ! n=1   threshold = 109678.764 cm^{-1}   (non-LTE)  [Lyman limit]
      IF (WAVENO .LT. H_SERIES_LIMITS(1)) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 1, 1)
      A = X * 2.0D0 * 1.0D0 * (BHYD(J,1) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,1) - EHVKT(J))

      EXIT levels  ! all levels processed
    END DO levels

    ! ------------------------------------------------------------------
    ! High-n explicit bound-free extension: levels 16 to NMAX_HLEVELS.
    !
    ! These levels are above the principal-quantum-number range covered
    ! by the explicit Karzas-Latter table (n <= 15), but XKARZAS handles
    ! n > 15 by rescaling the n=15 tabulation using the hydrogenic
    ! frequency scaling.  At low and moderate density (where w_n ~ 1
    ! for these levels), this loop adds the LTE bound-free contribution
    ! from levels 16-30 that would otherwise be missing entirely.
    !
    ! In the original F77 atlas12.for, this contribution was included
    ! via an analytic Boltzmann integral over n=9..infinity using pure
    ! Kramers cross sections (the "EXLIM/BOLTEX" term).  The F90 was
    ! restructured to use explicit Karzas-Latter levels through n=15
    ! plus the dissolved-fraction pseudo-continuum, which is more
    ! accurate for n=9..15 but inadvertently dropped the n>=16 LTE
    ! tail.  This loop restores that contribution.
    !
    ! Each level is weighted by w_n; the dissolved fraction (1-w_n) is
    ! handled by the pseudo-continuum loop below.
    ! ------------------------------------------------------------------
    DO n_hi = 16, NMAX_HLEVELS
      thresh_n = ELIM_HI / dble(n_hi)**2
      IF (WAVENO .LT. thresh_n) CYCLE              ! below this level's threshold
      w_n = occupation_prob(n_hi, XNE(J))
      IF (w_n .LT. 1.0D-6) CYCLE                   ! fully dissolved -- skip
      X = XKARZAS(FREQ, 1.0D0, n_hi, n_hi)
      E_n = ELIM_HI - RYDBERG_H / dble(n_hi)**2
      A = w_n * X * 2.0D0 * dble(n_hi)**2 * EXP(-E_n * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)
    END DO

    ! ------------------------------------------------------------------
    ! Pseudo-continuum from dissolved hydrogen levels.
    !
    ! For each level n where the occupation probability w_n < 1, the
    ! dissolved fraction (1 - w_n) contributes a smooth Kramers-like
    ! opacity at ALL frequencies.  This replaces the former fixed
    ! "dissolved levels n >= 16" Boltzmann integral.
    !
    ! Above threshold, this term restores the opacity reduced by the
    ! w_n weighting of the explicit bound-free levels, so the total
    ! above-threshold opacity is approximately unchanged.
    ! Below threshold, this is new opacity that smooths the region
    ! near series limits (the Balmer jump, Paschen jump, etc.).
    !
    ! Cross-section: Kramers formula sigma = 2.815e29 / (n^5 * nu^3)
    ! Population: g_n * exp(-E_n * hc/kT) = 2*n^2 * exp(-E_n * hc/kT)
    ! Combined: 2.815e29 * 2 / (n^3 * nu^3) * exp(-E_n * hc/kT)
    !
    ! References:
    !   Hubeny, I., Hummer, D.G. & Lanz, T. 1994, A&A 282, 151
    !   Tremblay, P.-E. & Bergeron, P. 2009, ApJ 696, 1755
    ! ------------------------------------------------------------------
    FREQ3_INV = 1.0D0 / (FREQ * FREQ * FREQ)
    DO I_PC = 7, NMAX_HLEVELS
      w_n = occupation_prob(I_PC, XNE(J))
      IF (1.0D0 - w_n .LT. 1.0D-6) CYCLE   ! fully bound, no dissolved fraction
      E_n = ELIM_HI - RYDBERG_H / dble(I_PC)**2
      A = (1.0D0 - w_n) * 2.815D29 * FREQ3_INV &
        * 2.0D0 / dble(I_PC)**3 * EXP(-E_n * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)
    END DO

    ! ------------------------------------------------------------------
    ! Apply population normalization and add free-free
    ! ------------------------------------------------------------------
    H = H * XNFP(J, 1) / RHO(J)
    S = S * XNFP(J, 1) / RHO(J)

    ! Free-free: kappa_ff = 3.6919e8 / sqrt(T) * g_ff / nu^3 * n_e * n_H+ * STIM / rho
    A = COEFF_FF / SQRT(T(J)) * COULFF(J,1) / FREQ * XNE(J) / FREQ &
      * XNFP(J, 2) / FREQ * STIM(J) / RHO(J)
    H = H + A
    S = S + A * BNU(J)

    AHYD(J) = H
    IF (H .GT. 0.0D0) SHYD(J) = S / H

  END DO

  RETURN

END SUBROUTINE HOP


!=========================================================================
! FUNCTION XKARZAS(FREQ, ZEFF2, N, L)
!
! Hydrogen bound-free photoionization cross-section from Karzas & Latter.
!
! Returns the bound-free cross-section (cm^2) for a hydrogen-like ion
! at frequency FREQ from level (N, L), scaled by 1/Z^2.
!
! Three regimes:
!   1. N <= 15, L >= N (or N > 6):  Level-averaged cross-section from
!      tabulated XN(29,15). Interpolated in log(sigma) vs log(freq).
!   2. N <= 6, L < N:  L-resolved cross-section from tabulated
!      XL(29,6,6). Uses the correct angular momentum substate.
!   3. N > 15:  Rescales the N=15 tabulation using hydrogenic frequency
!      scaling (nu_edge proportional to 1/n^2).
!
! Tables are stored in external data files (read on first call from
! DATADIR) to keep the source code clean. The tables use REAL*4
! precision as in the original Karzas & Latter tabulation.
!
! Reference: Karzas, W.J. and Latter, R. 1961, ApJS 6, 167-212
! Bug fixes: 24 Jul 2003 — five corrected values in X11-X13 data
!=========================================================================

FUNCTION XKARZAS(FREQ, ZEFF2, N, L)

  IMPLICIT NONE

  ! --- Arguments ---
  REAL(8),  INTENT(IN) :: FREQ    ! frequency (Hz)
  REAL(8),  INTENT(IN) :: ZEFF2   ! effective nuclear charge squared
  INTEGER, INTENT(IN) :: N       ! principal quantum number
  INTEGER, INTENT(IN) :: L       ! orbital quantum number

  ! --- Return value ---
  REAL(8) :: XKARZAS

  ! --- Tabulated cross-section data (read from files on first call) ---
  ! Tables use REAL*4 as in the original Karzas & Latter tabulation
  INTEGER, PARAMETER :: NPTS = 29   ! energy grid points per level
  INTEGER, PARAMETER :: NMAX = 15   ! max tabulated level for XN
  INTEGER, PARAMETER :: NLMAX = 6   ! max tabulated level for XL (l-resolved)

  REAL(4), SAVE :: FREQN(NPTS, NMAX)     ! log10(freq/Z^2) grid per level
  REAL(4), SAVE :: XN(NPTS, NMAX)        ! log10(sigma), level-averaged
  REAL(4), SAVE :: XL(NPTS, NLMAX, NLMAX) ! log10(sigma), l-resolved
  REAL(4), SAVE :: EKARZAS(NPTS)          ! energy grid for n>15 scaling
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! --- Local variables ---
  REAL(4)  :: FREQLG             ! log10(freq/Z^2)
  REAL(4)  :: X                  ! interpolated log10(sigma)
  REAL(4)  :: FREQN15(NPTS)     ! frequency grid for n>15 (constructed)
  INTEGER :: I

  ! --- First-call initialization: read tables from data files ---
  IF (.NOT. INITIALIZED) THEN
    CALL read_karzas_tables()
    INITIALIZED = .TRUE.
  END IF

  ! --- Compute cross-section ---
  FREQLG = REAL(log10(FREQ / ZEFF2))
  XKARZAS = 0.0D0

  IF (L .GE. N .OR. N .GT. NLMAX) THEN
    !-----------------------------------------------------------------
    ! Level-averaged cross-section (L >= N means use averaged tables,
    ! or N > 6 where only averaged tables are available)
    !-----------------------------------------------------------------
    IF (N .GT. NMAX) THEN
      ! N > 15: rescale from N=15 using hydrogenic frequency scaling
      FREQN15(NPTS) = log10(FREQ_RYDH / REAL(N)**2)
      IF (FREQLG .LT. FREQN15(NPTS)) RETURN
      DO I = 2, NPTS - 1
        FREQN15(I) = log10((EKARZAS(I) + 1.0 / REAL(N)**2) &
                   * FREQ_RYDH)
        IF (FREQLG .GT. FREQN15(I)) EXIT
      END DO
      IF (I .GT. NPTS - 1) I = NPTS
      X = (FREQLG - FREQN15(I)) / (FREQN15(I-1) - FREQN15(I)) &
        * (XN(I-1, NMAX) - XN(I, NMAX)) + XN(I, NMAX)
    ELSE
      ! N = 1-15: direct table lookup
      IF (FREQLG .LT. FREQN(NPTS, N)) RETURN
      DO I = 2, NPTS
        IF (FREQLG .GT. FREQN(I, N)) EXIT
      END DO
      X = (FREQLG - FREQN(I, N)) / (FREQN(I-1, N) - FREQN(I, N)) &
        * (XN(I-1, N) - XN(I, N)) + XN(I, N)
    END IF

  ELSE
    !-----------------------------------------------------------------
    ! L-resolved cross-section (N <= 6, L < N)
    !-----------------------------------------------------------------
    IF (FREQLG .LT. FREQN(NPTS, N)) RETURN
    DO I = 2, NPTS
      IF (FREQLG .GT. FREQN(I, N)) EXIT
    END DO
    X = (FREQLG - FREQN(I, N)) / (FREQN(I-1, N) - FREQN(I, N)) &
      * (XL(I-1, L+1, N) - XL(I, L+1, N)) + XL(I, L+1, N)
  END IF

  XKARZAS = exp(dble(X) * LN10) / ZEFF2

  RETURN

CONTAINS

  !-------------------------------------------------------------------
  ! Read Karzas & Latter cross-section tables from data files
  !-------------------------------------------------------------------
  SUBROUTINE read_karzas_tables()
    INTEGER :: I, J, K, IU, IOS
    CHARACTER(256) :: FILEPATH

    IU = 89  ! scratch unit (same convention as other data reads)

    ! Read EKARZAS (energy grid)
    FILEPATH = trim(DATADIR) // 'karzas_ekarzas.dat'
    OPEN(UNIT=IU, FILE=FILEPATH, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      STOP 'XKARZAS: cannot read Karzas-Latter data'
    END IF
    READ(IU, '(A)') FILEPATH  ! skip header
    READ(IU, '(A)') FILEPATH  ! skip header
    DO I = 1, NPTS
      READ(IU, *) EKARZAS(I)
    END DO
    CLOSE(IU)

    ! Read XN (level-averaged cross-sections, column-major)
    FILEPATH = trim(DATADIR) // 'karzas_xn.dat'
    OPEN(UNIT=IU, FILE=FILEPATH, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      STOP 'XKARZAS: cannot read Karzas-Latter data'
    END IF
    READ(IU, '(A)') FILEPATH  ! skip header
    READ(IU, '(A)') FILEPATH  ! skip header
    READ(IU, '(A)') FILEPATH  ! skip header
    DO J = 1, NMAX
      DO I = 1, NPTS
        READ(IU, *) XN(I, J)
      END DO
    END DO
    CLOSE(IU)

    ! Read FREQN (frequency grids, column-major)
    FILEPATH = trim(DATADIR) // 'karzas_freqn.dat'
    OPEN(UNIT=IU, FILE=FILEPATH, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      STOP 'XKARZAS: cannot read Karzas-Latter data'
    END IF
    READ(IU, '(A)') FILEPATH  ! skip header
    READ(IU, '(A)') FILEPATH  ! skip header
    DO J = 1, NMAX
      DO I = 1, NPTS
        READ(IU, *) FREQN(I, J)
      END DO
    END DO
    CLOSE(IU)

    ! Read XL (l-resolved cross-sections, column-major)
    FILEPATH = trim(DATADIR) // 'karzas_xl.dat'
    OPEN(UNIT=IU, FILE=FILEPATH, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      STOP 'XKARZAS: cannot read Karzas-Latter data'
    END IF
    READ(IU, '(A)') FILEPATH  ! skip header
    READ(IU, '(A)') FILEPATH  ! skip header
    READ(IU, '(A)') FILEPATH  ! skip header
    DO K = 1, NLMAX
      DO J = 1, NLMAX
        DO I = 1, NPTS
          READ(IU, *) XL(I, J, K)
        END DO
      END DO
    END DO
    CLOSE(IU)

  END SUBROUTINE read_karzas_tables

END FUNCTION XKARZAS

!=========================================================================
! FUNCTION COULX(N, FREQ, Z)
!
! Hydrogenic bound-free photoionization cross-section (Kramers formula).
!
! Returns sigma_bf (cm^2) for ionization from level N of a hydrogen-like
! ion with nuclear charge Z at frequency FREQ.
!
! The base formula is Kramers:
!   sigma = 2.815e29 / (nu^3 * n^5) * Z^4
! with polynomial Gaunt factor corrections for N=2-6 (fitted to
! Karzas & Latter tables), and the exact 1s Gaunt factor from
! COULBF1S for N=1. For N>6, the uncorrected Kramers formula is used.
!
! Returns 0 if FREQ is below the ionization threshold for level N.
!=========================================================================

FUNCTION COULX(N, FREQ, Z)

  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN) :: N       ! principal quantum number
  REAL(8),  INTENT(IN) :: FREQ    ! frequency (Hz)
  REAL(8),  INTENT(IN) :: Z       ! nuclear charge

  ! --- Return value ---
  REAL(8) :: COULX

  ! --- Gaunt factor polynomial coefficients for N=1-6 ---
  ! sigma = Kramers * (A + (B + C*(Z^2/nu)) * (Z^2/nu))
  REAL(8), PARAMETER :: A(6) = (/ 0.9916D0, 1.105D0, 1.101D0, &
                                  1.101D0, 1.102D0, 1.0986D0 /)
  REAL(8), PARAMETER :: B(6) = (/ 2.719D13, -2.375D14, -9.863D13, &
                                  -5.765D13, -3.909D13, -2.704D13 /)
  REAL(8), PARAMETER :: C(6) = (/ -2.268D30, 4.077D28, 1.035D28, &
                                  4.593D27, 2.371D27, 1.229D27 /)

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COULX'

  ! Check if frequency is above the ionization threshold: nu_edge = Z^2 * R_inf * c / n^2
  IF (FREQ .LT. Z * Z * FREQ_RYDH / dble(N)**2) THEN
    COULX = 0.0D0
    RETURN
  END IF

  ! Kramers cross-section
  COULX = 2.815D29 / FREQ**3 / dble(N)**5 * Z**4

  ! Apply Gaunt factor correction
  IF (N .LE. 6) THEN
    IF (N .EQ. 1) THEN
      ! Exact 1s Gaunt factor from tabulation
      COULX = COULX * COULBF1S(FREQ, Z)
    ELSE
      ! Polynomial fit for N=2-6
      COULX = COULX * (A(N) + (B(N) + C(N) * (Z * Z / FREQ)) * (Z * Z / FREQ))
    END IF
  END IF

  RETURN

END FUNCTION COULX

!=========================================================================
! FUNCTION COULBF1S(FREQ, Z)
!
! Bound-free Gaunt factor for the hydrogen 1s ground state.
!
! Returns the exact Gaunt factor g_bf(1s) by linear interpolation in a
! 151-point table covering log10(nu / nu_edge) = 0.00 to 3.00 in steps
! of 0.02. The tabulated values range from ~0.80 near threshold to
! ~0.20 at high energies (E/E_edge ~ 1000).
!
! Returns 0 if frequency is below the ionization threshold.
!=========================================================================

FUNCTION COULBF1S(FREQ, Z)

  IMPLICIT NONE

  ! --- Arguments ---
  REAL(8), INTENT(IN) :: FREQ     ! frequency (Hz)
  REAL(8), INTENT(IN) :: Z        ! nuclear charge

  ! --- Return value ---
  REAL(8) :: COULBF1S

  ! --- 1s Gaunt factor table ---
  ! 151 points at log10(nu/nu_edge) = 0.00, 0.02, 0.04, ..., 3.00
  INTEGER, PARAMETER :: NGAUNT = 151
  REAL(8), PARAMETER :: DLOG_STEP = 0.02D0

  REAL(8), PARAMETER :: GAUNT1S(151) = (/ &
    0.7973D0, 0.8094D0, 0.8212D0, 0.8328D0, 0.8439D0, &
    0.8548D0, 0.8653D0, 0.8754D0, 0.8852D0, 0.8946D0, &
    0.9035D0, 0.9120D0, 0.9201D0, 0.9278D0, 0.9351D0, &
    0.9420D0, 0.9484D0, 0.9544D0, 0.9601D0, 0.9653D0, &
    0.9702D0, 0.9745D0, 0.9785D0, 0.9820D0, 0.9852D0, &
    0.9879D0, 0.9903D0, 0.9922D0, 0.9938D0, 0.9949D0, &
    0.9957D0, 0.9960D0, 0.9960D0, 0.9957D0, 0.9949D0, &
    0.9938D0, 0.9923D0, 0.9905D0, 0.9884D0, 0.9859D0, &
    0.9832D0, 0.9801D0, 0.9767D0, 0.9730D0, 0.9688D0, &
    0.9645D0, 0.9598D0, 0.9550D0, 0.9499D0, 0.9445D0, &
    0.9389D0, 0.9330D0, 0.9269D0, 0.9206D0, 0.9140D0, &
    0.9071D0, 0.9001D0, 0.8930D0, 0.8856D0, 0.8781D0, &
    0.8705D0, 0.8627D0, 0.8546D0, 0.8464D0, 0.8381D0, &
    0.8298D0, 0.8213D0, 0.8128D0, 0.8042D0, 0.7954D0, &
    0.7866D0, 0.7777D0, 0.7685D0, 0.7593D0, 0.7502D0, &
    0.7410D0, 0.7318D0, 0.7226D0, 0.7134D0, 0.7042D0, &
    0.6951D0, 0.6859D0, 0.6767D0, 0.6675D0, 0.6584D0, &
    0.6492D0, 0.6401D0, 0.6310D0, 0.6219D0, 0.6129D0, &
    0.6039D0, 0.5948D0, 0.5859D0, 0.5769D0, 0.5680D0, &
    0.5590D0, 0.5502D0, 0.5413D0, 0.5324D0, 0.5236D0, &
    0.5148D0, 0.5063D0, 0.4979D0, 0.4896D0, 0.4814D0, &
    0.4733D0, 0.4652D0, 0.4572D0, 0.4493D0, 0.4415D0, &
    0.4337D0, 0.4261D0, 0.4185D0, 0.4110D0, 0.4035D0, &
    0.3962D0, 0.3889D0, 0.3818D0, 0.3749D0, 0.3680D0, &
    0.3611D0, 0.3544D0, 0.3478D0, 0.3413D0, 0.3348D0, &
    0.3285D0, 0.3222D0, 0.3160D0, 0.3099D0, 0.3039D0, &
    0.2980D0, 0.2923D0, 0.2866D0, 0.2810D0, 0.2755D0, &
    0.2701D0, 0.2648D0, 0.2595D0, 0.2544D0, 0.2493D0, &
    0.2443D0, 0.2394D0, 0.2345D0, 0.2298D0, 0.2251D0, &
    0.2205D0, 0.2160D0, 0.2115D0, 0.2072D0, 0.2029D0, &
    0.1987D0 /)

  ! --- Local variables ---
  REAL(8)  :: ELOG      ! log10(nu / nu_edge)
  INTEGER :: I         ! table index

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COULBF1S'

  ! Below ionization threshold
  IF (FREQ / Z**2 .LT. FREQ_RYDH) THEN
    COULBF1S = 0.0D0
    RETURN
  END IF

  ! Linear interpolation in log10(nu/nu_edge)
  ELOG = log10(FREQ / Z**2 / FREQ_RYDH)
  I = int(ELOG / DLOG_STEP)
  I = max(min(I + 1, NGAUNT - 1), 1)

  COULBF1S = GAUNT1S(I) + (GAUNT1S(I+1) - GAUNT1S(I)) / DLOG_STEP &
           * (ELOG - dble(I - 1) * DLOG_STEP)

  RETURN

END FUNCTION COULBF1S

!=========================================================================
! FUNCTION COULFF(J, NZ)
!
! Thermally averaged free-free (bremsstrahlung) Gaunt factor for a
! hydrogen-like ion of effective charge Z.
!
! Returns <g_ff(gamma^2, u)> by bilinear interpolation of the table
! published by van Hoof et al. (2014), MNRAS 444, 420.  The table
! covers canonical log10(gamma^2) in [-6, +10] and log10(u) in
! [-16, +13] at 0.2-dex spacing (81 x 146 = 11,826 points), computed
! to ~1e-5 relative precision using arbitrary-precision arithmetic.
!
! Physical variables:
!   gamma^2 = Z^2 Ry / kT      = Z^2 * 1.57877e5 K / T
!   u       = h*nu / kT
!
! AXIS CONVENTION NOTE:
!   The original Kurucz COULFF stored the lookup variable as
!     GAMLOG = log10(158000 * Z^4 / T) = 2 * log10(gamma^2)  (canonical)
!     HVKTLG = log10((hnu/kT)^2)       = 2 * log10(u)        (canonical)
!   (see legacy comment "GAMLOG=LOG10(158000*Z*Z/T)*2" in the original
!   Fortran 77 code).  The old 11x12 table labels "log gamma^2 = -6..+4"
!   and "log u = -8..+3" were therefore twice the canonical values.
!   Once decoded, the old table agreed with van Hoof (2014) to ~1%.
!   The new code computes canonical log gamma^2 and log u directly.
!
! Data file:  $ATLAS12/gauntff_vanhoof.dat  (ASCII; original copyright
!   notice preserved.  See van Hoof et al. 2014; data available at
!   http://data.nublado.org/gauntff/.)
!
! Arguments:
!   J  - depth index (for accessing T, TLOG, FREQ, FREQLG via module arrays)
!   NZ - charge index (1-6) for Z^2 lookup
!
! Callers: HOP (NZ=1, neutral H), HE2OP (NZ=2, He+).
!=========================================================================

FUNCTION COULFF(J, NZ)

  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN) :: J       ! depth index
  INTEGER, INTENT(IN) :: NZ      ! charge index (1=H, 2=He+, ..., 6)

  ! --- Return value ---
  REAL(8) :: COULFF

  ! --- log10(Z^2) for NZ = 1-6 (canonical gamma^2 scales as Z^2) ---
  REAL(8), PARAMETER :: Z2LOG(6) = (/ &
    0.0D0, 0.60206D0, 0.95424D0, 1.20412D0, 1.39794D0, 1.55630D0 /)

  ! --- van Hoof (2014) table grid parameters ---
  INTEGER, PARAMETER :: GFF_N_GAM2       = 81
  INTEGER, PARAMETER :: GFF_N_U          = 146
  REAL(8),  PARAMETER :: GFF_LOG_GAM2_MIN = -6.0D0
  REAL(8),  PARAMETER :: GFF_LOG_U_MIN    = -16.0D0
  REAL(8),  PARAMETER :: GFF_DLOG         =  0.2D0
  INTEGER, PARAMETER :: GFF_MAGIC        = 20140210

  ! --- Cached Gaunt factor table (stored as log10(g_ff) for interpolation) ---
  ! First index: log(gamma^2); second index: log(u).
  REAL(8), SAVE :: GFF_LOGTAB(GFF_N_GAM2, GFF_N_U)
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! --- Local variables ---
  REAL(8)  :: LGAM2          ! canonical log10(gamma^2)
  REAL(8)  :: LU             ! canonical log10(u)
  REAL(8)  :: XI, XJ         ! fractional grid indices
  REAL(8)  :: FI, FJ         ! fractional parts
  REAL(8)  :: LG00, LG01, LG10, LG11  ! log(g_ff) at 4 grid corners
  REAL(8)  :: LG             ! interpolated log(g_ff)
  INTEGER :: I0, J0         ! lower-corner grid indices (1-based Fortran)

  ! --- Derived physical constants (computed from mod_constants).
  !     log10(Ry_H / k) = log10(1.57881e5 K), with Ry in energy units.
  !     log10(k/h)      = log10(2.08366e10 Hz/K).
  REAL(8), PARAMETER :: LOG10_RYH_OVER_K = log10(FREQ_RYDH * HOVERK)
  REAL(8), PARAMETER :: LOG10_K_OVER_H   = -log10(HOVERK)

  ! --- First-call initialization ---
  IF (.NOT. INITIALIZED) THEN
    CALL read_gauntff_table()
    INITIALIZED = .TRUE.
  END IF

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING COULFF'

  ! --- Compute canonical log gamma^2 = log10(Z^2 * Ry / (k*T))
  !       = log10(Z^2) + log10(Ry/k) - log10(T)
  LGAM2 = Z2LOG(NZ) + LOG10_RYH_OVER_K - TLOG(J) / LN10

  ! --- Compute canonical log u = log10(h*nu / (k*T))
  !       = log10(nu) - log10(T) - log10(k/h)
  LU = (FREQLG - TLOG(J)) / LN10 - LOG10_K_OVER_H

  ! --- Fractional grid indices (1-based).  Clamp to valid interior [1, N-1]
  !     so we always have a 4-cell neighborhood for bilinear interpolation.
  XI = (LGAM2 - GFF_LOG_GAM2_MIN) / GFF_DLOG + 1.0D0
  XJ = (LU    - GFF_LOG_U_MIN)    / GFF_DLOG + 1.0D0
  I0 = max(1, min(int(XI), GFF_N_GAM2 - 1))
  J0 = max(1, min(int(XJ), GFF_N_U    - 1))
  FI = XI - dble(I0)
  FJ = XJ - dble(J0)
  ! Clamp fractions too, so out-of-range lookups use the edge value rather
  ! than extrapolating.  Canonical stellar conditions sit well inside the
  ! grid; clamping is a safety net for pathological T/nu combinations.
  FI = max(0.0D0, min(FI, 1.0D0))
  FJ = max(0.0D0, min(FJ, 1.0D0))

  ! --- Bilinear interpolation of log10(g_ff), following van Hoof's own
  !     interpolator convention (keeps the interpolated quantity smooth
  !     over many orders of magnitude at large |log u|).
  LG00 = GFF_LOGTAB(I0,     J0)
  LG01 = GFF_LOGTAB(I0,     J0 + 1)
  LG10 = GFF_LOGTAB(I0 + 1, J0)
  LG11 = GFF_LOGTAB(I0 + 1, J0 + 1)
  LG   = (1.0D0 - FI) * ((1.0D0 - FJ) * LG00 + FJ * LG01) &
       +          FI  * ((1.0D0 - FJ) * LG10 + FJ * LG11)

  COULFF = 10.0D0 ** LG

  RETURN

CONTAINS

  !-------------------------------------------------------------------
  ! Read the van Hoof (2014) Gaunt factor table from
  ! $DATADIR/gauntff_vanhoof.dat.
  !
  ! The file is the original ASCII-formatted table distributed at
  ! http://data.nublado.org/gauntff/ (BSD-style license preserved).
  ! It contains a commented header followed by two 146-row data
  ! tables: first the Gaunt factors, then their uncertainty estimates.
  ! Each data row has 81 values, one per log10(gamma^2) grid point.
  !
  ! We store log10(g_ff) in GFF_LOGTAB(I, J) where I indexes
  ! log10(gamma^2) and J indexes log10(u), so (I, J) order matches
  ! Fortran column-major layout (fast axis = gamma^2).
  !-------------------------------------------------------------------
  SUBROUTINE read_gauntff_table()
    INTEGER :: IU, IOS, I, J
    INTEGER :: MAGIC_IN, N_GAM2_IN, N_U_IN
    REAL(8)  :: LOG_GAM2_IN, LOG_U_IN, DLOG_IN
    REAL(8)  :: ROW(GFF_N_GAM2)
    CHARACTER(256)  :: FILEPATH, LINE
    CHARACTER(2048) :: DATALINE   ! big enough for 81 %e values (~1300 chars)

    IU = 89
    FILEPATH = trim(DATADIR) // 'gauntff_vanhoof.dat'
    OPEN(UNIT=IU, FILE=FILEPATH, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'COULFF: ERROR opening ', trim(FILEPATH)
      STOP 'COULFF: cannot read van Hoof Gaunt factor table'
    END IF

    ! --- Skip copyright header (lines starting with '#') and consume the
    !     5 numeric metadata lines (magic, grid dims, start values, step).
    CALL read_nonblank_line(IU, LINE);  READ(LINE, *) MAGIC_IN
    CALL read_nonblank_line(IU, LINE);  READ(LINE, *) N_GAM2_IN, N_U_IN
    CALL read_nonblank_line(IU, LINE);  READ(LINE, *) LOG_GAM2_IN
    CALL read_nonblank_line(IU, LINE);  READ(LINE, *) LOG_U_IN
    CALL read_nonblank_line(IU, LINE);  READ(LINE, *) DLOG_IN

    ! --- Validate header against expected constants ---
    IF (MAGIC_IN  .NE. GFF_MAGIC .OR. &
        N_GAM2_IN .NE. GFF_N_GAM2 .OR. &
        N_U_IN    .NE. GFF_N_U .OR. &
        abs(LOG_GAM2_IN - GFF_LOG_GAM2_MIN) .GT. 1.0D-6 .OR. &
        abs(LOG_U_IN    - GFF_LOG_U_MIN)    .GT. 1.0D-6 .OR. &
        abs(DLOG_IN     - GFF_DLOG)         .GT. 1.0D-6) THEN
      WRITE(6,*) 'COULFF: ERROR file header mismatch for ', trim(FILEPATH)
      WRITE(6,*) '   magic:  expect ', GFF_MAGIC,        ' got ', MAGIC_IN
      WRITE(6,*) '   n_gam2: expect ', GFF_N_GAM2,       ' got ', N_GAM2_IN
      WRITE(6,*) '   n_u:    expect ', GFF_N_U,          ' got ', N_U_IN
      WRITE(6,*) '   log_gam2_min: expect ', GFF_LOG_GAM2_MIN, ' got ', LOG_GAM2_IN
      WRITE(6,*) '   log_u_min:    expect ', GFF_LOG_U_MIN,    ' got ', LOG_U_IN
      WRITE(6,*) '   dlog:         expect ', GFF_DLOG,         ' got ', DLOG_IN
      STOP 'COULFF: van Hoof table header does not match compiled constants'
    END IF

    ! --- Skip blank lines and '#' comment lines between the metadata and
    !     the first data row, then read the 146 rows of Gaunt factor data.
    !
    !     Each row corresponds to a fixed log(u) value (starting at
    !     log(u) = -16) and sweeps log(gamma^2) from -6 to +10 across its
    !     81 entries.  On-disk layout is [U, GAM2] row-major.  We
    !     transpose on the fly into [GAM2, U] so Fortran column-major
    !     indexing puts the hot lookup axis (gamma^2) first.  We also
    !     take log10 of each value here so runtime interpolation is
    !     purely additive.
    !
    !     Each data line has 81 floating-point values occupying ~1300
    !     characters -- we need a large line buffer.
    DO J = 1, GFF_N_U
      CALL read_nonblank_line(IU, DATALINE)
      READ(DATALINE, *, IOSTAT=IOS) (ROW(I), I = 1, GFF_N_GAM2)
      IF (IOS .NE. 0) THEN
        WRITE(6,*) 'COULFF: error parsing row ', J, ' of Gaunt factor data'
        STOP 'COULFF: corrupted Gaunt factor data file'
      END IF
      DO I = 1, GFF_N_GAM2
        GFF_LOGTAB(I, J) = log10(ROW(I))
      END DO
    END DO

    CLOSE(IU)
    ! --- Uncertainty table (second half of file) is not read; it is only
    !     useful for sanity-checking the published data itself, not for
    !     runtime opacity calculation.

  END SUBROUTINE read_gauntff_table

  !-------------------------------------------------------------------
  ! Read one non-blank, non-comment line from unit IU into LINE.
  ! Treats '#' as a line-initial comment marker (per the van Hoof
  ! file format).  Also strips inline '# comment' tails so the caller
  ! can read numeric values with list-directed READ.
  !-------------------------------------------------------------------
  SUBROUTINE read_nonblank_line(IU, LINE)
    INTEGER, INTENT(IN) :: IU
    CHARACTER(*), INTENT(OUT) :: LINE
    INTEGER :: IOS, K
    DO
      READ(IU, '(A)', IOSTAT=IOS) LINE
      IF (IOS .NE. 0) THEN
        WRITE(6,*) 'COULFF: unexpected EOF reading gauntff_vanhoof.dat'
        STOP 'COULFF: corrupted Gaunt factor data file'
      END IF
      LINE = adjustl(LINE)
      IF (len_trim(LINE) .EQ. 0) CYCLE
      IF (LINE(1:1) .EQ. '#')    CYCLE
      K = index(LINE, '#')
      IF (K .GT. 0) LINE(K:) = ' '
      RETURN
    END DO
  END SUBROUTINE read_nonblank_line

END FUNCTION COULFF

!=========================================================================
! SUBROUTINE H2PLOP
!
! H2+ molecular ion opacity.
!
! Computes the bound-free and free-free absorption by the H2+ molecular
! ion at the current frequency for all depth points. H2+ forms when a
! neutral H atom captures a proton, and is important in cool stellar
! atmospheres (T < ~8000 K) at wavelengths longward of Lyman limit.
!
! The cross-section is expressed as:
!   sigma(nu) = exp(-E_s/kT + F_R)
! where:
!   F_R = polynomial fit in log10(freq) for the frequency dependence
!   E_s = polynomial fit in freq for the temperature dependence
!         (activation/dissociation energy in eV)
!
! The opacity per gram is then:
!   kappa = sigma * n(H,1s) * 2 * b(1) * n(H+)/U(H+) / rho * (1 - e^(-hv/kT))
!
! where the factor 2*b(1) accounts for the departure coefficient and
! statistical weight, and the stimulated emission correction is applied.
!
! Only contributes for frequencies below the Lyman limit (13.6 eV).
!=========================================================================

SUBROUTINE H2PLOP

  IMPLICIT NONE

  REAL(8)  :: FR        ! log cross-section frequency factor (polynomial in FREQLG)
  REAL(8)  :: ES        ! energy factor (polynomial in FREQ), in eV
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING H2PLOP'

  ! No contribution above Lyman limit
  IF (FREQ .GT. FREQ_RYDH) RETURN

  ! Frequency-dependent cross-section factor (polynomial in log10(freq))
  FR = -3.0233D3 + (3.7797D2 + (-1.82496D1 &
     + (3.9207D-1 - 3.1672D-3 * FREQLG) * FREQLG) * FREQLG) * FREQLG

  ! Temperature-dependent energy factor (polynomial in freq), in eV
  ES = -7.342D-3 + (-2.409D-15 + (1.028D-30 &
     + (-4.230D-46 + (1.224D-61 - 1.351D-77 * FREQ) * FREQ) &
     * FREQ) * FREQ) * FREQ

  AH2P = exp(-ES / TKEV + FR) &
          * XNFP(:, 1) * 2.0D0 * BHYD(:, 1) * XNFP(:, 2) &
          / RHO * STIM

  RETURN

END SUBROUTINE H2PLOP

!=========================================================================
! SUBROUTINE HMINOP
!
! H⁻ bound-free and free-free opacity.
!
! Computes the absorption coefficient and source function for the
! negative hydrogen ion at the current frequency for all depth points.
! H⁻ is the dominant continuum opacity source in solar-type stars
! (Teff ~ 4000-8000 K).
!
! Two contributions:
!
! 1. Bound-free (photodetachment): H⁻ → H + e⁻
!    Cross-section from Mathisen (1984), based on Wishart (1979) and
!    Broad & Reinhardt (1976). Tabulated at 85 wavelength points from
!    18 to 1643.91 nm (the photodetachment threshold at 0.7542 eV).
!    Includes non-LTE departure coefficient BMIN for H⁻ population.
!    Only contributes for lambda < 1643.91 nm (freq > 1.82365e14 Hz).
!
! 2. Free-free (inverse bremsstrahlung): H + e⁻ → H + e⁻ + photon
!    From Bell & Berrington (1987, J.Phys.B 20, 801). Tabulated on a
!    22-wavelength × 11-theta grid (theta = 5040/T). Interpolated in
!    log(cross-section) vs log(wavelength) at each theta, then
!    interpolated in theta at each depth point.
!
! The source function SHMIN accounts for non-LTE effects in the
! bound-free component via departure coefficient BMIN, while
! free-free emits in LTE (S = B_nu).
!
! Temperature-dependent quantities (H⁻ population XHMIN, theta) are
! cached and only recomputed when ITEMP changes.
!
! -------------------------------------------------------------------------
! MODERNIZATION REVIEW (2026):
!
! Both data sources here were reviewed against the modern literature as
! part of the F90 modernization, with the conclusion that NO UPGRADE is
! warranted.  The reasoning, for posterity:
!
! BOUND-FREE: Wishart (1979) / Broad & Reinhardt (1976), as tabulated by
!   Mathisen (1984), remains the community standard reference.  It is the
!   cross-section used by the 2024 Barklem & Amarsi non-LTE H⁻ study
!   (A&A 689, A100), which explicitly adopted Wishart (1979) after reviewing
!   alternatives.  McLaughlin et al. (2017, ApJ 842, 65) provide new
!   R-matrix calculations that add Feshbach+shape resonances between 10.92
!   and 14.35 eV (88-113 nm), but from their Fig. 1:
!     "The differences compared to the Wishart (1979) data are very small
!      in the visual and UV (photon energies up to 7 eV, roughly redward
!      of 1700 Å)."
!   The current 85-point table already captures the first resonance at
!   ~113 nm (10.97 eV) with 10 points reaching a peak of 95e-18 cm².  The
!   additional resonances at 11.5-14.35 eV fall below the shortest non-
!   resonance point of the current grid (~100 nm) and are only relevant
!   for the far-UV spectra of hot stars (Teff >~ 20,000 K).
!
!   John (1988, A&A 193, 189) is NOT a higher-accuracy calculation: it is
!   an analytic fit formula to the Wishart data.  Adopting John (1988)
!   would be a convenience-for-accuracy trade that cannot represent the
!   resonance structure at all, and is not an improvement.
!
! FREE-FREE: Bell & Berrington (1987) is an R-matrix calculation using 1s,
!   2s, 2p hydrogen states plus three pseudostates, quoted in the paper as
!   "currently the most accurate available" and still so as of this writing.
!   No modern replacement exists; Barklem & Amarsi (2024) also use this
!   source.  Again, John (1988) gives an analytic fit to these same
!   Bell & Berrington values -- a convenience formulation, not a new
!   calculation.
!
! NET RESULT: Both the bound-free and free-free cross-sections currently in
!   use here are the same data that modern 2024-era non-LTE studies adopt.
!   Any replacement would either be (a) the same data in a different form
!   (John 1988 fit formulas), which is lower accuracy, or (b) McLaughlin
!   (2017), which only changes things below 108 nm -- outside the regime
!   where ATLAS is typically used and which the current table already
!   covers adequately via Mathisen's tabulation of the first resonance.
!
! The original Kurucz code was well-chosen; no change needed.
! -------------------------------------------------------------------------
!
! References:
!   Wishart, A.W. 1979, MNRAS 187, 59P           [bound-free cross-section]
!   Broad, J.T. & Reinhardt, W.P. 1976, Phys.Rev.A 14, 2159
!                                                [resonance region of bf]
!   Mathisen, R. 1984, Inst.Theor.Astrophys.Oslo Pub.Series No. 1
!                                                [85-point tabulation used]
!   Bell, K.L. & Berrington, K.A. 1987, J.Phys.B 20, 801  [free-free]
!   Hotop, H. & Lineberger, W.C. 1985, J.Phys.Chem.Ref.Data 14, 731
!                                  [H⁻ electron affinity = 0.754209 eV]
!
!   McLaughlin, B.M. et al. 2017, ApJ 842, 65    [modern bf; for far-UV
!                                                 resonances only, not
!                                                 adopted here]
!   Barklem, P.S. & Amarsi, A.M. 2024, A&A 689, A100
!                                  [modern non-LTE H⁻ study using the
!                                   same Wishart/Bell-Berrington data]
!=========================================================================

SUBROUTINE HMINOP

  IMPLICIT NONE

  ! --- Bound-free cross-section table ---
  ! 85 points: wavelength (nm) and sigma (1e-18 cm^2)
  ! From Mathisen (1984) after Wishart (1979), Broad & Reinhardt (1976)
  INTEGER, PARAMETER :: NBF = 85

  REAL(8), PARAMETER :: WBF(85) = (/ &
    18.00D0, 19.60D0, 21.40D0, 23.60D0, 26.40D0, 29.80D0, 34.30D0, &
    40.40D0, 49.10D0, 62.60D0, 111.30D0, 112.10D0, 112.67D0, &
    112.95D0, 113.05D0, 113.10D0, 113.20D0, 113.23D0, 113.50D0, &
    114.40D0, 121.00D0, 139.00D0, 164.00D0, 175.00D0, 200.00D0, &
    225.00D0, 250.00D0, 275.00D0, 300.00D0, 325.00D0, 350.00D0, &
    375.00D0, 400.00D0, 425.00D0, 450.00D0, 475.00D0, 500.00D0, &
    525.00D0, 550.00D0, 575.00D0, 600.00D0, 625.00D0, 650.00D0, &
    675.00D0, 700.00D0, 725.00D0, 750.00D0, 775.00D0, 800.00D0, &
    825.00D0, 850.00D0, 875.00D0, 900.00D0, 925.00D0, 950.00D0, &
    975.00D0, 1000.00D0, 1025.00D0, 1050.00D0, 1075.00D0, &
    1100.00D0, 1125.00D0, 1150.00D0, 1175.00D0, 1200.00D0, &
    1225.00D0, 1250.00D0, 1275.00D0, 1300.00D0, 1325.00D0, &
    1350.00D0, 1375.00D0, 1400.00D0, 1425.00D0, 1450.00D0, &
    1475.00D0, 1500.00D0, 1525.00D0, 1550.00D0, 1575.00D0, &
    1600.00D0, 1610.00D0, 1620.00D0, 1630.00D0, 1643.91D0 /)

  REAL(8), PARAMETER :: BF(85) = (/ &
    0.067D0, 0.088D0, 0.117D0, 0.155D0, 0.206D0, 0.283D0, 0.414D0, &
    0.703D0, 1.24D0, 2.33D0, 11.60D0, 13.90D0, 24.30D0, 66.70D0, &
    95.00D0, 56.60D0, 20.00D0, 14.60D0, 8.50D0, 7.10D0, 5.43D0, &
    5.91D0, 7.29D0, 7.918D0, 9.453D0, 11.08D0, 12.75D0, 14.46D0, &
    16.19D0, 17.92D0, 19.65D0, 21.35D0, 23.02D0, 24.65D0, 26.24D0, &
    27.77D0, 29.23D0, 30.62D0, 31.94D0, 33.17D0, 34.32D0, 35.37D0, &
    36.32D0, 37.17D0, 37.91D0, 38.54D0, 39.07D0, 39.48D0, 39.77D0, &
    39.95D0, 40.01D0, 39.95D0, 39.77D0, 39.48D0, 39.06D0, 38.53D0, &
    37.89D0, 37.13D0, 36.25D0, 35.28D0, 34.19D0, 33.01D0, 31.72D0, &
    30.34D0, 28.87D0, 27.33D0, 25.71D0, 24.02D0, 22.26D0, 20.46D0, &
    18.62D0, 16.74D0, 14.85D0, 12.95D0, 11.07D0, 9.211D0, 7.407D0, &
    5.677D0, 4.052D0, 2.575D0, 1.302D0, 0.8697D0, 0.4974D0, &
    0.1989D0, 0.0D0 /)

  ! --- Free-free cross-section table ---
  ! Bell & Berrington (1987): FF(theta_index, wave_index)
  ! 11 theta values × 22 wavelength values (in units of 1e-26 cm^2)
  INTEGER, PARAMETER :: NFF_THETA = 11, NFF_WAVE = 22

  REAL(8), PARAMETER :: THETAFF(11) = (/ &
    0.5D0, 0.6D0, 0.8D0, 1.0D0, 1.2D0, 1.4D0, &
    1.6D0, 1.8D0, 2.0D0, 2.8D0, 3.6D0 /)

  REAL(8), PARAMETER :: WAVEK(22) = (/ &
    0.50D0, 0.40D0, 0.35D0, 0.30D0, 0.25D0, 0.20D0, 0.18D0, &
    0.16D0, 0.14D0, 0.12D0, 0.10D0, 0.09D0, 0.08D0, 0.07D0, &
    0.06D0, 0.05D0, 0.04D0, 0.03D0, 0.02D0, 0.01D0, 0.008D0, &
    0.006D0 /)

  REAL(8), PARAMETER :: FF(11, 22) = reshape( (/ &
    .0178, .0222, .0308, .0402, .0498, .0596, .0695, .0795, .0896, .131,  .172,  &
    .0228, .0280, .0388, .0499, .0614, .0732, .0851, .0972, .110,  .160,  .211,  &
    .0277, .0342, .0476, .0615, .0760, .0908, .105,  .121,  .136,  .199,  .262,  &
    .0364, .0447, .0616, .0789, .0966, .114,  .132,  .150,  .169,  .243,  .318,  &
    .0520, .0633, .0859, .108,  .131,  .154,  .178,  .201,  .225,  .321,  .418,  &
    .0791, .0959, .129,  .161,  .194,  .227,  .260,  .293,  .327,  .463,  .602,  &
    .0965, .117,  .157,  .195,  .234,  .272,  .311,  .351,  .390,  .549,  .711,  &
    .121,  .146,  .195,  .241,  .288,  .334,  .381,  .428,  .475,  .667,  .861,  &
    .154,  .188,  .249,  .309,  .367,  .424,  .482,  .539,  .597,  .830,  1.07,  &
    .208,  .250,  .332,  .409,  .484,  .557,  .630,  .702,  .774,  1.06,  1.36,  &
    .293,  .354,  .468,  .576,  .677,  .777,  .874,  .969,  1.06,  1.45,  1.83,  &
    .358,  .432,  .572,  .702,  .825,  .943,  1.06,  1.17,  1.28,  1.73,  2.17,  &
    .448,  .539,  .711,  .871,  1.02,  1.16,  1.29,  1.43,  1.57,  2.09,  2.60,  &
    .579,  .699,  .924,  1.13,  1.33,  1.51,  1.69,  1.86,  2.02,  2.67,  3.31,  &
    .781,  .940,  1.24,  1.52,  1.78,  2.02,  2.26,  2.48,  2.69,  3.52,  4.31,  &
    1.11,  1.34,  1.77,  2.17,  2.53,  2.87,  3.20,  3.51,  3.80,  4.92,  5.97,  &
    1.73,  2.08,  2.74,  3.37,  3.90,  4.50,  5.01,  5.50,  5.95,  7.59,  9.06,  &
    3.04,  3.65,  4.80,  5.86,  6.86,  7.79,  8.67,  9.50,  10.3,  13.2,  15.6,  &
    6.79,  8.16,  10.7,  13.1,  15.3,  17.4,  19.4,  21.2,  23.0,  29.5,  35.0,  &
    27.0,  32.4,  42.6,  51.9,  60.7,  68.9,  76.8,  84.2,  91.4,  117.,  140.,  &
    42.3,  50.6,  66.4,  80.8,  94.5,  107.,  120.,  131.,  142.,  183.,  219.,  &
    75.1,  90.0,  118.,  144.,  168.,  191.,  212.,  234.,  253.,  325.,  388.   &
    /), (/ NFF_THETA, NFF_WAVE /) )

  ! --- Cached arrays (recomputed when ITEMP changes) ---
  REAL(8),  SAVE :: XHMIN(kw)       ! H⁻ number density / rho factor
  REAL(8),  SAVE :: THETA(kw)       ! 5040/T at each depth
  INTEGER, SAVE :: ITEMP1 = 0      ! cached ITEMP for skip logic

  ! --- One-time initialization arrays ---
  REAL(8),  SAVE :: WFFLOG(NFF_WAVE)              ! log(wavelength) grid for ff
  REAL(8),  SAVE :: FFLOG(NFF_WAVE, NFF_THETA)    ! log(cross-section) table
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! --- Local variables ---
  REAL(8)  :: FFTT(NFF_THETA)  ! ff cross-section interpolated to current wavelength
  REAL(8)  :: FFTHETA(1)       ! ff cross-section interpolated to current depth's theta
  REAL(8)  :: WAVELOG(1)       ! log(wavelength) at current frequency
  REAL(8)  :: HMINBF_ARR(1)     ! bound-free cross-section (1e-18 cm^2)
  REAL(8)  :: WAVE_ARR(1)
  REAL(8)  :: HMINFF           ! free-free opacity per gram
  REAL(8)  :: H                ! bound-free opacity per gram
  REAL(8)  :: FFTLOG(1)        ! temporary for LINTER output
  INTEGER :: J, IWAVE, ITHETA, MAXWAVE

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HMINOP'

  !---------------------------------------------------------------------
  ! One-time initialization: precompute log(wavelength) and log(sigma)
  ! grids for free-free interpolation
  !---------------------------------------------------------------------
  IF (.NOT. INITIALIZED) THEN
    DO IWAVE = 1, NFF_WAVE
      ! 91.134 nm is the conversion factor from Bell & Berrington
      WFFLOG(IWAVE) = log(91.134D0 / WAVEK(IWAVE))
      DO ITHETA = 1, NFF_THETA
        FFLOG(IWAVE, ITHETA) = log(FF(ITHETA, IWAVE) * 1.0D-26)
      END DO
    END DO
    INITIALIZED = .TRUE.
  END IF

  !---------------------------------------------------------------------
  ! Recompute temperature-dependent quantities if T has changed
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    THETA = THETA_COEFF / T
    ! H⁻ population: Saha equation for H⁻ ↔ H + e⁻
    ! 0.754209 eV = H⁻ electron affinity (Hotop & Lineberger 1985)
    XHMIN = exp(HMINUS_EA / TKEV) &
             / (SAHA_PREFAC * T * sqrt(T)) &
             * BMIN * BHYD(:, 1) * XNFP(:, 1) * XNE
  END IF

  !---------------------------------------------------------------------
  ! Free-free: interpolate table to current wavelength, then to each
  ! depth's theta value
  !---------------------------------------------------------------------
  WAVELOG(1) = log(WAVE)
  DO ITHETA = 1, NFF_THETA
    CALL LINTER(WFFLOG, FFLOG(1, ITHETA), NFF_WAVE, WAVELOG(1), FFTLOG(1))
    FFTT(ITHETA) = exp(FFTLOG(1)) / THETAFF(ITHETA) * THETA_COEFF * KBOL
  END DO

  !---------------------------------------------------------------------
  ! Bound-free: interpolate cross-section table at current wavelength
  ! (only contributes for lambda < 1643.91 nm, i.e. freq > 1.82365e14)
  !---------------------------------------------------------------------
  HMINBF_ARR(1) = 0.0D0
  IF (FREQ .GT. 1.82365D14) THEN
    WAVE_ARR(1) = WAVE
    MAXWAVE = MAP1(WBF, BF, NBF, WAVE_ARR, HMINBF_ARR, 1)
  END IF

  !---------------------------------------------------------------------
  ! Assemble total opacity and source function at each depth
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    CALL LINTER(THETAFF, FFTT, NFF_THETA, THETA(J), FFTHETA(1))

    ! Free-free opacity per gram
    HMINFF = FFTHETA(1) * XNFP(J, 1) * 2.0D0 * BHYD(J, 1) * XNE(J) / RHO(J)

    ! Bound-free opacity per gram (with non-LTE and stimulated emission)
    H = HMINBF_ARR(1) * 1.0D-18 * (1.0D0 - EHVKT(J) / BMIN(J)) &
      * XHMIN(J) / RHO(J)

    AHMIN(J) = H + HMINFF

    ! Source function: bf weighted by non-LTE, ff in LTE
    SHMIN(J) = (H * BNU(J) * STIM(J) / (BMIN(J) - EHVKT(J)) &
             + HMINFF * BNU(J)) / AHMIN(J)
  END DO

  RETURN

END SUBROUTINE HMINOP

!=========================================================================
! SUBROUTINE LINTER
!
! Single-point linear interpolation in a monotonically increasing table.
!
! Given a table (XOLD, YOLD) of NOLD points with XOLD increasing,
! interpolates to find YNEW at the single point XNEW. Extrapolates
! using the nearest interval if XNEW is outside the table range.
!=========================================================================

SUBROUTINE LINTER(XOLD, YOLD, NOLD, XNEW, YNEW)

  IMPLICIT NONE

  ! --- Arguments ---
  INTEGER, INTENT(IN)  :: NOLD
  REAL(8),  INTENT(IN)  :: XOLD(NOLD), YOLD(NOLD)
  REAL(8),  INTENT(IN)  :: XNEW
  REAL(8),  INTENT(OUT) :: YNEW

  ! --- Local variables ---
  INTEGER :: I

  ! Find the bracketing interval
  DO I = 2, NOLD
    IF (XNEW .LT. XOLD(I)) EXIT
  END DO
  IF (I .GT. NOLD) I = NOLD

  ! Linear interpolation (or extrapolation at boundaries)
  YNEW = YOLD(I-1) + (YOLD(I) - YOLD(I-1)) &
       / (XOLD(I) - XOLD(I-1)) * (XNEW - XOLD(I-1))

  RETURN

END SUBROUTINE LINTER

!=========================================================================
! SUBROUTINE HRAYOP
!
! Hydrogen Rayleigh scattering opacity.
!
! Full Gavrila (1966) treatment matching atlas7lib.for.
! Reference: Gavrila, M. 1966, "Coherent scattering of light by
!   atomic hydrogen", JILA Report No. 86.
!
! The scattering amplitude f(ν) is tabulated across five frequency
! ranges, capturing the resonance structure near each Lyman line:
!
!   Range 1: ν < 0.01 ν_L       — extrapolation using f(1)²×(ν/ν_L)⁴
!   Range 2: ν ≤ 0.74 ν_L       — GAVRILAM (74 pts, step 0.01 ν_L)
!   Range 3: 0.77–0.885 ν_L     — GAVRILAMAB (27 pts, near Ly β)
!   Range 4: 0.890–0.936 ν_L    — GAVRILAMBC (24 pts, near Ly γ)
!   Range 5: 0.938–0.959 ν_L    — GAVRILAMCD (22 pts, near Ly δ)
!   Range 6: 0.961–1.000 ν_L    — constant (series limit)
!   Range 7: ν > ν_L            — GAVRILALYMANCONT (64 pts, continuum)
!
! Cross-section: σ = 6.65×10⁻²⁵ × f² (in ranges 1–5)
!            or: σ = 6.65×10⁻²⁵ × f   (in range 7, above Lyman limit)
!
! Gaps between ranges (0.74–0.77, 0.885–0.890, 0.936–0.938,
! 0.959–0.961 ν_L) have XSECT = 0.
!
! The opacity per gram:
!   κ_scat = σ × n(H,1s) × 2 × b(1) / ρ
!=========================================================================

SUBROUTINE HRAYOP

  IMPLICIT NONE

  ! FREQ_RYDH and SIGMA_THOMSON from mod_constants

  ! Range 1: scattering amplitude f, ν/ν_L = 0.01 to 0.74 by 0.01
  REAL(8), PARAMETER :: GAVRILAM(74) = (/ &
   -0.000113D0,  -0.000450D0,  -0.001014D0,  -0.001804D0,  -0.002823D0, &
   -0.004072D0,  -0.005553D0,  -0.007269D0,  -0.009223D0,  -0.011419D0, &
   -0.013861D0,  -0.016553D0,  -0.019500D0,  -0.022709D0,  -0.026185D0, &
   -0.029936D0,  -0.033968D0,  -0.038291D0,  -0.042913D0,  -0.047843D0, &
   -0.053093D0,  -0.058674D0,  -0.064599D0,  -0.070882D0,  -0.077537D0, &
   -0.084581D0,  -0.092031D0,  -0.099907D0,  -0.108230D0,  -0.117022D0, &
   -0.126308D0,  -0.136117D0,  -0.146477D0,  -0.157422D0,  -0.168987D0, &
   -0.181213D0,  -0.194143D0,  -0.207825D0,  -0.222313D0,  -0.237667D0, &
   -0.253953D0,  -0.271245D0,  -0.289626D0,  -0.309189D0,  -0.330041D0, &
   -0.352300D0,  -0.376103D0,  -0.401605D0,  -0.428985D0,  -0.458448D0, &
   -0.490235D0,  -0.524625D0,  -0.561947D0,  -0.602591D0,  -0.647023D0, &
   -0.695805D0,  -0.749619D0,  -0.809306D0,  -0.875910D0,  -0.950750D0, &
   -1.035515D0,  -1.132403D0,  -1.244337D0,  -1.375285D0,  -1.530787D0, &
   -1.718821D0,  -1.951320D0,  -2.246993D0,  -2.636960D0,  -3.177142D0, &
   -3.979234D0,  -5.303624D0,  -7.930999D0, -15.763602D0 /)

  ! Range 3: near Ly β, ν/ν_L = 0.755 to 0.885 by 0.005
  REAL(8), PARAMETER :: GAVRILAMAB(27) = (/ &
   31.008832D0,  15.382871D0,  10.160646D0,   7.538338D0,   5.955062D0, &
    4.890397D0,   4.121176D0,   3.535672D0,   3.071659D0,   2.691623D0, &
    2.371483D0,   2.094936D0,   1.850395D0,   1.629203D0,   1.424526D0, &
    1.230596D0,   1.042127D0,   0.853766D0,   0.659460D0,   0.451533D0, &
    0.219115D0,  -0.054939D0,  -0.400868D0,  -0.879559D0,  -1.637857D0, &
   -3.150374D0,  -8.326078D0 /)

  ! Range 4: near Ly γ, ν/ν_L = 0.890 to 0.936 by 0.002
  REAL(8), PARAMETER :: GAVRILAMBC(24) = (/ &
   32.260389D0,  11.880702D0,   7.418436D0,   5.442077D0,   4.313409D0, &
    3.573504D0,   3.043218D0,   2.637983D0,   2.312466D0,   2.039959D0, &
    1.803441D0,   1.591244D0,   1.394717D0,   1.206823D0,   1.021148D0, &
    0.831020D0,   0.628449D0,   0.402484D0,   0.136127D0,  -0.200462D0, &
   -0.667435D0,  -1.410661D0,  -2.906862D0,  -8.169314D0 /)

  ! Range 5: near Ly δ, ν/ν_L = 0.938 to 0.959 by 0.001
  REAL(8), PARAMETER :: GAVRILAMCD(22) = (/ &
   27.981406D0,   9.816495D0,   6.145775D0,   4.544224D0,   3.630968D0, &
    3.029081D0,   2.593248D0,   2.255265D0,   1.978565D0,   1.741426D0, &
    1.529699D0,   1.333240D0,   1.143898D0,   0.954154D0,   0.755875D0, &
    0.538760D0,   0.287687D0,  -0.022759D0,  -0.441666D0,  -1.081712D0, &
   -2.278530D0,  -5.705843D0 /)

  ! Range 7: above Lyman limit, correction factors
  REAL(8), PARAMETER :: GAVRILALYMANCONT(64) = (/ &
    2.667783D0,   2.526696D0,   2.408970D0,   2.308970D0,   2.222736D0, &
    2.147415D0,   2.080913D0,   2.021653D0,   1.968431D0,   1.920304D0, &
    1.876527D0,   1.799739D0,   1.734455D0,   1.678180D0,   1.629118D0, &
    1.585943D0,   1.547643D0,   1.513435D0,   1.482700D0,   1.454941D0, &
    1.429751D0,   1.406798D0,   1.385804D0,   1.366536D0,   1.348797D0, &
    1.332419D0,   1.317257D0,   1.303187D0,   1.290100D0,   1.277901D0, &
    1.266509D0,   1.255848D0,   1.245856D0,   1.236474D0,   1.227652D0, &
    1.219344D0,   1.190492D0,   1.167227D0,   1.148153D0,   1.132293D0, &
    1.118945D0,   1.107593D0,   1.097848D0,   1.089413D0,   1.082059D0, &
    1.075606D0,   1.069908D0,   1.064850D0,   1.060338D0,   1.056294D0, &
    1.052655D0,   1.038936D0,   1.030042D0,   1.023928D0,   1.019536D0, &
    1.016269D0,   1.011814D0,   1.008986D0,   1.007074D0,   1.005720D0, &
    1.004724D0,   1.003970D0,   1.003385D0,   1.003140D0 /)

  ! Frequency grid for Lyman continuum table (ν/ν_L)
  REAL(8), PARAMETER :: FGAVRILALYMANCONT(64) = (/ &
    1.00D0, 1.05D0, 1.10D0, 1.15D0, 1.20D0, 1.25D0, 1.30D0, 1.35D0, &
    1.40D0, 1.45D0, 1.5D0,  1.6D0,  1.7D0,  1.8D0,  1.9D0,  2.0D0,  &
    2.1D0,  2.2D0,  2.3D0,  2.4D0,  2.5D0,  2.6D0,  2.7D0,  2.8D0,  &
    2.9D0,  3.0D0,  3.1D0,  3.2D0,  3.3D0,  3.4D0,  3.5D0,  3.6D0,  &
    3.7D0,  3.8D0,  3.9D0,  4.0D0,  4.4D0,  4.8D0,  5.2D0,  5.6D0,  &
    6.0D0,  6.4D0,  6.8D0,  7.2D0,  7.6D0,  8.0D0,  8.4D0,  8.8D0,  &
    9.2D0,  9.6D0, 10.0D0, 12.0D0, 14.0D0, 16.0D0, 18.0D0, 20.0D0, &
   24.0D0, 28.0D0, 32.0D0, 36.0D0, 40.0D0, 44.0D0, 48.0D0, 50.0D0 /)

  ! Local variables
  REAL(8)  :: XSECT, G
  INTEGER :: I, J, IDUM

  ! Frequency step sizes for each range
  REAL(8), PARAMETER :: DFREQ1  = 0.01D0  * FREQ_RYDH  ! 0.01 × ν_L
  REAL(8), PARAMETER :: DFREQAB = 0.005D0 * FREQ_RYDH  ! 0.005 × ν_L
  REAL(8), PARAMETER :: DFREQBC = 0.002D0 * FREQ_RYDH  ! 0.002 × ν_L
  REAL(8), PARAMETER :: DFREQCD = 0.001D0 * FREQ_RYDH  ! 0.001 × ν_L

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HRAYOP'

  XSECT = 0.0D0

  IF (FREQ .LT. DFREQ1) THEN
    ! Below 0.01 ν_L: extrapolate using lowest table value
    XSECT = SIGMA_THOMSON * GAVRILAM(1)**2 * (FREQ / DFREQ1)**4

  ELSE IF (FREQ .LE. 0.74D0 * FREQ_RYDH) THEN
    ! Range 1: 0.01–0.74 ν_L, table step = 0.01 ν_L
    I = int(FREQ / DFREQ1)
    I = MIN(I + 1, 74)
    G = GAVRILAM(I-1) + (GAVRILAM(I) - GAVRILAM(I-1)) / DFREQ1 &
      * (FREQ - dble(I-1) * DFREQ1)
    XSECT = SIGMA_THOMSON * G**2

  ELSE IF (FREQ .LT. 0.77D0 * FREQ_RYDH) THEN
    ! Gap between range 1 and Ly β region
    XSECT = 0.0D0

  ELSE IF (FREQ .LE. 0.885D0 * FREQ_RYDH) THEN
    ! Range 3: 0.755–0.885 ν_L (near Ly β), step = 0.005 ν_L
    I = int((FREQ - 0.755D0 * FREQ_RYDH) / DFREQAB)
    I = I + 1
    I = MIN(I + 1, 27)
    G = GAVRILAMAB(I-1) + (GAVRILAMAB(I) - GAVRILAMAB(I-1)) / DFREQAB &
      * (FREQ - (0.755D0 * FREQ_RYDH + dble(I-1-1) * DFREQAB))
    XSECT = SIGMA_THOMSON * G**2

  ELSE IF (FREQ .LT. 0.890D0 * FREQ_RYDH) THEN
    ! Gap between Ly β and Ly γ regions
    XSECT = 0.0D0

  ELSE IF (FREQ .LE. 0.936D0 * FREQ_RYDH) THEN
    ! Range 4: 0.890–0.936 ν_L (near Ly γ), step = 0.002 ν_L
    I = int((FREQ - 0.890D0 * FREQ_RYDH) / DFREQBC)
    I = I + 1
    I = MIN(I + 1, 24)
    G = GAVRILAMBC(I-1) + (GAVRILAMBC(I) - GAVRILAMBC(I-1)) / DFREQBC &
      * (FREQ - (0.890D0 * FREQ_RYDH + dble(I-1-1) * DFREQBC))
    XSECT = SIGMA_THOMSON * G**2

  ELSE IF (FREQ .LT. 0.938D0 * FREQ_RYDH) THEN
    ! Gap between Ly γ and Ly δ regions
    XSECT = 0.0D0

  ELSE IF (FREQ .LE. 0.959D0 * FREQ_RYDH) THEN
    ! Range 5: 0.938–0.959 ν_L (near Ly δ), step = 0.001 ν_L
    I = int((FREQ - 0.938D0 * FREQ_RYDH) / DFREQCD)
    I = I + 1
    I = MIN(I + 1, 22)
    G = GAVRILAMCD(I-1) + (GAVRILAMCD(I) - GAVRILAMCD(I-1)) / DFREQCD &
      * (FREQ - (0.938D0 * FREQ_RYDH + dble(I-1-1) * DFREQCD))
    XSECT = SIGMA_THOMSON * G**2

  ELSE IF (FREQ .LT. 0.961D0 * FREQ_RYDH) THEN
    ! Gap between Ly δ and series limit
    XSECT = 0.0D0

  ELSE IF (FREQ .LE. FREQ_RYDH) THEN
    ! Near series limit: constant value
    XSECT = SIGMA_THOMSON * GAVRILALYMANCONT(1)

  ELSE
    ! Above Lyman limit: interpolate from continuum table
    BLOCK
      REAL(8) :: XNEW(1), FNEW(1)
      XNEW(1) = FREQ / FREQ_RYDH
      IDUM = MAP1(FGAVRILALYMANCONT, GAVRILALYMANCONT, 64, &
                  XNEW, FNEW, 1)
      XSECT = SIGMA_THOMSON * FNEW(1)
    END BLOCK

  END IF

  ! Apply to all depth points
  SIGH = XSECT * XNFP(:, 1) * 2.0D0 * BHYD(:, 1) / RHO

  RETURN

END SUBROUTINE HRAYOP

!=========================================================================
! SUBROUTINE HE1OP
!
! Neutral helium (He I) bound-free and free-free opacity.
!
! Computes the He I continuum absorption at the current frequency for
! all depth points. The opacity includes:
!
! 1. Bound-free from 10 resolved low-lying levels:
!      1s^2 1S (ground), 1s2s 3S, 1s2s 1S, 1s2p 3P, 1s2p 1P,
!      1s3s 3S, 1s3s 1S, 1s3p 3P, 1s3d 3D+1D, 1s3p 1P
!    Ground state uses the CROSSHE fit; 2s/2p states use dedicated
!    fits (HE12s3S, HE12s1S, HE12p3P, HE12p1P); 3-shell states
!    use hydrogenic XKARZAS with fitted effective charges.
!
! 2. Inner-shell ionization: when the photon can also eject the inner
!    1s electron from excited He I, leaving He II in n=2 or n=3.
!    These "two-electron" transitions add to the cross-sections of
!    the low-lying levels at sufficiently high frequencies.
!
! 3. High-n bound-free (n=4-27): hydrogenic with Z_eff^2 = 4 - 3/n^2,
!    only for freq > 1.25408e16 Hz (He I series limit from n=4).
!
! 4. Dissolved levels near the series limit (BOLTEX vs EXLIM).
!
! 5. Free-free: He II + e⁻ bremsstrahlung via COULFF.
!
! Temperature-dependent quantities are cached and only recomputed
! when ITEMP changes.
!=========================================================================

SUBROUTINE HE1OP

  IMPLICIT NONE

  ! --- Energy level data for 10 resolved low-lying states ---
  REAL(8), PARAMETER :: CHI(10) = (/ &
    0.0D0, 19.819D0, 20.615D0, 20.964D0, 21.217D0, &
    22.718D0, 22.920D0, 23.006D0, 23.073D0, 23.086D0 /)

  REAL(8), PARAMETER :: HEFREQ(10) = (/ &
    5.945209D15, 1.152844D15, 0.9603331D15, 0.8761076D15, &
    0.8147104D15, 0.4519048D15, 0.4030971D15, 0.3821191D15, &
    0.3660215D15, 0.3627891D15 /)

  REAL(8), PARAMETER :: G(10) = (/ &
    1.0D0, 3.0D0, 1.0D0, 9.0D0, 3.0D0, &
    3.0D0, 1.0D0, 9.0D0, 20.0D0, 3.0D0 /)

  ! --- Cached temperature-dependent arrays ---
  REAL(8),  SAVE :: BOLT(kw, 10)     ! Boltzmann factors for 10 levels
  REAL(8),  SAVE :: BOLTN(kw, 27)    ! Boltzmann factors for high-n levels
  REAL(8),  SAVE :: EXLIM(kw)        ! population at ionization limit (24.587 eV)
  REAL(8),  SAVE :: BOLTEX(kw)       ! population at dissolved limit (23.730 eV)
  REAL(8),  SAVE :: FREET(kw)        ! free-free factor: n_e * n(He+) / (rho * sqrt(T))
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Local variables ---
  REAL(8)  :: TRANS(10)     ! bound-free cross-sections for 10 levels
  REAL(8)  :: TRANSN(27)    ! bound-free cross-sections for high-n levels
  REAL(8)  :: FREQ3, CFREE, C
  REAL(8)  :: RYD, ELIM, FREQHE, ZEFF2
  REAL(8)  :: XR, EX, HE1
  INTEGER :: J, N, IMIN

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE1OP'

  !---------------------------------------------------------------------
  ! Recompute temperature-dependent quantities if T has changed
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    RYD = RYDBERG_HE * CLIGHT
    DO J = 1, NRHOX
      ! Boltzmann populations for 10 resolved levels
      DO N = 1, 10
        BOLT(J, N) = exp(-CHI(N) / TKEV(J)) * G(N) * XNFP(J, 3) / RHO(J)
      END DO
      ! High-n levels (n=4-27): hydrogenic with E_n = 24.587*(1-1/n^2) eV
      DO N = 4, 27
        BOLTN(J, N) = exp(-24.587D0 * (1.0D0 - 1.0D0 / dble(N)**2) / TKEV(J)) &
                    * 4.0D0 * dble(N)**2 * XNFP(J, 3) / RHO(J)
      END DO
      FREET(J) = XNE(J) * XNF(J, 4) / RHO(J) / sqrt(T(J))
      XR = XNFP(J, 3) * (4.0D0 / 2.0D0 / 13.595D0) * TKEV(J) / RHO(J)
      BOLTEX(J) = exp(-23.730D0 / TKEV(J)) * XR
      EXLIM(J) = exp(-24.587D0 / TKEV(J)) * XR
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Cross-sections at the current frequency
  !---------------------------------------------------------------------
  FREQ3 = FREQ**3
  CFREE = COEFF_FF / FREQ3
  C = 2.815D29 / FREQ3
  RYD = RYDBERG_HE * CLIGHT

  ! Find lowest level whose threshold is at or below FREQ
  IMIN = 0
  DO N = 1, 10
    IF (HEFREQ(N) .LE. FREQ) THEN
      IMIN = N
      EXIT
    END IF
  END DO

  ! Compute bound-free cross-sections for levels IMIN..10
  ! (HEFREQ is in decreasing order: ground state has the highest threshold.
  ! IMIN is the first level whose threshold <= FREQ. The computed GOTO in
  ! the original fell through from label 20+IMIN to 30, so all levels
  ! N = IMIN..10 are computed.)
  TRANS = 0.0D0
  IF (IMIN .GT. 0 .AND. IMIN .LE. 1)  TRANS(1)  = CROSSHE(FREQ)
  IF (IMIN .GT. 0 .AND. IMIN .LE. 2)  TRANS(2)  = HE12s3S(FREQ)
  IF (IMIN .GT. 0 .AND. IMIN .LE. 3)  TRANS(3)  = HE12s1S(FREQ)
  IF (IMIN .GT. 0 .AND. IMIN .LE. 4)  TRANS(4)  = HE12p3P(FREQ)
  IF (IMIN .GT. 0 .AND. IMIN .LE. 5)  TRANS(5)  = HE12p1P(FREQ)
  ! 1s3s 3S
  IF (IMIN .GT. 0 .AND. IMIN .LE. 6)  TRANS(6)  = XKARZAS(FREQ, 1.236439D0, 3, 0)
  ! 1s3s 1S
  IF (IMIN .GT. 0 .AND. IMIN .LE. 7)  TRANS(7)  = XKARZAS(FREQ, 1.102898D0, 3, 0)
  ! 1s3p 3P
  IF (IMIN .GT. 0 .AND. IMIN .LE. 8)  TRANS(8)  = XKARZAS(FREQ, 1.045499D0, 3, 1)
  ! 1s3d 3D+1D
  IF (IMIN .GT. 0 .AND. IMIN .LE. 9)  TRANS(9)  = XKARZAS(FREQ, 1.001427D0, 3, 2)
  ! 1s3p 1P
  IF (IMIN .GT. 0 .AND. IMIN .LE. 10) TRANS(10) = XKARZAS(FREQ, 0.9926D0, 3, 1)

  !---------------------------------------------------------------------
  ! Inner-shell ionization: He I excited state → He II n=2
  ! (Adds hydrogenic 1s cross-section at shifted threshold)
  !---------------------------------------------------------------------
  IF (IMIN .GE. 1) THEN
    ELIM = 527490.06D0
    ! 1s2p 1P → He II 2p
    FREQHE = (ELIM - 171135.000D0) * CLIGHT
    IF (FREQ .GE. FREQHE) THEN
      ZEFF2 = FREQHE / RYD
      TRANS(5) = TRANS(5) + XKARZAS(FREQ, ZEFF2, 1, 0)
      ! 1s2p 3P → He II 2p
      FREQHE = (ELIM - 169087.0D0) * CLIGHT
      IF (FREQ .GE. FREQHE) THEN
        ZEFF2 = FREQHE / RYD
        TRANS(4) = TRANS(4) + XKARZAS(FREQ, ZEFF2, 1, 0)
        ! 1s2s 1S → He II 2s
        FREQHE = (ELIM - 166277.546D0) * CLIGHT
        IF (FREQ .GE. FREQHE) THEN
          ZEFF2 = FREQHE / RYD
          TRANS(3) = TRANS(3) + XKARZAS(FREQ, ZEFF2, 1, 0)
          ! 1s2s 3S → He II 2s
          FREQHE = (ELIM - 159856.069D0) * CLIGHT
          IF (FREQ .GE. FREQHE) THEN
            ZEFF2 = FREQHE / RYD
            TRANS(2) = TRANS(2) + XKARZAS(FREQ, ZEFF2, 1, 0)
          END IF
        END IF
      END IF
    END IF

    !-------------------------------------------------------------------
    ! Inner-shell ionization: He I excited state → He II n=3
    !-------------------------------------------------------------------
    ELIM = 588451.59D0
    FREQHE = (ELIM - 186209.471D0) * CLIGHT
    IF (FREQ .GE. FREQHE) THEN
      ZEFF2 = FREQHE / RYD
      TRANS(10) = TRANS(10) + XKARZAS(FREQ, ZEFF2, 1, 0)
      FREQHE = (ELIM - 186101.0D0) * CLIGHT
      IF (FREQ .GE. FREQHE) THEN
        ZEFF2 = FREQHE / RYD
        TRANS(9) = TRANS(9) + XKARZAS(FREQ, ZEFF2, 1, 0)
        FREQHE = (ELIM - 185564.0D0) * CLIGHT
        IF (FREQ .GE. FREQHE) THEN
          ZEFF2 = FREQHE / RYD
          TRANS(8) = TRANS(8) + XKARZAS(FREQ, ZEFF2, 1, 0)
          FREQHE = (ELIM - 184864.0D0) * CLIGHT
          IF (FREQ .GE. FREQHE) THEN
            ZEFF2 = FREQHE / RYD
            TRANS(7) = TRANS(7) + XKARZAS(FREQ, ZEFF2, 1, 0)
            FREQHE = (ELIM - 183236.0D0) * CLIGHT
            IF (FREQ .GE. FREQHE) THEN
              ZEFF2 = FREQHE / RYD
              TRANS(6) = TRANS(6) + XKARZAS(FREQ, ZEFF2, 1, 0)
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF

  !---------------------------------------------------------------------
  ! High-n levels (n=4-27): hydrogenic with Z_eff^2 = 4 - 3/n^2
  !---------------------------------------------------------------------
  TRANSN = 0.0D0
  IF (FREQ .GE. 1.25408D16) THEN
    DO N = 4, 27
      ZEFF2 = 4.0D0 - 3.0D0 / dble(N)**2
      TRANSN(N) = XKARZAS(FREQ, ZEFF2, 1, 0)
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Assemble total opacity at each depth
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    ! Dissolved-level contribution
    EX = BOLTEX(J)
    IF (FREQ .LT. 2.055D14) EX = EXLIM(J) / EHVKT(J)
    HE1 = (EX - EXLIM(J)) * C

    ! Bound-free from resolved levels
    IF (IMIN .GT. 0) THEN
      DO N = IMIN, 10
        HE1 = HE1 + TRANS(N) * BOLT(J, N)
      END DO
    END IF

    ! High-n levels
    IF (FREQ .GE. 1.25408D16) THEN
      DO N = 4, 27
        HE1 = HE1 + TRANSN(N) * BOLTN(J, N)
      END DO
    END IF

    ! Total: bound + free-free (all LTE source function)
    !      AHE1BOUND(J) = HE1 * STIM(J)
    !      AHE1FREE(J) = (COULFF(J,1) * FREET(J) * CFREE) * STIM(J)
    AHE1(J) = (HE1 + COULFF(J, 1) * FREET(J) * CFREE) * STIM(J)
  END DO

  RETURN

END SUBROUTINE HE1OP

!=======================================================================
! CROSSHE: He I ground-state photoionization cross-section
!
! Experimental cross-section from Marr, G.V. and West, J.B.,
! Atomic Data and Nuclear Data Tables, vol 18, 497-508, 1976.
!
! Four wavelength regimes with different table spacings:
!   50–505 Å:  Δλ = 5 Å  (92 points, X505)
!   20–50  Å:  Δλ = 2 Å  (16 points, X50)
!   10–20  Å:  Δλ = 1 Å  (11 points, X20)
!    0–10  Å:  Δλ = 0.5 Å (21 points, X10)
! Threshold at 504 Å (ν = 5.945209 × 10¹⁵ Hz, He I ionization
! at 24.587 eV). Cross-sections in units of 10⁻¹⁸ cm².
!=======================================================================

FUNCTION CROSSHE(FREQ)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ
  REAL(8) :: CROSSHE

  REAL(8), PARAMETER :: X505(92) = (/ &
    7.58D0, 7.46D0, 7.33D0, 7.19D0, 7.06D0, 6.94D0, 6.81D0, &
    6.68D0, 6.55D0, 6.43D0, 6.30D0, 6.18D0, 6.05D0, 5.93D0, &
    5.81D0, 5.69D0, 5.57D0, 5.45D0, 5.33D0, 5.21D0, 5.10D0, &
    4.98D0, 4.87D0, 4.76D0, 4.64D0, 4.53D0, 4.42D0, 4.31D0, &
    4.20D0, 4.09D0, 4.00D0, 3.88D0, 3.78D0, 3.68D0, 3.57D0, &
    3.47D0, 3.37D0, 3.27D0, 3.18D0, 3.08D0, 2.98D0, 2.89D0, &
    2.80D0, 2.70D0, 2.61D0, 2.52D0, 2.44D0, 2.35D0, 2.26D0, &
    2.18D0, 2.10D0, 2.02D0, 1.94D0, 1.86D0, 1.78D0, 1.70D0, &
    1.63D0, 1.55D0, 1.48D0, 1.41D0, 1.34D0, 1.28D0, 1.21D0, &
    1.14D0, 1.08D0, 1.02D0, 0.961D0, 0.903D0, 0.847D0, 0.792D0, &
    0.738D0, 0.687D0, 0.637D0, 0.588D0, 0.542D0, 0.497D0, &
    0.454D0, 0.412D0, 0.373D0, 0.335D0, 0.299D0, 0.265D0, &
    0.233D0, 0.202D0, 0.174D0, 0.147D0, 0.123D0, 0.100D0, &
    0.0795D0, 0.0609D0, 0.0443D0, 0.0315D0 /)

  REAL(8), PARAMETER :: X50(16) = (/ &
    0.0315D0, 0.0282D0, 0.0250D0, 0.0220D0, 0.0193D0, 0.0168D0, &
    0.0145D0, 0.0124D0, 0.0105D0, 0.00885D0, 0.00736D0, &
    0.00604D0, 0.00489D0, 0.00389D0, 0.00303D0, 0.00231D0 /)

  REAL(8), PARAMETER :: X20(11) = (/ &
    0.00231D0, 0.00199D0, 0.00171D0, 0.00145D0, 0.00122D0, &
    0.00101D0, 0.000832D0, 0.000673D0, 0.000535D0, 0.000417D0, &
    0.000318D0 /)

  REAL(8), PARAMETER :: X10(21) = (/ &
    0.000318D0, 0.000274D0, 0.000235D0, 0.000200D0, 0.000168D0, &
    0.000139D0, 0.000115D0, 0.000093D0, 0.000074D0, 0.000057D0, &
    0.000044D0, 0.000032D0, 0.000023D0, 0.000016D0, 0.000010D0, &
    0.000006D0, 0.000003D0, 0.000001D0, 0.0000006D0, &
    0.0000003D0, 0.0D0 /)

  REAL(8)  :: WAVE
  INTEGER :: I

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING CROSSHE'

  CROSSHE = 0.0D0
  IF (FREQ .LT. 5.945209D15) RETURN

  WAVE = CLIGHT_ANGS / FREQ

  IF (WAVE .GT. 50.0D0) THEN
    ! 50–505 Å regime (Δλ = 5 Å)
    I = int(93.0D0 - (WAVE - 50.0D0) / 5.0D0)
    I = min(92, max(2, I))
    CROSSHE = ((WAVE - (92 - I) * 5.0D0 - 50.0D0) / 5.0D0 &
             * (X505(I-1) - X505(I)) + X505(I)) * 1.D-18
  ELSE IF (WAVE .GT. 20.0D0) THEN
    ! 20–50 Å regime (Δλ = 2 Å)
    I = int(17.0D0 - (WAVE - 20.0D0) / 2.0D0)
    I = min(16, max(2, I))
    CROSSHE = ((WAVE - (16 - I) * 2.0D0 - 20.0D0) / 2.0D0 &
             * (X50(I-1) - X50(I)) + X50(I)) * 1.D-18
  ELSE IF (WAVE .GT. 10.0D0) THEN
    ! 10–20 Å regime (Δλ = 1 Å)
    I = int(12.0D0 - (WAVE - 10.0D0) / 1.0D0)
    I = min(11, max(2, I))
    CROSSHE = ((WAVE - (11 - I) * 1.0D0 - 10.0D0) / 1.0D0 &
             * (X20(I-1) - X20(I)) + X20(I)) * 1.D-18
  ELSE
    ! 0–10 Å regime (Δλ = 0.5 Å)
    I = int(22.0D0 - WAVE / 0.5D0)
    I = min(21, max(2, I))
    CROSSHE = ((WAVE - (21 - I) * 0.5D0) / 0.5D0 &
             * (X10(I-1) - X10(I)) + X10(I)) * 1.D-18
  END IF
  RETURN

END FUNCTION CROSSHE

!=======================================================================
! HE111S: He I 1s² ¹S ground-state bound-free cross-section
!
! Photoionization cross-section for He I from the ground state,
! after Mathisen. 64-point table of σ (in units of 10⁻¹⁸ cm²) vs
! wavelength (Å), linearly interpolated. Threshold at λ = 504.3 Å
! (ν = 5.945209 × 10¹⁵ Hz, corresponding to He I ionization at
! 24.587 eV).
!=======================================================================

FUNCTION HE111S(FREQ)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ
  REAL(8) :: HE111S

  ! He I 1s² ¹S ground-state bound-free cross-section (after Mathisen)
  ! Linear interpolation in σ vs wavelength, 64-point table
  INTEGER, PARAMETER :: NP = 64
  REAL(8), PARAMETER :: W(64) = (/ &
    504.3D0, 501.5D0, 498.7D0, 493.3D0, 488.1D0, 480.3D0, 477.8D0, &
    454.0D0, 443.0D0, 395.0D0, 356.4D0, 348.2D0, 324.6D0, 302.0D0, &
    298.1D0, 275.6D0, 260.6D0, 256.2D0, 239.4D0, 224.6D0, 220.0D0, &
    215.0D0, 210.0D0, 205.0D0, 200.0D0, 195.0D0, 190.0D0, 185.0D0, &
    180.0D0, 175.0D0, 170.0D0, 165.0D0, 160.0D0, 155.0D0, 150.0D0, &
    145.0D0, 135.0D0, 130.0D0, 125.0D0, 120.0D0, 115.0D0, 110.0D0, &
    105.0D0, 100.0D0, 95.0D0, 90.0D0, 85.0D0, 80.0D0, 75.0D0, &
    70.0D0, 65.0D0, 60.0D0, 55.0D0, 50.0D0, 45.0D0, 40.0D0, &
    35.0D0, 30.0D0, 25.0D0, 20.0D0, 15.0D0, 10.0D0, 5.0D0, 0.0D0 /)
  REAL(8), PARAMETER :: X(64) = (/ &
    7.346D0, 7.317D0, 7.259D0, 7.143D0, 7.030D0, 6.857D0, 6.800D0, &
    6.284D0, 6.041D0, 4.977D0, 4.138D0, 3.961D0, 3.474D0, 3.025D0, &
    2.945D0, 2.522D0, 2.259D0, 2.179D0, 1.901D0, 1.684D0, 1.61D0, &
    1.53D0, 1.45D0, 1.38D0, 1.30D0, 1.22D0, 1.14D0, 1.08D0, &
    1.02D0, 0.961D0, 0.903D0, 0.847D0, 0.792D0, 0.738D0, 0.687D0, &
    0.637D0, 0.542D0, 0.497D0, 0.454D0, 0.412D0, 0.373D0, 0.335D0, &
    0.299D0, 0.265D0, 0.233D0, 0.202D0, 0.174D0, 0.147D0, 0.124D0, &
    0.103D0, 0.0840D0, 0.0676D0, 0.0535D0, 0.0414D0, 0.0311D0, &
    0.0266D0, 0.0158D0, 0.0104D0, 0.00637D0, 0.00349D0, 0.00161D0, &
    0.00054D0, 0.000083D0, 0.0D0 /)

  REAL(8)  :: WAVE
  INTEGER :: I

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE111S'

  HE111S = 0.0D0
  IF (FREQ .LT. 5.945209D15) RETURN

  WAVE = CLIGHT_ANGS / FREQ
  DO I = 2, NP
    IF (WAVE .GT. W(I)) EXIT
  END DO
  IF (I .GT. NP) I = NP

  HE111S = ((WAVE - W(I)) / (W(I-1) - W(I)) * (X(I-1) - X(I)) + X(I)) * 1.0D-18
  RETURN

END FUNCTION HE111S

!=======================================================================
! HE12S1S: He I 1s2s ¹S bound-free photoionization cross-section
!
! Tabulated cross-section for photoionization of He I from the
! 1s2s ¹S metastable state. 16-point table of log₁₀(σ) vs log₁₀(ν),
! linearly interpolated. Threshold at 32033.214 cm⁻¹.
! Above 2.4 Rydberg: Fano resonance profile formula with single
! autoionization resonance (¹S channel, q = 76.21).
!=======================================================================

FUNCTION HE12S1S(FREQ)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ
  REAL(8) :: HE12S1S

  ! He I 1s2s ¹S bound-free cross-section
  ! Table interpolation below 2.4 Rydberg, resonance formula above
  INTEGER, PARAMETER :: NP = 16
  REAL(8), PARAMETER :: FREQ1S(16) = (/ &
    15.947182D0, 15.913654D0, 15.877320D0, 15.837666D0, 15.794025D0, &
    15.745503D0, 15.690869D0, 15.628361D0, 15.555317D0, 15.467455D0, &
    15.357189D0, 15.289399D0, 15.251073D0, 15.209035D0, 15.162487D0, &
    14.982421D0 /)
  REAL(8), PARAMETER :: X1S(16) = (/ &
    -19.635557D0, -19.159345D0, -18.958474D0, -18.809535D0, &
    -18.676481D0, -18.546006D0, -18.410962D0, -18.264821D0, &
    -18.100205D0, -17.909165D0, -17.684370D0, -17.557867D0, &
    -17.490360D0, -17.417876D0, -17.349386D0, -17.084441D0 /)

  REAL(8)  :: FREQLG, X, EK, EPS, WAVNO
  INTEGER :: I

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE12S1S'

  HE12S1S = 0.0D0
  IF (FREQ .LT. 32033.214D0 * CLIGHT) RETURN

  IF (FREQ .GT. 2.4D0 * RYDBERG_HE * CLIGHT) THEN
    ! High-energy resonance formula
    WAVNO = FREQ / CLIGHT
    EK = (WAVNO - 32033.214D0) / RYDBERG_HE
    EPS = 2.0D0 * (EK - 2.612316D0) / 0.00322D0
    HE12S1S = 0.008175D0 * (484940.0D0 / WAVNO)**2.71D0 * 8.067D-18 &
            * (EPS + 76.21D0)**2 / (1.0D0 + EPS**2)
  ELSE
    ! Table interpolation in log₁₀(σ) vs log₁₀(ν)
    FREQLG = log10(FREQ)
    DO I = 2, NP
      IF (FREQLG .GT. FREQ1S(I)) EXIT
    END DO
    IF (I .GT. NP) I = NP
    X = (FREQLG - FREQ1S(I)) / (FREQ1S(I-1) - FREQ1S(I)) &
      * (X1S(I-1) - X1S(I)) + X1S(I)
    HE12S1S = exp(X * LN10)
  END IF
  RETURN

END FUNCTION HE12S1S

!=======================================================================
! HE12S3S: He I 1s2s ³S bound-free photoionization cross-section
!
! Tabulated cross-section for photoionization of He I from the
! 1s2s ³S metastable state. 16-point table of log₁₀(σ) vs log₁₀(ν),
! linearly interpolated. Threshold at 38454.691 cm⁻¹.
! Above 2.4 Rydberg: Fano resonance profile formula with single
! autoionization resonance (³S channel, q = -122.4).
!=======================================================================

FUNCTION HE12S3S(FREQ)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ
  REAL(8) :: HE12S3S

  ! He I 1s2s ³S bound-free cross-section
  ! Table interpolation below 2.4 Rydberg, resonance formula above
  INTEGER, PARAMETER :: NP = 16
  REAL(8), PARAMETER :: FREQ3S(16) = (/ &
    15.956523D0, 15.923736D0, 15.888271D0, 15.849649D0, 15.807255D0, &
    15.760271D0, 15.707580D0, 15.647601D0, 15.577992D0, 15.495055D0, &
    15.392451D0, 15.330345D0, 15.295609D0, 15.257851D0, 15.216496D0, &
    15.061770D0 /)
  REAL(8), PARAMETER :: X3S(16) = (/ &
    -18.426022D0, -18.610700D0, -18.593051D0, -18.543304D0, &
    -18.465513D0, -18.378707D0, -18.278574D0, -18.164329D0, &
    -18.033346D0, -17.882435D0, -17.705542D0, -17.605584D0, &
    -17.553459D0, -17.500667D0, -17.451318D0, -17.266686D0 /)

  REAL(8)  :: FREQLG, X, EK, EPS, WAVNO
  INTEGER :: I

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE12S3S'

  HE12S3S = 0.0D0
  IF (FREQ .LT. 38454.691D0 * CLIGHT) RETURN

  IF (FREQ .GT. 2.4D0 * RYDBERG_HE * CLIGHT) THEN
    ! High-energy resonance formula
    WAVNO = FREQ / CLIGHT
    EK = (WAVNO - 38454.691D0) / RYDBERG_HE
    EPS = 2.0D0 * (EK - 2.47898D0) / 0.000780D0
    HE12S3S = 0.01521D0 * (470310.0D0 / WAVNO)**3.12D0 * 8.067D-18 &
            * (EPS - 122.4D0)**2 / (1.0D0 + EPS**2)
  ELSE
    ! Table interpolation in log₁₀(σ) vs log₁₀(ν)
    FREQLG = log10(FREQ)
    DO I = 2, NP
      IF (FREQLG .GT. FREQ3S(I)) EXIT
    END DO
    IF (I .GT. NP) I = NP
    X = (FREQLG - FREQ3S(I)) / (FREQ3S(I-1) - FREQ3S(I)) &
      * (X3S(I-1) - X3S(I)) + X3S(I)
    HE12S3S = exp(X * LN10)
  END IF
  RETURN

END FUNCTION HE12S3S

!=======================================================================
! HE12P1P: He I 1s2p ¹P bound-free photoionization cross-section
!
! Tabulated cross-section for photoionization of He I from the
! 1s2p ¹P excited state. 16-point table of log₁₀(σ) vs log₁₀(ν),
! linearly interpolated. Threshold at 27175.76 cm⁻¹.
! Above 2.4 Rydberg: Fano resonance profile formula with two
! autoionization channels (¹S with q = -29.30, ¹D with q = 172.4).
!=======================================================================

FUNCTION HE12P1P(FREQ)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ
  REAL(8) :: HE12P1P

  ! He I 1s2p ¹P bound-free cross-section
  ! Table interpolation below 2.4 Rydberg, two-resonance formula above
  INTEGER, PARAMETER :: NP = 16
  REAL(8), PARAMETER :: FREQ1P(16) = (/ &
    15.939981D0, 15.905870D0, 15.868850D0, 15.828377D0, 15.783742D0, &
    15.733988D0, 15.677787D0, 15.613218D0, 15.537343D0, 15.445346D0, &
    15.328474D0, 15.255641D0, 15.214064D0, 15.168081D0, 15.116647D0, &
    14.911002D0 /)
  REAL(8), PARAMETER :: X1P(16) = (/ &
    -18.798876D0, -19.685922D0, -20.011664D0, -20.143030D0, &
    -20.091354D0, -19.908333D0, -19.656788D0, -19.367745D0, &
    -19.043016D0, -18.674484D0, -18.240861D0, -17.989700D0, &
    -17.852015D0, -17.702677D0, -17.525347D0, -16.816344D0 /)

  REAL(8)  :: FREQLG, X, EK, EPS1S, EPS1D, WAVNO
  INTEGER :: I

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE12P1P'

  HE12P1P = 0.0D0
  IF (FREQ .LT. 27175.76D0 * CLIGHT) RETURN

  IF (FREQ .GT. 2.4D0 * RYDBERG_HE * CLIGHT) THEN
    ! High-energy: two autoionization resonances (¹S and ¹D channels)
    WAVNO = FREQ / CLIGHT
    EK = (WAVNO - 27175.76D0) / RYDBERG_HE
    EPS1S = 2.0D0 * (EK - 2.446534D0) / 0.01037D0
    EPS1D = 2.0D0 * (EK - 2.59427D0) / 0.00538D0
    HE12P1P = 0.0009487D0 * (466750.0D0 / WAVNO)**3.69D0 * 8.067D-18 &
            * ((EPS1S - 29.30D0)**2 / (1.0D0 + EPS1S**2) &
            + (EPS1D + 172.4D0)**2 / (1.0D0 + EPS1D**2))
  ELSE
    ! Table interpolation in log₁₀(σ) vs log₁₀(ν)
    FREQLG = log10(FREQ)
    DO I = 2, NP
      IF (FREQLG .GT. FREQ1P(I)) EXIT
    END DO
    IF (I .GT. NP) I = NP
    X = (FREQLG - FREQ1P(I)) / (FREQ1P(I-1) - FREQ1P(I)) &
      * (X1P(I-1) - X1P(I)) + X1P(I)
    HE12P1P = exp(X * LN10)
  END IF
  RETURN

END FUNCTION HE12P1P

!=======================================================================
! HE12P3P: He I 1s2p ³P bound-free photoionization cross-section
!
! Tabulated cross-section for photoionization of He I from the
! 1s2p ³P excited state. 16-point table of log₁₀(σ) vs log₁₀(ν),
! linearly interpolated. Threshold at 29223.753 cm⁻¹.
!=======================================================================

FUNCTION HE12P3P(FREQ)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ
  REAL(8) :: HE12P3P

  ! He I 1s2p ³P bound-free cross-section
  ! Linear interpolation in log₁₀(σ) vs log₁₀(ν), 16-point table
  INTEGER, PARAMETER :: NP = 16
  REAL(8), PARAMETER :: FREQ3P(16) = (/ &
    15.943031D0, 15.909169D0, 15.872441D0, 15.832318D0, 15.788107D0, &
    15.738880D0, 15.683351D0, 15.619667D0, 15.545012D0, 15.454805D0, &
    15.340813D0, 15.270195D0, 15.230054D0, 15.185821D0, 15.136567D0, &
    14.942557D0 /)
  REAL(8), PARAMETER :: X3P(16) = (/ &
    -19.791021D0, -19.697886D0, -19.591421D0, -19.471855D0, &
    -19.337053D0, -19.183958D0, -19.009750D0, -18.807990D0, &
    -18.570571D0, -18.288361D0, -17.943476D0, -17.738737D0, &
    -17.624154D0, -17.497163D0, -17.403183D0, -17.032999D0 /)

  REAL(8)  :: FREQLG, X
  INTEGER :: I

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE12P3P'

  HE12P3P = 0.0D0
  IF (FREQ .LT. 29223.753D0 * CLIGHT) RETURN

  FREQLG = log10(FREQ)
  DO I = 2, NP
    IF (FREQLG .GT. FREQ3P(I)) EXIT
  END DO
  IF (I .GT. NP) I = NP

  X = (FREQLG - FREQ3P(I)) / (FREQ3P(I-1) - FREQ3P(I)) &
    * (X3P(I-1) - X3P(I)) + X3P(I)
  HE12P3P = exp(X * LN10)
  RETURN

END FUNCTION HE12P3P

!=======================================================================
! HE2OP: He II bound-free and free-free opacity
!
! Hydrogenic opacity for singly-ionized helium (Z=2). Frequencies
! are 4× hydrogen, ionization potential = 54.403 eV.
!
! Bound-free: sum over levels n=1–9 using XKARZAS cross-sections
!   with Z_eff² = 4, plus hydrogenic high-n contribution from the
!   series limit. Below 1.31522 × 10¹⁴ Hz the high-n term uses
!   EXLIM/EHVKT instead of the direct Boltzmann factor.
!
! Free-free: Kramers formula with Gaunt factor from COULFF(J,2).
!
! Temperature-dependent Boltzmann populations cached when ITEMP changes.
!=======================================================================

!==========================================================================
! SUBROUTINE HE2OP
!
! He II bound-free and free-free opacity.
!
! Full level-by-level version matching atlas7lib.for.  He II is
! hydrogenic (Z=2), so the structure parallels HOP but with:
!   - Ionization limit  = 4 × Rydberg = 438908.85 cm⁻¹ (54.403 eV)
!   - 4 × Rydberg const = 438889.068 cm⁻¹
!   - Z⁴ = 16 in the cross-section scaling
!   - Z² = 4 in the free-free scaling
!
! Structure:
!   n >= 10 : dissolved-level integral, LTE, S = B_nu
!   n = 7-9 : explicit hydrogenic (Kramers: ν³×Z⁴/n⁵), LTE
!   n = 1-6 : Gaunt-factor-corrected hydrogenic, with departure
!             coefficients.  In LTE, (BHE2-EHVKT) = STIM.
!   free-free: Kramers with Z²=4, Gaunt factor COULFF(J,2),
!              population n_e × n(He III).
!
! Note: line 5611 of atlas7lib.for marks a typo correction:
!   n=3 denominator is 243 (= 3⁵), not 343 (= 7³).
!==========================================================================

SUBROUTINE HE2OP

  IMPLICIT NONE

  REAL(8)  :: FREQ3, XNFPRHO
  REAL(8)  :: H, S, A, X
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HE2OP'

  FREQ3 = 2.815D29 / FREQ / FREQ / FREQ

  DO J = 1, NRHOX

    ! Population factor: He II (mode-11) / rho
    XNFPRHO = XNFP(J, 4) / RHO(J)

    ! --- n >= 10 to infinity: dissolved-level integral (LTE) ---
    ! Z⁴=16, gfactor=2, ionization limit = 438908.85 cm⁻¹
    ! Lower integration bound at n=10: E_10 = 438908.85 - 438908.85/100 = 434519.959
    H = FREQ3 * 16.0D0 * 2.0D0 / 2.0D0 / (438889.068D0 * HCKT(J)) &
      * (EXP(-MAX(434519.959D0, 438908.85D0 - WAVENO) * HCKT(J)) &
       - EXP(-438908.85D0 * HCKT(J))) * STIM(J) * XNFPRHO
    S = H * BNU(J)

    levels: DO

      ! --- n = 9 (threshold 5418.390 cm⁻¹) : hydrogenic ---
      IF (WAVENO .LT. 5418.390D0) EXIT levels
      X = FREQ3 / 59049.0D0 * 16.0D0
      A = X * 162.0D0 * EXP(-433490.46D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 8 (threshold 6857.660 cm⁻¹) : hydrogenic ---
      IF (WAVENO .LT. 6857.660D0) EXIT levels
      X = FREQ3 * 16.0D0 / 32768.0D0
      A = X * 128.0D0 * EXP(-432051.19D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 7 (threshold 8956.950 cm⁻¹) : hydrogenic ---
      IF (WAVENO .LT. 8956.950D0) EXIT levels
      X = FREQ3 * 16.0D0 / 16807.0D0
      A = X * 98.0D0 * EXP(-429951.90D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 6 (threshold 12191.437 cm⁻¹) : Gaunt-corrected, LTE ---
      IF (WAVENO .LT. 12191.437D0) EXIT levels
      X = FREQ3 * 16.0D0 / 7776.0D0 &
        * (1.0986D0 + (-2.704D13 + 1.229D27 / FREQ) / FREQ)
      A = X * 72.0D0 * EXP(-426717.413D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 5 (threshold 17555.715 cm⁻¹) : Gaunt-corrected, LTE ---
      IF (WAVENO .LT. 17555.715D0) EXIT levels
      X = FREQ3 * 16.0D0 / 3125.0D0 &
        * (1.102D0 + (-3.909D13 + 2.371D27 / FREQ) / FREQ)
      A = X * 50.0D0 * EXP(-421353.135D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 4 (threshold 27430.925 cm⁻¹) : Gaunt-corrected, LTE ---
      IF (WAVENO .LT. 27430.925D0) EXIT levels
      X = FREQ3 * 16.0D0 / 1024.0D0 &
        * (1.101D0 + (-5.765D13 + 4.593D27 / FREQ) / FREQ)
      A = X * 32.0D0 * EXP(-411477.925D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 3 (threshold 48766.491 cm⁻¹) : Gaunt-corrected, LTE ---
      ! Note: denominator is 243 (= 3⁵), not 343 (typo corrected in atlas7lib)
      IF (WAVENO .LT. 48766.491D0) EXIT levels
      X = FREQ3 * 16.0D0 / 243.0D0 &
        * (1.101D0 + (-9.863D13 + 1.035D28 / FREQ) / FREQ)
      A = X * 18.0D0 * EXP(-390142.359D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 2 (threshold 109726.529 cm⁻¹) : Gaunt-corrected, LTE ---
      IF (WAVENO .LT. 109726.529D0) EXIT levels
      X = FREQ3 * 16.0D0 / 32.0D0 &
        * (1.105D0 + (-2.375D14 + 4.077D28 / FREQ) / FREQ)
      A = X * 8.0D0 * EXP(-329182.321D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 1 (threshold 438908.850 cm⁻¹) : Gaunt-corrected, LTE ---
      IF (WAVENO .LT. 438908.850D0) EXIT levels
      X = FREQ3 * 16.0D0 / 1.0D0 &
        * (0.9916D0 + (2.719D13 - 2.268D30 / FREQ) / FREQ)
      A = X * 2.0D0 * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      EXIT levels
    END DO levels

    ! --- Free-free (bremsstrahlung, Z²=4) ---
    A = COEFF_FF * 4.0D0 / SQRT(T(J)) * COULFF(J, 2) / FREQ * XNE(J) &
      / FREQ * XNFP(J, 5) / FREQ * STIM(J) / RHO(J)
    H = H + A
    S = S + A * BNU(J)

    AHE2(J) = H
    SHE2(J) = BNU(J)
    IF (H .GT. 0.0D0) SHE2(J) = S / H

  END DO

  RETURN

END SUBROUTINE HE2OP

!=======================================================================
! HEMIOP
!=======================================================================

SUBROUTINE HEMIOP

  IMPLICIT NONE
  REAL(8)  :: A, B, C
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HEMIOP'
  ! He⁻ free-free opacity: polynomial fit in 1/freq and T
  A = 3.397D-46 + (-5.216D-31 + 7.039D-15 / FREQ) / FREQ
  B = -4.116D-42 + (1.067D-26 + 8.135D-11 / FREQ) / FREQ
  C = 5.081D-37 + (-8.724D-23 - 5.659D-8 / FREQ) / FREQ
  AHEMIN = (A * T + B + C / T) * XNE * XNFP(:, 3) / RHO
  RETURN

END SUBROUTINE HEMIOP

!=======================================================================
! HERAOP: He I Rayleigh scattering opacity
!
! Rayleigh scattering by ground-state neutral helium.  Cross-section
! is a two-parameter dispersion-corrected Rayleigh formula:
!
!   σ(λ) = (5.484e-14 / λ⁴) * [1 + (2.44e5 + 5.94e10/(λ²-2.90e5)) / λ²]²
!
! with λ in Å and σ in cm².  The formula has two distinct pieces:
!
! Leading (long-wavelength) term:
!   At λ → ∞, the bracketed dispersion factor → 1 and the cross section
!   reduces to 5.484e-14 / λ⁴.  This coefficient is exactly the textbook
!   low-frequency Rayleigh limit
!
!       σ(ω→0) = (8π/3)(ω/c)⁴ α²    =    (128 π⁵ / 3) α² / λ⁴
!
!   evaluated with the He ground-state static dipole polarizability
!   α_He = 0.20494e-24 cm³ (≈ 1.3832 a₀³).  Plugging in this α value
!   reproduces 5.484e-14 to 4 significant figures.  So the long-
!   wavelength behaviour is correct by construction and tied to the
!   fundamental polarizability of helium.
!
! Dispersion correction (short-wavelength factor):
!   The bracketed factor [α(ω)/α(0)]² accounts for the frequency
!   dependence of the helium polarizability as the photon energy
!   approaches the first excited bound states.  The two fitted
!   constants (2.44e5 and 5.94e10 Å², plus the pole at 2.90e5 Å²) are
!   NOT sums over physical He line positions: the pole at λ²=2.90e5 Å²
!   corresponds to λ ≈ 538 Å, which is between but not at the
!   He I 1s²→1snp resonances (584 Å, 537 Å, 522 Å, ...).  It behaves
!   like an effective oscillator-strength-weighted resonance wavelength.
!   The coefficients almost certainly trace to a compact fit of the
!   form used by Dalgarno in the early 1960s (Cauchy-moment expansion
!   of α(ω) from oscillator-strength sums, fitted to a two-pole Padé
!   form), but the precise source of these four numbers is not cited
!   in the original F77 Kurucz code and has not been identified.  The
!   most likely origin is Dalgarno (1962), "Spectral Reflectivity of
!   the Earth's Atmosphere III: The Scattering of Light by Atomic
!   Systems," Geophysical Corporation of America Report, which is
!   cited by later work (e.g. Rohrmann 2018 MNRAS 473, 457) as the
!   standard astrophysical He Rayleigh reference of that era but is
!   a technical report rather than a journal article.
!
! Validity and accuracy:
!   The leading 1/λ⁴ term is exact in the long-wavelength limit.  The
!   size of the dispersion correction itself (i.e., how much [α(ω)/α(0)]²
!   differs from unity) is about 2% at 5000 Å, 6% at 3000 Å, 13% at
!   2000 Å, and grows rapidly in the FUV as the 584 Å resonance is
!   approached.  The error of the fit relative to a modern calculation
!   (e.g. Rohrmann 2018) has not been quantified here but is expected
!   to be at the sub-percent level in the visible/near-IR and to grow
!   toward the UV.  He Rayleigh is a negligible opacity source in FGK
!   stars regardless of formula choice; it matters only in the UV of
!   A/B stars and in He-rich chemically peculiar stars.
!
! Frequency cap:
!   The evaluated frequency is capped at 5.15e15 Hz (λ ≈ 582 Å) to
!   prevent evaluation blueward of the He I 584 Å resonance, where
!   the Kurucz dispersion factor would go singular (the denominator
!   λ²-2.90e5 changes sign) and where the Rayleigh picture is no
!   longer physically valid anyway (real bound-bound absorption
!   dominates).  For ν > 5.15e15 Hz the cross section is frozen at
!   its λ = 582 Å value.
!
! Future work:
!   If increased precision in the FUV becomes necessary — e.g. for
!   detailed atmosphere modelling of He-rich CP stars, sdB stars, or
!   cool helium-atmosphere white dwarfs where He Rayleigh can be a
!   non-negligible opacity source — the current two-parameter fit
!   should be replaced with Rohrmann (2018, MNRAS 473, 457).  That
!   work evaluates the He polarizability using a comprehensive
!   oscillator-strength distribution (including doubly excited states
!   and the photoionization continuum), tabulates the full Rayleigh
!   cross section including the first two physical resonances, and
!   provides analytical fits valid for λ > 505 Å.  For warmer-gas
!   applications, Sneep & Ubachs (2005, JQSRT 92, 293) and Thalman
!   et al. (2014, JQSRT 147, 171) provide refractive-index-based
!   cross sections calibrated to modern lab measurements; these are
!   accurate in the visible/near-UV but do not extend into the FUV
!   resonance region.
!=======================================================================

SUBROUTINE HERAOP
  IMPLICIT NONE
  REAL(8)  :: W, WW, SIG
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HERAOP'

  ! Wavelength in Å, capped just redward of the He I 584 Å resonance
  ! to prevent the dispersion denominator from going singular.
  W  = CLIGHT_ANGS / min(FREQ, 5.15D15)
  WW = W**2                                   ! λ² [Å²]

  ! σ(λ) = (5.484e-14 / λ⁴) * [α(ω)/α(0)]²
  !   Leading coefficient 5.484e-14 = (128π⁵/3) * α_He(static)²
  !   Bracketed factor is the fitted two-pole dispersion correction.
  SIG = 5.484D-14 / WW / WW * (1.0D0 + (2.44D5 + 5.94D10 &
      / (WW - 2.90D5)) / WW)**2

  ! Apply to each depth point using neutral He number density
  ! (XNFP(J,3) is N(He I)); /RHO converts to mass opacity [cm²/g].
  SIGHE = SIG * XNFP(:, 3) / RHO
  RETURN

END SUBROUTINE HERAOP

!=========================================================================
! SUBROUTINE CONT_METAL_OPACITY_LEGACY
!
! Consolidated continuum-metal opacity using Kurucz's legacy hand-coded
! cross sections (Luo & Pradhan + Burke & Taylor Fanos for C I, Seaton
! formulas for Mg I / Al I / Si I / Fe I / Ca II, Peterson fits for
! N I / O I / C II / Mg II / Si II, and the 60-edge HOTOP modified-
! Seaton table).  Also includes H2-H2 collision-induced absorption
! and CH / OH molecular opacity.
!
! This routine replaces the former COOLOP + WARMOP + HOTOP wrappers.
! It exists so that a TOPbase-based variant (CONT_METAL_OPACITY_TOPBASE)
! can be dropped in as an alternative metal-continuum backend without
! disturbing the surrounding KAPP dispatch.
!
! Output: ACONT_METAL(J) in cm^2/g, stimulated-emission corrected,
! ready to sum into the total continuum absorption in KAPP.
!
! Structure:
!   1. Cool block (below Lyman limit only): neutral-metal bound-free
!      (C I, Mg I, Al I, Si I, Fe I), H2-H2 CIA, CH and OH molecular.
!      Gated by FREQ <= FREQ_RYDH because neutral-metal thresholds all
!      lie below the Lyman limit and H bound-free dominates above it.
!   2. Warm block (all frequencies): first-ion bound-free (C II, Mg II,
!      Si II) plus inline N I, O I, Ca II cross sections times
!      populations.
!   3. Hot block (all frequencies): ion free-free for C/N/O/Ne/Mg/Si/S/Fe
!      stages I-V with Coulomb Gaunt factors, plus the 60 modified-
!      Seaton bound-free edges from hotop.dat.  Accumulates into a
!      local array with the "< 1% of running total" skip optimization
!      preserved bit-for-bit from the original HOTOP, then applies
!      STIM/RHO before adding to ACONT_METAL.
!=========================================================================

SUBROUTINE CONT_METAL_OPACITY_LEGACY

  IMPLICIT NONE

  ! --- HOTOP persistent state (preserved across calls) ---
  INTEGER, PARAMETER :: NUM = 60    ! number of bound-free edges
  INTEGER, PARAMETER :: NPAR = 7    ! parameters per edge
  REAL(8),  SAVE :: A_HOT(NPAR, NUM)
  LOGICAL,  SAVE :: HOT_INITIALIZED = .FALSE.

  ! --- Locals ---
  REAL(8)  :: AHOT_TMP(kw)          ! running HOT total before STIM/RHO
  REAL(8)  :: AC2OP(kw)             ! C II b-f placeholder (disabled)
  REAL(8)  :: FREE, XSECT, XX, FRATIO
  INTEGER  :: I, J, ID
  CHARACTER(256) :: LINE

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING CONT_METAL_OPACITY_LEGACY'

  ! Zero output
  ACONT_METAL = 0.0D0

  !======================================================================
  ! (1) Cool block: below Lyman limit only.
  !     Neutrals C I, Mg I, Al I, Si I, Fe I + H2-H2 CIA + CH/OH molec.
  !======================================================================
  IF (FREQ .LE. FREQ_RYDH) THEN
    CALL C1OP
    CALL MG1OP
    CALL AL1OP
    CALL SI1OP
    CALL FE1OP
    CALL H2COLLOP(AH2COLL)        ! fills module-scope AH2COLL(kw)
    DO J = 1, NRHOX
      ACONT_METAL(J) = ACONT_METAL(J)                                    &
        + AC1(J) + AMG1(J) + AAL1(J) + ASI1(J) + AFE1(J)                 &
        + AH2COLL(J)                                                     &
        + (CHOP(J) * XNFP(J, 846) + OHOP(J) * XNFP(J, 848))              &
          * STIM(J) / RHO(J)
    END DO
  END IF

  !======================================================================
  ! (2) Warm block: all frequencies.
  !     First-ion bound-free C II / Mg II / Si II + inline N I / O I /
  !     Ca II cross sections * populations * STIM / RHO.
  !======================================================================
  CALL C2OP
  CALL MG2OP
  CALL SI2OP
  DO J = 1, NRHOX
    ACONT_METAL(J) = ACONT_METAL(J)                                      &
      + AC2(J) + AMG2(J) + ASI2(J)                                       &
      + (N1OP(J) * XNFP(J, 28) + O1OP(J) * XNFP(J, 36)                   &
       + CA2OP(J) * XNFP(J, 211)) * STIM(J) / RHO(J)
  END DO

  !======================================================================
  ! (3) Hot block: all frequencies.
  !     Ion free-free (C/N/O/Ne/Mg/Si/S/Fe I-V) + 60 bound-free edges
  !     from hotop.dat.  Accumulates into local AHOT_TMP to preserve
  !     the "< 1% running total" skip optimization in the edge loop,
  !     then applies STIM/RHO once before adding to ACONT_METAL.
  !======================================================================

  ! Read edge parameters from file on first call
  IF (.NOT. HOT_INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'hotop.dat', STATUS='OLD', ACTION='READ')
    DO
      READ(89, '(A)') LINE
      IF (LINE(1:1) .NE. '#') THEN
        BACKSPACE(89)
        EXIT
      END IF
    END DO
    DO I = 1, NUM
      READ(89, *) A_HOT(:, I)
    END DO
    CLOSE(89)
    HOT_INITIALIZED = .TRUE.
  END IF

  ! Free-free opacity: C,N,O,Ne,Mg,Si,S,Fe ions (stages I-V)
  DO J = 1, NRHOX
    FREE = COULFF(J,1) * 1.0D0                                           &
            * (XNF(J,22) + XNF(J,29) + XNF(J,37) + XNF(J,56)             &
             + XNF(J,79) + XNF(J,106) + XNF(J,137) + XNF(J,352))         &
         + COULFF(J,2) * 4.0D0                                           &
            * (XNF(J,23) + XNF(J,30) + XNF(J,38) + XNF(J,57)             &
             + XNF(J,80) + XNF(J,107) + XNF(J,138) + XNF(J,353))         &
         + COULFF(J,3) * 9.0D0                                           &
            * (XNF(J,24) + XNF(J,31) + XNF(J,39) + XNF(J,58)             &
             + XNF(J,81) + XNF(J,108) + XNF(J,139) + XNF(J,354))         &
         + COULFF(J,4) * 16.0D0                                          &
            * (XNF(J,25) + XNF(J,32) + XNF(J,40) + XNF(J,59)             &
             + XNF(J,82) + XNF(J,109) + XNF(J,140) + XNF(J,355))         &
         + COULFF(J,5) * 25.0D0                                          &
            * (XNF(J,26) + XNF(J,33) + XNF(J,41) + XNF(J,60)             &
             + XNF(J,83) + XNF(J,110) + XNF(J,141) + XNF(J,356))

    AHOT_TMP(J) = FREE * COEFF_FF / (FREQ * FREQ * FREQ) * XNE(J) / sqrt(T(J))
  END DO

  ! C II bound-free: disabled (handled in warm block per Bischoff 4 Jun 2003)
  ! The placeholder is kept for source-code parity with HOTOP.
  AC2OP = 0.0D0
  AHOT_TMP = AHOT_TMP + AC2OP

  ! 60 modified-Seaton bound-free edges.  The "< 1% of running total"
  ! skip condition is preserved exactly from the original HOTOP so that
  ! the refactor is bit-for-bit identical.
  DO I = 1, NUM
    IF (FREQ .LT. A_HOT(1, I)) CYCLE
    FRATIO = A_HOT(1, I) / FREQ
    XSECT = A_HOT(2, I) * (A_HOT(3, I) + FRATIO - A_HOT(3, I) * FRATIO)  &
          * sqrt(FRATIO ** int(A_HOT(4, I)))
    ID = int(A_HOT(7, I))
    DO J = 1, NRHOX
      XX = XSECT * XNFP(J, ID) * A_HOT(5, I)
      IF (XX .GT. AHOT_TMP(J) / 100.0D0) THEN
        AHOT_TMP(J) = AHOT_TMP(J) + XX / exp(A_HOT(6, I) / TKEV(J))
      END IF
    END DO
  END DO

  ! Apply stimulated emission and density normalization, sum into output.
  ! AHOT(J) (module-level) is kept populated for diagnostic consumers
  ! that may still reference it; nothing in the current pipeline reads
  ! AHOT after this point.
  AHOT_TMP = AHOT_TMP * STIM / RHO
  AHOT = AHOT_TMP
  ACONT_METAL = ACONT_METAL + AHOT_TMP

  RETURN

END SUBROUTINE CONT_METAL_OPACITY_LEGACY

!=========================================================================
! SUBROUTINE CONT_METAL_OPACITY_TOPBASE
!
! Dispatcher for metal continuous bound-free and free-free opacity.
!
! When USE_TOPBASE_MBF = .TRUE. (production default):
!   - MBF_TOPBASE supplies 28 species (Li I - Ca II, excluding Fe) from
!     TOPbase / Allende Prieto et al. (2003, ApJS 147, 363) processed
!     cross sections.
!   - MBF_HIGH_ION supplies high-ionization bound-free + Coulomb
!     free-free (C III+, N III+, O III+, Ne III+, Fe III-V) from
!     filtered hotop.dat.  hotop.dat contains NO iron bound-free edges
!     (only C/N/O/Ne), so no double-counting with felo can occur.
!   - FELO_OPACITY supplies Fe I and Fe II bound-free using OP / Iron
!     Project R-matrix data (Bautista 1997 for Fe I, Nahar & Pradhan
!     1994 for Fe II).
!
! When USE_TOPBASE_MBF = .FALSE. (legacy mode):
!   - Same MBF_TOPBASE and MBF_HIGH_ION.
!   - FE1OP supplies Fe I via the 48-level Kurucz Fano fit.
!   - Fe II has no dedicated bound-free treatment (known gap; Coulomb
!     free-free for Fe II is provided by MBF_HIGH_ION regardless).
!
! Output: ACONT_METAL(J) in cm^2/g, stimulated-emission corrected.
!
! Known limitation (TOPbase species only): the .xs files identify levels
! by LS term (ISLP) without resolving J.  Barklem & Collet (2016)
! partition functions used in Saha DO resolve J, so per-species opacity
! carries a factor U_TB/U_BC that is within a few percent for most
! species but reaches ~13% for Ar II, ~7% for Ne II, ~6% for Si II.
! J-resolved data would be needed for precision work on those species.
!=========================================================================

SUBROUTINE CONT_METAL_OPACITY_TOPBASE

  IMPLICIT NONE
  INTEGER :: J
  REAL(8) :: AFELO(kw)

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING CONT_METAL_OPACITY_TOPBASE'

  ! One-time initialization of TOPbase data from data/mbf/*.xs
  IF (.NOT. TB_INITIALIZED) CALL INIT_MBF_TOPBASE

  ! Zero output
  ACONT_METAL = 0.0D0

  !--- TOPbase metal bound-free (28 species) -----------------------------
  CALL MBF_TOPBASE

  !--- High-ionization bound-free + Coulomb free-free from hotop.dat -----
  CALL MBF_HIGH_ION

  !--- Fe I and Fe II bound-free -----------------------------------------
  IF (USE_TOPBASE_MBF) THEN
    CALL FELO_OPACITY(1, AFELO)                    ! Fe I
    ACONT_METAL = ACONT_METAL + AFELO
    CALL FELO_OPACITY(2, AFELO)                    ! Fe II
    ACONT_METAL = ACONT_METAL + AFELO
  ELSE
    CALL FE1OP
    ACONT_METAL = ACONT_METAL + AFE1
  END IF

  RETURN

END SUBROUTINE CONT_METAL_OPACITY_TOPBASE

!=========================================================================
! SUBROUTINE MBF_TOPBASE
!
! Runtime level-by-level photoionization opacity sum for the 30 TOPbase
! species.  For each species s with population n_s and partition function
! U_s (implicit in XNFP = n_s / U_s) the contribution is:
!
!     kappa_s(nu, T, rho) = XNFP(s) * sum_i [ g_i exp(-E_exc_i/kT) sigma_i(nu) ]
!                           * STIM / RHO
!
! (U_s cancels between Saha's n_s/U_s and the unnormalized level sum.)
!
! Performance notes:
!   - Levels whose Boltzmann factor at T=25000 K is below 1e-4 of the
!     ground-term contribution are dropped at INIT_MBF_TOPBASE time
!     (init-time truncation).  Typical 2-5x reduction in level count.
!   - Loop order is species -> level -> depth, so each level's log-grid
!     sigma(nu) interpolation is done once per frequency rather than once
!     per depth.
!   - Per-species quick reject: if the photon energy is below every
!     retained level's threshold in this species, the whole species is
!     skipped with one comparison.  Kills most species at visible/IR
!     wavelengths.
!   - Boltzmann factors g_i exp(-E_exc_i/kT_J) are cached per depth and
!     rebuilt only when T(J) changes.
!=========================================================================

SUBROUTINE MBF_TOPBASE

  IMPLICIT NONE

  INTEGER :: J, IS, IL, IFLAT, IFLAT_SPECIES_BASE
  REAL(8) :: TJ, KT_RY, SIGMA_NU, SIG_MB
  REAL(8) :: XI
  INTEGER :: K, NLS, NP_L
  REAL(8), PARAMETER :: RY_EV = 13.6056923D0          ! 1 Ry in eV
  REAL(8), PARAMETER :: RY_HZ = 3.289842D15           ! 1 Ry in Hz
  REAL(8) :: EPHOTON_RY   ! current FREQ in Ry
  REAL(8) :: LOG_E_PHOTON ! log10(EPHOTON_RY), hoisted

  ! --- Photon energy in Ry for current FREQ (computed once) ---
  EPHOTON_RY   = FREQ / RY_HZ
  LOG_E_PHOTON = log10(EPHOTON_RY)

  ! --- Rebuild Boltzmann cache per depth if T has changed ---
  ! kT in Ry: kT[eV] = T[K] * 8.617e-5;  kT[Ry] = kT[eV] / 13.6057
  DO J = 1, NRHOX
    IF (abs(T(J) - TB_T_LAST(J)) .LT. 1.0D-6) CYCLE
    TJ = T(J)
    KT_RY = TJ * 8.617333D-5 / RY_EV
    IFLAT = 0
    DO IS = 1, N_TB_SPECIES
      DO IL = 1, TB_SPECIES(IS)%NLEVELS
        IFLAT = IFLAT + 1
        TB_BOLTZ(IFLAT, J) = dble(TB_SPECIES(IS)%LEVELS(IL)%G_STAT) &
          * exp(-TB_SPECIES(IS)%LEVELS(IL)%E_EXC_RY / KT_RY)
      END DO
    END DO
    TB_T_LAST(J) = TJ
  END DO

  ! Precompute STIM/RHO per depth (invariant over species/level/frequency).
  BLOCK
    REAL(8) :: STIM_OVER_RHO(kw)

    DO J = 1, NRHOX
      STIM_OVER_RHO(J) = STIM(J) / RHO(J)
    END DO

  ! --- Species-by-species accumulation ---
  !
  ! Loop structure:
  !   species -> [cheap reject if FREQ below species' lowest threshold]
  !     level -> [compute sigma(nu) ONCE]
  !       depth -> accumulate Boltzmann * sigma
  !
  ! Moving the log-grid interpolation out of the depth loop reduces that
  ! work by a factor of NRHOX.  The per-species quick-reject eliminates
  ! all work for species whose lowest threshold is above the current
  ! photon energy -- true for most species across visible and IR.
  IFLAT_SPECIES_BASE = 0
  DO IS = 1, N_TB_SPECIES
    NLS = TB_SPECIES(IS)%NLEVELS

    ! Species-level quick reject: photon below every retained threshold
    IF (EPHOTON_RY .LT. TB_SPECIES(IS)%LOWEST_THR_RY) THEN
      IFLAT_SPECIES_BASE = IFLAT_SPECIES_BASE + NLS
      CYCLE
    END IF

    DO IL = 1, NLS
      ! Per-level threshold test
      IFLAT = IFLAT_SPECIES_BASE + IL
      IF (EPHOTON_RY .LT. abs(TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY)) CYCLE

      ! Interpolate sigma(nu) on this level's native log grid - ONCE
      XI = (LOG_E_PHOTON - TB_SPECIES(IS)%LEVELS(IL)%LOG_E_MIN) &
           / TB_SPECIES(IS)%LEVELS(IL)%LOG_E_STRIDE
      K = int(XI) + 1
      IF (K .LT. 1) CYCLE
      NP_L = TB_SPECIES(IS)%LEVELS(IL)%NP
      IF (K .GE. NP_L) THEN
        ! Above grid: Kramers tail sigma ~ (nu_top/nu)^3
        SIG_MB = TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB(NP_L) &
          * (TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY(NP_L) / EPHOTON_RY)**3
      ELSE
        SIG_MB = TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB(K) &
          + (TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB(K+1) &
             - TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB(K)) * (XI - dble(K-1))
      END IF
      ! Convert Mb to cm^2
      SIGMA_NU = SIG_MB * 1.0D-18

      ! Accumulate over depths.  sigma(nu) is the same for all depths.
      DO J = 1, NRHOX
        ACONT_METAL(J) = ACONT_METAL(J) &
          + XNFP(J, TB_SPECIES(IS)%XNFP_INDEX) &
          * TB_BOLTZ(IFLAT, J) * SIGMA_NU * STIM_OVER_RHO(J)
      END DO
    END DO

    IFLAT_SPECIES_BASE = IFLAT_SPECIES_BASE + NLS
  END DO
  END BLOCK

  RETURN

END SUBROUTINE MBF_TOPBASE

!=========================================================================
! SUBROUTINE MBF_HIGH_ION
!
! High-ionization bound-free (C III/IV, N III-V, O III-VI, Ne III-VI)
! from filtered hotop.dat, plus Coulomb free-free for C/N/O/Ne/Mg/Si/S/Fe
! ions (stages I-V).  On first call, hotop.dat is read and edges whose
! species_id is already covered by TOPbase (IDs 22, 29, 37, 55, 56) are
! skipped, leaving ~40 of the original 60 edges in HI_EDGE_DATA.
!=========================================================================

SUBROUTINE MBF_HIGH_ION

  IMPLICIT NONE
  INTEGER, PARAMETER :: NPAR = 7
  INTEGER, PARAMETER :: N_TB_COV = 5
  INTEGER, PARAMETER :: TB_COVERED(N_TB_COV) = (/22, 29, 37, 55, 56/)

  REAL(8) :: AHOT_TMP(kw)
  REAL(8) :: FREE, XSECT, XX, FRATIO
  REAL(8) :: ROW(NPAR)
  INTEGER :: I, J, ID, IREAD, IKEEP
  CHARACTER(256) :: LINE
  LOGICAL :: SKIP

  !--- First-call initialization of filtered hotop.dat ---
  IF (.NOT. HIGH_ION_INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'hotop.dat', STATUS='OLD', ACTION='READ')
    DO
      READ(89, '(A)') LINE
      IF (LINE(1:1) .NE. '#') THEN
        BACKSPACE(89)
        EXIT
      END IF
    END DO
    ! First pass: count how many edges we keep
    IKEEP = 0
    DO IREAD = 1, 60
      READ(89, *) ROW
      ID = int(ROW(7))
      SKIP = .FALSE.
      DO I = 1, N_TB_COV
        IF (ID .EQ. TB_COVERED(I)) SKIP = .TRUE.
      END DO
      IF (.NOT. SKIP) IKEEP = IKEEP + 1
    END DO
    N_HIGH_ION_EDGES = IKEEP
    ALLOCATE(HI_EDGE_DATA(NPAR, N_HIGH_ION_EDGES))

    ! Second pass: actually store them
    REWIND(89)
    DO
      READ(89, '(A)') LINE
      IF (LINE(1:1) .NE. '#') THEN
        BACKSPACE(89)
        EXIT
      END IF
    END DO
    IKEEP = 0
    DO IREAD = 1, 60
      READ(89, *) ROW
      ID = int(ROW(7))
      SKIP = .FALSE.
      DO I = 1, N_TB_COV
        IF (ID .EQ. TB_COVERED(I)) SKIP = .TRUE.
      END DO
      IF (.NOT. SKIP) THEN
        IKEEP = IKEEP + 1
        HI_EDGE_DATA(:, IKEEP) = ROW
      END IF
    END DO
    CLOSE(89)
    HIGH_ION_INITIALIZED = .TRUE.
    IF (IDEBUG .EQ. 1) THEN
      WRITE(6,'(A,I3,A,I3,A)') ' MBF_HIGH_ION: kept ', N_HIGH_ION_EDGES, &
        ' of 60 hotop.dat edges after filtering TOPbase-covered species'
    END IF
  END IF

  !--- Coulomb free-free: C/N/O/Ne/Mg/Si/S/Fe ions (stages I-V) ---
  DO J = 1, NRHOX
    FREE = COULFF(J,1) * 1.0D0                                           &
            * (XNF(J,22) + XNF(J,29) + XNF(J,37) + XNF(J,56)             &
             + XNF(J,79) + XNF(J,106) + XNF(J,137) + XNF(J,352))         &
         + COULFF(J,2) * 4.0D0                                           &
            * (XNF(J,23) + XNF(J,30) + XNF(J,38) + XNF(J,57)             &
             + XNF(J,80) + XNF(J,107) + XNF(J,138) + XNF(J,353))         &
         + COULFF(J,3) * 9.0D0                                           &
            * (XNF(J,24) + XNF(J,31) + XNF(J,39) + XNF(J,58)             &
             + XNF(J,81) + XNF(J,108) + XNF(J,139) + XNF(J,354))         &
         + COULFF(J,4) * 16.0D0                                          &
            * (XNF(J,25) + XNF(J,32) + XNF(J,40) + XNF(J,59)             &
             + XNF(J,82) + XNF(J,109) + XNF(J,140) + XNF(J,355))         &
         + COULFF(J,5) * 25.0D0                                          &
            * (XNF(J,26) + XNF(J,33) + XNF(J,41) + XNF(J,60)             &
             + XNF(J,83) + XNF(J,110) + XNF(J,141) + XNF(J,356))
    AHOT_TMP(J) = FREE * COEFF_FF / (FREQ * FREQ * FREQ) * XNE(J) / sqrt(T(J))
  END DO

  !--- Filtered bound-free edges (modified Seaton, with "< 1%" skip) ---
  DO I = 1, N_HIGH_ION_EDGES
    IF (FREQ .LT. HI_EDGE_DATA(1, I)) CYCLE
    FRATIO = HI_EDGE_DATA(1, I) / FREQ
    XSECT = HI_EDGE_DATA(2, I) * (HI_EDGE_DATA(3, I) + FRATIO &
          - HI_EDGE_DATA(3, I) * FRATIO)                      &
          * sqrt(FRATIO ** int(HI_EDGE_DATA(4, I)))
    ID = int(HI_EDGE_DATA(7, I))
    DO J = 1, NRHOX
      XX = XSECT * XNFP(J, ID) * HI_EDGE_DATA(5, I)
      IF (XX .GT. AHOT_TMP(J) / 100.0D0) THEN
        AHOT_TMP(J) = AHOT_TMP(J) + XX / exp(HI_EDGE_DATA(6, I) / TKEV(J))
      END IF
    END DO
  END DO

  !--- Apply STIM/RHO and accumulate into ACONT_METAL ---
  ACONT_METAL = ACONT_METAL + AHOT_TMP * STIM / RHO

  RETURN

END SUBROUTINE MBF_HIGH_ION

!=========================================================================
! SUBROUTINE INIT_FELO
!
! Lazy loader for Fe I / Fe II bound-free OP data.  On first call, reads
! op_fe1.dat and op_fe2.dat from DATADIR into the FELO_* module arrays.
! Safe to call repeatedly; only reads on first invocation.
!=========================================================================

SUBROUTINE INIT_FELO

  IMPLICIT NONE
  INTEGER :: ION
  CHARACTER(len=16), PARAMETER :: FN(FELO_N_ION) = &
    (/ 'op_fe1.dat      ', 'op_fe2.dat      ' /)
  INTEGER :: I
  REAL(8), PARAMETER :: CLIGHT_CGS = 2.99792458D10
  REAL(8) :: NU_MAX

  IF (FELO_INITIALIZED) RETURN

  DO ION = 1, FELO_N_ION
    CALL FELO_LOAD_ION(ION, trim(DATADIR)//trim(FN(ION)))
  END DO

  ! Excitation energies in cm^-1, referenced to the most-bound level
  ! of each ion (nu_th_max).  E_exc(i) = (nu_th_max - nu_th(i)) / c.
  DO ION = 1, FELO_N_ION
    NU_MAX = maxval(FELO_NUTH(1:FELO_NLEV(ION), ION))
    DO I = 1, FELO_NLEV(ION)
      FELO_E(I, ION) = (NU_MAX - FELO_NUTH(I, ION)) / CLIGHT_CGS
    END DO
  END DO

  FELO_INITIALIZED = .TRUE.

  RETURN

END SUBROUTINE INIT_FELO

!=========================================================================
! SUBROUTINE FELO_LOAD_ION
!
! Read a single op_feN.dat file into the FELO_* module arrays at slot ION.
!
! File format (ASCII):
!   - Header: lines starting with '#' or blank (skipped).
!   - Body (pure numeric + quoted labels, no '#' comments allowed):
!       NLEV
!       ilev  g  nu_th  'label'             (x NLEV)
!       ilev  ifancy  nfit                  (x NLEV, each followed by:)
!         ifancy == 2:  s0  alpha  beta  gamma
!         ifancy > 100: nfit pairs of x, logsig (free-format, may span
!                       multiple lines)
!
! List-directed reads (READ(LUN, *) ...) handle whitespace-separated
! values spanning multiple lines, including the quoted label (Fortran
! strips the quotes when the receiving variable is CHARACTER).
!=========================================================================

SUBROUTINE FELO_LOAD_ION(ION, PATH)

  IMPLICIT NONE
  INTEGER,          INTENT(IN) :: ION
  CHARACTER(len=*), INTENT(IN) :: PATH

  INTEGER            :: LUN, IOS, NLEV, I, K, ILEV, IFANCY, NFIT
  REAL(8)            :: G_VAL, NU_TH
  CHARACTER(len=32)  :: LABEL
  CHARACTER(len=256) :: LINE

  LUN = 89
  OPEN(LUN, FILE=PATH, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS .NE. 0) THEN
    WRITE(6,'(A,A)') ' FELO_LOAD_ION ERROR: cannot open ', trim(PATH)
    STOP 1
  END IF

  ! Skip header comment/blank lines, then back up so list-directed
  ! reads can consume the first data line cleanly.
  DO
    READ(LUN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) THEN
      WRITE(6,'(A)') ' FELO_LOAD_ION ERROR: unexpected EOF in header'
      STOP 1
    END IF
    LINE = adjustl(LINE)
    IF (len_trim(LINE) .EQ. 0)   CYCLE
    IF (LINE(1:1)      .EQ. '#') CYCLE
    BACKSPACE(LUN)
    EXIT
  END DO

  ! NLEV
  READ(LUN, *) NLEV
  IF (NLEV .GT. FELO_NLEV_MAX) THEN
    WRITE(6,'(A,I4,A,I4)') ' FELO_LOAD_ION ERROR: NLEV=', NLEV, &
          ' exceeds FELO_NLEV_MAX=', FELO_NLEV_MAX
    STOP 1
  END IF
  FELO_NLEV(ION) = NLEV

  ! Level records: ilev, g, nu_th, 'label'.
  DO I = 1, NLEV
    READ(LUN, *) ILEV, G_VAL, NU_TH, LABEL
    FELO_G   (I, ION) = G_VAL
    FELO_NUTH(I, ION) = NU_TH
  END DO

  ! Bound-free records
  DO I = 1, NLEV
    READ(LUN, *) ILEV, IFANCY, NFIT
    FELO_IFANCY(I, ION) = IFANCY
    FELO_NPTS  (I, ION) = NFIT

    IF (IFANCY .EQ. 2) THEN
      READ(LUN, *) FELO_SEATON(1, I, ION), FELO_SEATON(2, I, ION), &
                   FELO_SEATON(3, I, ION), FELO_SEATON(4, I, ION)
    ELSE IF (IFANCY .GT. 100) THEN
      IF (NFIT .GT. FELO_NPTS_MAX) THEN
        WRITE(6,'(A,I4,A,I4)') ' FELO_LOAD_ION ERROR: NFIT=', NFIT, &
              ' exceeds FELO_NPTS_MAX=', FELO_NPTS_MAX
        STOP 1
      END IF
      READ(LUN, *) (FELO_X     (K, I, ION), &
                    FELO_LOGSIG(K, I, ION), K = 1, NFIT)
    ELSE
      WRITE(6,'(A,I4)') ' FELO_LOAD_ION ERROR: unexpected IFANCY=', IFANCY
      STOP 1
    END IF
  END DO

  CLOSE(LUN)

  RETURN

END SUBROUTINE FELO_LOAD_ION

!=========================================================================
! FUNCTION FELO_SIGMA
!
! Photoionization cross-section [cm^2] at wavenumber WAVE_CM for level I
! of ion stage ION (1 = Fe I, 2 = Fe II).  Returns zero below threshold
! and above the tabulated range.
!
! Format per TLUSTY User's Guide sec. 11.2:
!     x       = log10(nu / nu_threshold)
!     sig[Mb] = 10^logsig(x), linearly interpolated
! The Peach/Seaton analytic form (IFANCY=2) carries a * 1e-18 factor,
! so its s0 is in Mb as well.
!=========================================================================

FUNCTION FELO_SIGMA(ION, I, WAVE_CM) RESULT(SIG)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ION, I
  REAL(8), INTENT(IN) :: WAVE_CM        ! cm^-1
  REAL(8)             :: SIG

  REAL(8), PARAMETER :: CLIGHT_CGS = 2.99792458D10
  REAL(8), PARAMETER :: MB_TO_CM2  = 1.0D-18
  REAL(8) :: NU, NU_TH, X, XMIN, XMAX, LOGS, FRAC, RATIO
  REAL(8) :: S0, ALPHA, BETA, GAMMA
  INTEGER :: NFIT, K

  SIG   = 0.0D0
  NU    = WAVE_CM * CLIGHT_CGS
  NU_TH = FELO_NUTH(I, ION)

  ! Peach/Seaton analytic form (IFANCY = 2).
  IF (FELO_IFANCY(I, ION) .EQ. 2) THEN
    IF (NU .LT. NU_TH) RETURN
    S0    = FELO_SEATON(1, I, ION)
    ALPHA = FELO_SEATON(2, I, ION)
    BETA  = FELO_SEATON(3, I, ION)
    GAMMA = FELO_SEATON(4, I, ION)
    RATIO = NU_TH / NU
    SIG   = S0 * (BETA + (1.0D0 - BETA) * RATIO**GAMMA) * RATIO**ALPHA
    SIG   = SIG * MB_TO_CM2
    RETURN
  END IF

  ! OP fit-point form.
  IF (NU .LE. 0.0D0) RETURN
  X = log10(NU / NU_TH)

  NFIT = FELO_NPTS(I, ION)
  IF (NFIT .LT. 2) RETURN
  XMIN = FELO_X(1,    I, ION)
  XMAX = FELO_X(NFIT, I, ION)

  ! Zero below threshold and outside tabulated range -- we do not
  ! extrapolate, since TLUSTY's manual does not specify a convention.
  IF (X .LT. max(XMIN, 0.0D0)) RETURN
  IF (X .GT. XMAX)             RETURN

  ! Linear interpolation in x, in log10(sigma) space.  Scan is O(N) but
  ! NFIT is small (<= ~120) and each level is visited once per frequency
  ! call, so a binary search is not worth the extra code.
  DO K = 1, NFIT - 1
    IF (X .LE. FELO_X(K + 1, I, ION)) EXIT
  END DO
  FRAC = (X - FELO_X(K, I, ION)) / &
         (FELO_X(K + 1, I, ION) - FELO_X(K, I, ION))
  LOGS = FELO_LOGSIG(K, I, ION) + FRAC * &
         (FELO_LOGSIG(K + 1, I, ION) - FELO_LOGSIG(K, I, ION))
  SIG  = (10.0D0**LOGS) * MB_TO_CM2

END FUNCTION FELO_SIGMA

!=========================================================================
! SUBROUTINE FELO_OPACITY
!
! Fe I or Fe II bound-free photoionization opacity at the current
! wavenumber (module-level WAVENO).  Fills A(1:NRHOX) in cm^2/g,
! stimulated-emission and density corrected:
!
!   A(j) = STIM(j) * XNFP(j, IDX) / RHO(j)
!        * SUM_i [ g_i * exp(-E_i * HCKT(j)) * sigma_i(WAVENO) ]
!
! where IDX = FELO_XNFP_IDX(ION) points to the Fe I or Fe II column of
! XNFP.  XNFP carries n_species / U(T) -- the B&C partition function is
! already applied by the populations layer, so we use raw Boltzmann
! factors here.
!
! References:
!   Fe I  : Bautista (1997), A&AS 122, 167  (Iron Project XX)
!   Fe II : Nahar & Pradhan (1994), J. Phys. B 27, 429  (OP XIX)
! Source data from TLUSTY model atoms (Hubeny & Lanz 2017,
! arXiv:1706.01937), reformatted into op_fe1.dat / op_fe2.dat.
!
! ION = 1 selects Fe I; ION = 2 selects Fe II.  Level data are lazily
! loaded on first call via INIT_FELO.
!=========================================================================

SUBROUTINE FELO_OPACITY(ION, A)

  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ION
  REAL(8), INTENT(INOUT) :: A(kw)

  INTEGER :: I, J, NLEV, IDX
  REAL(8) :: SIG_I

  IF (.NOT. FELO_INITIALIZED) CALL INIT_FELO

  NLEV = FELO_NLEV(ION)
  IDX  = FELO_XNFP_IDX(ION)

  A(1:NRHOX) = 0.0D0

  DO I = 1, NLEV
    SIG_I = FELO_SIGMA(ION, I, WAVENO)
    IF (SIG_I .LE. 0.0D0) CYCLE
    DO J = 1, NRHOX
      A(J) = A(J) + SIG_I * FELO_G(I, ION) * &
                    EXP(-FELO_E(I, ION) * HCKT(J))
    END DO
  END DO

  DO J = 1, NRHOX
    A(J) = A(J) * STIM(J) * XNFP(J, IDX) / RHO(J)
  END DO

  RETURN

END SUBROUTINE FELO_OPACITY

!=========================================================================
! SUBROUTINE INIT_MBF_TOPBASE
!
! First-call setup for the TOPbase/Allende-Prieto metal bound-free.
!
! Actions:
!   1. Define the species table (filename, XNFP index, NIST IP).
!   2. For each species, read data/mbf/<sym><ion>.xs via READ_XS_FILE.
!   3. Compute NIST shift = IP_NIST_Ry - |E_bind_ground_raw|, apply to
!      every level's E_BIND and to every level's photon energy grid.
!   4. Compute E_EXC_i = |E_bind_ground_shifted| - |E_bind_i_shifted|
!      relative to the ground-state threshold.
!   5. Precompute LOG_E_MIN, LOG_E_STRIDE for fast log-grid lookup.
!   6. Decode G_STAT from ISLP.
!   7. Allocate TB_BOLTZ cache.
!   8. If IDEBUG==1, print a PF comparison vs Barklem & Collet.
!=========================================================================

SUBROUTINE INIT_MBF_TOPBASE

  IMPLICIT NONE

  ! --- Species table (30 entries) ---
  ! Hardcoded arrays: filenames, XNFP indices, NIST IPs in Ry
  INTEGER, PARAMETER :: NSPEC = 30
  CHARACTER(len=16), PARAMETER :: SPEC_LABEL(NSPEC) = (/ &
    'LiI  ', 'LiII ', 'BeI  ', 'BeII ', 'BI   ', 'BII  ',                &
    'CI   ', 'CII  ', 'NI   ', 'NII  ', 'OI   ', 'OII  ',                &
    'FI   ', 'FII  ', 'NeI  ', 'NeII ', 'NaI  ', 'NaII ',                &
    'MgI  ', 'MgII ', 'AlI  ', 'AlII ', 'SiI  ', 'SiII ',                &
    'SI   ', 'SII  ', 'ArI  ', 'ArII ', 'CaI  ', 'CaII ' /)
  CHARACTER(len=8), PARAMETER :: FN(NSPEC) = (/ &
    'li1.xs  ', 'li2.xs  ', 'be1.xs  ', 'be2.xs  ', 'b1.xs   ', 'b2.xs   ',  &
    'c1.xs   ', 'c2.xs   ', 'n1.xs   ', 'n2.xs   ', 'o1.xs   ', 'o2.xs   ',  &
    'f1.xs   ', 'f2.xs   ', 'ne1.xs  ', 'ne2.xs  ', 'na1.xs  ', 'na2.xs  ',  &
    'mg1.xs  ', 'mg2.xs  ', 'al1.xs  ', 'al2.xs  ', 'si1.xs  ', 'si2.xs  ',  &
    's1.xs   ', 's2.xs   ', 'ar1.xs  ', 'ar2.xs  ', 'ca1.xs  ', 'ca2.xs  ' /)
  INTEGER, PARAMETER :: XNFP_IDX(NSPEC) = (/ &
      6,   7,  10,  11,  15,  16,   &    ! Li I/II, Be I/II, B I/II
     21,  22,  28,  29,  36,  37,   &    ! C I/II, N I/II, O I/II
     45,  46,  55,  56,  66,  67,   &    ! F I/II, Ne I/II, Na I/II
     78,  79,  91,  92, 105, 106,   &    ! Mg I/II, Al I/II, Si I/II
    136, 137, 171, 172, 210, 211 /)      ! S I/II, Ar I/II, Ca I/II
  ! Z and ionization stage, aligned with SPEC_LABEL order (for PF diagnostic)
  INTEGER, PARAMETER :: SPEC_Z(NSPEC) = (/ &
      3,  3,   4,  4,   5,  5,          &    ! Li, Be, B
      6,  6,   7,  7,   8,  8,          &    ! C, N, O
      9,  9,  10, 10,  11, 11,          &    ! F, Ne, Na
     12, 12,  13, 13,  14, 14,          &    ! Mg, Al, Si
     16, 16,  18, 18,  20, 20 /)             ! S, Ar, Ca
  INTEGER, PARAMETER :: SPEC_ION(NSPEC) = (/ &
      1,  2,   1,  2,   1,  2,          &
      1,  2,   1,  2,   1,  2,          &
      1,  2,   1,  2,   1,  2,          &
      1,  2,   1,  2,   1,  2,          &
      1,  2,   1,  2,   1,  2 /)
  REAL(8), PARAMETER :: IP_EV(NSPEC) = (/ &
       5.3917D0, 75.6400D0,  9.3227D0, 18.2112D0,  8.2980D0, 25.1548D0,    &
      11.2603D0, 24.3833D0, 14.5341D0, 29.6013D0, 13.6181D0, 35.1211D0,    &
      17.4228D0, 34.9708D0, 21.5645D0, 40.9630D0,  5.1391D0, 47.2864D0,    &
       7.6462D0, 15.0353D0,  5.9858D0, 18.8286D0,  8.1517D0, 16.3459D0,    &
      10.3600D0, 23.3379D0, 15.7596D0, 27.6297D0,  6.1132D0, 11.8717D0 /)
  REAL(8), PARAMETER :: RY_EV = 13.6056923D0

  INTEGER :: IS, IL, K, IP_MIN
  REAL(8) :: E_BIND_MIN_RAW, IP_RY, SHIFT_RY
  REAL(8) :: TCHECK(3) = (/ 5000.0D0, 8000.0D0, 12000.0D0 /)
  REAL(8) :: U_TB, U_BC_VAL, KT_RY
  INTEGER :: IS_NLEVELS, ITC, IZ, IZ_IPION
  CHARACTER(len=256) :: PATH

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING INIT_MBF_TOPBASE'

  ! --- Read each species' .xs file ---
  N_TB_SPECIES = NSPEC
  DO IS = 1, NSPEC
    TB_SPECIES(IS)%LABEL = SPEC_LABEL(IS)
    TB_SPECIES(IS)%XNFP_INDEX = XNFP_IDX(IS)
    TB_SPECIES(IS)%E_IP_NIST_RY = IP_EV(IS) / RY_EV

    PATH = trim(DATADIR) // 'mbf/' // trim(FN(IS))
    CALL READ_XS_FILE(trim(PATH), TB_SPECIES(IS)%LEVELS, IS_NLEVELS)
    TB_SPECIES(IS)%NLEVELS = IS_NLEVELS

    ! Identify ground state: most negative E_BIND
    E_BIND_MIN_RAW = TB_SPECIES(IS)%LEVELS(1)%E_BIND_RY
    IP_MIN = 1
    DO IL = 2, IS_NLEVELS
      IF (TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY .LT. E_BIND_MIN_RAW) THEN
        E_BIND_MIN_RAW = TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY
        IP_MIN = IL
      END IF
    END DO

    ! NIST shift: we want the ground state to be at -IP_NIST exactly.
    ! raw ground E_BIND = E_BIND_MIN_RAW (negative).  Shift so that
    ! shifted ground = -IP_NIST_RY.  SHIFT = -IP_NIST_RY - E_BIND_MIN_RAW.
    ! Apply this shift (additive) to every level's E_BIND.
    ! We do NOT shift photon energies: the cross-section sigma(E_photon)
    ! tables are with respect to absolute photon energy, which doesn't
    ! move; only the level threshold moves.
    IP_RY = TB_SPECIES(IS)%E_IP_NIST_RY
    SHIFT_RY = -IP_RY - E_BIND_MIN_RAW
    TB_SPECIES(IS)%E_SHIFT_RY = SHIFT_RY

    DO IL = 1, IS_NLEVELS
      TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY = &
        TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY + SHIFT_RY
      ! Excitation energy: relative to ground = |E_BIND_ground| - |E_BIND_this|
      !                  = -E_BIND_ground - (-E_BIND_this)
      !                  = E_BIND_this - E_BIND_ground   (both negative)
      !                  = (E_BIND_this + SHIFT) - (E_BIND_ground + SHIFT)
      !                  = E_BIND_this_raw - E_BIND_ground_raw  (SHIFT cancels)
      ! After shift, E_BIND_ground = -IP_NIST_RY.
      TB_SPECIES(IS)%LEVELS(IL)%E_EXC_RY = &
        TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY - (-IP_RY)

      ! Decode G_STAT from ISLP: ISLP = 100*(2S+1) + 10*L + parity
      ! g_stat = (2S+1)*(2L+1)
      K = TB_SPECIES(IS)%LEVELS(IL)%ISLP
      TB_SPECIES(IS)%LEVELS(IL)%G_STAT = (K / 100) * (2 * mod(K / 10, 10) + 1)

      ! Precompute log10 grid parameters for fast interpolation
      TB_SPECIES(IS)%LEVELS(IL)%LOG_E_MIN = &
        log10(TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY(1))
      IF (TB_SPECIES(IS)%LEVELS(IL)%NP .GE. 2) THEN
        TB_SPECIES(IS)%LEVELS(IL)%LOG_E_STRIDE = &
          log10(TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY(2)) &
          - TB_SPECIES(IS)%LEVELS(IL)%LOG_E_MIN
      ELSE
        TB_SPECIES(IS)%LEVELS(IL)%LOG_E_STRIDE = 1.0D-3
      END IF
    END DO

    ! Level truncation: drop levels whose Boltzmann factor at T_MAX is
    ! below TRUNC_TOL times the ground-term contribution.  These levels
    ! cannot contribute to opacity at any temperature this code will
    ! encounter, so we free their memory and skip them at runtime.
    BLOCK
      REAL(8), PARAMETER :: T_MAX_TRUNC = 25000.0D0
      REAL(8), PARAMETER :: TRUNC_TOL   = 1.0D-4
      REAL(8) :: KT_MAX, REF_WEIGHT, WEIGHT, E_BIND_ABS_MIN
      INTEGER :: NL_OLD, NL_NEW
      TYPE(TOPBASE_LEVEL), ALLOCATABLE :: TMP_LEVELS(:)

      KT_MAX = T_MAX_TRUNC * 8.617333D-5 / RY_EV
      REF_WEIGHT = dble(TB_SPECIES(IS)%LEVELS(IP_MIN)%G_STAT)

      NL_OLD = TB_SPECIES(IS)%NLEVELS
      NL_NEW = 0
      ALLOCATE(TMP_LEVELS(NL_OLD))
      DO IL = 1, NL_OLD
        WEIGHT = dble(TB_SPECIES(IS)%LEVELS(IL)%G_STAT) &
               * exp(-TB_SPECIES(IS)%LEVELS(IL)%E_EXC_RY / KT_MAX)
        IF (WEIGHT .GE. TRUNC_TOL * REF_WEIGHT) THEN
          NL_NEW = NL_NEW + 1
          TMP_LEVELS(NL_NEW)%ISLP        = TB_SPECIES(IS)%LEVELS(IL)%ISLP
          TMP_LEVELS(NL_NEW)%ILV         = TB_SPECIES(IS)%LEVELS(IL)%ILV
          TMP_LEVELS(NL_NEW)%G_STAT      = TB_SPECIES(IS)%LEVELS(IL)%G_STAT
          TMP_LEVELS(NL_NEW)%NP          = TB_SPECIES(IS)%LEVELS(IL)%NP
          TMP_LEVELS(NL_NEW)%E_BIND_RY   = TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY
          TMP_LEVELS(NL_NEW)%E_EXC_RY    = TB_SPECIES(IS)%LEVELS(IL)%E_EXC_RY
          TMP_LEVELS(NL_NEW)%LOG_E_MIN   = TB_SPECIES(IS)%LEVELS(IL)%LOG_E_MIN
          TMP_LEVELS(NL_NEW)%LOG_E_STRIDE= TB_SPECIES(IS)%LEVELS(IL)%LOG_E_STRIDE
          CALL MOVE_ALLOC(TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY, &
                          TMP_LEVELS(NL_NEW)%E_PHOTON_RY)
          CALL MOVE_ALLOC(TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB,    &
                          TMP_LEVELS(NL_NEW)%SIGMA_MB)
        ELSE
          IF (ALLOCATED(TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY)) &
            DEALLOCATE(TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY)
          IF (ALLOCATED(TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB)) &
            DEALLOCATE(TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB)
        END IF
      END DO

      ! Swap packed array into place
      DEALLOCATE(TB_SPECIES(IS)%LEVELS)
      ALLOCATE(TB_SPECIES(IS)%LEVELS(NL_NEW))
      DO IL = 1, NL_NEW
        TB_SPECIES(IS)%LEVELS(IL)%ISLP        = TMP_LEVELS(IL)%ISLP
        TB_SPECIES(IS)%LEVELS(IL)%ILV         = TMP_LEVELS(IL)%ILV
        TB_SPECIES(IS)%LEVELS(IL)%G_STAT      = TMP_LEVELS(IL)%G_STAT
        TB_SPECIES(IS)%LEVELS(IL)%NP          = TMP_LEVELS(IL)%NP
        TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY   = TMP_LEVELS(IL)%E_BIND_RY
        TB_SPECIES(IS)%LEVELS(IL)%E_EXC_RY    = TMP_LEVELS(IL)%E_EXC_RY
        TB_SPECIES(IS)%LEVELS(IL)%LOG_E_MIN   = TMP_LEVELS(IL)%LOG_E_MIN
        TB_SPECIES(IS)%LEVELS(IL)%LOG_E_STRIDE= TMP_LEVELS(IL)%LOG_E_STRIDE
        CALL MOVE_ALLOC(TMP_LEVELS(IL)%E_PHOTON_RY, &
                        TB_SPECIES(IS)%LEVELS(IL)%E_PHOTON_RY)
        CALL MOVE_ALLOC(TMP_LEVELS(IL)%SIGMA_MB,    &
                        TB_SPECIES(IS)%LEVELS(IL)%SIGMA_MB)
      END DO
      DEALLOCATE(TMP_LEVELS)
      TB_SPECIES(IS)%NLEVELS = NL_NEW

      ! Record the smallest photon threshold among retained levels
      ! (used for a cheap per-species reject at runtime).
      E_BIND_ABS_MIN = abs(TB_SPECIES(IS)%LEVELS(1)%E_BIND_RY)
      DO IL = 2, NL_NEW
        IF (abs(TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY) .LT. E_BIND_ABS_MIN) &
          E_BIND_ABS_MIN = abs(TB_SPECIES(IS)%LEVELS(IL)%E_BIND_RY)
      END DO
      TB_SPECIES(IS)%LOWEST_THR_RY = E_BIND_ABS_MIN

      IF (IDEBUG .EQ. 1) THEN
        WRITE(6,'(A,A,A,I4,A,I4,A,F7.4,A)') ' INIT_MBF_TOPBASE: ',    &
          trim(TB_SPECIES(IS)%LABEL), ' truncated ', NL_OLD,            &
          ' -> ', NL_NEW, ' levels; lowest threshold = ',                &
          E_BIND_ABS_MIN * RY_EV, ' eV'
      END IF
    END BLOCK

    IF (IDEBUG .EQ. 1) THEN
      WRITE(6,'(A,A,I3,A,F8.4,A,F8.4,A)') ' INIT_MBF_TOPBASE: ',          &
        trim(TB_SPECIES(IS)%LABEL), TB_SPECIES(IS)%NLEVELS,                 &
        ' levels, IP_NIST=', IP_EV(IS), ' eV, shift=',                      &
        SHIFT_RY * RY_EV, ' eV'
    END IF
  END DO

  ! --- Allocate Boltzmann cache ---
  N_TB_LEVELS_TOTAL = 0
  DO IS = 1, N_TB_SPECIES
    N_TB_LEVELS_TOTAL = N_TB_LEVELS_TOTAL + TB_SPECIES(IS)%NLEVELS
  END DO
  ALLOCATE(TB_BOLTZ(N_TB_LEVELS_TOTAL, kw))
  ALLOCATE(TB_T_LAST(kw))
  ALLOCATE(TB_LEVEL_SPECIES(N_TB_LEVELS_TOTAL))
  TB_BOLTZ = 0.0D0
  TB_T_LAST = -1.0D0   ! force first-call rebuild
  K = 0
  DO IS = 1, N_TB_SPECIES
    DO IL = 1, TB_SPECIES(IS)%NLEVELS
      K = K + 1
      TB_LEVEL_SPECIES(K) = IS
    END DO
  END DO

  ! --- Diagnostic: compare TOPbase-derived U(T) vs Barklem & Collet ---
  IF (IDEBUG .EQ. 1) THEN
    WRITE(6,'(A)') ' INIT_MBF_TOPBASE: U(T) comparison (TOPbase / B&C):'
    WRITE(6,'(A)') &
      '   Species   Z  ION        T=5000        T=8000       T=12000'
    DO IS = 1, N_TB_SPECIES
      IZ = SPEC_Z(IS)
      IZ_IPION = SPEC_ION(IS)
      WRITE(6,'(A11,I4,I5)', advance='no') &
        '  '//trim(TB_SPECIES(IS)%LABEL), IZ, IZ_IPION
      DO ITC = 1, 3
        KT_RY = TCHECK(ITC) * 8.617333D-5 / RY_EV
        U_TB = 0.0D0
        DO IL = 1, TB_SPECIES(IS)%NLEVELS
          U_TB = U_TB + dble(TB_SPECIES(IS)%LEVELS(IL)%G_STAT) &
               * exp(-TB_SPECIES(IS)%LEVELS(IL)%E_EXC_RY / KT_RY)
        END DO
        IF (IZ .GT. 0 .AND. IZ_IPION .GT. 0 .AND. IZ_IPION .LE. 3) THEN
          U_BC_VAL = U_BC(IZ, IZ_IPION, TCHECK(ITC))
          IF (U_BC_VAL .GT. 0.0D0) THEN
            WRITE(6,'(F14.4)', advance='no') U_TB / U_BC_VAL
          ELSE
            WRITE(6,'(A14)', advance='no') '         --'
          END IF
        ELSE
          WRITE(6,'(A14)', advance='no') '         --'
        END IF
      END DO
      WRITE(6,'(A)') ''
    END DO
  END IF

  TB_INITIALIZED = .TRUE.

  RETURN

END SUBROUTINE INIT_MBF_TOPBASE

!=========================================================================
! SUBROUTINE READ_XS_FILE
!
! Parse one Allende-Prieto .xs file into an array of TOPBASE_LEVEL.
!
! File format:
!   5-line header (with '===' banner)
!   per level:
!     header line: I NZ NE ISLP ILV  E(RYD)  NP
!     NP data lines:  E_photon(Ry)  sigma(Mb)
!
! Returns LEVELS(1:NLEVELS_OUT) allocated.
!=========================================================================

SUBROUTINE READ_XS_FILE(FNAME, LEVELS, NLEVELS_OUT)

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN)  :: FNAME
  TYPE(TOPBASE_LEVEL), ALLOCATABLE, INTENT(OUT) :: LEVELS(:)
  INTEGER, INTENT(OUT) :: NLEVELS_OUT

  INTEGER, PARAMETER :: MAX_LEVELS = 500
  TYPE(TOPBASE_LEVEL) :: TMP(MAX_LEVELS)
  INTEGER :: LUN, IOS, II, NZ, NE, ISLP, ILV, NP, K
  REAL(8) :: EBIND
  CHARACTER(len=256) :: LINE
  LOGICAL :: AT_LEVEL_HEADER

  LUN = 89
  OPEN(UNIT=LUN, FILE=FNAME, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS .NE. 0) THEN
    WRITE(6,*) 'READ_XS_FILE: cannot open ', trim(FNAME)
    STOP 'READ_XS_FILE: missing file'
  END IF

  ! Skip header lines: read until first data-looking line, then backspace
  AT_LEVEL_HEADER = .FALSE.
  DO
    READ(LUN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) EXIT
    IF (LINE(1:1) .EQ. ' ' .AND. index(LINE, '=') .EQ. 0 .AND. &
        index(LINE, 'I') .EQ. 0) THEN
      BACKSPACE(LUN)
      EXIT
    END IF
  END DO

  ! Read level blocks until EOF
  NLEVELS_OUT = 0
  DO
    READ(LUN, *, IOSTAT=IOS) II, NZ, NE, ISLP, ILV, EBIND, NP
    IF (IOS .NE. 0) EXIT
    NLEVELS_OUT = NLEVELS_OUT + 1
    IF (NLEVELS_OUT .GT. MAX_LEVELS) THEN
      WRITE(6,*) 'READ_XS_FILE: level overflow in ', trim(FNAME)
      STOP
    END IF
    TMP(NLEVELS_OUT)%ISLP = ISLP
    TMP(NLEVELS_OUT)%ILV = ILV
    TMP(NLEVELS_OUT)%E_BIND_RY = EBIND
    TMP(NLEVELS_OUT)%NP = NP
    ALLOCATE(TMP(NLEVELS_OUT)%E_PHOTON_RY(NP))
    ALLOCATE(TMP(NLEVELS_OUT)%SIGMA_MB(NP))
    DO K = 1, NP
      READ(LUN, *) TMP(NLEVELS_OUT)%E_PHOTON_RY(K), TMP(NLEVELS_OUT)%SIGMA_MB(K)
    END DO
  END DO
  CLOSE(LUN)

  ! Pack into output array
  ALLOCATE(LEVELS(NLEVELS_OUT))
  DO II = 1, NLEVELS_OUT
    LEVELS(II)%ISLP = TMP(II)%ISLP
    LEVELS(II)%ILV = TMP(II)%ILV
    LEVELS(II)%E_BIND_RY = TMP(II)%E_BIND_RY
    LEVELS(II)%NP = TMP(II)%NP
    CALL MOVE_ALLOC(TMP(II)%E_PHOTON_RY, LEVELS(II)%E_PHOTON_RY)
    CALL MOVE_ALLOC(TMP(II)%SIGMA_MB,    LEVELS(II)%SIGMA_MB)
  END DO

  RETURN

END SUBROUTINE READ_XS_FILE


!=========================================================================
! SUBROUTINE C1OP
!
! C I bound-free photoionization opacity.
!
! Full level-by-level version matching atlas7lib.for.  Computes opacity
! from 14 levels of neutral carbon with explicit cross-section formulas,
! including Luo & Pradhan (1989) background fits and Burke & Taylor
! (1979) Fano resonance profiles for the ground-term levels.
!
! Structure (evaluated from lowest to highest threshold):
!
!   Levels 13-9 (3p 1S,1D,3P,3S,3D): X=0 (cross-sections not available)
!   Level 8  (3p 1P):  sigma = 2.1E-18 * (nu0/nu)^1.5
!   Level 6  (3s 1P):  sigma = 1.54E-18 * (nu0/nu)^1.2
!   Level 5  (3s 3P):  sigma = 0.2E-18 * (nu0/nu)^1.2
!   Level 14 (2s2p3 3P, 4P limit): X=0
!   Level 3  (2p2 1S):  Luo & Pradhan + 1 Fano resonance
!   Level 7  (2s2p3 3D, 4P limit): sigma = 16E-18 * (nu0/nu)^3
!   Level 2  (2p2 1D):  Luo & Pradhan + 2 Fano resonances
!   Level 1  (2p2 3P):  Luo & Pradhan, 3 fine-structure x 2 limits
!   Level 4  (2s2p3 5S, 4P limit): sigma = 1E-18 * (nu0/nu)^3
!
! Ionization limits:
!   C II 2P_0.5 = 90820.42 cm-1
!   C II 2P_1.5 = 90883.84 cm-1
!   C II 4P     = 133856.20 cm-1
!
! In LTE, departure coefficients = 1, so (BC1-BC2*EHVKT) = STIM
! and the source function S/H = BNU.
!
! References:
!   Luo, D. & Pradhan, A.K. 1989, J.Phys.B 22, 3377
!   Burke, P.G. & Taylor, K.T. 1979, J.Phys.B 12, 2971
!=========================================================================

SUBROUTINE C1OP

  IMPLICIT NONE

  REAL(8), PARAMETER :: RYD = 109732.298D0

  REAL(8)  :: H, S, A, B, X, EPS
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING C1OP'

  DO J = 1, NRHOX
    H = 1.0D-30
    S = 0.0D0

    levels: DO

      ! --- Level 13: 3p 1S (E=73975.91, g=1, threshold 16886.790) X=0 ---
      IF (WAVENO .LT. 16886.790D0) EXIT levels

      ! --- Level 12: 3p 1D (E=72610.72, g=5, threshold 18251.980) X=0 ---
      IF (WAVENO .LT. 18251.980D0) EXIT levels

      ! --- Level 11: 3p 3P (E=71374.90, g=9, threshold 19487.800) X=0 ---
      IF (WAVENO .LT. 19487.800D0) EXIT levels

      ! --- Level 10: 3p 3S (E=70743.95, g=3, threshold 20118.750) X=0 ---
      IF (WAVENO .LT. 20118.750D0) EXIT levels

      ! --- Level 9: 3p 3D (E=69722.00, g=15, threshold 21140.700) X=0 ---
      IF (WAVENO .LT. 21140.700D0) EXIT levels

      ! --- Level 8: 3p 1P (E=68856.33, g=3, threshold 22006.370) ---
      IF (WAVENO .LT. 22006.370D0) EXIT levels
      X = 2.1D-18 * (22006.370D0 / WAVENO)**1.5D0
      A = X * 3.0D0 * EXP(-68856.33D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 6: 3s 1P (E=61981.82, g=3, threshold 28880.880) ---
      IF (WAVENO .LT. 28880.880D0) EXIT levels
      X = 1.54D-18 * (28880.880D0 / WAVENO)**1.2D0
      A = X * 3.0D0 * EXP(-61981.82D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 5: 3s 3P (E=60373.00, g=9, threshold 30489.700) ---
      IF (WAVENO .LT. 30489.700D0) EXIT levels
      X = 0.2D-18 * (30489.700D0 / WAVENO)**1.2D0
      A = X * 9.0D0 * EXP(-60373.00D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 14: 2s2p3 3P (E=75254.93, g=9, threshold 58601.270) X=0 ---
      !     Ionizes to 4P limit at 133856.20 cm-1
      IF (WAVENO .LT. 58601.270D0) EXIT levels

      ! --- Level 3: 2p2 1S (E=21648.02, g=1) ---
      !     Ionizes to 2P_0.5 at 90820.42 -> threshold 69172.400
      !     Luo & Pradhan background + Burke & Taylor Fano resonance
      IF (WAVENO .LT. 69172.400D0) EXIT levels
      X = 10.0D0**(-16.80D0 - (WAVENO - 69172.400D0) / 3.0D0 / RYD)
      EPS = (WAVENO - 97700.0D0) * 2.0D0 / 2743.0D0
      A = 68.0D-18;  B = 118.0D-18
      X = X + (A * EPS + B) / (EPS**2 + 1.0D0)
      X = X / 3.0D0
      A = X * 1.0D0 * EXP(-21648.02D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     Also ionizes to 2P_1.5 at 90883.84 -> threshold 69235.820
      IF (WAVENO .LT. 69235.820D0) EXIT levels
      A = A * 2.0D0
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 7: 2s2p3 3D (E=64088.85, g=15, threshold 69767.350) ---
      !     Ionizes to 4P limit at 133856.20
      IF (WAVENO .LT. 69767.350D0) EXIT levels
      X = 16.0D-18 * (69767.350D0 / WAVENO)**3
      A = X * 15.0D0 * EXP(-64088.85D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 2: 2p2 1D (E=10192.66, g=5) ---
      !     Ionizes to 2P_0.5 at 90820.42 -> threshold 80627.760
      !     Luo & Pradhan background + two Burke & Taylor Fano resonances
      IF (WAVENO .LT. 80627.760D0) EXIT levels
      X = 10.0D0**(-16.80D0 - (WAVENO - 80627.760D0) / 3.0D0 / RYD)
      EPS = (WAVENO - 93917.0D0) * 2.0D0 / 9230.0D0
      A = 22.0D-18;  B = 26.0D-18
      X = X + (A * EPS + B) / (EPS**2 + 1.0D0)
      EPS = (WAVENO - 111130.0D0) * 2.0D0 / 2743.0D0
      A = -10.5D-18;  B = 46.0D-18
      X = X + (A * EPS + B) / (EPS**2 + 1.0D0)
      X = X / 3.0D0
      A = X * 5.0D0 * EXP(-10192.66D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     Also ionizes to 2P_1.5 -> threshold 80691.180
      IF (WAVENO .LT. 80691.180D0) EXIT levels
      A = A * 2.0D0
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 1: 2p2 3P (ground term, 3 fine-structure components) ---
      !     Ionizes to 2P_0.5 at 90820.42
      !     3P_2 (E=43.42, g=5) -> threshold 90777.000
      IF (WAVENO .LT. 90777.000D0) EXIT levels
      X = 10.0D0**(-16.80D0 - (WAVENO - 90777.000D0) / 3.0D0 / RYD)
      X = X / 3.0D0
      A = X * 5.0D0 * EXP(-43.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_1 (E=16.42, g=3) -> threshold 90804.000
      IF (WAVENO .LT. 90804.000D0) EXIT levels
      A = X * 3.0D0 * EXP(-16.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_0 (E=0.00, g=1) -> threshold 90820.420
      IF (WAVENO .LT. 90820.420D0) EXIT levels
      A = X * 1.0D0 * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     Ionizes to 2P_1.5 at 90883.84 (cross-section x 2)
      !     3P_2 -> threshold 90840.420
      IF (WAVENO .LT. 90840.420D0) EXIT levels
      X = X * 2.0D0
      A = X * 5.0D0 * EXP(-43.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_1 -> threshold 90867.420
      IF (WAVENO .LT. 90867.420D0) EXIT levels
      A = X * 3.0D0 * EXP(-16.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_0 -> threshold 90883.840
      IF (WAVENO .LT. 90883.840D0) EXIT levels
      A = X * 1.0D0 * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 4: 2s2p3 5S (E=33735.20, g=5, threshold 100121.000) ---
      !     Ionizes to 4P limit at 133856.20
      IF (WAVENO .LT. 100121.000D0) EXIT levels
      X = 1.0D-18 * (100121.000D0 / WAVENO)**3
      A = X * 5.0D0 * EXP(-33735.20D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      EXIT levels
    END DO levels

    ! Scale by C I population / rho
    H = H * XNFP(J, 21) / RHO(J)
    S = S * XNFP(J, 21) / RHO(J)

    AC1(J) = H
    IF (H .GT. 0.0D0) SC1(J) = S / H

  END DO

  RETURN

END SUBROUTINE C1OP

!=======================================================================
! SEATON
!=======================================================================

FUNCTION SEATON(FREQ0, XSECT, POWER, A)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: FREQ0, XSECT, POWER, A
  REAL(8) :: SEATON

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING SEATON'
  ! Seaton (1958) photoionization cross-section formula:
  ! sigma(nu) = sigma_0 * [A + (1-A)*(nu_0/nu)] * (nu_0/nu)^POWER
  SEATON = XSECT * (A + (1.0D0 - A) * (FREQ0 / FREQ)) &
         * sqrt((FREQ0 / FREQ)**int(2.0D0 * POWER + 0.01D0))
  RETURN

END FUNCTION SEATON

!=========================================================================
! SUBROUTINE MG1OP
!
! Mg I bound-free photoionization opacity.
!
! Computes cross-sections for 15 energy levels of neutral magnesium,
! all ionizing to the Mg II 3s ²S ground state (ELIM = 61671.02 cm⁻¹):
!
!   Levels 1–2:  3s4f ³F, ¹F — hydrogenic (XKARZAS n=4, l=3)
!   Levels 3–4:  3s4d ³D, ¹D — hydrogenic (XKARZAS n=4, l=2)
!   Level  5:    3s4p ¹P     — hydrogenic (XKARZAS n=4, l=1)
!   Level  6:    3s3d ³D     — power-law fit σ₀(ν₀/ν)^2.7
!   Level  7:    3s4p ³P     — power-law fit σ₀(ν₀/ν)^2.8
!   Level  8:    3s3d ¹D     — power-law fit σ₀(ν₀/ν)^2.7
!   Level  9:    3s4s ¹S     — power-law fit σ₀(ν₀/ν)^2.6
!   Level 10:    3s4s ³S     — power-law fit σ₀(ν₀/ν)^2.6
!   Level 11:    3s3p ¹P     — two-term power-law fit
!   Levels 12–14: 3s3p ³P₂,₁,₀ — power-law with floor
!   Level 15:    3s² ¹S (ground) — steep power-law (ν₀/ν)^10
!
! Final assembly sums over all levels with Boltzmann populations plus
! a hydrogenic high-n (n≥5) dissolved-level contribution.
!
! Temperature-dependent Boltzmann factors are cached and only
! recomputed when ITEMP changes.
!
! Note: Castelli index correction (25 Sep 2002) on ground-state
! threshold test (level 15, not 13).
!=========================================================================

SUBROUTINE MG1OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLEV = 15

  ! --- Atomic data: energy levels (cm⁻¹) and statistical weights ---
  !
  ! Levels 1–2:   3s4f  (³F, ¹F)
  ! Levels 3–4:   3s4d  (³D, ¹D)
  ! Level  5:     3s4p  (¹P)
  ! Level  6:     3s3d  (³D)
  ! Level  7:     3s4p  (³P)
  ! Level  8:     3s3d  (¹D)
  ! Levels 9–10:  3s4s  (¹S, ³S)
  ! Level 11:     3s3p  (¹P)
  ! Levels 12–14: 3s3p  (³P₂, ³P₁, ³P₀)
  ! Level 15:     3s²   (¹S ground)

  ! NOTE: ELEV(10)=41197.043 may be a digit transposition for 41197.403
  !   (NIST 3s4s ³S₁); power-law threshold 20473.617 = ELIM - 41197.403
  ! NOTE: ELEV(12)=21919.178 may be a typo for 21911.178
  !   (original comment says 21911.178; power-law threshold 39759.842 = ELIM - 21911.178)
  REAL(8), PARAMETER :: ELEV(NLEV) = (/ &
    54676.710D0, 54676.438D0, 54192.284D0, 53134.642D0, 49346.729D0, &
    47957.034D0, 47847.797D0, 46403.065D0, 43503.333D0, 41197.043D0, &  
                                                                        
    35051.264D0, 21919.178D0, 21870.464D0, 21850.405D0,     0.000D0 /)

  REAL(8), PARAMETER :: GLEV(NLEV) = (/ &
    21.D0, 7.D0, 15.D0, 5.D0, 3.D0, 15.D0, 9.D0, 5.D0, &
     1.D0, 3.D0,  3.D0, 5.D0, 3.D0,  1.D0, 1.D0 /)

  ! Quantum numbers (n,l) for XKARZAS — only levels 1–5 are hydrogenic
  INTEGER, PARAMETER :: NQ(5) = (/ 4, 4, 4, 4, 4 /)
  INTEGER, PARAMETER :: LQ(5) = (/ 3, 3, 2, 2, 1 /)

  ! Ionization limit (cm⁻¹): Mg II 3s ²S
  REAL(8), PARAMETER :: ELIM = 61671.02D0
  REAL(8), PARAMETER :: RYD  = 109732.298D0

  ! Threshold wavenumbers for power-law fits (ELIM - ELEV)
  REAL(8), PARAMETER :: THR6  = 13713.986D0   ! 3s3d ³D
  REAL(8), PARAMETER :: THR7  = 13823.223D0   ! 3s4p ³P
  REAL(8), PARAMETER :: THR8  = 15267.955D0   ! 3s3d ¹D
  REAL(8), PARAMETER :: THR9  = 18167.687D0   ! 3s4s ¹S
  REAL(8), PARAMETER :: THR10 = 20473.617D0   ! 3s4s ³S
  REAL(8), PARAMETER :: THR11 = 26619.756D0   ! 3s3p ¹P
  REAL(8), PARAMETER :: THR12 = 39759.842D0   ! 3s3p ³P

  ! --- Persistent state ---
  REAL(8),  SAVE :: BOLT(NLEV, kw)
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Local variables ---
  REAL(8)  :: X(NLEV)
  REAL(8)  :: ZEFF2, FREQ3, H, RATIO
  INTEGER :: I, J, K

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING MG1OP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors when temperature structure changes
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      DO I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      END DO
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Compute cross-sections at the current frequency
  !---------------------------------------------------------------------
  X(:) = 0.0D0

  levels: DO

    ! Levels 1–5: hydrogenic (XKARZAS)
    DO I = 1, 5
      IF (WAVENO .LT. ELIM - ELEV(I)) EXIT levels
      ZEFF2 = 16.0D0 / RYD * (ELIM - ELEV(I))
      X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
    END DO

    ! Level 6: 3s3d ³D — power-law fit
    IF (WAVENO .LT. ELIM - ELEV(6)) EXIT levels
    X(6) = 25.0D-18 * (THR6 / WAVENO)**2.7D0

    ! Level 7: 3s4p ³P — power-law fit
    IF (WAVENO .LT. ELIM - ELEV(7)) EXIT levels
    X(7) = 33.8D-18 * (THR7 / WAVENO)**2.8D0

    ! Level 8: 3s3d ¹D — power-law fit
    IF (WAVENO .LT. ELIM - ELEV(8)) EXIT levels
    X(8) = 45.0D-18 * (THR8 / WAVENO)**2.7D0

    ! Level 9: 3s4s ¹S — power-law fit
    IF (WAVENO .LT. ELIM - ELEV(9)) EXIT levels
    X(9) = 0.43D-18 * (THR9 / WAVENO)**2.6D0

    ! Level 10: 3s4s ³S — power-law fit
    IF (WAVENO .LT. ELIM - ELEV(10)) EXIT levels
    X(10) = 2.1D-18 * (THR10 / WAVENO)**2.6D0

    ! Level 11: 3s3p ¹P — two-term power-law
    IF (WAVENO .LT. ELIM - ELEV(11)) EXIT levels
    RATIO = THR11 / WAVENO
    X(11) = 16.0D-18 * RATIO**2.1D0 - 7.8D-18 * RATIO**9.5D0

    ! Levels 12–14: 3s3p ³P₂,₁,₀ — power-law with floor
    DO I = 12, 14
      IF (WAVENO .LT. ELIM - ELEV(I)) EXIT levels
      RATIO = THR12 / WAVENO
      X(I) = max(20.0D-18 * RATIO**2.7D0, 40.0D-18 * RATIO**14.0D0)
    END DO

    ! Level 15: 3s² ¹S ground state — steep power-law
    ! (Castelli index correction 25 Sep 2002: test on level 15, not 13)
    IF (WAVENO .LT. ELIM - ELEV(15)) EXIT levels
    X(15) = 1.1D-18 * ((ELIM - ELEV(15)) / WAVENO)**10.0D0

    EXIT levels
  END DO levels

  !---------------------------------------------------------------------
  ! Assemble opacity over depth: levels + dissolved high-n contribution
  !---------------------------------------------------------------------
  FREQ3 = 2.815D29 / (FREQ * FREQ * FREQ)

  DO J = 1, NRHOX
    ! High-n dissolved levels (n >= 5 to infinity), GFACTOR = 2
    H = FREQ3 * 2.0D0 * 2.0D0 / 2.0D0 / (RYD * HCKT(J)) &
      * (exp(-max(ELIM - RYD / 25.0D0, ELIM - WAVENO) * HCKT(J)) &
       - exp(-ELIM * HCKT(J)))

    DO I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    END DO

    AMG1(J) = H * XNFP(J, 78) * STIM(J) / RHO(J)
  END DO

  RETURN

END SUBROUTINE MG1OP

!==========================================================================
! SUBROUTINE AL1OP
!
! Al I bound-free photoionization opacity.
!
! Full multi-level version matching atlas7lib.for.  Computes opacity
! from 10 energy levels of neutral aluminum, all ionizing to the
! Al II 3s² ¹S ground state (ELIM = 48278.37 cm⁻¹):
!
!   Level 1 (BAL1 idx 9): 4f ²F   E=41319.377  g=14  threshold  6958.993
!   Level 2 (BAL1 idx 8): 5p ²P   E=40275.903  g= 6  threshold  8002.467
!   Level 3 (BAL1 idx 7): 4d ²D   E=38932.139  g=10  threshold  9346.231
!   Level 4 (BAL1 idx 6): 5s ²S   E=37689.413  g= 2  threshold 10588.957
!   Level 5 (BAL1 idx 5): 4p ²P   E=32960.363  g= 6  threshold 15318.007
!   Level 6 (BAL1 idx 4): 3d ²D   E=32436.241  g=10  threshold 15842.129
!   Level 7 (BAL1 idx 2): 4s ²S   E=25347.756  g= 2  threshold 22930.614
!   Level 8 (BAL1 idx 1): 3p ²P₃/₂ E=  112.061 g= 4  threshold 48166.309
!   Level 9 (BAL1 idx 1): 3p ²P₁/₂ E=    0.000 g= 2  threshold 48278.370
!   Level 10(BAL1 idx 3): P2 ⁴P   E=29097.110  g=12  threshold 55903.260
!
! In LTE (all departure coefficients = 1), the stimulated emission
! correction simplifies to STIM(J) = 1 - exp(-hν/kT).
!
! Note: Level 1 (4f ²F) has cross-section X=0 in the Kurucz data,
! so it contributes nothing.  It is retained for structural fidelity.
!==========================================================================

SUBROUTINE AL1OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLEV = 10

  ! Threshold wavenumbers (cm⁻¹) for each level
  REAL(8), PARAMETER :: THR(NLEV) = (/ &
     6958.993D0,  8002.467D0,  9346.231D0, 10588.957D0, 15318.007D0, &
    15842.129D0, 22930.614D0, 48166.309D0, 48278.370D0, 55903.260D0 /)

  ! Level energies (cm⁻¹)
  REAL(8), PARAMETER :: ELEV(NLEV) = (/ &
    41319.377D0, 40275.903D0, 38932.139D0, 37689.413D0, 32960.363D0, &
    32436.241D0, 25347.756D0,   112.061D0,     0.000D0, 29097.110D0 /)

  ! Statistical weights
  REAL(8), PARAMETER :: GLEV(NLEV) = (/ &
    14.0D0,  6.0D0, 10.0D0,  2.0D0,  6.0D0, &
    10.0D0,  2.0D0,  4.0D0,  2.0D0, 12.0D0 /)

  REAL(8)  :: H, S, A, X
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING AL1OP'

  DO J = 1, NRHOX
    H = 1.0D-30
    S = 0.0D0

    levels: DO
      ! Level 1: 4f ²F (X=0, no contribution — retained for fidelity)
      IF (WAVENO .LT. THR(1)) EXIT levels
      ! X = 0, so A = 0, skip

      ! Level 2: 5p ²P
      IF (WAVENO .LT. THR(2)) EXIT levels
      X = 50.0D-18 * (THR(2) / WAVENO)**3
      A = X * GLEV(2) * EXP(-ELEV(2) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 3: 4d ²D
      IF (WAVENO .LT. THR(3)) EXIT levels
      X = 50.0D-18 * (THR(3) / WAVENO)**3
      A = X * GLEV(3) * EXP(-ELEV(3) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 4: 5s ²S
      IF (WAVENO .LT. THR(4)) EXIT levels
      X = 56.7D-18 * (THR(4) / WAVENO)**1.9D0
      A = X * GLEV(4) * EXP(-ELEV(4) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 5: 4p ²P
      IF (WAVENO .LT. THR(5)) EXIT levels
      X = 14.5D-18 * THR(5) / WAVENO
      A = X * GLEV(5) * EXP(-ELEV(5) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 6: 3d ²D
      IF (WAVENO .LT. THR(6)) EXIT levels
      X = 47.0D-18 * (THR(6) / WAVENO)**1.83D0
      A = X * GLEV(6) * EXP(-ELEV(6) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 7: 4s ²S
      IF (WAVENO .LT. THR(7)) EXIT levels
      X = 10.0D-18 * (THR(7) / WAVENO)**2
      A = X * GLEV(7) * EXP(-ELEV(7) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 8: 3p ²P₃/₂ (ground fine-structure)
      IF (WAVENO .LT. THR(8)) EXIT levels
      X = 65.0D-18 * (THR(8) / WAVENO)**5
      A = X * GLEV(8) * EXP(-ELEV(8) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 9: 3p ²P₁/₂ (ground state) — same cross-section formula
      IF (WAVENO .LT. THR(9)) EXIT levels
      A = X * GLEV(9) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 10: P2 ⁴P (inner-shell)
      IF (WAVENO .LT. THR(10)) EXIT levels
      X = 10.0D-18 * (THR(10) / WAVENO)**2
      A = X * GLEV(10) * EXP(-ELEV(10) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      EXIT levels
    END DO levels

    IF (H .GT. 0.0D0) SAL1(J) = S / H
    AAL1(J) = H * XNFP(J, 91) / RHO(J)
  END DO

  RETURN

END SUBROUTINE AL1OP

!==========================================================================
! SUBROUTINE SI1OP
!
! Si I bound-free photoionization opacity.
!
! Computes cross-sections for 33 energy levels of neutral silicon,
! organized into four blocks by ionization limit of the Si II parent:
!
!   Block 1 (levels 1-22): Excited states (4s,4p,4d,4f,3d configurations)
!     ionizing to Si II 3s² 3p ²P average (ELIM = 65939.18 cm⁻¹).
!     Cross-sections from hydrogenic XKARZAS with effective charges.
!
!   Block 2 (levels 23-27): Low-lying 3p² states (¹S, ¹D, ³P₂, ³P₁, ³P₀)
!     ionizing to Si II 3s² 3p ²P₁/₂ (ELIM = 65747.55 cm⁻¹).
!     Cross-sections from Nahar & Pradhan (1993, J.Phys.B 26, 1109)
!     fits with Fano resonance profiles.
!
!   Block 3 (levels 23-27): Same low-lying states, contribution from
!     ionization to Si II 3s² 3p ²P₃/₂ (ELIM = 65747.55 + 287.45 cm⁻¹).
!     Added to Block 2 cross-sections with 2/3 statistical weight.
!
!   Block 4 (levels 28-33): Inner-shell 3s3p³ configurations ionizing
!     to Si II 3s 3p² ⁴P₁/₂ (ELIM = 65747.5 + 42824.35 cm⁻¹).
!     Cross-sections from XKARZAS with n=3, l=1.
!
! Final assembly sums over all levels with Boltzmann populations plus
! a hydrogenic high-n (n≥5) contribution near the series limit.
!
! Temperature-dependent Boltzmann factors cached when ITEMP changes.
!==========================================================================

SUBROUTINE SI1OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLEV = 33

  ! Energy levels (cm⁻¹) and statistical weights for 33 Si I states
  REAL(8), PARAMETER :: ELEV(33) = (/ &
    59962.284D0, 59100.000D0, 59077.112D0, 58893.40D0, 58801.529D0, &  ! 4d³P, 4f, 4d³D, 4d¹F, 4d¹P
    58777.000D0, 57488.974D0, 56503.346D0, 54225.621D0, 53387.34D0, &  ! 4f, 4d³F, 4d¹D, 3d³D, 3d¹P
    53362.24D0,  51612.012D0, 50533.424D0, 50189.389D0, 49965.894D0, & ! 3d¹F, 4p¹S, 3d³P, 4p¹D, 3d³F
    49399.670D0, 49128.131D0, 48161.459D0, 47351.554D0, 47284.061D0, & ! 4p³S, 4p³P, 4p³D, 3d¹D, 4p¹P
    40991.884D0, 39859.920D0, 15394.370D0,  6298.850D0,   223.157D0, & ! 4s¹P, 4s³P, 3p²¹S, 3p²¹D, 3p²³P₂
    77.115D0,        0.000D0, 94000.000D0, 79664.000D0, 72000.000D0, & ! ³P₁, ³P₀, 3s3p³¹P, ³S, ¹D
    56698.738D0, 45303.310D0, 33326.053D0 /)                           ! ³P, ³D, ⁵S

  REAL(8), PARAMETER :: GLEV(33) = (/ &
    9.D0, 56.D0, 15.D0, 7.D0, 3.D0, 28.D0, 21.D0, 5.D0, 15.D0, &
    3.D0, 7.D0, 1.D0, 9.D0, 5.D0, 21.D0, 3.D0, 9.D0, 15.D0, &
    5.D0, 3.D0, 3.D0, 9.D0, 1.D0, 5.D0, 5.D0, 3.D0, 1.D0, &
    3.D0, 3.D0, 5.D0, 12.D0, 15.D0, 5.D0 /)

  ! L quantum numbers for XKARZAS calls (levels 1-22)
  ! d-orbital → l=2, f-orbital → l=3, p-orbital → l=1, s-orbital → l=0
  INTEGER, PARAMETER :: LQNUM(22) = (/ &
    2, 3, 2, 2, 2, 3, 2, 2,  &  ! levels 1-8 (n=4 shell)
    2, 2, 2, 1, 2, 1, 2, 1,  &  ! levels 9-16 (n=3d and n=4p)
    1, 1, 2, 1, 0, 0 /)          ! levels 17-22 (n=4p, 3d, 4s)

  ! Principal quantum numbers for XKARZAS calls (levels 1-22)
  INTEGER, PARAMETER :: NQNUM(22) = (/ &
    4, 4, 4, 4, 4, 4, 4, 4,  &  ! levels 1-8
    3, 3, 3, 4, 3, 4, 3, 4,  &  ! levels 9-16
    4, 4, 3, 4, 4, 4 /)          ! levels 17-22

  ! Z_eff² prefactor: n²/RYD for each level (levels 1-22)
  ! n=3 levels use 9/RYD, n=4 levels use 16/RYD
  REAL(8), PARAMETER :: RYD = 109732.298D0

  REAL(8), SAVE :: BOLT(NLEV, kw)
  INTEGER, SAVE :: ITEMP1 = 0

  REAL(8)  :: X(NLEV), FREQ3, ZEFF2, EPS, RESON1
  REAL(8)  :: H, GFACTOR, DEGEN, Z, ELIM_BLK
  INTEGER :: I, J, K

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING SI1OP'

  ! Recompute Boltzmann factors when temperature changes
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      DO I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      END DO
    END DO
  END IF

  Z = 1.0D0
  FREQ3 = 2.815D29 / FREQ / FREQ / FREQ * Z**4
  WAVENO = FREQ / CLIGHT

  DO I = 1, NLEV
    X(I) = 0.0D0
  END DO

  ! =====================================================================
  ! Block 1: Excited states → Si II 3s² 3p ²P average (levels 1-22)
  !   Hydrogenic cross-sections via XKARZAS
  ! =====================================================================
  ELIM_BLK = 65939.18D0

  DO I = 1, 22
    IF (WAVENO .LT. ELIM_BLK - ELEV(I)) EXIT
    ZEFF2 = dble(NQNUM(I))**2 / RYD * (ELIM_BLK - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQNUM(I), LQNUM(I))
  END DO

  ! =====================================================================
  ! Block 2: Low-lying 3p² states → Si II 3s² 3p ²P₁/₂ (levels 23-27)
  !   Nahar & Pradhan (1993) fits with Fano resonance profiles
  !   Statistical weight factor: 1/3 for ²P₁/₂ (g=2 of g_total=6)
  ! =====================================================================
  ELIM_BLK = 65747.55D0

  !  3s² 3p² ¹S (level 23)
  IF (WAVENO .GE. ELIM_BLK - ELEV(23)) THEN
    EPS = (WAVENO - 70000.0D0) * 2.0D0 / 6500.0D0
    ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
    RESON1 = (97.D-18 * EPS + 94.D-18) / (EPS**2 + 1.0D0)
!     EPS=(WAVENO-89700.)*2./75.
!     RESON2=900.E-18/(EPS**2+1.)
    X(23) = (37.D-18 * (50353.180D0 / WAVENO)**2.40D0 + RESON1) / 3.0D0
!     X(23)=46.E-18*(50353.180/WAVENO)**.5/3.

    !  3s² 3p² ¹D (level 24)
    IF (WAVENO .GE. ELIM_BLK - ELEV(24)) THEN
      ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     IF(WAVENO.LT.59448.700)GO TO 30
      EPS = (WAVENO - 78600.0D0) * 2.0D0 / 13000.0D0
      RESON1 = (-10.D-18 * EPS + 77.D-18) / (EPS**2 + 1.0D0)
!     EPS=(WAVENO-98000.)*2./400.
!     RESON2=(65.E-18*EPS+65.E-18)/(EPS**2+1.)
!     EPS=(WAVENO-99500.)*2./50.
!     RESON3=204.E-18/(EPS**2+1.)
      X(24) = (24.5D-18 * (59448.700D0 / WAVENO)**1.85D0 + RESON1) / 3.0D0
!     X(24)=35.E-18*(59448.700/WAVENO)**3/3.

      !  3s² 3p² ³P₂ (level 25)
      IF (WAVENO .GE. ELIM_BLK - ELEV(25)) THEN
        ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     RESON1=0.
!     RESON2=0.
!     RESON3=0.
!     IF(WAVENO.LT.100000.)THEN
!     EPS=(WAVENO-92200.)*2./950.
!     RESON1=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.92200.AND.WAVENO.LT.103600.)THEN
!     EPS=(WAVENO-100000.)*2./300.
!     RESON2=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.100000.AND.WAVENO.LT.105000.)THEN
!     EPS=(WAVENO-103600.)*2./150.
!     RESON3=20.E-18*EPS/(EPS**2+1.)
!     ENDIF
        IF (WAVENO .LE. 74000.0D0) THEN
          X(25) = 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 / 3.0D0
        ELSE
          X(25) = 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 / 3.0D0
        END IF
!     X(25)=X(25)+(RESON1+RESON2+RESON3)/3.
!     X(25)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)/3.

        !  3s² 3p² ³P₁ (level 26)
        IF (WAVENO .GE. ELIM_BLK - ELEV(26)) THEN
          ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     RESON1=0.
!     RESON2=0.
!     RESON3=0.
!     IF(WAVENO.LT.100000.)THEN
!     EPS=(WAVENO-92200.)*2./950.
!     RESON1=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.92200.AND.WAVENO.LT.103600.)THEN
!     EPS=(WAVENO-100000.)*2./300.
!     RESON2=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.100000.AND.WAVENO.LT.105000.)THEN
!     EPS=(WAVENO-103600.)*2./150.
!     RESON3=20.E-18*EPS/(EPS**2+1.)
!     ENDIF
          IF (WAVENO .LE. 74000.0D0) THEN
            X(26) = 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
          ELSE
            X(26) = 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
          END IF
!     X(26)=X(26)+(RESON1+RESON2+RESON3)*2./3.
!     X(26)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.

          !  3s² 3p² ³P₀ (level 27)
          IF (WAVENO .GE. ELIM_BLK - ELEV(27)) THEN
            ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     RESON1=0.
!     RESON2=0.
!     RESON3=0.
!     IF(WAVENO.LT.100000.)THEN
!     EPS=(WAVENO-92200.)*2./950.
!     RESON1=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.92200.AND.WAVENO.LT.103600.)THEN
!     EPS=(WAVENO-100000.)*2./300.
!     RESON2=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.100000.AND.WAVENO.LT.105000.)THEN
!     EPS=(WAVENO-103600.)*2./150.
!     RESON3=20.E-18*EPS/(EPS**2+1.)
!     ENDIF
            IF (WAVENO .LE. 74000.0D0) THEN
              X(27) = 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 / 3.0D0
            ELSE
              X(27) = 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 / 3.0D0
            END IF
!     X(27)=X(27)+(RESON1+RESON2+RESON3)/3.
!     X(27)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)/3.
          END IF  ! level 27
        END IF  ! level 26
      END IF  ! level 25
    END IF  ! level 24
  END IF  ! level 23

  ! =====================================================================
  ! Block 3: Same low-lying 3p² states → Si II 3s² 3p ²P₃/₂ (levels 23-27)
  !   Adds ²P₃/₂ contribution (statistical weight 2/3) to X(23)-X(27)
  ! =====================================================================
  ELIM_BLK = 65747.55D0 + 287.45D0

  !  3s² 3p² ¹S (level 23)
  IF (WAVENO .GE. ELIM_BLK - ELEV(23)) THEN
    EPS = (WAVENO - 70000.0D0) * 2.0D0 / 6500.0D0
    ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
    RESON1 = (97.D-18 * EPS + 94.D-18) / (EPS**2 + 1.0D0)
!     EPS=(WAVENO-89700.)*2./75.
!     RESON2=900.E-18/(EPS**2+1.)
    X(23) = X(23) + (37.D-18 * (50353.180D0 / WAVENO)**2.40D0 + RESON1) * 2.0D0 / 3.0D0
!     X(23)=X(23)+46.E-18*(50353.180/WAVENO)**.5*2./3.

    !  3s² 3p² ¹D (level 24)
    IF (WAVENO .GE. ELIM_BLK - ELEV(24)) THEN
      ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!      IF(WAVENO.LT.59448.700)GO TO 30
      EPS = (WAVENO - 78600.0D0) * 2.0D0 / 13000.0D0
      RESON1 = (-10.D-18 * EPS + 77.D-18) / (EPS**2 + 1.0D0)
!     EPS=(WAVENO-98000.)*2./400.
!     RESON2=(65.E-18*EPS+65.E-18)/(EPS**2+1.)
!     EPS=(WAVENO-99500.)*2./50.
!     RESON3=204.E-18/(EPS**2+1.)
      X(24) = X(24) + (24.5D-18 * (59448.700D0 / WAVENO)**1.85D0 + RESON1) * 2.0D0 / 3.0D0
!     X(24)=X(24)+35.E-18*(59448.700/WAVENO)**3*2./3.

      !  3s² 3p² ³P₂ (level 25)
      IF (WAVENO .GE. ELIM_BLK - ELEV(25)) THEN
        ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     RESON1=0.
!     RESON2=0.
!     RESON3=0.
!     IF(WAVENO.LT.100000.)THEN
!     EPS=(WAVENO-92200.)*2./950.
!     RESON1=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.92200.AND.WAVENO.LT.103600.)THEN
!     EPS=(WAVENO-100000.)*2./300.
!     RESON2=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.100000.AND.WAVENO.LT.105000.)THEN
!     EPS=(WAVENO-103600.)*2./150.
!     RESON3=20.E-18*EPS/(EPS**2+1.)
!     ENDIF
        IF (WAVENO .LE. 74000.0D0) THEN
          X(25) = X(25) + 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
        ELSE
          X(25) = X(25) + 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
        END IF
!     X(25)=X(25)+(RESON1+RESON2+RESON3)*2./3.
!     X(25)=X(25)+37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.

        !  3s² 3p² ³P₁ (level 26)
        IF (WAVENO .GE. ELIM_BLK - ELEV(26)) THEN
          ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     RESON1=0.
!     RESON2=0.
!     RESON3=0.
!     IF(WAVENO.LT.100000.)THEN
!     EPS=(WAVENO-92200.)*2./950.
!     RESON1=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.92200.AND.WAVENO.LT.103600.)THEN
!     EPS=(WAVENO-100000.)*2./300.
!     RESON2=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.100000.AND.WAVENO.LT.105000.)THEN
!     EPS=(WAVENO-103600.)*2./150.
!     RESON3=20.E-18*EPS/(EPS**2+1.)
!     ENDIF
          IF (WAVENO .LE. 74000.0D0) THEN
            X(26) = X(26) + 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
          ELSE
            X(26) = X(26) + 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
          END IF
!     X(26)=X(26)+(RESON1+RESON2+RESON3)*2./3.
!     X(26)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.

          !  3s² 3p² ³P₀ (level 27)
          IF (WAVENO .GE. ELIM_BLK - ELEV(27)) THEN
            ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
!     RESON1=0.
!     RESON2=0.
!     RESON3=0.
!     IF(WAVENO.LT.100000.)THEN
!     EPS=(WAVENO-92200.)*2./950.
!     RESON1=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.92200.AND.WAVENO.LT.103600.)THEN
!     EPS=(WAVENO-100000.)*2./300.
!     RESON2=30.E-18*EPS/(EPS**2+1.)
!     ENDIF
!     IF(WAVENO.GT.100000.AND.WAVENO.LT.105000.)THEN
!     EPS=(WAVENO-103600.)*2./150.
!     RESON3=20.E-18*EPS/(EPS**2+1.)
!     ENDIF
            IF (WAVENO .LE. 74000.0D0) THEN
              X(27) = X(27) + 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
            ELSE
              X(27) = X(27) + 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
            END IF
!     X(27)=X(27)+(RESON1+RESON2+RESON3)*2./3.
!     X(27)=X(27)+37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.
          END IF  ! level 27
        END IF  ! level 26
      END IF  ! level 25
    END IF  ! level 24
  END IF  ! level 23

  ! =====================================================================
  ! Block 4: Inner-shell states → Si II 3s 3p² ⁴P₁/₂ (levels 28-33)
  !   3s3p³ configurations, cross-sections from XKARZAS with n=3, l=1
  ! =====================================================================
  ELIM_BLK = 65747.5D0 + 42824.35D0

  ! NOTE: DEGEN=3 is set once and reused for levels 28-33, but
  ! GLEV varies (3,3,5,12,15,5). This appears to be a possible bug
  ! in the original Kurucz code where DEGEN should perhaps equal
  ! GLEV for each level. The original behavior is preserved here.
  DEGEN = 3.0D0

  DO I = 28, NLEV
    IF (WAVENO .LT. ELIM_BLK - ELEV(I)) EXIT
    ZEFF2 = 9.0D0 / RYD * (ELIM_BLK - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, 3, 1) * DEGEN
  END DO

  ! =====================================================================
  ! Final assembly: sum over all levels with Boltzmann populations
  ! =====================================================================
  ELIM_BLK = 65747.55D0
  GFACTOR = 6.0D0

  DO J = 1, NRHOX
    ! High-n hydrogenic contribution (n=5 to infinity)
    H = FREQ3 * GFACTOR * 2.0D0 / 2.0D0 / (RYD * Z**2 * HCKT(J)) &
      * (exp(-max(ELIM_BLK - RYD * Z**2 / 5.0D0**2, ELIM_BLK - WAVENO) * HCKT(J)) &
      - exp(-ELIM_BLK * HCKT(J)))
    ! Sum resolved levels
    DO I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    END DO
    ASI1(J) = H * XNFP(J, 105) * STIM(J) / RHO(J)
  END DO
  RETURN

END SUBROUTINE SI1OP

!==========================================================================
! SUBROUTINE FE1OP
!
! Fe I bound-free photoionization opacity.
!
! Full multi-level version matching atlas7lib.for.  Sums opacity from
! 48 energy levels of neutral iron using Fano-profile cross-sections:
!
!   σ(ν) = 3×10⁻¹⁸ / [1 + ((ν₀+3000-ν) / (ν₀×0.1))⁴]
!
! where ν₀ = WNO(i) is the threshold wavenumber for each level.
! Only contributes above 21000 cm⁻¹ (λ < 4762 Å).
!
! The source function SFE1 uses a fudge factor from the Si I ground
! state population (BSI1(J,1) in F77).  In LTE, SFE1 = BNU.
!==========================================================================

SUBROUTINE FE1OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLEV = 48

  REAL(8), PARAMETER :: G(48) = (/ &
    25.D0, 35.D0, 21.D0, 15.D0,  9.D0, 35.D0, 33.D0, 21.D0, &
    27.D0, 49.D0,  9.D0, 21.D0, 27.D0,  9.D0,  9.D0, 25.D0, &
    33.D0, 15.D0, 35.D0,  3.D0,  5.D0, 11.D0, 15.D0, 13.D0, &
    15.D0,  9.D0, 21.D0, 15.D0, 21.D0, 25.D0, 35.D0,  9.D0, &
     5.D0, 45.D0, 27.D0, 21.D0, 15.D0, 21.D0, 15.D0, 25.D0, &
    21.D0, 35.D0,  5.D0, 15.D0, 45.D0, 35.D0, 55.D0, 25.D0 /)

  REAL(8), PARAMETER :: E(48) = (/ &
      500.D0,  7500.D0, 12500.D0, 17500.D0, 19000.D0, 19500.D0, &
    19500.D0, 21000.D0, 22000.D0, 23000.D0, 23000.D0, 24000.D0, &
    24000.D0, 24500.D0, 24500.D0, 26000.D0, 26500.D0, 26500.D0, &
    27000.D0, 27500.D0, 28500.D0, 29000.D0, 29500.D0, 29500.D0, &
    29500.D0, 30000.D0, 31500.D0, 31500.D0, 33500.D0, 33500.D0, &
    34000.D0, 34500.D0, 34500.D0, 35000.D0, 35500.D0, 37000.D0, &
    37000.D0, 37000.D0, 38500.D0, 40000.D0, 40000.D0, 41000.D0, &
    41000.D0, 43000.D0, 43000.D0, 43000.D0, 43000.D0, 44000.D0 /)

  REAL(8), PARAMETER :: WNO(48) = (/ &
    63500.D0, 58500.D0, 53500.D0, 59500.D0, 45000.D0, 44500.D0, &
    44500.D0, 43000.D0, 58000.D0, 41000.D0, 54000.D0, 40000.D0, &
    40000.D0, 57500.D0, 55500.D0, 38000.D0, 57500.D0, 57500.D0, &
    37000.D0, 54500.D0, 53500.D0, 55000.D0, 34500.D0, 34500.D0, &
    34500.D0, 34000.D0, 32500.D0, 32500.D0, 32500.D0, 32500.D0, &
    32000.D0, 29500.D0, 29500.D0, 31000.D0, 30500.D0, 29000.D0, &
    27000.D0, 54000.D0, 27500.D0, 24000.D0, 47000.D0, 23000.D0, &
    44000.D0, 42000.D0, 42000.D0, 21000.D0, 42000.D0, 42000.D0 /)

  REAL(8)  :: XSECT
  INTEGER :: I, J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING FE1OP'

  ! Zero output
  AFE1 = 0.0D0

  ! No contribution below 21000 cm⁻¹
  IF (WAVENO .LT. 21000.0D0) RETURN

  ! Sum over all 48 levels
  DO I = 1, NLEV
    IF (WNO(I) .GT. WAVENO) CYCLE
    XSECT = 3.0D-18 / (1.0D0 + ((WNO(I) + 3000.0D0 - WAVENO) &
          / WNO(I) / 0.1D0)**4)
    AFE1 = AFE1 + XSECT * G(I) * EXP(-E(I) * HCKT)
  END DO

  ! Apply population, stimulated emission, and density factors
  AFE1 = AFE1 * STIM * XNFP(:, 351) / RHO
  SFE1 = BNU

  RETURN

END SUBROUTINE FE1OP

!=========================================================================
! FUNCTION CHOP(J)
!
! CH molecular bound-free + bound-bound opacity.
!
! Returns the CH contribution to the absorption coefficient at depth
! point J for the current frequency.  The cross-section is obtained by
! bilinear interpolation in a precomputed table CROSSCH(T, E):
!
!   Energy axis:  E = 2.0–10.5 eV in 0.1 eV steps  (indices 20–105)
!   Temperature axis:  T = 2000–9000 K in 500 K steps  (15 points)
!
! The table values are log10(σ·Q) where Q is the CH partition function.
! A separate partition function table PARTCH covers T = 1000–9000 K in
! 200 K steps (41 points) and is interpolated independently.
!
! The final opacity is:  CHOP = 10^(σ·Q interpolated) × Q(T)
!
! Frequency-dependent quantities (energy interpolation of CROSSCHT) are
! cached and only recomputed when FREQ changes between calls.
!
! Returns zero for T ≥ 9000 K or photon energy outside 2.0–10.5 eV.
!
! NOTE: The original code recomputes WAVENO = FREQ/c here, which is
! redundant with the module-level WAVENO already set by the caller.
! We use the module WAVENO directly instead.
!=========================================================================

FUNCTION CHOP(J)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: J
  REAL(8)  :: CHOP

  ! --- Cross-section × partition function table ---
  ! CROSSCH(IT, N): log10(σ·Q) for temperature index IT and energy index N
  !   IT = 1–15:  T = 2000, 2500, ..., 9000 K
  !   N  = 1–105: E = 0.1, 0.2, ..., 10.5 eV  (only N=20–105 used)
  ! Loaded from crossch.dat on first call.
  REAL(8),    SAVE :: CROSSCH(15,105)
  LOGICAL,   SAVE :: CROSSCH_LOADED = .FALSE.
  INTEGER :: II, JJ

  ! CH partition function table: T = 1000 to 9000 K in 200 K steps (41 points)
  REAL(8), PARAMETER :: PARTCH(41) = (/ &
    203.741D0,  249.643D0,  299.341D0,  353.477D0,  412.607D0,  477.237D0,  547.817D0, &
    624.786D0,  708.543D0,  799.463D0,  897.912D0, 1004.227D0, 1118.738D0, 1241.761D0, &
   1373.588D0, 1514.481D0, 1664.677D0, 1824.394D0, 1993.801D0, 2173.050D0, 2362.234D0, &
   2561.424D0, 2770.674D0, 2989.930D0, 3219.204D0, 3458.378D0, 3707.355D0, 3966.005D0, &
   4234.155D0, 4511.604D0, 4798.135D0, 5093.554D0, 5397.593D0, 5709.948D0, 6030.401D0, &
   6358.646D0, 6694.379D0, 7037.313D0, 7387.147D0, 7743.579D0, 8106.313D0 /)

  ! --- Cached frequency-interpolated cross-sections ---
  REAL(8),  SAVE :: CROSSCHT(15)     ! cross-section slice at current energy
  REAL(8),  SAVE :: FREQ1 = 0.0D0    ! cached frequency
  INTEGER, SAVE :: N_CACHED = 0     ! cached energy index

  ! --- Local variables ---
  REAL(8)  :: EVOLT, EN, TJ, PART, TN
  INTEGER :: N, IT_PART, IT_XSEC

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING CHOP'

  ! --- Load cross-section table from file on first call ---
  IF (.NOT. CROSSCH_LOADED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'crossch.dat', STATUS='OLD', ACTION='READ', IOSTAT=IT_XSEC)
    IF (IT_XSEC .NE. 0) THEN
      WRITE(6, '(A)') ' CHOP: ERROR opening ' // trim(DATADIR) // 'crossch.dat'
      STOP 'CHOP: crossch.dat not found'
    END IF
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')  ! skip header
    READ(89, '(A)') ; READ(89, '(A)')                     ! skip header
    DO JJ = 1, 105
      READ(89, *) (CROSSCH(II, JJ), II = 1, 15)
    END DO
    CLOSE(89)
    CROSSCH_LOADED = .TRUE.
  END IF

  CHOP = 0.0D0

  !---------------------------------------------------------------------
  ! Interpolate cross-section table in energy (cached on frequency)
  !---------------------------------------------------------------------
  IF (FREQ .NE. FREQ1) THEN
    FREQ1 = FREQ
    EVOLT = WAVENO / 8065.479D0
    N = int(EVOLT * 10.0D0)
    N_CACHED = N
    EN = dble(N) * 0.1D0
    IF (N .LT. 20 .OR. N .GE. 105) RETURN
    DO IT_XSEC = 1, 15
      CROSSCHT(IT_XSEC) = CROSSCH(IT_XSEC, N) &
        + (CROSSCH(IT_XSEC, N+1) - CROSSCH(IT_XSEC, N)) * (EVOLT - EN) / 0.1D0
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Temperature-dependent interpolation at depth point J
  !---------------------------------------------------------------------
  TJ = T(J)
  N = N_CACHED
  IF (TJ .GE. 9000.0D0) RETURN
  IF (N .LT. 20 .OR. N .GE. 105) RETURN

  ! Partition function interpolation (200 K grid, T = 1000–9000 K)
  IT_PART = max(int((TJ - 1000.0D0) / 200.0D0) + 1, 1)
  IT_PART = min(IT_PART, 40)     ! guard upper bound: IT_PART+1 <= 41
  TN = dble(IT_PART) * 200.0D0 + 800.0D0
  PART = PARTCH(IT_PART) + (PARTCH(IT_PART+1) - PARTCH(IT_PART)) * (TJ - TN) / 200.0D0

  ! Cross-section interpolation (500 K grid, T = 2000–9000 K)
  IT_XSEC = max(int((TJ - 2000.0D0) / 500.0D0) + 1, 1)
  IT_XSEC = min(IT_XSEC, 14)    ! guard upper bound: IT_XSEC+1 <= 15
  TN = dble(IT_XSEC) * 500.0D0 + 1500.0D0

  CHOP = exp((CROSSCHT(IT_XSEC) &
    + (CROSSCHT(IT_XSEC+1) - CROSSCHT(IT_XSEC)) * (TJ - TN) / 500.0D0) * LN10) * PART

  RETURN

END FUNCTION CHOP


!=========================================================================
! FUNCTION OHOP(J)
!
! OH molecular bound-free + bound-bound opacity.
!
! Returns the OH contribution to the absorption coefficient at depth
! point J for the current frequency.  The cross-section is obtained by
! bilinear interpolation in a precomputed table CROSSOH(T, E):
!
!   Energy axis:  E = 2.1-15.0 eV in 0.1 eV steps  (indices 1-130)
!   Temperature axis:  T = 2000-9000 K in 500 K steps  (15 points)
!
! The table values are log10(sigma*Q) where Q is the OH partition function.
! A separate partition function table PARTOH covers T = 1000-9000 K in
! 200 K steps (41 points) and is interpolated independently.
!
! The final opacity is:  OHOP = 10^(sigma*Q interpolated) * Q(T)
!
! Frequency-dependent quantities (energy interpolation of CROSSOHT) are
! cached and only recomputed when FREQ changes between calls.
!
! Returns zero for T >= 9000 K or photon energy outside 2.1-15.0 eV.
!
! NOTE: The original code recomputed WAVENO = FREQ/c here, which is
! redundant with the module-level WAVENO already set by the caller.
! We use the module WAVENO directly instead.
!=========================================================================

FUNCTION OHOP(J)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: J
  REAL(8)  :: OHOP

  ! --- Cross-section x partition function table ---
  ! CROSSOH(IT, N): log10(sigma*Q) for temperature index IT and energy index N
  !   IT = 1-15:  T = 2000, 2500, ..., 9000 K
  !   N  = 1-130: E = 2.1, 2.2, ..., 15.0 eV
  ! Loaded from crossoh.dat on first call.
  REAL(8),    SAVE :: CROSSOH(15,130)
  LOGICAL,   SAVE :: CROSSOH_LOADED = .FALSE.
  INTEGER :: II, JJ

  ! OH partition function table: T = 1000 to 9000 K in 200 K steps (41 points)
  REAL(8), PARAMETER :: PARTOH(41) = (/ &
      145.979D0,    178.033D0,    211.618D0,    247.053D0,    284.584D0,    324.398D0,    366.639D0, &
      411.425D0,    458.854D0,    509.012D0,    561.976D0,    617.823D0,    676.626D0,    738.448D0, &
      803.363D0,    871.437D0,    942.735D0,   1017.330D0,   1095.284D0,   1176.654D0,   1261.510D0, &
     1349.898D0,   1441.875D0,   1537.483D0,   1636.753D0,   1739.733D0,   1846.434D0,   1956.883D0, &
     2071.080D0,   2189.029D0,   2310.724D0,   2436.155D0,   2565.283D0,   2698.103D0,   2834.571D0, &
     2974.627D0,   3118.242D0,   3265.366D0,   3415.912D0,   3569.837D0,   3727.077D0 /)

  ! --- Cached frequency-interpolated cross-sections ---
  REAL(8),  SAVE :: CROSSOHT(15)     ! cross-section slice at current energy
  REAL(8),  SAVE :: FREQ1 = 0.0D0    ! cached frequency
  INTEGER, SAVE :: N_CACHED = 0     ! cached energy index

  ! --- Local variables ---
  REAL(8)  :: EVOLT, EN, TJ, PART, TN
  INTEGER :: N, IT_PART, IT_XSEC

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING OHOP'

  ! --- Load cross-section table from file on first call ---
  IF (.NOT. CROSSOH_LOADED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'crossoh.dat', STATUS='OLD', ACTION='READ', IOSTAT=IT_XSEC)
    IF (IT_XSEC .NE. 0) THEN
      WRITE(6, '(A)') ' OHOP: ERROR opening ' // trim(DATADIR) // 'crossoh.dat'
      STOP 'OHOP: crossoh.dat not found'
    END IF
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')  ! skip header
    READ(89, '(A)') ; READ(89, '(A)')                     ! skip header
    DO JJ = 1, 130
      READ(89, *) (CROSSOH(II, JJ), II = 1, 15)
    END DO
    CLOSE(89)
    CROSSOH_LOADED = .TRUE.
  END IF

  OHOP = 0.0D0

  !---------------------------------------------------------------------
  ! Interpolate cross-section table in energy (cached on frequency)
  ! OHOP energy index is offset: N = EVOLT*10 - 20  (E starts at 2.1 eV)
  !---------------------------------------------------------------------
  IF (FREQ .NE. FREQ1) THEN
    FREQ1 = FREQ
    EVOLT = WAVENO / 8065.479D0
    N = int(EVOLT * 10.0D0 - 20.0D0)
    N_CACHED = N
    EN = dble(N) * 0.1D0 + 2.0D0
    IF (N .LE. 0 .OR. N .GE. 130) RETURN
    DO IT_XSEC = 1, 15
      CROSSOHT(IT_XSEC) = CROSSOH(IT_XSEC, N) &
        + (CROSSOH(IT_XSEC, N+1) - CROSSOH(IT_XSEC, N)) * (EVOLT - EN) / 0.1D0
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Temperature-dependent interpolation at depth point J
  !---------------------------------------------------------------------
  TJ = T(J)
  N = N_CACHED
  IF (TJ .GE. 9000.0D0) RETURN
  IF (N .LE. 0 .OR. N .GE. 130) RETURN

  ! Partition function interpolation (200 K grid, T = 1000-9000 K)
  IT_PART = max(int((TJ - 1000.0D0) / 200.0D0) + 1, 1)
  IT_PART = min(IT_PART, 40)     ! guard upper bound: IT_PART+1 <= 41
  TN = dble(IT_PART) * 200.0D0 + 800.0D0
  PART = PARTOH(IT_PART) + (PARTOH(IT_PART+1) - PARTOH(IT_PART)) * (TJ - TN) / 200.0D0

  ! Cross-section interpolation (500 K grid, T = 2000-9000 K)
  IT_XSEC = max(int((TJ - 2000.0D0) / 500.0D0) + 1, 1)
  IT_XSEC = min(IT_XSEC, 14)    ! guard upper bound: IT_XSEC+1 <= 15
  TN = dble(IT_XSEC) * 500.0D0 + 1500.0D0

  OHOP = exp((CROSSOHT(IT_XSEC) &
    + (CROSSOHT(IT_XSEC+1) - CROSSOHT(IT_XSEC)) * (TJ - TN) / 500.0D0) * LN10) * PART

  RETURN

END FUNCTION OHOP

!=========================================================================
! SUBROUTINE H2COLLOP
!
! H2 collision-induced dipole opacity (H2-H2 and H2-He).
!
! Computes the collision-induced absorption (CIA) from molecular hydrogen
! colliding with H2 and He, using tabulated absorption coefficients from
! Borysow, Jorgensen & Zheng (1997, A&A 324, 185).
!
! Two tables are read from h2collop.dat on first call and interpolated
! bilinearly in (wavenumber, temperature):
!   H2HE(7,81): log10 of H2-He CIA coefficient
!   H2H2(7,81): log10 of H2-H2 CIA coefficient
!   Temperature axis:  T = 1000, 2000, ..., 7000 K  (7 points)
!   Wavenumber axis:   nu = 0, 250, ..., 20000 cm-1  (81 points)
!
! The H2 equilibrium number density XNH2 is computed from EQUILH2(T)
! and cached when the temperature structure (ITEMP) changes.
!
! Final opacity:
!   AH2COLL(J) = (sigma_H2He*n_He + sigma_H2H2*n_H2) * n_H2 / rho * stim
!
! NOTE (bug fix): The original Kurucz code had swapped weights in the
! temperature interpolation: f = f_low*DELT + f_high*(1-DELT) instead
! of the standard f = f_low*(1-DELT) + f_high*DELT.  This gave more
! weight to the farther grid point.  CORRECTED in this version.
!=========================================================================

SUBROUTINE H2COLLOP(AH2COLL)

  IMPLICIT NONE

  REAL(8), INTENT(OUT) :: AH2COLL(kw)

  ! --- CIA tables: log10 of absorption coefficient ---
  ! Read from h2collop.dat on first call
  ! Borysow, Jorgensen & Zheng (1997, A&A 324, 185)
  REAL(8), SAVE :: H2HE(7,81)
  REAL(8), SAVE :: H2H2(7,81)

  ! --- Persistent state ---
  INTEGER, SAVE :: ITEMP1 = 0
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! --- Local variables ---
  REAL(8)  :: H2HENU(7), H2H2NU(7)
  REAL(8)  :: DELNU, DELT, XH2H2, XH2HE
  INTEGER :: J, IT, NU, I
  CHARACTER(256) :: LINE

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING H2COLLOP'

  !---------------------------------------------------------------------
  ! Read CIA tables from file on first call
  !---------------------------------------------------------------------
  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'h2collop.dat', STATUS='OLD', ACTION='READ')
    ! Skip comment lines starting with #
    DO
      READ(89, '(A)') LINE
      IF (LINE(1:1) .NE. '#') THEN
        BACKSPACE(89)
        EXIT
      END IF
    END DO
    ! Read H2HE table (81 lines of 7 values)
    DO I = 1, 81
      READ(89, *) H2HE(:, I)
    END DO
    ! Read H2H2 table (81 lines of 7 values)
    DO I = 1, 81
      READ(89, *) H2H2(:, I)
    END DO
    CLOSE(89)
    INITIALIZED = .TRUE.
  END IF

  !---------------------------------------------------------------------
  ! Recompute H2 equilibrium number densities when T structure changes
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO J = 1, NRHOX
      XNH2(J) = 0.0D0
      IF (T(J) .LE. 20000.0D0) THEN
        XNH2(J) = (XNFP(J,1) * 2.0D0 * BHYD(J,1))**2 * EQUILH2(T(J))
      END IF
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Return zero for wavenumber > 20000 cm-1 (above table range)
  !---------------------------------------------------------------------
  IF (WAVENO .GT. 20000.0D0) THEN
    AH2COLL = 0.0D0
    RETURN
  END IF

  !---------------------------------------------------------------------
  ! Interpolate tables in wavenumber
  !---------------------------------------------------------------------
  NU = int(WAVENO / 250.0D0) + 1
  NU = min(NU, 80)
  DELNU = (WAVENO - 250.0D0 * dble(NU - 1)) / 250.0D0

  DO IT = 1, 7
    H2H2NU(IT) = H2H2(IT, NU) * (1.0D0 - DELNU) + H2H2(IT, NU+1) * DELNU
    H2HENU(IT) = H2HE(IT, NU) * (1.0D0 - DELNU) + H2HE(IT, NU+1) * DELNU
  END DO

  !---------------------------------------------------------------------
  ! Interpolate in temperature and assemble opacity at each depth
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    IT = int(T(J) / 1000.0D0)
    IT = max(1, min(6, IT))
    DELT = (T(J) - 1000.0D0 * dble(IT)) / 1000.0D0
    DELT = max(0.0D0, min(1.0D0, DELT))

    ! Temperature interpolation (weights corrected from original)
    XH2H2 = H2H2NU(IT) * (1.0D0 - DELT) + H2H2NU(IT+1) * DELT
    XH2HE = H2HENU(IT) * (1.0D0 - DELT) + H2HENU(IT+1) * DELT

    AH2COLL(J) = (10.0D0**XH2HE * XNF(J,3) + 10.0D0**XH2H2 * XNH2(J)) &
               * XNH2(J) / RHO(J) * STIM(J)
  END DO

  RETURN

END SUBROUTINE H2COLLOP

!=========================================================================
! SUBROUTINE C2OP
!
! C II bound-free opacity.
!
! 34 energy levels ionizing to four C III limits:
!
!   Block 1 (ELIM1 = 196664.70 cm-1, C III 2s2 1S0):
!     Levels 1-5:   2s2 5g/5f/5d/5p/5s, hydrogenic n=5, l=4..0
!     Levels 6-9:   2s2 4f/4d/4p/4s,    hydrogenic n=4, l=3..0
!     Levels 10-12: 2s2 3d/3p/3s,        hydrogenic n=3, l=2..0
!     Level 13:     2s2 2p 2P (handled in HOTOP, not computed here)
!
!   Block 2 (ELIM2 = 249031.76 cm-1, C III 2s2p 3P0):
!     Levels 14-19: 2s2p 3d states, hydrogenic n=3, l=2
!     Levels 20-25: 2s2p 3p states, hydrogenic n=3, l=1
!     Levels 26-27: 2s2p 3s states, hydrogenic n=3, l=0
!     Levels 28-31: 2s2p2 states  (handled in HOTOP, not computed here)
!
!   Block 3 (ELIM3 = 299016.74 cm-1, C III 2s2p 1P1):
!     No active levels.
!
!   Block 4 (ELIM4 = 334090.40 cm-1, C III 2p2 3P0):
!     Levels 32-34: 2p3 2P/2D/4S, hydrogenic n=2 l=1, DEGEN=3
!
! Dissolved high-n contributions:
!   Block 1: n >= 6 to infinity, GFACTOR = 1
!   Block 2: n >= 4 to infinity, GFACTOR = 9
!
! Species index: XNFP(J,22), Z = 2.
!
! NOTE: GLEV(9) = 1 and GLEV(12) = 1 for 2S terms where g = 2J+1 = 2
! may be typos (cf. GLEV(5) = 2 for the analogous 5s 2S).
! Preserved from original.
!=========================================================================

SUBROUTINE C2OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLEV = 34

  ! --- Atomic data: energy levels (cm-1) and statistical weights ---
  REAL(8), PARAMETER :: ELEV(NLEV) = (/ &
    179073.05D0, 178955.94D0, 178495.47D0, 175292.30D0, 173347.84D0, &  ! 5g,5f,5d,5p,5s
    168978.34D0, 168124.17D0, 162522.34D0, 157234.07D0,              &  ! 4f,4d,4p,4s
    145550.10D0, 131731.80D0, 116537.65D0,    42.28D0,               &  ! 3d,3p,3s, 2p(inactive)
    202188.07D0, 199965.31D0, 198856.92D0, 198431.96D0, 196572.80D0, &  ! 2s2p3d states
    195786.71D0, 190000.00D0, 188601.54D0, 186452.13D0, 184690.98D0, &  ! 2s2p3d/3p states
    182036.89D0, 181741.65D0, 177787.22D0, 167009.29D0,              &  ! 2s2p3p/3s states
    110651.76D0,  96493.74D0,  74931.11D0,  43035.80D0,              &  ! 2s2p2 (inactive)
    230407.20D0, 150464.60D0, 142027.10D0 /)                            ! 2p3 states

  REAL(8), PARAMETER :: GLEV(NLEV) = (/ &
    18.D0, 14.D0, 10.D0,  6.D0,  2.D0, 14.D0, 10.D0,  6.D0,  1.D0, &  ! NOTE: GLEV(9)=1 for 2S, expected 2
    10.D0,  6.D0,  1.D0,  3.D0,                                      &  ! NOTE: GLEV(12)=1 for 2S, expected 2
     6.D0, 10.D0, 12.D0, 10.D0, 20.D0, 28.D0,  2.D0, 10.D0, 12.D0, &
     4.D0,  6.D0, 20.D0,  6.D0, 12.D0,                               &
     6.D0,  2.D0, 10.D0, 12.D0,                                      &
     6.D0, 10.D0,  4.D0 /)

  ! Quantum numbers (n,l) for XKARZAS calls
  INTEGER, PARAMETER :: NQ(NLEV) = (/ &
    5,5,5,5,5,  4,4,4,4,  3,3,3, 0,  &  ! Block 1 (level 13 inactive)
    3,3,3,3,3,3,  3,3,3,3,3,3,  3,3,  &  ! Block 2
    0,0,0,0,  2,2,2 /)                    ! levels 28-31 inactive, Block 4

  INTEGER, PARAMETER :: LQ(NLEV) = (/ &
    4,3,2,1,0,  3,2,1,0,  2,1,0, 0,  &  ! Block 1
    2,2,2,2,2,2,  1,1,1,1,1,1,  0,0,  &  ! Block 2
    0,0,0,0,  1,1,1 /)                    ! Block 4

  ! Effective charge prefactor: n² for each level (Z²eff = NPRE/RYD*(ELIM-E))
  REAL(8), PARAMETER :: NPRE(NLEV) = (/ &
    25.D0,25.D0,25.D0,25.D0,25.D0,  16.D0,16.D0,16.D0,16.D0, &  ! Block 1
     9.D0, 9.D0, 9.D0,  0.D0,                                 &
     9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0,                     &  ! Block 2
     9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0,        &
     0.D0, 0.D0, 0.D0, 0.D0,                                  &  ! inactive
     4.D0, 4.D0, 4.D0 /)                                         ! Block 4

  ! Degeneracy multiplier for inner-shell levels
  REAL(8), PARAMETER :: DEGEN(NLEV) = (/ &
    1.D0,1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0,1.D0, 1.D0, &
    1.D0,1.D0,1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0,1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0, &
    1.D0,1.D0,1.D0,1.D0, &
    3.D0,3.D0,3.D0 /)

  ! Ionization limits (cm-1)
  REAL(8), PARAMETER :: ELIM1 = 196664.70D0                 ! C III 2s2 1S0
  REAL(8), PARAMETER :: ELIM2 = 196664.70D0 + 52367.06D0   ! C III 2s2p 3P0
  REAL(8), PARAMETER :: ELIM4 = 196664.70D0 + 137425.70D0  ! C III 2p2 3P0

  REAL(8), PARAMETER :: RYD = 109732.298D0
  REAL(8), PARAMETER :: Z4 = 16.0D0    ! Z**4 (Z=2 for C II)
  REAL(8), PARAMETER :: Z2 = 4.0D0     ! Z**2

  ! --- Persistent state ---
  REAL(8),  SAVE :: BOLT(NLEV, kw)
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Local variables ---
  REAL(8)  :: X(NLEV)
  REAL(8)  :: ELIM, ZEFF2, FREQ3, H
  INTEGER :: I, J, K

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING C2OP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors when temperature structure changes
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      DO I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      END DO
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Compute cross-sections at the current frequency
  !---------------------------------------------------------------------
  X(:) = 0.0D0

  ! --- Block 1: excited states -> C III 2s2 1S0 ---
  ELIM = ELIM1
  DO I = 1, 12
    IF (WAVENO .LT. ELIM - ELEV(I)) EXIT
    ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
  END DO
  ! Level 13 (2s2 2p 2P) handled in HOTOP

  ! --- Block 2: inner-shell states -> C III 2s2p 3P0 ---
  ELIM = ELIM2
  DO I = 14, 27
    IF (WAVENO .LT. ELIM - ELEV(I)) EXIT
    ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
  END DO
  ! Levels 28-31 (2s2p2 states) handled in HOTOP

  ! --- Block 4: 2p3 states -> C III 2p2 3P0 ---
  ELIM = ELIM4
  DO I = 32, 34
    IF (WAVENO .LT. ELIM - ELEV(I)) EXIT
    ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I)) * DEGEN(I)
  END DO

  !---------------------------------------------------------------------
  ! Assemble opacity over depth: levels + dissolved high-n contributions
  !---------------------------------------------------------------------
  FREQ3 = 2.815D29 / (FREQ * FREQ * FREQ) * Z4

  DO J = 1, NRHOX
    ! Dissolved levels, Block 1: n >= 6 to infinity, GFACTOR = 1
    H = FREQ3 * 1.0D0 * 2.0D0 / 2.0D0 / (RYD * Z2 * HCKT(J)) &
      * (exp(-max(ELIM1 - RYD * Z2 / 36.0D0, ELIM1 - WAVENO) * HCKT(J)) &
       - exp(-ELIM1 * HCKT(J)))

    ! Dissolved levels, Block 2: n >= 4 to infinity, GFACTOR = 9
    H = H + FREQ3 * 9.0D0 * 2.0D0 / 2.0D0 / (RYD * Z2 * HCKT(J)) &
      * (exp(-max(ELIM2 - RYD * Z2 / 16.0D0, ELIM2 - WAVENO) * HCKT(J)) &
       - exp(-ELIM2 * HCKT(J)))

    ! Sum bound-free over all active levels
    DO I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    END DO

    AC2(J) = H * XNFP(J, 22) * STIM(J) / RHO(J)
  END DO

  RETURN

END SUBROUTINE C2OP

!=======================================================================
! N1OP: N I bound-free opacity cross section
!=======================================================================

REAL(8) FUNCTION N1OP(J)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: J

  ! N I cross-section times partition function
  ! Three Seaton-formula edges: 853A (ground), 1020A, 1130A (excited)
  REAL(8), SAVE :: X853 = 0.0D0, X1020 = 0.0D0, X1130 = 0.0D0
  REAL(8), SAVE :: C1130(kw), C1020(kw)
  REAL(8), SAVE :: FREQ1 = 0.0D0
  INTEGER, SAVE :: ITEMP1 = 0
  INTEGER :: K

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING N1OP'

  ! Recompute Boltzmann factors when temperature changes
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      C1130(K) = 6.0D0 * exp(-3.575D0 / TKEV(K))
      C1020(K) = 10.0D0 * exp(-2.384D0 / TKEV(K))
    END DO
  END IF

  ! Recompute cross-sections when frequency changes
  IF (FREQ .NE. FREQ1) THEN
    X1130 = 0.0D0
    X1020 = 0.0D0
    X853 = 0.0D0
    IF (FREQ .GE. 3.517915D15) X853 = SEATON(3.517915D15, 1.142D-17, 2.0D0, 4.29D0)
    IF (FREQ .GE. 2.941534D15) X1020 = SEATON(2.941534D15, 4.41D-18, 1.5D0, 3.85D0)
    IF (FREQ .GE. 2.653317D15) X1130 = SEATON(2.653317D15, 4.2D-18, 1.5D0, 4.34D0)
    FREQ1 = FREQ
  END IF

  N1OP = X853 * 4.0D0 + X1020 * C1020(J) + X1130 * C1130(J)
  RETURN

END FUNCTION N1OP

!=======================================================================
! O1OP
!=======================================================================

FUNCTION O1OP(J)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: J
  REAL(8) :: O1OP

  ! O I bound-free cross-section times partition function (after Peach)
  REAL(8), SAVE :: X911 = 0.0D0
  REAL(8), SAVE :: FREQ1 = 0.0D0

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING O1OP'

  ! Recompute cross-section only when frequency changes
  IF (FREQ .NE. FREQ1) THEN
    X911 = 0.0D0
    IF (FREQ .GE. FREQ_RYDH) X911 = SEATON(FREQ_RYDH, 2.94D-18, 1.0D0, 2.66D0)
    FREQ1 = FREQ
  END IF

  O1OP = X911 * 9.0D0
  RETURN

END FUNCTION O1OP

!=========================================================================
! SUBROUTINE MG2OP
!
! Mg II bound-free opacity.
!
! 14 energy levels ionizing to Mg III 2p6 1S0 (ELIM = 121267.61 cm-1).
! Z = 2 (Mg III core charge), species index XNFP(J,79).
!
! Levels 1-2:  n=7,6 super-levels (all l combined), hydrogenic XKARZAS
!   Level 1: n=7 composite, GLEV=98=2*49, XKARZAS(7,7)
!   Level 2: n=6 composite, GLEV=72=2*36, XKARZAS(6,6)
!   NOTE: l=n in XKARZAS calls (super-level convention)
! Levels 3-7:  n=5 states (5g,5f,5d,5p,5s), Z2eff=25/RYD*(ELIM-E)
! Levels 8-9:  n=4 states (4f,4d), Z2eff=16/RYD*(ELIM-E)
! Level 10:    4p 2P, Z2eff=16, (4,1) [Stift 2009 correction from (4,2)]
! Level 11:    4s 2S, Z2eff=16, (4,0) [GLEV corrected from 22 to 2]
! Level 12:    3d 2D, Z2eff=9/RYD*(ELIM-E)
! Level 13:    3p 2P, Z2eff=9/RYD*(ELIM-E)
! Level 14:    3s 2S (ground), non-hydrogenic power-law fit
!
! Dissolved high-n: n >= 8 to infinity, GFACTOR = 2 (implicit in code).
!
! NOTE: Thresholds are NOT monotonically increasing (the 5s level has
! a higher threshold than 4f, and 4s higher than 3d).  The original code
! tests levels sequentially with GO TO 30 on failure, which means levels
! 8-14 are skipped if level 7 fails, even if level 8 is above threshold.
! This is a pre-existing behavior preserved here for fidelity.
!=========================================================================

SUBROUTINE MG2OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLEV = 14

  ! --- Atomic data: energy levels (cm-1) and statistical weights ---
  REAL(8), PARAMETER :: ELEV(NLEV) = (/ &
    112197.00D0, 108900.00D0, 103705.66D0, 103689.89D0, 103419.82D0, &  ! n=7,6, 5g,5f,5d
     97464.32D0,  92790.51D0,  93799.70D0,  93310.80D0,  80639.85D0, &  ! 5p,5s, 4f,4d, 4p
     69804.95D0,  71490.54D0,  35730.36D0,      0.00D0 /)                ! 4s, 3d, 3p, 3s

  ! GLEV(11) corrected from 22 to 2 (4s 2S: g=2J+1=2)
  REAL(8), PARAMETER :: GLEV(NLEV) = (/ &
    98.D0, 72.D0, 18.D0, 14.D0, 10.D0, 6.D0, 2.D0, &
    14.D0, 10.D0,  6.D0,  2.D0, 10.D0, 6.D0, 2.D0 /)

  ! Quantum numbers (n,l) for XKARZAS calls
  ! Levels 1-2: l=n convention for super-levels
  ! Level 10: l=1 (Stift 2009 correction from l=2)
  INTEGER, PARAMETER :: NQ(NLEV) = (/ 7,6, 5,5,5,5,5, 4,4,4,4, 3,3, 0 /)
  INTEGER, PARAMETER :: LQ(NLEV) = (/ 7,6, 4,3,2,1,0, 3,2,1,0, 2,1, 0 /)

  ! Effective charge prefactor: n2 for each level
  REAL(8), PARAMETER :: NPRE(NLEV) = (/ &
    49.D0, 36.D0, 25.D0, 25.D0, 25.D0, 25.D0, 25.D0, &
    16.D0, 16.D0, 16.D0, 16.D0,  9.D0,  9.D0,  0.D0 /)

  ! Ionization limit (cm-1): Mg III 2p6 1S0
  REAL(8), PARAMETER :: ELIM = 121267.61D0
  REAL(8), PARAMETER :: RYD  = 109732.298D0
  REAL(8), PARAMETER :: Z2   = 4.0D0     ! Z**2 (Z=2 for Mg II)
  REAL(8), PARAMETER :: Z4   = 16.0D0    ! Z**4

  ! --- Persistent state ---
  REAL(8),  SAVE :: BOLT(NLEV, kw)
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Local variables ---
  REAL(8)  :: X(NLEV)
  REAL(8)  :: ZEFF2, FREQ3, H, RATIO
  INTEGER :: I, J, K

  ! --- External functions ---

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING MG2OP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors when temperature structure changes
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      DO I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      END DO
    END DO
  END IF

  !---------------------------------------------------------------------
  ! Compute cross-sections at the current frequency
  !---------------------------------------------------------------------
  X(:) = 0.0D0

  levels: DO

    ! Levels 1-13: hydrogenic via XKARZAS
    DO I = 1, 13
      IF (WAVENO .LT. ELIM - ELEV(I)) EXIT levels
      ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
      X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
    END DO

    ! Level 14: 3s 2S ground state — non-hydrogenic power-law fit
    IF (WAVENO .LT. ELIM - ELEV(14)) EXIT levels
    RATIO = (ELIM - ELEV(14)) / WAVENO
    X(14) = 0.14D-18 * (6.700D0 * RATIO**4 - 5.700D0 * RATIO**5)

    EXIT levels
  END DO levels

  !---------------------------------------------------------------------
  ! Assemble opacity over depth: levels + dissolved high-n contribution
  !---------------------------------------------------------------------
  FREQ3 = 2.815D29 / (FREQ * FREQ * FREQ) * Z4

  DO J = 1, NRHOX
    ! Dissolved levels: n >= 8 to infinity
    ! GFACTOR = 2 (implicit: 2./2. = 1 in prefactor)
    H = FREQ3 * 2.0D0 / 2.0D0 / (RYD * Z2 * HCKT(J)) &
      * (exp(-max(ELIM - RYD * Z2 / 64.0D0, ELIM - WAVENO) * HCKT(J)) &
       - exp(-ELIM * HCKT(J)))

    DO I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    END DO

    AMG2(J) = H * XNFP(J, 79) * STIM(J) / RHO(J)
  END DO

  RETURN

END SUBROUTINE MG2OP

!==========================================================================
! SUBROUTINE SI2OP
!
! Si II bound-free photoionization opacity.
!
! Uses a pretabulated opacity strategy for efficiency: on first call,
! builds a 2D lookup table XTAB(200,51) of log(opacity) indexed by
! wavenumber (1000–200000 cm⁻¹, step 1000) and temperature (51 points,
! T = 10^(3.48 + k*0.02) for k=1..51, spanning ~3020–60250 K).
!
! 46 energy levels organized by Si III parent ion limit:
!
!   Levels 1–11:  3s² nl excited states → Si III 3s² ¹S
!                 (ELIM = 131838.4 cm⁻¹). XKARZAS cross-sections.
!
!   Levels 12–37: 3s3p nl states → Si III 3s3p ³P₀
!                 (ELIM = 184563.09 cm⁻¹). XKARZAS cross-sections.
!
!   Levels 38–41: 3s3p² states → Si III 3s3p ³P₀
!                 (ELIM = 184563.09 cm⁻¹). XKARZAS × 2 (degeneracy).
!
!   Levels 42–44: 3p³ states → Si III 3p² ³P₀
!                 (ELIM = 254052.92 cm⁻¹). XKARZAS × 3 (degeneracy).
!
!   Levels 45–46: High-n series limits (n≥6 for 3s², n≥5 for 3s3p).
!
! At runtime: bilinear interpolation in the pretabulated grid for
! wavenumber ≥ 12192.48 cm⁻¹; pure hydrogenic for lower frequencies.
!==========================================================================

SUBROUTINE SI2OP

  IMPLICIT NONE

  INTEGER, PARAMETER :: NWGRID = 200, NTGRID = 51
  REAL(8), PARAMETER :: RYD = 109732.298D0, Z = 2.0D0

  ! Ionization limits (cm⁻¹) for each level
  REAL(8), PARAMETER :: ELIMLEV(46) = (/ &
    131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, &  !  1-5
    131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, &  !  6-10
    131838.4D0, &                                                    !  11
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  12-15
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  16-19
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  20-23
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  24-27
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  28-31
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  32-35
    184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &           !  36-39
    184563.09D0, 184563.09D0, &                                      !  40-41
    254052.92D0, 254052.92D0, 254052.92D0, &                         !  42-44
    131838.4D0, 184563.09D0 /)                                       !  45-46

  ! Energy levels (cm⁻¹)
  REAL(8), PARAMETER :: ELEV(46) = (/ &
    114177.4D0, 113760.48D0, 112394.92D0, 103877.34D0, 97972.35D0, &  !  1-5:  3s²5g,5f,5d,5p,5s
    103556.36D0, 101024.09D0, 81231.57D0, 79348.67D0, 65500.73D0, &   !  6-10: 3s²4f,4d,4p,3d,4s
    191.55D0, &                                                        !  11:   3s²3p
    157396.6D0, 157188.8D0, 156838.9D0, 156836.9D0, 155663.4D0, &    !  12-16: 3s3p4f
    155593.7D0, 155555.0D0, 153523.1D0, 153147.2D0, 152977.0D0, &    !  17-21: 3s3p4f,4d,4p
    152480.7D0, 151245.1D0, 149905.6D0, 140696.0D0, 134905.34D0, &   !  22-26: 3s3p4d,4p
    134136.03D0, 132648.5D0, 132012.27D0, 131815.5D0, 126250.9D0, &  !  27-31: 3s3p4p,3d
    124595.5D0, 124373.8D0, 121541.76D0, 117058.95D0, 114415.54D0, & !  32-36: 3s3p3d,4s,3d
    108804.1D0, 83937.09D0, 76665.61D0, 55319.11D0, 43002.27D0, &    !  37-41: 3s3p3d,3p²
    143990.0D0, 135300.5D0, 123033.6D0, &                             !  42-44: 3p³
    119645.92D0, 167005.92D0 /)                                        !  45-46: high-n limits

  ! Threshold energies TLEV = ELIM - ELEV (cm⁻¹)
  REAL(8), PARAMETER :: TLEV(46) = (/ &
    17661.0D0, 18077.92D0, 19443.48D0, 27961.06D0, 33866.05D0, &
    28282.04D0, 30814.31D0, 50606.83D0, 52489.73D0, 66337.67D0, &
    131646.85D0, &
    27166.49D0, 27374.29D0, 27724.19D0, 27726.19D0, 28899.69D0, &
    28969.39D0, 29008.09D0, 31039.99D0, 31415.89D0, 31586.09D0, &
    32082.39D0, 33317.99D0, 34657.49D0, 43867.09D0, 49657.75D0, &
    50427.06D0, 51914.59D0, 52550.82D0, 52747.59D0, 58312.19D0, &
    59967.59D0, 60189.29D0, 63021.33D0, 67504.14D0, 70147.55D0, &
    75758.99D0, 100526.00D0, 107897.48D0, 129243.98D0, 141560.82D0, &
    110052.92D0, 118752.42D0, 131019.32D0, &
    12192.48D0, 17557.17D0 /)

  ! Statistical weights
  REAL(8), PARAMETER :: GLEV(46) = (/ &
    18.D0, 14.D0, 10.D0, 6.D0, 2.D0, 14.D0, 10.D0, 6.D0, 10.D0, &
    1.D0, 6.D0, 20.D0, 10.D0, 18.D0, 36.D0, 28.D0, 10.D0, 10.D0, &
    6.D0, 12.D0, 2.D0, 20.D0, 28.D0, 10.D0, 10.D0, 4.D0, 12.D0, &
    6.D0, 20.D0, 10.D0, 6.D0, 12.D0, 20.D0, 6.D0, 12.D0, 28.D0, &
    10.D0, 6.D0, 2.D0, 10.D0, 12.D0, 6.D0, 10.D0, 4.D0, &
    1.D0, 9.D0 /)

  ! Angular momentum quantum numbers
  INTEGER, PARAMETER :: LLEV(46) = (/ &
    4, 3, 2, 1, 0, 3, 2, 1, 2, 0, 1, &
    3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, &
    2, 2, 2, 0, 0, 2, 2, 1, 1, 1, 1, &
    1, 1, 1, 0, 0 /)

  ! Principal quantum numbers
  INTEGER, PARAMETER :: NQLEV(46) = (/ &
    5, 5, 5, 5, 5, 4, 4, 4, 3, 4, 3, &
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, &
    3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 6, 5 /)

  ! Pretabulated opacity grid and supporting arrays
  REAL(8),  SAVE :: XTAB(NWGRID, NTGRID) = 0.0D0
  REAL(8),  SAVE :: HCKTTAB(NTGRID), BOLT3s2(NTGRID), BOLT3s3p(NTGRID)
  REAL(8),  SAVE :: BOLT(44, NTGRID)
  REAL(8),  SAVE :: BOLTN(kw)
  INTEGER, SAVE :: INDEXT(kw)
  REAL(8),  SAVE :: TFRAC(kw)
  INTEGER, SAVE :: ITEMP1 = 0
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  REAL(8)  :: ZEFF2LEV(44)
  REAL(8)  :: X(46), FREQ3, WNOTAB, FREQTAB, TTAB, H, WFRAC, TLOG10
  INTEGER :: I, J, K, NU, IT, IW

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING SI2OP'

  ! ==================================================================
  ! One-time pretabulation of opacity grid
  ! ==================================================================
  IF (.NOT. INITIALIZED) THEN
    ! Compute effective charges for levels 1-44
    DO I = 1, 44
      ZEFF2LEV(I) = dble(NQLEV(I))**2 / RYD * TLEV(I)
    END DO

    ! Build temperature grid
    DO K = 1, NTGRID
      TTAB = 10.0D0**(3.48D0 + K * 0.02D0)
      HCKTTAB(K) = HCK / TTAB
      BOLT3s2(K) = exp(-ELIMLEV(1) * HCKTTAB(K))
      BOLT3s3p(K) = exp(-ELIMLEV(12) * HCKTTAB(K))
      DO I = 1, 44
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKTTAB(K))
      END DO
    END DO

    ! Build opacity table: loop over wavenumber grid
    DO NU = 1, NWGRID
      WNOTAB = dble(NU) * 1000.0D0
      FREQTAB = WNOTAB * CLIGHT
      FREQ3 = 2.815D29 / FREQTAB / FREQTAB / FREQTAB * Z**4

      DO I = 1, 46
        X(I) = 0.0D0
      END DO

      ! Levels 1-11: → Si III 3s²
      DO I = 1, 11
        IF (WNOTAB .LT. TLEV(I)) EXIT
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I))
      END DO

      ! Levels 12-37: → Si III 3s3p ³P₀
      DO I = 12, 37
        IF (WNOTAB .LT. TLEV(I)) EXIT
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I))
      END DO

      ! Levels 38-41: → Si III 3s3p ³P₀ (×2 degeneracy)
      DO I = 38, 41
        IF (WNOTAB .LT. TLEV(I)) EXIT
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I)) * 2.0D0
      END DO

      ! Levels 42-44: → Si III 3p² ³P₀ (×3 degeneracy)
      DO I = 42, 44
        IF (WNOTAB .LT. TLEV(I)) EXIT
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I)) * 3.0D0
      END DO

      ! Assemble opacity for each temperature point
      DO K = 1, NTGRID
        ! High-n series: 3s² n≥6
        H = FREQ3 * GLEV(45) * 2.0D0 / 2.0D0 / (RYD * Z**2 * HCKTTAB(K)) &
          * (exp(-max(ELEV(45), ELIMLEV(45) - WNOTAB) * HCKTTAB(K)) &
          - BOLT3s2(K))
        ! High-n series: 3s3p n≥5
        H = H + FREQ3 * GLEV(46) * 2.0D0 / 2.0D0 / (RYD * Z**2 * HCKTTAB(K)) &
          * (exp(-max(ELEV(46), ELIMLEV(46) - WNOTAB) * HCKTTAB(K)) &
          - BOLT3s3p(K))
        ! Sum resolved levels
        DO I = 1, 44
          H = H + X(I) * BOLT(I, K)
        END DO
        XTAB(NU, K) = log(H)
      END DO
    END DO

    INITIALIZED = .TRUE.
  END IF

  ! ==================================================================
  ! Temperature-dependent precomputation (when ITEMP changes)
  ! ==================================================================
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO J = 1, NRHOX
      BOLTN(J) = (GLEV(45) * exp(-ELIMLEV(45) * HCKT(J)) &
               + GLEV(46) * exp(-ELIMLEV(46) * HCKT(J))) &
               / (RYD * Z**2 * HCKT(J))
      TLOG10 = TLOG(J) / LN10
      IT = int((TLOG10 - 3.48D0) / 0.02D0)
      IT = max(min(IT, 50), 1)
      INDEXT(J) = IT
      TFRAC(J) = (TLOG10 - 3.48D0 - IT * 0.02D0) / 0.02D0
    END DO
  END IF

  ! ==================================================================
  ! Frequency evaluation: bilinear interpolation or low-freq fallback
  ! ==================================================================
  IF (WAVENO .GE. 12192.48D0) THEN
    ! Bilinear interpolation in pretabulated grid
    IW = int(WAVENO * 0.001D0)
    IW = max(min(IW, 199), 1)
    WFRAC = (WAVENO - IW * 1000.0D0) / 1000.0D0
    DO J = 1, NRHOX
      IT = INDEXT(J)
      H = (XTAB(IW, IT) * (1.0D0 - TFRAC(J)) &
        + XTAB(IW, IT+1) * TFRAC(J)) * (1.0D0 - WFRAC) &
        + (XTAB(IW+1, IT) * (1.0D0 - TFRAC(J)) &
        + XTAB(IW+1, IT+1) * TFRAC(J)) * WFRAC
      ASI2(J) = exp(H) * XNFP(J, 106) * STIM(J) / RHO(J)
    END DO
  ELSE
    ! Low frequency: only high-n hydrogenic contribution
    ! Si III 3s² (ELIM=131838.4, n≥6) + Si III 3s3p ³P₀ (ELIM=184563.09, n≥5)
    FREQ3 = 2.815D29 / FREQ / FREQ / FREQ * Z**4
    DO J = 1, NRHOX
      H = FREQ3 * (1.0D0 / EHVKT(J) - 1.0D0) * BOLTN(J)
      ASI2(J) = H * XNFP(J, 106) * STIM(J) / RHO(J)
    END DO
  END IF
  RETURN

END SUBROUTINE SI2OP

!=======================================================================
! CA2OP: Ca II bound-free photoionization cross-section
!
! Cross-section × partition function for Ca II (singly-ionized calcium).
! Three edges:
!   1044 Å (ground 4s ²S, g=2): ν⁻³ hydrogenic
!   1218 Å (3d ²D, g=10): ν⁻½ fit, Boltzmann factor exp(-1.697/kT)
!   1420 Å (4p ²P, g=6): Seaton formula, Boltzmann factor exp(-3.142/kT)
! Temperature- and frequency-cached for efficiency.
!=======================================================================

FUNCTION CA2OP(J)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: J
  REAL(8) :: CA2OP

  ! Ca II bound-free cross-section times partition function
  ! Three edges: 1044A (ground), 1218A, 1420A (excited)
  REAL(8), SAVE :: C1218(kw), C1420(kw)
  REAL(8), SAVE :: X1044 = 0.0D0, X1218 = 0.0D0, X1420 = 0.0D0
  REAL(8), SAVE :: FREQ1 = 0.0D0
  INTEGER, SAVE :: ITEMP1 = 0
  INTEGER :: K

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING CA2OP'

  ! Recompute Boltzmann factors when temperature changes
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      C1218(K) = 10.0D0 * exp(-1.697D0 / TKEV(K))
      C1420(K) = 6.0D0 * exp(-3.142D0 / TKEV(K))
    END DO
  END IF

  ! Recompute cross-sections when frequency changes
  IF (FREQ .NE. FREQ1) THEN
    X1420 = 0.0D0
    X1218 = 0.0D0
    X1044 = 0.0D0
    IF (FREQ .GE. 2.870454D15) X1044 = 5.4D-20 * (2.870454D15 / FREQ)**3
    IF (FREQ .GE. 2.460127D15) X1218 = 1.64D-17 * sqrt(2.460127D15 / FREQ)
    IF (FREQ .GE. 2.110779D15) X1420 = SEATON(2.110779D15, 4.13D-18, 3.0D0, 0.69D0)
    FREQ1 = FREQ
  END IF

  CA2OP = X1044 * 2.0D0 + X1218 * C1218(J) + X1420 * C1420(J)
  RETURN

END FUNCTION CA2OP

!=======================================================================
! ELECOP
!=======================================================================

SUBROUTINE ELECOP

  IMPLICIT NONE
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING ELECOP'
  ! Thomson electron scattering: sigma_T = 0.6653e-24 cm^2
  SIGEL = SIGMA_THOMSON * XNE / RHO
  RETURN

END SUBROUTINE ELECOP

!=======================================================================
! H2RAOP: H₂ Rayleigh scattering opacity
!
! Reference:
!   Dalgarno, A. & Williams, D. W. 1962, ApJ 136, 690,
!   "Rayleigh Scattering by Molecular Hydrogen."
!
! The Dalgarno & Williams (1962) cross-section for Rayleigh scattering
! by ground-state molecular hydrogen is a three-term expansion in
! inverse powers of wavelength:
!
!   σ(H₂) = (8.14e-13 / λ⁴) + (1.28e-6 / λ⁶) + (1.61 / λ⁸)     [cm²]
!
! with λ in Ångströms.  The leading 1/λ⁴ term is the standard Rayleigh
! low-frequency limit from the static polarizability; the 1/λ⁶ and
! 1/λ⁸ terms are dispersion corrections that become important in the
! blue/UV as the photon energy approaches the lowest H₂ electronic
! transitions.  The formula is a closed-form fit to a Kramers-Heisenberg
! sum over H₂ excited electronic states, evaluated at frequencies well
! below the Lyman and Werner band systems (~1110 Å and below).
!
! Validity regime:
!   The Dalgarno & Williams formula is accurate for λ >> 1110 Å, i.e.
!   well into the red of the H₂ Lyman/Werner electronic bands, where
!   the Rayleigh approximation (elastic scattering with the photon
!   energy far from any resonance) is clean.  Approaching the Lyman
!   bands from the red, real bound-bound and bound-free absorption
!   takes over and this pure-Rayleigh expression is no longer valid.
!   To prevent extrapolation into that regime, we cap the evaluated
!   frequency at 2.922e15 Hz (λ ≈ 1026 Å), which places the cutoff
!   just blueward of the visible/near-UV region where the formula is
!   reliable and just redward of the onset of H₂ electronic bands.
!   For ν > 2.922e15 Hz the cross section is frozen at its λ = 1026 Å
!   value.  In stellar atmosphere calculations this regime is of no
!   practical concern: the H₂ number density is negligible in layers
!   where FUV photons dominate.
!
! H₂ number density:
!   XNH2(J) is recomputed from chemical equilibrium via EQUILH2(T(J))
!   whenever ITEMP changes, and cached otherwise.  The formation
!   expression is
!
!      XNH2 = [2 * N(H I) * b(H I,n=1)]² * K_eq(T)
!
!   where N(H I) comes from XNFP(J,1), b(H I,n=1) = BHYD(J,1) is the
!   NLTE departure coefficient of the ground state (unity in LTE), the
!   factor of 2 converts proton number density to atomic H nucleus
!   density, and K_eq(T) = EQUILH2(T) is the H₂/H equilibrium constant.
!   Above T = 20,000 K the molecule is fully dissociated and XNH2 is
!   set to zero to avoid spurious contributions.
!
!   (The commented-out blocks below preserve two earlier inline
!   polynomial fits for the equilibrium constant that were superseded
!   by the external EQUILH2 function during modernization; they are
!   kept as a historical record of the F77 expressions.)
!
! Scaling:
!   SIGH2(J) = σ(λ) * XNH2(J) / RHO(J)
!
!   where the /RHO(J) factor converts number density to a per-unit-mass
!   opacity (cm²/g), matching the units convention of the other SIG*
!   arrays.  The /RHO(J) correction was added by K. Bischof on
!   14 Sep 2004; the original F77 omitted it (bug).
!=======================================================================

SUBROUTINE H2RAOP

  IMPLICIT NONE

  INTEGER, SAVE :: ITEMP1 = 0
  REAL(8)  :: W, WW, SIG
  INTEGER :: J

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING H2RAOP'

  ! Recompute H₂ number density when the temperature structure changes.
  ! XNH2 is a function of T(J) only (through EQUILH2) and the H I ground
  ! state population, so it is cached between calls at the same ITEMP.
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO J = 1, NRHOX
      XNH2(J) = 0.0D0
      IF (T(J) .LE. 20000.0D0) THEN
        ! Historical F77 inline polynomial fits for ln(K_eq(T)),
        ! superseded by EQUILH2(T) during modernization.  Preserved
        ! here only as a reference to the original expressions.
!C  11 XNH2(J)=(XNFP(J,1)*2.*BHYD(J,1))**2*EXP(4.477/TKEV(J)-4.6628E1+
!C    1(1.8031D-3+(-5.0239D-7+(8.1424D-11-5.0501D-15*T(J))*T(J))*T(J))*
!C    2T(J)-1.5*TLOG(J))/RHO(J)
!      XNH2(J)=(XNFP(J,1)*2.*BHYD(J,1))**2*EXP(4.478/TKEV(J)-4.64584E1+
!     1(1.63660D-3+(-4.93992D-7+(1.11822D-10+(-1.49567D-14+
!     2(1.06206D-18-3.08720D-23*T(J))*T(J))*T(J))*T(J))*T(J))*T(J)-
!     3 1.5*TLOG(J))
        XNH2(J) = (XNFP(J,1) * 2 * BHYD(J,1))**2 * EQUILH2(T(J))
      END IF
    END DO
  END IF

  ! Evaluate the Dalgarno & Williams (1962) cross-section.  Cap the
  ! frequency at 2.922e15 Hz (λ ≈ 1026 Å) to prevent extrapolation
  ! into the H₂ Lyman/Werner band region where the Rayleigh
  ! approximation breaks down.
  W  = CLIGHT_ANGS / min(FREQ, 2.922D15)   ! wavelength [Å], capped
  WW = W**2                                ! λ² [Å²]

  ! σ(λ) = (8.14e-13 + 1.28e-6/λ² + 1.61/λ⁴) / λ⁴   [cm²]
  ! Rewritten with WW to share the λ² computation between terms.
  SIG = (8.14D-13 + 1.28D-6 / WW + 1.61D0 / (WW * WW)) / (WW * WW)

  ! Apply to each depth point.  The /RHO(J) factor converts from
  ! cross-section-times-number-density [cm⁻¹] to mass opacity [cm²/g].
  SIGH2 = SIG * XNH2 / RHO
  RETURN

END SUBROUTINE H2RAOP

!=========================================================================
! SUBROUTINE HLINOP
!
! ***************** DEAD CODE *****************
!
! Hydrogen line opacity using Stark-broadened profiles.
!
! Computes the absorption coefficient AHLINE(J) and source function
! SHLINE(J) from hydrogen bound-bound transitions for lower levels
! N = 1 (Lyman), 2 (Balmer), 3 (Paschen), 4 (Brackett).
!
! For each depth point J, the line opacity is summed over a window of
! upper levels M around the level nearest to the current frequency.
! Each transition uses Stark-broadened profiles from STARK(N,M,J).
!
! Each upper level M is weighted by its occupation probability w_M
! from the Hummer & Mihalas (1988) formalism, which smoothly
! transitions levels from fully bound (w=1) to fully dissolved (w=0)
! based on the Holtsmark microfield distribution.  This replaces the
! former sharp cutoff at MLAST = 1100 / N_e^(2/15) from the
! Inglis-Teller formula.
!
! The dissolved fraction (1 - w_M) of each level's opacity is handled
! by HOP's pseudo-continuum loop (Kramers bound-free extrapolated to
! all frequencies), not by HLINOP.  This avoids mixing continuum
! opacity into the line opacity array (AHLINE), which would distort
! line cores.
!
! Note: this is distinct from the merged-continuum tapering used in
! SYNTHE's compute_line_opacity, which retains the empirical
! Inglis-Teller formula with a larger coefficient (1600 vs 1100)
! because the tapering serves a different purpose (opacity accounting
! near series limits, not level dissolution).
!
! Requires external functions: STARK(N,M,J), COULX(N,freq,Z)
!
! References:
!   Hummer, D.G. & Mihalas, D. 1988, ApJ 331, 794
!   Inglis, D.R. & Teller, E. 1939, ApJ 90, 439 (predecessor)
!=========================================================================

SUBROUTINE HLINOP

  IMPLICIT NONE

  REAL(8), PARAMETER :: NU_LYMAN = FREQ_RYDH   ! Lyman limit frequency (Hz)

  ! Maximum level for occupation probability computation.
  ! Beyond this, w_n is effectively 0 for any stellar density.
  INTEGER, PARAMETER :: NMAX_OCC = 80

  ! --- Persistent state ---
  REAL(8),  SAVE :: BOLT(kw, 4)            ! Boltzmann × population for n=1..4
  REAL(8),  SAVE :: W_OCC(kw, NMAX_OCC)    ! occupation probabilities per depth
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Local variables ---
  REAL(8)  :: H, S, A, BHYDJM, w_m
  INTEGER :: J, N, M, M1, M2, MFREQ

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HLINOP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors and occupation probabilities when T changes
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    DO J = 1, NRHOX
      ! Occupation probabilities for levels 2..NMAX_OCC
      W_OCC(J, 1) = 1.0D0
      DO N = 2, NMAX_OCC
        W_OCC(J, N) = occupation_prob(N, XNE(J))
      END DO
      ! Boltzmann factors for lower levels 1..4
      DO N = 1, 4
        BOLT(J, N) = exp(-(13.595D0 - 13.595D0 / dble(N)**2) / TKEV(J)) &
                   * 2.0D0 * dble(N)**2 * BHYD(J, N) * XNFP(J, 1) / RHO(J)
      END DO
    END DO
    ITEMP1 = ITEMP
  END IF

  !---------------------------------------------------------------------
  ! Determine lower level N from frequency
  !---------------------------------------------------------------------
  N = int(sqrt(NU_LYMAN / FREQ))
  IF (N .EQ. 0 .OR. N .GT. 4) RETURN

  ! Low-frequency cutoffs by series
  SELECT CASE (N)
  CASE (1)
    IF (FREQ .LT. 2.0D15) RETURN       ! Lyman: lambda > 1500 A
  CASE (2)
    IF (FREQ .LT. 4.44D14) RETURN       ! Balmer: lambda > 6756 A
  END SELECT

  !---------------------------------------------------------------------
  ! Upper level nearest to current frequency
  !---------------------------------------------------------------------
  MFREQ = nint(sqrt(NU_LYMAN / (NU_LYMAN / dble(N)**2 - FREQ)))

  !---------------------------------------------------------------------
  ! Sum Stark-broadened line opacity over depth
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    M1 = MFREQ
    M2 = M1 + 1
    M1 = max(M1, N + 1)
    H = 0.0D0
    S = 0.0D0

    IF (M1 .GT. 6) THEN
      ! High upper levels: check if effectively dissolved
      IF (M1 .LE. NMAX_OCC) THEN
        IF (W_OCC(J, M1) .LT. 1.0D-3) THEN
          ! Fully dissolved: no line opacity.
          ! The continuum opacity is provided by HOP's pseudo-continuum.
          AHLINE(J) = 0.0D0
          SHLINE(J) = BNU(J)
          CYCLE
        END IF
      ELSE
        ! Beyond tabulated range: fully dissolved
        AHLINE(J) = 0.0D0
        SHLINE(J) = BNU(J)
        CYCLE
      END IF
      M1 = M1 - 1
      M2 = M2 + 3
      ! Special case: add Paschen-alpha (3->4) if computing Brackett (N=4)
      IF (N .GE. 4 .AND. M1 .LE. 8) THEN
        H = STARK_MMM(3, 4, J) * (1.0D0 - EHVKT(J) * BHYD(J, 4) / BHYD(J, 3)) * BOLT(J, 3)
        S = H * BNU(J) * STIM(J) / (BHYD(J, 3) / BHYD(J, 4) - EHVKT(J))
      END IF
    END IF

    ! Sum over upper levels in window, weighted by occupation probability.
    ! Only the surviving fraction w_m contributes as line opacity.
    ! The dissolved fraction (1 - w_m) is handled by HOP's pseudo-continuum.
    DO M = M1, M2
      BHYDJM = 1.0D0
      IF (M .LE. 6) BHYDJM = BHYD(J, M)

      ! Occupation probability weight for upper level
      w_m = 1.0D0
      IF (M .GE. 2 .AND. M .LE. NMAX_OCC) w_m = W_OCC(J, M)
      IF (M .GT. NMAX_OCC) w_m = 0.0D0

      ! Line opacity (surviving fraction only)
      A = w_m * STARK_MMM(N, M, J) * (1.0D0 - EHVKT(J) * BHYDJM / BHYD(J, N)) * BOLT(J, N)
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J, N) / BHYDJM - EHVKT(J))
    END DO

    AHLINE(J) = H
    IF (H .GT. 0.0D0) THEN
      SHLINE(J) = S / H
    ELSE
      SHLINE(J) = BNU(J)
    END IF
  END DO

  RETURN

END SUBROUTINE HLINOP

!=========================================================================
! FUNCTION STARK(N, M, J)
!
! Stark-broadened hydrogen line profile for transition N -> M at depth J.
!
! Returns the absorption cross-section contribution from this hydrogen
! line at the current frequency FREQ.  Uses the Vidal-Cooper-Smith (VCS)
! unified Stark broadening theory with quasistatic ion and impact electron
! contributions.
!
! Key quantities:
!   F0(J) = 1.25e-9 * n_e^(2/3)  — Holtsmark normal field strength
!   KNM:  transition matrix element (tabulated for M-N <= 5, asymptotic beyond)
!   FNM:  oscillator strength proxy (tabulated for M-N <= 10, extrapolated)
!   BETA: normalized detuning = |FREQ - FREQ_NM| * c / (FREQ_NM^2 * F0 * KNM)
!
! Profile regimes:
!   BETA <= 20: PROF = 8/(80 + beta^3), core/intermediate wings
!   BETA >  20: PROF = 1.5/beta^(5/2), far wings with DIOI correction
!
! Arguments: N = lower level, M = upper level, J = depth index
!=========================================================================

FUNCTION STARK(N, M, J)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, M, J
  REAL(8)  :: STARK

  ! --- Transition matrix element table: KNMTAB(M-N, N) for M-N=1..5, N=1..4 ---
  REAL(8), PARAMETER :: KNMTAB(5,4) = reshape((/ &
    .000356D0, .000523D0, .00109D0, .00149D0, .00225D0, &
    .0125D0,   .0177D0,   .028D0,   .0348D0,  .0493D0, &
    .124D0,    .171D0,    .223D0,   .261D0,   .342D0, &
    .683D0,    .866D0,   1.02D0,   1.19D0,   1.46D0 /), shape(KNMTAB))

  ! --- Oscillator strength proxy: FSTARK(M-N, N) for M-N=1..10, N=1..4 ---
  REAL(8), PARAMETER :: FSTARK(10,4) = reshape((/ &
    .1387D0,  .07910D0, .02126D0,  .01394D0,  .006462D0, &
    .004814D0, .002779D0, .002216D0, .001443D0, .001201D0, &
    .3921D0,  .1193D0,  .03766D0,  .02209D0,  .01139D0, &
    .008036D0, .005007D0, .003850D0, .002658D0, .002151D0, &
    .6103D0,  .1506D0,  .04931D0,  .02768D0,  .01485D0, &
    .01023D0, .006588D0, .004996D0, .003542D0, .002838D0, &
    .8163D0,  .1788D0,  .05985D0,  .03189D0,  .01762D0, &
    .01196D0, .007825D0, .005882D0, .004233D0, .003375D0 /), shape(FSTARK))

  REAL(8), PARAMETER :: RYD = FREQ_RYDH       ! Rydberg frequency (Hz)
  ! PI from mod_constants; CLIGHT_ANGS replaces local CLIGHT (Å·Hz)
  REAL(8), PARAMETER :: A0 = 0.0265384D0       ! profile normalization constant

  ! --- Persistent state ---
  REAL(8),  SAVE :: F0(kw)       ! Holtsmark normal field strength
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Local variables ---
  REAL(8)  :: XN, XM, X, XX, NN, MM
  REAL(8)  :: KNM, FNM, FREQNM, DEL, DBETA, BETA
  REAL(8)  :: Y1, Y2, QSTAT, IMPACT, EXY2
  REAL(8)  :: PROF, RATIO, DIOI
  INTEGER :: K, MMINN

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING STARK'

  !---------------------------------------------------------------------
  ! Recompute Holtsmark field when temperature structure changes
  !---------------------------------------------------------------------
  IF (ITEMP .NE. ITEMP1) THEN
    DO K = 1, NRHOX
      F0(K) = 1.25D-9 * XNE(K)**(2.0D0/3.0D0)
    END DO
    ITEMP1 = ITEMP
  END IF

  !---------------------------------------------------------------------
  ! Transition properties
  !---------------------------------------------------------------------
  XN = dble(N)
  XM = dble(M)
  X  = XN / XM
  XX = X**2
  NN = XN * XN
  MM = XM * XM
  MMINN = M - N

  ! Transition matrix element KNM
  IF (MMINN .LE. 5) THEN
    KNM = KNMTAB(MMINN, N)
  ELSE
    KNM = 5.5D-5 * (NN * MM)**2 / (MM - NN)
  END IF

  ! Oscillator strength proxy FNM
  IF (MMINN .LE. 10) THEN
    FNM = FSTARK(MMINN, N)
  ELSE
    FNM = FSTARK(10, N) * ((20.0D0 * XN + 100.0D0) &
         / ((XN + 10.0D0) * XM * (1.0D0 - XX)))**3
  END IF

  !---------------------------------------------------------------------
  ! Line center, detuning, and normalized detuning
  !---------------------------------------------------------------------
  FREQNM = RYD * (1.0D0 / NN - 1.0D0 / MM)
  DEL = abs(FREQ - FREQNM)
  DBETA = CLIGHT_ANGS / FREQNM**2 / F0(J) / KNM
  BETA = DBETA * DEL

  !---------------------------------------------------------------------
  ! Quasistatic ion + impact electron broadening
  !---------------------------------------------------------------------
  Y1 = MM * DEL * HKT(J) / 2.0D0
  Y2 = (PI**2 / 2.0D0 / A0 / CLIGHT) * DEL**2 / XNE(J)
  QSTAT = 1.5D0 + 0.5D0 * (Y1**2 - 1.384D0) / (Y1**2 + 1.384D0)

  IMPACT = 0.0D0
  IF (Y1 .LT. 8.0D0 .AND. Y1 .LT. Y2) THEN
    EXY2 = 0.0D0
    IF (Y2 .LE. 8.0D0) EXY2 = EXINT(Y2)
    IMPACT = 1.438D0 * sqrt(Y1 * (1.0D0 - XX)) &
           * (0.4D0 * exp(-Y1) + EXINT(Y1) - 0.5D0 * EXY2)
  END IF

  !---------------------------------------------------------------------
  ! Profile and ratio assembly
  !---------------------------------------------------------------------
  IF (BETA .LE. 20.0D0) THEN
    ! Core and intermediate wings
    PROF  = 8.0D0 / (80.0D0 + BETA**3)
    RATIO = QSTAT + IMPACT
  ELSE
    ! Far wings with debye-shielding correction
    PROF = 1.5D0 / BETA / BETA / sqrt(BETA)
    DIOI = 6.28D0 * 1.48D-25 * (2.0D0 * MM * RYD / DEL) * XNE(J) &
         * (sqrt(2.0D0 * MM * RYD / DEL) * (1.3D0 * QSTAT + 0.30D0 * IMPACT) &
          - 3.9D0 * RYD * HKT(J))
    RATIO = QSTAT * min(1.0D0 + DIOI, 1.25D0) + IMPACT
  END IF

  STARK = A0 * FNM * PROF * DBETA * RATIO
  RETURN

CONTAINS

  !-----------------------------------------------------------------------
  ! Exponential integral E1(x) approximation (Abramowitz & Stegun 5.1.53)
  !-----------------------------------------------------------------------
  PURE REAL(8) FUNCTION EXINT(X)
    REAL(8), INTENT(IN) :: X
    EXINT = -log(X) - 0.57516D0 + (0.97996D0 + (-0.21654D0 &
          + (0.033572D0 + (-0.0029222D0 + 1.05439D-4 * X) * X) * X) * X) * X
  END FUNCTION EXINT

END FUNCTION STARK


!=======================================================================
! SUBROUTINE INIT_STARK_TABLES
!
! Read preprocessed Stehlé MMM hydrogen Stark broadening tables from
! binary data files.  Called once at startup.
!
! The data files are produced by the Python preprocessor
! (preprocess_stehle.py) from the original CDS VI/98A tables
! (Stehlé & Hutcheon 1999) and the Brackett tables from
! Stehlé & Fouquet (2010).
!
! Each file contains area-normalized Stark profiles I(Δα) for one
! hydrogen series (Lyman, Balmer, Paschen, or Brackett) on a grid
! of (electron density, temperature, reduced detuning Δα).
!
! References:
!   Stehlé, C. & Hutcheon, R. 1999, A&AS 140, 93
!   Stehlé, C. & Fouquet, S. 2010, Int. J. Spectrosc. 2010, 506346
!=======================================================================

SUBROUTINE INIT_STARK_TABLES

  IMPLICIT NONE

  CHARACTER(len=256) :: filepath
  CHARACTER(len=*), PARAMETER :: SERIES_FILES(4) = &
    (/ 'stehle_lyman.bin   ', 'stehle_balmer.bin  ', &
       'stehle_paschen.bin ', 'stehle_brackett.bin' /)
  INTEGER :: iseries, iu, n_lower, n_upper_min, n_upper_max
  INTEGER :: n_dens, n_temps, n_dalpha, n_trans, ios

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING INIT_STARK_TABLES'

  DO iseries = 1, NSTARK_SERIES

    ! Build file path
    filepath = trim(DATADIR) // '/' // trim(SERIES_FILES(iseries))

    ! Try to open the file
    iu = 200 + iseries
    OPEN(iu, FILE=trim(filepath), STATUS='old', FORM='unformatted', &
         ACCESS='sequential', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(A,A)') '  INIT_STARK_TABLES: file not found: ', trim(filepath)
      STEHLE_DATA(iseries)%loaded = .FALSE.
      CYCLE
    END IF

    ! Record 1: dimensions
    READ(iu) n_lower, n_upper_min, n_upper_max, n_dens, n_temps, n_dalpha
    n_trans = n_upper_max - n_upper_min + 1

    STEHLE_DATA(iseries)%n_lower     = n_lower
    STEHLE_DATA(iseries)%n_upper_min = n_upper_min
    STEHLE_DATA(iseries)%n_upper_max = n_upper_max
    STEHLE_DATA(iseries)%n_transitions = n_trans
    STEHLE_DATA(iseries)%n_dens      = n_dens

    ! Record 2: density grid
    READ(iu) STEHLE_DATA(iseries)%density_grid(1:n_dens)

    ! Record 3: temperature grid
    READ(iu) STEHLE_DATA(iseries)%temp_grid(1:n_temps)

    ! Record 4: log Δα grid
    READ(iu) STEHLE_DATA(iseries)%log_dalpha_grid(1:n_dalpha)

    ! Records 5-6: per-transition metadata
    ALLOCATE(STEHLE_DATA(iseries)%max_dens_idx(n_trans))
    ALLOCATE(STEHLE_DATA(iseries)%k_alpha(n_trans))
    READ(iu) STEHLE_DATA(iseries)%max_dens_idx(1:n_trans)
    READ(iu) STEHLE_DATA(iseries)%k_alpha(1:n_trans)

    ! Record 7: profile array
    ALLOCATE(STEHLE_DATA(iseries)%profiles(n_dalpha, n_temps, n_dens, n_trans))
    READ(iu) STEHLE_DATA(iseries)%profiles

    CLOSE(iu)
    STEHLE_DATA(iseries)%loaded = .TRUE.

  END DO

  STEHLE_TABLES_LOADED = .TRUE.

END SUBROUTINE INIT_STARK_TABLES


!=======================================================================
! FUNCTION hydrogen_f_value(N, M)
!
! Exact hydrogen oscillator strength f_{N→M} for bound-bound
! absorption.  Uses tabulated values for low transitions and the
! Kramers approximation with empirical correction for higher ones.
!
! These are the TRUE oscillator strengths (not FSTARK from the old
! Kurucz profile code, which uses a different normalization).
!=======================================================================

FUNCTION hydrogen_f_value(N, M)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, M
  REAL(8) :: hydrogen_f_value

  ! Tabulated gf = g_n × f_{N→M} for hydrogen
  ! Source: NIST Atomic Spectra Database (Wiese et al.)
  ! g_n = 2N²
  REAL(8), PARAMETER :: GF_LYMAN(29) = (/ &
    0.8324D0, 0.1580D0, 0.05798D0, 0.02787D0, 0.01551D0, &
    0.009466D0, 0.006158D0, 0.004220D0, 0.003014D0, 0.002225D0, &
    0.001688D0, 0.001312D0, 0.001038D0, 0.000835D0, 0.000681D0, &
    0.000562D0, 0.000469D0, 0.000395D0, 0.000336D0, 0.000288D0, &
    0.000249D0, 0.000216D0, 0.000189D0, 0.000166D0, 0.000146D0, &
    0.000129D0, 0.000115D0, 0.000103D0, 0.000092D0 /)

  REAL(8), PARAMETER :: GF_BALMER(28) = (/ &
    5.126D0, 0.9543D0, 0.3571D0, 0.1770D0, 0.1023D0, &
    0.06497D0, 0.04394D0, 0.03104D0, 0.02270D0, 0.01708D0, &
    0.01314D0, 0.01028D0, 0.00818D0, 0.00659D0, 0.00537D0, &
    0.00443D0, 0.00369D0, 0.00310D0, 0.00263D0, 0.00225D0, &
    0.00194D0, 0.00168D0, 0.00146D0, 0.00128D0, 0.00113D0, &
    0.00100D0, 0.000889D0, 0.000793D0 /)

  REAL(8), PARAMETER :: GF_PASCHEN(27) = (/ &
    15.16D0, 2.715D0, 1.001D0, 0.4937D0, 0.2850D0, &
    0.1813D0, 0.1232D0, 0.08804D0, 0.06526D0, 0.04975D0, &
    0.03882D0, 0.03089D0, 0.02498D0, 0.02050D0, 0.01703D0, &
    0.01430D0, 0.01212D0, 0.01036D0, 0.00893D0, 0.00775D0, &
    0.00677D0, 0.00595D0, 0.00526D0, 0.00467D0, 0.00416D0, &
    0.00373D0, 0.00335D0 /)

  REAL(8), PARAMETER :: GF_BRACKETT(3) = (/ &
    33.22D0, 5.731D0, 2.092D0 /)

  INTEGER :: idx
  REAL(8)  :: gf, xn, xm

  xn = dble(N)
  xm = dble(M)

  IF (M .LE. N) THEN
    hydrogen_f_value = 0.0D0
    RETURN
  END IF

  idx = M - N  ! for Lyman: idx = M-1; for Balmer: idx = M-2; etc.

  IF (N .EQ. 1 .AND. idx .LE. 29) THEN
    gf = GF_LYMAN(idx)
  ELSE IF (N .EQ. 2 .AND. idx .LE. 28) THEN
    gf = GF_BALMER(idx)
  ELSE IF (N .EQ. 3 .AND. idx .LE. 27) THEN
    gf = GF_PASCHEN(idx)
  ELSE IF (N .EQ. 4 .AND. idx .LE. 3) THEN
    gf = GF_BRACKETT(idx)
  ELSE
    ! Kramers approximation with empirical correction for higher n
    ! gf ≈ 1.96 × n^2 / (m^3 × (1/n^2 - 1/m^2)^3) × correction
    ! The correction factor accounts for the Gaunt factor
    gf = 1.9603D0 * xn**2 / (xm**3 * (1.0D0/xn**2 - 1.0D0/xm**2)**3)
    ! Apply empirical scale to match exact values at table boundary
    gf = gf * 0.80D0  ! approximate average Gaunt factor
  END IF

  hydrogen_f_value = gf / (2.0D0 * xn**2)

END FUNCTION hydrogen_f_value


!=======================================================================
! FUNCTION STARK_MMM(N, M, J)
!
! Stark-broadened hydrogen line cross-section from Stehlé MMM tables.
!
! Returns the same quantity as the old STARK(N,M,J) function: the
! line absorption cross-section σ(ν) at the current frequency FREQ
! for the hydrogen transition N → M at depth point J.
!
! Uses trilinear interpolation in (log N_e, T, log Δα) on the
! preprocessed Stehlé table grid.
!
! Falls back to the old analytic STARK function if:
!   - Tables not loaded
!   - Transition not in table range
!   - Density or temperature outside grid
!
! References:
!   Stehlé, C. & Hutcheon, R. 1999, A&AS 140, 93
!   Stehlé, C. & Fouquet, S. 2010, Int. J. Spectrosc. 2010, 506346
!=======================================================================

FUNCTION STARK_MMM(N, M, J)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, M, J
  REAL(8) :: STARK_MMM

  ! Physical constants
  REAL(8), PARAMETER :: A0 = 0.0265384D0         ! π e² / (m_e c) [cm²/s]

  ! Local variables
  REAL(8)  :: xn, xm, lambda0, freq_nm, del_freq
  REAL(8)  :: F0, dalpha, log_dalpha, log_ne, ne_eff
  REAL(8)  :: I_dalpha, f_nm
  REAL(8)  :: frac_d, frac_t, frac_a
  REAL(8)  :: v000, v001, v010, v011, v100, v101, v110, v111
  REAL(8)  :: v00, v01, v10, v11, v0, v1
  INTEGER :: iseries, itrans, id1, id2, it1, it2, ia1, ia2
  INTEGER :: nd, nt

  TYPE(stark_series_t), POINTER :: S

  ! ---------------------------------------------------------------
  ! Lazy initialization: load tables on first call
  ! ---------------------------------------------------------------
  IF (.NOT. STEHLE_TABLES_LOADED) THEN
    CALL INIT_STARK_TABLES
  END IF

  ! ---------------------------------------------------------------
  ! Determine which series this transition belongs to
  ! ---------------------------------------------------------------
  IF (N .LT. 1 .OR. N .GT. 4 .OR. M .LE. N) THEN
    STARK_MMM = STARK(N, M, J)
    RETURN
  END IF
  iseries = N

  ! Check if tables are loaded for this series
  IF (.NOT. STEHLE_TABLES_LOADED .OR. .NOT. STEHLE_DATA(iseries)%loaded) THEN
    STARK_MMM = STARK(N, M, J)
    RETURN
  END IF

  S => STEHLE_DATA(iseries)

  ! Check if this transition is in the table range
  IF (M .LT. S%n_upper_min .OR. M .GT. S%n_upper_max) THEN
    STARK_MMM = STARK(N, M, J)
    RETURN
  END IF

  itrans = M - S%n_upper_min + 1

  ! ---------------------------------------------------------------
  ! Density regime dispatch.
  !
  ! Four cases, matching synthe_module's STARK_MMM policy:
  !
  !   (a) XNE > density_grid(n_dens)  (above entire series table)
  !       Fall back to the analytic K-P STARK function.  The upper
  !       levels for the whole series are dissolving and the Stehle
  !       table has nothing to say at these densities.
  !
  !   (b) XNE > density_grid(max_dens_idx(itrans))  (within the grid
  !       but above the per-transition Inglis-Teller density)
  !       Stay in the Stehle path; cap XNE at the per-transition IT
  !       density for both the F0 calculation and the density
  !       bracketing below.  This gives the broadest tabulated
  !       profile for this specific transition and merges smoothly
  !       into the pseudo-continuum region.  The dissolved-fraction
  !       opacity is separately accounted for by the Hummer-Mihalas
  !       occupation-probability formalism in HOP (pseudo-continuum)
  !       and HLINOP (w_m weighting of the bound-bound sum).
  !
  !   (c) density_grid(1) <= XNE <= density_grid(max_dens_idx(itrans))
  !       Normal interpolation on the Stehle grid with ne_eff = XNE.
  !
  !   (d) XNE < density_grid(1)  (below entire grid ~ 1e10 cm^-3)
  !       Stay in the Stehle path; frac_d clamps to 0 (lowest row)
  !       below.  F0 uses the true small XNE so the far-wing
  !       asymptotic expression still degrades gracefully.  Falling
  !       back to K-P here was found to introduce profile
  !       discontinuities at shallow layers (see note further down).
  ! ---------------------------------------------------------------
  IF (XNE(J) .GT. S%density_grid(S%n_dens)) THEN
    STARK_MMM = STARK(N, M, J)
    RETURN
  END IF
  ne_eff = MIN(XNE(J), S%density_grid(S%max_dens_idx(itrans)))

  ! ---------------------------------------------------------------
  ! Compute detuning in Δα units (using capped density)
  ! ---------------------------------------------------------------
  xn = dble(N)
  xm = dble(M)
  lambda0 = RYD_ANG * (xn * xm)**2 / ((xm - xn) * (xm + xn))
  freq_nm = CLIGHT_ANGS / lambda0
  del_freq = abs(FREQ - freq_nm)

  F0 = 1.25D-9 * ne_eff**(2.0D0/3.0D0)
  IF (F0 .LE. 0.0D0) THEN
    STARK_MMM = STARK(N, M, J)
    RETURN
  END IF

  ! Δα = Δλ / F0, and Δλ = λ₀²/c × Δν (for small Δλ)
  dalpha = lambda0**2 / CLIGHT_ANGS * del_freq / F0

  ! ---------------------------------------------------------------
  ! Check bounds and find bracketing indices
  ! ---------------------------------------------------------------
  nd = S%n_dens
  nt = NSTARK_TEMPS

  ! Density bounds.
  ! Below the tabulated grid (ne < density_grid(1) ~ 1e10 cm^-3):
  ! stay in the Stehle path rather than falling back to K-P.  The
  ! frac_d clamp below pins the interpolation to the lowest density
  ! row, and at these densities Stark broadening is weak anyway so
  ! the profile is dominated by the Doppler convolution baked into
  ! the Stehle tables.  Falling back to K-P here was found (by the
  ! synthe_module work) to introduce a discontinuity in the emergent
  ! profile at the K-P nwid/hfwid boundary (~0.12 A from line centre
  ! for H18) for the shallowest solar layers.  The effect is also
  ! very relevant for M-dwarf atmospheres where most layers have
  ! Ne < 1e10 cm^-3.  Matches synthe_module policy.
  ! (Above-grid already handled further up by the K-P fallback.)
  log_ne = LOG10(ne_eff)

  ! Temperature bounds
  IF (T(J) .LT. S%temp_grid(1) * 0.5D0 .OR. &
      T(J) .GT. S%temp_grid(nt) * 2.0D0) THEN
    STARK_MMM = STARK(N, M, J)
    RETURN
  END IF

  ! ---------------------------------------------------------------
  ! Find density bracket (using capped ne_eff, not raw XNE)
  ! ---------------------------------------------------------------
  id1 = 1
  DO id1 = 1, nd - 1
    IF (ne_eff .LE. S%density_grid(id1 + 1)) EXIT
  END DO
  id1 = max(1, min(nd - 1, id1))
  id2 = id1 + 1
  frac_d = (log_ne - LOG10(S%density_grid(id1))) / &
           (LOG10(S%density_grid(id2)) - LOG10(S%density_grid(id1)))
  frac_d = max(0.0D0, min(1.0D0, frac_d))

  ! ---------------------------------------------------------------
  ! Find temperature bracket
  ! ---------------------------------------------------------------
  it1 = 1
  DO it1 = 1, nt - 1
    IF (T(J) .LE. S%temp_grid(it1 + 1)) EXIT
  END DO
  it1 = max(1, min(nt - 1, it1))
  it2 = it1 + 1
  frac_t = (T(J) - S%temp_grid(it1)) / &
           (S%temp_grid(it2) - S%temp_grid(it1))
  frac_t = max(0.0D0, min(1.0D0, frac_t))

  ! ---------------------------------------------------------------
  ! Find Δα bracket (in log space)
  ! ---------------------------------------------------------------
  IF (dalpha .LE. 0.0D0) THEN
    ! At line centre: use first grid point value
    log_dalpha = S%log_dalpha_grid(1)
    ia1 = 1
    ia2 = 1
    frac_a = 0.0D0
  ELSE
    log_dalpha = LOG10(dalpha)
    IF (log_dalpha .LE. S%log_dalpha_grid(1)) THEN
      ia1 = 1
      ia2 = 1
      frac_a = 0.0D0
    ELSE IF (log_dalpha .GE. S%log_dalpha_grid(NSTARK_DALPHA)) THEN
      ! Beyond table: use asymptotic wing K_alpha / Δα^2.5
      f_nm = hydrogen_f_value(N, M)
      I_dalpha = S%k_alpha(itrans) / dalpha**2.5D0
      STARK_MMM = A0 * f_nm * I_dalpha * lambda0**2 / (CLIGHT_ANGS * F0)
      RETURN
    ELSE
      ia1 = 1
      DO ia1 = 1, NSTARK_DALPHA - 1
        IF (log_dalpha .LE. S%log_dalpha_grid(ia1 + 1)) EXIT
      END DO
      ia1 = max(1, min(NSTARK_DALPHA - 1, ia1))
      ia2 = ia1 + 1
      frac_a = (log_dalpha - S%log_dalpha_grid(ia1)) / &
               (S%log_dalpha_grid(ia2) - S%log_dalpha_grid(ia1))
      frac_a = max(0.0D0, min(1.0D0, frac_a))
    END IF
  END IF

  ! ---------------------------------------------------------------
  ! Trilinear interpolation in (density, temperature, Δα)
  ! Interpolate in log(I) for better accuracy in the wings
  ! ---------------------------------------------------------------
  v000 = S%profiles(ia1, it1, id1, itrans)
  v001 = S%profiles(ia2, it1, id1, itrans)
  v010 = S%profiles(ia1, it2, id1, itrans)
  v011 = S%profiles(ia2, it2, id1, itrans)
  v100 = S%profiles(ia1, it1, id2, itrans)
  v101 = S%profiles(ia2, it1, id2, itrans)
  v110 = S%profiles(ia1, it2, id2, itrans)
  v111 = S%profiles(ia2, it2, id2, itrans)

  ! Guard against zero/negative values before taking log
  IF (v000 .LE. 0.0D0 .OR. v001 .LE. 0.0D0 .OR. &
      v010 .LE. 0.0D0 .OR. v011 .LE. 0.0D0 .OR. &
      v100 .LE. 0.0D0 .OR. v101 .LE. 0.0D0 .OR. &
      v110 .LE. 0.0D0 .OR. v111 .LE. 0.0D0) THEN
    ! Linear interpolation fallback if any vertex is zero
    v00 = v000 + frac_a * (v001 - v000)
    v01 = v010 + frac_a * (v011 - v010)
    v10 = v100 + frac_a * (v101 - v100)
    v11 = v110 + frac_a * (v111 - v110)
    v0 = v00 + frac_t * (v01 - v00)
    v1 = v10 + frac_t * (v11 - v10)
    I_dalpha = v0 + frac_d * (v1 - v0)
  ELSE
    ! Log-space interpolation for better wing accuracy
    v000 = LOG(v000); v001 = LOG(v001)
    v010 = LOG(v010); v011 = LOG(v011)
    v100 = LOG(v100); v101 = LOG(v101)
    v110 = LOG(v110); v111 = LOG(v111)

    v00 = v000 + frac_a * (v001 - v000)
    v01 = v010 + frac_a * (v011 - v010)
    v10 = v100 + frac_a * (v101 - v100)
    v11 = v110 + frac_a * (v111 - v110)
    v0 = v00 + frac_t * (v01 - v00)
    v1 = v10 + frac_t * (v11 - v10)
    I_dalpha = EXP(v0 + frac_d * (v1 - v0))
  END IF

  I_dalpha = max(I_dalpha, 0.0D0)

  ! ---------------------------------------------------------------
  ! Convert I(Δα) to cross-section
  !   σ = (πe²/m_ec) × f_nm × I(Δα) × λ₀² / (c × F₀)
  ! ---------------------------------------------------------------
  f_nm = hydrogen_f_value(N, M)
  STARK_MMM = A0 * f_nm * I_dalpha * lambda0**2 / (CLIGHT_ANGS * F0)

END FUNCTION STARK_MMM


!=========================================================================
! SUBROUTINE LINOP1
!
! First-iteration line opacity calculation.
!
! Reads packed integer line records from in-memory LINEDATA array
! converts to physical quantities, computes Voigt line profiles, and
! accumulates line opacity into XLINES(J,NU).
!
! Uses a coarse/fine depth grid strategy for efficiency:
!   Phase 1: Evaluate at every 8th depth point (coarse grid)
!   Phase 2: Fill intermediate depths only where neighboring coarse
!            points showed significant line opacity
!
! Two Voigt profile regimes:
!   ADAMP <= 0.20: pretabulated Voigt (H0TAB/H1TAB/H2TAB polynomial)
!                  with Lorentzian wing approximation for VVOIGT > 10
!   ADAMP >  0.20: full VOIGT function evaluation
!
! Wings extend <=100 frequency steps in each direction until the
! line opacity drops below the local continuum.
!=========================================================================

SUBROUTINE LINOP1

  IMPLICIT NONE

  ! Named constants (derived from mod_constants)
  REAL(8), PARAMETER :: CEN_PREFAC4 = 0.026538D0 / SQRTPI / CLIGHT_NMS
  REAL(8), PARAMETER :: FOURPI_C_INV = 1.0D0 / (FOURPI * CLIGHT_NMS)
  REAL(8), PARAMETER :: LORENTZ_PREFAC = INVSQRTPI
  REAL(8), PARAMETER :: ADAMP_THRESH = 0.20D0
  INTEGER, PARAMETER :: MAX_WING = 100

  ! Local variables
  REAL(8)  :: ELO, CGF, GAMMAR, GAMMAS, GAMMAW
  REAL(8)  :: ADAMP, CENTER, CV, DOPWAVE, VVOIGT, WLVAC4
  REAL(8)  :: TXNXN(kw)
  REAL(8)  :: RATIOLG, START, STOP
  INTEGER(4) :: LINEREC(4)
  INTEGER(4) :: IFJ(kw+1)
  INTEGER :: LINE, J, K, NU, IW, I, IV, NUCONT, IWLOLD, IFLINE, LINEUSED

  IF (IDEBUG .EQ. 1) WRITE(6, '(A)') ' RUNNING LINOP1'

  IFJ(1) = 0
  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)

  ! Initialize line opacity array
  DO NU = 1, NUMNU
    XLINES(:, NU) = 0.0
  END DO

  ! Precompute van der Waals broadening proxy
  TXNXN = (XNF(:, 1) + 0.42D0*XNF(:, 3) + 0.85D0*XNF(:, 841)) &
             * (T / 10000.D0)**0.3D0

  NUCONT = 1
  NU = 1
  IWLOLD = 0
  START = WAVESET(NULO) - 1.0D0
  STOP  = WAVESET(NUHI) + 1.0D0
  LINEUSED = 0

  !---------------------------------------------------------------------
  ! Main line loop (reads from in-memory LINEDATA array)
  !---------------------------------------------------------------------
  DO LINE = 1, NLINES_STORED
    LINEREC = LINEDATA(:, LINE)
    CALL UNPACK_LINEDATA(LINEREC)
    IF (mod(LINE, 100000) .EQ. 1 .AND. ITER .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW

    ! Check wavelength ordering
    ! This will occur e.g., when multiple line lists (lowlines+diatomics)
    ! are concatenated together.  This is not bad, it just slows down the 
    ! subsequent search.  If it happens a lot, it can be a performance hit
    IF (IWL .LT. IWLOLD) THEN
       !write(6, *) IWL, IWLOLD
       !write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
       NUCONT = 1
       NU = 1
    END IF

    ! Advance continuous opacity bin
    DO WHILE (IWL .GE. IWAVETAB(NUCONT))
      NUCONT = NUCONT + 1
      IF (NUCONT .GT. 344) THEN
        WRITE(6, '(A)') ' WARNING: NUCONT > 344, clamping'
        NUCONT = 344
        EXIT
      END IF
    END DO

    ! Bounds check on species index
    NELION = abs(IELION / 10)
    IF (NELION .LT. 1 .OR. NELION .GT. mion) THEN
      IF (IDEBUG .EQ. 1) WRITE(6, '(A,I6,A,I10)') &
        '  LINOP1: NELION=', NELION, ' OOB, LINE=', LINE
      IWLOLD = IWL
      CYCLE
    END IF

    ! Convert wavelength
    WLVAC = exp(IWL * RATIOLG)
    WLVAC4 = WLVAC
    IF (WLVAC .LT. START .OR. WLVAC .GT. STOP) THEN
      IWLOLD = IWL
      CYCLE
    END IF

    ! Advance frequency grid to match line position
    DO WHILE (WLVAC .GE. WAVESET(NU))
      NU = NU + 1
      IF (NU .GE. NUMNU) THEN
        ! Last point may miss blue line wings
        IWLOLD = IWL
        CYCLE
      END IF
    END DO

    ! Convert line parameters
    CGF = CEN_PREFAC4 * WLVAC4 * TABLOG(IGFLOG)
    ELO = TABLOG(IELO)
    ADAMP = 0.0D0
    IFLINE = 0

    !-------------------------------------------------------------------
    ! Phase 1: Coarse depth grid (every 8th depth)
    !-------------------------------------------------------------------
    DO J = 8, NRHOX, 8
      IFJ(J+1) = 0
      CENTER = CGF * XNFDOP(J, NELION)
      IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
      CENTER = CENTER * exp(-ELO * HCKT(J))
      IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
      IFJ(J+1) = 1
      IFLINE = 1

      ! Compute damping (once per line)
      IF (ADAMP .EQ. 0.0D0) THEN
        GAMMAR = TABLOG(IGR) * WLVAC4 * FOURPI_C_INV
        GAMMAS = TABLOG(IGS) * WLVAC4 * FOURPI_C_INV
        GAMMAW = TABLOG(IGW) * WLVAC4 * FOURPI_C_INV
      END IF

      ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
      DOPWAVE = DOPPLE(J, NELION) * WLVAC4

      IF (ADAMP .GT. ADAMP_THRESH) THEN
        ! --- Full Voigt regime ---
        ! Red wing
        DO IW = NU, min(NU + MAX_WING, NUMNU)
          CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
        ! Blue wing
        DO I = 1, MAX_WING
          IW = NU - I
          IF (IW .LE. 0) EXIT
          CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
      ELSE
        ! --- Pretabulated Voigt regime ---
        ! Red wing
        DO IW = NU, min(NU + MAX_WING, NUMNU)
          VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
          IF (VVOIGT .GT. 10.0D0) THEN
            CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
          ELSE
            IV = int(VVOIGT * 200.0D0 + 1.5D0)
            CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
          END IF
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
        ! Blue wing
        DO I = 1, MAX_WING
          IW = NU - I
          IF (IW .LE. 0) EXIT
          VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
          IF (VVOIGT .GT. 10.0D0) THEN
            CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
          ELSE
            IV = int(VVOIGT * 200.0D0 + 1.5D0)
            CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
          END IF
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
      END IF
    END DO  ! coarse depth J

    !-------------------------------------------------------------------
    ! Phase 2: Fine depth grid (intermediate points between coarse)
    !-------------------------------------------------------------------
    DO K = 8, NRHOX, 8
      IF (IFJ(K-7) + IFJ(K+1) .EQ. 0) CYCLE
      DO J = K-7, K-1
        CENTER = CGF * XNFDOP(J, NELION)
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
        CENTER = CENTER * exp(-ELO * HCKT(J))
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE

        ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
        DOPWAVE = DOPPLE(J, NELION) * WLVAC4

        IF (ADAMP .GT. ADAMP_THRESH) THEN
          ! --- Full Voigt regime ---
          ! Red wing
          DO IW = NU, min(NU + MAX_WING, NUMNU)
            CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
          ! Blue wing
          DO I = 1, MAX_WING
            IW = NU - I
            IF (IW .LE. 0) EXIT
            CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
        ELSE
          ! --- Pretabulated Voigt regime ---
          ! Red wing
          DO IW = NU, min(NU + MAX_WING, NUMNU)
            VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
            IF (VVOIGT .GT. 10.0D0) THEN
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            ELSE
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            END IF
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
          ! Blue wing
          DO I = 1, MAX_WING
            IW = NU - I
            IF (IW .LE. 0) EXIT
            VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
            IF (VVOIGT .GT. 10.0D0) THEN
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            ELSE
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            END IF
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
        END IF
      END DO  ! fine depth J
    END DO  ! coarse block K

    IF (IFLINE .EQ. 1) LINEUSED = LINEUSED + 1

    IWLOLD = IWL
  END DO  ! LINE

  RETURN

END SUBROUTINE LINOP1


!=========================================================================
! SUBROUTINE JOSH(IFSCAT, IFSURF)
!
! Feautrier radiative transfer solver.
!
! Solves for the source function S_nu, mean intensity J_nu, Eddington
! flux H_nu, and K-integral K_nu at the current frequency.
!
! Arguments:
!   IFSCAT = 1: solve integral equation for source function (scattering)
!   IFSCAT = 0: set S_nu = S_bar (no scattering)
!   IFSURF = 0: compute J and H at all depths
!   IFSURF = 1: compute surface flux only
!   IFSURF = 2: compute surface specific intensity at all angles
!
! Uses a fixed 51-point optical depth grid (XTAU8) with pre-computed
! operator matrices COEFJ and COEFH for the Feautrier solution.
! Deep layers (tau > 20) use a variable Eddington factor method on
! the physical depth grid.
!=========================================================================

SUBROUTINE JOSH(IFSCAT, IFSURF)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: IFSCAT, IFSURF

  ! Fixed quadrature grid parameters
  INTEGER, PARAMETER :: NXTAU = 51
  INTEGER, PARAMETER :: MAX_ITER = 50  ! convergence ceiling (both iterations)

  ! H (flux) quadrature weights on the 51-point grid
  REAL(8), PARAMETER :: CH_J(NXTAU) = (/ &
    7.15528131D-07, 1.49142693D-06, 1.52106577D-06, 2.98150826D-06, &
    5.33941056D-06, 9.13329677D-06, 1.61715943D-05, 2.97035986D-05, &
    5.33166603D-05, 9.11154202D-05, 1.61084638D-04, 2.95118050D-04, &
    5.27450291D-04, 8.67939554D-04, 1.61498412D-03, 2.50720908D-03, &
    3.20994272D-03, 5.61912498D-03, 8.60872678D-03, 1.04706492D-02, &
    1.33110350D-02, 1.62635669D-02, 1.90288834D-02, 2.18877215D-02, &
    2.36015432D-02, 2.10819542D-02, 1.80345085D-02, 1.64786074D-02, &
    1.49382707D-02, 1.19676525D-02, 9.90213640D-03, 8.17766134D-03, &
    6.11252524D-03, 4.84035723D-03, 3.06078210D-03, 2.40512565D-03, &
    2.01712688D-03, 1.33288081D-03, 7.83530239D-04, 4.31428343D-04, &
    1.76504589D-04, 4.75738016D-05, 1.65963702D-05, 5.41117970D-06, &
    2.08043571D-06, 7.11612643D-07, 8.08788982D-08, 1.95130507D-08, &
    4.33638281D-09, 8.87765583D-10, 3.90236420D-11 /)

  ! K-integral quadrature weights on the 51-point grid
  REAL(8), PARAMETER :: CK_J(NXTAU) = (/ &
    3.57771910D-07, 7.45730404D-07, 7.60575176D-07, 1.49091113D-06, &
    2.67016185D-06, 4.56793896D-06, 8.08956065D-06, 1.48632944D-05, &
    2.66928291D-05, 4.56529851D-05, 8.08134864D-05, 1.48363324D-04, &
    2.66052346D-04, 4.39771306D-04, 8.25088180D-04, 1.29440730D-03, &
    1.67680858D-03, 2.98973685D-03, 4.68314718D-03, 5.84855257D-03, &
    7.64854718D-03, 9.63155832D-03, 1.16419578D-02, 1.38551742D-02, &
    1.54840983D-02, 1.42877987D-02, 1.25930300D-02, 1.17983138D-02, &
    1.09717194D-02, 8.98320694D-03, 7.59950886D-03, 6.38808031D-03, &
    4.86854184D-03, 3.91568616D-03, 2.51398841D-03, 2.00142385D-03, &
    1.70069211D-03, 1.14058319D-03, 6.80292083D-04, 3.80097074D-04, &
    1.57705377D-04, 4.31706540D-05, 1.51795348D-05, 4.98576401D-06, &
    1.92979223D-06, 6.63957223D-07, 7.65236692D-08, 1.84933668D-08, &
    4.12596224D-09, 8.47334369D-10, 3.81791959D-11 /)

  ! Fixed optical depth grid (51 points, 0 to 20)
  REAL(8), PARAMETER :: XTAU8(NXTAU) = (/ &
    0.0D0, 0.0000032D0, 0.0000056D0, 0.00001D0, 0.000018D0, &
    0.000032D0, 0.000056D0, 0.0001D0, 0.00018D0, 0.00032D0, &
    0.00056D0, 0.001D0, 0.0018D0, 0.0032D0, 0.0056D0, &
    0.01D0, 0.016D0, 0.025D0, 0.042D0, 0.065D0, &
    0.096D0, 0.139D0, 0.196D0, 0.273D0, 0.375D0, &
    0.5D0, 0.63D0, 0.78D0, 0.95D0, 1.15D0, &
    1.35D0, 1.6D0, 1.85D0, 2.15D0, 2.45D0, &
    2.75D0, 3.15D0, 3.65D0, 4.25D0, 5.0D0, &
    6.0D0, 7.0D0, 8.0D0, 9.0D0, 10.0D0, &
    11.5D0, 13.0D0, 14.5D0, 16.0D0, 18.0D0, 20.0D0 /)

  ! Cached exp(-tau/mu) for surface intensity (persists between calls)
  REAL(8), SAVE :: EXTAU(NXTAU, 20) = 0.0D0

  ! Local variables
  REAL(8) :: XS(NXTAU), XSBAR(NXTAU), XALPHA(NXTAU), DIAG(NXTAU)
  REAL(8) :: XH(NXTAU), XJS(NXTAU)
  REAL(8) :: XSBAR8(NXTAU), XALPHA8(NXTAU), XS8(NXTAU)
  REAL(8) :: XH8(NXTAU), XJS8(NXTAU)
  REAL(8) :: A(kw), B(kw), C(kw), SNUBAR(kw), CTWO(kw), B2CT(kw), B2CT1(kw)
  REAL(8) :: DELXS, ERRORX, XK, ERROR, SNEW, EXNEW
  REAL(8) :: TANGLE, D, DDDDD, OLD, SUM_VAL
  INTEGER :: J, JJ, K, KK, L, M, MAXJ, MAXJ1, MU, N1, NM1, NMJ, MDUMMY
  INTEGER :: IFERR, IFNEG
  INTEGER :: output_mode

  IF (IDEBUG .EQ. 1) WRITE(6, '(A)') ' RUNNING JOSH'

  !---------------------------------------------------------------------
  ! Compute total opacity, scattering fraction, and thermal source
  !---------------------------------------------------------------------
  ABTOT = ACONT + ALINE + SIGMAC + SIGMAL
  ALPHA = (SIGMAC + SIGMAL) / ABTOT
  SNUBAR = (ACONT*SCONT + ALINE*SLINE) &
             / (ACONT + ALINE)
  CALL INTEG(RHOX, ABTOT, TAUNU, NRHOX, ABTOT(1)*RHOX(1))
  MAXJ = 0

  !---------------------------------------------------------------------
  ! Solve radiative transfer. Output mode depends on IFSURF:
  !   0 → full J, H, K solution
  !   1 → surface flux H(1) only
  !   2 → surface intensity SURFI(mu) via piecewise parabolic
  ! The solver block exits early when the output mode is determined.
  !---------------------------------------------------------------------
  output_mode = 0

  solver: DO   ! single-pass block for structured exit

    !-------------------------------------------------------------------
    ! No-scattering path: S_nu = S_bar
    !-------------------------------------------------------------------
    IF (IFSCAT .EQ. 0) THEN
      SNU = SNUBAR
      IF (IFSURF .EQ. 2) THEN
        output_mode = 70
        EXIT solver
      END IF
      MAXJ = MAP1(TAUNU, SNU, NRHOX, XTAU8, XS8, NXTAU)
      DO L = 1, NXTAU
        XS(L) = XS8(L)
      END DO
      IF (IFSURF .EQ. 1) THEN
        output_mode = 60
        EXIT solver
      END IF
      ALPHA = 0.0D0
    END IF

    !-------------------------------------------------------------------
    ! Scattering solution on 51-point grid (Lambda iteration)
    !-------------------------------------------------------------------
    IF (TAUNU(1) .GT. XTAU8(NXTAU)) MAXJ = 1

    IF (MAXJ .NE. 1) THEN
      MAXJ = MAP1(TAUNU, SNUBAR, NRHOX, XTAU8, XSBAR8, NXTAU)
      MAXJ = MAP1(TAUNU, ALPHA, NRHOX, XTAU8, XALPHA8, NXTAU)

      DO L = 1, NXTAU
        ! Clamp in case of bad interpolation
        XALPHA8(L) = max(XALPHA8(L), 0.0D0)
        XSBAR8(L)  = max(XSBAR8(L), 1.0D-38)
        XALPHA(L) = XALPHA8(L)
        XSBAR(L)  = XSBAR8(L)
        ! Extrapolate if tau grid starts above model surface
        IF (XTAU8(L) .LT. TAUNU(1)) THEN
          XSBAR8(L) = max(SNUBAR(1), 1.0D-38)
          XALPHA8(L) = max(ALPHA(1), 0.0D0)
          XSBAR(L) = XSBAR8(L)
          XALPHA(L) = XALPHA8(L)
        END IF
        XS(L) = XSBAR(L)
        DIAG(L) = 1.0D0 - XALPHA(L) * COEFJ(L, L)
        IF (abs(DIAG(L)) .LT. 1.0D-30) DIAG(L) = sign(1.0D-30, DIAG(L))
        XSBAR(L) = (1.0D0 - XALPHA(L)) * XSBAR(L)
      END DO

      ! Lambda iteration: Gauss-Seidel sweeps (max MAX_ITER iterations)
      BLOCK
        REAL(8)  :: WORST_ERR
        INTEGER :: WORST_K
        DO L = 1, MAX_ITER
          IFERR = 0
          WORST_ERR = 0.0D0
          WORST_K   = 0
          K = NXTAU + 1
          DO KK = 1, NXTAU
            K = K - 1
            DELXS = 0.0D0
            DO M = 1, NXTAU
              DELXS = DELXS + COEFJ(K, M) * XS(M)
            END DO
            DELXS = (DELXS * XALPHA(K) + XSBAR(K) - XS(K)) / DIAG(K)
            ERRORX = abs(DELXS / XS(K))
            IF (ERRORX .GT. 0.00001D0) IFERR = 1
            IF (ERRORX .GT. WORST_ERR) THEN
              WORST_ERR = ERRORX
              WORST_K   = K
            END IF
            XS(K) = max(XS(K) + DELXS, 1.0D-37)
          END DO
          IF (IFERR .EQ. 0) EXIT
        END DO
        IF (IFERR .NE. 0) &
          WRITE(6, '(A,I4,A,1PE12.4,A,I3,A,E9.2)') &
            ' JOSH WARNING: Lambda iter not converged in ', MAX_ITER, &
            ' sweeps  wave=', CLIGHT_NMS/FREQ, '  tau_pt=', WORST_K, '  err=', WORST_ERR
      END BLOCK

      ! Post-iteration dispatch
      IF (IFSURF .EQ. 1) THEN
        output_mode = 60
        EXIT solver
      END IF
      DO M = 1, NXTAU
        XS8(M) = XS(M)
      END DO
      IF (IFSURF .EQ. 2) THEN
        output_mode = 670
        EXIT solver
      END IF
      MDUMMY = MAP1(XTAU8, XS8, NXTAU, TAUNU, SNU, MAXJ)
    END IF  ! MAXJ /= 1

    !-------------------------------------------------------------------
    ! Deep atmosphere: variable Eddington factor on physical grid
    !-------------------------------------------------------------------
    IF (MAXJ .NE. NRHOX) THEN
      MAXJ1 = MAXJ + 1
      IF (MAXJ .EQ. 1) MAXJ1 = 1
      DO J = MAXJ1, NRHOX
        SNU(J) = SNUBAR(J)
      END DO
      M = max(MAXJ - 1, 1)
      NM1 = NRHOX - M + 1
      NMJ = NRHOX - MAXJ + 1

      ! Variable Eddington iteration (max MAX_ITER iterations)
      DO L = 1, MAX_ITER
        ERROR = 0.0D0
        IFNEG = 0

        ! Safety check: negative SNU → reset to Planck function
        DO J = M, NRHOX
          IF (SNU(J) .LE. 0.0D0) THEN
            IFNEG = 1
            DO JJ = M, NRHOX
              SNUBAR(JJ) = BNU(JJ)
              SNU(JJ) = BNU(JJ)
            END DO
            EXIT
          END IF
        END DO

        CALL DERIV(TAUNU(M), SNU(M), HNU(M), NM1)

        ! Safety check: negative HNU → reset to Planck function
        DO J = M, NRHOX
          IF (HNU(J) .LE. 0.0D0) THEN
            IFNEG = 1
            DO JJ = M, NRHOX
              SNUBAR(JJ) = BNU(JJ)
              SNU(JJ) = BNU(JJ)
            END DO
            CALL DERIV(TAUNU(M), SNU(M), HNU(M), NM1)
            EXIT
          END IF
        END DO

        DO J = M, NRHOX
          HNU(J) = HNU(J) / 3.0D0
        END DO
        CALL DERIV(TAUNU(MAXJ), HNU(MAXJ), JMINS(MAXJ), NMJ)
        DO J = MAXJ1, NRHOX
          IF (IFNEG .EQ. 1) JMINS(J) = 0.0D0
          JNU(J) = JMINS(J) + SNU(J)
          SNEW = (1.0D0 - ALPHA(J)) * SNUBAR(J) + ALPHA(J) * JNU(J)
          ERROR = abs(SNEW - SNU(J)) / SNEW + ERROR
          SNU(J) = SNEW
        END DO
        IF (ERROR .LT. 0.00001D0) EXIT
      END DO
      IF (ERROR .GE. 0.00001D0) &
        WRITE(6, '(A,I4,A,1PE10.3,A,0PF10.3)') ' JOSH WARNING: Eddington iteration did not converge in ', &
          MAX_ITER, ' sweeps, err=', ERROR, '  wave=', CLIGHT_NMS/FREQ
    END IF  ! MAXJ /= NRHOX

    !-------------------------------------------------------------------
    ! Post-solution dispatch
    !-------------------------------------------------------------------
    IF (IFSURF .EQ. 2) THEN
      output_mode = 70
      EXIT solver
    END IF
    IF (MAXJ .EQ. 1) THEN
      KNU(1) = JNU(1) / 3.0D0
      RETURN
    END IF

    ! Compute J, H, K from 51-point operator matrices
    DO L = 1, NXTAU
      XJS(L) = -XS(L)
      DO M = 1, NXTAU
        XJS(L) = XJS(L) + COEFJ(L, M) * XS(M)
      END DO
      XJS8(L) = XJS(L)
      XH(L) = 0.0D0
      DO M = 1, NXTAU
        XH(L) = XH(L) + COEFH(L, M) * XS(M)
      END DO
      XH8(L) = XH(L)
    END DO
    MDUMMY = MAP1(XTAU8, XJS8, NXTAU, TAUNU, JMINS, MAXJ)
    MDUMMY = MAP1(XTAU8, XH8, NXTAU, TAUNU, HNU, MAXJ)
    XK = 0.0D0
    DO M = 1, NXTAU
      XK = XK + CK_J(M) * XS(M)
    END DO
    KNU(1) = XK
    DO J = 1, MAXJ
      SNU(J)  = max(SNU(J), 1.0D-38)
      JNU(J) = max(JMINS(J) + SNU(J), 1.0D-38)
    END DO
    RETURN

  END DO solver

  !---------------------------------------------------------------------
  ! Output mode dispatch (reached via EXIT solver)
  !---------------------------------------------------------------------
  SELECT CASE (output_mode)

  CASE (60)
    ! Surface flux: H(1) = sum of CH_J weights times source function
    XH(1) = 0.0D0
    DO M = 1, NXTAU
      XH(1) = XH(1) + CH_J(M) * XS(M)
    END DO
    HNU(1) = XH(1)

  CASE (670)
    ! Surface intensity from 51-point grid (piecewise parabolic)
    CALL PARCOE(XS8, XTAU8, A, B, C, NXTAU)
    N1 = NXTAU - 1
    DO J = 1, NXTAU
      CTWO(J) = C(J) * 2.0D0
      B2CT(J) = B(J) + CTWO(J) * XTAU8(J)
    END DO
    DO J = 1, N1
      B2CT1(J) = B(J) + CTWO(J) * XTAU8(J+1)
    END DO
    ! Compute and cache exp(-tau/mu) if not yet done
    IF (EXTAU(1,1) .EQ. 0.0D0) THEN
      DO MU = 1, NMU
        DO J = 1, NXTAU
          TANGLE = XTAU8(J) / ANGLE(MU)
          IF (TANGLE .LT. 300.0D0) EXTAU(J, MU) = exp(-TANGLE)
        END DO
      END DO
    END IF
    DO MU = 1, NMU
      SURFI(MU) = 0.0D0
      DO J = 1, N1
        IF (EXTAU(J, MU) .EQ. 0.0D0) EXIT
        SURFI(MU) = SURFI(MU) &
          + EXTAU(J, MU) * (XS8(J) + (B2CT(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU)) &
          - EXTAU(J+1, MU) * (XS8(J+1) + (B2CT1(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU))
      END DO
      SURFI(MU) = SURFI(MU) &
        + EXTAU(NXTAU, MU) * (XS8(NXTAU) + (B2CT(NXTAU) + CTWO(NXTAU)*ANGLE(MU)) * ANGLE(MU))
    END DO

  CASE (70)
    ! Surface intensity from physical grid (piecewise parabolic)
    CALL PARCOE(SNU, TAUNU, A, B, C, NRHOX)
    N1 = NRHOX - 1
    CTWO = C * 2.0D0
    B2CT = B + CTWO * TAUNU
    DO J = 1, N1
      B2CT1(J) = B(J) + CTWO(J) * TAUNU(J+1)
    END DO
    DO MU = 1, NMU
      OLD = 1.0D0
      SUM_VAL = 0.0D0
      BLOCK
        LOGICAL :: tangle_done
        tangle_done = .FALSE.
        DO J = 1, N1
          TANGLE = TAUNU(J+1) / ANGLE(MU)
          EXNEW = exp(-TANGLE)
          D = TANGLE - TAUNU(J) / ANGLE(MU)
          IF (D .GT. 0.03D0) THEN
            SUM_VAL = SUM_VAL &
              + OLD * (SNU(J) + (B2CT(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU)) &
              - EXNEW * (SNU(J+1) + (B2CT1(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU))
            IF (TANGLE .GE. 300.0D0) THEN
              SURFI(MU) = SUM_VAL
              tangle_done = .TRUE.
              EXIT
            END IF
          ELSE
            ! Small optical depth increment: Taylor expansion for stability
            DDDDD = 1.0D0
            IF (D .GT. 0.001D0) &
              DDDDD = ((((D/9.0D0 + 1.0D0)*D/8.0D0 + 1.0D0)*D/7.0D0 + 1.0D0) &
                        *D/6.0D0 + 1.0D0)*D/5.0D0 + 1.0D0
            SUM_VAL = SUM_VAL + EXNEW * (SNU(J) + (SNU(J) + B2CT(J)*ANGLE(MU) &
              + (SNU(J) + (B2CT(J) + CTWO(J)*ANGLE(MU))*ANGLE(MU)) &
              * (DDDDD*D/4.0D0 + 1.0D0)*D/3.0D0)*D/2.0D0) * D
          END IF
          OLD = EXNEW
        END DO
        IF (.NOT. tangle_done) THEN
          SURFI(MU) = SUM_VAL &
            + OLD * (SNU(NRHOX) + (B2CT(NRHOX) + CTWO(NRHOX)*ANGLE(MU)) * ANGLE(MU))
        END IF
      END BLOCK
    END DO

  END SELECT
  RETURN

END SUBROUTINE JOSH

!=========================================================================
! SUBROUTINE BLOCKJ
!
! Load J-coefficient matrix for radiative transfer.
!
! Reads the 51x51 matrix of J-coefficients from blockj.dat on the
! first call. These coefficients are used in the formal solution of
! the radiative transfer equation to compute the mean intensity J
! from the source function at each depth point.
!
! The data file contains 2601 REAL*4 values stored in column-major
! order, which are read into the module array COEFJ(51,51).
!=========================================================================

SUBROUTINE BLOCKJ

  IMPLICIT NONE

  REAL(4)  :: CJ(2601)
  INTEGER :: I, IOS
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING BLOCKJ'

  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'blockj.dat', STATUS='OLD', &
         ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'BLOCKJ: ERROR opening ', trim(DATADIR)//'blockj.dat'
      STOP 'BLOCKJ: cannot read J-coefficient data'
    END IF
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')
    READ(89, *) (CJ(I), I=1, 2601)
    CLOSE(89)
    COEFJ = reshape(dble(CJ), (/ 51, 51 /))
    INITIALIZED = .TRUE.
  END IF

  RETURN

END SUBROUTINE BLOCKJ

!=========================================================================
! SUBROUTINE BLOCKH
!
! Load H-coefficient matrix for radiative transfer.
!
! Reads the 51x51 matrix of H-coefficients from blockh.dat on the
! first call. These coefficients are used in the formal solution of
! the radiative transfer equation to compute the Eddington flux H
! from the source function at each depth point.
!
! The data file contains 2601 REAL*4 values stored in column-major
! order, which are read into the module array COEFH(51,51).
!=========================================================================

SUBROUTINE BLOCKH

  IMPLICIT NONE

  REAL(4)  :: CH(2601)
  INTEGER :: I, IOS
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING BLOCKH'

  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'blockh.dat', STATUS='OLD', &
         ACTION='READ', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      WRITE(6,*) 'BLOCKH: ERROR opening ', trim(DATADIR)//'blockh.dat'
      STOP 'BLOCKH: cannot read H-coefficient data'
    END IF
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')
    READ(89, *) (CH(I), I=1, 2601)
    CLOSE(89)
    COEFH = reshape(dble(CH), (/ 51, 51 /))
    INITIALIZED = .TRUE.
  END IF

  RETURN

END SUBROUTINE BLOCKH

!=========================================================================
! SUBROUTINE IONPOTS
!
! Load ionization potentials and compute cumulative sums.
!
! Reads the array POTION(999) of individual ionization potentials from
! ionpots.dat on the first call. Then computes POTIONSUM, the
! cumulative sum of ionization potentials for each element/ion:
!
!   POTIONSUM(neutral) = 0
!   POTIONSUM(ion I)   = POTION(neutral)
!   POTIONSUM(ion II)  = POTION(neutral) + POTION(ion I)
!   etc.
!
! POTIONSUM gives the total energy needed to reach a given ionization
! stage, used in the Saha equation for ionization equilibrium.
!
! Two regimes:
!   IZ = 1-30  (H through Zn): all IZ+1 ionization stages tracked
!   IZ = 31-99 (Ga through Es): only 5 ionization stages (neutral + 4 ions)
!
! Total entries: 495 (light) + 345 (heavy) = 840
!=========================================================================

SUBROUTINE IONPOTS

  IMPLICIT NONE

  INTEGER :: IZ, ION, IOST
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING IONPOTS'

  ! Read ionization potentials on first call
  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'ionpots.dat', STATUS='OLD', &
         ACTION='READ', IOSTAT=IOST)
    IF (IOST .NE. 0) THEN
      WRITE(6,*) 'IONPOTS: ERROR opening ', trim(DATADIR)//'ionpots.dat'
      STOP 'IONPOTS: cannot read ionization potential data'
    END IF
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')
    READ(89, *) (POTION(IZ), IZ=1, 999)
    CLOSE(89)
    INITIALIZED = .TRUE.
  END IF

  ! Compute cumulative ionization potential sums
  NELION = 0

  ! Light elements (Z=1-30): all ionization stages
  DO IZ = 1, 30
    NELION = NELION + 1
    POTIONSUM(NELION) = 0.0D0       ! neutral atom
    DO ION = 2, IZ + 1
      NELION = NELION + 1
      POTIONSUM(NELION) = POTION(NELION - 1) + POTIONSUM(NELION - 1)
    END DO
  END DO

  ! Heavy elements (Z=31-99): 5 stages only (neutral + 4 ions)
  DO IZ = 31, 99
    NELION = NELION + 1
    POTIONSUM(NELION) = 0.0D0       ! neutral atom
    DO ION = 1, 4
      NELION = NELION + 1
      POTIONSUM(NELION) = POTION(NELION - 1) + POTIONSUM(NELION - 1)
    END DO
  END DO

  RETURN

END SUBROUTINE IONPOTS

!=========================================================================
! SUBROUTINE ISOTOPES
!
! Load isotope data and populate the ISOTOPE array.
!
! Reads the isotope mass and abundance fraction data from isotopes.dat
! on the first call. The file contains ISOION(20, 265): for each
! species index (1-265), 10 isotope masses and 10 isotope abundance
! fractions.
!
! The data is then copied into the ISOTOPE(10, 2, N) array indexed
! by ion state N:
!   ISOTOPE(I, 1, N) = isotope mass I     for the element of state N
!   ISOTOPE(I, 2, N) = isotope fraction I for the element of state N
!
! Since isotope data depends only on the element (not ionization
! stage), all ionization stages of the same element receive identical
! copies.
!
! Three regimes:
!   IZ = 1-30:   atoms with IZ+1 ionization stages each
!   IZ = 31-99:  atoms with 5 ionization stages each
!   IZ = 100-265: molecular species (1 state each)
!=========================================================================

SUBROUTINE ISOTOPES

  IMPLICIT NONE

  REAL(4),  SAVE :: ISOION(20, 265)
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  INTEGER :: IZ, ION, N, I, IOST

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING ISOTOPES'

  ! Read isotope data on first call
  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'isotopes.dat', STATUS='OLD', &
         ACTION='READ', IOSTAT=IOST)
    IF (IOST .NE. 0) THEN
      WRITE(6,*) 'ISOTOPES: ERROR opening ', trim(DATADIR)//'isotopes.dat'
      STOP 'ISOTOPES: cannot read isotope data'
    END IF
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')
    READ(89, *) ((ISOION(I, IZ), I=1,20), IZ=1,265)
    CLOSE(89)
    INITIALIZED = .TRUE.
  END IF

  ! Copy isotope data to all ionization stages
  N = 0

  ! Light elements (Z=1-30): IZ+1 ionization stages
  DO IZ = 1, 30
    DO ION = 1, IZ + 1
      N = N + 1
      DO I = 1, 10
        ISOTOPE(I, 1, N) = ISOION(I, IZ)
        ISOTOPE(I, 2, N) = ISOION(I + 10, IZ)
      END DO
    END DO
  END DO

  ! Heavy elements (Z=31-99): 5 ionization stages
  DO IZ = 31, 99
    DO ION = 1, 5
      N = N + 1
      DO I = 1, 10
        ISOTOPE(I, 1, N) = ISOION(I, IZ)
        ISOTOPE(I, 2, N) = ISOION(I + 10, IZ)
      END DO
    END DO
  END DO

  ! Molecular species (IZ=100-265): 1 state each
  DO IZ = 100, 265
    N = N + 1
    DO I = 1, 10
      ISOTOPE(I, 1, N) = ISOION(I, IZ)
      ISOTOPE(I, 2, N) = ISOION(I + 10, IZ)
    END DO
  END DO

  RETURN

END SUBROUTINE ISOTOPES

!=========================================================================
! SUBROUTINE KAPCONT
!
! Tabulate continuous opacity at 343 wavelength points.
!
! Computes the continuous (non-line) opacity at a fixed grid of
! wavelengths from 9.09 to 400000 nm, storing the results in
! TABCONT(J,NU) for later use in line opacity calculations.
!
! For each wavelength point:
!   1. Set FREQ, WAVENO, and thermodynamic quantities (EHVKT, STIM, BNU)
!   2. Call KAPP to compute continuous opacity ACONT and scattering SIGMAC
!   3. Store total opacity: TABCONT = (ACONT + SIGMAC) * 0.001 / STIM
!
! Also computes IWAVETAB (integer wavelength indices for fast lookup)
! and prints a diagnostic table.
!=========================================================================

SUBROUTINE KAPCONT

  IMPLICIT NONE

  INTEGER, PARAMETER :: NWAVE = 343

  ! --- Wavelength grid (nm): 343 points from 9.09 to 400000 nm ---
  REAL(8), PARAMETER :: WBIG(NWAVE) = (/ &
    9.09D0, 9.35D0, 9.61D0, 9.77D0, 9.96D0, 10.2D0, 10.38D0, 10.56D0, &
    10.77D0, 11.04D0, 11.4D0, 11.78D0, 12.13D0, 12.48D0, 12.71D0, 12.84D0, &
    13.05D0, 13.24D0, 13.39D0, 13.66D0, 13.98D0, 14.33D0, 14.72D0, 15.1D0, &
    15.52D0, 15.88D0, 16.2D0, 16.6D0, 17.03D0, 17.34D0, 17.68D0, 18.02D0, &
    18.17D0, 18.61D0, 19.1D0, 19.39D0, 19.84D0, 20.18D0, 20.5D0, 21.05D0, &
    21.62D0, 21.98D0, 22.3D0, 22.68D0, 23.0D0, 23.4D0, 24.0D0, 24.65D0, &
    25.24D0, 25.68D0, 26.0D0, 26.4D0, 26.85D0, 27.35D0, 27.85D0, 28.4D0, &
    29.0D0, 29.6D0, 30.1D0, 30.8D0, 31.8D0, 32.8D0, 33.8D0, 34.8D0, &
    35.7D0, 36.6D0, 37.5D0, 38.5D0, 39.5D0, 40.5D0, 41.4D0, 42.2D0, &
    43.0D0, 44.1D0, 45.1D0, 46.0D0, 47.0D0, 48.0D0, 49.0D0, 50.0D0, &
    50.6D0, 51.4D0, 53.0D0, 55.0D0, 56.7D0, 58.5D0, 60.5D0, 62.5D0, &
    64.5D0, 66.3D0, 68.0D0, 70.0D0, 71.6D0, 73.0D0, 75.0D0, 77.0D0, &
    79.0D0, 81.0D0, 83.0D0, 85.0D0, 87.0D0, 89.0D0, 90.6D0, 92.6D0, &
    96.0D0, 100.0D0, 104.0D0, 108.0D0, 111.5D0, 114.5D0, 118.0D0, 122.0D0, &
    126.0D0, 130.0D0, 134.0D0, 138.0D0, 142.0D0, 146.0D0, 150.0D0, 154.0D0, &
    160.0D0, 165.0D0, 169.0D0, 173.0D0, 177.5D0, 182.0D0, 186.0D0, 190.5D0, &
    195.0D0, 200.0D0, 204.5D0, 208.5D0, 212.5D0, 217.5D0, 222.5D0, 227.5D0, &
    232.5D0, 237.5D0, 242.5D0, 248.0D0, 253.0D0, 257.5D0, 262.5D0, 267.5D0, &
    272.5D0, 277.5D0, 282.5D0, 287.5D0, 295.0D0, 305.0D0, 315.0D0, 325.0D0, &
    335.0D0, 345.0D0, 355.0D0, 362.0D0, 367.0D0, 375.0D0, 385.0D0, 395.0D0, &
    405.0D0, 415.0D0, 425.0D0, 435.0D0, 455.0D0, 465.0D0, 475.0D0, 485.0D0, &
    495.0D0, 505.0D0, 515.0D0, 525.0D0, 535.0D0, 545.0D0, 555.0D0, 565.0D0, &
    575.0D0, 585.0D0, 595.0D0, 605.0D0, 615.0D0, 625.0D0, 635.0D0, 645.0D0, &
    655.0D0, 665.0D0, 675.0D0, 685.0D0, 695.0D0, 705.0D0, 715.0D0, 725.0D0, &
    735.0D0, 745.0D0, 755.0D0, 765.0D0, 775.0D0, 785.0D0, 795.0D0, 805.0D0, &
    815.0D0, 825.0D0, 835.0D0, 845.0D0, 855.0D0, 865.0D0, 875.0D0, 885.0D0, &
    895.0D0, 905.0D0, 915.0D0, 925.0D0, 935.0D0, 945.0D0, 955.0D0, 965.0D0, &
    975.0D0, 985.0D0, 995.0D0, 1012.5D0, 1037.5D0, 1062.5D0, 1087.5D0, 1112.5D0, &
    1137.5D0, 1162.5D0, 1187.5D0, 1212.5D0, 1237.5D0, 1262.5D0, 1287.5D0, 1312.5D0, &
    1337.5D0, 1362.5D0, 1387.5D0, 1412.5D0, 1442.0D0, 1467.0D0, 1487.5D0, 1512.5D0, &
    1537.5D0, 1562.5D0, 1587.5D0, 1620.0D0, 1660.0D0, 1700.0D0, 1740.0D0, 1780.0D0, &
    1820.0D0, 1860.0D0, 1900.0D0, 1940.0D0, 1980.0D0, 2025.0D0, 2075.0D0, 2125.0D0, &
    2175.0D0, 2225.0D0, 2265.0D0, 2290.0D0, 2325.0D0, 2375.0D0, 2425.0D0, 2475.0D0, &
    2525.0D0, 2575.0D0, 2625.0D0, 2675.0D0, 2725.0D0, 2775.0D0, 2825.0D0, 2875.0D0, &
    2925.0D0, 2975.0D0, 3025.0D0, 3075.0D0, 3125.0D0, 3175.0D0, 3240.0D0, 3340.0D0, &
    3450.0D0, 3550.0D0, 3650.0D0, 3750.0D0, 3850.0D0, 3950.0D0, 4050.0D0, 4150.0D0, &
    4250.0D0, 4350.0D0, 4450.0D0, 4550.0D0, 4650.0D0, 4750.0D0, 4850.0D0, 4950.0D0, &
    5050.0D0, 5150.0D0, 5250.0D0, 5350.0D0, 5450.0D0, 5550.0D0, 5650.0D0, 5750.0D0, &
    5850.0D0, 5950.0D0, 6050.0D0, 6150.0D0, 6250.0D0, 6350.0D0, 6500.0D0, 6700.0D0, &
    6900.0D0, 7100.0D0, 7300.0D0, 7500.0D0, 7700.0D0, 7900.0D0, 8100.0D0, 8300.0D0, &
    8500.0D0, 8700.0D0, 8900.0D0, 9100.0D0, 9300.0D0, 9500.0D0, 9700.0D0, 9900.0D0, &
    10000.0D0, 20000.0D0, 40000.0D0, 60000.0D0, 80000.0D0, 100000.0D0, 120000.0D0, 140000.0D0, &
    160000.0D0, 200000.0D0, 240000.0D0, 280000.0D0, 320000.0D0, 360000.0D0, 400000.0D0 /)

  ! --- Local variables ---
  REAL(8)  :: RATIOLG, FREQ15
  INTEGER :: NU, NU9, J, N

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING KAPCONT'

  !---------------------------------------------------------------------
  ! Copy wavelength grid to module array
  !---------------------------------------------------------------------
  WAVETAB(1:NWAVE) = WBIG

  !---------------------------------------------------------------------
  ! Wavelength-to-index conversion factor
  !---------------------------------------------------------------------
  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)
  IWAVETAB(NWAVE + 1) = 2**30    ! sentinel

  !---------------------------------------------------------------------
  ! Compute continuous opacity at each wavelength point
  !---------------------------------------------------------------------
  DO NU = 1, NWAVE
    WAVE = WAVETAB(NU)
    IWAVETAB(NU) = int(log(WAVE) / RATIOLG + 0.5D0)
    FREQ = CLIGHT_NMS / WAVETAB(NU)
    WAVENO = 1.0D7 / WAVETAB(NU)
    FREQLG = log(FREQ)

    ! Set up thermodynamic quantities at each depth
    FREQ15 = FREQ / 1.0D15
    ACONT = 1.0D10
    SCONT = 0.0D0
    EHVKT = exp(-FREQ * HKT)
    STIM = 1.0D0 - EHVKT
    BNU = BNU_PREFAC * FREQ15**3 * EHVKT / STIM

    ! Compute continuous opacity
    N = 0
    IF (WAVE .GT. WAVESET(NULO)) CALL KAPP

    ! Store tabulated opacity
    TABCONT(:, NU) = (ACONT + SIGMAC) * 0.001D0 / STIM
  END DO

  ! Pad last entry
  TABCONT(:, NWAVE + 1) = TABCONT(:, NWAVE)
  WAVETAB(NWAVE + 1) = WAVETAB(NWAVE)

  !---------------------------------------------------------------------
  ! Diagnostic printout
  !---------------------------------------------------------------------

  IF (IDEBUG .EQ. 1) THEN
     DO NU = 1, NWAVE, 10
        NU9 = min(NU + 9, NWAVE)
        WRITE(6, '(5X, 10F12.2)') (WAVETAB(N), N = NU, NU9)
        DO J = 1, NRHOX
           WRITE(6, '(I5, 1P10E12.3)') J, (TABCONT(J, N), N = NU, NU9)
        END DO
     END DO
  ENDIF

  RETURN

END SUBROUTINE KAPCONT

!=========================================================================
! SUBROUTINE SELECTLINES
!
! First-pass line selection filter: reads line databases and selects
! lines strong enough to affect the model atmosphere.
!
! Reads from 7 line databases in sequence:
!   (1) LOWLINES predicted (unit 11) — atomic predicted lines
!   (2) LOWLINES observed  (unit 111) — atomic observed lines
!   (3) HILINES            (unit 21) — atomic high-excitation lines
!   (4) DIATOMICS          (unit 31) — diatomic molecular lines
!   (5) TiO                (unit 41) — titanium oxide lines
!   (6) H2O                (unit 51) — water lines (special format)
!   (7) H3+               (unit 61) — trihydrogen cation lines
!
! Selection criteria (for each line):
!   1. Species must be present: XNFDOPMAX(NELION,NU) > threshold
!   2. Line-center/continuum ratio: CENRATIO = (pi*e^2/m_e/c/sqrt(pi))
!      * gf * max_J(XNFP*DOP/TABCONT) / freq > 1
!   3. Boltzmann excitation: CENRATIO * exp(-E_lo * hc/kT_max) > 1
!
! Selected lines are stored in the in-memory LINEDATA array.
!
! Special handling:
!   - Diatomics: isotope gf correction via MOLCODES/ISOX lookup
!   - TiO: isotope gf corrections for 46-50Ti via select case
!   - H2O: different record format (3 integers); isotope encoded in
!     signs of IELO/IGFLOG; radiation damping computed from frequency
!   - H3+: shares opacity slot NELION=895 with TiO
!=========================================================================

SUBROUTINE SELECTLINES

  IMPLICIT NONE

  INTEGER, PARAMETER :: NMOL = 46
  INTEGER, PARAMETER :: MAX_LINES = 500000000  ! max lines per database

  ! Line-center opacity prefactor: pi*e^2/(m_e*c*sqrt(pi)) in CGS
  REAL(8), PARAMETER :: CEN_PREFAC = 0.026538D0 / SQRTPI

  ! Molecular species codes (for DIATOMICS database lookup)
  INTEGER, PARAMETER :: MOLCODES(NMOL) = (/ &
    8410, 8411, 8460, 8461, 8470, 8471, 8480, 8481, 8482, 8510, &
    8511, 8512, 8530, 8531, 8532, 8580, 8581, 8582, 8583, 8584, &
    8620, 8621, 8622, 8623, 8640, 8641, 8642, 8643, 8680, 8681, &
    8682, 8690, 8691, 8692, 8693, 8700, 8701, 8702, 8703, 8704, &
    8705, 8890, 8891, 8892, 8896, 8960 /)

  ! Isotope gf corrections for diatomic molecules
  ! ISOX(IMOL): additive correction to IGFLOG for each isotopologue
  INTEGER, PARAMETER :: ISOX(NMOL) = (/ &
    0, -4469, -5, -1955, -2, -2444, -1, -3398, -2690, &            ! H2,HDD,CH,13CH,NH,15NH,OH,17OH,18OH
    -105, -996, -947, -35, -1331, -1516, &                          ! MgH isotopes, SiH isotopes
    -13, -2189, -2870, -1681, -4398, &                              ! CaH isotopes
    -1362, -77, -1022, -1626, &                                     ! CrH isotopes
    -1237, -38, -1658, -2553, &                                     ! FeH isotopes
    -10, -1960, -3910, &                                            ! C2 isotopes
    -7, -1957, -2449, -4399, &                                      ! CN isotopes
    -6, -1956, -3403, -5353, -2695, -4645, &                        ! CO isotopes
    -36, -1332, -1517, -2725, &                                     ! SiO isotopes
    -2 /)                                                           ! VO

  ! TiO isotope gf corrections: 46Ti..50Ti (ISO = |IELION| - 8949)
  INTEGER, PARAMETER :: TIO_ISOCORR(5) = (/ -1101, -1138, -131, -1259, -1272 /)

  ! H2O isotope gf corrections (ISO determined from signs of IELO/IGFLOG)
  INTEGER, PARAMETER :: H2O_ISOCORR(4) = (/ -1, -3398, -2690, -5000 /)

  ! --- Local variables ---
  REAL(8)  :: XNFDOPMAX(mion, 344)
  REAL(8)  :: CENRATIO, RATIOLG, GR, tablog8
  INTEGER(4) :: LINEREC(4)
  INTEGER :: NU, J, K, LINE
  INTEGER :: N12, N122, N22, N32, N42, N52, N62, N18
  INTEGER :: MOLCODE, MOLCODEOLD, KGFLOG, ISO, IMOL
  INTEGER :: LINEDATA_CAP, IOS

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING SELECTLINES'

  ! Allocate in-memory line storage (replaces fort.12)
  LINEDATA_CAP = 150000000   ! 150M lines initial (~2.4 GB)
  IF (allocated(LINEDATA)) DEALLOCATE(LINEDATA)
  ALLOCATE(LINEDATA(4, LINEDATA_CAP), STAT=IOS)
  IF (IOS .NE. 0) THEN
    WRITE(6, '(A,F6.2,A)') ' SELECTLINES: cannot allocate ', &
      LINEDATA_CAP * 16.0D0 / 1.0D9, ' GB for LINEDATA'
    STOP 'SELECTLINES: insufficient memory'
  END IF
  NLINES_STORED = 0

  !---------------------------------------------------------------------
  ! Compute maximum XNFDOP/TABCONT ratio over depth for each (ion, nu)
  !---------------------------------------------------------------------
  XNFDOPMAX = 0.0D0
  DO NU = 1, 344
    DO NELION = 1, MION
      DO J = 1, NRHOX
        XNFDOPMAX(NELION, NU) = max(XNFDOPMAX(NELION, NU), &
          XNFDOP(J, NELION) / TABCONT(J, NU))
      END DO
    END DO
  END DO

  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)
  N12 = 0; N122 = 0; N22 = 0; N32 = 0; N42 = 0; N52 = 0; N62 = 0
  NU = 1

  !=====================================================================
  ! (1) LOWLINES predicted (unit 11)
  !=====================================================================
  OPEN(UNIT=11, FILE=trim(DATADIR)//'lowlines_pl.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', &
       ACCESS='STREAM', ERR=669)

  DO LINE = 1, MAX_LINES
    READ(11, END=581) LINEREC
    CALL UNPACK_LINEDATA(LINEREC)
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW

    ! Advance wavelength bin to match line position
    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    END DO

    ! Apply selection filters
    NELION = abs(IELION / 10)
    IF (NELION .LT. 1 .OR. NELION .GT. mion) THEN
      IF (IDEBUG .EQ. 1) WRITE(6, '(A,I6,A,I10)') &
        '  SELECTLINES: NELION=', NELION, ' OOB, LINE=', LINE
      CYCLE
    END IF
    IF (XNFDOPMAX(NELION, NU) .LE. 1.0D-37) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    tablog8 = TABLOG(IELO)
    IF (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N12 = N12 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading LOWLINES predicted'
  STOP 'SELECTLINES: increase MAX_LINES'

  581 WRITE(6, '(I12,A)') N12, ' LINES FROM LOWLINES'
  CLOSE(UNIT=11)

  !=====================================================================
  ! (2) LOWLINES observed (unit 111)
  !=====================================================================
  OPEN(UNIT=111, FILE=trim(DATADIR)//'lowlines_obs.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', &
       ACCESS='STREAM', ERR=669)

  DO LINE = 1, MAX_LINES
    READ(111, END=5819) LINEREC
    CALL UNPACK_LINEDATA(LINEREC)

    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    END DO

    NELION = abs(IELION / 10)
    IF (XNFDOPMAX(NELION, NU) .LE. 1.0D-37) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    tablog8 = TABLOG(IELO)
    IF (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    WLVAC = exp(IWL * RATIOLG)
    N122 = N122 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading LOWLINES observed'
  STOP 'SELECTLINES: increase MAX_LINES'

  5819 WRITE(6, '(I12,A)') N122, ' LINES FROM LOWLINES observed'
  CLOSE(UNIT=111)

  !=====================================================================
  ! (3) HILINES (unit 21)
  !=====================================================================
  669 OPEN(UNIT=21, FILE=trim(DATADIR)//'hilines.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', &
           ACCESS='STREAM', ERR=769)
  NU = 1

  DO LINE = 1, MAX_LINES
    READ(21, END=681) LINEREC
    CALL UNPACK_LINEDATA(LINEREC)

    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    END DO

    NELION = abs(IELION / 10)
    IF (XNFDOPMAX(NELION, NU) .EQ. 0.0D0) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    tablog8 = TABLOG(IELO)
    IF (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N22 = N22 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading HILINES'
  STOP 'SELECTLINES: increase MAX_LINES'

  681 WRITE(6, '(I12,A)') N22, ' LINES FROM HILINES'
  CLOSE(UNIT=21)

  !=====================================================================
  ! (4) DIATOMICS (unit 31)
  !=====================================================================
  769 OPEN(UNIT=31, FILE=trim(DATADIR)//'diatomicspacksrt.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', ERR=869)
  NU = 1
  MOLCODEOLD = 0
  IMOL = 0

  DO LINE = 1, MAX_LINES
    READ(31, END=781) LINEREC
    CALL UNPACK_LINEDATA(LINEREC)

    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    END DO

    MOLCODE = abs(IELION)
    KGFLOG = IGFLOG
    IGS = 1

    ! Look up molecular species code (cache last match)
    IF (MOLCODE .NE. MOLCODEOLD) THEN
      MOLCODEOLD = MOLCODE
      IMOL = 0
      DO K = 1, NMOL
        IF (MOLCODE .EQ. MOLCODES(K)) THEN
          IMOL = K
          EXIT
        END IF
      END DO
      IF (IMOL .EQ. 0) THEN
        WRITE(6, '(9I12)') LINE, LINEREC
        STOP 1
      END IF
    END IF

    ! Apply isotope gf correction
    IGFLOG = max(KGFLOG + ISOX(IMOL), 1)

    NELION = abs(IELION / 10)
    IF (XNFDOPMAX(NELION, NU) .EQ. 0.0D0) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    tablog8 = TABLOG(IELO)
    IF (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    CALL PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N32 = N32 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading DIATOMICS'
  STOP 'SELECTLINES: increase MAX_LINES'

  781 WRITE(6, '(I12,A)') N32, ' LINES FROM DIATOMICS'
  CLOSE(UNIT=31)

  !=====================================================================
  ! (5) TiO (unit 41)
  !=====================================================================
  869 OPEN(UNIT=41, FILE=trim(DATADIR)//'schwenke.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', &
           ACCESS='STREAM', ERR=1869)
  NU = 1

  DO LINE = 1, MAX_LINES
    READ(41, END=881) LINEREC
    CALL UNPACK_LINEDATA(LINEREC)

    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    END DO

    KGFLOG = IGFLOG
    ISO = abs(IELION) - 8949

    ! Isotope gf correction for 46Ti-50Ti
    IF (ISO .GE. 1 .AND. ISO .LE. 5) THEN
      IGFLOG = max(KGFLOG + TIO_ISOCORR(ISO), 1)
    END IF

    NELION = 895
    IF (XNFDOPMAX(NELION, NU) .EQ. 0.0D0) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    tablog8 = TABLOG(IELO)
    IF (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    IGS = 1
    IGW = 9384
    CALL PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N42 = N42 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading TIOLINES'
  STOP 'SELECTLINES: increase MAX_LINES'

  881 WRITE(6, '(I12,A)') N42, ' LINES FROM TIOLINES'
  CLOSE(UNIT=41)

  !=====================================================================
  ! (6) H2O (unit 51) — special 3-integer record format
  !=====================================================================
  1869 OPEN(UNIT=51, FILE=trim(DATADIR)//'h2ofastfix.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', &
            ACCESS='STREAM', ERR=2869)
  NU = 1

  DO LINE = 1, MAX_LINES
    READ(51, END=1881) IWL, IELO, IGFLOG

    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      ! Radiation damping from frequency
      GAMMAR = 2.474D-22 * FREQ**2 * 0.001
      GR = log10(dble(GAMMAR))
      IGR = int(GR * 1000.0D0 + 16384.5D0)
      NU = NU + 1
    END DO

    ! Determine isotope from signs of IELO and IGFLOG
    IF (IELO .GT. 0 .AND. IGFLOG .GT. 0) THEN
      ISO = 1      ! 1H1H16O
    ELSE IF (IELO .GT. 0) THEN
      ISO = 2      ! 1H1H17O
    ELSE IF (IGFLOG .GT. 0) THEN
      ISO = 3      ! 1H1H18O
    ELSE
      ISO = 4      ! 1H2H16O
    END IF

    IELION = -(9399 + ISO)
    ELO = dble(abs(IELO))
    KGFLOG = abs(IGFLOG)
    IGFLOG = max(KGFLOG + H2O_ISOCORR(ISO), 1)

    NELION = 940
    IF (XNFDOPMAX(NELION, NU) .EQ. 0.0D0) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    IF (CENRATIO * exp(-ELO * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    IELO = int(log10(ELO) * 1000.0D0 + 16384.5D0)
    IGS = 1
    IGW = 9384
    CALL PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    IF (mod(LINE, 100000) .EQ. 1 .AND. IDEBUG .EQ. 1) &
      WRITE(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N52 = N52 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading H2OFAST'
  STOP 'SELECTLINES: increase MAX_LINES'

  1881 WRITE(6, '(I12,A)') N52, ' LINES FROM H2OFAST'
  CLOSE(UNIT=51)

  !=====================================================================
  ! (7) H3+ (unit 61)
  !=====================================================================
  2869 OPEN(UNIT=61, FILE=trim(DATADIR)//'h3plus.dat', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', &
            ACCESS='STREAM', ERR=1882)
  NU = 1

  DO LINE = 1, MAX_LINES
    READ(61, END=2881) LINEREC
    CALL UNPACK_LINEDATA(LINEREC)

    DO WHILE (IWL .GE. IWAVETAB(NU))
      FREQ = CLIGHT_NMS / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    END DO

    KGFLOG = IGFLOG
    IGFLOG = max(KGFLOG - 1272, 1)

    NELION = 895
    IF (XNFDOPMAX(NELION, NU) .EQ. 0.0D0) CYCLE
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    IF (CENRATIO .LT. 1.0D0) CYCLE
    tablog8 = TABLOG(IELO)
    IF (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) .LT. 1.0D0) CYCLE

    IGS = 1
    IGW = 9384
    CALL PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    IF (NLINES_STORED .GT. LINEDATA_CAP) THEN
      WRITE(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      STOP 'SELECTLINES: increase LINEDATA_CAP'
    END IF
    LINEDATA(:, NLINES_STORED) = LINEREC
    N62 = N62 + 1
  END DO
  WRITE(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading H3PLUS'
  STOP 'SELECTLINES: increase MAX_LINES'

  2881 WRITE(6, '(I10,A)') N62, ' LINES FROM H3PLUS'
  CLOSE(UNIT=61)

  !=====================================================================
  ! Summary
  !=====================================================================
  1882 N18 = N12 + N122 + N22 + N32 + N42 + N52 + N62
  WRITE(6, '(I12,A)') N18, ' LINES TOTAL'
  WRITE(6, '(I12,A,F6.2,A)') NLINES_STORED, ' LINES STORED IN MEMORY (', &
       NLINES_STORED * 16.0D0 / 1.0D9, ' GB)'
  IF (NLINES_STORED .EQ. 0) THEN
    STOP 'SELECTLINES ERROR: no lines found'
  END IF

  FLUSH(6)
  
END SUBROUTINE SELECTLINES

!=========================================================================
! FUNCTION MAP4(XOLD, FOLD, NOLD, XNEW, FNEW, NNEW)
!
! Piecewise quadratic interpolation with curvature-weighted blending.
!
! Remaps function values FOLD on grid XOLD (NOLD points) onto new grid
! XNEW (NNEW points), storing results in FNEW.  Both grids must be
! monotonically increasing.
!
! For each new point, two local quadratic fits are constructed:
!   - "backward" through points (L-2, L-1, L)
!   - "forward"  through points (L-1, L, L+1)
! where L is the upper bracket index.  These are blended using weights
! proportional to curvature: WT = |C_fwd| / (|C_fwd| + |C_bak|).
! The blended quadratic is evaluated as f(x) = A + (B + C*x)*x.
!
! Edge cases:
!   - Near start (L=2): linear interpolation/extrapolation
!   - Near end (L=NOLD): backward quadratic only
!   - Past end (L>NOLD): linear extrapolation from last two points
!   - Same interval as previous point: reuse coefficients
!   - Adjacent interval: shift (backward = previous forward)
!
! Returns: index of last interval used (LL-1).
!=========================================================================

FUNCTION MAP4(XOLD, FOLD, NOLD, XNEW, FNEW, NNEW)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NOLD, NNEW
  REAL(8),  INTENT(IN)  :: XOLD(NOLD), FOLD(NOLD)
  REAL(8),  INTENT(IN)  :: XNEW(NNEW)
  REAL(8),  INTENT(OUT) :: FNEW(NNEW)
  INTEGER :: MAP4

  ! --- Local variables ---
  REAL(8)  :: A, B, C, D
  REAL(8)  :: ABAC, BBAC, CBAC
  REAL(8)  :: AFOR, BFOR, CFOR
  REAL(8)  :: WT
  INTEGER :: K, L, LL, L1, L2

  L  = 2
  LL = 0
  A  = 0.0D0
  B  = 0.0D0
  C  = 0.0D0
  CBAC = 0.0D0; BBAC = 0.0D0; ABAC = 0.0D0
  CFOR = 0.0D0; BFOR = 0.0D0; AFOR = 0.0D0

  DO K = 1, NNEW

    !-------------------------------------------------------------------
    ! Bracket: find L such that XOLD(L-1) <= XNEW(K) < XOLD(L)
    !-------------------------------------------------------------------
    DO WHILE (L .LE. NOLD)
      IF (XNEW(K) .LT. XOLD(L)) EXIT
      L = L + 1
    END DO

    ! If same interval as last point, reuse coefficients
    IF (L .EQ. LL) THEN
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      CYCLE
    END IF

    !-------------------------------------------------------------------
    ! Past end or at start: linear interpolation/extrapolation
    !-------------------------------------------------------------------
    IF (L .GT. NOLD .OR. L .EQ. 2) THEN
      L = min(NOLD, L)
      C = 0.0D0
      B = (FOLD(L) - FOLD(L-1)) / (XOLD(L) - XOLD(L-1))
      A = FOLD(L) - XOLD(L) * B
      LL = L
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      CYCLE
    END IF

    !-------------------------------------------------------------------
    ! Interior point: construct backward and forward quadratics
    !-------------------------------------------------------------------
    L1 = L - 1

    ! Backward quadratic through (L-2, L-1, L)
    IF (L .LE. LL + 1 .AND. L .NE. 3) THEN
      ! Adjacent interval: shift (backward = previous forward)
      CBAC = CFOR
      BBAC = BFOR
      ABAC = AFOR
    ELSE
      ! Compute from scratch
      L2 = L - 2
      D = (FOLD(L1) - FOLD(L2)) / (XOLD(L1) - XOLD(L2))
      CBAC = FOLD(L) / ((XOLD(L) - XOLD(L1)) * (XOLD(L) - XOLD(L2))) &
           + (FOLD(L2) / (XOLD(L) - XOLD(L2)) &
            - FOLD(L1) / (XOLD(L) - XOLD(L1))) / (XOLD(L1) - XOLD(L2))
      BBAC = D - (XOLD(L1) + XOLD(L2)) * CBAC
      ABAC = FOLD(L2) - XOLD(L2) * D + XOLD(L1) * XOLD(L2) * CBAC
    END IF

    ! At last point: use backward quadratic only
    IF (L .GE. NOLD) THEN
      C = CBAC
      B = BBAC
      A = ABAC
      LL = L
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      CYCLE
    END IF

    ! Forward quadratic through (L-1, L, L+1)
    D = (FOLD(L) - FOLD(L1)) / (XOLD(L) - XOLD(L1))
    CFOR = FOLD(L+1) / ((XOLD(L+1) - XOLD(L)) * (XOLD(L+1) - XOLD(L1))) &
         + (FOLD(L1) / (XOLD(L+1) - XOLD(L1)) &
          - FOLD(L) / (XOLD(L+1) - XOLD(L))) / (XOLD(L) - XOLD(L1))
    BFOR = D - (XOLD(L) + XOLD(L1)) * CFOR
    AFOR = FOLD(L1) - XOLD(L1) * D + XOLD(L) * XOLD(L1) * CFOR

    ! Curvature-weighted blending
    WT = 0.0D0
    IF (abs(CFOR) .NE. 0.0D0) WT = abs(CFOR) / (abs(CFOR) + abs(CBAC))
    A = AFOR + WT * (ABAC - AFOR)
    B = BFOR + WT * (BBAC - BFOR)
    C = CFOR + WT * (CBAC - CFOR)
    LL = L

    FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
  END DO

  MAP4 = LL - 1
  RETURN

END FUNCTION MAP4

!=======================================================================
! TABVOIGT: Pretabulate Voigt profile components
!
! Builds lookup tables for fast Voigt profile evaluation during
! line opacity computation. Tabulates three arrays over N points
! at spacing 1/VSTEPS Doppler widths:
!   H0TAB(i) = exp(-v²)         — Doppler core (H₀)
!   H1TAB(i) = H₁(v)           — first Lorentzian correction
!                                  (interpolated from 81-pt table)
!   H2TAB(i) = (1-2v²)exp(-v²) — second correction (H₂)
! 100 steps per Doppler width gives ~2% accuracy.
!=======================================================================

SUBROUTINE TABVOIGT(VSTEPS, N)

  IMPLICIT NONE
  REAL(8),  INTENT(IN) :: VSTEPS
  INTEGER, INTENT(IN) :: N

  ! Pretabulate Voigt function components for fast line profile evaluation
  ! 100 steps per Doppler width gives ~2% accuracy
  INTEGER, PARAMETER :: NTAB = 81
  REAL(8), PARAMETER :: TABVI(81) = (/ &
    0.0D0, 0.1D0, 0.2D0, 0.3D0, 0.4D0, 0.5D0, 0.6D0, 0.7D0, &
    0.8D0, 0.9D0, 1.0D0, 1.1D0, 1.2D0, 1.3D0, 1.4D0, 1.5D0, &
    1.6D0, 1.7D0, 1.8D0, 1.9D0, 2.0D0, 2.1D0, 2.2D0, 2.3D0, &
    2.4D0, 2.5D0, 2.6D0, 2.7D0, 2.8D0, 2.9D0, 3.0D0, 3.1D0, &
    3.2D0, 3.3D0, 3.4D0, 3.5D0, 3.6D0, 3.7D0, 3.8D0, 3.9D0, &
    4.0D0, 4.2D0, 4.4D0, 4.6D0, 4.8D0, 5.0D0, 5.2D0, 5.4D0, &
    5.6D0, 5.8D0, 6.0D0, 6.2D0, 6.4D0, 6.6D0, 6.8D0, 7.0D0, &
    7.2D0, 7.4D0, 7.6D0, 7.8D0, 8.0D0, 8.2D0, 8.4D0, 8.6D0, &
    8.8D0, 9.0D0, 9.2D0, 9.4D0, 9.6D0, 9.8D0, 10.0D0, 10.2D0, &
    10.4D0, 10.6D0, 10.8D0, 11.0D0, 11.2D0, 11.4D0, 11.6D0, &
    11.8D0, 12.0D0 /)
  REAL(8), PARAMETER :: TABH1(81) = (/ &
    -1.12838D0, -1.10596D0, -1.04048D0, -0.93703D0, -0.80346D0, &
    -0.64945D0, -0.48552D0, -0.32192D0, -0.16772D0, -0.03012D0, &
    0.08594D0, 0.17789D0, 0.24537D0, 0.28981D0, 0.31394D0, &
    0.32130D0, 0.31573D0, 0.30094D0, 0.28027D0, 0.25648D0, &
    0.231726D0, 0.207528D0, 0.184882D0, 0.164341D0, 0.146128D0, &
    0.130236D0, 0.116515D0, 0.104739D0, 0.094653D0, 0.086005D0, &
    0.078565D0, 0.072129D0, 0.066526D0, 0.061615D0, 0.057281D0, &
    0.053430D0, 0.049988D0, 0.046894D0, 0.044098D0, 0.041561D0, &
    0.039250D0, 0.035195D0, 0.031762D0, 0.028824D0, 0.026288D0, &
    0.024081D0, 0.022146D0, 0.020441D0, 0.018929D0, 0.017582D0, &
    0.016375D0, 0.015291D0, 0.014312D0, 0.013426D0, 0.012620D0, &
    0.0118860D0, 0.0112145D0, 0.0105990D0, 0.0100332D0, &
    0.0095119D0, 0.0090306D0, 0.0085852D0, 0.0081722D0, &
    0.0077885D0, 0.0074314D0, 0.0070985D0, 0.0067875D0, &
    0.0064967D0, 0.0062243D0, 0.0059688D0, 0.0057287D0, &
    0.0055030D0, 0.0052903D0, 0.0050898D0, 0.0049006D0, &
    0.0047217D0, 0.0045526D0, 0.0043924D0, 0.0042405D0, &
    0.0040964D0, 0.0039595D0 /)

  INTEGER :: I, IDUM
  REAL(8)  :: VV

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING TABVOIGT'

  DO I = 1, N
    H0TAB(I) = dble(I - 1) / VSTEPS
  END DO
  IDUM = MAP4(TABVI, TABH1, NTAB, H0TAB, H1TAB, N)
  DO I = 1, N
    VV = (dble(I - 1) / VSTEPS)**2
    H0TAB(I) = exp(-VV)
    H2TAB(I) = H0TAB(I) - (VV + VV) * H0TAB(I)
  END DO
  RETURN

END SUBROUTINE TABVOIGT

!=========================================================================
! SUBROUTINE XLINOP
!
! Extended line opacity calculation for subsequent iterations.
!
! Reads pre-computed line data from unit 19 and accumulates line opacity
! into XLINES(J,NU).  Handles five line types plus merged continuum:
!
!   TYPE =  0: Normal line (Voigt profile, coarse/fine depth grid)
!   TYPE =  3: PRD line (treated as normal)
!   TYPE = -1: Hydrogen line (Stark-broadened via HPROF4)
!   TYPE =  1: Autoionizing line (Fano profile)
!   TYPE =  2: Coronal line (skipped)
!   else:      Merged continuum (flat opacity to dissolution limit)
!
! Normal lines use the same coarse/fine depth grid strategy as LINOP1
! but with a 2000-step wing limit and a blue-wing cutoff at the
! relevant continuum edge (WCON).
!=========================================================================

SUBROUTINE XLINOP

  IMPLICIT NONE

  ! Named constants
  REAL(8), PARAMETER :: LORENTZ_PREFAC = 0.5642D0  ! 1/sqrt(pi)
  REAL(8), PARAMETER :: ADAMP_THRESH = 0.20D0
  INTEGER, PARAMETER :: MAX_WING = 2000

  ! Ionization edges (cm^-1): CONTX(edge_index, species)
  ! 25 edges × 16 species. Column 1 = H, 4 = He I, 6 = C I, etc.
  ! Column 1 (H I) edges come from H_SERIES_LIMITS; 15 trailing zeros pad the column.
  REAL(8), PARAMETER :: CONTX_COL1_ZEROS(15) = 0.0D0     ! pads the H column
  REAL(8), PARAMETER :: CONTX(25,16) = reshape( (/ &
    H_SERIES_LIMITS(1:10), CONTX_COL1_ZEROS,                                              &
    0.0D0, 198310.76D0, 38454.691D0, 32033.214D0, 29223.753D0, 27175.76D0, 15073.868D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 438908.85D0, 109726.529D0, 48766.491D0, 27430.925D0, 17555.715D0, 12191.437D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 90883.84D0, 90867.42D0, 90840.42D0, 90820.42D0, 90804.0D0, &
    90777.0D0, 80691.18D0, 80627.76D0, 69235.82D0, 69172.4D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 61671.02D0, 39820.615D0, 39800.556D0, &
    39759.842D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 48278.37D0, &
    48166.309D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 66035.0D0, 65957.885D0, 65811.843D0, 65747.55D0, 65670.435D0, 65524.393D0, 59736.15D0, &
    59448.7D0, 50640.63D0, 50553.18D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), (/ 25, 16 /) )

  ! Local variables
  REAL(8)  :: ELO, CGF, GAMMAR, GAMMAS, GAMMAW
  REAL(8)  :: ADAMP, CENTER, CV, DOPWAVE, VVOIGT, WLVAC4, GF, G
  REAL(8)  :: XSECTG, CON, FRELIN, EPSIL, ASHORE, BSHORE
  REAL(8)  :: TXNXN(kw)
  REAL(8)  :: BOLTH(kw, 100), EH(100)
  REAL(8)  :: WCON, WMERGE, WSHIFT, EMERGE(kw), Z, WMAX
  REAL(8)  :: NSTARK, NMERGE, RATIOLG
  ! Temporaries matching the binary record layout of nltelines_obs.bin
  REAL(4)  :: ELO4, GF4, GAMMAR4, GAMMAS4, GAMMAW4
  INTEGER(4) :: IFJ(kw+1)
  INTEGER :: TYPE, NLAST
  INTEGER :: LINE, J, K, N, NU, IW, I, IV, NUCONT, IOS_RD
  INTEGER :: NBLO, NBUP, NCON, NELIONX, LIM

  IF (IDEBUG .EQ. 1) WRITE(6, '(A)') ' RUNNING XLINOP'
  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)

  !---------------------------------------------------------------------
  ! Precompute depth-dependent quantities
  !---------------------------------------------------------------------
  DO J = 1, NRHOX
    TXNXN(J) = (XNF(J,1) + 0.42D0*XNF(J,3) + 0.85D0*XNF(J,841)) &
               * (T(J) / 10000.D0)**0.3D0

    ! Stark dissolution quantum number (empirical)
    NSTARK = 1600.0D0 / XNE(J)**(2.0D0/15.0D0)
    NMERGE = NSTARK

    ! Hydrogen energy levels (cm^-1)
    EH(1) = 0.0D0
    EH(2) = 82259.105D0
    EH(3) = 97492.302D0
    EH(4) = 102823.893D0
    EH(5) = 105291.651D0
    EH(6) = 106632.160D0
    EH(7) = 107440.444D0
    EH(8) = 107965.051D0
    EH(9) = 108324.720D0
    EH(10) = 108581.988D0
    DO N = 11, 100
      EH(N) = ELIM_HI - RYDBERG_H / dble(N)**2
    END DO

    ! Hydrogen Boltzmann factors × number density
    DO NBLO = 1, 100
      BOLTH(J, NBLO) = exp(-EH(NBLO) * HCKT(J)) * XNFDOP(J, 1)
    END DO

    ! Merge energy: dissolved levels above this are continuum
    EMERGE(J) = RYDBERG_INF / NMERGE**2
  END DO

  !---------------------------------------------------------------------
  ! Main line loop (unit 19)
  !---------------------------------------------------------------------
  REWIND 19
  NUCONT = 1
  NU = 1
  IFJ(1) = 0    ! NOTE: was uninitialized in original code

  DO LINE = 1, 500000
    READ(19, IOSTAT=IOS_RD) WLVAC, ELO4, GF4, NBLO, NBUP, NELION, TYPE, &
                      NCON, NELIONX, GAMMAR4, GAMMAS4, GAMMAW4, IWL, LIM
    IF (IOS_RD .NE. 0) RETURN   ! end-of-file or read error
    ELO = dble(ELO4)
    GF  = dble(GF4)
    GAMMAR = dble(GAMMAR4)
    GAMMAS = dble(GAMMAS4)
    GAMMAW = dble(GAMMAW4)
    CGF = GF
    G = GF
    NLAST = TYPE
    ASHORE = GAMMAS
    BSHORE = GAMMAW

    IF (WLVAC .GT. WAVESET(NUHI)) RETURN
    WLVAC4 = WLVAC

    ! Advance continuous opacity bin
    DO WHILE (IWL .GE. IWAVETAB(NUCONT))
      NUCONT = NUCONT + 1
      IF (NUCONT .GT. 344) THEN
        WRITE(6, '(A)') ' WARNING: NUCONT > 344, clamping'
        NUCONT = 344
        EXIT
      END IF
    END DO

    ! Advance frequency grid
    DO WHILE (WLVAC .GE. WAVESET(NU))
      NU = NU + 1
      IF (NU .GE. NUMNU) RETURN
    END DO

    !-----------------------------------------------------------------
    ! Dispatch on line type
    !-----------------------------------------------------------------
    SELECT CASE (TYPE)

    CASE (2)
      ! CORONAL LINE — skip
      CYCLE

    CASE (0, 3)
      ! NORMAL LINE (and PRD treated as normal)
       IF (mod(LINE, 1000) .EQ. 1.AND.IDEBUG .EQ. 1) &
            WRITE(6, '(2I10,2F12.6,1P5E12.3,I10)') &
            LINE, NU, WLVAC, WAVESET(NU), CGF, ELO, GAMMAR, GAMMAS, GAMMAW, NELION

      ! --- Coarse depth grid (every 8th depth) ---
      DO J = 8, NRHOX, 8
        IFJ(J+1) = 0
        CENTER = CGF * XNFDOP(J, NELION)
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
        CENTER = CENTER * exp(-ELO * HCKT(J))
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
        IFJ(J+1) = 1

        ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
        DOPWAVE = DOPPLE(J, NELION) * WLVAC4

        ! Blue-wing cutoff at continuum edge
        WCON = 0.0D0
        IF (NCON .GT. 10) NCON = 0
        IF (NCON .GT. 0) WCON = 1.0D7 / (CONTX(NCON, NELIONX) - EMERGE(J))
        IF (WLVAC .LT. WCON) CYCLE

        IF (ADAMP .GT. ADAMP_THRESH) THEN
          ! Red wing — full Voigt
          DO IW = NU, min(NU + MAX_WING, NUMNU)
            CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
          ! Blue wing — full Voigt
          DO I = 1, MAX_WING
            IW = NU - I
            IF (IW .LE. 0) EXIT
            IF (WAVESET(IW) .LT. WCON) EXIT
            CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
        ELSE
          ! Red wing — pretabulated
          DO IW = NU, min(NU + MAX_WING, NUMNU)
            VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
            IF (VVOIGT .GT. 10.0D0) THEN
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            ELSE
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            END IF
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
          ! Blue wing — pretabulated
          DO I = 1, MAX_WING
            IW = NU - I
            IF (IW .LE. 0) EXIT
            IF (WAVESET(IW) .LT. WCON) EXIT
            VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
            IF (VVOIGT .GT. 10.0D0) THEN
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            ELSE
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            END IF
            XLINES(J, IW) = XLINES(J, IW) + CV
            IF (CV .LT. TABCONT(J, NUCONT)) EXIT
          END DO
        END IF
      END DO  ! coarse depth J

      ! --- Fine depth grid (intermediate points) ---
      DO K = 8, NRHOX, 8
        IF (IFJ(K-7) + IFJ(K+1) .EQ. 0) CYCLE
        DO J = K-7, K-1
          CENTER = CGF * XNFDOP(J, NELION)
          IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
          CENTER = CENTER * exp(-ELO * HCKT(J))
          IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE

          ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
          DOPWAVE = DOPPLE(J, NELION) * WLVAC4

          WCON = 0.0D0
          IF (NCON .GT. 10) NCON = 0
          IF (NCON .GT. 0) WCON = 1.0D7 / (CONTX(NCON, NELIONX) - EMERGE(J))
          IF (WLVAC .LT. WCON) CYCLE

          IF (ADAMP .GT. ADAMP_THRESH) THEN
            ! Red wing — full Voigt
            DO IW = NU, min(NU + MAX_WING, NUMNU)
              CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
              XLINES(J, IW) = XLINES(J, IW) + CV
              IF (CV .LT. TABCONT(J, NUCONT)) EXIT
            END DO
            ! Blue wing — full Voigt
            DO I = 1, MAX_WING
              IW = NU - I
              IF (IW .LE. 0) EXIT
              IF (WAVESET(IW) .LT. WCON) EXIT
              CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
              XLINES(J, IW) = XLINES(J, IW) + CV
              IF (CV .LT. TABCONT(J, NUCONT)) EXIT
            END DO
          ELSE
            ! Red wing — pretabulated
            DO IW = NU, min(NU + MAX_WING, NUMNU)
              VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
              IF (VVOIGT .GT. 10.0D0) THEN
                CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
              ELSE
                IV = int(VVOIGT * 200.0D0 + 1.5D0)
                CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
              END IF
              XLINES(J, IW) = XLINES(J, IW) + CV
              IF (CV .LT. TABCONT(J, NUCONT)) EXIT
            END DO
            ! Blue wing — pretabulated
            DO I = 1, MAX_WING
              IW = NU - I
              IF (IW .LE. 0) EXIT
              IF (WAVESET(IW) .LT. WCON) EXIT
              VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
              IF (VVOIGT .GT. 10.0D0) THEN
                CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
              ELSE
                IV = int(VVOIGT * 200.0D0 + 1.5D0)
                CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
              END IF
              XLINES(J, IW) = XLINES(J, IW) + CV
              IF (CV .LT. TABCONT(J, NUCONT)) EXIT
            END DO
          END IF
        END DO  ! fine depth J
      END DO  ! coarse block K

    CASE (-1)
      ! HYDROGEN LINE — Stark-broadened profile, all depths
      DO J = 1, NRHOX
        CENTER = CGF * BOLTH(J, NBLO)
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
        WCON = 1.0D7 / (CONTX(NCON, 1) - EMERGE(J))
        ! Red wing
        DO IW = NU, min(NU + MAX_WING, NUMNU)
          IF (WAVESET(IW) .LT. WCON) CYCLE
          CV = CENTER * HPROF4(NBLO, NBUP, J, WAVESET(IW) - WLVAC)
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
        ! Blue wing
        DO I = 1, MAX_WING
          IW = NU - I
          IF (IW .LE. 0) EXIT
          IF (WAVESET(IW) .LT. WCON) EXIT
          CV = CENTER * HPROF4(NBLO, NBUP, J, WAVESET(IW) - WLVAC)
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
      END DO

    CASE (1)
      ! AUTOIONIZING LINE — Fano profile, all depths
      FRELIN = CLIGHT_NMS / WLVAC
      DO J = 1, NRHOX
        CENTER = BSHORE * G * XNFP(J, NELION) / RHO(J)
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
        CENTER = CENTER * exp(-ELO * HCKT(J))
        IF (CENTER .LT. TABCONT(J, NUCONT)) CYCLE
        ! Red wing
        DO IW = NU, min(NU + MAX_WING, NUMNU)
          EPSIL = 2.0D0 * (CLIGHT_NMS / WAVESET(IW) - FRELIN) / GAMMAR
          CV = CENTER * (ASHORE * EPSIL + BSHORE) / (EPSIL**2 + 1.0D0) / BSHORE
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
        ! Blue wing
        DO I = 1, MAX_WING
          IW = NU - I
          IF (IW .LE. 0) EXIT
          EPSIL = 2.0D0 * (CLIGHT_NMS / WAVESET(IW) - FRELIN) / GAMMAR
          CV = CENTER * (ASHORE * EPSIL + BSHORE) / (EPSIL**2 + 1.0D0) / BSHORE
          XLINES(J, IW) = XLINES(J, IW) + CV
          IF (CV .LT. TABCONT(J, NUCONT)) EXIT
        END DO
      END DO

    CASE DEFAULT
      ! MERGED CONTINUUM — flat opacity from line to dissolution limit
      Z = 1.0D0
      IF (NELION .EQ. 4) Z = 2.0D0
      WSHIFT = 1.0D7 / (1.0D7 / WLVAC - RYDBERG_INF * Z**2 / dble(NLAST)**2)
      XSECTG = GF
      IF (IDEBUG .EQ. 1) &
           WRITE(6, '(2I10,2F12.6,1P5E12.3,I10)') &
           LINE, NU, WLVAC, WAVESET(NU), CGF, ELO, GAMMAR, GAMMAS, GAMMAW, NELION
      DO J = 1, NRHOX
        WMERGE = 1.0D7 / (1.0D7 / WLVAC - EMERGE(J) * Z**2)
        WMAX = max(WMERGE, WSHIFT)
        CON = XSECTG * XNFP(J, NELION) * exp(-ELO * HCKT(J)) / RHO(J)
        DO IW = NU, NU + 1000
          IF (WMAX .LT. WAVESET(IW)) EXIT
          XLINES(J, IW) = XLINES(J, IW) + CON
        END DO
      END DO

    END SELECT
  END DO  ! LINE

  RETURN

END SUBROUTINE XLINOP

!=======================================================================
! HFNM: Hydrogen oscillator strength f(n,m)
!
! Approximate oscillator strength for the hydrogen transition
! n → m (m > n) using the formula of Menzel & Pekeris with
! empirical corrections for Gaunt factor behavior at large and
! small Δn. Results cached: N-dependent quantities (GINF, GCA,
! FKN, WTC) recomputed only when N changes, and the full f-value
! only when M changes.
!=======================================================================

FUNCTION HFNM(N, M)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, M
  REAL(8) :: HFNM

  ! Hydrogen oscillator strength f(n,m) for transition n → m
  ! Uses approximate formula with caching for repeated N and M values
  REAL(8),  SAVE :: GINF, GCA, FKN, WTC, FNM
  INTEGER, SAVE :: NSTR = 0, MSTR = 0
  REAL(8)  :: XN, XM, XMN, XMN12, FK, WT

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HFNM'

  HFNM = 0.0D0
  IF (M .LE. N) RETURN

  ! Recompute N-dependent quantities if N changed
  IF (N .NE. NSTR) THEN
    XN = dble(N)
    GINF = 0.2027D0 / XN**0.71D0
    GCA = 0.124D0 / XN
    FKN = XN * 1.9603D0
    WTC = 0.45D0 - 2.4D0 / XN**3 * (XN - 1.0D0)
    NSTR = N
    MSTR = 0   ! force M recomputation
  END IF

  ! Recompute M-dependent quantities if M changed
  IF (M .NE. MSTR) THEN
    XM = dble(M)
    XMN = dble(M - N)
    FK = FKN * (XM / (XMN * (XM + dble(N))))**3
    XMN12 = XMN**1.2D0
    WT = (XMN12 - 1.0D0) / (XMN12 + WTC)
    FNM = FK * (1.0D0 - WT * GINF - (0.222D0 + GCA / XM) * (1.0D0 - WT))
    MSTR = M
  END IF

  HFNM = FNM
  RETURN

END FUNCTION HFNM

!=======================================================================
! VCSE1F: E₁(x) exponential integral for Stark line wing profiles
!
! Approximation to E_1(x) = ∫_x^∞ (e^{-t}/t) dt used in the
! Vidal-Cooper-Smith (VCS) theory of hydrogen line broadening.
! Three regimes optimized for speed:
!   x ≤ 0.01:  leading series terms -ln(x) - γ + x
!   x ≤ 1.0:   polynomial (Abramowitz & Stegun 5.1.53)
!   x ≤ 30:    rational approximation (A&S 5.1.56)
!   x > 30:    returns 0 (negligible)
!=======================================================================
FUNCTION VCSE1F(X)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: X
  REAL(8) :: VCSE1F

  ! E_1(x) exponential integral approximation (Abramowitz & Stegun style)
  ! Three regimes: small x (series), medium x (polynomial), large x (asymptotic)
  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING VCSE1F'

  IF (X .LE. 0.0D0) THEN
    VCSE1F = 0.0D0
  ELSE IF (X .LE. 0.01D0) THEN
    ! Small x: leading terms of series expansion
    VCSE1F = -log(X) - 0.577215D0 + X
  ELSE IF (X .LE. 1.0D0) THEN
    ! Medium x: polynomial approximation (Abramowitz & Stegun 5.1.53)
    VCSE1F = -log(X) - 0.57721566D0 + X * (0.99999193D0 &
           + X * (-0.24991055D0 + X * (0.05519968D0 &
           + X * (-0.00976004D0 + X * 0.00107857D0))))
  ELSE IF (X .LE. 30.0D0) THEN
    ! Large x: rational approximation (Abramowitz & Stegun 5.1.56)
    VCSE1F = (X * (X + 2.334733D0) + 0.25062D0) &
           / (X * (X + 3.330657D0) + 1.681534D0) / X * exp(-X)
  ELSE
    VCSE1F = 0.0D0
  END IF
  RETURN

END FUNCTION VCSE1F

!=========================================================================
! FUNCTION STARK_PROFILE(B, P, N, M)
!
! Quasistatic Stark profile wing function S(beta, P) for hydrogen lines.
! (Formerly SOFBET = "S of beta".)
!
! Generates the quasistatic ion microfield broadening profile as a
! function of normalized detuning B (= beta) and Debye shielding
! parameter P.  The alpha and beta lines of the first three series
! are explicitly tabulated; all other transitions use the Balmer-18
! (H2-18) profile as a generic template.
!
! Profiles are renormalized to full oscillator strength.
!
! Profile regimes:
!   B <= 25.12: bilinear interpolation in (P, beta) correction table,
!               blended core (B<=10) and wing (B>=8) approximate profiles
!   25.12 < B <= 500: asymptotic wing with C/D correction
!   B > 500: pure asymptotic wing (1.5/sqrt(B) + 27/B^2) / B^2
!
! Tables read from stark_profile.dat on first call.
!
! NOTE (bug fix): The EQUIVALENCE statements linking PROB1-7 to PROPBM
! and C1-7/D1-7 to C/D were lost in an earlier translation pass, leaving
! PROPBM/C/D uninitialized.  Fixed by reading all tables from file.
!=========================================================================

FUNCTION STARK_PROFILE(B, P, N, M)

  IMPLICIT NONE

  REAL(8),  INTENT(IN) :: B, P
  INTEGER, INTENT(IN) :: N, M
  REAL(8)  :: STARK_PROFILE

  ! --- Grid points ---
  REAL(8), PARAMETER :: PP(5) = (/ 0.0D0, 0.2D0, 0.4D0, 0.6D0, 0.8D0 /)
  REAL(8), PARAMETER :: BETAGRID(15) = (/ &
    1.0D0, 1.259D0, 1.585D0, 1.995D0, 2.512D0, 3.162D0, 3.981D0, &
    5.012D0, 6.310D0, 7.943D0, 10.0D0, 12.59D0, 15.85D0, 19.95D0, 25.12D0 /)

  ! --- Tables (read from file on first call) ---
  REAL(8),  SAVE :: PROPBM(5, 15, 7)   ! correction table
  REAL(8),  SAVE :: C(5, 7)             ! asymptotic correction coeff
  REAL(8),  SAVE :: D(5, 7)             ! asymptotic correction coeff
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! --- Local variables ---
  REAL(8)  :: CORR, B2, SB
  REAL(8)  :: WTPP, WTPM, WTBP, WTBM, CBP, CBM
  REAL(8)  :: PR1, PR2, WT, CC, DD
  INTEGER :: INDX, MMN, IM, IP, J, JM, JP
  INTEGER :: I, K
  CHARACTER(256) :: LINE

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING STARK_PROFILE'

  !---------------------------------------------------------------------
  ! Read tables from file on first call
  !---------------------------------------------------------------------
  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'stark_profile.dat', &
         STATUS='OLD', ACTION='READ')
    DO K = 1, 7
      ! Skip comment lines
      DO
        READ(89, '(A)') LINE
        IF (LINE(1:1) .NE. '#') THEN
          BACKSPACE(89)
          EXIT
        END IF
      END DO
      ! Read PROPBM block: 15 rows of 5 values
      DO I = 1, 15
        READ(89, *) PROPBM(:, I, K)
      END DO
      ! Skip C comment
      DO
        READ(89, '(A)') LINE
        IF (LINE(1:1) .NE. '#') THEN
          BACKSPACE(89)
          EXIT
        END IF
      END DO
      READ(89, *) C(:, K)
      ! Skip D comment
      DO
        READ(89, '(A)') LINE
        IF (LINE(1:1) .NE. '#') THEN
          BACKSPACE(89)
          EXIT
        END IF
      END DO
      READ(89, *) D(:, K)
    END DO
    CLOSE(89)
    INITIALIZED = .TRUE.
  END IF

  !---------------------------------------------------------------------
  ! Select profile index
  !---------------------------------------------------------------------
  INDX = 7                              ! default: generic (Balmer 18)
  MMN = M - N
  IF (N .LE. 3 .AND. MMN .LE. 2) INDX = 2 * (N - 1) + MMN

  CORR = 1.0D0
  B2 = B * B
  SB = sqrt(B)

  !---------------------------------------------------------------------
  ! B > 500: pure asymptotic wing
  !---------------------------------------------------------------------
  IF (B .GT. 500.0D0) THEN
    STARK_PROFILE = (1.5D0 / SB + 27.0D0 / B2) / B2 * CORR
    RETURN
  END IF

  ! Debye parameter interpolation weights
  IM = min(int(5.0D0 * P) + 1, 4)
  IP = IM + 1
  WTPP = 5.0D0 * (P - PP(IM))
  WTPM = 1.0D0 - WTPP

  !---------------------------------------------------------------------
  ! 25.12 < B <= 500: asymptotic with C/D correction
  !---------------------------------------------------------------------
  IF (B .GT. 25.12D0) THEN
    CC = C(IP, INDX) * WTPP + C(IM, INDX) * WTPM
    DD = D(IP, INDX) * WTPP + D(IM, INDX) * WTPM
    CORR = 1.0D0 + DD / (CC + B * SB)
    STARK_PROFILE = (1.5D0 / SB + 27.0D0 / B2) / B2 * CORR
    RETURN
  END IF

  !---------------------------------------------------------------------
  ! B <= 25.12: bilinear interpolation in correction table
  !---------------------------------------------------------------------
  ! Find beta grid interval
  JM = 1
  DO J = 2, 15
    IF (B .LE. BETAGRID(J)) THEN
      JM = J - 1
      EXIT
    END IF
  END DO
  JP = JM + 1

  WTBP = (B - BETAGRID(JM)) / (BETAGRID(JP) - BETAGRID(JM))
  WTBM = 1.0D0 - WTBP

  CBP = PROPBM(IP, JP, INDX) * WTPP + PROPBM(IM, JP, INDX) * WTPM
  CBM = PROPBM(IP, JM, INDX) * WTPP + PROPBM(IM, JM, INDX) * WTPM
  CORR = 1.0D0 + CBP * WTBP + CBM * WTBM

  ! Blended core/wing approximate profile
  PR1 = 0.0D0
  PR2 = 0.0D0
  WT = max(min(0.5D0 * (10.0D0 - B), 1.0D0), 0.0D0)
  IF (B .LE. 10.0D0) PR1 = 8.0D0 / (83.0D0 + (2.0D0 + 0.95D0 * B2) * B)
  IF (B .GE. 8.0D0)  PR2 = (1.5D0 / SB + 27.0D0 / B2) / B2

  STARK_PROFILE = (PR1 * WT + PR2 * (1.0D0 - WT)) * CORR
  RETURN

END FUNCTION STARK_PROFILE

!==========================================================================
! FUNCTION HPROF4
!
! Hydrogen line profile function including Stark, resonance, van der Waals,
! radiative, and Doppler broadening with fine structure.
!
! From Deane Peterson. Computes φ(Δλ) in Voigt function normalization.
!
! Three broadening regimes evaluated depending on whether the point
! is in the line core (IFCORE=1) or wing (IFCORE=0):
!
!   Core: only the dominant mechanism (Doppler, Lorentz, or Stark)
!   Wing: all three mechanisms summed
!
!   (1) Doppler core with resolved fine-structure components:
!       Exact patterns for alpha lines (Δn=1, n≤4),
!       m=∞ approximation for other lines (n≤4, m≤10),
!       single unresolved component otherwise.
!
!   (2) Lorentz wing: resonance + van der Waals + radiative damping.
!       Special treatment for Lyman alpha (n=1, m=2):
!       - Allard & Koester (1992) quasi-molecular H₂ satellite
!       - Allard (1997) H₂⁺ satellite cutoff profiles
!       - 4× enhanced resonance broadening near line center
!
!   (3) Stark wing: quasistatic ion microfield profile via STARK_PROFILE
!       with Griem electron impact corrections via EXPI.
!       Lyman alpha: half proton/half electron quasistatic,
!       with H₂⁺ quasi-molecular satellite.
!
! Temperature-dependent quantities cached when ITEMP changes.
! Line-dependent quantities cached when N,M change.
!
! Radiative damping: ASUM tables (switched from analytic, Aug 2009).
! Van der Waals: includes H₂ contribution (John Lester, Aug 2009).
!==========================================================================

REAL(8) FUNCTION HPROF4(N, M, J, DELW)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, M, J
  REAL(8),  INTENT(IN) :: DELW

  ! --- Saved depth-dependent arrays (recomputed when ITEMP changes) ---
  REAL(8),  SAVE :: PP(kw), FO(kw), GCON1(kw), GCON2(kw)
  REAL(8),  SAVE :: Y1B(kw), Y1S(kw), C1D(kw), C2D(kw)
  REAL(8),  SAVE :: T3NHE(kw), T3NH2(kw), EXP4492T(kw)
  REAL(8),  SAVE :: XNE4(kw), XNF4(kw), XNFP4(kw)
  INTEGER, SAVE :: ITEMP1 = 0

  ! --- Saved line-dependent quantities (recomputed when N,M change) ---
  INTEGER, SAVE :: N1 = 0, M1 = 0
  INTEGER, SAVE :: MMN, IFINS
  REAL(8),  SAVE :: XN, XN2, XM, XM2, XMN2, XM2MN2, GNM
  REAL(8),  SAVE :: XKNM, FREQNM, DBETA, WAVENM
  REAL(8),  SAVE :: C1CON, C2CON, RADAMP, RESONT, VDW, STARKC
  REAL(8),  SAVE :: Y1NUM, Y1WHT
  REAL(8),  SAVE :: FINEST(14), FINSWT(14)

  ! --- Stark pattern constants ---
  REAL(8), PARAMETER :: XKNMTB(4,3) = reshape((/ &
    0.0001716D0, 0.009019D0, 0.1001D0, 0.5820D0, &
    0.0005235D0, 0.01772D0,  0.171D0,  0.866D0, &
    0.0008912D0, 0.02507D0,  0.223D0,  1.02D0 /), (/4,3/))
  REAL(8), PARAMETER :: Y1WTM(2,2) = reshape((/ &
    1.D18, 1.D17, 1.D16, 1.D14 /), (/2,2/))
  ! RYDH (hydrogen Rydberg frequency) now FREQ_RYDH from mod_constants

  ! --- Fine structure for alpha lines (Δn=1) in freq × 10⁻⁷ ---
  REAL(8), PARAMETER :: STALPH(34) = (/ &
    -730.D0, 370.D0, 188.D0, 515.D0, 327.D0, 619.D0, &
    -772.D0, -473.D0, -369.D0, 120.D0, 256.D0, 162.D0, &
    285.D0, -161.D0, -38.3D0, 6.82D0, -174.D0, -147.D0, &
    -101.D0, -77.5D0, 55.D0, 126.D0, 75.D0, 139.D0, &
    -60.D0, 3.7D0, 27.D0, -69.D0, -42.D0, -18.D0, &
    -5.5D0, -9.1D0, -33.D0, -24.D0 /)
  REAL(8), PARAMETER :: STWTAL(34) = (/ &
    1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 2.D0, &
    3.D0, 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 4.D0, 6.D0, &
    1.D0, 2.D0, 3.D0, 4.D0, 1.D0, 2.D0, 1.D0, 2.D0, &
    1.D0, 4.D0, 6.D0, 1.D0, 7.D0, 6.D0, 4.D0, 4.D0, &
    4.D0, 5.D0 /)
  INTEGER, PARAMETER :: ISTAL(4) = (/ 1, 3, 10, 21 /)
  INTEGER, PARAMETER :: LNGHAL(4) = (/ 2, 7, 11, 14 /)

  ! --- Fine structure for m=∞ in freq × 10⁻⁷ ---
  REAL(8), PARAMETER :: STCOMP(5,4) = reshape((/ &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    468.D0, 576.D0, -522.D0, 0.D0, 0.D0, &
    260.D0, 290.D0, -33.D0, -140.D0, 0.0D0, &
    140.D0, 150.D0, 18.D0, -27.D0, -51.D0 /), (/5,4/))
  REAL(8), PARAMETER :: STCPWT(5,4) = reshape((/ &
    1.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    1.D0, 1.D0, 2.D0, 0.D0, 0.D0, &
    1.D0, 1.D0, 4.D0, 3.D0, 0.D0, &
    1.D0, 1.D0, 4.D0, 6.D0, 4.D0 /), (/5,4/))
  INTEGER, PARAMETER :: LNCOMP(4) = (/ 1, 3, 4, 5 /)

  ! --- Lyman alpha quasi-molecular satellite cutoffs (Allard 1997) ---
  ! H₂⁺ cutoff: Δν = -15000+100*(i-1) cm⁻¹, i=1..111 (to -4000)
  REAL(8), PARAMETER :: CUTOFFH2PLUS(111) = (/ &
    -15.14, -15.06, -14.97, -14.88, -14.80, -14.71, -14.62, -14.53, &
    -14.44, -14.36, -14.27, -14.18, -14.09, -14.01, -13.92, -13.83, &
    -13.74, -13.65, -13.57, -13.48, -13.39, -13.30, -13.21, -13.13, &
    -13.04, -12.95, -12.86, -12.77, -12.69, -12.60, -12.51, -12.40, &
    -12.29, -12.15, -12.02, -11.90, -11.76, -11.63, -11.53, -11.41, &
    -11.30, -11.22, -11.15, -11.09, -11.07, -11.06, -11.07, -11.09, &
    -11.12, -11.16, -11.19, -11.21, -11.24, -11.27, -11.30, -11.33, &
    -11.36, -11.39, -11.42, -11.45, -11.48, -11.48, -11.48, -11.48, &
    -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, &
    -11.48, -11.48, -11.48, -11.48, -11.41, -11.40, -11.39, -11.38, &
    -11.37, -11.36, -11.35, -11.34, -11.33, -11.32, -11.30, -11.29, &
    -11.28, -11.27, -11.27, -11.27, -11.26, -11.25, -11.24, -11.23, &
    -11.22, -11.21, -11.20, -11.19, -11.18, -11.17, -11.15, -11.14, &
    -11.13, -11.12, -11.11, -11.10, -11.09, -11.08, -11.07 /)

  ! H₂ cutoff: Δν = -22000+200*(i-1) cm⁻¹, i=1..91 (to -4000)
  REAL(8), PARAMETER :: CUTOFFH2(91) = (/ &
    -13.43, -13.32, -13.21, -13.10, -12.98, -12.86, -12.79, -12.72, &
    -12.65, -12.58, -12.51, -12.47, -12.45, -12.45, -12.48, -12.51, &
    -12.53, -12.56, -12.59, -12.62, -12.65, -12.69, -12.73, -12.77, &
    -12.81, -12.85, -12.87, -12.89, -12.90, -12.90, -12.90, -12.90, &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.89, -12.88, &
    -12.87, -12.86, -12.85, -12.84, -12.83, -12.81, -12.80, -12.79, &
    -12.78, -12.76, -12.74, -12.72, -12.70, -12.68, -12.65, -12.62, &
    -12.59, -12.56, -12.53 /)

  ! Radiative damping: summed A-values (switched from analytic, Aug 2009)
  REAL(8), PARAMETER :: ASUMLYMAN(100) = (/ &
    0.000E+00, 6.265E+08, 1.897E+08, 8.126E+07, 4.203E+07, &
    2.450E+07, 1.236E+07, 8.249E+06, 5.782E+06, 4.208E+06, &
    3.158E+06, 2.430E+06, 1.910E+06, 1.567E+06, 1.274E+06, &
    1.050E+06, 8.752E+05, 7.373E+05, 6.269E+05, 5.375E+05, &
    4.643E+05, 4.038E+05, 3.534E+05, 3.111E+05, 2.752E+05, &
    2.447E+05, 2.185E+05, 1.959E+05, 1.763E+05, 1.593E+05, &
    1.443E+05, 1.312E+05, 1.197E+05, 1.094E+05, 1.003E+05, &
    9.216E+04, 8.489E+04, 7.836E+04, 7.249E+04, 6.719E+04, &
    6.239E+04, 5.804E+04, 5.408E+04, 5.048E+04, 4.719E+04, &
    4.418E+04, 4.142E+04, 3.888E+04, 3.655E+04, 3.440E+04, &
    3.242E+04, 3.058E+04, 2.888E+04, 2.731E+04, 2.585E+04, &
    2.449E+04, 2.322E+04, 2.204E+04, 2.094E+04, 1.991E+04, &
    1.894E+04, 1.804E+04, 1.720E+04, 1.640E+04, 1.566E+04, &
    1.496E+04, 1.430E+04, 1.368E+04, 1.309E+04, 1.254E+04, &
    1.201E+04, 1.152E+04, 1.105E+04, 1.061E+04, 1.019E+04, &
    9.796E+03, 9.419E+03, 9.061E+03, 8.721E+03, 8.398E+03, &
    8.091E+03, 7.799E+03, 7.520E+03, 7.255E+03, 7.002E+03, &
    6.760E+03, 6.530E+03, 6.310E+03, 6.100E+03, 5.898E+03, &
    5.706E+03, 5.522E+03, 5.346E+03, 5.177E+03, 5.015E+03, &
    4.860E+03, 4.711E+03, 4.569E+03, 4.432E+03, 4.300E+03 /)

  REAL(8), PARAMETER :: ASUM(100) = (/ &
    0.000E+00, 4.696E+08, 9.980E+07, 3.017E+07, 1.155E+07, &
    5.189E+06, 2.616E+06, 1.437E+06, 8.444E+05, 5.234E+05, &
    3.389E+05, 2.275E+05, 1.575E+05, 1.120E+05, 8.142E+04, &
    6.040E+04, 4.560E+04, 3.496E+04, 2.719E+04, 2.141E+04, &
    1.711E+04, 1.377E+04, 1.119E+04, 9.166E+03, 7.572E+03, &
    6.341E+03, 5.338E+03, 4.523E+03, 3.854E+03, 3.302E+03, &
    2.844E+03, 2.460E+03, 2.138E+03, 1.866E+03, 1.635E+03, &
    1.438E+03, 1.269E+03, 1.124E+03, 9.983E+02, 8.894E+02, &
    7.947E+02, 7.120E+02, 6.396E+02, 5.759E+02, 5.198E+02, &
    4.703E+02, 4.263E+02, 3.873E+02, 3.526E+02, 3.215E+02, &
    2.938E+02, 2.689E+02, 2.465E+02, 2.264E+02, 2.082E+02, &
    1.918E+02, 1.769E+02, 1.634E+02, 1.512E+02, 1.400E+02, &
    1.298E+02, 1.206E+02, 1.121E+02, 1.043E+02, 9.720E+01, &
    9.066E+01, 8.465E+01, 7.912E+01, 7.403E+01, 6.933E+01, &
    6.498E+01, 6.097E+01, 5.725E+01, 5.381E+01, 5.061E+01, &
    4.765E+01, 4.489E+01, 4.232E+01, 3.994E+01, 3.771E+01, &
    3.563E+01, 3.369E+01, 3.188E+01, 3.019E+01, 2.860E+01, &
    2.712E+01, 2.572E+01, 2.442E+01, 2.319E+01, 2.204E+01, &
    2.096E+01, 1.994E+01, 1.898E+01, 1.808E+01, 1.722E+01, &
    1.642E+01, 1.566E+01, 1.495E+01, 1.427E+01, 1.363E+01 /)

  ! Local variables
  REAL(8)  :: DELstark, WL, FREQ4, DEL, DOP, HFWID
  REAL(8)  :: HWSTK, HWVDW, HWRAD, HWRES, HWLOR, HHW
  REAL(8)  :: HPROFLOR, HPROFRES, HPROFRAD, HPROFVDW
  REAL(8)  :: D, WTY1, Y1SCAL, C1, C2, G1, GNOT, BETA, Y1, Y2, GAM
  REAL(8)  :: PRQS, F, P1, FNS
  REAL(8)  :: CUTOFF, SPACING, CUTFREQ, FREQ15000, FREQ22000
  REAL(8)  :: BETA4000, PRQSP4000, CUTOFF4000
  REAL(8)  :: XNE16, TT4, T4, T43
  INTEGER :: I, K, NWID, IFCORE, IPOS, ICUT

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING HPROF4'

  ! ==================================================================
  ! Section 1: Temperature-dependent depth vectors
  ! ==================================================================
  IF (ITEMP .NE. ITEMP1) THEN
    ITEMP1 = ITEMP
    DO K = 1, NRHOX
      XNE4(K) = XNE(K)
      XNE16 = XNE(K)**(1.0D0/6.0D0)
      TT4 = T(K)
      PP(K) = XNE16 * 0.08989D0 / sqrt(TT4)
      FO(K) = XNE16**4 * 1.25D-9
      Y1B(K) = 2.0D0 / (1.0D0 + 0.012D0 / T(K) * sqrt(XNE(K) / T(K)))
      T4 = T(K) / 10000.0D0
      T43 = T4**0.3D0
      Y1S(K) = T43 / XNE16
      T3NHE(K) = T43 * XNFP(K, 3)
      ! 05aug2007
      T3NH2(K) = T43 * XNH2(K)
      C1D(K) = FO(K) * 78940.0D0 / TT4
      C2D(K) = FO(K)**2 / 5.96D-23 / XNE4(K)
      GCON1(K) = 0.2D0 + 0.09D0 * sqrt(T4) / (1.0D0 + XNE4(K) / 1.D13)
      GCON2(K) = 0.2D0 / (1.0D0 + XNE(K) / 1.D15)
      EXP4492T(K) = exp(-4492.0D0 / T(K))
    END DO
  END IF

  ! ==================================================================
  ! Section 2: Line-dependent constants (cached per N,M)
  ! ==================================================================
  IF (N .NE. N1 .OR. M .NE. M1) THEN
    N1 = N
    M1 = M
    MMN = M - N
    XN = dble(N)
    XN2 = XN * XN
    XM = dble(M)
    XM2 = XM * XM
    XMN2 = XM2 * XN2
    XM2MN2 = XM2 - XN2
    GNM = XM2MN2 / XMN2

    ! Stark pattern constant
    IF (MMN .LE. 3 .AND. N .LE. 4) THEN
      XKNM = XKNMTB(N, MMN)
    ELSE
      XKNM = 5.5D-5 / GNM * XMN2 / (1.0D0 + 0.13D0 / dble(MMN))
    END IF

    Y1NUM = 320.0D0
    IF (M .EQ. 2) Y1NUM = 550.0D0
    IF (M .EQ. 3) Y1NUM = 380.0D0

    Y1WHT = 1.D13
    IF (MMN .LE. 3) Y1WHT = 1.D14
    IF (MMN .LE. 2 .AND. N .LE. 2) Y1WHT = Y1WTM(N, MMN)

    FREQNM = FREQ_RYDH * GNM
    DBETA = CLIGHT_ANGS / FREQNM**2 / XKNM
    WAVENM = CLIGHT_ANGS / FREQNM
    C1CON = XKNM / WAVENM * GNM * XM2MN2
    C2CON = (XKNM / WAVENM)**2

    ! Radiative damping from ASUM tables (02aug2009)
!     RADAMP=1.389E9/XM**4.53/(1.+5./XM2/XM)
!     IF(N.NE.1)RADAMP=RADAMP+1.389E9/XN**4.53/(1.+5./XN2/XN)
    RADAMP = dble(ASUM(N)) + dble(ASUM(M))
    IF (N .EQ. 1) RADAMP = dble(ASUMLYMAN(M))
    RADAMP = RADAMP / FOURPI
    RADAMP = RADAMP / FREQNM

    ! Resonance broadening
    RESONT = HFNM(1, M) / XM / (1.0D0 - 1.0D0 / XM2)
    IF (N .NE. 1) RESONT = RESONT + HFNM(1, N) / XN / (1.0D0 - 1.0D0 / XN2)
!      FUDGE TO BASCHEK*2
!      RESONT=HFNM(1,M)/XM/(1.-1./XM2)*XM/3.*.791*2.
!      IF(N.NE.1)RESONT=RESONT+HFNM(1,N)/XN/(1.-1./XN2)*XN/3.*.791*2.
!     2 IS FOR CONVERTING XNFPH TO XNFH
!     RESONT=RESONT*5.593E-24/GNM*2.
!     RESONT=RESONT*5.593E-24/GNM
!     error in constant corrected 26nov95
    RESONT = RESONT * 3.92D-24 / GNM

    ! Van der Waals broadening
    VDW = 4.45D-26 / GNM * (XM2 * (7.0D0 * XM2 + 5.0D0))**0.4D0

    ! Stark broadening constant
    STARKC = 1.6678D-18 * FREQNM * XKNM

    ! Fine structure components
    IF (N .GT. 4 .OR. M .GT. 10) THEN
      ! Single unresolved component
      IFINS = 1
      FINEST(1) = 0.0D0
      FINSWT(1) = 1.0D0
    ELSE IF (MMN .NE. 1) THEN
      ! Non-alpha: use m=∞ structure
      IFINS = LNCOMP(N)
      DO I = 1, IFINS
        FINEST(I) = STCOMP(I, N) * 1.D7
        FINSWT(I) = STCPWT(I, N) / XN2
      END DO
    ELSE
      ! Alpha lines: exact pattern
      IFINS = LNGHAL(N)
      IPOS = ISTAL(N)
      DO I = 1, IFINS
        K = IPOS - 1 + I
        FINEST(I) = STALPH(K) * 1.D7
        FINSWT(I) = STWTAL(K) / XN2 / 3.0D0
      END DO
    END IF
  END IF  ! N,M changed

  ! ==================================================================
  ! Section 3: Profile evaluation at this depth and wavelength
  ! ==================================================================
  DELstark = -10.0D0 * DELW / WAVENM * FREQNM
  WL = WAVENM + DELW * 10.0D0
  FREQ4 = CLIGHT_ANGS / WL
  DEL = abs(FREQ4 - FREQNM)
  WL = WL / 10.0D0   ! convert to nm

  ! Half-widths (as Δν/ν)
  HWSTK = STARKC * FO(J)
  ! 05aug2009 John Lester: include H₂ van der Waals
  HWVDW = VDW * T3NHE(J) + 2.0D0 * VDW * T3NH2(J)
  HWRAD = RADAMP
  XNF4(J) = XNF(J, 1)
  HWRES = RESONT * XNF4(J)
  HWLOR = HWRES + HWVDW + HWRAD

  ! Determine dominant broadening mechanism
  ! NWID: 1=Doppler, 2=Lorentz, 3=Stark
  IF (DOPPLE(J, 1) .GE. HWSTK .AND. DOPPLE(J, 1) .GE. HWLOR) THEN
    NWID = 1
  ELSE IF (HWLOR .GE. HWSTK) THEN
    NWID = 2
  ELSE
    NWID = 3
  END IF

  HFWID = FREQNM * max(DOPPLE(J, 1), HWLOR, HWSTK)
  HPROF4 = 0.0D0
  IFCORE = 0
  IF (abs(DEL) .LE. HFWID) IFCORE = 1
  DOP = FREQNM * DOPPLE(J, 1)

  ! ------------------------------------------------------------------
  ! Doppler section: fine-structure resolved Gaussian core
  ! In wing (IFCORE=0): always computed. In core: only if NWID=1.
  ! ------------------------------------------------------------------
  IF (IFCORE .EQ. 0 .OR. NWID .EQ. 1) THEN
    DO I = 1, IFINS
      D = abs(FREQ4 - FREQNM - FINEST(I)) / DOP
!     IF(D.LE.7.)HPROF4=HPROF4+EXP(-D*D)/1.77245/DOP*FINSWT(I)
!     SAME NORMALIZATION AS VOIGT FUNCTION
      IF (D .LE. 7.0D0) HPROF4 = HPROF4 + exp(-D * D) * FINSWT(I)
    END DO
    IF (IFCORE .EQ. 1) RETURN
  END IF

  ! ------------------------------------------------------------------
  ! Lorentz section: resonance + van der Waals + radiative damping
  ! In wing (IFCORE=0): always computed. In core: only if NWID=2.
  ! ------------------------------------------------------------------
  IF (IFCORE .EQ. 0 .OR. NWID .EQ. 2) THEN
    IF (N .EQ. 1 .AND. M .EQ. 2) THEN
      ! === Lyman alpha special treatment ===
      ! Modify resonance broadening to match at 4000 cm⁻¹
      HWRES = HWRES * 4.0D0
      HWLOR = HWRES + HWVDW + HWRAD
      HHW = FREQNM * HWLOR

      IF (FREQ4 .GT. (82259.105D0 - 4000.0D0) * CLIGHT) THEN
        ! Near center: enhanced resonance Lorentzian
        ! error found by John Lester 31jul2009
        HPROFRES = HWRES * FREQNM / PI / (DEL**2 + HHW**2) &
                 * SQRTPI * DOP
      ELSE
        ! Far red wing: Allard & Koester (1992) H₂ satellite
        CUTOFF = 0.0D0
        IF (FREQ4 .GE. 50000.0D0 * CLIGHT) THEN
          ! Tabulated at 200 cm⁻¹ spacing
          SPACING = 200.0D0 * CLIGHT
          FREQ22000 = (82259.105D0 - 22000.0D0) * CLIGHT
          IF (FREQ4 .LT. FREQ22000) THEN
            CUTOFF = dble(CUTOFFH2(2) - CUTOFFH2(1)) / SPACING &
                   * (FREQ4 - FREQ22000) + dble(CUTOFFH2(1))
          ELSE
            ICUT = int((FREQ4 - FREQ22000) / SPACING)
            CUTFREQ = dble(ICUT) * SPACING + FREQ22000
            CUTOFF = dble(CUTOFFH2(ICUT+2) - CUTOFFH2(ICUT+1)) / SPACING &
                   * (FREQ4 - CUTFREQ) + dble(CUTOFFH2(ICUT+1))
          END IF
          XNFP4(J) = XNFP(J, 1)
          CUTOFF = (10.0D0**(CUTOFF - 14.0D0)) * XNFP4(J) * 2.0D0 &
                 / CLIGHT
        END IF
        HPROFRES = CUTOFF * SQRTPI * DOP
      END IF

      ! Radiative damping (cut off below Lyman edge to avoid double-counting
      ! with Rayleigh scattering in HRAYOP)
      HPROFRAD = HWRAD * FREQNM / PI / (DEL**2 + HHW**2) &
               * SQRTPI * DOP
!     CORRECTION TO LORENTZ PROFILE   ALLER P.164   NOT USED
!     HPROFRAD=HPROFRAD*4.*FREQ**2/(FREQ**2+FREQNM**2)
      IF (FREQ4 .LE. 2.463D15) HPROFRAD = 0.0D0

      ! Van der Waals (cut off below 60000 cm⁻¹ from line center)
      HPROFVDW = HWVDW * FREQNM / PI / (DEL**2 + HHW**2) &
               * SQRTPI * DOP
      IF (FREQ4 .LT. 1.8D15) HPROFVDW = 0.0D0

      HPROFLOR = HPROFRES + HPROFRAD + HPROFVDW
      HPROF4 = HPROF4 + HPROFLOR

    ELSE
      ! === Non-Lyman-alpha: simple combined Lorentzian ===
      HHW = FREQNM * HWLOR
      HPROFLOR = HHW / PI / (DEL**2 + HHW**2) * SQRTPI * DOP
      HPROF4 = HPROF4 + HPROFLOR
    END IF

    IF (IFCORE .EQ. 1) RETURN
  END IF

  ! ------------------------------------------------------------------
  ! Stark section: quasistatic ion microfield + electron impact
  ! Always reached in wing; in core only if NWID=3.
  ! ------------------------------------------------------------------
  WTY1 = 1.0D0 / (1.0D0 + XNE4(J) / Y1WHT)
  Y1SCAL = Y1NUM * Y1S(J) * WTY1 + Y1B(J) * (1.0D0 - WTY1)
  C1 = C1D(J) * C1CON * Y1SCAL
  C2 = C2D(J) * C2CON
  G1 = 6.77D0 * sqrt(C1)
  GNOT = G1 * max(0.0D0, 0.2114D0 + log(sqrt(C2) / C1)) &
       * (1.0D0 - GCON1(J) - GCON2(J))
  BETA = abs(DELstark) / FO(J) * DBETA
  Y1 = C1 * BETA
  Y2 = C2 * BETA**2
  GAM = GNOT

  IF (.NOT. (Y2 .LE. 1.D-4 .AND. Y1 .LE. 1.D-5)) THEN
!     GAM=G1*(.5*EXP(-MIN(80.,Y1))+VCSE1F(Y1)-.5*VCSE1F(Y2))*
    GAM = G1 * (0.5D0 * exp(-min(80.D0, Y1)) + EXPI(1, Y1) &
        - 0.5D0 * EXPI(1, Y2)) &
        * (1.0D0 - GCON1(J) / (1.0D0 + (90.0D0 * Y1)**3) &
        - GCON2(J) / (1.0D0 + 2000.0D0 * Y1))
    IF (GAM .LE. 1.D-20) GAM = 0.0D0
  END IF

  PRQS = STARK_PROFILE(BETA, PP(J), N, M)

  IF (M .LE. 2) THEN
    ! Lyman alpha Stark: assume quasistatic is half protons, half electrons
    PRQS = PRQS * 0.5D0
    CUTOFF = 0.0D0

    ! H₂⁺ quasi-molecular satellite (Allard 1997)
    IF (FREQ4 .GE. (82259.105D0 - 20000.0D0) * CLIGHT) THEN
      IF (FREQ4 .GT. (82259.105D0 - 4000.0D0) * CLIGHT) THEN
        ! Near center: use ratio method
        BETA4000 = 4000.0D0 * CLIGHT / FO(J) * DBETA
        PRQSP4000 = STARK_PROFILE(BETA4000, PP(J), N, M) * 0.5D0 / FO(J) * DBETA
        CUTOFF4000 = (10.0D0**(-11.07D0 - 14.0D0)) / CLIGHT * XNF(J, 2)
        HPROF4 = HPROF4 + CUTOFF4000 / PRQSP4000 * PRQS / FO(J) * DBETA &
               * SQRTPI * DOP
      ELSE
        ! Interpolate H₂⁺ cutoff table (100 cm⁻¹ spacing)
        FREQ15000 = (82259.105D0 - 15000.0D0) * CLIGHT
        SPACING = 100.0D0 * CLIGHT
        IF (FREQ4 .LT. FREQ15000) THEN
          CUTOFF = dble(CUTOFFH2PLUS(2) - CUTOFFH2PLUS(1)) / SPACING &
                 * (FREQ4 - FREQ15000) + dble(CUTOFFH2PLUS(1))
        ELSE
          ICUT = int((FREQ4 - FREQ15000) / SPACING)
          CUTFREQ = dble(ICUT) * SPACING + FREQ15000
          CUTOFF = dble(CUTOFFH2PLUS(ICUT+2) - CUTOFFH2PLUS(ICUT+1)) &
                 / SPACING * (FREQ4 - CUTFREQ) + dble(CUTOFFH2PLUS(ICUT+1))
        END IF
        CUTOFF = (10.0D0**(CUTOFF - 14.0D0)) / CLIGHT * XNF(J, 2)
        HPROF4 = HPROF4 + CUTOFF * SQRTPI * DOP
      END IF
    END IF
  END IF

  ! Final Stark assembly
  F = 0.0D0
  IF (GAM .GT. 0.0D0) F = GAM / PI / (GAM**2 + BETA**2)
  P1 = (0.9D0 * Y1)**2
  FNS = (P1 + 0.03D0 * sqrt(Y1)) / (P1 + 1.0D0)
  ! Same normalization as Voigt function
  HPROF4 = HPROF4 + (PRQS * (1.0D0 + FNS) + F) / FO(J) * DBETA &
         * SQRTPI * DOP
  RETURN

END FUNCTION HPROF4

!=======================================================================
! FUNCTION VOIGT
!
! Fast Voigt profile H(a,v) using pretabulated coefficients.
!
!   H(a,v) = (a/pi) * integral[-inf,+inf] exp(-t^2)/(a^2+(v-t)^2) dt
!
! where a = Lorentz damping parameter (gamma / 4*pi*DopplerWidth)
!       v = frequency offset in Doppler widths
!
! Algorithm uses four regions:
!   1. Small a (a < 0.2), v <= 10:
!      Taylor expansion in a: H ≈ H0 + H1*a + H2*a^2
!      using pretabulated H0TAB, H1TAB, H2TAB from TABVOIGT
!
!   2. Small a (a < 0.2), v > 10:
!      Lorentz wing: H ≈ a / (sqrt(pi) * v^2)
!
!   3. Moderate a (0.2 <= a <= 1.4, a+v <= 3.2):
!      4th-order Horner polynomial in a, built from H0/H1/H2 tables
!      with empirical correction factor
!
!   4. Large a or v (a > 1.4 or a+v > 3.2):
!      Asymptotic Lorentz expansion: H ≈ a/(sqrt(2*pi)*U) * correction
!      where U = sqrt(2) * (a^2 + v^2)
!
! Tables are filled by TABVOIGT at 200 steps per Doppler width.
!-----------------------------------------------------------------------

REAL(8) FUNCTION VOIGT(V, A)
  
  IMPLICIT NONE

  REAL(8), INTENT(IN) :: V, A

  INTEGER :: IV
  REAL(8)  :: VV, AA, U, AAU, VVU, UU
  REAL(8)  :: HH1, HH2, HH3, HH4

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING VOIGT'

  ! Table index: 200 steps per Doppler width, offset by 1
  IV = INT(V * 200.d0 + 1.5d0)

  IF (A .LT. 0.2d0) THEN
    !-----------------------------------------------------------------
    ! Small damping: Taylor expansion or Lorentz wing
    !-----------------------------------------------------------------
    IF (V .LE. 10.d0) THEN
      ! Quadratic Taylor expansion in a using pretabulated coefficients
      VOIGT = (H2TAB(IV)*A + H1TAB(IV))*A + H0TAB(IV)
    ELSE
      ! Far wing: Lorentz profile (a << v)
      VOIGT = 0.5642d0 * A / V**2
    ENDIF

  ELSE IF (A .LE. 1.4d0 .AND. A + V .LE. 3.2d0) THEN
    !-----------------------------------------------------------------
    ! Moderate damping: 4th-order Horner polynomial in a
    ! with empirical correction for normalization
    !-----------------------------------------------------------------
    VV = V * V
    HH1 = H1TAB(IV) + H0TAB(IV) * 1.12838d0
    HH2 = H2TAB(IV) + HH1 * 1.12838d0 - H0TAB(IV)
    HH3 = (1.d0 - H2TAB(IV)) * 0.37613d0 - HH1 * 0.66667d0*VV &
           + HH2 * 1.12838d0
    HH4 = (3.d0*HH3 - HH1) * 0.37613d0 + H0TAB(IV) * 0.66667d0*VV*VV

    VOIGT = ((((HH4*A + HH3)*A + HH2)*A + HH1)*A + H0TAB(IV)) &
          * (((-.122727278d0*A + .532770573d0)*A - .96284325d0)*A &
             + .979895032d0)

  ELSE
    !-----------------------------------------------------------------
    ! Large damping or far wing: asymptotic Lorentz with correction
    ! H ~ a / (sqrt(2*pi) * U) * [1 + correction terms / U^2]
    !-----------------------------------------------------------------
    AA = A * A
    VV = V * V
    U  = (AA + VV) * 1.4142d0   ! sqrt(2) * (a^2 + v^2)

    VOIGT = A * 0.79788d0 / U   ! a / sqrt(2*pi) / (a^2+v^2)

    IF (A .LE. 100.d0) THEN
      ! Higher-order correction terms
      AAU = AA / U
      VVU = VV / U
      UU  = U * U
      VOIGT = ((((AAU - 10.d0*VVU)*AAU*3.d0 + 15.d0*VVU*VVU) &
               + 3.d0*VV - AA) / UU + 1.d0) * VOIGT
    ENDIF
  ENDIF

END FUNCTION VOIGT

!=======================================================================
! EXPI: Exponential integral E_n(x)
!=======================================================================

FUNCTION EXPI(N, X)
!-----------------------------------------------------------------------
! Exponential integral E_n(x) for positive arguments.
!
!   E_n(x) = integral from 1 to infinity of exp(-x*t) / t^n dt
!
! Algorithm:
!   E_1(x) computed via rational approximations from Cody & Thacher,
!   Mathematics of Computation, 22, 641 (1968), in three ranges:
!     0 < x <= 1 : polynomial/polynomial - ln(x)
!     1 < x <= 4 : polynomial/polynomial * exp(-x)
!         x > 4  : asymptotic rational approximation * exp(-x)/x
!         x <= 0 : returns 0
!
!   E_n(x) for n > 1 via the recurrence:
!     E_n(x) = [exp(-x) - x * E_{n-1}(x)] / (n-1)
!
! Caching: E_1 result is cached for repeated calls with the same x
! but different n (common in radiative transfer calculations).
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N
  REAL(8),  INTENT(IN) :: X
  REAL(8) :: EXPI

  ! Cached values from previous call
  REAL(8), SAVE :: X_prev  = -1.D20
  REAL(8), SAVE :: EX_save = 0.d0   ! exp(-x)
  REAL(8), SAVE :: E1_save = 0.d0   ! E_1(x)

  ! Local variables
  REAL(8) :: EX, E1
  INTEGER :: I

  ! Cody & Thacher rational approximation coefficients
  ! Range 0 < x <= 1: E_1(x) = P(x)/Q(x) - ln(x)
  REAL(8), PARAMETER :: A0 = -44178.5471728217d0
  REAL(8), PARAMETER :: A1 =  57721.7247139444d0
  REAL(8), PARAMETER :: A2 =   9938.31388962037d0
  REAL(8), PARAMETER :: A3 =   1842.11088668000d0
  REAL(8), PARAMETER :: A4 =    101.093806161906d0
  REAL(8), PARAMETER :: A5 =      5.03416184097568d0
  REAL(8), PARAMETER :: B0 =  76537.3323337614d0
  REAL(8), PARAMETER :: B1 =  32597.1881290275d0
  REAL(8), PARAMETER :: B2 =   6106.10794245759d0
  REAL(8), PARAMETER :: B3 =    635.419418378382d0
  REAL(8), PARAMETER :: B4 =     37.2298352833327d0

  ! Range 1 < x <= 4: E_1(x) = exp(-x) * P(x)/Q(x)
  REAL(8), PARAMETER :: C0 = 4.65627107975096D-7
  REAL(8), PARAMETER :: C1 = 0.999979577051595d0
  REAL(8), PARAMETER :: C2 = 9.04161556946329d0
  REAL(8), PARAMETER :: C3 = 24.3784088791317d0
  REAL(8), PARAMETER :: C4 = 23.0192559391333d0
  REAL(8), PARAMETER :: C5 = 6.90522522784444d0
  REAL(8), PARAMETER :: C6 = 0.430967839469389d0
  REAL(8), PARAMETER :: D1 = 10.0411643829054d0
  REAL(8), PARAMETER :: D2 = 32.4264210695138d0
  REAL(8), PARAMETER :: D3 = 41.2807841891424d0
  REAL(8), PARAMETER :: D4 = 20.4494785013794d0
  REAL(8), PARAMETER :: D5 = 3.31909213593302d0
  REAL(8), PARAMETER :: D6 = 0.103400130404874d0

  ! Range x > 4: E_1(x) = exp(-x)/x * [1 + P(1/x)/Q(1/x)]
  REAL(8), PARAMETER :: E0 = -0.999999999998447d0
  REAL(8), PARAMETER :: E1C = -26.6271060431811d0   ! renamed from E1 to avoid conflict
  REAL(8), PARAMETER :: E2 = -241.055827097015d0
  REAL(8), PARAMETER :: E3 = -895.927957772937d0
  REAL(8), PARAMETER :: E4 = -1298.85688746484d0
  REAL(8), PARAMETER :: E5 = -545.374158883133d0
  REAL(8), PARAMETER :: E6 = -5.66575206533869d0
  REAL(8), PARAMETER :: F1 = 28.6271060422192d0
  REAL(8), PARAMETER :: F2 = 292.310039388533d0
  REAL(8), PARAMETER :: F3 = 1332.78537748257d0
  REAL(8), PARAMETER :: F4 = 2777.61949509163d0
  REAL(8), PARAMETER :: F5 = 2404.01713225909d0
  REAL(8), PARAMETER :: F6 = 631.657483280800d0

  !=====================================================================
  ! Compute E_1(x) (use cache if x unchanged)
  !=====================================================================
  IF (X .NE. X_prev) THEN
    EX = EXP(-X)
    X_prev  = X
    EX_save = EX

    IF (X .GT. 4.d0) THEN
      ! Asymptotic range: E_1 ~ exp(-x)/x * rational(1/x)
      E1 = (EX + EX*(E0 + (E1C + (E2 + (E3 + (E4 + (E5 + E6/X)/X)/X)/X)/X)/X) &
           / (X + F1 + (F2 + (F3 + (F4 + (F5 + F6/X)/X)/X)/X)/X)) / X

    ELSE IF (X .GT. 1.d0) THEN
      ! Mid-range rational approximation
      E1 = EX * (C6 + (C5 + (C4 + (C3 + (C2 + (C1 + C0*X)*X)*X)*X)*X)*X) &
              / (D6 + (D5 + (D4 + (D3 + (D2 + (D1 + X)*X)*X)*X)*X)*X)

    ELSE IF (X .GT. 0.d0) THEN
      ! Small-x: rational approximation minus logarithmic singularity
      E1 = (A0 + (A1 + (A2 + (A3 + (A4 + A5*X)*X)*X)*X)*X) &
         / (B0 + (B1 + (B2 + (B3 + (B4 + X)*X)*X)*X)*X) - LOG(X)

    ELSE
      ! x <= 0: return 0
      E1 = 0.d0
    ENDIF

    E1_save = E1
  ENDIF

  !=====================================================================
  ! Return E_1 or apply recurrence for E_n (n > 1)
  !=====================================================================
  EXPI = E1_save
  IF (N .EQ. 1) RETURN

  ! Recurrence: E_n(x) = [exp(-x) - x * E_{n-1}(x)] / (n-1)
  DO I = 1, N - 1
    EXPI = (EX_save - X * EXPI) / dble(I)
  END DO

END FUNCTION EXPI

!=========================================================================
! SUBROUTINE PFIRON(NELEM, ION, TLOG8, POTLOW8, PF)
!
!   Compute partition functions for iron-group elements (Ca through Ni,
!   Z = 20..28) by interpolation in a precomputed table.
!
!   The table PFTAB(7, 56, 10, 9) is indexed by:
!     dim 1: lowering potential bin (7 values: 500..32000 cm^-1)
!     dim 2: temperature bin (56 values covering log T = 3.32..5.25)
!     dim 3: ionization stage (1..10)
!     dim 4: element (1..9, mapped as NELEM - 19)
!
!   Temperature grid:
!     log T = 3.32..3.70 in steps of 0.02  (IT = 2..21)
!     log T = 3.70..4.00 in steps of 0.03  (IT = 21..31)
!     log T = 4.00..5.25 in steps of 0.05  (IT = 31..56)
!
!   Interpolation is bilinear in (temperature, potential lowering).
!
!   Arguments:
!     NELEM    — atomic number (20..28)
!     ION      — ionization stage (1..10; outside this range returns PF=0)
!     TLOG8    — log10(T) (REAL*8)
!     POTLOW8  — potential lowering in cm^-1 (REAL*8)
!     PF       — output: LINEAR partition function (not log10)
!
!   Data source: pfiron.dat (read on first call from DATADIR)
!
!   2026 modernization: an optional Barklem & Collet (2016) multiplicative
!   hybrid is applied at the end of this routine, rescaling the Kurucz
!   table to B&C NIST-based values over the photospheric T range.  See
!   the BC_RATIO precompute and runtime application below for details.
!=========================================================================

SUBROUTINE PFIRON(NELEM, ION, TLOG8, POTLOW8, PF)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NELEM, ION
  REAL(8),  INTENT(IN) :: TLOG8, POTLOW8
  REAL(8),  INTENT(OUT) :: PF

  ! --- Table dimensions ---
  INTEGER, PARAMETER :: NPOT = 7       ! potential lowering bins
  INTEGER, PARAMETER :: NTEMP = 56     ! temperature bins
  INTEGER, PARAMETER :: NION_TAB = 10  ! max ionization stages in table
  INTEGER, PARAMETER :: NELEM_TAB = 9  ! elements: Ca(20)..Ni(28)

  ! --- Persistent table data (read once from file) ---
  REAL(8),  SAVE :: PFTAB(NPOT, NTEMP, NION_TAB, NELEM_TAB)
  REAL(8),  SAVE :: POTLO(NPOT)
  REAL(8),  SAVE :: POTLOLOG(NPOT)
  LOGICAL, SAVE :: INITIALIZED = .FALSE.

  ! --- Precomputed B&C hybrid cache (for ION in 1..3 only) ---
  ! BC_RATIO(IT, ION, IELEM) is a LINEAR multiplicative rescaling applied
  ! at runtime as:   PF = PF_kurucz * BC_RATIO
  ! so a value of 1.0 is a no-op.
  !
  ! At init time, BC_RATIO is populated as:
  !   IT corresponding to T_IT <= T_LOW_BC :
  !       BC_RATIO = U_BC(T_IT) / PFTAB(1, IT, ION, IELEM)       (full hybrid)
  !   IT corresponding to T_IT >= T_HIGH_BC :
  !       BC_RATIO = 1                                            (pure Kurucz)
  !   IT in the blend window (T_LOW_BC < T_IT < T_HIGH_BC) :
  !       BC_RATIO = 1 + (1 - w) * (R - 1)  where R = U_BC/PFTAB,
  !                                        w = smoothstep((T-TLOW)/(THIGH-TLOW))
  !
  ! The smoothstep taper makes BC_RATIO a C^1-continuous function of T_IT
  ! at IT=31 (T=10000 K), matching the PFGROUND_HYBRID dispatcher policy.
  ! This prevents a spurious discontinuity in d(PF)/dT that would otherwise
  ! trip up ATLAS's adiabatic-gradient finite differences.
  REAL(8),  SAVE :: BC_RATIO(NTEMP, 3, NELEM_TAB) = 1.0d0

  ! --- Local variables ---
  REAL(8)  :: TLOG, POTLOW, F, P
  INTEGER :: IT, LOW, I, J, K, L, IELEM

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING PFIRON'

  ! --- Read table on first call ---
  IF (.NOT. INITIALIZED) THEN
    OPEN(UNIT=89, FILE=trim(DATADIR)//'pfiron.dat', STATUS='OLD', ACTION='READ')
    READ(89, '(A)') ; READ(89, '(A)') ; READ(89, '(A)')
    READ(89, *) ((((PFTAB(I,J,K,L), I=1,NPOT), J=1,NTEMP), K=1,NION_TAB), L=1,NELEM_TAB)
    CLOSE(89)
    POTLO = (/ 500.D0, 1000.D0, 2000.D0, 4000.D0, 8000.D0, 16000.D0, 32000.D0 /)
    POTLOLOG = (/ 2.69897D0, 3.D0, 3.30103D0, 3.60206D0, 3.90309D0, 4.20412D0, 4.50515D0 /)

    ! --- Precompute linear B&C-to-Kurucz ratio at each PFIRON T grid point.
    !     Populate IT_init = 2..31 (log T = 3.32..4.00, T ~ 2000..10000 K).
    !     Below T_LOW_BC (9000 K) apply the full ratio.  Above T_HIGH_BC
    !     (10000 K) apply no correction.  In between, smoothstep-taper
    !     the ratio so BC_RATIO is C^1-continuous at IT=31.  IT >= 32
    !     remains at the default 1.0 (pure Kurucz).
    !
    !     Safeguards:
    !       - U_BC must be positive and finite (< 1e6).
    !       - PFTAB(1, ...) must be a sensible partition function
    !         (>= 0.99 and < 1e6; PFTAB starts at 1.0 for most species).
    !       - Raw ratio must be in [0.5, 2.0] (factor of 2 from Kurucz;
    !         anything more extreme is almost certainly bad data).
    !     If any guard fails, leave BC_RATIO at its default 1.0 -- the
    !     hybrid becomes a pure no-op for that (IT, ION, IELEM) slot.
    BLOCK
      REAL(8) :: TLOG_IT, T_IT, U_BC_val, RATIO, x, w, taper
      INTEGER :: IT_init, ION_init
      REAL(8), PARAMETER :: T_LOW_BC  =  9000.0d0
      REAL(8), PARAMETER :: T_HIGH_BC = 10000.0d0
      DO IELEM = 1, NELEM_TAB
        DO ION_init = 1, 3
          DO IT_init = 2, 31
            IF (IT_init .LE. 21) THEN
              TLOG_IT = 3.32D0 + (IT_init - 2) * 0.02D0
            ELSE
              TLOG_IT = 3.70D0 + (IT_init - 21) * 0.03D0
            END IF
            T_IT = 10.0D0**TLOG_IT
            ! Skip the top endpoint: at T_IT = 10000 K U_BC returns
            ! the sentinel anyway, and we want BC_RATIO = 1.0 there.
            IF (T_IT .GE. T_HIGH_BC) CYCLE
            U_BC_val = U_BC(IELEM + 19, ION_init, T_IT)
            IF (U_BC_val .GT. 0.0D0 .AND. U_BC_val .LT. 1.0D6 .AND. &
                PFTAB(1, IT_init, ION_init, IELEM) .GE. 0.99D0 .AND. &
                PFTAB(1, IT_init, ION_init, IELEM) .LT. 1.0D6) THEN
              RATIO = U_BC_val / PFTAB(1, IT_init, ION_init, IELEM)
              IF (RATIO .GE. 0.5D0 .AND. RATIO .LE. 2.0D0) THEN
                ! Compute the smoothstep taper.  taper = 1 below T_LOW_BC,
                ! 0 above T_HIGH_BC, smooth in between.
                IF (T_IT .LE. T_LOW_BC) THEN
                  taper = 1.0D0
                ELSE
                  x = (T_IT - T_LOW_BC) / (T_HIGH_BC - T_LOW_BC)
                  w = x * x * (3.0D0 - 2.0D0 * x)   ! smoothstep
                  taper = 1.0D0 - w
                END IF
                ! Attenuate the ratio toward 1.0 by (1 - taper).
                BC_RATIO(IT_init, ION_init, IELEM) = &
                  1.0D0 + taper * (RATIO - 1.0D0)
              END IF
            END IF
          END DO
        END DO
      END DO
    END BLOCK

    INITIALIZED = .TRUE.
  END IF

  ! --- Bounds check: ION must be within table range ---
  !     For higher ionization stages (ION > NION_TAB), the table has no data.
  !     Return PF = 0 (matching pristine behavior; callers of PFSAHA never
  !     hit this path in practice because NION in 1..10 covers all stages).
  IF (ION .LT. 1 .OR. ION .GT. NION_TAB) THEN
    IF (IDEBUG .EQ. 1) WRITE(6,'(A,I4,A,I4,A)') &
      '  PFIRON: ION =', ION, ' out of table range (1..', NION_TAB, '), returning PF=0'
    PF = 0.0D0
    RETURN
  END IF

  ! --- Bounds check: NELEM must map to valid table index ---
  IELEM = NELEM - 19
  IF (IELEM .LT. 1 .OR. IELEM .GT. NELEM_TAB) THEN
    WRITE(6,'(A,I4,A)') ' PFIRON: NELEM =', NELEM, ' out of range (20..28), returning PF=0'
    PF = 0.0D0
    RETURN
  END IF

  TLOG = TLOG8
  POTLOW = POTLOW8

  ! --- Determine temperature bin IT and interpolation fraction F ---
  IF (TLOG .LE. 3.7D0) THEN
    ! Low-T grid: log T = 3.32..3.70, step = 0.02
    IT = int((TLOG - 3.32D0) / 0.02D0) + 2
    IT = max(IT, 2)
    F = (TLOG - (IT - 2) * 0.02D0 - 3.32D0) / 0.02D0
  ELSE IF (TLOG .LE. 4.0D0) THEN
    ! Mid-T grid: log T = 3.70..4.00, step = 0.03
    IT = int((TLOG - 3.7D0) / 0.03D0) + 21
    F = (TLOG - (IT - 21) * 0.03D0 - 3.7D0) / 0.03D0
  ELSE
    ! High-T grid: log T = 4.00..5.25, step = 0.05
    IT = int((TLOG - 4.0D0) / 0.05D0) + 31
    IT = min(IT, NTEMP)
    F = (TLOG - (IT - 31) * 0.05D0 - 4.0D0) / 0.05D0
  END IF

  ! Safety clamp: ensure IT-1 >= 1 and IT <= NTEMP
  IT = max(IT, 2)
  IT = min(IT, NTEMP)
  F = max(0.0D0, min(1.0D0, F))

  ! --- Determine potential lowering bin LOW ---
  LOW = 1
  IF (POTLOW .GE. POTLO(1)) THEN
    DO LOW = 2, NPOT
      IF (POTLOW .LT. POTLO(LOW)) EXIT
    END DO
    IF (LOW .GT. NPOT) LOW = NPOT
  END IF

  IF (LOW .GT. 1 .AND. POTLOW .GE. POTLO(1) .AND. POTLOW .LT. POTLO(LOW)) THEN
    ! Two-bin interpolation in potential lowering
    P = (log10(POTLOW) - POTLOLOG(LOW-1)) / 0.30103D0
    PF = P * (F * PFTAB(LOW, IT, ION, IELEM) &
            + (1.0D0 - F) * PFTAB(LOW, IT-1, ION, IELEM)) &
       + (1.0D0 - P) * (F * PFTAB(LOW-1, IT, ION, IELEM) &
            + (1.0D0 - F) * PFTAB(LOW-1, IT-1, ION, IELEM))
  ELSE
    ! Single-bin interpolation (POTLOW below first bin or above last)
    PF = F * PFTAB(LOW, IT, ION, IELEM) &
       + (1.0D0 - F) * PFTAB(LOW, IT-1, ION, IELEM)
  END IF

  ! --- B&C hybrid rescaling (LINEAR multiplicative correction).
  !     Kurucz PFIRON returns a linear partition function.  Multiply by
  !     the precomputed BC_RATIO array (T-interpolated) to anchor against
  !     modern B&C (2016) NIST-based values for ION in 1..3 and
  !     T <= 10000 K.  Above 10000 K (IT >= 32) BC_RATIO is 1.0, so the
  !     hybrid is a no-op and PFIRON returns the pristine Kurucz value.
  !     Across the 9000..10000 K blend window BC_RATIO is smoothstep-
  !     tapered toward 1.0 to keep PF continuous and C^1 at IT=31,
  !     which ATLAS's convergence and finite-difference routines require.
  IF (ION .LE. 3 .AND. IT .GE. 2) THEN
    PF = PF * (F * BC_RATIO(IT, ION, IELEM) &
            + (1.0D0 - F) * BC_RATIO(IT-1, ION, IELEM))
  END IF

END SUBROUTINE PFIRON

!=======================================================================
! FUNCTION PFGROUND_KURUCZ(NELION, T)
!
!   Kurucz's original ground-state partition function routine, preserved
!   verbatim from the pristine ATLAS12 source.  Returns statistical
!   weight of the ground configuration (with temperature-dependent
!   corrections for certain species where Kurucz hand-coded multi-level
!   Boltzmann sums) for each (Z, ION) combination, identified via
!     NELION = (Z-1)*6 + ION,   ION in 1..6 (1=neutral, 2=singly ionized).
!
!   Used as the high-T fallback by PFGROUND_HYBRID() -- above 10,000 K B&C
!   is undefined and this routine takes over.  Also used as the fallback
!   for species outside B&C's Z=1..92, ION=1..3 coverage, at any T.
!
!   Original data from Kurucz, with corrections by J. Laird (2004) and
!   K. Bischof.  Historical bug fix: missing EXP() in F V (NELION=54).
!
!   Iron-group elements (Z=20..28, NELION=115..168) return 1.0 because
!   their partition functions are handled via PFIRON.
!=======================================================================

FUNCTION PFGROUND_KURUCZ(NELION, T)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)  :: NELION
  REAL(8),  INTENT(IN)  :: T
  REAL(8)               :: PFGROUND_KURUCZ

  ! Local constants
  ! HCK now from mod_constants (updated to CODATA 2018)

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING PFGROUND_KURUCZ'

  ! Default: bare ground state
  PFGROUND_KURUCZ = 1.0d0

  !---------------------------------------------------------------------
  ! Iron group (Z=20-28, Ca-Ni): partition functions from PFIRON tables.
  ! NELION = 115..168 for these elements.
  !---------------------------------------------------------------------
  IF (NELION .GE. 115 .AND. NELION .LE. 168) RETURN

  SELECT CASE (NELION)

  ! --- Z= 1  H  ---
  CASE (  1); PFGROUND_KURUCZ = 2.                    ! H  I
  CASE (  2); PFGROUND_KURUCZ = 1.                    ! H  II
  CASE (  3); PFGROUND_KURUCZ = 1.                    ! H  III
  CASE (  4); PFGROUND_KURUCZ = 1.                    ! H  IV
  CASE (  5); PFGROUND_KURUCZ = 1.                    ! H  V
  CASE (  6); PFGROUND_KURUCZ = 1.                    ! H  VI

  ! --- Z= 2  He ---
  CASE (  7); PFGROUND_KURUCZ = 1.                    ! He I
  CASE (  8); PFGROUND_KURUCZ = 2.                    ! He II
  CASE (  9); PFGROUND_KURUCZ = 1.                    ! He III
  CASE ( 10); PFGROUND_KURUCZ = 1.                    ! He IV
  CASE ( 11); PFGROUND_KURUCZ = 1.                    ! He V
  CASE ( 12); PFGROUND_KURUCZ = 1.                    ! He VI

  ! --- Z= 3  Li ---
  CASE ( 13); PFGROUND_KURUCZ = 2.                    ! Li I
  CASE ( 14); PFGROUND_KURUCZ = 1.                    ! Li II
  CASE ( 15); PFGROUND_KURUCZ = 2.                    ! Li III
  CASE ( 16); PFGROUND_KURUCZ = 1.                    ! Li IV
  CASE ( 17); PFGROUND_KURUCZ = 1.                    ! Li V
  CASE ( 18); PFGROUND_KURUCZ = 1.                    ! Li VI

  ! --- Z= 4  Be ---
  CASE ( 19); PFGROUND_KURUCZ = 1.                    ! Be I
  CASE ( 20); PFGROUND_KURUCZ = 2.                    ! Be II
  CASE ( 21); PFGROUND_KURUCZ = 1.                    ! Be III
  CASE ( 22); PFGROUND_KURUCZ = 2.                    ! Be IV
  CASE ( 23); PFGROUND_KURUCZ = 1.                    ! Be V
  CASE ( 24); PFGROUND_KURUCZ = 1.                    ! Be VI

  ! --- Z= 5  B  ---
  CASE ( 25)  ! B  I
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*15.254)
  CASE ( 26)  ! B  II
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*37336.7)+3.*EXP(-HCK/T*37342.4)+ 5.*EXP(-HCK/T* 37358.3)+3.*EXP(-HCK/T*73396.60)
  CASE ( 27)  ! B  III
    PFGROUND_KURUCZ = 2.+2.*EXP(-HCK/T*48358.40)+4.*EXP(-HCK/T*48392.50)
  CASE ( 28); PFGROUND_KURUCZ = 1.                    ! B  IV
  CASE ( 29); PFGROUND_KURUCZ = 2.                    ! B  V
  CASE ( 30); PFGROUND_KURUCZ = 1.                    ! B  VI

  ! --- Z= 6  C  ---
  CASE ( 31)  ! C  I
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*16.40)+5.*EXP(-HCK/T*43.40)+ 5.*EXP(-HCK/T*10192.63)
  CASE ( 32)  ! C  II
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*63.42)
  CASE ( 33); PFGROUND_KURUCZ = 1.                    ! C  III
  CASE ( 34)  ! C  IV
    PFGROUND_KURUCZ = 2.+2.*EXP(-HCK/T*64484.0)+4.*EXP(-HCK/T*64591.7)
  CASE ( 35); PFGROUND_KURUCZ = 1.                    ! C  V
  CASE ( 36); PFGROUND_KURUCZ = 2.                    ! C  VI

  ! --- Z= 7  N  ---
  CASE ( 37); PFGROUND_KURUCZ = 4.                    ! N  I
  CASE ( 38)  ! N  II
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*48.7)+5.*EXP(-HCK/T*130.8)+ 5.*EXP(-HCK/T*15316.2)
  CASE ( 39)  ! N  III
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*174.4)
  CASE ( 40); PFGROUND_KURUCZ = 1.                    ! N  IV
  CASE ( 41)  ! N  V
    PFGROUND_KURUCZ = 2.+2.*EXP(-HCK/T*80463.2)+4.*EXP(-HCK/T*80721.9)
  CASE ( 42); PFGROUND_KURUCZ = 1.                    ! N  VI

  ! --- Z= 8  O  ---
  CASE ( 43)  ! O  I
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*158.265)+EXP(-HCK/T*226.977)
  CASE ( 44)  ! O  II
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*26810.55)+4.*EXP(-HCK/T*26830.57)
  CASE ( 45)  ! O  III
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*113.178)+5.*EXP(-HCK/T*306.174)+ 5.*EXP(-HCK/T*20273.27)+1.*EXP(-HCK/T*43185.74)+ &
      5.*EXP(-HCK/T*60324.79)
  CASE ( 46)  ! O  IV
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*385.9)+2.*EXP(-HCK/T*71439.8)+ 4.*EXP(-HCK/T*71570.1)+6.*EXP(-HCK/T*71755.5)
  CASE ( 47)  ! O  V
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*81942.5)+3.*EXP(-HCK/T*82078.6)+ 5.*EXP(-HCK/T*82385.3)
  CASE ( 48)  ! O  VI
    PFGROUND_KURUCZ = 2.+2.*EXP(-HCK/T*96375.0)+4.*EXP(-HCK/T*96907.5)

  ! --- Z= 9  F  ---
  CASE ( 49)  ! F  I
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*404.1)
  CASE ( 50)  ! F  II
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*341.0)+EXP(-HCK/T*489.9)+ 5.*EXP(-HCK/T*20873.4)+1.*EXP(-HCK/T*44918.1)
  CASE ( 51)  ! F  III
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*34087.4)+4.*EXP(-HCK/T*34123.2)+ 4.*EXP(-HCK/T*51561.4)+2.*EXP(-HCK/T*51560.6)
  CASE ( 52)  ! F  IV
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*225.2)+5.*EXP(-HCK/T*612.2)+ 5.*EXP(-HCK/T*25238.2)+1.*EXP(-HCK/T*53541.2)+ 5.*EXP(-HCK/T*74194.7)
  CASE ( 53)  ! F  V
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*744.5)+ 2.*EXP(-HCK/T*85790.2)+4.*EXP(-HCK/T*86043.5)+ 6.*EXP(-HCK/T*86407.0)
  CASE ( 54)  ! F  VI
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*96590.)+3.*EXP(-HCK/T*96850.)+ 5.*EXP(-HCK/T*97427.)

  ! --- Z=10  Ne ---
  CASE ( 55); PFGROUND_KURUCZ = 1.                    ! Ne I
  CASE ( 56)  ! Ne II
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*780.45)
  CASE ( 57)  ! Ne III
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*642.9)+EXP(-HCK/T*920.4)+ 4.*EXP(-HCK/T*96907.5)+5.*EXP(-HCK/T*25840.8)+ 1.*EXP(-HCK/T*55750.6)
  CASE ( 58)  ! Ne IV
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*41234.6)+4.*EXP(-HCK/T*41279.5)+ 2.*EXP(-HCK/T*62434.6)+4.*EXP(-HCK/T*62441.3)
  CASE ( 59)  ! Ne V
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*414.)+5.*EXP(-HCK/T*1112.)+ 5.*EXP(-HCK/T*30291.5)+1.*EXP(-HCK/T*63913.6)+ 5.*EXP(-HCK/T*88360.)
  CASE ( 60)  ! Ne VI
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*1310.)+2.*EXP(-HCK/T*100261.)+ 4.*EXP(-HCK/T*100704.)+6.*EXP(-HCK/T*101347.)

  ! --- Z=11  Na ---
  CASE ( 61); PFGROUND_KURUCZ = 2.                    ! Na I
  CASE ( 62); PFGROUND_KURUCZ = 1.                    ! Na II
  CASE ( 63)  ! Na III
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*780.45)
  CASE ( 64)  ! Na IV
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*642.9)+EXP(-HCK/T*920.4)+ 5.*EXP(-HCK/T*30839.8)+1.*EXP(-HCK/T*66496.)
  CASE ( 65)  ! Na V
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*48330.)+4.*EXP(-HCK/T*48366.)+ 2.*EXP(-HCK/T*73218.)+4.*EXP(-HCK/T*73255.)
  CASE ( 66)  ! Na VI
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*414.)+5.*EXP(-HCK/T*1112.)+ 5.*EXP(-HCK/T*35498.)+1.*EXP(-HCK/T*74414.)

  ! --- Z=12  Mg ---
  CASE ( 67); PFGROUND_KURUCZ = 1.                    ! Mg I
  CASE ( 68); PFGROUND_KURUCZ = 2.                    ! Mg II
  CASE ( 69); PFGROUND_KURUCZ = 1.                    ! Mg III
  CASE ( 70)  ! Mg IV
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*2238.)
  CASE ( 71)  ! Mg V
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*1782.1)+EXP(-HCK/T*2521.8)+ 5.*EXP(-HCK/T*35926.)+1.*EXP(-HCK/T*77279.)
  CASE ( 72)  ! Mg VI
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*55356.)+4.*EXP(-HCK/T*55372.8)+ 2.*EXP(-HCK/T*83920.0)+4.*EXP(-HCK/T*84028.4)

  ! --- Z=13  Al ---
  CASE ( 73)  ! Al I
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*112.061)
  CASE ( 74); PFGROUND_KURUCZ = 1.                    ! Al II
  CASE ( 75); PFGROUND_KURUCZ = 2.                    ! Al III
  CASE ( 76); PFGROUND_KURUCZ = 1.                    ! Al IV
  CASE ( 77)  ! Al V
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*3442.)
  CASE ( 78)  ! Al VI
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*2732.)+EXP(-HCK/T*3829.)+ 5.*EXP(-HCK/T*41167.)+1.*EXP(-HCK/T*88213.)

  ! --- Z=14  Si ---
  CASE ( 79)  ! Si I
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*77.115)+5.*EXP(-HCK/T*223.157)+ 5.*EXP(-HCK/T*6298.850)
  CASE ( 80)  ! Si II
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*287.32)+2.*EXP(-HCK/T*42824.35)+ 4.*EXP(-HCK/T*42932.68)+6.*EXP(-HCK/T*43107.97)
  CASE ( 81)  ! Si III
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*52724.69)+3.*EXP(-HCK/T*52853.28)+ 5.*EXP(-HCK/T*53115.01)+3.*EXP(-HCK/T*82884.41)
  CASE ( 82); PFGROUND_KURUCZ = 2.                    ! Si IV
  CASE ( 83); PFGROUND_KURUCZ = 1.                    ! Si V
  CASE ( 84)  ! Si VI
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*5090.)

  ! --- Z=15  P  ---
  CASE ( 85); PFGROUND_KURUCZ = 4.                    ! P  I
  CASE ( 86)  ! P  II
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*164.90)+5.*EXP(-HCK/T*469.12)+ 5.*EXP(-HCK/T*8882.31)+1.*EXP(-HCK/T*21575.63)
  CASE ( 87)  ! P  III
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*559.14)+2.*EXP(-HCK/T*56021.67)+ 4.*EXP(-HCK/T*57125.98)+6.*EXP(-HCK/T*57454.00)+ &
      4.*EXP(-HCK/T*74916.85)+6.*EXP(-HCK/T*74945.86)
  CASE ( 88)  ! P  IV
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*67918.03)+3.*EXP(-HCK/T*68146.48)+ 5.*EXP(-HCK/T*68615.17)
  CASE ( 89)  ! P  V
    PFGROUND_KURUCZ = 2.+2.*EXP(-HCK/T*88651.87)+4.*EXP(-HCK/T*89447.25)
  CASE ( 90); PFGROUND_KURUCZ = 1.                    ! P  VI

  ! --- Z=16  S  ---
  CASE ( 91)  ! S  I
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*396.055)+EXP(-HCK/T*573.640)+ 5.*EXP(-HCK/T*9238.609)
  CASE ( 92)  ! S  II
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*14852.94)+6.*EXP(-HCK/T*14884.73)+ 2.*EXP(-HCK/T*24524.83)+4.*EXP(-HCK/T*24571.54)
  CASE ( 93)  ! S  III
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*298.69)+5.*EXP(-HCK/T*833.08)+ 5.*EXP(-HCK/T*11322.7)+1.*EXP(-HCK/T*27161.0)
  CASE ( 94)  ! S  IV
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*951.43)+2.*EXP(-HCK/T*71184.1)+ 4.*EXP(-HCK/T*71528.7)+6.*EXP(-HCK/T*72074.4)+ &
      4.*EXP(-HCK/T*94103.1)+6.*EXP(-HCK/T*94150.4)
  CASE ( 95)  ! S  V
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*83024.0)+3.*EXP(-HCK/T*83393.5)+ 5.*EXP(-HCK/T*84155.2)
  CASE ( 96); PFGROUND_KURUCZ = 2.                    ! S  VI

  ! --- Z=17  Cl ---
  CASE ( 97)  ! Cl I
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*882.36)
  CASE ( 98)  ! Cl II
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*696.1)+EXP(-HCK/T*996.4)+ 5.*EXP(-HCK/T*11653.58)+1.*EXP(-HCK/T*27878.02)
  CASE ( 99)  ! Cl III
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*18053.)+6.*EXP(-HCK/T*18118.6)+ 2.*EXP(-HCK/T*29812.)+4.*EXP(-HCK/T*29907.)
  CASE (100)  ! Cl IV
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*491.)+5.*EXP(-HCK/T*1341.)+ 5.*EXP(-HCK/T*13767.6)+1.*EXP(-HCK/T*32547.8)+ 5.*EXP(-HCK/T*65000.)
  CASE (101)  ! Cl V
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*1490.8)+2.*EXP(-HCK/T*86000.)+ 4.*EXP(-HCK/T*86538.)+6.*EXP(-HCK/T*87381.)
  CASE (102)  ! Cl VI
    PFGROUND_KURUCZ = 1.+1.*EXP(-HCK/T*97405.)+3.*EXP(-HCK/T*97958.)+ 5.*EXP(-HCK/T*99123.)

  ! --- Z=18  Ar ---
  CASE (103); PFGROUND_KURUCZ = 1.                    ! Ar I
  CASE (104)  ! Ar II
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*1431.41)
  CASE (105)  ! Ar III
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*1112.1)+EXP(-HCK/T*1570.2)+ 5.*EXP(-HCK/T*14010.004)+1.*EXP(-HCK/T*33265.724)
  CASE (106)  ! Ar IV
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*21090.4)+6.*EXP(-HCK/T*21219.3)+ 2.*EXP(-HCK/T*34855.5)+4.*EXP(-HCK/T*35032.6)
  CASE (107)  ! Ar V
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*765.)+5.*EXP(-HCK/T*2030.)+ 5.*EXP(-HCK/T*16298.9)+1.*EXP(-HCK/T*37912.0)+ 5.*EXP(-HCK/T*84100.0)
  CASE (108)  ! Ar VI
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*2208.)

  ! --- Z=19  K  ---
  CASE (109); PFGROUND_KURUCZ = 2.                    ! K  I
  CASE (110); PFGROUND_KURUCZ = 1.                    ! K  II
  CASE (111)  ! K  III
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*2166.)
  CASE (112)  ! K  IV
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*1673.)+EXP(-HCK/T*2325.)+ 5.*EXP(-HCK/T*16384.1)+EXP(-HCK/T*38546.3)
  CASE (113)  ! K  V
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*24012.5)+6.*EXP(-HCK/T*24249.6)+ 2.*EXP(-HCK/T*39758.1)+4.*EXP(-HCK/T*40080.2)
  CASE (114)  ! K  VI
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*1132.)+5.*EXP(-HCK/T*2924.)+ 5.*EXP(-HCK/T*18977.8)+1.*EXP(-HCK/T*43358.8)

  ! --- Z=29  Cu ---
  CASE (169); PFGROUND_KURUCZ = 2.                    ! Cu I
  CASE (170); PFGROUND_KURUCZ = 1.                    ! Cu II
  CASE (171)  ! Cu III
    PFGROUND_KURUCZ = 6.+4.*EXP(-HCK/T*2071.8)
  CASE (172); PFGROUND_KURUCZ = 1.                    ! Cu IV
  CASE (173); PFGROUND_KURUCZ = 1.                    ! Cu V
  CASE (174); PFGROUND_KURUCZ = 1.                    ! Cu VI

  ! --- Z=30  Zn ---
  CASE (175); PFGROUND_KURUCZ = 1.                    ! Zn I
  CASE (176); PFGROUND_KURUCZ = 2.                    ! Zn II
  CASE (177); PFGROUND_KURUCZ = 1.                    ! Zn III
  CASE (178); PFGROUND_KURUCZ = 1.                    ! Zn IV
  CASE (179); PFGROUND_KURUCZ = 1.                    ! Zn V
  CASE (180); PFGROUND_KURUCZ = 1.                    ! Zn VI

  ! --- Z=31  Ga ---
  CASE (181)  ! Ga I
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*826.19)
  CASE (182); PFGROUND_KURUCZ = 1.                    ! Ga II
  CASE (183); PFGROUND_KURUCZ = 2.                    ! Ga III
  CASE (184); PFGROUND_KURUCZ = 1.                    ! Ga IV
  CASE (185); PFGROUND_KURUCZ = 1.                    ! Ga V
  CASE (186); PFGROUND_KURUCZ = 1.                    ! Ga VI

  ! --- Z=32  Ge ---
  CASE (187)  ! Ge I
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*557.134)+5.*EXP(-HCK/T*1409.961)+ 5.*EXP(-HCK/T*7125.299)
  CASE (188)  ! Ge II
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*1767.356)
  CASE (189); PFGROUND_KURUCZ = 1.                    ! Ge III
  CASE (190); PFGROUND_KURUCZ = 1.                    ! Ge IV
  CASE (191); PFGROUND_KURUCZ = 1.                    ! Ge V
  CASE (192); PFGROUND_KURUCZ = 1.                    ! Ge VI

  ! --- Z=33  As ---
  CASE (193); PFGROUND_KURUCZ = 4.                    ! As I
  CASE (194)  ! As II
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*1061.)+5.*EXP(-HCK/T*2538.)
  CASE (195)  ! As III
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*2940.)
  CASE (196); PFGROUND_KURUCZ = 1.                    ! As IV
  CASE (197); PFGROUND_KURUCZ = 1.                    ! As V
  CASE (198); PFGROUND_KURUCZ = 1.                    ! As VI

  ! --- Z=34  Se ---
  CASE (199)  ! Se I
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*1989.49)+EXP(-HCK/T*2534.35)+ 5.*EXP(-HCK/T*9576.149)+1.*EXP(-HCK/T*22446.202)
  CASE (200)  ! Se II
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*13168.2)+6.*EXP(-HCK/T*13784.4)+ 2.*EXP(-HCK/T*23038.3)+4.*EXP(-HCK/T*23894.8)
  CASE (201)  ! Se III
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*1741.)+5.*EXP(-HCK/T*3937.)+ 5.*EXP(-HCK/T*13032.)+1.*EXP(-HCK/T*28430.)
  CASE (202); PFGROUND_KURUCZ = 1.                    ! Se IV
  CASE (203); PFGROUND_KURUCZ = 1.                    ! Se V
  CASE (204); PFGROUND_KURUCZ = 1.                    ! Se VI

  ! --- Z=35  Br ---
  CASE (205)  ! Br I
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*3685.24)
  CASE (206)  ! Br II
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*3136.4)+EXP(-HCK/T*3837.5)+ 5.*EXP(-HCK/T*12089.1)
  CASE (207)  ! Br III
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*15042.0)+6.*EXP(-HCK/T*16301.0)+ 2.*EXP(-HCK/T*26915.0)+4.*EXP(-HCK/T*28579.0)
  CASE (208); PFGROUND_KURUCZ = 1.                    ! Br IV
  CASE (209); PFGROUND_KURUCZ = 1.                    ! Br V
  CASE (210); PFGROUND_KURUCZ = 1.                    ! Br VI

  ! --- Z=36  Kr ---
  CASE (211); PFGROUND_KURUCZ = 1.                    ! Kr I
  CASE (212)  ! Kr II
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*5371.)
  CASE (213)  ! Kr III
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*3136.4)+EXP(-HCK/T*3837.5)+ 5.*EXP(-HCK/T*14644.3)+1.*EXP(-HCK/T*33079.6)
  CASE (214); PFGROUND_KURUCZ = 1.                    ! Kr IV
  CASE (215); PFGROUND_KURUCZ = 1.                    ! Kr V
  CASE (216); PFGROUND_KURUCZ = 1.                    ! Kr VI

  ! --- Z=37  Rb ---
  CASE (217); PFGROUND_KURUCZ = 2.                    ! Rb I
  CASE (218); PFGROUND_KURUCZ = 1.                    ! Rb II
  CASE (219)  ! Rb III
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*7380.)
  CASE (220); PFGROUND_KURUCZ = 1.                    ! Rb IV
  CASE (221); PFGROUND_KURUCZ = 1.                    ! Rb V
  CASE (222); PFGROUND_KURUCZ = 1.                    ! Rb VI

  ! --- Z=38  Sr ---
  CASE (223); PFGROUND_KURUCZ = 1.                    ! Sr I
  CASE (224)  ! Sr II
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*14555.50)+6.*EXP(-HCK/T*14836.24)
  CASE (225); PFGROUND_KURUCZ = 1.                    ! Sr III
  CASE (226); PFGROUND_KURUCZ = 1.                    ! Sr IV
  CASE (227); PFGROUND_KURUCZ = 1.                    ! Sr V
  CASE (228); PFGROUND_KURUCZ = 1.                    ! Sr VI

  ! --- Z=39  Y  ---
  CASE (229)  ! Y  I
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*530.36)
  CASE (230)  ! Y  II
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*840.198)+5.*EXP(-HCK/T*1045.076)+ 7.*EXP(-HCK/T*1449.752)+5.*EXP(-HCK/T*3296.280)+ &
      5.*EXP(-HCK/T*8003.126)+7.*EXP(-HCK/T*8328.039)+ 9.*EXP(-HCK/T*8743.322)
  CASE (231)  ! Y  III
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*724.15)+2.*EXP(-HCK/T*7467.10)
  CASE (232); PFGROUND_KURUCZ = 1.                    ! Y  IV
  CASE (233); PFGROUND_KURUCZ = 1.                    ! Y  V
  CASE (234); PFGROUND_KURUCZ = 1.                    ! Y  VI

  ! --- Z=40  Zr ---
  CASE (235)  ! Zr I
    PFGROUND_KURUCZ = 5.+7.*EXP(-HCK/T*570.41)+9.*EXP(-HCK/T*1240.84)+ 1.*EXP(-HCK/T*4196.85)+3.*EXP(-HCK/T*4376.28)+ &
      5.*EXP(-HCK/T*4186.11)+3.*EXP(-HCK/T*4870.53)+ 5.*EXP(-HCK/T*5023.41)+7.*EXP(-HCK/T*5249.07)+ 9.*EXP(-HCK/T*5540.54)+ &
      11.*EXP(-HCK/T*5888.93)+ 5.*EXP(-HCK/T*5101.68)+9.*EXP(-HCK/T*8057.30)
  CASE (236)  ! Zr II
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*314.67)+8.*EXP(-HCK/T*763.44)+ 10.*EXP(-HCK/T*1322.91)+4.*EXP(-HCK/T*2572.21)+ &
      6.*EXP(-HCK/T*2895.00)+8.*EXP(-HCK/T*3299.58)+ 10.*EXP(-HCK/T*3757.63)+4.*EXP(-HCK/T*4247.97)+ 6.*EXP(-HCK/T*4505.30)+ &
      2.*EXP(-HCK/T*5723.78)+ 4.*EXP(-HCK/T*6111.16)+6.*EXP(-HCK/T*5752.55)+ 8.*EXP(-HCK/T*6467.10)+2.*EXP(-HCK/T*7512.61)+ &
      4.*EXP(-HCK/T*7736.05)+6.*EXP(-HCK/T*8058.27)+ 8.*EXP(-HCK/T*7837.49)+10.*EXP(-HCK/T*8152.57)+ 2.*EXP(-HCK/T*9553.13)+ &
      4.*EXP(-HCK/T*9742.80)+ 6.*EXP(-HCK/T*9968.75)
  CASE (237)  ! Zr III
    PFGROUND_KURUCZ = 5.+7.*EXP(-HCK/T*681.2)+9.*EXP(-HCK/T*1485.8)+ 5.*EXP(-HCK/T*5742.8)+1.*EXP(-HCK/T*8062.7)+ &
      3.*EXP(-HCK/T*8327.0)+5.*EXP(-HCK/T*8839.7)+ 9.*EXP(-HCK/T*11049.9)+1.*EXP(-HCK/T*23974.9)+ 3.*EXP(-HCK/T*18400.8)+ &
      5.*EXP(-HCK/T*18804.7)+ 7.*EXP(-HCK/T*19535.3)+5.*EXP(-HCK/T*25066.9)+ 1.*EXP(-HCK/T*36473.7)
  CASE (238); PFGROUND_KURUCZ = 1.                    ! Zr IV
  CASE (239); PFGROUND_KURUCZ = 1.                    ! Zr V
  CASE (240); PFGROUND_KURUCZ = 1.                    ! Zr VI

  ! --- Z=41  Nb ---
  CASE (241)  ! Nb I
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*154.19)+6.*EXP(-HCK/T*391.99)+ 8.*EXP(-HCK/T*695.25)+10.*EXP(-HCK/T*1050.26)+ &
      4.*EXP(-HCK/T*1142.79)+6.*EXP(-HCK/T*1586.90)+ 8.*EXP(-HCK/T*2154.11)+10.*EXP(-HCK/T*2805.36)+ 2.*EXP(-HCK/T*4998.17)+ &
      4.*EXP(-HCK/T*5297.92)+ 6.*EXP(-HCK/T*5965.45)+2.*EXP(-HCK/T*8410.90)+ 4.*EXP(-HCK/T*8705.32)+6.*EXP(-HCK/T*9043.14)+ &
      8.*EXP(-HCK/T*9497.52)+8.*EXP(-HCK/T*8827.00)+ 10.*EXP(-HCK/T*9328.88)+4.*EXP(-HCK/T*9439.08)+ 6.*EXP(-HCK/T*10237.51)
  CASE (242)  ! Nb II
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*158.99)+5.*EXP(-HCK/T*438.38)+ 7.*EXP(-HCK/T*801.38)+9.*EXP(-HCK/T*1224.87)+ &
      3.*EXP(-HCK/T*2356.76)+5.*EXP(-HCK/T*2629.07)+ 7.*EXP(-HCK/T*3029.57)+9.*EXP(-HCK/T*3542.50)+ 11.*EXP(-HCK/T*4146.00)+ &
      1.*EXP(-HCK/T*5562.26)+ 3.*EXP(-HCK/T*6192.33)+5.*EXP(-HCK/T*7261.33)+ 5.*EXP(-HCK/T*7505.78)+7.*EXP(-HCK/T*7900.65)+ &
      9.*EXP(-HCK/T*8320.40)+9.*EXP(-HCK/T*9509.67)+ 11.*EXP(-HCK/T*9812.56)+13.*EXP(-HCK/T*10186.41)
  CASE (243)  ! Nb III
    PFGROUND_KURUCZ = 4.+6.*EXP(-HCK/T*515.8)+8.*EXP(-HCK/T*1176.6)+ 10.*EXP(-HCK/T*1939.0)+ 2.*EXP(-HCK/T*8664.3)+ &
      4.*EXP(-HCK/T*8607.5)+ 6.*EXP(-HCK/T*9593.7)+8.*EXP(-HCK/T*9236.1)+ 10.*EXP(-HCK/T*9804.5)+4.*EXP(-HCK/T*10912.2)+ &
      6.*EXP(-HCK/T*13094.0)+10.*EXP(-HCK/T*12916.0)+ 12.*EXP(-HCK/T*13263.8)+6.*EXP(-HCK/T*19975.0)+ 8.*EXP(-HCK/T*19861.0)+ &
      4.*EXP(-HCK/T*25220.2)+ 6.*EXP(-HCK/T*25735.2)+8.*EXP(-HCK/T*26463.7)+ 10.*EXP(-HCK/T*27373.5)
  CASE (244); PFGROUND_KURUCZ = 1.                    ! Nb IV
  CASE (245); PFGROUND_KURUCZ = 1.                    ! Nb V
  CASE (246); PFGROUND_KURUCZ = 1.                    ! Nb VI

  ! --- Z=42  Mo ---
  CASE (247); PFGROUND_KURUCZ = 7.                    ! Mo I
  CASE (248); PFGROUND_KURUCZ = 6.                    ! Mo II
  CASE (249)  ! Mo III
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*243.10)+5.*EXP(-HCK/T*669.60)+ 7.*EXP(-HCK/T*1225.20)+9.*EXP(-HCK/T*1873.80)
  CASE (250); PFGROUND_KURUCZ = 1.                    ! Mo IV
  CASE (251); PFGROUND_KURUCZ = 1.                    ! Mo V
  CASE (252); PFGROUND_KURUCZ = 1.                    ! Mo VI

  ! --- Z=43  Tc ---
  CASE (253)  ! Tc I
    PFGROUND_KURUCZ = 6.   +10.*EXP(-HCK/T*2572.89)+8.*EXP(-HCK/T*3250.91)+ 6.*EXP(-HCK/T*3700.54)+4.*EXP(-HCK/T*4002.57)+ &
      2.*EXP(-HCK/T*4178.75)
  CASE (254)  ! Tc II
    PFGROUND_KURUCZ = 7.  +9.*EXP(-HCK/T*3461.27)+7.*EXP(-HCK/T*4217.17)+ 5.*EXP(-HCK/T*4669.22)+3.*EXP(-HCK/T*4961.14)+ &
      1.*EXP(-HCK/T*5100.98)
  CASE (255); PFGROUND_KURUCZ = 6.                    ! Tc III
  CASE (256); PFGROUND_KURUCZ = 1.                    ! Tc IV
  CASE (257); PFGROUND_KURUCZ = 1.                    ! Tc V
  CASE (258); PFGROUND_KURUCZ = 1.                    ! Tc VI

  ! --- Z=44  Ru ---
  CASE (259)  ! Ru I
    PFGROUND_KURUCZ = 11.+9.*EXP(-HCK/T*1190.64)+7.*EXP(-HCK/T*2091.54)+ 5.*EXP(-HCK/T*2713.24)+3.*EXP(-HCK/T*3105.49)+ &
      9.*EXP(-HCK/T*6545.03)+7.*EXP(-HCK/T*8084.12)+ 5.*EXP(-HCK/T*9183.66)+9.*EXP(-HCK/T*7483.07)+ 7.*EXP(-HCK/T*8575.42)+ &
      5.*EXP(-HCK/T*9057.64)+ 3.*EXP(-HCK/T*9072.98)+1.*EXP(-HCK/T*8492.37)+ 7.*EXP(-HCK/T*8770.93)+5.*EXP(-HCK/T*8043.69)+ &
      3.*EXP(-HCK/T*9620.29)+9.*EXP(-HCK/T*9120.63)
  CASE (260)  ! Ru II
    PFGROUND_KURUCZ = 10.+8.*EXP(-HCK/T*1523.1)+6.*EXP(-HCK/T*2493.9)+ 4.*EXP(-HCK/T*3104.2)+6.*EXP(-HCK/T*8256.7)+ &
      4.*EXP(-HCK/T*8477.4)+2.*EXP(-HCK/T*9373.4)+ 10.*EXP(-HCK/T*9151.6)
  CASE (261)  ! Ru III
    PFGROUND_KURUCZ = 9.+7.*EXP(-HCK/T*1158.8)+5.*EXP(-HCK/T*1826.3)+ 3.*EXP(-HCK/T*2266.3)+EXP(-HCK/T*2476.0)+ 7.*EXP(-HCK/T*27162.8)+ &
      5.*EXP(-HCK/T*41111.7)
  CASE (262); PFGROUND_KURUCZ = 1.                    ! Ru IV
  CASE (263); PFGROUND_KURUCZ = 1.                    ! Ru V
  CASE (264); PFGROUND_KURUCZ = 1.                    ! Ru VI

  ! --- Z=45  Rh ---
  CASE (265)  ! Rh I
    PFGROUND_KURUCZ = 10.+8.*EXP(-HCK/T*1529.97)+6.*EXP(-HCK/T*2598.03)+ 4.*EXP(-HCK/T*3472.68)+6.*EXP(-HCK/T*3309.86)+ &
      4.*EXP(-HCK/T*5657.97)+8.*EXP(-HCK/T*5690.97)+ 6.*EXP(-HCK/T*7791.23)+6.*EXP(-HCK/T*9221.22)
  CASE (266)  ! Rh II
    PFGROUND_KURUCZ = 9.+7.*EXP(-HCK/T*2401.3)+5.*EXP(-HCK/T*3580.7)+ 5.*EXP(-HCK/T*8164.4)+1.*EXP(-HCK/T*10760.8)+ &
      3.*EXP(-HCK/T*10515.0)+5.*EXP(-HCK/T*11643.7)+ 9.*EXP(-HCK/T*14855.4)+11.*EXP(-HCK/T*16884.8)+ 9.*EXP(-HCK/T*18540.4)+ &
      7.*EXP(-HCK/T*19792.4)
  CASE (267)  ! Rh III
    PFGROUND_KURUCZ = 10.+8.*EXP(-HCK/T*2147.8)+6.*EXP(-HCK/T*3485.7)+ 4.*EXP(-HCK/T*4322.0)+6.*EXP(-HCK/T*11062.3)+ &
      4.*EXP(-HCK/T*10997.1)+2.*EXP(-HCK/T*12469.8)+ 10.*EXP(-HCK/T*14044.0)+8.*EXP(-HCK/T*15256.8)+ 4.*EXP(-HCK/T*16870.7)+ &
      2.*EXP(-HCK/T*18303.7)+ 12.*EXP(-HCK/T*19490.2)+6.*EXP(-HCK/T*19528.5)
  CASE (268); PFGROUND_KURUCZ = 1.                    ! Rh IV
  CASE (269); PFGROUND_KURUCZ = 1.                    ! Rh V
  CASE (270); PFGROUND_KURUCZ = 1.                    ! Rh VI

  ! --- Z=46  Pd ---
  CASE (271)  ! Pd I
    PFGROUND_KURUCZ = 1.+7.*EXP(-HCK/T*6564.11)+5.*EXP(-HCK/T*7754.99)
  CASE (272)  ! Pd II
    PFGROUND_KURUCZ = 6.+4.*EXP(-HCK/T*3539.2)
  CASE (273)  ! Pd III
    PFGROUND_KURUCZ = 9.+7.*EXP(-HCK/T*3229.3)+5.*EXP(-HCK/T*4687.5)+ 5.*EXP(-HCK/T*10229.3)+3.*EXP(-HCK/T*13468.9)+ &
      1.*EXP(-HCK/T*13697.5)+5.*EXP(-HCK/T*14634.4)+ 9.*EXP(-HCK/T*17879.3)
  CASE (274); PFGROUND_KURUCZ = 1.                    ! Pd IV
  CASE (275); PFGROUND_KURUCZ = 1.                    ! Pd V
  CASE (276); PFGROUND_KURUCZ = 1.                    ! Pd VI

  ! --- Z=47  Ag ---
  CASE (277); PFGROUND_KURUCZ = 2.                    ! Ag I
  CASE (278); PFGROUND_KURUCZ = 1.                    ! Ag II
  CASE (279)  ! Ag III
    PFGROUND_KURUCZ = 6.+4.*EXP(-HCK/T*4607.)
  CASE (280); PFGROUND_KURUCZ = 1.                    ! Ag IV
  CASE (281); PFGROUND_KURUCZ = 1.                    ! Ag V
  CASE (282); PFGROUND_KURUCZ = 1.                    ! Ag VI

  ! --- Z=48  Cd ---
  CASE (283); PFGROUND_KURUCZ = 1.                    ! Cd I
  CASE (284); PFGROUND_KURUCZ = 2.                    ! Cd II
  CASE (285); PFGROUND_KURUCZ = 1.                    ! Cd III
  CASE (286); PFGROUND_KURUCZ = 1.                    ! Cd IV
  CASE (287); PFGROUND_KURUCZ = 1.                    ! Cd V
  CASE (288); PFGROUND_KURUCZ = 1.                    ! Cd VI

  ! --- Z=49  In ---
  CASE (289)  ! In I
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*2212.598)
  CASE (290); PFGROUND_KURUCZ = 1.                    ! In II
  CASE (291); PFGROUND_KURUCZ = 2.                    ! In III
  CASE (292); PFGROUND_KURUCZ = 1.                    ! In IV
  CASE (293); PFGROUND_KURUCZ = 1                    ! In V
  CASE (294); PFGROUND_KURUCZ = 1.                    ! In VI

  ! --- Z=50  Sn ---
  CASE (295)  ! Sn I
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*1691.8)+5.*EXP(-HCK/T*3427.7)+ 5.*EXP(-HCK/T*6513.0)
  CASE (296)  ! Sn II
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*4251.4)
  CASE (297); PFGROUND_KURUCZ = 1.                    ! Sn III
  CASE (298); PFGROUND_KURUCZ = 1.                    ! Sn IV
  CASE (299); PFGROUND_KURUCZ = 1.                    ! Sn V
  CASE (300); PFGROUND_KURUCZ = 1.                    ! Sn VI

  ! --- Z=51  Sb ---
  CASE (301)  ! Sb I
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*8512.1)+6.*EXP(-HCK/T*9854.1)
  CASE (302)  ! Sb II
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*3055.0)+5.*EXP(-HCK/T*5659.0)
  CASE (303)  ! Sb III
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*6576.)
  CASE (304); PFGROUND_KURUCZ = 1.                    ! Sb IV
  CASE (305); PFGROUND_KURUCZ = 1.                    ! Sb V
  CASE (306); PFGROUND_KURUCZ = 1.                    ! Sb VI

  ! --- Z=52  Te ---
  CASE (307)  ! Te I
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*4750.712)+EXP(-HCK/T*4706.5)
  CASE (308)  ! Te II
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*10222.385)+6.*EXP(-HCK/T*12421.854)+ 2.*EXP(-HCK/T*20546.591)+4.*EXP(-HCK/T*24032.2)
  CASE (309)  ! Te III
    PFGROUND_KURUCZ = 1.+3.*EXP(-HCK/T*4756.5)+5.*EXP(-HCK/T*8166.9)+ 5.*EXP(-HCK/T*17358.)
  CASE (310); PFGROUND_KURUCZ = 1.                    ! Te IV
  CASE (311); PFGROUND_KURUCZ = 1.                    ! Te V
  CASE (312); PFGROUND_KURUCZ = 1.                    ! Te VI

  ! --- Z=53  I  ---
  CASE (313)  ! I  I
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*7063.15)
  CASE (314)  ! I  II
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*7087.0)+EXP(-HCK/T*6447.9)+ 5.*EXP(-HCK/T*13727.2)+1.*EXP(-HCK/T*29501.3)
  CASE (315)  ! I  III
    PFGROUND_KURUCZ = 4.+4.*EXP(-HCK/T*11711.2)+6.*EXP(-HCK/T*14901.9)+ 2.*EXP(-HCK/T*24299.3)+4.*EXP(-HCK/T*29636.8)
  CASE (316); PFGROUND_KURUCZ = 1.                    ! I  IV
  CASE (317); PFGROUND_KURUCZ = 1.                    ! I  V
  CASE (318); PFGROUND_KURUCZ = 1.                    ! I  VI

  ! --- Z=54  Xe ---
  CASE (319); PFGROUND_KURUCZ = 1.                    ! Xe I
  CASE (320)  ! Xe II
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*10537.01)
  CASE (321)  ! Xe III
    PFGROUND_KURUCZ = 5.+3.*EXP(-HCK/T*9794.36)+EXP(-HCK/T*8130.08)+ 5.*EXP(-HCK/T*17098.73)+1.*EXP(-HCK/T*36102.94)
  CASE (322); PFGROUND_KURUCZ = 1.                    ! Xe IV
  CASE (323); PFGROUND_KURUCZ = 1.                    ! Xe V
  CASE (324); PFGROUND_KURUCZ = 1.                    ! Xe VI

  ! --- Z=55  Cs ---
  CASE (325); PFGROUND_KURUCZ = 2.                    ! Cs I
  CASE (326); PFGROUND_KURUCZ = 1.                    ! Cs II
  CASE (327)  ! Cs III
    PFGROUND_KURUCZ = 4.+2.*EXP(-HCK/T*13884.)
  CASE (328); PFGROUND_KURUCZ = 1.                    ! Cs IV
  CASE (329); PFGROUND_KURUCZ = 1.                    ! Cs V
  CASE (330); PFGROUND_KURUCZ = 1.                    ! Cs VI

  ! --- Z=56  Ba ---
  CASE (331); PFGROUND_KURUCZ = 1.                    ! Ba I
  CASE (332)  ! Ba II
    PFGROUND_KURUCZ = 2.+4.*EXP(-HCK/T*4873.852)+6.*EXP(-HCK/T*5674.807)
  CASE (333); PFGROUND_KURUCZ = 1.                    ! Ba III
  CASE (334); PFGROUND_KURUCZ = 1.                    ! Ba IV
  CASE (335); PFGROUND_KURUCZ = 1.                    ! Ba V
  CASE (336); PFGROUND_KURUCZ = 1.                    ! Ba VI

  ! Default for all other NELION values
  CASE DEFAULT
    PFGROUND_KURUCZ = 1.0d0

  END SELECT

END FUNCTION PFGROUND_KURUCZ

!=======================================================================
! FUNCTION PFGROUND_HYBRID(NELION, T)
!
!   Atomic partition function for the ion identified by
!     NELION = (Z - 1) * 6 + ION,  ION in 1..6
!   at temperature T [K].  Returns a LINEAR partition function (not log10).
!
!   Hybrid dispatcher over two data sources:
!     T <= 9000 K          : Barklem & Collet (2016), via U_BC().
!     9000 K < T < 10000 K : smoothstep blend of U_BC and PFGROUND_KURUCZ.
!     T >= 10000 K         : Kurucz hand-coded values, via PFGROUND_KURUCZ().
!
!   The blend region is C^1-continuous: both value and first derivative
!   match at the blend endpoints, which keeps the ATLAS atmosphere model
!   iterating smoothly across the transition.  (A hard switch at 10000 K
!   would create a small-but-real discontinuity in PF, which ENERGY_DENSITY
!   would amplify into a large spurious adiabatic-gradient contribution
!   via its finite-difference d(ln U)/d(ln T) calculation.)
!
!   For species outside B&C coverage (ION > 3, or (Z, ION) not tabulated,
!   or T > 10000 K) the output is pristine Kurucz -- identical to the
!   original ATLAS12.  For species inside B&C coverage at T < 9000 K,
!   the output reflects modern NIST-level-summation partition functions.
!
!   Iron group (Z=20..28, NELION=115..168): returns 1.0 and defers to
!   PFIRON, matching pristine behavior.
!
!   Callers: PFSAHA uses this as a low-T safety floor via
!     PART(ION) = MAX(PFGROUND_HYBRID(...), PART(ION))
!   in the NNN-polynomial branch.
!=======================================================================

FUNCTION PFGROUND_HYBRID(NELION, T)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)  :: NELION
  REAL(8),  INTENT(IN)  :: T
  REAL(8)               :: PFGROUND_HYBRID

  ! Locals
  INTEGER :: IZ, ION
  REAL(8)  :: U_bc_val, U_kurucz_val, x, w

  ! Blend window [K].  Chosen so that T_HIGH_BC matches the top of the
  ! B&C tabulated grid (10000 K) -- no extrapolation ever.
  REAL(8), PARAMETER :: T_LOW_BC  =  9000.0d0
  REAL(8), PARAMETER :: T_HIGH_BC = 10000.0d0

  IF (IDEBUG .EQ. 1) WRITE(6,'(A)') ' RUNNING PFGROUND_HYBRID'

  ! --- Iron group deferred to PFIRON; short-circuit ---
  IF (NELION .GE. 115 .AND. NELION .LE. 168) THEN
    PFGROUND_HYBRID = 1.0d0
    RETURN
  END IF

  ! --- Above the blend window: pure Kurucz ---
  IF (T .GE. T_HIGH_BC) THEN
    PFGROUND_HYBRID = PFGROUND_KURUCZ(NELION, T)
    RETURN
  END IF

  ! --- Decode NELION into (IZ, ION); out-of-range ION -> pure Kurucz ---
  IZ  = (NELION - 1) / 6 + 1
  ION = NELION - (IZ - 1) * 6
  IF (ION .LT. 1 .OR. ION .GT. 3) THEN
    PFGROUND_HYBRID = PFGROUND_KURUCZ(NELION, T)
    RETURN
  END IF

  ! --- Ask B&C for its value; sentinel means no B&C data for this species ---
  U_bc_val = U_BC(IZ, ION, T)
  IF (U_bc_val .LE. 0.0d0) THEN
    ! No B&C value available; use Kurucz everywhere, including the blend
    ! window.  (U_BC returns <0 for out-of-coverage species regardless of
    ! T, so this branch is entered for Z>92 or for those (Z, ION) combos
    ! B&C doesn't tabulate.)
    PFGROUND_HYBRID = PFGROUND_KURUCZ(NELION, T)
    RETURN
  END IF

  ! --- Below the blend window: pure B&C ---
  IF (T .LE. T_LOW_BC) THEN
    PFGROUND_HYBRID = U_bc_val
    RETURN
  END IF

  ! --- Blend region: smoothstep from B&C (w=0) at T_LOW_BC to Kurucz
  !     (w=1) at T_HIGH_BC.  w = 3x^2 - 2x^3 with x in [0,1]:
  !       w(0) = 0, w(1) = 1
  !       w'(0) = 0, w'(1) = 0
  !     gives C^1 continuity at both blend endpoints.
  U_kurucz_val = PFGROUND_KURUCZ(NELION, T)
  x = (T - T_LOW_BC) / (T_HIGH_BC - T_LOW_BC)
  w = x * x * (3.0d0 - 2.0d0 * x)
  PFGROUND_HYBRID = (1.0d0 - w) * U_bc_val + w * U_kurucz_val

END FUNCTION PFGROUND_HYBRID



!=======================================================================
! FUNCTION PARTFNH2(T)
!
! H2 partition function by linear interpolation in temperature.
!
! Reads a two-column text file (T, Q) on first call. The file
! must have comment lines starting with '#' and data lines with
! T(K) and Q(H2) in free format, sorted by ascending T.
! The table is assumed to be on a uniform 100K grid starting at 100K.
!
! Default file: 'partfnh2_table.dat' (can be replaced with updated
! data, e.g. from Barklem & Collet 2016).
!
! Original hardcoded table: Kurucz, R.L. 1985, CfA preprint no. 2162.
!=======================================================================

FUNCTION PARTFNH2(T)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: T
  REAL(8) :: PARTFNH2

  INTEGER, PARAMETER :: NMAX = 500   ! max table entries
  REAL(8),  SAVE :: PF(NMAX)          ! partition function values
  REAL(8),  SAVE :: TSTEP = 100.0d0   ! temperature step (K)
  REAL(8),  SAVE :: TSTART = 100.0d0  ! first temperature in table
  INTEGER, SAVE :: NPF = 0           ! number of entries loaded
  LOGICAL, SAVE :: INITIALIZED = .FALSE.
  
  INTEGER :: N, IOS, LUN
  REAL(8)  :: frac, TDUM, QDUM
  CHARACTER(len=256) :: LINE

  ! Read table from file on first call
  IF (.NOT. INITIALIZED) THEN
    INITIALIZED = .TRUE.
    NPF = 0
    LUN = 89
    OPEN(LUN, FILE=trim(DATADIR)//'partfnh2.dat', STATUS='old', IOSTAT=IOS)
    IF (IOS .NE. 0) THEN
      STOP ' PARTFNH2 ERROR: cannot open '//trim(DATADIR)//'partfnh2.dat'
   ENDIF
   DO WHILE (NPF .LT. NMAX)
      READ(LUN, '(A)', IOSTAT=IOS) LINE
      IF (IOS .NE. 0) EXIT
      ! Skip comment and blank lines
      LINE = adjustl(LINE)
      IF (LINE(1:1) .EQ. '#' .OR. len_trim(LINE) .EQ. 0) CYCLE
      READ(LINE, *, IOSTAT=IOS) TDUM, QDUM
      IF (IOS .NE. 0) CYCLE
      NPF = NPF + 1
      PF(NPF) = QDUM
      IF (NPF .EQ. 1) TSTART = TDUM
      IF (NPF .EQ. 2) TSTEP = TDUM - TSTART
   END DO
   CLOSE(LUN)
  END IF

  ! Fallback if no table loaded
  IF (NPF .EQ. 0) THEN
    PARTFNH2 = max(T / 100.0d0, 0.5d0)
    RETURN
  END IF

  ! Table index from temperature
  N = INT((T - TSTART) / TSTEP) + 1
  N = MIN(NPF - 1, MAX(1, N))

  ! Linear interpolation
  frac = (T - (TSTART + (N-1) * TSTEP)) / TSTEP
  frac = MIN(1.0d0, MAX(0.0d0, frac))
  PARTFNH2 = PF(N) + (PF(N+1) - PF(N)) * frac

END FUNCTION PARTFNH2

!=======================================================================
! FUNCTION EQUILH2(T)
!
! H2 equilibrium constant K(T) for the reaction H + H ⇌ H2.
!
! Returns K(T) such that n(H2) = n(H)^2 * K(T), where:
!
!   K(T) = Q_H2(T) * (2^1.5 / 4) / (2*pi*m_H2*k*T/h^2)^1.5
!          * exp(D0 * h*c / (k*T))
!
! Q_H2(T) is the internal partition function from PARTFNH2.
! D0(H2) = 36118.11 cm^-1 is the dissociation energy.
!
! Reference: Kurucz, R.L. 1985, CfA preprint no. 2162.
! Revised 8 Nov 2005: extended to 20000 K for white dwarfs.
!=======================================================================

FUNCTION EQUILH2(T)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: T
  REAL(8) :: EQUILH2

  ! External function

  ! Physical constants from mod_constants; only local spectroscopic data here
  REAL(8), PARAMETER :: m_H  = 1.008d0 * AMU        ! H atom mass [g]

  ! H2 dissociation energy
  REAL(8), PARAMETER :: D0_cm = 36118.11d0           ! D0(H2) [cm^-1]
  REAL(8), PARAMETER :: D0_over_kT_coeff = D0_cm * HCK
  !                                     = 51967.8 K  (D0/k)

  ! Translational partition function prefactor (T-independent part)
  ! = 2^1.5 / 4 / (2*pi*m_H*k/h^2)^1.5
  ! where m_H appears because the reduced mass for equal-mass dissociation
  ! products cancels to give the atomic H mass in the Saha-like expression.
  REAL(8), PARAMETER :: trans_prefactor = 2.d0**1.5d0 / 4.d0 &
    / (2.d0 * PI * m_H * KBOL / HPLANCK**2)**1.5d0

  REAL(8) :: Q_H2

  ! Internal partition function
  Q_H2 = PARTFNH2(T)

  ! Equilibrium constant
  EQUILH2 = Q_H2 * trans_prefactor / T**1.5d0 * EXP(D0_over_kT_coeff / T)

END FUNCTION EQUILH2


!=======================================================================
! FUNCTION holtsmark_Q(beta)
!
! Holtsmark cumulative microfield distribution Q(beta).
!
! Returns the probability that the electric microfield strength
! F < beta * F_0, where F_0 is the Holtsmark normal field strength.
! Uses linear interpolation on the pre-tabulated Q_HOLTSMARK grid
! (uniform in log10(beta), stored in mod_constants).
!
! For beta < 0.01: returns 0 (field always below critical → bound)
! For beta > 50:   returns 1 (field always below critical → bound)
!
! Note on sign convention: large beta_crit means the critical field
! is much larger than the typical field, so Q ~ 1 and the level
! survives (is bound).  Small beta_crit means the level is dissolved.
!
! Reference: Holtsmark, J. 1919, Ann. Phys. 363, 577
!=======================================================================

FUNCTION holtsmark_Q(beta)

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: beta
  REAL(8) :: holtsmark_Q

  REAL(8)  :: log_beta, idx_f, frac
  INTEGER :: idx

  IF (beta .LE. 0.01D0) THEN
    holtsmark_Q = 0.0D0
    RETURN
  END IF
  IF (beta .GE. 50.0D0) THEN
    holtsmark_Q = 1.0D0
    RETURN
  END IF

  log_beta = LOG10(beta)
  idx_f = (log_beta - LOG_BETA_MIN) / LOG_BETA_STEP
  idx = INT(idx_f)
  idx = MAX(0, MIN(NQ_HOLTSMARK - 2, idx))
  frac = idx_f - DBLE(idx)

  holtsmark_Q = Q_HOLTSMARK(idx + 1) &
              + frac * (Q_HOLTSMARK(idx + 2) - Q_HOLTSMARK(idx + 1))

END FUNCTION holtsmark_Q


!=======================================================================
! FUNCTION occupation_prob(n, xne)
!
! Occupation probability for hydrogen level n in the Hummer & Mihalas
! (1988) formalism.
!
! Returns w_n = Q(beta_crit), where:
!   beta_crit = BETA_COEFF_HM88 / (n^5 * N_e^(2/3))
!   Q = Holtsmark cumulative microfield distribution
!
! w_n = 1 means the level is fully bound (low density or low n).
! w_n = 0 means the level is fully dissolved (high density or high n).
!
! For hydrogen lines, the line opacity should be multiplied by w_n
! for the upper level.  The "missing" opacity (1 - w_n) appears as
! pseudo-continuum opacity (handled separately in HOP).
!
! Arguments:
!   n   — principal quantum number (integer, >= 1)
!   xne — electron number density [cm^-3]
!
! References:
!   Hummer, D.G. & Mihalas, D. 1988, ApJ 331, 794
!   Nayfonov, A. et al. 1999, ApJ 526, 451
!=======================================================================

FUNCTION occupation_prob(n, xne)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(8),  INTENT(IN) :: xne
  REAL(8) :: occupation_prob

  REAL(8) :: beta

  IF (n .LE. 1 .OR. xne .LE. 0.0D0) THEN
    occupation_prob = 1.0D0
    RETURN
  END IF

  beta = BETA_COEFF_HM88 / (DBLE(n)**5 * xne**(2.0D0/3.0D0))
  occupation_prob = holtsmark_Q(beta)

END FUNCTION occupation_prob


END MODULE mod_atlas_data

