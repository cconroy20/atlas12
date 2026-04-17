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
!    mod_parameters          Compile-time constants (kw, mion, maxmol, ...)
!    mod_atlas_data          All shared state (replaces 58 COMMON blocks)
!
!  MAIN PROGRAM
!    ATLAS12                 Iteration driver: read input, iterate atmosphere
!
!  MODEL OUTPUT
!    PUTOUT                  Write model atmosphere, fluxes, intensities
!
!  TEMPERATURE CORRECTION & STRUCTURE
!    TCORR                   Temperature correction via flux conservation
!    STATEQ                  Statistical equilibrium (H NLTE rates)
!    RADIAP                  Radiative acceleration iteration
!    ROSS                    Rosseland mean opacity computation
!    CONVEC                  Convective flux by mixing-length theory
!    COMPUTE_PTURB           Turbulent pressure computation
!    VTURB_VARYDEPTH         Standard microturbulent velocity profile
!    HIGH                    High-temperature opacity flag setting
!    TTAUP                   Starting model: T-tau-P relation
!
!  EQUATION OF STATE
!    NELECT                  Electron density by charge conservation
!    COMPUTE_ONE_POP         Single-species population (Saha-Boltzmann)
!    COMPUTE_ALL_POPS        Populations for all species
!    PFSAHA                  Partition functions and Saha ionization
!    PFGROUND                Ground-state partition functions
!    PFIRON                  Iron-group partition function interpolation
!    MOLEC                   Single-molecule equilibrium
!    NMOLEC                  Full molecular equilibrium solution
!    READMOL                 Read molecular equilibrium input data
!    EQUILH2                 H2 equilibrium constant
!    PARTFNH2                H2 molecular partition function
!
!  RADIATIVE TRANSFER
!    JOSH                    Source function solver (Feautrier method)
!    ENERGY_DENSITY          Radiation energy density integration
!    BLOCKJ                  J-coefficient matrix (loaded from file)
!    BLOCKH                  H-coefficient matrix (loaded from file)
!
!  CONTINUOUS OPACITY  — Hydrogen
!    HOP                     H I bound-free + free-free
!    XKARZAS                 Karzas-Latter bound-free cross sections
!    COULX                   Hydrogenic photoionization cross section
!    COULBF1S                Ground-state Coulomb bound-free
!    COULFF                  Coulomb free-free gaunt factor
!    HMINOP                  H- bound-free + free-free
!    H2PLOP                  H2+ opacity
!    H2RAOP                  H2 Rayleigh scattering
!    H2COLLOP                H2 collision-induced absorption
!    HRAYOP                  H I Rayleigh scattering
!
!  CONTINUOUS OPACITY  — Helium
!    HE1OP                   He I bound-free (with satellites)
!    CROSSHE                 He I photoionization cross section
!    HE111S, HE12S1S, ...    He I level-specific cross sections
!    HE2OP                   He II bound-free + free-free
!    HEMIOP                  He- free-free
!    HERAOP                  He I Rayleigh scattering
!
!  CONTINUOUS OPACITY  — Metals
!    C1OP, C2OP              Carbon  I–II
!    N1OP                    Nitrogen I
!    O1OP                    Oxygen   I 
!    MG1OP, MG2OP            Magnesium I–II
!    AL1OP                   Aluminum I
!    SI1OP, SI2OP            Silicon  I–II
!    CA2OP                   Calcium  II
!    FE1OP                   Iron I
!    CHOP                    CH molecular opacity
!    OHOP                    OH molecular opacity
!
!  CONTINUOUS OPACITY  — Assembly
!    KAPP                    Total continuous opacity from all sources
!    KAPCONT                 Continuous opacity tabulation
!    COOLOP                  Cool-star opacity assembly (T < 7000 K)
!    WARMOP                  Lukewarm opacity (7000–12000 K)
!    HOTOP                   Hot-star opacity (T > 12000 K)
!    ELECOP                  Electron scattering
!    XCONOP                  Tabulated continuum opacity lookup
!    ROSSTAB                 Rosseland opacity table interpolation
!
!  LINE OPACITY
!    SELECTLINES             Line selection for current frequency
!    LINOP1                  Line opacity computation (single frequency)
!    XLINOP                  Extended line opacity computation
!    IONPOTS                 Load ionization potentials
!    ISOTOPES                Load isotope data
!
!  LINE PROFILES
!    HLINOP                  Hydrogen line opacity (Balmer, Lyman, ...)
!    HPROF4                  Hydrogen Stark + resonance profile
!    STARK                   Stark broadening half-width
!    VOIGT                   Voigt function H(a,v)
!    TABVOIGT                Voigt profile tabulation
!    HFNM                    Hydrogen oscillator strength f(n,m)
!    VCSE1F                  VCS electric field distribution
!    STARK_PROFILE           Stark profile wing function (formerly SOFBET)
!    (FASTE1 removed — replaced by direct EXPI calls)
!    EXPI                    Exponential integral E_n(x)
!
!  NUMERICAL UTILITIES
!    DERIV                   Differentiation by parabolic interpolation
!    INTEG                   Integration using parabolic coefficients
!    PARCOE                  Parabolic interpolation coefficients
!    MAP1                    Linear interpolation/remapping
!    MAP4                    Cubic interpolation/remapping
!    SOLVIT                  LU decomposition linear equation solver
!    LINTER                  Linear interpolation utility
!
!  I/O & PARSING
!    READIN                  Read and parse all input control cards
!    FREEFF                  Free-format floating point reader
!    FREEFR                  Free-format real number reader
!    NEXTWORD                 Free-format word reader (returns character string)
!    DUMP_ARRAY              Debug array print utility
!
!=========================================================================

!=========================================================================
! mod_parameters: Compile-time constants for ATLAS12
!=========================================================================

module mod_parameters

  implicit none
  integer, parameter :: kw = 80        ! Number of atmospheric depth points
  integer, parameter :: mion = 1006    ! Number of ion species
  integer, parameter :: maxmol = 200   ! Maximum number of molecules
  integer, parameter :: max1 = maxmol + 1
  integer, parameter :: maxeq = 35     ! Maximum number of equilibrium equations
  integer, parameter :: maxloc = 3 * maxmol

  ! Debug flag: set to 1 to print subroutine entry tracing to unit 6
  integer, parameter :: IDEBUG = 0

end module mod_parameters

!=========================================================================
! mod_constants: Physical constants and derived quantities
!
!   Fundamental constants: CODATA 2018 / 2019 SI redefinition.
!   k_B and h are exact by definition; c is exact by the definition
!   of the metre; sigma_T and m_u carry experimental uncertainties
!   negligible for our purposes.
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
!     CODATA 2018 (Tiesinga et al. 2021, Rev. Mod. Phys. 93, 025010).
!
!   Previously used (CODATA 1963):
!     k_B = 1.38054E-16,  h = 6.6256E-27,  m_u = 1.660E-24
!=========================================================================

module mod_constants

  implicit none

  ! =====================================================================
  !  FUNDAMENTAL CONSTANTS  (CGS-Gaussian units)
  !
  !  Source: CODATA 2018 (Tiesinga et al. 2021, RMP 93, 025010)
  !  Values marked (exact) are exact under the 2019 SI redefinition.
  ! =====================================================================

  real*8, parameter :: PI        = 3.14159265358979323846D0
  real*8, parameter :: FOURPI    = 4.0D0 * PI               ! 12.566370614...
  real*8, parameter :: SQRTPI    = 1.7724538509055159D0      ! sqrt(pi)
  real*8, parameter :: INVSQRTPI = 0.5641895835477563D0      ! 1/sqrt(pi)

  ! Boltzmann constant [erg K⁻¹]  (exact)
  real*8, parameter :: KBOL = 1.380649D-16

  ! Planck constant [erg s]  (exact)
  real*8, parameter :: HPLANCK = 6.62607015D-27

  ! Speed of light [cm s⁻¹]  (exact)
  real*8, parameter :: CLIGHT = 2.99792458D10

  ! Speed of light in other wavelength·frequency unit systems
  real*8, parameter :: CLIGHT_NM  = 2.99792458D17   ! [nm Hz]   freq = CLIGHT_NM / lambda_nm
  real*8, parameter :: CLIGHT_ANG = 2.99792458D18   ! [Å Hz]    freq = CLIGHT_ANG / lambda_Ang
  real*8, parameter :: CLIGHT_KMS = 2.99792458D5    ! [km s⁻¹]

  ! Unified atomic mass unit [g]
  ! CODATA 2018: 1.66053906660(50) × 10⁻²⁴ g
  real*8, parameter :: AMU = 1.66053907D-24

  ! Electron mass [g]
  ! CODATA 2018: 9.1093837015(28) × 10⁻²⁸ g
  real*8, parameter :: EMASS = 9.1093837D-28

  ! Elementary charge [esu = statcoulomb]
  ! CODATA 2018: 4.80320451(10) × 10⁻¹⁰ esu
  real*8, parameter :: ECHARGE = 4.80320451D-10

  ! Stefan-Boltzmann constant [erg cm⁻² s⁻¹ K⁻⁴]
  ! σ = 2π⁵k⁴/(15h³c²)
  ! CODATA 2018: 5.670374419 × 10⁻⁵
  real*8, parameter :: SIGMA_SB = 5.670374419D-5

  ! Thomson scattering cross section [cm²]
  ! σ_T = 8π/3 × (e²/m_e c²)²
  ! CODATA 2018: 6.6524587321(60) × 10⁻²⁵ cm²
  real*8, parameter :: SIGMA_THOMSON = 6.6524587D-25

  ! Electron volt [erg]
  ! CODATA 2018: 1.602176634 × 10⁻¹² erg  (exact)
  real*8, parameter :: EV_ERG = 1.602176634D-12

  ! =====================================================================
  !  DERIVED CONSTANTS
  !
  !  Pre-computed combinations that appear frequently in the code.
  !  All derived from the fundamental constants above.
  ! =====================================================================

  ! hc/k  [cm K]  — exponent factor: exp(-E_cm * HCK / T)
  real*8, parameter :: HCK = HPLANCK * CLIGHT / KBOL
  !                        = 1.43877688... (vs old 1.43879 from CODATA 1963)

  ! h/k  [s K]  — for computing HKT = h/(kT) = HOVERK / T(K)
  ! then hν/kT = HKT * freq
  real*8, parameter :: HOVERK = HPLANCK / KBOL

  ! k/eV  [K⁻¹]  — so TKEV = T * KBOL_EV
  real*8, parameter :: KBOL_EV = KBOL / EV_ERG
  !                            = 8.617333262... × 10⁻⁵ eV/K

  ! Reciprocal temperature coefficient: θ = THETA_COEFF / T
  ! where θ is the conventional log₁₀ Boltzmann parameter
  ! THETA_COEFF = log₁₀(e) / k_B(eV) = 5039.81...
  real*8, parameter :: THETA_COEFF = 0.4342944819032518D0 / KBOL_EV

  ! 4π/c  [cm⁻¹ s]  — radiation energy density factor
  real*8, parameter :: FOURPI_OVER_C = FOURPI / CLIGHT

  ! 4σ/π  [erg cm⁻² s⁻¹ K⁻⁴]  — appears in convection
  real*8, parameter :: FOUR_SIGMA_OVER_PI = 4.0D0 * SIGMA_SB / PI

  ! Planck function prefactor: 2h/c² × (10¹⁵)³  [cgs]
  ! For B_ν = BNU_PREFAC * (ν/10¹⁵)³ × exp(-hν/kT) / (1 - exp(-hν/kT))
  real*8, parameter :: BNU_PREFAC = 2.0D0 * HPLANCK / CLIGHT**2 * 1.0D45
  !                               = 1.47449...E-2  (vs old 1.47439E-2)

  ! Saha equation prefactor: 2 × (2π m_e k / h²)^(3/2)  [cm⁻³ K⁻³/²]
  ! Usage: n_e × n⁺/n₀ = SAHA_PREFAC × T^(3/2) × (U⁺/U₀) × exp(-χ/kT)
  real*8, parameter :: SAHA_PREFAC = 2.0D0 &
    * (2.0D0 * PI * EMASS * KBOL / HPLANCK**2)**1.5D0
  !                               = 2.4150...E15  (vs old 2.4148E15)

  ! Free-free opacity coefficient  [cgs]
  ! κ_ff = COEFF_FF / ν³ × g_ff × n_e × n_ion / (ρ √T)
  ! COEFF_FF = 32π² e⁶ / (3√3 m_e c h) × √(2π / (3 m_e k))
  !          = 3.6919...E8
  ! The full derivation is given by, e.g., Mihalas (1978) eq. 4-64.
  real*8, parameter :: COEFF_FF = 3.69196D8

  ! =====================================================================
  !  SPECTROSCOPIC CONSTANTS
  !
  !  These are measured physical quantities, not fundamental constants.
  !  Updated values and references are given where relevant.
  ! =====================================================================

  ! Rydberg constant for infinite mass [cm⁻¹]
  ! CODATA 2018: R_∞ = 109737.31568160(21)
  ! For hydrogen: R_H = R_∞ × m_p/(m_e + m_p) = 109677.5774...
  real*8, parameter :: RYDBERG_INF = 109737.31568D0   ! R_∞ [cm⁻¹]
  real*8, parameter :: RYDBERG_H   = 109677.576D0     ! R_H [cm⁻¹]

  ! H I ionization limit [cm⁻¹]  (observed, including quantum defect)
  real*8, parameter :: ELIM_HI = 109678.764D0

  ! Hydrogen ionization frequency [Hz]:  ν_∞ = R_H × c
  real*8, parameter :: FREQ_RYDH = RYDBERG_H * CLIGHT
  !                              = 3.288055...E15  (vs old 3.28805E15)

  ! He I Rydberg [cm⁻¹]:  R_He = R_∞ × M_He/(m_e + M_He) ≈ 109722.267
  real*8, parameter :: RYDBERG_HE = 109722.267D0

  ! H⁻ electron affinity [eV]
  ! Lykke, Murray & Lineberger (1991, Phys. Rev. A 43, 6104):
  !   E_A(H⁻) = 6082.99 ± 0.15 cm⁻¹ = 0.75419(2) eV
  ! Previously used: 0.754209 eV (Hotop & Lineberger 1985)
  real*8, parameter :: HMINUS_EA = 0.75419D0   ! [eV]


  ! =====================================================================
  !  HOLTSMARK MICROFIELD DISTRIBUTION TABLE
  !
  !  Q(beta) = probability that the electric microfield F < beta * F_0
  !  at a test charge surrounded by randomly distributed perturber ions.
  !  F_0 = 2.6031 * Z_p * e * (4*pi*N_ion/3)^(2/3) is the Holtsmark
  !  normal field strength.
  !
  !  Used for the Hummer & Mihalas (1988, ApJ 331, 794) occupation
  !  probability formalism: the probability that hydrogen level n
  !  survives dissolution is w_n = Q(beta_crit), where
  !    beta_crit = BETA_COEFF_HM88 / (n^5 * N_e^(2/3))
  !
  !  Table: 150 points on a uniform log10(beta) grid from 0.01 to 50.
  !  For beta < 0.01, Q ~ 0 (level fully bound).
  !  For beta > 50,   Q ~ 1 (level fully dissolved / always survives).
  !
  !  References:
  !    Holtsmark, J. 1919, Ann. Phys. 363, 577
  !    Hummer, D.G. & Mihalas, D. 1988, ApJ 331, 794
  !    Nayfonov, A. et al. 1999, ApJ 526, 451
  ! =====================================================================
  integer, parameter :: NQ_HOLTSMARK = 150
  real*8,  parameter :: LOG_BETA_MIN  = -2.0000000000000000D+00
  real*8,  parameter :: LOG_BETA_MAX  = 1.6989700043360187D+00
  real*8,  parameter :: LOG_BETA_STEP = 2.4825302042523617D-02

  real*8,  parameter :: Q_HOLTSMARK(150) = [ &
    1.414671303101461D-07, &
    1.679306576215498D-07, &
    1.993444994252022D-07, &
    2.366346434794731D-07, &
    2.809002805410456D-07, &
    3.334461984807552D-07, &
    3.958212340900693D-07, &
    4.698639150426803D-07, &
    5.577566360795190D-07, &
    6.620899645927724D-07, &
    7.859389687582573D-07, &
    9.329538149374408D-07, &
    1.107467300593689D-06, &
    1.314622486716282D-06, &
    1.560524184270007D-06, &
    1.852418749731333D-06, &
    2.198907475766682D-06, &
    2.610199848759794D-06, &
    3.098414113872032D-06, &
    3.677933974566057D-06, &
    4.365831897218437D-06, &
    5.182371440131865D-06, &
    6.151603336163127D-06, &
    7.302072795786128D-06, &
    8.667658741262937D-06, &
    1.028856952546623D-05, &
    1.221252424026525D-05, &
    1.449615410836965D-05, &
    1.720666483126171D-05, &
    2.042380831343906D-05, &
    2.424222111026057D-05, &
    2.877419750061379D-05, &
    3.415297755655442D-05, &
    4.053664530978627D-05, &
    4.811274949652059D-05, &
    5.710377986116004D-05, &
    6.777365615449882D-05, &
    8.043541539930427D-05, &
    9.546031643883215D-05, &
    1.132886200658068D-04, &
    1.344423491071419D-04, &
    1.595403868046069D-04, &
    1.893163349213092D-04, &
    2.246396266107989D-04, &
    2.665404747620448D-04, &
    3.162393359894821D-04, &
    3.751816855286531D-04, &
    4.450790309992992D-04, &
    5.279572453545271D-04, &
    6.262134733849464D-04, &
    7.426830638017916D-04, &
    8.807182017925553D-04, &
    1.044280166083684D-03, &
    1.238047410111971D-03, &
    1.467541967671182D-03, &
    1.739277006102283D-03, &
    2.060928688577077D-03, &
    2.441535851096500D-03, &
    2.891731333794821D-03, &
    3.424009106941961D-03, &
    4.053031566870220D-03, &
    4.795981500493777D-03, &
    5.672963167595142D-03, &
    6.707456645870922D-03, &
    7.926828918285955D-03, &
    9.362904019250405D-03, &
    1.105259272456321D-02, &
    1.303857956034104D-02, &
    1.537006106786488D-02, &
    1.810352400480589D-02, &
    2.130354516765621D-02, &
    2.504358545069891D-02, &
    2.940673929734126D-02, &
    3.448638660516894D-02, &
    4.038667732307838D-02, &
    4.722275959985061D-02, &
    5.512064100258095D-02, &
    6.421655023415057D-02, &
    7.465564600958080D-02, &
    8.658990347532045D-02, &
    1.001750012765402D-01, &
    1.155660400501043D-01, &
    1.329119530555128D-01, &
    1.523485300474136D-01, &
    1.739900743830450D-01, &
    1.979198568471881D-01, &
    2.241797192777529D-01, &
    2.527594102861463D-01, &
    2.835864859104609D-01, &
    3.165178467055500D-01, &
    3.513341604156768D-01, &
    3.877384740815520D-01, &
    4.253601858765986D-01, &
    4.637651729352414D-01, &
    5.024722401913855D-01, &
    5.409752086892451D-01, &
    5.787690154213567D-01, &
    6.153773349921458D-01, &
    6.503786770405923D-01, &
    6.834278586173648D-01, &
    7.142702953716408D-01, &
    7.427476403869326D-01, &
    7.687946986070335D-01, &
    7.924289205182540D-01, &
    8.137347922132799D-01, &
    8.328458692792989D-01, &
    8.499270198328847D-01, &
    8.651588016097912D-01, &
    8.787250560031417D-01, &
    8.908040093261350D-01, &
    9.015625910920573D-01, &
    9.111533622086631D-01, &
    9.197133565854327D-01, &
    9.273641961484557D-01, &
    9.342129682962188D-01, &
    9.403534976718694D-01, &
    9.458677646981941D-01, &
    9.508273225009358D-01, &
    9.552946289851881D-01, &
    9.593242514728686D-01, &
    9.629639386838167D-01, &
    9.662555465715050D-01, &
    9.692358425479507D-01, &
    9.719372021216177D-01, &
    9.743882111568676D-01, &
    9.766141278943968D-01, &
    9.786373698849190D-01, &
    9.804778420731699D-01, &
    9.821532811516551D-01, &
    9.836795130160503D-01, &
    9.850706740272015D-01, &
    9.863394492399232D-01, &
    9.874972091274948D-01, &
    9.885541671623754D-01, &
    9.895195463975206D-01, &
    9.904015671946884D-01, &
    9.912077851485863D-01, &
    9.919449707318418D-01, &
    9.926192297837247D-01, &
    9.932361028219120D-01, &
    9.938006248722731D-01, &
    9.943173688128598D-01, &
    9.947904964619882D-01, &
    9.952237635494758D-01, &
    9.956206036156809D-01, &
    9.959840685797559D-01, &
    9.963172215001466D-01, &
    9.966224531295311D-01, &
    9.969021568157113D-01, &
    9.971585495409699D-01  ]

  ! beta_crit(n, N_e) = BETA_COEFF_HM88 / (n^5 * N_e^(2/3))
  ! Derived from Inglis-Teller critical field F_crit = Ry/(3*n^5*e*a_0)
  ! and Holtsmark normal field F_0 = 2.6031*e*(4*pi*N_e/3)^(2/3)
  real*8,  parameter :: BETA_COEFF_HM88 = 8.798905203085208D+14

end module mod_constants

!=========================================================================
! mod_atlas_data: Shared arrays replacing all COMMON blocks
!   Original ATLAS12 used 58 COMMON blocks for shared state.
!   All are consolidated here as module variables with SAVE attribute.
!=========================================================================

module mod_atlas_data

  use mod_parameters
  use mod_constants
  implicit none
  save

  ! --- Rosseland optical depth scale ---
  ! COMMON /ABROSS/
  real*8  :: ABROSS(kw), TAUROS(kw)

  ! --- Total absorption and scattering ---
  ! COMMON /ABTOT/
  real*8  :: ABTOT(kw), ALPHA(kw)

  ! --- Cool-star continuous opacities (C I, Mg I, Al I, Si I, Fe I) ---
  ! COMMON /ACOOL/
  real*8  :: AC1(kw), AMG1(kw), AAL1(kw), ASI1(kw), AFE1(kw)
  ! Per-species source functions: computed by C1OP/AL1OP/FE1OP/HE2OP
  ! but NOT consumed by KAPP, which follows the atlas12.for convention
  ! of weighting all metal/He continua by BNU.  Kept for potential
  ! future per-species source function treatment (e.g. NLTE departure
  ! coefficient hooks).
  real*8  :: SAL1(kw), SFE1(kw), SHE2(kw), SC1(kw)

  ! --- Lukewarm opacities (C II, N I, O I, Mg II, Si II, Ca II) ---
  ! COMMON /ALUKE/
  real*8  :: AC2(kw), AN1(kw), AO1(kw), AMG2(kw), ASI2(kw), ACA2(kw)

  ! --- Convection parameters ---
  ! COMMON /CONV/
  real*8  :: DLTDLP(kw) = 0.0d0, HEATCP(kw) = 0.0d0, DLRDLT(kw) = 0.0d0, VELSND(kw) = 0.0d0
  real*8  :: GRDADB(kw) = 0.0d0, HSCALE(kw) = 0.0d0, FLXCNV(kw) = 0.0d0, VCONV(kw) = 0.0d0
  real*8  :: MIXLTH = 2.0d0, OVERWT = 0.0d0
  real*8  :: FLXCNV0(kw), FLXCNV1(kw)
  integer :: IFCONV = 1, NCONV = 30

  ! --- NLTE departure coefficients ---
  ! COMMON /DEPART/
  real*8  :: BHYD(kw, 6) = 1.0d0, BMIN(kw) = 1.0d0
  integer :: NLTEON = 0

  ! --- Electron density ---
  ! COMMON /EDENS/
  real*8  :: EDENS(kw)
  integer :: IFEDNS

  ! --- Element abundances, atomic masses, labels ---
  ! COMMON /ELEM/
  real*8       :: YABUND(99)
  ! Default solar abundances: Anders & Grevesse (1989, Geochim. Cosmochim.
  ! Acta, 53, 197).  H and He are fractional number densities (H+He=1);
  ! Z >= 3 are log10(N_Z / N_total).  Elements with -20.00 have no
  ! astrophysically relevant abundance.
  real*8       :: ABUND(99) = (/ &
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
  real*8       :: ATMASS(99) = (/ 1.008D0, 4.003D0, &
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
  character(2) :: ELEM(99) = (/ 'H ','He', &
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
  real*8  :: FLUX = 0.0d0, FLXERR(kw) = 0.0d0, FLXDRV(kw) = 0.0d0, FLXRAD(kw)

  ! --- Free-format input parser state ---
  ! COMMON /FREE/
  integer :: NUMCOL, LETCOL, LAST, MORE, IFFAIL, MAXPOW

  ! --- Current frequency point ---
  ! COMMON /FREQ/
  real*8  :: FREQ, FREQLG, EHVKT(kw), STIM(kw), BNU(kw), WAVE, WAVENO

  ! --- Frequency set for opacity distribution functions ---
  ! COMMON /FRESET/
  real*8  :: WAVESET(30000), RCOSET(30000)
  integer :: NULO, NUHI, NUMNU = 0, NUSTEP

  ! --- Hydrogen bound-free opacity tables ---
  ! COMMON /H1TAB/
  real*8  :: H0TAB(2001), H1TAB(2001), H2TAB(2001)

  ! --- Geometric height scale ---
  ! COMMON /HEIGHT/
  real*8  :: HEIGHT(kw)

  ! --- Control flags ---
  ! COMMON /IF/
  integer :: IFCORR = 1, IFPRES = 1, IFSURF = 0, IFSCAT = 1, IFMOL = 1, IFREADLINES = 1

  ! ROSSTAB interpolation mode:
  !   1 = original bilinear (4-quadrant nearest neighbor)
  !   2 = Shepard (K-nearest, inverse-distance weighted, p=3)
  integer :: IROSSTAB = 1

  ! --- Molecular equilibrium ---
  ! COMMON /IFEQUA/
  real*8  :: XNMOLCODE(maxmol), EQUIL(6, maxmol)
  integer :: IFEQUA(101), KCOMPS(maxloc), LOCJ(max1), IDEQUA(maxeq)
  integer :: NEQUA, NEQUA1, NEQNEQ, NUMMOL, NLOC

  ! --- Opacity control flags ---
  ! COMMON /IFOP/
  integer :: IFOP(20) = (/1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0/)

  ! --- Population control flag ---
  ! COMMON /IFPOP/
  integer :: IFPOP

  ! --- Line data integer header words ---
  ! COMMON /IIIIIII/ (renamed LINEREC) — packed binary line-data record
  integer    :: IWL
  integer*2  :: IELION, IELO, IGFLOG, IGR, IGS, IGW

  ! --- Isotope fractional abundances ---
  ! COMMON /ISOTOPE/
  real*8  :: ISOTOPE(10, 2, mion)

  ! --- Iteration control ---
  ! COMMON /ITER/
  integer :: ITER, ifprnt(60) = 2, ifpnch(60) = 0, NUMITS = 0

  ! --- Title and header data ---
  ! COMMON /JUNK/
  character(1) :: TITLE(74) = ' '
  real*8  :: XSCALE = 1.0d0
  character(4) :: WLTE = 'LTE '
  character(256) :: DATADIR    ! Path to data files (from $ATLAS12)
  integer :: INPUTDATA

  ! --- Flag set if SYNTHE is running ---
  integer :: IFSYNTHE = 0

  ! --- Radiative transfer J-coefficient matrix ---
  ! COMMON /MATXJ/
  real*8  :: COEFJ(51, 51)

  ! --- Radiative transfer H-coefficient matrix ---
  ! COMMON /MATXH/
  real*8  :: COEFH(51, 51)

  ! --- Angular quadrature ---
  ! COMMON /MUS/
  real*8  :: ANGLE(20), SURFI(20) = 0.0d0
  integer :: NMU = 1

  ! --- Individual continuous opacity sources ---
  ! COMMON /OPS/
  real*8  :: AHYD(kw), AH2P(kw), AHMIN(kw), SIGH(kw)
  real*8  :: AHE1(kw), AHE2(kw), AHEMIN(kw), SIGHE(kw)
  real*8  :: ACOOL(kw), ALUKE(kw), AHOT(kw)
  real*8  :: SIGEL(kw), SIGH2(kw), AHLINE(kw), ALINES(kw), SIGLIN(kw)
  real*8  :: AXLINE(kw), SIGXL(kw), AXCONT(kw), SIGX(kw)
  real*8  :: SHYD(kw), SHMIN(kw), SHLINE(kw), SXLINE(kw), SXCONT(kw)

  ! --- Total opacity ---
  ! COMMON /OPTOT/
  real*8  :: ACONT(kw), SCONT(kw), ALINE(kw), SLINE(kw)
  real*8  :: SIGMAC(kw), SIGMAL(kw)

  ! --- Ionization potentials ---
  ! COMMON /POTION/
  real*8  :: POTION(999), POTIONSUM(999)

  ! --- Total pressure ---
  ! COMMON /PTOTAL/
  real*8  :: PTOTAL(kw)

  ! --- Output control ---
  ! COMMON /PUT/
  real*8  :: PUT
  integer :: IPUT

  ! --- Zero-point pressure and radiation field ---
  ! COMMON /PZERO/
  real*8  :: PZERO, PCON, PRADK0, PTURB0
  real*8  :: KNU(kw), PRADK(kw), RADEN(kw)

  ! --- Radiation pressure ---
  ! COMMON /RAD/
  real*8  :: ACCRAD(kw) = 0.0d0, PRAD(kw) = 0.0d0

  ! --- Column mass depth scale ---
  ! COMMON /RHOX/
  real*8  :: RHOX(kw)
  integer :: NRHOX = 0

  ! --- Mean intensity diagnostic ---
  ! COMMON /RR/
  real*8  :: RJMINSNU(kw), RDIAGJNU(kw)

  ! --- Thermodynamic state (pressure, electron density, etc.) ---
  ! COMMON /STATE/
  real*8  :: P(kw), XNE(kw), XNATOM(kw), RHO(kw), CHARGESQ(kw)

  ! --- Depth integration step parameters ---
  ! COMMON /STEPLG/
  real*8  :: STEPLG = 0.125d0, TAU1LG = -6.875d0
  integer :: KRHOX = 0

  ! --- Continuous opacity table ---
  ! COMMON /TABCONT/
  real*8  :: TABCONT(kw, 344)
  real*8  :: WAVETAB(344)
  integer :: IWAVETAB(344)

  ! --- Log lookup table (replaced array with function) ---
  ! COMMON /TABLOG/ removed — now computed inline via TABLOG()

  ! --- Optical depth, source function, flux moments ---
  ! COMMON /TAUSHJ/
  real*8  :: TAUNU(kw), SNU(kw), HNU(kw), JNU(kw), JMINS(kw)

  ! --- Standard optical depth scale ---
  ! COMMON /TAUSTD/
  real*8  :: TAUSTD(kw)

  ! --- Effective temperature and gravity ---
  ! COMMON /TEFF/
  real*8  :: TEFF = 0.0d0, GRAV = 0.0d0, GLOG

  ! --- Temperature and derived quantities ---
  ! COMMON /TEMP/
  real*8  :: T(kw), TKEV(kw), TK(kw), HKT(kw), HCKT(kw), TLOG(kw)
  integer :: ITEMP

  ! --- Temperature smoothing ---
  ! COMMON /TSMOOTH/
  integer :: J1SMOOTH = 0, J2SMOOTH = 0
  real*8  :: WTJM1 = 0.3d0, WTJ = 0.4d0, WTJP1 = 0.3d0, TSMOOTH(kw)

  ! --- Turbulent pressure ---
  ! COMMON /TURBPR/
  real*8  :: VTURB(kw) = 0.0d0, PTURB(kw) = 0.0d0, TRBFDG = 0.0d0, TRBCON = 0.0d0, TRBPOW = 0.0d0, TRBSND = 0.0d0
  integer :: IFTURB = 0

  ! --- Wavelength grid control ---
  ! COMMON /WAVEY/
  real*8  :: WBEGIN, DELTAW
  integer :: IFWAVE = 0

  ! --- Current line data ---
  ! COMMON /WWWWWWW/ (renamed LINEPARAM) — packed line-parameter record
  real*8  :: WLVAC
  integer :: NELION
  real*8  :: CGF, ELO, GAMMAR, GAMMAS, GAMMAW

  ! --- Depth-dependent abundances ---
  ! COMMON /XABUND/
  real*8  :: XABUND(kw, 99), WTMOLE(kw), XRELATIVE(99) = 0.D0

  ! --- Isotope fractions and masses ---
  ! COMMON /XISO/
  real*8  :: XISO(10, mion), AMASSISO(10, mion)

  ! --- Line opacity distribution ---
  ! COMMON /XLINES/
  real*8  :: XLINES(kw, 30000)

  ! --- Ion number densities ---
  ! COMMON /XNF/
  real*8  :: XNF(kw, mion), XNFP(kw, mion), XNH2(kw)

  ! --- Doppler-broadened line opacity ---
  ! COMMON /XNFDOP/
  real*8  :: XNFDOP(kw, mion), DOPPLE(kw, mion)

  ! --- Molecular number densities ---
  ! COMMON /XNMOL/
  real*8  :: XNMOL(kw, maxmol), XNFPMOL(kw, maxmol)

  ! Flag: set to 1 after MOLEC has read molecular data from INPUTDATA
  ! (or when NMOLEC has populated the arrays directly, skipping the read).
  integer :: MOLEC_IREAD = 0

  ! --- Saved number densities for equilibrium ---
  ! COMMON /XNSAVE/
  real*8  :: XNSAVE(kw, maxeq)

  ! --- In-memory line data storage (replaces fort.12 I/O) ---
  integer*4, allocatable :: LINEDATA(:,:)   ! (4, NLINES_STORED)
  integer :: NLINES_STORED = 0

  ! --- Stehlé MMM hydrogen Stark broadening tables ---
  ! Preprocessed from Stehlé & Hutcheon (1999) and Stehlé & Fouquet (2010).
  ! Loaded by INIT_STARK_TABLES; used by STARK_MMM.
  integer, parameter :: NSTARK_SERIES = 4         ! Ly, Ba, Pa, Br
  integer, parameter :: NSTARK_DALPHA = 60        ! Δα grid points
  integer, parameter :: NSTARK_TEMPS  = 10        ! temperature grid points
  integer, parameter :: NSTARK_DENS_MAX = 20      ! max density grid points

  type :: stark_series_t
    integer :: n_lower                             ! lower quantum number
    integer :: n_upper_min, n_upper_max            ! upper quantum number range
    integer :: n_transitions                       ! = n_upper_max - n_upper_min + 1
    integer :: n_dens                              ! actual number of density points
    real*8  :: density_grid(NSTARK_DENS_MAX)       ! Ne [cm^-3]
    real*8  :: temp_grid(NSTARK_TEMPS)             ! T [K]
    real*8  :: log_dalpha_grid(NSTARK_DALPHA)      ! log10(Δα)
    integer, allocatable :: max_dens_idx(:)        ! (n_transitions) highest valid density
    real*8,  allocatable :: k_alpha(:)             ! (n_transitions) asymptotic wing constant
    real*8,  allocatable :: profiles(:,:,:,:)      ! (NSTARK_DALPHA, NSTARK_TEMPS, n_dens, n_transitions)
    logical :: loaded = .false.
  end type stark_series_t

  type(stark_series_t), target :: STEHLE_DATA(NSTARK_SERIES)
  logical :: STEHLE_TABLES_LOADED = .false.

contains

  ! Unpack LINEREC(4) record into module line-data variables
  ! Replaces the removed EQUIVALENCE (LINEREC(1),IWL)
  subroutine UNPACK_LINEDATA(III)
    integer*4, intent(in) :: III(4)
    integer*2 :: pair(2)
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
  end subroutine UNPACK_LINEDATA

  ! Pack module line-data variables back into LINEREC(4) record
  ! For use before WRITE(12) LINEREC when variables have been modified
  subroutine PACK_LINEDATA(III)
    integer*4, intent(out) :: III(4)
    integer*2 :: pair(2)
    III(1) = IWL
    pair(1) = IELION; pair(2) = IELO
    III(2) = TRANSFER(pair, III(2))
    pair(1) = IGFLOG; pair(2) = IGR
    III(3) = TRANSFER(pair, III(3))
    pair(1) = IGS; pair(2) = IGW
    III(4) = TRANSFER(pair, III(4))
  end subroutine PACK_LINEDATA

  ! Pack variables into LINEPARAM(4) matching F77 COMMON /LINEPARAM/ layout:
  ! LINEPARAM(1) = WLVAC (REAL*8)
  ! LINEPARAM(2) = [NELION (INTEGER*4), CGF (REAL*4)]
  ! LINEPARAM(3) = [ELO (REAL*4), GAMMAR (REAL*4)]
  ! LINEPARAM(4) = [GAMMAS (REAL*4), GAMMAW (REAL*4)]
  subroutine PACK_LINEPARAM(W, WLVAC_in, NELION_in, CGF_in, ELO_in, GAMMAR_in, GAMMAS_in, GAMMAW_in)
    real*8, intent(out) :: W(4)
    real*8, intent(in) :: WLVAC_in
    integer, intent(in) :: NELION_in
    real*4, intent(in) :: CGF_in, ELO_in, GAMMAR_in, GAMMAS_in, GAMMAW_in
    real*4 :: pair(2)
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
  end subroutine PACK_LINEPARAM

  ! Unpack LINEPARAM(4) into variables matching F77 COMMON /LINEPARAM/ layout
  subroutine UNPACK_LINEPARAM(W, WLVAC_out, NELION_out, CGF_out, ELO_out, GAMMAR_out, GAMMAS_out, GAMMAW_out)
    real*8, intent(in) :: W(4)
    real*8, intent(out) :: WLVAC_out
    integer, intent(out) :: NELION_out
    real*4, intent(out) :: CGF_out, ELO_out, GAMMAR_out, GAMMAS_out, GAMMAW_out
    real*4 :: pair(2)
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
  end subroutine UNPACK_LINEPARAM

  elemental real*8 function TABLOG(I)
    integer*2, intent(in) :: I
    TABLOG = 10.D0**((I-16384)*.001D0)
  end function TABLOG


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

  implicit none

  integer, intent(in) :: MODE

  ! --- Local variables ---
  real*8  :: HSURF, TAUEND
  real*8  :: HLAM, HLAMLG, HLAMMG, HNULG, HNUMG
  real*8  :: FLXCNVRATIO(kw)
  real*8  :: SURFIN(20)       ! Note: SURFIN accumulation was disabled in original
                               ! code (opacity-sampling replaced frequency groups).
                               ! SURFI from JOSH holds the actual surface intensity.
  integer :: J, I, IZ, MU, JTAU1

  ! Persistent state across calls (NU and IFHEAD survive between MODE=1 and MODE=2-4)
  integer, save :: NU = 0
  integer, save :: IFHEAD = 0

  ! =================================================================

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING PUTOUT'
  select case (MODE)

  ! =================================================================
  ! MODE 1 — Write column headings to output and punch files
  ! =================================================================
  case (1)

    if (IFPRNT(ITER) == 0) return
    IFHEAD = 0
    NU = NULO - NUSTEP

    ! Write header to surface flux file (unit 8)
    if (IFPNCH(ITER) < 2) return
    write(8, '("TEFF ",F7.0,"  GRAVITY",F8.4,1X,A4 / "TITLE ",74A1)') &
      TEFF, GLOG, WLTE, TITLE

    ! For intensity mode (IFSURF=2), also write angle grid
    if (IFSURF /= 2) return
    write(8, '(I3," ANGLES",10F7.4 / 10X,10F7.4)') NMU, (ANGLE(MU), MU=1,NMU)

  ! =================================================================
  ! MODE 2,3,4 — Print flux/intensity at current wavelength point
  ! =================================================================
  case (2, 3, 4)

    NU = NU + NUSTEP
    HSURF = HNU(1)
    TAUEND = log10(TAUNU(NRHOX))
    if (HSURF <= 0.0D0) HSURF = 1.0D-50

    ! --- Detailed frequency-by-frequency output (print level >= 2) ---
    if (IFPRNT(ITER) >= 2.AND.IDEBUG == 1) then

      ! Find first depth where optical depth exceeds unity
      JTAU1 = NRHOX
      do J = 1, NRHOX
        if (TAUNU(J) > 1.0D0) then
          JTAU1 = J
          exit
        end if
      end do
      TAUEND = log10(TAUNU(NRHOX))

      ! --- Flux output format (IFSURF = 0 or 1) ---
      if (IFSURF == 0 .or. IFSURF == 1) then
        if (IFHEAD == 0) then
          write(6, 101)
101       format('1', ///// 10X, 'WAVE', 7X, 'HLAMBDA', 7X, 'LOG H', 7X, 'MAG', &
            10X, 'FREQUENCY', 8X, 'HNU', 10X, 'LOG H', 7X, 'MAG', 2X, 'TAUONE', &
            ' TAUNU')
        end if
        IFHEAD = 1

        ! Convert H_nu to H_lambda and compute magnitudes
        HLAM   = HSURF * FREQ / WAVE
        HNULG  = log10(HSURF)
        HLAMLG = log10(HLAM)
        HLAMMG = -2.5D0 * HLAMLG
        HNUMG  = -2.5D0 * HNULG

        write(6, 401) NU, WAVE, HLAM, HLAMLG, HLAMMG, &
                       FREQ, HSURF, HNULG, HNUMG, JTAU1, TAUEND, NU
401     format(I6, F10.3, 1PE13.4, 0PF12.5, F10.3, 1PE20.6, E13.4, &
          0PF12.5, F10.3, I6, F6.2, I6)

        ! Write tau(nu) profile to diagnostic file (unit 50)
        write(50, "(I6,F12.3,(100F10.4))") NU, WAVE, &
          (log10(TAUNU(J)), J=1,NRHOX)
      end if

      ! --- Intensity output format (IFSURF = 2) ---
      if (IFSURF == 2) then
        if (IFHEAD == 0) then
          write(6, 102)
102       format('1', ///// 10X, 'WAVE', 5X, 'FREQUENCY', 3X, 'TAUONE TAUNU', &
            5('   MU  INTENSITY '))
        end if
        IFHEAD = 1

        write(6, 406) NU, WAVE, FREQ, JTAU1, TAUEND, &
                       (ANGLE(MU), SURFIN(MU), MU=1,NMU)
406     format(I6, F9.3, 1PE15.6, I6, 0PF6.2, 5(0PF7.4, 1PE10.3) / &
          (42X, 5(0PF7.4, 1PE10.3)))
      end if

    end if  ! IFPRNT >= 2

    ! --- Write to punch file (unit 8) ---
    if (IFPNCH(ITER) < 2) return
    if (IFSURF > 2) return

    if (IFSURF == 2) then
      ! Intensity mode: write surface intensity at all angles
      write(8, 416) NU, WAVE, FREQ, (SURFIN(MU), MU=1,NMU)
416   format('INTENSITY', I5, F9.2, 1PE15.6 / (1P8E10.3))
      if (NU == NUHI) write(8, 416)
    else
      ! Flux mode: write wavelength (in nm) and Eddington flux H_nu
      write(8, '(F13.2,E13.4)') WAVE * 10.0D0, HSURF
    end if

  ! =================================================================
  ! MODE 5 — End-of-iteration summary and model punch output
  ! =================================================================
  case (5)

    ! --- Print iteration summary tables to screen ---
    if (IFPRNT(ITER) /= 0) then

       if (IDEBUG == 1) then
          ! Convection and structure parameters vs. depth
          write(6, 501) (J, RHOX(J), PTOTAL(J), PTURB(J), GRDADB(J), DLTDLP(J), &
               VELSND(J), DLRDLT(J), HEATCP(J), HSCALE(J), VCONV(J), FLXCNV(J), &
               J=1,NRHOX)
501       format(///, '       RHOX       PTOTAL     PTURB      GRDADB', &
               '     DLTDLP     VELSND     DLRDLT     HEATCP     HSCALE     VCONV', &
               '     FLXCNV            ', / (I3, 1P11E11.3))
          
          ! Number densities: H I, H II, He I, He II, He III + turbulent velocity
          write(6, 503) (J, XNATOM(J), RADEN(J), PRADK(J), XNFP(J,1), XNFP(J,2), &
               XNFP(J,3), XNFP(J,4), XNFP(J,5), VTURB(J), &
               FLXCNV0(J), FLXCNV1(J), J=1,NRHOX)
503       format(///, '      XNATOM      RADEN      PRADK     XNFPH1', &
               '    XNFPH2     XNFPHE1    XNFPHE2    XNFPHE3     VTURB', &
               / (I3, 1P11E11.3))
       
       endif

      ! Compute convective flux fraction at each depth
      do J = 1, NRHOX
        if (IFCORR == 0) FLXRAD(J) = FLUX - FLXCNV(J)
        FLXCNVRATIO(J) = FLXCNV(J) / (FLXCNV(J) + FLXRAD(J))
      end do

      ! Full model atmosphere table to unit 66
      write(66, 542) (J, RHOX(J), T(J), P(J), XNE(J), RHO(J), ABROSS(J), &
        HEIGHT(J), TAUROS(J), FLXCNVRATIO(J), ACCRAD(J), FLXERR(J), FLXDRV(J), &
        J=1,NRHOX)
542   format('0', 35X, 'ELECTRON', 11X, &
        'ROSSELAND    HEIGHT   ROSSELAND   FRACTION  RADIATIVE', &
        '         PERCENT FLUX', / &
        '       RHOX      TEMP    PRESSURE    NUMBER', &
        '    DENSITY      MEAN       (KM)      DEPTH', &
        '    CONV FLUX  ACCELERATION', &
        '   ERROR   DERIV', / (I3, 1PE10.3, 0PF9.1, 1P8E11.3, 2P2E11.3))

    end if  ! IFPRNT

    ! --- Write model to punch file (unit 7) ---
    if (IFPNCH(ITER) == 0) return

    ! Model parameters and abundances
    write(7, '("TEFF ",F7.0,"  GRAVITY",F8.4,1X,A4)') TEFF, GLOG, WLTE
    write(7, '("TITLE ",74A1)') TITLE

    ! Abundance table with relative offsets
    write(7, 553) ELEM(1), ABUND(1), ELEM(2), ABUND(2), &
      (IZ, ELEM(IZ), ABUND(IZ), XRELATIVE(IZ), IZ=3,99)
553 format(' ABUNDANCE TABLE' / '    1', A2, F10.6, '       2', A2, F10.6 &
      / (5(I5, A2, F7.3, F6.3)))

    ! Model structure: column mass, T, P, N_e, kappa_Ross, g_rad, v_turb,
    !                  convective flux, convective velocity
    write(7, 554) NRHOX, (RHOX(J), T(J), P(J), XNE(J), ABROSS(J), ACCRAD(J), &
      VTURB(J), FLXCNV(J), VCONV(J), J=1,NRHOX)
554 format('READ DECK6', I3, &
      '     RHOX         T         P       XNE', &
      '     ABROSS    ACCRAD     VTURB    FLXCNV     VCONV' &
      / (13X, 1PE12.5, 0PF10.2, 1P7E10.3))

    ! Surface radiation pressure constant
    write(7, '("PRADK",1PE11.4)') PRADK0

    ! NLTE departure coefficients (if NLTE mode is on)
    if (NLTEON /= 0) then
      write(7, 556) NRHOX, (RHOX(J), (BHYD(J,I), I=1,6), BMIN(J), J=1,NRHOX)
556   format('READ DEPARTURE COEFFICIENTS', I3, ' RHOX  BHYD 1-6  BMIN' &
        / (1PE11.4, 0P7F9.4))
    end if

    ! Signal completion
    write(7, '("BEGIN",20X,"ITERATION ",I3," COMPLETED")') ITER
    close(UNIT=7)

  end select

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

  implicit none

  ! --- Arguments ---
  integer, intent(in)  :: MODE
  real*8,  intent(in)  :: RCOWT   ! quadrature weight for current frequency

  ! --- Named constants ---
  ! SIGMA_SB and FOURPI now from mod_constants
  real*8, parameter :: EULER_M1  = 0.922784335098467D0  ! 1 - gamma_Euler + ln(2)
  real*8, parameter :: RDIAGJ_FLOOR = 1.0D-30    ! prevent division by zero
  real*8, parameter :: ABROSS_FLOOR = 1.0D-30    ! prevent division by zero
  real*8, parameter :: DEL_FLOOR    = 1.0D-10    ! prevent division by zero in superadiabatic excess

  ! --- Persistent locals (accumulate across MODE 1→2→3 calls) ---
  real*8, save :: RJMINS(kw)   ! integrated opacity-weighted (J - S)
  real*8, save :: RDABH(kw)    ! integrated flux divergence term
  real*8, save :: RDIAGJ(kw)   ! integrated diagonal Lambda-operator response
  real*8, save :: OLDT1(kw)    ! previous iteration's total correction (for damping)

  ! --- MODE 2 locals ---
  real*8  :: DABTOT(kw)        ! d(kappa_tot)/d(RHOX)
  real*8  :: TERM1, TERM2, D, EX, DIAGJ, DBDT

  ! --- MODE 3 locals: derivatives and gradients ---
  real*8  :: DTDRHX(kw)        ! dT/d(RHOX)
  real*8  :: DDLT(kw)          ! d(nabla)/d(RHOX)
  real*8  :: DABROS(kw)        ! d(kappa_Ross)/d(RHOX)

  ! --- MODE 3 locals: convection ---
  real*8  :: CNVFLX(kw)        ! local copy of convective flux (smoothed)
  real*8  :: SMOOTH(kw)        ! smoothing buffer for convective flux
  real*8  :: DDEL(kw)          ! convective response factor (1 + D/(D+DEL))/DEL
  real*8  :: HRATIO(kw)        ! convective-to-total flux ratio
  real*8  :: DEL, VCO, FLUXCO, TAUB, CNVFL

  ! --- MODE 3 locals: flux correction ---
  real*8  :: CODRHX(kw)        ! integrand for G(tau) integrating factor
  real*8  :: G(kw)             ! integrating factor exp(integral of CODRHX)
  real*8  :: GFLUX(kw)         ! G * (flux error) / (flux response)
  real*8  :: DTAU(kw)          ! integrated tau correction
  real*8  :: DTFLUX(kw)        ! flux-constancy temperature correction

  ! --- MODE 3 locals: Lambda correction ---
  real*8  :: DTLAMB(kw)        ! Lambda-iteration temperature correction
  real*8  :: DTSURF(kw)        ! surface boundary correction
  real*8  :: DTCONV(kw)        ! convective flux correction (Crivellari-Simonneau)
  real*8  :: T1(kw)            ! total correction
  real*8  :: TEFF25            ! clamp: Teff/25
  real*8  :: DTSUR             ! surface correction magnitude
  real*8  :: DUM(kw), TINTEG(kw), TAV
  real*8  :: TONE_ARR(1), TTWO_ARR(1), XNEW_TMP(1)

  ! --- MODE 3 locals: convective adiabatic sweep ---
  real*8  :: T_SWEEP(kw)       ! corrected T for adiabatic integration
  real*8  :: DT_RAD            ! radiative correction at a layer
  real*8  :: DT_ADIAB          ! adiabatic sweep correction at a layer
  real*8  :: DLNP              ! Δln P between adjacent layers
  real*8  :: T_ADIAB           ! adiabatic extrapolation temperature
  real*8  :: FCONV_RATIO       ! convective flux fraction at a depth
  integer :: JANCHOR           ! shallowest convective layer (anchor)

  ! --- MODE 3 locals: RHOX correction / remapping ---
  real*8  :: TPLUS(kw), TNEW1(kw), TNEW2(kw), PRDNEW(kw)
  real*8  :: AB1(kw), PTOT1(kw), P1(kw)
  real*8  :: AB2(kw), PTOT2(kw), P2(kw)
  real*8  :: PPP(kw), RRR(kw), DRHOX(kw)
  real*8  :: REMAP(kw, 10)     ! buffer for MAP1 remapping

  ! --- Scalar temps ---
  integer :: J, I, K, IDUM, IFUDGE, JCONV
  real*8  :: ABROSS_safe

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING TCORR'

  !=====================================================================
  ! MODE 1: Zero frequency-integral accumulators
  !=====================================================================
  if (MODE == 1) then

    do J = 1, NRHOX
      RJMINS(J) = 0.0D0
      RDABH(J)  = 0.0D0
      RDIAGJ(J) = 0.0D0
      FLXRAD(J) = 0.0D0
    end do
    return

  !=====================================================================
  ! MODE 2: Accumulate frequency integrals at current wavelength
  !=====================================================================
  else if (MODE == 2) then

    ! --- Opacity gradient term: d(ln kappa)/d(RHOX) * H_nu ---
    call DERIV(RHOX, ABTOT, DABTOT, NRHOX)
    do J = 1, NRHOX
      RDABH(J)  = RDABH(J) + DABTOT(J) / ABTOT(J) * HNU(J) * RCOWT
      RJMINSNU(J) = ABTOT(J) * JMINS(J) * RCOWT
      RJMINS(J) = RJMINS(J) + RJMINSNU(J)
      FLXRAD(J) = FLXRAD(J) + HNU(J) * RCOWT
    end do

    ! --- Diagonal Lambda operator: sum kappa * (Lambda_diag - 1) * (1-eps) * dB/dT ---
    !     Lambda_diag is computed from E_3 exponential integrals of the
    !     monochromatic optical depth increments between adjacent layers.
    TERM2 = 0.0D0
    do J = 1, NRHOX
      TERM1 = TERM2

      ! Optical depth increment to next layer
      if (J /= NRHOX) then
        D = TAUNU(J+1) - TAUNU(J)
      end if
      D = max(1.0D-10, D)

      if (D <= 0.01D0) then
        ! Small optical depth: Taylor series for E_3 integral
        ! TERM2 ≈ (1 - gamma + ln2 - lnD)*D/4 + D^2/12 - D^3/96 + D^4/720
        TERM2 = (EULER_M1 - log(D)) * D / 4.0D0 &
              + D**2 / 12.0D0 - D**3 / 96.0D0 + D**4 / 720.0D0
      else
        ! Standard E_3 path
        EX = 0.0D0
        if (D < 10.0D0) EX = EXPI(3, D)
        ! Cool-star stability patch: suppress E_3 for narrow Dtau range
        if (TEFF <= 4250.0D0 .and. D > 0.005D0 .and. D < 0.02D0) EX = 0.0D0
        TERM2 = 0.5D0 * (D + EX - 0.5D0) / D
      end if

      DIAGJ = TERM1 + TERM2

      ! dB_nu/dT = B_nu * h*nu / (kT^2) / (1 - e^{-h*nu/kT})
      if (NUMNU == 1) then
        DBDT = FLUX * 16.0D0 / T(J)
      else
        ! Guard: STIM = 1 - exp(-hv/kT) can vanish in Rayleigh-Jeans limit
        DBDT = BNU(J) * FREQ * HKT(J) / T(J) / max(STIM(J), 1.0D-30)
      end if

      ! Accumulate: kappa * (Lambda_diag - 1) / (1 - eps*Lambda_diag) * (1-eps) * dB/dT
      RDIAGJNU(J) = ABTOT(J) * (DIAGJ - 1.0D0) &
                   / (1.0D0 - ALPHA(J) * DIAGJ) &
                   * (1.0D0 - ALPHA(J)) * DBDT * RCOWT
      RDIAGJ(J) = RDIAGJ(J) + RDIAGJNU(J)
    end do
    return

  end if

  !=====================================================================
  ! MODE 3: Compute and apply temperature corrections
  !=====================================================================

  ! --- Compute needed derivatives ---
  call DERIV(RHOX, T, DTDRHX, NRHOX)
  call DERIV(RHOX, DLTDLP, DDLT, NRHOX)
  call DERIV(RHOX, ABROSS, DABROS, NRHOX)

  !---------------------------------------------------------------------
  ! (A) Prepare smoothed convective flux
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    CNVFLX(J) = 0.0D0
    if (IFCONV == 1 .and. J >= 3) CNVFLX(J) = FLXCNV(J)
  end do

  ! 1-2-1 smoothing filter on convective flux (interior points)
  do J = 2, NRHOX - 1
    SMOOTH(J) = 0.25D0 * CNVFLX(J-1) + 0.50D0 * CNVFLX(J) + 0.25D0 * CNVFLX(J+1)
  end do
  ! Asymmetric boundary kernel at deepest layer: 75-25 split
  SMOOTH(NRHOX) = 0.75D0 * CNVFLX(NRHOX) + 0.25D0 * CNVFLX(NRHOX-1)
  do J = 2, NRHOX
    CNVFLX(J) = SMOOTH(J)
    FLXCNV(J) = CNVFLX(J)
  end do

  !---------------------------------------------------------------------
  ! (B) Build integrating factor for flux correction
  !     CODRHX = d(ln G)/d(RHOX) incorporates opacity gradient and
  !     convective flux response to temperature perturbations
  !---------------------------------------------------------------------

  ! Initialize DDEL for all layers (safe default for non-convective case)
  DDEL(:) = 1.0D0

  do J = 1, NRHOX
    ABROSS_safe = max(ABROSS(J), ABROSS_FLOOR)
    RDABH(J) = RDABH(J) - FLXRAD(J) * DABROS(J) / ABROSS_safe

    ! Convective response: DEL = superadiabatic excess, D = radiative leak parameter
    DEL  = 1.0D0
    D    = 0.0D0

    if (CNVFLX(J) > 0.0D0 .and. FLXCNV0(J) > 0.0D0) then
      DEL = DLTDLP(J) - GRDADB(J)
      DEL = max(DEL, DEL_FLOOR)

      VCO = 0.5D0 * MIXLTH * sqrt(max(-0.5D0 * PTOTAL(J) / RHO(J) * DLRDLT(J), 0.0D0))
      FLUXCO = 0.5D0 * RHO(J) * HEATCP(J) * T(J) * MIXLTH / FOURPI

      if (MIXLTH > 0.0D0 .and. VCO > 0.0D0) then
        D = 8.0D0 * SIGMA_SB * T(J)**4 &
          / (ABROSS_safe * HSCALE(J) * RHO(J)) / (FLUXCO * FOURPI) / VCO
      end if

      TAUB = ABROSS_safe * RHO(J) * MIXLTH * HSCALE(J)
      D = D * TAUB**2 / (2.0D0 + TAUB**2)
      D = D**2 / 2.0D0
      DDEL(J) = (1.0D0 + D / (D + DEL)) / DEL
    end if

    ! Only include convective coupling when it's significant
    CNVFL = 0.0D0
    if (max(FLXRAD(J), 1.0D-30) > 0.0D0) then
      if (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) > 1.0D-3 .and. &
          FLXCNV0(J) / max(FLXRAD(J), 1.0D-30) > 1.0D-3) then
        CNVFL = CNVFLX(J)
      end if
    end if

    CODRHX(J) = (RDABH(J) &
               + CNVFL * (DTDRHX(J) / T(J) * (1.0D0 - 9.0D0 * D / (D + DEL)) &
                        + 1.5D0 * DDLT(J) / DEL * (1.0D0 + D / (D + DEL)))) &
              / (FLXRAD(J) + CNVFLX(J) * 1.5D0 * DLTDLP(J) / DEL &
                           * (1.0D0 + D / (D + DEL)))
  end do

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
  call INTEG(RHOX, CODRHX, G, NRHOX, 0.0D0)
  do J = 1, NRHOX
    G(J) = exp(G(J))
    GFLUX(J) = G(J) * (FLXRAD(J) + CNVFLX(J) - FLUX) &
             / (FLXRAD(J) + CNVFLX(J) * 1.5D0 * DLTDLP(J) * DDEL(J))

    ! Suppress GFLUX in convection-dominated layers (with hysteresis)
    if (IFCONV == 1 .and. &
        (CNVFLX(J) > 0.0D0 .or. FLXCNV0(J) > 0.0D0)) then
      FCONV_RATIO = max(CNVFLX(J), FLXCNV0(J)) &
                  / (max(CNVFLX(J), FLXCNV0(J)) + max(FLXRAD(J), 1.0D-30))
      GFLUX(J) = GFLUX(J) * (1.0D0 - FCONV_RATIO)
    end if
  end do

  call INTEG(TAUROS, GFLUX, DTAU, NRHOX, 0.0D0)
  do J = 1, NRHOX
    DTAU(J) = DTAU(J) / G(J)
    ! Clamp tau correction to ±tau/3 to prevent overshooting
    DTAU(J) = max(-TAUROS(J) / 3.0D0, min(TAUROS(J) / 3.0D0, DTAU(J)))
    ABROSS_safe = max(ABROSS(J), ABROSS_FLOOR)
    DTFLUX(J) = -DTAU(J) * DTDRHX(J) / ABROSS_safe
  end do

  TEFF25 = TEFF / 25.0D0

  !---------------------------------------------------------------------
  ! (C2) Suppress DTFLUX in convective layers.
  !      The Avrett-Krook integral handles radiative flux errors only;
  !      in convective layers, DTFLUX is attenuated by (1 - f_conv).
  !      The convective correction (adiabatic sweep) is applied later,
  !      after DTLAMB and DTSURF have been computed.
  !---------------------------------------------------------------------
  DTCONV(:) = 0.0D0

  if (IFCONV == 1) then
    do J = 1, NRHOX
      if (CNVFLX(J) > 0.0D0 .or. FLXCNV0(J) > 0.0D0) then
        FCONV_RATIO = max(CNVFLX(J), FLXCNV0(J)) &
                    / (max(CNVFLX(J), FLXCNV0(J)) + max(FLXRAD(J), 1.0D-30))
        DTFLUX(J) = DTFLUX(J) * (1.0D0 - FCONV_RATIO)
      end if
    end do
  end if

  !---------------------------------------------------------------------
  ! (D) DTLAMB: Lambda-iteration correction (optically thin layers)
  !     Uses thermal imbalance (J - S) and diagonal Lambda operator
  !     to drive local radiative balance.
  !---------------------------------------------------------------------

  ! Flux error as percentage
  do J = 1, NRHOX
    FLXERR(J) = (FLXRAD(J) + CNVFLX(J) - FLUX) / FLUX * 100.0D0
  end do
  call DERIV(TAUROS, FLXERR, FLXDRV, NRHOX)

  ! Find first layer where convection becomes significant
  ! (used to safely bound the DTLAMB backward-damping reach)
  JCONV = NRHOX + 1
  do J = 1, NRHOX
    if (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) >= 1.0D-5 .or. TAUROS(J) >= 1.0D0) then
      JCONV = J
      exit
    end if
  end do

  do J = 1, NRHOX
    ! In radiative layers, replace flux derivative with thermal imbalance
    if (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) < 1.0D-5) then
      ABROSS_safe = max(ABROSS(J), ABROSS_FLOOR)
      FLXDRV(J) = RJMINS(J) / ABROSS_safe / FLUX * 100.0D0
    end if

    ! Lambda correction: DT = -(dF/dtau) / (dRDIAGJ) * kappa_Ross
    DTLAMB(J) = -FLXDRV(J) * FLUX / 100.0D0 &
              / max(abs(RDIAGJ(J)), RDIAGJ_FLOOR) * sign(1.0D0, RDIAGJ(J)) &
              * max(ABROSS(J), ABROSS_FLOOR)

    ! Zero DTLAMB in convective / optically thick layers, and
    ! damp preceding layers to ensure smooth transition
    if (CNVFLX(J) / max(FLXRAD(J), 1.0D-30) >= 1.0D-5 .or. TAUROS(J) >= 1.0D0) then
      DTLAMB(J) = 0.0D0
      ! Damp 5 preceding layers (with bounds check)
      do K = 1, 5
        if (J - K >= 1) DTLAMB(J - K) = DTLAMB(J - K) / 2.0D0
      end do
    end if

    ! Clamp to ±Teff/25
    DTLAMB(J) = max(-TEFF25, min(TEFF25, DTLAMB(J)))
  end do

  !---------------------------------------------------------------------
  ! (E) DTSURF: surface boundary correction
  !     Uniform shift from flux error at tau=0, adjusted to not fight
  !     the integral of (DTFLUX + DTLAMB) between tau=0.1 and tau=2.
  !---------------------------------------------------------------------
  DTSUR = (FLUX - FLXRAD(1)) / FLUX * 0.25D0 * T(1)
  DTSUR = max(-TEFF25, min(TEFF25, DTSUR))

  do J = 1, NRHOX
    DUM(J) = DTFLUX(J) + DTLAMB(J)
  end do
  call INTEG(TAUROS, DUM, TINTEG, NRHOX, 0.0D0)
  XNEW_TMP(1) = 0.1D0
  IDUM = MAP1(TAUROS, TINTEG, NRHOX, XNEW_TMP, TONE_ARR, 1)
  XNEW_TMP(1) = 2.0D0
  IDUM = MAP1(TAUROS, TINTEG, NRHOX, XNEW_TMP, TTWO_ARR, 1)

  TAV = (TTWO_ARR(1) - TONE_ARR(1)) / 2.0D0
  if (DTSUR * TAV <= 0.0D0) TAV = 0.0D0
  if (abs(TAV) > abs(DTSUR)) TAV = DTSUR
  DTSUR = DTSUR - TAV

  do J = 1, NRHOX
    DTSURF(J) = DTSUR
    HRATIO(J) = CNVFLX(J) / (CNVFLX(J) + max(FLXRAD(J), 1.0D-30))
  end do

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

  if (IFCONV == 1) then

    ! --- Find anchor: shallowest layer with convective flux ---
    JANCHOR = 0
    do J = 1, NRHOX
      if (CNVFLX(J) > 0.0D0 .or. FLXCNV0(J) > 0.0D0) then
        JANCHOR = J
        exit
      end if
    end do

    ! --- Adiabatic sweep from anchor downward ---
    if (JANCHOR > 0) then
      ! Anchor gets normal radiative corrections only (DTCONV=0 there).
      ! T_SWEEP tracks the corrected temperature for integration.
      T_SWEEP(JANCHOR) = T(JANCHOR) + DTFLUX(JANCHOR) &
                       + DTLAMB(JANCHOR) + DTSURF(JANCHOR)

      do J = JANCHOR + 1, NRHOX
        ! Convective fraction at this layer (with hysteresis)
        if (CNVFLX(J) > 0.0D0 .or. FLXCNV0(J) > 0.0D0) then
          FCONV_RATIO = max(CNVFLX(J), FLXCNV0(J)) &
                      / (max(CNVFLX(J), FLXCNV0(J)) + max(FLXRAD(J), 1.0D-30))
        else
          FCONV_RATIO = 0.0D0
        end if

        if (FCONV_RATIO > 0.0D0 .and. PTOTAL(J) > 0.0D0 &
            .and. PTOTAL(J-1) > 0.0D0) then
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
        else
          ! Non-convective layer below anchor: no DTCONV, pass through
          T_SWEEP(J) = T(J) + DTFLUX(J) + DTLAMB(J) + DTSURF(J)
        end if
      end do
    end if
  end if

  !---------------------------------------------------------------------
  ! (F) Total correction and iteration damping
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    T1(J) = DTFLUX(J) + DTLAMB(J) + DTSURF(J) + DTCONV(J)
  end do

  ! Iteration damping: compare current correction T1 against previous
  ! iteration's correction OLDT1 to detect convergence behavior.
  !   Same sign, shrinking  → accelerate by 1.25x (monotone convergence)
  !   Same sign, growing    → cap at previous magnitude (prevent runaway)
  !   Sign flip             → damp by 0.5x (oscillation control)
  ! Skip damping on first iteration only (no history).
  do J = 1, NRHOX
    if (ITER == 1) then
      ! No damping — accept raw correction
    else if (OLDT1(J) * T1(J) > 0.0D0) then
      ! Same sign as previous iteration
      if (abs(T1(J)) < abs(OLDT1(J))) then
        ! Shrinking: accelerate by 25%
        T1(J) = T1(J) * 1.25D0
      else
        ! Growing: cap at previous magnitude to prevent runaway
        T1(J) = sign(abs(OLDT1(J)), T1(J))
      end if
    else if (OLDT1(J) * T1(J) < 0.0D0) then
      ! Sign flip: damp by 50%
      T1(J) = T1(J) * 0.5D0
    end if
    OLDT1(J) = T1(J)
  end do

  ! Diagnostic output (after damping, so T1 reflects what is actually applied)
  if (IFPRNT(ITER) /= 0) then
    write(67, 100) (J, log10(max(TAUROS(J),1.0D-30)), T(J), DTLAMB(J), &
                    DTSURF(J), DTFLUX(J), DTCONV(J), T1(J), HRATIO(J), &
                    FLXERR(J), FLXDRV(J), DLTDLP(J), GRDADB(J), J=1,NRHOX)
100 format('0', 2X, 'lgTAUROS', 6X, 'T', 6X, &
      'DTLAMB   DTSURF   DTFLUX   DTCONV', 5X, &
      'T1   CONV/TOTAL      ERROR     DERIV   NABLA  NABLA_AD', / &
      (I3, F8.3, F10.1, 5F9.1, 1X, 1PE11.3, 1X, &
       0P, 2F10.3, 2F8.4))
    flush(67)
  end if

  !---------------------------------------------------------------------
  ! (G) Compute RHOX correction to maintain constant TAUROS grid
  !     Uses finite-difference: run TTAUP with T and T+DT, compare
  !     total pressures to infer needed RHOX adjustment.
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    TPLUS(J)  = T(J) + T1(J)
    TAUSTD(J) = 10.0D0**(TAU1LG + (J - 1) * STEPLG)
  end do

  IDUM = MAP1(TAUROS, T,    NRHOX, TAUSTD, TNEW1,  NRHOX)
  IDUM = MAP1(TAUROS, PRAD, NRHOX, TAUSTD, PRDNEW, NRHOX)
  call TTAUP(TNEW1, TAUSTD, AB1, PTOT1, P1, PRDNEW, PTURB, VTURB, GRAV, NRHOX)

  IDUM = MAP1(TAUROS, TPLUS, NRHOX, TAUSTD, TNEW2, NRHOX)
  call TTAUP(TNEW2, TAUSTD, AB2, PTOT2, P2, PRDNEW, PTURB, VTURB, GRAV, NRHOX)

  do J = 1, NRHOX
    PPP(J) = (PTOT2(J) - PTOT1(J)) / PTOT1(J)
  end do
  IDUM = MAP1(TAUSTD, PPP, NRHOX, TAUROS, RRR, NRHOX)
  do J = 1, NRHOX
    DRHOX(J) = RRR(J) * RHOX(J)
  end do

  !---------------------------------------------------------------------
  ! (H) Apply temperature correction
  !---------------------------------------------------------------------
  do J = 1, NRHOX
     !apply damping to the temperature correction
    ! T1(J) = sign(min(abs(T1(J)), 0.02D0 * T(J)), T1(J))
     T(J) = T(J) + T1(J)
  end do

  ! Optional smoothing
  if (J1SMOOTH > 0) then
    do J = J1SMOOTH, J2SMOOTH
      TSMOOTH(J) = WTJM1 * T(J-1) + WTJ * T(J) + WTJP1 * T(J+1)
    end do
    do J = J1SMOOTH, J2SMOOTH
      T(J) = TSMOOTH(J)
    end do
  end if

  ! Force monotonicity: T must increase inward
  do I = 2, NRHOX
    J = NRHOX + 1 - I
    T(J) = min(T(J), T(J+1) - 1.0D0)
  end do

  ! Floor temperature options
  ! various opacity tables only extend to 2000K
  do J = 1, NRHOX
     T(J) = max(T(J), 2000.0D0)
    ! T(J) = max(T(J),0.7*TEFF)
  end do

  !---------------------------------------------------------------------
  ! (I) Recompute thermodynamic quantities from corrected T
  !---------------------------------------------------------------------
  IFUDGE = 0

  do J = 1, NRHOX
    TK(J)   = KBOL * T(J)
    HKT(J)  = HPLANCK / TK(J)
    HCKT(J) = HKT(J) * CLIGHT
    TKEV(J) = 8.6171D-5 * T(J)
    TLOG(J) = log(T(J))
  end do

  if (IFUDGE == 1) return

  !---------------------------------------------------------------------
  ! (J) Apply RHOX correction and remap atmosphere onto standard
  !     Rosseland optical depth grid TAUSTD
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    RHOX(J) = RHOX(J) + DRHOX(J)
  end do

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
  do J = 1, NRHOX
    if (TAUROS(1) > TAUSTD(J) .and. REMAP(J,5) < 0.0D0) then
      REMAP(J,5) = ABROSS(1)
    end if
  end do

  ! Copy remapped values back to module arrays
  do J = 1, NRHOX
    RHOX(J)   = REMAP(J, 1)
    T(J)      = REMAP(J, 2)
    TK(J)     = KBOL * T(J)
    HKT(J)    = HPLANCK / TK(J)
    HCKT(J)   = HKT(J) * CLIGHT
    TKEV(J)   = 8.6171D-5 * T(J)
    TLOG(J)   = log(T(J))
    P(J)      = REMAP(J, 3)
    XNE(J)    = REMAP(J, 4)
    ABROSS(J) = REMAP(J, 5)
    PRAD(J)   = REMAP(J, 6)
    PRADK(J)  = PRAD(J) + PRADK0
    VTURB(J)  = REMAP(J, 7)
    BMIN(J)   = REMAP(J, 8)
    PTURB(J)  = REMAP(J, 9)
  end do

  ! Remap hydrogen NLTE departure coefficients
  do I = 1, 6
    IDUM = MAP1(TAUROS, BHYD(1,I), NRHOX, TAUSTD, REMAP(1,1), NRHOX)
    do J = 1, NRHOX
      BHYD(J,I) = REMAP(J, 1)
    end do
  end do

  ! Replace old TAUROS with the standard grid
  do J = 1, NRHOX
    TAUROS(J) = TAUSTD(J)
  end do

  return

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
  
  implicit none

  ! --- Arguments ---
  integer, intent(in) :: MODE
  real*8,  intent(in) :: RCOWT    ! quadrature weight for current frequency

  ! --- Number of bound levels in the model atom ---
  integer, parameter :: NLEV = 6   ! bound levels solved explicitly
  integer, parameter :: NATOM = 8  ! total levels (6 bound + H+ + H-)

  ! --- Persistent locals (accumulate across MODE 1 → 2 → 3) ---
  real*8, save :: QRADIK(kw, NLEV)  ! radiative photoionization rates (level i → continuum)
  real*8, save :: QRADKI(kw, NLEV)  ! radiative recombination rates (continuum → level i)
  real*8, save :: DQRAD(kw, NLEV)   ! dQ_recomb/dT (for first-order T correction)
  real*8, save :: QRDHMK(kw)        ! H- radiative detachment rate (H- → H + e)
  real*8, save :: QRDKHM(kw)        ! H- radiative attachment rate (H + e → H-)
  real*8, save :: DQRD(kw)          ! dQ_Hminus_attach/dT
  real*8, save :: TOLD(kw)          ! temperature at MODE=1 (for dT correction)

  ! --- MODE 2 locals ---
  real*8  :: HCONT(NLEV)         ! hydrogen bound-free cross-sections at current freq
  real*8  :: HMINBF              ! H- bound-free cross-section at current freq
  real*8  :: RFRWT               ! 4*pi / (h*nu) * RCOWT — rate prefactor
  real*8  :: HVC                 ! 2*h*nu*(nu/c)^2 — stimulated emission correction
  real*8  :: RJ, RJE, RJEDT     ! rate building blocks at each depth

  ! --- MODE 3 locals ---
  real*8  :: QCOLL(NATOM, NATOM) ! electron collision rate matrix
  real*8  :: A(NLEV, NLEV)       ! rate equation matrix (destroyed by SOLVIT)
  real*8  :: RIGHT(NLEV)         ! RHS vector → solution vector
  integer :: IPIVOT(NLEV)        ! pivot scratch for SOLVIT
  real*8  :: DT                  ! temperature change since MODE 1
  real*8  :: THETA               ! 5040/T (excitation temperature parameter)
  real*8  :: TH                  ! ionization potential / kT = 13.595 / T_eV
  real*8  :: Y, Z_lev            ! real-valued level indices
  real*8  :: GIK, X0, Q          ! collision rate intermediates
  real*8  :: QELECT              ! H- collisional detachment by electrons
  real*8  :: QASSOC              ! H- associative detachment (H + H → H2 + e)
  real*8  :: QCHARG              ! H- charge exchange (H- + H+ → 2H)
  real*8  :: DENOM               ! denominator for BMIN (with guard)

  ! --- Loop indices ---
  integer :: I, J, K, L, N

  ! --- Oscillator strengths for bound-bound transitions ---
  ! F(I,K) = oscillator strength for transition I → K
  ! From Peterson (1968); upper triangle only (F is symmetric by convention here)
  real*8, parameter :: F(NATOM, NATOM) = reshape( [ &
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

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING STATEQ'

  !=====================================================================
  ! MODE 1: Zero radiative rate accumulators and save current temperature
  !=====================================================================
  if (MODE == 1) then

    do I = 1, NLEV
      do J = 1, NRHOX
        QRADIK(J, I) = 0.0D0
        QRADKI(J, I) = 0.0D0
        DQRAD(J, I)  = 0.0D0
      end do
    end do
    do J = 1, NRHOX
      TOLD(J)   = T(J)
      QRDHMK(J) = 0.0D0
      QRDKHM(J) = 0.0D0
      DQRD(J)   = 0.0D0
    end do
    return

  !=====================================================================
  ! MODE 2: Accumulate radiative rates at current frequency
  !=====================================================================
  else if (MODE == 2) then

    ! Rate prefactor: (4*pi / h*nu) * quadrature weight
    RFRWT = FOURPI / HPLANCK * RCOWT / FREQ

    ! Stimulated emission term: 2*h*nu * (nu/c)^2
    HVC = 2.0D0 * HPLANCK * FREQ * (FREQ / CLIGHT)**2

    ! Hydrogen bound-free cross-sections for levels 2-6
    do N = 2, NLEV
      HCONT(N) = COULX(N, FREQ, 1.0D0)
    end do

    ! H- bound-free cross-section (polynomial fit)
    HMINBF = 0.0D0
    if (FREQ > 1.8259D14 .and. FREQ < 2.111D14) then
      HMINBF = 3.695D-16 + (-1.251D-1 + 1.052D13 / FREQ) / FREQ
    end if
    if (FREQ >= 2.111D14) then
      HMINBF = 6.801D-20 + (5.358D-3 + (1.481D13 &
             + (-5.519D27 + 4.808D41 / FREQ) / FREQ) / FREQ) / FREQ
    end if

    ! Accumulate rates at each depth
    do J = 1, NRHOX
      ! RJ  = (4*pi*J_nu / h*nu) * weight — photoionization driver
      ! RJE = (4*pi*e^{-hv/kT} * (J_nu + 2hv^3/c^2) / h*nu) * weight — recombination
      ! RJEDT = d(RJE)/dT — temperature derivative for MODE 3 correction
      RJ    = RFRWT * JNU(J)
      RJE   = RFRWT * EHVKT(J) * (JNU(J) + HVC)
      RJEDT = RJE * HKT(J) * FREQ / T(J)

      ! Hydrogen levels 2-6 (level 1 has no photoionization at these frequencies)
      do I = 2, NLEV
        QRADIK(J, I) = QRADIK(J, I) + HCONT(I) * RJ
        QRADKI(J, I) = QRADKI(J, I) + HCONT(I) * RJE
        DQRAD(J, I)  = DQRAD(J, I)  + HCONT(I) * RJEDT
      end do

      ! H- rates
      QRDHMK(J) = QRDHMK(J) + HMINBF * RJ
      QRDKHM(J) = QRDKHM(J) + HMINBF * RJE
      DQRD(J)   = DQRD(J)   + HMINBF * RJEDT
    end do
    return

  end if

  !=====================================================================
  ! MODE 3: Compute collision rates, solve rate equations, update BHYD/BMIN
  !=====================================================================

  !-------------------------------------------------------------------
  ! (A) H- departure coefficient: detailed balance of formation/destruction
  !-------------------------------------------------------------------
  if (IFPRNT(ITER) > 0) then
    write(6, 201)
201 format('1', /////, 36X, 'HMINUS STATISTICAL EQUILIBRIUM' / &
      10X, 'RHOX', 7X, 'QELECT', 6X, 'QASSOC', 6X, 'QCHARG', &
      6X, 'QRDKHM', 6X, 'QRDHMK', 7X, 'BMIN')
  end if

  do J = 1, NRHOX
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

    if (IFPRNT(ITER) > 0) then
      write(6, 211) J, RHOX(J), QELECT, QASSOC, QCHARG, &
                     QRDKHM(J), QRDHMK(J), BMIN(J)
211   format(I5, 1P6E12.3, 0PF10.4)
    end if
  end do

  !-------------------------------------------------------------------
  ! (B) Hydrogen 6-level atom: build and solve rate equations
  !-------------------------------------------------------------------
  if (IFPRNT(ITER) > 0) then
    write(6, 31)
31  format('1', /////, 30X, &
      'STATISTICAL EQUILIBRIUM RATES    RATE=SIGN(ALOG10(MAX(ABS(RATE*1.E20),1.)),RATE)', /, &
      '0 RAD   1-K   K-1   2-K   K-2   3-K   K-3   4-K   K-4   5-K', &
      '   K-5   6-K   K-6   COLL  1-K   2-K   3-K   4-K   5-K   6-K   5-8', &
      '   6-8  ', /, &
      '  COLL  1-2   1-3   1-4   1-5   1-6   1-7   2-3   2-4   2-5', &
      '   2-6   2-7   3-4   3-5   3-6   3-7   4-5   4-6   4-7   5-6   5-7', &
      '   6-7  ')
  end if

  do J = 1, NRHOX
    DT = T(J) - TOLD(J)
    TH = 13.595D0 / TKEV(J)    ! chi_H / kT

    !--- Electron collision rates ---
    ! Diagonal: bound-free collisional ionization
    ! Off-diagonal: bound-bound collisional excitation/de-excitation
    do I = 1, NATOM
      Y = dble(I)

      ! Bound-free collisional ionization rate (Seaton-like)
      QCOLL(I, I) = 2.2D-8 * Y**3 / sqrt(TH) * exp(-TH / Y**2) * XNE(J)

      if (I < NATOM) then
        do K = I + 1, NATOM
          Z_lev = dble(K)
          GIK = 1.0D0 / Y**2 - 1.0D0 / Z_lev**2
          X0  = TH * GIK

          ! Burke, Ormonde & Whitaker cross-section with E_1 and E_5 integrals
          Q = 2.186D-10 * F(I, K) / GIK**2 * X0 * sqrt(T(J)) &
            * (EXPI(1, X0) + 0.148D0 * X0 * EXPI(5, X0))

          QCOLL(I, K) = Q * XNE(J)
          ! Detailed balance: reverse rate includes Boltzmann factor
          QCOLL(K, I) = QCOLL(I, K) * (Y / Z_lev)**2 * exp(X0)
        end do
      end if
    end do

    !--- Assemble 6x6 rate matrix ---
    ! A(I,I) = total destruction rate for level I (radiative + collisional)
    ! A(I,K) = -collisional rate I → K (off-diagonal, K ≠ I)
    ! RIGHT(I) = total formation rate for level I from continuum/H-
    do I = 1, NLEV
      ! Start diagonal with radiative ionization rate
      A(I, I) = QRADIK(J, I)

      ! First-order T correction to recombination rate
      QRADKI(J, I) = QRADKI(J, I) + DQRAD(J, I) * DT

      ! RHS: recombination + collisional ionization + coupling to levels 7,8
      RIGHT(I) = QRADKI(J, I) + QCOLL(I, I) + QCOLL(I, 7) + QCOLL(I, 8)

      ! Add all collisional rates to diagonal (total destruction)
      do K = 1, NATOM
        A(I, I) = A(I, I) + QCOLL(I, K)
      end do

      ! Off-diagonal: collisional coupling between bound levels
      if (I < NLEV) then
        do K = I + 1, NLEV
          A(I, K) = -QCOLL(I, K)
          A(K, I) = -QCOLL(K, I)
        end do
      end if
    end do

    !--- Solve for departure coefficients ---
    call SOLVIT(A, NLEV, RIGHT, IPIVOT)
    do L = 1, NLEV
      BHYD(J, L) = RIGHT(L)
    end do

    !--- Diagnostic output (log-scaled rates) ---
    if (IFPRNT(ITER) > 1) then
      ! Convert rates to sign-preserving log scale for display
      do I = 1, NLEV
        QRADKI(J, I) = sign(log10(max(abs(QRADKI(J, I) * 1.0D20), 1.0D0)), &
                             QRADKI(J, I))
        QRADIK(J, I) = sign(log10(max(abs(QRADIK(J, I) * 1.0D20), 1.0D0)), &
                             QRADIK(J, I))
      end do
      do I = 1, NATOM
        do K = 1, NATOM
          QCOLL(I, K) = sign(log10(max(abs(QCOLL(I, K) * 1.0D20), 1.0D0)), &
                              QCOLL(I, K))
        end do
      end do

      write(6, 100) J, (QRADIK(J, I), QRADKI(J, I), I=1,NLEV), &
                        (QCOLL(I, I), I=1,NLEV), QCOLL(5, 8), QCOLL(6, 8)
100   format('0', I5, 12F6.2, 6X, 8F6.2)
      write(6, 110) (QCOLL(1, K), K=2,7), (QCOLL(2, K), K=3,7), &
                    (QCOLL(3, K), K=4,7), (QCOLL(4, K), K=5,7), &
                    (QCOLL(5, K), K=6,7), QCOLL(6, 7)
110   format(6X, 21F6.2)
    end if

  end do  ! depth loop

  !-------------------------------------------------------------------
  ! (C) Print final departure coefficients
  !-------------------------------------------------------------------
  write(6, 170) (J, RHOX(J), (BHYD(J, I), I=1,NLEV), J=1,NRHOX)
170 format('1', /////, 30X, 'STATISTICAL EQUILIBRIUM FOR HYDROGEN' / &
      15X, 'RHOX', 10X, 'B1', 8X, 'B2', 8X, 'B3', 8X, 'B4', 8X, 'B5', 8X, 'B6', / &
      (8X, I2, 1PE11.4, 1X, 0P, 6F10.4))

  return

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
  
  implicit none

  ! --- Arguments ---
  integer, intent(in) :: MODE
  real*8,  intent(in) :: RCOWT    ! quadrature weight for current frequency

  ! --- Named constants: FOURPI_OVER_C now from mod_constants ---

  ! --- Persistent local: frequency-integrated flux (for error scaling) ---
  real*8, save :: H_TOTAL(kw)    ! integrated Eddington flux at each depth

  ! --- Local variables ---
  real*8  :: ERRORMAX            ! max flux ratio (for stability scaling)
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING RADIAP'

  !=====================================================================
  ! MODE 1: Zero frequency-integral accumulators
  !=====================================================================
  if (MODE == 1) then

    do J = 1, NRHOX
      H_TOTAL(J) = 0.0D0
      RADEN(J)   = 0.0D0
      ACCRAD(J)  = 0.0D0
    end do
    PRADK0 = 0.0D0
    return

  !=====================================================================
  ! MODE 2: Accumulate frequency integrals at current wavelength
  !=====================================================================
  else if (MODE == 2) then

    do J = 1, NRHOX
      ! Mean radiation energy density: J_nu * weight
      RADEN(J) = RADEN(J) + JNU(J) * RCOWT

      ! Total Eddington flux: H_nu * weight
      H_TOTAL(J) = H_TOTAL(J) + HNU(J) * RCOWT

      ! Radiative acceleration: kappa_nu * H_nu * weight
      ACCRAD(J) = ACCRAD(J) + ABTOT(J) * HNU(J) * RCOWT
    end do

    ! Surface K-integral (second moment, = radiation pressure at tau=0)
    PRADK0 = PRADK0 + KNU(1) * RCOWT
    return

  end if

  !=====================================================================
  ! MODE 3: Convert to physical units and integrate for P_rad
  !=====================================================================

  ! Convert from frequency sums to physical units: multiply by 4*pi/c
  ! This converts:
  !   RADEN:  sum(J_nu * dnu) → (4*pi/c) * J = radiation energy density
  !   ACCRAD: sum(kappa * H_nu * dnu) → (4*pi/c) * kappa * H = radiative acceleration
  ERRORMAX = 0.0D0
  do J = 1, NRHOX
    RADEN(J)  = RADEN(J) * FOURPI_OVER_C
    ACCRAD(J) = ACCRAD(J) * FOURPI_OVER_C

    ! Stability safeguard: if the integrated flux exceeds the target
    ! (FLUX = sigma*Teff^4 / 4*pi), the radiative acceleration is
    ! artificially scaled down to prevent runaway radiation pressure
    ! in early iterations with large flux errors.
    if (H_TOTAL(J) / FLUX > 1.0D0) then
      ACCRAD(J) = ACCRAD(J) * FLUX / H_TOTAL(J)
    end if
    ERRORMAX = max(ERRORMAX, H_TOTAL(J) / FLUX)
  end do

  ! Same scaling for surface K-integral
  PRADK0 = PRADK0 * FOURPI_OVER_C
  if (ERRORMAX > 1.0D0) PRADK0 = PRADK0 / ERRORMAX

  ! Integrate radiative acceleration inward to get radiation pressure:
  !   P_rad(RHOX) = integral from 0 to RHOX of ACCRAD * dRHOX
  ! with boundary value ACCRAD(1)*RHOX(1) (linear extrapolation to surface)
  call INTEG(RHOX, ACCRAD, PRAD, NRHOX, ACCRAD(1) * RHOX(1))

  ! Total radiation pressure = depth-integrated + surface boundary term
  do J = 1, NRHOX
    PRADK(J) = PRAD(J) + PRADK0
  end do

  return

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

  implicit none

  integer, intent(in) :: MODE
  real*8,  intent(in) :: RCOWT

  ! --- Local variables ---
  real*8  :: DBDT
  integer :: J

  ! FOUR_SIGMA_OVER_PI now from mod_constants

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING ROSS'
  select case (MODE)

  ! --- MODE 1: Zero the Rosseland mean accumulator ---
  case (1)
    do J = 1, NRHOX
      ABROSS(J) = 0.0D0
    end do

  ! --- MODE 2: Accumulate (dB_nu/dT) / kappa_nu at current frequency ---
  case (2)
    do J = 1, NRHOX
      ! dB_nu/dT = B_nu * (h*nu/kT) / T / (1 - e^{-h*nu/kT})
      !          = BNU * FREQ * HKT / T / STIM
      DBDT = BNU(J) * FREQ * HKT(J) / T(J) / STIM(J)
      ! Single-frequency fallback: use total dB/dT = (4*sigma/pi)*T^3
      if (NUMNU == 1) DBDT = FOUR_SIGMA_OVER_PI * T(J)**3
      ABROSS(J) = ABROSS(J) + DBDT / ABTOT(J) * RCOWT
    end do

  ! --- MODE 3: Finalize — compute kappa_Ross and tau_Ross scale ---
  case (3)
    ! Convert accumulated sum to Rosseland mean opacity:
    !   kappa_Ross = (4*sigma/pi)*T^3 / sum[ (dB/dT)/kappa_nu * dnu ]
    do J = 1, NRHOX
      ABROSS(J) = FOUR_SIGMA_OVER_PI * T(J)**3 / ABROSS(J)
    end do
    ! Integrate kappa_Ross over column mass to get Rosseland optical depth
    call INTEG(RHOX, ABROSS, TAUROS, NRHOX, ABROSS(1) * RHOX(1))

  end select

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

  implicit none

  integer, intent(in)  :: N
  real*8,  intent(in)  :: X(*), F(*)
  real*8,  intent(out) :: DFDX(*)

  ! --- Local variables ---
  real*8  :: S, SCALE, D_RIGHT, D_LEFT, T_RIGHT, T_LEFT
  integer :: J

  ! --- Endpoint derivatives: one-sided finite differences ---
  DFDX(1) = (F(2) - F(1)) / (X(2) - X(1))
  DFDX(N) = (F(N) - F(N-1)) / (X(N) - X(N-1))
  if (N == 2) return

  ! --- Sign of the X spacing (handles decreasing grids) ---
  S = abs(X(2) - X(1)) / (X(2) - X(1))

  ! --- Interior points: tangent half-angle averaging ---
  do J = 2, N - 1
    ! Normalization scale to prevent overflow
    SCALE = max(abs(F(J-1)), abs(F(J)), abs(F(J+1))) / abs(X(J))
    if (SCALE == 0.0D0) SCALE = 1.0D0

    ! Scaled right and left finite-difference slopes
    D_RIGHT = (F(J+1) - F(J)) / (X(J+1) - X(J)) / SCALE
    D_LEFT  = (F(J) - F(J-1)) / (X(J) - X(J-1)) / SCALE

    ! Convert slopes to half-angle tangents: t = d / (S*sqrt(1+d^2) + 1)
    T_RIGHT = D_RIGHT / (S * sqrt(1.0D0 + D_RIGHT**2) + 1.0D0)
    T_LEFT  = D_LEFT  / (S * sqrt(1.0D0 + D_LEFT**2)  + 1.0D0)

    ! Tangent addition formula: tan(a+b) = (tan(a)+tan(b))/(1-tan(a)*tan(b))
    DFDX(J) = (T_RIGHT + T_LEFT) / (1.0D0 - T_RIGHT * T_LEFT) * SCALE
  end do

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

  implicit none

  integer, intent(in)  :: N
  real*8,  intent(in)  :: X(*), F(*), START
  real*8,  intent(out) :: FINT(*)

  ! --- Local variables ---
  real*8  :: A(kw), B(kw), C(kw)
  real*8  :: XLO, XHI
  integer :: I

  ! Fit piecewise parabolas to F(X)
  call PARCOE(F, X, A, B, C, N)

  ! Integrate by summing the analytic parabolic integrals
  FINT(1) = START
  do I = 1, N - 1
    XLO = X(I)
    XHI = X(I+1)
    FINT(I+1) = FINT(I) + (A(I) + B(I) * 0.5D0 * (XHI + XLO) &
                + C(I) / 3.0D0 * ((XHI + XLO) * XHI + XLO * XLO)) &
                * (XHI - XLO)
  end do

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

  implicit none

  integer, intent(in)  :: N
  real*8,  intent(in)  :: F(*), X(*)
  real*8,  intent(out) :: A(*), B(*), C(*)

  ! --- Local variables ---
  real*8  :: D, WT
  integer :: J

  ! --- Step 1: Endpoint linear fits (C = 0) ---
  C(1) = 0.0D0
  B(1) = (F(2) - F(1)) / (X(2) - X(1))
  A(1) = F(1) - X(1) * B(1)

  C(N) = 0.0D0
  B(N) = (F(N) - F(N-1)) / (X(N) - X(N-1))
  A(N) = F(N) - X(N) * B(N)

  if (N == 2) return

  ! --- N = 3: use linear interpolation on both intervals ---
  !     Copy last-interval linear fit to the middle point.
  !     (Returning early prevents the post-loop code from accessing
  !     out-of-bounds F(4) and X(4).)
  if (N == 3) then
    A(2) = A(3)
    B(2) = B(3)
    C(2) = C(3)
    return
  end if

  ! --- Step 2: Parabolic fits at interior points (N >= 4) ---
  !     For each interior point J, fit a parabola through
  !     (X_{J-1}, F_{J-1}), (X_J, F_J), (X_{J+1}, F_{J+1}).
  !     C(J) is the second divided difference (curvature term).
  do J = 2, N - 1
    D = (F(J) - F(J-1)) / (X(J) - X(J-1))
    C(J) = F(J+1) / ((X(J+1) - X(J)) * (X(J+1) - X(J-1))) &
         - F(J)   / ((X(J) - X(J-1))   * (X(J+1) - X(J)))   &
         + F(J-1) / ((X(J) - X(J-1))   * (X(J+1) - X(J-1)))
    B(J) = D - (X(J) + X(J-1)) * C(J)
    A(J) = F(J-1) - X(J-1) * D + X(J) * X(J-1) * C(J)
  end do

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
  do J = 2, N - 1
    if (C(J) /= 0.0D0) then
      WT = abs(C(J+1)) / (abs(C(J+1)) + abs(C(J)))
      A(J) = A(J+1) + WT * (A(J) - A(J+1))
      B(J) = B(J+1) + WT * (B(J) - B(J+1))
      C(J) = C(J+1) + WT * (C(J) - C(J+1))
    end if
  end do

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

  implicit none

  integer, intent(in)  :: NOLD, NNEW
  real*8,  intent(in)  :: XOLD(*), FOLD(*), XNEW(*)
  real*8,  intent(out) :: FNEW(*)
  integer              :: MAP1

  ! --- Local variables ---
  real*8  :: A, B, C                   ! Current interpolation coefficients
  real*8  :: ABAC, BBAC, CBAC          ! Backward parabola coefficients
  real*8  :: AFOR, BFOR, CFOR          ! Forward parabola coefficients
  real*8  :: D, WT
  integer :: K, L, LL, L1, L2
  logical :: past_right_edge           ! True if XNEW(K) is beyond XOLD(NOLD)

  L  = 2     ! Pointer into XOLD: XNEW(K) is in [XOLD(L-1), XOLD(L)]
  LL = 0     ! Previous L value (0 = no coefficients computed yet)

  do K = 1, NNEW

    ! -----------------------------------------------------------------
    ! Step 1: Advance L to bracket XNEW(K) in [XOLD(L-1), XOLD(L)]
    ! -----------------------------------------------------------------
    past_right_edge = .false.
    do while (XNEW(K) >= XOLD(L))
      L = L + 1
      if (L > NOLD) then
        ! Beyond right edge of old grid — use boundary linear fit
        past_right_edge = .true.
        exit
      end if
    end do

    ! -----------------------------------------------------------------
    ! Step 2: If still in the same interval, reuse existing coefficients
    ! -----------------------------------------------------------------
    if (.not. past_right_edge .and. L == LL) then
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      cycle
    end if

    ! -----------------------------------------------------------------
    ! Step 3: Boundary cases — use linear interpolation
    !         (L=2,3: too few points to the left for a parabola;
    !          past right edge: linear extrapolation from last two points)
    ! -----------------------------------------------------------------
    if (L <= 3 .or. past_right_edge) then
      L = min(NOLD, L)
      if (.not. past_right_edge .and. L == LL) then
        FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
        cycle
      end if
      C = 0.0D0
      B = (FOLD(L) - FOLD(L-1)) / (XOLD(L) - XOLD(L-1))
      A = FOLD(L) - XOLD(L) * B
      LL = L
      FNEW(K) = A + B * XNEW(K)
      cycle
    end if

    ! -----------------------------------------------------------------
    ! Step 4: Interior — compute backward parabola through (L-2, L-1, L)
    ! -----------------------------------------------------------------
    L1 = L - 1

    ! Check if we can reuse the previous forward parabola as backward
    if (L == LL + 1 .and. L /= 3 .and. L /= 4) then
      CBAC = CFOR
      BBAC = BFOR
      ABAC = AFOR
    else
      ! Compute fresh backward parabola
      L2 = L - 2
      D = (FOLD(L1) - FOLD(L2)) / (XOLD(L1) - XOLD(L2))
      CBAC = FOLD(L)  / ((XOLD(L) - XOLD(L1)) * (XOLD(L) - XOLD(L2))) &
           + (FOLD(L2) / (XOLD(L) - XOLD(L2)) &
           -  FOLD(L1) / (XOLD(L) - XOLD(L1))) &
           / (XOLD(L1) - XOLD(L2))
      BBAC = D - (XOLD(L1) + XOLD(L2)) * CBAC
      ABAC = FOLD(L2) - XOLD(L2) * D + XOLD(L1) * XOLD(L2) * CBAC
    end if

    ! -----------------------------------------------------------------
    ! Step 5: At right edge (L=NOLD) — use backward parabola only
    ! -----------------------------------------------------------------
    if (L >= NOLD) then
      A = ABAC;  B = BBAC;  C = CBAC
      LL = L
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      cycle
    end if

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
    if (abs(CFOR) /= 0.0D0) WT = abs(CFOR) / (abs(CFOR) + abs(CBAC))
    A = AFOR + WT * (ABAC - AFOR)
    B = BFOR + WT * (BBAC - BFOR)
    C = CFOR + WT * (CBAC - CFOR)
    LL = L

    ! -----------------------------------------------------------------
    ! Evaluate interpolant at XNEW(K)
    ! -----------------------------------------------------------------
    FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)

  end do

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

  implicit none

  integer, intent(in)    :: N
  real*8,  intent(inout) :: A(N,N), B(N)
  integer, intent(out)   :: IPIVOT(N)

  ! --- Local variables ---
  real*8  :: PIVOT, T, C
  integer :: I, J, K, M

  ! Threshold for treating a pivot as effectively zero.
  ! Machine epsilon for REAL*8 is ~2.2e-16; we use a generous threshold
  ! relative to the column's maximum element to catch near-singular cases.
  real*8, parameter :: PIVOT_TOL = 1.0D-30

  ! =================================================================
  ! Phase 1: LU factorization with partial pivoting
  ! =================================================================

  do I = 1, N - 1

    ! --- (a) Find pivot: row with largest |A(K,I)| for K >= I ---
    M = I
    do K = I + 1, N
      if (abs(A(K,I)) > abs(A(M,I))) M = K
    end do
    IPIVOT(I) = M

    ! --- (b) Swap rows I and M for columns I+1..N ---
    !     Column I is not swapped; its entries are handled implicitly
    !     through the multiplier computation below.
    if (M /= I) then
      do K = I + 1, N
        T = A(I,K)
        A(I,K) = A(M,K)
        A(M,K) = T
      end do
    end if

    ! --- (c) Store reciprocal pivot on diagonal ---
    if (abs(A(M,I)) < PIVOT_TOL) then
      write(6,'(A,I4,A,1PE12.4)') &
        ' SOLVIT: near-singular matrix at column ', I, &
        ', pivot = ', A(M,I)
      A(M,I) = sign(PIVOT_TOL, A(M,I))
    end if
    PIVOT = 1.0D0 / A(M,I)
    A(M,I) = A(I,I)
    A(I,I) = PIVOT

    ! --- (d) Compute multipliers L(K,I) = A(K,I) / pivot ---
    do K = I + 1, N
      A(K,I) = A(K,I) * PIVOT
    end do

    ! --- (e) Eliminate: update submatrix A(I+1:N, I+1:N) ---
    do J = I + 1, N
      C = A(I,J)
      if (C == 0.0D0) cycle
      do K = I + 1, N
        A(K,J) = A(K,J) - A(K,I) * C
      end do
    end do

  end do

  ! Reciprocal of the last pivot
  if (abs(A(N,N)) < PIVOT_TOL) then
    write(6,'(A,I4,A,1PE12.4)') &
      ' SOLVIT: near-singular matrix at column ', N, &
      ', pivot = ', A(N,N)
    A(N,N) = sign(PIVOT_TOL, A(N,N))
  end if
  A(N,N) = 1.0D0 / A(N,N)

  ! =================================================================
  ! Phase 2: Forward substitution — apply row swaps and L^{-1} to B
  ! =================================================================

  do I = 1, N - 1
    M = IPIVOT(I)
    if (M /= I) then
      T = B(M)
      B(M) = B(I)
      B(I) = T
    end if
    C = B(I)
    do K = I + 1, N
      B(K) = B(K) - A(K,I) * C
    end do
  end do

  ! =================================================================
  ! Phase 3: Back substitution — apply U^{-1} to B
  ! =================================================================
  !   A(J,J) stores 1/U(J,J); A(K,J) for K<J stores U(K,J).

  do J = N, 2, -1
    B(J) = B(J) * A(J,J)
    C = B(J)
    do K = 1, J - 1
      B(K) = B(K) - A(K,J) * C
    end do
  end do
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

  implicit none

  character(*), intent(in) :: LABEL
  real*8,       intent(in) :: ARR(*)
  integer,      intent(in) :: N

  integer :: I

  write(6, '(A6,1P10E12.4 / (7X,10E12.4))') LABEL, (ARR(I), I=1,N)

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

  implicit none

  real*8, intent(in) :: TEMP, PRESSURE, V
  real*8             :: ROSSTAB

  ! --- Table storage (persistent across calls) ---
  integer, parameter :: MAXTAB = kw * 60

  real*8,  save :: ROSS(MAXTAB)     ! log10(kappa_Ross)
  real*8,  save :: TABT(MAXTAB)     ! normalized log10(T)
  real*8,  save :: TABP(MAXTAB)     ! normalized log10(P)
  real*8,  save :: ZEROT, ZEROP     ! normalization offsets
  real*8,  save :: SLOPET, SLOPEP   ! normalization ranges
  integer, save :: NROSS = 0        ! number of entries in table

  ! --- Shepard parameters ---
  integer, parameter :: KFIT = 12          ! nearest neighbors
  real*8,  parameter :: SHEP_PHALF = 1.5D0 ! p/2 where p=3
  real*8,  parameter :: SHEP_EPS = 1.0D-12 ! softening

  ! --- Local variables ---
  real*8  :: TEMPLOG, PRESSLOG      ! normalized query coordinates
  real*8  :: DT, DP, RADIUS2        ! displacement and distance
  real*8  :: R, RWT, W              ! interpolated result, weights

  ! Bilinear variables
  real*8  :: RMIN_PP, RMIN_PM, RMIN_MP, RMIN_MM
  integer :: IDX_PP,  IDX_PM,  IDX_MP,  IDX_MM
  real*8  :: DIST_PP, DIST_PM, DIST_MP, DIST_MM
  real*8  :: TPP, PPP, RPP, TPM, PPM, RPM
  real*8  :: TMP, PMP, RMP, TMM, PMM, RMM
  real*8  :: R_HIPRESS, R_LOPRESS
  real*8  :: P_HIPRESS, P_LOPRESS

  ! Shepard variables
  real*8  :: DIST2(KFIT)
  integer :: KIDX(KFIT)
  real*8  :: DMAX2
  integer :: KMAX, KFOUND

  integer :: I, J, K

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING ROSSTAB'

  ! =================================================================
  ! MODE 1: STORE — append current kappa_Ross profile to table
  ! =================================================================

  if (TEMP <= 0.0D0) then

    if (NROSS == 0) then
      ZEROT  = log10(T(1))
      ZEROP  = log10(P(1))
      SLOPET = log10(T(NRHOX)) - ZEROT
      SLOPEP = log10(P(NRHOX)) - ZEROP
    end if

    do J = 1, NRHOX
      NROSS = NROSS + 1
      if (NROSS > MAXTAB) then
        write(6,*) 'ROSSTAB ERROR: table overflow, NROSS =', NROSS
        NROSS = MAXTAB
        exit
      end if
      TABT(NROSS) = (log10(T(J)) - ZEROT) / SLOPET
      TABP(NROSS) = (log10(P(J)) - ZEROP) / SLOPEP
      ROSS(NROSS) = log10(ABROSS(J))
      if (IDEBUG == 1) &
        write(6, '(" ROSSTAB",I5,F10.1,F10.5,1PE12.3,0PF10.5,F10.5)') &
        NROSS, T(J), TABT(NROSS), P(J), TABP(NROSS), ROSS(NROSS)
    end do

    ROSSTAB = 0.0D0
    return
  end if

  ! =================================================================
  ! MODE 2: LOOKUP
  ! =================================================================

  TEMPLOG  = (log10(TEMP) - ZEROT) / SLOPET
  PRESSLOG = (log10(PRESSURE) - ZEROP) / SLOPEP

  ! -----------------------------------------------------------------
  ! IROSSTAB = 1: Original bilinear (4-quadrant nearest neighbor)
  ! -----------------------------------------------------------------
  if (IROSSTAB == 1) then

    RMIN_PP = 1.0D30;  IDX_PP = 0
    RMIN_PM = 1.0D30;  IDX_PM = 0
    RMIN_MP = 1.0D30;  IDX_MP = 0
    RMIN_MM = 1.0D30;  IDX_MM = 0

    do I = 1, NROSS
      DT = TABT(I) - TEMPLOG
      DP = TABP(I) - PRESSLOG
      RADIUS2 = DT**2 + DP**2

      if (DT >= 0.0D0) then
        if (DP >= 0.0D0) then
          if (RADIUS2 < RMIN_PP) then
            RMIN_PP = RADIUS2;  IDX_PP = I
          end if
        else
          if (RADIUS2 < RMIN_PM) then
            RMIN_PM = RADIUS2;  IDX_PM = I
          end if
        end if
      else
        if (DP >= 0.0D0) then
          if (RADIUS2 < RMIN_MP) then
            RMIN_MP = RADIUS2;  IDX_MP = I
          end if
        else
          if (RADIUS2 < RMIN_MM) then
            RMIN_MM = RADIUS2;  IDX_MM = I
          end if
        end if
      end if
    end do

    if (IDX_PP /= 0 .and. IDX_PM /= 0 .and. &
        IDX_MP /= 0 .and. IDX_MM /= 0) then

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

    else

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

    end if

  ! -----------------------------------------------------------------
  ! IROSSTAB = 2: Shepard (K-nearest, p=3)
  ! -----------------------------------------------------------------
  else

    KFOUND = 0
    DMAX2  = 0.0D0
    KMAX   = 1

    do I = 1, NROSS
      DT = TABT(I) - TEMPLOG
      DP = TABP(I) - PRESSLOG
      RADIUS2 = DT * DT + DP * DP

      if (KFOUND < KFIT) then
        KFOUND = KFOUND + 1
        DIST2(KFOUND) = RADIUS2
        KIDX(KFOUND)  = I
        if (RADIUS2 > DMAX2) then
          DMAX2 = RADIUS2
          KMAX  = KFOUND
        end if
      else if (RADIUS2 < DMAX2) then
        DIST2(KMAX) = RADIUS2
        KIDX(KMAX)  = I
        DMAX2 = DIST2(1)
        KMAX  = 1
        do K = 2, KFIT
          if (DIST2(K) > DMAX2) then
            DMAX2 = DIST2(K)
            KMAX  = K
          end if
        end do
      end if
    end do

    R   = 0.0D0
    RWT = 0.0D0
    do K = 1, KFOUND
      W   = 1.0D0 / (DIST2(K) + SHEP_EPS) ** SHEP_PHALF
      R   = R   + W * ROSS(KIDX(K))
      RWT = RWT + W
    end do
    R = R / RWT

  end if

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

  implicit none

  ! --- Arguments ---
  integer, intent(in)    :: NUMTAU
  real*8,  intent(in)    :: T_in(NUMTAU)
  real*8,  intent(in)    :: TAU(NUMTAU)
  real*8,  intent(out)   :: ABSTD(NUMTAU)
  real*8,  intent(out)   :: PTOTAL(NUMTAU)
  real*8,  intent(out)   :: P_out(NUMTAU)
  real*8,  intent(in)    :: PRAD_in(NUMTAU)
  real*8,  intent(in)    :: PTURB_in(NUMTAU)
  real*8,  intent(in)    :: VTURB_in(NUMTAU)
  real*8,  intent(in)    :: GRAV_in

  ! --- Tuning parameters ---
  integer, parameter :: MAX_CORRECTOR  = 1000
  real*8,  parameter :: CONV_TOL       = 5.0D-5   ! convergence in ln(P)
  real*8,  parameter :: DPLOG_MAX      = 10.0D0   ! max d(ln P) per step
  real*8,  parameter :: OPACITY_REBOOT = 2.0D0    ! opacity ratio triggering surface redo
  integer, parameter :: RELAX_TRIGGER  = 4        ! stall count before under-relaxation
  real*8,  parameter :: RELAX_FACTOR   = 0.3D0    ! relaxation weight when oscillating

  ! --- Local variables ---
  real*8  :: DLGTAU              ! log spacing: ln(TAU(2)/TAU(1))
  real*8  :: PLOG                ! current ln(P_total) estimate
  real*8  :: PNEW               ! corrected ln(P_total)
  real*8  :: DPLOG              ! d(ln P)/d(ln tau) * DLGTAU at current point
  real*8  :: ERROR              ! |PNEW - PLOG| convergence measure
  real*8  :: PREV_ERROR         ! previous corrector error
  real*8  :: ABSTD1_INIT        ! initial surface opacity guess (for bootstrap)
  real*8  :: RELAX              ! relaxation weight for current step
  integer :: N_STALL            ! count of stalled corrector steps

  ! Adams-Bashforth history (previous 4 points)
  real*8  :: PLOG1, PLOG2, PLOG3, PLOG4
  real*8  :: DPLOG1, DPLOG2, DPLOG3

  integer :: J, N
  logical :: CONVERGED

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING TTAUP'

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
  if (PRAD_in(1) > 0.0D0) then
    ABSTD(1) = min(0.1D0, GRAV_in * TAU(1) / PRAD_in(1) / 2.0D0)
  end if
  ABSTD1_INIT = ABSTD(1)

  !---------------------------------------------------------------------
  ! March inward through the atmosphere
  !---------------------------------------------------------------------
  do J = 1, NUMTAU

    ! --- Predictor: initial guess for ln(P_total) ---
    if (J == 1) then
      ! Surface: P_total = g * tau / kappa (constant-opacity approximation)
      PLOG = log(GRAV_in / ABSTD(1) * TAU(1))
    else if (J <= 4) then
      ! Linear extrapolation from previous point
      PLOG = PLOG1 + DPLOG1
    else
      ! 4th-order Adams-Bashforth predictor
      PLOG = (3.0D0 * PLOG4 + 8.0D0 * DPLOG1 &
            - 4.0D0 * DPLOG2 + 8.0D0 * DPLOG3) / 3.0D0
    end if

    ! --- Corrector iterations ---
    ERROR      = 1.0D0
    PREV_ERROR = 1.0D0
    N          = 0
    CONVERGED  = .false.
    N_STALL    = 0

    do
      ! Convert ln(P_total) → total pressure → gas pressure → opacity
      PTOTAL(J) = exp(PLOG)

      ! Gas pressure: subtract radiation and turbulent pressure changes
      ! relative to their surface values
      P_out(J) = PTOTAL(J) + (PRAD_in(1) - PRAD_in(J)) &
                            + (PTURB_in(1) - PTURB_in(J))

      if (P_out(J) <= 0.0D0) then
        ! Gas pressure went negative — corrector guess for P_total is too low.
        ! Floor P_gas to a small fraction of P_total so the opacity lookup
        ! returns a large kappa, which drives DPLOG upward and self-corrects.
        P_out(J) = PTOTAL(J) * 1.0D-4
        if (P_out(J) <= 0.0D0) P_out(J) = 1.0D-10
      end if

      ! Look up Rosseland opacity from the interpolation table
      ABSTD(J) = ROSSTAB(T_in(J), P_out(J), VTURB_in(J))

      ! Hydrostatic equilibrium derivative
      DPLOG = GRAV_in / ABSTD(J) * TAU(J) / PTOTAL(J) * DLGTAU

      ! Clamp to prevent unphysical pressure jumps
      DPLOG = max(0.0D0, min(DPLOG_MAX, DPLOG))

      N = N + 1
      if (N == 1) cycle    ! first pass: go straight to corrector

      ! --- Corrector: update ln(P) ---
      if (J == 1) then
        PNEW = log(GRAV_in / ABSTD(1) * TAU(1))
      else if (J <= 4) then
        ! Trapezoidal corrector
        PNEW = (PLOG + 2.0D0 * PLOG1 + DPLOG + DPLOG1) / 3.0D0
      else
        ! 4th-order Adams-Moulton corrector
        PNEW = (126.0D0 * PLOG1 - 14.0D0 * PLOG3 + 9.0D0 * PLOG4 &
              + 42.0D0 * DPLOG + 108.0D0 * DPLOG1 &
              - 54.0D0 * DPLOG2 + 24.0D0 * DPLOG3) / 121.0D0
      end if

      PREV_ERROR = ERROR
      ERROR = abs(PNEW - PLOG)

      ! Detect stalling: error not decreasing meaningfully
      if (N > 3 .and. ERROR > 0.9D0 * PREV_ERROR) then
        N_STALL = N_STALL + 1
      else
        N_STALL = max(N_STALL - 1, 0)
      end if

      ! Choose relaxation weight:
      !   Normal: average old and new (0.5), same as original code
      !   Oscillating: bias toward old value to damp oscillation
      if (N_STALL >= RELAX_TRIGGER) then
        RELAX = RELAX_FACTOR
      else
        RELAX = 0.5D0
      end if
      PLOG = (1.0D0 - RELAX) * PLOG + RELAX * PNEW

      if (ERROR <= CONV_TOL) then
        CONVERGED = .true.
        exit
      end if
      if (N > MAX_CORRECTOR) exit
    end do

    ! Warn on convergence failure
    if (.not. CONVERGED) then
      write(6, '(A,I4,A,ES10.3,A,I5)') &
        ' TTAUP WARNING: corrector did not converge at J=', J, &
        ', error=', ERROR, ', iter=', N
    end if

    ! --- Surface bootstrap ---
    ! If the converged opacity at J=1 differs substantially from the
    ! initial blind guess, redo J=1 with the better opacity value.
    ! This prevents a poor surface guess from corrupting the entire
    ! Adams-Bashforth history.
    if (J == 1 .and. (ABSTD(1) / ABSTD1_INIT > OPACITY_REBOOT .or. &
                       ABSTD1_INIT / ABSTD(1) > OPACITY_REBOOT)) then
      ABSTD1_INIT = ABSTD(1)   ! prevent infinite re-bootstrapping
      PLOG = log(GRAV_in / ABSTD(1) * TAU(1))
      PTOTAL(1) = exp(PLOG)
      P_out(1)  = PTOTAL(1)    ! at surface: PRAD(1)-PRAD(1) = 0, same for PTURB
      if (P_out(1) > 0.0D0) then
        ABSTD(1) = ROSSTAB(T_in(1), P_out(1), VTURB_in(1))
        PLOG     = log(GRAV_in / ABSTD(1) * TAU(1))
        PTOTAL(1) = exp(PLOG)
        P_out(1)  = PTOTAL(1)
      end if
      DPLOG = GRAV_in / ABSTD(1) * TAU(1) / PTOTAL(1) * DLGTAU
      DPLOG = max(0.0D0, min(DPLOG_MAX, DPLOG))
    end if

    ! --- Shift Adams-Bashforth history ---
    PLOG4  = PLOG3
    PLOG3  = PLOG2
    PLOG2  = PLOG1
    PLOG1  = PLOG
    DPLOG3 = DPLOG2
    DPLOG2 = DPLOG1
    DPLOG1 = DPLOG

  end do

  return

END SUBROUTINE TTAUP



!=======================================================================
! READIN: Read and parse all input data
!=======================================================================

SUBROUTINE READIN(MODE)

  implicit none

  integer, intent(in) :: MODE
!     MODE=1  COMPUTE A MODEL (ATLAS12: reads input_model.dat + stdin)
!     MODE=20 READ A MODEL FOR SYNTHESIS (SYNTHE: caller opens file on unit 5)

  ! Local scalars
  real*8  :: XSCALELOG
  integer :: I, IZ
  integer :: J, MU
  integer :: IOS_CARD
  logical :: FIRST_KW
  character(20) :: KEYWORD

  character(1) :: CARD(132)

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING READIN'

  if (MODE == 1) then
    !-------------------------------------------------------------------
    ! ATLAS12 path: read model from input_model.dat, then stdin overrides
    !-------------------------------------------------------------------
    OPEN(UNIT=3, FILE='input_model.dat', STATUS='OLD', ACTION='READ')
    INPUTDATA = 3
    IFPRES = 1
    IFCORR = 1
    call READMOL
  else
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
    call READMOL
  end if

  LAST=133
  MAXPOW=38

  !---------------------------------------------------------------------
  ! Main card-reading loop
  !---------------------------------------------------------------------
  card_loop: do
    MORE = 0
    LETCOL = 1
    READ(INPUTDATA, '(132A1)', IOSTAT=IOS_CARD) CARD
    if (IOS_CARD /= 0) exit card_loop

    !-------------------------------------------------------------------
    ! Keyword parsing loop (multiple keywords per card)
    !-------------------------------------------------------------------
    KEYWORD = NEXTWORD(CARD)
    NUMCOL = LETCOL
    FIRST_KW = .true.

    keyword_loop: do
      if (.not. FIRST_KW) then
        ! Get next keyword from same card
        LETCOL = MAX(LETCOL, NUMCOL)
        MORE = 1
        KEYWORD = NEXTWORD(CARD)
        if (IFFAIL == 1) exit keyword_loop
        MORE = 0
        NUMCOL = LETCOL
      end if
      FIRST_KW = .false.

      !-----------------------------------------------------------------
      ! TEFF
      !-----------------------------------------------------------------
      if (trim(KEYWORD) == 'TEFF') then
        TEFF = FREEFF(CARD)
        FLUX = SIGMA_SB / FOURPI * TEFF**4
        cycle keyword_loop

      !-----------------------------------------------------------------
      ! GRAVITY
      !-----------------------------------------------------------------
      else if (trim(KEYWORD) == 'GRAVITY') then
        GRAV = FREEFF(CARD)
        if (GRAV < 10.0D0) GRAV = 10.0D0**(GRAV)
        GLOG = LOG10(GRAV)
        cycle keyword_loop

      !-----------------------------------------------------------------
      ! ABUNDANCE (sub-keywords: SCALE, CHANGE, ABSOLUTE, RELATIVE, TABLE)
      !-----------------------------------------------------------------
      else if (trim(KEYWORD) == 'ABUNDANCE') then
        KEYWORD = NEXTWORD(CARD)
        ! SCALE
        if (trim(KEYWORD) == 'SCALE') then
          NUMCOL = LETCOL
          XSCALE = FREEFF(CARD)
          if (XSCALE > 0.0D0) XSCALELOG = LOG10(XSCALE)
          do IZ = 3, 99
            XRELATIVE(IZ) = XSCALELOG
          end do
          XSCALE = 1.0D0
          cycle keyword_loop
        ! CHANGE
        else if (trim(KEYWORD) == 'CHANGE') then
          MORE = 1
          do
            IZ = FREEFF(CARD)
            if (IFFAIL == 1) exit keyword_loop
            ABUND(IZ) = FREEFF(CARD)
            if (IZ > 2 .and. ABUND(IZ) > 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
          end do
        ! ABSOLUTE
        else if (trim(KEYWORD) == 'ABSOLUTE') then
          MORE = 1
          do
            IZ = FREEFF(CARD)
            if (IFFAIL == 1) exit keyword_loop
            ABUND(IZ) = FREEFF(CARD)
            if (IZ > 2 .and. ABUND(IZ) > 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
            XRELATIVE(IZ) = 0.0D0
          end do
        ! RELATIVE
        else if (trim(KEYWORD) == 'RELATIVE') then
          MORE = 1
          do
            IZ = FREEFF(CARD)
            if (IFFAIL == 1) exit keyword_loop
            XRELATIVE(IZ) = FREEFF(CARD)
          end do
        ! TABLE
        else if (trim(KEYWORD) == 'TABLE') then
          READ(INPUTDATA, '(7X,F10.6,10X,F10.6/(5(7X,F7.3,F6.3)))') &
            ABUND(1), ABUND(2), (ABUND(IZ), XRELATIVE(IZ), IZ=3,99)
          do IZ = 3, 99
            if (ABUND(IZ) > 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
          end do
          XSCALE = 1.0D0
          exit keyword_loop
        else
          write(6, '(A,A)') ' ABUNDANCE: unknown sub-keyword: ', trim(KEYWORD)
          stop 1
        end if

      !-----------------------------------------------------------------
      ! READ (sub-keyword: DECK/DECK6)
      !-----------------------------------------------------------------
      else if (trim(KEYWORD) == 'READ') then
        KEYWORD = NEXTWORD(CARD)
        NUMCOL = LETCOL

        ! DECK / DECK6
        if (trim(KEYWORD) == 'DECK' .or. trim(KEYWORD) == 'DECK6') then
          NRHOX = FREEFF(CARD)
          do J = 1, NRHOX
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
          end do
          if (RHOX(1) < 0.0D0) then
            do J = 1, NRHOX
              RHOX(J) = 10.0D0**(RHOX(J))
            end do
          end if
          PRADK0 = 0.0D0
          PTURB0 = PTURB(1)
          PCON = 0.0D0
          PZERO = PCON + PRADK0 + PTURB0
          call INTEG(RHOX, ABROSS, TAUROS, NRHOX, ABROSS(1)*RHOX(1))
          if (trim(KEYWORD) == 'DECK6') then
            ! DECK6: read additional PRADK0 card
            READ(INPUTDATA, '(132A1)') CARD
            NUMCOL = 1
            PRADK0 = FREEFF(CARD)
            do J = 1, NRHOX
              ACCRAD(J) = PRAD(J)
            end do
            call INTEG(RHOX, ACCRAD, PRAD, NRHOX, ACCRAD(1)*RHOX(1))
            do J = 1, NRHOX
              PRADK(J) = PRAD(J) + PRADK0
            end do
          end if
          exit keyword_loop

        else
          write(6, '(A,A)') ' READ: unknown sub-keyword: ', trim(KEYWORD)
          stop 1
        end if

      !-----------------------------------------------------------------
      ! BEGIN — end of model file; exit card loop for finalization
      !-----------------------------------------------------------------
      else if (trim(KEYWORD) == 'BEGIN') then
        if (INPUTDATA == 3) CLOSE(UNIT=3)
        exit card_loop

      !-----------------------------------------------------------------
      ! TURBULENCE (sub-keywords: ON, OFF)
      !-----------------------------------------------------------------
      else if (trim(KEYWORD) == 'TURBULENCE') then
        KEYWORD = NEXTWORD(CARD)
        if (trim(KEYWORD) == 'ON') then
          IFTURB = 1
          NUMCOL = LETCOL
          TRBFDG = FREEFF(CARD)
          TRBPOW = FREEFF(CARD)
          TRBSND = FREEFF(CARD)
          TRBCON = FREEFF(CARD)
        else if (trim(KEYWORD) == 'OFF') then
          IFTURB = 0
          TRBFDG = 0.0D0
          TRBPOW = 0.0D0
          TRBSND = 0.0D0
          TRBCON = 0.0D0
        else
          write(6, '(A,A)') ' TURBULENCE: unknown sub-keyword: ', trim(KEYWORD)
          stop 1
        end if
        cycle keyword_loop

      !-----------------------------------------------------------------
      ! SURFACE (sub-keywords: INTENSITY, FLUX, OFF)
      !-----------------------------------------------------------------
      else if (trim(KEYWORD) == 'SURFACE') then
        KEYWORD = NEXTWORD(CARD)
        if (trim(KEYWORD) == 'INTENSITY') then
          NMU = FREEFF(CARD)
          do MU = 1, NMU
            ANGLE(MU) = FREEFR(CARD)
          end do
          IFSURF = 2
        else if (trim(KEYWORD) == 'FLUX') then
          IFSURF = 1
        else if (trim(KEYWORD) == 'OFF') then
          IFSURF = 0
        else
          write(6, '(A,A)') ' SURFACE: unknown sub-keyword: ', trim(KEYWORD)
          stop 1
        end if
        cycle keyword_loop

      !-----------------------------------------------------------------
      ! Unrecognized keyword — skip to next card
      !-----------------------------------------------------------------
      else
        exit keyword_loop
      end if

    end do keyword_loop
  end do card_loop

  ! ===================================================================
  ! FINALIZATION — compute derived quantities
  ! ===================================================================
  if (ABUND(1) < 0.0D0) ABUND(1) = 10.0D0**ABUND(1)
  if (ABUND(2) < 0.0D0) ABUND(2) = 10.0D0**ABUND(2)
  do IZ = 3, 99
    if (ABUND(IZ) > 0.0D0) ABUND(IZ) = LOG10(ABUND(IZ))
  end do
  ! Write abbreviated list of abundances
  WRITE(6, '(/" TEFF",F7.0,"   LOGG",F8.4/' // &
       '"   1",A2,F10.6,"     2",A2,F10.6/(5(I4,A2,F7.2,F5.2)))') &
       TEFF, GLOG, ELEM(1), ABUND(1), ELEM(2), ABUND(2), &
     (IZ, ELEM(IZ), ABUND(IZ), XRELATIVE(IZ), IZ=3,32)
  do J = 1, NRHOX
    do IZ = 3, 99
      XABUND(J,IZ) = 10.0D0**(ABUND(IZ) + XRELATIVE(IZ))
    end do
    XABUND(J,1) = ABUND(1)
    XABUND(J,2) = ABUND(2)
    WTMOLE(J) = 0.0D0
    do IZ = 1, 99
      WTMOLE(J) = WTMOLE(J) + XABUND(J,IZ) * ATMASS(IZ)
    end do
  end do
  YABUND(1) = ABUND(1)
  YABUND(2) = ABUND(2)
  do IZ = 3, 99
    YABUND(IZ) = ABUND(IZ) + XRELATIVE(IZ)
  end do
  do J = 1, NRHOX
    TK(J) = KBOL * T(J)
    HKT(J) = HPLANCK / TK(J)
    HCKT(J) = HKT(J) * CLIGHT
    TKEV(J) = 8.6171D-5 * T(J)
    TLOG(J) = LOG(T(J))
    XNATOM(J) = P(J) / TK(J) - XNE(J)
    RHO(J) = XNATOM(J) * WTMOLE(J) * AMU
    if (IFTURB > 0) PTURB(J) = 0.5D0 * RHO(J) * VTURB(J)**2
    CHARGESQ(J) = XNE(J)
  end do
  WRITE(6, '(/" H1",I2," H2PLUS",I2," HMINUS",I2," HRAY",I2,' // &
    '" HE1",I2," HE2",I2," HEMINUS",I2," HERAY",I2," COOL",I2,' // &
    '" LUKE",I2/" HOT",I2," ELECTRON",I2," H2RAY",I2," HLINES",' // &
    'I2," LINES",I2," LINESCAT",I2," XLINES",I2," XLSCAT",I2,' // &
    '" XCONT",I2," XSCAT",I2)') IFOP
  WRITE(6, '(" IFCORR",I2,"  IFPRES",I2,"  IFSURF",I2,' // &
    '"  IFSCAT",I2,"  IFCONV",I2,"  MIXLTH",F6.2,"  IFMOL",I2/' // &
    '" IFTURB",I2,"  TRBFDG",F6.2,"  TRBPOW",F6.2,' // &
    '"  TRBSND",F6.2,"  TRBCON",F6.2)') &
    IFCORR, IFPRES, IFSURF, IFSCAT, IFCONV, MIXLTH, IFMOL, &
    IFTURB, TRBFDG, TRBPOW, TRBSND, TRBCON
  write(6, *)
  WRITE(6, '(" NUMITS",I3,"  IFPRNT",*(I2))') NUMITS, IFPRNT(1:NUMITS)
  WRITE(6, '(10X,"  IFPNCH",*(I2))') IFPNCH(1:NUMITS)

  ! Auto-generate title from mixing length
  block
    character(74) :: TITLEBUF
    write(TITLEBUF, '(A,F5.2)') 'ATLAS12 l/H=', MIXLTH
    do I = 1, 74
      TITLE(I) = TITLEBUF(I:I)
    end do
  end block

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

  implicit none

  real*8, intent(in) :: TEFF_NEW, LOGG_NEW

  ! --- Tau grid parameters ---
  ! NDEPTHS: number of depth points in the output model
  ! TAU1LG and STEPLG are module-level variables (mod_atlas_data),
  ! also used by TCORR to rebuild the grid during iterations.
  ! The deepest point is at log10(tau) = TAU1LG + (NDEPTHS-1)*STEPLG
  integer, parameter :: NDEPTHS = 80

  ! Local work arrays for interpolation
  real*8  :: DUM1(kw), DUM2(kw), DUM3(kw), DUM4(kw)
  real*8  :: DUM5(kw), DUM6(kw), DUM7(kw), DUM8(kw)
  real*8  :: GNEW
  integer :: I, J, IDUM

  ! Convert logg to linear gravity
  GNEW = 10.0D0**LOGG_NEW

  ! --- Step 1: Build new uniform log-tau grid ---
  do J = 1, NDEPTHS
    TAUSTD(J) = 10.0D0**(TAU1LG + (J-1)*STEPLG)
  end do

  ! --- Step 2: Integrate ABROSS to get TAUROS from current model ---
  call INTEG(RHOX, ABROSS, TAUROS, NRHOX, ABROSS(1)*RHOX(1))

  ! --- Step 3: Interpolate all structure variables onto new grid ---
  IDUM = MAP1(TAUROS, RHOX,   NRHOX, TAUSTD, DUM1, NDEPTHS)
  IDUM = MAP1(TAUROS, T,      NRHOX, TAUSTD, DUM2, NDEPTHS)
  IDUM = MAP1(TAUROS, P,      NRHOX, TAUSTD, DUM3, NDEPTHS)
  IDUM = MAP1(TAUROS, XNE,    NRHOX, TAUSTD, DUM4, NDEPTHS)
  IDUM = MAP1(TAUROS, ABROSS, NRHOX, TAUSTD, DUM5, NDEPTHS)
  IDUM = MAP1(TAUROS, PRAD,   NRHOX, TAUSTD, DUM6, NDEPTHS)
  IDUM = MAP1(TAUROS, VTURB,  NRHOX, TAUSTD, DUM7, NDEPTHS)
  IDUM = MAP1(TAUROS, BMIN,   NRHOX, TAUSTD, DUM8, NDEPTHS)
  do J = 1, NDEPTHS
    RHOX(J)   = DUM1(J)
    T(J)      = DUM2(J)
    P(J)      = DUM3(J)
    XNE(J)    = DUM4(J)
    ABROSS(J) = DUM5(J)
    PRAD(J)   = DUM6(J)
    PRADK(J)  = PRAD(J) + PRADK0
    VTURB(J)  = DUM7(J)
    BMIN(J)   = DUM8(J)
  end do
  do I = 1, 6
    IDUM = MAP1(TAUROS, BHYD(1,I), NRHOX, TAUSTD, DUM1, NDEPTHS)
    do J = 1, NDEPTHS
      BHYD(J,I) = DUM1(J)
    end do
  end do
  NRHOX = NDEPTHS

  ! --- Step 4: Rescale to new Teff/logg if different ---
  if (TEFF_NEW == 0.0D0) return
  if (TEFF_NEW == TEFF .and. GNEW == GRAV) return
  if (TEFF_NEW < TEFF+1.0D0 .and. TEFF_NEW > TEFF-1.0D0 .and. &
      GNEW < GRAV*1.001D0 .and. GNEW > GRAV*0.999D0) return

  do J = 1, NRHOX
    TAUROS(J) = TAUSTD(J)
    T(J)      = T(J) * TEFF_NEW / TEFF
    PTURB(J)  = 0.0D0
    PRADK(J)  = PRADK(J) * (TEFF_NEW/TEFF)**4
    PRAD(J)   = PRAD(J)  * (TEFF_NEW/TEFF)**4
  end do
  PRADK0 = PRADK0 * (TEFF_NEW/TEFF)**4
  PZERO  = PCON + PRADK0 + PTURB0
  TEFF   = TEFF_NEW
  FLUX   = SIGMA_SB / FOURPI * TEFF**4
  GRAV   = GNEW
  GLOG   = LOG10(GRAV)
  do J = 1, NRHOX
    PTOTAL(J) = P(J) + PRAD(J) + PTURB(J)
    RHOX(J)   = PTOTAL(J) / GRAV
    PTOTAL(J) = PTOTAL(J) + PZERO
  end do

END SUBROUTINE SCALE_MODEL

!=======================================================================
! FREEFR: Free-format real number reader
!=======================================================================

FUNCTION FREEFR(CARD)

  implicit none
  character(1), intent(inout) :: CARD(*)
  real*8 :: FREEFR
  integer :: I, L

  MORE = 1
  FREEFR = FREEFF(CARD)
  if (IFFAIL == 0) return
  L = LAST - 1
  read(INPUTDATA, 1) (CARD(I), I=1, L)
    1 format(132A1)
  NUMCOL = 1
  FREEFR = FREEFF(CARD)
  return

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

  implicit none

  character(1), intent(in) :: CARD(*)
  real*8 :: FREEFF

  ! Digit lookup table: 0-9
  character(1), parameter :: DIGIT(10) = (/ &
    '0','1','2','3','4','5','6','7','8','9' /)

  ! Parser states
  integer, parameter :: ST_INTEGER  = 1  ! before decimal point
  integer, parameter :: ST_FRACTION = 2  ! after decimal point
  integer, parameter :: ST_EXPSIGN  = 3  ! after E (expecting sign or digit)
  integer, parameter :: ST_EXPDIGIT = 4  ! exponent digits

  character(1) :: C
  real*8  :: ANSWER, ASIGN
  integer :: I, N, NCOL, NPT, NPOWER, ISIGN, IF0, STATE

  ! Initialize
  IFFAIL = 0
  if (NUMCOL > LAST) then
    IFFAIL = 1
    if (MORE > 0) then
      FREEFF = 0.0D0
      return
    end if
    write(6, '(A/(1X,131A1))') '1FREEFF HAS READ OFF THE END', &
      (CARD(I), I = 1, LAST)
    stop 1
  end if

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
  scan_loop: do NCOL = NUMCOL, LAST
    C = CARD(NCOL)

    select case (STATE)

    !-----------------------------------------------------------------
    ! State 1: Integer part (before decimal point)
    !-----------------------------------------------------------------
    case (ST_INTEGER)
      if (C == ' ') then
        if (N == 0) then
          call freeff_reset()
          cycle scan_loop
        end if
        ! Blank after digits — return integer value
        FREEFF = ANSWER * ASIGN
        NUMCOL = NCOL + 1
        return
      end if
      ! Check for digit
      do I = 1, 10
        if (C == DIGIT(I)) then
          N = N + 1
          ANSWER = 10.0D0 * ANSWER + dble(I - 1)
          cycle scan_loop
        end if
      end do
      ! Check for decimal point
      if (C == '.') then
        STATE = ST_FRACTION
        cycle scan_loop
      end if
      ! Check for comma (delimiter)
      if (C == ',') then
        if (N == 0) then
          ! Reset — comma before any digits
          call freeff_reset()
          cycle scan_loop
        end if
        FREEFF = ANSWER * ASIGN
        NUMCOL = NCOL + 1
        return
      end if
      ! Check for minus sign
      if (C == '-') then
        if (N == 0) then
          ASIGN = -1.0D0
          cycle scan_loop
        end if
        ! Minus after digits — reset
        call freeff_reset()
        cycle scan_loop
      end if
      ! Unrecognized character — reset
      call freeff_reset()

    !-----------------------------------------------------------------
    ! State 2: Fractional part (after decimal point)
    !-----------------------------------------------------------------
    case (ST_FRACTION)
      ! Check for digit
      do I = 1, 10
        if (C == DIGIT(I)) then
          N = N + 1
          NPT = NPT + 1
          ANSWER = 10.0D0 * ANSWER + dble(I - 1)
          cycle scan_loop
        end if
      end do
      if (C == 'E') then
        STATE = ST_EXPSIGN
        cycle scan_loop
      end if
      if (C == '-') then
        ! Minus after decimal digits — treat as exponent sign
        ISIGN = -1
        NPOWER = 0
        STATE = ST_EXPDIGIT
        cycle scan_loop
      end if
      if (C == '+') then
        NPOWER = 0
        STATE = ST_EXPDIGIT
        cycle scan_loop
      end if
      if (C == ' ' .or. C == ',') then
        if (N == 0) then
          call freeff_reset()
          cycle scan_loop
        end if
        FREEFF = ANSWER * ASIGN / 10.0D0**NPT
        NUMCOL = NCOL + 1
        return
      end if
      ! Unrecognized — reset
      call freeff_reset()

    !-----------------------------------------------------------------
    ! State 3: After E (expecting sign or first exponent digit)
    !-----------------------------------------------------------------
    case (ST_EXPSIGN)
      ! Check for digit
      do I = 1, 10
        if (C == DIGIT(I)) then
          NPOWER = I - 1
          IF0 = 1
          STATE = ST_EXPDIGIT
          cycle scan_loop
        end if
      end do
      if (C == ' ' .or. C == '+') then
        NPOWER = 0
        STATE = ST_EXPDIGIT
        cycle scan_loop
      end if
      if (C == '-') then
        ISIGN = -1
        NPOWER = 0
        STATE = ST_EXPDIGIT
        cycle scan_loop
      end if
      ! Unrecognized — reset
      call freeff_reset()

    !-----------------------------------------------------------------
    ! State 4: Exponent digits
    !-----------------------------------------------------------------
    case (ST_EXPDIGIT)
      ! Check for digit
      do I = 1, 10
        if (C == DIGIT(I)) then
          NPOWER = 10 * NPOWER + I - 1
          IF0 = 1
          if (NPOWER >= MAXPOW) then
            call freeff_reset()
            cycle scan_loop
          end if
          cycle scan_loop
        end if
      end do
      if (C == ',' .or. C == ' ') then
        if (IF0 == 0) then
          call freeff_reset()
          cycle scan_loop
        end if
        FREEFF = ANSWER * ASIGN * 10.0D0**(ISIGN * NPOWER - NPT)
        NUMCOL = NCOL + 1
        return
      end if
      ! Unrecognized — reset
      call freeff_reset()

    end select
  end do scan_loop

  !---------------------------------------------------------------------
  ! Fell off end of card
  !---------------------------------------------------------------------
  NUMCOL = LAST + 1
  IFFAIL = 1
  if (MORE > 0) then
    FREEFF = 0.0D0
    return
  end if
  write(6, '(A/(1X,131A1))') '1FREEFF HAS READ OFF THE END', &
    (CARD(I), I = 1, LAST)
  stop 1

contains

  subroutine freeff_reset()
    ! Reset parser to initial searching state
    ASIGN  = 1.0D0
    ANSWER = 0.0D0
    ISIGN  = 1
    NPT    = 0
    IF0    = 0
    N      = 0
    STATE  = ST_INTEGER
  end subroutine freeff_reset

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

  implicit none

  character(1), intent(in) :: CARD(*)
  character(20) :: NEXTWORD

  character(1) :: C
  integer :: I, NCOL, N
  logical :: IN_WORD, SKIP_E, IS_ALPHA

  ! Initialize
  NEXTWORD = ' '
  IFFAIL = 0
  N = 0
  IN_WORD = .false.

  if (LETCOL > LAST) then
    IFFAIL = 1
    if (MORE > 0) return
    write(6, '(A/(1X,131A1))') 'NEXTWORD HAS READ OFF THE END', &
      (CARD(I), I = 1, LAST)
    stop 1
  end if

  !---------------------------------------------------------------------
  ! Main scan loop
  !---------------------------------------------------------------------
  scan_loop: do NCOL = LETCOL, LAST
    C = CARD(NCOL)

    if (.not. IN_WORD) then
      !--- State: searching for word start ---
      if (C == ' ') cycle scan_loop

      ! Check if C is a letter (A-Z)
      IS_ALPHA = (C >= 'A' .and. C <= 'Z')
      if (.not. IS_ALPHA) then
        ! Not a letter — reset and continue searching
        N = 0
        cycle scan_loop
      end if

      ! Special case: is this an "E" in scientific notation?
      if (C == 'E' .and. NCOL > 1) then
        SKIP_E = .false.
        ! Check preceding character: digit or "."
        if ((CARD(NCOL-1) >= '0' .and. CARD(NCOL-1) <= '9') .or. &
            CARD(NCOL-1) == '.') then
          ! Check following character: digit or blank
          if (CARD(NCOL+1) == ' ' .or. &
              (CARD(NCOL+1) >= '0' .and. CARD(NCOL+1) <= '9')) then
            SKIP_E = .true.
          end if
        end if
        if (SKIP_E) then
          N = 0
          cycle scan_loop
        end if
      end if

      ! Start a new word
      N = 1
      NEXTWORD(1:1) = C
      IN_WORD = .true.

    else
      !--- State: inside a word ---
      if (C == ' ' .or. C == '=' .or. C == ',') then
        ! Word delimiter — return the word
        LETCOL = NCOL + 1
        return
      end if

      ! Check if C is alphanumeric (A-Z or 0-9)
      IS_ALPHA = (C >= 'A' .and. C <= 'Z') .or. &
                 (C >= '0' .and. C <= '9')
      if (.not. IS_ALPHA) then
        ! Not alphanumeric — abandon word, reset to searching
        NEXTWORD = ' '
        N = 0
        IN_WORD = .false.
        cycle scan_loop
      end if

      ! Accumulate character into word (max 20 chars)
      N = N + 1
      if (N <= 20) then
        NEXTWORD(N:N) = C
      end if
    end if
  end do scan_loop

  !---------------------------------------------------------------------
  ! Fell off end of card without completing a word
  !---------------------------------------------------------------------
  LETCOL = LAST + 1
  IFFAIL = 1
  if (MORE > 0) return
  write(6, '(A/(1X,131A1))') 'NEXTWORD HAS READ OFF THE END', &
    (CARD(I), I = 1, LAST)
  stop 1

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
!   Each line: F18.2, F7.3, 5E11.4
!     Column 1-18:  species code (e.g. 60808.00)
!     Column 19-25: first equilibrium constant
!     Column 26-80: five more equilibrium constants
!   Terminated by a line with code = 0.0
!=========================================================================

SUBROUTINE READMOL

  implicit none

  ! --- Code extraction table ---
  ! XCODE(I) is the place value for the I-th pair of digits in the
  ! species code. Dividing the code by XCODE and truncating extracts
  ! the atomic number encoded at that position.
  real*8, parameter :: XCODE(8) = &
    (/ 1.0D14, 1.0D12, 1.0D10, 1.0D8, 1.0D6, 1.0D4, 1.0D2, 1.0D0 /)

  ! --- Local variables ---
  real*8  :: C                   ! species code read from file
  real*8  :: E1, E2, E3, E4, E5, E6  ! equilibrium constants
  real*8  :: X                   ! remaining code after extracting components
  integer :: JMOL               ! molecule counter
  integer :: KLOC               ! position in KCOMPS array
  integer :: ID                  ! atomic number extracted from code
  integer :: ION                 ! ionization state from decimal part
  integer :: II                  ! starting position in XCODE for decoding
  integer :: I, IEQUA, IOS

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING READMOL'

  open(UNIT=2, FILE=trim(DATADIR)//'molecules.dat', &
       STATUS='OLD', ACTION='READ', iostat=IOS)
    if (IOS /= 0) then
      stop 'molecules.dat file not found'
    end if

  ! --- Initialize: no elements require equations yet ---
  IFEQUA(:) = 0

  ! --- Read molecules and build component table ---
  KLOC    = 1
  LOCJ(1) = 1

  do JMOL = 1, maxmol

    if (KLOC > MAXLOC) then
      write(6, '(A)') ' READMOL ERROR: TOO MANY COMPONENTS'
      close(UNIT=2)
      stop 1
    end if

    read(2, '(F18.2, F7.3, 5E11.4)') C, E1, E2, E3, E4, E5, E6

    ! Code = 0 terminates the list
    if (C == 0.0D0) exit

    if (IDEBUG == 1) write(6, '(I5, F18.2, F7.3, 1P5E11.4)') &
         JMOL, C, E1, E2, E3, E4, E5, E6

    ! --- Decode species code into atomic components ---
    ! Find the leading non-zero pair of digits
    II = 0
    do I = 1, 8
      if (C >= XCODE(I)) then
        II = I
        exit
      end if
    end do
    if (II == 0) then
      write(6, '(A,F18.2)') ' READMOL ERROR: invalid species code ', C
      close(UNIT=2)
      stop 1
    end if

    ! Extract atomic numbers from successive digit pairs
    X = C
    do I = II, 8
      ID = int(X / XCODE(I) + 0.5D0)
      X  = X - dble(ID) * XCODE(I)
      if (ID == 0) ID = 100       ! code 00 → neutral atom (element 100)
      IFEQUA(ID) = 1              ! flag: this element needs an equation
      KCOMPS(KLOC) = ID
      KLOC = KLOC + 1
    end do

    ! Extract ionization state from decimal part: .NN → NN electrons removed
    ! Each removed electron adds component 101 (free electron)
    ION = int(X * 100.0D0 + 0.5D0)
    if (ION >= 1) then
      IFEQUA(100) = 1    ! neutral atom appears in ion balance
      IFEQUA(101) = 1    ! electron appears in ion balance
      do I = 1, ION
        KCOMPS(KLOC) = 101
        KLOC = KLOC + 1
      end do
    end if

    ! Store molecule data
    LOCJ(JMOL + 1) = KLOC
    XNMOLCODE(JMOL) = C
    EQUIL(1, JMOL) = E1
    EQUIL(2, JMOL) = E2
    EQUIL(3, JMOL) = E3
    EQUIL(4, JMOL) = E4
    EQUIL(5, JMOL) = E5
    EQUIL(6, JMOL) = E6

  end do

  NUMMOL = JMOL - 1
  if (NUMMOL >= maxmol .and. C /= 0.0D0) then
    write(6, '(A,I5)') ' READMOL ERROR: molecule count exceeds maxmol =', maxmol
    close(UNIT=2)
    stop 1
  end if
  NLOC   = KLOC - 1

  ! --- Assign equation numbers to each element ---
  ! Equation 1 is for total particle number (XNATOM).
  ! Subsequent equations are assigned in order of element number.
  ! If any ions are present, the last equation is charge conservation
  ! (variable NEQUA = XNE, variable NEQUA1 = 1/XNE).
  IEQUA = 1
  do I = 1, 100
    if (IFEQUA(I) == 0) cycle
    IEQUA = IEQUA + 1
    IFEQUA(I) = IEQUA        ! element I → equation number IEQUA
    IDEQUA(IEQUA) = I        ! equation IEQUA → element I
  end do

  NEQUA  = IEQUA
  NEQUA1 = NEQUA + 1
  IFEQUA(101) = NEQUA1       ! electron → variable index NEQUA+1
  NEQNEQ = NEQUA**2

  ! --- Remap KCOMPS from element IDs to equation indices ---
  ! After this, KCOMPS(K) gives the equation number for the K-th
  ! component, ready for direct use in building the rate matrix.
  do I = 1, NLOC
    ID = KCOMPS(I)
    KCOMPS(I) = IFEQUA(ID)
  end do

  write(6,*)
  write(6, '(A,I4,A,I4)') ' MOLECULES  USED', NUMMOL, '  MAX', MAXMOL
  write(6, '(A,I4,A,I4)') ' COMPONENTS USED', NLOC,   '  MAX', MAXLOC
  write(6, '(A,I4,A,I4)') ' EQUATIONS  USED', NEQUA,  '  MAX', MAXEQ

  close(UNIT=2)
  return

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

  implicit none

  ! --- Arguments ---
  real*8,  intent(in)    :: CODE
  integer, intent(in)    :: MODE
  real*8,  intent(inout) :: NUMBER(kw, *)

  ! --- Persistent: tracks whether electron density is current ---
  integer, save :: ITEMP_PREV = 0

  ! --- Local variables ---
  integer :: IZ       ! atomic number
  integer :: NION     ! number of ionization stages requested
  integer :: NNNN     ! number of stages to convert to number density
  integer :: J, ION

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COMPUTE_ONE_POP'

  !=====================================================================
  ! Branch: molecular equilibrium ON or OFF
  !=====================================================================
  if (IFMOL == 1) then

    ! --- Molecular equilibrium path ---
    ! NMOLEC solves the full molecular + ionization equilibrium system,
    ! including electron density, in one shot.
    if (IFPRES == 1 .and. ITEMP /= ITEMP_PREV) call NMOLEC(MODE)
    ITEMP_PREV = ITEMP

    if (CODE == 0.0D0) return

    ! Look up requested species from molecular equilibrium tables
    call MOLEC(CODE, MODE, NUMBER)

  else

    ! --- Atomic-only path (no molecules) ---
    ! NELECT solves for electron density from Saha ionization alone.
    if (IFPRES == 1 .and. ITEMP /= ITEMP_PREV) call NELECT
    ITEMP_PREV = ITEMP

    if (CODE == 0.0D0) return

    if (CODE >= 100.0D0) then
      ! Molecular species requested but molecules are off
      write(6, '(A)') ' COMPUTE_ONE_POP ERROR: molecular species requested but IFMOL=0'
      stop 1
    end if

    ! --- Decode species code ---
    IZ   = int(CODE)                              ! atomic number
    NION = int((CODE - dble(IZ)) * 100.0D0 + 1.5D0)  ! number of ion stages

    ! --- Compute ionization fractions at each depth via Saha equation ---
    do J = 1, NRHOX
      call PFSAHA(J, IZ, NION, MODE, NUMBER)

      ! Convert fractions to absolute number densities:
      !   n_ion = fraction * n_total_atoms * abundance_of_element
      ! For MODE < 10: only the first stage (ground ionization)
      ! For MODE >= 10: all NION stages
      NNNN = 1
      if (MODE >= 10) NNNN = NION

      do ION = 1, NNNN
        NUMBER(J, ION) = NUMBER(J, ION) * XNATOM(J) * XABUND(J, IZ)
      end do
    end do

  end if

  return

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

  implicit none

  ! --- Temperature perturbation for finite differences ---
  real*8, parameter :: DELTA_P     = 1.001D0     ! T+ multiplier
  real*8, parameter :: DELTA_M     = 0.999D0     ! T- multiplier
  real*8, parameter :: RATIO_PM    = 0.999D0 / 1.001D0  ! T- / T+ (to go from T+ to T-)
  real*8, parameter :: DLNUDT_SCALE = 1000.0D0   ! 2 / (ln(1.001) - ln(0.999))

  ! --- Number of ionization stages per element (for PFSAHA calls) ---
  integer, parameter :: NION_Z(30) = (/ &
    2, 3, 4, 5, 5,                        &  ! Z = 1-5
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,     &  ! Z = 6-16
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  &  ! Z = 17-28
    3, 3 /)                                   ! Z = 29-30
  ! Z = 31-99: all use NION = 3

  ! --- Local arrays ---
  real*8  :: PFPLUS(mion)    ! partition functions at T+
  real*8  :: PFMINUS(mion)   ! partition functions at T-

  ! --- Local variables ---
  real*8  :: XNTOT           ! total particle density (atoms + electrons)
  real*8  :: UPF, DLNU       ! partition function sum and log-derivative
  integer :: J, IION, NOFF

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING ENERGY_DENSITY'
  if (IFEDNS == 0) return

  do J = 1, NRHOX

    ! --- Compute partition functions at T+ = T * 1.001 ---
    call perturb_T(J, DELTA_P)
    call compute_partfcns(J, PFPLUS)

    ! --- Compute partition functions at T- = T * 0.999 ---
    ! T is currently at T*1.001, so multiply by 0.999/1.001
    call perturb_T(J, RATIO_PM)
    call compute_partfcns(J, PFMINUS)

    ! --- Restore original temperature ---
    ! T is currently at T*0.999, so divide by 0.999
    call perturb_T(J, 1.0D0 / DELTA_M)

    ! --- Assemble energy density ---
    ! Kinetic: (3/2) n_total kT
    XNTOT = XNE(J) + XNATOM(J)
    EDENS(J) = 1.5D0 * XNTOT * TK(J)

    ! Ionization + excitation for each ion stage (IION = 1..840)
    do IION = 1, 840
      UPF  = PFPLUS(IION) + PFMINUS(IION) + 1.0D-30
      DLNU = (PFPLUS(IION) - PFMINUS(IION)) / UPF * DLNUDT_SCALE

      EDENS(J) = EDENS(J) + XNF(J, IION) * TK(J) &
               * (POTIONSUM(IION) * HCKT(J) + DLNU)
    end do

    ! Convert to energy per unit mass
    EDENS(J) = EDENS(J) / RHO(J)

  end do

  return

contains

  !---------------------------------------------------------------------
  ! Multiply temperature and related arrays at depth J by factor FAC.
  !---------------------------------------------------------------------
  subroutine perturb_T(J_in, FAC)
    integer, intent(in) :: J_in
    real*8,  intent(in) :: FAC
    real*8 :: FAC_INV
    FAC_INV = 1.0D0 / FAC
    T(J_in)    = T(J_in) * FAC
    TK(J_in)   = TK(J_in) * FAC
    TKEV(J_in) = TKEV(J_in) * FAC
    HCKT(J_in) = HCKT(J_in) * FAC_INV
    HKT(J_in)  = HKT(J_in) * FAC_INV
    TLOG(J_in) = log(T(J_in))
  end subroutine perturb_T

  !---------------------------------------------------------------------
  ! Compute partition functions for all elements (Z=1-99) at depth J
  ! via PFSAHA MODE 5. Results stored in PF(:) at the standard NELION
  ! offsets: IZ*(IZ+1)/2 for Z=1-30, 496+(IZ-31)*5 for Z=31-99.
  !---------------------------------------------------------------------
  subroutine compute_partfcns(J_in, PF)
    integer, intent(in)  :: J_in
    real*8,  intent(out) :: PF(mion)
    integer :: IZ_loc

    do IZ_loc = 1, 30
      NOFF = IZ_loc * (IZ_loc + 1) / 2
      call PFSAHA(J_in, IZ_loc, NION_Z(IZ_loc), 5, PF(NOFF))
    end do
    do IZ_loc = 31, 99
      NOFF = 496 + (IZ_loc - 31) * 5
      call PFSAHA(J_in, IZ_loc, 3, 5, PF(NOFF))
    end do
  end subroutine compute_partfcns

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

  implicit none

  ! --- Convergence parameters ---
  integer, parameter :: MAX_ITER_XNE = 200
  real*8,  parameter :: XNE_TOL      = 1.0D-4   ! relative convergence

  ! --- Number of ionization stages per element for Saha computation ---
  ! These differ slightly from COMPUTE_ALL_POPS because NELECT needs
  ! enough stages to capture charge balance accurately.
  integer, parameter :: NION_NELECT(30) = (/ &
    2, 3, 4, 4, 4,                        &  ! Z = 1-5
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,     &  ! Z = 6-16
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  &  ! Z = 17-28
    3, 3 /)                                   ! Z = 29-30
  ! Z = 31-99: all use NION = 3

  ! --- Local variables ---
  real*8  :: XNTOT           ! total particle density = P / kT
  real*8  :: XNENEW          ! new electron density estimate
  real*8  :: CHARGESQUARE    ! sum of n_ion * q^2
  real*8  :: ERROR           ! relative change in XNE
  real*8  :: DEBYE           ! Debye shielding length
  integer :: J, IZ, ION, IION, NOFF, ITER_XNE

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING NELECT'

  !---------------------------------------------------------------------
  ! Loop over all depth points
  !---------------------------------------------------------------------
  do J = 1, NRHOX

    ! Total particle density from ideal gas law
    XNTOT = P(J) / TK(J)
    XNATOM(J) = XNTOT - XNE(J)

    !-------------------------------------------------------------------
    ! Iterate to converge electron density
    !-------------------------------------------------------------------
    do ITER_XNE = 1, MAX_ITER_XNE

      ! Zero all ion populations
      XNF(J, 1:MION) = 0.0D0

      XNENEW = 0.0D0
      CHARGESQUARE = 0.0D0

      ! --- Compute Saha ionization for elements Z=1-30 ---
      do IZ = 1, 30
        NOFF = IZ * (IZ + 1) / 2
        call PFSAHA(J, IZ, NION_NELECT(IZ), 12, XNF(1, NOFF))
      end do

      ! --- Accumulate electron density and charge sum for Z=1-30 ---
      ! Each element IZ has IZ+1 ion stages in the NELION layout
      IION = 0
      do IZ = 1, 30
        do ION = 1, IZ + 1
          IION = IION + 1
          ! Convert fraction to number density
          XNF(J, IION) = XNF(J, IION) * XNATOM(J) * XABUND(J, IZ)
          ! Charge = ION - 1 (neutral = 0, singly ionized = 1, etc.)
          CHARGESQUARE = CHARGESQUARE + XNF(J, IION) * dble((ION - 1)**2)
          XNENEW = XNENEW + XNF(J, IION) * dble(ION - 1)
        end do
      end do

      ! --- Elements Z=31-99: 3 stages computed, 5 slots allocated ---
      do IZ = 31, 99
        NOFF = 496 + (IZ - 31) * 5
        call PFSAHA(J, IZ, 3, 12, XNF(1, NOFF))
        do ION = 1, 5
          IION = IION + 1
          XNF(J, IION) = XNF(J, IION) * XNATOM(J) * XABUND(J, IZ)
          CHARGESQUARE = CHARGESQUARE + XNF(J, IION) * dble((ION - 1)**2)
          XNENEW = XNENEW + XNF(J, IION) * dble(ION - 1)
        end do
      end do

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

      if (ERROR < XNE_TOL) exit

    end do  ! ITER_XNE

    if (ITER_XNE > MAX_ITER_XNE) then
      write(6, '(A,I4,A,ES10.3)') &
        ' NELECT WARNING: XNE did not converge at J=', J, ', error=', ERROR
    end if

    ! Mass density from atom count and mean molecular weight
    RHO(J) = XNATOM(J) * WTMOLE(J) * AMU

  end do  ! depth loop

  ! --- Debug output ---
  if (IDEBUG == 1) then
    write(6, '(3X,4X,A,A,A,A,A,A,A,A)') &
      'RHOX', '     T    ', '     P     ', '     XNE    ', &
      '    XNATOM  ', '    WTMOLE  ', '     RHO    ', '    CHARGESQ'
    do J = 1, NRHOX
      write(6, '(I3,1PE15.7,0PF10.1,1P6E12.3)') &
        J, RHOX(J), T(J), P(J), XNE(J), XNATOM(J), WTMOLE(J), &
        RHO(J), CHARGESQ(J)
    end do
  end if

  return

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

  implicit none

  ! Arguments
  integer, intent(in)    :: J       ! Depth point index
  integer, intent(in)    :: IZ      ! Atomic number
  integer, intent(in)    :: NION    ! Number of ionization stages requested
  integer, intent(in)    :: MODE    ! Output mode selector
  real*8,  intent(inout) :: ANSWER(kw, *)

  ! External function

  ! NNN interpolation table (loaded from file on first call)
  integer, save :: NNN(6, 365)
  logical, save :: INITIALIZED = .false.

  ! Local arrays
  real*8  :: F(31), IP(31), PART(31), POTLO(31)
  integer :: LOCZ(29)

  ! Scale factors for NNN unpacking
  real*8 :: SCALE(4)

  ! Detailed energy level tables for special elements (cm^-1)
  real*8 :: EHYD(6), GHYD(6)
  real*8 :: EHE1(29), GHE1(29), EHE2(6), GHE2(6)
  real*8 :: EB1(7), GB1(7)
  real*8 :: EC1(14), GC1(14), EC2(6), GC2(6)
  real*8 :: EO1(13), GO1(13)
  real*8 :: ENA1(8), GNA1(8)
  real*8 :: EMG1(11), GMG1(11), EMG2(6), GMG2(6)
  real*8 :: EAL1(9), GAL1(9)
  real*8 :: ESI1(11), GSI1(11), ESI2(6), GSI2(6)
  real*8 :: EK1(8), GK1(8)
  real*8 :: ECA1(8), GCA1(8), ECA2(5), GCA2(5)

  ! Local scalars
  real*8  :: Z_ion, TV, DEBYE, POTLOW, T2000, DT, P1, P2, PMIN
  real*8  :: G, D1, D2, CF, B_dep
  integer :: MODE1, N, NIONS, NION2, ION, I, L, IT, NNN100, INDEX
  integer :: K1, K2, K3, KSCALE, KP1
  integer :: jj  ! loop variable for data read

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
  if (IDEBUG == 1) write(6,'(A)') ' RUNNING PFSAHA'
  if (.not. INITIALIZED) then
    open(unit=89, file=TRIM(DATADIR)//'pfsaha.dat', status='OLD', action='READ')
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')  ! skip 3 header lines
    read(89, *) ((NNN(I,jj), I=1,6), jj=1,365)
    close(89)
    ! Ensure ionization potentials are loaded (PFSAHA needs POTION for IP)
    call IONPOTS
    INITIALIZED = .true.
  endif

  !=====================================================================
  ! Setup: determine NNN row range and number of ions for this element
  !=====================================================================
  MODE1 = MODE
  if (MODE1 > 10) MODE1 = MODE1 - 10

  ! Debye-Hückel lowering of ionization potential (volts per unit Zeff)
  DEBYE  = SQRT(TK(J) / FOURPI / ECHARGE**2 / CHARGESQ(J))
  POTLOW = MIN(1.D0, 1.44D-7 / DEBYE)
  TV     = TKEV(J)

  ! Locate element in NNN table
  if (IZ <= 28) then
    N     = LOCZ(IZ)
    NIONS = LOCZ(IZ+1) - N
  else
    N     = 3*IZ + 54
    NIONS = 3
  endif

  ! Override for C and N (extended tables in NNN67 block)
  if (IZ == 6) then; N = 354; NIONS = 6; endif
  if (IZ == 7) then; N = 360; NIONS = 6; endif

  ! Iron group elements have 10 ionization stages in PFIRON
  if (IZ >= 20 .and. IZ < 29) NIONS = 10

  NION2 = MIN(NION + 2, NIONS)
  N = N - 1   ! will be incremented at start of loop

  !=====================================================================
  ! Phase 1: Compute partition functions PART(ION) for each ion stage
  !=====================================================================
  do ION = 1, NION2
    Z_ion = dble(ION)
    POTLO(ION) = POTLOW * Z_ion
    N = N + 1

    ! Ionization potential from POTION table
    NNN100 = NNN(6,N) / 100
    if (IZ <= 30) then
      INDEX = (IZ*(IZ+1))/2 + ION - 1
    else
      INDEX = IZ*5 + 341 + ION - 1
    endif
    IP(ION) = POTION(INDEX) / 8065.479d0
    if (IP(ION) == 0.d0) IP(ION) = POTION(INDEX-1) / 8065.479d0

    ! Iron group: use PFIRON for partition functions
    if (IZ >= 20 .and. IZ < 29) then
      call PFIRON(IZ, ION, TLOG(J)/2.30258509299405d0, &
                  POTLO(ION)*8065.479d0, PART(ION))
      cycle  ! next ION
    endif

    ! Statistical weight of highest included level
    G = NNN(6,N) - NNN100*100

    !-----------------------------------------------------------------
    ! Detailed level summation for special elements
    ! B_dep = NLTE departure coefficient (always 1.0 in LTE mode)
    !-----------------------------------------------------------------
    B_dep = 1.0d0

    select case (N)

    case (1)  ! ---- Hydrogen I ----
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
      do I = 2, 6
        PART(1) = PART(1) + occupation_prob(I, XNE(J)) &
                * GHYD(I) * B_dep * EXP(-EHYD(I)*HCKT(J))
      end do
      ! Levels n=7 and above: hydrogenic energies, weighted by w_n.
      ! Sum truncates naturally where w_n -> 0 or Boltzmann factor -> 0.
      do I = 7, 80
        D1 = occupation_prob(I, XNE(J))
        if (D1 < 1.0d-6) exit   ! remaining levels fully dissolved
        PART(1) = PART(1) + D1 &
                * 2.0d0 * dble(I)**2 * EXP(-(ELIM_HI - RYDBERG_H/dble(I)**2)*HCKT(J))
      end do
      cycle

    case (3)  ! ---- Helium I ----
      PART(1) = B_dep
      do I = 2, 29
        PART(1) = PART(1) + GHE1(I) * B_dep * EXP(-EHE1(I)*HCKT(J))
      end do
      D1 = RYDBERG_H / 5.5d0 / 5.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (4)  ! ---- Helium II ----
      PART(2) = 2.d0 * B_dep
      do I = 2, 6
        PART(2) = PART(2) + GHE2(I) * B_dep * EXP(-EHE2(I)*HCKT(J))
      end do
      D1 = 4.d0 * RYDBERG_HE / 6.5d0 / 6.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (14)  ! ---- Boron I ----
      PART(1) = B_dep * (2.d0 + 4.d0*EXP(-15.25d0*HCKT(J)))
      do I = 2, 7
        PART(1) = PART(1) + GB1(I) * B_dep * EXP(-EB1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 6.d0*EXP(-57786.80d0*HCKT(J)) &
        + 10.d0*EXP(-59989.d0*HCKT(J)) + 14.d0*EXP(-60031.03d0*HCKT(J)) &
        + 2.d0*EXP(-63561.d0*HCKT(J))
      G = 2.d0
      D1 = 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (45)  ! ---- Sodium I ----
      PART(1) = B_dep * 2.d0
      do I = 2, 8
        PART(1) = PART(1) + GNA1(I) * B_dep * EXP(-ENA1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 10.d0*EXP(-34548.745d0*HCKT(J)) &
        + 14.d0*EXP(-34586.96d0*HCKT(J))
      G = 2.d0
      D1 = 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (51)  ! ---- Magnesium I ----
      PART(1) = B_dep
      do I = 2, 11
        PART(1) = PART(1) + GMG1(I) * B_dep * EXP(-EMG1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 5.d0*EXP(-53134.d0*HCKT(J)) &
        + 15.d0*EXP(-54192.d0*HCKT(J)) + 28.d0*EXP(-54676.d0*HCKT(J)) &
        + 9.d0*EXP(-57853.d0*HCKT(J))
      G = 4.d0
      D1 = 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (52)  ! ---- Magnesium II ----
      PART(2) = B_dep * 2.d0
      do I = 2, 6
        PART(2) = PART(2) + GMG2(I) * B_dep * EXP(-EMG2(I)*HCKT(J))
      end do
      PART(2) = PART(2) + 10.d0*EXP(-93310.80d0*HCKT(J)) &
        + 14.d0*EXP(-93799.70d0*HCKT(J)) + 6.d0*EXP(-97464.32d0*HCKT(J)) &
        + 10.d0*EXP(-103419.82d0*HCKT(J)) + 14.d0*EXP(-103689.89d0*HCKT(J)) &
        + 18.d0*EXP(-103705.66d0*HCKT(J))
      G = 2.d0
      D1 = 4.d0 * 109734.83d0 / 5.5d0 / 5.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (57)  ! ---- Aluminum I ----
      PART(1) = B_dep * (2.d0 + 4.d0*EXP(-112.061d0*HCKT(J)))
      do I = 2, 9
        PART(1) = PART(1) + GAL1(I) * B_dep * EXP(-EAL1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 10.d0*EXP(-42235.d0*HCKT(J)) &
        + 14.d0*EXP(-43831.d0*HCKT(J))
      G = 2.d0
      D1 = 109735.08d0 / 5.5d0 / 5.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (63)  ! ---- Silicon I ----
      PART(1) = B_dep * (1.d0 + 3.d0*EXP(-77.115d0*HCKT(J)) &
        + 5.d0*EXP(-223.157d0*HCKT(J)))
      do I = 2, 11
        PART(1) = PART(1) + GSI1(I) * B_dep * EXP(-ESI1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 76.d0*EXP(-53000.d0*HCKT(J)) &
        + 71.d0*EXP(-57000.d0*HCKT(J)) + 191.d0*EXP(-60000.d0*HCKT(J)) &
        + 240.d0*EXP(-62000.d0*HCKT(J)) + 251.d0*EXP(-63000.d0*HCKT(J)) &
        + 300.d0*EXP(-65000.d0*HCKT(J))
      cycle  ! Si I has no high-level correction (goes direct to 18)

    case (64)  ! ---- Silicon II ----
      PART(2) = B_dep * (2.d0 + 4.d0*EXP(-287.32d0*HCKT(J)))
      do I = 2, 6
        PART(2) = PART(2) + GSI2(I) * B_dep * EXP(-ESI2(I)*HCKT(J))
      end do
      PART(2) = PART(2) + 6.d0*EXP(-81231.59d0*HCKT(J)) &
        + 6.d0*EXP(-83937.08d0*HCKT(J)) + 10.d0*EXP(-101024.09d0*HCKT(J)) &
        + 14.d0*EXP(-103556.35d0*HCKT(J)) + 10.d0*EXP(-108800.d0*HCKT(J)) &
        + 42.d0*EXP(-115000.d0*HCKT(J)) + 6.d0*EXP(-121000.d0*HCKT(J)) &
        + 38.d0*EXP(-125000.d0*HCKT(J)) + 34.d0*EXP(-132000.d0*HCKT(J))
      G = 2.d0
      D1 = 4.d0 * 109734.83d0 / 4.5d0 / 4.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (91)  ! ---- Potassium I ----
      PART(1) = B_dep * 2.d0
      do I = 2, 8
        PART(1) = PART(1) + GK1(I) * B_dep * EXP(-EK1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 10.d0*EXP(-27397.077d0*HCKT(J)) &
        + 14.d0*EXP(-28127.85d0*HCKT(J))
      G = 2.d0
      D1 = 109734.83d0 / 5.5d0 / 5.5d0 * HCKT(J)
      call pfsaha_highlevels(PART(ION), G, IP(ION), Z_ion, TV, &
                             POTLO(ION), D1)
      cycle

    case (354)  ! ---- Carbon I (extended, NNN67) ----
      PART(1) = B_dep * (1.d0 + 3.d0*EXP(-16.42d0*HCKT(J)) &
        + 5.d0*EXP(-43.42d0*HCKT(J)))
      do I = 2, 14
        PART(1) = PART(1) + GC1(I) * B_dep * EXP(-EC1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 108.d0*EXP(-80000.d0*HCKT(J)) &
        + 189.d0*EXP(-84000.d0*HCKT(J)) + 247.d0*EXP(-87000.d0*HCKT(J)) &
        + 231.d0*EXP(-88000.d0*HCKT(J)) + 190.d0*EXP(-89000.d0*HCKT(J)) &
        + 300.d0*EXP(-90000.d0*HCKT(J))
      cycle  ! C I has no high-level correction

    case (355)  ! ---- Carbon II (extended) ----
      PART(2) = B_dep * (2.d0 + 4.d0*EXP(-63.42d0*HCKT(J)))
      do I = 2, 6
        PART(2) = PART(2) + GC2(I) * B_dep * EXP(-EC2(I)*HCKT(J))
      end do
      PART(2) = PART(2) + 6.d0*EXP(-131731.80d0*HCKT(J)) &
        + 4.d0*EXP(-142027.1d0*HCKT(J)) + 10.d0*EXP(-145550.13d0*HCKT(J)) &
        + 10.d0*EXP(-150463.62d0*HCKT(J)) + 2.d0*EXP(-157234.07d0*HCKT(J)) &
        + 6.d0*EXP(-162500.d0*HCKT(J)) + 42.d0*EXP(-168000.d0*HCKT(J)) &
        + 56.d0*EXP(-178000.d0*HCKT(J)) + 102.d0*EXP(-183000.d0*HCKT(J)) &
        + 400.d0*EXP(-188000.d0*HCKT(J))
      cycle  ! C II has no high-level correction

    case (367)  ! ---- Oxygen I (extended) ----
      PART(1) = B_dep * (5.d0 + 3.d0*EXP(-158.265d0*HCKT(J)) &
        + EXP(-226.977d0*HCKT(J)))
      do I = 2, 13
        PART(1) = PART(1) + GO1(I) * B_dep * EXP(-EO1(I)*HCKT(J))
      end do
      PART(1) = PART(1) + 15.d0*EXP(-101140.d0*HCKT(J)) &
        + 131.d0*EXP(-103000.d0*HCKT(J)) + 128.d0*EXP(-105000.d0*HCKT(J)) &
        + 600.d0*EXP(-107000.d0*HCKT(J))
      cycle  ! O I has no high-level correction

    case default
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

      if (MOD(IT,2) == 0) then
        ! Even IT: interpolate between K3 and next column K1
        P1 = K3 * SCALE(KSCALE)
        K1 = NNN(I+1,N) / 100000
        KSCALE = MOD(NNN(I+1,N), 10)
        P2 = K1 * SCALE(KSCALE)
      else
        ! Odd IT: interpolate between K1 and K3
        P1 = K1 * SCALE(KSCALE)
        P2 = K3 * SCALE(KSCALE)
        if (DT < 0.d0 .and. KSCALE <= 1) then
          KP1 = int(P1)
          if (KP1 == INT(P2 + 0.5d0)) PMIN = dble(KP1)
        endif
      endif

      PART(ION) = MAX(PMIN, P1 + (P2 - P1)*DT)

      ! Patch (14 Jun 2004, J. Laird): ensure PFGROUND >= PFSAHA at low T
      if (T(J) < T2000*2.0d0) then
        PART(ION) = MAX(PFGROUND((IZ-1)*6 + ION, T(J)), PART(ION))
        cycle  ! skip high-level correction at low T
      endif

      ! High-level (Rydberg series) correction
      if (G == 0.d0 .or. POTLO(ION) < 0.1d0 .or. T(J) < T2000*4.d0) cycle
      if (T(J) > T2000*11.d0) TV = (T2000*11.d0) * 8.6171D-5
      D1 = 0.1d0 / TV
      D2 = POTLO(ION) / TV
      PART(ION) = PART(ION) + G*EXP(-IP(ION)/TV) * &
        ( SQRT(13.595d0*Z_ion*Z_ion/TV/D2)**3 * &
          (1.d0/3.d0 + (1.d0 - (0.5d0 + (1.d0/18.d0 + D2/120.d0)*D2)*D2)*D2) &
        - SQRT(13.595d0*Z_ion*Z_ion/TV/D1)**3 * &
          (1.d0/3.d0 + (1.d0 - (0.5d0 + (1.d0/18.d0 + D1/120.d0)*D1)*D1)*D1) )
      TV = TKEV(J)
      cycle

    end select
  end do  ! ION

  !=====================================================================
  ! Phase 2: Saha equation and output
  !=====================================================================
  if (MODE1 == 5) then
    ! MODE 5: partition functions + cumulative ionization potentials
    ANSWER(32,1) = 0.d0
    do ION = 1, NION
      ANSWER(ION,1)    = PART(ION)
      ANSWER(ION+32,1) = IP(ION) + ANSWER(ION+31,1)
    end do
    return
  end if

  if (MODE1 /= 3) then
    ! Compute Saha ionization fractions
    N = N - NION2
    CF = SAHA_PREFAC * T(J) * SQRT(T(J)) / XNE(J)

    do ION = 2, NION2
      N = N + 1
      F(ION) = CF * PART(ION) / PART(ION-1) * EXP(-(IP(ION-1) - POTLO(ION-1))/TV)
    end do

    ! Normalize fractions
    F(1) = 1.d0
    L = NION2 + 1
    do ION = 2, NION2
      L = L - 1
      F(1) = 1.d0 + F(L)*F(1)
    end do
    F(1) = 1.d0 / F(1)

    do ION = 2, NION2
      F(ION) = F(ION-1) * F(ION)
    end do
  end if

  !---------------------------------------------------------------------
  ! Fill ANSWER based on MODE
  !---------------------------------------------------------------------
  if (MODE >= 10) then
    ! Multi-ion output modes (MODE = 11, 12, 13, 14)
    select case (MODE1)
    case (1)
      do ION = 1, NION
        ANSWER(J,ION) = F(ION) / PART(ION)
      end do
    case (2)
      do ION = 1, NION
        ANSWER(J,ION) = F(ION)
      end do
    case (3)
      do ION = 1, NION
        ANSWER(J,ION) = PART(ION)
      end do
    case (4)
      ANSWER(J,1) = 0.d0
      do ION = 2, NION2
        ANSWER(J,1) = ANSWER(J,1) + F(ION) * dble(ION-1)
      end do
    end select
  else
    ! Single-ion output modes (MODE = 1, 2, 3, 4)
    select case (MODE1)
    case (1);  ANSWER(J,1) = F(NION) / PART(NION)
    case (2);  ANSWER(J,1) = F(NION)
    case (3);  ANSWER(J,1) = PART(NION)
    case (4)
      ANSWER(J,1) = 0.d0
      do ION = 2, NION2
        ANSWER(J,1) = ANSWER(J,1) + F(ION) * dble(ION-1)
      end do
    end select
  endif
  return

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
  implicit none
  real*8, intent(inout) :: PART_ION
  real*8, intent(in)    :: G, IP_ION, Z_ion, TV_in, POTLO_ION, D1
  real*8 :: D2

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

  implicit none

  ! Arguments
  real*8,  intent(in)    :: CODOUT
  integer, intent(in)    :: MODE
  real*8,  intent(inout) :: NUMBER(kw, 1)

  ! Local variables
  integer :: JMOL, J, NN, ION, ID, I, II
  real*8  :: C
  logical :: found

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING MOLEC'

  !=====================================================================
  ! Read molecular input data (first call only)
  !=====================================================================
  if (IFPOP /= 2 .and. MOLEC_IREAD /= 1 .and. IFPRES /= 1) then
    read(INPUTDATA, '(I5)') NUMMOL

    do JMOL = 1, NUMMOL
      read(INPUTDATA, '(F20.2)') XNMOLCODE(JMOL)
      read(INPUTDATA, '(1P8E10.3)') (XNMOL(J,JMOL), J=1,NRHOX)
      if (IDEBUG == 1) write(6, '(F20.2/(1P8E10.3))') XNMOLCODE(JMOL), (XNMOL(J,JMOL), J=1,NRHOX)
    end do

    read(INPUTDATA, '(1P8E10.3)')  ! skip header line
    read(INPUTDATA, '(1P8E10.3)') (XNATOM(J), RHO(J), J=1,NRHOX)
    if (IDEBUG == 1) write(6, '(1P8E10.3)') (XNATOM(J), RHO(J), J=1,NRHOX)

    read(INPUTDATA, '(1P8E10.3)')  ! skip header line
    read(INPUTDATA, '(1P8E10.3)') (XNE(J), J=1,NRHOX)
    if (IDEBUG == 1) write(6, '(1P8E10.3)') (XNE(J), J=1,NRHOX)

    MOLEC_IREAD = 1
  endif

  !=====================================================================
  ! Path 1: Molecular species lookup (CODOUT >= 100)
  !=====================================================================
  if (CODOUT >= 100.d0) then
    do JMOL = 1, NUMMOL
      if (XNMOLCODE(JMOL) == CODOUT) then
        do J = 1, NRHOX
          if (MODE == 1 .or. MODE == 11) NUMBER(J,1) = XNFPMOL(J,JMOL)
          if (MODE == 2 .or. MODE == 12) NUMBER(J,1) = XNMOL(J,JMOL)
        end do
        return
      endif
    end do
    ! Species not in molecular equilibrium table: return zero populations
    do J = 1, NRHOX
      NUMBER(J,1) = 0.d0
    end do
    return
  endif

  !=====================================================================
  ! Path 2: Atomic element lookup (CODOUT < 100)
  !=====================================================================
  C = CODOUT
  NN = 1
  if (MODE == 11 .or. MODE == 12) NN = INT((C - AINT(C))*100.d0 + 1.5d0)

  do I = 1, NN
    ION = NN - I + 1

    ! Search molecule table for approximate code match
    found = .false.
    do JMOL = 1, NUMMOL
      if (XNMOLCODE(JMOL) + 0.001d0 > C .and. &
          XNMOLCODE(JMOL) - 0.001d0 < C) then
        do J = 1, NRHOX
          if (MODE == 1 .or. MODE == 11) NUMBER(J,ION) = XNFPMOL(J,JMOL)
          if (MODE == 2 .or. MODE == 12) NUMBER(J,ION) = XNMOL(J,JMOL)
        end do
        found = .true.
        exit  ! exit JMOL loop
      endif
    end do

    if (.not. found) then
      ! Try matching by integer element code
      ID = INT(CODOUT)
      found = .false.
      do JMOL = 1, NUMMOL
        if (INT(XNMOLCODE(JMOL)) == ID) then
          ! Element exists but this specific ion not tabulated: zero it
          do J = 1, NRHOX
            NUMBER(J,ION) = 0.d0
          end do
          found = .true.
          exit  ! exit JMOL loop
        endif
      end do

      if (.not. found) then
        ! Element not in molecule table at all: use Saha equation
        ION = INT((CODOUT - dble(ID))*100.d0 + 1.5d0)
        NN = ION
        if (MODE == 1) NN = 1
        do J = 1, NRHOX
          call PFSAHA(J, ID, ION, MODE, NUMBER)
          do II = 1, NN
            NUMBER(J,II) = NUMBER(J,II) * XNATOM(J) * XABUND(J,ID)
          end do
        end do
        return  ! exit entire subroutine (matches original GO TO 400)
      endif
    endif

    C = C - 0.01d0
  end do  ! I

  return

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

  implicit none

  ! --- Arguments ---
  integer, intent(in) :: MODE

  ! --- Local arrays ---
  real*8  :: EQUILJ(maxmol)          ! equilibrium constants at current depth
  real*8  :: XNZ(kw, maxeq)          ! converged number densities per equation
  real*8  :: EQ(maxeq)               ! residual vector
  real*8  :: XN(maxeq)               ! current number density estimates
  real*8  :: XAB(maxeq)              ! element abundances
  real*8  :: DEQ(maxeq * maxeq)      ! Jacobian matrix (flattened)
  real*8  :: EQOLD(maxeq)            ! previous residuals (for oscillation detection)
  integer :: IPIVOT(maxeq)           ! pivot array for SOLVIT

  real*8  :: FRAC(kw, 6)             ! scratch for PFSAHA output
  real*8  :: PFP(61), PFM(61)        ! partition functions at T+/T- (atoms)
  real*8  :: PFPLUS(kw), PFMIN(kw)   ! partition functions at T+/T- (molecules)

  ! For electron contribution diagnostic
  real*8  :: E_CONTRIB(kw, maxmol)
  real*8  :: XE_CONTRIB(kw, maxmol)
  integer :: IDZ(maxmol), NION_MOL(maxmol)

  ! --- Local scalars ---
  real*8  :: XNTOT          ! total particle density = P / kT
  real*8  :: RATIO          ! pressure ratio for depth extrapolation
  real*8  :: TERM, D        ! molecule contribution to Jacobian
  real*8  :: X, XNEQ, XN100, SCALE
  real*8  :: AMASS          ! molecular mass accumulator
  integer :: J, K, M, KK, K1, MK
  integer :: JMOL, NCOMP, ION, ID
  integer :: LOCK, LOCJ1, LOCJ2, LOCM, NEQUAK
  integer :: IFERR
  integer :: JMOL1, JMOL10, NN, NZ = 0, IZ, IZ1, IZ10

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING NMOLEC'
  
  NEQUA1 = NEQUA + 1
  NEQNEQ = NEQUA**2

  !=====================================================================
  ! Initialize abundances (constant with depth in this version)
  !=====================================================================
  XAB(:) = 0.0D0
  J = 1
  do K = 2, NEQUA
    ID = IDEQUA(K)
    if (ID < 100) XAB(K) = max(XABUND(J, ID), 1.0D-20)
  end do
  if (ID == 100) XAB(NEQUA) = 0.0D0

  ! Initial guess for number densities at the surface
  XNTOT = P(1) / TK(1)
  XN(1) = XNTOT / 2.0D0
  if (T(1) < 4000.0D0) XN(1) = XNTOT
  X = XN(1) / 10.0D0
  do K = 2, NEQUA
    XN(K) = X * XAB(K)
  end do
  if (ID == 100) XN(NEQUA) = X
  XNE(1) = X
  
  !=====================================================================
  ! Main depth loop: solve molecular equilibrium at each depth
  !=====================================================================
  do J = 1, NRHOX
     
    XNTOT = P(J) / TK(J)

    ! Extrapolate from previous depth as initial guess
    if (J > 1) then
      RATIO = P(J) / P(J - 1)
      XNE(J) = XNE(J - 1) * RATIO
      do K = 1, NEQUA
        XN(K) = XN(K) * RATIO
      end do
    end if

    ! If computing energy density, restore saved values from previous
    ! normal (non-EDENS) run as the starting point
    if (IFEDNS /= 0) then
      do K = 1, NEQUA
        XN(K) = XNSAVE(J, K)
      end do
    end if

    !-------------------------------------------------------------------
    ! Compute equilibrium constants for all molecules at this depth
    !-------------------------------------------------------------------
    do JMOL = 1, NUMMOL
      NCOMP = LOCJ(JMOL + 1) - LOCJ(JMOL)

      if (EQUIL(1, JMOL) /= 0.0D0) then
        ! True molecule: has equilibrium constant polynomial
        ION = int((XNMOLCODE(JMOL) - aint(XNMOLCODE(JMOL))) * 100.0D0 + 0.5D0)
        EQUILJ(JMOL) = 0.0D0

        if (XNMOLCODE(JMOL) == 101.0D0) then
          ! H2: special equilibrium function (Lester 2005 update)
          if (T(J) <= 20000.0D0) EQUILJ(JMOL) = EQUILH2(T(J))
        else
          ! General molecule: polynomial equilibrium constant
          if (T(J) <= 10000.0D0) then
            EQUILJ(JMOL) = exp(EQUIL(1, JMOL) / TKEV(J) - EQUIL(2, JMOL) &
              + (EQUIL(3, JMOL) + (-EQUIL(4, JMOL) + (EQUIL(5, JMOL) &
              - EQUIL(6, JMOL) * T(J)) * T(J)) * T(J)) * T(J) &
              - 1.5D0 * (NCOMP - ION - ION - 1) * TLOG(J))
          end if
        end if

      else if (NCOMP == 1) then
        ! Single atom: trivial equilibrium
        EQUILJ(JMOL) = 1.0D0

      else
        ! Multi-stage atom: get ionization equilibrium from Saha
        ID = int(XNMOLCODE(JMOL))
        ION = NCOMP - 1
        call PFSAHA(J, ID, NCOMP, 12, FRAC)
        EQUILJ(JMOL) = FRAC(J, NCOMP) / FRAC(J, 1) * XNE(J)**ION
      end if
    end do

    !-------------------------------------------------------------------
    ! Newton-Raphson iteration for chemical equilibrium
    !-------------------------------------------------------------------
    EQOLD(:) = 0.0D0

    newton_loop: do

      ! --- Build Jacobian (DEQ) and residual (EQ) ---
      DEQ(1:NEQNEQ) = 0.0D0

      ! Equation 1: total particle conservation
      EQ(1) = -XNTOT
      K1 = 1
      KK = 1
      do K = 2, NEQUA
        EQ(1) = EQ(1) + XN(K)
        K1 = K1 + NEQUA
        DEQ(K1) = 1.0D0              ! d(eq1)/d(XN_k) = 1
        EQ(K) = XN(K) - XAB(K) * XN(1)  ! element conservation
        KK = KK + NEQUA1
        DEQ(KK) = 1.0D0              ! d(eq_k)/d(XN_k) = 1
        DEQ(K) = -XAB(K)             ! d(eq_k)/d(XN_1) = -A_k
      end do

      ! Charge conservation (if ions present)
      if (IDEQUA(NEQUA) >= 100) then
        EQ(NEQUA) = -XN(NEQUA)
        DEQ(NEQNEQ) = -1.0D0
      end if

      ! --- Add molecular contributions ---
      do JMOL = 1, NUMMOL
        NCOMP = LOCJ(JMOL + 1) - LOCJ(JMOL)
        if (NCOMP == 1) cycle         ! single atoms don't contribute

        ! Compute molecule number density: TERM = K * product(n_k)
        TERM = EQUILJ(JMOL)
        LOCJ1 = LOCJ(JMOL)
        LOCJ2 = LOCJ(JMOL + 1) - 1
        do LOCK = LOCJ1, LOCJ2
          K = KCOMPS(LOCK)
          if (K == NEQUA1) then
            TERM = TERM / XN(NEQUA)   ! electron: divide by n_e
          else
            TERM = TERM * XN(K)        ! atom: multiply by n_k
          end if
        end do

        ! Add to total particle count
        EQ(1) = EQ(1) + TERM

        ! Add to element conservation and Jacobian
        do LOCK = LOCJ1, LOCJ2
          K = KCOMPS(LOCK)
          if (K >= NEQUA1) then
            K = NEQUA
            D = -TERM / XN(K)          ! electron: d(TERM)/d(n_e) = -TERM/n_e
          else
            D = TERM / XN(K)           ! atom: d(TERM)/d(n_k) = TERM/n_k
          end if
          EQ(K) = EQ(K) + TERM
          NEQUAK = NEQUA * K - NEQUA
          K1 = NEQUAK + 1
          DEQ(K1) = DEQ(K1) + D
          do LOCM = LOCJ1, LOCJ2
            M = KCOMPS(LOCM)
            if (M == NEQUA1) M = NEQUA
            MK = M + NEQUAK
            DEQ(MK) = DEQ(MK) + D
          end do
        end do

        ! Correction to charge equation for negative ions
        ! (species whose last component is element 100 = neutral atom)
        K = KCOMPS(LOCJ2)
        if (IDEQUA(K) == 100) then
          do LOCK = LOCJ1, LOCJ2
            K = KCOMPS(LOCK)
            D = TERM / XN(K)
            if (K == NEQUA) EQ(K) = EQ(K) - TERM - TERM
            NEQUAK = NEQUA * K - NEQUA
            do LOCM = LOCJ1, LOCJ2
              M = KCOMPS(LOCM)
              if (M == NEQUA) then
                MK = M + NEQUAK
                DEQ(MK) = DEQ(MK) - D - D
              end if
            end do
          end do
        end if

      end do  ! JMOL

      ! --- Solve linearized system ---
      call SOLVIT(DEQ, NEQUA, EQ, IPIVOT)

      ! --- Apply corrections with damping ---
      IFERR = 0
      SCALE = 100.0D0
      do K = 1, NEQUA
        RATIO = abs(EQ(K) / XN(K))
        if (RATIO > 1.0D-4) IFERR = 1

        ! Reduce correction if oscillating (sign flip)
        if (EQOLD(K) * EQ(K) < 0.0D0) EQ(K) = EQ(K) * 0.69D0

        XNEQ = XN(K) - EQ(K)
        XN100 = XN(K) / 100.0D0

        if (abs(XNEQ) >= XN100) then
          ! Normal correction: accept new value (force positive)
          XN(K) = abs(XNEQ)
        else
          ! Correction too large relative to current value: scale down
          XN(K) = XN(K) / SCALE
          if (EQOLD(K) * EQ(K) < 0.0D0) SCALE = sqrt(SCALE)
        end if

        EQOLD(K) = EQ(K)
      end do

      if (IFERR == 0) exit newton_loop

    end do newton_loop

    !-------------------------------------------------------------------
    ! Store converged solution
    !-------------------------------------------------------------------
    do K = 1, NEQUA
      XNZ(J, K) = XN(K)
    end do
    XNATOM(J) = XN(1)
    RHO(J) = XNATOM(J) * WTMOLE(J) * AMU
    if (IDEQUA(NEQUA) == 100) XNE(J) = XN(NEQUA)

    ! Compute molecular number densities from converged XN
    do JMOL = 1, NUMMOL
      XNMOL(J, JMOL) = EQUILJ(JMOL)
      LOCJ1 = LOCJ(JMOL)
      LOCJ2 = LOCJ(JMOL + 1) - 1
      do LOCK = LOCJ1, LOCJ2
        K = KCOMPS(LOCK)
        if (K == NEQUA1) then
          XNMOL(J, JMOL) = XNMOL(J, JMOL) / XN(NEQUA)
        else
          XNMOL(J, JMOL) = XNMOL(J, JMOL) * XN(K)
        end if
      end do
    end do

  end do  ! depth loop J

  !=====================================================================
  ! If computing energy density (IFEDNS=1), jump to that section
  !=====================================================================
  if (IFEDNS == 1) then
    call nmolec_energy_density()
    return
  end if

  
  !=====================================================================
  ! Save solution for next call / EDENS restart
  !=====================================================================

  do K = 1, NEQUA
    do J = 1, NRHOX
      XNSAVE(J, K) = XNZ(J, K)
    end do
  end do

  ! --- Diagnostic output ---
  if ((ITER == 1 .or. ITER == NUMITS) .and. IDEBUG == 1) then
    write(6, '(10X,"RHOX",9X,"T",11X,"P",10X,"XNE",8X,"XNATOM",8X,"RHO")') 
    do J = 1, NRHOX
      write(6, '(I5,1P6E12.3)') J, RHOX(J), T(J), P(J), XNE(J), XNATOM(J), RHO(J)
    end do
  end if

  if (ITER == NUMITS .and. IFSYNTHE == 1) then
    do JMOL1 = 1, NUMMOL, 10
      JMOL10 = min(JMOL1 + 9, NUMMOL)
      write(35, '(49X,"MOLECULAR NUMBER DENSITIES"/5X,10F12.2)') &
        (XNMOLCODE(JMOL), JMOL = JMOL1, JMOL10)
      do J = 1, NRHOX
        write(35, '(I5,1P10E12.3)') J, (XNMOL(J, JMOL), JMOL = JMOL1, JMOL10)
      end do
    end do
  end if

  !=====================================================================
  ! MODE 1 or 11: compute n / partition_function for opacity
  !=====================================================================
  if (MODE /= 2 .and. MODE /= 12) then

  ! Convert element number densities to n/U
  do K = 2, NEQUA
    ID = IDEQUA(K)
    if (ID == 100) then
      ! Electrons: divide by electron partition function (= 2)
      ! and translational PF: (2*pi*m_e*kT/h^2)^(3/2) = 2.4148e15 * T^(3/2)
      do J = 1, NRHOX
        XNZ(J, K) = XNZ(J, K) / SAHA_PREFAC / T(J) / sqrt(T(J))
      end do
    else
      ! Atoms: divide by partition function and translational PF
      do J = 1, NRHOX
        call PFSAHA(J, ID, 1, 3, FRAC)
        XNZ(J, K) = XNZ(J, K) / FRAC(J, 1) / 1.8786D20 &
                   / sqrt((ATMASS(ID) * T(J))**3)
      end do
    end if
  end do

  ! Compute n/U for molecules and track element IDs for electron diagnostics
  IZ = 1
  IDZ(IZ) = 1
  do JMOL = 1, NUMMOL
    NCOMP = LOCJ(JMOL + 1) - LOCJ(JMOL)

    if (EQUIL(1, JMOL) /= 0.0D0) then
      ! True molecule: n_mol/U = exp(D0/kT) * product(n_k/U_k) * translational_PF
      AMASS = 0.0D0
      LOCJ1 = LOCJ(JMOL)
      LOCJ2 = LOCJ(JMOL + 1) - 1

      do J = 1, NRHOX
        XNFPMOL(J, JMOL) = exp(EQUIL(1, JMOL) / TKEV(J))
      end do

      do LOCK = LOCJ1, LOCJ2
        K = KCOMPS(LOCK)
        if (K == NEQUA1) then
          do J = 1, NRHOX
            XNFPMOL(J, JMOL) = XNFPMOL(J, JMOL) / XNZ(J, NEQUA)
          end do
        else
          ID = IDEQUA(K)
          if (ID < 100) AMASS = AMASS + ATMASS(ID)
          do J = 1, NRHOX
            XNFPMOL(J, JMOL) = XNFPMOL(J, JMOL) * XNZ(J, K)
          end do
        end if
      end do

      do J = 1, NRHOX
        XNFPMOL(J, JMOL) = XNFPMOL(J, JMOL) * 1.8786D20 * sqrt((AMASS * T(J))**3)
      end do

    else
      ! Atom/ion: use PFSAHA to get n/U
      ID = int(XNMOLCODE(JMOL))
      if (ID /= IDZ(IZ)) then
        if (JMOL >= 3) NION_MOL(IZ) = LOCJ(JMOL - 1) - LOCJ(JMOL - 2)
        IZ = IZ + 1
        IDZ(IZ) = ID
      end if

      do J = 1, NRHOX
        call PFSAHA(J, ID, NCOMP, 3, FRAC)
        XNFPMOL(J, JMOL) = XNMOL(J, JMOL) / FRAC(J, 1)
      end do
    end if
  end do

  ! Compute electron contribution per element
  NZ = IZ
  do IZ = 1, NZ
    do J = 1, NRHOX
      call PFSAHA(J, IDZ(IZ), NION_MOL(IZ), 4, E_CONTRIB)
      XE_CONTRIB(J, IZ) = E_CONTRIB(J, 1) * XNATOM(J) &
                         * XABUND(J, IDZ(IZ)) / XNE(J)
    end do
  end do

  end if  ! MODE /= 2 and MODE /= 12

  ! --- Print n/U and electron contribution diagnostics ---
  if (ITER == NUMITS.AND. IDEBUG == 1) then
    NN = ((NUMMOL / 10) + 1) * 10
    do JMOL1 = 1, NN, 10
      JMOL10 = JMOL1 + 9
      write(6, '("1",40X,"NUMBER DENSITIES / PARTITION FUNCTIONS"/5X,10F12.2/(I5,1P10E12.3))') &
        (XNMOLCODE(JMOL), JMOL = JMOL1, JMOL10), &
        (J, (XNFPMOL(J, JMOL), JMOL = JMOL1, JMOL10), J = 1, NRHOX)
    end do

    do IZ1 = 1, NZ, 10
      IZ10 = IZ1 + 9
      write(6, '("1",40X,"ELECTRON CONTRIBUTION Ne(elem)/Ne_tot"/5X,10I12/(I5,1P10E12.3))') &
        (IDZ(IZ), IZ = IZ1, IZ10), &
        (J, (XE_CONTRIB(J, IZ), IZ = IZ1, IZ10), J = 1, NRHOX)
    end do
  end if

  return

contains

  !-------------------------------------------------------------------
  ! Compute molecular contribution to energy density (IFEDNS=1 path).
  ! Called from CONVEC via COMPUTE_ONE_POP when energy density needed.
  !-------------------------------------------------------------------
  subroutine nmolec_energy_density()
    integer :: J_e, JMOL_e, NCOMP_e, ION_e, ID_e
    integer :: LOCK_e, LOCJ1_e, LOCJ2_e, K_e
    real*8  :: TPLUS_e, TMINUS_e

    ! Initialize with kinetic energy
    do J_e = 1, NRHOX
      XNTOT = P(J_e) / TK(J_e)
      EDENS(J_e) = 1.5D0 * XNTOT * TK(J_e)
    end do

    do JMOL_e = 1, NUMMOL
      NCOMP_e = LOCJ(JMOL_e + 1) - LOCJ(JMOL_e)

      if (EQUIL(1, JMOL_e) /= 0.0D0) then
        !---------------------------------------------------------------
        ! Molecules: partition function derivative by finite differences
        !---------------------------------------------------------------
        do J_e = 1, NRHOX
          if (XNMOL(J_e, JMOL_e) <= 0.0D0) cycle
          TPLUS_e  = T(J_e) * 1.001D0
          TMINUS_e = T(J_e) * 0.999D0
          PFMIN(J_e)  = 0.0D0
          PFPLUS(J_e) = 0.0D0

          if (XNMOLCODE(JMOL_e) == 101.0D0) then
            ! H2: use dedicated partition function
            PFMIN(J_e)  = PARTFNH2(TMINUS_e)
            PFPLUS(J_e) = PARTFNH2(TPLUS_e)
          else
            ! General molecule: PF from equilibrium constant polynomial
            if (T(J_e) <= 10000.0D0) then
              PFPLUS(J_e) = exp(-EQUIL(2, JMOL_e) &
                + (EQUIL(3, JMOL_e) + (-EQUIL(4, JMOL_e) &
                + (EQUIL(5, JMOL_e) - EQUIL(6, JMOL_e) * TPLUS_e) &
                * TPLUS_e) * TPLUS_e) * TPLUS_e) + 1.0D-30
              PFMIN(J_e) = exp(-EQUIL(2, JMOL_e) &
                + (EQUIL(3, JMOL_e) + (-EQUIL(4, JMOL_e) &
                + (EQUIL(5, JMOL_e) - EQUIL(6, JMOL_e) * TMINUS_e) &
                * TMINUS_e) * TMINUS_e) * TMINUS_e) + 1.0D-30
            end if
          end if
        end do

        ! For non-H2 molecules: multiply PF by component atom PFs at T+/T-
        if (XNMOLCODE(JMOL_e) /= 101.0D0) then
          LOCJ1_e = LOCJ(JMOL_e)
          LOCJ2_e = LOCJ(JMOL_e + 1) - 1
          do LOCK_e = LOCJ1_e, LOCJ2_e
            K_e = KCOMPS(LOCK_e)
            if (K_e == NEQUA) cycle    ! electron: skip
            if (K_e > NEQUA) exit      ! past last equation
            ID_e = IDEQUA(K_e)
            do J_e = 1, NRHOX
              if (XNMOL(J_e, JMOL_e) <= 0.0D0) cycle
              ! PF at T+
              T(J_e) = T(J_e) * 1.001D0
              TK(J_e) = TK(J_e) * 1.001D0
              TKEV(J_e) = TKEV(J_e) * 1.001D0
              call PFSAHA(J_e, ID_e, 1, 3, FRAC)
              PFPLUS(J_e) = PFPLUS(J_e) * FRAC(J_e, 1)
              ! PF at T-
              T(J_e) = T(J_e) / 1.001D0 * 0.999D0
              TK(J_e) = TK(J_e) / 1.001D0 * 0.999D0
              TKEV(J_e) = TKEV(J_e) / 1.001D0 * 0.999D0
              call PFSAHA(J_e, ID_e, 1, 3, FRAC)
              PFMIN(J_e) = PFMIN(J_e) * FRAC(J_e, 1)
              ! Restore T
              T(J_e) = T(J_e) / 0.999D0
              TK(J_e) = TK(J_e) / 0.999D0
              TKEV(J_e) = TKEV(J_e) / 0.999D0
            end do
          end do
        end if

        ! Accumulate energy density contribution
        if (XNMOLCODE(JMOL_e) == 101.0D0) then
          ! H2: dissociation energy = 36118.11 cm^-1
          do J_e = 1, NRHOX
            if (XNMOL(J_e, JMOL_e) <= 0.0D0) cycle
            EDENS(J_e) = EDENS(J_e) + XNMOL(J_e, JMOL_e) * TK(J_e) &
              * (-36118.11D0 * HCKT(J_e) &
              + (PFPLUS(J_e) - PFMIN(J_e)) / (PFPLUS(J_e) + PFMIN(J_e) + 1.0D-30) &
              * 1000.0D0)
          end do
        else
          ! General molecule
          do J_e = 1, NRHOX
            if (XNMOL(J_e, JMOL_e) <= 0.0D0) cycle
            EDENS(J_e) = EDENS(J_e) + XNMOL(J_e, JMOL_e) * TK(J_e) &
              * (-EQUIL(1, JMOL_e) / TKEV(J_e) &
              + (PFPLUS(J_e) - PFMIN(J_e)) / (PFPLUS(J_e) + PFMIN(J_e) + 1.0D-30) &
              * 1000.0D0)
          end do
        end if

      else
        !---------------------------------------------------------------
        ! Atoms: partition function derivative by finite differences
        !---------------------------------------------------------------
        ID_e = int(XNMOLCODE(JMOL_e))
        do J_e = 1, NRHOX
          if (XNMOL(J_e, JMOL_e) <= 0.0D0) cycle
          ! PF at T+
          T(J_e) = T(J_e) * 1.001D0
          TK(J_e) = TK(J_e) * 1.001D0
          TKEV(J_e) = TKEV(J_e) * 1.001D0
          call PFSAHA(J_e, ID_e, NCOMP_e, 5, PFP)
          ! PF at T-
          T(J_e) = T(J_e) / 1.001D0 * 0.999D0
          TK(J_e) = TK(J_e) / 1.001D0 * 0.999D0
          TKEV(J_e) = TKEV(J_e) / 1.001D0 * 0.999D0
          call PFSAHA(J_e, ID_e, NCOMP_e, 5, PFM)
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
        end do

      end if
    end do  ! JMOL_e

    ! Convert to energy per unit mass
    do J_e = 1, NRHOX
      EDENS(J_e) = EDENS(J_e) / RHO(J_e)
    end do

  end subroutine nmolec_energy_density

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

  implicit none

  ! --- Species code tables for Z=1-30 ---
  ! CODE = IZ + NION/100:  IZ = atomic number, NION = number of ion stages
  ! MODE=12 (XNF): for line opacity — fewer stages needed
  ! MODE=11 (XNFP): for continuum opacity — more stages for transition metals

  ! Number of ion stages per element for MODE=12 (XNF, line opacity)
  integer, parameter :: NION_XNF(30) = (/ &
    1, 2, 3, 3, 3,                        &  ! Z = 1-5  (H, He, Li, Be, B)
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,     &  ! Z = 6-16
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  &  ! Z = 17-28
    2, 2 /)                                   ! Z = 29-30

  ! Number of ion stages per element for MODE=11 (XNFP, continuum opacity)
  integer, parameter :: NION_XNFP(30) = (/ &
    1, 2, 3, 3, 3,                         &  ! Z = 1-5
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,      &  ! Z = 6-16
    5, 4, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9,   &  ! Z = 17-28
    2, 2 /)                                    ! Z = 29-30

  ! NELION offsets: IZ*(IZ+1)/2 for Z=1-30
  integer, parameter :: NOFF(30) = (/ &
      1,   3,   6,  10,  15,  21,  28,  36,  45,  55, &
     66,  78,  91, 105, 120, 136, 153, 171, 190, 210, &
    231, 253, 276, 300, 325, 351, 378, 406, 435, 465 /)

  ! --- Molecular species for IFMOL=1 path ---
  ! CODE → NELION mapping for molecules computed by NMOLEC
  ! H2=841, CH=846, NH=847, OH=848, MgH=851, SiH=853,
  ! CaH=858, CrH=862, FeH=864, C2=868, CN=869, CO=870,
  ! SiO=889, TiO=895, VO=896, H2O=940
  integer, parameter :: NMOL_SPECIES = 16
  real*8,  parameter :: MOL_CODE(16) = (/ &
    101.0D0, 106.0D0, 107.0D0, 108.0D0, 112.0D0, 114.0D0, &
    120.0D0, 124.0D0, 126.0D0, 606.0D0, 607.0D0, 608.0D0, &
    814.0D0, 822.0D0, 823.0D0, 10108.0D0 /)
  integer, parameter :: MOL_NELION(16) = (/ &
    841, 846, 847, 848, 851, 853, &
    858, 862, 864, 868, 869, 870, &
    889, 895, 896, 940 /)

  ! --- Local variables ---
  real*8  :: CODE
  integer :: IZ, J, IMOL

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COMPUTE_ALL_POPS'

  !=====================================================================
  ! Pass 1: XNF (MODE=12) — number densities for line opacity
  !=====================================================================
  do IZ = 1, 30
    CODE = dble(IZ) + dble(NION_XNF(IZ)) / 100.0D0
    call COMPUTE_ONE_POP(CODE, 12, XNF(1, NOFF(IZ)))
  end do

  !=====================================================================
  ! Pass 2: XNFP (MODE=11) — n/partition function for continuum opacity
  !=====================================================================
  do IZ = 1, 30
    CODE = dble(IZ) + dble(NION_XNFP(IZ)) / 100.0D0
    call COMPUTE_ONE_POP(CODE, 11, XNFP(1, NOFF(IZ)))
  end do

  !=====================================================================
  ! Elements Z=31-99: 2 stages for both passes
  !=====================================================================
  do IZ = 31, 99
    CODE = dble(IZ) + 0.02D0
    call COMPUTE_ONE_POP(CODE, 11, XNFP(1, 496 + (IZ - 31) * 5))
    call COMPUTE_ONE_POP(CODE, 12, XNF(1, 496 + (IZ - 31) * 5))
  end do

  !=====================================================================
  ! Approximate H2 and CO for cool stars when molecules are off
  ! (analytic equilibrium fits, valid for T < 9000 K)
  !=====================================================================
  do J = 1, NRHOX
    if (T(J) > 9000.0D0) cycle

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
  end do

  !=====================================================================
  ! Molecular populations from NMOLEC (when molecular equilibrium is on)
  !=====================================================================
  if (IFMOL == 0) return

  do IMOL = 1, NMOL_SPECIES
    call COMPUTE_ONE_POP(MOL_CODE(IMOL), 1, XNFP(1, MOL_NELION(IMOL)))
  end do

  return

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

  implicit none

  ! --- Local arrays: finite-difference energy density and density ---
  real*8 :: EDENS1(kw), EDENS2(kw)   ! E at T+, T-
  real*8 :: EDENS3(kw), EDENS4(kw)   ! E at P+, P-
  real*8 :: RHO1(kw),   RHO2(kw)     ! ρ at T+, T-
  real*8 :: RHO3(kw),   RHO4(kw)     ! ρ at P+, P-
  real*8 :: SAVXNE(kw),  SAVXNA(kw),  SAVRHO(kw)  ! saved state
  real*8 :: DILUT(kw)                 ! dilution factor 1 - exp(-τ_Ross)

  real*8 :: DTDRHX(kw)    ! dT/d(RHOX) from DERIV
  real*8 :: ABCONV(kw)    ! convective opacity (harmonic mean at T±ΔT)
  real*8 :: DELTAT(kw)    ! temperature excess of convective element
  real*8 :: ROSST(kw)     ! Rosseland opacity at actual T, P
  real*8 :: CNVINT(kw)    ! integrated convective flux (for overshooting)
  real*8 :: DELHGT(kw)    ! overshooting height increment

  ! --- Local scalars: thermodynamic derivatives ---
  real*8  :: DEDT          ! (dE/dT)_P
  real*8  :: DRDT          ! (dρ/dT)_P
  real*8  :: DEDPG         ! (dE/dP)_T
  real*8  :: DRDPG         ! (dρ/dP)_T
  real*8  :: DPDPG         ! dP_total/dP_gas (= 1 when ignoring P_turb)
  real*8  :: DPDT          ! dP_rad/dT ≈ 4*P_rad/T * dilution
  real*8  :: HEATCV        ! specific heat at constant volume

  ! --- Local scalars: mixing length theory ---
  real*8  :: DEL           ! superadiabatic excess ∇ - ∇_ad
  real*8  :: VCO           ! convective velocity scale
  real*8  :: FLUXCO        ! convective flux coefficient
  real*8  :: D             ! radiative damping parameter
  real*8  :: TAUB          ! optical thickness of convective bubble
  real*8  :: DDD           ! discriminant for cubic equation
  real*8  :: DELTA         ! convective efficiency parameter
  real*8  :: TERM, UP, DOWN  ! series expansion variables
  real*8  :: DPLUS, DMINUS ! opacity ratios at T±ΔT
  real*8  :: OLDDELT       ! previous iteration's DELTAT
  integer :: ITS30         ! max opacity iterations (30 or 1)
  integer :: IDELTAT       ! opacity iteration counter

  ! --- Local scalars: overshooting ---
  real*8  :: WTCNV         ! overshooting weight
  real*8  :: CNV1(1), CNV2(1)    ! interpolated integrated fluxes
  real*8  :: XNEW_CNV(1)
  integer :: M             ! MAP1 return value

  ! --- Other locals ---
  integer :: J
  integer :: JTOP, JBOT, JA, JB  ! gap-filling indices
  real*8  :: WGHT                 ! interpolation weight

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING CONVEC'

  call DERIV(RHOX, T, DTDRHX, NRHOX)

  !=====================================================================
  ! Finite-difference thermodynamic derivatives
  ! Compute E and ρ at T±0.1% and P±0.1% to get dE/dT, dρ/dT, etc.
  !=====================================================================
  IFEDNS = 1

  ! Save current state
  do J = 1, NRHOX
    DILUT(J) = 1.0D0 - exp(-TAUROS(J))
    ! DILUT(J) = PRAD(J) / PRADK(J)
    SAVXNE(J) = XNE(J)
    SAVXNA(J) = XNATOM(J)
    SAVRHO(J) = RHO(J)
  end do

  ! --- Perturbation 1: T + 0.1% ---
  do J = 1, NRHOX
    TLOG(J) = TLOG(J) + 0.0009995003D0
    T(J)    = T(J) * 1.001D0
    TK(J)   = TK(J) * 1.001D0
    HKT(J)  = HKT(J) / 1.001D0
    HCKT(J) = HCKT(J) / 1.001D0
    TKEV(J) = TKEV(J) * 1.001D0
  end do
  ITEMP = ITEMP + 1
  call COMPUTE_ONE_POP(0.0D0, 1, XNE)

  do J = 1, NRHOX
    ! 3*PRADK is approximately RADEN (radiation energy density).
    ! PRADK is used because it can be reconstructed from model decks
    ! whereas RADEN cannot. Rigorously the radiation field should be
    ! recalculated.
    EDENS1(J) = EDENS(J) + 3.0D0 * PRADK(J) / RHO(J) &
              * (1.0D0 + DILUT(J) * (1.001D0**4 - 1.0D0))
    RHO1(J) = RHO(J)
  end do

  ! --- Perturbation 2: T - 0.1% ---
  do J = 1, NRHOX
    TLOG(J) = TLOG(J) - 0.0009995003D0 - 0.0010005003D0
    T(J)    = T(J) / 1.001D0 * 0.999D0
    TK(J)   = TK(J) / 1.001D0 * 0.999D0
    HKT(J)  = HKT(J) * 1.001D0 / 0.999D0
    HCKT(J) = HCKT(J) * 1.001D0 / 0.999D0
    TKEV(J) = TKEV(J) / 1.001D0 * 0.999D0
  end do
  ITEMP = ITEMP + 1
  call COMPUTE_ONE_POP(0.0D0, 1, XNE)

  do J = 1, NRHOX
    EDENS2(J) = EDENS(J) + 3.0D0 * PRADK(J) / RHO(J) &
              * (1.0D0 + DILUT(J) * (0.999D0**4 - 1.0D0))
    RHO2(J) = RHO(J)
  end do

  ! --- Perturbation 3: P + 0.1% (restore T first) ---
  do J = 1, NRHOX
    TLOG(J) = TLOG(J) + 0.0010005003D0
    T(J)    = T(J) / 0.999D0
    TK(J)   = TK(J) / 0.999D0
    HKT(J)  = HKT(J) * 0.999D0
    HCKT(J) = HCKT(J) * 0.999D0
    TKEV(J) = TKEV(J) / 0.999D0
    P(J)    = P(J) * 1.001D0
  end do
  ITEMP = ITEMP + 1
  call COMPUTE_ONE_POP(0.0D0, 1, XNE)

  do J = 1, NRHOX
    EDENS3(J) = EDENS(J) + 3.0D0 * PRADK(J) / RHO(J)
    RHO3(J) = RHO(J)
    P(J) = P(J) / 1.001D0 * 0.999D0
  end do

  ! --- Perturbation 4: P - 0.1% ---
  ITEMP = ITEMP + 1
  call COMPUTE_ONE_POP(0.0D0, 1, XNE)

  do J = 1, NRHOX
    EDENS4(J) = EDENS(J) + 3.0D0 * PRADK(J) / RHO(J)
    RHO4(J) = RHO(J)
    ! Restore saved state
    XNE(J)    = SAVXNE(J)
    XNATOM(J) = SAVXNA(J)
    RHO(J)    = SAVRHO(J)
    P(J)      = P(J) / 0.999D0
  end do

  !=====================================================================
  ! Compute thermodynamic quantities and convective flux at each depth
  !=====================================================================
  do J = 1, NRHOX

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
    if (HEATCV <= 0.0D0) then
      VELSND(J) = 0.0D0
    else
      VELSND(J) = sqrt(max(HEATCP(J) / HEATCV * DPDPG / DRDPG, 0.0D0))
    end if

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
    if (MIXLTH == 0.0D0) cycle
    if (J < 4) cycle

    DEL = DLTDLP(J) - GRDADB(J)
    if (DEL < 0.0D0) cycle    ! subadiabatic: no convection

    VCO = 0.5D0 * MIXLTH &
        * sqrt(max(-0.5D0 * PTOTAL(J) / RHO(J) * DLRDLT(J), 0.0D0))
    if (VCO == 0.0D0) cycle

    FLUXCO = 0.5D0 * RHO(J) * HEATCP(J) * T(J) * MIXLTH / FOURPI
    ROSST(J) = ROSSTAB(T(J), P(J), VTURB(J))
    OLDDELT = 0.0D0

    ! --- Iterate on the convective opacity ---
    ITS30 = 30
    if (IFCONV == 0) ITS30 = 1

    do IDELTAT = 1, ITS30
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
      if (DDD >= 0.5D0) then
        DELTA = (1.0D0 - sqrt(1.0D0 - DDD)) / DDD
      else
        DELTA = 0.5D0
        TERM  = 0.5D0
        UP    = -1.0D0
        DOWN  = 2.0D0
        do
          UP   = UP + 2.0D0
          DOWN = DOWN + 2.0D0
          TERM = UP / DOWN * DDD * TERM
          DELTA = DELTA + TERM
          if (TERM <= 1.0D-6) exit
        end do
      end if
      DELTA = DELTA * DEL**2 / (D + DEL)

      ! Convective velocity and flux
      VCONV(J)  = VCO * sqrt(DELTA)
      FLXCNV(J) = FLUXCO * VCONV(J) * DELTA
      FLXCNV(J) = max(FLXCNV(J), 0.0D0)

      ! Temperature excess of convective element
      DELTAT(J) = T(J) * MIXLTH * DELTA
      DELTAT(J) = min(DELTAT(J), T(J) * 0.15D0)
      DELTAT(J) = DELTAT(J) * 0.7D0 + OLDDELT * 0.3D0

      if (DELTAT(J) < OLDDELT + 0.5D0 .and. &
          DELTAT(J) > OLDDELT - 0.5D0) exit
      OLDDELT = DELTAT(J)
    end do  ! IDELTAT
    if (IDELTAT > ITS30) then
      write(6,'(A,I3,A,F8.1,A,F8.1)') &
        ' CONVEC WARNING: opacity iteration did not converge at layer ', &
        J, '  DELTAT=', DELTAT(J), '  OLDDELT=', OLDDELT
    end if

  end do  ! J depth loop

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
  do J = 1, NRHOX
    FLXCNV0(J) = FLXCNV(J)
  end do
  do J = 2, NRHOX - 1
    FLXCNV(J) = 0.25D0 * FLXCNV0(J - 1) + 0.50D0 * FLXCNV0(J) &
              + 0.25D0 * FLXCNV0(J + 1)
  end do
  ! Asymmetric boundary kernel at deepest layer: 75-25 split
  FLXCNV(NRHOX) = 0.75D0 * FLXCNV0(NRHOX) + 0.25D0 * FLXCNV0(NRHOX - 1)

  ! Fill radiative gaps inside the convection zone.
  ! If layers above and below are both convective, the layer in between
  ! must be too — a radiative pocket inside a convection zone is unphysical.
  ! Find the shallowest and deepest convective layers, then fill any
  ! zero-flux layers between them by interpolating from the boundaries.
  JTOP = 0
  JBOT = 0
  do J = 1, NRHOX
    if (FLXCNV(J) > 0.0D0) then
      if (JTOP == 0) JTOP = J
      JBOT = J
    end if
  end do
  if (JTOP > 0 .and. JBOT > JTOP + 1) then
    do J = JTOP + 1, JBOT - 1
      if (FLXCNV(J) == 0.0D0) then
        ! Linear interpolation in log between nearest convective neighbors
        ! Find nearest convective layer above
        JA = J - 1
        do while (JA > JTOP .and. FLXCNV(JA) == 0.0D0)
          JA = JA - 1
        end do
        ! Find nearest convective layer below
        JB = J + 1
        do while (JB < JBOT .and. FLXCNV(JB) == 0.0D0)
          JB = JB + 1
        end do
        if (FLXCNV(JA) > 0.0D0 .and. FLXCNV(JB) > 0.0D0) then
          WGHT = dble(J - JA) / dble(JB - JA)
          FLXCNV(J) = FLXCNV(JA) * (1.0D0 - WGHT) + FLXCNV(JB) * WGHT
        end if
      end if
    end do
  end if

  do J = 1, NRHOX
    FLXCNV0(J) = FLXCNV(J)
  end do

  !=====================================================================
  ! Overshooting: extend convective flux above the formal boundary
  ! Assumes overshooting by 0.5 H_P if convection is strong,
  ! none if weak. Setting OVERWT=0 turns off overshooting entirely.
  !=====================================================================
  if (OVERWT > 0.0D0) then
    ! Find maximum convective-to-total flux ratio
    !     WTCNV = MIN(FLXCNV(NRHOX)/FLUX, 1.D0) * OVERWT
    ! Correction from Fiorella Castelli:
    WTCNV = 0.0D0
    do J = 1, NRHOX
      WTCNV = max(WTCNV, FLXCNV(J) / FLUX)
    end do
    WTCNV = min(WTCNV, 1.0D0) * OVERWT

    do J = 1, NRHOX
      !      DELHGT(J) = MIN(HSCALE(J)*MIXLTH*0.5D-5, HEIGHT(NRHOX)-HEIGHT(J),
      DELHGT(J) = min(HSCALE(J) * 0.5D-5 * WTCNV, &
                      HEIGHT(NRHOX) - HEIGHT(J), &
                      HEIGHT(J) - HEIGHT(1))
      !      WRITE(6,775) J, HEIGHT(J), DELHGT(J), CNVINT(J)
      FLXCNV0(J) = FLXCNV(J)
      FLXCNV1(J) = 0.0D0
    end do

    call INTEG(HEIGHT, FLXCNV, CNVINT, NRHOX, 0.0D0)

    do J = NRHOX / 2, NRHOX - 1
      if (DELHGT(J) == 0.0D0) cycle
      XNEW_CNV(1) = HEIGHT(J) - DELHGT(J)
      M = MAP1(HEIGHT, CNVINT, NRHOX, XNEW_CNV, CNV1, 1)
      XNEW_CNV(1) = HEIGHT(J) + DELHGT(J)
      M = MAP1(HEIGHT, CNVINT, NRHOX, XNEW_CNV, CNV2, 1)
      FLXCNV1(J) = FLXCNV1(J) + (CNV2(1) - CNV1(1)) / DELHGT(J) / 2.0D0
    end do

    do J = 1, NRHOX
      FLXCNV(J) = max(FLXCNV0(J), FLXCNV1(J))
    end do
  end if

  ! Patch to remove numerical artifacts in outermost layers
  do J = 1, NCONV
    ! DO 7779 J=1,NRHOX/3
    FLXCNV(J) = 0.0D0
  end do

  return

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

  implicit none

  real*8 :: RHOINV(kw)
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COMPUTE_HEIGHT'

  do J = 1, NRHOX
    RHOINV(J) = 1.0D-5 / RHO(J)
  end do

  call INTEG(RHOX, RHOINV, HEIGHT, NRHOX, 0.0D0)

  return

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

  implicit none

  ! --- Arguments ---
  real*8, intent(in) :: VNEW

  ! --- Avrett Solar Model C microturbulence profile ---
  ! 30-point tabulation of v_turb (cm/s) vs log10(tau_Ross)
  integer, parameter :: NSOLAR = 30
  real*8, parameter :: VTURB_SOLAR_MAX = 1.83D5  ! peak value (cm/s)

  real*8, parameter :: VSTANDARD(30) = (/ &
    0.50D5, 0.50D5, 0.50D5, 0.51D5, 0.52D5, 0.55D5, 0.63D5, &
    0.80D5, 0.90D5, 1.00D5, 1.10D5, 1.20D5, 1.30D5, 1.40D5, &
    1.46D5, 1.52D5, 1.56D5, 1.60D5, 1.64D5, 1.68D5, 1.71D5, &
    1.74D5, 1.76D5, 1.78D5, 1.80D5, 1.81D5, 1.82D5, 1.83D5, &
    1.83D5, 1.83D5 /)

  real*8, parameter :: TAUSTANDARD(30) = (/ &
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
  integer, parameter :: NG_TAB = 13, NT_TAB = 25
  real*8, parameter :: VMAXSTD(NG_TAB, NT_TAB) = reshape( (/ &
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
  real*8  :: VMAX           ! amplitude scale for the profile (cm/s)
  real*8  :: TAULOG(kw)     ! log10(tau_Ross) at each depth
  real*8  :: DELG, DELT     ! bilinear interpolation weights
  integer :: IG, IT         ! table indices
  integer :: J, MM

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING VTURB_VARYDEPTH'

  !---------------------------------------------------------------------
  ! Determine VMAX: either from table or from user input
  !---------------------------------------------------------------------
  if (VNEW == -99.0D5) then
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
  else
    VMAX = abs(VNEW)
  end if

  !---------------------------------------------------------------------
  ! Build depth-dependent profile by interpolating the solar template
  ! onto the model's optical depth grid, then scaling by VMAX
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    TAULOG(J) = log10(TAUSTD(J))
  end do

  MM = MAP1(TAUSTANDARD, VSTANDARD, NSOLAR, TAULOG, VTURB, NRHOX)

  do J = 1, NRHOX
    VTURB(J) = VTURB(J) * VMAX / VTURB_SOLAR_MAX
  end do

  return

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

  implicit none

  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COMPUTE_PTURB'

  do J = 1, NRHOX
    VTURB(J) = (TRBFDG * RHO(J)**TRBPOW &
              + TRBSND * VELSND(J) / 1.0D5 &
              + TRBCON) * 1.0D5
    PTURB(J) = RHO(J) * VTURB(J)**2 * 0.5D0
  end do

  return

END SUBROUTINE COMPUTE_PTURB

!=========================================================================
! SUBROUTINE KAPP
!
! Assemble total monochromatic opacity from all continuum and line sources.
!
! At the current frequency (set by BNU, FREQ, etc. before calling KAPP),
! this routine:
!   1. Zeros all per-source opacity arrays
!   2. Calls each enabled opacity source (controlled by IFOP flags):
!        IFOP(1)  HOP     — hydrogen bound-free and free-free
!        IFOP(2)  H2PLOP  — H2+ opacity
!        IFOP(3)  HMINOP  — H- bound-free and free-free
!        IFOP(4)  HRAYOP  — hydrogen Rayleigh scattering
!        IFOP(5)  HE1OP   — neutral helium bound-free
!        IFOP(6)  HE2OP   — ionized helium bound-free and free-free
!        IFOP(7)  HEMIOP  — He- free-free
!        IFOP(8)  HERAOP  — helium Rayleigh scattering
!        IFOP(9)  COOLOP  — metals (Mg, Al, Si, Ca, Fe, etc.)
!        IFOP(10) WARMOP  — additional metal opacities
!        IFOP(11) HOTOP   — high-ionization opacities
!        IFOP(12) ELECOP  — electron (Thomson) scattering
!        IFOP(13) H2RAOP  — H2 Rayleigh scattering
!        IFOP(14) HLINOP  — hydrogen line opacity (Stark broadened)
!        IFOP(19) XCONOP  — extra continuum opacity (user-defined)
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

  implicit none

  integer :: J
  real*8  :: A         ! sum of non-H, non-Hminus, non-extra continuum absorption

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING KAPP'

  ! Zero all per-source opacity arrays
  do J = 1, NRHOX
    AHYD(J)   = 0.0D0
    AH2P(J)   = 0.0D0
    AHMIN(J)  = 0.0D0
    SIGH(J)   = 0.0D0
    AHE1(J)   = 0.0D0
    AHE2(J)   = 0.0D0
    AHEMIN(J) = 0.0D0
    SIGHE(J)  = 0.0D0
    ACOOL(J)  = 0.0D0
    ALUKE(J)  = 0.0D0
    AHOT(J)   = 0.0D0
    SIGEL(J)  = 0.0D0
    SIGH2(J)  = 0.0D0
    AHLINE(J) = 0.0D0
    ALINES(J) = 0.0D0
    SIGLIN(J) = 0.0D0
    AXLINE(J) = 0.0D0
    SIGXL(J)  = 0.0D0
    AXCONT(J) = 0.0D0
    SIGX(J)   = 0.0D0
    SHYD(J)   = 0.0D0
    SHMIN(J)  = 0.0D0
    SHLINE(J) = 0.0D0
    SXLINE(J) = 0.0D0
    SXCONT(J) = 0.0D0
    SAL1(J)   = 0.0D0
    SFE1(J)   = 0.0D0
    SHE2(J)   = 0.0D0
    SC1(J)    = 0.0D0
  end do

  ! Call each enabled opacity source
  if (IFOP(1)  == 1) call HOP
  if (IFOP(2)  == 1) call H2PLOP
  if (IFOP(3)  == 1) call HMINOP
  if (IFOP(4)  == 1) call HRAYOP
  if (IFOP(5)  == 1) call HE1OP
  if (IFOP(6)  == 1) call HE2OP
  if (IFOP(7)  == 1) call HEMIOP
  if (IFOP(8)  == 1) call HERAOP
  if (IFOP(9)  == 1) call COOLOP
  if (IFOP(10) == 1) call WARMOP
  if (IFOP(11) == 1) call HOTOP
  if (IFOP(12) == 1) call ELECOP
  if (IFOP(13) == 1) call H2RAOP
  if (IFOP(14) == 1) call HLINOP
  if (IFOP(19) == 1) call XCONOP

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
  ! Note on the absorption sum: ACONT here uses the individual metal
  ! opacity arrays (AC1, AMG1, AAL1, ASI1, AFE1) directly rather than
  ! the lumped ACOOL bundle from atlas12.for, because the F90 broke
  ! the metal continua out into separate routines during translation.
  ! The numerical sum is identical -- just a regrouping.
  do J = 1, NRHOX
    ! Sources weighted by BNU in the source function
    A = AH2P(J) + AHE1(J) + AHE2(J) + AHEMIN(J) + ALUKE(J) + AHOT(J) &
      + AC1(J)  + AMG1(J) + AAL1(J) + ASI1(J)   + AFE1(J)

    ! Total continuum absorption
    ACONT(J) = A + AHYD(J) + AHMIN(J) + AXCONT(J)

    ! Continuum source function (atlas12.for form)
    SCONT(J) = BNU(J)
    if (ACONT(J) > 0.0D0) then
      SCONT(J) = (A * BNU(J) + AHYD(J) * SHYD(J) &
               + AHMIN(J) * SHMIN(J) + AXCONT(J) * SXCONT(J)) / ACONT(J)
    end if

    ! Line absorption and source function
    ALINE(J) = AHLINE(J) + ALINES(J) + AXLINE(J)
    SLINE(J) = BNU(J)
    if (ALINE(J) > 0.0D0) then
      SLINE(J) = (AHLINE(J) * SHLINE(J) + ALINES(J) * BNU(J) &
               + AXLINE(J) * SXLINE(J)) / ALINE(J)
    end if

    ! Scattering: continuum and line
    SIGMAC(J) = SIGH(J) + SIGHE(J) + SIGEL(J) + SIGH2(J) + SIGX(J)
    SIGMAL(J) = SIGLIN(J) + SIGXL(J)
  end do

  return

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

  implicit none

  ! Maximum hydrogen level included in HOP.  The explicit-level loop runs
  ! n=1..15 (with full Karzas-Latter cross sections); a high-n extension
  ! loop covers n=16..NMAX_HLEVELS using XKARZAS's hydrogenic rescaling
  ! from the n=15 table; and the dissolved-fraction pseudo-continuum
  ! covers n=7..NMAX_HLEVELS.  All three loops share the same upper bound
  ! so the bookkeeping stays consistent.
  integer, parameter :: NMAX_HLEVELS = 30

  real*8  :: X            ! bound-free cross-section from XKARZAS
  real*8  :: A            ! single-level opacity contribution
  real*8  :: H, S         ! running sums for opacity and source function
  real*8  :: w_n          ! occupation probability for current level
  real*8  :: E_n          ! excitation energy for the current level [cm^-1]
  real*8  :: thresh_n     ! ionization threshold of the current level [cm^-1]
  real*8  :: FREQ3_INV    ! 1 / FREQ^3 (precomputed for pseudo-continuum)
  integer :: J, I_PC, n_hi

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HOP'

  do J = 1, NRHOX

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
    levels: do

      ! n=15  threshold =    487.456 cm^{-1}
      IF (WAVENO < 487.456D0) EXIT levels
      w_n = occupation_prob(15, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 15, 15)
      A = w_n * X * 450.0D0 * EXP(-109191.313D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=14  threshold =    559.579 cm^{-1}
      IF (WAVENO < 559.579D0) EXIT levels
      w_n = occupation_prob(14, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 14, 14)
      A = w_n * X * 392.0D0 * EXP(-109119.188D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=13  threshold =    648.980 cm^{-1}
      IF (WAVENO < 648.980D0) EXIT levels
      w_n = occupation_prob(13, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 13, 13)
      A = w_n * X * 338.0D0 * EXP(-109029.789D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=12  threshold =    761.649 cm^{-1}
      IF (WAVENO < 761.649D0) EXIT levels
      w_n = occupation_prob(12, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 12, 12)
      A = w_n * X * 288.0D0 * EXP(-108917.117D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=11  threshold =    906.426 cm^{-1}
      IF (WAVENO < 906.426D0) EXIT levels
      w_n = occupation_prob(11, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 11, 11)
      A = w_n * X * 242.0D0 * EXP(-108772.336D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=10  threshold =   1096.776 cm^{-1}
      IF (WAVENO < 1096.776D0) EXIT levels
      w_n = occupation_prob(10, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 10, 10)
      A = w_n * X * 200.0D0 * EXP(-108581.992D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=9   threshold =   1354.044 cm^{-1}
      IF (WAVENO < 1354.044D0) EXIT levels
      w_n = occupation_prob(9, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 9, 9)
      A = w_n * X * 162.0D0 * EXP(-108324.719D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=8   threshold =   1713.713 cm^{-1}
      IF (WAVENO < 1713.713D0) EXIT levels
      w_n = occupation_prob(8, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 8, 8)
      A = w_n * X * 128.0D0 * EXP(-107965.051D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=7   threshold =   2238.320 cm^{-1}
      IF (WAVENO < 2238.320D0) EXIT levels
      w_n = occupation_prob(7, XNE(J))
      X = XKARZAS(FREQ, 1.0D0, 7, 7)
      A = w_n * X * 98.0D0 * EXP(-107440.444D0 * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)

      ! n=6   threshold =   3046.604 cm^{-1}   (non-LTE)
      IF (WAVENO < 3046.604D0) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 6, 6)
      A = X * 72.0D0 * EXP(-106632.160D0 * HCKT(J)) * (BHYD(J,6) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,6) - EHVKT(J))

      ! n=5   threshold =   4387.113 cm^{-1}   (non-LTE)
      IF (WAVENO < 4387.113D0) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 5, 5)
      A = X * 50.0D0 * EXP(-105291.651D0 * HCKT(J)) * (BHYD(J,5) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,5) - EHVKT(J))

      ! n=4   threshold =   6854.871 cm^{-1}   (non-LTE)
      IF (WAVENO < 6854.871D0) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 4, 4)
      A = X * 32.0D0 * EXP(-102823.893D0 * HCKT(J)) * (BHYD(J,4) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,4) - EHVKT(J))

      ! n=3   threshold =  12186.462 cm^{-1}   (non-LTE)
      IF (WAVENO < 12186.462D0) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 3, 3)
      A = X * 18.0D0 * EXP(-97492.302D0 * HCKT(J)) * (BHYD(J,3) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,3) - EHVKT(J))

      ! n=2   threshold =  27419.659 cm^{-1}   (non-LTE)  [Balmer limit]
      IF (WAVENO < 27419.659D0) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 2, 2)
      A = X * 8.0D0 * EXP(-82259.105D0 * HCKT(J)) * (BHYD(J,2) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,2) - EHVKT(J))

      ! n=1   threshold = 109678.764 cm^{-1}   (non-LTE)  [Lyman limit]
      IF (WAVENO < 109678.764D0) EXIT levels
      X = XKARZAS(FREQ, 1.0D0, 1, 1)
      A = X * 2.0D0 * 1.0D0 * (BHYD(J,1) - EHVKT(J))
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J,1) - EHVKT(J))

      EXIT levels  ! all levels processed
    end do levels

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
    do n_hi = 16, NMAX_HLEVELS
      thresh_n = ELIM_HI / dble(n_hi)**2
      if (WAVENO < thresh_n) cycle              ! below this level's threshold
      w_n = occupation_prob(n_hi, XNE(J))
      if (w_n < 1.0D-6) cycle                   ! fully dissolved -- skip
      X = XKARZAS(FREQ, 1.0D0, n_hi, n_hi)
      E_n = ELIM_HI - RYDBERG_H / dble(n_hi)**2
      A = w_n * X * 2.0D0 * dble(n_hi)**2 * EXP(-E_n * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)
    end do

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
    do I_PC = 7, NMAX_HLEVELS
      w_n = occupation_prob(I_PC, XNE(J))
      if (1.0D0 - w_n < 1.0D-6) cycle   ! fully bound, no dissolved fraction
      E_n = ELIM_HI - RYDBERG_H / dble(I_PC)**2
      A = (1.0D0 - w_n) * 2.815D29 * FREQ3_INV &
        * 2.0D0 / dble(I_PC)**3 * EXP(-E_n * HCKT(J)) * STIM(J)
      H = H + A
      S = S + A * BNU(J)
    end do

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
    IF (H > 0.0D0) SHYD(J) = S / H

  end do

  return

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

  implicit none

  ! --- Arguments ---
  real*8,  intent(in) :: FREQ    ! frequency (Hz)
  real*8,  intent(in) :: ZEFF2   ! effective nuclear charge squared
  integer, intent(in) :: N       ! principal quantum number
  integer, intent(in) :: L       ! orbital quantum number

  ! --- Return value ---
  real*8 :: XKARZAS

  ! --- Tabulated cross-section data (read from files on first call) ---
  ! Tables use REAL*4 as in the original Karzas & Latter tabulation
  integer, parameter :: NPTS = 29   ! energy grid points per level
  integer, parameter :: NMAX = 15   ! max tabulated level for XN
  integer, parameter :: NLMAX = 6   ! max tabulated level for XL (l-resolved)

  real*4, save :: FREQN(NPTS, NMAX)     ! log10(freq/Z^2) grid per level
  real*4, save :: XN(NPTS, NMAX)        ! log10(sigma), level-averaged
  real*4, save :: XL(NPTS, NLMAX, NLMAX) ! log10(sigma), l-resolved
  real*4, save :: EKARZAS(NPTS)          ! energy grid for n>15 scaling
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*4  :: FREQLG             ! log10(freq/Z^2)
  real*4  :: X                  ! interpolated log10(sigma)
  real*4  :: FREQN15(NPTS)     ! frequency grid for n>15 (constructed)
  real*8, parameter :: LN10 = 2.30258509299405D0
  integer :: I

  ! --- First-call initialization: read tables from data files ---
  if (.not. INITIALIZED) then
    call read_karzas_tables()
    INITIALIZED = .true.
  end if

  ! --- Compute cross-section ---
  FREQLG = real(log10(FREQ / ZEFF2))
  XKARZAS = 0.0D0

  if (L >= N .or. N > NLMAX) then
    !-----------------------------------------------------------------
    ! Level-averaged cross-section (L >= N means use averaged tables,
    ! or N > 6 where only averaged tables are available)
    !-----------------------------------------------------------------
    if (N > NMAX) then
      ! N > 15: rescale from N=15 using hydrogenic frequency scaling
      FREQN15(NPTS) = log10(FREQ_RYDH / real(N)**2)
      if (FREQLG < FREQN15(NPTS)) return
      do I = 2, NPTS - 1
        FREQN15(I) = log10((EKARZAS(I) + 1.0 / real(N)**2) &
                   * FREQ_RYDH)
        if (FREQLG > FREQN15(I)) exit
      end do
      if (I > NPTS - 1) I = NPTS
      X = (FREQLG - FREQN15(I)) / (FREQN15(I-1) - FREQN15(I)) &
        * (XN(I-1, NMAX) - XN(I, NMAX)) + XN(I, NMAX)
    else
      ! N = 1-15: direct table lookup
      if (FREQLG < FREQN(NPTS, N)) return
      do I = 2, NPTS
        if (FREQLG > FREQN(I, N)) exit
      end do
      X = (FREQLG - FREQN(I, N)) / (FREQN(I-1, N) - FREQN(I, N)) &
        * (XN(I-1, N) - XN(I, N)) + XN(I, N)
    end if

  else
    !-----------------------------------------------------------------
    ! L-resolved cross-section (N <= 6, L < N)
    !-----------------------------------------------------------------
    if (FREQLG < FREQN(NPTS, N)) return
    do I = 2, NPTS
      if (FREQLG > FREQN(I, N)) exit
    end do
    X = (FREQLG - FREQN(I, N)) / (FREQN(I-1, N) - FREQN(I, N)) &
      * (XL(I-1, L+1, N) - XL(I, L+1, N)) + XL(I, L+1, N)
  end if

  XKARZAS = exp(dble(X) * LN10) / ZEFF2

  return

contains

  !-------------------------------------------------------------------
  ! Read Karzas & Latter cross-section tables from data files
  !-------------------------------------------------------------------
  subroutine read_karzas_tables()
    integer :: I, J, K, IU, IOS
    character(256) :: FILEPATH

    IU = 89  ! scratch unit (same convention as other data reads)

    ! Read EKARZAS (energy grid)
    FILEPATH = trim(DATADIR) // 'karzas_ekarzas.dat'
    open(unit=IU, file=FILEPATH, status='OLD', action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      stop 'XKARZAS: cannot read Karzas-Latter data'
    end if
    read(IU, '(A)') FILEPATH  ! skip header
    read(IU, '(A)') FILEPATH  ! skip header
    do I = 1, NPTS
      read(IU, *) EKARZAS(I)
    end do
    close(IU)

    ! Read XN (level-averaged cross-sections, column-major)
    FILEPATH = trim(DATADIR) // 'karzas_xn.dat'
    open(unit=IU, file=FILEPATH, status='OLD', action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      stop 'XKARZAS: cannot read Karzas-Latter data'
    end if
    read(IU, '(A)') FILEPATH  ! skip header
    read(IU, '(A)') FILEPATH  ! skip header
    read(IU, '(A)') FILEPATH  ! skip header
    do J = 1, NMAX
      do I = 1, NPTS
        read(IU, *) XN(I, J)
      end do
    end do
    close(IU)

    ! Read FREQN (frequency grids, column-major)
    FILEPATH = trim(DATADIR) // 'karzas_freqn.dat'
    open(unit=IU, file=FILEPATH, status='OLD', action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      stop 'XKARZAS: cannot read Karzas-Latter data'
    end if
    read(IU, '(A)') FILEPATH  ! skip header
    read(IU, '(A)') FILEPATH  ! skip header
    do J = 1, NMAX
      do I = 1, NPTS
        read(IU, *) FREQN(I, J)
      end do
    end do
    close(IU)

    ! Read XL (l-resolved cross-sections, column-major)
    FILEPATH = trim(DATADIR) // 'karzas_xl.dat'
    open(unit=IU, file=FILEPATH, status='OLD', action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'XKARZAS: ERROR opening ', trim(FILEPATH)
      stop 'XKARZAS: cannot read Karzas-Latter data'
    end if
    read(IU, '(A)') FILEPATH  ! skip header
    read(IU, '(A)') FILEPATH  ! skip header
    read(IU, '(A)') FILEPATH  ! skip header
    do K = 1, NLMAX
      do J = 1, NLMAX
        do I = 1, NPTS
          read(IU, *) XL(I, J, K)
        end do
      end do
    end do
    close(IU)

  end subroutine read_karzas_tables

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

  implicit none

  ! --- Arguments ---
  integer, intent(in) :: N       ! principal quantum number
  real*8,  intent(in) :: FREQ    ! frequency (Hz)
  real*8,  intent(in) :: Z       ! nuclear charge

  ! --- Return value ---
  real*8 :: COULX

  ! --- Gaunt factor polynomial coefficients for N=1-6 ---
  ! sigma = Kramers * (A + (B + C*(Z^2/nu)) * (Z^2/nu))
  real*8, parameter :: A(6) = (/ 0.9916D0, 1.105D0, 1.101D0, &
                                  1.101D0, 1.102D0, 1.0986D0 /)
  real*8, parameter :: B(6) = (/ 2.719D13, -2.375D14, -9.863D13, &
                                  -5.765D13, -3.909D13, -2.704D13 /)
  real*8, parameter :: C(6) = (/ -2.268D30, 4.077D28, 1.035D28, &
                                  4.593D27, 2.371D27, 1.229D27 /)

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COULX'

  ! Check if frequency is above the ionization threshold: nu_edge = Z^2 * R_inf * c / n^2
  if (FREQ < Z * Z * FREQ_RYDH / dble(N)**2) then
    COULX = 0.0D0
    return
  end if

  ! Kramers cross-section
  COULX = 2.815D29 / FREQ**3 / dble(N)**5 * Z**4

  ! Apply Gaunt factor correction
  if (N <= 6) then
    if (N == 1) then
      ! Exact 1s Gaunt factor from tabulation
      COULX = COULX * COULBF1S(FREQ, Z)
    else
      ! Polynomial fit for N=2-6
      COULX = COULX * (A(N) + (B(N) + C(N) * (Z * Z / FREQ)) * (Z * Z / FREQ))
    end if
  end if

  return

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

  implicit none

  ! --- Arguments ---
  real*8, intent(in) :: FREQ     ! frequency (Hz)
  real*8, intent(in) :: Z        ! nuclear charge

  ! --- Return value ---
  real*8 :: COULBF1S

  ! --- 1s Gaunt factor table ---
  ! 151 points at log10(nu/nu_edge) = 0.00, 0.02, 0.04, ..., 3.00
  integer, parameter :: NGAUNT = 151
  real*8, parameter :: DLOG_STEP = 0.02D0

  real*8, parameter :: GAUNT1S(151) = (/ &
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
  real*8  :: ELOG      ! log10(nu / nu_edge)
  integer :: I         ! table index

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COULBF1S'

  ! Below ionization threshold
  if (FREQ / Z**2 < FREQ_RYDH) then
    COULBF1S = 0.0D0
    return
  end if

  ! Linear interpolation in log10(nu/nu_edge)
  ELOG = log10(FREQ / Z**2 / FREQ_RYDH)
  I = int(ELOG / DLOG_STEP)
  I = max(min(I + 1, NGAUNT - 1), 1)

  COULBF1S = GAUNT1S(I) + (GAUNT1S(I+1) - GAUNT1S(I)) / DLOG_STEP &
           * (ELOG - dble(I - 1) * DLOG_STEP)

  return

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

  implicit none

  ! --- Arguments ---
  integer, intent(in) :: J       ! depth index
  integer, intent(in) :: NZ      ! charge index (1=H, 2=He+, ..., 6)

  ! --- Return value ---
  real*8 :: COULFF

  ! --- log10(Z^2) for NZ = 1-6 (canonical gamma^2 scales as Z^2) ---
  real*8, parameter :: Z2LOG(6) = (/ &
    0.0D0, 0.60206D0, 0.95424D0, 1.20412D0, 1.39794D0, 1.55630D0 /)

  ! --- van Hoof (2014) table grid parameters ---
  integer, parameter :: GFF_N_GAM2       = 81
  integer, parameter :: GFF_N_U          = 146
  real*8,  parameter :: GFF_LOG_GAM2_MIN = -6.0D0
  real*8,  parameter :: GFF_LOG_U_MIN    = -16.0D0
  real*8,  parameter :: GFF_DLOG         =  0.2D0
  integer, parameter :: GFF_MAGIC        = 20140210

  ! --- Cached Gaunt factor table (stored as log10(g_ff) for interpolation) ---
  ! First index: log(gamma^2); second index: log(u).
  real*8, save :: GFF_LOGTAB(GFF_N_GAM2, GFF_N_U)
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*8  :: LGAM2          ! canonical log10(gamma^2)
  real*8  :: LU             ! canonical log10(u)
  real*8  :: XI, XJ         ! fractional grid indices
  real*8  :: FI, FJ         ! fractional parts
  real*8  :: LG00, LG01, LG10, LG11  ! log(g_ff) at 4 grid corners
  real*8  :: LG             ! interpolated log(g_ff)
  integer :: I0, J0         ! lower-corner grid indices (1-based Fortran)

  ! --- Derived physical constants (computed from mod_constants).
  !     log10(Ry_H / k) = log10(1.57881e5 K), with Ry in energy units.
  !     log10(k/h)      = log10(2.08366e10 Hz/K).
  real*8, parameter :: LOG10_RYH_OVER_K = log10(FREQ_RYDH * HOVERK)
  real*8, parameter :: LOG10_K_OVER_H   = -log10(HOVERK)
  real*8, parameter :: LN10             = 2.30258509299405D0

  ! --- First-call initialization ---
  if (.not. INITIALIZED) then
    call read_gauntff_table()
    INITIALIZED = .true.
  end if

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COULFF'

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

  return

contains

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
  subroutine read_gauntff_table()
    integer :: IU, IOS, I, J
    integer :: MAGIC_IN, N_GAM2_IN, N_U_IN
    real*8  :: LOG_GAM2_IN, LOG_U_IN, DLOG_IN
    real*8  :: ROW(GFF_N_GAM2)
    character(256)  :: FILEPATH, LINE
    character(2048) :: DATALINE   ! big enough for 81 %e values (~1300 chars)

    IU = 89
    FILEPATH = trim(DATADIR) // 'gauntff_vanhoof.dat'
    open(unit=IU, file=FILEPATH, status='OLD', action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'COULFF: ERROR opening ', trim(FILEPATH)
      stop 'COULFF: cannot read van Hoof Gaunt factor table'
    end if

    ! --- Skip copyright header (lines starting with '#') and consume the
    !     5 numeric metadata lines (magic, grid dims, start values, step).
    call read_nonblank_line(IU, LINE);  read(LINE, *) MAGIC_IN
    call read_nonblank_line(IU, LINE);  read(LINE, *) N_GAM2_IN, N_U_IN
    call read_nonblank_line(IU, LINE);  read(LINE, *) LOG_GAM2_IN
    call read_nonblank_line(IU, LINE);  read(LINE, *) LOG_U_IN
    call read_nonblank_line(IU, LINE);  read(LINE, *) DLOG_IN

    ! --- Validate header against expected constants ---
    if (MAGIC_IN  /= GFF_MAGIC .or. &
        N_GAM2_IN /= GFF_N_GAM2 .or. &
        N_U_IN    /= GFF_N_U .or. &
        abs(LOG_GAM2_IN - GFF_LOG_GAM2_MIN) > 1.0D-6 .or. &
        abs(LOG_U_IN    - GFF_LOG_U_MIN)    > 1.0D-6 .or. &
        abs(DLOG_IN     - GFF_DLOG)         > 1.0D-6) then
      write(6,*) 'COULFF: ERROR file header mismatch for ', trim(FILEPATH)
      write(6,*) '   magic:  expect ', GFF_MAGIC,        ' got ', MAGIC_IN
      write(6,*) '   n_gam2: expect ', GFF_N_GAM2,       ' got ', N_GAM2_IN
      write(6,*) '   n_u:    expect ', GFF_N_U,          ' got ', N_U_IN
      write(6,*) '   log_gam2_min: expect ', GFF_LOG_GAM2_MIN, ' got ', LOG_GAM2_IN
      write(6,*) '   log_u_min:    expect ', GFF_LOG_U_MIN,    ' got ', LOG_U_IN
      write(6,*) '   dlog:         expect ', GFF_DLOG,         ' got ', DLOG_IN
      stop 'COULFF: van Hoof table header does not match compiled constants'
    end if

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
    do J = 1, GFF_N_U
      call read_nonblank_line(IU, DATALINE)
      read(DATALINE, *, iostat=IOS) (ROW(I), I = 1, GFF_N_GAM2)
      if (IOS /= 0) then
        write(6,*) 'COULFF: error parsing row ', J, ' of Gaunt factor data'
        stop 'COULFF: corrupted Gaunt factor data file'
      end if
      do I = 1, GFF_N_GAM2
        GFF_LOGTAB(I, J) = log10(ROW(I))
      end do
    end do

    close(IU)
    ! --- Uncertainty table (second half of file) is not read; it is only
    !     useful for sanity-checking the published data itself, not for
    !     runtime opacity calculation.

  end subroutine read_gauntff_table

  !-------------------------------------------------------------------
  ! Read one non-blank, non-comment line from unit IU into LINE.
  ! Treats '#' as a line-initial comment marker (per the van Hoof
  ! file format).  Also strips inline '# comment' tails so the caller
  ! can read numeric values with list-directed READ.
  !-------------------------------------------------------------------
  subroutine read_nonblank_line(IU, LINE)
    integer, intent(in) :: IU
    character(*), intent(out) :: LINE
    integer :: IOS, K
    do
      read(IU, '(A)', iostat=IOS) LINE
      if (IOS /= 0) then
        write(6,*) 'COULFF: unexpected EOF reading gauntff_vanhoof.dat'
        stop 'COULFF: corrupted Gaunt factor data file'
      end if
      LINE = adjustl(LINE)
      if (len_trim(LINE) == 0) cycle
      if (LINE(1:1) == '#')    cycle
      K = index(LINE, '#')
      if (K > 0) LINE(K:) = ' '
      return
    end do
  end subroutine read_nonblank_line

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

  implicit none

  real*8  :: FR        ! log cross-section frequency factor (polynomial in FREQLG)
  real*8  :: ES        ! energy factor (polynomial in FREQ), in eV
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING H2PLOP'

  ! No contribution above Lyman limit
  if (FREQ > FREQ_RYDH) return

  ! Frequency-dependent cross-section factor (polynomial in log10(freq))
  FR = -3.0233D3 + (3.7797D2 + (-1.82496D1 &
     + (3.9207D-1 - 3.1672D-3 * FREQLG) * FREQLG) * FREQLG) * FREQLG

  ! Temperature-dependent energy factor (polynomial in freq), in eV
  ES = -7.342D-3 + (-2.409D-15 + (1.028D-30 &
     + (-4.230D-46 + (1.224D-61 - 1.351D-77 * FREQ) * FREQ) &
     * FREQ) * FREQ) * FREQ

  do J = 1, NRHOX
    AH2P(J) = exp(-ES / TKEV(J) + FR) &
            * XNFP(J, 1) * 2.0D0 * BHYD(J, 1) * XNFP(J, 2) &
            / RHO(J) * STIM(J)
  end do

  return

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

  implicit none

  ! --- Bound-free cross-section table ---
  ! 85 points: wavelength (nm) and sigma (1e-18 cm^2)
  ! From Mathisen (1984) after Wishart (1979), Broad & Reinhardt (1976)
  integer, parameter :: NBF = 85

  real*8, parameter :: WBF(85) = (/ &
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

  real*8, parameter :: BF(85) = (/ &
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
  integer, parameter :: NFF_THETA = 11, NFF_WAVE = 22

  real*8, parameter :: THETAFF(11) = (/ &
    0.5D0, 0.6D0, 0.8D0, 1.0D0, 1.2D0, 1.4D0, &
    1.6D0, 1.8D0, 2.0D0, 2.8D0, 3.6D0 /)

  real*8, parameter :: WAVEK(22) = (/ &
    0.50D0, 0.40D0, 0.35D0, 0.30D0, 0.25D0, 0.20D0, 0.18D0, &
    0.16D0, 0.14D0, 0.12D0, 0.10D0, 0.09D0, 0.08D0, 0.07D0, &
    0.06D0, 0.05D0, 0.04D0, 0.03D0, 0.02D0, 0.01D0, 0.008D0, &
    0.006D0 /)

  real*8, parameter :: FF(11, 22) = reshape( (/ &
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
  real*8,  save :: XHMIN(kw)       ! H⁻ number density / rho factor
  real*8,  save :: THETA(kw)       ! 5040/T at each depth
  integer, save :: ITEMP1 = 0      ! cached ITEMP for skip logic

  ! --- One-time initialization arrays ---
  real*8,  save :: WFFLOG(NFF_WAVE)              ! log(wavelength) grid for ff
  real*8,  save :: FFLOG(NFF_WAVE, NFF_THETA)    ! log(cross-section) table
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*8  :: FFTT(NFF_THETA)  ! ff cross-section interpolated to current wavelength
  real*8  :: FFTHETA(1)       ! ff cross-section interpolated to current depth's theta
  real*8  :: WAVELOG(1)       ! log(wavelength) at current frequency
  real*8  :: HMINBF_ARR(1)     ! bound-free cross-section (1e-18 cm^2)
  real*8  :: WAVE_ARR(1)
  real*8  :: HMINFF           ! free-free opacity per gram
  real*8  :: H                ! bound-free opacity per gram
  real*8  :: FFTLOG(1)        ! temporary for LINTER output
  integer :: J, IWAVE, ITHETA, MAXWAVE

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HMINOP'

  !---------------------------------------------------------------------
  ! One-time initialization: precompute log(wavelength) and log(sigma)
  ! grids for free-free interpolation
  !---------------------------------------------------------------------
  if (.not. INITIALIZED) then
    do IWAVE = 1, NFF_WAVE
      ! 91.134 nm is the conversion factor from Bell & Berrington
      WFFLOG(IWAVE) = log(91.134D0 / WAVEK(IWAVE))
      do ITHETA = 1, NFF_THETA
        FFLOG(IWAVE, ITHETA) = log(FF(ITHETA, IWAVE) * 1.0D-26)
      end do
    end do
    INITIALIZED = .true.
  end if

  !---------------------------------------------------------------------
  ! Recompute temperature-dependent quantities if T has changed
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do J = 1, NRHOX
      THETA(J) = THETA_COEFF / T(J)
      ! H⁻ population: Saha equation for H⁻ ↔ H + e⁻
      ! 0.754209 eV = H⁻ electron affinity (Hotop & Lineberger 1985)
      XHMIN(J) = exp(HMINUS_EA / TKEV(J)) &
               / (SAHA_PREFAC * T(J) * sqrt(T(J))) &
               * BMIN(J) * BHYD(J, 1) * XNFP(J, 1) * XNE(J)
    end do
  end if

  !---------------------------------------------------------------------
  ! Free-free: interpolate table to current wavelength, then to each
  ! depth's theta value
  !---------------------------------------------------------------------
  WAVELOG(1) = log(WAVE)
  do ITHETA = 1, NFF_THETA
    call LINTER(WFFLOG, FFLOG(1, ITHETA), NFF_WAVE, WAVELOG(1), FFTLOG(1))
    FFTT(ITHETA) = exp(FFTLOG(1)) / THETAFF(ITHETA) * THETA_COEFF * KBOL
  end do

  !---------------------------------------------------------------------
  ! Bound-free: interpolate cross-section table at current wavelength
  ! (only contributes for lambda < 1643.91 nm, i.e. freq > 1.82365e14)
  !---------------------------------------------------------------------
  HMINBF_ARR(1) = 0.0D0
  if (FREQ > 1.82365D14) then
    WAVE_ARR(1) = WAVE
    MAXWAVE = MAP1(WBF, BF, NBF, WAVE_ARR, HMINBF_ARR, 1)
  end if

  !---------------------------------------------------------------------
  ! Assemble total opacity and source function at each depth
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    call LINTER(THETAFF, FFTT, NFF_THETA, THETA(J), FFTHETA(1))

    ! Free-free opacity per gram
    HMINFF = FFTHETA(1) * XNFP(J, 1) * 2.0D0 * BHYD(J, 1) * XNE(J) / RHO(J)

    ! Bound-free opacity per gram (with non-LTE and stimulated emission)
    H = HMINBF_ARR(1) * 1.0D-18 * (1.0D0 - EHVKT(J) / BMIN(J)) &
      * XHMIN(J) / RHO(J)

    AHMIN(J) = H + HMINFF

    ! Source function: bf weighted by non-LTE, ff in LTE
    SHMIN(J) = (H * BNU(J) * STIM(J) / (BMIN(J) - EHVKT(J)) &
             + HMINFF * BNU(J)) / AHMIN(J)
  end do

  return

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

  implicit none

  ! --- Arguments ---
  integer, intent(in)  :: NOLD
  real*8,  intent(in)  :: XOLD(NOLD), YOLD(NOLD)
  real*8,  intent(in)  :: XNEW
  real*8,  intent(out) :: YNEW

  ! --- Local variables ---
  integer :: I

  ! Find the bracketing interval
  do I = 2, NOLD
    if (XNEW < XOLD(I)) exit
  end do
  if (I > NOLD) I = NOLD

  ! Linear interpolation (or extrapolation at boundaries)
  YNEW = YOLD(I-1) + (YOLD(I) - YOLD(I-1)) &
       / (XOLD(I) - XOLD(I-1)) * (XNEW - XOLD(I-1))

  return

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

  implicit none

  ! FREQ_RYDH and SIGMA_THOMSON from mod_constants

  ! Range 1: scattering amplitude f, ν/ν_L = 0.01 to 0.74 by 0.01
  real*8, parameter :: GAVRILAM(74) = (/ &
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
  real*8, parameter :: GAVRILAMAB(27) = (/ &
   31.008832D0,  15.382871D0,  10.160646D0,   7.538338D0,   5.955062D0, &
    4.890397D0,   4.121176D0,   3.535672D0,   3.071659D0,   2.691623D0, &
    2.371483D0,   2.094936D0,   1.850395D0,   1.629203D0,   1.424526D0, &
    1.230596D0,   1.042127D0,   0.853766D0,   0.659460D0,   0.451533D0, &
    0.219115D0,  -0.054939D0,  -0.400868D0,  -0.879559D0,  -1.637857D0, &
   -3.150374D0,  -8.326078D0 /)

  ! Range 4: near Ly γ, ν/ν_L = 0.890 to 0.936 by 0.002
  real*8, parameter :: GAVRILAMBC(24) = (/ &
   32.260389D0,  11.880702D0,   7.418436D0,   5.442077D0,   4.313409D0, &
    3.573504D0,   3.043218D0,   2.637983D0,   2.312466D0,   2.039959D0, &
    1.803441D0,   1.591244D0,   1.394717D0,   1.206823D0,   1.021148D0, &
    0.831020D0,   0.628449D0,   0.402484D0,   0.136127D0,  -0.200462D0, &
   -0.667435D0,  -1.410661D0,  -2.906862D0,  -8.169314D0 /)

  ! Range 5: near Ly δ, ν/ν_L = 0.938 to 0.959 by 0.001
  real*8, parameter :: GAVRILAMCD(22) = (/ &
   27.981406D0,   9.816495D0,   6.145775D0,   4.544224D0,   3.630968D0, &
    3.029081D0,   2.593248D0,   2.255265D0,   1.978565D0,   1.741426D0, &
    1.529699D0,   1.333240D0,   1.143898D0,   0.954154D0,   0.755875D0, &
    0.538760D0,   0.287687D0,  -0.022759D0,  -0.441666D0,  -1.081712D0, &
   -2.278530D0,  -5.705843D0 /)

  ! Range 7: above Lyman limit, correction factors
  real*8, parameter :: GAVRILALYMANCONT(64) = (/ &
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
  real*8, parameter :: FGAVRILALYMANCONT(64) = (/ &
    1.00D0, 1.05D0, 1.10D0, 1.15D0, 1.20D0, 1.25D0, 1.30D0, 1.35D0, &
    1.40D0, 1.45D0, 1.5D0,  1.6D0,  1.7D0,  1.8D0,  1.9D0,  2.0D0,  &
    2.1D0,  2.2D0,  2.3D0,  2.4D0,  2.5D0,  2.6D0,  2.7D0,  2.8D0,  &
    2.9D0,  3.0D0,  3.1D0,  3.2D0,  3.3D0,  3.4D0,  3.5D0,  3.6D0,  &
    3.7D0,  3.8D0,  3.9D0,  4.0D0,  4.4D0,  4.8D0,  5.2D0,  5.6D0,  &
    6.0D0,  6.4D0,  6.8D0,  7.2D0,  7.6D0,  8.0D0,  8.4D0,  8.8D0,  &
    9.2D0,  9.6D0, 10.0D0, 12.0D0, 14.0D0, 16.0D0, 18.0D0, 20.0D0, &
   24.0D0, 28.0D0, 32.0D0, 36.0D0, 40.0D0, 44.0D0, 48.0D0, 50.0D0 /)

  ! Local variables
  real*8  :: XSECT, G
  integer :: I, J, IDUM

  ! Frequency step sizes for each range
  real*8, parameter :: DFREQ1  = 0.01D0  * FREQ_RYDH  ! 0.01 × ν_L
  real*8, parameter :: DFREQAB = 0.005D0 * FREQ_RYDH  ! 0.005 × ν_L
  real*8, parameter :: DFREQBC = 0.002D0 * FREQ_RYDH  ! 0.002 × ν_L
  real*8, parameter :: DFREQCD = 0.001D0 * FREQ_RYDH  ! 0.001 × ν_L

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HRAYOP'

  XSECT = 0.0D0

  if (FREQ < DFREQ1) then
    ! Below 0.01 ν_L: extrapolate using lowest table value
    XSECT = SIGMA_THOMSON * GAVRILAM(1)**2 * (FREQ / DFREQ1)**4

  else if (FREQ <= 0.74D0 * FREQ_RYDH) then
    ! Range 1: 0.01–0.74 ν_L, table step = 0.01 ν_L
    I = int(FREQ / DFREQ1)
    I = MIN(I + 1, 74)
    G = GAVRILAM(I-1) + (GAVRILAM(I) - GAVRILAM(I-1)) / DFREQ1 &
      * (FREQ - dble(I-1) * DFREQ1)
    XSECT = SIGMA_THOMSON * G**2

  else if (FREQ < 0.77D0 * FREQ_RYDH) then
    ! Gap between range 1 and Ly β region
    XSECT = 0.0D0

  else if (FREQ <= 0.885D0 * FREQ_RYDH) then
    ! Range 3: 0.755–0.885 ν_L (near Ly β), step = 0.005 ν_L
    I = int((FREQ - 0.755D0 * FREQ_RYDH) / DFREQAB)
    I = I + 1
    I = MIN(I + 1, 27)
    G = GAVRILAMAB(I-1) + (GAVRILAMAB(I) - GAVRILAMAB(I-1)) / DFREQAB &
      * (FREQ - (0.755D0 * FREQ_RYDH + dble(I-1-1) * DFREQAB))
    XSECT = SIGMA_THOMSON * G**2

  else if (FREQ < 0.890D0 * FREQ_RYDH) then
    ! Gap between Ly β and Ly γ regions
    XSECT = 0.0D0

  else if (FREQ <= 0.936D0 * FREQ_RYDH) then
    ! Range 4: 0.890–0.936 ν_L (near Ly γ), step = 0.002 ν_L
    I = int((FREQ - 0.890D0 * FREQ_RYDH) / DFREQBC)
    I = I + 1
    I = MIN(I + 1, 24)
    G = GAVRILAMBC(I-1) + (GAVRILAMBC(I) - GAVRILAMBC(I-1)) / DFREQBC &
      * (FREQ - (0.890D0 * FREQ_RYDH + dble(I-1-1) * DFREQBC))
    XSECT = SIGMA_THOMSON * G**2

  else if (FREQ < 0.938D0 * FREQ_RYDH) then
    ! Gap between Ly γ and Ly δ regions
    XSECT = 0.0D0

  else if (FREQ <= 0.959D0 * FREQ_RYDH) then
    ! Range 5: 0.938–0.959 ν_L (near Ly δ), step = 0.001 ν_L
    I = int((FREQ - 0.938D0 * FREQ_RYDH) / DFREQCD)
    I = I + 1
    I = MIN(I + 1, 22)
    G = GAVRILAMCD(I-1) + (GAVRILAMCD(I) - GAVRILAMCD(I-1)) / DFREQCD &
      * (FREQ - (0.938D0 * FREQ_RYDH + dble(I-1-1) * DFREQCD))
    XSECT = SIGMA_THOMSON * G**2

  else if (FREQ < 0.961D0 * FREQ_RYDH) then
    ! Gap between Ly δ and series limit
    XSECT = 0.0D0

  else if (FREQ <= FREQ_RYDH) then
    ! Near series limit: constant value
    XSECT = SIGMA_THOMSON * GAVRILALYMANCONT(1)

  else
    ! Above Lyman limit: interpolate from continuum table
    BLOCK
      real*8 :: XNEW(1), FNEW(1)
      XNEW(1) = FREQ / FREQ_RYDH
      IDUM = MAP1(FGAVRILALYMANCONT, GAVRILALYMANCONT, 64, &
                  XNEW, FNEW, 1)
      XSECT = SIGMA_THOMSON * FNEW(1)
    END BLOCK

  end if

  ! Apply to all depth points
  do J = 1, NRHOX
    SIGH(J) = XSECT * XNFP(J, 1) * 2.0D0 * BHYD(J, 1) / RHO(J)
  end do

  return

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

  implicit none

  ! --- Energy level data for 10 resolved low-lying states ---
  real*8, parameter :: CHI(10) = (/ &
    0.0D0, 19.819D0, 20.615D0, 20.964D0, 21.217D0, &
    22.718D0, 22.920D0, 23.006D0, 23.073D0, 23.086D0 /)

  real*8, parameter :: HEFREQ(10) = (/ &
    5.945209D15, 1.152844D15, 0.9603331D15, 0.8761076D15, &
    0.8147104D15, 0.4519048D15, 0.4030971D15, 0.3821191D15, &
    0.3660215D15, 0.3627891D15 /)

  real*8, parameter :: G(10) = (/ &
    1.0D0, 3.0D0, 1.0D0, 9.0D0, 3.0D0, &
    3.0D0, 1.0D0, 9.0D0, 20.0D0, 3.0D0 /)

  ! --- Cached temperature-dependent arrays ---
  real*8,  save :: BOLT(kw, 10)     ! Boltzmann factors for 10 levels
  real*8,  save :: BOLTN(kw, 27)    ! Boltzmann factors for high-n levels
  real*8,  save :: EXLIM(kw)        ! population at ionization limit (24.587 eV)
  real*8,  save :: BOLTEX(kw)       ! population at dissolved limit (23.730 eV)
  real*8,  save :: FREET(kw)        ! free-free factor: n_e * n(He+) / (rho * sqrt(T))
  integer, save :: ITEMP1 = 0

  ! --- Local variables ---
  real*8  :: TRANS(10)     ! bound-free cross-sections for 10 levels
  real*8  :: TRANSN(27)    ! bound-free cross-sections for high-n levels
  real*8  :: FREQ3, CFREE, C
  real*8  :: RYD, ELIM, FREQHE, ZEFF2
  real*8  :: XR, EX, HE1
  integer :: J, N, IMIN

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE1OP'

  !---------------------------------------------------------------------
  ! Recompute temperature-dependent quantities if T has changed
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    RYD = RYDBERG_HE * CLIGHT
    do J = 1, NRHOX
      ! Boltzmann populations for 10 resolved levels
      do N = 1, 10
        BOLT(J, N) = exp(-CHI(N) / TKEV(J)) * G(N) * XNFP(J, 3) / RHO(J)
      end do
      ! High-n levels (n=4-27): hydrogenic with E_n = 24.587*(1-1/n^2) eV
      do N = 4, 27
        BOLTN(J, N) = exp(-24.587D0 * (1.0D0 - 1.0D0 / dble(N)**2) / TKEV(J)) &
                    * 4.0D0 * dble(N)**2 * XNFP(J, 3) / RHO(J)
      end do
      FREET(J) = XNE(J) * XNF(J, 4) / RHO(J) / sqrt(T(J))
      XR = XNFP(J, 3) * (4.0D0 / 2.0D0 / 13.595D0) * TKEV(J) / RHO(J)
      BOLTEX(J) = exp(-23.730D0 / TKEV(J)) * XR
      EXLIM(J) = exp(-24.587D0 / TKEV(J)) * XR
    end do
  end if

  !---------------------------------------------------------------------
  ! Cross-sections at the current frequency
  !---------------------------------------------------------------------
  FREQ3 = FREQ**3
  CFREE = COEFF_FF / FREQ3
  C = 2.815D29 / FREQ3
  RYD = RYDBERG_HE * CLIGHT

  ! Find lowest level whose threshold is at or below FREQ
  IMIN = 0
  do N = 1, 10
    if (HEFREQ(N) <= FREQ) then
      IMIN = N
      exit
    end if
  end do

  ! Compute bound-free cross-sections for levels IMIN..10
  ! (HEFREQ is in decreasing order: ground state has the highest threshold.
  ! IMIN is the first level whose threshold <= FREQ. The computed GOTO in
  ! the original fell through from label 20+IMIN to 30, so all levels
  ! N = IMIN..10 are computed.)
  TRANS = 0.0D0
  if (IMIN > 0 .and. IMIN <= 1)  TRANS(1)  = CROSSHE(FREQ)
  if (IMIN > 0 .and. IMIN <= 2)  TRANS(2)  = HE12s3S(FREQ)
  if (IMIN > 0 .and. IMIN <= 3)  TRANS(3)  = HE12s1S(FREQ)
  if (IMIN > 0 .and. IMIN <= 4)  TRANS(4)  = HE12p3P(FREQ)
  if (IMIN > 0 .and. IMIN <= 5)  TRANS(5)  = HE12p1P(FREQ)
  ! 1s3s 3S
  if (IMIN > 0 .and. IMIN <= 6)  TRANS(6)  = XKARZAS(FREQ, 1.236439D0, 3, 0)
  ! 1s3s 1S
  if (IMIN > 0 .and. IMIN <= 7)  TRANS(7)  = XKARZAS(FREQ, 1.102898D0, 3, 0)
  ! 1s3p 3P
  if (IMIN > 0 .and. IMIN <= 8)  TRANS(8)  = XKARZAS(FREQ, 1.045499D0, 3, 1)
  ! 1s3d 3D+1D
  if (IMIN > 0 .and. IMIN <= 9)  TRANS(9)  = XKARZAS(FREQ, 1.001427D0, 3, 2)
  ! 1s3p 1P
  if (IMIN > 0 .and. IMIN <= 10) TRANS(10) = XKARZAS(FREQ, 0.9926D0, 3, 1)

  !---------------------------------------------------------------------
  ! Inner-shell ionization: He I excited state → He II n=2
  ! (Adds hydrogenic 1s cross-section at shifted threshold)
  !---------------------------------------------------------------------
  if (IMIN >= 1) then
    ELIM = 527490.06D0
    ! 1s2p 1P → He II 2p
    FREQHE = (ELIM - 171135.000D0) * CLIGHT
    if (FREQ >= FREQHE) then
      ZEFF2 = FREQHE / RYD
      TRANS(5) = TRANS(5) + XKARZAS(FREQ, ZEFF2, 1, 0)
      ! 1s2p 3P → He II 2p
      FREQHE = (ELIM - 169087.0D0) * CLIGHT
      if (FREQ >= FREQHE) then
        ZEFF2 = FREQHE / RYD
        TRANS(4) = TRANS(4) + XKARZAS(FREQ, ZEFF2, 1, 0)
        ! 1s2s 1S → He II 2s
        FREQHE = (ELIM - 166277.546D0) * CLIGHT
        if (FREQ >= FREQHE) then
          ZEFF2 = FREQHE / RYD
          TRANS(3) = TRANS(3) + XKARZAS(FREQ, ZEFF2, 1, 0)
          ! 1s2s 3S → He II 2s
          FREQHE = (ELIM - 159856.069D0) * CLIGHT
          if (FREQ >= FREQHE) then
            ZEFF2 = FREQHE / RYD
            TRANS(2) = TRANS(2) + XKARZAS(FREQ, ZEFF2, 1, 0)
          end if
        end if
      end if
    end if

    !-------------------------------------------------------------------
    ! Inner-shell ionization: He I excited state → He II n=3
    !-------------------------------------------------------------------
    ELIM = 588451.59D0
    FREQHE = (ELIM - 186209.471D0) * CLIGHT
    if (FREQ >= FREQHE) then
      ZEFF2 = FREQHE / RYD
      TRANS(10) = TRANS(10) + XKARZAS(FREQ, ZEFF2, 1, 0)
      FREQHE = (ELIM - 186101.0D0) * CLIGHT
      if (FREQ >= FREQHE) then
        ZEFF2 = FREQHE / RYD
        TRANS(9) = TRANS(9) + XKARZAS(FREQ, ZEFF2, 1, 0)
        FREQHE = (ELIM - 185564.0D0) * CLIGHT
        if (FREQ >= FREQHE) then
          ZEFF2 = FREQHE / RYD
          TRANS(8) = TRANS(8) + XKARZAS(FREQ, ZEFF2, 1, 0)
          FREQHE = (ELIM - 184864.0D0) * CLIGHT
          if (FREQ >= FREQHE) then
            ZEFF2 = FREQHE / RYD
            TRANS(7) = TRANS(7) + XKARZAS(FREQ, ZEFF2, 1, 0)
            FREQHE = (ELIM - 183236.0D0) * CLIGHT
            if (FREQ >= FREQHE) then
              ZEFF2 = FREQHE / RYD
              TRANS(6) = TRANS(6) + XKARZAS(FREQ, ZEFF2, 1, 0)
            end if
          end if
        end if
      end if
    end if
  end if

  !---------------------------------------------------------------------
  ! High-n levels (n=4-27): hydrogenic with Z_eff^2 = 4 - 3/n^2
  !---------------------------------------------------------------------
  TRANSN = 0.0D0
  if (FREQ >= 1.25408D16) then
    do N = 4, 27
      ZEFF2 = 4.0D0 - 3.0D0 / dble(N)**2
      TRANSN(N) = XKARZAS(FREQ, ZEFF2, 1, 0)
    end do
  end if

  !---------------------------------------------------------------------
  ! Assemble total opacity at each depth
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    ! Dissolved-level contribution
    EX = BOLTEX(J)
    if (FREQ < 2.055D14) EX = EXLIM(J) / EHVKT(J)
    HE1 = (EX - EXLIM(J)) * C

    ! Bound-free from resolved levels
    if (IMIN > 0) then
      do N = IMIN, 10
        HE1 = HE1 + TRANS(N) * BOLT(J, N)
      end do
    end if

    ! High-n levels
    if (FREQ >= 1.25408D16) then
      do N = 4, 27
        HE1 = HE1 + TRANSN(N) * BOLTN(J, N)
      end do
    end if

    ! Total: bound + free-free (all LTE source function)
    !      AHE1BOUND(J) = HE1 * STIM(J)
    !      AHE1FREE(J) = (COULFF(J,1) * FREET(J) * CFREE) * STIM(J)
    AHE1(J) = (HE1 + COULFF(J, 1) * FREET(J) * CFREE) * STIM(J)
  end do

  return

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

  implicit none
  real*8, intent(in) :: FREQ
  real*8 :: CROSSHE

  real*8, parameter :: X505(92) = (/ &
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

  real*8, parameter :: X50(16) = (/ &
    0.0315D0, 0.0282D0, 0.0250D0, 0.0220D0, 0.0193D0, 0.0168D0, &
    0.0145D0, 0.0124D0, 0.0105D0, 0.00885D0, 0.00736D0, &
    0.00604D0, 0.00489D0, 0.00389D0, 0.00303D0, 0.00231D0 /)

  real*8, parameter :: X20(11) = (/ &
    0.00231D0, 0.00199D0, 0.00171D0, 0.00145D0, 0.00122D0, &
    0.00101D0, 0.000832D0, 0.000673D0, 0.000535D0, 0.000417D0, &
    0.000318D0 /)

  real*8, parameter :: X10(21) = (/ &
    0.000318D0, 0.000274D0, 0.000235D0, 0.000200D0, 0.000168D0, &
    0.000139D0, 0.000115D0, 0.000093D0, 0.000074D0, 0.000057D0, &
    0.000044D0, 0.000032D0, 0.000023D0, 0.000016D0, 0.000010D0, &
    0.000006D0, 0.000003D0, 0.000001D0, 0.0000006D0, &
    0.0000003D0, 0.0D0 /)

  real*8  :: WAVE
  integer :: I

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING CROSSHE'

  CROSSHE = 0.0D0
  if (FREQ < 5.945209D15) return

  WAVE = CLIGHT_ANG / FREQ

  if (WAVE > 50.0D0) then
    ! 50–505 Å regime (Δλ = 5 Å)
    I = int(93.0D0 - (WAVE - 50.0D0) / 5.0D0)
    I = min(92, max(2, I))
    CROSSHE = ((WAVE - (92 - I) * 5.0D0 - 50.0D0) / 5.0D0 &
             * (X505(I-1) - X505(I)) + X505(I)) * 1.D-18
  else if (WAVE > 20.0D0) then
    ! 20–50 Å regime (Δλ = 2 Å)
    I = int(17.0D0 - (WAVE - 20.0D0) / 2.0D0)
    I = min(16, max(2, I))
    CROSSHE = ((WAVE - (16 - I) * 2.0D0 - 20.0D0) / 2.0D0 &
             * (X50(I-1) - X50(I)) + X50(I)) * 1.D-18
  else if (WAVE > 10.0D0) then
    ! 10–20 Å regime (Δλ = 1 Å)
    I = int(12.0D0 - (WAVE - 10.0D0) / 1.0D0)
    I = min(11, max(2, I))
    CROSSHE = ((WAVE - (11 - I) * 1.0D0 - 10.0D0) / 1.0D0 &
             * (X20(I-1) - X20(I)) + X20(I)) * 1.D-18
  else
    ! 0–10 Å regime (Δλ = 0.5 Å)
    I = int(22.0D0 - WAVE / 0.5D0)
    I = min(21, max(2, I))
    CROSSHE = ((WAVE - (21 - I) * 0.5D0) / 0.5D0 &
             * (X10(I-1) - X10(I)) + X10(I)) * 1.D-18
  end if
  return

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

  implicit none
  real*8, intent(in) :: FREQ
  real*8 :: HE111S

  ! He I 1s² ¹S ground-state bound-free cross-section (after Mathisen)
  ! Linear interpolation in σ vs wavelength, 64-point table
  integer, parameter :: NP = 64
  real*8, parameter :: W(64) = (/ &
    504.3D0, 501.5D0, 498.7D0, 493.3D0, 488.1D0, 480.3D0, 477.8D0, &
    454.0D0, 443.0D0, 395.0D0, 356.4D0, 348.2D0, 324.6D0, 302.0D0, &
    298.1D0, 275.6D0, 260.6D0, 256.2D0, 239.4D0, 224.6D0, 220.0D0, &
    215.0D0, 210.0D0, 205.0D0, 200.0D0, 195.0D0, 190.0D0, 185.0D0, &
    180.0D0, 175.0D0, 170.0D0, 165.0D0, 160.0D0, 155.0D0, 150.0D0, &
    145.0D0, 135.0D0, 130.0D0, 125.0D0, 120.0D0, 115.0D0, 110.0D0, &
    105.0D0, 100.0D0, 95.0D0, 90.0D0, 85.0D0, 80.0D0, 75.0D0, &
    70.0D0, 65.0D0, 60.0D0, 55.0D0, 50.0D0, 45.0D0, 40.0D0, &
    35.0D0, 30.0D0, 25.0D0, 20.0D0, 15.0D0, 10.0D0, 5.0D0, 0.0D0 /)
  real*8, parameter :: X(64) = (/ &
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

  real*8  :: WAVE
  integer :: I

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE111S'

  HE111S = 0.0D0
  if (FREQ < 5.945209D15) return

  WAVE = CLIGHT_ANG / FREQ
  do I = 2, NP
    if (WAVE > W(I)) exit
  end do
  if (I > NP) I = NP

  HE111S = ((WAVE - W(I)) / (W(I-1) - W(I)) * (X(I-1) - X(I)) + X(I)) * 1.0D-18
  return

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

  implicit none
  real*8, intent(in) :: FREQ
  real*8 :: HE12S1S

  ! He I 1s2s ¹S bound-free cross-section
  ! Table interpolation below 2.4 Rydberg, resonance formula above
  integer, parameter :: NP = 16
  real*8, parameter :: FREQ1S(16) = (/ &
    15.947182D0, 15.913654D0, 15.877320D0, 15.837666D0, 15.794025D0, &
    15.745503D0, 15.690869D0, 15.628361D0, 15.555317D0, 15.467455D0, &
    15.357189D0, 15.289399D0, 15.251073D0, 15.209035D0, 15.162487D0, &
    14.982421D0 /)
  real*8, parameter :: X1S(16) = (/ &
    -19.635557D0, -19.159345D0, -18.958474D0, -18.809535D0, &
    -18.676481D0, -18.546006D0, -18.410962D0, -18.264821D0, &
    -18.100205D0, -17.909165D0, -17.684370D0, -17.557867D0, &
    -17.490360D0, -17.417876D0, -17.349386D0, -17.084441D0 /)

  real*8  :: FREQLG, X, EK, EPS, WAVNO
  integer :: I

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE12S1S'

  HE12S1S = 0.0D0
  if (FREQ < 32033.214D0 * CLIGHT) return

  if (FREQ > 2.4D0 * RYDBERG_HE * CLIGHT) then
    ! High-energy resonance formula
    WAVNO = FREQ / CLIGHT
    EK = (WAVNO - 32033.214D0) / RYDBERG_HE
    EPS = 2.0D0 * (EK - 2.612316D0) / 0.00322D0
    HE12S1S = 0.008175D0 * (484940.0D0 / WAVNO)**2.71D0 * 8.067D-18 &
            * (EPS + 76.21D0)**2 / (1.0D0 + EPS**2)
  else
    ! Table interpolation in log₁₀(σ) vs log₁₀(ν)
    FREQLG = log10(FREQ)
    do I = 2, NP
      if (FREQLG > FREQ1S(I)) exit
    end do
    if (I > NP) I = NP
    X = (FREQLG - FREQ1S(I)) / (FREQ1S(I-1) - FREQ1S(I)) &
      * (X1S(I-1) - X1S(I)) + X1S(I)
    HE12S1S = exp(X * 2.30258509299405D0)
  end if
  return

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

  implicit none
  real*8, intent(in) :: FREQ
  real*8 :: HE12S3S

  ! He I 1s2s ³S bound-free cross-section
  ! Table interpolation below 2.4 Rydberg, resonance formula above
  integer, parameter :: NP = 16
  real*8, parameter :: FREQ3S(16) = (/ &
    15.956523D0, 15.923736D0, 15.888271D0, 15.849649D0, 15.807255D0, &
    15.760271D0, 15.707580D0, 15.647601D0, 15.577992D0, 15.495055D0, &
    15.392451D0, 15.330345D0, 15.295609D0, 15.257851D0, 15.216496D0, &
    15.061770D0 /)
  real*8, parameter :: X3S(16) = (/ &
    -18.426022D0, -18.610700D0, -18.593051D0, -18.543304D0, &
    -18.465513D0, -18.378707D0, -18.278574D0, -18.164329D0, &
    -18.033346D0, -17.882435D0, -17.705542D0, -17.605584D0, &
    -17.553459D0, -17.500667D0, -17.451318D0, -17.266686D0 /)

  real*8  :: FREQLG, X, EK, EPS, WAVNO
  integer :: I

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE12S3S'

  HE12S3S = 0.0D0
  if (FREQ < 38454.691D0 * CLIGHT) return

  if (FREQ > 2.4D0 * RYDBERG_HE * CLIGHT) then
    ! High-energy resonance formula
    WAVNO = FREQ / CLIGHT
    EK = (WAVNO - 38454.691D0) / RYDBERG_HE
    EPS = 2.0D0 * (EK - 2.47898D0) / 0.000780D0
    HE12S3S = 0.01521D0 * (470310.0D0 / WAVNO)**3.12D0 * 8.067D-18 &
            * (EPS - 122.4D0)**2 / (1.0D0 + EPS**2)
  else
    ! Table interpolation in log₁₀(σ) vs log₁₀(ν)
    FREQLG = log10(FREQ)
    do I = 2, NP
      if (FREQLG > FREQ3S(I)) exit
    end do
    if (I > NP) I = NP
    X = (FREQLG - FREQ3S(I)) / (FREQ3S(I-1) - FREQ3S(I)) &
      * (X3S(I-1) - X3S(I)) + X3S(I)
    HE12S3S = exp(X * 2.30258509299405D0)
  end if
  return

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

  implicit none
  real*8, intent(in) :: FREQ
  real*8 :: HE12P1P

  ! He I 1s2p ¹P bound-free cross-section
  ! Table interpolation below 2.4 Rydberg, two-resonance formula above
  integer, parameter :: NP = 16
  real*8, parameter :: FREQ1P(16) = (/ &
    15.939981D0, 15.905870D0, 15.868850D0, 15.828377D0, 15.783742D0, &
    15.733988D0, 15.677787D0, 15.613218D0, 15.537343D0, 15.445346D0, &
    15.328474D0, 15.255641D0, 15.214064D0, 15.168081D0, 15.116647D0, &
    14.911002D0 /)
  real*8, parameter :: X1P(16) = (/ &
    -18.798876D0, -19.685922D0, -20.011664D0, -20.143030D0, &
    -20.091354D0, -19.908333D0, -19.656788D0, -19.367745D0, &
    -19.043016D0, -18.674484D0, -18.240861D0, -17.989700D0, &
    -17.852015D0, -17.702677D0, -17.525347D0, -16.816344D0 /)

  real*8  :: FREQLG, X, EK, EPS1S, EPS1D, WAVNO
  integer :: I

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE12P1P'

  HE12P1P = 0.0D0
  if (FREQ < 27175.76D0 * CLIGHT) return

  if (FREQ > 2.4D0 * RYDBERG_HE * CLIGHT) then
    ! High-energy: two autoionization resonances (¹S and ¹D channels)
    WAVNO = FREQ / CLIGHT
    EK = (WAVNO - 27175.76D0) / RYDBERG_HE
    EPS1S = 2.0D0 * (EK - 2.446534D0) / 0.01037D0
    EPS1D = 2.0D0 * (EK - 2.59427D0) / 0.00538D0
    HE12P1P = 0.0009487D0 * (466750.0D0 / WAVNO)**3.69D0 * 8.067D-18 &
            * ((EPS1S - 29.30D0)**2 / (1.0D0 + EPS1S**2) &
            + (EPS1D + 172.4D0)**2 / (1.0D0 + EPS1D**2))
  else
    ! Table interpolation in log₁₀(σ) vs log₁₀(ν)
    FREQLG = log10(FREQ)
    do I = 2, NP
      if (FREQLG > FREQ1P(I)) exit
    end do
    if (I > NP) I = NP
    X = (FREQLG - FREQ1P(I)) / (FREQ1P(I-1) - FREQ1P(I)) &
      * (X1P(I-1) - X1P(I)) + X1P(I)
    HE12P1P = exp(X * 2.30258509299405D0)
  end if
  return

END FUNCTION HE12P1P

!=======================================================================
! HE12P3P: He I 1s2p ³P bound-free photoionization cross-section
!
! Tabulated cross-section for photoionization of He I from the
! 1s2p ³P excited state. 16-point table of log₁₀(σ) vs log₁₀(ν),
! linearly interpolated. Threshold at 29223.753 cm⁻¹.
!=======================================================================

FUNCTION HE12P3P(FREQ)

  implicit none
  real*8, intent(in) :: FREQ
  real*8 :: HE12P3P

  ! He I 1s2p ³P bound-free cross-section
  ! Linear interpolation in log₁₀(σ) vs log₁₀(ν), 16-point table
  integer, parameter :: NP = 16
  real*8, parameter :: FREQ3P(16) = (/ &
    15.943031D0, 15.909169D0, 15.872441D0, 15.832318D0, 15.788107D0, &
    15.738880D0, 15.683351D0, 15.619667D0, 15.545012D0, 15.454805D0, &
    15.340813D0, 15.270195D0, 15.230054D0, 15.185821D0, 15.136567D0, &
    14.942557D0 /)
  real*8, parameter :: X3P(16) = (/ &
    -19.791021D0, -19.697886D0, -19.591421D0, -19.471855D0, &
    -19.337053D0, -19.183958D0, -19.009750D0, -18.807990D0, &
    -18.570571D0, -18.288361D0, -17.943476D0, -17.738737D0, &
    -17.624154D0, -17.497163D0, -17.403183D0, -17.032999D0 /)

  real*8  :: FREQLG, X
  integer :: I

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE12P3P'

  HE12P3P = 0.0D0
  if (FREQ < 29223.753D0 * CLIGHT) return

  FREQLG = log10(FREQ)
  do I = 2, NP
    if (FREQLG > FREQ3P(I)) exit
  end do
  if (I > NP) I = NP

  X = (FREQLG - FREQ3P(I)) / (FREQ3P(I-1) - FREQ3P(I)) &
    * (X3P(I-1) - X3P(I)) + X3P(I)
  HE12P3P = exp(X * 2.30258509299405D0)
  return

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

  implicit none

  real*8  :: FREQ3, XNFPRHO
  real*8  :: H, S, A, X
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HE2OP'

  FREQ3 = 2.815D29 / FREQ / FREQ / FREQ

  do J = 1, NRHOX

    ! Population factor: He II (mode-11) / rho
    XNFPRHO = XNFP(J, 4) / RHO(J)

    ! --- n >= 10 to infinity: dissolved-level integral (LTE) ---
    ! Z⁴=16, gfactor=2, ionization limit = 438908.85 cm⁻¹
    ! Lower integration bound at n=10: E_10 = 438908.85 - 438908.85/100 = 434519.959
    H = FREQ3 * 16.0D0 * 2.0D0 / 2.0D0 / (438889.068D0 * HCKT(J)) &
      * (EXP(-MAX(434519.959D0, 438908.85D0 - WAVENO) * HCKT(J)) &
       - EXP(-438908.85D0 * HCKT(J))) * STIM(J) * XNFPRHO
    S = H * BNU(J)

    levels: do

      ! --- n = 9 (threshold 5418.390 cm⁻¹) : hydrogenic ---
      if (WAVENO < 5418.390D0) exit levels
      X = FREQ3 / 59049.0D0 * 16.0D0
      A = X * 162.0D0 * EXP(-433490.46D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 8 (threshold 6857.660 cm⁻¹) : hydrogenic ---
      if (WAVENO < 6857.660D0) exit levels
      X = FREQ3 * 16.0D0 / 32768.0D0
      A = X * 128.0D0 * EXP(-432051.19D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 7 (threshold 8956.950 cm⁻¹) : hydrogenic ---
      if (WAVENO < 8956.950D0) exit levels
      X = FREQ3 * 16.0D0 / 16807.0D0
      A = X * 98.0D0 * EXP(-429951.90D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 6 (threshold 12191.437 cm⁻¹) : Gaunt-corrected, LTE ---
      if (WAVENO < 12191.437D0) exit levels
      X = FREQ3 * 16.0D0 / 7776.0D0 &
        * (1.0986D0 + (-2.704D13 + 1.229D27 / FREQ) / FREQ)
      A = X * 72.0D0 * EXP(-426717.413D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 5 (threshold 17555.715 cm⁻¹) : Gaunt-corrected, LTE ---
      if (WAVENO < 17555.715D0) exit levels
      X = FREQ3 * 16.0D0 / 3125.0D0 &
        * (1.102D0 + (-3.909D13 + 2.371D27 / FREQ) / FREQ)
      A = X * 50.0D0 * EXP(-421353.135D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 4 (threshold 27430.925 cm⁻¹) : Gaunt-corrected, LTE ---
      if (WAVENO < 27430.925D0) exit levels
      X = FREQ3 * 16.0D0 / 1024.0D0 &
        * (1.101D0 + (-5.765D13 + 4.593D27 / FREQ) / FREQ)
      A = X * 32.0D0 * EXP(-411477.925D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 3 (threshold 48766.491 cm⁻¹) : Gaunt-corrected, LTE ---
      ! Note: denominator is 243 (= 3⁵), not 343 (typo corrected in atlas7lib)
      if (WAVENO < 48766.491D0) exit levels
      X = FREQ3 * 16.0D0 / 243.0D0 &
        * (1.101D0 + (-9.863D13 + 1.035D28 / FREQ) / FREQ)
      A = X * 18.0D0 * EXP(-390142.359D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 2 (threshold 109726.529 cm⁻¹) : Gaunt-corrected, LTE ---
      if (WAVENO < 109726.529D0) exit levels
      X = FREQ3 * 16.0D0 / 32.0D0 &
        * (1.105D0 + (-2.375D14 + 4.077D28 / FREQ) / FREQ)
      A = X * 8.0D0 * EXP(-329182.321D0 * HCKT(J)) * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      ! --- n = 1 (threshold 438908.850 cm⁻¹) : Gaunt-corrected, LTE ---
      if (WAVENO < 438908.850D0) exit levels
      X = FREQ3 * 16.0D0 / 1.0D0 &
        * (0.9916D0 + (2.719D13 - 2.268D30 / FREQ) / FREQ)
      A = X * 2.0D0 * STIM(J) * XNFPRHO
      H = H + A;  S = S + A * BNU(J)

      exit levels
    end do levels

    ! --- Free-free (bremsstrahlung, Z²=4) ---
    A = COEFF_FF * 4.0D0 / SQRT(T(J)) * COULFF(J, 2) / FREQ * XNE(J) &
      / FREQ * XNFP(J, 5) / FREQ * STIM(J) / RHO(J)
    H = H + A
    S = S + A * BNU(J)

    AHE2(J) = H
    SHE2(J) = BNU(J)
    if (H > 0.0D0) SHE2(J) = S / H

  end do

  return

END SUBROUTINE HE2OP

!=======================================================================
! HEMIOP
!=======================================================================

SUBROUTINE HEMIOP

  implicit none
  real*8  :: A, B, C
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HEMIOP'
  ! He⁻ free-free opacity: polynomial fit in 1/freq and T
  A = 3.397D-46 + (-5.216D-31 + 7.039D-15 / FREQ) / FREQ
  B = -4.116D-42 + (1.067D-26 + 8.135D-11 / FREQ) / FREQ
  C = 5.081D-37 + (-8.724D-23 - 5.659D-8 / FREQ) / FREQ
  do J = 1, NRHOX
    AHEMIN(J) = (A * T(J) + B + C / T(J)) * XNE(J) * XNFP(J, 3) / RHO(J)
  end do
  return

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
  implicit none
  real*8  :: W, WW, SIG
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HERAOP'

  ! Wavelength in Å, capped just redward of the He I 584 Å resonance
  ! to prevent the dispersion denominator from going singular.
  W  = CLIGHT_ANG / min(FREQ, 5.15D15)
  WW = W**2                                   ! λ² [Å²]

  ! σ(λ) = (5.484e-14 / λ⁴) * [α(ω)/α(0)]²
  !   Leading coefficient 5.484e-14 = (128π⁵/3) * α_He(static)²
  !   Bracketed factor is the fitted two-pole dispersion correction.
  SIG = 5.484D-14 / WW / WW * (1.0D0 + (2.44D5 + 5.94D10 &
      / (WW - 2.90D5)) / WW)**2

  ! Apply to each depth point using neutral He number density
  ! (XNFP(J,3) is N(He I)); /RHO converts to mass opacity [cm²/g].
  do J = 1, NRHOX
    SIGHE(J) = SIG * XNFP(J, 3) / RHO(J)
  end do
  return

END SUBROUTINE HERAOP

!=======================================================================
! COOLOP: Cool-star opacity assembly (T < 7000 K)
!=======================================================================

SUBROUTINE COOLOP

  implicit none
  integer :: J
  real*8  :: AH2COLL(kw)

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING COOLOP'
  ! Cool-star opacities: only contribute below Lyman limit
  if (FREQ > FREQ_RYDH) return
  call C1OP
  call MG1OP
  call AL1OP
  call SI1OP
  call FE1OP
  call H2COLLOP(AH2COLL)
  do J = 1, NRHOX
    ACOOL(J) = AC1(J) + AMG1(J) + AAL1(J) + ASI1(J) + AFE1(J) &
             + (CHOP(J) * XNFP(J, 846) + OHOP(J) * XNFP(J, 848)) &
             * STIM(J) / RHO(J) + AH2COLL(J)
  end do
  return

END SUBROUTINE COOLOP

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

  implicit none

  real*8, parameter :: RYD = 109732.298D0

  real*8  :: H, S, A, B, X, EPS
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING C1OP'

  do J = 1, NRHOX
    H = 1.0D-30
    S = 0.0D0

    levels: do

      ! --- Level 13: 3p 1S (E=73975.91, g=1, threshold 16886.790) X=0 ---
      if (WAVENO < 16886.790D0) exit levels

      ! --- Level 12: 3p 1D (E=72610.72, g=5, threshold 18251.980) X=0 ---
      if (WAVENO < 18251.980D0) exit levels

      ! --- Level 11: 3p 3P (E=71374.90, g=9, threshold 19487.800) X=0 ---
      if (WAVENO < 19487.800D0) exit levels

      ! --- Level 10: 3p 3S (E=70743.95, g=3, threshold 20118.750) X=0 ---
      if (WAVENO < 20118.750D0) exit levels

      ! --- Level 9: 3p 3D (E=69722.00, g=15, threshold 21140.700) X=0 ---
      if (WAVENO < 21140.700D0) exit levels

      ! --- Level 8: 3p 1P (E=68856.33, g=3, threshold 22006.370) ---
      if (WAVENO < 22006.370D0) exit levels
      X = 2.1D-18 * (22006.370D0 / WAVENO)**1.5D0
      A = X * 3.0D0 * EXP(-68856.33D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 6: 3s 1P (E=61981.82, g=3, threshold 28880.880) ---
      if (WAVENO < 28880.880D0) exit levels
      X = 1.54D-18 * (28880.880D0 / WAVENO)**1.2D0
      A = X * 3.0D0 * EXP(-61981.82D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 5: 3s 3P (E=60373.00, g=9, threshold 30489.700) ---
      if (WAVENO < 30489.700D0) exit levels
      X = 0.2D-18 * (30489.700D0 / WAVENO)**1.2D0
      A = X * 9.0D0 * EXP(-60373.00D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 14: 2s2p3 3P (E=75254.93, g=9, threshold 58601.270) X=0 ---
      !     Ionizes to 4P limit at 133856.20 cm-1
      if (WAVENO < 58601.270D0) exit levels

      ! --- Level 3: 2p2 1S (E=21648.02, g=1) ---
      !     Ionizes to 2P_0.5 at 90820.42 -> threshold 69172.400
      !     Luo & Pradhan background + Burke & Taylor Fano resonance
      if (WAVENO < 69172.400D0) exit levels
      X = 10.0D0**(-16.80D0 - (WAVENO - 69172.400D0) / 3.0D0 / RYD)
      EPS = (WAVENO - 97700.0D0) * 2.0D0 / 2743.0D0
      A = 68.0D-18;  B = 118.0D-18
      X = X + (A * EPS + B) / (EPS**2 + 1.0D0)
      X = X / 3.0D0
      A = X * 1.0D0 * EXP(-21648.02D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     Also ionizes to 2P_1.5 at 90883.84 -> threshold 69235.820
      if (WAVENO < 69235.820D0) exit levels
      A = A * 2.0D0
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 7: 2s2p3 3D (E=64088.85, g=15, threshold 69767.350) ---
      !     Ionizes to 4P limit at 133856.20
      if (WAVENO < 69767.350D0) exit levels
      X = 16.0D-18 * (69767.350D0 / WAVENO)**3
      A = X * 15.0D0 * EXP(-64088.85D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 2: 2p2 1D (E=10192.66, g=5) ---
      !     Ionizes to 2P_0.5 at 90820.42 -> threshold 80627.760
      !     Luo & Pradhan background + two Burke & Taylor Fano resonances
      if (WAVENO < 80627.760D0) exit levels
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
      if (WAVENO < 80691.180D0) exit levels
      A = A * 2.0D0
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 1: 2p2 3P (ground term, 3 fine-structure components) ---
      !     Ionizes to 2P_0.5 at 90820.42
      !     3P_2 (E=43.42, g=5) -> threshold 90777.000
      if (WAVENO < 90777.000D0) exit levels
      X = 10.0D0**(-16.80D0 - (WAVENO - 90777.000D0) / 3.0D0 / RYD)
      X = X / 3.0D0
      A = X * 5.0D0 * EXP(-43.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_1 (E=16.42, g=3) -> threshold 90804.000
      if (WAVENO < 90804.000D0) exit levels
      A = X * 3.0D0 * EXP(-16.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_0 (E=0.00, g=1) -> threshold 90820.420
      if (WAVENO < 90820.420D0) exit levels
      A = X * 1.0D0 * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     Ionizes to 2P_1.5 at 90883.84 (cross-section x 2)
      !     3P_2 -> threshold 90840.420
      if (WAVENO < 90840.420D0) exit levels
      X = X * 2.0D0
      A = X * 5.0D0 * EXP(-43.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_1 -> threshold 90867.420
      if (WAVENO < 90867.420D0) exit levels
      A = X * 3.0D0 * EXP(-16.42D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      !     3P_0 -> threshold 90883.840
      if (WAVENO < 90883.840D0) exit levels
      A = X * 1.0D0 * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! --- Level 4: 2s2p3 5S (E=33735.20, g=5, threshold 100121.000) ---
      !     Ionizes to 4P limit at 133856.20
      if (WAVENO < 100121.000D0) exit levels
      X = 1.0D-18 * (100121.000D0 / WAVENO)**3
      A = X * 5.0D0 * EXP(-33735.20D0 * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      exit levels
    end do levels

    ! Scale by C I population / rho
    H = H * XNFP(J, 21) / RHO(J)
    S = S * XNFP(J, 21) / RHO(J)

    AC1(J) = H
    if (H > 0.0D0) SC1(J) = S / H

  end do

  return

END SUBROUTINE C1OP

!=======================================================================
! SEATON
!=======================================================================

FUNCTION SEATON(FREQ0, XSECT, POWER, A)

  implicit none
  real*8, intent(in) :: FREQ0, XSECT, POWER, A
  real*8 :: SEATON

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING SEATON'
  ! Seaton (1958) photoionization cross-section formula:
  ! sigma(nu) = sigma_0 * [A + (1-A)*(nu_0/nu)] * (nu_0/nu)^POWER
  SEATON = XSECT * (A + (1.0D0 - A) * (FREQ0 / FREQ)) &
         * sqrt((FREQ0 / FREQ)**int(2.0D0 * POWER + 0.01D0))
  return

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

  implicit none

  integer, parameter :: NLEV = 15

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
  real*8, parameter :: ELEV(NLEV) = (/ &
    54676.710D0, 54676.438D0, 54192.284D0, 53134.642D0, 49346.729D0, &
    47957.034D0, 47847.797D0, 46403.065D0, 43503.333D0, 41197.043D0, &  
                                                                        
    35051.264D0, 21919.178D0, 21870.464D0, 21850.405D0,     0.000D0 /)

  real*8, parameter :: GLEV(NLEV) = (/ &
    21.D0, 7.D0, 15.D0, 5.D0, 3.D0, 15.D0, 9.D0, 5.D0, &
     1.D0, 3.D0,  3.D0, 5.D0, 3.D0,  1.D0, 1.D0 /)

  ! Quantum numbers (n,l) for XKARZAS — only levels 1–5 are hydrogenic
  integer, parameter :: NQ(5) = (/ 4, 4, 4, 4, 4 /)
  integer, parameter :: LQ(5) = (/ 3, 3, 2, 2, 1 /)

  ! Ionization limit (cm⁻¹): Mg II 3s ²S
  real*8, parameter :: ELIM = 61671.02D0
  real*8, parameter :: RYD  = 109732.298D0

  ! Threshold wavenumbers for power-law fits (ELIM - ELEV)
  real*8, parameter :: THR6  = 13713.986D0   ! 3s3d ³D
  real*8, parameter :: THR7  = 13823.223D0   ! 3s4p ³P
  real*8, parameter :: THR8  = 15267.955D0   ! 3s3d ¹D
  real*8, parameter :: THR9  = 18167.687D0   ! 3s4s ¹S
  real*8, parameter :: THR10 = 20473.617D0   ! 3s4s ³S
  real*8, parameter :: THR11 = 26619.756D0   ! 3s3p ¹P
  real*8, parameter :: THR12 = 39759.842D0   ! 3s3p ³P

  ! --- Persistent state ---
  real*8,  save :: BOLT(NLEV, kw)
  integer, save :: ITEMP1 = 0

  ! --- Local variables ---
  real*8  :: X(NLEV)
  real*8  :: ZEFF2, FREQ3, H, RATIO
  integer :: I, J, K

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING MG1OP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors when temperature structure changes
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      do I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      end do
    end do
  end if

  !---------------------------------------------------------------------
  ! Compute cross-sections at the current frequency
  !---------------------------------------------------------------------
  X(:) = 0.0D0

  levels: do

    ! Levels 1–5: hydrogenic (XKARZAS)
    do I = 1, 5
      if (WAVENO < ELIM - ELEV(I)) exit levels
      ZEFF2 = 16.0D0 / RYD * (ELIM - ELEV(I))
      X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
    end do

    ! Level 6: 3s3d ³D — power-law fit
    if (WAVENO < ELIM - ELEV(6)) exit levels
    X(6) = 25.0D-18 * (THR6 / WAVENO)**2.7D0

    ! Level 7: 3s4p ³P — power-law fit
    if (WAVENO < ELIM - ELEV(7)) exit levels
    X(7) = 33.8D-18 * (THR7 / WAVENO)**2.8D0

    ! Level 8: 3s3d ¹D — power-law fit
    if (WAVENO < ELIM - ELEV(8)) exit levels
    X(8) = 45.0D-18 * (THR8 / WAVENO)**2.7D0

    ! Level 9: 3s4s ¹S — power-law fit
    if (WAVENO < ELIM - ELEV(9)) exit levels
    X(9) = 0.43D-18 * (THR9 / WAVENO)**2.6D0

    ! Level 10: 3s4s ³S — power-law fit
    if (WAVENO < ELIM - ELEV(10)) exit levels
    X(10) = 2.1D-18 * (THR10 / WAVENO)**2.6D0

    ! Level 11: 3s3p ¹P — two-term power-law
    if (WAVENO < ELIM - ELEV(11)) exit levels
    RATIO = THR11 / WAVENO
    X(11) = 16.0D-18 * RATIO**2.1D0 - 7.8D-18 * RATIO**9.5D0

    ! Levels 12–14: 3s3p ³P₂,₁,₀ — power-law with floor
    do I = 12, 14
      if (WAVENO < ELIM - ELEV(I)) exit levels
      RATIO = THR12 / WAVENO
      X(I) = max(20.0D-18 * RATIO**2.7D0, 40.0D-18 * RATIO**14.0D0)
    end do

    ! Level 15: 3s² ¹S ground state — steep power-law
    ! (Castelli index correction 25 Sep 2002: test on level 15, not 13)
    if (WAVENO < ELIM - ELEV(15)) exit levels
    X(15) = 1.1D-18 * ((ELIM - ELEV(15)) / WAVENO)**10.0D0

    exit levels
  end do levels

  !---------------------------------------------------------------------
  ! Assemble opacity over depth: levels + dissolved high-n contribution
  !---------------------------------------------------------------------
  FREQ3 = 2.815D29 / (FREQ * FREQ * FREQ)

  do J = 1, NRHOX
    ! High-n dissolved levels (n >= 5 to infinity), GFACTOR = 2
    H = FREQ3 * 2.0D0 * 2.0D0 / 2.0D0 / (RYD * HCKT(J)) &
      * (exp(-max(ELIM - RYD / 25.0D0, ELIM - WAVENO) * HCKT(J)) &
       - exp(-ELIM * HCKT(J)))

    do I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    end do

    AMG1(J) = H * XNFP(J, 78) * STIM(J) / RHO(J)
  end do

  return

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

  implicit none

  integer, parameter :: NLEV = 10

  ! Threshold wavenumbers (cm⁻¹) for each level
  real*8, parameter :: THR(NLEV) = (/ &
     6958.993D0,  8002.467D0,  9346.231D0, 10588.957D0, 15318.007D0, &
    15842.129D0, 22930.614D0, 48166.309D0, 48278.370D0, 55903.260D0 /)

  ! Level energies (cm⁻¹)
  real*8, parameter :: ELEV(NLEV) = (/ &
    41319.377D0, 40275.903D0, 38932.139D0, 37689.413D0, 32960.363D0, &
    32436.241D0, 25347.756D0,   112.061D0,     0.000D0, 29097.110D0 /)

  ! Statistical weights
  real*8, parameter :: GLEV(NLEV) = (/ &
    14.0D0,  6.0D0, 10.0D0,  2.0D0,  6.0D0, &
    10.0D0,  2.0D0,  4.0D0,  2.0D0, 12.0D0 /)

  real*8  :: H, S, A, X
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING AL1OP'

  do J = 1, NRHOX
    H = 1.0D-30
    S = 0.0D0

    levels: do
      ! Level 1: 4f ²F (X=0, no contribution — retained for fidelity)
      if (WAVENO < THR(1)) exit levels
      ! X = 0, so A = 0, skip

      ! Level 2: 5p ²P
      if (WAVENO < THR(2)) exit levels
      X = 50.0D-18 * (THR(2) / WAVENO)**3
      A = X * GLEV(2) * EXP(-ELEV(2) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 3: 4d ²D
      if (WAVENO < THR(3)) exit levels
      X = 50.0D-18 * (THR(3) / WAVENO)**3
      A = X * GLEV(3) * EXP(-ELEV(3) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 4: 5s ²S
      if (WAVENO < THR(4)) exit levels
      X = 56.7D-18 * (THR(4) / WAVENO)**1.9D0
      A = X * GLEV(4) * EXP(-ELEV(4) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 5: 4p ²P
      if (WAVENO < THR(5)) exit levels
      X = 14.5D-18 * THR(5) / WAVENO
      A = X * GLEV(5) * EXP(-ELEV(5) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 6: 3d ²D
      if (WAVENO < THR(6)) exit levels
      X = 47.0D-18 * (THR(6) / WAVENO)**1.83D0
      A = X * GLEV(6) * EXP(-ELEV(6) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 7: 4s ²S
      if (WAVENO < THR(7)) exit levels
      X = 10.0D-18 * (THR(7) / WAVENO)**2
      A = X * GLEV(7) * EXP(-ELEV(7) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 8: 3p ²P₃/₂ (ground fine-structure)
      if (WAVENO < THR(8)) exit levels
      X = 65.0D-18 * (THR(8) / WAVENO)**5
      A = X * GLEV(8) * EXP(-ELEV(8) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 9: 3p ²P₁/₂ (ground state) — same cross-section formula
      if (WAVENO < THR(9)) exit levels
      A = X * GLEV(9) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      ! Level 10: P2 ⁴P (inner-shell)
      if (WAVENO < THR(10)) exit levels
      X = 10.0D-18 * (THR(10) / WAVENO)**2
      A = X * GLEV(10) * EXP(-ELEV(10) * HCKT(J)) * STIM(J)
      H = H + A;  S = S + A * BNU(J)

      exit levels
    end do levels

    if (H > 0.0D0) SAL1(J) = S / H
    AAL1(J) = H * XNFP(J, 91) / RHO(J)
  end do

  return

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

  implicit none

  integer, parameter :: NLEV = 33

  ! Energy levels (cm⁻¹) and statistical weights for 33 Si I states
  real*8, parameter :: ELEV(33) = (/ &
    59962.284D0, 59100.000D0, 59077.112D0, 58893.40D0, 58801.529D0, &  ! 4d³P, 4f, 4d³D, 4d¹F, 4d¹P
    58777.000D0, 57488.974D0, 56503.346D0, 54225.621D0, 53387.34D0, &  ! 4f, 4d³F, 4d¹D, 3d³D, 3d¹P
    53362.24D0,  51612.012D0, 50533.424D0, 50189.389D0, 49965.894D0, & ! 3d¹F, 4p¹S, 3d³P, 4p¹D, 3d³F
    49399.670D0, 49128.131D0, 48161.459D0, 47351.554D0, 47284.061D0, & ! 4p³S, 4p³P, 4p³D, 3d¹D, 4p¹P
    40991.884D0, 39859.920D0, 15394.370D0,  6298.850D0,   223.157D0, & ! 4s¹P, 4s³P, 3p²¹S, 3p²¹D, 3p²³P₂
    77.115D0,        0.000D0, 94000.000D0, 79664.000D0, 72000.000D0, & ! ³P₁, ³P₀, 3s3p³¹P, ³S, ¹D
    56698.738D0, 45303.310D0, 33326.053D0 /)                           ! ³P, ³D, ⁵S

  real*8, parameter :: GLEV(33) = (/ &
    9.D0, 56.D0, 15.D0, 7.D0, 3.D0, 28.D0, 21.D0, 5.D0, 15.D0, &
    3.D0, 7.D0, 1.D0, 9.D0, 5.D0, 21.D0, 3.D0, 9.D0, 15.D0, &
    5.D0, 3.D0, 3.D0, 9.D0, 1.D0, 5.D0, 5.D0, 3.D0, 1.D0, &
    3.D0, 3.D0, 5.D0, 12.D0, 15.D0, 5.D0 /)

  ! L quantum numbers for XKARZAS calls (levels 1-22)
  ! d-orbital → l=2, f-orbital → l=3, p-orbital → l=1, s-orbital → l=0
  integer, parameter :: LQNUM(22) = (/ &
    2, 3, 2, 2, 2, 3, 2, 2,  &  ! levels 1-8 (n=4 shell)
    2, 2, 2, 1, 2, 1, 2, 1,  &  ! levels 9-16 (n=3d and n=4p)
    1, 1, 2, 1, 0, 0 /)          ! levels 17-22 (n=4p, 3d, 4s)

  ! Principal quantum numbers for XKARZAS calls (levels 1-22)
  integer, parameter :: NQNUM(22) = (/ &
    4, 4, 4, 4, 4, 4, 4, 4,  &  ! levels 1-8
    3, 3, 3, 4, 3, 4, 3, 4,  &  ! levels 9-16
    4, 4, 3, 4, 4, 4 /)          ! levels 17-22

  ! Z_eff² prefactor: n²/RYD for each level (levels 1-22)
  ! n=3 levels use 9/RYD, n=4 levels use 16/RYD
  real*8, parameter :: RYD = 109732.298D0

  real*8, save :: BOLT(NLEV, kw)
  integer, save :: ITEMP1 = 0

  real*8  :: X(NLEV), FREQ3, ZEFF2, EPS, RESON1
  real*8  :: H, GFACTOR, DEGEN, Z, ELIM_BLK
  integer :: I, J, K

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING SI1OP'

  ! Recompute Boltzmann factors when temperature changes
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      do I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      end do
    end do
  end if

  Z = 1.0D0
  FREQ3 = 2.815D29 / FREQ / FREQ / FREQ * Z**4
  WAVENO = FREQ / CLIGHT

  do I = 1, NLEV
    X(I) = 0.0D0
  end do

  ! =====================================================================
  ! Block 1: Excited states → Si II 3s² 3p ²P average (levels 1-22)
  !   Hydrogenic cross-sections via XKARZAS
  ! =====================================================================
  ELIM_BLK = 65939.18D0

  do I = 1, 22
    if (WAVENO < ELIM_BLK - ELEV(I)) exit
    ZEFF2 = dble(NQNUM(I))**2 / RYD * (ELIM_BLK - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQNUM(I), LQNUM(I))
  end do

  ! =====================================================================
  ! Block 2: Low-lying 3p² states → Si II 3s² 3p ²P₁/₂ (levels 23-27)
  !   Nahar & Pradhan (1993) fits with Fano resonance profiles
  !   Statistical weight factor: 1/3 for ²P₁/₂ (g=2 of g_total=6)
  ! =====================================================================
  ELIM_BLK = 65747.55D0

  !  3s² 3p² ¹S (level 23)
  if (WAVENO >= ELIM_BLK - ELEV(23)) then
    EPS = (WAVENO - 70000.0D0) * 2.0D0 / 6500.0D0
    ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
    RESON1 = (97.D-18 * EPS + 94.D-18) / (EPS**2 + 1.0D0)
!     EPS=(WAVENO-89700.)*2./75.
!     RESON2=900.E-18/(EPS**2+1.)
    X(23) = (37.D-18 * (50353.180D0 / WAVENO)**2.40D0 + RESON1) / 3.0D0
!     X(23)=46.E-18*(50353.180/WAVENO)**.5/3.

    !  3s² 3p² ¹D (level 24)
    if (WAVENO >= ELIM_BLK - ELEV(24)) then
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
      if (WAVENO >= ELIM_BLK - ELEV(25)) then
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
        if (WAVENO <= 74000.0D0) then
          X(25) = 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 / 3.0D0
        else
          X(25) = 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 / 3.0D0
        end if
!     X(25)=X(25)+(RESON1+RESON2+RESON3)/3.
!     X(25)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)/3.

        !  3s² 3p² ³P₁ (level 26)
        if (WAVENO >= ELIM_BLK - ELEV(26)) then
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
          if (WAVENO <= 74000.0D0) then
            X(26) = 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
          else
            X(26) = 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
          end if
!     X(26)=X(26)+(RESON1+RESON2+RESON3)*2./3.
!     X(26)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.

          !  3s² 3p² ³P₀ (level 27)
          if (WAVENO >= ELIM_BLK - ELEV(27)) then
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
            if (WAVENO <= 74000.0D0) then
              X(27) = 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 / 3.0D0
            else
              X(27) = 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 / 3.0D0
            end if
!     X(27)=X(27)+(RESON1+RESON2+RESON3)/3.
!     X(27)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)/3.
          end if  ! level 27
        end if  ! level 26
      end if  ! level 25
    end if  ! level 24
  end if  ! level 23

  ! =====================================================================
  ! Block 3: Same low-lying 3p² states → Si II 3s² 3p ²P₃/₂ (levels 23-27)
  !   Adds ²P₃/₂ contribution (statistical weight 2/3) to X(23)-X(27)
  ! =====================================================================
  ELIM_BLK = 65747.55D0 + 287.45D0

  !  3s² 3p² ¹S (level 23)
  if (WAVENO >= ELIM_BLK - ELEV(23)) then
    EPS = (WAVENO - 70000.0D0) * 2.0D0 / 6500.0D0
    ! fits to Nahar & Pradhan, J.Phys.B 26, 1109-1127, 1993.
    RESON1 = (97.D-18 * EPS + 94.D-18) / (EPS**2 + 1.0D0)
!     EPS=(WAVENO-89700.)*2./75.
!     RESON2=900.E-18/(EPS**2+1.)
    X(23) = X(23) + (37.D-18 * (50353.180D0 / WAVENO)**2.40D0 + RESON1) * 2.0D0 / 3.0D0
!     X(23)=X(23)+46.E-18*(50353.180/WAVENO)**.5*2./3.

    !  3s² 3p² ¹D (level 24)
    if (WAVENO >= ELIM_BLK - ELEV(24)) then
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
      if (WAVENO >= ELIM_BLK - ELEV(25)) then
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
        if (WAVENO <= 74000.0D0) then
          X(25) = X(25) + 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
        else
          X(25) = X(25) + 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
        end if
!     X(25)=X(25)+(RESON1+RESON2+RESON3)*2./3.
!     X(25)=X(25)+37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.

        !  3s² 3p² ³P₁ (level 26)
        if (WAVENO >= ELIM_BLK - ELEV(26)) then
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
          if (WAVENO <= 74000.0D0) then
            X(26) = X(26) + 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
          else
            X(26) = X(26) + 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
          end if
!     X(26)=X(26)+(RESON1+RESON2+RESON3)*2./3.
!     X(26)=37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.

          !  3s² 3p² ³P₀ (level 27)
          if (WAVENO >= ELIM_BLK - ELEV(27)) then
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
            if (WAVENO <= 74000.0D0) then
              X(27) = X(27) + 72.D-18 * (65524.393D0 / WAVENO)**1.90D0 * 2.0D0 / 3.0D0
            else
              X(27) = X(27) + 93.D-18 * (65524.393D0 / WAVENO)**4.00D0 * 2.0D0 / 3.0D0
            end if
!     X(27)=X(27)+(RESON1+RESON2+RESON3)*2./3.
!     X(27)=X(27)+37.E-18*MIN(1.D0,(74000./WAVENO)**5)*2./3.
          end if  ! level 27
        end if  ! level 26
      end if  ! level 25
    end if  ! level 24
  end if  ! level 23

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

  do I = 28, NLEV
    if (WAVENO < ELIM_BLK - ELEV(I)) exit
    ZEFF2 = 9.0D0 / RYD * (ELIM_BLK - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, 3, 1) * DEGEN
  end do

  ! =====================================================================
  ! Final assembly: sum over all levels with Boltzmann populations
  ! =====================================================================
  ELIM_BLK = 65747.55D0
  GFACTOR = 6.0D0

  do J = 1, NRHOX
    ! High-n hydrogenic contribution (n=5 to infinity)
    H = FREQ3 * GFACTOR * 2.0D0 / 2.0D0 / (RYD * Z**2 * HCKT(J)) &
      * (exp(-max(ELIM_BLK - RYD * Z**2 / 5.0D0**2, ELIM_BLK - WAVENO) * HCKT(J)) &
      - exp(-ELIM_BLK * HCKT(J)))
    ! Sum resolved levels
    do I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    end do
    ASI1(J) = H * XNFP(J, 105) * STIM(J) / RHO(J)
  end do
  return

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

  implicit none

  integer, parameter :: NLEV = 48

  real*8, parameter :: G(48) = (/ &
    25.D0, 35.D0, 21.D0, 15.D0,  9.D0, 35.D0, 33.D0, 21.D0, &
    27.D0, 49.D0,  9.D0, 21.D0, 27.D0,  9.D0,  9.D0, 25.D0, &
    33.D0, 15.D0, 35.D0,  3.D0,  5.D0, 11.D0, 15.D0, 13.D0, &
    15.D0,  9.D0, 21.D0, 15.D0, 21.D0, 25.D0, 35.D0,  9.D0, &
     5.D0, 45.D0, 27.D0, 21.D0, 15.D0, 21.D0, 15.D0, 25.D0, &
    21.D0, 35.D0,  5.D0, 15.D0, 45.D0, 35.D0, 55.D0, 25.D0 /)

  real*8, parameter :: E(48) = (/ &
      500.D0,  7500.D0, 12500.D0, 17500.D0, 19000.D0, 19500.D0, &
    19500.D0, 21000.D0, 22000.D0, 23000.D0, 23000.D0, 24000.D0, &
    24000.D0, 24500.D0, 24500.D0, 26000.D0, 26500.D0, 26500.D0, &
    27000.D0, 27500.D0, 28500.D0, 29000.D0, 29500.D0, 29500.D0, &
    29500.D0, 30000.D0, 31500.D0, 31500.D0, 33500.D0, 33500.D0, &
    34000.D0, 34500.D0, 34500.D0, 35000.D0, 35500.D0, 37000.D0, &
    37000.D0, 37000.D0, 38500.D0, 40000.D0, 40000.D0, 41000.D0, &
    41000.D0, 43000.D0, 43000.D0, 43000.D0, 43000.D0, 44000.D0 /)

  real*8, parameter :: WNO(48) = (/ &
    63500.D0, 58500.D0, 53500.D0, 59500.D0, 45000.D0, 44500.D0, &
    44500.D0, 43000.D0, 58000.D0, 41000.D0, 54000.D0, 40000.D0, &
    40000.D0, 57500.D0, 55500.D0, 38000.D0, 57500.D0, 57500.D0, &
    37000.D0, 54500.D0, 53500.D0, 55000.D0, 34500.D0, 34500.D0, &
    34500.D0, 34000.D0, 32500.D0, 32500.D0, 32500.D0, 32500.D0, &
    32000.D0, 29500.D0, 29500.D0, 31000.D0, 30500.D0, 29000.D0, &
    27000.D0, 54000.D0, 27500.D0, 24000.D0, 47000.D0, 23000.D0, &
    44000.D0, 42000.D0, 42000.D0, 21000.D0, 42000.D0, 42000.D0 /)

  real*8  :: XSECT
  integer :: I, J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING FE1OP'

  ! Zero output
  do J = 1, NRHOX
    AFE1(J) = 0.0D0
  end do

  ! No contribution below 21000 cm⁻¹
  if (WAVENO < 21000.0D0) return

  ! Sum over all 48 levels
  do I = 1, NLEV
    if (WNO(I) > WAVENO) cycle
    XSECT = 3.0D-18 / (1.0D0 + ((WNO(I) + 3000.0D0 - WAVENO) &
          / WNO(I) / 0.1D0)**4)
    do J = 1, NRHOX
      AFE1(J) = AFE1(J) + XSECT * G(I) * EXP(-E(I) * HCKT(J))
    end do
  end do

  ! Apply population, stimulated emission, and density factors
  do J = 1, NRHOX
    AFE1(J) = AFE1(J) * STIM(J) * XNFP(J, 351) / RHO(J)
    SFE1(J) = BNU(J)
  end do

  return

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

  implicit none

  integer, intent(in) :: J
  real*8  :: CHOP

  ! --- Cross-section × partition function table ---
  ! CROSSCH(IT, N): log10(σ·Q) for temperature index IT and energy index N
  !   IT = 1–15:  T = 2000, 2500, ..., 9000 K
  !   N  = 1–105: E = 0.1, 0.2, ..., 10.5 eV  (only N=20–105 used)
  ! Loaded from crossch.dat on first call.
  real*8,    save :: CROSSCH(15,105)
  logical,   save :: CROSSCH_LOADED = .false.
  integer :: II, JJ

  real*8, parameter :: LN10 = 2.30258509299405D0

  ! CH partition function table: T = 1000 to 9000 K in 200 K steps (41 points)
  real*8, parameter :: PARTCH(41) = (/ &
    203.741D0,  249.643D0,  299.341D0,  353.477D0,  412.607D0,  477.237D0,  547.817D0, &
    624.786D0,  708.543D0,  799.463D0,  897.912D0, 1004.227D0, 1118.738D0, 1241.761D0, &
   1373.588D0, 1514.481D0, 1664.677D0, 1824.394D0, 1993.801D0, 2173.050D0, 2362.234D0, &
   2561.424D0, 2770.674D0, 2989.930D0, 3219.204D0, 3458.378D0, 3707.355D0, 3966.005D0, &
   4234.155D0, 4511.604D0, 4798.135D0, 5093.554D0, 5397.593D0, 5709.948D0, 6030.401D0, &
   6358.646D0, 6694.379D0, 7037.313D0, 7387.147D0, 7743.579D0, 8106.313D0 /)

  ! --- Cached frequency-interpolated cross-sections ---
  real*8,  save :: CROSSCHT(15)     ! cross-section slice at current energy
  real*8,  save :: FREQ1 = 0.0D0    ! cached frequency
  integer, save :: N_CACHED = 0     ! cached energy index

  ! --- Local variables ---
  real*8  :: EVOLT, EN, TJ, PART, TN
  integer :: N, IT_PART, IT_XSEC

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING CHOP'

  ! --- Load cross-section table from file on first call ---
  if (.not. CROSSCH_LOADED) then
    open(unit=89, file=trim(DATADIR)//'crossch.dat', status='OLD', action='READ', iostat=IT_XSEC)
    if (IT_XSEC /= 0) then
      write(6, '(A)') ' CHOP: ERROR opening ' // trim(DATADIR) // 'crossch.dat'
      stop 'CHOP: crossch.dat not found'
    end if
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')  ! skip header
    read(89, '(A)') ; read(89, '(A)')                     ! skip header
    do JJ = 1, 105
      read(89, *) (CROSSCH(II, JJ), II = 1, 15)
    end do
    close(89)
    CROSSCH_LOADED = .true.
  end if

  CHOP = 0.0D0

  !---------------------------------------------------------------------
  ! Interpolate cross-section table in energy (cached on frequency)
  !---------------------------------------------------------------------
  if (FREQ /= FREQ1) then
    FREQ1 = FREQ
    EVOLT = WAVENO / 8065.479D0
    N = int(EVOLT * 10.0D0)
    N_CACHED = N
    EN = dble(N) * 0.1D0
    if (N < 20 .or. N >= 105) return
    do IT_XSEC = 1, 15
      CROSSCHT(IT_XSEC) = CROSSCH(IT_XSEC, N) &
        + (CROSSCH(IT_XSEC, N+1) - CROSSCH(IT_XSEC, N)) * (EVOLT - EN) / 0.1D0
    end do
  end if

  !---------------------------------------------------------------------
  ! Temperature-dependent interpolation at depth point J
  !---------------------------------------------------------------------
  TJ = T(J)
  N = N_CACHED
  if (TJ >= 9000.0D0) return
  if (N < 20 .or. N >= 105) return

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

  return

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

  implicit none

  integer, intent(in) :: J
  real*8  :: OHOP

  ! --- Cross-section x partition function table ---
  ! CROSSOH(IT, N): log10(sigma*Q) for temperature index IT and energy index N
  !   IT = 1-15:  T = 2000, 2500, ..., 9000 K
  !   N  = 1-130: E = 2.1, 2.2, ..., 15.0 eV
  ! Loaded from crossoh.dat on first call.
  real*8,    save :: CROSSOH(15,130)
  logical,   save :: CROSSOH_LOADED = .false.
  integer :: II, JJ

  real*8, parameter :: LN10 = 2.30258509299405D0

  ! OH partition function table: T = 1000 to 9000 K in 200 K steps (41 points)
  real*8, parameter :: PARTOH(41) = (/ &
      145.979D0,    178.033D0,    211.618D0,    247.053D0,    284.584D0,    324.398D0,    366.639D0, &
      411.425D0,    458.854D0,    509.012D0,    561.976D0,    617.823D0,    676.626D0,    738.448D0, &
      803.363D0,    871.437D0,    942.735D0,   1017.330D0,   1095.284D0,   1176.654D0,   1261.510D0, &
     1349.898D0,   1441.875D0,   1537.483D0,   1636.753D0,   1739.733D0,   1846.434D0,   1956.883D0, &
     2071.080D0,   2189.029D0,   2310.724D0,   2436.155D0,   2565.283D0,   2698.103D0,   2834.571D0, &
     2974.627D0,   3118.242D0,   3265.366D0,   3415.912D0,   3569.837D0,   3727.077D0 /)

  ! --- Cached frequency-interpolated cross-sections ---
  real*8,  save :: CROSSOHT(15)     ! cross-section slice at current energy
  real*8,  save :: FREQ1 = 0.0D0    ! cached frequency
  integer, save :: N_CACHED = 0     ! cached energy index

  ! --- Local variables ---
  real*8  :: EVOLT, EN, TJ, PART, TN
  integer :: N, IT_PART, IT_XSEC

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING OHOP'

  ! --- Load cross-section table from file on first call ---
  if (.not. CROSSOH_LOADED) then
    open(unit=89, file=trim(DATADIR)//'crossoh.dat', status='OLD', action='READ', iostat=IT_XSEC)
    if (IT_XSEC /= 0) then
      write(6, '(A)') ' OHOP: ERROR opening ' // trim(DATADIR) // 'crossoh.dat'
      stop 'OHOP: crossoh.dat not found'
    end if
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')  ! skip header
    read(89, '(A)') ; read(89, '(A)')                     ! skip header
    do JJ = 1, 130
      read(89, *) (CROSSOH(II, JJ), II = 1, 15)
    end do
    close(89)
    CROSSOH_LOADED = .true.
  end if

  OHOP = 0.0D0

  !---------------------------------------------------------------------
  ! Interpolate cross-section table in energy (cached on frequency)
  ! OHOP energy index is offset: N = EVOLT*10 - 20  (E starts at 2.1 eV)
  !---------------------------------------------------------------------
  if (FREQ /= FREQ1) then
    FREQ1 = FREQ
    EVOLT = WAVENO / 8065.479D0
    N = int(EVOLT * 10.0D0 - 20.0D0)
    N_CACHED = N
    EN = dble(N) * 0.1D0 + 2.0D0
    if (N <= 0 .or. N >= 130) return
    do IT_XSEC = 1, 15
      CROSSOHT(IT_XSEC) = CROSSOH(IT_XSEC, N) &
        + (CROSSOH(IT_XSEC, N+1) - CROSSOH(IT_XSEC, N)) * (EVOLT - EN) / 0.1D0
    end do
  end if

  !---------------------------------------------------------------------
  ! Temperature-dependent interpolation at depth point J
  !---------------------------------------------------------------------
  TJ = T(J)
  N = N_CACHED
  if (TJ >= 9000.0D0) return
  if (N <= 0 .or. N >= 130) return

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

  return

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

  implicit none

  real*8, intent(out) :: AH2COLL(kw)

  ! --- CIA tables: log10 of absorption coefficient ---
  ! Read from h2collop.dat on first call
  ! Borysow, Jorgensen & Zheng (1997, A&A 324, 185)
  real*8, save :: H2HE(7,81)
  real*8, save :: H2H2(7,81)

  ! --- Persistent state ---
  integer, save :: ITEMP1 = 0
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*8  :: H2HENU(7), H2H2NU(7)
  real*8  :: DELNU, DELT, XH2H2, XH2HE
  integer :: J, IT, NU, I
  character(256) :: LINE

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING H2COLLOP'

  !---------------------------------------------------------------------
  ! Read CIA tables from file on first call
  !---------------------------------------------------------------------
  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'h2collop.dat', status='OLD', action='READ')
    ! Skip comment lines starting with #
    do
      read(89, '(A)') LINE
      if (LINE(1:1) /= '#') then
        backspace(89)
        exit
      end if
    end do
    ! Read H2HE table (81 lines of 7 values)
    do I = 1, 81
      read(89, *) H2HE(:, I)
    end do
    ! Read H2H2 table (81 lines of 7 values)
    do I = 1, 81
      read(89, *) H2H2(:, I)
    end do
    close(89)
    INITIALIZED = .true.
  end if

  !---------------------------------------------------------------------
  ! Recompute H2 equilibrium number densities when T structure changes
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do J = 1, NRHOX
      XNH2(J) = 0.0D0
      if (T(J) <= 20000.0D0) then
        XNH2(J) = (XNFP(J,1) * 2.0D0 * BHYD(J,1))**2 * EQUILH2(T(J))
      end if
    end do
  end if

  !---------------------------------------------------------------------
  ! Return zero for wavenumber > 20000 cm-1 (above table range)
  !---------------------------------------------------------------------
  if (WAVENO > 20000.0D0) then
    do J = 1, NRHOX
      AH2COLL(J) = 0.0D0
    end do
    return
  end if

  !---------------------------------------------------------------------
  ! Interpolate tables in wavenumber
  !---------------------------------------------------------------------
  NU = int(WAVENO / 250.0D0) + 1
  NU = min(NU, 80)
  DELNU = (WAVENO - 250.0D0 * dble(NU - 1)) / 250.0D0

  do IT = 1, 7
    H2H2NU(IT) = H2H2(IT, NU) * (1.0D0 - DELNU) + H2H2(IT, NU+1) * DELNU
    H2HENU(IT) = H2HE(IT, NU) * (1.0D0 - DELNU) + H2HE(IT, NU+1) * DELNU
  end do

  !---------------------------------------------------------------------
  ! Interpolate in temperature and assemble opacity at each depth
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    IT = int(T(J) / 1000.0D0)
    IT = max(1, min(6, IT))
    DELT = (T(J) - 1000.0D0 * dble(IT)) / 1000.0D0
    DELT = max(0.0D0, min(1.0D0, DELT))

    ! Temperature interpolation (weights corrected from original)
    XH2H2 = H2H2NU(IT) * (1.0D0 - DELT) + H2H2NU(IT+1) * DELT
    XH2HE = H2HENU(IT) * (1.0D0 - DELT) + H2HENU(IT+1) * DELT

    AH2COLL(J) = (10.0D0**XH2HE * XNF(J,3) + 10.0D0**XH2H2 * XNH2(J)) &
               * XNH2(J) / RHO(J) * STIM(J)
  end do

  return

END SUBROUTINE H2COLLOP

!=======================================================================
! WARMOP: Lukewarm opacity assembly (7000-12000 K)
!=======================================================================

SUBROUTINE WARMOP
  implicit none
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING WARMOP'
  ! Lukewarm metal opacities: N I, O I, Ca II + C II, Mg II, Si II
  call C2OP
  call MG2OP
  call SI2OP
  do J = 1, NRHOX
    ALUKE(J) = (N1OP(J) * XNFP(J, 28) + O1OP(J) * XNFP(J, 36) &
             + CA2OP(J) * XNFP(J, 211)) * STIM(J) / RHO(J) &
             + AC2(J) + AMG2(J) + ASI2(J)
  end do
  return

END SUBROUTINE WARMOP

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

  implicit none

  integer, parameter :: NLEV = 34

  ! --- Atomic data: energy levels (cm-1) and statistical weights ---
  real*8, parameter :: ELEV(NLEV) = (/ &
    179073.05D0, 178955.94D0, 178495.47D0, 175292.30D0, 173347.84D0, &  ! 5g,5f,5d,5p,5s
    168978.34D0, 168124.17D0, 162522.34D0, 157234.07D0,              &  ! 4f,4d,4p,4s
    145550.10D0, 131731.80D0, 116537.65D0,    42.28D0,               &  ! 3d,3p,3s, 2p(inactive)
    202188.07D0, 199965.31D0, 198856.92D0, 198431.96D0, 196572.80D0, &  ! 2s2p3d states
    195786.71D0, 190000.00D0, 188601.54D0, 186452.13D0, 184690.98D0, &  ! 2s2p3d/3p states
    182036.89D0, 181741.65D0, 177787.22D0, 167009.29D0,              &  ! 2s2p3p/3s states
    110651.76D0,  96493.74D0,  74931.11D0,  43035.80D0,              &  ! 2s2p2 (inactive)
    230407.20D0, 150464.60D0, 142027.10D0 /)                            ! 2p3 states

  real*8, parameter :: GLEV(NLEV) = (/ &
    18.D0, 14.D0, 10.D0,  6.D0,  2.D0, 14.D0, 10.D0,  6.D0,  1.D0, &  ! NOTE: GLEV(9)=1 for 2S, expected 2
    10.D0,  6.D0,  1.D0,  3.D0,                                      &  ! NOTE: GLEV(12)=1 for 2S, expected 2
     6.D0, 10.D0, 12.D0, 10.D0, 20.D0, 28.D0,  2.D0, 10.D0, 12.D0, &
     4.D0,  6.D0, 20.D0,  6.D0, 12.D0,                               &
     6.D0,  2.D0, 10.D0, 12.D0,                                      &
     6.D0, 10.D0,  4.D0 /)

  ! Quantum numbers (n,l) for XKARZAS calls
  integer, parameter :: NQ(NLEV) = (/ &
    5,5,5,5,5,  4,4,4,4,  3,3,3, 0,  &  ! Block 1 (level 13 inactive)
    3,3,3,3,3,3,  3,3,3,3,3,3,  3,3,  &  ! Block 2
    0,0,0,0,  2,2,2 /)                    ! levels 28-31 inactive, Block 4

  integer, parameter :: LQ(NLEV) = (/ &
    4,3,2,1,0,  3,2,1,0,  2,1,0, 0,  &  ! Block 1
    2,2,2,2,2,2,  1,1,1,1,1,1,  0,0,  &  ! Block 2
    0,0,0,0,  1,1,1 /)                    ! Block 4

  ! Effective charge prefactor: n² for each level (Z²eff = NPRE/RYD*(ELIM-E))
  real*8, parameter :: NPRE(NLEV) = (/ &
    25.D0,25.D0,25.D0,25.D0,25.D0,  16.D0,16.D0,16.D0,16.D0, &  ! Block 1
     9.D0, 9.D0, 9.D0,  0.D0,                                 &
     9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0,                     &  ! Block 2
     9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0, 9.D0,        &
     0.D0, 0.D0, 0.D0, 0.D0,                                  &  ! inactive
     4.D0, 4.D0, 4.D0 /)                                         ! Block 4

  ! Degeneracy multiplier for inner-shell levels
  real*8, parameter :: DEGEN(NLEV) = (/ &
    1.D0,1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0,1.D0, 1.D0, &
    1.D0,1.D0,1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0,1.D0,1.D0,1.D0,1.D0, 1.D0,1.D0, &
    1.D0,1.D0,1.D0,1.D0, &
    3.D0,3.D0,3.D0 /)

  ! Ionization limits (cm-1)
  real*8, parameter :: ELIM1 = 196664.70D0                 ! C III 2s2 1S0
  real*8, parameter :: ELIM2 = 196664.70D0 + 52367.06D0   ! C III 2s2p 3P0
  real*8, parameter :: ELIM4 = 196664.70D0 + 137425.70D0  ! C III 2p2 3P0

  real*8, parameter :: RYD = 109732.298D0
  real*8, parameter :: Z4 = 16.0D0    ! Z**4 (Z=2 for C II)
  real*8, parameter :: Z2 = 4.0D0     ! Z**2

  ! --- Persistent state ---
  real*8,  save :: BOLT(NLEV, kw)
  integer, save :: ITEMP1 = 0

  ! --- Local variables ---
  real*8  :: X(NLEV)
  real*8  :: ELIM, ZEFF2, FREQ3, H
  integer :: I, J, K

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING C2OP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors when temperature structure changes
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      do I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      end do
    end do
  end if

  !---------------------------------------------------------------------
  ! Compute cross-sections at the current frequency
  !---------------------------------------------------------------------
  X(:) = 0.0D0

  ! --- Block 1: excited states -> C III 2s2 1S0 ---
  ELIM = ELIM1
  do I = 1, 12
    if (WAVENO < ELIM - ELEV(I)) exit
    ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
  end do
  ! Level 13 (2s2 2p 2P) handled in HOTOP

  ! --- Block 2: inner-shell states -> C III 2s2p 3P0 ---
  ELIM = ELIM2
  do I = 14, 27
    if (WAVENO < ELIM - ELEV(I)) exit
    ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
  end do
  ! Levels 28-31 (2s2p2 states) handled in HOTOP

  ! --- Block 4: 2p3 states -> C III 2p2 3P0 ---
  ELIM = ELIM4
  do I = 32, 34
    if (WAVENO < ELIM - ELEV(I)) exit
    ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
    X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I)) * DEGEN(I)
  end do

  !---------------------------------------------------------------------
  ! Assemble opacity over depth: levels + dissolved high-n contributions
  !---------------------------------------------------------------------
  FREQ3 = 2.815D29 / (FREQ * FREQ * FREQ) * Z4

  do J = 1, NRHOX
    ! Dissolved levels, Block 1: n >= 6 to infinity, GFACTOR = 1
    H = FREQ3 * 1.0D0 * 2.0D0 / 2.0D0 / (RYD * Z2 * HCKT(J)) &
      * (exp(-max(ELIM1 - RYD * Z2 / 36.0D0, ELIM1 - WAVENO) * HCKT(J)) &
       - exp(-ELIM1 * HCKT(J)))

    ! Dissolved levels, Block 2: n >= 4 to infinity, GFACTOR = 9
    H = H + FREQ3 * 9.0D0 * 2.0D0 / 2.0D0 / (RYD * Z2 * HCKT(J)) &
      * (exp(-max(ELIM2 - RYD * Z2 / 16.0D0, ELIM2 - WAVENO) * HCKT(J)) &
       - exp(-ELIM2 * HCKT(J)))

    ! Sum bound-free over all active levels
    do I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    end do

    AC2(J) = H * XNFP(J, 22) * STIM(J) / RHO(J)
  end do

  return

END SUBROUTINE C2OP

!=======================================================================
! N1OP: N I bound-free opacity cross section
!=======================================================================

REAL*8 FUNCTION N1OP(J)

  implicit none
  integer, intent(in) :: J

  ! N I cross-section times partition function
  ! Three Seaton-formula edges: 853A (ground), 1020A, 1130A (excited)
  real*8, save :: X853 = 0.0D0, X1020 = 0.0D0, X1130 = 0.0D0
  real*8, save :: C1130(kw), C1020(kw)
  real*8, save :: FREQ1 = 0.0D0
  integer, save :: ITEMP1 = 0
  integer :: K

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING N1OP'

  ! Recompute Boltzmann factors when temperature changes
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      C1130(K) = 6.0D0 * exp(-3.575D0 / TKEV(K))
      C1020(K) = 10.0D0 * exp(-2.384D0 / TKEV(K))
    end do
  end if

  ! Recompute cross-sections when frequency changes
  if (FREQ /= FREQ1) then
    X1130 = 0.0D0
    X1020 = 0.0D0
    X853 = 0.0D0
    if (FREQ >= 3.517915D15) X853 = SEATON(3.517915D15, 1.142D-17, 2.0D0, 4.29D0)
    if (FREQ >= 2.941534D15) X1020 = SEATON(2.941534D15, 4.41D-18, 1.5D0, 3.85D0)
    if (FREQ >= 2.653317D15) X1130 = SEATON(2.653317D15, 4.2D-18, 1.5D0, 4.34D0)
    FREQ1 = FREQ
  end if

  N1OP = X853 * 4.0D0 + X1020 * C1020(J) + X1130 * C1130(J)
  return

END FUNCTION N1OP

!=======================================================================
! O1OP
!=======================================================================

FUNCTION O1OP(J)

  implicit none
  integer, intent(in) :: J
  real*8 :: O1OP

  ! O I bound-free cross-section times partition function (after Peach)
  real*8, save :: X911 = 0.0D0
  real*8, save :: FREQ1 = 0.0D0

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING O1OP'

  ! Recompute cross-section only when frequency changes
  if (FREQ /= FREQ1) then
    X911 = 0.0D0
    if (FREQ >= FREQ_RYDH) X911 = SEATON(FREQ_RYDH, 2.94D-18, 1.0D0, 2.66D0)
    FREQ1 = FREQ
  end if

  O1OP = X911 * 9.0D0
  return

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

  implicit none

  integer, parameter :: NLEV = 14

  ! --- Atomic data: energy levels (cm-1) and statistical weights ---
  real*8, parameter :: ELEV(NLEV) = (/ &
    112197.00D0, 108900.00D0, 103705.66D0, 103689.89D0, 103419.82D0, &  ! n=7,6, 5g,5f,5d
     97464.32D0,  92790.51D0,  93799.70D0,  93310.80D0,  80639.85D0, &  ! 5p,5s, 4f,4d, 4p
     69804.95D0,  71490.54D0,  35730.36D0,      0.00D0 /)                ! 4s, 3d, 3p, 3s

  ! GLEV(11) corrected from 22 to 2 (4s 2S: g=2J+1=2)
  real*8, parameter :: GLEV(NLEV) = (/ &
    98.D0, 72.D0, 18.D0, 14.D0, 10.D0, 6.D0, 2.D0, &
    14.D0, 10.D0,  6.D0,  2.D0, 10.D0, 6.D0, 2.D0 /)

  ! Quantum numbers (n,l) for XKARZAS calls
  ! Levels 1-2: l=n convention for super-levels
  ! Level 10: l=1 (Stift 2009 correction from l=2)
  integer, parameter :: NQ(NLEV) = (/ 7,6, 5,5,5,5,5, 4,4,4,4, 3,3, 0 /)
  integer, parameter :: LQ(NLEV) = (/ 7,6, 4,3,2,1,0, 3,2,1,0, 2,1, 0 /)

  ! Effective charge prefactor: n2 for each level
  real*8, parameter :: NPRE(NLEV) = (/ &
    49.D0, 36.D0, 25.D0, 25.D0, 25.D0, 25.D0, 25.D0, &
    16.D0, 16.D0, 16.D0, 16.D0,  9.D0,  9.D0,  0.D0 /)

  ! Ionization limit (cm-1): Mg III 2p6 1S0
  real*8, parameter :: ELIM = 121267.61D0
  real*8, parameter :: RYD  = 109732.298D0
  real*8, parameter :: Z2   = 4.0D0     ! Z**2 (Z=2 for Mg II)
  real*8, parameter :: Z4   = 16.0D0    ! Z**4

  ! --- Persistent state ---
  real*8,  save :: BOLT(NLEV, kw)
  integer, save :: ITEMP1 = 0

  ! --- Local variables ---
  real*8  :: X(NLEV)
  real*8  :: ZEFF2, FREQ3, H, RATIO
  integer :: I, J, K

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING MG2OP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors when temperature structure changes
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      do I = 1, NLEV
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKT(K))
      end do
    end do
  end if

  !---------------------------------------------------------------------
  ! Compute cross-sections at the current frequency
  !---------------------------------------------------------------------
  X(:) = 0.0D0

  levels: do

    ! Levels 1-13: hydrogenic via XKARZAS
    do I = 1, 13
      if (WAVENO < ELIM - ELEV(I)) exit levels
      ZEFF2 = NPRE(I) / RYD * (ELIM - ELEV(I))
      X(I) = XKARZAS(FREQ, ZEFF2, NQ(I), LQ(I))
    end do

    ! Level 14: 3s 2S ground state — non-hydrogenic power-law fit
    if (WAVENO < ELIM - ELEV(14)) exit levels
    RATIO = (ELIM - ELEV(14)) / WAVENO
    X(14) = 0.14D-18 * (6.700D0 * RATIO**4 - 5.700D0 * RATIO**5)

    exit levels
  end do levels

  !---------------------------------------------------------------------
  ! Assemble opacity over depth: levels + dissolved high-n contribution
  !---------------------------------------------------------------------
  FREQ3 = 2.815D29 / (FREQ * FREQ * FREQ) * Z4

  do J = 1, NRHOX
    ! Dissolved levels: n >= 8 to infinity
    ! GFACTOR = 2 (implicit: 2./2. = 1 in prefactor)
    H = FREQ3 * 2.0D0 / 2.0D0 / (RYD * Z2 * HCKT(J)) &
      * (exp(-max(ELIM - RYD * Z2 / 64.0D0, ELIM - WAVENO) * HCKT(J)) &
       - exp(-ELIM * HCKT(J)))

    do I = 1, NLEV
      H = H + X(I) * BOLT(I, J)
    end do

    AMG2(J) = H * XNFP(J, 79) * STIM(J) / RHO(J)
  end do

  return

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

  implicit none

  integer, parameter :: NWGRID = 200, NTGRID = 51
  real*8, parameter :: RYD = 109732.298D0, Z = 2.0D0

  ! Ionization limits (cm⁻¹) for each level
  real*8, parameter :: ELIMLEV(46) = (/ &
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
  real*8, parameter :: ELEV(46) = (/ &
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
  real*8, parameter :: TLEV(46) = (/ &
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
  real*8, parameter :: GLEV(46) = (/ &
    18.D0, 14.D0, 10.D0, 6.D0, 2.D0, 14.D0, 10.D0, 6.D0, 10.D0, &
    1.D0, 6.D0, 20.D0, 10.D0, 18.D0, 36.D0, 28.D0, 10.D0, 10.D0, &
    6.D0, 12.D0, 2.D0, 20.D0, 28.D0, 10.D0, 10.D0, 4.D0, 12.D0, &
    6.D0, 20.D0, 10.D0, 6.D0, 12.D0, 20.D0, 6.D0, 12.D0, 28.D0, &
    10.D0, 6.D0, 2.D0, 10.D0, 12.D0, 6.D0, 10.D0, 4.D0, &
    1.D0, 9.D0 /)

  ! Angular momentum quantum numbers
  integer, parameter :: LLEV(46) = (/ &
    4, 3, 2, 1, 0, 3, 2, 1, 2, 0, 1, &
    3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, &
    2, 2, 2, 0, 0, 2, 2, 1, 1, 1, 1, &
    1, 1, 1, 0, 0 /)

  ! Principal quantum numbers
  integer, parameter :: NQLEV(46) = (/ &
    5, 5, 5, 5, 5, 4, 4, 4, 3, 4, 3, &
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, &
    3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 6, 5 /)

  ! Pretabulated opacity grid and supporting arrays
  real*8,  save :: XTAB(NWGRID, NTGRID) = 0.0D0
  real*8,  save :: HCKTTAB(NTGRID), BOLT3s2(NTGRID), BOLT3s3p(NTGRID)
  real*8,  save :: BOLT(44, NTGRID)
  real*8,  save :: BOLTN(kw)
  integer, save :: INDEXT(kw)
  real*8,  save :: TFRAC(kw)
  integer, save :: ITEMP1 = 0
  logical, save :: INITIALIZED = .false.

  real*8  :: ZEFF2LEV(44)
  real*8  :: X(46), FREQ3, WNOTAB, FREQTAB, TTAB, H, WFRAC, TLOG10
  integer :: I, J, K, NU, IT, IW

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING SI2OP'

  ! ==================================================================
  ! One-time pretabulation of opacity grid
  ! ==================================================================
  if (.not. INITIALIZED) then
    ! Compute effective charges for levels 1-44
    do I = 1, 44
      ZEFF2LEV(I) = dble(NQLEV(I))**2 / RYD * TLEV(I)
    end do

    ! Build temperature grid
    do K = 1, NTGRID
      TTAB = 10.0D0**(3.48D0 + K * 0.02D0)
      HCKTTAB(K) = HCK / TTAB
      BOLT3s2(K) = exp(-ELIMLEV(1) * HCKTTAB(K))
      BOLT3s3p(K) = exp(-ELIMLEV(12) * HCKTTAB(K))
      do I = 1, 44
        BOLT(I, K) = GLEV(I) * exp(-ELEV(I) * HCKTTAB(K))
      end do
    end do

    ! Build opacity table: loop over wavenumber grid
    do NU = 1, NWGRID
      WNOTAB = dble(NU) * 1000.0D0
      FREQTAB = WNOTAB * CLIGHT
      FREQ3 = 2.815D29 / FREQTAB / FREQTAB / FREQTAB * Z**4

      do I = 1, 46
        X(I) = 0.0D0
      end do

      ! Levels 1-11: → Si III 3s²
      do I = 1, 11
        if (WNOTAB < TLEV(I)) exit
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I))
      end do

      ! Levels 12-37: → Si III 3s3p ³P₀
      do I = 12, 37
        if (WNOTAB < TLEV(I)) exit
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I))
      end do

      ! Levels 38-41: → Si III 3s3p ³P₀ (×2 degeneracy)
      do I = 38, 41
        if (WNOTAB < TLEV(I)) exit
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I)) * 2.0D0
      end do

      ! Levels 42-44: → Si III 3p² ³P₀ (×3 degeneracy)
      do I = 42, 44
        if (WNOTAB < TLEV(I)) exit
        X(I) = XKARZAS(FREQTAB, ZEFF2LEV(I), NQLEV(I), LLEV(I)) * 3.0D0
      end do

      ! Assemble opacity for each temperature point
      do K = 1, NTGRID
        ! High-n series: 3s² n≥6
        H = FREQ3 * GLEV(45) * 2.0D0 / 2.0D0 / (RYD * Z**2 * HCKTTAB(K)) &
          * (exp(-max(ELEV(45), ELIMLEV(45) - WNOTAB) * HCKTTAB(K)) &
          - BOLT3s2(K))
        ! High-n series: 3s3p n≥5
        H = H + FREQ3 * GLEV(46) * 2.0D0 / 2.0D0 / (RYD * Z**2 * HCKTTAB(K)) &
          * (exp(-max(ELEV(46), ELIMLEV(46) - WNOTAB) * HCKTTAB(K)) &
          - BOLT3s3p(K))
        ! Sum resolved levels
        do I = 1, 44
          H = H + X(I) * BOLT(I, K)
        end do
        XTAB(NU, K) = log(H)
      end do
    end do

    INITIALIZED = .true.
  end if

  ! ==================================================================
  ! Temperature-dependent precomputation (when ITEMP changes)
  ! ==================================================================
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do J = 1, NRHOX
      BOLTN(J) = (GLEV(45) * exp(-ELIMLEV(45) * HCKT(J)) &
               + GLEV(46) * exp(-ELIMLEV(46) * HCKT(J))) &
               / (RYD * Z**2 * HCKT(J))
      TLOG10 = TLOG(J) / 2.30258509299405D0
      IT = int((TLOG10 - 3.48D0) / 0.02D0)
      IT = max(min(IT, 50), 1)
      INDEXT(J) = IT
      TFRAC(J) = (TLOG10 - 3.48D0 - IT * 0.02D0) / 0.02D0
    end do
  end if

  ! ==================================================================
  ! Frequency evaluation: bilinear interpolation or low-freq fallback
  ! ==================================================================
  if (WAVENO >= 12192.48D0) then
    ! Bilinear interpolation in pretabulated grid
    IW = int(WAVENO * 0.001D0)
    IW = max(min(IW, 199), 1)
    WFRAC = (WAVENO - IW * 1000.0D0) / 1000.0D0
    do J = 1, NRHOX
      IT = INDEXT(J)
      H = (XTAB(IW, IT) * (1.0D0 - TFRAC(J)) &
        + XTAB(IW, IT+1) * TFRAC(J)) * (1.0D0 - WFRAC) &
        + (XTAB(IW+1, IT) * (1.0D0 - TFRAC(J)) &
        + XTAB(IW+1, IT+1) * TFRAC(J)) * WFRAC
      ASI2(J) = exp(H) * XNFP(J, 106) * STIM(J) / RHO(J)
    end do
  else
    ! Low frequency: only high-n hydrogenic contribution
    ! Si III 3s² (ELIM=131838.4, n≥6) + Si III 3s3p ³P₀ (ELIM=184563.09, n≥5)
    FREQ3 = 2.815D29 / FREQ / FREQ / FREQ * Z**4
    do J = 1, NRHOX
      H = FREQ3 * (1.0D0 / EHVKT(J) - 1.0D0) * BOLTN(J)
      ASI2(J) = H * XNFP(J, 106) * STIM(J) / RHO(J)
    end do
  end if
  return

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

  implicit none
  integer, intent(in) :: J
  real*8 :: CA2OP

  ! Ca II bound-free cross-section times partition function
  ! Three edges: 1044A (ground), 1218A, 1420A (excited)
  real*8, save :: C1218(kw), C1420(kw)
  real*8, save :: X1044 = 0.0D0, X1218 = 0.0D0, X1420 = 0.0D0
  real*8, save :: FREQ1 = 0.0D0
  integer, save :: ITEMP1 = 0
  integer :: K

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING CA2OP'

  ! Recompute Boltzmann factors when temperature changes
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      C1218(K) = 10.0D0 * exp(-1.697D0 / TKEV(K))
      C1420(K) = 6.0D0 * exp(-3.142D0 / TKEV(K))
    end do
  end if

  ! Recompute cross-sections when frequency changes
  if (FREQ /= FREQ1) then
    X1420 = 0.0D0
    X1218 = 0.0D0
    X1044 = 0.0D0
    if (FREQ >= 2.870454D15) X1044 = 5.4D-20 * (2.870454D15 / FREQ)**3
    if (FREQ >= 2.460127D15) X1218 = 1.64D-17 * sqrt(2.460127D15 / FREQ)
    if (FREQ >= 2.110779D15) X1420 = SEATON(2.110779D15, 4.13D-18, 3.0D0, 0.69D0)
    FREQ1 = FREQ
  end if

  CA2OP = X1044 * 2.0D0 + X1218 * C1218(J) + X1420 * C1420(J)
  return

END FUNCTION CA2OP

!=========================================================================
! SUBROUTINE HOTOP
!
! Hot-star opacity assembly (significant for T > ~12000 K).
!
! Combines three contributions into AHOT(J):
!
!   1. Free-free opacity from C, N, O, Ne, Mg, Si, S, Fe (ions I-V)
!      using Coulomb free-free Gaunt factors COULFF(J,Z) for Z=1..5.
!      Formula: sigma_ff = 3.6919e8 * Z^2 * g_ff * n_e / (nu^3 * sqrt(T))
!
!   2. C II bound-free (AC2OP) — currently disabled (handled in LUKE
!      per Bischoff 2003); placeholder zeroed array preserved.
!
!   3. 60 bound-free edges from various species, read from hotop.dat.
!      Each edge has 7 parameters: freq0, sigma0, s, n_power, g, E_exc, id
!      Cross-section (modified Seaton):
!        sigma = sigma0 * (s + freq0/nu - s*freq0/nu) * sqrt((freq0/nu)^n)
!      Applied with stat weight g and Boltzmann factor exp(-E_exc/kT).
!      Optimization: edges contributing < 1% of running total are skipped.
!=========================================================================

SUBROUTINE HOTOP

  implicit none

  integer, parameter :: NUM = 60    ! number of bound-free edges
  integer, parameter :: NPAR = 7    ! parameters per edge

  ! COEFF_FF replaced by COEFF_FF from mod_constants

  ! --- Edge parameters (read from file on first call) ---
  ! A(1:7,I): freq0, sigma0, s, n_power, g_stat, E_exc, species_id
  real*8,  save :: A(NPAR, NUM)
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*8  :: AC2OP(kw)
  real*8  :: FREE, XSECT, XX, FRATIO
  integer :: I, J, ID
  character(256) :: LINE

  ! --- External functions ---

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HOTOP'

  !---------------------------------------------------------------------
  ! Read edge parameters from file on first call
  !---------------------------------------------------------------------
  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'hotop.dat', status='OLD', action='READ')
    do
      read(89, '(A)') LINE
      if (LINE(1:1) /= '#') then
        backspace(89)
        exit
      end if
    end do
    do I = 1, NUM
      read(89, *) A(:, I)
    end do
    close(89)
    INITIALIZED = .true.
  end if

  !---------------------------------------------------------------------
  ! Free-free opacity: C,N,O,Ne,Mg,Si,S,Fe ions (stages I-V)
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    FREE = COULFF(J,1) * 1.0D0 &
            * (XNF(J,22) + XNF(J,29) + XNF(J,37) + XNF(J,56) &
             + XNF(J,79) + XNF(J,106) + XNF(J,137) + XNF(J,352)) &
         + COULFF(J,2) * 4.0D0 &
            * (XNF(J,23) + XNF(J,30) + XNF(J,38) + XNF(J,57) &
             + XNF(J,80) + XNF(J,107) + XNF(J,138) + XNF(J,353)) &
         + COULFF(J,3) * 9.0D0 &
            * (XNF(J,24) + XNF(J,31) + XNF(J,39) + XNF(J,58) &
             + XNF(J,81) + XNF(J,108) + XNF(J,139) + XNF(J,354)) &
         + COULFF(J,4) * 16.0D0 &
            * (XNF(J,25) + XNF(J,32) + XNF(J,40) + XNF(J,59) &
             + XNF(J,82) + XNF(J,109) + XNF(J,140) + XNF(J,355)) &
         + COULFF(J,5) * 25.0D0 &
            * (XNF(J,26) + XNF(J,33) + XNF(J,41) + XNF(J,60) &
             + XNF(J,83) + XNF(J,110) + XNF(J,141) + XNF(J,356))

    AHOT(J) = FREE * COEFF_FF / (FREQ * FREQ * FREQ) * XNE(J) / sqrt(T(J))
  end do

  !---------------------------------------------------------------------
  ! C II bound-free — disabled (handled in LUKE per Bischoff 4 Jun 2003)
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    AC2OP(J) = 0.0D0
  end do
  ! CALL C2OP(AC2OP)
  do J = 1, NRHOX
    AHOT(J) = AHOT(J) + AC2OP(J)
  end do

  !---------------------------------------------------------------------
  ! Bound-free edges (60 modified-Seaton edges from various species)
  !---------------------------------------------------------------------
  do I = 1, NUM
    if (FREQ < A(1, I)) cycle
    FRATIO = A(1, I) / FREQ
    XSECT = A(2, I) * (A(3, I) + FRATIO - A(3, I) * FRATIO) &
          * sqrt(FRATIO ** int(A(4, I)))
    ID = int(A(7, I))
    do J = 1, NRHOX
      XX = XSECT * XNFP(J, ID) * A(5, I)
      if (XX > AHOT(J) / 100.0D0) then
        AHOT(J) = AHOT(J) + XX / exp(A(6, I) / TKEV(J))
      end if
    end do
  end do

  !---------------------------------------------------------------------
  ! Apply stimulated emission and normalize by density
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    AHOT(J) = AHOT(J) * STIM(J) / RHO(J)
  end do

  return

END SUBROUTINE HOTOP

!=======================================================================
! ELECOP
!=======================================================================

SUBROUTINE ELECOP

  implicit none
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING ELECOP'
  ! Thomson electron scattering: sigma_T = 0.6653e-24 cm^2
  do J = 1, NRHOX
    SIGEL(J) = SIGMA_THOMSON * XNE(J) / RHO(J)
  end do
  return

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

  implicit none

  integer, save :: ITEMP1 = 0
  real*8  :: W, WW, SIG
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING H2RAOP'

  ! Recompute H₂ number density when the temperature structure changes.
  ! XNH2 is a function of T(J) only (through EQUILH2) and the H I ground
  ! state population, so it is cached between calls at the same ITEMP.
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do J = 1, NRHOX
      XNH2(J) = 0.0D0
      if (T(J) <= 20000.0D0) then
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
      end if
    end do
  end if

  ! Evaluate the Dalgarno & Williams (1962) cross-section.  Cap the
  ! frequency at 2.922e15 Hz (λ ≈ 1026 Å) to prevent extrapolation
  ! into the H₂ Lyman/Werner band region where the Rayleigh
  ! approximation breaks down.
  W  = CLIGHT_ANG / min(FREQ, 2.922D15)   ! wavelength [Å], capped
  WW = W**2                                ! λ² [Å²]

  ! σ(λ) = (8.14e-13 + 1.28e-6/λ² + 1.61/λ⁴) / λ⁴   [cm²]
  ! Rewritten with WW to share the λ² computation between terms.
  SIG = (8.14D-13 + 1.28D-6 / WW + 1.61D0 / (WW * WW)) / (WW * WW)

  ! Apply to each depth point.  The /RHO(J) factor converts from
  ! cross-section-times-number-density [cm⁻¹] to mass opacity [cm²/g].
  do J = 1, NRHOX
    SIGH2(J) = SIG * XNH2(J) / RHO(J)
  end do
  return

END SUBROUTINE H2RAOP

!=========================================================================
! SUBROUTINE HLINOP
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

  implicit none

  real*8, parameter :: NU_LYMAN = FREQ_RYDH   ! Lyman limit frequency (Hz)

  ! Maximum level for occupation probability computation.
  ! Beyond this, w_n is effectively 0 for any stellar density.
  integer, parameter :: NMAX_OCC = 80

  ! --- Persistent state ---
  real*8,  save :: BOLT(kw, 4)            ! Boltzmann × population for n=1..4
  real*8,  save :: W_OCC(kw, NMAX_OCC)    ! occupation probabilities per depth
  integer, save :: ITEMP1 = 0

  ! --- Local variables ---
  real*8  :: H, S, A, BHYDJM, w_m
  integer :: J, N, M, M1, M2, MFREQ

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HLINOP'

  !---------------------------------------------------------------------
  ! Recompute Boltzmann factors and occupation probabilities when T changes
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    do J = 1, NRHOX
      ! Occupation probabilities for levels 2..NMAX_OCC
      W_OCC(J, 1) = 1.0D0
      do N = 2, NMAX_OCC
        W_OCC(J, N) = occupation_prob(N, XNE(J))
      end do
      ! Boltzmann factors for lower levels 1..4
      do N = 1, 4
        BOLT(J, N) = exp(-(13.595D0 - 13.595D0 / dble(N)**2) / TKEV(J)) &
                   * 2.0D0 * dble(N)**2 * BHYD(J, N) * XNFP(J, 1) / RHO(J)
      end do
    end do
    ITEMP1 = ITEMP
  end if

  !---------------------------------------------------------------------
  ! Determine lower level N from frequency
  !---------------------------------------------------------------------
  N = int(sqrt(NU_LYMAN / FREQ))
  if (N == 0 .or. N > 4) return

  ! Low-frequency cutoffs by series
  select case (N)
  case (1)
    if (FREQ < 2.0D15) return       ! Lyman: lambda > 1500 A
  case (2)
    if (FREQ < 4.44D14) return       ! Balmer: lambda > 6756 A
  end select

  !---------------------------------------------------------------------
  ! Upper level nearest to current frequency
  !---------------------------------------------------------------------
  MFREQ = nint(sqrt(NU_LYMAN / (NU_LYMAN / dble(N)**2 - FREQ)))

  !---------------------------------------------------------------------
  ! Sum Stark-broadened line opacity over depth
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    M1 = MFREQ
    M2 = M1 + 1
    M1 = max(M1, N + 1)
    H = 0.0D0
    S = 0.0D0

    if (M1 > 6) then
      ! High upper levels: check if effectively dissolved
      if (M1 <= NMAX_OCC) then
        if (W_OCC(J, M1) < 1.0D-3) then
          ! Fully dissolved: no line opacity.
          ! The continuum opacity is provided by HOP's pseudo-continuum.
          AHLINE(J) = 0.0D0
          SHLINE(J) = BNU(J)
          cycle
        end if
      else
        ! Beyond tabulated range: fully dissolved
        AHLINE(J) = 0.0D0
        SHLINE(J) = BNU(J)
        cycle
      end if
      M1 = M1 - 1
      M2 = M2 + 3
      ! Special case: add Paschen-alpha (3->4) if computing Brackett (N=4)
      if (N >= 4 .and. M1 <= 8) then
        H = STARK_MMM(3, 4, J) * (1.0D0 - EHVKT(J) * BHYD(J, 4) / BHYD(J, 3)) * BOLT(J, 3)
        S = H * BNU(J) * STIM(J) / (BHYD(J, 3) / BHYD(J, 4) - EHVKT(J))
      end if
    end if

    ! Sum over upper levels in window, weighted by occupation probability.
    ! Only the surviving fraction w_m contributes as line opacity.
    ! The dissolved fraction (1 - w_m) is handled by HOP's pseudo-continuum.
    do M = M1, M2
      BHYDJM = 1.0D0
      if (M <= 6) BHYDJM = BHYD(J, M)

      ! Occupation probability weight for upper level
      w_m = 1.0D0
      if (M >= 2 .and. M <= NMAX_OCC) w_m = W_OCC(J, M)
      if (M > NMAX_OCC) w_m = 0.0D0

      ! Line opacity (surviving fraction only)
      A = w_m * STARK_MMM(N, M, J) * (1.0D0 - EHVKT(J) * BHYDJM / BHYD(J, N)) * BOLT(J, N)
      H = H + A
      S = S + A * BNU(J) * STIM(J) / (BHYD(J, N) / BHYDJM - EHVKT(J))
    end do

    AHLINE(J) = H
    if (H > 0.0D0) then
      SHLINE(J) = S / H
    else
      SHLINE(J) = BNU(J)
    end if
  end do

  return

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

  implicit none

  integer, intent(in) :: N, M, J
  real*8  :: STARK

  ! --- Transition matrix element table: KNMTAB(M-N, N) for M-N=1..5, N=1..4 ---
  real*8, parameter :: KNMTAB(5,4) = reshape((/ &
    .000356D0, .000523D0, .00109D0, .00149D0, .00225D0, &
    .0125D0,   .0177D0,   .028D0,   .0348D0,  .0493D0, &
    .124D0,    .171D0,    .223D0,   .261D0,   .342D0, &
    .683D0,    .866D0,   1.02D0,   1.19D0,   1.46D0 /), shape(KNMTAB))

  ! --- Oscillator strength proxy: FSTARK(M-N, N) for M-N=1..10, N=1..4 ---
  real*8, parameter :: FSTARK(10,4) = reshape((/ &
    .1387D0,  .07910D0, .02126D0,  .01394D0,  .006462D0, &
    .004814D0, .002779D0, .002216D0, .001443D0, .001201D0, &
    .3921D0,  .1193D0,  .03766D0,  .02209D0,  .01139D0, &
    .008036D0, .005007D0, .003850D0, .002658D0, .002151D0, &
    .6103D0,  .1506D0,  .04931D0,  .02768D0,  .01485D0, &
    .01023D0, .006588D0, .004996D0, .003542D0, .002838D0, &
    .8163D0,  .1788D0,  .05985D0,  .03189D0,  .01762D0, &
    .01196D0, .007825D0, .005882D0, .004233D0, .003375D0 /), shape(FSTARK))

  real*8, parameter :: RYD = FREQ_RYDH       ! Rydberg frequency (Hz)
  ! PI from mod_constants; CLIGHT_ANG replaces local CLIGHT (Å·Hz)
  real*8, parameter :: A0 = 0.0265384D0       ! profile normalization constant

  ! --- Persistent state ---
  real*8,  save :: F0(kw)       ! Holtsmark normal field strength
  integer, save :: ITEMP1 = 0

  ! --- Local variables ---
  real*8  :: XN, XM, X, XX, NN, MM
  real*8  :: KNM, FNM, FREQNM, DEL, DBETA, BETA
  real*8  :: Y1, Y2, QSTAT, IMPACT, EXY2
  real*8  :: PROF, RATIO, DIOI
  integer :: K, MMINN

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING STARK'

  !---------------------------------------------------------------------
  ! Recompute Holtsmark field when temperature structure changes
  !---------------------------------------------------------------------
  if (ITEMP /= ITEMP1) then
    do K = 1, NRHOX
      F0(K) = 1.25D-9 * XNE(K)**0.6666667D0
    end do
    ITEMP1 = ITEMP
  end if

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
  if (MMINN <= 5) then
    KNM = KNMTAB(MMINN, N)
  else
    KNM = 5.5D-5 * (NN * MM)**2 / (MM - NN)
  end if

  ! Oscillator strength proxy FNM
  if (MMINN <= 10) then
    FNM = FSTARK(MMINN, N)
  else
    FNM = FSTARK(10, N) * ((20.0D0 * XN + 100.0D0) &
         / ((XN + 10.0D0) * XM * (1.0D0 - XX)))**3
  end if

  !---------------------------------------------------------------------
  ! Line center, detuning, and normalized detuning
  !---------------------------------------------------------------------
  FREQNM = RYD * (1.0D0 / NN - 1.0D0 / MM)
  DEL = abs(FREQ - FREQNM)
  DBETA = CLIGHT_ANG / FREQNM**2 / F0(J) / KNM
  BETA = DBETA * DEL

  !---------------------------------------------------------------------
  ! Quasistatic ion + impact electron broadening
  !---------------------------------------------------------------------
  Y1 = MM * DEL * HKT(J) / 2.0D0
  Y2 = (PI**2 / 2.0D0 / A0 / CLIGHT) * DEL**2 / XNE(J)
  QSTAT = 1.5D0 + 0.5D0 * (Y1**2 - 1.384D0) / (Y1**2 + 1.384D0)

  IMPACT = 0.0D0
  if (Y1 < 8.0D0 .and. Y1 < Y2) then
    EXY2 = 0.0D0
    if (Y2 <= 8.0D0) EXY2 = EXINT(Y2)
    IMPACT = 1.438D0 * sqrt(Y1 * (1.0D0 - XX)) &
           * (0.4D0 * exp(-Y1) + EXINT(Y1) - 0.5D0 * EXY2)
  end if

  !---------------------------------------------------------------------
  ! Profile and ratio assembly
  !---------------------------------------------------------------------
  if (BETA <= 20.0D0) then
    ! Core and intermediate wings
    PROF  = 8.0D0 / (80.0D0 + BETA**3)
    RATIO = QSTAT + IMPACT
  else
    ! Far wings with debye-shielding correction
    PROF = 1.5D0 / BETA / BETA / sqrt(BETA)
    DIOI = 6.28D0 * 1.48D-25 * (2.0D0 * MM * RYD / DEL) * XNE(J) &
         * (sqrt(2.0D0 * MM * RYD / DEL) * (1.3D0 * QSTAT + 0.30D0 * IMPACT) &
          - 3.9D0 * RYD * HKT(J))
    RATIO = QSTAT * min(1.0D0 + DIOI, 1.25D0) + IMPACT
  end if

  STARK = A0 * FNM * PROF * DBETA * RATIO
  return

contains

  !-----------------------------------------------------------------------
  ! Exponential integral E1(x) approximation (Abramowitz & Stegun 5.1.53)
  !-----------------------------------------------------------------------
  pure real*8 function EXINT(X)
    real*8, intent(in) :: X
    EXINT = -log(X) - 0.57516D0 + (0.97996D0 + (-0.21654D0 &
          + (0.033572D0 + (-0.0029222D0 + 1.05439D-4 * X) * X) * X) * X) * X
  end function EXINT

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

  implicit none

  character(len=256) :: filepath
  character(len=*), parameter :: SERIES_FILES(4) = &
    (/ 'stehle_lyman.bin   ', 'stehle_balmer.bin  ', &
       'stehle_paschen.bin ', 'stehle_brackett.bin' /)
  integer :: iseries, iu, n_lower, n_upper_min, n_upper_max
  integer :: n_dens, n_temps, n_dalpha, n_trans, ios

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING INIT_STARK_TABLES'

  do iseries = 1, NSTARK_SERIES

    ! Build file path
    filepath = trim(DATADIR) // '/' // trim(SERIES_FILES(iseries))

    ! Try to open the file
    iu = 200 + iseries
    open(iu, file=trim(filepath), status='old', form='unformatted', &
         access='sequential', iostat=ios)
    if (ios /= 0) then
      write(6,'(A,A)') '  INIT_STARK_TABLES: file not found: ', trim(filepath)
      STEHLE_DATA(iseries)%loaded = .false.
      cycle
    end if

    ! Record 1: dimensions
    read(iu) n_lower, n_upper_min, n_upper_max, n_dens, n_temps, n_dalpha
    n_trans = n_upper_max - n_upper_min + 1

    STEHLE_DATA(iseries)%n_lower     = n_lower
    STEHLE_DATA(iseries)%n_upper_min = n_upper_min
    STEHLE_DATA(iseries)%n_upper_max = n_upper_max
    STEHLE_DATA(iseries)%n_transitions = n_trans
    STEHLE_DATA(iseries)%n_dens      = n_dens

    ! Record 2: density grid
    read(iu) STEHLE_DATA(iseries)%density_grid(1:n_dens)

    ! Record 3: temperature grid
    read(iu) STEHLE_DATA(iseries)%temp_grid(1:n_temps)

    ! Record 4: log Δα grid
    read(iu) STEHLE_DATA(iseries)%log_dalpha_grid(1:n_dalpha)

    ! Records 5-6: per-transition metadata
    allocate(STEHLE_DATA(iseries)%max_dens_idx(n_trans))
    allocate(STEHLE_DATA(iseries)%k_alpha(n_trans))
    read(iu) STEHLE_DATA(iseries)%max_dens_idx(1:n_trans)
    read(iu) STEHLE_DATA(iseries)%k_alpha(1:n_trans)

    ! Record 7: profile array
    allocate(STEHLE_DATA(iseries)%profiles(n_dalpha, n_temps, n_dens, n_trans))
    read(iu) STEHLE_DATA(iseries)%profiles

    close(iu)
    STEHLE_DATA(iseries)%loaded = .true.

  end do

  STEHLE_TABLES_LOADED = .true.

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

  implicit none

  integer, intent(in) :: N, M
  real*8 :: hydrogen_f_value

  ! Tabulated gf = g_n × f_{N→M} for hydrogen
  ! Source: NIST Atomic Spectra Database (Wiese et al.)
  ! g_n = 2N²
  real*8, parameter :: GF_LYMAN(29) = (/ &
    0.8324D0, 0.1580D0, 0.05798D0, 0.02787D0, 0.01551D0, &
    0.009466D0, 0.006158D0, 0.004220D0, 0.003014D0, 0.002225D0, &
    0.001688D0, 0.001312D0, 0.001038D0, 0.000835D0, 0.000681D0, &
    0.000562D0, 0.000469D0, 0.000395D0, 0.000336D0, 0.000288D0, &
    0.000249D0, 0.000216D0, 0.000189D0, 0.000166D0, 0.000146D0, &
    0.000129D0, 0.000115D0, 0.000103D0, 0.000092D0 /)

  real*8, parameter :: GF_BALMER(28) = (/ &
    5.126D0, 0.9543D0, 0.3571D0, 0.1770D0, 0.1023D0, &
    0.06497D0, 0.04394D0, 0.03104D0, 0.02270D0, 0.01708D0, &
    0.01314D0, 0.01028D0, 0.00818D0, 0.00659D0, 0.00537D0, &
    0.00443D0, 0.00369D0, 0.00310D0, 0.00263D0, 0.00225D0, &
    0.00194D0, 0.00168D0, 0.00146D0, 0.00128D0, 0.00113D0, &
    0.00100D0, 0.000889D0, 0.000793D0 /)

  real*8, parameter :: GF_PASCHEN(27) = (/ &
    15.16D0, 2.715D0, 1.001D0, 0.4937D0, 0.2850D0, &
    0.1813D0, 0.1232D0, 0.08804D0, 0.06526D0, 0.04975D0, &
    0.03882D0, 0.03089D0, 0.02498D0, 0.02050D0, 0.01703D0, &
    0.01430D0, 0.01212D0, 0.01036D0, 0.00893D0, 0.00775D0, &
    0.00677D0, 0.00595D0, 0.00526D0, 0.00467D0, 0.00416D0, &
    0.00373D0, 0.00335D0 /)

  real*8, parameter :: GF_BRACKETT(3) = (/ &
    33.22D0, 5.731D0, 2.092D0 /)

  integer :: idx
  real*8  :: gf, xn, xm

  xn = dble(N)
  xm = dble(M)

  if (M <= N) then
    hydrogen_f_value = 0.0D0
    return
  end if

  idx = M - N  ! for Lyman: idx = M-1; for Balmer: idx = M-2; etc.

  if (N == 1 .and. idx <= 29) then
    gf = GF_LYMAN(idx)
  else if (N == 2 .and. idx <= 28) then
    gf = GF_BALMER(idx)
  else if (N == 3 .and. idx <= 27) then
    gf = GF_PASCHEN(idx)
  else if (N == 4 .and. idx <= 3) then
    gf = GF_BRACKETT(idx)
  else
    ! Kramers approximation with empirical correction for higher n
    ! gf ≈ 1.96 × n^2 / (m^3 × (1/n^2 - 1/m^2)^3) × correction
    ! The correction factor accounts for the Gaunt factor
    gf = 1.9603D0 * xn**2 / (xm**3 * (1.0D0/xn**2 - 1.0D0/xm**2)**3)
    ! Apply empirical scale to match exact values at table boundary
    gf = gf * 0.80D0  ! approximate average Gaunt factor
  end if

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

  implicit none

  integer, intent(in) :: N, M, J
  real*8 :: STARK_MMM

  ! Physical constants
  real*8, parameter :: A0 = 0.0265384D0         ! π e² / (m_e c) [cm²/s]
  real*8, parameter :: CLIGHT_ANG = 2.9979D18   ! speed of light [Å/s]
  real*8, parameter :: RYD_ANG = 911.7633455D0  ! Rydberg wavelength [Å]

  ! Local variables
  real*8  :: xn, xm, lambda0, freq_nm, del_freq
  real*8  :: F0, dalpha, log_dalpha, log_ne
  real*8  :: I_dalpha, f_nm
  real*8  :: frac_d, frac_t, frac_a
  real*8  :: v000, v001, v010, v011, v100, v101, v110, v111
  real*8  :: v00, v01, v10, v11, v0, v1
  integer :: iseries, itrans, id1, id2, it1, it2, ia1, ia2
  integer :: nd, nt

  type(stark_series_t), pointer :: S

  ! ---------------------------------------------------------------
  ! Lazy initialization: load tables on first call
  ! ---------------------------------------------------------------
  if (.not. STEHLE_TABLES_LOADED) then
    call INIT_STARK_TABLES
  end if

  ! ---------------------------------------------------------------
  ! Determine which series this transition belongs to
  ! ---------------------------------------------------------------
  if (N < 1 .or. N > 4 .or. M <= N) then
    STARK_MMM = STARK(N, M, J)
    return
  end if
  iseries = N

  ! Check if tables are loaded for this series
  if (.not. STEHLE_TABLES_LOADED .or. .not. STEHLE_DATA(iseries)%loaded) then
    STARK_MMM = STARK(N, M, J)
    return
  end if

  S => STEHLE_DATA(iseries)

  ! Check if this transition is in the table range
  if (M < S%n_upper_min .or. M > S%n_upper_max) then
    STARK_MMM = STARK(N, M, J)
    return
  end if

  itrans = M - S%n_upper_min + 1

  ! ---------------------------------------------------------------
  ! Compute detuning in Δα units
  ! ---------------------------------------------------------------
  xn = dble(N)
  xm = dble(M)
  lambda0 = RYD_ANG * (xn * xm)**2 / ((xm - xn) * (xm + xn))
  freq_nm = CLIGHT_ANG / lambda0
  del_freq = abs(FREQ - freq_nm)

  F0 = 1.25D-9 * XNE(J)**(2.0D0/3.0D0)
  if (F0 <= 0.0D0) then
    STARK_MMM = STARK(N, M, J)
    return
  end if

  ! Δα = Δλ / F0, and Δλ = λ₀²/c × Δν (for small Δλ)
  dalpha = lambda0**2 / CLIGHT_ANG * del_freq / F0

  ! ---------------------------------------------------------------
  ! Check bounds and find bracketing indices
  ! ---------------------------------------------------------------
  nd = S%n_dens
  nt = NSTARK_TEMPS

  ! Density bounds.
  ! Below the tabulated grid: lines are well-defined but the table doesn't
  ! reach this density; fall back to the analytic profile.
  ! Above the tabulated grid: we are past the Inglis-Teller limit for every
  ! transition in this series, so the upper level has dissolved into the
  ! continuum and there is no bound-bound opacity to add. The dissolved
  ! oscillator strength is accounted for on the b-f side via the
  ! Hummer-Mihalas occupation probability formalism in HOP.
  log_ne = LOG10(XNE(J))
  if (XNE(J) < S%density_grid(1)) then
    STARK_MMM = STARK(N, M, J)
    return
  end if
  if (XNE(J) > S%density_grid(nd)) then
    STARK_MMM = 0.0D0
    return
  end if

  ! Per-transition Inglis-Teller limit: the upper level m has dissolved
  ! at this density for this specific transition. Return zero opacity;
  ! the f-strength has moved to the pseudo-continuum (handled in HOP).
  if (XNE(J) > S%density_grid(S%max_dens_idx(itrans))) then
    STARK_MMM = 0.0D0
    return
  end if

  ! Temperature bounds
  if (T(J) < S%temp_grid(1) * 0.5D0 .or. &
      T(J) > S%temp_grid(nt) * 2.0D0) then
    STARK_MMM = STARK(N, M, J)
    return
  end if

  ! ---------------------------------------------------------------
  ! Find density bracket
  ! ---------------------------------------------------------------
  id1 = 1
  do id1 = 1, nd - 1
    if (XNE(J) <= S%density_grid(id1 + 1)) exit
  end do
  id1 = max(1, min(nd - 1, id1))
  id2 = id1 + 1
  frac_d = (log_ne - LOG10(S%density_grid(id1))) / &
           (LOG10(S%density_grid(id2)) - LOG10(S%density_grid(id1)))
  frac_d = max(0.0D0, min(1.0D0, frac_d))

  ! ---------------------------------------------------------------
  ! Find temperature bracket
  ! ---------------------------------------------------------------
  it1 = 1
  do it1 = 1, nt - 1
    if (T(J) <= S%temp_grid(it1 + 1)) exit
  end do
  it1 = max(1, min(nt - 1, it1))
  it2 = it1 + 1
  frac_t = (T(J) - S%temp_grid(it1)) / &
           (S%temp_grid(it2) - S%temp_grid(it1))
  frac_t = max(0.0D0, min(1.0D0, frac_t))

  ! ---------------------------------------------------------------
  ! Find Δα bracket (in log space)
  ! ---------------------------------------------------------------
  if (dalpha <= 0.0D0) then
    ! At line centre: use first grid point value
    log_dalpha = S%log_dalpha_grid(1)
    ia1 = 1
    ia2 = 1
    frac_a = 0.0D0
  else
    log_dalpha = LOG10(dalpha)
    if (log_dalpha <= S%log_dalpha_grid(1)) then
      ia1 = 1
      ia2 = 1
      frac_a = 0.0D0
    else if (log_dalpha >= S%log_dalpha_grid(NSTARK_DALPHA)) then
      ! Beyond table: use asymptotic wing K_alpha / Δα^2.5
      f_nm = hydrogen_f_value(N, M)
      I_dalpha = S%k_alpha(itrans) / dalpha**2.5D0
      STARK_MMM = A0 * f_nm * I_dalpha * lambda0**2 / (CLIGHT_ANG * F0)
      return
    else
      ia1 = 1
      do ia1 = 1, NSTARK_DALPHA - 1
        if (log_dalpha <= S%log_dalpha_grid(ia1 + 1)) exit
      end do
      ia1 = max(1, min(NSTARK_DALPHA - 1, ia1))
      ia2 = ia1 + 1
      frac_a = (log_dalpha - S%log_dalpha_grid(ia1)) / &
               (S%log_dalpha_grid(ia2) - S%log_dalpha_grid(ia1))
      frac_a = max(0.0D0, min(1.0D0, frac_a))
    end if
  end if

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
  if (v000 <= 0.0D0 .or. v001 <= 0.0D0 .or. &
      v010 <= 0.0D0 .or. v011 <= 0.0D0 .or. &
      v100 <= 0.0D0 .or. v101 <= 0.0D0 .or. &
      v110 <= 0.0D0 .or. v111 <= 0.0D0) then
    ! Linear interpolation fallback if any vertex is zero
    v00 = v000 + frac_a * (v001 - v000)
    v01 = v010 + frac_a * (v011 - v010)
    v10 = v100 + frac_a * (v101 - v100)
    v11 = v110 + frac_a * (v111 - v110)
    v0 = v00 + frac_t * (v01 - v00)
    v1 = v10 + frac_t * (v11 - v10)
    I_dalpha = v0 + frac_d * (v1 - v0)
  else
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
  end if

  I_dalpha = max(I_dalpha, 0.0D0)

  ! ---------------------------------------------------------------
  ! Convert I(Δα) to cross-section
  !   σ = (πe²/m_ec) × f_nm × I(Δα) × λ₀² / (c × F₀)
  ! ---------------------------------------------------------------
  f_nm = hydrogen_f_value(N, M)
  STARK_MMM = A0 * f_nm * I_dalpha * lambda0**2 / (CLIGHT_ANG * F0)

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

  implicit none

  ! Named constants (derived from mod_constants)
  real*8, parameter :: CEN_PREFAC4 = 0.026538D0 / SQRTPI / CLIGHT_NM
  real*8, parameter :: FOURPI_C_INV = 1.0D0 / (FOURPI * CLIGHT_NM)
  real*8, parameter :: LORENTZ_PREFAC = INVSQRTPI
  real*8, parameter :: ADAMP_THRESH = 0.20D0
  integer, parameter :: MAX_WING = 100

  ! Local variables
  real*8  :: ELO, CGF, GAMMAR, GAMMAS, GAMMAW
  real*8  :: ADAMP, CENTER, CV, DOPWAVE, VVOIGT, WLVAC4
  real*8  :: TXNXN(kw)
  real*8  :: RATIOLG, START, STOP
  integer*4 :: LINEREC(4)
  integer*4 :: IFJ(kw+1)
  integer :: LINE, J, K, NU, IW, I, IV, NUCONT, IWLOLD, IFLINE, LINEUSED

  if (IDEBUG == 1) write(6, '(A)') ' RUNNING LINOP1'

  IFJ(1) = 0
  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)

  ! Initialize line opacity array
  do NU = 1, NUMNU
    do J = 1, NRHOX
      XLINES(J, NU) = 0.0
    end do
  end do

  ! Precompute van der Waals broadening proxy
  do J = 1, NRHOX
    TXNXN(J) = (XNF(J,1) + 0.42D0*XNF(J,3) + 0.85D0*XNF(J,841)) &
               * (T(J) / 10000.D0)**0.3D0
  end do

  NUCONT = 1
  NU = 1
  IWLOLD = 0
  START = WAVESET(NULO) - 1.0D0
  STOP  = WAVESET(NUHI) + 1.0D0
  LINEUSED = 0

  !---------------------------------------------------------------------
  ! Main line loop (reads from in-memory LINEDATA array)
  !---------------------------------------------------------------------
  do LINE = 1, NLINES_STORED
    LINEREC = LINEDATA(:, LINE)
    call UNPACK_LINEDATA(LINEREC)
    if (mod(LINE, 100000) == 1 .and. ITER == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW

    ! Check wavelength ordering
    ! This will occur e.g., when multiple line lists (lowlines+diatomics)
    ! are concatenated together.  This is not bad, it just slows down the 
    ! subsequent search.  If it happens a lot, it can be a performance hit
    if (IWL < IWLOLD) then
       !write(6, *) IWL, IWLOLD
       !write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
       NUCONT = 1
       NU = 1
    end if

    ! Advance continuous opacity bin
    do while (IWL >= IWAVETAB(NUCONT))
      NUCONT = NUCONT + 1
      if (NUCONT > 344) then
        write(6, '(A)') ' WARNING: NUCONT > 344, clamping'
        NUCONT = 344
        exit
      end if
    end do

    ! Bounds check on species index
    NELION = abs(IELION / 10)
    if (NELION < 1 .or. NELION > mion) then
      if (IDEBUG == 1) write(6, '(A,I6,A,I10)') &
        '  LINOP1: NELION=', NELION, ' OOB, LINE=', LINE
      IWLOLD = IWL
      cycle
    end if

    ! Convert wavelength
    WLVAC = exp(IWL * RATIOLG)
    WLVAC4 = WLVAC
    if (WLVAC < START .or. WLVAC > STOP) then
      IWLOLD = IWL
      cycle
    end if

    ! Advance frequency grid to match line position
    do while (WLVAC >= WAVESET(NU))
      NU = NU + 1
      if (NU >= NUMNU) then
        ! Last point may miss blue line wings
        IWLOLD = IWL
        cycle
      end if
    end do

    ! Convert line parameters
    CGF = CEN_PREFAC4 * WLVAC4 * TABLOG(IGFLOG)
    ELO = TABLOG(IELO)
    ADAMP = 0.0D0
    IFLINE = 0

    !-------------------------------------------------------------------
    ! Phase 1: Coarse depth grid (every 8th depth)
    !-------------------------------------------------------------------
    do J = 8, NRHOX, 8
      IFJ(J+1) = 0
      CENTER = CGF * XNFDOP(J, NELION)
      if (CENTER < TABCONT(J, NUCONT)) cycle
      CENTER = CENTER * exp(-ELO * HCKT(J))
      if (CENTER < TABCONT(J, NUCONT)) cycle
      IFJ(J+1) = 1
      IFLINE = 1

      ! Compute damping (once per line)
      if (ADAMP == 0.0D0) then
        GAMMAR = TABLOG(IGR) * WLVAC4 * FOURPI_C_INV
        GAMMAS = TABLOG(IGS) * WLVAC4 * FOURPI_C_INV
        GAMMAW = TABLOG(IGW) * WLVAC4 * FOURPI_C_INV
      end if

      ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
      DOPWAVE = DOPPLE(J, NELION) * WLVAC4

      if (ADAMP > ADAMP_THRESH) then
        ! --- Full Voigt regime ---
        ! Red wing
        do IW = NU, min(NU + MAX_WING, NUMNU)
          CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
        ! Blue wing
        do I = 1, MAX_WING
          IW = NU - I
          if (IW <= 0) exit
          CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
      else
        ! --- Pretabulated Voigt regime ---
        ! Red wing
        do IW = NU, min(NU + MAX_WING, NUMNU)
          VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
          if (VVOIGT > 10.0D0) then
            CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
          else
            IV = int(VVOIGT * 200.0D0 + 1.5D0)
            CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
          end if
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
        ! Blue wing
        do I = 1, MAX_WING
          IW = NU - I
          if (IW <= 0) exit
          VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
          if (VVOIGT > 10.0D0) then
            CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
          else
            IV = int(VVOIGT * 200.0D0 + 1.5D0)
            CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
          end if
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
      end if
    end do  ! coarse depth J

    !-------------------------------------------------------------------
    ! Phase 2: Fine depth grid (intermediate points between coarse)
    !-------------------------------------------------------------------
    do K = 8, NRHOX, 8
      if (IFJ(K-7) + IFJ(K+1) == 0) cycle
      do J = K-7, K-1
        CENTER = CGF * XNFDOP(J, NELION)
        if (CENTER < TABCONT(J, NUCONT)) cycle
        CENTER = CENTER * exp(-ELO * HCKT(J))
        if (CENTER < TABCONT(J, NUCONT)) cycle

        ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
        DOPWAVE = DOPPLE(J, NELION) * WLVAC4

        if (ADAMP > ADAMP_THRESH) then
          ! --- Full Voigt regime ---
          ! Red wing
          do IW = NU, min(NU + MAX_WING, NUMNU)
            CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
          ! Blue wing
          do I = 1, MAX_WING
            IW = NU - I
            if (IW <= 0) exit
            CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
        else
          ! --- Pretabulated Voigt regime ---
          ! Red wing
          do IW = NU, min(NU + MAX_WING, NUMNU)
            VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
            if (VVOIGT > 10.0D0) then
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            else
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            end if
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
          ! Blue wing
          do I = 1, MAX_WING
            IW = NU - I
            if (IW <= 0) exit
            VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
            if (VVOIGT > 10.0D0) then
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            else
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            end if
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
        end if
      end do  ! fine depth J
    end do  ! coarse block K

    if (IFLINE == 1) LINEUSED = LINEUSED + 1

    IWLOLD = IWL
  end do  ! LINE

  return

END SUBROUTINE LINOP1


!=======================================================================
! XCONOP
!=======================================================================

SUBROUTINE XCONOP

  implicit none
  integer :: J

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING XCONOP'
  ! Extra user-defined continuum opacity from tabulated Rosseland mean
  ! Source function = (sigma/pi) * T^4 (Stefan-Boltzmann)
  do J = 1, NRHOX
    AXCONT(J) = ROSSTAB(T(J), P(J), VTURB(J))
    SXCONT(J) = SIGMA_SB / FOURPI * T(J)**4 * 4.0D0
  end do
  return

END SUBROUTINE XCONOP

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

  implicit none

  integer, intent(in) :: IFSCAT, IFSURF

  ! Fixed quadrature grid parameters
  integer, parameter :: NXTAU = 51
  integer, parameter :: MAX_ITER = 50  ! convergence ceiling (both iterations)

  ! H (flux) quadrature weights on the 51-point grid
  real*8, parameter :: CH_J(NXTAU) = (/ &
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
  real*8, parameter :: CK_J(NXTAU) = (/ &
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
  real*8, parameter :: XTAU8(NXTAU) = (/ &
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
  real*8, save :: EXTAU(NXTAU, 20) = 0.0D0

  ! Local variables
  real*8 :: XS(NXTAU), XSBAR(NXTAU), XALPHA(NXTAU), DIAG(NXTAU)
  real*8 :: XH(NXTAU), XJS(NXTAU)
  real*8 :: XSBAR8(NXTAU), XALPHA8(NXTAU), XS8(NXTAU)
  real*8 :: XH8(NXTAU), XJS8(NXTAU)
  real*8 :: A(kw), B(kw), C(kw), SNUBAR(kw), CTWO(kw), B2CT(kw), B2CT1(kw)
  real*8 :: DELXS, ERRORX, XK, ERROR, SNEW, EXNEW
  real*8 :: TANGLE, D, DDDDD, OLD, SUM_VAL
  integer :: J, JJ, K, KK, L, M, MAXJ, MAXJ1, MU, N1, NM1, NMJ, MDUMMY
  integer :: IFERR, IFNEG
  integer :: output_mode

  if (IDEBUG == 1) write(6, '(A)') ' RUNNING JOSH'

  !---------------------------------------------------------------------
  ! Compute total opacity, scattering fraction, and thermal source
  !---------------------------------------------------------------------
  do J = 1, NRHOX
    ABTOT(J) = ACONT(J) + ALINE(J) + SIGMAC(J) + SIGMAL(J)
    ALPHA(J) = (SIGMAC(J) + SIGMAL(J)) / ABTOT(J)
    SNUBAR(J) = (ACONT(J)*SCONT(J) + ALINE(J)*SLINE(J)) &
               / (ACONT(J) + ALINE(J))
  end do
  call INTEG(RHOX, ABTOT, TAUNU, NRHOX, ABTOT(1)*RHOX(1))
  MAXJ = 0

  !---------------------------------------------------------------------
  ! Solve radiative transfer. Output mode depends on IFSURF:
  !   0 → full J, H, K solution
  !   1 → surface flux H(1) only
  !   2 → surface intensity SURFI(mu) via piecewise parabolic
  ! The solver block exits early when the output mode is determined.
  !---------------------------------------------------------------------
  output_mode = 0

  solver: do   ! single-pass block for structured exit

    !-------------------------------------------------------------------
    ! No-scattering path: S_nu = S_bar
    !-------------------------------------------------------------------
    if (IFSCAT == 0) then
      do J = 1, NRHOX
        SNU(J) = SNUBAR(J)
      end do
      if (IFSURF == 2) then
        output_mode = 70
        exit solver
      end if
      MAXJ = MAP1(TAUNU, SNU, NRHOX, XTAU8, XS8, NXTAU)
      do L = 1, NXTAU
        XS(L) = XS8(L)
      end do
      if (IFSURF == 1) then
        output_mode = 60
        exit solver
      end if
      do J = 1, NRHOX
        ALPHA(J) = 0.0D0
      end do
    end if

    !-------------------------------------------------------------------
    ! Scattering solution on 51-point grid (Lambda iteration)
    !-------------------------------------------------------------------
    if (TAUNU(1) > XTAU8(NXTAU)) MAXJ = 1

    if (MAXJ /= 1) then
      MAXJ = MAP1(TAUNU, SNUBAR, NRHOX, XTAU8, XSBAR8, NXTAU)
      MAXJ = MAP1(TAUNU, ALPHA, NRHOX, XTAU8, XALPHA8, NXTAU)

      do L = 1, NXTAU
        ! Clamp in case of bad interpolation
        XALPHA8(L) = max(XALPHA8(L), 0.0D0)
        XSBAR8(L)  = max(XSBAR8(L), 1.0D-38)
        XALPHA(L) = XALPHA8(L)
        XSBAR(L)  = XSBAR8(L)
        ! Extrapolate if tau grid starts above model surface
        if (XTAU8(L) < TAUNU(1)) then
          XSBAR8(L) = max(SNUBAR(1), 1.0D-38)
          XALPHA8(L) = max(ALPHA(1), 0.0D0)
          XSBAR(L) = XSBAR8(L)
          XALPHA(L) = XALPHA8(L)
        end if
        XS(L) = XSBAR(L)
        DIAG(L) = 1.0D0 - XALPHA(L) * COEFJ(L, L)
        if (abs(DIAG(L)) < 1.0D-30) DIAG(L) = sign(1.0D-30, DIAG(L))
        XSBAR(L) = (1.0D0 - XALPHA(L)) * XSBAR(L)
      end do

      ! Lambda iteration: Gauss-Seidel sweeps (max MAX_ITER iterations)
      BLOCK
        real*8  :: WORST_ERR
        integer :: WORST_K
        do L = 1, MAX_ITER
          IFERR = 0
          WORST_ERR = 0.0D0
          WORST_K   = 0
          K = NXTAU + 1
          do KK = 1, NXTAU
            K = K - 1
            DELXS = 0.0D0
            do M = 1, NXTAU
              DELXS = DELXS + COEFJ(K, M) * XS(M)
            end do
            DELXS = (DELXS * XALPHA(K) + XSBAR(K) - XS(K)) / DIAG(K)
            ERRORX = abs(DELXS / XS(K))
            if (ERRORX > 0.00001D0) IFERR = 1
            if (ERRORX > WORST_ERR) then
              WORST_ERR = ERRORX
              WORST_K   = K
            end if
            XS(K) = max(XS(K) + DELXS, 1.0D-37)
          end do
          if (IFERR == 0) exit
        end do
        if (IFERR /= 0) &
          write(6, '(A,I4,A,1PE12.4,A,I3,A,E9.2)') &
            ' JOSH WARNING: Lambda iter not converged in ', MAX_ITER, &
            ' sweeps  wave=', CLIGHT_NM/FREQ, '  tau_pt=', WORST_K, '  err=', WORST_ERR
      END BLOCK

      ! Post-iteration dispatch
      if (IFSURF == 1) then
        output_mode = 60
        exit solver
      end if
      do M = 1, NXTAU
        XS8(M) = XS(M)
      end do
      if (IFSURF == 2) then
        output_mode = 670
        exit solver
      end if
      MDUMMY = MAP1(XTAU8, XS8, NXTAU, TAUNU, SNU, MAXJ)
    end if  ! MAXJ /= 1

    !-------------------------------------------------------------------
    ! Deep atmosphere: variable Eddington factor on physical grid
    !-------------------------------------------------------------------
    if (MAXJ /= NRHOX) then
      MAXJ1 = MAXJ + 1
      if (MAXJ == 1) MAXJ1 = 1
      do J = MAXJ1, NRHOX
        SNU(J) = SNUBAR(J)
      end do
      M = max(MAXJ - 1, 1)
      NM1 = NRHOX - M + 1
      NMJ = NRHOX - MAXJ + 1

      ! Variable Eddington iteration (max MAX_ITER iterations)
      do L = 1, MAX_ITER
        ERROR = 0.0D0
        IFNEG = 0

        ! Safety check: negative SNU → reset to Planck function
        do J = M, NRHOX
          if (SNU(J) <= 0.0D0) then
            IFNEG = 1
            do JJ = M, NRHOX
              SNUBAR(JJ) = BNU(JJ)
              SNU(JJ) = BNU(JJ)
            end do
            exit
          end if
        end do

        call DERIV(TAUNU(M), SNU(M), HNU(M), NM1)

        ! Safety check: negative HNU → reset to Planck function
        do J = M, NRHOX
          if (HNU(J) <= 0.0D0) then
            IFNEG = 1
            do JJ = M, NRHOX
              SNUBAR(JJ) = BNU(JJ)
              SNU(JJ) = BNU(JJ)
            end do
            call DERIV(TAUNU(M), SNU(M), HNU(M), NM1)
            exit
          end if
        end do

        do J = M, NRHOX
          HNU(J) = HNU(J) / 3.0D0
        end do
        call DERIV(TAUNU(MAXJ), HNU(MAXJ), JMINS(MAXJ), NMJ)
        do J = MAXJ1, NRHOX
          if (IFNEG == 1) JMINS(J) = 0.0D0
          JNU(J) = JMINS(J) + SNU(J)
          SNEW = (1.0D0 - ALPHA(J)) * SNUBAR(J) + ALPHA(J) * JNU(J)
          ERROR = abs(SNEW - SNU(J)) / SNEW + ERROR
          SNU(J) = SNEW
        end do
        if (ERROR < 0.00001D0) exit
      end do
      if (ERROR >= 0.00001D0) &
        write(6, '(A,I4,A,1PE10.3,A,0PF10.3)') ' JOSH WARNING: Eddington iteration did not converge in ', &
          MAX_ITER, ' sweeps, err=', ERROR, '  wave=', CLIGHT_NM/FREQ
    end if  ! MAXJ /= NRHOX

    !-------------------------------------------------------------------
    ! Post-solution dispatch
    !-------------------------------------------------------------------
    if (IFSURF == 2) then
      output_mode = 70
      exit solver
    end if
    if (MAXJ == 1) then
      KNU(1) = JNU(1) / 3.0D0
      return
    end if

    ! Compute J, H, K from 51-point operator matrices
    do L = 1, NXTAU
      XJS(L) = -XS(L)
      do M = 1, NXTAU
        XJS(L) = XJS(L) + COEFJ(L, M) * XS(M)
      end do
      XJS8(L) = XJS(L)
      XH(L) = 0.0D0
      do M = 1, NXTAU
        XH(L) = XH(L) + COEFH(L, M) * XS(M)
      end do
      XH8(L) = XH(L)
    end do
    MDUMMY = MAP1(XTAU8, XJS8, NXTAU, TAUNU, JMINS, MAXJ)
    MDUMMY = MAP1(XTAU8, XH8, NXTAU, TAUNU, HNU, MAXJ)
    XK = 0.0D0
    do M = 1, NXTAU
      XK = XK + CK_J(M) * XS(M)
    end do
    KNU(1) = XK
    do J = 1, MAXJ
      SNU(J)  = max(SNU(J), 1.0D-38)
      JNU(J) = max(JMINS(J) + SNU(J), 1.0D-38)
    end do
    return

  end do solver

  !---------------------------------------------------------------------
  ! Output mode dispatch (reached via EXIT solver)
  !---------------------------------------------------------------------
  select case (output_mode)

  case (60)
    ! Surface flux: H(1) = sum of CH_J weights times source function
    XH(1) = 0.0D0
    do M = 1, NXTAU
      XH(1) = XH(1) + CH_J(M) * XS(M)
    end do
    HNU(1) = XH(1)

  case (670)
    ! Surface intensity from 51-point grid (piecewise parabolic)
    call PARCOE(XS8, XTAU8, A, B, C, NXTAU)
    N1 = NXTAU - 1
    do J = 1, NXTAU
      CTWO(J) = C(J) * 2.0D0
      B2CT(J) = B(J) + CTWO(J) * XTAU8(J)
    end do
    do J = 1, N1
      B2CT1(J) = B(J) + CTWO(J) * XTAU8(J+1)
    end do
    ! Compute and cache exp(-tau/mu) if not yet done
    if (EXTAU(1,1) == 0.0D0) then
      do MU = 1, NMU
        do J = 1, NXTAU
          TANGLE = XTAU8(J) / ANGLE(MU)
          if (TANGLE < 300.0D0) EXTAU(J, MU) = exp(-TANGLE)
        end do
      end do
    end if
    do MU = 1, NMU
      SURFI(MU) = 0.0D0
      do J = 1, N1
        if (EXTAU(J, MU) == 0.0D0) exit
        SURFI(MU) = SURFI(MU) &
          + EXTAU(J, MU) * (XS8(J) + (B2CT(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU)) &
          - EXTAU(J+1, MU) * (XS8(J+1) + (B2CT1(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU))
      end do
      SURFI(MU) = SURFI(MU) &
        + EXTAU(NXTAU, MU) * (XS8(NXTAU) + (B2CT(NXTAU) + CTWO(NXTAU)*ANGLE(MU)) * ANGLE(MU))
    end do

  case (70)
    ! Surface intensity from physical grid (piecewise parabolic)
    call PARCOE(SNU, TAUNU, A, B, C, NRHOX)
    N1 = NRHOX - 1
    do J = 1, NRHOX
      CTWO(J) = C(J) * 2.0D0
      B2CT(J) = B(J) + CTWO(J) * TAUNU(J)
    end do
    do J = 1, N1
      B2CT1(J) = B(J) + CTWO(J) * TAUNU(J+1)
    end do
    do MU = 1, NMU
      OLD = 1.0D0
      SUM_VAL = 0.0D0
      BLOCK
        logical :: tangle_done
        tangle_done = .false.
        do J = 1, N1
          TANGLE = TAUNU(J+1) / ANGLE(MU)
          EXNEW = exp(-TANGLE)
          D = TANGLE - TAUNU(J) / ANGLE(MU)
          if (D > 0.03D0) then
            SUM_VAL = SUM_VAL &
              + OLD * (SNU(J) + (B2CT(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU)) &
              - EXNEW * (SNU(J+1) + (B2CT1(J) + CTWO(J)*ANGLE(MU)) * ANGLE(MU))
            if (TANGLE >= 300.0D0) then
              SURFI(MU) = SUM_VAL
              tangle_done = .true.
              exit
            end if
          else
            ! Small optical depth increment: Taylor expansion for stability
            DDDDD = 1.0D0
            if (D > 0.001D0) &
              DDDDD = ((((D/9.0D0 + 1.0D0)*D/8.0D0 + 1.0D0)*D/7.0D0 + 1.0D0) &
                        *D/6.0D0 + 1.0D0)*D/5.0D0 + 1.0D0
            SUM_VAL = SUM_VAL + EXNEW * (SNU(J) + (SNU(J) + B2CT(J)*ANGLE(MU) &
              + (SNU(J) + (B2CT(J) + CTWO(J)*ANGLE(MU))*ANGLE(MU)) &
              * (DDDDD*D/4.0D0 + 1.0D0)*D/3.0D0)*D/2.0D0) * D
          end if
          OLD = EXNEW
        end do
        if (.not. tangle_done) then
          SURFI(MU) = SUM_VAL &
            + OLD * (SNU(NRHOX) + (B2CT(NRHOX) + CTWO(NRHOX)*ANGLE(MU)) * ANGLE(MU))
        end if
      END BLOCK
    end do

  end select
  return

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

  implicit none

  real*4  :: CJ(2601)
  integer :: I, IOS
  logical, save :: INITIALIZED = .false.

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING BLOCKJ'

  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'blockj.dat', status='OLD', &
         action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'BLOCKJ: ERROR opening ', trim(DATADIR)//'blockj.dat'
      stop 'BLOCKJ: cannot read J-coefficient data'
    end if
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')
    read(89, *) (CJ(I), I=1, 2601)
    close(89)
    COEFJ = reshape(dble(CJ), (/ 51, 51 /))
    INITIALIZED = .true.
  end if

  return

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

  implicit none

  real*4  :: CH(2601)
  integer :: I, IOS
  logical, save :: INITIALIZED = .false.

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING BLOCKH'

  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'blockh.dat', status='OLD', &
         action='READ', iostat=IOS)
    if (IOS /= 0) then
      write(6,*) 'BLOCKH: ERROR opening ', trim(DATADIR)//'blockh.dat'
      stop 'BLOCKH: cannot read H-coefficient data'
    end if
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')
    read(89, *) (CH(I), I=1, 2601)
    close(89)
    COEFH = reshape(dble(CH), (/ 51, 51 /))
    INITIALIZED = .true.
  end if

  return

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

  implicit none

  integer :: IZ, ION, IOST
  logical, save :: INITIALIZED = .false.

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING IONPOTS'

  ! Read ionization potentials on first call
  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'ionpots.dat', status='OLD', &
         action='READ', iostat=IOST)
    if (IOST /= 0) then
      write(6,*) 'IONPOTS: ERROR opening ', trim(DATADIR)//'ionpots.dat'
      stop 'IONPOTS: cannot read ionization potential data'
    end if
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')
    read(89, *) (POTION(IZ), IZ=1, 999)
    close(89)
    INITIALIZED = .true.
  end if

  ! Compute cumulative ionization potential sums
  NELION = 0

  ! Light elements (Z=1-30): all ionization stages
  do IZ = 1, 30
    NELION = NELION + 1
    POTIONSUM(NELION) = 0.0D0       ! neutral atom
    do ION = 2, IZ + 1
      NELION = NELION + 1
      POTIONSUM(NELION) = POTION(NELION - 1) + POTIONSUM(NELION - 1)
    end do
  end do

  ! Heavy elements (Z=31-99): 5 stages only (neutral + 4 ions)
  do IZ = 31, 99
    NELION = NELION + 1
    POTIONSUM(NELION) = 0.0D0       ! neutral atom
    do ION = 1, 4
      NELION = NELION + 1
      POTIONSUM(NELION) = POTION(NELION - 1) + POTIONSUM(NELION - 1)
    end do
  end do

  return

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

  implicit none

  real*4,  save :: ISOION(20, 265)
  logical, save :: INITIALIZED = .false.

  integer :: IZ, ION, N, I, IOST

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING ISOTOPES'

  ! Read isotope data on first call
  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'isotopes.dat', status='OLD', &
         action='READ', iostat=IOST)
    if (IOST /= 0) then
      write(6,*) 'ISOTOPES: ERROR opening ', trim(DATADIR)//'isotopes.dat'
      stop 'ISOTOPES: cannot read isotope data'
    end if
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')
    read(89, *) ((ISOION(I, IZ), I=1,20), IZ=1,265)
    close(89)
    INITIALIZED = .true.
  end if

  ! Copy isotope data to all ionization stages
  N = 0

  ! Light elements (Z=1-30): IZ+1 ionization stages
  do IZ = 1, 30
    do ION = 1, IZ + 1
      N = N + 1
      do I = 1, 10
        ISOTOPE(I, 1, N) = ISOION(I, IZ)
        ISOTOPE(I, 2, N) = ISOION(I + 10, IZ)
      end do
    end do
  end do

  ! Heavy elements (Z=31-99): 5 ionization stages
  do IZ = 31, 99
    do ION = 1, 5
      N = N + 1
      do I = 1, 10
        ISOTOPE(I, 1, N) = ISOION(I, IZ)
        ISOTOPE(I, 2, N) = ISOION(I + 10, IZ)
      end do
    end do
  end do

  ! Molecular species (IZ=100-265): 1 state each
  do IZ = 100, 265
    N = N + 1
    do I = 1, 10
      ISOTOPE(I, 1, N) = ISOION(I, IZ)
      ISOTOPE(I, 2, N) = ISOION(I + 10, IZ)
    end do
  end do

  return

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

  implicit none

  integer, parameter :: NWAVE = 343

  ! --- Wavelength grid (nm): 343 points from 9.09 to 400000 nm ---
  real*8, parameter :: WBIG(NWAVE) = (/ &
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
  real*8  :: RATIOLG, FREQ15
  integer :: NU, NU9, J, N

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING KAPCONT'

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
  do NU = 1, NWAVE
    WAVE = WAVETAB(NU)
    IWAVETAB(NU) = int(log(WAVE) / RATIOLG + 0.5D0)
    FREQ = CLIGHT_NM / WAVETAB(NU)
    WAVENO = 1.0D7 / WAVETAB(NU)
    FREQLG = log(FREQ)

    ! Set up thermodynamic quantities at each depth
    FREQ15 = FREQ / 1.0D15
    do J = 1, NRHOX
      ACONT(J) = 1.0D10
      SCONT(J) = 0.0D0
      EHVKT(J) = exp(-FREQ * HKT(J))
      STIM(J) = 1.0D0 - EHVKT(J)
      BNU(J) = BNU_PREFAC * FREQ15**3 * EHVKT(J) / STIM(J)
    end do

    ! Compute continuous opacity
    N = 0
    if (WAVE > WAVESET(NULO)) call KAPP

    ! Store tabulated opacity
    do J = 1, NRHOX
      TABCONT(J, NU) = (ACONT(J) + SIGMAC(J)) * 0.001D0 / STIM(J)
    end do
  end do

  ! Pad last entry
  do J = 1, NRHOX
    TABCONT(J, NWAVE + 1) = TABCONT(J, NWAVE)
  end do
  WAVETAB(NWAVE + 1) = WAVETAB(NWAVE)

  !---------------------------------------------------------------------
  ! Diagnostic printout
  !---------------------------------------------------------------------

  if (IDEBUG == 1) then
     do NU = 1, NWAVE, 10
        NU9 = min(NU + 9, NWAVE)
        write(6, '(5X, 10F12.2)') (WAVETAB(N), N = NU, NU9)
        do J = 1, NRHOX
           write(6, '(I5, 1P10E12.3)') J, (TABCONT(J, N), N = NU, NU9)
        end do
     end do
  endif

  return

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

  implicit none

  integer, parameter :: NMOL = 46
  integer, parameter :: MAX_LINES = 500000000  ! max lines per database

  ! Line-center opacity prefactor: pi*e^2/(m_e*c*sqrt(pi)) in CGS
  real*8, parameter :: CEN_PREFAC = 0.026538D0 / SQRTPI

  ! Molecular species codes (for DIATOMICS database lookup)
  integer, parameter :: MOLCODES(NMOL) = (/ &
    8410, 8411, 8460, 8461, 8470, 8471, 8480, 8481, 8482, 8510, &
    8511, 8512, 8530, 8531, 8532, 8580, 8581, 8582, 8583, 8584, &
    8620, 8621, 8622, 8623, 8640, 8641, 8642, 8643, 8680, 8681, &
    8682, 8690, 8691, 8692, 8693, 8700, 8701, 8702, 8703, 8704, &
    8705, 8890, 8891, 8892, 8896, 8960 /)

  ! Isotope gf corrections for diatomic molecules
  ! ISOX(IMOL): additive correction to IGFLOG for each isotopologue
  integer, parameter :: ISOX(NMOL) = (/ &
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
  integer, parameter :: TIO_ISOCORR(5) = (/ -1101, -1138, -131, -1259, -1272 /)

  ! H2O isotope gf corrections (ISO determined from signs of IELO/IGFLOG)
  integer, parameter :: H2O_ISOCORR(4) = (/ -1, -3398, -2690, -5000 /)

  ! --- Local variables ---
  real*8  :: XNFDOPMAX(mion, 344)
  real*8  :: CENRATIO, RATIOLG, GR, tablog8
  integer*4 :: LINEREC(4)
  integer :: NU, J, K, LINE
  integer :: N12, N122, N22, N32, N42, N52, N62, N18
  integer :: MOLCODE, MOLCODEOLD, KGFLOG, ISO, IMOL
  integer :: LINEDATA_CAP, IOS

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING SELECTLINES'

  ! Allocate in-memory line storage (replaces fort.12)
  LINEDATA_CAP = 150000000   ! 150M lines initial (~2.4 GB)
  if (allocated(LINEDATA)) deallocate(LINEDATA)
  allocate(LINEDATA(4, LINEDATA_CAP), stat=IOS)
  if (IOS /= 0) then
    write(6, '(A,F6.2,A)') ' SELECTLINES: cannot allocate ', &
      LINEDATA_CAP * 16.0D0 / 1.0D9, ' GB for LINEDATA'
    stop 'SELECTLINES: insufficient memory'
  end if
  NLINES_STORED = 0

  !---------------------------------------------------------------------
  ! Compute maximum XNFDOP/TABCONT ratio over depth for each (ion, nu)
  !---------------------------------------------------------------------
  XNFDOPMAX = 0.0D0
  do NU = 1, 344
    do NELION = 1, MION
      do J = 1, NRHOX
        XNFDOPMAX(NELION, NU) = max(XNFDOPMAX(NELION, NU), &
          XNFDOP(J, NELION) / TABCONT(J, NU))
      end do
    end do
  end do

  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)
  N12 = 0; N122 = 0; N22 = 0; N32 = 0; N42 = 0; N52 = 0; N62 = 0
  NU = 1

  !=====================================================================
  ! (1) LOWLINES predicted (unit 11)
  !=====================================================================
  open(unit=11, file=trim(DATADIR)//'lowlines_pl.bin', &
       status='OLD', form='UNFORMATTED', action='READ', &
       access='STREAM', err=669)

  do LINE = 1, MAX_LINES
    read(11, end=581) LINEREC
    call UNPACK_LINEDATA(LINEREC)
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW

    ! Advance wavelength bin to match line position
    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    end do

    ! Apply selection filters
    NELION = abs(IELION / 10)
    if (NELION < 1 .or. NELION > mion) then
      if (IDEBUG == 1) write(6, '(A,I6,A,I10)') &
        '  SELECTLINES: NELION=', NELION, ' OOB, LINE=', LINE
      cycle
    end if
    if (XNFDOPMAX(NELION, NU) <= 1.0D-37) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    tablog8 = TABLOG(IELO)
    if (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) < 1.0D0) cycle

    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N12 = N12 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading LOWLINES predicted'
  stop 'SELECTLINES: increase MAX_LINES'

  581 write(6, '(I12,A)') N12, ' LINES FROM LOWLINES'
  close(unit=11)

  !=====================================================================
  ! (2) LOWLINES observed (unit 111)
  !=====================================================================
  open(unit=111, file=trim(DATADIR)//'lowlines_obs.bin', &
       status='OLD', form='UNFORMATTED', action='READ', &
       access='STREAM', err=669)

  do LINE = 1, MAX_LINES
    read(111, end=5819) LINEREC
    call UNPACK_LINEDATA(LINEREC)

    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    end do

    NELION = abs(IELION / 10)
    if (XNFDOPMAX(NELION, NU) <= 1.0D-37) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    tablog8 = TABLOG(IELO)
    if (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) < 1.0D0) cycle

    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    WLVAC = exp(IWL * RATIOLG)
    N122 = N122 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading LOWLINES observed'
  stop 'SELECTLINES: increase MAX_LINES'

  5819 write(6, '(I12,A)') N122, ' LINES FROM LOWLINES observed'
  close(unit=111)

  !=====================================================================
  ! (3) HILINES (unit 21)
  !=====================================================================
  669 open(unit=21, file=trim(DATADIR)//'hilines.bin', &
       status='OLD', form='UNFORMATTED', action='READ', &
           access='STREAM', err=769)
  NU = 1

  do LINE = 1, MAX_LINES
    read(21, end=681) LINEREC
    call UNPACK_LINEDATA(LINEREC)

    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    end do

    NELION = abs(IELION / 10)
    if (XNFDOPMAX(NELION, NU) == 0.0D0) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    tablog8 = TABLOG(IELO)
    if (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) < 1.0D0) cycle

    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N22 = N22 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading HILINES'
  stop 'SELECTLINES: increase MAX_LINES'

  681 write(6, '(I12,A)') N22, ' LINES FROM HILINES'
  close(unit=21)

  !=====================================================================
  ! (4) DIATOMICS (unit 31)
  !=====================================================================
  769 open(unit=31, file=trim(DATADIR)//'diatomicspacksrt.bin', &
       status='OLD', form='UNFORMATTED', action='READ', err=869)
  NU = 1
  MOLCODEOLD = 0
  IMOL = 0

  do LINE = 1, MAX_LINES
    read(31, end=781) LINEREC
    call UNPACK_LINEDATA(LINEREC)

    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    end do

    MOLCODE = abs(IELION)
    KGFLOG = IGFLOG
    IGS = 1

    ! Look up molecular species code (cache last match)
    if (MOLCODE /= MOLCODEOLD) then
      MOLCODEOLD = MOLCODE
      IMOL = 0
      do K = 1, NMOL
        if (MOLCODE == MOLCODES(K)) then
          IMOL = K
          exit
        end if
      end do
      if (IMOL == 0) then
        write(6, '(9I12)') LINE, LINEREC
        stop 1
      end if
    end if

    ! Apply isotope gf correction
    IGFLOG = max(KGFLOG + ISOX(IMOL), 1)

    NELION = abs(IELION / 10)
    if (XNFDOPMAX(NELION, NU) == 0.0D0) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    tablog8 = TABLOG(IELO)
    if (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) < 1.0D0) cycle

    call PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N32 = N32 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading DIATOMICS'
  stop 'SELECTLINES: increase MAX_LINES'

  781 write(6, '(I12,A)') N32, ' LINES FROM DIATOMICS'
  close(unit=31)

  !=====================================================================
  ! (5) TiO (unit 41)
  !=====================================================================
  869 open(unit=41, file=trim(DATADIR)//'schwenke.bin', &
       status='OLD', form='UNFORMATTED', action='READ', &
           access='STREAM', err=1869)
  NU = 1

  do LINE = 1, MAX_LINES
    read(41, end=881) LINEREC
    call UNPACK_LINEDATA(LINEREC)

    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    end do

    KGFLOG = IGFLOG
    ISO = abs(IELION) - 8949

    ! Isotope gf correction for 46Ti-50Ti
    if (ISO >= 1 .and. ISO <= 5) then
      IGFLOG = max(KGFLOG + TIO_ISOCORR(ISO), 1)
    end if

    NELION = 895
    if (XNFDOPMAX(NELION, NU) == 0.0D0) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    tablog8 = TABLOG(IELO)
    if (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) < 1.0D0) cycle

    IGS = 1
    IGW = 9384
    call PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N42 = N42 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading TIOLINES'
  stop 'SELECTLINES: increase MAX_LINES'

  881 write(6, '(I12,A)') N42, ' LINES FROM TIOLINES'
  close(unit=41)

  !=====================================================================
  ! (6) H2O (unit 51) — special 3-integer record format
  !=====================================================================
  1869 open(unit=51, file=trim(DATADIR)//'h2ofastfix.bin', &
       status='OLD', form='UNFORMATTED', action='READ', &
            access='STREAM', err=2869)
  NU = 1

  do LINE = 1, MAX_LINES
    read(51, end=1881) IWL, IELO, IGFLOG

    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      ! Radiation damping from frequency
      GAMMAR = 2.474D-22 * FREQ**2 * 0.001
      GR = log10(dble(GAMMAR))
      IGR = int(GR * 1000.0D0 + 16384.5D0)
      NU = NU + 1
    end do

    ! Determine isotope from signs of IELO and IGFLOG
    if (IELO > 0 .and. IGFLOG > 0) then
      ISO = 1      ! 1H1H16O
    else if (IELO > 0) then
      ISO = 2      ! 1H1H17O
    else if (IGFLOG > 0) then
      ISO = 3      ! 1H1H18O
    else
      ISO = 4      ! 1H2H16O
    end if

    IELION = -(9399 + ISO)
    ELO = dble(abs(IELO))
    KGFLOG = abs(IGFLOG)
    IGFLOG = max(KGFLOG + H2O_ISOCORR(ISO), 1)

    NELION = 940
    if (XNFDOPMAX(NELION, NU) == 0.0D0) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    if (CENRATIO * exp(-ELO * HCKT(NRHOX)) < 1.0D0) cycle

    IELO = int(log10(ELO) * 1000.0D0 + 16384.5D0)
    IGS = 1
    IGW = 9384
    call PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    if (mod(LINE, 100000) == 1 .and. IDEBUG == 1) &
      write(6, '(8I15)') LINE, IWL, IELION, IELO, IGFLOG, IGR, IGS, IGW
    N52 = N52 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading H2OFAST'
  stop 'SELECTLINES: increase MAX_LINES'

  1881 write(6, '(I12,A)') N52, ' LINES FROM H2OFAST'
  close(unit=51)

  !=====================================================================
  ! (7) H3+ (unit 61)
  !=====================================================================
  2869 open(unit=61, file=trim(DATADIR)//'h3plus.dat', &
       status='OLD', form='UNFORMATTED', action='READ', &
            access='STREAM', err=1882)
  NU = 1

  do LINE = 1, MAX_LINES
    read(61, end=2881) LINEREC
    call UNPACK_LINEDATA(LINEREC)

    do while (IWL >= IWAVETAB(NU))
      FREQ = CLIGHT_NM / WAVETAB(NU)
      ! (FREQ removed — using FREQ directly)
      NU = NU + 1
    end do

    KGFLOG = IGFLOG
    IGFLOG = max(KGFLOG - 1272, 1)

    NELION = 895
    if (XNFDOPMAX(NELION, NU) == 0.0D0) cycle
    CENRATIO = CEN_PREFAC * TABLOG(IGFLOG) * XNFDOPMAX(NELION, NU) / FREQ
    if (CENRATIO < 1.0D0) cycle
    tablog8 = TABLOG(IELO)
    if (CENRATIO * exp(-tablog8 * HCKT(NRHOX)) < 1.0D0) cycle

    IGS = 1
    IGW = 9384
    call PACK_LINEDATA(LINEREC)
    NLINES_STORED = NLINES_STORED + 1
    if (NLINES_STORED > LINEDATA_CAP) then
      write(6, '(A,I12)') ' SELECTLINES: LINEDATA overflow at ', NLINES_STORED
      stop 'SELECTLINES: increase LINEDATA_CAP'
    end if
    LINEDATA(:, NLINES_STORED) = LINEREC
    N62 = N62 + 1
  end do
  write(6, '(A,I12,A)') ' FATAL: MAX_LINES (', MAX_LINES, ') exhausted reading H3PLUS'
  stop 'SELECTLINES: increase MAX_LINES'

  2881 write(6, '(I10,A)') N62, ' LINES FROM H3PLUS'
  close(unit=61)

  !=====================================================================
  ! Summary
  !=====================================================================
  1882 N18 = N12 + N122 + N22 + N32 + N42 + N52 + N62
  write(6, '(I12,A)') N18, ' LINES TOTAL'
  write(6, '(I12,A,F6.2,A)') NLINES_STORED, ' LINES STORED IN MEMORY (', &
       NLINES_STORED * 16.0D0 / 1.0D9, ' GB)'
  if (NLINES_STORED == 0) then
    stop 'SELECTLINES ERROR: no lines found'
  end if

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

  implicit none

  integer, intent(in)  :: NOLD, NNEW
  real*8,  intent(in)  :: XOLD(NOLD), FOLD(NOLD)
  real*8,  intent(in)  :: XNEW(NNEW)
  real*8,  intent(out) :: FNEW(NNEW)
  integer :: MAP4

  ! --- Local variables ---
  real*8  :: A, B, C, D
  real*8  :: ABAC, BBAC, CBAC
  real*8  :: AFOR, BFOR, CFOR
  real*8  :: WT
  integer :: K, L, LL, L1, L2

  L  = 2
  LL = 0
  A  = 0.0D0
  B  = 0.0D0
  C  = 0.0D0
  CBAC = 0.0D0; BBAC = 0.0D0; ABAC = 0.0D0
  CFOR = 0.0D0; BFOR = 0.0D0; AFOR = 0.0D0

  do K = 1, NNEW

    !-------------------------------------------------------------------
    ! Bracket: find L such that XOLD(L-1) <= XNEW(K) < XOLD(L)
    !-------------------------------------------------------------------
    do while (L <= NOLD)
      if (XNEW(K) < XOLD(L)) exit
      L = L + 1
    end do

    ! If same interval as last point, reuse coefficients
    if (L == LL) then
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      cycle
    end if

    !-------------------------------------------------------------------
    ! Past end or at start: linear interpolation/extrapolation
    !-------------------------------------------------------------------
    if (L > NOLD .or. L == 2) then
      L = min(NOLD, L)
      C = 0.0D0
      B = (FOLD(L) - FOLD(L-1)) / (XOLD(L) - XOLD(L-1))
      A = FOLD(L) - XOLD(L) * B
      LL = L
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      cycle
    end if

    !-------------------------------------------------------------------
    ! Interior point: construct backward and forward quadratics
    !-------------------------------------------------------------------
    L1 = L - 1

    ! Backward quadratic through (L-2, L-1, L)
    if (L <= LL + 1 .and. L /= 3) then
      ! Adjacent interval: shift (backward = previous forward)
      CBAC = CFOR
      BBAC = BFOR
      ABAC = AFOR
    else
      ! Compute from scratch
      L2 = L - 2
      D = (FOLD(L1) - FOLD(L2)) / (XOLD(L1) - XOLD(L2))
      CBAC = FOLD(L) / ((XOLD(L) - XOLD(L1)) * (XOLD(L) - XOLD(L2))) &
           + (FOLD(L2) / (XOLD(L) - XOLD(L2)) &
            - FOLD(L1) / (XOLD(L) - XOLD(L1))) / (XOLD(L1) - XOLD(L2))
      BBAC = D - (XOLD(L1) + XOLD(L2)) * CBAC
      ABAC = FOLD(L2) - XOLD(L2) * D + XOLD(L1) * XOLD(L2) * CBAC
    end if

    ! At last point: use backward quadratic only
    if (L >= NOLD) then
      C = CBAC
      B = BBAC
      A = ABAC
      LL = L
      FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
      cycle
    end if

    ! Forward quadratic through (L-1, L, L+1)
    D = (FOLD(L) - FOLD(L1)) / (XOLD(L) - XOLD(L1))
    CFOR = FOLD(L+1) / ((XOLD(L+1) - XOLD(L)) * (XOLD(L+1) - XOLD(L1))) &
         + (FOLD(L1) / (XOLD(L+1) - XOLD(L1)) &
          - FOLD(L) / (XOLD(L+1) - XOLD(L))) / (XOLD(L) - XOLD(L1))
    BFOR = D - (XOLD(L) + XOLD(L1)) * CFOR
    AFOR = FOLD(L1) - XOLD(L1) * D + XOLD(L) * XOLD(L1) * CFOR

    ! Curvature-weighted blending
    WT = 0.0D0
    if (abs(CFOR) /= 0.0D0) WT = abs(CFOR) / (abs(CFOR) + abs(CBAC))
    A = AFOR + WT * (ABAC - AFOR)
    B = BFOR + WT * (BBAC - BFOR)
    C = CFOR + WT * (CBAC - CFOR)
    LL = L

    FNEW(K) = A + (B + C * XNEW(K)) * XNEW(K)
  end do

  MAP4 = LL - 1
  return

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

  implicit none
  real*8,  intent(in) :: VSTEPS
  integer, intent(in) :: N

  ! Pretabulate Voigt function components for fast line profile evaluation
  ! 100 steps per Doppler width gives ~2% accuracy
  integer, parameter :: NTAB = 81
  real*8, parameter :: TABVI(81) = (/ &
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
  real*8, parameter :: TABH1(81) = (/ &
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

  integer :: I, IDUM
  real*8  :: VV

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING TABVOIGT'

  do I = 1, N
    H0TAB(I) = dble(I - 1) / VSTEPS
  end do
  IDUM = MAP4(TABVI, TABH1, NTAB, H0TAB, H1TAB, N)
  do I = 1, N
    VV = (dble(I - 1) / VSTEPS)**2
    H0TAB(I) = exp(-VV)
    H2TAB(I) = H0TAB(I) - (VV + VV) * H0TAB(I)
  end do
  return

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

  implicit none

  ! Named constants
  real*8, parameter :: LORENTZ_PREFAC = 0.5642D0  ! 1/sqrt(pi)
  real*8, parameter :: ADAMP_THRESH = 0.20D0
  integer, parameter :: MAX_WING = 2000

  ! Ionization edges (cm^-1): CONTX(edge_index, species)
  ! 25 edges × 16 species. Column 1 = H, 4 = He I, 6 = C I, etc.
  real*8, parameter :: CONTX(25,16) = reshape( (/ &
    109678.764D0, 27419.659D0, 12186.462D0, 6854.871D0, 4387.113D0, 3046.604D0, 2238.32D0, 1713.711D0, &
    1354.044D0, 1096.776D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
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
  real*8  :: ELO, CGF, GAMMAR, GAMMAS, GAMMAW
  real*8  :: ADAMP, CENTER, CV, DOPWAVE, VVOIGT, WLVAC4, GF, G
  real*8  :: XSECTG, CON, FRELIN, EPSIL, ASHORE, BSHORE
  real*8  :: TXNXN(kw)
  real*8  :: BOLTH(kw, 100), EH(100)
  real*8  :: WCON, WMERGE, WSHIFT, EMERGE(kw), Z, WMAX
  real*8  :: NSTARK, NMERGE, RATIOLG
  ! Temporaries matching the binary record layout of nltelines_obs.bin
  real*4  :: ELO4, GF4, GAMMAR4, GAMMAS4, GAMMAW4
  integer*4 :: IFJ(kw+1)
  integer :: TYPE, NLAST
  integer :: LINE, J, K, N, NU, IW, I, IV, NUCONT, IOS_RD
  integer :: NBLO, NBUP, NCON, NELIONX, LIM

  if (IDEBUG == 1) write(6, '(A)') ' RUNNING XLINOP'
  RATIOLG = log(1.0D0 + 1.0D0 / 2000000.0D0)

  !---------------------------------------------------------------------
  ! Zero XLINES if LINOP1 was not called
  !---------------------------------------------------------------------
  if (IFOP(15) == 0) then
    do NU = 1, NUMNU
      do J = 1, NRHOX
        XLINES(J, NU) = 0.0
      end do
    end do
  end if

  !---------------------------------------------------------------------
  ! Precompute depth-dependent quantities
  !---------------------------------------------------------------------
  do J = 1, NRHOX
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
    do N = 11, 100
      EH(N) = ELIM_HI - RYDBERG_H / dble(N)**2
    end do

    ! Hydrogen Boltzmann factors × number density
    do NBLO = 1, 100
      BOLTH(J, NBLO) = exp(-EH(NBLO) * HCKT(J)) * XNFDOP(J, 1)
    end do

    ! Merge energy: dissolved levels above this are continuum
    EMERGE(J) = 109737.312D0 / NMERGE**2
  end do

  !---------------------------------------------------------------------
  ! Main line loop (unit 19)
  !---------------------------------------------------------------------
  rewind 19
  NUCONT = 1
  NU = 1
  IFJ(1) = 0    ! NOTE: was uninitialized in original code

  do LINE = 1, 500000
    read(19, iostat=IOS_RD) WLVAC, ELO4, GF4, NBLO, NBUP, NELION, TYPE, &
                      NCON, NELIONX, GAMMAR4, GAMMAS4, GAMMAW4, IWL, LIM
    if (IOS_RD /= 0) return   ! end-of-file or read error
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

    if (WLVAC > WAVESET(NUHI)) return
    WLVAC4 = WLVAC

    ! Advance continuous opacity bin
    do while (IWL >= IWAVETAB(NUCONT))
      NUCONT = NUCONT + 1
      if (NUCONT > 344) then
        write(6, '(A)') ' WARNING: NUCONT > 344, clamping'
        NUCONT = 344
        exit
      end if
    end do

    ! Advance frequency grid
    do while (WLVAC >= WAVESET(NU))
      NU = NU + 1
      if (NU >= NUMNU) return
    end do

    !-----------------------------------------------------------------
    ! Dispatch on line type
    !-----------------------------------------------------------------
    select case (TYPE)

    case (2)
      ! CORONAL LINE — skip
      cycle

    case (0, 3)
      ! NORMAL LINE (and PRD treated as normal)
       if (mod(LINE, 1000) == 1.AND.IDEBUG == 1) &
            write(6, '(2I10,2F12.6,1P5E12.3,I10)') &
            LINE, NU, WLVAC, WAVESET(NU), CGF, ELO, GAMMAR, GAMMAS, GAMMAW, NELION

      ! --- Coarse depth grid (every 8th depth) ---
      do J = 8, NRHOX, 8
        IFJ(J+1) = 0
        CENTER = CGF * XNFDOP(J, NELION)
        if (CENTER < TABCONT(J, NUCONT)) cycle
        CENTER = CENTER * exp(-ELO * HCKT(J))
        if (CENTER < TABCONT(J, NUCONT)) cycle
        IFJ(J+1) = 1

        ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
        DOPWAVE = DOPPLE(J, NELION) * WLVAC4

        ! Blue-wing cutoff at continuum edge
        WCON = 0.0D0
        if (NCON > 10) NCON = 0
        if (NCON > 0) WCON = 1.0D7 / (CONTX(NCON, NELIONX) - EMERGE(J))
        if (WLVAC < WCON) cycle

        if (ADAMP > ADAMP_THRESH) then
          ! Red wing — full Voigt
          do IW = NU, min(NU + MAX_WING, NUMNU)
            CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
          ! Blue wing — full Voigt
          do I = 1, MAX_WING
            IW = NU - I
            if (IW <= 0) exit
            if (WAVESET(IW) < WCON) exit
            CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
        else
          ! Red wing — pretabulated
          do IW = NU, min(NU + MAX_WING, NUMNU)
            VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
            if (VVOIGT > 10.0D0) then
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            else
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            end if
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
          ! Blue wing — pretabulated
          do I = 1, MAX_WING
            IW = NU - I
            if (IW <= 0) exit
            if (WAVESET(IW) < WCON) exit
            VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
            if (VVOIGT > 10.0D0) then
              CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
            else
              IV = int(VVOIGT * 200.0D0 + 1.5D0)
              CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
            end if
            XLINES(J, IW) = XLINES(J, IW) + CV
            if (CV < TABCONT(J, NUCONT)) exit
          end do
        end if
      end do  ! coarse depth J

      ! --- Fine depth grid (intermediate points) ---
      do K = 8, NRHOX, 8
        if (IFJ(K-7) + IFJ(K+1) == 0) cycle
        do J = K-7, K-1
          CENTER = CGF * XNFDOP(J, NELION)
          if (CENTER < TABCONT(J, NUCONT)) cycle
          CENTER = CENTER * exp(-ELO * HCKT(J))
          if (CENTER < TABCONT(J, NUCONT)) cycle

          ADAMP = (GAMMAR + GAMMAS*XNE(J) + GAMMAW*TXNXN(J)) / DOPPLE(J, NELION)
          DOPWAVE = DOPPLE(J, NELION) * WLVAC4

          WCON = 0.0D0
          if (NCON > 10) NCON = 0
          if (NCON > 0) WCON = 1.0D7 / (CONTX(NCON, NELIONX) - EMERGE(J))
          if (WLVAC < WCON) cycle

          if (ADAMP > ADAMP_THRESH) then
            ! Red wing — full Voigt
            do IW = NU, min(NU + MAX_WING, NUMNU)
              CV = CENTER * VOIGT((WAVESET(IW) - WLVAC) / DOPWAVE, ADAMP)
              XLINES(J, IW) = XLINES(J, IW) + CV
              if (CV < TABCONT(J, NUCONT)) exit
            end do
            ! Blue wing — full Voigt
            do I = 1, MAX_WING
              IW = NU - I
              if (IW <= 0) exit
              if (WAVESET(IW) < WCON) exit
              CV = CENTER * VOIGT((WLVAC - WAVESET(IW)) / DOPWAVE, ADAMP)
              XLINES(J, IW) = XLINES(J, IW) + CV
              if (CV < TABCONT(J, NUCONT)) exit
            end do
          else
            ! Red wing — pretabulated
            do IW = NU, min(NU + MAX_WING, NUMNU)
              VVOIGT = (WAVESET(IW) - WLVAC) / DOPWAVE
              if (VVOIGT > 10.0D0) then
                CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
              else
                IV = int(VVOIGT * 200.0D0 + 1.5D0)
                CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
              end if
              XLINES(J, IW) = XLINES(J, IW) + CV
              if (CV < TABCONT(J, NUCONT)) exit
            end do
            ! Blue wing — pretabulated
            do I = 1, MAX_WING
              IW = NU - I
              if (IW <= 0) exit
              if (WAVESET(IW) < WCON) exit
              VVOIGT = (WLVAC - WAVESET(IW)) / DOPWAVE
              if (VVOIGT > 10.0D0) then
                CV = CENTER * LORENTZ_PREFAC * ADAMP / VVOIGT**2
              else
                IV = int(VVOIGT * 200.0D0 + 1.5D0)
                CV = CENTER * ((H2TAB(IV)*ADAMP + H1TAB(IV))*ADAMP + H0TAB(IV))
              end if
              XLINES(J, IW) = XLINES(J, IW) + CV
              if (CV < TABCONT(J, NUCONT)) exit
            end do
          end if
        end do  ! fine depth J
      end do  ! coarse block K

    case (-1)
      ! HYDROGEN LINE — Stark-broadened profile, all depths
      do J = 1, NRHOX
        CENTER = CGF * BOLTH(J, NBLO)
        if (CENTER < TABCONT(J, NUCONT)) cycle
        WCON = 1.0D7 / (CONTX(NCON, 1) - EMERGE(J))
        ! Red wing
        do IW = NU, min(NU + MAX_WING, NUMNU)
          if (WAVESET(IW) < WCON) cycle
          CV = CENTER * HPROF4(NBLO, NBUP, J, WAVESET(IW) - WLVAC)
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
        ! Blue wing
        do I = 1, MAX_WING
          IW = NU - I
          if (IW <= 0) exit
          if (WAVESET(IW) < WCON) exit
          CV = CENTER * HPROF4(NBLO, NBUP, J, WAVESET(IW) - WLVAC)
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
      end do

    case (1)
      ! AUTOIONIZING LINE — Fano profile, all depths
      FRELIN = CLIGHT_NM / WLVAC
      do J = 1, NRHOX
        CENTER = BSHORE * G * XNFP(J, NELION) / RHO(J)
        if (CENTER < TABCONT(J, NUCONT)) cycle
        CENTER = CENTER * exp(-ELO * HCKT(J))
        if (CENTER < TABCONT(J, NUCONT)) cycle
        ! Red wing
        do IW = NU, min(NU + MAX_WING, NUMNU)
          EPSIL = 2.0D0 * (CLIGHT_NM / WAVESET(IW) - FRELIN) / GAMMAR
          CV = CENTER * (ASHORE * EPSIL + BSHORE) / (EPSIL**2 + 1.0D0) / BSHORE
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
        ! Blue wing
        do I = 1, MAX_WING
          IW = NU - I
          if (IW <= 0) exit
          EPSIL = 2.0D0 * (CLIGHT_NM / WAVESET(IW) - FRELIN) / GAMMAR
          CV = CENTER * (ASHORE * EPSIL + BSHORE) / (EPSIL**2 + 1.0D0) / BSHORE
          XLINES(J, IW) = XLINES(J, IW) + CV
          if (CV < TABCONT(J, NUCONT)) exit
        end do
      end do

    case default
      ! MERGED CONTINUUM — flat opacity from line to dissolution limit
      Z = 1.0D0
      if (NELION == 4) Z = 2.0D0
      WSHIFT = 1.0D7 / (1.0D7 / WLVAC - 109737.312D0 * Z**2 / dble(NLAST)**2)
      XSECTG = GF
      IF (IDEBUG == 1) &
           write(6, '(2I10,2F12.6,1P5E12.3,I10)') &
           LINE, NU, WLVAC, WAVESET(NU), CGF, ELO, GAMMAR, GAMMAS, GAMMAW, NELION
      do J = 1, NRHOX
        WMERGE = 1.0D7 / (1.0D7 / WLVAC - EMERGE(J) * Z**2)
        WMAX = max(WMERGE, WSHIFT)
        CON = XSECTG * XNFP(J, NELION) * exp(-ELO * HCKT(J)) / RHO(J)
        do IW = NU, NU + 1000
          if (WMAX < WAVESET(IW)) exit
          XLINES(J, IW) = XLINES(J, IW) + CON
        end do
      end do

    end select
  end do  ! LINE

  return

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

  implicit none
  integer, intent(in) :: N, M
  real*8 :: HFNM

  ! Hydrogen oscillator strength f(n,m) for transition n → m
  ! Uses approximate formula with caching for repeated N and M values
  real*8,  save :: GINF, GCA, FKN, WTC, FNM
  integer, save :: NSTR = 0, MSTR = 0
  real*8  :: XN, XM, XMN, XMN12, FK, WT

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HFNM'

  HFNM = 0.0D0
  if (M <= N) return

  ! Recompute N-dependent quantities if N changed
  if (N /= NSTR) then
    XN = dble(N)
    GINF = 0.2027D0 / XN**0.71D0
    GCA = 0.124D0 / XN
    FKN = XN * 1.9603D0
    WTC = 0.45D0 - 2.4D0 / XN**3 * (XN - 1.0D0)
    NSTR = N
    MSTR = 0   ! force M recomputation
  end if

  ! Recompute M-dependent quantities if M changed
  if (M /= MSTR) then
    XM = dble(M)
    XMN = dble(M - N)
    FK = FKN * (XM / (XMN * (XM + dble(N))))**3
    XMN12 = XMN**1.2D0
    WT = (XMN12 - 1.0D0) / (XMN12 + WTC)
    FNM = FK * (1.0D0 - WT * GINF - (0.222D0 + GCA / XM) * (1.0D0 - WT))
    MSTR = M
  end if

  HFNM = FNM
  return

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

  implicit none
  real*8, intent(in) :: X
  real*8 :: VCSE1F

  ! E_1(x) exponential integral approximation (Abramowitz & Stegun style)
  ! Three regimes: small x (series), medium x (polynomial), large x (asymptotic)
  if (IDEBUG == 1) write(6,'(A)') ' RUNNING VCSE1F'

  if (X <= 0.0D0) then
    VCSE1F = 0.0D0
  else if (X <= 0.01D0) then
    ! Small x: leading terms of series expansion
    VCSE1F = -log(X) - 0.577215D0 + X
  else if (X <= 1.0D0) then
    ! Medium x: polynomial approximation (Abramowitz & Stegun 5.1.53)
    VCSE1F = -log(X) - 0.57721566D0 + X * (0.99999193D0 &
           + X * (-0.24991055D0 + X * (0.05519968D0 &
           + X * (-0.00976004D0 + X * 0.00107857D0))))
  else if (X <= 30.0D0) then
    ! Large x: rational approximation (Abramowitz & Stegun 5.1.56)
    VCSE1F = (X * (X + 2.334733D0) + 0.25062D0) &
           / (X * (X + 3.330657D0) + 1.681534D0) / X * exp(-X)
  else
    VCSE1F = 0.0D0
  end if
  return

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

  implicit none

  real*8,  intent(in) :: B, P
  integer, intent(in) :: N, M
  real*8  :: STARK_PROFILE

  ! --- Grid points ---
  real*8, parameter :: PP(5) = (/ 0.0D0, 0.2D0, 0.4D0, 0.6D0, 0.8D0 /)
  real*8, parameter :: BETAGRID(15) = (/ &
    1.0D0, 1.259D0, 1.585D0, 1.995D0, 2.512D0, 3.162D0, 3.981D0, &
    5.012D0, 6.310D0, 7.943D0, 10.0D0, 12.59D0, 15.85D0, 19.95D0, 25.12D0 /)

  ! --- Tables (read from file on first call) ---
  real*8,  save :: PROPBM(5, 15, 7)   ! correction table
  real*8,  save :: C(5, 7)             ! asymptotic correction coeff
  real*8,  save :: D(5, 7)             ! asymptotic correction coeff
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*8  :: CORR, B2, SB
  real*8  :: WTPP, WTPM, WTBP, WTBM, CBP, CBM
  real*8  :: PR1, PR2, WT, CC, DD
  integer :: INDX, MMN, IM, IP, J, JM, JP
  integer :: I, K
  character(256) :: LINE

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING STARK_PROFILE'

  !---------------------------------------------------------------------
  ! Read tables from file on first call
  !---------------------------------------------------------------------
  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'stark_profile.dat', &
         status='OLD', action='READ')
    do K = 1, 7
      ! Skip comment lines
      do
        read(89, '(A)') LINE
        if (LINE(1:1) /= '#') then
          backspace(89)
          exit
        end if
      end do
      ! Read PROPBM block: 15 rows of 5 values
      do I = 1, 15
        read(89, *) PROPBM(:, I, K)
      end do
      ! Skip C comment
      do
        read(89, '(A)') LINE
        if (LINE(1:1) /= '#') then
          backspace(89)
          exit
        end if
      end do
      read(89, *) C(:, K)
      ! Skip D comment
      do
        read(89, '(A)') LINE
        if (LINE(1:1) /= '#') then
          backspace(89)
          exit
        end if
      end do
      read(89, *) D(:, K)
    end do
    close(89)
    INITIALIZED = .true.
  end if

  !---------------------------------------------------------------------
  ! Select profile index
  !---------------------------------------------------------------------
  INDX = 7                              ! default: generic (Balmer 18)
  MMN = M - N
  if (N <= 3 .and. MMN <= 2) INDX = 2 * (N - 1) + MMN

  CORR = 1.0D0
  B2 = B * B
  SB = sqrt(B)

  !---------------------------------------------------------------------
  ! B > 500: pure asymptotic wing
  !---------------------------------------------------------------------
  if (B > 500.0D0) then
    STARK_PROFILE = (1.5D0 / SB + 27.0D0 / B2) / B2 * CORR
    return
  end if

  ! Debye parameter interpolation weights
  IM = min(int(5.0D0 * P) + 1, 4)
  IP = IM + 1
  WTPP = 5.0D0 * (P - PP(IM))
  WTPM = 1.0D0 - WTPP

  !---------------------------------------------------------------------
  ! 25.12 < B <= 500: asymptotic with C/D correction
  !---------------------------------------------------------------------
  if (B > 25.12D0) then
    CC = C(IP, INDX) * WTPP + C(IM, INDX) * WTPM
    DD = D(IP, INDX) * WTPP + D(IM, INDX) * WTPM
    CORR = 1.0D0 + DD / (CC + B * SB)
    STARK_PROFILE = (1.5D0 / SB + 27.0D0 / B2) / B2 * CORR
    return
  end if

  !---------------------------------------------------------------------
  ! B <= 25.12: bilinear interpolation in correction table
  !---------------------------------------------------------------------
  ! Find beta grid interval
  JM = 1
  do J = 2, 15
    if (B <= BETAGRID(J)) then
      JM = J - 1
      exit
    end if
  end do
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
  if (B <= 10.0D0) PR1 = 8.0D0 / (83.0D0 + (2.0D0 + 0.95D0 * B2) * B)
  if (B >= 8.0D0)  PR2 = (1.5D0 / SB + 27.0D0 / B2) / B2

  STARK_PROFILE = (PR1 * WT + PR2 * (1.0D0 - WT)) * CORR
  return

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

REAL*8 FUNCTION HPROF4(N, M, J, DELW)

  implicit none
  integer, intent(in) :: N, M, J
  real*8,  intent(in) :: DELW

  ! --- Saved depth-dependent arrays (recomputed when ITEMP changes) ---
  real*8,  save :: PP(kw), FO(kw), GCON1(kw), GCON2(kw)
  real*8,  save :: Y1B(kw), Y1S(kw), C1D(kw), C2D(kw)
  real*8,  save :: T3NHE(kw), T3NH2(kw), EXP4492T(kw)
  real*8,  save :: XNE4(kw), XNF4(kw), XNFP4(kw)
  integer, save :: ITEMP1 = 0

  ! --- Saved line-dependent quantities (recomputed when N,M change) ---
  integer, save :: N1 = 0, M1 = 0
  integer, save :: MMN, IFINS
  real*8,  save :: XN, XN2, XM, XM2, XMN2, XM2MN2, GNM
  real*8,  save :: XKNM, FREQNM, DBETA, WAVENM
  real*8,  save :: C1CON, C2CON, RADAMP, RESONT, VDW, STARKC
  real*8,  save :: Y1NUM, Y1WHT
  real*8,  save :: FINEST(14), FINSWT(14)

  ! --- Stark pattern constants ---
  real*8, parameter :: XKNMTB(4,3) = reshape((/ &
    0.0001716D0, 0.009019D0, 0.1001D0, 0.5820D0, &
    0.0005235D0, 0.01772D0,  0.171D0,  0.866D0, &
    0.0008912D0, 0.02507D0,  0.223D0,  1.02D0 /), (/4,3/))
  real*8, parameter :: Y1WTM(2,2) = reshape((/ &
    1.D18, 1.D17, 1.D16, 1.D14 /), (/2,2/))
  ! RYDH (hydrogen Rydberg frequency) now FREQ_RYDH from mod_constants

  ! --- Fine structure for alpha lines (Δn=1) in freq × 10⁻⁷ ---
  real*8, parameter :: STALPH(34) = (/ &
    -730.D0, 370.D0, 188.D0, 515.D0, 327.D0, 619.D0, &
    -772.D0, -473.D0, -369.D0, 120.D0, 256.D0, 162.D0, &
    285.D0, -161.D0, -38.3D0, 6.82D0, -174.D0, -147.D0, &
    -101.D0, -77.5D0, 55.D0, 126.D0, 75.D0, 139.D0, &
    -60.D0, 3.7D0, 27.D0, -69.D0, -42.D0, -18.D0, &
    -5.5D0, -9.1D0, -33.D0, -24.D0 /)
  real*8, parameter :: STWTAL(34) = (/ &
    1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 2.D0, &
    3.D0, 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 4.D0, 6.D0, &
    1.D0, 2.D0, 3.D0, 4.D0, 1.D0, 2.D0, 1.D0, 2.D0, &
    1.D0, 4.D0, 6.D0, 1.D0, 7.D0, 6.D0, 4.D0, 4.D0, &
    4.D0, 5.D0 /)
  integer, parameter :: ISTAL(4) = (/ 1, 3, 10, 21 /)
  integer, parameter :: LNGHAL(4) = (/ 2, 7, 11, 14 /)

  ! --- Fine structure for m=∞ in freq × 10⁻⁷ ---
  real*8, parameter :: STCOMP(5,4) = reshape((/ &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    468.D0, 576.D0, -522.D0, 0.D0, 0.D0, &
    260.D0, 290.D0, -33.D0, -140.D0, 0.0D0, &
    140.D0, 150.D0, 18.D0, -27.D0, -51.D0 /), (/5,4/))
  real*8, parameter :: STCPWT(5,4) = reshape((/ &
    1.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    1.D0, 1.D0, 2.D0, 0.D0, 0.D0, &
    1.D0, 1.D0, 4.D0, 3.D0, 0.D0, &
    1.D0, 1.D0, 4.D0, 6.D0, 4.D0 /), (/5,4/))
  integer, parameter :: LNCOMP(4) = (/ 1, 3, 4, 5 /)

  ! --- Lyman alpha quasi-molecular satellite cutoffs (Allard 1997) ---
  ! H₂⁺ cutoff: Δν = -15000+100*(i-1) cm⁻¹, i=1..111 (to -4000)
  real*8, parameter :: CUTOFFH2PLUS(111) = (/ &
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
  real*8, parameter :: CUTOFFH2(91) = (/ &
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
  real*8, parameter :: ASUMLYMAN(100) = (/ &
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

  real*8, parameter :: ASUM(100) = (/ &
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
  real*8  :: DELstark, WL, FREQ4, DEL, DOP, HFWID
  real*8  :: HWSTK, HWVDW, HWRAD, HWRES, HWLOR, HHW
  real*8  :: HPROFLOR, HPROFRES, HPROFRAD, HPROFVDW
  real*8  :: D, WTY1, Y1SCAL, C1, C2, G1, GNOT, BETA, Y1, Y2, GAM
  real*8  :: PRQS, F, P1, FNS
  real*8  :: CUTOFF, SPACING, CUTFREQ, FREQ15000, FREQ22000
  real*8  :: BETA4000, PRQSP4000, CUTOFF4000
  real*8  :: XNE16, TT4, T4, T43
  integer :: I, K, NWID, IFCORE, IPOS, ICUT

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING HPROF4'

  ! ==================================================================
  ! Section 1: Temperature-dependent depth vectors
  ! ==================================================================
  if (ITEMP /= ITEMP1) then
    ITEMP1 = ITEMP
    do K = 1, NRHOX
      XNE4(K) = XNE(K)
      XNE16 = XNE(K)**0.1666667D0
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
    end do
  end if

  ! ==================================================================
  ! Section 2: Line-dependent constants (cached per N,M)
  ! ==================================================================
  if (N /= N1 .or. M /= M1) then
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
    if (MMN <= 3 .and. N <= 4) then
      XKNM = XKNMTB(N, MMN)
    else
      XKNM = 5.5D-5 / GNM * XMN2 / (1.0D0 + 0.13D0 / dble(MMN))
    end if

    Y1NUM = 320.0D0
    if (M == 2) Y1NUM = 550.0D0
    if (M == 3) Y1NUM = 380.0D0

    Y1WHT = 1.D13
    if (MMN <= 3) Y1WHT = 1.D14
    if (MMN <= 2 .and. N <= 2) Y1WHT = Y1WTM(N, MMN)

    FREQNM = FREQ_RYDH * GNM
    DBETA = CLIGHT_ANG / FREQNM**2 / XKNM
    WAVENM = CLIGHT_ANG / FREQNM
    C1CON = XKNM / WAVENM * GNM * XM2MN2
    C2CON = (XKNM / WAVENM)**2

    ! Radiative damping from ASUM tables (02aug2009)
!     RADAMP=1.389E9/XM**4.53/(1.+5./XM2/XM)
!     IF(N.NE.1)RADAMP=RADAMP+1.389E9/XN**4.53/(1.+5./XN2/XN)
    RADAMP = dble(ASUM(N)) + dble(ASUM(M))
    if (N == 1) RADAMP = dble(ASUMLYMAN(M))
    RADAMP = RADAMP / FOURPI
    RADAMP = RADAMP / FREQNM

    ! Resonance broadening
    RESONT = HFNM(1, M) / XM / (1.0D0 - 1.0D0 / XM2)
    if (N /= 1) RESONT = RESONT + HFNM(1, N) / XN / (1.0D0 - 1.0D0 / XN2)
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
    if (N > 4 .or. M > 10) then
      ! Single unresolved component
      IFINS = 1
      FINEST(1) = 0.0D0
      FINSWT(1) = 1.0D0
    else if (MMN /= 1) then
      ! Non-alpha: use m=∞ structure
      IFINS = LNCOMP(N)
      do I = 1, IFINS
        FINEST(I) = STCOMP(I, N) * 1.D7
        FINSWT(I) = STCPWT(I, N) / XN2
      end do
    else
      ! Alpha lines: exact pattern
      IFINS = LNGHAL(N)
      IPOS = ISTAL(N)
      do I = 1, IFINS
        K = IPOS - 1 + I
        FINEST(I) = STALPH(K) * 1.D7
        FINSWT(I) = STWTAL(K) / XN2 / 3.0D0
      end do
    end if
  end if  ! N,M changed

  ! ==================================================================
  ! Section 3: Profile evaluation at this depth and wavelength
  ! ==================================================================
  DELstark = -10.0D0 * DELW / WAVENM * FREQNM
  WL = WAVENM + DELW * 10.0D0
  FREQ4 = CLIGHT_ANG / WL
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
  if (DOPPLE(J, 1) >= HWSTK .and. DOPPLE(J, 1) >= HWLOR) then
    NWID = 1
  else if (HWLOR >= HWSTK) then
    NWID = 2
  else
    NWID = 3
  end if

  HFWID = FREQNM * max(DOPPLE(J, 1), HWLOR, HWSTK)
  HPROF4 = 0.0D0
  IFCORE = 0
  if (abs(DEL) <= HFWID) IFCORE = 1
  DOP = FREQNM * DOPPLE(J, 1)

  ! ------------------------------------------------------------------
  ! Doppler section: fine-structure resolved Gaussian core
  ! In wing (IFCORE=0): always computed. In core: only if NWID=1.
  ! ------------------------------------------------------------------
  if (IFCORE == 0 .or. NWID == 1) then
    do I = 1, IFINS
      D = abs(FREQ4 - FREQNM - FINEST(I)) / DOP
!     IF(D.LE.7.)HPROF4=HPROF4+EXP(-D*D)/1.77245/DOP*FINSWT(I)
!     SAME NORMALIZATION AS VOIGT FUNCTION
      if (D <= 7.0D0) HPROF4 = HPROF4 + exp(-D * D) * FINSWT(I)
    end do
    if (IFCORE == 1) return
  end if

  ! ------------------------------------------------------------------
  ! Lorentz section: resonance + van der Waals + radiative damping
  ! In wing (IFCORE=0): always computed. In core: only if NWID=2.
  ! ------------------------------------------------------------------
  if (IFCORE == 0 .or. NWID == 2) then
    if (N == 1 .and. M == 2) then
      ! === Lyman alpha special treatment ===
      ! Modify resonance broadening to match at 4000 cm⁻¹
      HWRES = HWRES * 4.0D0
      HWLOR = HWRES + HWVDW + HWRAD
      HHW = FREQNM * HWLOR

      if (FREQ4 > (82259.105D0 - 4000.0D0) * CLIGHT) then
        ! Near center: enhanced resonance Lorentzian
        ! error found by John Lester 31jul2009
        HPROFRES = HWRES * FREQNM / PI / (DEL**2 + HHW**2) &
                 * SQRTPI * DOP
      else
        ! Far red wing: Allard & Koester (1992) H₂ satellite
        CUTOFF = 0.0D0
        if (FREQ4 >= 50000.0D0 * CLIGHT) then
          ! Tabulated at 200 cm⁻¹ spacing
          SPACING = 200.0D0 * CLIGHT
          FREQ22000 = (82259.105D0 - 22000.0D0) * CLIGHT
          if (FREQ4 < FREQ22000) then
            CUTOFF = dble(CUTOFFH2(2) - CUTOFFH2(1)) / SPACING &
                   * (FREQ4 - FREQ22000) + dble(CUTOFFH2(1))
          else
            ICUT = int((FREQ4 - FREQ22000) / SPACING)
            CUTFREQ = dble(ICUT) * SPACING + FREQ22000
            CUTOFF = dble(CUTOFFH2(ICUT+2) - CUTOFFH2(ICUT+1)) / SPACING &
                   * (FREQ4 - CUTFREQ) + dble(CUTOFFH2(ICUT+1))
          end if
          XNFP4(J) = XNFP(J, 1)
          CUTOFF = (10.0D0**(CUTOFF - 14.0D0)) * XNFP4(J) * 2.0D0 &
                 / CLIGHT
        end if
        HPROFRES = CUTOFF * SQRTPI * DOP
      end if

      ! Radiative damping (cut off below Lyman edge to avoid double-counting
      ! with Rayleigh scattering in HRAYOP)
      HPROFRAD = HWRAD * FREQNM / PI / (DEL**2 + HHW**2) &
               * SQRTPI * DOP
!     CORRECTION TO LORENTZ PROFILE   ALLER P.164   NOT USED
!     HPROFRAD=HPROFRAD*4.*FREQ**2/(FREQ**2+FREQNM**2)
      if (FREQ4 <= 2.463D15) HPROFRAD = 0.0D0

      ! Van der Waals (cut off below 60000 cm⁻¹ from line center)
      HPROFVDW = HWVDW * FREQNM / PI / (DEL**2 + HHW**2) &
               * SQRTPI * DOP
      if (FREQ4 < 1.8D15) HPROFVDW = 0.0D0

      HPROFLOR = HPROFRES + HPROFRAD + HPROFVDW
      HPROF4 = HPROF4 + HPROFLOR

    else
      ! === Non-Lyman-alpha: simple combined Lorentzian ===
      HHW = FREQNM * HWLOR
      HPROFLOR = HHW / PI / (DEL**2 + HHW**2) * SQRTPI * DOP
      HPROF4 = HPROF4 + HPROFLOR
    end if

    if (IFCORE == 1) return
  end if

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

  if (.not. (Y2 <= 1.D-4 .and. Y1 <= 1.D-5)) then
!     GAM=G1*(.5*EXP(-MIN(80.,Y1))+VCSE1F(Y1)-.5*VCSE1F(Y2))*
    GAM = G1 * (0.5D0 * exp(-min(80.D0, Y1)) + EXPI(1, Y1) &
        - 0.5D0 * EXPI(1, Y2)) &
        * (1.0D0 - GCON1(J) / (1.0D0 + (90.0D0 * Y1)**3) &
        - GCON2(J) / (1.0D0 + 2000.0D0 * Y1))
    if (GAM <= 1.D-20) GAM = 0.0D0
  end if

  PRQS = STARK_PROFILE(BETA, PP(J), N, M)

  if (M <= 2) then
    ! Lyman alpha Stark: assume quasistatic is half protons, half electrons
    PRQS = PRQS * 0.5D0
    CUTOFF = 0.0D0

    ! H₂⁺ quasi-molecular satellite (Allard 1997)
    if (FREQ4 >= (82259.105D0 - 20000.0D0) * CLIGHT) then
      if (FREQ4 > (82259.105D0 - 4000.0D0) * CLIGHT) then
        ! Near center: use ratio method
        BETA4000 = 4000.0D0 * CLIGHT / FO(J) * DBETA
        PRQSP4000 = STARK_PROFILE(BETA4000, PP(J), N, M) * 0.5D0 / FO(J) * DBETA
        CUTOFF4000 = (10.0D0**(-11.07D0 - 14.0D0)) / CLIGHT * XNF(J, 2)
        HPROF4 = HPROF4 + CUTOFF4000 / PRQSP4000 * PRQS / FO(J) * DBETA &
               * SQRTPI * DOP
      else
        ! Interpolate H₂⁺ cutoff table (100 cm⁻¹ spacing)
        FREQ15000 = (82259.105D0 - 15000.0D0) * CLIGHT
        SPACING = 100.0D0 * CLIGHT
        if (FREQ4 < FREQ15000) then
          CUTOFF = dble(CUTOFFH2PLUS(2) - CUTOFFH2PLUS(1)) / SPACING &
                 * (FREQ4 - FREQ15000) + dble(CUTOFFH2PLUS(1))
        else
          ICUT = int((FREQ4 - FREQ15000) / SPACING)
          CUTFREQ = dble(ICUT) * SPACING + FREQ15000
          CUTOFF = dble(CUTOFFH2PLUS(ICUT+2) - CUTOFFH2PLUS(ICUT+1)) &
                 / SPACING * (FREQ4 - CUTFREQ) + dble(CUTOFFH2PLUS(ICUT+1))
        end if
        CUTOFF = (10.0D0**(CUTOFF - 14.0D0)) / CLIGHT * XNF(J, 2)
        HPROF4 = HPROF4 + CUTOFF * SQRTPI * DOP
      end if
    end if
  end if

  ! Final Stark assembly
  F = 0.0D0
  if (GAM > 0.0D0) F = GAM / PI / (GAM**2 + BETA**2)
  P1 = (0.9D0 * Y1)**2
  FNS = (P1 + 0.03D0 * sqrt(Y1)) / (P1 + 1.0D0)
  ! Same normalization as Voigt function
  HPROF4 = HPROF4 + (PRQS * (1.0D0 + FNS) + F) / FO(J) * DBETA &
         * SQRTPI * DOP
  return

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

REAL*8 FUNCTION VOIGT(V, A)
  
  implicit none

  real*8, intent(in) :: V, A

  integer :: IV
  real*8  :: VV, AA, U, AAU, VVU, UU
  real*8  :: HH1, HH2, HH3, HH4

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING VOIGT'

  ! Table index: 200 steps per Doppler width, offset by 1
  IV = INT(V * 200.d0 + 1.5d0)

  if (A < 0.2d0) then
    !-----------------------------------------------------------------
    ! Small damping: Taylor expansion or Lorentz wing
    !-----------------------------------------------------------------
    if (V <= 10.d0) then
      ! Quadratic Taylor expansion in a using pretabulated coefficients
      VOIGT = (H2TAB(IV)*A + H1TAB(IV))*A + H0TAB(IV)
    else
      ! Far wing: Lorentz profile (a << v)
      VOIGT = 0.5642d0 * A / V**2
    endif

  else if (A <= 1.4d0 .and. A + V <= 3.2d0) then
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

  else
    !-----------------------------------------------------------------
    ! Large damping or far wing: asymptotic Lorentz with correction
    ! H ~ a / (sqrt(2*pi) * U) * [1 + correction terms / U^2]
    !-----------------------------------------------------------------
    AA = A * A
    VV = V * V
    U  = (AA + VV) * 1.4142d0   ! sqrt(2) * (a^2 + v^2)

    VOIGT = A * 0.79788d0 / U   ! a / sqrt(2*pi) / (a^2+v^2)

    if (A <= 100.d0) then
      ! Higher-order correction terms
      AAU = AA / U
      VVU = VV / U
      UU  = U * U
      VOIGT = ((((AAU - 10.d0*VVU)*AAU*3.d0 + 15.d0*VVU*VVU) &
               + 3.d0*VV - AA) / UU + 1.d0) * VOIGT
    endif
  endif

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

  implicit none

  integer, intent(in) :: N
  real*8,  intent(in) :: X
  real*8 :: EXPI

  ! Cached values from previous call
  real*8, save :: X_prev  = -1.D20
  real*8, save :: EX_save = 0.d0   ! exp(-x)
  real*8, save :: E1_save = 0.d0   ! E_1(x)

  ! Local variables
  real*8 :: EX, E1
  integer :: I

  ! Cody & Thacher rational approximation coefficients
  ! Range 0 < x <= 1: E_1(x) = P(x)/Q(x) - ln(x)
  real*8, parameter :: A0 = -44178.5471728217d0
  real*8, parameter :: A1 =  57721.7247139444d0
  real*8, parameter :: A2 =   9938.31388962037d0
  real*8, parameter :: A3 =   1842.11088668000d0
  real*8, parameter :: A4 =    101.093806161906d0
  real*8, parameter :: A5 =      5.03416184097568d0
  real*8, parameter :: B0 =  76537.3323337614d0
  real*8, parameter :: B1 =  32597.1881290275d0
  real*8, parameter :: B2 =   6106.10794245759d0
  real*8, parameter :: B3 =    635.419418378382d0
  real*8, parameter :: B4 =     37.2298352833327d0

  ! Range 1 < x <= 4: E_1(x) = exp(-x) * P(x)/Q(x)
  real*8, parameter :: C0 = 4.65627107975096D-7
  real*8, parameter :: C1 = 0.999979577051595d0
  real*8, parameter :: C2 = 9.04161556946329d0
  real*8, parameter :: C3 = 24.3784088791317d0
  real*8, parameter :: C4 = 23.0192559391333d0
  real*8, parameter :: C5 = 6.90522522784444d0
  real*8, parameter :: C6 = 0.430967839469389d0
  real*8, parameter :: D1 = 10.0411643829054d0
  real*8, parameter :: D2 = 32.4264210695138d0
  real*8, parameter :: D3 = 41.2807841891424d0
  real*8, parameter :: D4 = 20.4494785013794d0
  real*8, parameter :: D5 = 3.31909213593302d0
  real*8, parameter :: D6 = 0.103400130404874d0

  ! Range x > 4: E_1(x) = exp(-x)/x * [1 + P(1/x)/Q(1/x)]
  real*8, parameter :: E0 = -0.999999999998447d0
  real*8, parameter :: E1C = -26.6271060431811d0   ! renamed from E1 to avoid conflict
  real*8, parameter :: E2 = -241.055827097015d0
  real*8, parameter :: E3 = -895.927957772937d0
  real*8, parameter :: E4 = -1298.85688746484d0
  real*8, parameter :: E5 = -545.374158883133d0
  real*8, parameter :: E6 = -5.66575206533869d0
  real*8, parameter :: F1 = 28.6271060422192d0
  real*8, parameter :: F2 = 292.310039388533d0
  real*8, parameter :: F3 = 1332.78537748257d0
  real*8, parameter :: F4 = 2777.61949509163d0
  real*8, parameter :: F5 = 2404.01713225909d0
  real*8, parameter :: F6 = 631.657483280800d0

  !=====================================================================
  ! Compute E_1(x) (use cache if x unchanged)
  !=====================================================================
  if (X /= X_prev) then
    EX = EXP(-X)
    X_prev  = X
    EX_save = EX

    if (X > 4.d0) then
      ! Asymptotic range: E_1 ~ exp(-x)/x * rational(1/x)
      E1 = (EX + EX*(E0 + (E1C + (E2 + (E3 + (E4 + (E5 + E6/X)/X)/X)/X)/X)/X) &
           / (X + F1 + (F2 + (F3 + (F4 + (F5 + F6/X)/X)/X)/X)/X)) / X

    else if (X > 1.d0) then
      ! Mid-range rational approximation
      E1 = EX * (C6 + (C5 + (C4 + (C3 + (C2 + (C1 + C0*X)*X)*X)*X)*X)*X) &
              / (D6 + (D5 + (D4 + (D3 + (D2 + (D1 + X)*X)*X)*X)*X)*X)

    else if (X > 0.d0) then
      ! Small-x: rational approximation minus logarithmic singularity
      E1 = (A0 + (A1 + (A2 + (A3 + (A4 + A5*X)*X)*X)*X)*X) &
         / (B0 + (B1 + (B2 + (B3 + (B4 + X)*X)*X)*X)*X) - LOG(X)

    else
      ! x <= 0: return 0
      E1 = 0.d0
    endif

    E1_save = E1
  endif

  !=====================================================================
  ! Return E_1 or apply recurrence for E_n (n > 1)
  !=====================================================================
  EXPI = E1_save
  if (N == 1) return

  ! Recurrence: E_n(x) = [exp(-x) - x * E_{n-1}(x)] / (n-1)
  do I = 1, N - 1
    EXPI = (EX_save - X * EXPI) / dble(I)
  end do

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
!     ION      — ionization stage (1..10; higher stages return PF=1)
!     TLOG8    — log10(T) (REAL*8)
!     POTLOW8  — potential lowering in cm^-1 (REAL*8)
!     PF       — output: log10(partition function)
!
!   Data source: pfiron.dat (read on first call from DATADIR)
!=========================================================================

SUBROUTINE PFIRON(NELEM, ION, TLOG8, POTLOW8, PF)

  implicit none

  integer, intent(in) :: NELEM, ION
  real*8,  intent(in) :: TLOG8, POTLOW8
  real*8,  intent(out) :: PF

  ! --- Table dimensions ---
  integer, parameter :: NPOT = 7       ! potential lowering bins
  integer, parameter :: NTEMP = 56     ! temperature bins
  integer, parameter :: NION_TAB = 10  ! max ionization stages in table
  integer, parameter :: NELEM_TAB = 9  ! elements: Ca(20)..Ni(28)

  ! --- Persistent table data (read once from file) ---
  real*8,  save :: PFTAB(NPOT, NTEMP, NION_TAB, NELEM_TAB)
  real*8,  save :: POTLO(NPOT)
  real*8,  save :: POTLOLOG(NPOT)
  logical, save :: INITIALIZED = .false.

  ! --- Local variables ---
  real*8  :: TLOG, POTLOW, F, P
  integer :: IT, LOW, I, J, K, L, IELEM

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING PFIRON'

  ! --- Read table on first call ---
  if (.not. INITIALIZED) then
    open(unit=89, file=trim(DATADIR)//'pfiron.dat', status='OLD', action='READ')
    read(89, '(A)') ; read(89, '(A)') ; read(89, '(A)')
    read(89, *) ((((PFTAB(I,J,K,L), I=1,NPOT), J=1,NTEMP), K=1,NION_TAB), L=1,NELEM_TAB)
    close(89)
    POTLO = (/ 500.D0, 1000.D0, 2000.D0, 4000.D0, 8000.D0, 16000.D0, 32000.D0 /)
    POTLOLOG = (/ 2.69897D0, 3.D0, 3.30103D0, 3.60206D0, 3.90309D0, 4.20412D0, 4.50515D0 /)
    INITIALIZED = .true.
  end if

  ! --- Bounds check: ION must be within table range ---
  !     For higher ionization stages (ION > 10), the table has no data.
  !     Return PF = 0 (i.e., log10(partition function) = 0 → PF = 1).
  if (ION < 1 .or. ION > NION_TAB) then
    if (IDEBUG == 1) write(6,'(A,I4,A,I4,A)') &
      '  PFIRON: ION =', ION, ' out of table range (1..', NION_TAB, '), returning PF=0'
    PF = 0.0D0
    return
  end if

  ! --- Bounds check: NELEM must map to valid table index ---
  IELEM = NELEM - 19
  if (IELEM < 1 .or. IELEM > NELEM_TAB) then
    write(6,'(A,I4,A)') ' PFIRON: NELEM =', NELEM, ' out of range (20..28), returning PF=0'
    PF = 0.0D0
    return
  end if

  TLOG = TLOG8
  POTLOW = POTLOW8

  ! --- Determine temperature bin IT and interpolation fraction F ---
  if (TLOG <= 3.7D0) then
    ! Low-T grid: log T = 3.32..3.70, step = 0.02
    IT = int((TLOG - 3.32D0) / 0.02D0) + 2
    IT = max(IT, 2)
    F = (TLOG - (IT - 2) * 0.02D0 - 3.32D0) / 0.02D0
  else if (TLOG <= 4.0D0) then
    ! Mid-T grid: log T = 3.70..4.00, step = 0.03
    IT = int((TLOG - 3.7D0) / 0.03D0) + 21
    F = (TLOG - (IT - 21) * 0.03D0 - 3.7D0) / 0.03D0
  else
    ! High-T grid: log T = 4.00..5.25, step = 0.05
    IT = int((TLOG - 4.0D0) / 0.05D0) + 31
    IT = min(IT, NTEMP)
    F = (TLOG - (IT - 31) * 0.05D0 - 4.0D0) / 0.05D0
  end if

  ! Safety clamp: ensure IT-1 >= 1 and IT <= NTEMP
  IT = max(IT, 2)
  IT = min(IT, NTEMP)
  F = max(0.0D0, min(1.0D0, F))

  ! --- Determine potential lowering bin LOW ---
  LOW = 1
  if (POTLOW >= POTLO(1)) then
    do LOW = 2, NPOT
      if (POTLOW < POTLO(LOW)) exit
    end do
    if (LOW > NPOT) LOW = NPOT
  end if

  if (LOW > 1 .and. POTLOW >= POTLO(1) .and. POTLOW < POTLO(LOW)) then
    ! Two-bin interpolation in potential lowering
    P = (log10(POTLOW) - POTLOLOG(LOW-1)) / 0.30103D0
    PF = P * (F * PFTAB(LOW, IT, ION, IELEM) &
            + (1.0D0 - F) * PFTAB(LOW, IT-1, ION, IELEM)) &
       + (1.0D0 - P) * (F * PFTAB(LOW-1, IT, ION, IELEM) &
            + (1.0D0 - F) * PFTAB(LOW-1, IT-1, ION, IELEM))
  else
    ! Single-bin interpolation (POTLOW below first bin or above last)
    PF = F * PFTAB(LOW, IT, ION, IELEM) &
       + (1.0D0 - F) * PFTAB(LOW, IT-1, ION, IELEM)
  end if

END SUBROUTINE PFIRON

!=======================================================================
! PFGROUND: Ground-state partition functions
!   Returns statistical weight of ground configuration for each ion.
!   NELION = (Z-1)*6 + ionization_stage  (1=neutral, 2=singly ionized, etc.)
!   Covers H-K (Z=1-19), Cu-Ba (Z=29-56). Iron group (Z=20-28) uses PFIRON.
!   Original data from Kurucz, with corrections by J. Laird (2004) and
!   K. Bischof. Bug fix: added missing EXP() in F V (NELION=54).
!
!=======================================================================

FUNCTION PFGROUND(NELION, T)

  implicit none

  ! Arguments
  integer, intent(in)  :: NELION
  real*8,  intent(in)  :: T
  real*8               :: PFGROUND

  ! Local constants
  ! HCK now from mod_constants (updated to CODATA 2018)

  if (IDEBUG == 1) write(6,'(A)') ' RUNNING PFGROUND'

  ! Default: bare ground state
  PFGROUND = 1.0d0

  !---------------------------------------------------------------------
  ! Iron group (Z=20-28, Ca-Ni): partition functions from PFIRON tables.
  ! NELION = 115..168 for these elements.
  !---------------------------------------------------------------------
  if (NELION >= 115 .and. NELION <= 168) return

  select case (NELION)

  ! --- Z= 1  H  ---
  case (  1); PFGROUND = 2.                    ! H  I
  case (  2); PFGROUND = 1.                    ! H  II
  case (  3); PFGROUND = 1.                    ! H  III
  case (  4); PFGROUND = 1.                    ! H  IV
  case (  5); PFGROUND = 1.                    ! H  V
  case (  6); PFGROUND = 1.                    ! H  VI

  ! --- Z= 2  He ---
  case (  7); PFGROUND = 1.                    ! He I
  case (  8); PFGROUND = 2.                    ! He II
  case (  9); PFGROUND = 1.                    ! He III
  case ( 10); PFGROUND = 1.                    ! He IV
  case ( 11); PFGROUND = 1.                    ! He V
  case ( 12); PFGROUND = 1.                    ! He VI

  ! --- Z= 3  Li ---
  case ( 13); PFGROUND = 2.                    ! Li I
  case ( 14); PFGROUND = 1.                    ! Li II
  case ( 15); PFGROUND = 2.                    ! Li III
  case ( 16); PFGROUND = 1.                    ! Li IV
  case ( 17); PFGROUND = 1.                    ! Li V
  case ( 18); PFGROUND = 1.                    ! Li VI

  ! --- Z= 4  Be ---
  case ( 19); PFGROUND = 1.                    ! Be I
  case ( 20); PFGROUND = 2.                    ! Be II
  case ( 21); PFGROUND = 1.                    ! Be III
  case ( 22); PFGROUND = 2.                    ! Be IV
  case ( 23); PFGROUND = 1.                    ! Be V
  case ( 24); PFGROUND = 1.                    ! Be VI

  ! --- Z= 5  B  ---
  case ( 25)  ! B  I
    PFGROUND = 2.+4.*EXP(-HCK/T*15.254)
  case ( 26)  ! B  II
    PFGROUND = 1.+1.*EXP(-HCK/T*37336.7)+3.*EXP(-HCK/T*37342.4)+ 5.*EXP(-HCK/T* 37358.3)+3.*EXP(-HCK/T*73396.60)
  case ( 27)  ! B  III
    PFGROUND = 2.+2.*EXP(-HCK/T*48358.40)+4.*EXP(-HCK/T*48392.50)
  case ( 28); PFGROUND = 1.                    ! B  IV
  case ( 29); PFGROUND = 2.                    ! B  V
  case ( 30); PFGROUND = 1.                    ! B  VI

  ! --- Z= 6  C  ---
  case ( 31)  ! C  I
    PFGROUND = 1.+3.*EXP(-HCK/T*16.40)+5.*EXP(-HCK/T*43.40)+ 5.*EXP(-HCK/T*10192.63)
  case ( 32)  ! C  II
    PFGROUND = 2.+4.*EXP(-HCK/T*63.42)
  case ( 33); PFGROUND = 1.                    ! C  III
  case ( 34)  ! C  IV
    PFGROUND = 2.+2.*EXP(-HCK/T*64484.0)+4.*EXP(-HCK/T*64591.7)
  case ( 35); PFGROUND = 1.                    ! C  V
  case ( 36); PFGROUND = 2.                    ! C  VI

  ! --- Z= 7  N  ---
  case ( 37); PFGROUND = 4.                    ! N  I
  case ( 38)  ! N  II
    PFGROUND = 1.+3.*EXP(-HCK/T*48.7)+5.*EXP(-HCK/T*130.8)+ 5.*EXP(-HCK/T*15316.2)
  case ( 39)  ! N  III
    PFGROUND = 2.+4.*EXP(-HCK/T*174.4)
  case ( 40); PFGROUND = 1.                    ! N  IV
  case ( 41)  ! N  V
    PFGROUND = 2.+2.*EXP(-HCK/T*80463.2)+4.*EXP(-HCK/T*80721.9)
  case ( 42); PFGROUND = 1.                    ! N  VI

  ! --- Z= 8  O  ---
  case ( 43)  ! O  I
    PFGROUND = 5.+3.*EXP(-HCK/T*158.265)+EXP(-HCK/T*226.977)
  case ( 44)  ! O  II
    PFGROUND = 4.+6.*EXP(-HCK/T*26810.55)+4.*EXP(-HCK/T*26830.57)
  case ( 45)  ! O  III
    PFGROUND = 1.+3.*EXP(-HCK/T*113.178)+5.*EXP(-HCK/T*306.174)+ 5.*EXP(-HCK/T*20273.27)+1.*EXP(-HCK/T*43185.74)+ &
      5.*EXP(-HCK/T*60324.79)
  case ( 46)  ! O  IV
    PFGROUND = 2.+4.*EXP(-HCK/T*385.9)+2.*EXP(-HCK/T*71439.8)+ 4.*EXP(-HCK/T*71570.1)+6.*EXP(-HCK/T*71755.5)
  case ( 47)  ! O  V
    PFGROUND = 1.+1.*EXP(-HCK/T*81942.5)+3.*EXP(-HCK/T*82078.6)+ 5.*EXP(-HCK/T*82385.3)
  case ( 48)  ! O  VI
    PFGROUND = 2.+2.*EXP(-HCK/T*96375.0)+4.*EXP(-HCK/T*96907.5)

  ! --- Z= 9  F  ---
  case ( 49)  ! F  I
    PFGROUND = 4.+2.*EXP(-HCK/T*404.1)
  case ( 50)  ! F  II
    PFGROUND = 5.+3.*EXP(-HCK/T*341.0)+EXP(-HCK/T*489.9)+ 5.*EXP(-HCK/T*20873.4)+1.*EXP(-HCK/T*44918.1)
  case ( 51)  ! F  III
    PFGROUND = 4.+6.*EXP(-HCK/T*34087.4)+4.*EXP(-HCK/T*34123.2)+ 4.*EXP(-HCK/T*51561.4)+2.*EXP(-HCK/T*51560.6)
  case ( 52)  ! F  IV
    PFGROUND = 1.+3.*EXP(-HCK/T*225.2)+5.*EXP(-HCK/T*612.2)+ 5.*EXP(-HCK/T*25238.2)+1.*EXP(-HCK/T*53541.2)+ 5.*EXP(-HCK/T*74194.7)
  case ( 53)  ! F  V
    PFGROUND = 2.+4.*EXP(-HCK/T*744.5)+ 2.*EXP(-HCK/T*85790.2)+4.*EXP(-HCK/T*86043.5)+ 6.*EXP(-HCK/T*86407.0)
  case ( 54)  ! F  VI
    PFGROUND = 1.+1.*EXP(-HCK/T*96590.)+3.*EXP(-HCK/T*96850.)+ 5.*EXP(-HCK/T*97427.)

  ! --- Z=10  Ne ---
  case ( 55); PFGROUND = 1.                    ! Ne I
  case ( 56)  ! Ne II
    PFGROUND = 4.+2.*EXP(-HCK/T*780.45)
  case ( 57)  ! Ne III
    PFGROUND = 5.+3.*EXP(-HCK/T*642.9)+EXP(-HCK/T*920.4)+ 4.*EXP(-HCK/T*96907.5)+5.*EXP(-HCK/T*25840.8)+ 1.*EXP(-HCK/T*55750.6)
  case ( 58)  ! Ne IV
    PFGROUND = 4.+6.*EXP(-HCK/T*41234.6)+4.*EXP(-HCK/T*41279.5)+ 2.*EXP(-HCK/T*62434.6)+4.*EXP(-HCK/T*62441.3)
  case ( 59)  ! Ne V
    PFGROUND = 1.+3.*EXP(-HCK/T*414.)+5.*EXP(-HCK/T*1112.)+ 5.*EXP(-HCK/T*30291.5)+1.*EXP(-HCK/T*63913.6)+ 5.*EXP(-HCK/T*88360.)
  case ( 60)  ! Ne VI
    PFGROUND = 2.+4.*EXP(-HCK/T*1310.)+2.*EXP(-HCK/T*100261.)+ 4.*EXP(-HCK/T*100704.)+6.*EXP(-HCK/T*101347.)

  ! --- Z=11  Na ---
  case ( 61); PFGROUND = 2.                    ! Na I
  case ( 62); PFGROUND = 1.                    ! Na II
  case ( 63)  ! Na III
    PFGROUND = 4.+2.*EXP(-HCK/T*780.45)
  case ( 64)  ! Na IV
    PFGROUND = 5.+3.*EXP(-HCK/T*642.9)+EXP(-HCK/T*920.4)+ 5.*EXP(-HCK/T*30839.8)+1.*EXP(-HCK/T*66496.)
  case ( 65)  ! Na V
    PFGROUND = 4.+6.*EXP(-HCK/T*48330.)+4.*EXP(-HCK/T*48366.)+ 2.*EXP(-HCK/T*73218.)+4.*EXP(-HCK/T*73255.)
  case ( 66)  ! Na VI
    PFGROUND = 1.+3.*EXP(-HCK/T*414.)+5.*EXP(-HCK/T*1112.)+ 5.*EXP(-HCK/T*35498.)+1.*EXP(-HCK/T*74414.)

  ! --- Z=12  Mg ---
  case ( 67); PFGROUND = 1.                    ! Mg I
  case ( 68); PFGROUND = 2.                    ! Mg II
  case ( 69); PFGROUND = 1.                    ! Mg III
  case ( 70)  ! Mg IV
    PFGROUND = 4.+2.*EXP(-HCK/T*2238.)
  case ( 71)  ! Mg V
    PFGROUND = 5.+3.*EXP(-HCK/T*1782.1)+EXP(-HCK/T*2521.8)+ 5.*EXP(-HCK/T*35926.)+1.*EXP(-HCK/T*77279.)
  case ( 72)  ! Mg VI
    PFGROUND = 4.+6.*EXP(-HCK/T*55356.)+4.*EXP(-HCK/T*55372.8)+ 2.*EXP(-HCK/T*83920.0)+4.*EXP(-HCK/T*84028.4)

  ! --- Z=13  Al ---
  case ( 73)  ! Al I
    PFGROUND = 2.+4.*EXP(-HCK/T*112.061)
  case ( 74); PFGROUND = 1.                    ! Al II
  case ( 75); PFGROUND = 2.                    ! Al III
  case ( 76); PFGROUND = 1.                    ! Al IV
  case ( 77)  ! Al V
    PFGROUND = 4.+2.*EXP(-HCK/T*3442.)
  case ( 78)  ! Al VI
    PFGROUND = 5.+3.*EXP(-HCK/T*2732.)+EXP(-HCK/T*3829.)+ 5.*EXP(-HCK/T*41167.)+1.*EXP(-HCK/T*88213.)

  ! --- Z=14  Si ---
  case ( 79)  ! Si I
    PFGROUND = 1.+3.*EXP(-HCK/T*77.115)+5.*EXP(-HCK/T*223.157)+ 5.*EXP(-HCK/T*6298.850)
  case ( 80)  ! Si II
    PFGROUND = 2.+4.*EXP(-HCK/T*287.32)+2.*EXP(-HCK/T*42824.35)+ 4.*EXP(-HCK/T*42932.68)+6.*EXP(-HCK/T*43107.97)
  case ( 81)  ! Si III
    PFGROUND = 1.+1.*EXP(-HCK/T*52724.69)+3.*EXP(-HCK/T*52853.28)+ 5.*EXP(-HCK/T*53115.01)+3.*EXP(-HCK/T*82884.41)
  case ( 82); PFGROUND = 2.                    ! Si IV
  case ( 83); PFGROUND = 1.                    ! Si V
  case ( 84)  ! Si VI
    PFGROUND = 4.+2.*EXP(-HCK/T*5090.)

  ! --- Z=15  P  ---
  case ( 85); PFGROUND = 4.                    ! P  I
  case ( 86)  ! P  II
    PFGROUND = 1.+3.*EXP(-HCK/T*164.90)+5.*EXP(-HCK/T*469.12)+ 5.*EXP(-HCK/T*8882.31)+1.*EXP(-HCK/T*21575.63)
  case ( 87)  ! P  III
    PFGROUND = 2.+4.*EXP(-HCK/T*559.14)+2.*EXP(-HCK/T*56021.67)+ 4.*EXP(-HCK/T*57125.98)+6.*EXP(-HCK/T*57454.00)+ &
      4.*EXP(-HCK/T*74916.85)+6.*EXP(-HCK/T*74945.86)
  case ( 88)  ! P  IV
    PFGROUND = 1.+1.*EXP(-HCK/T*67918.03)+3.*EXP(-HCK/T*68146.48)+ 5.*EXP(-HCK/T*68615.17)
  case ( 89)  ! P  V
    PFGROUND = 2.+2.*EXP(-HCK/T*88651.87)+4.*EXP(-HCK/T*89447.25)
  case ( 90); PFGROUND = 1.                    ! P  VI

  ! --- Z=16  S  ---
  case ( 91)  ! S  I
    PFGROUND = 5.+3.*EXP(-HCK/T*396.055)+EXP(-HCK/T*573.640)+ 5.*EXP(-HCK/T*9238.609)
  case ( 92)  ! S  II
    PFGROUND = 4.+4.*EXP(-HCK/T*14852.94)+6.*EXP(-HCK/T*14884.73)+ 2.*EXP(-HCK/T*24524.83)+4.*EXP(-HCK/T*24571.54)
  case ( 93)  ! S  III
    PFGROUND = 1.+3.*EXP(-HCK/T*298.69)+5.*EXP(-HCK/T*833.08)+ 5.*EXP(-HCK/T*11322.7)+1.*EXP(-HCK/T*27161.0)
  case ( 94)  ! S  IV
    PFGROUND = 2.+4.*EXP(-HCK/T*951.43)+2.*EXP(-HCK/T*71184.1)+ 4.*EXP(-HCK/T*71528.7)+6.*EXP(-HCK/T*72074.4)+ &
      4.*EXP(-HCK/T*94103.1)+6.*EXP(-HCK/T*94150.4)
  case ( 95)  ! S  V
    PFGROUND = 1.+1.*EXP(-HCK/T*83024.0)+3.*EXP(-HCK/T*83393.5)+ 5.*EXP(-HCK/T*84155.2)
  case ( 96); PFGROUND = 2.                    ! S  VI

  ! --- Z=17  Cl ---
  case ( 97)  ! Cl I
    PFGROUND = 4.+2.*EXP(-HCK/T*882.36)
  case ( 98)  ! Cl II
    PFGROUND = 5.+3.*EXP(-HCK/T*696.1)+EXP(-HCK/T*996.4)+ 5.*EXP(-HCK/T*11653.58)+1.*EXP(-HCK/T*27878.02)
  case ( 99)  ! Cl III
    PFGROUND = 4.+4.*EXP(-HCK/T*18053.)+6.*EXP(-HCK/T*18118.6)+ 2.*EXP(-HCK/T*29812.)+4.*EXP(-HCK/T*29907.)
  case (100)  ! Cl IV
    PFGROUND = 1.+3.*EXP(-HCK/T*491.)+5.*EXP(-HCK/T*1341.)+ 5.*EXP(-HCK/T*13767.6)+1.*EXP(-HCK/T*32547.8)+ 5.*EXP(-HCK/T*65000.)
  case (101)  ! Cl V
    PFGROUND = 2.+4.*EXP(-HCK/T*1490.8)+2.*EXP(-HCK/T*86000.)+ 4.*EXP(-HCK/T*86538.)+6.*EXP(-HCK/T*87381.)
  case (102)  ! Cl VI
    PFGROUND = 1.+1.*EXP(-HCK/T*97405.)+3.*EXP(-HCK/T*97958.)+ 5.*EXP(-HCK/T*99123.)

  ! --- Z=18  Ar ---
  case (103); PFGROUND = 1.                    ! Ar I
  case (104)  ! Ar II
    PFGROUND = 4.+2.*EXP(-HCK/T*1431.41)
  case (105)  ! Ar III
    PFGROUND = 5.+3.*EXP(-HCK/T*1112.1)+EXP(-HCK/T*1570.2)+ 5.*EXP(-HCK/T*14010.004)+1.*EXP(-HCK/T*33265.724)
  case (106)  ! Ar IV
    PFGROUND = 4.+4.*EXP(-HCK/T*21090.4)+6.*EXP(-HCK/T*21219.3)+ 2.*EXP(-HCK/T*34855.5)+4.*EXP(-HCK/T*35032.6)
  case (107)  ! Ar V
    PFGROUND = 1.+3.*EXP(-HCK/T*765.)+5.*EXP(-HCK/T*2030.)+ 5.*EXP(-HCK/T*16298.9)+1.*EXP(-HCK/T*37912.0)+ 5.*EXP(-HCK/T*84100.0)
  case (108)  ! Ar VI
    PFGROUND = 2.+4.*EXP(-HCK/T*2208.)

  ! --- Z=19  K  ---
  case (109); PFGROUND = 2.                    ! K  I
  case (110); PFGROUND = 1.                    ! K  II
  case (111)  ! K  III
    PFGROUND = 4.+2.*EXP(-HCK/T*2166.)
  case (112)  ! K  IV
    PFGROUND = 5.+3.*EXP(-HCK/T*1673.)+EXP(-HCK/T*2325.)+ 5.*EXP(-HCK/T*16384.1)+EXP(-HCK/T*38546.3)
  case (113)  ! K  V
    PFGROUND = 4.+4.*EXP(-HCK/T*24012.5)+6.*EXP(-HCK/T*24249.6)+ 2.*EXP(-HCK/T*39758.1)+4.*EXP(-HCK/T*40080.2)
  case (114)  ! K  VI
    PFGROUND = 1.+3.*EXP(-HCK/T*1132.)+5.*EXP(-HCK/T*2924.)+ 5.*EXP(-HCK/T*18977.8)+1.*EXP(-HCK/T*43358.8)

  ! --- Z=29  Cu ---
  case (169); PFGROUND = 2.                    ! Cu I
  case (170); PFGROUND = 1.                    ! Cu II
  case (171)  ! Cu III
    PFGROUND = 6.+4.*EXP(-HCK/T*2071.8)
  case (172); PFGROUND = 1.                    ! Cu IV
  case (173); PFGROUND = 1.                    ! Cu V
  case (174); PFGROUND = 1.                    ! Cu VI

  ! --- Z=30  Zn ---
  case (175); PFGROUND = 1.                    ! Zn I
  case (176); PFGROUND = 2.                    ! Zn II
  case (177); PFGROUND = 1.                    ! Zn III
  case (178); PFGROUND = 1.                    ! Zn IV
  case (179); PFGROUND = 1.                    ! Zn V
  case (180); PFGROUND = 1.                    ! Zn VI

  ! --- Z=31  Ga ---
  case (181)  ! Ga I
    PFGROUND = 2.+4.*EXP(-HCK/T*826.19)
  case (182); PFGROUND = 1.                    ! Ga II
  case (183); PFGROUND = 2.                    ! Ga III
  case (184); PFGROUND = 1.                    ! Ga IV
  case (185); PFGROUND = 1.                    ! Ga V
  case (186); PFGROUND = 1.                    ! Ga VI

  ! --- Z=32  Ge ---
  case (187)  ! Ge I
    PFGROUND = 1.+3.*EXP(-HCK/T*557.134)+5.*EXP(-HCK/T*1409.961)+ 5.*EXP(-HCK/T*7125.299)
  case (188)  ! Ge II
    PFGROUND = 2.+4.*EXP(-HCK/T*1767.356)
  case (189); PFGROUND = 1.                    ! Ge III
  case (190); PFGROUND = 1.                    ! Ge IV
  case (191); PFGROUND = 1.                    ! Ge V
  case (192); PFGROUND = 1.                    ! Ge VI

  ! --- Z=33  As ---
  case (193); PFGROUND = 4.                    ! As I
  case (194)  ! As II
    PFGROUND = 1.+3.*EXP(-HCK/T*1061.)+5.*EXP(-HCK/T*2538.)
  case (195)  ! As III
    PFGROUND = 2.+4.*EXP(-HCK/T*2940.)
  case (196); PFGROUND = 1.                    ! As IV
  case (197); PFGROUND = 1.                    ! As V
  case (198); PFGROUND = 1.                    ! As VI

  ! --- Z=34  Se ---
  case (199)  ! Se I
    PFGROUND = 5.+3.*EXP(-HCK/T*1989.49)+EXP(-HCK/T*2534.35)+ 5.*EXP(-HCK/T*9576.149)+1.*EXP(-HCK/T*22446.202)
  case (200)  ! Se II
    PFGROUND = 4.+4.*EXP(-HCK/T*13168.2)+6.*EXP(-HCK/T*13784.4)+ 2.*EXP(-HCK/T*23038.3)+4.*EXP(-HCK/T*23894.8)
  case (201)  ! Se III
    PFGROUND = 1.+3.*EXP(-HCK/T*1741.)+5.*EXP(-HCK/T*3937.)+ 5.*EXP(-HCK/T*13032.)+1.*EXP(-HCK/T*28430.)
  case (202); PFGROUND = 1.                    ! Se IV
  case (203); PFGROUND = 1.                    ! Se V
  case (204); PFGROUND = 1.                    ! Se VI

  ! --- Z=35  Br ---
  case (205)  ! Br I
    PFGROUND = 4.+2.*EXP(-HCK/T*3685.24)
  case (206)  ! Br II
    PFGROUND = 5.+3.*EXP(-HCK/T*3136.4)+EXP(-HCK/T*3837.5)+ 5.*EXP(-HCK/T*12089.1)
  case (207)  ! Br III
    PFGROUND = 4.+4.*EXP(-HCK/T*15042.0)+6.*EXP(-HCK/T*16301.0)+ 2.*EXP(-HCK/T*26915.0)+4.*EXP(-HCK/T*28579.0)
  case (208); PFGROUND = 1.                    ! Br IV
  case (209); PFGROUND = 1.                    ! Br V
  case (210); PFGROUND = 1.                    ! Br VI

  ! --- Z=36  Kr ---
  case (211); PFGROUND = 1.                    ! Kr I
  case (212)  ! Kr II
    PFGROUND = 4.+2.*EXP(-HCK/T*5371.)
  case (213)  ! Kr III
    PFGROUND = 5.+3.*EXP(-HCK/T*3136.4)+EXP(-HCK/T*3837.5)+ 5.*EXP(-HCK/T*14644.3)+1.*EXP(-HCK/T*33079.6)
  case (214); PFGROUND = 1.                    ! Kr IV
  case (215); PFGROUND = 1.                    ! Kr V
  case (216); PFGROUND = 1.                    ! Kr VI

  ! --- Z=37  Rb ---
  case (217); PFGROUND = 2.                    ! Rb I
  case (218); PFGROUND = 1.                    ! Rb II
  case (219)  ! Rb III
    PFGROUND = 4.+2.*EXP(-HCK/T*7380.)
  case (220); PFGROUND = 1.                    ! Rb IV
  case (221); PFGROUND = 1.                    ! Rb V
  case (222); PFGROUND = 1.                    ! Rb VI

  ! --- Z=38  Sr ---
  case (223); PFGROUND = 1.                    ! Sr I
  case (224)  ! Sr II
    PFGROUND = 2.+4.*EXP(-HCK/T*14555.50)+6.*EXP(-HCK/T*14836.24)
  case (225); PFGROUND = 1.                    ! Sr III
  case (226); PFGROUND = 1.                    ! Sr IV
  case (227); PFGROUND = 1.                    ! Sr V
  case (228); PFGROUND = 1.                    ! Sr VI

  ! --- Z=39  Y  ---
  case (229)  ! Y  I
    PFGROUND = 4.+6.*EXP(-HCK/T*530.36)
  case (230)  ! Y  II
    PFGROUND = 1.+3.*EXP(-HCK/T*840.198)+5.*EXP(-HCK/T*1045.076)+ 7.*EXP(-HCK/T*1449.752)+5.*EXP(-HCK/T*3296.280)+ &
      5.*EXP(-HCK/T*8003.126)+7.*EXP(-HCK/T*8328.039)+ 9.*EXP(-HCK/T*8743.322)
  case (231)  ! Y  III
    PFGROUND = 4.+6.*EXP(-HCK/T*724.15)+2.*EXP(-HCK/T*7467.10)
  case (232); PFGROUND = 1.                    ! Y  IV
  case (233); PFGROUND = 1.                    ! Y  V
  case (234); PFGROUND = 1.                    ! Y  VI

  ! --- Z=40  Zr ---
  case (235)  ! Zr I
    PFGROUND = 5.+7.*EXP(-HCK/T*570.41)+9.*EXP(-HCK/T*1240.84)+ 1.*EXP(-HCK/T*4196.85)+3.*EXP(-HCK/T*4376.28)+ &
      5.*EXP(-HCK/T*4186.11)+3.*EXP(-HCK/T*4870.53)+ 5.*EXP(-HCK/T*5023.41)+7.*EXP(-HCK/T*5249.07)+ 9.*EXP(-HCK/T*5540.54)+ &
      11.*EXP(-HCK/T*5888.93)+ 5.*EXP(-HCK/T*5101.68)+9.*EXP(-HCK/T*8057.30)
  case (236)  ! Zr II
    PFGROUND = 4.+6.*EXP(-HCK/T*314.67)+8.*EXP(-HCK/T*763.44)+ 10.*EXP(-HCK/T*1322.91)+4.*EXP(-HCK/T*2572.21)+ &
      6.*EXP(-HCK/T*2895.00)+8.*EXP(-HCK/T*3299.58)+ 10.*EXP(-HCK/T*3757.63)+4.*EXP(-HCK/T*4247.97)+ 6.*EXP(-HCK/T*4505.30)+ &
      2.*EXP(-HCK/T*5723.78)+ 4.*EXP(-HCK/T*6111.16)+6.*EXP(-HCK/T*5752.55)+ 8.*EXP(-HCK/T*6467.10)+2.*EXP(-HCK/T*7512.61)+ &
      4.*EXP(-HCK/T*7736.05)+6.*EXP(-HCK/T*8058.27)+ 8.*EXP(-HCK/T*7837.49)+10.*EXP(-HCK/T*8152.57)+ 2.*EXP(-HCK/T*9553.13)+ &
      4.*EXP(-HCK/T*9742.80)+ 6.*EXP(-HCK/T*9968.75)
  case (237)  ! Zr III
    PFGROUND = 5.+7.*EXP(-HCK/T*681.2)+9.*EXP(-HCK/T*1485.8)+ 5.*EXP(-HCK/T*5742.8)+1.*EXP(-HCK/T*8062.7)+ &
      3.*EXP(-HCK/T*8327.0)+5.*EXP(-HCK/T*8839.7)+ 9.*EXP(-HCK/T*11049.9)+1.*EXP(-HCK/T*23974.9)+ 3.*EXP(-HCK/T*18400.8)+ &
      5.*EXP(-HCK/T*18804.7)+ 7.*EXP(-HCK/T*19535.3)+5.*EXP(-HCK/T*25066.9)+ 1.*EXP(-HCK/T*36473.7)
  case (238); PFGROUND = 1.                    ! Zr IV
  case (239); PFGROUND = 1.                    ! Zr V
  case (240); PFGROUND = 1.                    ! Zr VI

  ! --- Z=41  Nb ---
  case (241)  ! Nb I
    PFGROUND = 2.+4.*EXP(-HCK/T*154.19)+6.*EXP(-HCK/T*391.99)+ 8.*EXP(-HCK/T*695.25)+10.*EXP(-HCK/T*1050.26)+ &
      4.*EXP(-HCK/T*1142.79)+6.*EXP(-HCK/T*1586.90)+ 8.*EXP(-HCK/T*2154.11)+10.*EXP(-HCK/T*2805.36)+ 2.*EXP(-HCK/T*4998.17)+ &
      4.*EXP(-HCK/T*5297.92)+ 6.*EXP(-HCK/T*5965.45)+2.*EXP(-HCK/T*8410.90)+ 4.*EXP(-HCK/T*8705.32)+6.*EXP(-HCK/T*9043.14)+ &
      8.*EXP(-HCK/T*9497.52)+8.*EXP(-HCK/T*8827.00)+ 10.*EXP(-HCK/T*9328.88)+4.*EXP(-HCK/T*9439.08)+ 6.*EXP(-HCK/T*10237.51)
  case (242)  ! Nb II
    PFGROUND = 1.+3.*EXP(-HCK/T*158.99)+5.*EXP(-HCK/T*438.38)+ 7.*EXP(-HCK/T*801.38)+9.*EXP(-HCK/T*1224.87)+ &
      3.*EXP(-HCK/T*2356.76)+5.*EXP(-HCK/T*2629.07)+ 7.*EXP(-HCK/T*3029.57)+9.*EXP(-HCK/T*3542.50)+ 11.*EXP(-HCK/T*4146.00)+ &
      1.*EXP(-HCK/T*5562.26)+ 3.*EXP(-HCK/T*6192.33)+5.*EXP(-HCK/T*7261.33)+ 5.*EXP(-HCK/T*7505.78)+7.*EXP(-HCK/T*7900.65)+ &
      9.*EXP(-HCK/T*8320.40)+9.*EXP(-HCK/T*9509.67)+ 11.*EXP(-HCK/T*9812.56)+13.*EXP(-HCK/T*10186.41)
  case (243)  ! Nb III
    PFGROUND = 4.+6.*EXP(-HCK/T*515.8)+8.*EXP(-HCK/T*1176.6)+ 10.*EXP(-HCK/T*1939.0)+ 2.*EXP(-HCK/T*8664.3)+ &
      4.*EXP(-HCK/T*8607.5)+ 6.*EXP(-HCK/T*9593.7)+8.*EXP(-HCK/T*9236.1)+ 10.*EXP(-HCK/T*9804.5)+4.*EXP(-HCK/T*10912.2)+ &
      6.*EXP(-HCK/T*13094.0)+10.*EXP(-HCK/T*12916.0)+ 12.*EXP(-HCK/T*13263.8)+6.*EXP(-HCK/T*19975.0)+ 8.*EXP(-HCK/T*19861.0)+ &
      4.*EXP(-HCK/T*25220.2)+ 6.*EXP(-HCK/T*25735.2)+8.*EXP(-HCK/T*26463.7)+ 10.*EXP(-HCK/T*27373.5)
  case (244); PFGROUND = 1.                    ! Nb IV
  case (245); PFGROUND = 1.                    ! Nb V
  case (246); PFGROUND = 1.                    ! Nb VI

  ! --- Z=42  Mo ---
  case (247); PFGROUND = 7.                    ! Mo I
  case (248); PFGROUND = 6.                    ! Mo II
  case (249)  ! Mo III
    PFGROUND = 1.+3.*EXP(-HCK/T*243.10)+5.*EXP(-HCK/T*669.60)+ 7.*EXP(-HCK/T*1225.20)+9.*EXP(-HCK/T*1873.80)
  case (250); PFGROUND = 1.                    ! Mo IV
  case (251); PFGROUND = 1.                    ! Mo V
  case (252); PFGROUND = 1.                    ! Mo VI

  ! --- Z=43  Tc ---
  case (253)  ! Tc I
    PFGROUND = 6.   +10.*EXP(-HCK/T*2572.89)+8.*EXP(-HCK/T*3250.91)+ 6.*EXP(-HCK/T*3700.54)+4.*EXP(-HCK/T*4002.57)+ &
      2.*EXP(-HCK/T*4178.75)
  case (254)  ! Tc II
    PFGROUND = 7.  +9.*EXP(-HCK/T*3461.27)+7.*EXP(-HCK/T*4217.17)+ 5.*EXP(-HCK/T*4669.22)+3.*EXP(-HCK/T*4961.14)+ &
      1.*EXP(-HCK/T*5100.98)
  case (255); PFGROUND = 6.                    ! Tc III
  case (256); PFGROUND = 1.                    ! Tc IV
  case (257); PFGROUND = 1.                    ! Tc V
  case (258); PFGROUND = 1.                    ! Tc VI

  ! --- Z=44  Ru ---
  case (259)  ! Ru I
    PFGROUND = 11.+9.*EXP(-HCK/T*1190.64)+7.*EXP(-HCK/T*2091.54)+ 5.*EXP(-HCK/T*2713.24)+3.*EXP(-HCK/T*3105.49)+ &
      9.*EXP(-HCK/T*6545.03)+7.*EXP(-HCK/T*8084.12)+ 5.*EXP(-HCK/T*9183.66)+9.*EXP(-HCK/T*7483.07)+ 7.*EXP(-HCK/T*8575.42)+ &
      5.*EXP(-HCK/T*9057.64)+ 3.*EXP(-HCK/T*9072.98)+1.*EXP(-HCK/T*8492.37)+ 7.*EXP(-HCK/T*8770.93)+5.*EXP(-HCK/T*8043.69)+ &
      3.*EXP(-HCK/T*9620.29)+9.*EXP(-HCK/T*9120.63)
  case (260)  ! Ru II
    PFGROUND = 10.+8.*EXP(-HCK/T*1523.1)+6.*EXP(-HCK/T*2493.9)+ 4.*EXP(-HCK/T*3104.2)+6.*EXP(-HCK/T*8256.7)+ &
      4.*EXP(-HCK/T*8477.4)+2.*EXP(-HCK/T*9373.4)+ 10.*EXP(-HCK/T*9151.6)
  case (261)  ! Ru III
    PFGROUND = 9.+7.*EXP(-HCK/T*1158.8)+5.*EXP(-HCK/T*1826.3)+ 3.*EXP(-HCK/T*2266.3)+EXP(-HCK/T*2476.0)+ 7.*EXP(-HCK/T*27162.8)+ &
      5.*EXP(-HCK/T*41111.7)
  case (262); PFGROUND = 1.                    ! Ru IV
  case (263); PFGROUND = 1.                    ! Ru V
  case (264); PFGROUND = 1.                    ! Ru VI

  ! --- Z=45  Rh ---
  case (265)  ! Rh I
    PFGROUND = 10.+8.*EXP(-HCK/T*1529.97)+6.*EXP(-HCK/T*2598.03)+ 4.*EXP(-HCK/T*3472.68)+6.*EXP(-HCK/T*3309.86)+ &
      4.*EXP(-HCK/T*5657.97)+8.*EXP(-HCK/T*5690.97)+ 6.*EXP(-HCK/T*7791.23)+6.*EXP(-HCK/T*9221.22)
  case (266)  ! Rh II
    PFGROUND = 9.+7.*EXP(-HCK/T*2401.3)+5.*EXP(-HCK/T*3580.7)+ 5.*EXP(-HCK/T*8164.4)+1.*EXP(-HCK/T*10760.8)+ &
      3.*EXP(-HCK/T*10515.0)+5.*EXP(-HCK/T*11643.7)+ 9.*EXP(-HCK/T*14855.4)+11.*EXP(-HCK/T*16884.8)+ 9.*EXP(-HCK/T*18540.4)+ &
      7.*EXP(-HCK/T*19792.4)
  case (267)  ! Rh III
    PFGROUND = 10.+8.*EXP(-HCK/T*2147.8)+6.*EXP(-HCK/T*3485.7)+ 4.*EXP(-HCK/T*4322.0)+6.*EXP(-HCK/T*11062.3)+ &
      4.*EXP(-HCK/T*10997.1)+2.*EXP(-HCK/T*12469.8)+ 10.*EXP(-HCK/T*14044.0)+8.*EXP(-HCK/T*15256.8)+ 4.*EXP(-HCK/T*16870.7)+ &
      2.*EXP(-HCK/T*18303.7)+ 12.*EXP(-HCK/T*19490.2)+6.*EXP(-HCK/T*19528.5)
  case (268); PFGROUND = 1.                    ! Rh IV
  case (269); PFGROUND = 1.                    ! Rh V
  case (270); PFGROUND = 1.                    ! Rh VI

  ! --- Z=46  Pd ---
  case (271)  ! Pd I
    PFGROUND = 1.+7.*EXP(-HCK/T*6564.11)+5.*EXP(-HCK/T*7754.99)
  case (272)  ! Pd II
    PFGROUND = 6.+4.*EXP(-HCK/T*3539.2)
  case (273)  ! Pd III
    PFGROUND = 9.+7.*EXP(-HCK/T*3229.3)+5.*EXP(-HCK/T*4687.5)+ 5.*EXP(-HCK/T*10229.3)+3.*EXP(-HCK/T*13468.9)+ &
      1.*EXP(-HCK/T*13697.5)+5.*EXP(-HCK/T*14634.4)+ 9.*EXP(-HCK/T*17879.3)
  case (274); PFGROUND = 1.                    ! Pd IV
  case (275); PFGROUND = 1.                    ! Pd V
  case (276); PFGROUND = 1.                    ! Pd VI

  ! --- Z=47  Ag ---
  case (277); PFGROUND = 2.                    ! Ag I
  case (278); PFGROUND = 1.                    ! Ag II
  case (279)  ! Ag III
    PFGROUND = 6.+4.*EXP(-HCK/T*4607.)
  case (280); PFGROUND = 1.                    ! Ag IV
  case (281); PFGROUND = 1.                    ! Ag V
  case (282); PFGROUND = 1.                    ! Ag VI

  ! --- Z=48  Cd ---
  case (283); PFGROUND = 1.                    ! Cd I
  case (284); PFGROUND = 2.                    ! Cd II
  case (285); PFGROUND = 1.                    ! Cd III
  case (286); PFGROUND = 1.                    ! Cd IV
  case (287); PFGROUND = 1.                    ! Cd V
  case (288); PFGROUND = 1.                    ! Cd VI

  ! --- Z=49  In ---
  case (289)  ! In I
    PFGROUND = 2.+4.*EXP(-HCK/T*2212.598)
  case (290); PFGROUND = 1.                    ! In II
  case (291); PFGROUND = 2.                    ! In III
  case (292); PFGROUND = 1.                    ! In IV
  case (293); PFGROUND = 1                    ! In V
  case (294); PFGROUND = 1.                    ! In VI

  ! --- Z=50  Sn ---
  case (295)  ! Sn I
    PFGROUND = 1.+3.*EXP(-HCK/T*1691.8)+5.*EXP(-HCK/T*3427.7)+ 5.*EXP(-HCK/T*6513.0)
  case (296)  ! Sn II
    PFGROUND = 2.+4.*EXP(-HCK/T*4251.4)
  case (297); PFGROUND = 1.                    ! Sn III
  case (298); PFGROUND = 1.                    ! Sn IV
  case (299); PFGROUND = 1.                    ! Sn V
  case (300); PFGROUND = 1.                    ! Sn VI

  ! --- Z=51  Sb ---
  case (301)  ! Sb I
    PFGROUND = 4.+4.*EXP(-HCK/T*8512.1)+6.*EXP(-HCK/T*9854.1)
  case (302)  ! Sb II
    PFGROUND = 1.+3.*EXP(-HCK/T*3055.0)+5.*EXP(-HCK/T*5659.0)
  case (303)  ! Sb III
    PFGROUND = 2.+4.*EXP(-HCK/T*6576.)
  case (304); PFGROUND = 1.                    ! Sb IV
  case (305); PFGROUND = 1.                    ! Sb V
  case (306); PFGROUND = 1.                    ! Sb VI

  ! --- Z=52  Te ---
  case (307)  ! Te I
    PFGROUND = 5.+3.*EXP(-HCK/T*4750.712)+EXP(-HCK/T*4706.5)
  case (308)  ! Te II
    PFGROUND = 4.+4.*EXP(-HCK/T*10222.385)+6.*EXP(-HCK/T*12421.854)+ 2.*EXP(-HCK/T*20546.591)+4.*EXP(-HCK/T*24032.2)
  case (309)  ! Te III
    PFGROUND = 1.+3.*EXP(-HCK/T*4756.5)+5.*EXP(-HCK/T*8166.9)+ 5.*EXP(-HCK/T*17358.)
  case (310); PFGROUND = 1.                    ! Te IV
  case (311); PFGROUND = 1.                    ! Te V
  case (312); PFGROUND = 1.                    ! Te VI

  ! --- Z=53  I  ---
  case (313)  ! I  I
    PFGROUND = 4.+2.*EXP(-HCK/T*7063.15)
  case (314)  ! I  II
    PFGROUND = 5.+3.*EXP(-HCK/T*7087.0)+EXP(-HCK/T*6447.9)+ 5.*EXP(-HCK/T*13727.2)+1.*EXP(-HCK/T*29501.3)
  case (315)  ! I  III
    PFGROUND = 4.+4.*EXP(-HCK/T*11711.2)+6.*EXP(-HCK/T*14901.9)+ 2.*EXP(-HCK/T*24299.3)+4.*EXP(-HCK/T*29636.8)
  case (316); PFGROUND = 1.                    ! I  IV
  case (317); PFGROUND = 1.                    ! I  V
  case (318); PFGROUND = 1.                    ! I  VI

  ! --- Z=54  Xe ---
  case (319); PFGROUND = 1.                    ! Xe I
  case (320)  ! Xe II
    PFGROUND = 4.+2.*EXP(-HCK/T*10537.01)
  case (321)  ! Xe III
    PFGROUND = 5.+3.*EXP(-HCK/T*9794.36)+EXP(-HCK/T*8130.08)+ 5.*EXP(-HCK/T*17098.73)+1.*EXP(-HCK/T*36102.94)
  case (322); PFGROUND = 1.                    ! Xe IV
  case (323); PFGROUND = 1.                    ! Xe V
  case (324); PFGROUND = 1.                    ! Xe VI

  ! --- Z=55  Cs ---
  case (325); PFGROUND = 2.                    ! Cs I
  case (326); PFGROUND = 1.                    ! Cs II
  case (327)  ! Cs III
    PFGROUND = 4.+2.*EXP(-HCK/T*13884.)
  case (328); PFGROUND = 1.                    ! Cs IV
  case (329); PFGROUND = 1.                    ! Cs V
  case (330); PFGROUND = 1.                    ! Cs VI

  ! --- Z=56  Ba ---
  case (331); PFGROUND = 1.                    ! Ba I
  case (332)  ! Ba II
    PFGROUND = 2.+4.*EXP(-HCK/T*4873.852)+6.*EXP(-HCK/T*5674.807)
  case (333); PFGROUND = 1.                    ! Ba III
  case (334); PFGROUND = 1.                    ! Ba IV
  case (335); PFGROUND = 1.                    ! Ba V
  case (336); PFGROUND = 1.                    ! Ba VI

  ! Default for all other NELION values
  case default
    PFGROUND = 1.0d0

  end select

END FUNCTION PFGROUND

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

  implicit none

  real*8, intent(in) :: T
  real*8 :: PARTFNH2

  integer, parameter :: NMAX = 500   ! max table entries
  real*8,  save :: PF(NMAX)          ! partition function values
  real*8,  save :: TSTEP = 100.0d0   ! temperature step (K)
  real*8,  save :: TSTART = 100.0d0  ! first temperature in table
  integer, save :: NPF = 0           ! number of entries loaded
  logical, save :: INITIALIZED = .false.
  
  integer :: N, IOS, LUN
  real*8  :: frac, TDUM, QDUM
  character(len=256) :: LINE

  ! Read table from file on first call
  if (.not. INITIALIZED) then
    INITIALIZED = .true.
    NPF = 0
    LUN = 89
    open(LUN, file=trim(DATADIR)//'partfnh2.dat', status='old', iostat=IOS)
    if (IOS /= 0) then
      stop ' PARTFNH2 ERROR: cannot open '//trim(DATADIR)//'partfnh2.dat'
   endif
   do while (NPF < NMAX)
      read(LUN, '(A)', iostat=IOS) LINE
      if (IOS /= 0) exit
      ! Skip comment and blank lines
      LINE = adjustl(LINE)
      if (LINE(1:1) == '#' .or. len_trim(LINE) == 0) cycle
      read(LINE, *, iostat=IOS) TDUM, QDUM
      if (IOS /= 0) cycle
      NPF = NPF + 1
      PF(NPF) = QDUM
      if (NPF == 1) TSTART = TDUM
      if (NPF == 2) TSTEP = TDUM - TSTART
   end do
   close(LUN)
  end if

  ! Fallback if no table loaded
  if (NPF == 0) then
    PARTFNH2 = max(T / 100.0d0, 0.5d0)
    return
  end if

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

  implicit none

  real*8, intent(in) :: T
  real*8 :: EQUILH2

  ! External function

  ! Physical constants from mod_constants; only local spectroscopic data here
  real*8, parameter :: m_H  = 1.008d0 * AMU        ! H atom mass [g]

  ! H2 dissociation energy
  real*8, parameter :: D0_cm = 36118.11d0           ! D0(H2) [cm^-1]
  real*8, parameter :: D0_over_kT_coeff = D0_cm * HCK
  !                                     = 51967.8 K  (D0/k)

  ! Translational partition function prefactor (T-independent part)
  ! = 2^1.5 / 4 / (2*pi*m_H*k/h^2)^1.5
  ! where m_H appears because the reduced mass for equal-mass dissociation
  ! products cancels to give the atomic H mass in the Saha-like expression.
  real*8, parameter :: trans_prefactor = 2.d0**1.5d0 / 4.d0 &
    / (2.d0 * PI * m_H * KBOL / HPLANCK**2)**1.5d0

  real*8 :: Q_H2

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

  implicit none

  real*8, intent(in) :: beta
  real*8 :: holtsmark_Q

  real*8  :: log_beta, idx_f, frac
  integer :: idx

  if (beta <= 0.01D0) then
    holtsmark_Q = 0.0D0
    return
  end if
  if (beta >= 50.0D0) then
    holtsmark_Q = 1.0D0
    return
  end if

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

  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: xne
  real*8 :: occupation_prob

  real*8 :: beta

  if (n <= 1 .or. xne <= 0.0D0) then
    occupation_prob = 1.0D0
    return
  end if

  beta = BETA_COEFF_HM88 / (DBLE(n)**5 * xne**(2.0D0/3.0D0))
  occupation_prob = holtsmark_Q(beta)

END FUNCTION occupation_prob


end module mod_atlas_data

