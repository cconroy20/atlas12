!=========================================================================
! PROGRAM ATLAS12
!
!   Opacity-sampling stellar atmosphere code.
!   Robert L. Kurucz, Harvard-Smithsonian Center for Astrophysics.
!
!   This program iteratively solves for the structure of a plane-parallel,
!   LTE (or NLTE) stellar atmosphere in radiative and hydrostatic
!   equilibrium using opacity sampling over ~30,000 wavelength points.
!
!   The main loop structure is:
!     1. Read input model parameters (READIN)
!     2. For each iteration (1..NUMITS):
!        a. Solve hydrostatic equilibrium for gas pressure
!        b. Compute populations and Doppler widths for all species
!        c. Build opacity tables (first iteration only)
!        d. Compute line opacities (LINOP1, XLINOP)
!        e. Loop over all frequencies:
!           - Compute continuous + line opacity at this frequency
!           - Solve radiative transfer (Feautrier method via JOSH)
!           - Accumulate flux integrals for T-correction
!        f. Compute Rosseland mean, convection, T-correction
!        g. Output results
!     3. Loop back to read next model
!
!=========================================================================

PROGRAM ATLAS12

  USE mod_atlas_data
  IMPLICIT NONE

  ! --- Local variables ---
  REAL(8)        :: VSTEPS
  REAL(8)        :: EXCESS, XMAX
  ! NB: FREQ15, found_negative, RCOWT used to live here; they moved to
  !     RUN_RT_PASS along with the frequency loop body.
  REAL(8)        :: RX
  INTEGER        :: I, J, NU, ITERAT
  INTEGER        :: NUCI, NULYMAN, NUHEI, NUHEII, NUSTART
  INTEGER(8)     :: CLOCK_START, CLOCK_END, CLOCK_RATE
  INTEGER(8)     :: CLOCK_TOTAL_START

  ! --- Environment variable for data directory ---
  CHARACTER(256) :: ENVVAL
  INTEGER        :: ENVLEN, ENVSTAT

  ! --- Command-line output file base name and options ---
  CHARACTER(256) :: OUTBASE, ARGBUF
  CHARACTER(256) :: ABUND_FILE
  CHARACTER(256) :: ABUND_LINE
  CHARACTER(32)  :: CMD_SOLAR
  INTEGER        :: NARGS, IEQPOS, ISTAT, IPOSARG
  REAL(8)        :: VTURB_KMS
  INTEGER        :: ICZC_CLI      ! czc= CLI value (0 = off, nonzero = on)
  REAL(8)        :: CMD_TEFF, CMD_LOGG
  REAL(8)        :: CMD_ZSCALE, CMD_HEABND
  REAL(8)        :: Z_TOTAL
  INTEGER        :: IZ_ABUND, IOS_ABUND
  REAL(8)        :: ABUND_VAL
  LOGICAL        :: have_abund_overrides, have_overrides
  LOGICAL        :: file_exists

  ! --- INITIALIZATION PHASE ---------------------------------------------

  CALL SYSTEM_CLOCK(CLOCK_TOTAL_START, CLOCK_RATE)

  ! --- Parse command-line arguments ---
  !     Usage: atlas12.exe <input_atm> [basename] [numit=N] [vturb=X] [mlt=X]
  !            [teff=X] [logg=X] [solar=name] [zscale=X] [heabnd=X] [abund=file]
  !     Pass --help (or -h, help) to print usage and exit.
  !     input_atm  : input atmosphere model file (REQUIRED, first positional)
  !     basename   : output file base name (default 'mystar')
  !     Output files: <basename>.atm, .flux, .iter
  !       (.iter includes the former .tcorr correction diagnostics; the
  !        separate .taunu and .tcorr files are no longer written)
  !     numit=N    : number of iterations (default 30)
  !     vturb=X    : microturbulence in km/s (default: from model)
  !     mlt=X      : mixing length parameter (default 2.0)
  !     teff=X     : rescale model to this Teff (default: from model)
  !     logg=X     : rescale model to this log g (default: from model)
  !     solar=name : solar reference abundance pattern (default: from
  !                  model file; compiled-in default is ag89).  Replaces
  !                  the pattern the model file carries; zscale/abund/
  !                  heabnd overrides then apply on top of it.
  !     zscale=X   : metal abundance scale factor (default: no scaling)
  !     heabnd=X   : He number fraction Y (default: from model);
  !                  H is computed as X = 1 - Y - Z for consistency
  !     abund=file : file with individual element overrides (Z  log_abund)
  !     czc=0|1    : deep-CZ temperature constructor off/on (default on)
  !     smooth=0|1 : interior 1-2-1 FLXCNV smoothing off/on (default on)
  !     rosstab=N  : Rosseland-table interpolation 1|2|3 (default 1)
  !     quad=0|1   : INTEG quadrature, 0 = legacy blended-parabola,
  !                  1 = Steffen monotone cubic (default 0)
  
  OUTBASE    = 'mystar'
  ABUND_FILE = ''
  CMD_SOLAR  = ''
  INPUT_MODEL_FILE = ''
  NUMITS     = 30
  VTURB_KMS  = -1.0D0    ! sentinel: not set
  CMD_TEFF   = -1.0D0    ! sentinel: not set
  CMD_LOGG   = -99.0D0   ! sentinel: not set
  CMD_ZSCALE = -1.0D0    ! sentinel: not set
  CMD_HEABND = -1.0D0    ! sentinel: not set

  NARGS = COMMAND_ARGUMENT_COUNT()

  ! First pass: scan for help triggers anywhere on the command line
  DO I = 1, NARGS
    CALL GET_COMMAND_ARGUMENT(I, ARGBUF)
    SELECT CASE (TRIM(ARGBUF))
    CASE ('--help', '-h', 'help')
      CALL print_usage()
      CALL EXIT(0)
    END SELECT
  END DO

  ! No args at all -> show usage and exit with error
  IF (NARGS .EQ. 0) THEN
    CALL print_usage()
    CALL EXIT(1)
  END IF

  ! Sequential parse: positional args fill input_atm then basename, in order
  IPOSARG = 0
  DO I = 1, NARGS
    CALL GET_COMMAND_ARGUMENT(I, ARGBUF)
    IEQPOS = INDEX(ARGBUF, '=')
    IF (IEQPOS .EQ. 0) THEN
      ! positional argument
      IPOSARG = IPOSARG + 1
      SELECT CASE (IPOSARG)
      CASE (1)
        INPUT_MODEL_FILE = ARGBUF
      CASE (2)
        OUTBASE = ARGBUF
      CASE DEFAULT
        WRITE(6, '(A)') ' WARNING: extra positional argument ignored: '//TRIM(ARGBUF)
      END SELECT
      CYCLE
    END IF

    ! key=value: split once, then dispatch on key
    BLOCK
      CHARACTER(LEN=32)  :: key
      CHARACTER(LEN=256) :: val
      key = ARGBUF(1:IEQPOS-1)
      val = ARGBUF(IEQPOS+1:)
      ISTAT = 0
      SELECT CASE (TRIM(key))
      CASE ('numit');  READ(val, *, IOSTAT=ISTAT) NUMITS
      CASE ('vturb');  READ(val, *, IOSTAT=ISTAT) VTURB_KMS
      CASE ('czc');    READ(val, *, IOSTAT=ISTAT) ICZC_CLI
                       IF (ISTAT .EQ. 0) USE_CZ_CONSTRUCTOR = ICZC_CLI .NE. 0
      CASE ('smooth'); READ(val, *, IOSTAT=ISTAT) ICZC_CLI
                       IF (ISTAT .EQ. 0) USE_FLXCNV_SMOOTH = ICZC_CLI .NE. 0
      CASE ('rosstab'); READ(val, *, IOSTAT=ISTAT) IROSSTAB
                       IF (ISTAT .EQ. 0 .AND. (IROSSTAB .LT. 1 .OR. IROSSTAB .GT. 3)) THEN
                         WRITE(6, '(A)') ' ERROR: rosstab must be 1, 2, or 3'
                         CALL EXIT(1)
                       END IF
      CASE ('quad');   READ(val, *, IOSTAT=ISTAT) IQUAD
                       IF (ISTAT .EQ. 0 .AND. (IQUAD .LT. 0 .OR. IQUAD .GT. 1)) THEN
                         WRITE(6, '(A)') ' ERROR: quad must be 0 or 1'
                         CALL EXIT(1)
                       END IF
      CASE ('mlt');    READ(val, *, IOSTAT=ISTAT) MIXLTH
      CASE ('teff');   READ(val, *, IOSTAT=ISTAT) CMD_TEFF
      CASE ('logg');   READ(val, *, IOSTAT=ISTAT) CMD_LOGG
      CASE ('zscale'); READ(val, *, IOSTAT=ISTAT) CMD_ZSCALE
      CASE ('heabnd'); READ(val, *, IOSTAT=ISTAT) CMD_HEABND
      CASE ('solar')
        ! Normalize to lowercase, then validate the name now (fail fast).
        ! This also pre-loads ABUND, but READIN may overwrite it from the
        ! model file, so the selection is re-applied after READIN below.
        DO J = 1, LEN_TRIM(val)
          IF (val(J:J) .GE. 'A' .AND. val(J:J) .LE. 'Z') &
            val(J:J) = CHAR(ICHAR(val(J:J)) + 32)
        END DO
        CMD_SOLAR = TRIM(val)
        CALL SET_SOLAR_ABUND(TRIM(CMD_SOLAR), ISTAT)
        IF (ISTAT .NE. 0) THEN
          WRITE(6, '(A)') ' ERROR: unknown solar scale: '//TRIM(CMD_SOLAR)
          WRITE(6, '(A,*(1X,A))') '   valid names:', &
            (TRIM(SOLAR_SCALE_LIST(J)), J = 1, SIZE(SOLAR_SCALE_LIST))
          CALL EXIT(1)
        END IF
      CASE ('abund');  ABUND_FILE = TRIM(val)
      CASE DEFAULT
        WRITE(6, '(A)') ' WARNING: unknown argument: '//TRIM(ARGBUF)
      END SELECT
      IF (ISTAT .NE. 0) THEN
        WRITE(6, '(A,A,A)') ' ERROR: invalid ', TRIM(key), ' value: '//TRIM(ARGBUF)
        CALL EXIT(1)
      END IF
    END BLOCK
  END DO

  ! Validate: input atmosphere file is required
  IF (LEN_TRIM(INPUT_MODEL_FILE) .EQ. 0) THEN
    WRITE(6, '(A)') ' ERROR: missing required argument <input_atm>'
    WRITE(6, '(A)') ''
    CALL print_usage()
    CALL EXIT(1)
  END IF

  ! Validate: input atmosphere file must exist and be readable
  INQUIRE(FILE=TRIM(INPUT_MODEL_FILE), EXIST=file_exists)
  IF (.NOT. file_exists) THEN
    WRITE(6, '(A)') ' ERROR: input atmosphere file not found: '//TRIM(INPUT_MODEL_FILE)
    CALL EXIT(1)
  END IF

  ! Are any CLI overrides active?  (Controls whether we re-apply abundance
  ! math and/or recompute derived quantities after READIN.)
  have_abund_overrides = CMD_ZSCALE .GT. 0.0D0 .OR. CMD_HEABND .GT. 0.0D0 &
                         .OR. LEN_TRIM(ABUND_FILE) .GT. 0 &
                         .OR. LEN_TRIM(CMD_SOLAR) .GT. 0
  have_overrides       = have_abund_overrides .OR. VTURB_KMS .GE. 0.0D0 &
                         .OR. CMD_TEFF .GT. 0.0D0 .OR. CMD_LOGG .GT. -98.0D0

  ! Initialize IFPNCH: punch=2 only on last iteration
  ! Initialize IFPRNT: print=1 except last iteration=3
  DO I = 1, NUMITS
    IFPNCH(I) = 0
    IFPRNT(I) = 1
  END DO
  IFPNCH(NUMITS) = 2
  IFPRNT(NUMITS) = 3

  OPEN(UNIT=7,  FILE=TRIM(OUTBASE)//'.atm',   STATUS='REPLACE')
  OPEN(UNIT=8,  FILE=TRIM(OUTBASE)//'.flux',  STATUS='REPLACE')
  OPEN(UNIT=66, FILE=TRIM(OUTBASE)//'.iter',  STATUS='REPLACE')

  ! --- Locate data files via $ATLAS12 environment variable ---
  CALL GET_ENVIRONMENT_VARIABLE('ATLAS12', ENVVAL, ENVLEN, ENVSTAT)
  IF (ENVSTAT .EQ. 0 .AND. ENVLEN .GT. 0) THEN
    DATADIR = TRIM(ENVVAL)
    IF (DATADIR(ENVLEN:ENVLEN) .NE. '/') DATADIR = TRIM(DATADIR) // '/'
    DATADIR = TRIM(DATADIR) // 'data/'
  ELSE
    DATADIR = 'data/'
  END IF

  ! Point the B&C partition function module at the same data directory.
  CALL set_bc_data_dir(TRIM(DATADIR))

  ! --- Pre-tabulate Voigt profile H(a,v) at 200 steps per Doppler width ---
  VSTEPS = 200.0
  CALL TABVOIGT(VSTEPS, 2001)

  ! --- Standard input unit ---
  INPUTDATA = 5

  ! --- Load Feautrier coefficient matrices from external data files ---
  CALL BLOCKJ
  CALL BLOCKH

  ! --- Load ionization potentials for all species ---
  CALL IONPOTS

  ! --- Open optional pre-computed line data file (unit 19) ---
  !     ERR branch: if file doesn't exist, just continue
  OPEN(UNIT=19, FILE=TRIM(DATADIR)//'nltelinobsat12.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', ERR=10)
10 continue
  ITEMP = 0

  ! --- Read and compute a single model -----------------------------------
  CALL READIN(1)

    ! --- ABUNDANCE OVERRIDES (applied after READIN reads model file) ----
    ! Order: (1) solar= replaces the reference pattern, (2) zscale
    ! shifts metals, (3) abund= file overrides individual elements,
    ! (4) H is recomputed as X = 1 - Y - Z.
    ! At this point ABUND(1:2) are number fractions and ABUND(3:99)
    ! are log10(number fraction), as set by READIN finalization.

    ! --- Step 1: Select solar reference pattern (replaces the table the
    !     model file carried; XRELATIVE offsets survive, so an [M/H]
    !     scaling is preserved across a change of zero-point) ---
    IF (LEN_TRIM(CMD_SOLAR) .GT. 0) THEN
      CALL SET_SOLAR_ABUND(TRIM(CMD_SOLAR), ISTAT)
    END IF

    ! --- Step 2: Apply metallicity scaling (shifts all metals) ---
    IF (CMD_ZSCALE .GT. 0.0D0) THEN
      DO IZ_ABUND = 3, 99
        XRELATIVE(IZ_ABUND) = LOG10(CMD_ZSCALE)
      END DO
    END IF

    ! --- Step 3: Read individual element overrides from file ---
    !     File format: one element per line, two columns:
    !       Z   log10(number_fraction)
    !     e.g.:  6  -3.52
    !           26  -4.54
    !     Lines starting with '#' or '!' are comments.
    IF (LEN_TRIM(ABUND_FILE) .GT. 0) THEN
      OPEN(UNIT=4, FILE=TRIM(ABUND_FILE), STATUS='OLD', ACTION='READ', &
           IOSTAT=IOS_ABUND)
      IF (IOS_ABUND .NE. 0) THEN
        WRITE(6, '(A,A)') ' ERROR: cannot open abundance file: ', &
              TRIM(ABUND_FILE)
        CALL EXIT(1)
      END IF
      DO
        READ(4, '(A)', IOSTAT=IOS_ABUND) ABUND_LINE
        IF (IOS_ABUND .NE. 0) EXIT
        ABUND_LINE = ADJUSTL(ABUND_LINE)
        IF (LEN_TRIM(ABUND_LINE) .EQ. 0) CYCLE
        IF (ABUND_LINE(1:1) .EQ. '#' .OR. ABUND_LINE(1:1) .EQ. '!') CYCLE
        READ(ABUND_LINE, *, IOSTAT=IOS_ABUND) IZ_ABUND, ABUND_VAL
        IF (IOS_ABUND .NE. 0) CYCLE
        IF (IZ_ABUND .LT. 1 .OR. IZ_ABUND .GT. 99) CYCLE
        ABUND(IZ_ABUND) = ABUND_VAL
        IF (IZ_ABUND .GT. 2) XRELATIVE(IZ_ABUND) = 0.0D0
      END DO
      CLOSE(UNIT=4)
    END IF

    ! --- Step 4: Apply He override and recompute H = 1 - Y - Z ---
    IF (have_abund_overrides) THEN
      ! Compute total metal number fraction Z
      Z_TOTAL = 0.0D0
      DO IZ_ABUND = 3, 99
        Z_TOTAL = Z_TOTAL + 10.0D0**(ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND))
      END DO
      ! Set He: from command line or keep model value
      IF (CMD_HEABND .GT. 0.0D0) ABUND(2) = CMD_HEABND
      ! Compute H = 1 - Y - Z
      ABUND(1) = 1.0D0 - ABUND(2) - Z_TOTAL
      IF (ABUND(1) .LT. 0.0D0) THEN
        WRITE(6, '(A,F10.6)') ' ERROR: H abundance is negative: ', ABUND(1)
        WRITE(6, '(A,F10.6,A,F10.6)') '   Y = ', ABUND(2), '  Z = ', Z_TOTAL
        CALL EXIT(1)
      END IF
      WRITE(6, '(A,F10.6,A,F10.6,A,F10.6)') &
        ' Number fractions: X(H)=', ABUND(1), '  Y(He)=', ABUND(2), &
        '  Z(metals)=', Z_TOTAL
      ! Convert to mass fractions for display
      BLOCK
        REAL(8) :: MU_ABN, X_MASS, Y_MASS, Z_MASS
        MU_ABN = ABUND(1)*ATMASS(1) + ABUND(2)*ATMASS(2)
        DO IZ_ABUND = 3, 99
          MU_ABN = MU_ABN + 10.0D0**(ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND)) &
                   * ATMASS(IZ_ABUND)
        END DO
        X_MASS = ABUND(1) * ATMASS(1) / MU_ABN
        Y_MASS = ABUND(2) * ATMASS(2) / MU_ABN
        Z_MASS = 1.0D0 - X_MASS - Y_MASS
        WRITE(6, '(A,F10.6,A,F10.6,A,F10.6)') &
          ' Mass fractions:   X(H)=', X_MASS, '  Y(He)=', Y_MASS, &
          '  Z(metals)=', Z_MASS
      END BLOCK

      ! Recompute abundance-dependent arrays
      DO J = 1, NRHOX
        DO IZ_ABUND = 3, 99
          XABUND(J,IZ_ABUND) = 10.0D0**(ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND))
        END DO
        XABUND(J,1) = ABUND(1)
        XABUND(J,2) = ABUND(2)
        WTMOLE(J) = 0.0D0
        DO IZ_ABUND = 1, 99
          WTMOLE(J) = WTMOLE(J) + XABUND(J,IZ_ABUND) * ATMASS(IZ_ABUND)
        END DO
      END DO
      YABUND(1) = ABUND(1)
      YABUND(2) = ABUND(2)
      DO IZ_ABUND = 3, 99
        YABUND(IZ_ABUND) = ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND)
      END DO
    END IF

    ! --- Regrid and rescale model to target Teff/logg if specified ---
    !     If neither teff= nor logg= given, use values from the model file.
    !     If only one is given, the other defaults to the model file value.
    IF (CMD_TEFF .GT. 0.0D0 .OR. CMD_LOGG .GT. -98.0D0) THEN
      IF (CMD_TEFF .LT. 0.0D0) CMD_TEFF = TEFF
      IF (CMD_LOGG .LT. -98.0D0) CMD_LOGG = GLOG
      CALL SCALE_MODEL(CMD_TEFF, CMD_LOGG)
    END IF

    ! --- Apply command-line VTURB override (km/s → cm/s) if specified ---
    IF (VTURB_KMS .GE. 0.0D0) THEN
      VTURB = VTURB_KMS * 1.0D5
    END IF

    ! --- Recompute derived quantities after any overrides ---
    !     READIN's finalization computed these from the original model;
    !     if we changed abundances, regridded, or changed VTURB, they need updating.
    IF (have_overrides) THEN
        TK   = KBOL * T
        HKT  = HPLANCK / TK
        HCKT = HKT * CLIGHT
        TKEV = KBOL_EV * T
        TLOG = LOG(T)
        XNATOM = P / TK - XNE
        RHO  = XNATOM * WTMOLE * AMU
        IF (IFTURB .GT. 0) PTURB = 0.5D0 * RHO * VTURB**2
        ! Approximate CHARGESQ = sum(n_i * Z_i^2) + n_e for Debye shielding.
        ! Assuming mostly singly ionized gas, the ion sum ~ n_e, so
        ! CHARGESQ ~ 2*n_e.  NELECT recomputes this self-consistently.
        CHARGESQ = XNE * 2.0D0
    END IF

    ! --- Echo run parameters (resolved values after CLI/model merge) ---
    ! Mirrors SYNTHE's input-echo style.  All values are final at this
    ! point; the (CLI override) tag flags parameters whose value came from
    ! the command line rather than the model file.
    WRITE(6,*) ''
    WRITE(6,'(A,A)')     '  Input model      = ', TRIM(INPUT_MODEL_FILE)
    WRITE(6,'(A,A)')     '  Output basename  = ', TRIM(OUTBASE)
    WRITE(6,'(A,I5)')    '  numit            = ', NUMITS
    WRITE(6,'(A,F5.2)')  '  mlt              = ', MIXLTH
    WRITE(6,'(A,F5.2)')  '  vturb (km/s)     = ', VTURB(1) * 1.0D-5
    WRITE(6,'(A,I5)')    '  teff (K)         = ', INT(TEFF)
    WRITE(6,'(A,F5.2)')  '  logg             = ', GLOG
    IF (LEN_TRIM(CMD_SOLAR) .GT. 0) &
      WRITE(6,'(A,A,A)')    '  solar scale      = ', TRIM(CMD_SOLAR), '   (CLI override)'
    IF (IQUAD .NE. 0) &
      WRITE(6,'(A,I2,A)')  '  quad             = ', IQUAD, '   (CLI override)'
    IF (CMD_ZSCALE .GT. 0.0D0) &
      WRITE(6,'(A,F5.2,A)') '  zscale           = ', CMD_ZSCALE, '   (CLI override)'
    IF (CMD_HEABND .GT. 0.0D0) &
      WRITE(6,'(A,F5.2,A)') '  He abundance     = ', CMD_HEABND, '   (CLI override)'
    IF (LEN_TRIM(ABUND_FILE) .GT. 0) &
      WRITE(6,'(A,A)')       '  Abundance file   = ', TRIM(ABUND_FILE)
    WRITE(6,*) ''

    ! --- Zero out number densities and Doppler widths for all species ---
    XNF    = 0.0D0
    XNFP   = 0.0D0
    DOPPLE = 0.0D0

    ! --- Load isotope mass fractions ---
    CALL ISOTOPES

    ! --- Find dominant isotope mass for each ion species ---
    !     AMASSISO(1,NELION) = mass of the most abundant isotope
    DO NELION = 1, MION
      XMAX = 0.0D0
      DO I = 1, 10
        IF (ISOTOPE(I, 2, NELION) .GT. XMAX) THEN
          AMASSISO(1, NELION) = ISOTOPE(I, 1, NELION)
          XMAX = ISOTOPE(I, 2, NELION)
        END IF
      END DO
      ! Fallback: if no isotope fractions given, use first listed mass
      IF (XMAX .LE. 0.0D0 .AND. ISOTOPE(1, 1, NELION) .GT. 0.0) THEN
        AMASSISO(1, NELION) = ISOTOPE(1, 1, NELION)
      END IF
    END DO

    ! --- Write depth/Teff summary to unit 66 ---
    WRITE(66, '(I3,1X,F8.1)') NRHOX, TEFF

    ! --- SET UP WAVELENGTH GRID -----------------------------------------
    !   The wavelength grid is logarithmically spaced:
    !     WAVESET(NU) = 10^(1.0 + 0.0001*(NU + NUSTART - 1))  [nm]
    !   NUSTART shifts the grid blueward for hotter stars to capture
    !   the ionization edges of H, He I, He II.

    NULO   = 1
    NUHI   = 30000
    NUSTEP = 1

    ! Frequency indices for key ionization edges
    NUCI     = 11601    ! C I ionization edge
    NULYMAN  = 9599     ! H Lyman limit (91.2 nm)
    NUHEI    = 7027     ! He I ionization edge (50.4 nm)
    NUHEII   = 3577     ! He II ionization edge (22.8 nm)

    ! Select starting wavelength index based on Teff
    NUSTART = 1
    IF (TEFF .LT. 30000.0D0) NUSTART = NUHEII
    IF (TEFF .LT. 13000.0D0) NUSTART = NUHEI
    IF (TEFF .LT. 7250.0D0)  NUSTART = NULYMAN
    IF (TEFF .LT. 4500.0D0)  NUSTART = NUCI

    ! Build the wavelength array
    NUMNU = NUHI
    DO NU = NULO, NUHI, NUSTEP
      WAVESET(NU) = 10.0D0 ** (1.0D0 + 0.0001D0 * (NU + NUSTART - 1))
    END DO

    ! --- COMPUTE FREQUENCY INTEGRATION COEFFICIENTS (trapezoidal rule) ---
    !   RCOSET(NU) = integration weight in frequency space for each
    !   wavelength point, used in flux/correction integrals.
    !   Boundary conditions: flux = 0 at blue and red limits.

    ! Blue edge: assume flux = 0 at WAVESET(0)
    RCOSET(1) = (CLIGHT_NMS / WAVESET(1) &
               - CLIGHT_NMS / WAVESET(1 + NUSTEP)) * 1.5D0

    ! Interior points: centered differences
    DO NU = NULO + NUSTEP, NUMNU - NUSTEP, NUSTEP
      RCOSET(NU) = (CLIGHT_NMS / WAVESET(NU - NUSTEP) &
                  -  CLIGHT_NMS / WAVESET(NU + NUSTEP)) * 0.5D0
    END DO

    ! Red edge: assume flux = 0 at infinite wavelength (freq = 0)
    RCOSET(NUMNU) = (CLIGHT_NMS / WAVESET(NUMNU - NUSTEP) &
                   + CLIGHT_NMS / WAVESET(NUMNU)) * 0.25D0


    ! --- ITERATION LOOP — converge the atmospheric structure ------------

    iteration_loop: DO ITERAT = 1, NUMITS
      ITER = ITERAT
      CALL SYSTEM_CLOCK(CLOCK_START, CLOCK_RATE)

      ! ITEMP tracks temperature changes — incrementing tells subroutines
      ! that T has been updated and they should recompute T-dependent quantities
      ITEMP = ITEMP + ITER

      ! ---------------------------------------------------------------
      ! HYDROSTATIC EQUILIBRIUM
      ! ---------------------------------------------------------------
      IF (IFPRES .NE. 0) THEN

         IF (ITEMP .NE. 1) THEN
            
          ! Integrate the equation of hydrostatic equilibrium:
          !   P_gas = g * RHOX - P_rad - P_turb - P_con
          DO J = 1, NRHOX
            P(J) = GRAV * RHOX(J) - PRAD(J) - PTURB(J) - PCON
            IF (P(J) .LE. 0.0D0) THEN
               P(J) = MAX(GRAV * RHOX(J) * 1.0D-4, 1.0D-10)
            END IF
          END DO
        END IF

        ! Boundary pressure: P0 = P_continuous + P_rad(surface) + P_turb(surface)
        PZERO = PCON + PRADK0 + PTURB0

        ! Compute charge-squared sum and total pressure at each depth
        DO J = 1, NRHOX
          CHARGESQ(J) = XNE(J) * 2.0D0
          EXCESS = 2.0D0 * XNE(J) - P(J) / TK(J)
          ! Allowance for doubly ionized helium
          IF (EXCESS .GT. 0.0D0) CHARGESQ(J) = CHARGESQ(J) + 2.0D0 * EXCESS
          PTOTAL(J) = GRAV * RHOX(J) + PZERO
        END DO

        ! --- Compute populations ---
        IFEDNS = 0
        CALL COMPUTE_ONE_POP(0.D0, 1, XNE)
        CALL COMPUTE_ALL_POPS

      END IF  ! IFPRES

      ! --- Radiation energy density (if needed) ---
      IF (IFEDNS .EQ. 1) CALL ENERGY_DENSITY

      ! ---------------------------------------------------------------
      ! DOPPLER WIDTHS for all species at all depths
      ! ---------------------------------------------------------------
      !   Doppler width = sqrt(2kT/m + v_turb^2) / c
      !   XNFDOP = (number density / partition function) / (Doppler width * density)
      !          = line opacity coefficient per unit Doppler width

      DO J = 1, NRHOX
        DO NELION = 1, MION - 1
          IF (AMASSISO(1, NELION) .LE. 0.0D0) CYCLE
          DOPPLE(J, NELION) = SQRT(2.0D0 * TK(J) / AMASSISO(1, NELION) / AMU &
                                   + VTURB(J)**2) / CLIGHT
          XNFDOP(J, NELION) = XNFP(J, NELION) / DOPPLE(J, NELION) / RHO(J)
        END DO
      END DO

      ! ---------------------------------------------------------------
      ! OPACITY TABLES (first iteration of first model only)
      ! ---------------------------------------------------------------
      IF (ITER .EQ. 1 .AND. ITEMP .EQ. 1) CALL KAPCONT
      IF (ITER .EQ. 1 .AND. IFREADLINES .EQ. 1 .AND. ITEMP .EQ. 1) CALL SELECTLINES

      ! ---------------------------------------------------------------
      ! RADIATIVE TRANSFER PASS  (line opacities + frequency loop)
      ! ---------------------------------------------------------------
      !   Loops over all wavelength points.  At each frequency:
      !     1. Planck function and stimulated emission correction
      !     2. Continuous opacity (KAPP)
      !     3. OS line opacity (XLINES)
      !     4. Transfer solve (JOSH)
      !     5. Accumulate flux moments for T-correction etc.
      !   Factored out to atlas12_modules.RUN_RT_PASS so the same loop
      !   can be re-run with a perturbed T to build a numerical Jacobian
      !   for the Stage-2 hybrid Newton temperature correction.
      CALL RUN_RT_PASS(.TRUE.)

      ! For surface-only mode, skip the iteration finishing steps
      IF (IFSURF .LE. 0) THEN

        ! ---------------------------------------------------------------
        ! FINISH ITERATION — Rosseland mean, convection, T-correction
        ! ---------------------------------------------------------------
        CALL ROSS(3, 0.0D0)
        RX = ROSSTAB(0.D0, 0.D0, 0.D0)
        CALL RADIAP(3, 0.0D0)
        CALL COMPUTE_HEIGHT
        IF (IFPRES .EQ. 1 .AND. IFCONV .EQ. 1) CALL CONVEC(.FALSE.)
        IF (IFCORR .EQ. 1)  CALL TCORR(3, 0.0D0)
        IF (NLTEON .EQ. 1)  CALL STATEQ(3, 0.0D0)
        IF (IFTURB .EQ. 1)  CALL COMPUTE_PTURB
        CALL PUTOUT(5)

        CALL SYSTEM_CLOCK(CLOCK_END)
        WRITE(6, '(A,I2,A,I3,A)') ' Iteration ', ITERAT, ' completed in ', &
             INT(DBLE(CLOCK_END - CLOCK_START) / DBLE(CLOCK_RATE)), ' seconds'
        FLUSH(6)

      END IF

    END DO iteration_loop

  CALL SYSTEM_CLOCK(CLOCK_END)
  CALL report_elapsed(REAL(CLOCK_END - CLOCK_TOTAL_START, 8) / REAL(CLOCK_RATE, 8))

CONTAINS

  ! ------------------------------------------------------------------------
  !  report_elapsed(seconds)
  !
  !  Print wall-clock elapsed time in a human-readable form:
  !    < 1 min : SS.SSs
  !    < 1 hr  : Mm SSs
  !    >= 1 hr : Hh MMm SSs
  ! ------------------------------------------------------------------------
  SUBROUTINE report_elapsed(seconds)
    REAL(8), INTENT(IN) :: seconds
    INTEGER :: h, m, s

    IF (seconds .LT. 60.0D0) THEN
      WRITE(6,'(/A,F6.2,A)') ' Elapsed: ', seconds, 's'
    ELSE IF (seconds .LT. 3600.0D0) THEN
      m = INT(seconds / 60.0D0)
      s = NINT(seconds - 60.0D0 * m)
      WRITE(6,'(/A,I0,A,I2.2,A)') ' Elapsed: ', m, 'm ', s, 's'
    ELSE
      h = INT(seconds / 3600.0D0)
      m = INT((seconds - 3600.0D0 * h) / 60.0D0)
      s = NINT(seconds - 3600.0D0 * h - 60.0D0 * m)
      WRITE(6,'(/A,I0,A,I2.2,A,I2.2,A)') ' Elapsed: ', h, 'h ', m, 'm ', s, 's'
    END IF

  END SUBROUTINE report_elapsed

  ! ------------------------------------------------------------------------
  !  print_usage()
  !
  !  Print command-line usage to stderr-equivalent (unit 6).  Called when
  !  the user passes --help/-h/help, when no arguments are given, or when
  !  the required <input_atm> argument is missing.
  ! ------------------------------------------------------------------------
  SUBROUTINE print_usage()
    WRITE(6, '(A)') 'Usage: atlas12.exe <input_atm> [basename] [key=value ...]'
    WRITE(6, '(A)') ''
    WRITE(6, '(A)') 'Required:'
    WRITE(6, '(A)') '  input_atm    Input atmosphere model file'
    WRITE(6, '(A)') ''
    WRITE(6, '(A)') 'Optional positional:'
    WRITE(6, '(A)') '  basename     Output file base name (default ''mystar'')'
    WRITE(6, '(A)') '               Outputs: <basename>.atm, .flux, .iter'
    WRITE(6, '(A)') ''
    WRITE(6, '(A)') 'Options (key=value):'
    WRITE(6, '(A)') '  numit=N      Number of iterations (default 30)'
    WRITE(6, '(A)') '  vturb=X      Microturbulence in km/s (default: from model)'
    WRITE(6, '(A)') '  mlt=X        Mixing length parameter (default 2.0)'
    WRITE(6, '(A)') '  teff=X       Rescale model to this Teff (default: from model)'
    WRITE(6, '(A)') '  logg=X       Rescale model to this log g (default: from model)'
    WRITE(6, '(A)') '  solar=NAME   Solar reference abundance pattern: ag89, agss09, berg25'
    WRITE(6, '(A)') '               (default: from model file; replaces its table)'
    WRITE(6, '(A)') '  zscale=X     Metal abundance scale factor (default: no scaling)'
    WRITE(6, '(A)') '  heabnd=X     He number fraction Y; H = 1 - Y - Z (default: from model)'
    WRITE(6, '(A)') '  abund=file   File with individual element overrides (Z log_abund)'
    WRITE(6, '(A)') '  czc=0|1      Deep-CZ temperature constructor off/on (default on)'
    WRITE(6, '(A)') '  smooth=0|1   Interior 1-2-1 FLXCNV smoothing off/on (default on)'
    WRITE(6, '(A)') '  rosstab=N    Rosseland-table interpolation: 1=bilinear 2=Shepard 3=MLS'
    WRITE(6, '(A)') '  quad=0|1     INTEG quadrature: 0=blended-parabola 1=Steffen cubic (default 0)'
    WRITE(6, '(A)') ''
    WRITE(6, '(A)') 'Help:'
    WRITE(6, '(A)') '  --help, -h, help    Print this message and exit'
  END SUBROUTINE print_usage

END PROGRAM ATLAS12
