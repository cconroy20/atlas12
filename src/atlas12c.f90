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
!   Change log (selected):
!     17sep2015  Convective flux smoothing enabled
!     16sep2015  Enforcing temperature minimum of 0.75*Teff
!     27jan2015  Reformatted fort.8 output (surface flux) — C. Conroy
!     21jan2015  Changed to 80 depths — C. Conroy
!     04nov2009  Bugs in COMPUTE_ALL_POPS for Be, B, Ar — M. Stift
!     22oct2009  Oxygen bug in PFSAHA — M. Stift
!     05aug2009  van der Waals broadening by H2 in HPROF
!     See source header for full history.
!=========================================================================

PROGRAM ATLAS12

  use mod_atlas_data
  implicit none

  ! --- Physical constants (CGS) ---
  !   k_B  = 1.38054e-16  erg/K
  !   h    = 6.6256e-27   erg*s
  !   c    = 2.99792458e10 cm/s
  !   e    = 1.60210e-19   C  (not used directly here)
  !   amu  = 1.660e-24     g

  ! --- Named constants for clarity ---
  real*8, parameter  :: CLIGHT_ANGFREQ = 2.99792458D17  ! c in Angstrom/s (= c * 1e8)
  real*8, parameter  :: AMU_GRAMS      = 1.660D-24      ! atomic mass unit in grams
  real*8, parameter  :: CLIGHT_CGS     = 2.99792458D10  ! c in cm/s
  real*8, parameter  :: PLANCK_PREFAC  = 1.47439D-2     ! 2*h/c^2 in CGS-frequency units

  ! --- Local variables ---
  real*8             :: VSTEPS, RCOSUM, RCOWT
  real*8             :: EXCESS, XMAX, FREQ15, DOPPLE8
  real*8             :: RX
  real*4             :: stim4(kw)
  integer            :: I, J, NU, ITERAT
  integer            :: NUCI, NULYMAN, NUHEI, NUHEII, NUSTART
  logical            :: found_negative
  integer*8          :: CLOCK_START, CLOCK_END, CLOCK_RATE

  ! --- Environment variable for data directory ---
  character(256)     :: ENVVAL
  integer            :: ENVLEN, ENVSTAT

  ! --- Command-line output file base name and options ---
  character(256)     :: OUTBASE, ARGBUF
  character(256)     :: ABUND_FILE
  integer            :: NARGS, IEQPOS, ISTAT
  real*8             :: VTURB_KMS
  real*8             :: CMD_TEFF, CMD_LOGG
  real*8             :: CMD_ZSCALE, CMD_HEABND
  real*8             :: Z_TOTAL
  integer            :: IZ_ABUND, IOS_ABUND
  real*8             :: ABUND_VAL

  ! =====================================================================
  ! INITIALIZATION PHASE
  ! =====================================================================

  ! --- Parse command-line arguments ---
  !     Usage: atlas12.exe [basename] [numit=N] [vturb=X] [mlt=X] [teff=X]
  !            [logg=X] [zscale=X] [heabnd=X] [abund=file]
  !     Output files: <basename>.atm, .flux, .taunu, .iter, .tcorr
  !     basename   : output file base name (default 'mystar')
  !     numit=N    : number of iterations (default 30)
  !     vturb=X    : microturbulence in km/s (default: from model)
  !     mlt=X      : mixing length parameter (default 2.0)
  !     teff=X     : rescale model to this Teff (default: from model)
  !     logg=X     : rescale model to this log g (default: from model)
  !     zscale=X   : metal abundance scale factor (default: no scaling)
  !     heabnd=X   : He number fraction Y (default: from model);
  !                  H is computed as X = 1 - Y - Z for consistency
  !     abund=file : file with individual element overrides (Z  log_abund)
  OUTBASE = 'mystar'
  ABUND_FILE = ''
  NUMITS = 30
  VTURB_KMS  = -1.0D0   ! sentinel: not set
  CMD_TEFF   = -1.0D0   ! sentinel: not set
  CMD_LOGG   = -99.0D0  ! sentinel: not set
  CMD_ZSCALE = -1.0D0   ! sentinel: not set
  CMD_HEABND = -1.0D0   ! sentinel: not set
  NARGS = command_argument_count()
  do I = 1, NARGS
    call get_command_argument(I, ARGBUF)
    IEQPOS = INDEX(ARGBUF, '=')
    if (IEQPOS > 0) then
      ! keyword=value argument
      if (ARGBUF(1:IEQPOS) == 'numit=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) NUMITS
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid numit value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'vturb=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) VTURB_KMS
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid vturb value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'mlt=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) MIXLTH
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid mlt value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'teff=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) CMD_TEFF
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid teff value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'logg=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) CMD_LOGG
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid logg value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'zscale=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) CMD_ZSCALE
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid zscale value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'heabnd=') then
        read(ARGBUF(IEQPOS+1:), *, IOSTAT=ISTAT) CMD_HEABND
        if (ISTAT /= 0) then
          write(6, '(A)') ' ERROR: invalid heabnd value: '//trim(ARGBUF)
          stop 1
        end if
      else if (ARGBUF(1:IEQPOS) == 'abund=') then
        ABUND_FILE = ARGBUF(IEQPOS+1:)
      else
        write(6, '(A)') ' WARNING: unknown argument: '//trim(ARGBUF)
      end if
    else
      ! positional argument = basename
      OUTBASE = ARGBUF
    end if
  end do
  ! Initialize IFPNCH: punch=2 only on last iteration
  ! Initialize IFPRNT: print=1 except last iteration=3
  do I = 1, NUMITS
    IFPNCH(I) = 0
    IFPRNT(I) = 1
  end do
  IFPNCH(NUMITS) = 2
  IFPRNT(NUMITS) = 3

  open(unit=7,  file=trim(OUTBASE)//'.atm',   status='REPLACE')
  open(unit=8,  file=trim(OUTBASE)//'.flux',  status='REPLACE')
  open(unit=50, file=trim(OUTBASE)//'.taunu', status='REPLACE')
  open(unit=66, file=trim(OUTBASE)//'.iter',  status='REPLACE')
  open(unit=67, file=trim(OUTBASE)//'.tcorr', status='REPLACE')

  ! --- Locate data files via $ATLAS12 environment variable ---
  call GET_ENVIRONMENT_VARIABLE('ATLAS12', ENVVAL, ENVLEN, ENVSTAT)
  if (ENVSTAT == 0 .and. ENVLEN > 0) then
    DATADIR = trim(ENVVAL)
    if (DATADIR(ENVLEN:ENVLEN) /= '/') DATADIR = trim(DATADIR) // '/'
    DATADIR = TRIM(DATADIR) // 'data/'
  else
    DATADIR = 'data/'
  end if

  ! --- Pre-tabulate Voigt profile H(a,v) at 200 steps per Doppler width ---
  VSTEPS = 200.0
  call TABVOIGT(VSTEPS, 2001)

  ! --- Standard input unit ---
  INPUTDATA = 5

  ! --- Load Feautrier coefficient matrices from external data files ---
  call BLOCKJ
  call BLOCKH

  ! --- Load ionization potentials for all species ---
  call IONPOTS

  ! --- Open optional pre-computed line data file (unit 19) ---
  !     ERR branch: if file doesn't exist, just continue
  open(UNIT=19, FILE=trim(DATADIR)//'nltelines_obs.bin', &
       STATUS='OLD', FORM='UNFORMATTED', ACTION='READ', ERR=10)
10 continue
  ITEMP = 0

  ! =====================================================================
  ! =====================================================================
  ! Read and compute a single model
  ! =====================================================================

  call READIN(1)

    ! =================================================================
    ! ABUNDANCE OVERRIDES (applied after READIN reads model file)
    !
    ! Order: (1) zscale shifts metals, (2) abund= file overrides
    ! individual elements, (3) H is recomputed as X = 1 - Y - Z.
    ! At this point ABUND(1:2) are number fractions and ABUND(3:99)
    ! are log10(number fraction), as set by READIN finalization.
    ! =================================================================

    ! --- Step 1: Apply metallicity scaling (shifts all metals) ---
    if (CMD_ZSCALE > 0.0D0) then
      do IZ_ABUND = 3, 99
        XRELATIVE(IZ_ABUND) = LOG10(CMD_ZSCALE)
      end do
    end if

    ! --- Step 2: Read individual element overrides from file ---
    !     File format: one element per line, two columns:
    !       Z   log10(number_fraction)
    !     e.g.:  6  -3.52
    !           26  -4.54
    !     Lines starting with '#' or '!' are comments.
    if (len_trim(ABUND_FILE) > 0) then
      open(UNIT=4, FILE=trim(ABUND_FILE), STATUS='OLD', ACTION='READ', &
           IOSTAT=IOS_ABUND)
      if (IOS_ABUND /= 0) then
        write(6, '(A,A)') ' ERROR: cannot open abundance file: ', &
              trim(ABUND_FILE)
        stop 1
      end if
      do
        read(4, *, IOSTAT=IOS_ABUND) IZ_ABUND, ABUND_VAL
        if (IOS_ABUND /= 0) exit
        if (IZ_ABUND < 1 .or. IZ_ABUND > 99) cycle
        ABUND(IZ_ABUND) = ABUND_VAL
        if (IZ_ABUND > 2) XRELATIVE(IZ_ABUND) = 0.0D0
      end do
      close(UNIT=4)
    end if

    ! --- Step 3: Apply He override and recompute H = 1 - Y - Z ---
    if (CMD_ZSCALE > 0.0D0 .or. CMD_HEABND > 0.0D0 .or. &
        len_trim(ABUND_FILE) > 0) then
      ! Compute total metal number fraction Z
      Z_TOTAL = 0.0D0
      do IZ_ABUND = 3, 99
        Z_TOTAL = Z_TOTAL + 10.0D0**(ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND))
      end do
      ! Set He: from command line or keep model value
      if (CMD_HEABND > 0.0D0) ABUND(2) = CMD_HEABND
      ! Compute H = 1 - Y - Z
      ABUND(1) = 1.0D0 - ABUND(2) - Z_TOTAL
      if (ABUND(1) < 0.0D0) then
        write(6, '(A,F10.6)') ' ERROR: H abundance is negative: ', ABUND(1)
        write(6, '(A,F10.6,A,F10.6)') '   Y = ', ABUND(2), '  Z = ', Z_TOTAL
        stop 1
      end if
      write(6, '(A,F10.6,A,F10.6,A,F10.6)') &
        ' Number fractions: X(H)=', ABUND(1), '  Y(He)=', ABUND(2), &
        '  Z(metals)=', Z_TOTAL
      ! Convert to mass fractions for display
      block
        real*8 :: MU_ABN, X_MASS, Y_MASS, Z_MASS
        MU_ABN = ABUND(1)*ATMASS(1) + ABUND(2)*ATMASS(2)
        do IZ_ABUND = 3, 99
          MU_ABN = MU_ABN + 10.0D0**(ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND)) &
                   * ATMASS(IZ_ABUND)
        end do
        X_MASS = ABUND(1) * ATMASS(1) / MU_ABN
        Y_MASS = ABUND(2) * ATMASS(2) / MU_ABN
        Z_MASS = 1.0D0 - X_MASS - Y_MASS
        write(6, '(A,F10.6,A,F10.6,A,F10.6)') &
          ' Mass fractions:   X(H)=', X_MASS, '  Y(He)=', Y_MASS, &
          '  Z(metals)=', Z_MASS
      end block

      ! Recompute abundance-dependent arrays
      do J = 1, NRHOX
        do IZ_ABUND = 3, 99
          XABUND(J,IZ_ABUND) = 10.0D0**(ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND))
        end do
        XABUND(J,1) = ABUND(1)
        XABUND(J,2) = ABUND(2)
        WTMOLE(J) = 0.0D0
        do IZ_ABUND = 1, 99
          WTMOLE(J) = WTMOLE(J) + XABUND(J,IZ_ABUND) * ATMASS(IZ_ABUND)
        end do
      end do
      YABUND(1) = ABUND(1)
      YABUND(2) = ABUND(2)
      do IZ_ABUND = 3, 99
        YABUND(IZ_ABUND) = ABUND(IZ_ABUND) + XRELATIVE(IZ_ABUND)
      end do
    end if

    ! --- Regrid and rescale model to target Teff/logg if specified ---
    !     If neither teff= nor logg= given, use values from the model file.
    !     If only one is given, the other defaults to the model file value.
    if (CMD_TEFF > 0.0D0 .or. CMD_LOGG > -98.0D0) then
      if (CMD_TEFF < 0.0D0) CMD_TEFF = TEFF
      if (CMD_LOGG < -98.0D0) CMD_LOGG = GLOG
      call SCALE_MODEL(CMD_TEFF, CMD_LOGG)
    end if

    ! --- Apply command-line VTURB override (km/s → cm/s) if specified ---
    if (VTURB_KMS >= 0.0D0) then
      do J = 1, NRHOX
        VTURB(J) = VTURB_KMS * 1.0D5
      end do
    end if

    ! --- Recompute derived quantities after any overrides ---
    !     READIN's finalization computed these from the original model;
    !     if we changed abundances, regridded, or changed VTURB, they need updating.
    if (CMD_TEFF > 0.0D0 .or. CMD_LOGG > -98.0D0 .or. VTURB_KMS >= 0.0D0 .or. &
        CMD_ZSCALE > 0.0D0 .or. CMD_HEABND > 0.0D0 .or. &
        len_trim(ABUND_FILE) > 0) then
      do J = 1, NRHOX
        TK(J) = 1.38054D-16 * T(J)
        HKT(J) = 6.6256D-27 / TK(J)
        HCKT(J) = HKT(J) * 2.99792458D10
        TKEV(J) = 8.6171D-5 * T(J)
        TLOG(J) = LOG(T(J))
        XNATOM(J) = P(J) / TK(J) - XNE(J)
        RHO(J) = XNATOM(J) * WTMOLE(J) * 1.660D-24
        if (IFTURB > 0) PTURB(J) = 0.5D0 * RHO(J) * VTURB(J)**2
        CHARGESQ(J) = XNE(J)
      end do
    end if

    ! --- Zero out number densities and Doppler widths for all species ---
    do NELION = 1, MION
      do J = 1, NRHOX
        XNF(J, NELION)  = 0.0D0
        XNFP(J, NELION) = 0.0D0
        DOPPLE(J, NELION) = 0.0
      end do
    end do

    ! --- Load isotope mass fractions ---
    call ISOTOPES

    ! --- Find dominant isotope mass for each ion species ---
    !     AMASSISO(1,NELION) = mass of the most abundant isotope
    do NELION = 1, MION
      XMAX = 0.0D0
      do I = 1, 10
        if (ISOTOPE(I, 2, NELION) > XMAX) then
          AMASSISO(1, NELION) = ISOTOPE(I, 1, NELION)
          XMAX = ISOTOPE(I, 2, NELION)
        end if
      end do
      ! Fallback: if no isotope fractions given, use first listed mass
      if (XMAX == 0.0 .and. ISOTOPE(1, 1, NELION) > 0.0) then
        AMASSISO(1, NELION) = ISOTOPE(1, 1, NELION)
      end if
    end do

    ! --- Write depth/Teff summary to unit 66 ---
    write(66, '(I3,1X,F8.1)') NRHOX, TEFF

    ! =================================================================
    ! SET UP WAVELENGTH GRID
    ! =================================================================
    !   The wavelength grid is logarithmically spaced:
    !     WAVESET(NU) = 10^(1.0 + 0.0001*(NU + NUSTART - 1))  [Angstroms]
    !   NUSTART shifts the grid blueward for hotter stars to capture
    !   the ionization edges of H, He I, He II.

    NULO   = 1
    NUHI   = 30000
    NUSTEP = 1

    ! Frequency indices for key ionization edges
    NUCI     = 11601    ! C I ionization edge
    NULYMAN  = 9599     ! H Lyman limit (912 A)
    NUHEI    = 7027     ! He I ionization edge (504 A)
    NUHEII   = 3577     ! He II ionization edge (228 A)

    ! Select starting wavelength index based on Teff
    NUSTART = 1
    if (TEFF < 30000.0D0) NUSTART = NUHEII
    if (TEFF < 13000.0D0) NUSTART = NUHEI
    if (TEFF < 7250.0D0)  NUSTART = NULYMAN
    if (TEFF < 4500.0D0)  NUSTART = NUCI

    ! Build the wavelength array
    do NU = NULO, NUHI, NUSTEP
      NUMNU = NU
      WAVESET(NU) = 10.0D0 ** (1.0D0 + 0.0001D0 * (NU + NUSTART - 1))
    end do

    ! =================================================================
    ! COMPUTE FREQUENCY INTEGRATION COEFFICIENTS (trapezoidal rule)
    ! =================================================================
    !   RCOSET(NU) = integration weight in frequency space for each
    !   wavelength point, used in flux/correction integrals.
    !   Boundary conditions: flux = 0 at blue and red limits.

    ! Blue edge: assume flux = 0 at WAVESET(0)
    RCOSET(1) = (CLIGHT_ANGFREQ / WAVESET(1) &
               - CLIGHT_ANGFREQ / WAVESET(1 + NUSTEP)) * 1.5D0
    RCOSUM = RCOSET(1)

    ! Interior points: centered differences
    do NU = NULO + NUSTEP, NUMNU - NUSTEP, NUSTEP
      RCOSET(NU) = (CLIGHT_ANGFREQ / WAVESET(NU - NUSTEP) &
                  -  CLIGHT_ANGFREQ / WAVESET(NU + NUSTEP)) * 0.5D0
      RCOSUM = RCOSUM + RCOSET(NU)
    end do

    ! Red edge: assume flux = 0 at infinite wavelength (freq = 0)
    RCOSET(NUMNU) = (CLIGHT_ANGFREQ / WAVESET(NUMNU - NUSTEP) &
                   + CLIGHT_ANGFREQ / WAVESET(NUMNU)) * 0.25D0
    RCOSUM = RCOSUM + RCOSET(NUMNU)


    ! =================================================================
    ! ITERATION LOOP — converge the atmospheric structure
    ! =================================================================

    iteration_loop: do ITERAT = 1, NUMITS
      ITER = ITERAT
      call system_clock(CLOCK_START, CLOCK_RATE)

      ! ITEMP tracks temperature changes — incrementing tells subroutines
      ! that T has been updated and they should recompute T-dependent quantities
      ITEMP = ITEMP + ITER

      ! ---------------------------------------------------------------
      ! HYDROSTATIC EQUILIBRIUM
      ! ---------------------------------------------------------------
      if (IFPRES /= 0) then

        if (ITEMP /= 1) then
          ! Integrate the equation of hydrostatic equilibrium:
          !   P_gas = g * RHOX - P_rad - P_turb - P_con
          do J = 1, NRHOX
            P(J) = GRAV * RHOX(J) - PRAD(J) - PTURB(J) - PCON
            if (P(J) <= 0.0D0) then
               P(J) = max(GRAV * RHOX(J) * 1.0D-4, 1.0D-10)
            end if
          end do
        end if

        ! Boundary pressure: P0 = P_continuous + P_rad(surface) + P_turb(surface)
        PZERO = PCON + PRADK0 + PTURB0

        ! Compute charge-squared sum and total pressure at each depth
        do J = 1, NRHOX
          CHARGESQ(J) = XNE(J) * 2.0D0
          EXCESS = 2.0D0 * XNE(J) - P(J) / TK(J)
          ! Allowance for doubly ionized helium
          if (EXCESS > 0.0D0) CHARGESQ(J) = CHARGESQ(J) + 2.0D0 * EXCESS
          PTOTAL(J) = GRAV * RHOX(J) + PZERO
        end do

        ! --- Compute populations ---
        IFEDNS = 0
        call COMPUTE_ONE_POP(0.D0, 1, XNE)
        call COMPUTE_ALL_POPS

      end if  ! IFPRES

      ! --- Radiation energy density (if needed) ---
      if (IFEDNS == 1) call ENERGY_DENSITY

      ! ---------------------------------------------------------------
      ! DOPPLER WIDTHS for all species at all depths
      ! ---------------------------------------------------------------
      !   Doppler width = sqrt(2kT/m + v_turb^2) / c
      !   XNFDOP = (number density / partition function) / (Doppler width * density)
      !          = line opacity coefficient per unit Doppler width

      do J = 1, NRHOX
        do NELION = 1, MION - 1
          DOPPLE(J, NELION) = sqrt(2.0D0 * TK(J) / AMASSISO(1, NELION) / AMU_GRAMS &
                                   + VTURB(J)**2) / CLIGHT_CGS
          DOPPLE8 = DOPPLE(J, NELION)
          XNFDOP(J, NELION) = XNFP(J, NELION) / DOPPLE8 / RHO(J)
        end do
      end do

      ! ---------------------------------------------------------------
      ! OPACITY TABLES (first iteration of first model only)
      ! ---------------------------------------------------------------
      if (ITER == 1 .and. ITEMP == 1) call KAPCONT
      if (ITER == 1 .and. IFREADLINES == 1 .and. ITEMP == 1) call SELECTLINES

      ! --- Compute line opacities ---
      if (IFOP(15) == 1) call LINOP1
      if (IFOP(17) == 1) call XLINOP

      ! ---------------------------------------------------------------
      ! INITIALIZE FREQUENCY INTEGRALS (mode=1: erase accumulators)
      ! ---------------------------------------------------------------
      if (IFCORR == 1) call TCORR(1, 0.0D0)
      call ROSS(1, 0.0D0)
      call RADIAP(1, 0.0D0)
      if (NLTEON == 1) call STATEQ(1, 0.0D0)

      ! ---------------------------------------------------------------
      ! FREQUENCY INTEGRATION — the heart of the calculation
      ! ---------------------------------------------------------------
      !   Loop over all wavelength points. At each frequency:
      !   1. Compute Planck function and stimulated emission correction
      !   2. Assemble total continuous opacity (KAPP)
      !   3. Add line opacity from pre-computed tables
      !   4. Solve the transfer equation (JOSH)
      !   5. Accumulate flux moments for T-correction, Rosseland mean, etc.

      call PUTOUT(1)

      frequency_loop: do NU = NULO, NUHI, NUSTEP

        ! Set current wavelength and frequency
        WAVE   = WAVESET(NU)
        FREQ   = CLIGHT_ANGFREQ / WAVE
        WAVENO = 1.0D7 / WAVE
        RCOWT  = RCOSET(NU)
        FREQLG = log(FREQ)

        ! Compute Planck function B_nu(T) at each depth
        FREQ15 = FREQ / 1.0D15
        do J = 1, NRHOX
          EHVKT(J) = exp(-FREQ * HKT(J))
          STIM(J)  = 1.0D0 - EHVKT(J)
          BNU(J)   = PLANCK_PREFAC * FREQ15**3 * EHVKT(J) / STIM(J)
        end do

        ! Compute continuous opacity at this frequency
        call KAPP

        ! Add pre-computed line opacity (corrected for stimulated emission)
        do J = 1, NRHOX
          stim4(J) = stim(J)
          ALINES(J) = XLINES(J, NU) * stim4(J)
          ALINE(J)  = ALINE(J) + ALINES(J)
        end do

        ! Solve the radiative transfer equation (Feautrier method)
        call JOSH(IFSCAT, IFSURF)

        ! --- Check for unphysical negative values in the solution ---
        found_negative = .false.
        do J = 1, NRHOX
          if (SNU(J) < 0.0D0 .or. JNU(J) < 0.0D0 .or. HNU(J) < 0.0D0) then
            found_negative = .true.
            exit
          end if
        end do

        if (found_negative) then
           ! Floor negative values to small positive number to continue
           do J = 1, NRHOX
              JNU(J) = max(JNU(J), 1.0D-99)
              SNU(J) = max(SNU(J), 1.0D-99)
              HNU(J) = max(HNU(J), 1.0D-99)
           end do
        end if

        ! --- Accumulate frequency integrals (mode=2) ---
        if (IFSURF == 0) then
          if (IFCORR == 1) call TCORR(2, RCOWT)
          call RADIAP(2, RCOWT)
          call ROSS(2, RCOWT)
          if (NLTEON == 1) call STATEQ(2, RCOWT)
        end if

        call PUTOUT(4)

      end do frequency_loop

      ! For surface-only mode, skip the iteration finishing steps
      if (IFSURF <= 0) then

        ! ---------------------------------------------------------------
        ! FINISH ITERATION — Rosseland mean, convection, T-correction
        ! ---------------------------------------------------------------
        call ROSS(3, 0.0D0)
        RX = ROSSTAB(0.D0, 0.D0, 0.D0)
        call RADIAP(3, 0.0D0)
        call COMPUTE_HEIGHT
        if (IFPRES == 1 .and. IFCONV == 1) call CONVEC
        if (IFCORR == 1)  call TCORR(3, 0.0D0)
        if (NLTEON == 1)  call STATEQ(3, 0.0D0)
        if (IFTURB == 1)  call COMPUTE_PTURB
        call PUTOUT(5)

        call system_clock(CLOCK_END)
        write(6, '(A,I3,A,F8.1,A)') ' ITERATION ', ITERAT, ' completed in ', &
             dble(CLOCK_END - CLOCK_START) / dble(CLOCK_RATE), ' seconds'
        write(6,*)
        FLUSH(6)

      end if

    end do iteration_loop

END PROGRAM ATLAS12
