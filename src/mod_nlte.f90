! ==============================================================================
! mod_nlte.f90 — departure coefficients for selected lines in SYNTHE
!
! WHAT THIS IS FOR
! ----------------
! SYNTHE is an LTE synthesizer.  Line opacity is built from Boltzmann/Saha
! populations, and the transfer step (JOSH) is handed ONE line source
! function SLINE per depth per wavelength, which is B_nu.  That is wrong for
! resonance lines of minority species, most sharply for the alkalis: Olander
! et al. (2021, A&A 649, A103) find NLTE corrections of -0.06 to -0.37 dex on
! the K I resonance lines in M dwarfs, driven by photon losses that
! overpopulate the lower level and deepen the line.  Na I D is the same kind
! of transition.  Kurucz anticipated the correction -- gfall reserves
! NBLO/NBUP columns, read_gfall routes any record carrying them into the
! nlte_lines pool, and CODEX whitelists 17 species including Na I -- but the
! preprocessor that filled them (RNLTE) was never part of this port, so the
! hook has always been dead.  This module is the live end of it.
!
! THE PHYSICS, AND WHY ONLY TWO NUMBERS PER LINE PER DEPTH ARE NEEDED
! ------------------------------------------------------------------
! With departure coefficients b_l = n_l/n_l(LTE) and b_u = n_u/n_u(LTE), and
! x = h nu_0 / kT at line centre,
!
!     kappa = kappa_LTE * b_l * [1 - (b_u/b_l) e^-x] / [1 - e^-x]      (1)
!     S_l   = (2h nu^3/c^2) / [ (b_l/b_u) e^x - 1 ]                    (2)
!
! Both stimulated-emission corrections are exact, not Wien-limit.  Writing
! beta = b_u/b_l, (1) and (2) multiply to
!
!     kappa * S_l = b_l beta * kappa_LTE * (2h nu^3/c^2) e^-x
!                 = b_u * [ kappa_LTE (1 - e^-x) ] * B_nu              (3)
!
! -- the line EMISSIVITY is simply b_u times its LTE value.  That identity is
! what makes this cheap to retrofit.  SYNTHE accumulates a single summed line
! opacity per (lambda, depth); the opacity-weighted mean source function that
! JOSH wants is
!
!     SLINE = sum_i kappa_i S_i / sum_i kappa_i
!           = B_nu * [ sum_i kappa_i + sum_i kappa_i (r_i - 1) ] / sum_i kappa_i
!
! with r_i = S_i/B_nu.  The second sum runs only over lines that are NOT in
! LTE, because r_i - 1 vanishes for every other line in the window.  So the
! whole retrofit is one extra accumulator carrying the DEVIATION of the
! emissivity from LTE, and the tens of millions of LTE lines never touch it.
!
! Concretely, for a tagged line this module returns two scalars per depth:
!
!     fkappa = b_l [1 - beta e^-x] / [1 - e^-x]     opacity scale, eq. (1)
!     fdev   = (b_u - fkappa) / fkappa              deviation per unit of the
!                                                   ALREADY-SCALED opacity
!
! so that SYNTHE multiplies the line's kappa0 by fkappa as usual and adds
! profile*fdev into the deviation accumulator using the same profile array
! and the same index arithmetic.  Note nu_0 is used across the whole profile;
! over a line width e^-x is constant to ~1e-5.
!
! MODES
! -----
! Set NLTE_MODE in mod_parameters (developer flag, not a CLI option).
!
!   0  pure LTE.  nlte_on stays .FALSE., nothing here allocates or runs.
!   1  null/identity test: b_l = NLTE_TEST_BLO, b_u = NLTE_TEST_BUP at every
!      depth, applied to every tagged component.  With BLO = BUP = 1 the
!      arithmetic below collapses to fkappa = 1.0 and fdev = 0.0 EXACTLY --
!      not to within rounding -- because 1 - beta e^-x and 1 - e^-x are then
!      the same expression and y/y is exactly 1 in IEEE arithmetic, so mode 1
!      reproduces mode 0 bit for bit while still exercising every line of the
!      new path.  With BLO = BUP = c it scales the tagged lines' opacity by
!      exactly c, S_l still = B_nu, i.e. identical to shifting their log gf by
!      log10(c) -- a check of the kappa path against an answer known in
!      advance.  With BUP < BLO the source function drops below B_nu and the
!      deviation accumulator carries a nonzero number for the first time.
!   2  departure coefficients from a <model>.dep sidecar, produced from a
!      published grid by tools/nlte_make_dep.py and interpolated here onto
!      the model's own layers.  See read_dep_file for the format, the
!      interpolation, and why the grid reader lives in Python rather than
!      here.  Data source: the Amarsi et al. (2020) MARCS grids (3000-8000 K,
!      Na included), as repackaged by Gerber et al. (2023) for Turbospectrum.
!
! Which transitions are eligible is declared in mod_mklinelist
! (NLTE_TR_CODE/ELO/EUP); read_gfall tags the matching lte_lines entries
! unconditionally, exactly as it records ICA4227.
! ==============================================================================

MODULE mod_nlte

  USE mod_parameters, ONLY: kw, NLTE_MODE, NLTE_TEST_SHAPE, &
                            NLTE_TEST_BLO, NLTE_TEST_BUP
  USE mod_atlas_data, ONLY: DATADIR
  USE mod_mklinelist, ONLY: NLTE_NTRANS, N_NLTE_TAGGED, &
                            NLTE_TAG_IDX, NLTE_TAG_TRANS, &
                            NLTE_TR_LEVLO, NLTE_TR_LEVUP

  IMPLICIT NONE
  PUBLIC

  ! .TRUE. when NLTE_MODE /= 0 AND at least one tagged component is in the
  ! synthesis window.  Everything on the hot path tests this first, so a
  ! production run pays one logical compare per line.
  LOGICAL, SAVE :: nlte_on = .FALSE.

  ! Departure coefficients, per eligible transition per depth.  Initialised
  ! to 1 so that any depth or transition left unfilled is plain LTE.
  REAL(8), SAVE :: blo(NLTE_NTRANS, kw) = 1.0D0
  REAL(8), SAVE :: bup(NLTE_NTRANS, kw) = 1.0D0

  ! --- STRUCTURED TEST PROFILE (NLTE_MODE = 1, NLTE_TEST_SHAPE = 1) --------
  !
  ! Per transition, b ramps linearly in log10(column mass) between a DEEP
  ! value at the bottom of the atmosphere and a TOP value at the surface:
  !
  !     b(j) = b_deep + (b_top - b_deep) * f(j)
  !     f(j) = [log10(rhox_bottom) - log10(rhox_j)]
  !            / [log10(rhox_bottom) - log10(rhox_top)]      f: 1 top, 0 bottom
  !
  ! The numbers are invented, but the SHAPE is deliberate on three counts.
  ! (i) Both b_l and b_u sit at 1 deep down, where collisions enforce LTE --
  ! so the wings, which form deep, stay near their LTE values and the change
  ! is confined to the core.  A run that moves the far wings as much as the
  ! core has its depth axis reversed.  (ii) b_l /= b_u AND b_l /= 1, which
  ! the constant shape never arranges: with b_l = b_u the two are symmetric
  ! and a swap of their roles is invisible, and with b_l = 1 a b_l that was
  ! silently ignored would look correct.  (iii) The two transitions get
  ! DIFFERENT values, so coefficients attributed to the wrong line show up
  ! as D1 and D2 responding identically.
  !
  ! Ordering matches the transition table in mod_mklinelist: 1 = Na I D2,
  ! 2 = Na I D1.  S_l/B_nu at the surface is b_u/b_l = 0.50 for D2 and 0.73
  ! for D1 -- distinct, and both in the photon-loss direction.
  REAL(8), PARAMETER :: BLO_DEEP(NLTE_NTRANS) = [1.00D0, 1.00D0]
  REAL(8), PARAMETER :: BLO_TOP (NLTE_NTRANS) = [1.20D0, 1.10D0]
  REAL(8), PARAMETER :: BUP_DEEP(NLTE_NTRANS) = [1.00D0, 1.00D0]
  REAL(8), PARAMETER :: BUP_TOP (NLTE_NTRANS) = [0.60D0, 0.80D0]

  ! --- NLTE_MODE = 3: the single runtime grid file -------------------------
  ! Built by tools/nlte_build_runtime.py from the master store.  ONE file:
  ! axes, then the parameter -> record index, then the records themselves, so
  ! an index can never be paired with the wrong data, and there is one thing
  ! to install.  Layout (little-endian):
  !     int32   nT,nG,nF,nV,nD, nlev, ndep, nrec
  !     float32 the five axis value lists
  !     int32   level_id(nlev)     original atom level numbers
  !     int32   rec(...)           1-based record, 0 = absent, Fortran order
  !     int8    geom(...)          0 absent, 1 plane-parallel, 2 spherical
  !     float32 data(nrec, ndep + nlev*ndep)   tau then b(level,depth)
  ! Records are read by stream access at an explicit byte POS, so the
  ! arithmetic is visible and does not depend on what a compiler calls a
  ! "record length".  Dimensions come from the header rather than being
  ! compiled in, so a regenerated file with different levels or depths is
  ! picked up without touching this source.
  !
  ! Lives flat in $ATLAS12/data/nlte/ under a name carrying the element and
  ! the grid's release date, so several elements coexist without
  ! subdirectories.  $NLTE_GRID overrides with an absolute path.
  CHARACTER(LEN=*), PARAMETER :: NLTE_GRID_NAME = &
      'nlte/Na_MARCS_Jul-14-2023.nlte'

  INTEGER, SAVE :: gnT=0, gnG=0, gnF=0, gnV=0, gnD=0
  INTEGER, SAVE :: gnlev=0, gndep=0, gnrec=0
  INTEGER, SAVE :: grecw=0                  ! floats per record
  INTEGER(8), SAVE :: gdata0=0              ! byte offset of the data block
  INTEGER, ALLOCATABLE, SAVE :: glevid(:)   ! atom level number per slot
  ! Slot within a record for each transition's lower/upper level, resolved
  ! once from glevid: the transition table names ATOM levels, while which of
  ! them a given runtime file carries is a property of that file.
  INTEGER, SAVE :: gslot_lo(NLTE_NTRANS)=0, gslot_up(NLTE_NTRANS)=0
  REAL(4),    ALLOCATABLE, SAVE :: gT(:), gG(:), gF(:), gV(:), gDn(:)
  INTEGER(4), ALLOCATABLE, SAVE :: grec(:,:,:,:,:)
  INTEGER(1), ALLOCATABLE, SAVE :: ggeo(:,:,:,:,:)
  CHARACTER(LEN=512), SAVE :: gprefix = ''
  ! Set from the run by nlte_set_params, used to pick the cell.
  REAL(8), SAVE :: m_teff=0, m_logg=0, m_feh=0, m_vturb=0, m_ana=0

  ! Solar reference for turning the model's absolute abundances into [Fe/H].
  ! C3K models fold metallicity into the ABUNDANCE CHANGE values and leave
  ! ABUNDANCE SCALE at 1, so [Fe/H] cannot be read off the scale factor; it
  ! has to be A(Fe)_model - A(Fe)_sun.  The C3K grids are built on Grevesse &
  ! Sauval (1998) (see scripts/atlas12.sh), and the MARCS grid's own [Fe/H]
  ! axis is on a comparable scale, so GS98 is used here.  Different solar
  ! scales differ by a few 0.01 dex, which shifts the metallicity lookup
  ! slightly; A(Na) itself is absolute and unaffected.
  REAL(8), PARAMETER :: A_FE_SUN = 7.50D0     ! GS98 log eps(Fe)
  ! Reported after interpolation: how much of the 32-corner weight was
  ! actually available, and whether any corner came from a spherical model.
  REAL(8), SAVE :: g_wfrac = 0.0D0
  INTEGER, SAVE :: g_nsph = 0, g_ncorn = 0

  ! Cursor into NLTE_TAG_IDX for nlte_tag_at(); see the contract there.
  INTEGER, SAVE :: scan_ptr = 1

  ! --- Diagnostics, dumped to <model>.nlte when nlte_on --------------------
  ! Depth structure kept for the dump, and what nlte_factors ACTUALLY
  ! computed for each (transition, depth).  Recording the factors rather
  ! than only blo/bup is the point: it makes the dump able to catch an
  ! indexing error inside nlte_factors, not just one in fill_departures.
  INTEGER, SAVE :: nrhox_sv = 0
  REAL(8), SAVE :: rhox_sv(kw) = 0.0D0
  REAL(8), SAVE :: temp_sv(kw) = 0.0D0
  ! Continuum tau_5000 on this model's layers.  The published grids are
  ! tabulated against tau_5000, so this is what lets their b be placed on
  ! our layers without going through the MARCS model they were solved on.
  REAL(8), SAVE :: tau5_sv(kw) = 0.0D0
  REAL(8), SAVE :: fkap_sv(NLTE_NTRANS, kw) = 0.0D0
  REAL(8), SAVE :: fdev_sv(NLTE_NTRANS, kw) = 0.0D0
  REAL(8), SAVE :: xhnu_sv(NLTE_NTRANS, kw) = 0.0D0
  LOGICAL, SAVE :: seen_sv(NLTE_NTRANS, kw) = .FALSE.

  ! (transition, depth) pairs where b_u/b_l implied a population inversion
  ! and the line was forced back to LTE.  Flagged rather than counted so a
  ! transition evaluated once per hyperfine component is not counted ten
  ! times over; reported and marked in the dump.
  LOGICAL, SAVE :: inverted_sv(NLTE_NTRANS, kw) = .FALSE.

  ! Depths where mode 2 had to hold an endpoint of the .dep table because
  ! the grid did not reach that far.  Flagged in the dump so that a grid
  ! that simply does not cover the model cannot be mistaken for a broken
  ! interpolation.
  LOGICAL, SAVE :: clamped_sv(kw) = .FALSE.

CONTAINS

  ! --------------------------------------------------------------------------
  !  nlte_init(nrhox) -- decide whether NLTE is active and fill blo/bup.
  !
  !  Call once, after run_mklinelist (which does the tagging) and after the
  !  depth structure is known.  A mode that is switched on but finds no
  !  tagged line in the window is reported and then treated as LTE -- that
  !  is a legitimate outcome (synthesise 4000-4500 A and Na D is simply not
  !  there), unlike CA4227_MODE, which names one line and so hard-fails.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_init(nrhox, rhox, temp, tau5, depfile)

    INTEGER,          INTENT(IN) :: nrhox
    REAL(8),          INTENT(IN) :: rhox(:)   ! column mass [g/cm^2], 1=surface
    REAL(8),          INTENT(IN) :: temp(:)   ! temperature [K]
    REAL(8),          INTENT(IN) :: tau5(:)   ! continuum tau_5000, same layers
    CHARACTER(LEN=*), INTENT(IN) :: depfile   ! <model>.dep, used by mode 2
    INTEGER :: i

    nlte_on     = .FALSE.
    inverted_sv = .FALSE.
    seen_sv     = .FALSE.
    clamped_sv  = .FALSE.
    blo         = 1.0D0
    bup         = 1.0D0

    nrhox_sv              = nrhox
    rhox_sv(1:nrhox)      = rhox(1:nrhox)
    temp_sv(1:nrhox)      = temp(1:nrhox)
    tau5_sv(1:nrhox)      = tau5(1:nrhox)

    IF (NLTE_MODE .EQ. 0) RETURN

    IF (N_NLTE_TAGGED .EQ. 0) THEN
      WRITE(6,'(A,I0,A)') ' NLTE: mode ', NLTE_MODE, &
        ' is set but no eligible transition lies in the synthesis' // &
        ' window; continuing in LTE.'
      RETURN
    END IF

    CALL fill_departures(nrhox, depfile)

    nlte_on = .TRUE.

    WRITE(6,'(A,I0,A,I0,A,I0,A)') ' NLTE: mode ', NLTE_MODE, ', ', &
      N_NLTE_TAGGED, ' tagged line components over ', NLTE_NTRANS, &
      ' transitions'
    IF (NLTE_MODE .EQ. 1 .AND. NLTE_TEST_SHAPE .EQ. 0) &
      WRITE(6,'(A,F10.6,A,F10.6,A)') ' NLTE: constant test, b_lo = ', &
        NLTE_TEST_BLO, ', b_up = ', NLTE_TEST_BUP, ' at all depths'
    IF (NLTE_MODE .EQ. 1 .AND. NLTE_TEST_SHAPE .EQ. 1) THEN
      WRITE(6,'(A)') ' NLTE: structured test, b ramped in log10(column mass)'
      DO i = 1, NLTE_NTRANS
        WRITE(6,'(A,I0,A,F6.3,A,F6.3,A,F6.3,A,F6.3)') &
          ' NLTE:   transition ', i, &
          ':  b_lo ', BLO_DEEP(i), ' ->', BLO_TOP(i), &
          ',  b_up ', BUP_DEEP(i), ' ->', BUP_TOP(i)
      END DO
    END IF
    DO i = 1, N_NLTE_TAGGED
      IF (i .LE. 12) &
        WRITE(6,'(A,I9,A,I0)') ' NLTE:   lte_lines index', NLTE_TAG_IDX(i), &
          '  transition ', NLTE_TAG_TRANS(i)
    END DO
    IF (N_NLTE_TAGGED .GT. 12) &
      WRITE(6,'(A,I0,A)') ' NLTE:   (', N_NLTE_TAGGED - 12, ' more)'

  END SUBROUTINE nlte_init


  ! --------------------------------------------------------------------------
  !  fill_departures(nrhox) -- populate blo/bup.
  !
  !  This is the ONLY routine that needs replacing to move from the null
  !  test to real departure coefficients: everything downstream consumes
  !  blo/bup and nothing else.
  ! --------------------------------------------------------------------------
  SUBROUTINE fill_departures(nrhox, depfile)

    INTEGER,          INTENT(IN) :: nrhox
    CHARACTER(LEN=*), INTENT(IN) :: depfile
    INTEGER :: k, j
    REAL(8) :: lr_top, lr_bot, span, f

    SELECT CASE (NLTE_MODE)

    CASE (1)

      IF (NLTE_TEST_SHAPE .EQ. 0) THEN
        ! Constant b at every depth and every transition.
        DO j = 1, nrhox
          DO k = 1, NLTE_NTRANS
            blo(k,j) = NLTE_TEST_BLO
            bup(k,j) = NLTE_TEST_BUP
          END DO
        END DO

      ELSE
        ! Structured: linear ramp in log10(column mass), per transition.
        ! Column mass is the same variable mode 2 must interpolate its
        ! grid on, so exercising it here is not incidental.
        IF (nrhox .LT. 2 .OR. rhox_sv(1) .LE. 0.0D0 .OR. &
            rhox_sv(nrhox) .LE. 0.0D0) THEN
          WRITE(6,'(A)') ' ERROR: NLTE structured test needs a positive,' // &
                         ' multi-layer column-mass scale.'
          CALL EXIT(1)
        END IF
        lr_top = LOG10(rhox_sv(1))
        lr_bot = LOG10(rhox_sv(nrhox))
        span   = lr_bot - lr_top
        IF (ABS(span) .LT. 1.0D-8) THEN
          WRITE(6,'(A)') ' ERROR: NLTE structured test: column mass does' // &
                         ' not span a range.'
          CALL EXIT(1)
        END IF
        DO j = 1, nrhox
          f = (lr_bot - LOG10(rhox_sv(j))) / span    ! 1 at surface, 0 at base
          DO k = 1, NLTE_NTRANS
            blo(k,j) = BLO_DEEP(k) + (BLO_TOP(k) - BLO_DEEP(k)) * f
            bup(k,j) = BUP_DEEP(k) + (BUP_TOP(k) - BUP_DEEP(k)) * f
          END DO
        END DO
      END IF

    CASE (2)
      ! Departure coefficients from an external grid, pre-extracted for
      ! this model into a <model>.dep sidecar and interpolated here onto
      ! the model's own depth scale.
      CALL read_dep_file(depfile, nrhox)

    CASE (3)
      ! Straight from the local extracted grid, interpolated in stellar
      ! parameters here.  No sidecar, no per-model preparation step.
      CALL interp_grid(nrhox)

    CASE DEFAULT
      WRITE(6,'(A,I0)') ' ERROR: unknown NLTE_MODE = ', NLTE_MODE
      CALL EXIT(1)

    END SELECT

  END SUBROUTINE fill_departures



  ! --------------------------------------------------------------------------
  !  read_grid_index() -- load <prefix>.idx once.
  !
  !  The index is cell -> record number over the five real axes
  !  (Teff, log g, [Fe/H], v_turb, dNa), written by tools/nlte_build_index.py.
  !  [alpha/Fe] is deliberately absent: the MARCS '_st_' models tie it to
  !  [Fe/H] by the standard relation, so it is a label rather than a
  !  dimension, and an alpha-enhanced model necessarily gets b computed at
  !  the standard alpha for its metallicity.
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  !  nlte_set_params -- tell mode 3 which model this is.
  !
  !  Abundances come in as the code stores them: ABUND(IZ) is log10(N_Z/N_tot)
  !  for Z >= 3 with XRELATIVE(IZ) an additive log offset, while ABUND(1) is
  !  the hydrogen number FRACTION.  So
  !        A(X) = ABUND(X) + XRELATIVE(X) - log10(ABUND(H)) + 12
  !  and [Fe/H] follows by subtracting the solar reference above.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_set_params(teff, logg, vturb, abund_h, abund_na, abund_fe)
    REAL(8), INTENT(IN) :: teff, logg, vturb, abund_h, abund_na, abund_fe
    m_teff  = teff
    m_logg  = logg
    m_vturb = vturb
    m_ana   = abund_na - LOG10(abund_h) + 12.0D0
    m_feh   = (abund_fe - LOG10(abund_h) + 12.0D0) - A_FE_SUN
  END SUBROUTINE nlte_set_params


  SUBROUTINE read_grid_index()

    INTEGER :: u, hdr(8), n, k, i
    LOGICAL :: ex
    CHARACTER(LEN=512) :: fn
    INTEGER :: envlen, envstat
    CHARACTER(LEN=512) :: envval

    IF (gnT .GT. 0) RETURN                      ! already loaded

    CALL GET_ENVIRONMENT_VARIABLE('NLTE_GRID', envval, envlen, envstat)
    IF (envstat .EQ. 0 .AND. envlen .GT. 0) THEN
      fn = TRIM(envval)
    ELSE
      fn = TRIM(DATADIR) // NLTE_GRID_NAME
    END IF
    gprefix = fn

    INQUIRE(FILE=fn, EXIST=ex)
    IF (.NOT. ex) THEN
      WRITE(6,'(A,A)') ' ERROR: NLTE_MODE = 3 needs the runtime grid file;' &
                       // ' not found: ', TRIM(fn)
      WRITE(6,'(A)')   '        Build it with tools/nlte_extract_grid.py,' &
                       // ' nlte_build_index.py and nlte_build_runtime.py'
      WRITE(6,'(A)')   '        into $ATLAS12/data/nlte/, or point' &
                       // ' $NLTE_GRID at it.'
      CALL EXIT(1)
    END IF

    OPEN(NEWUNIT=u, FILE=fn, ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ')
    READ(u) hdr
    gnT = hdr(1);  gnG = hdr(2);  gnF = hdr(3);  gnV = hdr(4);  gnD = hdr(5)
    gnlev = hdr(6);  gndep = hdr(7);  gnrec = hdr(8)
    grecw = gndep + gnlev*gndep

    ALLOCATE(gT(gnT), gG(gnG), gF(gnF), gV(gnV), gDn(gnD), glevid(gnlev))
    READ(u) gT
    READ(u) gG
    READ(u) gF
    READ(u) gV
    READ(u) gDn
    READ(u) glevid
    ALLOCATE(grec(gnT,gnG,gnF,gnV,gnD), ggeo(gnT,gnG,gnF,gnV,gnD))
    READ(u) grec
    READ(u) ggeo
    INQUIRE(UNIT=u, POS=gdata0)        ! first byte of the data block
    CLOSE(u)

    ! Resolve each transition's atom levels to slots in the record.  The
    ! transition table names ATOM levels; which of them this file happens to
    ! carry is a property of the file, so a runtime file built with a
    ! narrower level set must fail loudly rather than read the wrong slot.
    DO k = 1, NLTE_NTRANS
      gslot_lo(k) = 0;  gslot_up(k) = 0
      DO i = 1, gnlev
        IF (glevid(i) .EQ. NLTE_TR_LEVLO(k)) gslot_lo(k) = i
        IF (glevid(i) .EQ. NLTE_TR_LEVUP(k)) gslot_up(k) = i
      END DO
      IF (gslot_lo(k) .EQ. 0 .OR. gslot_up(k) .EQ. 0) THEN
        WRITE(6,'(A,I0,A,I0,A,I0,A)') ' ERROR: transition ', k, ' needs atom' &
          // ' levels ', NLTE_TR_LEVLO(k), ' and ', NLTE_TR_LEVUP(k), &
          ' but the runtime grid file does not carry both.'
        WRITE(6,'(A)') '        Rebuild it with a wider --levels set.'
        CALL EXIT(1)
      END IF
    END DO

    n = COUNT(grec .GT. 0)
    WRITE(6,'(A,A)') ' NLTE: grid  ', TRIM(fn)
    WRITE(6,'(A,I0,4(A,I0),3(A,I0),A)') &
      ' NLTE:   axes ', gnT,'x',gnG,'x',gnF,'x',gnV,'x',gnD, &
      ', ', gnlev, ' levels, ', gnrec, ' records, ', &
      NINT(100.0*REAL(n)/REAL(SIZE(grec))), ' % of cells filled'

  END SUBROUTINE read_grid_index


  ! --------------------------------------------------------------------------
  !  bracket(ax, n, x, i, w) -- locate x on a sorted axis.
  !
  !  Returns the lower node i and the weight w of node i+1, both CLAMPED to
  !  the axis ends: outside the grid the edge value is held, never
  !  extrapolated, exactly as read_dep_file does in depth.
  ! --------------------------------------------------------------------------
  SUBROUTINE bracket(ax, n, x, i, w)
    REAL(4), INTENT(IN)  :: ax(:)
    INTEGER, INTENT(IN)  :: n
    REAL(8), INTENT(IN)  :: x
    INTEGER, INTENT(OUT) :: i
    REAL(8), INTENT(OUT) :: w
    INTEGER :: lo, hi, mid
    IF (n .EQ. 1) THEN
      i = 1;  w = 0.0D0;  RETURN
    END IF
    IF (x .LE. DBLE(ax(1)))  THEN
      i = 1;    w = 0.0D0;  RETURN
    END IF
    IF (x .GE. DBLE(ax(n))) THEN
      i = n-1;  w = 1.0D0;  RETURN
    END IF
    lo = 1;  hi = n
    DO WHILE (hi - lo .GT. 1)
      mid = (lo + hi)/2
      IF (DBLE(ax(mid)) .GT. x) THEN
        hi = mid
      ELSE
        lo = mid
      END IF
    END DO
    i = lo
    w = (x - DBLE(ax(lo))) / (DBLE(ax(lo+1)) - DBLE(ax(lo)))
  END SUBROUTINE bracket


  ! --------------------------------------------------------------------------
  !  interp_grid(nrhox) -- NLTE_MODE = 3.
  !
  !  Multilinear interpolation of b over the five axes, evaluated on this
  !  model's own layers.
  !
  !  TWO THINGS MAKE THIS MORE THAN A WEIGHTED SUM.
  !
  !  Every corner record carries its OWN tau_5000 grid -- they are different
  !  MARCS atmospheres -- so b cannot be combined index by index.  Each corner
  !  is first put onto THIS model's tau_5000 (computed by run_xnfpelsyn), and
  !  only then combined.  Interpolation is in log b against log tau, with the
  !  corner's endpoint value held outside its range for the same reason as in
  !  depth: b turns over sharply at the top and a linear continuation there is
  !  confident nonsense.
  !
  !  And corners can be MISSING: the grid is 41 per cent filled, because it is
  !  HR-diagram shaped.  Weight from absent corners is dropped and the rest
  !  renormalised, which is a real approximation near the edges of the
  !  populated region -- so the surviving weight fraction is recorded and
  !  reported rather than silently absorbed.
  ! --------------------------------------------------------------------------
  SUBROUTINE interp_grid(nrhox)

    INTEGER, INTENT(IN) :: nrhox

    INTEGER :: iT,iG,iF,iV,iD, cT,cG,cF,cV,cD, k, j, u, r
    REAL(8) :: wT,wG,wF,wV,wD, wgt, wsum
    REAL(4), ALLOCATABLE :: buf(:)
    REAL(8), ALLOCATABLE :: ltau(:), lb(:)
    REAL(8) :: lt, t
    REAL(8) :: acc_lo(NLTE_NTRANS,kw), acc_up(NLTE_NTRANS,kw)
    INTEGER :: lo, hi, mid, lev
    LOGICAL :: ex

    CALL read_grid_index()

    CALL bracket(gT,  gnT, m_teff,          iT, wT)
    CALL bracket(gG,  gnG, m_logg,          iG, wG)
    CALL bracket(gF,  gnF, m_feh,           iF, wF)
    CALL bracket(gV,  gnV, m_vturb,         iV, wV)
    CALL bracket(gDn, gnD, m_ana - m_feh,   iD, wD)

    acc_lo = 0.0D0;  acc_up = 0.0D0
    wsum   = 0.0D0
    g_nsph = 0;  g_ncorn = 0

    ALLOCATE(buf(grecw), ltau(gndep), lb(gndep))
    OPEN(NEWUNIT=u, FILE=TRIM(gprefix), ACCESS='STREAM', &
         FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')

    DO cT = 0, MIN(1, gnT-1)
     DO cG = 0, MIN(1, gnG-1)
      DO cF = 0, MIN(1, gnF-1)
       DO cV = 0, MIN(1, gnV-1)
        DO cD = 0, MIN(1, gnD-1)
          wgt = (MERGE(wT, 1.0D0-wT, cT.EQ.1)) * (MERGE(wG, 1.0D0-wG, cG.EQ.1)) &
              * (MERGE(wF, 1.0D0-wF, cF.EQ.1)) * (MERGE(wV, 1.0D0-wV, cV.EQ.1)) &
              * (MERGE(wD, 1.0D0-wD, cD.EQ.1))
          IF (wgt .LE. 0.0D0) CYCLE
          r = grec(iT+cT, iG+cG, iF+cF, iV+cV, iD+cD)
          IF (r .LE. 0) CYCLE                    ! absent corner
          g_ncorn = g_ncorn + 1
          IF (ggeo(iT+cT, iG+cG, iF+cF, iV+cV, iD+cD) .EQ. 2) g_nsph = g_nsph + 1

          READ(u, POS=gdata0 + INT(r-1, 8)*INT(grecw, 8)*4_8) buf
          DO j = 1, gndep
            ltau(j) = LOG10(MAX(DBLE(buf(j)), 1.0D-99))
          END DO

          DO k = 1, NLTE_NTRANS
            DO lev = 1, 2
              IF (lev .EQ. 1) THEN
                DO j = 1, gndep
                  lb(j) = LOG10(MAX(DBLE(buf(gndep + (gslot_lo(k)-1)*gndep + j)), 1.0D-99))
                END DO
              ELSE
                DO j = 1, gndep
                  lb(j) = LOG10(MAX(DBLE(buf(gndep + (gslot_up(k)-1)*gndep + j)), 1.0D-99))
                END DO
              END IF
              DO j = 1, nrhox
                lt = LOG10(MAX(tau5_sv(j), 1.0D-99))
                ! Mark layers outside the corner's own tau range.  Mode 2
                ! reports this and it is how the grid's truncation became
                ! visible at all; mode 3 must not lose it.
                IF (lt .LT. ltau(1) .OR. lt .GT. ltau(gndep)) clamped_sv(j) = .TRUE.
                IF (lt .LE. ltau(1)) THEN
                  t = lb(1)
                ELSE IF (lt .GE. ltau(gndep)) THEN
                  t = lb(gndep)
                ELSE
                  lo = 1;  hi = gndep
                  DO WHILE (hi - lo .GT. 1)
                    mid = (lo + hi)/2
                    IF (ltau(mid) .GT. lt) THEN
                      hi = mid
                    ELSE
                      lo = mid
                    END IF
                  END DO
                  t = lb(lo) + (lt - ltau(lo)) / (ltau(lo+1) - ltau(lo)) &
                             * (lb(lo+1) - lb(lo))
                END IF
                IF (lev .EQ. 1) THEN
                  acc_lo(k,j) = acc_lo(k,j) + wgt*t
                ELSE
                  acc_up(k,j) = acc_up(k,j) + wgt*t
                END IF
              END DO
            END DO
          END DO
          wsum = wsum + wgt
        END DO
       END DO
      END DO
     END DO
    END DO
    CLOSE(u)
    DEALLOCATE(buf, ltau, lb)

    IF (wsum .LE. 0.0D0) THEN
      WRITE(6,'(A)') ' ERROR: NLTE_MODE = 3 found no populated grid cell'  &
                     // ' near this model.'
      WRITE(6,'(A,F8.1,A,F6.2,A,F6.2,A,F5.1,A,F6.2)') &
        '        Teff=', m_teff, ' logg=', m_logg, ' [Fe/H]=', m_feh, &
        ' vturb=', m_vturb, ' A(Na)=', m_ana
      CALL EXIT(1)
    END IF

    g_wfrac = wsum
    DO k = 1, NLTE_NTRANS
      DO j = 1, nrhox
        blo(k,j) = 10.0D0**(acc_lo(k,j)/wsum)
        bup(k,j) = 10.0D0**(acc_up(k,j)/wsum)
      END DO
    END DO

    WRITE(6,'(A,F7.1,A,F5.2,A,F6.2,A,F4.1,A,F5.2)') &
      ' NLTE: model  Teff=', m_teff, '  logg=', m_logg, '  [Fe/H]=', m_feh, &
      '  vturb=', m_vturb, '  A(Na)=', m_ana
    WRITE(6,'(A,I3,A,F6.3,A,I0,A)') ' NLTE:   ', g_ncorn, &
      ' corners used, weight fraction ', g_wfrac, ', ', g_nsph, ' spherical'
    IF (COUNT(clamped_sv(1:nrhox)) .GT. 0) &
      WRITE(6,'(A,I0,A,I0,A)') ' NLTE: WARNING -- grid does not cover the' // &
        ' model in depth: ', COUNT(clamped_sv(1:nrhox)), ' of ', nrhox, &
        ' layers hold an endpoint b (not extrapolated).'
    IF (g_wfrac .LT. 0.999D0) &
      WRITE(6,'(A,F6.3,A)') ' NLTE: WARNING -- only ', g_wfrac, &
        ' of the interpolation weight was populated; the rest was dropped' &
        // ' and the remainder renormalised.'
    IF (g_nsph .GT. 0 .AND. g_nsph .LT. g_ncorn) &
      WRITE(6,'(A)') ' NLTE: WARNING -- corners mix spherical and' &
        // ' plane-parallel MARCS geometry.'

  END SUBROUTINE interp_grid


  ! --------------------------------------------------------------------------
  !  read_dep_file(depfile, nrhox) -- NLTE_MODE = 2.
  !
  !  Reads <model>.dep and interpolates it onto this model's depth scale.
  !
  !  WHY A SIDECAR RATHER THAN THE GRID ITSELF
  !  -----------------------------------------
  !  The published departure-coefficient grids (Amarsi et al. 2020, as
  !  repackaged by Gerber et al. 2023) are multi-gigabyte binaries indexed
  !  by (Teff, log g, [Fe/H], [alpha/Fe]) over thousands of MARCS models,
  !  with their own auxiliary index files and their own level ordering per
  !  model atom.  Putting a reader for that inside SYNTHE would bury the
  !  physics under file format, and would have to be rewritten whenever the
  !  grid is reissued.  Instead tools/nlte_make_dep.py does the grid I/O,
  !  the stellar-parameter interpolation, and the model-atom level -> our
  !  transition mapping, and emits one small text file per model.  What is
  !  left here is the part that genuinely belongs in the synthesiser:
  !  putting b onto this atmosphere's layers.
  !
  !  FILE FORMAT (free-form, '#' comments and blank lines ignored)
  !  ------------------------------------------------------------
  !    NTRANS  <m>          must equal NLTE_NTRANS; transition order must
  !                         match NLTE_TR_* in mod_mklinelist
  !    NDEPTH  <n>          number of table rows, n >= 2
  !    DEPTHVAR LOGRHOX     optional, but if present must say LOGRHOX.  The
  !                         published grids are tabulated against tau_5000,
  !                         which SYNTHE does not carry; the converter maps
  !                         them onto column mass using the MARCS model the
  !                         grid was computed on.  The keyword exists so that
  !                         a tau-based table cannot be fed in silently and
  !                         read as if its first column were log column mass.
  !    then n rows of  log10(rhox[g/cm^2])  b_lo(1) b_up(1) ... b_lo(m) b_up(m)
  !    rows ordered by INCREASING log10(rhox), i.e. surface first
  !
  !  INTERPOLATION
  !  -------------
  !  Linear in log10(b) against log10(column mass).  Column mass, not tau,
  !  because the grids are computed on MARCS structures and an equal-tau
  !  mapping between codes whose continuum opacities differ misplaces the
  !  layers -- the same trap that once faked a 440 K offset in a cross-code
  !  temperature comparison here.  log b rather than b because b is a ratio
  !  that can run over decades, and because interpolating it linearly can
  !  cross zero on a steep gradient while log b cannot.  All b in the file
  !  must therefore be strictly positive; a non-positive value is an error,
  !  not something to quietly floor.
  !
  !  OUTSIDE THE FILE'S RANGE the endpoint value is HELD, never
  !  extrapolated: departure coefficients turn over sharply at the top of
  !  the atmosphere and a linear continuation there produces confident
  !  nonsense.  Clamped layers are counted, reported, and flagged in the
  !  .nlte dump, because "my grid did not cover the model" and "my
  !  interpolation is wrong" are the two failure modes that otherwise look
  !  identical.
  ! --------------------------------------------------------------------------
  SUBROUTINE read_dep_file(depfile, nrhox)

    CHARACTER(LEN=*), INTENT(IN) :: depfile
    INTEGER,          INTENT(IN) :: nrhox

    INTEGER, PARAMETER :: MAXDEP = 500     ! rows allowed in the sidecar

    REAL(8) :: lrx(MAXDEP)
    REAL(8) :: blo_f(NLTE_NTRANS, MAXDEP), bup_f(NLTE_NTRANS, MAXDEP)
    REAL(8) :: row(1 + 2*NLTE_NTRANS)
    REAL(8) :: lr, t, lb, rtest
    INTEGER :: u, ios, ntr, ndep, n, k, j, i, lo, hi, mid
    INTEGER :: nclamp_top, nclamp_bot
    CHARACTER(LEN=512) :: line
    CHARACTER(LEN=32)  :: tag, dvar
    LOGICAL :: exists

    INQUIRE(FILE=depfile, EXIST=exists)
    IF (.NOT. exists) THEN
      WRITE(6,'(A,A)') ' ERROR: NLTE_MODE = 2 needs a departure-coefficient' // &
                       ' sidecar; not found: ', TRIM(depfile)
      WRITE(6,'(A)')   '        Build one with tools/nlte_make_dep.py.'
      CALL EXIT(1)
    END IF

    OPEN(NEWUNIT=u, FILE=depfile, STATUS='OLD', ACTION='READ')

    ntr  = -1
    ndep = -1
    n    =  0

    DO
      READ(u,'(A)',IOSTAT=ios) line
      IF (ios .NE. 0) EXIT
      line = ADJUSTL(line)
      IF (LEN_TRIM(line) .EQ. 0) CYCLE
      IF (line(1:1) .EQ. '#')    CYCLE

      ! Keyword or table row?  Decided by whether the first token parses as
      ! a number, not by position: an earlier version kept a "still in the
      ! header" flag that cleared as soon as NTRANS and NDEPTH had both been
      ! seen, so a DEPTHVAR line written after them was read as data and
      ! reported as a malformed row.
      READ(line,*,IOSTAT=ios) rtest
      IF (ios .NE. 0) THEN
        ! Read the keyword ALONE first, then re-read the line with the type
        ! that keyword's argument actually has.  Demanding an integer up
        ! front misreports "DEPTHVAR LOGTAU5000" as a malformed NTRANS.
        READ(line,*,IOSTAT=ios) tag
        IF (ios .NE. 0) THEN
          WRITE(6,'(A,A)') ' ERROR: .dep: cannot parse line: ', TRIM(line)
          CALL EXIT(1)
        END IF

        IF (TRIM(tag) .EQ. 'DEPTHVAR') THEN
          READ(line,*,IOSTAT=ios) tag, dvar
          IF (ios .NE. 0 .OR. TRIM(dvar) .NE. 'LOGRHOX') THEN
            WRITE(6,'(A,A)') ' ERROR: .dep DEPTHVAR must be LOGRHOX, got: ', &
                             TRIM(line)
            WRITE(6,'(A)')   '        The depth column is log10 column mass' // &
                             ' [g/cm^2], not an optical depth.'
            CALL EXIT(1)
          END IF
          CYCLE
        END IF

        READ(line,*,IOSTAT=ios) tag, i
        IF (ios .NE. 0) THEN
          WRITE(6,'(A,A)') ' ERROR: .dep: ', TRIM(tag) // &
                           ' needs an integer argument.'
          CALL EXIT(1)
        END IF
        SELECT CASE (TRIM(tag))
        CASE ('NTRANS');  ntr  = i
        CASE ('NDEPTH');  ndep = i
        CASE DEFAULT
          WRITE(6,'(A,A)') ' ERROR: .dep: unknown keyword ', TRIM(tag)
          CALL EXIT(1)
        END SELECT
        IF (ntr .GE. 0 .AND. ntr .NE. NLTE_NTRANS) THEN
          WRITE(6,'(A,I0,A,I0)') ' ERROR: .dep declares NTRANS = ', ntr, &
            ' but this build has NLTE_NTRANS = ', NLTE_NTRANS
          CALL EXIT(1)
        END IF
        IF (ndep .GT. MAXDEP) THEN
          WRITE(6,'(A,I0,A,I0)') ' ERROR: .dep declares NDEPTH = ', ndep, &
            ', above the compiled limit ', MAXDEP
          CALL EXIT(1)
        END IF
        CYCLE
      END IF

      ! Table body.
      IF (ntr .LT. 0 .OR. ndep .LT. 0) THEN
        WRITE(6,'(A)') ' ERROR: .dep has data before its NTRANS/NDEPTH' // &
                       ' headers.'
        CALL EXIT(1)
      END IF
      IF (n .GE. ndep) THEN
        WRITE(6,'(A,I0,A)') ' ERROR: .dep has more rows than the declared' // &
                            ' NDEPTH = ', ndep, '.'
        CALL EXIT(1)
      END IF
      READ(line,*,IOSTAT=ios) row
      IF (ios .NE. 0) THEN
        WRITE(6,'(A,I0,A,I0,A)') ' ERROR: .dep row ', n+1, ' needs ', &
          1 + 2*NLTE_NTRANS, ' numbers: log10(rhox) then b_lo b_up per transition.'
        CALL EXIT(1)
      END IF
      n = n + 1
      lrx(n) = row(1)
      DO k = 1, NLTE_NTRANS
        blo_f(k,n) = row(2*k)
        bup_f(k,n) = row(2*k + 1)
        IF (blo_f(k,n) .LE. 0.0D0 .OR. bup_f(k,n) .LE. 0.0D0) THEN
          WRITE(6,'(A,I0,A,I0)') ' ERROR: .dep has a non-positive departure' // &
            ' coefficient at row ', n, ', transition ', k
          CALL EXIT(1)
        END IF
      END DO
      IF (n .GE. 2) THEN
        IF (lrx(n) .LE. lrx(n-1)) THEN
          WRITE(6,'(A,I0,A)') ' ERROR: .dep log10(rhox) must increase' // &
            ' (surface first); row ', n, ' does not.'
          CALL EXIT(1)
        END IF
      END IF
    END DO
    CLOSE(u)

    IF (ntr .LT. 0 .OR. ndep .LT. 0) THEN
      WRITE(6,'(A)') ' ERROR: .dep is missing an NTRANS or NDEPTH header.'
      CALL EXIT(1)
    END IF
    IF (n .NE. ndep) THEN
      WRITE(6,'(A,I0,A,I0)') ' ERROR: .dep declared NDEPTH = ', ndep, &
        ' but supplied ', n
      CALL EXIT(1)
    END IF
    IF (n .LT. 2) THEN
      WRITE(6,'(A)') ' ERROR: .dep needs at least two depth rows.'
      CALL EXIT(1)
    END IF

    ! --- interpolate onto this model's layers -------------------------------
    nclamp_top = 0
    nclamp_bot = 0
    DO j = 1, nrhox
      lr = LOG10(rhox_sv(j))

      IF (lr .LE. lrx(1)) THEN
        lo = 1;  hi = 1;  t = 0.0D0
        IF (lr .LT. lrx(1)) THEN
          nclamp_top     = nclamp_top + 1
          clamped_sv(j)  = .TRUE.
        END IF
      ELSE IF (lr .GE. lrx(n)) THEN
        lo = n;  hi = n;  t = 0.0D0
        IF (lr .GT. lrx(n)) THEN
          nclamp_bot     = nclamp_bot + 1
          clamped_sv(j)  = .TRUE.
        END IF
      ELSE
        lo = 1;  hi = n
        DO WHILE (hi - lo .GT. 1)
          mid = (lo + hi) / 2
          IF (lrx(mid) .GT. lr) THEN
            hi = mid
          ELSE
            lo = mid
          END IF
        END DO
        t = (lr - lrx(lo)) / (lrx(hi) - lrx(lo))
      END IF

      ! t == 0 covers both clamped ends and an exact hit on a table node.
      ! Take the tabulated value straight through in that case rather than
      ! round-tripping it via 10**log10(b), which is not the identity in
      ! floating point.  A sidecar sampled on the model's own layers then
      ! reproduces itself exactly, which is what makes it usable as a
      ! reference against the mode-1 harness.
      IF (t .EQ. 0.0D0) THEN
        DO k = 1, NLTE_NTRANS
          blo(k,j) = blo_f(k,lo)
          bup(k,j) = bup_f(k,lo)
        END DO
      ELSE
        DO k = 1, NLTE_NTRANS
          lb = LOG10(blo_f(k,lo)) + t*(LOG10(blo_f(k,hi)) - LOG10(blo_f(k,lo)))
          blo(k,j) = 10.0D0**lb
          lb = LOG10(bup_f(k,lo)) + t*(LOG10(bup_f(k,hi)) - LOG10(bup_f(k,lo)))
          bup(k,j) = 10.0D0**lb
        END DO
      END IF
    END DO

    WRITE(6,'(A,A)')       ' NLTE: departures from ', TRIM(depfile)
    WRITE(6,'(A,I0,A,F8.4,A,F8.4)') ' NLTE:   ', n, &
      ' rows, log10(rhox) ', lrx(1), ' to ', lrx(n)
    WRITE(6,'(A,F8.4,A,F8.4)') ' NLTE:   model spans log10(rhox) ', &
      LOG10(rhox_sv(1)), ' to ', LOG10(rhox_sv(nrhox))
    IF (nclamp_top + nclamp_bot .GT. 0) THEN
      WRITE(6,'(A,I0,A,I0,A)') ' NLTE: WARNING -- grid does not cover the' // &
        ' model: ', nclamp_top, ' layers above its top, ', nclamp_bot, &
        ' below its base; endpoint b held (not extrapolated).'
    END IF

  END SUBROUTINE read_dep_file


  ! --------------------------------------------------------------------------
  !  nlte_scan_reset() / nlte_tag_at(iline)
  !
  !  Membership test for the inner line loop.  NLTE_TAG_IDX is sorted
  !  ascending (read_gfall appends as n_lte increases), and SYNTHE's line
  !  loop visits iline = 1, 2, 3, ..., so a single advancing cursor answers
  !  "is this line tagged?" in two integer compares instead of a search.
  !
  !  CONTRACT: call nlte_scan_reset() before each pass over lte_lines(), and
  !  call nlte_tag_at() with strictly increasing iline within a pass.  The
  !  cursor advances only on a hit, so CYCLEs later in the loop body are
  !  harmless -- but the call must sit ABOVE them.
  !
  !  Returns the transition index (1..NLTE_NTRANS), or 0 if not tagged.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_scan_reset()
    scan_ptr = 1
  END SUBROUTINE nlte_scan_reset

  INTEGER FUNCTION nlte_tag_at(iline)
    INTEGER, INTENT(IN) :: iline
    nlte_tag_at = 0
    IF (scan_ptr .GT. N_NLTE_TAGGED)       RETURN
    IF (NLTE_TAG_IDX(scan_ptr) .NE. iline) RETURN
    nlte_tag_at = NLTE_TAG_TRANS(scan_ptr)
    scan_ptr    = scan_ptr + 1
  END FUNCTION nlte_tag_at


  ! --------------------------------------------------------------------------
  !  nlte_factors(k, j, x, fkappa, fdev)
  !
  !  The two scalars of the header comment, for transition k at depth j.
  !
  !    x      = h nu_0 / kT at line centre (dimensionless), INPUT
  !    fkappa = kappa / kappa_LTE                            eq. (1)
  !    fdev   = (b_u - fkappa) / fkappa
  !
  !  fdev is expressed per unit of the ALREADY-SCALED opacity so the caller
  !  can reuse its profile array, which was built with the scaled kappa0:
  !  the deviation to accumulate at each grid point is just profile*fdev.
  !
  !  Population inversion (b_u/b_l > e^x) makes eq. (1) non-positive -- the
  !  line masing rather than absorbing.  Nothing downstream in SYNTHE can
  !  represent that, so the line is forced back to LTE at that depth and
  !  the occurrence counted.  With sensible photospheric departure
  !  coefficients it never triggers; if it does, the input grid is wrong.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_factors(k, j, x, fkappa, fdev)

    INTEGER, INTENT(IN)  :: k, j
    REAL(8), INTENT(IN)  :: x
    REAL(8), INTENT(OUT) :: fkappa, fdev

    REAL(8) :: bl, bu, ex, stim_lte, stim_nlte

    bl = blo(k,j)
    bu = bup(k,j)

    ex        = EXP(-x)
    stim_lte  = 1.0D0 - ex
    stim_nlte = 1.0D0 - (bu/bl)*ex

    IF (stim_nlte .LE. 0.0D0 .OR. stim_lte .LE. 0.0D0 .OR. bl .LE. 0.0D0) THEN
      inverted_sv(k,j) = .TRUE.
      fkappa = 1.0D0
      fdev   = 0.0D0
    ELSE
      fkappa = bl * stim_nlte / stim_lte
      fdev   = (bu - fkappa) / fkappa
    END IF

    ! Record for the dump.  Written on every call rather than only the
    ! first: the hyperfine components of one transition differ in x only
    ! in the seventh digit, so which one lands here is immaterial, and an
    ! unconditional store is cheaper than a test.
    seen_sv(k,j) = .TRUE.
    fkap_sv(k,j) = fkappa
    fdev_sv(k,j) = fdev
    xhnu_sv(k,j) = x

  END SUBROUTINE nlte_factors


  ! --------------------------------------------------------------------------
  !  nlte_dump(basename) -- write <basename>.nlte
  !
  !  One row per (transition, depth) carrying the model's depth scale, the
  !  departure coefficients that were installed, and the factors
  !  nlte_factors actually returned.  Small (a few hundred rows), so it is
  !  written unconditionally whenever NLTE is active rather than gated on
  !  more_output -- if NLTE is on at all, this is a developer run.
  !
  !  This is the file that makes the depth mapping observable.  Checking
  !  the dumped b against the profile it was supposed to come from -- the
  !  analytic ramp in the structured test, the input grid in mode 2 --
  !  tests the (transition, depth) indexing directly, rather than trying
  !  to infer it from an emergent flux that depends on everything at once.
  !  S_l/B_nu is written out too, from the same b and x the opacity loop
  !  used, so the source-function half is checkable independently of the
  !  opacity half.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_dump(basename)

    CHARACTER(LEN=*), INTENT(IN) :: basename

    INTEGER :: u, k, j
    REAL(8) :: ex, sratio

    IF (.NOT. nlte_on) RETURN

    OPEN(NEWUNIT=u, FILE=TRIM(basename)//'.nlte', STATUS='REPLACE', &
         ACTION='WRITE')
    WRITE(u,'(A,I0,A,I0)') '# NLTE_MODE=', NLTE_MODE, &
                           '  NLTE_TEST_SHAPE=', NLTE_TEST_SHAPE
    WRITE(u,'(A)') '# One row per (transition, depth).  inv=1 marks a depth' // &
                   ' where b_u/b_l implied a'
    WRITE(u,'(A)') '# population inversion and the line was forced back to' // &
                   ' LTE (fkappa=1, fdev=0).'
    WRITE(u,'(A)') '# seen=0 means the opacity loop never reached this' // &
                   ' transition at this depth.'
    WRITE(u,'(A)') '# clm=1 (mode 2 only) marks a depth outside the .dep' // &
                   ' table, where the endpoint'
    WRITE(u,'(A)') '# value was held rather than extrapolated.'
    WRITE(u,'(A)') '#'
    WRITE(u,'(A)') '# Values are written at full double precision on' // &
                   ' purpose: this file exists to be'
    WRITE(u,'(A)') '# differenced against the profile it came from, and a' // &
                   ' 5-digit dump would put a'
    WRITE(u,'(A)') '# 1e-6 floor under every such comparison.'
    WRITE(u,'(A)') '#'
    WRITE(u,'(A)') '#  k    j                   rhox            T' // &
                   '                   b_lo                   b_up' // &
                   '                x=hv/kT                 fkappa' // &
                   '                   fdev               S_l/B_nu' // &
                   '  seen  inv  clm'
    DO k = 1, NLTE_NTRANS
      DO j = 1, nrhox_sv
        ! S_l/B_nu from the same b and x the opacity loop used.  Equal to
        ! 1 + fdev/fkappa by construction; written explicitly because it is
        ! the quantity with physical meaning.
        ex = EXP(-xhnu_sv(k,j))
        IF (xhnu_sv(k,j) .GT. 0.0D0 .AND. bup(k,j) .GT. 0.0D0) THEN
          sratio = (1.0D0 - ex) / (blo(k,j)/bup(k,j) - ex)
        ELSE
          sratio = 1.0D0
        END IF
        WRITE(u,'(I4,1X,I4,1X,ES22.14,1X,F14.6,1X,F12.4,6(1X,ES22.14),3(1X,I5))') &
          k, j, rhox_sv(j), &
          MERGE(LOG10(MAX(tau5_sv(j),1.0D-99)), -99.0D0, tau5_sv(j) > 0.0D0), &
          temp_sv(j), blo(k,j), bup(k,j), xhnu_sv(k,j), &
          fkap_sv(k,j), fdev_sv(k,j), sratio, &
          MERGE(1, 0, seen_sv(k,j)), MERGE(1, 0, inverted_sv(k,j)), &
          MERGE(1, 0, clamped_sv(j))
      END DO
    END DO
    CLOSE(u)
    WRITE(6,'(A,A)') ' NLTE diagnostics = ', TRIM(basename)//'.nlte'

  END SUBROUTINE nlte_dump


  ! --------------------------------------------------------------------------
  !  nlte_report() -- one line of end-of-run diagnostics.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_report()
    INTEGER :: k, j, ninv, nseen
    IF (.NOT. nlte_on) RETURN
    ninv  = 0
    nseen = 0
    DO k = 1, NLTE_NTRANS
      DO j = 1, nrhox_sv
        IF (seen_sv(k,j))     nseen = nseen + 1
        IF (inverted_sv(k,j)) ninv  = ninv  + 1
      END DO
    END DO
    WRITE(6,'(A,I0,A,I0,A)') ' NLTE: departure coefficients applied at ', &
      nseen, ' of ', NLTE_NTRANS*nrhox_sv, ' (transition, depth) pairs'
    IF (ninv .GT. 0) &
      WRITE(6,'(A,I0,A)') ' NLTE: WARNING -- ', ninv, &
        ' (transition, depth) pairs implied a population inversion' // &
        ' and were forced to LTE'
  END SUBROUTINE nlte_report

END MODULE mod_nlte
