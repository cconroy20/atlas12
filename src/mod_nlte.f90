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
!   0  OFF -- pure LTE.  nlte_on stays .FALSE., nothing here allocates or
!      runs, and the line loop pays one logical compare per line.
!   1  ON -- departure coefficients interpolated from the local grid file over
!      Teff, log g, [Fe/H], v_turb and Na abundance, on this model's own
!      layers.  See interp_grid.  Data: Amarsi et al. (2020) on MARCS,
!      repackaged by Gerber et al. (2023) for Turbospectrum, reduced to one
!      file by the tools/nlte_* chain.
!
! Which transitions are eligible is declared in mod_mklinelist
! (NLTE_TR_CODE/ELO/EUP); read_gfall tags the matching lte_lines entries
! unconditionally, exactly as it records ICA4227.
! ==============================================================================

MODULE mod_nlte

  USE mod_parameters, ONLY: kw, NLTE_MODE
  USE mod_atlas_data, ONLY: DATADIR
  USE mod_mklinelist, ONLY: NLTE_NTRANS, N_NLTE_TAGGED, NLTE_NELEM, &
                            NLTE_EL_SYM, NLTE_EL_Z, NLTE_TR_EL, &
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

  ! --- The runtime grid files, one per element (NLTE_MODE = 1) -------------
  ! Built by tools/nlte_build_runtime.py from the master store.  ONE file per
  ! element: axes, then the parameter -> record index, then the records
  ! themselves, so an index can never be paired with the wrong data, and there
  ! is one thing to install per element.  Layout (little-endian):
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
  ! ONE FILE PER ELEMENT rather than one merged file: the elements come from
  ! separate published grids with their own level sets, release dates and
  ! record counts, and a synthesis that needs only Na D should not have to
  ! install Ca.  Only the elements that actually own a tagged transition in
  ! the window are opened -- see elem_used below.
  !
  ! They live flat in $ATLAS12/data/nlte/ under names carrying the element and
  ! the grid's release date.  Each has its own environment override.
  CHARACTER(LEN=32), PARAMETER :: NLTE_GRID_NAME(NLTE_NELEM) = &
      [ 'nlte/Na_MARCS_Jul-14-2023.nlte  ', &
        'nlte/Mg_MARCS_Nov-13-2024.nlte  ', &
        'nlte/Ca_MARCS_Jun-02-2021.nlte  ', &
        'nlte/Fe_MARCS_May-07-2021.nlte  ' ]
  CHARACTER(LEN=16), PARAMETER :: NLTE_GRID_ENV(NLTE_NELEM) = &
      [ 'NLTE_GRID_NA    ', 'NLTE_GRID_MG    ', 'NLTE_GRID_CA    ', &
        'NLTE_GRID_FE    ' ]

  ! Everything about one element's grid, so the five axes of one element can
  ! never be read with another's record index.
  TYPE :: grid_t
    CHARACTER(LEN=512)      :: fn    = ''
    INTEGER                 :: nT=0, nG=0, nF=0, nV=0, nD=0
    INTEGER                 :: nlev=0, ndep=0, nrec=0
    INTEGER                 :: recw=0            ! floats per record
    INTEGER(8)              :: data0=0           ! byte offset of data block
    INTEGER,    ALLOCATABLE :: levid(:)          ! atom level number per slot
    REAL(4),    ALLOCATABLE :: aT(:), aG(:), aF(:), aV(:), aD(:)
    INTEGER(4), ALLOCATABLE :: rec(:,:,:,:,:)
    INTEGER(1), ALLOCATABLE :: geo(:,:,:,:,:)
    REAL(8)                 :: ax    = 0.0D0     ! A(X) of THIS model
    ! Reported after interpolation: how much of the 32-corner weight was
    ! actually available, and whether any corner came from a spherical model.
    REAL(8)                 :: wfrac = 0.0D0
    INTEGER                 :: ncorn = 0, nsph = 0
  END TYPE grid_t

  TYPE(grid_t), SAVE :: gr(NLTE_NELEM)

  ! --- Production policy: never abort on missing DATA -------------------
  ! A grid run must not die because one model lands in a hole in a published
  ! grid.  Two things can go wrong that are properties of the data rather
  ! than of the installation, and each has a defined, reported fallback:
  !
  !   too hot        The model is hotter than the grid's last Teff node.
  !                  The grids stop at 8000 K, and a hotter model would
  !                  otherwise be handed the 8000 K departure coefficients
  !                  by the axis clamp -- silently, and badly wrong, since
  !                  Teff is the axis along which these species' departures
  !                  grow fastest and an A star is nothing like a 8000 K
  !                  one.  So this case, and only this case, reverts to LTE.
  !
  !                  Every OTHER axis clamps instead, which is the right
  !                  trade: log g 5.6 against an axis ending at 5.5, or an
  !                  abundance just off the end of the ladder, is a small
  !                  extrapolation of a slowly varying quantity, whereas
  !                  refusing it would punch holes in a production grid for
  !                  no physical reason.  Clamping is always REPORTED, per
  !                  axis and with the amount, so it can be audited.
  !
  !   empty cell     Every interpolation corner is absent.  The published
  !                  abundance ladders have holes -- at 4500 K/1.5/[Fe/H]=-2
  !                  the Mg ladder covers only [Mg/Fe] = -1.05..-0.45, so a
  !                  scaled-solar model finds nothing.  The fallback is the
  !                  nearest POPULATED point on the abundance axis alone,
  !                  reported with its offset in dex, because b depends far
  !                  more weakly on A(X) than on Teff or log g.  If even
  !                  that fails, the element reverts to LTE.
  !
  ! Reverting to LTE is exact, not approximate: b = 1 makes the factors
  ! identically 1 by the same algebra as the null test, so those lines come
  ! out bit-identical to an LTE run.  Every fallback is reported, and every
  ! model prints one greppable NLTE: STATUS line so a grid can be audited
  ! afterwards rather than trusted.
  !
  ! Set NLTE_STRICT to make the data cases fatal again while debugging.
  LOGICAL, PARAMETER :: NLTE_STRICT = .FALSE.

  ! Per-element outcome for this model, for the status line.
  LOGICAL, SAVE :: elem_lte(NLTE_NELEM) = .FALSE.
  CHARACTER(LEN=24), SAVE :: elem_stat(NLTE_NELEM) = 'unused'

  ! Which elements own a tagged transition in this window; only those are
  ! loaded.  Set by nlte_init from the tag list.
  LOGICAL, SAVE :: elem_used(NLTE_NELEM) = .FALSE.

  ! Slot within a record for each transition's lower/upper level, resolved
  ! once from levid: the transition table names ATOM levels, while which of
  ! them a given runtime file carries is a property of that file.
  INTEGER, SAVE :: gslot_lo(NLTE_NTRANS)=0, gslot_up(NLTE_NTRANS)=0

  ! Set from the run by nlte_set_params, used to pick the cell.
  REAL(8), SAVE :: m_teff=0, m_logg=0, m_feh=0, m_vturb=0

  ! Solar reference for turning the model's absolute abundances into [Fe/H].
  ! C3K models fold metallicity into the ABUNDANCE CHANGE values and leave
  ! ABUNDANCE SCALE at 1, so [Fe/H] cannot be read off the scale factor; it
  ! has to be A(Fe)_model - A(Fe)_sun.  The C3K grids are built on Grevesse &
  ! Sauval (1998) (see scripts/atlas12.sh), and the MARCS grids' own [Fe/H]
  ! axis is on a comparable scale, so GS98 is used here.  Different solar
  ! scales differ by a few 0.01 dex, which shifts the metallicity lookup
  ! slightly; A(X) itself is absolute and unaffected.
  REAL(8), PARAMETER :: A_FE_SUN = 7.50D0     ! GS98 log eps(Fe)

  ! Cursor into NLTE_TAG_IDX for nlte_tag_at(); see the contract there.
  INTEGER, SAVE :: scan_ptr = 1

  ! --- Diagnostics, dumped to <model>.nlte when nlte_on --------------------
  ! Depth structure kept for the dump, and what nlte_factors ACTUALLY
  ! computed for each (transition, depth).  Recording the factors rather
  ! than only blo/bup is the point: it makes the dump able to catch an
  ! indexing error inside nlte_factors, not just one in the grid read.
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

  ! Developer switch for bounding the truncation caveat.  The published
  ! grids' MARCS atmospheres stop short of our surfaces, so the outermost
  ! layers hold the grid's endpoint b rather than a computed one -- and the
  ! cores of strong lines form exactly there.  Setting this .TRUE. forces
  ! b = 1 (pure LTE) in every held layer instead, which is the opposite
  ! extreme: the truth lies between the two, so the pair brackets how much
  ! of a result rests on extrapolated data.  Not a physical option -- it
  ! exists to be A/B'd against the default.
  LOGICAL, PARAMETER :: NLTE_HELD_TO_LTE = .FALSE.

  ! Depths where the grid's tau range had to be held at an endpoint because
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
  SUBROUTINE nlte_init(nrhox, rhox, temp, tau5)

    INTEGER, INTENT(IN) :: nrhox
    REAL(8), INTENT(IN) :: rhox(:)   ! column mass [g/cm^2], 1 = surface
    REAL(8), INTENT(IN) :: temp(:)   ! temperature [K]
    REAL(8), INTENT(IN) :: tau5(:)   ! continuum tau_5000, same layers
    INTEGER :: i
    CHARACTER(LEN=160) :: statline

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
      WRITE(6,'(A)') ' NLTE: on, but no eligible transition lies in the' // &
        ' synthesis window; continuing in LTE.'
      RETURN
    END IF

    ! Which grids this window actually needs.  A synthesis covering only
    ! Na D must not require the Ca file to be installed.
    elem_used = .FALSE.
    DO i = 1, N_NLTE_TAGGED
      elem_used(NLTE_TR_EL(NLTE_TAG_TRANS(i))) = .TRUE.
    END DO

    CALL interp_grid(nrhox)

    ! One greppable line per model, so a grid run can be audited afterwards
    ! instead of trusted.  Every element that owns a tagged line in this
    ! window appears with what actually happened to it.
    statline = ''
    DO i = 1, NLTE_NELEM
      IF (.NOT. elem_used(i)) CYCLE
      statline = TRIM(statline) // ' ' // TRIM(NLTE_EL_SYM(i)) // '=' // &
                 TRIM(elem_stat(i))
    END DO
    WRITE(6,'(A,F7.1,A,F5.2,A,F6.2,A,A)') ' NLTE: STATUS  Teff=', m_teff, &
      ' logg=', m_logg, ' feh=', m_feh, ' |', TRIM(statline)

    ! If every element we needed fell back to LTE there is nothing left to
    ! apply; say so and switch the hot-path test off rather than running the
    ! machinery to multiply by one.
    IF (ALL(elem_lte .OR. .NOT. elem_used)) THEN
      WRITE(6,'(A)') ' NLTE: every element reverted to LTE for this model;' &
        // ' continuing in LTE.'
      nlte_on = .FALSE.
      RETURN
    END IF

    nlte_on = .TRUE.

    WRITE(6,'(A,I0,A,I0,A)') ' NLTE: on -- ', N_NLTE_TAGGED, &
      ' tagged line components over ', NLTE_NTRANS, ' transitions'
    DO i = 1, N_NLTE_TAGGED
      IF (i .LE. 12) &
        WRITE(6,'(A,I9,A,I0)') ' NLTE:   lte_lines index', NLTE_TAG_IDX(i), &
          '  transition ', NLTE_TAG_TRANS(i)
    END DO
    IF (N_NLTE_TAGGED .GT. 12) &
      WRITE(6,'(A,I0,A)') ' NLTE:   (', N_NLTE_TAGGED - 12, ' more)'

  END SUBROUTINE nlte_init




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
  !  nlte_set_params -- tell the grid lookup which model this is.
  !
  !  Abundances come in as the code stores them: ABUND(IZ) is log10(N_Z/N_tot)
  !  for Z >= 3 with XRELATIVE(IZ) an additive log offset, while ABUND(1) is
  !  the hydrogen number FRACTION.  So
  !        A(X) = ABUND(X) + XRELATIVE(X) - log10(ABUND(H)) + 12
  !  and [Fe/H] follows by subtracting the solar reference above.
  !
  !  The whole vectors are passed rather than one abundance per element:
  !  which elements have grids is a property of this module, and threading
  !  an extra argument through synthe.f90 for each new one is exactly the
  !  kind of edit that gets forgotten.
  ! --------------------------------------------------------------------------
  SUBROUTINE nlte_set_params(teff, logg, vturb, abund, xrel)
    REAL(8), INTENT(IN) :: teff, logg, vturb
    REAL(8), INTENT(IN) :: abund(:)   ! as the code stores them, index = Z
    REAL(8), INTENT(IN) :: xrel(:)    ! additive log offsets, index = Z
    INTEGER :: e, z
    m_teff  = teff
    m_logg  = logg
    m_vturb = vturb
    m_feh   = a_of(26) - A_FE_SUN
    DO e = 1, NLTE_NELEM
      z = NLTE_EL_Z(e)
      gr(e)%ax = a_of(z)
    END DO
  CONTAINS
    REAL(8) FUNCTION a_of(z)
      INTEGER, INTENT(IN) :: z
      a_of = abund(z) + xrel(z) - LOG10(abund(1)) + 12.0D0
    END FUNCTION a_of
  END SUBROUTINE nlte_set_params


  SUBROUTINE read_grid_index(e)

    INTEGER, INTENT(IN) :: e

    INTEGER :: u, hdr(8), n, k, i
    LOGICAL :: ex
    CHARACTER(LEN=512) :: fn
    INTEGER :: envlen, envstat
    CHARACTER(LEN=512) :: envval
    INTEGER(8) :: fsize, fwant

    IF (gr(e)%nT .GT. 0) RETURN                 ! already loaded

    CALL GET_ENVIRONMENT_VARIABLE(TRIM(NLTE_GRID_ENV(e)), envval, envlen, &
                                  envstat)
    IF (envstat .EQ. 0 .AND. envlen .GT. 0) THEN
      fn = TRIM(envval)
    ELSE
      fn = TRIM(DATADIR) // TRIM(NLTE_GRID_NAME(e))
    END IF
    gr(e)%fn = fn

    INQUIRE(FILE=fn, EXIST=ex)
    IF (.NOT. ex) THEN
      WRITE(6,'(A,A,A)') ' ERROR: NLTE needs the runtime grid file for ', &
                         NLTE_EL_SYM(e), '; not found:'
      WRITE(6,'(A,A)')   '        ', TRIM(fn)
      WRITE(6,'(A)')     '        Build it with tools/nlte_extract_grid.py,' &
                         // ' nlte_build_index.py and nlte_build_runtime.py'
      WRITE(6,'(A,A,A)') '        into $ATLAS12/data/nlte/, or point $', &
                         TRIM(NLTE_GRID_ENV(e)), ' at it.'
      CALL EXIT(1)
    END IF

    OPEN(NEWUNIT=u, FILE=fn, ACCESS='STREAM', FORM='UNFORMATTED', &
         STATUS='OLD', ACTION='READ')
    READ(u) hdr
    gr(e)%nT = hdr(1);  gr(e)%nG = hdr(2);  gr(e)%nF = hdr(3)
    gr(e)%nV = hdr(4);  gr(e)%nD = hdr(5)
    gr(e)%nlev = hdr(6);  gr(e)%ndep = hdr(7);  gr(e)%nrec = hdr(8)
    gr(e)%recw = gr(e)%ndep + gr(e)%nlev*gr(e)%ndep

    ALLOCATE(gr(e)%aT(gr(e)%nT), gr(e)%aG(gr(e)%nG), gr(e)%aF(gr(e)%nF), &
             gr(e)%aV(gr(e)%nV), gr(e)%aD(gr(e)%nD), gr(e)%levid(gr(e)%nlev))
    READ(u) gr(e)%aT
    READ(u) gr(e)%aG
    READ(u) gr(e)%aF
    READ(u) gr(e)%aV
    READ(u) gr(e)%aD
    READ(u) gr(e)%levid
    ALLOCATE(gr(e)%rec(gr(e)%nT,gr(e)%nG,gr(e)%nF,gr(e)%nV,gr(e)%nD), &
             gr(e)%geo(gr(e)%nT,gr(e)%nG,gr(e)%nF,gr(e)%nV,gr(e)%nD))
    READ(u) gr(e)%rec
    READ(u) gr(e)%geo
    INQUIRE(UNIT=u, POS=gr(e)%data0)   ! first byte of the data block
    CLOSE(u)

    ! Does the file actually hold what its header claims?  A MISSING grid
    ! announces itself; a TRUNCATED one -- an interrupted copy of a 787 MB
    ! file -- does not, and would either read garbage or die deep inside the
    ! interpolation with an unhelpful message.  The header fixes the size
    ! exactly, so check it once, here, where the message can still be useful.
    INQUIRE(FILE=fn, SIZE=fsize)
    fwant = gr(e)%data0 - 1_8 + INT(gr(e)%nrec, 8)*INT(gr(e)%recw, 8)*4_8
    IF (fsize .NE. fwant) THEN
      WRITE(6,'(A,A,A)') ' ERROR: the NLTE grid file for ', NLTE_EL_SYM(e), &
        ' is not the size its own header implies:'
      WRITE(6,'(A,A)')   '        ', TRIM(fn)
      WRITE(6,'(A,I0,A,I0,A)') '        header implies ', fwant, &
        ' bytes, file is ', fsize, ' bytes.'
      WRITE(6,'(A)') '        Most likely an interrupted copy or a partial' &
        // ' download; re-install it.'
      CALL EXIT(1)
    END IF

    ! Resolve this element's transitions to slots in the record.  The
    ! transition table names ATOM levels; which of them this file happens to
    ! carry is a property of the file, so a runtime file built with a
    ! narrower level set must fail loudly rather than read the wrong slot.
    DO k = 1, NLTE_NTRANS
      IF (NLTE_TR_EL(k) .NE. e) CYCLE
      gslot_lo(k) = 0;  gslot_up(k) = 0
      DO i = 1, gr(e)%nlev
        IF (gr(e)%levid(i) .EQ. NLTE_TR_LEVLO(k)) gslot_lo(k) = i
        IF (gr(e)%levid(i) .EQ. NLTE_TR_LEVUP(k)) gslot_up(k) = i
      END DO
      IF (gslot_lo(k) .EQ. 0 .OR. gslot_up(k) .EQ. 0) THEN
        WRITE(6,'(A,I0,A,I0,A,I0,A)') ' ERROR: transition ', k, ' needs atom' &
          // ' levels ', NLTE_TR_LEVLO(k), ' and ', NLTE_TR_LEVUP(k), &
          ' but the runtime grid file does not carry both.'
        WRITE(6,'(A)') '        Rebuild it with a wider --levels set.'
        CALL EXIT(1)
      END IF
    END DO

    n = COUNT(gr(e)%rec .GT. 0)
    WRITE(6,'(A,A,A,A)') ' NLTE: grid  ', NLTE_EL_SYM(e), '  ', TRIM(fn)
    WRITE(6,'(A,I0,4(A,I0),3(A,I0),A)') &
      ' NLTE:   axes ', gr(e)%nT,'x',gr(e)%nG,'x',gr(e)%nF,'x',gr(e)%nV, &
      'x',gr(e)%nD, ', ', gr(e)%nlev, ' levels, ', gr(e)%nrec, ' records, ', &
      NINT(100.0*REAL(n)/REAL(SIZE(gr(e)%rec))), ' % of cells filled'

  END SUBROUTINE read_grid_index


  ! --------------------------------------------------------------------------
  !  bracket(ax, n, x, i, w) -- locate x on a sorted axis.
  !
  !  Returns the lower node i and the weight w of node i+1, both CLAMPED to
  !  the axis ends: outside the grid the edge value is held, never
  !  extrapolated, as in depth.
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
  !  interp_grid(nrhox) -- install b for this model.
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

    INTEGER :: iT,iG,iF,iV,iD, cT,cG,cF,cV,cD, k, j, u, r, e, ic, nc
    REAL(8) :: wT,wG,wF,wV,wD, wgt, wsum
    REAL(4), ALLOCATABLE :: buf(:)
    REAL(8), ALLOCATABLE :: ltau(:), lb(:)
    REAL(8) :: lt, t, dxoff
    REAL(8) :: acc_lo(NLTE_NTRANS,kw), acc_up(NLTE_NTRANS,kw)
    INTEGER :: lo, hi, mid, lev, nd
    INTEGER :: crec(32), cgeo(32), jT,jG,jF,jV,jD, off, sgn
    CHARACTER(LEN=96) :: clampnote
    REAL(8) :: cwgt(32)
    LOGICAL :: found

    acc_lo = 0.0D0;  acc_up = 0.0D0
    elem_lte  = .FALSE.
    elem_stat = 'unused'

    ! One element at a time: each has its own axes, its own record index and
    ! its own abundance ladder, and only the elements with a tagged
    ! transition in this window are touched at all.
    DO e = 1, NLTE_NELEM

      IF (.NOT. elem_used(e)) CYCLE

      CALL read_grid_index(e)

      ! --- hotter than the grid?  that one is not clampable ---------------
      IF (m_teff .GT. DBLE(gr(e)%aT(gr(e)%nT))) THEN
        CALL data_problem(e, 'too-hot', &
          ' NLTE: ' // NLTE_EL_SYM(e) // ' -- model is hotter than the' &
          // ' grid''s last Teff node; this element reverts to LTE rather' &
          // ' than clamp (departures grow fastest along Teff).')
        CYCLE
      END IF

      ! --- every other axis clamps, but says so ---------------------------
      clampnote = ''
      CALL note_clamp(gr(e)%aT, gr(e)%nT, m_teff,           'Teff', clampnote)
      CALL note_clamp(gr(e)%aG, gr(e)%nG, m_logg,           'logg', clampnote)
      CALL note_clamp(gr(e)%aF, gr(e)%nF, m_feh,            'feh',  clampnote)
      CALL note_clamp(gr(e)%aV, gr(e)%nV, m_vturb,          'vtb',  clampnote)
      CALL note_clamp(gr(e)%aD, gr(e)%nD, gr(e)%ax - m_feh, 'dX',   clampnote)
      IF (LEN_TRIM(clampnote) .GT. 0) &
        WRITE(6,'(A,A,A,A)') ' NLTE: WARNING -- ', NLTE_EL_SYM(e), &
          ': outside the grid, held at the edge for', TRIM(clampnote)

      CALL bracket(gr(e)%aT, gr(e)%nT, m_teff,             iT, wT)
      CALL bracket(gr(e)%aG, gr(e)%nG, m_logg,             iG, wG)
      CALL bracket(gr(e)%aF, gr(e)%nF, m_feh,              iF, wF)
      CALL bracket(gr(e)%aV, gr(e)%nV, m_vturb,            iV, wV)
      CALL bracket(gr(e)%aD, gr(e)%nD, gr(e)%ax - m_feh,   iD, wD)

      ! --- PHASE 1: collect the populated corners and their weights -------
      nc = 0
      dxoff = 0.0D0
      DO cT = 0, MIN(1, gr(e)%nT-1)
       DO cG = 0, MIN(1, gr(e)%nG-1)
        DO cF = 0, MIN(1, gr(e)%nF-1)
         DO cV = 0, MIN(1, gr(e)%nV-1)
          DO cD = 0, MIN(1, gr(e)%nD-1)
            wgt = (MERGE(wT, 1.0D0-wT, cT.EQ.1)) * (MERGE(wG, 1.0D0-wG, cG.EQ.1)) &
                * (MERGE(wF, 1.0D0-wF, cF.EQ.1)) * (MERGE(wV, 1.0D0-wV, cV.EQ.1)) &
                * (MERGE(wD, 1.0D0-wD, cD.EQ.1))
            IF (wgt .LE. 0.0D0) CYCLE
            r = gr(e)%rec(iT+cT, iG+cG, iF+cF, iV+cV, iD+cD)
            IF (r .LE. 0) CYCLE                    ! absent corner
            nc = nc + 1
            crec(nc) = r
            cwgt(nc) = wgt
            cgeo(nc) = gr(e)%geo(iT+cT, iG+cG, iF+cF, iV+cV, iD+cD)
          END DO
         END DO
        END DO
       END DO
      END DO

      ! --- PHASE 2: nothing populated?  step along the abundance axis -----
      IF (nc .EQ. 0) THEN
        jT = MIN(iT + NINT(wT), gr(e)%nT)
        jG = MIN(iG + NINT(wG), gr(e)%nG)
        jF = MIN(iF + NINT(wF), gr(e)%nF)
        jV = MIN(iV + NINT(wV), gr(e)%nV)
        found = .FALSE.
        DO off = 0, gr(e)%nD - 1
          DO sgn = 1, -1, -2
            jD = iD + sgn*off
            IF (jD .LT. 1 .OR. jD .GT. gr(e)%nD) CYCLE
            IF (gr(e)%rec(jT,jG,jF,jV,jD) .GT. 0) THEN
              nc = 1
              crec(1) = gr(e)%rec(jT,jG,jF,jV,jD)
              cwgt(1) = 1.0D0
              cgeo(1) = gr(e)%geo(jT,jG,jF,jV,jD)
              dxoff = DBLE(gr(e)%aD(jD)) - (gr(e)%ax - m_feh)
              found = .TRUE.
              EXIT
            END IF
            IF (off .EQ. 0) EXIT
          END DO
          IF (found) EXIT
        END DO
        IF (.NOT. found) THEN
          CALL data_problem(e, 'no-cell', &
            ' NLTE: ' // NLTE_EL_SYM(e) // ' -- no populated grid cell' &
            // ' anywhere on the abundance axis here; reverts to LTE.')
          CYCLE
        END IF
        WRITE(6,'(A,A,A,F6.2,A)') ' NLTE: WARNING -- ', NLTE_EL_SYM(e), &
          ': the grid has no cell at this abundance; using the nearest' &
          // ' populated one, offset ', dxoff, ' dex.'
      END IF

      ! --- PHASE 3: read the corners and put them on our depth scale ------
      wsum = 0.0D0
      gr(e)%nsph = 0;  gr(e)%ncorn = 0
      nd = gr(e)%ndep

      ALLOCATE(buf(gr(e)%recw), ltau(nd), lb(nd))
      OPEN(NEWUNIT=u, FILE=TRIM(gr(e)%fn), ACCESS='STREAM', &
           FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')

      DO ic = 1, nc
        r   = crec(ic)
        wgt = cwgt(ic)
        gr(e)%ncorn = gr(e)%ncorn + 1
        IF (cgeo(ic) .EQ. 2) gr(e)%nsph = gr(e)%nsph + 1

        READ(u, POS=gr(e)%data0 + INT(r-1, 8)*INT(gr(e)%recw, 8)*4_8) buf
        DO j = 1, nd
          ltau(j) = LOG10(MAX(DBLE(buf(j)), 1.0D-99))
        END DO

        DO k = 1, NLTE_NTRANS
          IF (NLTE_TR_EL(k) .NE. e) CYCLE
          DO lev = 1, 2
            IF (lev .EQ. 1) THEN
              DO j = 1, nd
                lb(j) = LOG10(MAX(DBLE(buf(nd + (gslot_lo(k)-1)*nd + j)), 1.0D-99))
              END DO
            ELSE
              DO j = 1, nd
                lb(j) = LOG10(MAX(DBLE(buf(nd + (gslot_up(k)-1)*nd + j)), 1.0D-99))
              END DO
            END IF
            DO j = 1, nrhox
              lt = LOG10(MAX(tau5_sv(j), 1.0D-99))
              ! Mark layers outside the corner's own tau range.  This is
              ! reported and is how the grid's truncation became visible
              ! at all; it must not be lost.
              IF (lt .LT. ltau(1) .OR. lt .GT. ltau(nd)) clamped_sv(j) = .TRUE.
              IF (lt .LE. ltau(1)) THEN
                t = lb(1)
              ELSE IF (lt .GE. ltau(nd)) THEN
                t = lb(nd)
              ELSE
                lo = 1;  hi = nd
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
      CLOSE(u)
      DEALLOCATE(buf, ltau, lb)

      gr(e)%wfrac = wsum
      DO k = 1, NLTE_NTRANS
        IF (NLTE_TR_EL(k) .NE. e) CYCLE
        DO j = 1, nrhox
          blo(k,j) = 10.0D0**(acc_lo(k,j)/wsum)
          bup(k,j) = 10.0D0**(acc_up(k,j)/wsum)
          ! Bounding experiment: drop the held layers to LTE instead of
          ! holding the grid's edge value.  See NLTE_HELD_TO_LTE above.
          IF (NLTE_HELD_TO_LTE .AND. clamped_sv(j)) THEN
            blo(k,j) = 1.0D0
            bup(k,j) = 1.0D0
          END IF
        END DO
      END DO

      IF (dxoff .NE. 0.0D0) THEN
        WRITE(elem_stat(e),'(A,SP,F0.2,A)') 'nearest-dX(', dxoff, ')'
      ELSE
        WRITE(elem_stat(e),'(A,F0.3,A)') 'ok(w=', gr(e)%wfrac, ')'
      END IF

      WRITE(6,'(A,A,A,F7.1,A,F5.2,A,F6.2,A,F4.1,A,A,A,F5.2)') &
        ' NLTE: ', NLTE_EL_SYM(e), ' model  Teff=', m_teff, '  logg=', m_logg, &
        '  [Fe/H]=', m_feh, '  vturb=', m_vturb, '  A(', NLTE_EL_SYM(e), &
        ')=', gr(e)%ax
      WRITE(6,'(A,I3,A,F6.3,A,I0,A)') ' NLTE:   ', gr(e)%ncorn, &
        ' corners used, weight fraction ', gr(e)%wfrac, ', ', gr(e)%nsph, &
        ' spherical'
      IF (gr(e)%wfrac .LT. 0.999D0) &
        WRITE(6,'(A,A,A,F6.3,A)') ' NLTE: WARNING -- ', NLTE_EL_SYM(e), &
          ': only ', gr(e)%wfrac, ' of the interpolation weight was' &
          // ' populated; the rest was dropped and the remainder renormalised.'
      IF (gr(e)%nsph .GT. 0 .AND. gr(e)%nsph .LT. gr(e)%ncorn) &
        WRITE(6,'(A,A,A)') ' NLTE: WARNING -- ', NLTE_EL_SYM(e), &
          ': corners mix spherical and plane-parallel MARCS geometry.'

    END DO

    IF (COUNT(clamped_sv(1:nrhox)) .GT. 0) &
      WRITE(6,'(A,I0,A,I0,A)') ' NLTE: WARNING -- grid does not cover the' // &
        ' model in depth: ', COUNT(clamped_sv(1:nrhox)), ' of ', nrhox, &
        ' layers hold an endpoint b (not extrapolated).'

  END SUBROUTINE interp_grid


  ! --------------------------------------------------------------------------
  !  note_clamp(ax, n, x, name, note) -- record that x fell off this axis.
  !
  !  bracket() clamps silently.  Silence is the problem, not the clamping:
  !  holding the edge value is a reasonable thing to do for log g, [Fe/H],
  !  v_turb and abundance, and an unreasonable one for Teff above the grid,
  !  which is handled separately.  This appends "<axis><signed amount>" to a
  !  note so that whatever clamping did happen appears in the log.
  ! --------------------------------------------------------------------------
  SUBROUTINE note_clamp(ax, n, x, name, note)
    REAL(4),          INTENT(IN)    :: ax(:)
    INTEGER,          INTENT(IN)    :: n
    REAL(8),          INTENT(IN)    :: x
    CHARACTER(LEN=*), INTENT(IN)    :: name
    CHARACTER(LEN=*), INTENT(INOUT) :: note
    REAL(8) :: d
    CHARACTER(LEN=32) :: piece
    IF (n .EQ. 1) RETURN              ! a collapsed axis constrains nothing
    d = 0.0D0
    IF (x .LT. DBLE(ax(1))) d = x - DBLE(ax(1))
    IF (x .GT. DBLE(ax(n))) d = x - DBLE(ax(n))
    IF (d .EQ. 0.0D0) RETURN
    WRITE(piece,'(1X,A,SP,G0.3)') TRIM(name), d
    note = TRIM(note) // TRIM(piece)
  END SUBROUTINE note_clamp


  ! --------------------------------------------------------------------------
  !  data_problem(e, tag, msg) -- this element cannot be done for this model.
  !
  !  Reverts it to LTE and records why.  b stays 1, which makes the opacity
  !  and source-function factors identically 1 -- so those lines come out
  !  bit-identical to an LTE run, not approximately so.  Fatal only under
  !  NLTE_STRICT, which exists for debugging a grid build.
  ! --------------------------------------------------------------------------
  SUBROUTINE data_problem(e, tag, msg)
    INTEGER,          INTENT(IN) :: e
    CHARACTER(LEN=*), INTENT(IN) :: tag, msg
    WRITE(6,'(A)') TRIM(msg)
    WRITE(6,'(A,F8.1,A,F6.2,A,F6.2,A,F5.1,A,A,A,F6.2)') &
      '        Teff=', m_teff, ' logg=', m_logg, ' [Fe/H]=', m_feh, &
      ' vturb=', m_vturb, ' A(', NLTE_EL_SYM(e), ')=', gr(e)%ax
    IF (NLTE_STRICT) THEN
      WRITE(6,'(A)') '        NLTE_STRICT is set: stopping.'
      CALL EXIT(1)
    END IF
    elem_lte(e)  = .TRUE.
    elem_stat(e) = 'LTE:' // tag
  END SUBROUTINE data_problem



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
  !  the input grid --
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
    WRITE(u,'(A,I0)') '# NLTE_MODE=', NLTE_MODE
    WRITE(u,'(A)') '# One row per (transition, depth).  inv=1 marks a depth' // &
                   ' where b_u/b_l implied a'
    WRITE(u,'(A)') '# population inversion and the line was forced back to' // &
                   ' LTE (fkappa=1, fdev=0).'
    WRITE(u,'(A)') '# seen=0 means the opacity loop never reached this' // &
                   ' transition at this depth.'
    WRITE(u,'(A)') '# clm=1 marks a depth outside the grid tau range,' // &
                   ' where the endpoint'
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
