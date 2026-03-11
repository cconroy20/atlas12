! synthe_spectrv.f90 -- Merged XNFPELSYN + SYNTHE + SPECTRV main program
! Requires: atlas12_modules.f90, synthe_module.f90 (compile in that order)
!
! ============================================================================
!  PROGRAM SYNTHE_SPECTRV
!
!  Integrates the three formerly separate programs XNFPELSYN, SYNTHE, and
!  SPECTRV into a single executable, eliminating all intermediate file I/O:
!
!    xnfpelsyn  wrote fort.10 (model atmosphere + continuum opacities)
!    synthe     read  fort.10; wrote fort.9 (opacity vectors + line data)
!    spectrv    read  fort.9 and fort.10; wrote fort.7 (emergent spectrum)
!
!  All three data flows are now in-memory:
!    run_xnfpelsyn() fills module arrays (continuum opacities, ion pops)
!    asynth_sv(nrhox) carries one wavelength point's opacity vector
!    lindat8_sv/lindat4_sv/alinec_sv carry line-centre records
!
!  ATLAS library routines (JOSH, READIN, KAPP, etc.) are accessed via
!  USE mod_atlas_data from atlas12_modules.f90.  All former COMMON block
!  interfaces have been replaced with module variable access.
!
!  Unit assignments:
!    5   input : fort.5 model atmosphere control cards (read by READIN)
!   11   output: ASCII spectrum (wavelength in Angstroms, flux, continuum flux)
!   12   input : preprocessed LTE line data (unformatted sequential)
!   13   scratch: merged line archive (deleted after use)
!   14   input : preprocessed LTE line data (fort.14, deleted after merge)
!   15   scratch: (ILINE, KAPCEN) pairs for line-centre tracking (deleted)
!   16   output: line identification records (formatted)
!   17   input : continua.dat continuum edge frequency list (from DATADIR)
!   33   output: wavelength / flux / continuum / optical-depth table
!   19   input : NLTE line data from RNLTE
!   20   input : additional NLTE line data
!   25   input : spectrv.input run parameters (rhoxj, r1, r101, fudge exponents)
!   93   input : run parameters from SYNBEG (unformatted, deleted after read)
!
!  ASCII spectrum output format (fort.11):
!    One line per wavelength point:
!      col  1-11 : wavelength (Angstroms), F11.4
!      col 12-26 : flux or specific intensity, E15.6
!      col 27-41 : continuum flux, E15.6
!    Replaces the standalone SYNTOASCANG post-processing step.
!
!  Units REMOVED relative to the three originals:
!    8   -- was: n9 count SYNTHE->SPECTRV             (now: nlines_sv)
!    9   -- was: opacity vectors + line data           (now: in-memory)
!   10   -- was: fort.10 model atmosphere file         (now: module arrays)
!    7   -- was: binary spectrum for PLOTSY              (now: ASCII fort.11)
!   14*  -- was: direct-access opacity matrix scratch  (now: in-memory arrays)
!   28   -- was: XNFPEL diagnostic dump               (removed)
!   29   -- was: opacity spectrum diagnostic           (removed)
!   (syntoascanga post-processing step eliminated; ASCII output now inline)
! ============================================================================
PROGRAM synthe_spectrv
  USE synthe_module
  USE mod_atlas_data, only: &
    ! Procedures
    JOSH, READIN, KAPP, &
    ! Data directory path
    DATADIR, &
    ! Renamed imports (collide with synthe_module or local names)
    bhyd_a => BHYD, hckt_a => HCKT, hkt_a => HKT, itemp_a => ITEMP, &
    nrhox_a => NRHOX, rhox_a => RHOX, vturb_a => VTURB, &
    ahline_a => AHLINE, shline_a => SHLINE, &
    teff_a => TEFF, glog_a => GLOG, title_atlas => TITLE, &
    ! Direct imports (no name collision)
    ACONT, ALINE, ANGLE, AXLINE, BNU, DELTAW, EHVKT, FREQ, FREQLG, &
    HNU, IFPRES, IFSCAT, IFSURF, NMU, NUHI, NULO, NUMNU, &
    PRDDOP, PRDPOW, SCONT, SIGMAC, SIGMAL, SIGPRD, SLINE, STIM, &
    SURFI, SXLINE, TAUNU
  IMPLICIT NONE


  ! ATLAS variables are now accessed via USE mod_atlas_data above.
  ! (Formerly declared via COMMON blocks; removed in F90 modernisation.)

  ! --- LINDAT variables (local, for binary I/O and line identification) ---
  !
  ! lindat8_c(14) and lindat4_c(28) are the I/O buffers for unformatted
  ! READ/WRITE of line data records.  The named scalars below are unpacked
  ! from these arrays after each READ.
  !
  ! The original F77 code used EQUIVALENCE to overlay the scalars on the
  ! arrays, but this only works reliably in COMMON blocks.  With local
  ! variables, gfortran does not guarantee contiguous layout, so the
  ! overlay silently fails for all fields except the first.
  !
  ! lindat8 layout (14 words, REAL*8):
  !   1: wl   2: e   3: ep   4-5: label(2)   6-7: labelp(2)
  !   8-9: other1(2)   10-11: other2(2)   12: wlvac   13: center   14: concen
  !
  ! lindat4 layout (28 words, REAL*4 / INTEGER):
  !   1: nelion(I)  2: nblo(I)  3: nbup(I)  4: iso1(I)  5: iso2(I)
  !   6: gammar  7: gammas  8: gammaw  9: ref  10: x1  11: x2
  !   12: gflog  13: xj  14: xjp  15: code  16: elo  17: gf
  !   18: gs  19: gr  20: gw
  !   21: dwl  22: dgflog  23: dgammar  24: dgammas  25: dgammaw
  !   26: extra1  27: extra2  28: extra3
  REAL(8)  :: wl_c, e_c, ep_c, wlvac_c, center_c, concen_c
  REAL(8)  :: label_c(2), labelp_c(2), other1_c(2), other2_c(2)
  INTEGER  :: nelion_c, nblo_c, nbup_c, iso1_c, iso2_c
  REAL(4)  :: gammar_c, gammas_c, gammaw_c, ref_c, x1_c, x2_c
  REAL(4)  :: gflog_c, xj_c, xjp_c, code_c, elo_c, gf_c, gs_c, gr_c, gw_c
  REAL(4)  :: dwl_c, dgflog_c, dgammar_c, dgammas_c, dgammaw_c
  REAL(4)  :: extra1_c, extra2_c, extra3_c
  REAL(4)  :: alinec_c(kw)
  REAL(8)  :: lindat8_c(14)
  REAL(4)  :: lindat4_c(28)

  ! --- Variables formerly in ATLAS COMMON blocks, not in mod_atlas_data ---
  REAL(8) :: bone(kw)
  REAL(8) :: wlbeg_loc, wlend_loc, wbegin_loc
  REAL(8) :: air_loc, vt_loc
  INTEGER :: ifvac_loc, n10_loc, nmu2_loc

  ! ==========================================================================
  !  LOCAL PARAMETERS (aliases for synthe_module constants)
  ! ==========================================================================
  INTEGER, PARAMETER :: kw_p      = kw
  INTEGER, PARAMETER :: mw_p      = mw
  INTEGER, PARAMETER :: mw6_p     = mw6
  INTEGER, PARAMETER :: maxbuff_p = MAXBUFF
  INTEGER, PARAMETER :: maxlin_p  = MAXLIN

  ! ==========================================================================
  !  SYNTHE LOCAL ARRAYS
  ! ==========================================================================
  INTEGER(4) :: line_flag(maxlin_p)

  ! In-memory opacity matrices (replace unit 14 direct-access scratch file)
  REAL(4), ALLOCATABLE :: opacity_matrix(:,:)   ! (length, kw) wavelength × depth
  REAL(4), ALLOCATABLE :: linecen_matrix(:,:)   ! (n9, kw) line × depth

  ! --- Run parameters (unit 93) ---
  INTEGER   :: nlines_in, ifvac, ifnlte, n19, ifpred, linout
  REAL(4)   :: turbv, cutoff
  REAL(4)   :: deckj(7, kw_p)
  REAL(4)   :: velshift_arr(kw_p)
  REAL(4)   :: hfield(kw_p)

  ! --- Model atmosphere header (unit 10) ---
  INTEGER   :: nedge, ncon_tape
  REAL(8)   :: teff, glog
  REAL(8)   :: title(74)
  REAL(8)   :: frqedg(377), wledge(377), cmedge(377), confrq(1131)
  REAL(8)   :: deledge(377), halfedge(377)
  REAL(8)   :: idmol(mm), momass(mm)   ! mm=mw-39=100 molecules

  ! --- REAL*8 atmosphere read buffers ---
  REAL(8)   :: qt(kw_p), qtkev(kw_p), qtk(kw_p), qhkt(kw_p)
  REAL(8)   :: qtlog(kw_p), qhckt(kw_p)
  REAL(8)   :: qp(kw_p), qxne(kw_p), qxnatom(kw_p), qrho(kw_p)
  REAL(8)   :: qrhox(kw_p), qvturb(kw_p)
  REAL(8)   :: qxnfh(kw_p), qxnfhe(kw_p,2), qxnfh2(kw_p)
  REAL(8)   :: qdopple(mw6_p), qxnfpel(mw6_p)
  REAL(8)   :: qablog(1131)
  REAL(8)   :: qcongf

  REAL(4)   :: ablog(3, 377)

  REAL(4)   :: asynth(kw_p), alinec(kw_p)
  REAL(4)   :: asyncont(kw_p), alinecont(kw_p)
  INTEGER   :: mlinej(kw_p)
  INTEGER   :: iftp(kw_p), abmin_arr(kw_p)

  REAL(4)   :: hfac(kw_p), hefac(kw_p), h2fac(kw_p)

  REAL(8)   :: lindat8(14)
  REAL(4)   :: lindat4(28)

  ! ==========================================================================
  !  SYNTHE SCALAR WORK VARIABLES
  ! ==========================================================================
  INTEGER   :: i, j, l, n, nu, iedge, nbuff_i
  INTEGER   :: n12, nlines, iline, ilines, n191
  INTEGER   :: n9, ncen
  INTEGER   :: nvshift
  INTEGER   :: minred, maxred, maxblue, minblue_i, ibuff, istep
  INTEGER   :: n10dop, nstep, n1, maxstep
  INTEGER   :: nbuff_s, congf_nel
  INTEGER   :: maxline, i9, nelem_i
  REAL(4)   :: kappa0_s, kapmin_s, kapcen_s, kappa_s
  REAL(4)   :: adamp_s
  REAL(4)   :: congf_s, alpha_s, v2_s
  REAL(4)   :: gamrf, gamsf, gamwf
  REAL(4)   :: vsteps, tabstep, tabi, dvoigt
  REAL(4)   :: x_wing
  REAL(8)   :: wave8, freq8
  REAL(4)   :: elo_s, dopple_nel

  ! ==========================================================================
  !  SPECTRV LOCAL VARIABLES
  ! ==========================================================================

  ! Continuum opacity read buffers (unit 10 SPECTRV pass)
  REAL(8)   :: qcontabs(3 * (medge-1))
  REAL(8)   :: qcontscat(3 * (medge-1))

  ! Spectrum synthesis output
  REAL(8)   :: q(41)

  ! ASCII plot array
  CHARACTER(LEN=1) :: aplot(101)

  ! SPECTRV run parameters
  REAL(8)   :: rhoxj, r1, r101, ph1, pc1, psi1
  REAL(8)   :: slope, origin

  ! Wavelength loop variables (SPECTRV half)
  REAL(8)   :: wave_sv   ! current synthesis wavelength (nm)
  REAL(8)   :: wavold    ! previous line-centre wavelength (for bracket reset)
  REAL(8)   :: slinec    ! line source function contribution
  REAL(8)   :: resid     ! residual intensity

  ! Integer loop/index variables (SPECTRV half)
  INTEGER   :: mu, iresid, iplot, nlines_rd, n910
  INTEGER   :: iedge_sv  ! continuum edge bracket for SPECTRV (separate from SYNTHE's iedge)

  ! Flags
  INTEGER   :: ifvac_sv   ! local copy of ifvac for SPECTRV /TRASH/ and title

  ! ==========================================================================
  !  INITIALISE DATA DIRECTORY
  !
  !  Locate data files via $ATLAS12 environment variable.
  !  If unset, defaults to ./data/
  ! ==========================================================================
  BLOCK
    CHARACTER(256) :: envval
    INTEGER :: envlen, envstat
    CALL GET_ENVIRONMENT_VARIABLE('ATLAS12', envval, envlen, envstat)
    IF (envstat == 0 .AND. envlen > 0) THEN
      DATADIR = TRIM(envval)
      IF (DATADIR(envlen:envlen) /= '/') DATADIR = TRIM(DATADIR) // '/'
      DATADIR = TRIM(DATADIR) // 'data/'
    ELSE
      DATADIR = 'data/'
    END IF
    WRITE(6, '(A,A)') ' DATADIR = ', TRIM(DATADIR)
  END BLOCK

  ! ==========================================================================
  !  SECTION 1.  READ RUN PARAMETERS FROM UNIT 93
  ! ==========================================================================
  OPEN(UNIT=93, FORM='UNFORMATTED')
  READ(93) nlines_in, length, ifvac, ifnlte, n19, turbv, deckj, ifpred, &
       wlbeg, wlend, resolu, ratio, ratiolg, cutoff, linout
  WRITE(6,*) 'SYNTHE INPUT (UNIT 93): '
  WRITE(6,*) 'nlines,length,ifvac,ifnlte,n19,ifpred,', &
             'wlbeg,wlend,resolu,cutoff,linout'
  WRITE(6,'(I8,I6,1x,I1,1x,I1,1x,I5,1x,I1,1x,F7.1,1x,F7.1,1x,F9.1,1x,ES12.3,1x,I8)') &
       nlines_in, length, ifvac, ifnlte, n19, ifpred, &
       wlbeg, wlend, resolu, cutoff, linout
  CLOSE(UNIT=93, STATUS='DELETE')

  ixwlbeg = INT(LOG(wlbeg) / ratiolg)
  wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  IF (wbegin < wlbeg) THEN
    ixwlbeg = ixwlbeg + 1
    wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  END IF

  ! ==========================================================================
  !  SECTION 2.  OPEN FILES
  ! ==========================================================================
  ! Unit 10 (fort.10) eliminated -- data now in module arrays from run_xnfpelsyn()
  OPEN(UNIT=12, STATUS='OLD',     FORM='UNFORMATTED', POSITION='APPEND')
  OPEN(UNIT=13, STATUS='SCRATCH', FORM='UNFORMATTED')
  OPEN(UNIT=14, STATUS='OLD',     FORM='UNFORMATTED', POSITION='APPEND')
  OPEN(UNIT=19, STATUS='OLD',     FORM='UNFORMATTED', POSITION='APPEND')
  OPEN(UNIT=20, STATUS='OLD',     FORM='UNFORMATTED', POSITION='APPEND')

  ! ==========================================================================
  !  SECTION 3.  MERGE LINE ARCHIVES INTO UNIT 13
  ! ==========================================================================
  IF (linout >= 0) THEN
    IF (n19 > 0) THEN
      REWIND 20
      DO i = 1, n19
        READ(20)  lindat8, lindat4
        WRITE(13) lindat8, lindat4
      END DO
    END IF
    IF (nlines_in > 0) THEN
      REWIND 14
      DO i = 1, nlines_in
        READ(14)  lindat8, lindat4
        WRITE(13) lindat8, lindat4
      END DO
    END IF
  END IF
  CLOSE(UNIT=20, STATUS='DELETE')
  CLOSE(UNIT=14, STATUS='DELETE')

  ! ==========================================================================
  !  SECTION 4.  INITIALISE VOIGT LOOKUP TABLES
  ! ==========================================================================
  vsteps = 200.0
  CALL init_voigt_tables(vsteps, NTAB_VOIGT)

  ! ==========================================================================
  !  SPECTRV INITIALISATION  (formerly the opening sections of spectrv.f90)
  !
  !  These steps are done once, before the depth/opacity loop.
  !  READIN populates the mod_atlas_data atmosphere state, then
  !  run_xnfpelsyn computes continuum opacities and ion populations.
  ! ==========================================================================

  ! --- Read run parameters from spectrv.input ---
  BLOCK
    INTEGER :: ios25
    OPEN(UNIT=25, FILE=trim(DATADIR)//'spectrv.input', STATUS='OLD', &
         ACTION='READ', IOSTAT=ios25)
    IF (ios25 /= 0) THEN
      WRITE(6,'(A,A)') ' ERROR: cannot open ', trim(DATADIR)//'spectrv.input'
      STOP 1
    END IF
  END BLOCK
  READ(25, '(8F10.5)') rhoxj, r1, r101, ph1, pc1, psi1, PRDDOP, PRDPOW
  CLOSE(UNIT=25)
  slope  = 100.0D0 / (r101 - r1)
  origin = 1.5D0 - r1 * slope
  ! Store in module so process_wavelength_point / process_linecen_record can use them
  rhoxj_sv  = rhoxj;  r1_sv  = r1;  r101_sv = r101
  ph1_sv    = ph1;    pc1_sv = pc1; psi1_sv = psi1
  slope_sv  = slope;  origin_sv = origin

  ! --- Run READIN then run_xnfpelsyn ---
  ! IFPRES is set by READIN from the model cards (e.g. "PRESSURE OFF").
  itemp_a  = 1
  OPEN(UNIT=5,  FILE='fort.5',  STATUS='OLD', ACTION='READ')
  OPEN(UNIT=17, FILE=trim(DATADIR)//'continua.dat', STATUS='OLD', ACTION='READ')
  CALL readin(20)
  ! Keep unit 5 open: MOLEC reads molecular data from INPUTDATA (=5)
  ! on its first call during run_xnfpelsyn.
  CALL run_xnfpelsyn()
  CLOSE(UNIT=5)
  CLOSE(UNIT=17)

  ! ==========================================================================
  !  SECTIONS 5+6.  POPULATE SYNTHE WORKING ARRAYS FROM MODULE ARRAYS
  !
  !  READIN + run_xnfpelsyn have now populated the atlas module state and
  !  the synthe_module xf_* arrays.  Copy into the local working arrays.
  ! ==========================================================================

  ! --- Section 5: edge grid and header ---
  nrhox  = nrhox_a
  teff   = teff_a
  glog   = glog_a
  DO i = 1, 74
    title(i) = DBLE(ICHAR(title_atlas(i)))
  END DO
  WRITE(6,'(A,F10.1,3X,A,F6.3,3X,74A1)') ' TEFF=', teff, 'GRAV=', glog, &
        (CHAR(INT(title(i))), i=1,74)

  nedge  = nedge_m
  DO iedge = 1, nedge_m
    frqedg(iedge) = frqedg_m(iedge)
    wledge(iedge) = ABS(wledge_m(iedge))
    cmedge(iedge) = cmedge_m(iedge)
  END DO
  DO iedge = 1, mm
    idmol(iedge)  = idmol_m(iedge)
    momass(iedge) = momass_m(iedge)
  END DO
  ncon_tape = numnu_m
  DO nu = 1, numnu_m
    confrq(nu) = freqset_m(nu)
  END DO

  DO iedge = 2, nedge
    halfedge(iedge-1) = (wledge(iedge-1) + wledge(iedge)) * 0.5D0
    deledge(iedge-1)  = (wledge(iedge) - wledge(iedge-1))**2 * 0.5D0
  END DO

  ! --- Section 6: depth structure ---
  itemp = 1
  DO j = 1, nrhox
    t(j)       = REAL(xf_t(j))
    tkev(j)    = REAL(xf_tkev(j))
    tk(j)      = REAL(xf_tk(j))
    hkt(j)     = xf_hkt(j)
    tlog(j)    = REAL(xf_tlog(j))
    hckt(j)    = xf_hckt(j)
    p(j)       = REAL(xf_p(j))
    xne(j)     = REAL(xf_xne(j))
    xnatom(j)  = REAL(xf_xnatom(j))
    rho(j)     = REAL(xf_rho(j))
    rhox(j)    = REAL(xf_rhox(j))
    vturb(j)   = REAL(xf_vturb(j))
    xnfh(j)    = REAL(xf_xnfh(j))
    xnfhe(j,1) = REAL(xf_xnfhe(j,1))
    xnfhe(j,2) = REAL(xf_xnfhe(j,2))
    xnfh2(j)   = REAL(xf_xnfh2(j))
  END DO

  DO j = 1, kw_p
    velshift_arr(j) = deckj(1,j)
    hfield(j)       = deckj(2,j)
  END DO

  ! --- Unpack continuum opacity tables ---
  DO j = 1, nrhox_a
    nu = 0
    DO iedge = 1, nedge_m - 1
      nu = nu + 1;  contabs_sv(1,iedge,j)  = contabs_m(nu,j)
                    contscat_sv(1,iedge,j) = contscat_m(nu,j)
      nu = nu + 1;  contabs_sv(2,iedge,j)  = contabs_m(nu,j)
                    contscat_sv(2,iedge,j) = contscat_m(nu,j)
      nu = nu + 1;  contabs_sv(3,iedge,j)  = contabs_m(nu,j)
                    contscat_sv(3,iedge,j) = contscat_m(nu,j)
    END DO
  END DO

  ! --- SPECTRV Section 4: depth-dependent work arrays ---
  DO j = 1, nrhox
    bone(j)   = 1.0D0
    bfudge_sv(j) = bhyd_gs(j)**ph1 * (bc1_gs(j)/bc2_gs(j))**pc1 * &
                   (bsi1_gs(j)/bsi2_gs(j))**psi1
    fscat_sv(j)  = 0.0D0
    IF (rhoxj /= 0.0D0 .AND. rhox_a(j)/rhoxj < 100.0D0) &
      fscat_sv(j) = EXP(-rhox_a(j) / rhoxj)
  END DO

  ! --- SPECTRV Section 5: populate /TRASH/ and wavelength grid (replaces unit-9 header read) ---
  ifvac_sv      = ifvac
  IFPRES      = 0
  itemp_a       = 1
  wlbeg_loc       = wlbeg
  wlend_loc       = wlend
  vt_loc          = DBLE(turbv)
  air_loc         = DBLE(ifvac)    ! legacy: air_loc /= 0 means vacuum wavelengths
  ifvac_loc       = ifvac
  n10_loc         = 0              ! NLTE line centres (not currently used)

  ! Store 'A' or 'V' flag in title(74)
  BLOCK
    INTEGER(8) :: itmp
    itmp = INT(ICHAR('A'), 8)
    IF (ifvac == 1) itmp = INT(ICHAR('V'), 8)
    title(74) = TRANSFER(itmp, title(74))
  END BLOCK

  WRITE(6, '(A,F10.4,A,F10.1,A,F10.4,A,I6,A,I2)') &
    ' WLBEG=', wlbeg, '   RESOLU=', resolu, '   WLEND=', wlend, &
    '   LENGTH=', length, '   NRHOX=', nrhox

  ! NOTE: turbv is NOT added to vturb_a here.  In the original code,
  ! XNFPELSYN computed Doppler widths using only the model's vturb.
  ! The extra turbv from fort.93 is added later (line ~526) when building
  ! qdopple in the per-depth SYNTHE loop, matching the original flow.

  ! Build log-wavelength grid for SPECTRV (/FRESET/ COMMON)
  ratio     = 1.0D0 + 1.0D0 / resolu
  ratiolg   = LOG(ratio)
  ixwlbeg   = INT(LOG(wlbeg) / ratiolg)
  wbegin    = EXP(DBLE(ixwlbeg) * ratiolg)
  IF (wbegin < wlbeg) THEN
    ixwlbeg = ixwlbeg + 1
    wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  END IF
  wbegin_loc  = wbegin
  DELTAW  = resolu
  NULO    = 1
  NUHI    = length
  NUMNU   = length
  wlbeg_loc   = wbegin
  wlend_loc   = wbegin * ratio**(NUHI - 1)
  nmu2_loc    = NMU * 2

  ! Initialise SPECTRV ASCII plot array and edge bracket
  aplot   = ' '
  iedge_sv = 1

  ! ==========================================================================
  !  SECTION 7.  MAIN DEPTH LOOP  (SYNTHE opacity accumulation)
  ! ==========================================================================
  ALLOCATE(opacity_matrix(length, nrhox))
  opacity_matrix = 0.0

  ilines   = 0
  n12      = nlines_in
  nlines   = nlines_in + n19
  WRITE(6,'(A,I10,A,I10,A,I10)') ' TOTAL LINES: n12=', n12, '  n19=', n19, '  nlines=', nlines

  depth_loop: DO j = 1, nrhox

    REWIND 12
    buffer(1:length) = 0.0

    ! Read continuum opacity (total = abs+scat) from module array instead of fort.10
    ! continall_m(nu,j) holds the same values as the former CONTINALL(nu,j) record.
    DO nu = 1, numnu_m
      qablog(nu) = continall_m(nu, j)
    END DO

    nu = 0
    DO iedge = 1, nedge-1
      nu = nu + 1;  ablog(1,iedge) = REAL(qablog(nu))
      nu = nu + 1;  ablog(2,iedge) = REAL(qablog(nu))
      nu = nu + 1;  ablog(3,iedge) = REAL(qablog(nu))
    END DO

    iedge = 1
    DO nbuff_i = 1, length
      wave8 = wbegin * ratio**(nbuff_i-1)
      DO WHILE (wave8 >= wledge(iedge+1) .AND. iedge < nedge-1)
        iedge = iedge + 1
      END DO
      continuum(nbuff_i) = REAL( &
        ((wave8 - halfedge(iedge))*(wave8 - wledge(iedge+1))*ablog(1,iedge) + &
         (wledge(iedge) - wave8)*(wave8 - wledge(iedge+1))*2.0D0*ablog(2,iedge) + &
         (wave8 - wledge(iedge))*(wave8 - halfedge(iedge))*ablog(3,iedge)) / &
        deledge(iedge) )
    END DO
    DO nbuff_i = 1, length
      continuum(nbuff_i) = 10.0**continuum(nbuff_i)
    END DO

    ! Read ion populations and Doppler widths from module array instead of fort.10.
    ! xnfpel_m / dopple_m are (6, mw, kw) -- unpack into flat mw6 arrays.
    DO nelem_i = 1, mw
      DO i = 1, 6
        qxnfpel((nelem_i-1)*6 + i) = xnfpel_m(i, nelem_i, j)
        qdopple((nelem_i-1)*6 + i) = dopple_m (i, nelem_i, j)
      END DO
    END DO

    xnfph(j,1)   = REAL(qxnfpel(1))
    xnfph(j,2)   = REAL(qxnfpel(2))
    xnfphe(j,1)  = REAL(qxnfpel(7))
    xnfphe(j,2)  = REAL(qxnfpel(8))
    xnfphe(j,3)  = REAL(qxnfpel(9))

    DO i = 1, mw6_p
      xnfpel(i) = 0.0D0
      IF (qxnfpel(i) < 1.0D25) xnfpel(i) = qxnfpel(i)
    END DO

    DO i = 1, mw6_p
      qdopple(i)  = SQRT(qdopple(i)**2 + (DBLE(turbv)/CLIGHT_KMS)**2)
      dopple(i)   = qdopple(i)
      xnfpel(i)   = xnfpel(i) / rho(j)
      xnfdop(i)   = qxnfpel(i) / xf_rho(j) / qdopple(i)
    END DO

    txnxn(j) = (xnfh(j) + 0.42*xnfhe(j,1) + 0.85*xnfh2(j)) * &
               (t(j)/10000.0)**0.3

    nvshift = INT(resolu * DBLE(velshift_arr(j)) / CLIGHT_KMS + 0.5D0)

    mlines = 0
    IF (n19 > 0) CALL compute_line_opacity(j, n19, cutoff, velshift_arr(j), ifvac, linout)

    IF (n12 > 0) THEN
    n191    = n19 + 1
    alpha_s = 0.0

    DO iline = n191, nlines
      READ(12) nbuff_s, congf_s, congf_nel, elo_s, gamrf, gamsf, gamwf

      kappa0_s = congf_s * REAL(xnfdop(congf_nel))
      kapmin_s = continuum(MIN(MAX(nbuff_s,1),length)) * cutoff

      IF (kappa0_s < kapmin_s) CYCLE

      kappa0_s = kappa0_s * REAL(EXP(-elo_s * hckt(j)))
      IF (kappa0_s < kapmin_s) CYCLE

      IF (alpha_s /= 0.0) THEN
        nelem_i = INT(congf_nel/6) + 1
        v2_s    = (1.0 - alpha_s) / 2.0
        hfac(j)  = (t(j)/10000.0)**v2_s
        hefac(j) = 0.628*(2.0991E-4*t(j)*(0.25 + 1.008/atmass(nelem_i)))**v2_s
        h2fac(j) = 1.08 *(2.0991E-4*t(j)*(0.50 + 1.008/atmass(nelem_i)))**v2_s
        txnxn(j) = xnfh(j)*hfac(j) + xnfhe(j,1)*hefac(j) + xnfh2(j)*h2fac(j)
      END IF

      adamp_s = REAL((gamrf + gamsf*xne(j) + gamwf*txnxn(j)) / dopple(congf_nel))
      nbuff_s = nbuff_s + nvshift

      n10dop     = INT(10.0D0 * dopple(congf_nel) * DBLE(resolu))
      dopple_nel = REAL(dopple(congf_nel))

      centre_on_grid: IF (nbuff_s >= 1 .AND. nbuff_s <= length) THEN
        mlines = mlines + 1
        IF (adamp_s < VOIGT_APPROX_LIMIT) THEN
          kapcen_s = kappa0_s * (1.0 - 1.128*adamp_s)
        ELSE
          kapcen_s = kappa0_s * voigt_profile(0.0, adamp_s)
        END IF
        IF (linout >= 0) WRITE(15) iline, kapcen_s
        buffer(nbuff_s) = buffer(nbuff_s) + kapcen_s
      END IF centre_on_grid

      IF (adamp_s < VOIGT_APPROX_LIMIT) THEN
        tabstep = vsteps / (dopple_nel * REAL(resolu))
        tabi    = 1.5
        DO nstep = 1, n10dop
          tabi = tabi + tabstep
          profile(nstep) = kappa0_s * (h0tab(INT(tabi)) + adamp_s*h1tab(INT(tabi)))
          IF (profile(nstep) < kapmin_s) EXIT
        END DO
      ELSE
        dvoigt = 1.0 / dopple_nel / REAL(resolu)
        DO nstep = 1, n10dop
          profile(nstep) = kappa0_s * voigt_profile(REAL(nstep)*dvoigt, adamp_s)
          IF (profile(nstep) < kapmin_s) EXIT
        END DO
      END IF

      IF (nstep > n10dop) THEN
        x_wing  = profile(n10dop) * REAL(n10dop)**2
        maxstep = INT(SQRT(x_wing/kapmin_s)) + 1
        maxstep = MIN(maxstep, maxprof)
        n1 = n10dop + 1
        DO nstep = n1, maxstep
          profile(nstep) = x_wing / REAL(nstep)**2
        END DO
        nstep = maxstep
      END IF

      IF (nbuff_s+nstep < 1 .OR. nbuff_s-nstep > length) CYCLE

      IF (nbuff_s < length) THEN
        maxred    = MIN(length - nbuff_s, nstep)
        minblue_i = MAX(1, 1 - nbuff_s)
        DO istep = minblue_i, maxred
          buffer(nbuff_s + istep) = buffer(nbuff_s + istep) + profile(istep)
        END DO
        IF (nbuff_s <= 1) CYCLE
      END IF

      maxblue   = MIN(nbuff_s - 1, nstep)
      minblue_i = MAX(1, nbuff_s - length)
      DO istep = minblue_i, maxblue
        buffer(nbuff_s - istep) = buffer(nbuff_s - istep) + profile(istep)
      END DO

    END DO   ! iline loop

    END IF   ! n12 > 0

    ! Store opacity buffer into in-memory matrix (replaces unit 14 write)
    opacity_matrix(1:length, j) = buffer(1:length)

    WRITE(6,'(30X,2I10,A)') j, mlines, ' LINES USED'
    mlinej(j) = mlines
    ilines    = ilines + mlines

  END DO depth_loop

  ! ==========================================================================
  !  SECTION 8.  SPECTRV WAVELENGTH LOOP
  !
  !  For each wavelength point, assemble the opacity vector from the
  !  in-memory opacity_matrix, apply stimulated-emission correction,
  !  and call process_wavelength_point for the radiative transfer.
  ! ==========================================================================
  n9 = 0

  ! Open ASCII spectrum output (replaces standalone syntoascanga post-processing)
  OPEN(UNIT=11, FILE='fort.11', STATUS='REPLACE', ACTION='WRITE')
  BLOCK
    REAL(8) :: wend_loc, vstep_loc
    wend_loc  = wbegin * ratio**(length-1)
    vstep_loc = CLIGHT_KMS / resolu
  END BLOCK

  DO n9 = 1, length
    wave8 = wbegin * ratio**(n9-1)
    freq8 = CLIGHT_NM_HZ / wave8
    DO j = 1, nrhox
      asynth(j) = opacity_matrix(n9, j) * (1.0 - EXP(-freq8*hkt(j)))
    END DO

    ! --- SPECTRV wavelength loop body ---
    CALL process_wavelength_point(n9, wave8, asynth)

  END DO
  WRITE(6,'(I10,A)') length, ' OPACITY VECTORS PROCESSED'

  DEALLOCATE(opacity_matrix)

  ! ==========================================================================
  !  SECTION 9.  IDENTIFY USED LINES
  ! ==========================================================================
  n9 = 0
  identify_lines: IF (nlines > 0 .AND. linout >= 0 .AND. ilines > 0) THEN

    line_flag(1:nlines) = 0

    REWIND 15
    DO i = 1, ilines
      READ(15) iline
      line_flag(iline) = 1
    END DO

    ! First pass: count used lines
    REWIND 13
    DO i = 1, nlines
      READ(13) lindat8, lindat4
      IF (line_flag(i) == 0) CYCLE
      line_flag(i) = 0
      n9 = n9 + 1
      line_flag(i) = n9
    END DO

  END IF identify_lines

  nlines_sv = n9
  CLOSE(UNIT=12, STATUS='DELETE')

  IF (ilines == 0 .OR. linout < 0) THEN
    CLOSE(UNIT=13, STATUS='DELETE')
    CLOSE(UNIT=11)
    STOP
  END IF

  ! Second pass: store compacted lindat records in memory
  BLOCK
    REAL(8), ALLOCATABLE :: lindat8_arr(:,:)
    REAL(4), ALLOCATABLE :: lindat4_arr(:,:)
    INTEGER :: k9

    ALLOCATE(lindat8_arr(14, n9), lindat4_arr(28, n9))
    REWIND 13
    k9 = 0
    DO i = 1, nlines
      READ(13) lindat8, lindat4
      IF (line_flag(i) == 0) CYCLE
      k9 = k9 + 1
      lindat8_arr(:, k9) = lindat8
      lindat4_arr(:, k9) = lindat4
    END DO
    CLOSE(UNIT=13, STATUS='DELETE')

  ! ==========================================================================
  !  SECTION 10.  LINE-CENTRE OPACITIES AND IDENTIFICATION RECORDS
  !
  !  Build in-memory line-centre opacity matrix from unit 15 (iline, kapcen)
  !  pairs, then loop over lines calling process_linecen_record for each.
  ! ==========================================================================
  REWIND 15

  ALLOCATE(linecen_matrix(n9, nrhox))
  linecen_matrix = 0.0

  DO j = 1, nrhox
    maxline = mlinej(j)
    IF (maxline == 0) CYCLE
    DO l = 1, maxline
      READ(15) iline, kapcen_s
      i9 = line_flag(iline)
      IF (i9 == 0) CYCLE
      linecen_matrix(i9, j) = kapcen_s
    END DO
  END DO

  IF (n9 > maxlin_p) THEN
    WRITE(6,'(A)') ' TOO MANY LINES TO TRANSPOSE'
    STOP
  END IF

  ncen    = 0

  n910 = n9 + n10_loc

  DO i = 1, n9
    ncen = ncen + 1
    lindat8 = lindat8_arr(:, i)
    lindat4 = lindat4_arr(:, i)
    wl        = lindat8(1);   e         = lindat8(2);   ep        = lindat8(3)
    label(1)  = lindat8(4);   label(2)  = lindat8(5)
    labelp(1) = lindat8(6);   labelp(2) = lindat8(7)
    other1(1) = lindat8(8);   other1(2) = lindat8(9)
    other2(1) = lindat8(10);  other2(2) = lindat8(11)
    wlvac  = lindat8(12)
    center = lindat8(13)
    concen = lindat8(14)
    freq8  = CLIGHT_NM_HZ / wlvac
    DO j = 1, nrhox
      alinec(j) = linecen_matrix(i, j) * (1.0 - EXP(-freq8*hkt(j)))
    END DO

    CALL process_linecen_record(lindat8, lindat4, alinec)
  END DO

  DEALLOCATE(linecen_matrix)
  DEALLOCATE(lindat8_arr, lindat4_arr)
  END BLOCK

  CLOSE(UNIT=15, STATUS='DELETE')
  WRITE(6,'(I10,A)') ncen, ' LINE CENTER RECORDS PROCESSED'

  ! --- SPECTRV Section 8: NLTE line centres (unit 20) ---
  IF (n10_loc > 0) THEN
    REWIND 20
    wavold   = 0.0D0
    iedge_sv = 1

    nlte_loop: DO iline = 1, n10_loc
      READ(20) lindat8_c, lindat4_c
      CALL unpack_lindat()
      wave_sv = wlvac_c
      IF (wave_sv < wavold) iedge_sv = 1
      wavold = wave_sv

      CALL setup_opacity_sv(wave_sv)

      CALL josh(1, IFSURF)
      IF (IFSURF == 1) concen_c = HNU(1)
      IF (IFSURF == 2) concen_c = SURFI(1)

      CALL josh(1, IFSURF)
      IF (IFSURF == 1) center_c = HNU(1)
      IF (IFSURF == 2) center_c = SURFI(1)
      concen_c = surf_sv(1)
      resid    = center_c / concen_c

      WRITE(6, '(1H0,F10.4,F7.3,F5.1,F12.3,F5.1,F12.3,F9.2,1X,A8,A2,A8,A2,' // &
                'F12.4,F9.3,1P2E11.3,' // &
                '/1X,0PF10.4,I4,F6.2,F6.2,F6.2,A4,I2,I2,I3,F7.3,I3,F7.3,1X,' // &
                'A8,A2,A8,A2,F7.4,F7.3,3F6.2)') &
        wl_c, gflog_c, xj_c, e_c, xjp_c, ep_c, code_c, label_c, labelp_c, &
        wlvac_c, resid, center_c, concen_c, &
        wl_c, nelion_c, gr_c, gs_c, gw_c, ref_c, nblo_c, nbup_c, &
        iso1_c, x1_c, iso2_c, x2_c, &
        other1_c, other2_c, dwl_c, dgflog_c, dgammar_c, dgammas_c, dgammaw_c

    END DO nlte_loop
    CLOSE(UNIT=20, STATUS='DELETE')
  END IF

  CLOSE(UNIT=11)

  STOP

CONTAINS

  ! ==========================================================================
  !  SETUP_OPACITY_SV(wave)
  !
  !  SPECTRV continuum and frequency setup for one wavelength point.
  !  Advances iedge_sv to bracket WAVE, computes quadratic Lagrange weights,
  !  fills ACONT / SIGMAC from contabs_sv / contscat_sv, and sets all
  !  frequency-dependent ATLAS COMMON variables (FREQ, EHVKT, STIM,
  !  BNU, etc.).
  !
  !  Replaces setup_opacity_point() in standalone spectrv.f90.
  !  Uses iedge_sv (module-level bracket, separate from SYNTHE's iedge).
  ! ==========================================================================
  SUBROUTINE setup_opacity_sv(wave)
    REAL(8), INTENT(IN) :: wave

    REAL(8) :: freq15, c1, c2, c3
    INTEGER :: jj

    DO WHILE (wave >= wledge(iedge_sv + 1))
      iedge_sv = iedge_sv + 1
    END DO

    c1 = (wave - halfedge(iedge_sv)) * (wave - wledge(iedge_sv+1)) / deledge(iedge_sv)
    c2 = (wledge(iedge_sv) - wave)   * (wave - wledge(iedge_sv+1)) * 2.0D0 / deledge(iedge_sv)
    c3 = (wave - wledge(iedge_sv))   * (wave - halfedge(iedge_sv)) / deledge(iedge_sv)

    DO jj = 1, nrhox
      ACONT(jj)  = 10.0D0**(c1*contabs_sv(1,iedge_sv,jj)  + &
                               c2*contabs_sv(2,iedge_sv,jj)  + &
                               c3*contabs_sv(3,iedge_sv,jj))
      SIGMAC(jj) = 10.0D0**(c1*contscat_sv(1,iedge_sv,jj) + &
                               c2*contscat_sv(2,iedge_sv,jj) + &
                               c3*contscat_sv(3,iedge_sv,jj))
    END DO

    FREQ   = CLIGHT_NM_HZ / wave
    freq15   = FREQ / 1.0D15
    FREQLG = LOG(FREQ)
    DO jj = 1, nrhox
      EHVKT(jj)  = EXP(-FREQ * hkt_a(jj))
      STIM(jj)   = 1.0D0 - EHVKT(jj)
      BNU(jj)    = PLANCK_PREFACTOR * freq15**3 * EHVKT(jj) / STIM(jj)
      SIGPRD(jj) = 0.0D0
      ALINE(jj)  = 0.0D0
      SIGMAL(jj) = 0.0D0
      SLINE(jj)  = BNU(jj)
      SCONT(jj)  = BNU(jj)
    END DO

  END SUBROUTINE setup_opacity_sv


  ! ==========================================================================
  !  PROCESS_WAVELENGTH_POINT(nu_idx, wave_in, asyn)
  !
  !  Executes the SPECTRV wavelength loop body for one wavelength point.
  !  Replaces the READ(9) asynth + loop body that was in the standalone
  !  SPECTRV wavelength loop.
  !
  !  Arguments:
  !    nu_idx  : wavelength point index (1..LENGTH)
  !    wave_in : wavelength in Å (from SYNTHE grid; converted to nm internally)
  !    asyn    : stimulated-emission-corrected line opacity vector (nrhox)
  ! ==========================================================================
  SUBROUTINE process_wavelength_point(nu_idx, wave_in, asyn)
    INTEGER,  INTENT(IN) :: nu_idx
    REAL(8),  INTENT(IN) :: wave_in
    REAL(4),  INTENT(IN) :: asyn(kw)

    REAL(8) :: wave_nm, q_loc(41), slinec_loc
    INTEGER :: jj, mu_loc, iresid_loc, iplot_loc, ii

    ! wave_in is in nm (from the SYNTHE wavelength grid)
    wave_nm = wave_in

    CALL setup_opacity_sv(wave_nm)

    ! First JOSH call: pure continuum
    CALL josh(IFSCAT, IFSURF)
    DO mu_loc = 1, NMU
      IF (IFSURF == 1) surf_sv(mu_loc) = HNU(1)
      IF (IFSURF == 2) surf_sv(mu_loc) = SURFI(mu_loc)
    END DO

    IF (nu_idx == 1) WRITE(6, '(A)') ' '

    ! Add line opacity and compute blended source function
    DO jj = 1, nrhox
      ALINE(jj)  = DBLE(asyn(jj)) * (1.0D0 - fscat_sv(jj))
      SLINE(jj)  = BNU(jj) * STIM(jj) / (bfudge_sv(jj) - EHVKT(jj))
      SIGMAL(jj) = DBLE(asyn(jj)) * fscat_sv(jj)
    END DO

    ! Second JOSH call: continuum + line
    CALL josh(1, IFSURF)

    WRITE(33, '(F12.4,1P2E12.4,/(10E12.4))') wave_nm, HNU(1), surf_sv(1), &
         (TAUNU(jj), jj=1,nrhox)

    DO mu_loc = 1, NMU
      IF (IFSURF == 1) resid = HNU(1)    / surf_sv(mu_loc)
      IF (IFSURF == 2) resid = SURFI(mu_loc) / surf_sv(mu_loc)

      IF (nu_idx <= ABS(linout) .AND. mu_loc == 1) THEN
        iresid_loc   = NINT(resid * 1000.0D0)
        iplot_loc    = INT(resid * slope_sv + origin_sv)
        iplot_loc    = MAX(1, MIN(101, iplot_loc))
        aplot(iplot_loc) = 'X'
        WRITE(6, '(1X,I5,F11.4,I7,101A1)') nu_idx, wave_nm, iresid_loc, aplot
        aplot(iplot_loc) = ' '
      END IF

      q_loc(mu_loc)        = resid * surf_sv(mu_loc)
      q_loc(mu_loc + NMU) = surf_sv(mu_loc)
    END DO

    ! Write ASCII spectrum line: wavelength (Angstroms), flux, continuum flux
    WRITE(11, '(F11.4,2E15.6)') wave_nm * 10.0D0, q_loc(1), q_loc(NMU + 1)

  END SUBROUTINE process_wavelength_point


  ! ==========================================================================
  !  PROCESS_LINECEN_RECORD(ld8, ld4, alin)
  !
  !  Executes the SPECTRV line-centre loop body for one line identification
  !  record.  Replaces the READ(9) lindat8, lindat4, alinec + loop body in
  !  standalone SPECTRV section 7.
  !
  !  Arguments:
  !    ld8  : REAL*8 line parameter array (14 words)
  !    ld4  : REAL*4 line parameter array (28 words)
  !    alin : depth vector of line-centre opacities (nrhox)
  ! ==========================================================================
  SUBROUTINE process_linecen_record(ld8, ld4, alin)
    REAL(8), INTENT(IN) :: ld8(14)
    REAL(4), INTENT(IN) :: ld4(28)
    REAL(4), INTENT(IN) :: alin(kw)

    REAL(8) :: wave_nm, freq_loc, slinec_loc, concen_loc, center_loc
    INTEGER :: jj, mu_loc

    ! Copy into local arrays and unpack named scalars
    lindat8_c = ld8
    lindat4_c = ld4
    CALL unpack_lindat()
    ! alinec_c is in /LINDAT/ -- populate for line-centre processing
    alinec_c(1:nrhox) = alin(1:nrhox)

    wave_nm = wlvac_c

    IF (wave_nm < wavold) iedge_sv = 1
    wavold = wave_nm

    CALL setup_opacity_sv(wave_nm)

    ! Continuum-only JOSH
    CALL josh(1, IFSURF)
    IF (IFSURF == 1) concen_loc = HNU(1)
    IF (IFSURF == 2) concen_loc = SURFI(1)

    ! Add line + H line + extra line opacity
    freq_loc = CLIGHT_NM_HZ / wave_nm
    DO jj = 1, nrhox
      slinec_loc = BNU(jj) * STIM(jj) / (bfudge_sv(jj) - EHVKT(jj))
      ALINE(jj) = ahline_a(jj) + DBLE(alin(jj)) * (1.0D0 - fscat_sv(jj)) &
                    + AXLINE(jj)
      SLINE(jj) = BNU(jj)
      IF (ALINE(jj) > 0.0D0) &
        SLINE(jj) = (ahline_a(jj)*shline_a(jj) &
                      + DBLE(alin(jj))*slinec_loc*(1.0D0 - fscat_sv(jj)) &
                      + AXLINE(jj)*SXLINE(jj)) / ALINE(jj)
      SIGMAL(jj) = DBLE(alin(jj)) * fscat_sv(jj)
    END DO

    ! Line-centre JOSH
    CALL josh(1, IFSURF)
    IF (IFSURF == 1) center_loc = HNU(1)
    IF (IFSURF == 2) center_loc = SURFI(1)
    resid = center_loc / concen_loc

    ! Write line identification record to unit 16
    WRITE(16, '(F10.4,F7.3,F5.1,F12.3,F5.1,F12.3,F9.2,A8,A2,A8,A2,' // &
               '/0PF10.4,I4,F6.2,F6.2,F6.2,A4,I2,I2,I3,F7.3,I3,F7.3,' // &
               'A8,A2,A8,A2,F7.4,F7.3,3F6.2,F10.4)') &
      wl_c, gflog_c, xj_c, e_c, xjp_c, ep_c, code_c, label_c, labelp_c, &
      wl_c, nelion_c, gr_c, gs_c, gw_c, ref_c, nblo_c, nbup_c, &
      iso1_c, x1_c, iso2_c, x2_c, &
      other1_c, other2_c, dwl_c, dgflog_c, dgammar_c, dgammas_c, dgammaw_c, wl_c

  END SUBROUTINE process_linecen_record


  ! ==========================================================================
  !  UNPACK_LINDAT -- extract named scalars from lindat8_c / lindat4_c
  !
  !  Replaces the broken EQUIVALENCE overlay.  Must be called after every
  !  READ or assignment into lindat8_c / lindat4_c.
  ! ==========================================================================
  SUBROUTINE unpack_lindat()
    ! REAL(8) fields from lindat8_c
    wl_c        = lindat8_c(1)
    e_c         = lindat8_c(2)
    ep_c        = lindat8_c(3)
    label_c(1)  = lindat8_c(4)
    label_c(2)  = lindat8_c(5)
    labelp_c(1) = lindat8_c(6)
    labelp_c(2) = lindat8_c(7)
    other1_c(1) = lindat8_c(8)
    other1_c(2) = lindat8_c(9)
    other2_c(1) = lindat8_c(10)
    other2_c(2) = lindat8_c(11)
    wlvac_c     = lindat8_c(12)
    center_c    = lindat8_c(13)
    concen_c    = lindat8_c(14)
    ! INTEGER and REAL(4) fields from lindat4_c
    nelion_c  = TRANSFER(lindat4_c(1),  nelion_c)
    nblo_c    = TRANSFER(lindat4_c(2),  nblo_c)
    nbup_c    = TRANSFER(lindat4_c(3),  nbup_c)
    iso1_c    = TRANSFER(lindat4_c(4),  iso1_c)
    iso2_c    = TRANSFER(lindat4_c(5),  iso2_c)
    gammar_c  = lindat4_c(6)
    gammas_c  = lindat4_c(7)
    gammaw_c  = lindat4_c(8)
    ref_c     = lindat4_c(9)
    x1_c      = lindat4_c(10)
    x2_c      = lindat4_c(11)
    gflog_c   = lindat4_c(12)
    xj_c      = lindat4_c(13)
    xjp_c     = lindat4_c(14)
    code_c    = lindat4_c(15)
    elo_c     = lindat4_c(16)
    gf_c      = lindat4_c(17)
    gs_c      = lindat4_c(18)
    gr_c      = lindat4_c(19)
    gw_c      = lindat4_c(20)
    dwl_c     = lindat4_c(21)
    dgflog_c  = lindat4_c(22)
    dgammar_c = lindat4_c(23)
    dgammas_c = lindat4_c(24)
    dgammaw_c = lindat4_c(25)
    extra1_c  = lindat4_c(26)
    extra2_c  = lindat4_c(27)
    extra3_c  = lindat4_c(28)
  END SUBROUTINE unpack_lindat


END PROGRAM synthe_spectrv
