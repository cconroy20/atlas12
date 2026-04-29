! ==============================================================================
! mod_mklinelist.f90 — Kurucz line list preprocessor
!
! Public interface:
!   call run_mklinelist(wlbeg, wlend, resolu, teff, lines_list_path, datadir)
!
! After the call, lte_lines(:) and nlte_lines(:) are allocated and populated;
! nlines_lte and nlines_nlte hold their sizes.
!
! Line list sources are specified in a plain-text file (lines.list by default):
!   gfall   /path/to/gfall.dat
!   predict /path/to/predict.bin
!   mol     /path/to/molecule1.dat
!   mol     /path/to/molecule2.dat
!   h2o     /path/to/h2o.bin
! Blank lines and lines beginning with # are ignored.  Entry order does not
! matter; readers always run in canonical order: gfall -> predict -> mol(s)
! -> h2o.  The file must exist; a hard stop is issued if it is not found.
! ==============================================================================

MODULE mod_mklinelist
  IMPLICIT NONE
  PRIVATE

  ! --- Public types -----------------------------------------------------

  !  LTE line record — mirrors the fort.12 binary record layout:
  !    NBUFF(I4)  CGF(R4)  NELION(I4)  ELO(R4)
  !    GAMRF(R4)  GAMSF(R4)  GAMWF(R4)          [28 bytes]
  TYPE, PUBLIC :: lte_line_t
    INTEGER  :: nbuff    ! grid-index offset from ixwlbeg
    REAL(4)  :: cgf      ! line-strength factor (units: per (cm^-1 g^-1 cm^2))
    INTEGER  :: nelion   ! species index (nelem-1)*6 + ion_stage
    REAL(4)  :: elo      ! lower level energy (cm^-1)
    REAL(4)  :: gamrf    ! radiative  damping / 4πν
    REAL(4)  :: gamsf    ! Stark      damping / 4πν
    REAL(4)  :: gamwf    ! van der Waals damping / 4πν
  END TYPE lte_line_t

  !  NLTE / complex-profile line record — mirrors the fort.19 binary record:
  !    WLVAC(R8)  ELO(R4)  GF(R4)  NBLO(I4)  NBUP(I4)
  !    NELION(I4)  TYPE(I4)  NCON(I4)  NELIONX(I4)
  !    GAMMAR(R4)  GAMMAS(R4)  GAMMAW(R4)  NBUFF(I4)  LIM(I4)  [60 bytes]
  TYPE, PUBLIC :: nlte_line_t
    REAL(8)  :: wlvac    ! vacuum wavelength (nm)
    REAL(4)  :: elo      ! lower level energy (cm^-1)
    REAL(4)  :: gf       ! oscillator strength (CGF for most types; raw GF for
                         !   autoionizing TYPE=1 and continuum TYPE>3)
    INTEGER  :: nblo     ! lower NLTE level index
    INTEGER  :: nbup     ! upper NLTE level index
    INTEGER  :: nelion   ! species index
    INTEGER  :: itype    ! line type flag (see read_gfall dispatch table)
    INTEGER  :: ncon     ! number of continuum edges (isotope field reuse)
    INTEGER  :: nelionx  ! species index for departure coefficients
    REAL(4)  :: gammar   ! radiative  damping / 4πν
    REAL(4)  :: gammas   ! Stark      damping / 4πν
    REAL(4)  :: gammaw   ! van der Waals damping / 4πν
    INTEGER  :: nbuff    ! grid-index offset from ixwlbeg
    INTEGER  :: lim      ! wing extent index (0=widest, 6=narrowest)
  END TYPE nlte_line_t

  ! --- Public module arrays and counters --------------------------------

  TYPE(lte_line_t),  ALLOCATABLE, PUBLIC :: lte_lines(:)
  TYPE(nlte_line_t), ALLOCATABLE, PUBLIC :: nlte_lines(:)
  INTEGER,                        PUBLIC :: nlines_lte  = 0
  INTEGER,                        PUBLIC :: nlines_nlte = 0

  ! Set VERBOSE = 1 before calling run_mklinelist to get progress output on
  ! stdout.  Errors and warnings are always printed regardless of setting.
  INTEGER, PUBLIC :: VERBOSE = 0

  PUBLIC :: run_mklinelist

  ! --- Module-level ionisation potential table (replaces COMMON /potion/)
  REAL(8), SAVE :: potion(999)

  ! --- Teff gate for cool-atmosphere molecular lines (H2O, TiO) ---------
  ! Species skipped when Teff > this, even if listed in lines.list.
  REAL(8), PARAMETER :: TEFF_COOL_LIMIT = 5000.0D0

CONTAINS

  ! ============================================================================
  !  RUN_MKLINELIST — top-level entry point called from SYNTHE
  !
  !  Reads lines.list, dispatches to readers in canonical order, assembles
  !  the lte_lines and nlte_lines module arrays.
  ! ============================================================================
  SUBROUTINE run_mklinelist(wlbeg, wlend, resolu, teff, lines_list_path, datadir)
    REAL(8),          INTENT(IN) :: wlbeg, wlend, resolu, teff
    CHARACTER(LEN=*), INTENT(IN) :: lines_list_path
    CHARACTER(LEN=*), INTENT(IN) :: datadir

    ! synbeg_params — grid scalars passed to every reader
    REAL(8) :: ratio, ratiolg
    INTEGER :: ixwlbeg
    REAL(8) :: wbegin

    ! File paths parsed from lines.list
    CHARACTER(LEN=512) :: gfall_file, predict_file, h2o_file
    CHARACTER(LEN=512), SAVE :: mol_files(256)
    INTEGER :: nmol

    ! Temporary per-reader arrays (LTE)
    TYPE(lte_line_t),  ALLOCATABLE :: lte_gfall(:),   lte_predict(:)
    TYPE(lte_line_t),  ALLOCATABLE :: lte_mol(:),     lte_h2o(:)
    INTEGER :: n_lte_gfall, n_lte_predict, n_lte_mol, n_lte_h2o

    ! Temporary per-reader arrays (NLTE)
    TYPE(nlte_line_t), ALLOCATABLE :: nlte_gfall(:)
    INTEGER :: n_nlte_gfall

    INTEGER :: imol, ioff
    INTEGER :: n_lte_mol_before   ! per-file delta for mol entries

    ! Progress-report format strings
    CHARACTER(LEN=*), PARAMETER :: ROW_FMT  = '(a4,2x,i12,2x,i12,2x,a)'
    CHARACTER(LEN=*), PARAMETER :: SKIP_FMT = '(a4,2x,a12,2x,a12,2x,a,a)'

    CALL ionpots()

    ! --- Log-wavelength grid: wbegin is first grid point >= wlbeg ---------
    ratio   = 1.0D0 + 1.0D0 / resolu
    ratiolg = LOG(ratio)
    ixwlbeg = INT(LOG(wlbeg) / ratiolg)
    wbegin  = EXP(ixwlbeg * ratiolg)
    IF (wbegin .LT. wlbeg) THEN
      ixwlbeg = ixwlbeg + 1
      wbegin  = EXP(ixwlbeg * ratiolg)
    END IF

    ! --- Parse lines.list and init per-reader output arrays --------------
    gfall_file   = ''
    predict_file = ''
    h2o_file     = ''
    nmol         = 0
    CALL parse_lines_list(lines_list_path, datadir, &
                          gfall_file, predict_file, h2o_file, &
                          mol_files, nmol)

    n_lte_gfall   = 0
    n_lte_predict = 0
    n_lte_mol     = 0
    n_lte_h2o     = 0
    n_nlte_gfall  = 0
    ALLOCATE(lte_gfall(0))
    ALLOCATE(lte_predict(0))
    ALLOCATE(lte_mol(0))
    ALLOCATE(lte_h2o(0))
    ALLOCATE(nlte_gfall(0))

    IF (VERBOSE .EQ. 1) THEN
      WRITE(6,'(/a)') '=== mklinelist ==='
      WRITE(6,'(a4,2x,a12,2x,a12,2x,a)') 'src', 'n_lte', 'n_nlte', 'file'
    END IF

    ! --- Canonical reader order: gfall -> predict -> mol(s) -> h2o -------

    IF (gfall_file .NE. '') THEN
      CALL read_gfall(gfall_file, wlbeg, wlend, ratiolg, ixwlbeg, &
                      lte_gfall, n_lte_gfall, nlte_gfall, n_nlte_gfall)
      IF (VERBOSE .EQ. 1) &
        WRITE(6,ROW_FMT) 'gf', n_lte_gfall, n_nlte_gfall, TRIM(basename(gfall_file))
    END IF

    IF (predict_file .NE. '') THEN
      CALL read_predict(predict_file, wlbeg, wlend, ratiolg, ixwlbeg, &
                        lte_predict, n_lte_predict)
      IF (VERBOSE .EQ. 1) &
        WRITE(6,ROW_FMT) 'pr', n_lte_predict, 0, TRIM(basename(predict_file))
    END IF

    DO imol = 1, nmol
      IF (is_tio_file(mol_files(imol)) .AND. teff .GT. TEFF_COOL_LIMIT) THEN
        IF (VERBOSE .EQ. 1) &
          WRITE(6,SKIP_FMT) 'mol', 'skip', '', &
            'Teff > 5000 K: ', TRIM(basename(mol_files(imol)))
        CYCLE
      END IF
      n_lte_mol_before = n_lte_mol
      CALL read_molec(mol_files(imol), wlbeg, wlend, ratiolg, ixwlbeg, &
                      lte_mol, n_lte_mol)
      IF (VERBOSE .EQ. 1) &
        WRITE(6,ROW_FMT) 'mol', n_lte_mol - n_lte_mol_before, 0, &
          TRIM(basename(mol_files(imol)))
    END DO

    IF (h2o_file .NE. '') THEN
      IF (teff .GT. TEFF_COOL_LIMIT) THEN
        IF (VERBOSE .EQ. 1) &
          WRITE(6,SKIP_FMT) 'h2o', 'skip', '', &
            'Teff > 5000 K: ', TRIM(basename(h2o_file))
      ELSE
        CALL read_h2o(h2o_file, wlbeg, wlend, ratiolg, ixwlbeg, &
                      lte_h2o, n_lte_h2o)
        IF (VERBOSE .EQ. 1) &
          WRITE(6,ROW_FMT) 'h2o', n_lte_h2o, 0, TRIM(basename(h2o_file))
      END IF
    END IF

    ! --- Assemble module arrays from per-reader temporaries --------------
    nlines_nlte = n_nlte_gfall
    nlines_lte  = n_lte_gfall + n_lte_predict + n_lte_mol + n_lte_h2o

    ALLOCATE(nlte_lines(nlines_nlte))
    IF (nlines_nlte .GT. 0) &
      nlte_lines(1:nlines_nlte) = nlte_gfall(1:nlines_nlte)

    ALLOCATE(lte_lines(nlines_lte))
    ioff = 0
    IF (n_lte_gfall .GT. 0) THEN
      lte_lines(ioff+1 : ioff+n_lte_gfall) = lte_gfall(1:n_lte_gfall)
      ioff = ioff + n_lte_gfall
    END IF
    IF (n_lte_predict .GT. 0) THEN
      lte_lines(ioff+1 : ioff+n_lte_predict) = lte_predict(1:n_lte_predict)
      ioff = ioff + n_lte_predict
    END IF
    IF (n_lte_mol .GT. 0) THEN
      lte_lines(ioff+1 : ioff+n_lte_mol) = lte_mol(1:n_lte_mol)
      ioff = ioff + n_lte_mol
    END IF
    IF (n_lte_h2o .GT. 0) THEN
      lte_lines(ioff+1 : ioff+n_lte_h2o) = lte_h2o(1:n_lte_h2o)
    END IF

    DEALLOCATE(lte_gfall, lte_predict, lte_mol, lte_h2o, nlte_gfall)

    IF (VERBOSE .EQ. 1) THEN
      WRITE(6,'(a4,2x,a12,2x,a12)') '---', '------------', '------------'
      WRITE(6,'(a4,2x,i12,2x,i12)') 'tot', nlines_lte, nlines_nlte
      WRITE(6,'(a)') ''
    END IF

  END SUBROUTINE run_mklinelist


  ! ============================================================================
  !  PARSE_LINES_LIST — read lines.list and populate file path variables
  ! ============================================================================
  SUBROUTINE parse_lines_list(listfile, datadir, gfall_file, predict_file, h2o_file, &
                               mol_files, nmol)
    CHARACTER(LEN=*),   INTENT(IN)  :: listfile
    CHARACTER(LEN=*),   INTENT(IN)  :: datadir
    CHARACTER(LEN=512), INTENT(OUT) :: gfall_file, predict_file, h2o_file
    CHARACTER(LEN=512), INTENT(OUT) :: mol_files(256)
    INTEGER,            INTENT(OUT) :: nmol

    INTEGER, PARAMETER :: LU = 50
    CHARACTER(LEN=1024) :: line
    CHARACTER(LEN=32)   :: keyword
    CHARACTER(LEN=512)  :: filepath
    INTEGER :: ios

    gfall_file   = ''
    predict_file = ''
    h2o_file     = ''
    nmol         = 0

    OPEN(UNIT=LU, FILE=listfile, STATUS='old', ACTION='read', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') 'ERROR: lines.list not found: ', TRIM(listfile)
      CALL EXIT(1)
    END IF

    DO
      READ(LU, '(a)', IOSTAT=ios) line
      IF (ios .NE. 0) EXIT
      line = ADJUSTL(line)
      IF (line .EQ. '')            CYCLE
      IF (line(1:1) .EQ. '#')      CYCLE

      ! Read keyword only via list-directed (stops at first whitespace).
      ! List-directed read of the full line would truncate filepath at '/'
      ! since Fortran treats '/' as a record terminator in list-directed I/O.
      READ(line, *, IOSTAT=ios) keyword
      IF (ios .NE. 0) THEN
        WRITE(6,'(a,a)') 'WARNING: ignoring malformed line in lines.list: ', TRIM(line)
        CYCLE
      END IF
      ! Extract filepath as everything after the keyword, stripping any
      ! leading whitespace (spaces or tabs) between keyword and filename.
      ! adjustl only strips spaces, so we scan manually.
      BLOCK
        INTEGER :: i, klen
        klen = LEN_TRIM(keyword)
        i    = klen + 1
        DO WHILE (i .LE. LEN_TRIM(line))
          IF (line(i:i) .EQ. ' ' .OR. line(i:i) .EQ. CHAR(9)) THEN
            i = i + 1
          ELSE
            EXIT
          END IF
        END DO
        filepath = TRIM(datadir) // line(i:LEN_TRIM(line))
      END BLOCK

      SELECT CASE (TRIM(keyword))
      CASE ('gfall')
        gfall_file = TRIM(filepath)
      CASE ('predict')
        predict_file = TRIM(filepath)
      CASE ('mol')
        IF (nmol .LT. 256) THEN
          nmol = nmol + 1
          mol_files(nmol) = TRIM(filepath)
        ELSE
          WRITE(6,'(a)') 'WARNING: more than 256 mol entries in lines.list; extras ignored'
        END IF
      CASE ('h2o')
        h2o_file = TRIM(filepath)
      CASE DEFAULT
        WRITE(6,'(a,a)') 'WARNING: unrecognised keyword in lines.list: ', TRIM(keyword)
      END SELECT
    END DO

    CLOSE(LU)

  END SUBROUTINE parse_lines_list


  ! ============================================================================
  !  READ_GFALL — read Kurucz gfall ASCII atomic/ionic line list
  !               (translated from rgfall.for)
  !
  !  LTE plain Voigt lines go into lte_arr.
  !  Complex-profile lines (H, He, autoionizing, PRD, continuum edges, NLTE
  !  departure coefficient lines) go into nlte_arr.
  ! ============================================================================
  SUBROUTINE read_gfall(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                         lte_arr, n_lte, nlte_arr, n_nlte)
    CHARACTER(LEN=*),              INTENT(IN)    :: filename
    REAL(8),                       INTENT(IN)    :: wlbeg, wlend, ratiolg
    INTEGER,                       INTENT(IN)    :: ixwlbeg
    TYPE(lte_line_t),  ALLOCATABLE, INTENT(INOUT) :: lte_arr(:)
    INTEGER,                       INTENT(INOUT) :: n_lte
    TYPE(nlte_line_t), ALLOCATABLE, INTENT(INOUT) :: nlte_arr(:)
    INTEGER,                       INTENT(INOUT) :: n_nlte

    ! --- LINDAT fields ---
    REAL(8)  :: wl, e, ep, label(2), labelp(2), other1(2), other2(2)
    REAL(8)  :: wlvac
    REAL(4)  :: nelion_r
    REAL(4)  :: ref, gflog, xj, xjp, code, gf, gs, gr, gw
    REAL(4)  :: dwl, dgflog, dgammar, dgammas, dgammaw, dwliso
    INTEGER  :: nblo, nbup, iso1, iso2, isoshift
    REAL(4)  :: x1, x2

    REAL(8)  :: lindat8(14)
    REAL(4)  :: lindat4(28)
    EQUIVALENCE (lindat8(1), wl),      (lindat4(1), nelion_r)

    CHARACTER(LEN=10) :: cother1, cother2
    EQUIVALENCE (cother1, other1(1)), (cother2, other2(1))
    CHARACTER(LEN=3)  :: auto_flag
    CHARACTER(LEN=6)  :: ixfixfp

    ! CODEX: species eligible for NLTE departure coefficients (values are
    ! Kurucz species codes × 100, so 1.00 -> 100, 2.01 -> 201, etc.)
    INTEGER, PARAMETER :: CODEX(17) = [ &
      100, 200, 201, 600, 601, 1200, 1201, 1300, 1301, &
      1400, 1401, 2000, 2001, 800, 1100, 500, 1900]

    REAL(8), PARAMETER :: DELLIM(7) = &
      [100.D0, 30.D0, 10.D0, 3.D0, 1.D0, 0.3D0, 0.1D0]

    ! --- working variables ---
    INTEGER  :: ios, nelem, icharge, linesize, lim, ncon, nelionx, itype, ic
    INTEGER  :: ixwl, nbuff, ishift, ishiftp, icode_r
    REAL(8)  :: delfactor, eshift, eshiftp, frelin, cgf, frq4pi
    REAL(8)  :: effnsq, zeff, rsqup, rsqlo, eup
    REAL(8)  :: gammar_d, gammas_d, gammaw_d, elo_d
    INTEGER  :: nelion_i

    ! --- dynamic growth buffers ---
    INTEGER, PARAMETER :: CHUNK = 100000
    TYPE(lte_line_t),  ALLOCATABLE :: lte_buf(:)
    TYPE(nlte_line_t), ALLOCATABLE :: nlte_buf(:)
    INTEGER :: lte_cap, nlte_cap

    lte_cap  = CHUNK
    nlte_cap = CHUNK
    ALLOCATE(lte_buf(lte_cap))
    ALLOCATE(nlte_buf(nlte_cap))
    n_lte  = 0
    n_nlte = 0

    OPEN(UNIT=11, FILE=filename, STATUS='old', ACTION='read', &
         FORM='formatted', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') 'ERROR: cannot open gfall file: ', TRIM(filename)
      CALL EXIT(1)
    END IF

    delfactor = 1.0D0
    IF (wlbeg .GT. 500.0D0) delfactor = wlbeg / 500.0D0

    dwl       = 0.0
    dgflog    = 0.0
    dgammar   = 0.0
    dgammas   = 0.0
    dgammaw   = 0.0
    dwliso    = 0.0
    other1(2) = 0.0D0
    other2(1) = 0.0D0
    other2(2) = 0.0D0

    DO
      READ(11, &
        '(F11.4,F7.3,F6.2,F12.3,F5.1,1X,A8,A2,F12.3,F5.1,1X,A8,A2,' // &
        'F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,A8,A2,A8,A2,2I5,I6)', &
        IOSTAT=ios) &
        wl, gflog, code, e, xj, label, ep, xjp, labelp, &
        gr, gs, gw, ref, nblo, nbup, iso1, x1, iso2, x2, &
        other1, other2, lim, lim, isoshift
      IF (ios .LT. 0) EXIT
      IF (ios .GT. 0) THEN
        WRITE(6,'(a)') ' WARNING: read_gfall: read error, stopping line read'
        EXIT
      END IF

      READ(cother1, '(2I5)', IOSTAT=ios) ishift, ishiftp
      IF (ios .NE. 0) THEN
        ishift  = 0
        ishiftp = 0
      END IF
      READ(cother2, '(A6,I1,A3)', IOSTAT=ios) ixfixfp, linesize, auto_flag
      IF (ios .NE. 0) THEN
        linesize  = 0
        auto_flag = '   '
      END IF

      eshift  = ishift  * 0.001D0
      eshiftp = ishiftp * 0.001D0
      dwliso  = REAL(isoshift, 4) * 0.0001     ! milli-Angstrom -> nm

      wlvac = ABS(wl) + dwl
      IF (TRANSFER(labelp(1), ' ') .EQ. 'CONTINUU') &
        wlvac = 1.D7 / ABS(ABS(ep) + eshiftp - (ABS(e) + eshift)) + dwl + dwliso

      ! Guard: skip any record with non-physical wlvac (e.g. continuum edge
      ! where the energy level difference is zero giving wlvac=Inf or NaN)
      IF (.NOT. (wlvac .GT. 0.0D0 .AND. wlvac .LT. 1.0D8)) CYCLE

      IF (wlvac .GT. wlend + dellim(1)) EXIT

      ixwl  = INT(LOG(wlvac) / ratiolg + 0.5D0)
      nbuff = ixwl - ixwlbeg + 1

      icode_r = NINT(code * 100.0)

      lim = MIN(8 - linesize, 7)
      IF (icode_r .EQ. 100) lim = 1

      IF (wlvac .LT. wlbeg - dellim(lim) * delfactor) CYCLE
      IF (wlvac .GT. wlend + dellim(lim) * delfactor) CYCLE
      IF (auto_flag .EQ. 'COR') CYCLE

      ! --- oscillator strength and energy ---
      gf    = 10.0**(gflog + dgflog + x1 + x2)
      elo_d = DBLE(MIN(ABS(e), ABS(ep)))

      ! --- damping constants ---
      gammar_d = 10.0D0**(gr + dgammar)
      gammas_d = 10.0D0**(gs + dgammas)
      gammaw_d = 10.0D0**(gw + dgammaw)
      IF (auto_flag .EQ. 'AUT' .AND. gs .GT. 0.0) &
        gammas_d = -10.0D0**(-gs + dgammas)

      IF (ABS(gr) .LT. 1.0e-6) gammar_d = 2.223D13 / wlvac**2

      ! --- NELION ---
      nelem    = INT(code)
      icharge  = NINT((code - REAL(nelem,4)) * 100.0)
      zeff     = icharge + 1
      nelion_i = nelem*6 - 6 + INT(zeff)
      IF (nelem .GT. 19 .AND. nelem .LT. 29 .AND. icharge .GT. 5) &
        nelion_i = 6*(nelem + icharge*10 - 30) - 1

      ! --- default Stark width ---
      IF (ABS(gs) .LT. 1.0e-6) THEN
        IF (code .LT. 100.0) THEN
          eup    = DBLE(MAX(ABS(e), ABS(ep)))
          effnsq = 25.0D0
          CALL ionpot_index(nelem, icharge, eup, effnsq, zeff)
          gammas_d = 1.0D-8 * effnsq**2 * SQRT(effnsq)
          gs = REAL(LOG10(gammas_d), 4)
        ELSE
          gammas_d = 1.0D-5
        END IF
      END IF

      ! --- default van der Waals width ---
      IF (ABS(gw) .LT. 1.0e-6) THEN
        IF (code .LT. 100.0) THEN
          eup    = DBLE(MAX(ABS(e), ABS(ep)))
          effnsq = 25.0D0
          CALL ionpot_index(nelem, icharge, eup, effnsq, zeff)
          effnsq = MIN(effnsq, 1000.0D0)
          rsqup  = 2.5D0 * (effnsq / zeff)**2
          CALL ionpot_index_lo(nelem, icharge, elo_d, effnsq, zeff)
          effnsq = MIN(effnsq, 1000.0D0)
          rsqlo  = 2.5D0 * (effnsq / zeff)**2
          IF (code - zeff + 1.0D0 .GT. 20.0D0 .AND. code - zeff + 1.0D0 .LT. 29.0D0) THEN
            rsqup = (45.0D0 - (code - zeff + 1.0D0)) / zeff
            rsqlo = 0.0D0
          END IF
          BLOCK
            CHARACTER(LEN=8) :: clabelp_vdw
            clabelp_vdw = TRANSFER(labelp(1), clabelp_vdw)
            IF (clabelp_vdw .EQ. 'CONTINUU') rsqlo = 0.0D0
          END BLOCK
          IF (rsqup .LT. rsqlo) rsqup = 2.0D0 * rsqlo
          gammaw_d = 4.5D-9 * (rsqup - rsqlo)**0.4D0
        ELSE
          gammaw_d = 1.D-7 / zeff
        END IF
      END IF

      ! --- line type dispatch ---
      itype = 0
      IF (icode_r .EQ. 100)                     itype = -1
      IF (icode_r .EQ. 100 .AND. iso1 .EQ. 2)    itype = -2
      IF (icode_r .EQ. 200)                     itype = -3
      IF (icode_r .EQ. 200 .AND. iso1 .EQ. 3)    itype = -4
      IF (icode_r .EQ. 201)                     itype = -6
      IF (icode_r .EQ. 201 .AND. iso1 .EQ. 3)    itype = -6
      IF (auto_flag .EQ. 'AUT')                 itype =  1
      IF (auto_flag .EQ. 'PRD')                 itype =  3

      BLOCK
        CHARACTER(LEN=8) :: clabelp
        clabelp = TRANSFER(labelp(1), clabelp)
        IF (clabelp .EQ. 'CONTINUU') THEN
          itype = NINT(xjp)
          gf    = gf * (xj + xj + 1.0)
        END IF
      END BLOCK

      ncon = 0
      IF (iso1 .EQ. 0 .AND. iso2 .GT. 0) ncon = iso2

      ! --- CGF conversion and damping normalisation ---
      IF (itype .NE. 1 .AND. itype .LE. 3) THEN
        frelin   = 2.99792458D17 / wlvac
        cgf      = 0.026538D0 / 1.77245D0 * DBLE(gf) / frelin
        frq4pi   = 12.5664D0 * frelin
        IF (itype .EQ. 2) THEN
          gammar_d = DBLE(gr)
        ELSE
          gammar_d = gammar_d / frq4pi
          gammas_d = gammas_d / frq4pi
          gammaw_d = gammaw_d / frq4pi
        END IF
      END IF

      nbup = ABS(nbup)
      nblo = ABS(nblo)
      nelionx = 0

      ! --- append to appropriate buffer ---
      IF (itype .NE. 0) THEN
        ! complex-profile line -> nlte buffer
        IF (nblo + nbup .NE. 0) THEN
          DO ic = 1, 17
            IF (icode_r .EQ. codex(ic)) THEN
              nelionx = ic
              EXIT
            END IF
          END DO
          IF (nelionx .EQ. 0) THEN
            WRITE(6,'(a,f10.2)') ' WARNING: nlte line has unknown CODE: ', code
            CYCLE
          END IF
        END IF
        IF (n_nlte .EQ. nlte_cap) CALL grow_nlte(nlte_buf, nlte_cap)
        n_nlte = n_nlte + 1
        IF (itype .EQ. 1 .OR. itype .GT. 3) THEN
          nlte_buf(n_nlte) = nlte_line_t(wlvac, REAL(elo_d,4), REAL(gf,4), &
            nblo, nbup, nelion_i, itype, ncon, nelionx, &
            REAL(gammar_d,4), REAL(gammas_d,4), REAL(gammaw_d,4), nbuff, lim)
        ELSE
          nlte_buf(n_nlte) = nlte_line_t(wlvac, REAL(elo_d,4), REAL(cgf,4), &
            nblo, nbup, nelion_i, itype, ncon, nelionx, &
            REAL(gammar_d,4), REAL(gammas_d,4), REAL(gammaw_d,4), nbuff, lim)
        END IF

      ELSE IF (nblo + nbup .NE. 0) THEN
        ! TYPE=0 with departure coefficients -> nlte buffer
        DO ic = 1, 17
          IF (icode_r .EQ. codex(ic)) THEN
            nelionx = ic
            EXIT
          END IF
        END DO
        IF (nelionx .EQ. 0) THEN
          WRITE(6,'(a,f10.2)') ' WARNING: nlte line has unknown CODE: ', code
          CYCLE
        END IF
        IF (n_nlte .EQ. nlte_cap) CALL grow_nlte(nlte_buf, nlte_cap)
        n_nlte = n_nlte + 1
        nlte_buf(n_nlte) = nlte_line_t(wlvac, REAL(elo_d,4), REAL(cgf,4), &
          nblo, nbup, nelion_i, itype, ncon, nelionx, &
          REAL(gammar_d,4), REAL(gammas_d,4), REAL(gammaw_d,4), nbuff, lim)

      ELSE
        ! plain Voigt line -> lte buffer
        IF (n_lte .EQ. lte_cap) CALL grow_lte(lte_buf, lte_cap)
        n_lte = n_lte + 1
        lte_buf(n_lte) = lte_line_t(nbuff, REAL(cgf,4), nelion_i, &
          REAL(elo_d,4), REAL(gammar_d,4), REAL(gammas_d,4), REAL(gammaw_d,4))
      END IF

    END DO

    CLOSE(11)

    IF (ALLOCATED(lte_arr))  DEALLOCATE(lte_arr)
    IF (ALLOCATED(nlte_arr)) DEALLOCATE(nlte_arr)
    ALLOCATE(lte_arr(n_lte),   SOURCE=lte_buf(1:n_lte))
    ALLOCATE(nlte_arr(n_nlte), SOURCE=nlte_buf(1:n_nlte))
    DEALLOCATE(lte_buf, nlte_buf)

  END SUBROUTINE read_gfall


  ! ============================================================================
  !  READ_PREDICT — read Kurucz predicted-wavelength binary line list
  !                 (translated from rpredict.for)
  ! ============================================================================
  SUBROUTINE read_predict(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                           lte_arr, n_lte)
    CHARACTER(LEN=*),             INTENT(IN)    :: filename
    REAL(8),                      INTENT(IN)    :: wlbeg, wlend, ratiolg
    INTEGER,                      INTENT(IN)    :: ixwlbeg
    TYPE(lte_line_t), ALLOCATABLE, INTENT(INOUT) :: lte_arr(:)
    INTEGER,                      INTENT(INOUT) :: n_lte

    INTEGER(4) :: iiiiiii(4)
    INTEGER(4) :: iwl
    INTEGER(2) :: iwords(8)
    EQUIVALENCE (iiiiiii(1), iwl)
    EQUIVALENCE (iiiiiii(1), iwords(1))

    INTEGER :: nelionold(1005)
    INTEGER :: nelionolda(209), nelionoldb(286), nelionoldc(95)
    INTEGER :: nelionoldd(95),  nelionolde(95),  nelionoldf(60)
    INTEGER :: nelionoldg(165)
    EQUIVALENCE (nelionold(  1), nelionolda(1))
    EQUIVALENCE (nelionold(210), nelionoldb(1))
    EQUIVALENCE (nelionold(496), nelionoldc(1))
    EQUIVALENCE (nelionold(591), nelionoldd(1))
    EQUIVALENCE (nelionold(686), nelionolde(1))
    EQUIVALENCE (nelionold(781), nelionoldf(1))
    EQUIVALENCE (nelionold(841), nelionoldg(1))

    DATA nelionolda/ &
      1,  2, &
      7,  8,  9, &
     13, 14, 15, 16, &
     19, 20, 21, 22, 23, &
     25, 26, 27, 28, 29, 30, &
     31, 32, 33, 34, 35, 36,  0, &
     37, 38, 39, 40, 41, 42,  0,  0, &
     43, 44, 45, 46, 47, 48,  0,  0,  0, &
     49, 50, 51, 52, 53, 54,  0,  0,  0,  0, &
     55, 56, 57, 58, 59, 60,  0,  0,  0,  0,  0, &
     61, 62, 63, 64, 65, 66,  0,  0,  0,  0,  0,  0, &
     67, 68, 69, 70, 71, 72,  0,  0,  0,  0,  0,  0,  0, &
     73, 74, 75, 76, 77, 78,  0,  0,  0,  0,  0,  0,  0,  0, &
     79, 80, 81, 82, 83, 84,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     85, 86, 87, 88, 89, 90,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     91, 92, 93, 94, 95, 96,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     97, 98, 99,100,101,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
    103,104,105,106,107,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
    109,110,111,112,113,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
    DATA nelionoldb/ &
    115,116,117,118,119,120,299,359,419,479, 11*0, &
    121,122,123,124,125,126,305,365,425,485, 12*0, &
    127,128,129,130,131,132,311,371,431,491, 13*0, &
    133,134,135,136,137,138,317,377,437,497, 14*0, &
    139,140,141,142,143,144,323,383,443,503, 15*0, &
    145,146,147,148,149,150,329,389,449,509, 16*0, &
    151,152,153,154,155,156,335,395,455,515, 17*0, &
    157,158,159,160,161,162,341,401,461,521, 18*0, &
    163,164,165,166,167,168,347,407,467,527, 19*0, &
    169,170,171,             27*0, &
    175,176,177,             28*0/
    DATA nelionoldc/ &
    181,182,183, 0, 0, 187,188,189, 0, 0, 193,194,195, 0, 0, &
    199,200,201, 0, 0, 205,206,207, 0, 0, 211,212,213, 0, 0, &
    217,218,219, 0, 0, 223,224,225, 0, 0, 229,230,231, 0, 0, &
    235,236,237, 0, 0, 241,242,243, 0, 0, 247,248,249, 0, 0, &
    253,254,255, 0, 0, 259,260,261, 0, 0, 265,266,267, 0, 0, &
    271,272,273, 0, 0, 277,278,279, 0, 0, 283,284,285, 0, 0, &
    289,290,291, 0, 0/
    DATA nelionoldd/ &
    295,296,297, 0, 0, 301,302,303, 0, 0, 307,308,309, 0, 0, &
    313,314,315, 0, 0, 319,320,321, 0, 0, 325,326,327, 0, 0, &
    331,332,333, 0, 0, 337,338,339, 0, 0, 343,344,345, 0, 0, &
    349,350,351, 0, 0, 355,356,357, 0, 0, 361,362,363, 0, 0, &
    367,368,369, 0, 0, 373,374,375, 0, 0, 379,380,381, 0, 0, &
    385,386,387, 0, 0, 391,392,393, 0, 0, 397,398,399, 0, 0, &
    403,404,405, 0, 0/
    DATA nelionolde/ &
    409,410,411, 0, 0, 415,416,417, 0, 0, 421,422,423, 0, 0, &
    427,428,429, 0, 0, 433,434,435, 0, 0, 439,440,441, 0, 0, &
    445,446,447, 0, 0, 451,452,453, 0, 0, 457,458,459, 0, 0, &
    463,464,465, 0, 0, 469,470,471, 0, 0, 475,476,477, 0, 0, &
    481,482,483, 0, 0, 487,488,489, 0, 0, 493,494,495, 0, 0, &
    499,500,501, 0, 0, 505,506,507, 0, 0, 511,512,513, 0, 0, &
    517,518,519, 0, 0/
    DATA nelionoldf/ &
    523,524,525, 0, 0, 529,530,531, 0, 0, 535,536,537, 0, 0, &
    541,542,543, 0, 0, 547,548,549, 0, 0, 553,554,555, 0, 0, &
    559,560,561, 0, 0, 565,566,567, 0, 0, 571,572,573, 0, 0, &
    577,578,579, 0, 0, 583,584,585, 0, 0, 589,590,591, 0, 0/
    DATA nelionoldg/ &
    240,  0,378,384,390,246,252,258,396,  0, &
    300,306,312,402,336,408,  0,342,414,420, &
    426,432,438,444,558,564,570,264,270,276, &
      0,  0,  0,  0,282,288,  0,  0,  0,  0, &
      0,492,498,294,  0,  0,318,324,330,504, &
    348,510,354,360,366,372,516,522,528,576, &
    582,588,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,450,456,462,  0,  0, &
    468,474,480,  0,  0,  0,486,  0,  0,  0, &
    594,  0,  0,  0,  0,  0,  0,  0,  0,534, &
    540,546,  0,  0,552,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0/

    REAL(4), SAVE :: tablog(32768)
    REAL(8) :: wlvac_d, elo_d, gf_d, freq, congf, frq4pi
    REAL(8) :: gamrf, gamsf, gamwf, ratiolog_r2e6
    INTEGER :: ios, i, istart, istop, limitblue, limitred, newlimit
    INTEGER :: lengthfile, nstart, iwl1, nelionnew, nelion_i, nbuff
    INTEGER :: ixwl

    INTEGER, PARAMETER :: CHUNK = 200000
    TYPE(lte_line_t), ALLOCATABLE :: lte_buf(:)
    INTEGER :: lte_cap

    DO i = 1, 32768
      tablog(i) = 10.0**((i - 16384) * 0.001)
    END DO

    ratiolog_r2e6 = LOG(1.0D0 + 1.0D0/2000000.0D0)

    istart = INT(LOG(wlbeg - 1.0D0) / ratiolog_r2e6 + 0.5D0)
    istop  = INT(LOG(wlend + 1.0D0) / ratiolog_r2e6 + 0.5D0)

    OPEN(UNIT=11, FILE=filename, STATUS='old', FORM='unformatted', &
         ACCESS='direct', RECL=16, IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') 'ERROR: cannot open predict file: ', TRIM(filename)
      CALL EXIT(1)
    END IF

    READ(11, REC=1) iwl1
    IF (iwl1 .GT. istop) THEN
      IF (VERBOSE .EQ. 1) &
        WRITE(6,'(a)') '  No predicted lines in window (file starts past wlend).'
      CLOSE(11)
      IF (ALLOCATED(lte_arr)) DEALLOCATE(lte_arr)
      ALLOCATE(lte_arr(0))
      n_lte = 0
      RETURN
    END IF

    BLOCK
      INTEGER(8) :: fsize
      INQUIRE(FILE=filename, SIZE=fsize)
      lengthfile = INT(fsize / 16)
    END BLOCK

    READ(11, REC=lengthfile) iwl
    IF (iwl .LT. istart) THEN
      IF (VERBOSE .EQ. 1) &
        WRITE(6,'(a)') '  No predicted lines in window (file ends before wlbeg).'
      CLOSE(11)
      IF (ALLOCATED(lte_arr)) DEALLOCATE(lte_arr)
      ALLOCATE(lte_arr(0))
      n_lte = 0
      RETURN
    END IF

    ! binary search for first record >= istart
    limitblue = 1
    limitred  = lengthfile
    DO WHILE (limitred - limitblue .GT. 1)
      newlimit = (limitred + limitblue) / 2
      READ(11, REC=newlimit) iwl
      IF (iwl .LT. istart) THEN
        limitblue = newlimit
      ELSE
        limitred  = newlimit
      END IF
    END DO
    nstart = limitred

    lte_cap = CHUNK
    ALLOCATE(lte_buf(lte_cap))
    n_lte = 0

    DO i = nstart, lengthfile
      READ(11, REC=i, IOSTAT=ios) iiiiiii
      IF (ios .LT. 0) EXIT
      IF (ios .GT. 0) THEN
        WRITE(6,'(a,i0)') ' ERROR: read_predict: read error at record ', i
        EXIT
      END IF
      IF (iwl .GT. istop) EXIT

      nelionnew = ABS(iwords(3)) / 10
      IF (nelionnew .LT. 1 .OR. nelionnew .GT. 1005) CYCLE
      nelion_i  = nelionold(nelionnew)
      IF (nelion_i .EQ. 0) CYCLE

      wlvac_d = EXP(iwl * ratiolog_r2e6)
      freq    = 2.99792458D17 / wlvac_d
      gf_d    = tablog(iwords(5))
      congf   = 0.01502D0 * gf_d / freq
      elo_d   = tablog(iwords(4))
      frq4pi  = freq * 12.5664D0
      gamrf   = tablog(iwords(6)) / frq4pi
      gamsf   = tablog(iwords(7)) / frq4pi
      gamwf   = tablog(iwords(8)) / frq4pi

      ixwl  = INT(LOG(wlvac_d) / ratiolg + 0.5D0)
      nbuff = ixwl - ixwlbeg + 1

      IF (n_lte .EQ. lte_cap) CALL grow_lte(lte_buf, lte_cap)
      n_lte = n_lte + 1
      lte_buf(n_lte) = lte_line_t(nbuff, REAL(congf,4), nelion_i, &
        REAL(elo_d,4), REAL(gamrf,4), REAL(gamsf,4), REAL(gamwf,4))
    END DO

    CLOSE(11)

    IF (ALLOCATED(lte_arr)) DEALLOCATE(lte_arr)
    ALLOCATE(lte_arr(n_lte), SOURCE=lte_buf(1:n_lte))
    DEALLOCATE(lte_buf)

  END SUBROUTINE read_predict


  ! ============================================================================
  !  READ_MOLEC — dispatcher: binary sidecar if present, else ASCII
  ! ============================================================================
  SUBROUTINE read_molec(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                         lte_arr, n_lte)
    CHARACTER(LEN=*),             INTENT(IN)    :: filename
    REAL(8),                      INTENT(IN)    :: wlbeg, wlend, ratiolg
    INTEGER,                      INTENT(IN)    :: ixwlbeg
    TYPE(lte_line_t), ALLOCATABLE, INTENT(INOUT) :: lte_arr(:)
    INTEGER,                      INTENT(INOUT) :: n_lte

    CHARACTER(LEN=512) :: binfile
    LOGICAL :: have_bin
    INTEGER :: dot

    ! Build sidecar path: replace extension with .bin
    binfile = TRIM(filename)
    dot     = INDEX(binfile, '.', BACK=.TRUE.)
    IF (dot .GT. 0) THEN
      binfile = binfile(:dot-1) // '.bin'
    ELSE
      binfile = TRIM(binfile) // '.bin'
    END IF
    INQUIRE(FILE=binfile, EXIST=have_bin)

    IF (have_bin) THEN
      CALL read_molec_bin(filename, binfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                          lte_arr, n_lte)
    ELSE
      CALL read_molec_ascii(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                            lte_arr, n_lte)
    END IF

  END SUBROUTINE read_molec


  ! ============================================================================
  !  READ_MOLEC_BIN — binary sidecar path with O(log N) window seek
  ! ============================================================================
  SUBROUTINE read_molec_bin(ascfile, binfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                             lte_arr, n_lte)
    CHARACTER(LEN=*),             INTENT(IN)    :: ascfile, binfile
    REAL(8),                      INTENT(IN)    :: wlbeg, wlend, ratiolg
    INTEGER,                      INTENT(IN)    :: ixwlbeg
    TYPE(lte_line_t), ALLOCATABLE, INTENT(INOUT) :: lte_arr(:)
    INTEGER,                      INTENT(INOUT) :: n_lte

    INTEGER, PARAMETER :: MAGIC    = INT(z'4D4F4C32')
    INTEGER, PARAMETER :: RECBYTES = 32

    INTEGER(4) :: hdr(8), buf(8)
    REAL(4)    :: wlvac, gflog, e_r4, ep_r4
    INTEGER(4) :: icode, iso, loggr, labelp_x
    INTEGER    :: nrec, ios, ixwl, nbuff, lo, hi, mid, istart
    REAL(8)    :: wlvac_d, elo, gf, freq, congf, frq4pi
    REAL(8)    :: gammar, gammas, gammaw, gamrf, gamsf, gamwf
    REAL(4)    :: fudge, x1, x2
    INTEGER    :: nelion, iso1, iso2

    INTEGER, PARAMETER :: CHUNK = 200000
    TYPE(lte_line_t), ALLOCATABLE :: lte_buf(:)
    ! n_lte is a running total across molecule files; preserve existing contents
    TYPE(lte_line_t), ALLOCATABLE :: lte_old(:)
    INTEGER :: n_old, n_new, lte_cap

    OPEN(UNIT=21, FILE=binfile, STATUS='old', FORM='unformatted', &
         ACCESS='direct', RECL=RECBYTES, IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') '  WARNING: cannot open binary sidecar, falling back: ', &
        TRIM(binfile)
      CALL read_molec_ascii(ascfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                            lte_arr, n_lte)
      RETURN
    END IF

    READ(21, REC=1) hdr
    IF (hdr(1) .NE. MAGIC) THEN
      WRITE(6,'(a)') '  WARNING: binary sidecar wrong magic, falling back'
      CLOSE(21)
      CALL read_molec_ascii(ascfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                            lte_arr, n_lte)
      RETURN
    END IF
    nrec = hdr(2)

    ! Binary search: first record with wlvac >= wlbeg - 1.0
    lo = 1
    hi = nrec
    DO WHILE (lo .LT. hi)
      mid = (lo + hi) / 2
      READ(21, REC=mid+1, IOSTAT=ios) buf
      IF (ios .NE. 0) THEN
        WRITE(6,'(a,i0)') ' ERROR: read_molec_bin seek error at record ', mid+1
        CLOSE(21)
        RETURN
      END IF
      wlvac = TRANSFER(buf(1), wlvac)
      IF (wlvac .LT. REAL(wlbeg - 1.0D0, 4)) THEN
        lo = mid + 1
      ELSE
        hi = mid
      END IF
    END DO
    istart = lo

    lte_cap = CHUNK
    ALLOCATE(lte_buf(lte_cap))
    n_new = 0

    DO mid = istart, nrec
      READ(21, REC=mid+1, IOSTAT=ios) buf
      IF (ios .LT. 0) EXIT
      IF (ios .GT. 0) THEN
        WRITE(6,'(a,i0)') ' ERROR: read_molec_bin read error at record ', mid+1
        EXIT
      END IF

      wlvac    = TRANSFER(buf(1), wlvac)
      gflog    = TRANSFER(buf(2), gflog)
      e_r4     = TRANSFER(buf(3), e_r4)
      ep_r4    = TRANSFER(buf(4), ep_r4)
      icode    = buf(5)
      iso      = buf(6)
      loggr    = buf(7)
      labelp_x = buf(8)

      wlvac_d = DBLE(wlvac)
      IF (wlvac_d .GT. wlend + 1.0D0) EXIT

      CALL molec_dispatch(iso, icode, nelion, iso1, iso2, x1, x2, fudge)
      IF (nelion .EQ. 0) CYCLE

      gf  = EXP((DBLE(gflog) + x1 + x2 + fudge) * 2.30258509299405D0)
      elo = MIN(ABS(DBLE(e_r4)), ABS(DBLE(ep_r4)))

      ixwl   = INT(LOG(wlvac_d) / ratiolg + 0.5D0)
      nbuff  = ixwl - ixwlbeg + 1
      freq   = 2.99792458D17 / wlvac_d
      congf  = 0.026538D0 / 1.77245D0 * gf / freq
      frq4pi = freq * 12.5664D0

      IF (loggr .EQ. 0) THEN
        gammar = 2.223D13 / wlvac_d**2
      ELSE
        gammar = 10.0D0**(loggr * 0.01D0)
      END IF

      IF (labelp_x .EQ. 1) THEN
        gammas = 3.0D-8
        gammaw = 1.0D-8
      ELSE
        gammas = 3.0D-5
        gammaw = 1.0D-7
      END IF

      gamrf = gammar / frq4pi
      gamsf = gammas / frq4pi
      gamwf = gammaw / frq4pi

      IF (n_new .EQ. lte_cap) CALL grow_lte(lte_buf, lte_cap)
      n_new = n_new + 1
      lte_buf(n_new) = lte_line_t(nbuff, REAL(congf,4), nelion, &
        REAL(elo,4), REAL(gamrf,4), REAL(gamsf,4), REAL(gamwf,4))
    END DO

    CLOSE(21)

    ! Append to existing lte_arr (which may already have lines from prior mol files)
    n_old = n_lte
    IF (n_old .GT. 0 .AND. ALLOCATED(lte_arr)) THEN
      ALLOCATE(lte_old(n_old), SOURCE=lte_arr(1:n_old))
    END IF
    IF (ALLOCATED(lte_arr)) DEALLOCATE(lte_arr)
    n_lte = n_old + n_new
    ALLOCATE(lte_arr(n_lte))
    IF (n_old .GT. 0) lte_arr(1:n_old) = lte_old(1:n_old)
    IF (n_new .GT. 0) lte_arr(n_old+1:n_lte) = lte_buf(1:n_new)
    IF (ALLOCATED(lte_old)) DEALLOCATE(lte_old)
    DEALLOCATE(lte_buf)

  END SUBROUTINE read_molec_bin


  ! ============================================================================
  !  READ_MOLEC_ASCII — sequential ASCII scan with window guards
  ! ============================================================================
  SUBROUTINE read_molec_ascii(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                               lte_arr, n_lte)
    CHARACTER(LEN=*),             INTENT(IN)    :: filename
    REAL(8),                      INTENT(IN)    :: wlbeg, wlend, ratiolg
    INTEGER,                      INTENT(IN)    :: ixwlbeg
    TYPE(lte_line_t), ALLOCATABLE, INTENT(INOUT) :: lte_arr(:)
    INTEGER,                      INTENT(INOUT) :: n_lte

    REAL(8)          :: wl, e, ep, wlvac, elo, gf, freq, congf, frq4pi
    REAL(8)          :: gammar, gammas, gammaw, gamrf, gamsf, gamwf
    REAL(4)          :: gflog, xj, xjp, fudge, x1, x2
    INTEGER          :: loggr, icode, iso, nelion, iso1, iso2
    INTEGER          :: ios, ixwl, nbuff
    CHARACTER(LEN=8) :: clabel, clabelp

    INTEGER, PARAMETER :: CHUNK = 200000
    TYPE(lte_line_t), ALLOCATABLE :: lte_buf(:)
    TYPE(lte_line_t), ALLOCATABLE :: lte_old(:)
    INTEGER :: n_old, n_new, lte_cap

    OPEN(UNIT=11, FILE=filename, STATUS='old', ACTION='read', &
         FORM='formatted', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') 'ERROR: cannot open molecule file: ', TRIM(filename)
      CALL EXIT(1)
    END IF

    lte_cap = CHUNK
    ALLOCATE(lte_buf(lte_cap))
    n_new = 0

    DO
      READ(11, '(F10.4,F7.3,F5.1,F10.3,F5.1,F11.3,I4,A8,A8,I2,I4)', &
           IOSTAT=ios) wl, gflog, xj, e, xjp, ep, icode, clabel, clabelp, iso, loggr
      IF (ios .LT. 0) EXIT
      IF (ios .GT. 0) THEN
        WRITE(6,'(a,i0)') ' ERROR: read_molec_ascii: read error near line ', n_new
        EXIT
      END IF

      IF (ABS(wl) .GT. wlend + 2.0D0) EXIT

      wlvac = 1.0D7 / ABS(ABS(ep) - ABS(e))
      IF (wlvac .LT. wlbeg - 1.0D0) CYCLE
      IF (wlvac .GT. wlend + 1.0D0) CYCLE

      CALL molec_dispatch(iso, icode, nelion, iso1, iso2, x1, x2, fudge)
      IF (nelion .EQ. 0) CYCLE

      gf  = EXP((gflog + x1 + x2 + fudge) * 2.30258509299405D0)
      elo = MIN(ABS(e), ABS(ep))

      ixwl   = INT(LOG(wlvac) / ratiolg + 0.5D0)
      nbuff  = ixwl - ixwlbeg + 1
      freq   = 2.99792458D17 / wlvac
      congf  = 0.026538D0 / 1.77245D0 * gf / freq
      frq4pi = freq * 12.5664D0

      IF (loggr .EQ. 0) THEN
        gammar = 2.223D13 / wlvac**2
      ELSE
        gammar = 10.0D0**(loggr * 0.01D0)
      END IF

      IF (clabelp(1:1) .EQ. 'X') THEN
        gammas = 3.0D-8
        gammaw = 1.0D-8
      ELSE
        gammas = 3.0D-5
        gammaw = 1.0D-7
      END IF

      gamrf = gammar / frq4pi
      gamsf = gammas / frq4pi
      gamwf = gammaw / frq4pi

      IF (n_new .EQ. lte_cap) CALL grow_lte(lte_buf, lte_cap)
      n_new = n_new + 1
      lte_buf(n_new) = lte_line_t(nbuff, REAL(congf,4), nelion, &
        REAL(elo,4), REAL(gamrf,4), REAL(gamsf,4), REAL(gamwf,4))
    END DO

    CLOSE(11)

    n_old = n_lte
    IF (n_old .GT. 0 .AND. ALLOCATED(lte_arr)) THEN
      ALLOCATE(lte_old(n_old), SOURCE=lte_arr(1:n_old))
    END IF
    IF (ALLOCATED(lte_arr)) DEALLOCATE(lte_arr)
    n_lte = n_old + n_new
    ALLOCATE(lte_arr(n_lte))
    IF (n_old .GT. 0) lte_arr(1:n_old) = lte_old(1:n_old)
    IF (n_new .GT. 0) lte_arr(n_old+1:n_lte) = lte_buf(1:n_new)
    IF (ALLOCATED(lte_old)) DEALLOCATE(lte_old)
    DEALLOCATE(lte_buf)

  END SUBROUTINE read_molec_ascii


  ! ============================================================================
  !  READ_H2O — read Partridge-Schwenke H2O binary line list
  !              (translated from rh2ofast.for)
  ! ============================================================================
  SUBROUTINE read_h2o(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                       lte_arr, n_lte)
    CHARACTER(LEN=*),             INTENT(IN)    :: filename
    REAL(8),                      INTENT(IN)    :: wlbeg, wlend, ratiolg
    INTEGER,                      INTENT(IN)    :: ixwlbeg
    TYPE(lte_line_t), ALLOCATABLE, INTENT(INOUT) :: lte_arr(:)
    INTEGER,                      INTENT(INOUT) :: n_lte

    INTEGER(4) :: irec(2)
    INTEGER(4) :: iwl
    INTEGER(2) :: ielo_i
    EQUIVALENCE (irec(1), iwl)
    EQUIVALENCE (irec(2), ielo_i)

    REAL(4), SAVE :: tablog(32768)
    REAL(4), PARAMETER :: xiso(4)  = [ 0.9976,  0.0004,  0.0020, 0.00001]

    REAL(8) :: wlvac, wlvac1, freq, congf, ratiolog_r2e6, frq4pi
    REAL(8) :: gammar, gamrf, gamsf, gamwf
    INTEGER :: ios, istart, i, iso, ixwl, nbuff
    INTEGER :: limitblue, limitred, newlimit, lengthfile
    INTEGER(2) :: igflog_loc, ielo_loc

    INTEGER, PARAMETER :: CHUNK = 500000
    TYPE(lte_line_t), ALLOCATABLE :: lte_buf(:)
    INTEGER :: n_new, lte_cap

    DO i = 1, 32768
      tablog(i) = 10.0**((i - 16384) * 0.001)
    END DO

    ratiolog_r2e6 = LOG(1.0D0 + 1.0D0/2000000.0D0)

    OPEN(UNIT=11, FILE=filename, STATUS='old', FORM='unformatted', &
         ACCESS='direct', RECL=8, IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') 'ERROR: cannot open H2O file: ', TRIM(filename)
      CALL EXIT(1)
    END IF

    READ(11, REC=1) irec
    wlvac = EXP(iwl * ratiolog_r2e6)
    IF (wlvac .GT. wlend + 1.0D0) THEN
      CLOSE(11)
      IF (.NOT. ALLOCATED(lte_arr)) ALLOCATE(lte_arr(0))
      RETURN
    END IF

    BLOCK
      INTEGER(8) :: fsize
      INQUIRE(FILE=filename, SIZE=fsize)
      lengthfile = INT(fsize / 8)
    END BLOCK

    READ(11, REC=lengthfile) irec
    wlvac1 = EXP(iwl * ratiolog_r2e6)
    IF (wlbeg - 1.0D0 .GT. wlvac1) THEN
      CLOSE(11)
      IF (.NOT. ALLOCATED(lte_arr)) ALLOCATE(lte_arr(0))
      RETURN
    END IF

    limitblue = 1
    limitred  = lengthfile
    DO WHILE (limitred - limitblue .GT. 1)
      newlimit = (limitred + limitblue) / 2
      READ(11, REC=newlimit) irec
      wlvac = EXP(iwl * ratiolog_r2e6)
      IF (wlvac .LT. wlbeg - 1.0D0) THEN
        limitblue = newlimit
      ELSE
        limitred  = newlimit
      END IF
    END DO
    istart = newlimit

    lte_cap = CHUNK
    ALLOCATE(lte_buf(lte_cap))
    n_new  = 0

    DO i = istart, lengthfile
      READ(11, REC=i, IOSTAT=ios) irec
      IF (ios .LT. 0) EXIT
      IF (ios .GT. 0) THEN
        WRITE(6,'(a,i0)') ' ERROR: read_h2o: read error at record ', i
        EXIT
      END IF

      ielo_loc   = INT(IBITS(irec(2),  0, 16), 2)
      igflog_loc = INT(IBITS(irec(2), 16, 16), 2)

      wlvac = EXP(iwl * ratiolog_r2e6)
      IF (wlvac .GT. wlend + 1.0D0) EXIT

      freq  = 2.99792458D17 / wlvac
      ixwl  = INT(LOG(wlvac) / ratiolg + 0.5D0)
      nbuff = ixwl - ixwlbeg + 1

      IF      (ielo_loc .GT. 0 .AND. igflog_loc .GT. 0)  THEN
        iso = 1
      ELSE IF (ielo_loc .GT. 0 .AND. igflog_loc .LE. 0) THEN
        iso = 2
      ELSE IF (ielo_loc .LE. 0 .AND. igflog_loc .GT. 0) THEN
        iso = 3
      ELSE
        iso = 4
      END IF

      congf  = 0.01502D0 * tablog(ABS(igflog_loc)) / freq * xiso(iso)
      frq4pi = freq * 12.5664D0
      gammar = 2.223D13 / wlvac**2 * 0.001D0
      gamrf  = gammar / frq4pi
      gamsf  = tablog(1)    / frq4pi
      gamwf  = tablog(9384) / frq4pi

      IF (n_new .EQ. lte_cap) CALL grow_lte(lte_buf, lte_cap)
      n_new = n_new + 1
      lte_buf(n_new) = lte_line_t(nbuff, REAL(congf,4), 534, &
        REAL(ABS(ielo_loc),4), REAL(gamrf,4), REAL(gamsf,4), REAL(gamwf,4))
    END DO

    CLOSE(11)

    IF (ALLOCATED(lte_arr)) DEALLOCATE(lte_arr)
    ALLOCATE(lte_arr(n_new), SOURCE=lte_buf(1:n_new))
    n_lte = n_new
    DEALLOCATE(lte_buf)

  END SUBROUTINE read_h2o


  ! ============================================================================
  !  MOLEC_DISPATCH — map ISO isotope index to NELION and abundance corrections
  ! ============================================================================
  SUBROUTINE molec_dispatch(iso, icode, nelion, iso1, iso2, x1, x2, fudge)
    INTEGER, INTENT(IN)  :: iso, icode
    INTEGER, INTENT(OUT) :: nelion, iso1, iso2
    REAL(4), INTENT(OUT) :: x1, x2, fudge

    nelion = 0; iso1 = 0; iso2 = 0; x1 = 0.0; x2 = 0.0; fudge = 0.0

    SELECT CASE (iso)
    CASE (1);   nelion=240; iso1=1;  iso2=1;  x1= 0.0;    x2=-5.0
    CASE (2);   nelion=240; iso1=1;  iso2=2;  x1= 0.0;    x2=-4.469
    CASE (12)
      SELECT CASE (icode)
      CASE (606); nelion=264; iso1=12; iso2=12; x1=-.005; x2=-.005
      CASE (608); nelion=276; iso1=12; iso2=16; x1=-.005; x2=-.001
      CASE (106); nelion=246; iso1=1;  iso2=12; x1= 0.0;  x2=-.005
      CASE DEFAULT
                  nelion=270; iso1=12; iso2=14; x1=-.005; x2=-.002
      END SELECT
    CASE (13)
      SELECT CASE (icode)
      CASE (606); nelion=264; iso1=12; iso2=13; x1=-.005;  x2=-1.955
      CASE (608); nelion=276; iso1=13; iso2=16; x1=-1.955; x2=-.001
      CASE (106); nelion=246; iso1=1;  iso2=13; x1= 0.0;   x2=-1.955
      CASE DEFAULT
                  nelion=270; iso1=13; iso2=14; x1=-1.955; x2=-.002
      END SELECT
    CASE (14);  nelion=252; iso1=1;  iso2=14; x1= 0.0;   x2=-.002
    CASE (15)
      IF (icode .EQ. 607) THEN
                  nelion=270; iso1=12; iso2=15; x1=-.005; x2=-2.444
      ELSE;       nelion=252; iso1=1;  iso2=15; x1= 0.0;  x2=-2.444
      END IF
    CASE (16)
      IF (icode .EQ. 813) THEN
                  nelion=324; iso1=27; iso2=16; x1= 0.0;  x2=-.001
      ELSE;       nelion=258; iso1=1;  iso2=16; x1= 0.0;  x2=-.001
      END IF
    CASE (17)
      IF (icode .EQ. 813) THEN
                  nelion=324; iso1=27; iso2=17; x1= 0.0;  x2=-3.398
      ELSE;       nelion=276; iso1=12; iso2=17; x1=-.005; x2=-3.398
      END IF
    CASE (18)
      SELECT CASE (icode)
      CASE (814); nelion=330; iso1=28; iso2=18; x1=-.035; x2=-2.690
      CASE (608); nelion=276; iso1=12; iso2=18; x1=-.005; x2=-2.690
      CASE (813); nelion=324; iso1=27; iso2=18; x1= 0.0;  x2=-2.690
      CASE DEFAULT
                  nelion=258; iso1=1;  iso2=18; x1= 0.0;  x2=-2.690
      END SELECT
    CASE (23);  nelion=492; iso1=1;  iso2=23; x1= 0.0;   x2= 0.0
    CASE (24)
      IF (icode .EQ. 812) THEN
                  nelion=318; iso1=16; iso2=24; x1=-.001; x2=-.102
      ELSE;       nelion=300; iso1=1;  iso2=24; x1= 0.0;  x2=-.105
      END IF
    CASE (25)
      IF (icode .EQ. 812) THEN
                  nelion=318; iso1=16; iso2=25; x1=-.001; x2=-1.000
      ELSE;       nelion=300; iso1=1;  iso2=25; x1= 0.0;  x2=-.996
      END IF
    CASE (26)
      IF (icode .EQ. 812) THEN
                  nelion=318; iso1=16; iso2=26; x1=-.001; x2=-.958
      ELSE;       nelion=300; iso1=1;  iso2=26; x1= 0.0;  x2=-.947
      END IF
    CASE (28)
      IF (icode .EQ. 814) THEN
                  nelion=330; iso1=28; iso2=16; x1=-.035; x2=-.001
      ELSE;       nelion=312; iso1=1;  iso2=28; x1= 0.0;  x2=-.035
      END IF
    CASE (29)
      IF (icode .EQ. 814) THEN
                  nelion=330; iso1=29; iso2=16; x1=-1.328; x2=-.001
      ELSE;       nelion=312; iso1=1;  iso2=29; x1= 0.0;   x2=-1.331
      END IF
    CASE (30)
      IF (icode .EQ. 814) THEN
                  nelion=330; iso1=30; iso2=16; x1=-1.510; x2=-.001
      ELSE;       nelion=312; iso1=1;  iso2=30; x1= 0.0;   x2=-1.516
      END IF
    CASE (33);  nelion=264; iso1=13; iso2=13; x1=-1.955; x2=-1.955
    CASE (39);  nelion=498; iso1=39; iso2=1;  x1=-.030;  x2= 0.0
    CASE (40)
      IF (icode .EQ. 820) THEN
                  nelion=354; iso1=40; iso2=16; x1=-.013; x2=-.001
      ELSE;       nelion=342; iso1=40; iso2=1;  x1=-.013; x2= 0.0
      END IF
    CASE (41);  nelion=498; iso1=41; iso2=1;  x1=-1.172; x2= 0.0
    CASE (42);  nelion=342; iso1=42; iso2=1;  x1=-2.189; x2= 0.0
    CASE (43);  nelion=342; iso1=43; iso2=1;  x1=-2.870; x2= 0.0
    CASE (44);  nelion=342; iso1=44; iso2=1;  x1=-1.681; x2= 0.0
    CASE (46)
      IF (icode .EQ. 120) THEN
                  nelion=342; iso1=46; iso2=1;  x1=-4.398; x2= 0.0
      ELSE;       nelion=366; iso1=16; iso2=46; x1= 0.0;   x2=-1.101
      END IF
    CASE (47);  nelion=366; iso1=16; iso2=47; x1= 0.0;   x2=-1.138
    CASE (48)
      IF (icode .EQ. 120) THEN
                  nelion=342; iso1=48; iso2=1;  x1=-2.728; x2= 0.0
      ELSE;       nelion=366; iso1=16; iso2=48; x1= 0.0;   x2=-.131
      END IF
    CASE (49);  nelion=366; iso1=16; iso2=49; x1= 0.0;   x2=-1.259
    CASE (50)
      IF (icode .EQ. 124) THEN
                  nelion=432; iso1=50; iso2=1;  x1=-1.362; x2= 0.0
      ELSE;       nelion=366; iso1=16; iso2=50; x1= 0.0;   x2=-1.272
      END IF
    CASE (51);  nelion=372; iso1=16; iso2=51; x1= 0.0;   x2=-.001
    CASE (52);  nelion=432; iso1=52; iso2=1;  x1=-.077;  x2= 0.0
    CASE (53);  nelion=432; iso1=53; iso2=1;  x1=-1.022; x2= 0.0
    CASE (54)
      IF (icode .EQ. 124) THEN
                  nelion=432; iso1=54; iso2=1;  x1=-1.626; x2= 0.0
      ELSE;       nelion=444; iso1=54; iso2=1;  x1=-1.237; x2= 0.0
      END IF
    CASE (56);  nelion=444; iso1=56; iso2=1;  x1=-.038;  x2= 0.0
    CASE (57);  nelion=444; iso1=57; iso2=1;  x1=-1.658; x2= 0.0
    CASE (58);  nelion=444; iso1=58; iso2=1;  x1=-2.553; x2= 0.0
    CASE DEFAULT
      nelion = 0
    END SELECT
  END SUBROUTINE molec_dispatch


  ! ============================================================================
  !  GROW_LTE / GROW_NLTE — double buffer capacity
  ! ============================================================================
  SUBROUTINE grow_lte(buf, cap)
    TYPE(lte_line_t), ALLOCATABLE, INTENT(INOUT) :: buf(:)
    INTEGER,                       INTENT(INOUT) :: cap
    TYPE(lte_line_t), ALLOCATABLE :: tmp(:)
    ALLOCATE(tmp(cap*2))
    tmp(1:cap) = buf(1:cap)
    CALL MOVE_ALLOC(tmp, buf)
    cap = cap * 2
  END SUBROUTINE grow_lte

  SUBROUTINE grow_nlte(buf, cap)
    TYPE(nlte_line_t), ALLOCATABLE, INTENT(INOUT) :: buf(:)
    INTEGER,                        INTENT(INOUT) :: cap
    TYPE(nlte_line_t), ALLOCATABLE :: tmp(:)
    ALLOCATE(tmp(cap*2))
    tmp(1:cap) = buf(1:cap)
    CALL MOVE_ALLOC(tmp, buf)
    cap = cap * 2
  END SUBROUTINE grow_nlte


  ! ============================================================================
  !  IONPOT_INDEX — upper-level effective principal quantum number
  ! ============================================================================
  SUBROUTINE ionpot_index(nelem, icharge, eup, effnsq, zeff)
    INTEGER, INTENT(IN)    :: nelem, icharge
    REAL(8), INTENT(IN)    :: eup, zeff
    REAL(8), INTENT(INOUT) :: effnsq   ! in: safety default; out: computed value
    INTEGER :: idx
    REAL(8) :: deleup

    IF (nelem .LE. 30) THEN
      idx = nelem*(nelem+1)/2 + icharge
    ELSE
      idx = nelem*5 + 341 + icharge
    END IF
    deleup = potion(idx) - eup
    IF (deleup .GT. 0.0D0) THEN
      effnsq = 109737.31D0 * zeff**2 / deleup
    END IF
    ! else: retain the caller-supplied safety default (25.0)
  END SUBROUTINE ionpot_index


  ! ============================================================================
  !  IONPOT_INDEX_LO — lower-level effective principal quantum number
  ! ============================================================================
  SUBROUTINE ionpot_index_lo(nelem, icharge, elo, effnsq, zeff)
    INTEGER, INTENT(IN)  :: nelem, icharge
    REAL(8), INTENT(IN)  :: elo, zeff
    REAL(8), INTENT(OUT) :: effnsq
    INTEGER :: idx
    REAL(8) :: delelo

    IF (nelem .LE. 30) THEN
      idx = nelem*(nelem+1)/2 + icharge
    ELSE
      idx = nelem*5 + 341 + icharge
    END IF
    delelo = potion(idx) - elo
    effnsq = 109737.31D0 * zeff**2 / delelo
    effnsq = MIN(effnsq, 1000.0D0)
  END SUBROUTINE ionpot_index_lo



  ! ============================================================================
  !  BASENAME — return the filename portion of a path (after last '/')
  ! ============================================================================
  pure FUNCTION basename(path) RESULT(name)
    CHARACTER(LEN=*), INTENT(IN) :: path
    CHARACTER(LEN=512)           :: name
    INTEGER :: i
    i = INDEX(path, '/', BACK=.TRUE.)
    IF (i .GT. 0) THEN
      name = path(i+1:)
    ELSE
      name = path
    END IF
  END FUNCTION basename


  ! ============================================================================
  !  IS_TIO_FILE — does the basename begin with "tio" (case-insensitive)?
  !
  !  Used to gate TiO line lists out of hot-star runs.  Kurucz TiO filenames
  !  lead with the species symbol (tio*.asc, tiopred.bin); the prefix match
  !  avoids false positives like "ratio" / "station".
  ! ============================================================================
  pure FUNCTION is_tio_file(path) RESULT(isTiO)
    CHARACTER(LEN=*), INTENT(IN) :: path
    LOGICAL                      :: isTiO
    CHARACTER(LEN=512) :: bname
    CHARACTER(LEN=3)   :: prefix
    INTEGER :: k, c

    bname = basename(path)
    prefix = bname(1:3)
    ! Lowercase in place (ASCII only; Kurucz filenames are ASCII).
    DO k = 1, 3
      c = IACHAR(prefix(k:k))
      IF (c .GE. IACHAR('A') .AND. c .LE. IACHAR('Z')) THEN
        prefix(k:k) = ACHAR(c + 32)
      END IF
    END DO
    isTiO = (prefix .EQ. 'tio')
  END FUNCTION is_tio_file


  ! ============================================================================
  !  IONPOTS — fill the ionisation potential table (cm^-1)
  ! ============================================================================
  SUBROUTINE ionpots()
    potion = 0.0D0
    ! H, He
    potion(  1) = 109678.772D0
    potion(  3) = 198310.666D0; potion(  4) = 438908.879D0
    ! Li
    potion(  6) =  43487.114D0; potion(  7) = 610078.526D0
    potion(  8) = 987661.014D0
    ! Be
    potion( 10) =  75192.640D0; potion( 11) = 146882.86D0
    potion( 12) = 1241256.600D0; potion( 13) = 1756018.822D0
    ! B
    potion( 15) =  66928.040D0; potion( 16) = 202887.40D0
    potion( 17) = 305930.80D0;  potion( 18) = 2091972.D0
    potion( 19) = 2744107.936D0
    ! C
    potion( 21) =  90820.42D0;  potion( 22) = 196674.D0
    potion( 23) = 386241.0D0;   potion( 24) = 520175.8D0
    potion( 25) = 3162423.30D0; potion( 26) = 3952061.670D0
    ! N
    potion( 28) = 117225.70D0;  potion( 29) = 238750.20D0
    potion( 30) = 382672.D0;    potion( 31) = 624866.D0
    potion( 32) = 789537.D0;    potion( 33) = 4452723.30D0
    potion( 34) = 5380089.80D0
    ! O
    potion( 36) = 109837.02D0;  potion( 37) = 283270.90D0
    potion( 38) = 443085.0D0;   potion( 39) = 624382.0D0
    potion( 40) = 918657.D0;    potion( 41) = 1114004.D0
    potion( 42) = 5963073.00D0; potion( 43) = 7028394.70D0
    ! F
    potion( 45) = 140524.50D0;  potion( 46) = 282058.6D0
    potion( 47) = 505774.0D0;   potion( 48) = 703110.D0
    potion( 49) = 921480.D0;    potion( 50) = 1267606.0D0
    potion( 51) = 1493632.D0;   potion( 52) = 7693706.60D0
    potion( 53) = 8897242.50D0
    ! Ne
    potion( 55) = 173929.750D0; potion( 56) = 330388.60D0
    potion( 57) = 511544.D0;    potion( 58) = 783890.D0
    potion( 59) = 1018250.D0;   potion( 60) = 1273820.D0
    potion( 61) = 1671750.D0;   potion( 62) = 1928447.D0
    potion( 63) = 9644840.7D0;  potion( 64) = 10986877.20D0
    ! Na
    potion( 66) =  41449.451D0; potion( 67) = 381390.2D0
    potion( 68) = 577654.D0;    potion( 69) = 797970.D0
    potion( 70) = 1116300.D0;   potion( 71) = 1389100.D0
    potion( 72) = 1681700.D0;   potion( 73) = 2130850.D0
    potion( 74) = 2418500.D0;   potion( 75) = 11817106.70D0
    potion( 76) = 13297680.0D0
    ! Mg
    potion( 78) =  61671.050D0; potion( 79) = 121267.61D0
    potion( 80) = 646402.D0;    potion( 81) = 881285.D0
    potion( 82) = 1139900.D0;   potion( 83) = 1506300.D0
    potion( 84) = 1814900.D0;   potion( 85) = 2144820.D0
    potion( 86) = 2645400.D0;   potion( 87) = 2964000.D0
    potion( 88) = 14209914.7D0; potion( 89) = 15829950.D0
    ! Al
    potion( 91) =  48278.48D0;  potion( 92) = 151862.50D0
    potion( 93) = 229445.70D0;  potion( 94) = 967804.D0
    potion( 95) = 1240684.D0;   potion( 96) = 1536400.D0
    potion( 97) = 1949900.D0;   potion( 98) = 2295800.D0
    potion( 99) = 2663300.D0;   potion(100) = 3215300.D0
    potion(101) = 3565010.D0;   potion(102) = 16824539.3D0
    potion(103) = 18584143.0D0
    ! Si
    potion(105) =  65747.76D0;  potion(106) = 131838.14D0
    potion(107) = 270139.30D0;  potion(108) = 364093.10D0
    potion(109) = 1345070.D0;   potion(110) = 1655590.D0
    potion(111) = 1986700.D0;   potion(112) = 2449200.D0
    potion(113) = 2831800.D0;   potion(114) = 3237400.D0
    potion(115) = 3840600.D0;   potion(116) = 4221630.D0
    potion(117) = 19661038.9D0; potion(118) = 21560631.0D0
    ! P
    potion(120) =  84580.83D0;  potion(121) = 159451.70D0
    potion(122) = 243600.70D0;  potion(123) = 414922.8D0
    potion(124) = 524462.9D0;   potion(125) = 1777890.D0
    potion(126) = 2125800.D0;   potion(127) = 2497100.D0
    potion(128) = 3002900.D0;   potion(129) = 3423000.D0
    potion(130) = 3867000.D0;   potion(131) = 4521700.D0
    potion(132) = 4934020.D0;   potion(133) = 22719901.6D0
    potion(134) = 24759942.D0
    ! S
    potion(136) =  83559.1D0;   potion(137) = 188232.7D0
    potion(138) = 281100.D0;    potion(139) = 380870.D0
    potion(140) = 585514.D0;    potion(141) = 710195.D0
    potion(142) = 2266050.D0;   potion(143) = 2651900.D0
    potion(144) = 3063600.D0;   potion(145) = 3611300.D0
    potion(146) = 4069500.D0;   potion(147) = 4552200.D0
    potion(148) = 5258400.D0;   potion(149) = 5702290.D0
    potion(150) = 26001545.1D0; potion(151) = 28182526.D0
    ! Cl
    potion(153) = 104591.00D0;  potion(154) = 192070.0D0
    potion(155) = 321000.D0;    potion(156) = 429400.D0
    potion(157) = 545800.D0;    potion(158) = 781900.D0
    potion(159) = 921096.D0;    potion(160) = 2809280.D0
    potion(161) = 3233080.D0;   potion(162) = 3683000.D0
    potion(163) = 4274000.D0;   potion(164) = 4771400.D0
    potion(165) = 5293400.D0;   potion(166) = 6051000.D0
    potion(167) = 6526620.D0;   potion(168) = 29506532.5D0
    potion(169) = 31828983.D0
    ! Ar
    potion(171) = 127109.842D0; potion(172) = 222848.3D0
    potion(173) = 328550.D0;    potion(174) = 480600.D0
    potion(175) = 603700.D0;    potion(176) = 736300.D0
    potion(177) = 1003400.D0;   potion(178) = 1157056.D0
    potion(179) = 3408500.D0;   potion(180) = 3869500.D0
    potion(181) = 4359000.D0;   potion(182) = 4992000.D0
    potion(183) = 5528700.D0;   potion(184) = 6090500.D0
    potion(185) = 6899800.D0;   potion(186) = 7407190.D0
    potion(187) = 33235410.D0;  potion(188) = 35699895.D0
    ! K
    potion(190) =  35009.814D0; potion(191) = 255072.8D0
    potion(192) = 369427.D0;    potion(193) = 491330.D0
    potion(194) = 666700.D0;    potion(195) = 802000.D0
    potion(196) = 948200.D0;    potion(197) = 1249100.D0
    potion(198) = 1418063.D0;   potion(199) = 4062400.D0
    potion(200) = 4562000.D0;   potion(201) = 5090000.D0
    potion(202) = 5764000.D0;   potion(203) = 6342000.D0
    potion(204) = 6943800.D0;   potion(205) = 7805000.D0
    potion(206) = 8344140.D0;   potion(207) = 37189176.0D0
    potion(208) = 39795784.D0
    ! Ca
    potion(210) =  49305.924D0; potion(211) =  95751.870D0
    potion(212) = 410642.3D0;   potion(213) = 542595.D0
    potion(214) = 680200.D0;    potion(215) = 877400.D0
    potion(216) = 1026000.D0;   potion(217) = 1187600.D0
    potion(218) = 1520600.D0;   potion(219) = 1704050.D0
    potion(220) = 4771600.D0;   potion(221) = 5309000.D0
    potion(222) = 5877000.D0;   potion(223) = 6591000.D0
    potion(224) = 7210000.D0;   potion(225) = 7853000.D0
    potion(226) = 8766000.D0;   potion(227) = 9337690.D0
    potion(228) = 41367028.D0;  potion(229) = 44117409.D0
    ! Sc
    potion(231) =  52922.00D0;  potion(232) = 103237.1D0
    potion(233) = 199677.37D0;  potion(234) = 592732.D0
    potion(235) = 741600.D0;    potion(236) = 892700.D0
    potion(237) = 1113000.D0;   potion(238) = 1275000.D0
    potion(239) = 1452000.D0;   potion(240) = 1816200.D0
    potion(241) = 2014760.D0;   potion(242) = 5543900.D0
    potion(243) = 6111000.D0;   potion(244) = 6720000.D0
    potion(245) = 7473000.D0;   potion(246) = 8135000.D0
    potion(247) = 8820000.D0;   potion(248) = 9784000.D0
    potion(249) = 10388070.D0;  potion(250) = 45771185.D0
    potion(251) = 48665510.D0
    ! Ti
    potion(253) =  55072.50D0;  potion(254) = 109494.D0
    potion(255) = 221735.6D0;   potion(256) = 348973.3D0
    potion(257) = 800900.D0;    potion(258) = 964100.D0
    potion(259) = 1134700.D0;   potion(260) = 1375000.D0
    potion(261) = 1549000.D0;   potion(262) = 1741500.D0
    potion(263) = 2137900.D0;   potion(264) = 2351110.D0
    potion(265) = 6353000.D0;   potion(266) = 6969000.D0
    potion(267) = 7618000.D0;   potion(268) = 8408000.D0
    potion(269) = 9116000.D0;   potion(270) = 9842000.D0
    potion(271) = 10859000.D0;  potion(272) = 11495470.D0
    potion(273) = 50401766.D0;  potion(274) = 53440740.D0
    ! V
    potion(276) =  54411.67D0;  potion(277) = 117900.D0
    potion(278) = 236410.D0;    potion(279) = 376730.D0
    potion(280) = 526532.0D0;   potion(281) = 1033400.D0
    potion(282) = 1215700.D0;   potion(283) = 1399800.D0
    potion(284) = 1661000.D0;   potion(285) = 1859000.D0
    potion(286) = 2055000.D0;   potion(287) = 2488200.D0
    potion(288) = 2712230.D0;   potion(289) = 7227000.D0
    potion(290) = 7882000.D0;   potion(291) = 8573000.D0
    potion(292) = 9398000.D0;   potion(293) = 10153000.D0
    potion(294) = 10922000.D0;  potion(295) = 11991000.D0
    potion(296) = 12660130.D0;  potion(297) = 55259549.D0
    potion(298) = 58443920.D0
    ! Cr
    potion(300) =  54575.6D0;   potion(301) = 132971.02D0
    potion(302) = 249700.D0;    potion(303) = 396500.D0
    potion(304) = 560200.D0;    potion(305) = 731020.D0
    potion(306) = 1292800.D0;   potion(307) = 1490200.D0
    potion(308) = 1690100.D0;   potion(309) = 1972000.D0
    potion(310) = 2184000.D0;   potion(311) = 2393000.D0
    potion(312) = 2860500.D0;   potion(313) = 3098480.D0
    potion(314) = 8159000.D0;   potion(315) = 8850000.D0
    potion(316) = 9582000.D0;   potion(317) = 10443000.D0
    potion(318) = 11247000.D0;  potion(319) = 12059000.D0
    potion(320) = 13180000.D0;  potion(321) = 13882280.D0
    potion(322) = 60345293.D0;  potion(323) = 63675850.D0
    ! Mn
    potion(325) =  59959.4D0;   potion(326) = 126145.00D0
    potion(327) = 271550.D0;    potion(328) = 413000.D0
    potion(329) = 584000.D0;    potion(330) = 771100.D0
    potion(331) = 961440.D0;    potion(332) = 1577000.D0
    potion(333) = 1789600.D0;   potion(334) = 2005400.D0
    potion(335) = 2308000.D0;   potion(336) = 2536000.D0
    potion(337) = 2771000.D0;   potion(338) = 3250000.D0
    potion(339) = 3509900.D0;   potion(340) = 9144000.D0
    potion(341) = 9873000.D0;   potion(342) = 10649000.D0
    potion(343) = 11541000.D0;  potion(344) = 12398000.D0
    potion(345) = 13253000.D0;  potion(346) = 14427000.D0
    potion(347) = 15162200.D0;  potion(348) = 65659877.D0
    potion(349) = 69137430.D0
    ! Fe
    potion(351) =  63737.704D0; potion(352) = 130655.40D0
    potion(353) = 247220.D0;    potion(354) = 442900.D0
    potion(355) = 604900.D0;    potion(356) = 798370.D0
    potion(357) = 1008000.D0;   potion(358) = 1218380.D0
    potion(359) = 1884000.D0;   potion(360) = 2114000.D0
    potion(361) = 2346000.D0;   potion(362) = 2668000.D0
    potion(363) = 2912000.D0;   potion(364) = 3163000.D0
    potion(365) = 3680000.D0;   potion(366) = 3946570.D0
    potion(367) = 10184000.D0;  potion(368) = 10951000.D0
    potion(369) = 11770000.D0;  potion(370) = 12708000.D0
    potion(371) = 13607000.D0;  potion(372) = 14505000.D0
    potion(373) = 15731000.D0;  potion(374) = 16500160.D0
    potion(375) = 71204137.D0;  potion(376) = 74829550.D0
  END SUBROUTINE ionpots

END MODULE mod_mklinelist
