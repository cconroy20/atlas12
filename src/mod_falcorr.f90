! ==============================================================================
! mod_falcorr.f90 — per-line wavelength / log gf / vdW corrections from FAL
!
! Public interface:
!   call falcorr_load(path)                       ! once, before run_mklinelist
!   found = falcorr_lookup(src, icode, wlvac, elo, gflog, iso1, iso2, &
!                          dwl, dgf, dgw)         ! from the line readers
!   call falcorr_report(wlbeg, wlend)             ! once, after run_mklinelist
!   falcorr_strict                                ! .FALSE.: unmatched lines warn
!                                                 !          instead of stopping
!
! FAL (FALv3, falfit.jl) fits a correction (dlambda, dlog gf, dlog gamma_vdW) to
! each strong atomic line by matching Korg syntheses to six benchmark stars.
! It fits them against the synthe2026 master line list, which is a full-
! precision dump of THIS code's readers (korg/ai/mkll/mkmasterll.f90).  So each
! corrected line can be found again here by the quantities read_gfall and
! read_predict already compute, without any shared row index:
!
!   code   species code x 100 (2600 = Fe I, 2601 = Fe II)
!   wlvac  vacuum wavelength [nm], BEFORE any correction
!   elo    REAL(4) lower-level energy [cm^-1]
!   gflog  raw REAL(4) log gf, isotope fractions NOT folded in (gfall column;
!          for gfpred, log10 of the tablog value)
!   iso1, iso2
!
! cgf is deliberately not part of the key: read_gfall forms gf in REAL(4) and
! mkmasterll in REAL(8), so the two differ in the last bits.
!
! The catalog is built by the FAL project (Phillip Cargile) with its
! code/postproc/mkfalcorr.py.  Plain text:
!   '#' comment lines, then one line per corrected line,
!     wlvac_nm code elo_cm gflog iso1 iso2 dwl_nm dloggf dloggw [ignored ...]
!   sorted by wlvac_nm.  The corrections are applied exactly as falfit does:
!     wlvac += dwl_nm,  log gf += dloggf,  gamma_vdW *= 10**dloggw
!
! Only the atomic readers (gfall, gfpred) are hooked.  A catalog may carry
! molecular rows (code >= 100, mkfalcorr.py --include-mol); they are loaded and
! counted in the report but never applied, and do not trip the unmatched check.
!
! Nothing here runs unless falcorr= is given: with no catalog loaded the lookup
! returns .FALSE. on its first test, so the synthesis is bit-identical.
! ==============================================================================

MODULE mod_falcorr
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: falcorr_load, falcorr_lookup, falcorr_report
  PUBLIC :: falcorr_active, falcorr_strict
  PUBLIC :: FC_SRC_GFALL, FC_SRC_PREDICT

  ! Identity tolerances.  The tightest real separation they must resolve is
  ! the wavelength step of the sources (1e-4 nm gfall, 2.5e-4 nm gfpred) and
  ! the 0.001 dex step of gfall's F7.3 / gfpred's tablog.  Keep in lockstep
  ! with TOL_* in mkfalcorr.py.
  REAL(8), PARAMETER :: TOL_WL    = 1.0D-6   ! nm
  REAL(8), PARAMETER :: TOL_ELO   = 1.0D-2   ! cm^-1
  REAL(8), PARAMETER :: TOL_GFLOG = 5.0D-4   ! dex

  ! Which reader matched a catalog line, for the report.
  INTEGER, PARAMETER :: FC_SRC_GFALL = 1, FC_SRC_PREDICT = 2, FC_NSRC = 2
  CHARACTER(LEN=6), PARAMETER :: FC_SRC_NAME(FC_NSRC) = ['gfall ', 'gfpred']

  LOGICAL :: falcorr_active = .FALSE.
  LOGICAL :: falcorr_strict = .TRUE.
  INTEGER :: nfc = 0

  REAL(8), ALLOCATABLE :: fc_wl(:), fc_elo(:), fc_gflog(:)
  REAL(8), ALLOCATABLE :: fc_dwl(:), fc_dgf(:), fc_dgw(:)
  INTEGER, ALLOCATABLE :: fc_code(:), fc_iso1(:), fc_iso2(:)
  INTEGER, ALLOCATABLE :: fc_hits(:)   ! times each entry was matched
  INTEGER, ALLOCATABLE :: fc_src(:)    ! reader of the first match, 0 = none

CONTAINS

  ! ============================================================================
  !  FALCORR_LOAD — read the catalog; hard stop on any malformed or unsorted row
  ! ============================================================================
  SUBROUTINE falcorr_load(path)
    CHARACTER(LEN=*), INTENT(IN) :: path
    CHARACTER(LEN=512) :: line
    INTEGER :: u, ios, n, pass, i, nmol
    REAL(8) :: wl, code, elo, gflog, dwl, dgf, dgw
    INTEGER :: iso1, iso2

    OPEN(NEWUNIT=u, FILE=TRIM(path), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(6,'(a,a)') ' ERROR: cannot open falcorr catalog: ', TRIM(path)
      CALL EXIT(1)
    END IF

    ! Two passes: count, then fill.
    DO pass = 1, 2
      REWIND(u)
      n = 0
      DO
        READ(u, '(A)', IOSTAT=ios) line
        IF (ios .NE. 0) EXIT
        IF (LEN_TRIM(line) .EQ. 0) CYCLE
        IF (line(1:1) .EQ. '#') THEN
          IF (pass .EQ. 2 .AND. (INDEX(line, '# tag') .EQ. 1 .OR. &
              INDEX(line, '# master') .EQ. 1 .OR. INDEX(line, '# built') .EQ. 1)) &
            WRITE(6,'(a,a)') '   ', TRIM(line(3:))
          CYCLE
        END IF
        READ(line, *, IOSTAT=ios) wl, code, elo, gflog, iso1, iso2, dwl, dgf, dgw
        IF (ios .NE. 0) THEN
          WRITE(6,'(a,a)') ' ERROR: malformed falcorr line: ', TRIM(line)
          CALL EXIT(1)
        END IF
        n = n + 1
        IF (pass .EQ. 2) THEN
          fc_wl(n)    = wl
          fc_code(n)  = NINT(code * 100.0D0)
          fc_elo(n)   = elo
          fc_gflog(n) = gflog
          fc_iso1(n)  = iso1
          fc_iso2(n)  = iso2
          fc_dwl(n)   = dwl
          fc_dgf(n)   = dgf
          fc_dgw(n)   = dgw
        END IF
      END DO
      IF (pass .EQ. 1) THEN
        nfc = n
        ALLOCATE(fc_wl(n), fc_code(n), fc_elo(n), fc_gflog(n), fc_iso1(n), &
                 fc_iso2(n), fc_dwl(n), fc_dgf(n), fc_dgw(n), fc_hits(n), fc_src(n))
        fc_hits = 0
        fc_src  = 0
        WRITE(6,'(a,a)') ' FALCORR catalog  = ', TRIM(path)
      END IF
    END DO
    CLOSE(u)

    ! falcorr_lookup bisects on wavelength.
    DO i = 2, nfc
      IF (fc_wl(i) .LT. fc_wl(i-1)) THEN
        WRITE(6,'(a,i0)') ' ERROR: falcorr catalog not sorted by wavelength at line ', i
        CALL EXIT(1)
      END IF
    END DO

    falcorr_active = nfc .GT. 0
    IF (nfc .GT. 0) THEN
      nmol = COUNT(fc_code .GE. 10000)
      WRITE(6,'(a,i0,a,i0,a,i0,a)') '   ', nfc, ' catalog lines (', nfc - nmol, &
        ' atomic, ', nmol, ' molecular)'
      WRITE(6,'(a,f9.3,a,f9.3,a)') '   covering ', fc_wl(1), ' - ', fc_wl(nfc), ' nm'
      IF (nmol .GT. 0) WRITE(6,'(a)') &
        '   molecular rows are reported but NOT applied (no molecular reader hook)'
      IF (.NOT. falcorr_strict) WRITE(6,'(a)') &
        '   falcorr_strict = no: unmatched catalog lines will warn, not stop'
    ELSE
      WRITE(6,'(a)') '   catalog is empty; no corrections will be applied'
    END IF
  END SUBROUTINE falcorr_load


  ! ============================================================================
  !  FALCORR_LOOKUP — corrections for one line, .FALSE. (and zeros) if none
  ! ============================================================================
  LOGICAL FUNCTION falcorr_lookup(src, icode, wlvac, elo, gflog, iso1, iso2, &
                                  dwl, dgf, dgw)
    INTEGER, INTENT(IN)  :: src, icode, iso1, iso2
    REAL(8), INTENT(IN)  :: wlvac, elo, gflog
    REAL(8), INTENT(OUT) :: dwl, dgf, dgw
    INTEGER :: lo, hi, mid, i

    falcorr_lookup = .FALSE.
    dwl = 0.0D0
    dgf = 0.0D0
    dgw = 0.0D0
    IF (nfc .EQ. 0) RETURN
    IF (wlvac .LT. fc_wl(1) - TOL_WL .OR. wlvac .GT. fc_wl(nfc) + TOL_WL) RETURN

    ! first entry with fc_wl >= wlvac - TOL_WL
    lo = 1
    hi = nfc
    DO WHILE (lo .LT. hi)
      mid = (lo + hi) / 2
      IF (fc_wl(mid) .LT. wlvac - TOL_WL) THEN
        lo = mid + 1
      ELSE
        hi = mid
      END IF
    END DO

    DO i = lo, nfc
      IF (fc_wl(i) .GT. wlvac + TOL_WL) EXIT
      IF (fc_code(i) .NE. icode) CYCLE
      IF (fc_iso1(i) .NE. iso1 .OR. fc_iso2(i) .NE. iso2) CYCLE
      IF (ABS(fc_elo(i) - elo) .GT. TOL_ELO) CYCLE
      IF (ABS(fc_gflog(i) - gflog) .GT. TOL_GFLOG) CYCLE
      dwl = fc_dwl(i)
      dgf = fc_dgf(i)
      dgw = fc_dgw(i)
      fc_hits(i) = fc_hits(i) + 1
      IF (fc_src(i) .EQ. 0) fc_src(i) = src
      falcorr_lookup = .TRUE.
      RETURN
    END DO
  END FUNCTION falcorr_lookup


  ! ============================================================================
  !  FALCORR_REPORT — what was matched, what was applied, what was missed
  !
  !  Every atomic catalog line inside [wlbeg, wlend] must have been read by
  !  gfall or gfpred, both of which read past the window edges.  One that was
  !  not means the catalog was built against a different line-list generation
  !  than the one lines.list points at, and the corrections that DID match may
  !  be on the wrong lines too.  That stops the run unless falcorr_strict=no.
  ! ============================================================================
  SUBROUTINE falcorr_report(wlbeg, wlend)
    REAL(8), INTENT(IN) :: wlbeg, wlend
    INTEGER :: i, k, nin, nmol_in, nhit, nmiss, nmulti, nsynth, nshown
    INTEGER :: nsrc(FC_NSRC), nstrong, nweak, nvdw
    REAL(8) :: sabs_dwl, max_dwl, min_dgf, max_dgf, sabs_dgf, min_dgw, max_dgw
    CHARACTER(LEN=*), PARAMETER :: ROW = '(3x,a,t40,i9,a)'

    IF (.NOT. falcorr_active) RETURN

    nin = 0; nmol_in = 0; nhit = 0; nmulti = 0; nsynth = 0; nsrc = 0
    nstrong = 0; nweak = 0; nvdw = 0
    sabs_dwl = 0.0D0; max_dwl = 0.0D0; sabs_dgf = 0.0D0
    min_dgf = HUGE(1.0D0); max_dgf = -HUGE(1.0D0)
    min_dgw = HUGE(1.0D0); max_dgw = -HUGE(1.0D0)
    DO i = 1, nfc
      IF (fc_wl(i) .LT. wlbeg .OR. fc_wl(i) .GT. wlend) CYCLE
      IF (fc_code(i) .GE. 10000) THEN
        nmol_in = nmol_in + 1
        CYCLE
      END IF
      nin = nin + 1
      IF (fc_hits(i) .EQ. 0) CYCLE
      nhit   = nhit + 1
      nsynth = nsynth + fc_hits(i)
      IF (fc_hits(i) .GT. 1) nmulti = nmulti + 1
      nsrc(fc_src(i)) = nsrc(fc_src(i)) + 1
      sabs_dwl = sabs_dwl + ABS(fc_dwl(i))
      max_dwl  = MAX(max_dwl, ABS(fc_dwl(i)))
      sabs_dgf = sabs_dgf + ABS(fc_dgf(i))
      min_dgf  = MIN(min_dgf, fc_dgf(i))
      max_dgf  = MAX(max_dgf, fc_dgf(i))
      IF (fc_dgf(i) .GT. 0.0D0) nstrong = nstrong + 1
      IF (fc_dgf(i) .LT. 0.0D0) nweak = nweak + 1
      IF (fc_dgw(i) .NE. 0.0D0) THEN
        nvdw = nvdw + 1
        min_dgw = MIN(min_dgw, fc_dgw(i))
        max_dgw = MAX(max_dgw, fc_dgw(i))
      END IF
    END DO
    nmiss = nin - nhit

    WRITE(6,'(a)') ' FALCORR summary ---------------------------------------------'
    WRITE(6,ROW) 'catalog lines, all wavelengths', nfc, ''
    WRITE(6,ROW) 'atomic lines in synthesis window', nin, ''
    WRITE(6,ROW) '  matched', nhit, ''
    DO k = 1, FC_NSRC
      WRITE(6,ROW) '    by ' // TRIM(FC_SRC_NAME(k)) // ' reader', nsrc(k), ''
    END DO
    WRITE(6,ROW) '  NOT matched', nmiss, MERGE('   <-- see below', '                ', nmiss .GT. 0)
    WRITE(6,ROW) 'line-list records corrected', nsynth, &
      '   (incl. duplicate source records)'
    WRITE(6,ROW) '  catalog lines matching >1 record', nmulti, ''
    WRITE(6,ROW) 'molecular in window, not applied', nmol_in, ''
    IF (nhit .GT. 0) THEN
      WRITE(6,'(3x,a,t40,f9.3,a,f8.3,a)') 'applied |dlambda| [mA]  mean / max', &
        1.0D4 * sabs_dwl / nhit, ' /', 1.0D4 * max_dwl, ''
      WRITE(6,'(3x,a,t40,f9.3,a,f8.3,a,f8.3)') 'applied dlog gf  mean|.| / min / max', &
        sabs_dgf / nhit, ' /', min_dgf, ' /', max_dgf
      WRITE(6,'(3x,a,t40,i9,a,i8,a)') '  weakened / strengthened', nweak, ' /', nstrong, ''
      IF (nvdw .GT. 0) THEN
        WRITE(6,'(3x,a,t40,i9,a,f8.3,a,f8.3)') 'applied dlog gvdW  n / min / max', &
          nvdw, ' /', min_dgw, ' /', max_dgw
      ELSE
        WRITE(6,ROW) 'applied dlog gvdW  n', 0, ''
      END IF
    END IF
    WRITE(6,'(a)') ' -------------------------------------------------------------'
    IF (nmiss .EQ. 0) RETURN

    WRITE(6,'(a,i0,a,i0,a)') ' ', nmiss, ' of ', nin, &
      ' atomic catalog lines in the synthesis window were not found in the'
    WRITE(6,'(a)') ' line list.  The catalog was probably built against a different'
    WRITE(6,'(a)') ' line-list generation than lines.list points at.  First few:'
    WRITE(6,'(a)') '        wlvac_nm      code    elo_cm      gflog  iso1 iso2     dloggf'
    nshown = 0
    DO i = 1, nfc
      IF (fc_wl(i) .LT. wlbeg .OR. fc_wl(i) .GT. wlend) CYCLE
      IF (fc_code(i) .GE. 10000 .OR. fc_hits(i) .GT. 0) CYCLE
      WRITE(6,'(8x,f14.7,i6,f12.3,f9.3,2i5,f11.4)') fc_wl(i), fc_code(i), fc_elo(i), &
        fc_gflog(i), fc_iso1(i), fc_iso2(i), fc_dgf(i)
      nshown = nshown + 1
      IF (nshown .GE. 10) EXIT
    END DO
    IF (falcorr_strict) THEN
      WRITE(6,'(a)') ' ERROR: stopping.  Rebuild the catalog against this line list,'
      WRITE(6,'(a)') '        or pass falcorr_strict=no to continue without these lines.'
      CALL EXIT(1)
    END IF
    WRITE(6,'(a)') ' WARNING: continuing (falcorr_strict=no); the lines above are uncorrected.'
  END SUBROUTINE falcorr_report

END MODULE mod_falcorr
