! ============================================================================
!  PROGRAM SYNTHE
!
!  Merged XNFPELSYN + SYNTHE + SPECTRV spectral-synthesis driver.  All
!  inter-stage I/O is in-memory via module arrays; see synthe_module and
!  mod_mklinelist.  ATLAS library routines are accessed through
!  mod_atlas_data (atlas12_modules.f90).
!
!  Units:
!     5  input  : model atmosphere cards (CLI argument, read by READIN)
!    11  output : ASCII spectrum <model>.spec
!    17  input  : continua.dat continuum edge frequency list (from DATADIR)
!    33  output : <model>.linform (wavelength / flux / continuum / tau table)
!    35  output : <model>.mol molecular number density table
!
!  <model>.spec format, one line per wavelength point:
!    cols  1-11 : wavelength (Angstroms),        F11.4
!    cols 12-26 : flux or specific intensity,    E15.6
!    cols 27-41 : continuum flux,                E15.6
! ============================================================================

PROGRAM SYNTHE
  
  USE synthe_module
  USE mod_mklinelist, only: run_mklinelist, lte_lines, nlines_lte, nlines_nlte
  USE mod_atlas_data, only: &
    JOSH, READIN, BLOCKJ, BLOCKH, set_bc_data_dir, &
    DATADIR, IFSYNTHE, &
    ! Renamed to avoid collision with module-level names in synthe_module
    hkt_a => HKT, itemp_a => ITEMP, &
    nrhox_a => NRHOX, rhox_a => RHOX, &
    ! Direct imports
    ACONT, ALINE, BNU, DELTAW, EHVKT, FREQ, FREQLG, &
    HNU, IFSCAT, IFSURF, NMU, NUHI, NULO, NUMNU, &
    SCONT, SIGMAC, SIGMAL, SLINE, STIM, &
    SURFI, TAUNU, TEFF
  IMPLICIT NONE

  ! --- Source-function fudge parameters ----------------------------------
  !
  ! LINE_SCAT_RHOX_SCALE -- column-mass scale height (g/cm^2) for an
  ! empirical split of line opacity between thermal absorption (S = B_nu)
  ! and coherent scattering (S = J_nu).  Per depth:
  !
  !     f_scat(j) = exp(-rhox(j) / LINE_SCAT_RHOX_SCALE)      (0 if scale=0)
  !     ALINE (j) = alpha_line(j) * (1 - f_scat(j))
  !     SIGMAL(j) = alpha_line(j) *      f_scat(j)
  !
  ! Deep in the photosphere (rhox >> scale) f_scat -> 0, recovering pure
  ! LTE line opacity.  High up (rhox << scale) f_scat -> 1, driving the
  ! source function toward J_nu and mimicking scattering-dominated cores.
  ! Typical "on" value: 0.1-1.0 g/cm^2; ship value is 0 (pure LTE).  Set
  ! nonzero here and recompile to enable.
  !
  ! PH1, PC1, PSI1 -- exponents for a multiplicative NLTE-like fudge on
  ! the LTE line source function, using ground-state Boltzmann populations
  ! of H I, C I/II, Si I/II (bhyd_gs, bc1_gs/bc2_gs, bsi1_gs/bsi2_gs):
  !
  !     bfudge(j) = bhyd_gs**PH1 * (bc1_gs/bc2_gs)**PC1 * (bsi1_gs/bsi2_gs)**PSI1
  !     SLINE (j) = BNU(j) * STIM(j) / (bfudge(j) - EHVKT(j))
  !
  ! All zeros => bfudge = 1, standard LTE source function.  Nonzero only
  ! for tuning experiments.
  ! -----------------------------------------------------------------------
  REAL(8), PARAMETER :: LINE_SCAT_RHOX_SCALE = 0.0D0
  REAL(8), PARAMETER :: PH1                  = 0.0D0
  REAL(8), PARAMETER :: PC1                  = 0.0D0
  REAL(8), PARAMETER :: PSI1                 = 0.0D0

  CHARACTER(LEN=*), PARAMETER :: USAGE = &
    ' Usage: synthe_spectrv.exe <model_file> wlbeg=<nm> wlend=<nm> [resolu=<R>] [turbv=<kms>]'

  ! --- Fixed run parameters (formerly read from fort.93 / SYNBEG) -------
  !   ifvac  = 1    : vacuum wavelengths throughout
  !   cutoff = 1e-3 : line-wing truncation threshold (drop when local
  !                   opacity < cutoff * continuum)
  INTEGER, PARAMETER :: IFVAC  = 1
  REAL(4), PARAMETER :: CUTOFF = 1.0E-3

  ! --- Run parameters set from CLI --------------------------------------
  REAL(4)   :: turbv

  ! --- Model atmosphere header ------------------------------------------
  INTEGER   :: nedge
  REAL(8)   :: wledge(377), deledge(377), halfedge(377)

  ! --- Per-depth atmosphere read buffers --------------------------------
  REAL(8)   :: qdopple(mw6), qxnfpel(mw6)
  REAL(8)   :: qablog(1131)
  REAL(8)   :: ablog(3, 377)

  ! --- In-memory opacity matrix, (length, nrhox) ------------------------
  REAL(4), ALLOCATABLE :: opacity_matrix(:,:)

  ! --- Scalar work variables --------------------------------------------
  INTEGER   :: i, j, nu, iedge, nbuff_i, eqpos
  INTEGER   :: iline, iwave
  INTEGER   :: maxred, maxblue, minblue_i, istep
  INTEGER   :: n10dop, nstep, n1, maxstep
  INTEGER   :: nbuff_s, congf_nel, nelem_i
  REAL(4)   :: kappa0_s, kapmin_s, kapcen_s
  REAL(4)   :: adamp_s, congf_s, dvoigt, x_wing
  REAL(4)   :: gamrf, gamsf, gamwf
  REAL(4)   :: elo_s, dopple_nel
  REAL(8)   :: wave8, freq8, resid
  REAL(4)   :: asynth(kw)

  ! --- CLI / filename scratch -------------------------------------------
  CHARACTER(LEN=512) :: model_file, spec_file, linform_file, mol_file, model_base
  CHARACTER(LEN=64)  :: tmparg
  CHARACTER(256)     :: envval
  INTEGER            :: envlen, envstat, dotpos

  INTEGER   :: iedge_sv   ! SPECTRV continuum edge bracket (!= SYNTHE's iedge)

  ! --- Wall-clock timing ------------------------------------------------
  INTEGER(8) :: clock_start, clock_end, clock_rate

  CALL SYSTEM_CLOCK(clock_start, clock_rate)

  IFSYNTHE = 1

  ! --- Data directory (via $ATLAS12; defaults to ./data) ----------------
  CALL GET_ENVIRONMENT_VARIABLE('ATLAS12', envval, envlen, envstat)
  IF (envstat == 0 .AND. envlen > 0) THEN
    DATADIR = TRIM(envval)
    IF (DATADIR(envlen:envlen) /= '/') DATADIR = TRIM(DATADIR) // '/'
    DATADIR = TRIM(DATADIR) // 'data/'
  ELSE
    DATADIR = 'data/'
  END IF
  CALL set_bc_data_dir(TRIM(DATADIR))   ! point B&C partition-fn reader at the same dir

  ! --- Parse command-line arguments -------------------------------------
  !   Arg 1       : model atmosphere filename (positional, required)
  !   wlbeg=<nm>  : start wavelength (required)
  !   wlend=<nm>  : end wavelength   (required)
  !   resolu=<R>  : resolving power  (optional, default 300000)
  !   turbv=<kms> : extra microturbulence (optional, default 0.0)
  ! Output file basenames are model_file with extension stripped.
  ! -----------------------------------------------------------------------
  IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
    WRITE(6,'(A)') ' ERROR: expected model filename as first argument'
    WRITE(6,'(A)') USAGE
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, model_file)

  resolu = 300000.0D0
  turbv  = 0.0
  wlbeg  = 0.0D0    ! sentinel
  wlend  = 0.0D0    ! sentinel

  DO i = 2, COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(i, tmparg)
    eqpos = INDEX(tmparg, '=')
    IF (eqpos < 2) THEN
      WRITE(6,'(A,A)') ' ERROR: unrecognised argument (expected key=value): ', TRIM(tmparg)
      STOP 1
    END IF
    SELECT CASE (tmparg(1:eqpos-1))
      CASE ('wlbeg');  READ(tmparg(eqpos+1:), *) wlbeg
      CASE ('wlend');  READ(tmparg(eqpos+1:), *) wlend
      CASE ('resolu'); READ(tmparg(eqpos+1:), *) resolu
      CASE ('turbv');  READ(tmparg(eqpos+1:), *) turbv
      CASE DEFAULT
        WRITE(6,'(A,A)') ' ERROR: unknown keyword argument: ', TRIM(tmparg)
        STOP 1
    END SELECT
  END DO

  IF (wlbeg <= 0.0D0 .OR. wlend <= wlbeg) THEN
    WRITE(6,'(A)') ' ERROR: require 0 < wlbeg < wlend'
    WRITE(6,'(A)') USAGE
    STOP 1
  END IF

  WRITE(6,'(A,A)')     ' Input model      = ', TRIM(model_file)
  WRITE(6,'(A,F10.3)') ' wlbeg (nm)       = ', wlbeg
  WRITE(6,'(A,F10.3)') ' wlend (nm)       = ', wlend
  WRITE(6,'(A,F10.1)') ' resolu           = ', resolu
  WRITE(6,'(A,F8.4)')  ' turbv (km/s)     = ', turbv

  ! Derive output filenames: strip leading path, strip last extension
  dotpos = INDEX(TRIM(model_file), '/', BACK=.TRUE.)
  IF (dotpos > 0) THEN
    model_base = model_file(dotpos+1:)
  ELSE
    model_base = TRIM(model_file)
  END IF
  dotpos = INDEX(TRIM(model_base), '.', BACK=.TRUE.)
  IF (dotpos > 1) model_base = model_base(1:dotpos-1)
  spec_file    = TRIM(model_base) // '.spec'
  mol_file     = TRIM(model_base) // '.mol'
  linform_file = TRIM(model_base) // '.linform'
  WRITE(6,'(A,A)') ' Spectrum output  = ', TRIM(spec_file)

  ! --- Wavelength grid ---------------------------------------------------
  !   length: number of log-lambda steps from wlbeg to wlend at resolu
  !   (~636080 for 300-2500 nm at R=300000).
  !   wbegin is snapped upward to the first grid point >= wlbeg.
  ratio   = 1.0D0 + 1.0D0 / resolu
  ratiolg = LOG(ratio)
  length  = INT(LOG(wlend / wlbeg) / ratiolg)
  ixwlbeg = INT(LOG(wlbeg) / ratiolg)
  wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  IF (wbegin < wlbeg) THEN
    ixwlbeg = ixwlbeg + 1
    wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  END IF

  ! --- READIN + run_xnfpelsyn ------------------------------------------
  ! run_mklinelist runs after these so TEFF is set (gates cool-star
  ! line lists: H2O, TiO).  The two stages are otherwise independent.
  itemp_a = 1
  OPEN(UNIT=5,  FILE=TRIM(model_file),              STATUS='OLD', ACTION='READ')
  OPEN(UNIT=17, FILE=TRIM(DATADIR)//'continua.dat', STATUS='OLD', ACTION='READ')
  CALL readin(20)
  ! Keep unit 5 open: MOLEC reads from INPUTDATA(=5) on first call below.
  CALL run_xnfpelsyn()
  CLOSE(UNIT=5)
  CLOSE(UNIT=17)

  ! --- Build line lists in memory --------------------------------------
  ! Populates lte_lines(:), nlte_lines(:), nlines_lte, nlines_nlte in
  ! mod_mklinelist.  TEFF > TEFF_COOL_LIMIT skips H2O/TiO lists.
  CALL run_mklinelist(wlbeg, wlend, resolu, TEFF, &
                      TRIM(DATADIR) // 'lines.list', DATADIR)

  WRITE(6,'(A,I9,A,I9,A,I9)') ' Lines:  LTE =', nlines_lte, &
       '   NLTE =', nlines_nlte, '   total =', nlines_lte + nlines_nlte

  OPEN(UNIT=35, FILE=TRIM(mol_file), STATUS='REPLACE', ACTION='WRITE')

  ! --- Copy module state into local working arrays ---------------------
  nrhox = nrhox_a

  nedge = nedge_m
  DO iedge = 1, nedge_m
    wledge(iedge) = ABS(wledge_m(iedge))
  END DO

  DO iedge = 2, nedge
    halfedge(iedge-1) = (wledge(iedge-1) + wledge(iedge)) * 0.5D0
    deledge(iedge-1)  = (wledge(iedge) - wledge(iedge-1))**2 * 0.5D0
  END DO

  ! --- Depth structure ---
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

  ! --- SPECTRV depth-dependent work arrays (see fudge-parameter block) --
  DO j = 1, nrhox
    bfudge_sv(j) = bhyd_gs(j)**PH1 * (bc1_gs(j)/bc2_gs(j))**PC1 * &
                   (bsi1_gs(j)/bsi2_gs(j))**PSI1
    fscat_sv(j)  = 0.0D0
    IF (LINE_SCAT_RHOX_SCALE /= 0.0D0) &
      fscat_sv(j) = EXP(-rhox_a(j) / LINE_SCAT_RHOX_SCALE)
  END DO

  itemp_a = 1

  ! NB: CLI turbv is NOT folded into vturb_a here; it is added per depth
  ! inside the depth loop at qdopple = SQRT(qdopple^2 + (turbv/c)^2).

  DELTAW = resolu
  NULO   = 1
  NUHI   = length
  NUMNU  = length

  iedge_sv = 1

  ! JOSH operator matrices (COEFJ, COEFH).  F90 BLOCKJ/BLOCKH read + cache
  ! from file; call once here instead of the F77 per-call BLOCKJH path.
  CALL BLOCKJ
  CALL BLOCKH

  ! --- Main depth loop: SYNTHE opacity accumulation ---------------------
  ALLOCATE(opacity_matrix(length, nrhox))
  opacity_matrix = 0.0

  depth_loop: DO j = 1, nrhox

    buffer(1:length) = 0.0

    ! Unpack total continuum opacity (log10) into (3, iedge) quadratic
    ! interpolation basis over wavelength, then exponentiate.
    DO nu = 1, numnu_m
      qablog(nu) = continall_m(nu, j)
    END DO
    nu = 0
    DO iedge = 1, nedge-1
      nu = nu + 1;  ablog(1,iedge) = qablog(nu)
      nu = nu + 1;  ablog(2,iedge) = qablog(nu)
      nu = nu + 1;  ablog(3,iedge) = qablog(nu)
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

    ! Unpack ion populations and Doppler widths from (6, mw, kw) into flat mw6.
    DO nelem_i = 1, mw
      DO i = 1, 6
        qxnfpel((nelem_i-1)*6 + i) = xnfpel_m(i, nelem_i, j)
        qdopple((nelem_i-1)*6 + i) = dopple_m (i, nelem_i, j)
      END DO
    END DO

    xnfph (j,1) = REAL(qxnfpel(1))
    xnfph (j,2) = REAL(qxnfpel(2))
    xnfphe(j,1) = REAL(qxnfpel(7))
    xnfphe(j,2) = REAL(qxnfpel(8))
    xnfphe(j,3) = REAL(qxnfpel(9))

    DO i = 1, mw6
      xnfpel(i) = 0.0D0
      IF (qxnfpel(i) < 1.0D25) xnfpel(i) = qxnfpel(i)
    END DO

    DO i = 1, mw6
      qdopple(i) = SQRT(qdopple(i)**2 + (DBLE(turbv)/CLIGHT_KMS)**2)
      dopple(i)  = qdopple(i)
      xnfpel(i)  = xnfpel(i) / rho(j)
      IF (qdopple(i) > 0.0D0) THEN
        xnfdop(i) = qxnfpel(i) / xf_rho(j) / qdopple(i)
      ELSE
        xnfdop(i) = 0.0D0
      END IF
    END DO

    txnxn(j) = REAL( (xnfh(j) + 0.42D0*xnfhe(j,1) + 0.85D0*xnfh2(j)) * &
                     (t(j)/10000.0D0)**0.3D0 )

    ! NLTE / complex-profile lines (velshift argument reserved for future use)
    IF (nlines_nlte > 0) CALL compute_line_opacity(j, nlines_nlte, CUTOFF, 0.0, IFVAC)

    ! LTE metal lines: Voigt core on the wavelength grid + r^-2 far-wing tail
    IF (nlines_lte > 0) THEN
       DO iline = 1, nlines_lte
          nbuff_s   = lte_lines(iline)%nbuff
          congf_s   = lte_lines(iline)%cgf
          congf_nel = lte_lines(iline)%nelion
          elo_s     = lte_lines(iline)%elo
          gamrf     = lte_lines(iline)%gamrf
          gamsf     = lte_lines(iline)%gamsf
          gamwf     = lte_lines(iline)%gamwf

          kappa0_s = congf_s * REAL(xnfdop(congf_nel))
          kapmin_s = continuum(MIN(MAX(nbuff_s,1),length)) * CUTOFF
          IF (kappa0_s < kapmin_s) CYCLE
          kappa0_s = kappa0_s * REAL(EXP(-elo_s * hckt(j)))
          IF (kappa0_s < kapmin_s) CYCLE

          adamp_s    = REAL((gamrf + gamsf*xne(j) + gamwf*txnxn(j)) / dopple(congf_nel))
          n10dop     = INT(10.0D0 * dopple(congf_nel) * DBLE(resolu))
          dopple_nel = REAL(dopple(congf_nel))

          centre_on_grid: IF (nbuff_s >= 1 .AND. nbuff_s <= length) THEN
             kapcen_s = kappa0_s * voigt_profile(0.0, adamp_s)
             buffer(nbuff_s) = buffer(nbuff_s) + kapcen_s
          END IF centre_on_grid

          ! Voigt profile out to 10 Doppler widths; stop when below cutoff.
          dvoigt = 1.0 / dopple_nel / REAL(resolu)
          DO nstep = 1, n10dop
             profile(nstep) = kappa0_s * voigt_profile(REAL(nstep)*dvoigt, adamp_s)
             IF (profile(nstep) < kapmin_s) EXIT
          END DO

          ! Far wing: extend as x_wing / nstep^2 (Lorentzian asymptote).
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

          ! Red wing (plus +0 lines that fell just off-grid at centre).
          IF (nbuff_s < length) THEN
             maxred    = MIN(length - nbuff_s, nstep)
             minblue_i = MAX(1, 1 - nbuff_s)
             DO istep = minblue_i, maxred
                buffer(nbuff_s + istep) = buffer(nbuff_s + istep) + profile(istep)
             END DO
             IF (nbuff_s <= 1) CYCLE
          END IF

          ! Blue wing.
          maxblue   = MIN(nbuff_s - 1, nstep)
          minblue_i = MAX(1, nbuff_s - length)
          DO istep = minblue_i, maxblue
             buffer(nbuff_s - istep) = buffer(nbuff_s - istep) + profile(istep)
          END DO
       END DO
    END IF

    opacity_matrix(1:length, j) = buffer(1:length)

  END DO depth_loop

  ! --- SPECTRV wavelength loop: radiative transfer per wavelength point -
  OPEN(UNIT=11, FILE=TRIM(spec_file),    STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=33, FILE=TRIM(linform_file), STATUS='REPLACE', ACTION='WRITE')

  DO iwave = 1, length
    wave8 = wbegin * ratio**(iwave-1)
    freq8 = CLIGHT_NM_HZ / wave8
    DO j = 1, nrhox
      asynth(j) = REAL( opacity_matrix(iwave, j) * (1.0D0 - EXP(-freq8*hkt(j))) )
    END DO
    CALL process_wavelength_point(iwave, wave8, asynth)
  END DO

  DEALLOCATE(opacity_matrix)
  CLOSE(UNIT=11)
  CLOSE(UNIT=33)
  CLOSE(UNIT=35)

  CALL SYSTEM_CLOCK(clock_end)
  CALL report_elapsed(REAL(clock_end - clock_start, 8) / REAL(clock_rate, 8))

  STOP

CONTAINS

  ! ------------------------------------------------------------------------
  !  setup_opacity_sv(wave)
  !
  !  Per-wavelength continuum + frequency setup for SPECTRV.  Advances
  !  iedge_sv (separate from the depth-loop iedge) to bracket wave,
  !  interpolates ACONT/SIGMAC from the log-space tables via quadratic
  !  Lagrange weights, and sets FREQ, FREQLG, EHVKT, STIM, BNU, plus
  !  initialised line slots (ALINE=SIGMAL=0, SLINE=SCONT=BNU).
  ! ------------------------------------------------------------------------
  SUBROUTINE setup_opacity_sv(wave)
    REAL(8), INTENT(IN) :: wave

    REAL(8) :: freq15, c1, c2, c3
    INTEGER :: jj

    DO WHILE (wave >= wledge(iedge_sv + 1))
      iedge_sv = iedge_sv + 1
    END DO

    c1 = (wave - halfedge(iedge_sv)) * (wave - wledge(iedge_sv+1))         / deledge(iedge_sv)
    c2 = (wledge(iedge_sv) - wave)   * (wave - wledge(iedge_sv+1)) * 2.0D0 / deledge(iedge_sv)
    c3 = (wave - wledge(iedge_sv))   * (wave - halfedge(iedge_sv))         / deledge(iedge_sv)

    DO jj = 1, nrhox
      ACONT(jj)  = 10.0D0**(c1*contabs_sv (1,iedge_sv,jj) + &
                            c2*contabs_sv (2,iedge_sv,jj) + &
                            c3*contabs_sv (3,iedge_sv,jj))
      SIGMAC(jj) = 10.0D0**(c1*contscat_sv(1,iedge_sv,jj) + &
                            c2*contscat_sv(2,iedge_sv,jj) + &
                            c3*contscat_sv(3,iedge_sv,jj))
    END DO

    FREQ   = CLIGHT_NM_HZ / wave
    freq15 = FREQ / 1.0D15
    FREQLG = LOG(FREQ)
    DO jj = 1, nrhox
      EHVKT(jj)  = EXP(-FREQ * hkt_a(jj))
      STIM(jj)   = 1.0D0 - EHVKT(jj)
      BNU(jj)    = PLANCK_PREFACTOR * freq15**3 * EHVKT(jj) / STIM(jj)
      ALINE(jj)  = 0.0D0
      SIGMAL(jj) = 0.0D0
      SLINE(jj)  = BNU(jj)
      SCONT(jj)  = BNU(jj)
    END DO

  END SUBROUTINE setup_opacity_sv


  ! ------------------------------------------------------------------------
  !  process_wavelength_point(nu_idx, wave_in, asyn)
  !
  !  Radiative transfer at one wavelength point.  Calls JOSH twice (pure
  !  continuum, then continuum+line) and writes one record each to the
  !  .spec (flux vs wavelength) and .linform (tau table) outputs.
  !
  !    nu_idx  : wavelength point index (1..LENGTH)
  !    wave_in : wavelength in nm
  !    asyn    : stimulated-emission-corrected line opacity vector (nrhox)
  ! ------------------------------------------------------------------------
  SUBROUTINE process_wavelength_point(nu_idx, wave_in, asyn)
    INTEGER, INTENT(IN) :: nu_idx
    REAL(8), INTENT(IN) :: wave_in
    REAL(4), INTENT(IN) :: asyn(kw)

    REAL(8) :: q_loc(41)
    INTEGER :: jj, mu_loc

    CALL setup_opacity_sv(wave_in)

    ! Pure-continuum JOSH -> surf_sv (used as baseline for each mu)
    CALL josh(IFSCAT, IFSURF)
    DO mu_loc = 1, NMU
      IF (IFSURF == 1) surf_sv(mu_loc) = HNU(1)
      IF (IFSURF == 2) surf_sv(mu_loc) = SURFI(mu_loc)
    END DO

    IF (nu_idx == 1) WRITE(6,'(A)') ' '

    ! Blend line opacity + fudge source function, then JOSH again.
    DO jj = 1, nrhox
      ALINE(jj)  = DBLE(asyn(jj)) * (1.0D0 - fscat_sv(jj))
      SLINE(jj)  = BNU(jj) * STIM(jj) / (bfudge_sv(jj) - EHVKT(jj))
      SIGMAL(jj) = DBLE(asyn(jj)) * fscat_sv(jj)
    END DO
    CALL josh(1, IFSURF)

    WRITE(33, '(F12.4,1P2E12.4,/(10E12.4))') wave_in, HNU(1), surf_sv(1), &
         (TAUNU(jj), jj=1,nrhox)

    DO mu_loc = 1, NMU
      IF (IFSURF == 1) resid = HNU(1)          / surf_sv(mu_loc)
      IF (IFSURF == 2) resid = SURFI(mu_loc)   / surf_sv(mu_loc)
      q_loc(mu_loc)       = resid * surf_sv(mu_loc)
      q_loc(mu_loc + NMU) = surf_sv(mu_loc)
    END DO

    ! .spec record: wavelength (A), flux, continuum flux
    WRITE(11, '(F11.4,2E15.6)') wave_in * 10.0D0, q_loc(1), q_loc(NMU + 1)

  END SUBROUTINE process_wavelength_point


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

    IF (seconds < 60.0D0) THEN
      WRITE(6,'(/A,F6.2,A)') ' Elapsed: ', seconds, 's'
    ELSE IF (seconds < 3600.0D0) THEN
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


END PROGRAM SYNTHE
