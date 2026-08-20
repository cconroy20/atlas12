! ============================================================================
!  PROGRAM SYNTHE
!
!  Merged XNFPELSYN + SYNTHE + SPECTRV spectral-synthesis driver.  All
!  inter-stage I/O is in-memory via module arrays; see synthe_module and
!  mod_mklinelist.  ATLAS library routines are accessed through
!  mod_atlas_data (atlas12_modules.f90).
!
!  Example usage:
!
!    ./synthe.exe model.dat wlbeg=400 wlend=700
!
!    ./synthe.exe model.dat wlbeg=400 wlend=700 \
!                         resolu=500000 turbv=1.5 more_output=yes
!
!  Arguments:
!    <model_file>       (required, positional) model atmosphere file
!    wlbeg=<nm>         (required) start wavelength in nm
!    wlend=<nm>         (required) end wavelength in nm
!    resolu=<R>         (optional) resolving power, default 300000.
!                       Values below 300000 are REFUSED: the grid must
!                       resolve the ~0.7 km/s molecular thermal width.
!    turbv=<kms>        (optional) microturbulence in km/s.  Any value >= 0
!                       REPLACES the model atmosphere's microturbulence at all
!                       depths -- INCLUDING turbv=0, which forces a genuine
!                       zero.  If omitted, the per-layer value from the input
!                       model is used.  (Before 2026-08-10 the test was > 0,
!                       so turbv=0 silently did nothing.)
!    more_output=<v>    (optional) if yes/true/1/on, also emit the
!                       <model>.mol and <model>.linform files.  Default no.
!
!  Units:
!     5  input  : model atmosphere cards (CLI argument, read by READIN)
!    11  output : ASCII spectrum <model>.spec
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
  USE mod_mklinelist, only: run_mklinelist, lte_lines, nlte_lines, nlines_lte, &
                            nlines_nlte, ICA4227
  USE mod_parameters, only: LINE_CUTOFF, &
                            CA4227_MODE, CA4227_PP, CA4227_PS, CA4227_PX, &
                            CA4227_DLMAX, CA4227_WLA
  USE mod_nlte,       only: nlte_on, nlte_init, nlte_report, nlte_dump, &
                            nlte_scan_reset, nlte_tag_at, nlte_factors, &
                            nlte_set_params
  USE mod_atlas_data, only: &
    JOSH, READIN, BLOCKJ, BLOCKH, set_bc_data_dir, &
    DATADIR, IFSYNTHE, IFMOLOUT, DUMP_CONTINUUM, TURBV_UNSET, &
    ! Renamed to avoid collision with module-level names in synthe_module
    hkt_a => HKT, itemp_a => ITEMP, &
    nrhox_a => NRHOX, rhox_a => RHOX, &
    ! Direct imports
    ACONT, ALINE, BNU, DELTAW, EHVKT, FREQ, FREQLG, &
    HNU, IFSCAT, IFSURF, NMU, NUHI, NULO, NUMNU, &
    SCONT, SIGMAC, SIGMAL, SLINE, STIM, &
    SURFI, TAUNU, TEFF, GLOG, ABUND, XRELATIVE
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
    ' Usage: synthe_spectrv.exe <model_file> wlbeg=<nm> wlend=<nm> ' // &
    '[resolu=<R>] [turbv=<kms>] [more_output=yes|no]'

  ! --- Fixed run parameters (formerly read from fort.93 / SYNBEG) -------
  !   ifvac = 1 : vacuum wavelengths throughout
  ! Line-rejection threshold uses LINE_CUTOFF from mod_parameters; that
  ! single constant keeps SYNTHE in lockstep with ATLAS's SELECTLINES /
  ! TABCONT cutoff (see KAPCONT header).
  INTEGER, PARAMETER :: IFVAC  = 1

  ! --- Run parameters set from CLI --------------------------------------
  REAL(4)   :: turbv
  LOGICAL   :: more_output   ! print .mol (molecular densities) and .linform (tau tables)

  ! --- Model atmosphere header ------------------------------------------
  INTEGER   :: nedge
  REAL(8)   :: wledge(me), deledge(me), halfedge(me)

  ! --- Per-depth atmosphere read buffers --------------------------------
  REAL(8)   :: qdopple(mw6), qxnfpel(mw6)
  REAL(8)   :: qablog(nf)
  REAL(8)   :: ablog(3, me)

  ! --- In-memory opacity matrix, (length, nrhox) ------------------------
  REAL(4), ALLOCATABLE :: opacity_matrix(:,:)

  ! --- NLTE emissivity deviation, (length, nrhox) -----------------------
  ! Companion to opacity_matrix holding  sum_i kappa_i * (r_i - 1)  where
  ! r_i = S_i/B_nu, so that the opacity-weighted line source function is
  !   SLINE = BNU * (opacity_matrix + dev_matrix) / opacity_matrix.
  ! Only lines with departure coefficients contribute, so this stays
  ! identically zero away from them and the LTE line loop never writes to
  ! it.  Allocated ONLY when nlte_on -- it is the same size as
  ! opacity_matrix (about 100 MB for a full optical window at R = 3e5),
  ! and there is no reason to pay that in production.  dev_buffer is its
  ! per-depth accumulator, the twin of BUFFER.
  REAL(4), ALLOCATABLE :: dev_matrix(:,:)
  REAL(4), ALLOCATABLE :: dev_buffer(:)
  REAL(8)              :: sratio_sv(kw)   ! SLINE/BNU per depth, this lambda
  REAL(8)              :: fkappa_s, fdev_s, xhnukt
  REAL(4)              :: devfac_s
  INTEGER              :: k_nlte

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
  LOGICAL   :: is_ca4227
  INTEGER   :: ca4227_nmax
  REAL(4)   :: xv_s, prof_s
  REAL(8)   :: hmod_s
  REAL(8)   :: wave8, freq8, resid
  REAL(4)   :: asynth(kw)

  ! --- CLI / filename scratch -------------------------------------------
  CHARACTER(LEN=512) :: model_file, spec_file, linform_file, mol_file, model_base
  CHARACTER(LEN=512) :: cont_file
  ! Minimum supported resolving power (see the check below).
  REAL(8), PARAMETER :: RESOLU_MIN = 300000.0D0
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
  IF (envstat .EQ. 0 .AND. envlen .GT. 0) THEN
    DATADIR = TRIM(envval)
    IF (DATADIR(envlen:envlen) .NE. '/') DATADIR = TRIM(DATADIR) // '/'
    DATADIR = TRIM(DATADIR) // 'data/'
  ELSE
    DATADIR = 'data/'
  END IF
  CALL set_bc_data_dir(TRIM(DATADIR))   ! point B&C partition-fn reader at the same dir

  ! --- Parse command-line arguments -------------------------------------
  !   Arg 1            : model atmosphere filename (positional, required)
  !   wlbeg=<nm>       : start wavelength (required)
  !   wlend=<nm>       : end wavelength   (required)
  !   resolu=<R>       : resolving power (optional, default and MINIMUM
  !                      300000 -- lower values are refused, see below)
  !   turbv=<kms>      : microturbulence; any >=0 value replaces the model's
  !                      at all depths (0 forces zero); omitted keeps the model
  !   more_output=<v>  : if yes/true/1, also emit <model>.mol (molecular
  !                      number densities) and <model>.linform (tau table).
  !                      Accepted true values: yes, true, 1, on.
  !                      Accepted false values: no, false, 0, off.
  !                      Default: no.
  ! Output file basenames are model_file with extension stripped.
  ! -----------------------------------------------------------------------
  IF (COMMAND_ARGUMENT_COUNT() .LT. 1) THEN
    WRITE(6,'(A)') ' ERROR: expected model filename as first argument'
    WRITE(6,'(A)') USAGE
    CALL EXIT(1)
  END IF
  CALL GET_COMMAND_ARGUMENT(1, model_file)

  resolu      = 300000.0D0
  turbv       = TURBV_UNSET   ! negative sentinel: "not given on the CLI"
  wlbeg       = 0.0D0    ! sentinel
  wlend       = 0.0D0    ! sentinel
  more_output = .FALSE.

  DO i = 2, COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(i, tmparg)
    eqpos = INDEX(tmparg, '=')
    IF (eqpos .LT. 2) THEN
      WRITE(6,'(A,A)') ' ERROR: unrecognised argument (expected key=value): ', TRIM(tmparg)
      CALL EXIT(1)
    END IF
    SELECT CASE (tmparg(1:eqpos-1))
      CASE ('wlbeg');  READ(tmparg(eqpos+1:), *) wlbeg
      CASE ('wlend');  READ(tmparg(eqpos+1:), *) wlend
      CASE ('resolu'); READ(tmparg(eqpos+1:), *) resolu
      CASE ('turbv');  READ(tmparg(eqpos+1:), *) turbv
      CASE ('more_output')
        SELECT CASE (TRIM(ADJUSTL(tmparg(eqpos+1:))))
          CASE ('yes', 'YES', 'true', 'TRUE', '1', 'on', 'ON', 'y', 'Y')
            more_output = .TRUE.
          CASE ('no', 'NO', 'false', 'FALSE', '0', 'off', 'OFF', 'n', 'N')
            more_output = .FALSE.
          CASE DEFAULT
            WRITE(6,'(A,A)') ' ERROR: more_output expects yes/no, got: ', &
              TRIM(tmparg(eqpos+1:))
            CALL EXIT(1)
        END SELECT
      CASE DEFAULT
        WRITE(6,'(A,A)') ' ERROR: unknown keyword argument: ', TRIM(tmparg)
        CALL EXIT(1)
    END SELECT
  END DO

  IF (wlbeg .LE. 0.0D0 .OR. wlend .LE. wlbeg) THEN
    WRITE(6,'(A)') ' ERROR: require 0 < wlbeg < wlend'
    WRITE(6,'(A)') USAGE
    CALL EXIT(1)
  END IF

  ! Hard floor on the synthesis resolution.  The computation grid must
  ! resolve the thermal width of the molecular haze, which is ~0.7 km/s in
  ! a cool photosphere.  A convergence ladder against R = 500,000 (ratios
  ! after smoothing to the observed resolution) gives 0.898 at R = 50,000,
  ! 0.986 at 100,000 and 0.9987 at 300,000: below the floor the smoothed
  ! optical spectrum is systematically OVER-absorbed, by ~10% in the median
  ! and up to 17% in the strong TiO bands, because a 6 km/s grid undersamples
  ! lines it then integrates.  That is a silent bias, not noise, and it is
  ! wrong in a direction that looks like real absorption -- it produced a
  ! long-standing spurious offset before it was diagnosed.  Refuse rather
  ! than let anyone reintroduce it.
  IF (resolu .LT. RESOLU_MIN) THEN
    WRITE(6,'(A,F10.1)') ' ERROR: resolu below the supported floor: ', resolu
    WRITE(6,'(A,F10.1)') '        minimum is ', RESOLU_MIN
    WRITE(6,'(A)') '        Coarser grids undersample the molecular haze and'
    WRITE(6,'(A)') '        over-absorb the smoothed spectrum by ~10% (17% in'
    WRITE(6,'(A)') '        TiO bands).  See CHANGELOG, synthesis-resolution'
    WRITE(6,'(A)') '        convergence.  Raise resolu; do not lower the floor.'
    CALL EXIT(1)
  END IF

  ! Propagate more_output to the ATLAS-level flag read by NMOLEC
  IFMOLOUT = MERGE(1, 0, more_output)

  WRITE(6,*) ''
  WRITE(6,'(A,A)')     ' Input model      = ', TRIM(model_file)
  WRITE(6,'(A,F10.3)') ' wlbeg (nm)       = ', wlbeg
  WRITE(6,'(A,F10.3)') ' wlend (nm)       = ', wlend
  WRITE(6,'(A,F10.1)') ' resolu           = ', resolu
  IF (CA4227_MODE .LT. 0) THEN
    WRITE(6,'(A)')     ' Ca I 4227        = SUPPRESSED (CA4227_MODE < 0)'
  ELSE IF (CA4227_MODE .GT. 0) THEN
    WRITE(6,'(A,3(1X,1PE9.2),A,0PF7.1)') &
      ' Ca I 4227        = modified profile, (pp,ps,px) =', &
      CA4227_PP, CA4227_PS, CA4227_PX, ', |dlam| <', CA4227_DLMAX
  END IF
  ! Report the OVERRIDE, not a bare number that reads like the adopted
  ! microturbulence.  The value actually used is echoed after the model is
  ! read (it lives per-layer in the atmosphere until then).
  IF (turbv .GE. 0.0) THEN
    WRITE(6,'(A,F8.4)')  ' turbv override   = ', turbv
  ELSE
    WRITE(6,'(A)')       ' turbv override   = none (model VTURB)'
  END IF

  ! Derive output filenames: strip leading path, strip last extension
  dotpos = INDEX(TRIM(model_file), '/', BACK=.TRUE.)
  IF (dotpos .GT. 0) THEN
    model_base = model_file(dotpos+1:)
  ELSE
    model_base = TRIM(model_file)
  END IF
  dotpos = INDEX(TRIM(model_base), '.', BACK=.TRUE.)
  IF (dotpos .GT. 1) model_base = model_base(1:dotpos-1)
  spec_file    = TRIM(model_base) // '.spec'
  mol_file     = TRIM(model_base) // '.mol'
  linform_file = TRIM(model_base) // '.linform'
  cont_file    = TRIM(model_base) // '.cont'
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
  IF (wbegin .LT. wlbeg) THEN
    ixwlbeg = ixwlbeg + 1
    wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  END IF

  ! --- READIN + run_xnfpelsyn ------------------------------------------
  ! run_mklinelist runs after these so TEFF is set (gates cool-star
  ! line lists: H2O, TiO).  The two stages are otherwise independent.
  itemp_a = 1
  OPEN(UNIT=5,  FILE=TRIM(model_file), STATUS='OLD',     ACTION='READ')
  IF (more_output) &
    OPEN(UNIT=35, FILE=TRIM(mol_file), STATUS='REPLACE', ACTION='WRITE')
  ! Per-source continuum opacity decomposition, filled inside
  ! run_xnfpelsyn's KAPP loop.  All opacities cm^2/g.  Gated by the
  ! DUMP_CONTINUUM developer flag in mod_atlas_data, not by more_output.
  IF (DUMP_CONTINUUM) THEN
    OPEN(UNIT=36, FILE=TRIM(cont_file), STATUS='REPLACE', ACTION='WRITE')
    WRITE(36,'(A)') '# nu wave_nm j T rho ' // &
      'aHyd aH2plus aHmin aH2min aH2coll aHe1 aHe2 aHemin aMetal ACONT ' // &
      'sigH sigHe sigEl sigH2 sigX SIGMAC'
  END IF
  CALL readin(20)
  ! Keep unit 5 open: MOLEC reads from INPUTDATA(=5) on first call below.
  CALL run_xnfpelsyn(turbv)
  CLOSE(UNIT=5)

  ! --- Build line lists in memory --------------------------------------
  ! Populates lte_lines(:), nlte_lines(:), nlines_lte, nlines_nlte in
  ! mod_mklinelist.  TEFF > TEFF_COOL_LIMIT skips H2O/TiO lists.
  WRITE(6,*) ''
  WRITE(6,*) 'Assembling line list...'
  CALL run_mklinelist(wlbeg, wlend, resolu, TEFF, &
                      TRIM(DATADIR) // 'lines.list', DATADIR)

  WRITE(6,'(A,I9,A,I9,A,I9)') ' Lines:  LTE =', nlines_lte, &
       '   NLTE =', nlines_nlte, '   total =', nlines_lte + nlines_nlte

  ! A non-default Ca I 4227 mode that silently found no line to act on
  ! would look exactly like a null result.  Fail loudly instead.
  IF (CA4227_MODE .NE. 0) THEN
    IF (ICA4227 .EQ. 0) THEN
      WRITE(6,'(A)') ' ERROR: CA4227_MODE is set but the Ca I 4227 line is' // &
                     ' not in the synthesis window.'
      CALL EXIT(1)
    END IF
    WRITE(6,'(A,I9)') ' Ca I 4227 resonance line at LTE index', ICA4227
  END IF

  ! Optional: dump the assembled "used lines" to an ASCII file for external
  ! linelist validation (e.g., cross-checking the Korg master linelist).
  IF (more_output) CALL dump_used_lines()

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

  ! Model parameters for the NLTE_MODE = 3 grid lookup.  vturb is taken at
  ! the layer nearest tau_5000 = 1, i.e. where the lines form, rather than
  ! averaged over a scale on which it may vary.
  eqpos = 1
  DO j = 2, nrhox
    IF (ABS(LOG10(MAX(tau5000_sv(j),1.0D-99))) .LT. &
        ABS(LOG10(MAX(tau5000_sv(eqpos),1.0D-99)))) eqpos = j
  END DO
  CALL nlte_set_params(TEFF, GLOG, xf_vturb(eqpos)/1.0D5, &
                       ABUND(1), ABUND(11) + XRELATIVE(11), &
                       ABUND(26) + XRELATIVE(26))

  ! Departure coefficients for tagged transitions (NLTE_MODE developer
  ! flag).  Sets nlte_on; a no-op when NLTE_MODE = 0.  Placed here, after
  ! the depth structure is unpacked, because the b profile is a function
  ! of column mass -- for the structured test by construction, and for an
  ! external grid because that is the variable it must be interpolated on.
  CALL nlte_init(nrhox, xf_rhox(1:nrhox), xf_t(1:nrhox), &
                 tau5000_sv(1:nrhox))

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
    IF (LINE_SCAT_RHOX_SCALE .NE. 0.0D0) &
      fscat_sv(j) = EXP(-rhox_a(j) / LINE_SCAT_RHOX_SCALE)
  END DO

  itemp_a = 1

  ! NB: a CLI turbv >= 0 REPLACES the model atmosphere's microturbulence at
  ! all depths; this is applied inside run_xnfpelsyn (which substitutes turbv
  ! into VTURB before building the Doppler widths).  Omitting it (the negative
  ! sentinel) leaves the per-layer model value untouched.

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
  WRITE(6,*) ''
  WRITE(6,*) 'Accumulating opacity...'
  ALLOCATE(opacity_matrix(length, nrhox))
  opacity_matrix = 0.0
  IF (nlte_on) THEN
    ALLOCATE(dev_matrix(length, nrhox), dev_buffer(length))
    dev_matrix = 0.0
  END IF

  depth_loop: DO j = 1, nrhox

    buffer(1:length) = 0.0
    IF (nlte_on) dev_buffer(1:length) = 0.0

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
      DO WHILE (wave8 .GE. wledge(iedge+1) .AND. iedge .LT. nedge-1)
        iedge = iedge + 1
      END DO
      continuum(nbuff_i) = REAL( &
        ((wave8 - halfedge(iedge))*(wave8 - wledge(iedge+1))*ablog(1,iedge) + &
         (wledge(iedge) - wave8)*(wave8 - wledge(iedge+1))*2.0D0*ablog(2,iedge) + &
         (wave8 - wledge(iedge))*(wave8 - halfedge(iedge))*ablog(3,iedge)) / &
        deledge(iedge) )
      ! Defensive floor on the LOG10 interpolation: a quadratic over
      ! tightly-spaced edges can extrapolate the log-opacity to
      ! arbitrarily large negative values within an interval, which
      ! would crater kapmin_s in the line loop below.  Floor at
      ! log10 = -15 is far below any physically reasonable continuum
      ! opacity (~1e-10 cm^2/g at worst) yet large enough to keep
      ! kapmin_s well above the integer-overflow regime in the
      ! Voigt-wing extent calculation maxstep = SQRT(x_wing/kapmin_s).
      IF (continuum(nbuff_i) .LT. -15.0) continuum(nbuff_i) = -15.0
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
      IF (qxnfpel(i) .LT. 1.0D25) xnfpel(i) = qxnfpel(i)
    END DO

    ! qdopple already carries the effective microturbulence: run_xnfpelsyn
    ! substituted the CLI turbv into VTURB when it was positive, so both the
    ! thermal and turbulent terms are baked into the Doppler widths here.
    DO i = 1, mw6
      dopple(i)  = qdopple(i)
      xnfpel(i)  = xnfpel(i) / rho(j)
      IF (qdopple(i) .GT. 0.0D0) THEN
        xnfdop(i) = qxnfpel(i) / xf_rho(j) / qdopple(i)
      ELSE
        xnfdop(i) = 0.0D0
      END IF
    END DO

    txnxn(j) = REAL( (xnfh(j) + 0.42D0*xnfhe(j,1) + 0.85D0*xnfh2(j)) * &
                     (t(j)/10000.0D0)**0.3D0 )

    ! NLTE / complex-profile lines (velshift argument reserved for future use)
    IF (nlines_nlte .GT. 0) CALL compute_line_opacity(j, nlines_nlte, LINE_CUTOFF, 0.0, IFVAC)

    ! LTE metal lines: Voigt core on the wavelength grid + r^-2 far-wing tail
    IF (nlines_lte .GT. 0) THEN
       CALL nlte_scan_reset()
       DO iline = 1, nlines_lte
          ! NLTE membership.  Must be tested before any CYCLE below, since
          ! the cursor in nlte_tag_at only advances on a hit and expects to
          ! see every iline in order.  k_nlte = 0 for the overwhelming
          ! majority of lines; when NLTE is off it costs one compare.
          k_nlte = 0
          IF (nlte_on) k_nlte = nlte_tag_at(iline)

          nbuff_s   = lte_lines(iline)%nbuff
          congf_s   = lte_lines(iline)%cgf
          congf_nel = lte_lines(iline)%nelion
          elo_s     = lte_lines(iline)%elo
          gamrf     = lte_lines(iline)%gamrf
          gamsf     = lte_lines(iline)%gamsf
          gamwf     = lte_lines(iline)%gamwf

          ! Ca I 4227: developer switch for the resonance-line profile.
          ! ICA4227 is 0 unless that exact line is in the window.
          is_ca4227 = (CA4227_MODE .NE. 0) .AND. (iline .EQ. ICA4227)
          IF (is_ca4227 .AND. CA4227_MODE .LT. 0) CYCLE

          kappa0_s = congf_s * REAL(xnfdop(congf_nel))

          ! Departure coefficients: rescale the opacity and record the
          ! per-unit-opacity emissivity deviation.  Applied BEFORE the
          ! cutoff tests so a line weakened out of relevance by NLTE is
          ! dropped on its true strength.  devfac_s is 0 in LTE, so the
          ! deviation accumulator stays untouched for every other line.
          devfac_s = 0.0
          IF (k_nlte .GT. 0) THEN
             xhnukt = CLIGHT_NMS / (wbegin * ratio**(nbuff_s-1)) * hkt(j)
             CALL nlte_factors(k_nlte, j, xhnukt, fkappa_s, fdev_s)
             kappa0_s = kappa0_s * REAL(fkappa_s, 4)
             devfac_s = REAL(fdev_s, 4)
          END IF

          kapmin_s = continuum(MIN(MAX(nbuff_s,1),length)) * REAL(LINE_CUTOFF, 4)
          IF (kappa0_s .LT. kapmin_s) CYCLE
          kappa0_s = kappa0_s * REAL(EXP(-elo_s * hckt(j)))
          IF (kappa0_s .LT. kapmin_s) CYCLE

          adamp_s    = REAL((gamrf + gamsf*xne(j) + gamwf*txnxn(j)) / dopple(congf_nel))
          n10dop     = INT(10.0D0 * dopple(congf_nel) * DBLE(resolu))
          dopple_nel = REAL(dopple(congf_nel))

          centre_on_grid: IF (nbuff_s .GE. 1 .AND. nbuff_s .LE. length) THEN
             kapcen_s = kappa0_s * voigt_profile(0.0, adamp_s)
             buffer(nbuff_s) = buffer(nbuff_s) + kapcen_s
             IF (k_nlte .GT. 0) &
               dev_buffer(nbuff_s) = dev_buffer(nbuff_s) + kapcen_s*devfac_s
          END IF centre_on_grid

          dvoigt = 1.0 / dopple_nel / REAL(resolu)

          IF (is_ca4227) THEN
             ! Non-Voigt resonance profile: MAX(Voigt, Jones+23 form), run
             ! out to a declared detuning rather than to the wing cutoff
             ! (see CA4227_* in mod_parameters for why the cutoff cannot be
             ! trusted to terminate an x^-PX wing).  Both terms of the
             ! modified form decrease monotonically, so an early EXIT below
             ! kapmin_s is safe.
             ca4227_nmax = MIN(maxprof, &
                               INT(CA4227_DLMAX * resolu / CA4227_WLA))
             DO nstep = 1, ca4227_nmax
                xv_s   = REAL(nstep) * dvoigt
                hmod_s = EXP(-CA4227_PP * DBLE(xv_s)**2) + &
                         CA4227_PS * DBLE(adamp_s) / DBLE(xv_s)**CA4227_PX
                prof_s = kappa0_s * &
                         MAX(voigt_profile(xv_s, adamp_s), REAL(hmod_s,4))
                profile(nstep) = prof_s
                IF (prof_s .LT. kapmin_s) EXIT
             END DO
             nstep = MIN(nstep, ca4227_nmax)

          ELSE
             ! Voigt profile out to 10 Doppler widths; stop when below cutoff.
             DO nstep = 1, n10dop
                profile(nstep) = kappa0_s * voigt_profile(REAL(nstep)*dvoigt, adamp_s)
                IF (profile(nstep) .LT. kapmin_s) EXIT
             END DO

             ! Far wing: extend as x_wing / nstep^2 (Lorentzian asymptote).
             IF (nstep .GT. n10dop) THEN
                x_wing  = profile(n10dop) * REAL(n10dop)**2
                maxstep = INT(SQRT(x_wing/kapmin_s)) + 1
                maxstep = MIN(maxstep, maxprof)
                n1 = n10dop + 1
                DO nstep = n1, maxstep
                   profile(nstep) = x_wing / REAL(nstep)**2
                END DO
                nstep = maxstep
             END IF
          END IF

          IF (nbuff_s+nstep .LT. 1 .OR. nbuff_s-nstep .GT. length) CYCLE

          ! Red wing (plus +0 lines that fell just off-grid at centre).
          IF (nbuff_s .LT. length) THEN
             maxred    = MIN(length - nbuff_s, nstep)
             minblue_i = MAX(1, 1 - nbuff_s)
             DO istep = minblue_i, maxred
                buffer(nbuff_s + istep) = buffer(nbuff_s + istep) + profile(istep)
             END DO
             IF (k_nlte .GT. 0) THEN
                DO istep = minblue_i, maxred
                   dev_buffer(nbuff_s + istep) = dev_buffer(nbuff_s + istep) &
                                               + profile(istep)*devfac_s
                END DO
             END IF
             IF (nbuff_s .LE. 1) CYCLE
          END IF

          ! Blue wing.
          maxblue   = MIN(nbuff_s - 1, nstep)
          minblue_i = MAX(1, nbuff_s - length)
          DO istep = minblue_i, maxblue
             buffer(nbuff_s - istep) = buffer(nbuff_s - istep) + profile(istep)
          END DO
          IF (k_nlte .GT. 0) THEN
             DO istep = minblue_i, maxblue
                dev_buffer(nbuff_s - istep) = dev_buffer(nbuff_s - istep) &
                                            + profile(istep)*devfac_s
             END DO
          END IF
       END DO
    END IF

    opacity_matrix(1:length, j) = buffer(1:length)
    IF (nlte_on) dev_matrix(1:length, j) = dev_buffer(1:length)

  END DO depth_loop

  ! --- Wavelength loop: radiative transfer per wavelength point ---
  WRITE(6,*) 'Synthesizing spectrum...'

  OPEN(UNIT=11, FILE=TRIM(spec_file),    STATUS='REPLACE', ACTION='WRITE')
  IF (more_output) &
    OPEN(UNIT=33, FILE=TRIM(linform_file), STATUS='REPLACE', ACTION='WRITE')

  DO iwave = 1, length
    wave8 = wbegin * ratio**(iwave-1)
    freq8 = CLIGHT_NMS / wave8
    DO j = 1, nrhox
      asynth(j) = REAL( opacity_matrix(iwave, j) * (1.0D0 - EXP(-freq8*hkt(j))) )
    END DO
    ! Opacity-weighted line source function in units of B_nu.  The
    ! stimulated-emission factor is common to numerator and denominator and
    ! cancels, so this can be formed from the raw (pre-STIM) accumulators.
    ! Where no line has departure coefficients the deviation is identically
    ! zero and sratio is exactly 1.
    sratio_sv(1:nrhox) = 1.0D0
    IF (nlte_on) THEN
      DO j = 1, nrhox
        IF (opacity_matrix(iwave, j) .GT. 0.0) &
          sratio_sv(j) = 1.0D0 + DBLE(dev_matrix(iwave, j)) &
                               / DBLE(opacity_matrix(iwave, j))
      END DO
    END IF
    CALL process_wavelength_point(wave8, asynth, sratio_sv)
  END DO

  IF (nlte_on) DEALLOCATE(dev_matrix, dev_buffer)
  CALL nlte_report()
  CALL nlte_dump(TRIM(model_base))
  DEALLOCATE(opacity_matrix)
  CLOSE(UNIT=11)
  IF (more_output) THEN
    CLOSE(UNIT=33)
    CLOSE(UNIT=35)
  END IF

  CALL SYSTEM_CLOCK(clock_end)
  CALL report_elapsed(REAL(clock_end - clock_start, 8) / REAL(clock_rate, 8))

  CALL EXIT(0)

CONTAINS

  ! ------------------------------------------------------------------------
  !  dump_used_lines()
  !
  !  Emit the assembled line list (exactly what SYNTHE selected for this
  !  window) to an ASCII file <model>.lines, gated on more_output.  This is
  !  a general "what lines are in this synthesis" dump; it is not tied to any
  !  particular downstream code.
  !
  !  Columns: kind  wl_vac_ang  code  nelion  elo_cm  strength  g_rad  g_stark  g_vdw
  !
  !  Design notes:
  !   - Both the raw integer nelion AND the decoded Kurucz species code are
  !     written.  For LTE lines the code is the true species captured by the
  !     reader (Z.ion for atoms, molecule code for molecules, incl. the FeH
  !     156->126 remap); for NLTE lines (all atomic) it is decoded from
  !     nelion via nelem=(nelion-1)/6+1, ion=mod(nelion-1,6).
  !   - LTE line centers are reconstructed from the log-lambda grid index
  !     nbuff using SYNTHE's own convention  wave = wbegin*ratio**(nbuff-1)
  !     (see the opacity loop in synthe_module); NLTE lines carry wlvac.
  !   - Wavelengths are vacuum, written in Angstrom (nm*10).
  ! ------------------------------------------------------------------------
  SUBROUTINE dump_used_lines()
    INTEGER            :: iL, u, ne, io
    REAL(8)            :: wlang
    REAL(4)            :: code_nlte
    CHARACTER(LEN=512) :: lines_file

    lines_file = TRIM(model_base) // '.lines'
    OPEN(NEWUNIT=u, FILE=TRIM(lines_file), STATUS='REPLACE', ACTION='WRITE')
    WRITE(u,'(A)') '# kind  wl_vac_ang  code  nelion  elo_cm  strength  g_rad  g_stark  g_vdw'
    WRITE(u,'(A,I10,A,I10)') '# nlines_lte=', nlines_lte, &
                             '  nlines_nlte=', nlines_nlte
    DO iL = 1, nlines_lte
      wlang = wbegin * ratio**(lte_lines(iL)%nbuff - 1) * 10.0D0
      WRITE(u,'(A4,1X,F13.4,1X,F9.2,1X,I8,1X,ES14.6,1X,ES13.5,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
        'LTE ', wlang, lte_lines(iL)%code, lte_lines(iL)%nelion, lte_lines(iL)%elo, &
        lte_lines(iL)%cgf, lte_lines(iL)%gamrf, lte_lines(iL)%gamsf, lte_lines(iL)%gamwf
    END DO
    DO iL = 1, nlines_nlte
      wlang = nlte_lines(iL)%wlvac * 10.0D0
      ne = (nlte_lines(iL)%nelion - 1)/6 + 1
      io = MOD(nlte_lines(iL)%nelion - 1, 6)
      code_nlte = REAL(ne,4) + REAL(io,4)/100.0
      WRITE(u,'(A4,1X,F13.4,1X,F9.2,1X,I8,1X,ES14.6,1X,ES13.5,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
        'NLTE', wlang, code_nlte, nlte_lines(iL)%nelion, nlte_lines(iL)%elo, &
        nlte_lines(iL)%gf, nlte_lines(iL)%gammar, nlte_lines(iL)%gammas, nlte_lines(iL)%gammaw
    END DO
    CLOSE(u)
    WRITE(6,'(A,A)') ' Line list dump   = ', TRIM(lines_file)
  END SUBROUTINE dump_used_lines

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

    DO WHILE (wave .GE. wledge(iedge_sv + 1))
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

    FREQ   = CLIGHT_NMS / wave
    freq15 = FREQ / 1.0D15
    FREQLG = LOG(FREQ)
    DO jj = 1, nrhox
      EHVKT(jj)  = EXP(-FREQ * hkt_a(jj))
      STIM(jj)   = 1.0D0 - EHVKT(jj)
      BNU(jj)    = BNU_PREFAC * freq15**3 * EHVKT(jj) / STIM(jj)
      ALINE(jj)  = 0.0D0
      SIGMAL(jj) = 0.0D0
      SLINE(jj)  = BNU(jj)
      SCONT(jj)  = BNU(jj)
    END DO

  END SUBROUTINE setup_opacity_sv


  ! ------------------------------------------------------------------------
  !  process_wavelength_point(wave_in, asyn)
  !
  !  Radiative transfer at one wavelength point.  Calls JOSH twice (pure
  !  continuum, then continuum+line) and writes one record each to the
  !  .spec (flux vs wavelength) and .linform (tau table) outputs.
  !
  !    wave_in : wavelength in nm
  !    asyn    : stimulated-emission-corrected line opacity vector (nrhox)
  !    sratio  : opacity-weighted line source function in units of B_nu
  !              (all 1 in LTE; see the NLTE block in the wavelength loop)
  ! ------------------------------------------------------------------------
  SUBROUTINE process_wavelength_point(wave_in, asyn, sratio)
    REAL(8), INTENT(IN) :: wave_in
    REAL(4), INTENT(IN) :: asyn(kw)
    REAL(8), INTENT(IN) :: sratio(kw)

    REAL(8) :: q_loc(41)
    INTEGER :: jj, mu_loc

    CALL setup_opacity_sv(wave_in)

    ! Pure-continuum JOSH -> surf_sv (used as baseline for each mu)
    CALL josh(IFSCAT, IFSURF)
    DO mu_loc = 1, NMU
      IF (IFSURF .EQ. 1) surf_sv(mu_loc) = HNU(1)
      IF (IFSURF .EQ. 2) surf_sv(mu_loc) = SURFI(mu_loc)
    END DO

    ! Blend line opacity + fudge source function, then JOSH again.
    DO jj = 1, nrhox
      ALINE(jj)  = DBLE(asyn(jj)) * (1.0D0 - fscat_sv(jj))
      SLINE(jj)  = BNU(jj) * STIM(jj) / (bfudge_sv(jj) - EHVKT(jj))
      SIGMAL(jj) = DBLE(asyn(jj)) * fscat_sv(jj)
    END DO

    ! NLTE source function, applied in a SEPARATE guarded loop rather than
    ! as a fourth factor in the expression above.  sratio is exactly 1 in
    ! LTE, so folding it in is mathematically a no-op -- but it is not a
    ! CODEGEN no-op: the extra operand changed how gfortran contracted the
    ! neighbouring multiply-add at -O3 -march=native, which moved 2 of 2544
    ! points of a 3000 K spectrum by one unit in the last printed digit.
    ! Harmless in itself, but it costs the exact mode-0 == pre-change
    ! equality that makes the null test able to prove anything.  Leaving
    ! the LTE expression textually untouched keeps that guarantee.
    IF (nlte_on) THEN
      DO jj = 1, nrhox
        SLINE(jj) = SLINE(jj) * sratio(jj)
      END DO
    END IF
    CALL josh(1, IFSURF)

    IF (more_output) &
      WRITE(33, '(F12.4,1P2E12.4,/(10E12.4))') wave_in, HNU(1), surf_sv(1), &
           (TAUNU(jj), jj=1,nrhox)

    DO mu_loc = 1, NMU
      IF (IFSURF .EQ. 1) resid = HNU(1)          / surf_sv(mu_loc)
      IF (IFSURF .EQ. 2) resid = SURFI(mu_loc)   / surf_sv(mu_loc)
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

    IF (seconds .LT. 60.0D0) THEN
      WRITE(6,'(/A,F6.2,A)') ' Runtime: ', seconds, 's'
    ELSE IF (seconds .LT. 3600.0D0) THEN
      m = INT(seconds / 60.0D0)
      s = NINT(seconds - 60.0D0 * m)
      WRITE(6,'(/A,I0,A,I2.2,A)') ' Runtime: ', m, 'm ', s, 's'
    ELSE
      h = INT(seconds / 3600.0D0)
      m = INT((seconds - 3600.0D0 * h) / 60.0D0)
      s = NINT(seconds - 3600.0D0 * h - 60.0D0 * m)
      WRITE(6,'(/A,I0,A,I2.2,A,I2.2,A)') ' Runtime: ', h, 'h ', m, 'm ', s, 's'
    END IF

  END SUBROUTINE report_elapsed


END PROGRAM SYNTHE
