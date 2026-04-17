! ============================================================================
!  PROGRAM SYNTHE
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
!    5   input : model atmosphere control cards (command-line arg, read by READIN)
!   11   output: ASCII spectrum (wavelength in Angstroms, flux, continuum flux)
!   17   input : continua.dat continuum edge frequency list (from DATADIR)
!   33   output: wavelength / flux / continuum / optical-depth table
!   93   input : run parameters from SYNBEG (unformatted, deleted after read)
!
!  Units REMOVED relative to the original three-program pipeline:
!   12   -- was: preprocessed LTE  line data (fort.12)  now: lte_lines()  in mod_mklinelist
!   19   -- was: preprocessed NLTE line data (fort.19)  now: nlte_lines() in mod_mklinelist
!   13   scratch: merged line archive           (deleted; was unit 13)
!   15   scratch: (ILINE, KAPCEN) pairs         (now: in-memory)
!
!  ASCII spectrum output format (<model>.spec):
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
!    7   -- was: binary spectrum for PLOTSY              (now: <model>.spec)
!   14*  -- was: direct-access opacity matrix scratch  (now: in-memory arrays)
!   25   -- was: spectrv.input run parameters         (now: hardcoded PARAMETERs)
!   28   -- was: XNFPEL diagnostic dump               (removed)
!   29   -- was: opacity spectrum diagnostic           (removed)
!   (syntoascanga post-processing step eliminated; ASCII output now inline)
! ============================================================================

PROGRAM SYNTHE
  
  USE synthe_module
  USE mod_mklinelist, only: run_mklinelist, lte_lines, nlines_lte, nlines_nlte
  USE mod_atlas_data, only: &
    ! Procedures
    JOSH, READIN, KAPP, BLOCKJ, BLOCKH, &
    ! Data directory path
    DATADIR, &
    ! flag to broadcast if synthe is running
    IFSYNTHE, & 
    ! Renamed imports (collide with synthe_module or local names)
    hkt_a => HKT, itemp_a => ITEMP, &
    nrhox_a => NRHOX, rhox_a => RHOX, &
    ahline_a => AHLINE, shline_a => SHLINE, &
    teff_a => TEFF, glog_a => GLOG, title_atlas => TITLE, &
    ! Direct imports (no name collision)
    ACONT, ALINE, AXLINE, BNU, DELTAW, EHVKT, FREQ, FREQLG, &
    HNU, IFPRES, IFSCAT, IFSURF, NMU, NUHI, NULO, NUMNU, &
    SCONT, SIGMAC, SIGMAL, SLINE, STIM, &
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
  INTEGER :: ifvac_loc, nmu2_loc

  ! ==========================================================================
  !  LOCAL PARAMETERS (aliases for synthe_module constants)
  ! ==========================================================================
  INTEGER, PARAMETER :: kw_p      = kw
  INTEGER, PARAMETER :: mw_p      = mw
  INTEGER, PARAMETER :: mw6_p     = mw6
  INTEGER, PARAMETER :: maxbuff_p = MAXBUFF
  INTEGER, PARAMETER :: maxlin_p  = MAXLIN

  ! ==========================================================================
  !  SOURCE-FUNCTION FUDGE PARAMETERS
  !
  !  These four numbers were historically read from the file 'spectrv.input'
  !  in the F77 SYNTHE/SPECTRV pipeline.  In every standard distribution
  !  the file shipped with all values set to zero, and they have been
  !  hardcoded here as PARAMETERs to eliminate the unnecessary file I/O.
  !  If you ever need to experiment with them, edit the values below and
  !  recompile.
  !
  !  ----------------------------------------------------------------------
  !  LINE_SCAT_RHOX_SCALE  -- column-mass scale (g/cm²) for the empirical
  !  line-scattering treatment that mimics NLTE source-function divergence
  !  in the upper photosphere.
  !
  !  Background.  In strict LTE, line opacity is treated as pure thermal
  !  absorption coupled to S = B_ν.  In real stars, the line source
  !  function in the upper atmosphere drops below B_ν because the line
  !  photons scatter (and J_ν < B_ν there) rather than thermalizing.
  !  Solving this rigorously requires multilevel NLTE statistical
  !  equilibrium; SYNTHE does not.  Instead, Kurucz provided a hack:
  !  smoothly reassign a fraction of the line opacity from the absorption
  !  channel into the coherent-scattering channel as one moves up in the
  !  atmosphere, with the fraction controlled by a single column-mass
  !  scale height.  Per depth point j, with column mass rhox(j):
  !
  !      f_scat(j) = exp( -rhox(j) / LINE_SCAT_RHOX_SCALE )
  !
  !  with f_scat=0 forced when LINE_SCAT_RHOX_SCALE = 0.  Then for every
  !  line at every wavelength:
  !
  !      ALINE(j)  = alpha_line(j) * (1 - f_scat(j))   ! absorption (S=B*fudge)
  !      SIGMAL(j) = alpha_line(j) *      f_scat(j)    ! scattering (S=J)
  !
  !  Deep in the photosphere rhox >> scale, f_scat → 0, all line opacity is
  !  thermal absorption (LTE).  High up rhox << scale, f_scat → 1, all
  !  line opacity becomes coherent scattering and the source function is
  !  driven toward J_ν, deepening the line cores in the way real
  !  scattering-dominated transitions do.  A typical "on" value would be
  !  ~0.1-1.0 g/cm² (roughly the column mass at the temperature minimum
  !  of a solar-type star), but the standard distribution ships with this
  !  set to zero, leaving SYNTHE in pure-LTE mode for line cores.  Modern
  !  workflows that actually need accurate cores of NLTE-sensitive species
  !  use precomputed departure-coefficient grids (or RH for chromospheric
  !  resonance lines) instead of relying on this fudge.
  !
  !  ----------------------------------------------------------------------
  !  PH1, PC1, PSI1  -- exponents for an ad-hoc NLTE source-function
  !  correction applied as a multiplicative fudge to the LTE line source
  !  function.  Per depth point j, with ground-state Boltzmann populations
  !  of H I, C I, C II, Si I, and Si II (computed in PFSAHA and stored in
  !  bhyd_gs/bc1_gs/bc2_gs/bsi1_gs/bsi2_gs):
  !
  !      bfudge(j) = bhyd_gs(j)**PH1
  !                * (bc1_gs(j)/bc2_gs(j))**PC1
  !                * (bsi1_gs(j)/bsi2_gs(j))**PSI1
  !
  !      SLINE(j) = BNU(j) * STIM(j) / (bfudge(j) - EHVKT(j))
  !
  !  This is the photospheric analogue of using H, C, and Si departure
  !  coefficients as proxies for the dominant electron-donor non-LTE state
  !  in the line-forming layers.  All zeros means bfudge=1, recovering the
  !  pure-LTE source function.  This is the standard configuration; the
  !  exponents are nonzero only in tuning experiments.
  ! ==========================================================================
  REAL(8), PARAMETER :: LINE_SCAT_RHOX_SCALE = 0.0D0
  REAL(8), PARAMETER :: PH1                  = 0.0D0
  REAL(8), PARAMETER :: PC1                  = 0.0D0
  REAL(8), PARAMETER :: PSI1                 = 0.0D0

  ! ==========================================================================
  !  SYNTHE LOCAL ARRAYS
  ! ==========================================================================
  INTEGER(4) :: line_flag(maxlin_p)

  ! In-memory opacity matrices (replace unit 14 direct-access scratch file)
  REAL(4), ALLOCATABLE :: opacity_matrix(:,:)   ! (length, kw) wavelength × depth
  REAL(4), ALLOCATABLE :: linecen_matrix(:,:)   ! (n9, kw) line × depth

  ! In-memory merged line archive (replaces unit 13 scratch file)
  REAL(8), ALLOCATABLE :: merged_lindat8(:,:)   ! (14, nlines)
  REAL(4), ALLOCATABLE :: merged_lindat4(:,:)   ! (28, nlines)

  ! --- Run parameters (unit 93) ---
  INTEGER   :: nlines_in, ifvac, n19, linout
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
  REAL(8)   :: qdopple(mw6_p), qxnfpel(mw6_p)
  REAL(8)   :: qablog(1131)

  REAL(8)   :: ablog(3, 377)

  REAL(4)   :: asynth(kw_p), alinec(kw_p)
  INTEGER   :: mlinej(kw_p)

  REAL(4)   :: hfac(kw_p), hefac(kw_p), h2fac(kw_p)

  REAL(8)   :: lindat8(14)
  REAL(4)   :: lindat4(28)

  ! ==========================================================================
  !  SYNTHE SCALAR WORK VARIABLES
  ! ==========================================================================
  INTEGER   :: i, j, l, nu, iedge, nbuff_i, ios, eqpos
  INTEGER   :: n12, nlines, iline, ilines, n191
  INTEGER   :: n9, ncen
  INTEGER   :: nvshift
  INTEGER   :: maxred, maxblue, minblue_i, istep
  INTEGER   :: n10dop, nstep, n1, maxstep
  INTEGER   :: nbuff_s, congf_nel
  INTEGER   :: maxline, i9, nelem_i
  REAL(4)   :: kappa0_s, kapmin_s, kapcen_s
  REAL(4)   :: adamp_s
  REAL(4)   :: congf_s, alpha_s, v2_s
  REAL(4)   :: gamrf, gamsf, gamwf
  REAL(4)   :: dvoigt
  REAL(4)   :: x_wing
  REAL(8)   :: wave8, freq8
  REAL(4)   :: elo_s, dopple_nel

  ! ==========================================================================
  !  SPECTRV LOCAL VARIABLES
  ! ==========================================================================

  ! Command-line model filename and derived outputs
  CHARACTER(LEN=512) :: model_file, spec_file, linform_file, mol_file, model_base
  CHARACTER(LEN=64)  :: tmparg   ! scratch buffer for optional numeric arguments

  ! Wavelength loop variables (SPECTRV half)
  REAL(8)   :: wave_sv   ! current synthesis wavelength (nm)
  REAL(8)   :: wavold    ! previous line-centre wavelength (for bracket reset)
  REAL(8)   :: resid     ! residual intensity

  ! Integer loop/index variables (SPECTRV half)
  INTEGER   :: iedge_sv  ! continuum edge bracket for SPECTRV (separate from SYNTHE's iedge)

  ! Flags
  INTEGER   :: ifvac_sv   ! local copy of ifvac for SPECTRV /TRASH/ and title

  ! --- Variables formerly scoped in BLOCK constructs (hoisted to program level) ---
  CHARACTER(256) :: envval
  INTEGER        :: envlen, envstat
  INTEGER        :: dotpos
  INTEGER(8)     :: itmp
  INTEGER        :: jptr, k9
  REAL(8), ALLOCATABLE :: lindat8_arr(:,:)
  REAL(4), ALLOCATABLE :: lindat4_arr(:,:)

  ! ==========================================================================

  IFSYNTHE = 1
  
  ! ==========================================================================
  !  INITIALISE DATA DIRECTORY
  !
  !  Locate data files via $ATLAS12 environment variable.
  !  If unset, defaults to ./data/
  ! ==========================================================================
  CALL GET_ENVIRONMENT_VARIABLE('ATLAS12', envval, envlen, envstat)
  IF (envstat == 0 .AND. envlen > 0) THEN
    DATADIR = TRIM(envval)
    IF (DATADIR(envlen:envlen) /= '/') DATADIR = TRIM(DATADIR) // '/'
    DATADIR = TRIM(DATADIR) // 'data/'
  ELSE
    DATADIR = 'data/'
  END IF
  WRITE(6, '(/A,A)') ' DATADIR = ', TRIM(DATADIR)

  ! ==========================================================================
  !  PARSE COMMAND-LINE ARGUMENTS
  !
  !  Usage:  synthe_spectrv.exe <model_file> wlbeg=<nm> wlend=<nm> [resolu=<R>] [turbv=<kms>]
  !
  !  Argument 1  : model atmosphere filename (required, positional)
  !  wlbeg=<nm>  : start wavelength in nm (required)
  !  wlend=<nm>  : end wavelength in nm (required)
  !  resolu=<R>  : resolving power lambda/delta-lambda (optional, default 300000)
  !  turbv=<kms> : extra microturbulence in km/s (optional, default 0.0)
  !
  !  Keyword arguments may appear in any order after the model filename.
  !  The model filename extension is stripped and replaced with .spec / .linform
  !  for the output files.
  ! ==========================================================================
  IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
    WRITE(6, '(A)') ' ERROR: expected model filename as first argument'
    WRITE(6, '(A)') ' Usage: synthe_spectrv.exe <model_file> wlbeg=<nm> wlend=<nm> [resolu=<R>] [turbv=<kms>]'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, model_file)

  ! Set defaults for optional keyword args
  resolu = 300000.0D0
  turbv  = 0.0
  wlbeg  = 0.0D0    ! sentinel: must be set by keyword arg
  wlend  = 0.0D0    ! sentinel: must be set by keyword arg

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

  IF (wlbeg <= 0.0D0) THEN
    WRITE(6,'(A)') ' ERROR: wlbeg not specified or invalid'
    WRITE(6,'(A)') ' Usage: synthe_spectrv.exe <model_file> wlbeg=<nm> wlend=<nm> [resolu=<R>] [turbv=<kms>]'
    STOP 1
  END IF
  IF (wlend <= wlbeg) THEN
    WRITE(6,'(A)') ' ERROR: wlend not specified or <= wlbeg'
    WRITE(6,'(A)') ' Usage: synthe_spectrv.exe <model_file> wlbeg=<nm> wlend=<nm> [resolu=<R>] [turbv=<kms>]'
    STOP 1
  END IF

  WRITE(6, '(A,A)')      ' Input model      = ', TRIM(model_file)
  WRITE(6, '(A,F10.3)')  ' wlbeg (nm)       = ', wlbeg
  WRITE(6, '(A,F10.3)')  ' wlend (nm)       = ', wlend
  WRITE(6, '(A,F10.1)')  ' resolu           = ', resolu
  WRITE(6, '(A,F8.4)')   ' turbv (km/s)     = ', turbv

  ! Strip any leading directory from model_file for output filenames
  ! e.g. /path/to/sun.atm -> sun.atm
  dotpos = INDEX(TRIM(model_file), '/', BACK=.TRUE.)
  IF (dotpos > 0) THEN
    model_base = model_file(dotpos+1:)
  ELSE
    model_base = TRIM(model_file)
  END IF

  ! Strip extension: find last '.' and replace with .spec / .linform
  dotpos = INDEX(TRIM(model_base), '.', BACK=.TRUE.)
  IF (dotpos > 1) THEN
    spec_file    = model_base(1:dotpos-1) // '.spec'
    mol_file     = model_base(1:dotpos-1) // '.mol'
    linform_file = model_base(1:dotpos-1) // '.linform'
  ELSE
    spec_file    = TRIM(model_base) // '.spec'
    mol_file     = TRIM(model_base) // '.mol'
    linform_file = TRIM(model_base) // '.linform'
  END IF
  WRITE(6, '(A,A)') ' Spectrum output  = ', TRIM(spec_file)
  WRITE(6, '(A,A)') ' Linform output   = ', TRIM(linform_file)
  WRITE(6, '(A,A)') ' Mol nden output  = ', TRIM(mol_file)
  ! ==========================================================================
 ! OPEN(UNIT=93, FORM='UNFORMATTED')
 ! READ(93) nlines_in, length, ifvac, ifnlte, n19, turbv, deckj, ifpred, &
 !      wlbeg, wlend, resolu, ratio, ratiolg, cutoff, linout
 ! CLOSE(UNIT=93, STATUS='DELETE')

  ! deckj(1,j) = per-depth velocity shift (km/s); deckj(2,j) = magnetic field strength.
  ! Columns 3-7 are unused.  Default both to zero (no shift, no field).
  ! The fort.93 value is discarded here; when fort.93 is eventually removed,
  ! delete deckj from the READ above and this block becomes the only initialisation.
  deckj = 0.0

  ! ifvac=1: vacuum wavelengths throughout (sets 'V' flag in title, air_loc, ifvac_loc).
  ! ifnlte: read from fort.93 but never referenced anywhere downstream -- dead variable.
  ! ifpred: read from fort.93 but never referenced anywhere downstream -- dead variable.
  ! turbv, wlbeg, wlend, resolu: set from command-line args before the READ above; restored here
  !   because the READ overwrites them with whatever SYNBEG wrote.
  ! ratio, ratiolg: derived from resolu (also recomputed later in the SPECTRV section, but
  !   needed earlier at the ixwlbeg calculation below).
  ! length: derived as the number of log-wavelength steps from wlbeg to wlend at the given resolu.
  !   Verified: INT(log(wlend/wlbeg)/log(1+1/resolu)) = 636080 for 300-2500nm at R=300000.
  ! When fort.93 is eventually removed, delete all nine from the READ above and keep only these.
  ifvac  = 1
  resolu = 300000.0D0  ! re-apply defaults before keyword re-parse
  turbv  = 0.0
  DO i = 2, COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(i, tmparg)
    eqpos = INDEX(tmparg, '=')
    IF (eqpos < 2) CYCLE
    SELECT CASE (tmparg(1:eqpos-1))
      CASE ('wlbeg');  READ(tmparg(eqpos+1:), *) wlbeg
      CASE ('wlend');  READ(tmparg(eqpos+1:), *) wlend
      CASE ('resolu'); READ(tmparg(eqpos+1:), *) resolu
      CASE ('turbv');  READ(tmparg(eqpos+1:), *) turbv
    END SELECT
  END DO
  ratio   = 1.0D0 + 1.0D0 / resolu
  ratiolg = LOG(ratio)
  length  = INT(LOG(wlend / wlbeg) / ratiolg)

  ! cutoff=1e-3: line wing truncation threshold (kapmin = continuum * cutoff);
  !   lines whose wing opacity drops below this fraction of the local continuum
  !   are dropped.  Value confirmed from fort.93.
  ! linout=-30: negative value suppresses line identification output entirely
  !   (all linout>=0 branches are skipped).  Value confirmed from fort.93.
  ! When fort.93 is eventually removed, delete both from the READ above.
  cutoff = 1.0E-3
  linout = -30

  
  !nlines_in = 100
 
  
  ixwlbeg = INT(LOG(wlbeg) / ratiolg)
  wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  IF (wbegin < wlbeg) THEN
    ixwlbeg = ixwlbeg + 1
    wbegin  = EXP(DBLE(ixwlbeg) * ratiolg)
  END IF

  ! ==========================================================================
  !  SPECTRV INITIALISATION  (formerly the opening sections of spectrv.f90)
  !
  !  These steps are done once, before the depth/opacity loop.
  !  READIN populates the mod_atlas_data atmosphere state, then
  !  run_xnfpelsyn computes continuum opacities and ion populations.
  !
  !  This block runs *before* run_mklinelist so that teff_a is available
  !  to gate temperature-dependent line lists (H2O, TiO) out of hot-star
  !  runs.  The two stages are otherwise independent: line-list assembly
  !  does not touch atmosphere state, and the atmosphere initialisation
  !  does not touch the line lists.
  !
  !  The source-function fudge parameters (LINE_SCAT_RHOX_SCALE, PH1, PC1,
  !  PSI1) used to be read here from 'spectrv.input', but they're now
  !  hardcoded as PARAMETERs near the top of this file (see comments there).
  ! ==========================================================================

  ! --- Run READIN then run_xnfpelsyn ---
  ! IFPRES is set by READIN from the model cards (e.g. "PRESSURE OFF").
  itemp_a  = 1
  OPEN(UNIT=5,  FILE=TRIM(model_file),  STATUS='OLD', ACTION='READ')
  OPEN(UNIT=17, FILE=trim(DATADIR)//'continua.dat', STATUS='OLD', ACTION='READ')
  CALL readin(20)
  ! Keep unit 5 open: MOLEC reads molecular data from INPUTDATA (=5)
  ! on its first call during run_xnfpelsyn.
  CALL run_xnfpelsyn()
  CLOSE(UNIT=5)
  CLOSE(UNIT=17)

  ! ==========================================================================
  !  SECTION 2.  BUILD LINE LISTS IN MEMORY
  !
  !  run_mklinelist reads lines.list, dispatches to the appropriate readers
  !  (gfall, predict, mol, h2o) in canonical order, and populates the
  !  mod_mklinelist module arrays lte_lines(:) and nlte_lines(:).
  !  nlines_lte and nlines_nlte replace the former fort.12/fort.19 counts.
  !  fort.12 and fort.19 are never opened or written.
  !
  !  teff_a (set by readin above) gates H2O and TiO line lists: they are
  !  skipped when teff_a > TEFF_COOL_LIMIT (5000 K), even if listed in
  !  lines.list.
  ! ==========================================================================
  CALL run_mklinelist(wlbeg, wlend, resolu, teff_a, &
                      TRIM(DATADIR) // 'lines.list', DATADIR)

  nlines_in = nlines_lte
  n19       = nlines_nlte

  WRITE(6,'(A,I9)') ' nlines_in (LTE,  lte_lines)  = ', nlines_in
  WRITE(6,'(A,I9)') ' n19       (NLTE, nlte_lines) = ', n19


  OPEN(UNIT=35, FILE=TRIM(mol_file), STATUS='REPLACE', ACTION='WRITE')

  ! ==========================================================================
  !  SECTION 3.  MERGE LINE ARCHIVES INTO MEMORY
  !
  !  Line identification output (linout>=0) is permanently disabled (linout=-30);
  !  fort.14 and fort.20 are never opened or used.
  ! ==========================================================================

  ! ==========================================================================
  !  SECTION 4.  VOIGT FUNCTION
  !  The Weideman (1994) algorithm uses precomputed PARAMETER coefficients
  !  in synthe_module — no runtime initialisation is needed.
  ! ==========================================================================

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
  ! See the parameter declarations near the top of this file for the
  ! physical meaning of LINE_SCAT_RHOX_SCALE, PH1, PC1, PSI1.
  DO j = 1, nrhox
    bone(j)   = 1.0D0
    bfudge_sv(j) = bhyd_gs(j)**PH1 * (bc1_gs(j)/bc2_gs(j))**PC1 * &
                   (bsi1_gs(j)/bsi2_gs(j))**PSI1
    fscat_sv(j)  = 0.0D0
    IF (LINE_SCAT_RHOX_SCALE /= 0.0D0 .AND. &
        rhox_a(j)/LINE_SCAT_RHOX_SCALE < 100.0D0) &
      fscat_sv(j) = EXP(-rhox_a(j) / LINE_SCAT_RHOX_SCALE)
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

  ! Store 'A' or 'V' flag in title(74)
  itmp = INT(ICHAR('A'), 8)
  IF (ifvac == 1) itmp = INT(ICHAR('V'), 8)
  title(74) = TRANSFER(itmp, title(74))

  WRITE(6, '(/A,F6.1,A,I8,A,F6.1,A,I6,A,I2)') &
    ' WLBEG=', wlbeg, 'nm,   RESOLU=', int(resolu), ',   WLEND=', wlend, &
    'nm,   LENGTH=', length, ',   NRHOX=', nrhox

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

  iedge_sv = 1

  ! Initialise JOSH operator matrices (COEFJ, COEFH).
  ! The F77 code called BLOCKJH at the start of every JOSH call, which
  ! loaded the 51x51 Lambda-operator matrices from DATA statements.
  ! The F90 BLOCKJ/BLOCKH read from external files and cache the result,
  ! so they only need to be called once before the first JOSH call.
  CALL BLOCKJ
  CALL BLOCKH

  ! ==========================================================================
  !  SECTION 7.  MAIN DEPTH LOOP  (SYNTHE opacity accumulation)
  ! ==========================================================================
  ALLOCATE(opacity_matrix(length, nrhox))
  opacity_matrix = 0.0

  ilines   = 0
  n12      = nlines_in
  nlines   = nlines_in + n19
  WRITE(6,'(A,I9)') ' TOTAL LINES=', nlines

  depth_loop: DO j = 1, nrhox

    buffer(1:length) = 0.0

    ! Read continuum opacity (total = abs+scat) from module array instead of fort.10
    ! continall_m(nu,j) holds the same values as the former CONTINALL(nu,j) record.
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
      continuum(nbuff_i) = &
        ((wave8 - halfedge(iedge))*(wave8 - wledge(iedge+1))*ablog(1,iedge) + &
         (wledge(iedge) - wave8)*(wave8 - wledge(iedge+1))*2.0D0*ablog(2,iedge) + &
         (wave8 - wledge(iedge))*(wave8 - halfedge(iedge))*ablog(3,iedge)) / &
        deledge(iedge)
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
      if (qdopple(i) > 0.0D0) then
        xnfdop(i) = qxnfpel(i) / xf_rho(j) / qdopple(i)
      else
        xnfdop(i) = 0.0D0
      end if
    END DO
    
    txnxn(j) = (xnfh(j) + 0.42*xnfhe(j,1) + 0.85*xnfh2(j)) * &
               (t(j)/10000.0)**0.3

    nvshift = INT(resolu * DBLE(velshift_arr(j)) / CLIGHT_KMS + 0.5D0)


    mlines = 0
    !n19 is the number of NLTE lines
    IF (n19 > 0) CALL compute_line_opacity(j, n19, cutoff, velshift_arr(j), ifvac, linout)
    
    IF (n12 > 0) THEN
       
       n191    = n19 + 1
       alpha_s = 0.0
       
       DO iline = n191, nlines

          
          ! Access LTE line directly from in-memory array (replaces READ(12))
          nbuff_s  = lte_lines(iline - n19)%nbuff
          congf_s  = lte_lines(iline - n19)%cgf
          congf_nel= lte_lines(iline - n19)%nelion
          elo_s    = lte_lines(iline - n19)%elo
          gamrf    = lte_lines(iline - n19)%gamrf
          gamsf    = lte_lines(iline - n19)%gamsf
          gamwf    = lte_lines(iline - n19)%gamwf
          
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
             kapcen_s = kappa0_s * voigt_profile(0.0, adamp_s)
             IF (linout >= 0) CALL journal_append(iline, kapcen_s)
             buffer(nbuff_s) = buffer(nbuff_s) + kapcen_s
          END IF centre_on_grid
          
          dvoigt = 1.0 / dopple_nel / REAL(resolu)
          DO nstep = 1, n10dop
             profile(nstep) = kappa0_s * voigt_profile(REAL(nstep)*dvoigt, adamp_s)
             IF (profile(nstep) < kapmin_s) EXIT
          END DO
          
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
    
    mlinej(j) = mlines
    ilines    = ilines + mlines
    
  END DO depth_loop
  
  ! fort.19 eliminated -- NLTE lines now in nlte_lines() from mod_mklinelist

  ! ==========================================================================
  !  SECTION 8.  SPECTRV WAVELENGTH LOOP
  !
  !  For each wavelength point, assemble the opacity vector from the
  !  in-memory opacity_matrix, apply stimulated-emission correction,
  !  and call process_wavelength_point for the radiative transfer.
  ! ==========================================================================
  n9 = 0

  ! Open ASCII spectrum output (replaces standalone syntoascanga post-processing)
  OPEN(UNIT=11, FILE=TRIM(spec_file), STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=33, FILE=TRIM(linform_file), STATUS='REPLACE', ACTION='WRITE')

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

    DO i = 1, ilines
      line_flag(journal_iline(i)) = 1
    END DO

    ! First pass: count used lines and assign compacted indices
    DO i = 1, nlines
      IF (line_flag(i) == 0) CYCLE
      line_flag(i) = 0
      n9 = n9 + 1
      line_flag(i) = n9
    END DO

  END IF identify_lines

  nlines_sv = n9
  ! fort.12 eliminated -- LTE lines now in lte_lines() from mod_mklinelist

  IF (ilines == 0 .OR. linout < 0) THEN
    IF (ALLOCATED(merged_lindat8)) DEALLOCATE(merged_lindat8)
    IF (ALLOCATED(merged_lindat4)) DEALLOCATE(merged_lindat4)
    IF (ALLOCATED(journal_iline))  DEALLOCATE(journal_iline)
    IF (ALLOCATED(journal_kapcen)) DEALLOCATE(journal_kapcen)
    CLOSE(UNIT=11)
    CLOSE(UNIT=33)
    STOP
  END IF

  ! Second pass: store compacted lindat records in memory
  ALLOCATE(lindat8_arr(14, n9), lindat4_arr(28, n9))
  k9 = 0
  DO i = 1, nlines
    IF (line_flag(i) == 0) CYCLE
    k9 = k9 + 1
    lindat8_arr(:, k9) = merged_lindat8(:, i)
    lindat4_arr(:, k9) = merged_lindat4(:, i)
  END DO
  DEALLOCATE(merged_lindat8, merged_lindat4)

  ! ==========================================================================
  !  SECTION 10.  LINE-CENTRE OPACITIES AND IDENTIFICATION RECORDS
  !
  !  Build in-memory line-centre opacity matrix from journal_iline/journal_kapcen
  !  pairs, then loop over lines calling process_linecen_record for each.
  ! ==========================================================================
  ALLOCATE(linecen_matrix(n9, nrhox))
  linecen_matrix = 0.0

  jptr = 0
  DO j = 1, nrhox
    maxline = mlinej(j)
    IF (maxline == 0) CYCLE
    DO l = 1, maxline
      jptr = jptr + 1
      i9 = line_flag(journal_iline(jptr))
      IF (i9 == 0) CYCLE
      linecen_matrix(i9, j) = journal_kapcen(jptr)
    END DO
  END DO

  DEALLOCATE(journal_iline, journal_kapcen)
  journal_count = 0
  journal_size  = 0

  IF (n9 > maxlin_p) THEN
    WRITE(6,'(A)') ' TOO MANY LINES TO TRANSPOSE'
    STOP
  END IF

  ncen    = 0


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

  WRITE(6,'(I10,A)') ncen, ' LINE CENTER RECORDS PROCESSED'


  CLOSE(UNIT=11)
  CLOSE(UNIT=33)
  CLOSE(UNIT=35)

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

    REAL(8) :: wave_nm, q_loc(41)
    INTEGER :: jj, mu_loc

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
    INTEGER :: jj

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


END PROGRAM SYNTHE
